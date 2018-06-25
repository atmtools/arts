
/* Copyright (C) 2000-2012
   Stefan Buehler   <sbuehler@ltu.se>
   Patrick Eriksson <patrick.eriksson@chalmers.se>
   Axel von Engeln  <engeln@uni-bremen.de>
   Thomas Kuhn      <tkuhn@uni-bremen.de>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

// 

/**
   \file   m_abs.cc

   Stuff related to the calculation of absorption coefficients.

   \author Stefan Buehler
   \date   2001-03-12
*/
#include <cmath>
#include <algorithm>
#include "arts.h"
#include "matpackI.h"
#include "array.h"
#include "messages.h"
#include "file.h"
#include "auto_md.h"
#include "math_funcs.h"
#include "global_data.h"
#include "physics_funcs.h"
#include "absorption.h"
#include "continua.h"
#include "check_input.h"
#include "montecarlo.h"
#include "m_xml.h"
#include "optproperties.h"
#include "parameters.h"
#include "rte.h"
#include "xml_io.h"
#include "jacobian.h"

#ifdef ENABLE_NETCDF
#include <netcdf.h>
#include "nc_io.h"
#endif


extern const Numeric ELECTRON_CHARGE;
extern const Numeric ELECTRON_MASS;
extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric VACUUM_PERMITTIVITY;





/* Workspace method: Doxygen documentation will be auto-generated */
void AbsInputFromRteScalars(// WS Output:
                            Vector&        abs_p,
                            Vector&        abs_t,
                            Matrix&        abs_nlte,
                            Matrix&        abs_vmrs,
                            // WS Input:
                            const Numeric& rtp_pressure,
                            const Numeric& rtp_temperature,
                            const Vector&  rtp_nlte,
                            const Vector&  rtp_vmr,
                            const Verbosity&)
{
  // Prepare abs_p:
  abs_p.resize(1);
  abs_p = rtp_pressure;

  // Prepare abs_t:
  abs_t.resize(1);
  abs_t = rtp_temperature;
  
  // Prepare abs_nlte:
  abs_nlte.resize(rtp_nlte.nelem(),1);
  abs_nlte = rtp_nlte;

  // Prepare abs_vmrs:
  abs_vmrs.resize(rtp_vmr.nelem(),1);
  abs_vmrs = rtp_vmr;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetEmpty(// WS Output:
                                   ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                   // WS Input:
                                   const ArrayOfArrayOfSpeciesTag& tgs,
                                   const Verbosity&)
{
  // Make abs_lines_per_species the right size:
  abs_lines_per_species.resize( tgs.nelem() );
  
  for (Index i=0; i<tgs.nelem(); ++i)
    {
      abs_lines_per_species[i].resize(0);
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesArtscat4FromArtscat3(// WS Output:
                                   ArrayOfLineRecord& abs_lines,
                                   // Verbosity object:
                                   const Verbosity&)
{
    String fail_msg;
    bool failed = false;

    // Loop over all lines, use member function to do conversion.
#pragma omp parallel for      \
  if (!arts_omp_in_parallel()  \
      && abs_lines.nelem() >= arts_omp_get_max_threads())
    for ( Index i=0; i<abs_lines.nelem(); ++i )
    {
        try
        {
            abs_lines[i].ARTSCAT4FromARTSCAT3();
        }
        catch (const std::runtime_error &e)
        {
#pragma omp critical (abs_linesArtscat4FromArtscat3_fail)
            { fail_msg = e.what(); failed = true; }
        }
    }

    if (failed) throw runtime_error(fail_msg);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesArtscat5FromArtscat34(// WS Output:
                                    ArrayOfLineRecord& abs_lines,
                                    // Verbosity object:
                                    const Verbosity&)
{
    String fail_msg;
    bool failed = false;

    // Loop over all lines, use member function to do conversion.
#pragma omp parallel for      \
  if (!arts_omp_in_parallel()  \
      && abs_lines.nelem() >= arts_omp_get_max_threads())
    for ( Index i=0; i<abs_lines.nelem(); ++i )
    {
        try
        {
            if (abs_lines[i].Version() == 3)
                abs_lines[i].ARTSCAT5FromARTSCAT3();
            else if (abs_lines[i].Version() == 4)
                abs_lines[i].ARTSCAT5FromARTSCAT4();
            
        }
        catch (const std::runtime_error &e)
        {
#pragma omp critical (abs_linesArtscat4FromArtscat3_fail)
            { fail_msg = e.what(); failed = true; }
        }
    }

    if (failed) throw runtime_error(fail_msg);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReadFromLBLRTM(
                             // WS Output:
                             ArrayOfLineRecord& abs_lines,
                             // Control Parameters:
                             const String& filename,
                             const Numeric& fmin,
                             const Numeric& fmax,
                             const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  ifstream is;

  // Reset lines in case it already existed:
  abs_lines.resize(0);

  out2 << "  Reading file: " << filename << "\n";
  open_input_file(is, filename);

  bool go_on = true;
  while ( go_on )
    {
      LineRecord lr;
      if ( lr.ReadFromLBLRTMStream(is, verbosity) )
        {
          // If we are here the read function has reached eof and has
          // returned no data.
          go_on = false;
        }
      else
        {
          if ( fmin <= lr.F() )
            {
              if ( lr.F() <= fmax )
                abs_lines.push_back(std::move(lr));
              else
                go_on = false;
            }
        }
    }
  out2 << "  Read " << abs_lines.nelem() << " lines.\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReadFromHitranPre2004(
                             // WS Output:
                             ArrayOfLineRecord& abs_lines,
                             // Control Parameters:
                             const String& filename,
                             const Numeric& fmin,
                             const Numeric& fmax,
                             const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  ifstream is;

  // Reset lines in case it already existed:
  abs_lines.resize(0);

  out2 << "  Reading file: " << filename << "\n";
  open_input_file(is, filename);

  bool go_on = true;
  while ( go_on )
    {
      LineRecord lr;
      if ( lr.ReadFromHitran2001Stream(is, verbosity) )
        {
          // If we are here the read function has reached eof and has
          // returned no data.
          go_on = false;
        }
      else
        {
          if ( fmin <= lr.F() )
            {
              if ( lr.F() <= fmax )
                abs_lines.push_back(std::move(lr));
              else
                go_on = false;
            }
        }
    }
  out2 << "  Read " << abs_lines.nelem() << " lines.\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReadFromHitran(// WS Output:
                             ArrayOfLineRecord& abs_lines,
                             // Control Parameters:
                             const String& filename,
                             const Numeric& fmin,
                             const Numeric& fmax,
                             const Verbosity& verbosity)
{
    CREATE_OUT2;
    
    ifstream is;
    
    // Reset lines in case it already existed:
    abs_lines.resize(0);
    
    String filename_lower = filename;
    filename_lower.tolower();
    ArrayOfString splitted_fname;
    filename_lower.split(splitted_fname, "/");
    if (splitted_fname.nelem())
    {
        filename_lower = splitted_fname[splitted_fname.nelem()-1];
    }
    else
    {
        throw std::runtime_error("Catalog filename is empty");
    }
#ifdef USE_HITRAN2008
    if (filename_lower.nelem() < 8
        || (filename_lower.substr(0, 8) != "hitran08"
            && filename_lower.substr(0, 8) != "hitran04"))
    {
        ostringstream os;
        os << "'" << filename << "'\n"
        << "does not appear to be a HITRAN 2008 catalogue. The catalog filename\n"
        << "name must start with HITRAN08. This version of arts was compiled with\n"
        << "support only for HITRAN 2008. To switch back to the latest HITRAN\n"
        << "run 'cmake -DWITH_HITRAN2008=0 ..' and recompile arts";
        throw std::runtime_error(os.str());
    }
#else
    if (filename_lower.nelem() < 10
        || filename_lower.substr(0, 10) != "hitran2012")
    {
        ostringstream os;
        os << "'" << filename << "'\n"
        << "does not appear to be a HITRAN 2012 catalogue. The catalog filename\n"
        << "must start with HITRAN2012. If you intend to use a HITRAN 2008 catalog\n"
        << "run 'cmake -DWITH_HITRAN2008=1 ..' and recompile arts";
        throw std::runtime_error(os.str());
    }
#endif

    out2 << "  Reading file: " << filename << "\n";
    open_input_file(is, filename);
    
    bool go_on = true;
    Index linenr = 0;
    while ( go_on )
    {
        linenr++;
        try {
            LineRecord lr;
            if ( lr.ReadFromHitran2004Stream(is, verbosity, fmin) )
            {
                // If we are here the read function has reached eof and has
                // returned no data.
                go_on = false;
            }
            else
            {
                if ( fmin <= lr.F() )
                {
                    if ( lr.F() <= fmax )
                        abs_lines.push_back(std::move(lr));
                    else
                        go_on = false;
                }
            }
        }
        catch (const std::runtime_error &e)
        {
            ostringstream os;
            os << e.what() << "\n";
            os << "Error parsing line " << linenr << " from catalog.\n";
            throw runtime_error(os.str());
        }
    }
    out2 << "  Read " << abs_lines.nelem() << " lines.\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReadFromMytran2(// WS Output:
                              ArrayOfLineRecord& abs_lines,
                              // Control Parameters:
                              const String& filename,
                              const Numeric& fmin,
                              const Numeric& fmax,
                              const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  ifstream is;

  // Reset lines in case it already existed:
  abs_lines.resize(0);

  out2 << "  Reading file: " << filename << "\n";
  open_input_file(is, filename);

  bool go_on = true;
  while ( go_on )
    {
      LineRecord lr;
      if ( lr.ReadFromMytran2Stream(is, verbosity) )
        {
          // If we are here the read function has reached eof and has
          // returned no data.
          go_on = false;
        }
      else
        {
          // lines are not necessarily frequency sorted 
          if ( fmin <= lr.F() )
            if ( lr.F() <= fmax )
              abs_lines.push_back(std::move(lr));
        }
    }
  out2 << "  Read " << abs_lines.nelem() << " lines.\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReadFromJpl(// WS Output:
                          ArrayOfLineRecord& abs_lines,
                          // Control Parameters:
                          const String& filename,
                          const Numeric& fmin,
                          const Numeric& fmax,
                          const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  ifstream is;

  // Reset lines in case it already existed:
  abs_lines.resize(0);

  out2 << "  Reading file: " << filename << "\n";
  open_input_file(is, filename);

  bool go_on = true;
  while ( go_on )
    {
      LineRecord lr;
      if ( lr.ReadFromJplStream(is, verbosity) )
        {
          // If we are here the read function has reached eof and has
          // returned no data.
          go_on = false;
        }
      else
        {
          // we expect lines to be sorted
          if ( fmin <= lr.F() )
            {
              if ( lr.F() <= fmax )
                abs_lines.push_back(std::move(lr));
              else
                go_on = false;
            }
        }
    }
  out2 << "  Read " << abs_lines.nelem() << " lines.\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReadFromArts(// WS Output:
                           ArrayOfLineRecord& abs_lines,
                           // Control Parameters:
                           const String& filename,
                           const Numeric& fmin,
                           const Numeric& fmax,
                           const Verbosity& verbosity)
{
  xml_read_arts_catalogue_from_file (filename, abs_lines, fmin, fmax, verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesWriteToSplitArtscat(// WS Input:
                                              const String& output_file_format,
                                              const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                              // Control Parameters:
                                              const String& basename,
                                              const Verbosity& verbosity)
{
  using global_data::species_data;

  for (ArrayOfArrayOfLineRecord::const_iterator it = abs_lines_per_species.begin();
       it != abs_lines_per_species.end();
       it++)
  {
    if (it->nelem())
    {
      String species_filename = basename;
      if (basename.length() && basename[basename.length()-1] != '/')
        species_filename += ".";
      
      species_filename += species_data[(*it)[0].Species()].Name() + ".xml";
      WriteXML(output_file_format, *it, species_filename, 0,
               "", "", "", verbosity);
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReadFromSplitArtscat(// WS Output:
                                   ArrayOfLineRecord& abs_lines,
                                   const ArrayOfArrayOfSpeciesTag& abs_species,
                                   // Control Parameters:
                                   const String& basename,
                                   const Numeric& fmin,
                                   const Numeric& fmax,
                                   const Verbosity& verbosity)
{
  CREATE_OUT2;
  using global_data::species_data;
 
  // Build a set of species indices. Duplicates are ignored.
  set<Index> unique_species;
  for (ArrayOfArrayOfSpeciesTag::const_iterator asp = abs_species.begin();
       asp != abs_species.end(); asp++)
    for (ArrayOfSpeciesTag::const_iterator sp = asp->begin();
         sp != asp->end(); sp++)
    {
      // Of the four different types of SpeciesTags, only PLAIN and ZEEMAN
      // actually require explicit spectral lines.
      if (sp->Type()==SpeciesTag::TYPE_PLAIN ||
          sp->Type()==SpeciesTag::TYPE_ZEEMAN) {
          unique_species.insert(sp->Species());
      }
      
//      if (sp->Isotopologue() == species_data[sp->Species()].Isotopologue().nelem()
//          || !species_data[sp->Species()].Isotopologue()[sp->Isotopologue()].isContinuum())
//          {
//            unique_species.insert(sp->Species());
//          }
    }
  
  // Read catalogs for each identified species and put them all into
  // abs_lines.
  abs_lines.resize(0);
  for (set<Index>::const_iterator it = unique_species.begin();
       it != unique_species.end(); it++)
  {
    ArrayOfLineRecord more_abs_lines;
    String tmpbasename = basename;
    if (basename.length() && basename[basename.length()-1] != '/')
    {
       tmpbasename += '.';
    }

    abs_linesReadFromArts(more_abs_lines, tmpbasename + (species_data[*it].Name()) + ".xml",
                          fmin, fmax, verbosity);
    abs_lines.insert(abs_lines.end(), more_abs_lines.begin(), more_abs_lines.end());
  }
  
  out2 << "  Read " << abs_lines.nelem() << " lines in total.\n";
}


//! Obsolete old ARTS catalogue reading function
/*!
  This function is for the old ARTS catalogue format without XML
  header. It is no longer used, but left here for historical reasons. 
  
  \author Stefan Buehler

  \param abs_lines Line data.
  \param filename Name of catalogue file.
  \param fmin Minimum frequency.
  \param fmax Maximum frequency. */
void abs_linesReadFromArtsObsolete(// WS Output:
                                   ArrayOfLineRecord& abs_lines,
                                   // Control Parameters:
                                   const String& filename,
                                   const Numeric& fmin,
                                   const Numeric& fmax,
                                   const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  // The input stream:
  ifstream is;

  // We will use this line record to read in the line records in the
  // file one after another:
  LineRecord lr;

  // Reset lines in case it already existed:
  abs_lines.resize(0);

  out2 << "  Reading file: " << filename << "\n";
  open_input_file(is, filename);

  // Get version tag and check that it corresponds to the current version.
  {
    String v;
    is >> v;
    if ( v!=lr.VersionString() )
      {
        ostringstream os;
        
        // If what we read is the version String, it should have at elast 9 characters.
        if ( 9 <= v.nelem() )
          {
            if ( "ARTSCAT" == v.substr(0,7) )
            {
              os << "The ARTS line file you are trying contains a version tag "
                 << "different from the current version.\n"
                 << "Tag in file:     " << v << "\n"
                 << "Current version: " << lr.Version();
              throw runtime_error(os.str());
            }
           }

            os << "The ARTS line file you are trying to read does not contain a valid version tag.\n"
            << "Probably it was created with an older version of ARTS that used different units.";
            throw runtime_error(os.str());
      }
  }

  bool go_on = true;
  while ( go_on )
    {
      if ( lr.ReadFromArtscat3Stream(is, verbosity) )
        {
          // If we are here the read function has reached eof and has
          // returned no data.
          go_on = false;
        }
      else
        {
          // lines are not necessarily frequency sorted 
          if ( fmin <= lr.F() && lr.F() <= fmax )
            {
              abs_lines.push_back(std::move(lr));
            }
        }
    }
  out2 << "  Read " << abs_lines.nelem() << " lines.\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesReadFromCatalogues(// WS Output:
                                             ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                             // WS Input:
                                             const ArrayOfArrayOfSpeciesTag& tgs,
                                             // Control Parameters:
                                             const ArrayOfString& filenames,
                                             const ArrayOfString& formats,
                                             const Vector& fmin,
                                             const Vector& fmax,
                                             const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  const Index n_tg   = tgs.nelem();     // # tag groups
  const Index n_cat = filenames.nelem();        // # tag Catalogues

  // Check that dimensions of the keyword parameters are consistent
  // (must all be the same). 

  if ( n_cat != formats.nelem() ||
       n_cat != fmin.nelem() ||
       n_cat != fmax.nelem() )
    {
      ostringstream os;
      os << "abs_lines_per_speciesReadFromCatalogues: All keyword\n"
         << "parameters must get the same number of arguments.\n"
         << "You specified:\n"
         << "filenames: " << n_cat         << "\n"
         << "formats:   " << formats.nelem() << "\n"
         << "fmin:      " << fmin.nelem()    << "\n"
         << "fmax:      " << fmax.nelem();
      throw runtime_error(os.str());
    }
  
  // Furthermore, the dimension must be
  // smaller than or equal to the number of tag groups.

  if ( n_cat > n_tg )
    {
      ostringstream os;
      os << "abs_lines_per_speciesReadFromCatalogues: You specified more\n"
         << "catalugues than you have tag groups.\n"
         << "You specified:\n"
         << "Catalogues: " << n_cat << "\n"
         << "tgs: " << n_tg;
      throw runtime_error(os.str());
    }

  // There must be at least one tag group and at least one catalogue:

  if ( n_cat < 1 ||
       n_tg   < 1 )
    {
      ostringstream os;
      os << "abs_lines_per_speciesReadFromCatalogues: You must have at\n"
         << "least one catalogue and at least one tag group.\n"
         << "You specified:\n"
         << "Catalogues: " << n_cat << "\n"
         << "tgs: " << n_tg;
      throw runtime_error(os.str());
    }

  // There can be repetitions in the keyword paramters. We want to read
  // and process each catalogue only once, so we'll compile a set of
  // real catalogues, along with a data structure that tells us which
  // tag groups should use this catalogue.

  ArrayOfString real_filenames{filenames[0]};
  ArrayOfString real_formats{formats[0]};
  ArrayOfNumeric real_fmin{fmin[0]};
  ArrayOfNumeric real_fmax{fmax[0]};

  Array< ArrayOfIndex > real_tgs(1);
  real_tgs[0].resize(1);
  real_tgs[0][0] = 0;

  // The last specified catalogue, to which we should assign all
  // remaining lines. Index of this one in real_ arrays.
  Index last_cat = 0;           

  for ( Index i=1; i<n_tg; ++i )
    {
      // Is there a catalogue specified?
      if ( n_cat > i )
        {
          // Yes, there is a catalogue.

          // Has this been specified before?
          // We use the STL find algorithm to look for the catalogue
          // name in the real_catalogues. Find returns an iterator, so
          // to get an index we have to take the difference to
          // .begin(). 
          const Index that_cat = find( real_filenames.begin(),
                                        real_filenames.end(),
                                        filenames[i] ) - real_filenames.begin();
          if ( that_cat < real_filenames.nelem() )
            {
              // Yes, it has been specified before
              // ==> Assign to that catalogue
              real_tgs[that_cat].push_back(i);

              // Verify, that format, fmin, and fmax are consistent:
              if ( formats[i] != real_formats[that_cat] ||
                   fmin[i]    != real_fmin[that_cat]    ||
                   fmax[i]    != real_fmax[that_cat]       )
                {
                  ostringstream os;
                  os << "abs_lines_per_speciesReadFromCatalogues: If you specify the\n"
                     << "same catalogue repeatedly, format, fmin, and fmax must be\n"
                     << "consistent. There is an inconsistency between\n"
                     << "catalogue " << that_cat << " and " << i << ".";
                  throw runtime_error(os.str());
                }
            }
          else
            {
              // No, it has not been specified before.
              // ==> Add an entry to real_tgs and the other real_ variables:
              real_tgs.push_back({i});

              real_filenames.push_back( filenames[i] );
              real_formats.push_back  ( formats[i]   );  
              real_fmin.push_back     ( fmin[i]      );     
              real_fmax.push_back     ( fmax[i]      );

              last_cat = i;     // assign remainder of lines to this
                                // catalogue, if there is no other catalogue.
            }
        }
      else
        {
          // No, there is no catalogue.
          // ==> Assign to the last catalogue
          real_tgs[last_cat].push_back(i);
        }
    }

  Index n_real_cat = real_filenames.nelem(); // # real catalogues to read

  // Some output to low priority stream:
  out3 << "  Catalogues to read and tgs for which these will be used:\n";
  for ( Index i=0; i<n_real_cat; ++i )
    {
      out3 << "  " << real_filenames[i] << ":";
      for ( Index s=0; s<real_tgs[i].nelem(); ++s )
        out3 << " " << real_tgs[i][s];
      out3 << "\n";
    }

  // Make abs_lines_per_species the right size:
  abs_lines_per_species.resize( tgs.nelem() );

  // Loop through the catalogues to read:
  for ( Index i=0; i<n_real_cat; ++i )
    {
      ArrayOfLineRecord   abs_lines;

      // Read catalogue:

      if ( "HITRAN96"==real_formats[i] )
        {
          abs_linesReadFromHitranPre2004( abs_lines, real_filenames[i], real_fmin[i], real_fmax[i], verbosity );
        }
      else if ( "HITRAN04"==real_formats[i] )
        {
          abs_linesReadFromHitran( abs_lines, real_filenames[i], real_fmin[i], real_fmax[i], verbosity );
        }
      else if ( "MYTRAN2"==real_formats[i] )
        {
          abs_linesReadFromMytran2( abs_lines, real_filenames[i], real_fmin[i], real_fmax[i], verbosity );
        }
      else if ( "JPL"==real_formats[i] )
        {
          abs_linesReadFromJpl( abs_lines, real_filenames[i], real_fmin[i], real_fmax[i], verbosity );
        }
      else if ( "ARTS"==real_formats[i] )
        {
          abs_linesReadFromArts( abs_lines, real_filenames[i], real_fmin[i], real_fmax[i], verbosity );
        }
      else
        {
          ostringstream os;
          os << "abs_lines_per_speciesReadFromCatalogues: You specified the\n"
             << "format `" << real_formats[i] << "', which is unknown.\n"
             << "Allowd formats are: HITRAN96, HITRAN04, MYTRAN2, JPL, ARTS.";
          throw runtime_error(os.str());
        }

      // We need to make subset tgs for the groups that should
      // be read from this catalogue.
      ArrayOfArrayOfSpeciesTag  these_tgs(real_tgs[i].nelem());
      for ( Index s=0; s<real_tgs[i].nelem(); ++s )
        {
          these_tgs[s] = tgs[real_tgs[i][s]];
        }

      // Create these_abs_lines_per_species:
      ArrayOfArrayOfLineRecord these_abs_lines_per_species;
      abs_lines_per_speciesCreateFromLines( these_abs_lines_per_species, abs_lines, these_tgs, verbosity );

      // Put these lines in the right place in abs_lines_per_species:
      for ( Index s=0; s<real_tgs[i].nelem(); ++s )
        {
          abs_lines_per_species[real_tgs[i][s]] = std::move(these_abs_lines_per_species[s]);
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesCreateFromLines(// WS Output:
                                          ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                          // WS Input:
                                          const ArrayOfLineRecord& abs_lines,
                                          const ArrayOfArrayOfSpeciesTag& tgs,
                                          const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // The species lookup data:
  using global_data::species_data;

  // As a safety feature, we will watch out for the case that a
  // species is included in the calculation, but not all lines are
  // used. For this we need an array to flag the used species:
  
  // For some weird reason, Arrays of bool do not work, although all
  // other types seem to work fine. So in this single case, I'll use
  // the stl vector directly. The other place where this is done is in
  // the function executor in main.cc.
  // FIXME: Fix this when Array<bool> works.
  std::vector<bool> species_used (species_data.nelem(),false);
      
  // Make abs_lines_per_species the right size:
  abs_lines_per_species = ArrayOfArrayOfLineRecord(tgs.nelem());

  // Unfortunately, MTL conatains a bug that leads to all elements of
  // the outer Array of an Array<Array>> pointing to the same data
  // after creation. So we need to fix this explicitly:
  for ( Index i=0; i<abs_lines_per_species.nelem(); ++i )
    abs_lines_per_species[i] = ArrayOfLineRecord();

  // Loop all lines in the input line list:
  for ( Index i=0; i<abs_lines.nelem(); ++i )
    {
      // Get a convenient reference to the current line:
      const LineRecord& this_line = abs_lines[i];

      // We have to test each tag group in turn (in the order in which
      // they appear in the controlfile). The line is assigned to the
      // first tag group that fits.

      // The flag found is used to break the for loops when the right
      bool found = false;

      // We need to define j here, since we need the value outside the
      // for loop:
      Index j;

      // Loop the tag groups:
      for ( j=0; j<tgs.nelem() && !found ; ++j ) 
        {
          // A tag group can contain several tags:
          for ( Index k=0; k<tgs[j].nelem() && !found; ++k )
            {
              // Get a reference to the current tag (not really
              // necessary, but makes for nicer notation below):
              const SpeciesTag& this_tag = tgs[j][k];

              // Now we will test different attributes of the line
              // against matching attributes of the tag. If any
              // attribute does not match, we continue with the next tag
              // in the tag group. (Exception: Species, see just below.)

              // Test species. If this attribute does not match we don`t
              // have to test the other tags in this group, since all
              // tags must belong to the same species.
              if ( this_tag.Species() != this_line.Species() ) break;

              // Test isotopologue. The isotopologue can either match directly, or
              // the Isotopologue of the tag can be one larger than the
              // number of isotopologues, which means `all'. Test the second
              // condition first, since this will probably be more often
              // used.
              if ( this_tag.Isotopologue() != this_line.SpeciesData().Isotopologue().nelem() )
                if ( this_tag.Isotopologue() != this_line.Isotopologue() )
                  continue;

            // Test frequncy range we take both the lower (Lf) and the
            // upper (Uf) border frequency to include the `equal' case.
            // Both Lf and Uf can also be negative, which means `no limit'

            // Take the lower limit first:
              if ( this_tag.Lf() >= 0 )
                if ( this_tag.Lf() > this_line.F() )
                  continue;

            // Then the upper limit:
              if ( this_tag.Uf() >= 0 )
                if ( this_tag.Uf() < this_line.F() )
                  continue;

            // When we get here, this_tag has survived all tests. That
            // means it matches the line perfectly!
              found = true;
            }
        }

      // If a matching tag was found, this line can be used in the
      // calculation. Add it to the line list for this tag group.
      if (found)
        {
          // We have to use j-1 here, since j was still increased by
          // one after the matching tag has been found.
          abs_lines_per_species[j-1].push_back(this_line);

          // Flag this species as used, if not already done:
          if ( !species_used[this_line.Species()] )
            species_used[this_line.Species()] = true;
        }
      else
        {
          // Safety feature: Issue a warning messages if the lines for a
          // species are only partly covered by tags.
          if ( species_used[this_line.Species()] )
            {
              CREATE_OUT2;
              out2 << "  Your tags include other lines of species "
                   << this_line.SpeciesData().Name()
                   << ",\n"
                   << "  Why do you not include line "
                   << i
                   << " (at "
                   << this_line.F()
                   << " Hz)?\n";
            }
        }
   }

  // Write some information to the lowest priority output stream.
  for (Index i=0; i<tgs.nelem(); ++i)
    {
        out3 << "  " << i << ":";

        for (Index s=0; s<tgs[i].nelem(); ++s)
          out3 << " " << tgs[i][s].Name();

        out3 << ": " << abs_lines_per_species[i].nelem() << " lines\n";
    }

}

void abs_lines_per_speciesAddMirrorLines(
         ArrayOfArrayOfLineRecord& abs_lines_per_species,
   const Numeric& max_f,
   const Verbosity&)
{
  // We will simply append the mirror lines after the original
  // lines. This way we don't have to make a backup copy of the
  // original lines. 

  for ( Index i=0; i<abs_lines_per_species.nelem(); ++i )
    {
      // Get a reference to the current list of lines to save typing:
      ArrayOfLineRecord& ll = abs_lines_per_species[i];
      
      // It is important that we determine the size of ll *before*
      // we start the loop. After all, we are adding elements. And
      // we cerainly don't want to continue looping the newly added
      // elements, we want to loop only the original elements.
      //
      const Index n = ll.nelem();

      // For increased efficiency, reserve the necessary space:
      //
      Index nnew = n;
      //
      if( max_f >= 0 )
        { 
          nnew = 0;
          for ( Index j=0; j<n; ++j )
            {
              if( ll[j].F() <= max_f )
                { nnew += 1; }
            }
        }
      //
      ll.reserve( n + nnew );
      
      // Loop through all lines of this tag group:
      for ( Index j=0; j<n; ++j )
        {          
          if( max_f < 0  ||  ll[j].F() <= max_f )
            {
              LineRecord dummy = ll[j];
              dummy.setF( -dummy.F() );
              ll.push_back(std::move(dummy));
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesCompact(// WS Output:
                                  ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                  // WS Input:
                                  const ArrayOfLineshapeSpec& abs_lineshape,
                                  const Vector& f_grid,
                                  const Verbosity&)
{
  // Make sure abs_lines_per_species and abs_lineshape have the same
  // dimension. As a special case, abs_lineshape can have length 1, in
  // which case it is valid for all species.
  if ( (abs_lines_per_species.nelem() != abs_lineshape.nelem()) &&
       (abs_lineshape.nelem()         != 1) ) 
    {
      ostringstream os;
      os << "Dimension of abs_lines_per_species does\n"
         << "not match that of abs_lineshape.\n"
         << "(Dimension of abs_lineshape must be 1 or\n"
         << "equal to abs_lines_per_species.)";
      throw runtime_error(os.str());
    }
  
  // Make sure that the frequency grid is properly sorted:
  for ( Index s=0; s<f_grid.nelem()-1; ++s )
    {
      if ( f_grid[s+1] <= f_grid[s] )
        {
          ostringstream os;
          os << "The frequency grid f_grid is not properly sorted.\n"
             << "f_grid[" << s << "] = " << f_grid[s] << "\n"
             << "f_grid[" << s+1 << "] = " << f_grid[s+1];
          throw runtime_error(os.str());
        }
    }

  // Cycle through all tag groups:
  for ( Index i=0; i<abs_lines_per_species.nelem(); ++i )
    {
      // Get cutoff frequency of this tag group.
      // We have to also handle the special case that there is only
      // one lineshape for all species.
      Numeric cutoff;
      if (1==abs_lineshape.nelem())
        cutoff = abs_lineshape[0].Cutoff();
      else
        cutoff = abs_lineshape[i].Cutoff();

      // Check whether cutoff is defined:
      if ( cutoff != -1)
        {
          // Get a reference to the current list of lines to save typing:
          ArrayOfLineRecord& ll = abs_lines_per_species[i];

          // Calculate the borders:
          Numeric upp = f_grid[f_grid.nelem()-1] + cutoff;
          Numeric low = f_grid[0] - cutoff;

          // Cycle through all lines within this tag group. 
          Array<ArrayOfLineRecord::iterator> keep;
          for ( ArrayOfLineRecord::iterator j=ll.begin(); j<ll.end(); ++j )
            {
              // Center frequency:
              const Numeric F0 = j->F();

              if ( ( F0 >= low) && ( F0 <= upp) )
                {
                  // Build list of elements which should be kept
                  keep.push_back (j);
                  // The original implementation just erased the non-wanted
                  // elements from the array. Problem: On every erase the
                  // whole array will be copied which actually kills
                  // performance.
                }
            }

          // If next comparison is false, some elements have to be removed
          if (keep.nelem () != ll.nelem ())
            {
              ArrayOfLineRecord nll;
              // Copy all elements that should be kept to a new Array
              for (Array<ArrayOfLineRecord::iterator>::iterator j
                   = keep.begin(); j != keep.end(); j++)
                {
                  nll.push_back (std::move(**j));
                }
              // Overwrite the old array with the new one
              ll.resize (nll.nelem ());
              ll = std::move(nll);
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesDefineAllInScenario(// WS Output:
                                    ArrayOfArrayOfSpeciesTag& tgs,
                                    Index& propmat_clearsky_agenda_checked,
                                    Index& abs_xsec_agenda_checked,
                                    // Control Parameters:
                                    const String& basename,
                                    const Verbosity& verbosity)
{
  CREATE_OUT2;

  // Invalidate agenda check flags
  propmat_clearsky_agenda_checked = false;
  abs_xsec_agenda_checked = false;

  // Species lookup data:
  using global_data::species_data;

  // We want to make lists of included and excluded species:
  ArrayOfString included(0), excluded(0);

  tgs.resize(0);

  for ( Index i=0; i<species_data.nelem(); ++i )
    {
      const String specname = species_data[i].Name();
      
      String filename = basename;
      if (basename.length() && basename[basename.length()-1] != '/')
        filename += ".";
      filename += specname;

      try {
          find_xml_file(filename, verbosity);
          // Add to included list:
          included.push_back(specname);

          // Convert name of species to a SpeciesTag object:
          SpeciesTag this_tag(specname);

          // Create Array of SpeciesTags with length 1
          // (our tag group has only one tag):
          ArrayOfSpeciesTag this_group(1);
          this_group[0] = this_tag;

          // Add this tag group to tgs:
          tgs.push_back(this_group);
      }
      catch (const std::runtime_error &e)
      {
          // The file for the species could not be found.
          excluded.push_back(specname);
      }
    }
  
  // Some nice output:
  out2 << "  Included Species (" << included.nelem() << "):\n";
  for ( Index i=0; i<included.nelem(); ++i )
    out2 << "     " << included[i] << "\n";

  out2 << "  Excluded Species (" << excluded.nelem() << "):\n";
  for ( Index i=0; i<excluded.nelem(); ++i )
    out2 << "     " << excluded[i] << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lineshapeDefine(// WS Output:
                         ArrayOfLineshapeSpec&    abs_lineshape,
                         // WS Input:
                         const String&            shape,
                         const String&            normalizationfactor,
                         const Numeric&           cutoff,
                         const Verbosity&         verbosity)
{
  CREATE_OUT2;
  
  // Make lineshape and normalization factor data visible:
  using global_data::lineshape_data;
  using global_data::lineshape_norm_data;


  // generate the right number of elements
  //  Index tag_sz = tgs.nelem();
  // We generate only 1 copy of the lineshape settings. Absorption
  // routines check for this case and use it for all species.
  Index tag_sz = 1;
  abs_lineshape.resize(tag_sz);

  // Is this lineshape available?
  Index found0=-1;
  for ( Index i=0; i<lineshape_data.nelem() && (found0 == -1) ; ++i )
    {
      const String& str = lineshape_data[i].Name();
      if (str == shape) 
        {
          out2 << "  Selected lineshape: " << str << "\n";
          found0=i;
        }
    }

  // Is this normalization to the lineshape available?
  Index found1=-1;
  for ( Index i=0; i<lineshape_norm_data.nelem() && (found1 == -1); ++i )
    {
      const String& str = lineshape_norm_data[i].Name();
      if (str == normalizationfactor) 
        {
          out2 << "  Selected normalization factor  : " << normalizationfactor << "\n";

          if ( (cutoff != -1) && (cutoff < 0.0) )
            throw runtime_error("  Cutoff must be -1 or positive.");
          out2 << "  Selected cutoff frequency [Hz] : " << cutoff << "\n";
          found1=i;
        }
    }


  // did we find the lineshape and normalization factor?
  if (found0 == -1)
    throw runtime_error("Selected lineshape not available.");
  if (found1 == -1)
    throw runtime_error("Selected normalization to lineshape not available.");


  // now set the lineshape  
  for (Index i=0; i<tag_sz; i++)
    {
      abs_lineshape[i].SetInd_ls( found0 );
      abs_lineshape[i].SetInd_lsn( found1 );
      abs_lineshape[i].SetCutoff( cutoff );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lineshape_per_tgDefine(// WS Output:
                                ArrayOfLineshapeSpec& abs_lineshape,
                                // WS Input:
                                const ArrayOfArrayOfSpeciesTag&      tgs,
                                const ArrayOfString&  shape,
                                const ArrayOfString&  normalizationfactor,
                                const Vector&         cutoff,
                                const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  // Make lineshape and normalization factor data visible:
  using global_data::lineshape_data;
  using global_data::lineshape_norm_data;

  // check that the number of elements are equal
  Index tg_sz = tgs.nelem();
  if ( (tg_sz != shape.nelem()) ||
       (tg_sz != normalizationfactor.nelem()) || 
       (tg_sz != cutoff.nelem()) )
    {
      ostringstream os;
      os << "abs_lineshape_per_tgDefine: number of elements does\n"
         << "not match the number of tag groups defined.";
      throw runtime_error(os.str());
    }
      

  // generate the right number of elements
  abs_lineshape.resize(tg_sz);

  // Is this lineshape available?
  for (Index k=0; k<tg_sz; ++k)
    {
      Index found0=-1;
      for ( Index i=0; i<lineshape_data.nelem() && (found0 == -1); ++i )
        {
          const String& str = lineshape_data[i].Name();
          if (str == shape[k]) 
            {
              out2 << "  Tag Group: [";
              for (Index s=0; s<tgs[k].nelem()-1; ++s)
                out2 << tgs[k][s].Name() << ", "; 
              out2 << tgs[k][tgs[k].nelem()-1].Name() << "]\n";
              out2 << "  Selected lineshape: " << str << "\n";
              found0=i;
            }
        }

      // Is this normalization to the lineshape available?
      Index found1=-1;
      for ( Index i=0; i<lineshape_norm_data.nelem() && (found1 == -1); ++i )
        {
          const String& str = lineshape_norm_data[i].Name();
          if (str == normalizationfactor[k]) 
            {
              out2 << "  Selected normalization factor: " << normalizationfactor[k] << "\n";
              if ( (cutoff[k] != -1) && (cutoff[k] < 0.0) )
                throw runtime_error("  Cutoff must be -1 or positive.");
              out2 << "  Selected cutoff frequency    : " << cutoff[k] << "\n";
              found1=i;
            }
        }


      // did we find the lineshape and normalization factor?
      if (found0 == -1)
        {
          ostringstream os;
          os << "Selected lineshape not available: "<< shape[k] <<"\n";
          throw runtime_error(os.str());
        }
      if (found1 == -1)
        {
          ostringstream os;
          os << "Selected normalization to lineshape not available: "<< 
            normalizationfactor[k] <<"\n";
          throw runtime_error(os.str());
        }

      // now set the lineshape variables 
      abs_lineshape[k].SetInd_ls( found0 );
      abs_lineshape[k].SetInd_lsn( found1 );
      abs_lineshape[k].SetCutoff( cutoff[k] );
    }
}


//! abs_h2oSet.
/*!
 Sets abs_h2o to the profile of the first tag group containing
 water.
 
 This is necessary, because for example *abs_coefCalc* requires abs_h2o
 to contain the water vapour profile(the reason for this is the
 calculation of oxygen line broadening requires water vapour profile).
 Then this function can be used to copy the profile of the first tag
 group of water.
 
 \author Stefan Buehler
 
 \param[out] abs_h2o    WS Output
 \param[in]     abs_species WS Input
 \param[in]     abs_vmrs WS Input
 */
void abs_h2oSet(Vector&          abs_h2o,
                const ArrayOfArrayOfSpeciesTag& abs_species,
                const Matrix&    abs_vmrs,
                const Verbosity&)
{
  const Index h2o_index 
    = find_first_species_tg( abs_species,
                             species_index_from_species_name("H2O") );

  abs_h2o.resize( abs_vmrs.ncols() );
  if ( h2o_index < 0 )
    abs_h2o = -99;
  else
    abs_h2o = abs_vmrs(h2o_index,Range(joker));   
}


//!  abs_n2Set.
/*!
 Sets abs_n2 to the profile of the first tag group containing
 molecular nitrogen. See *abs_h2oSet* for more details.
 
 \author Stefan Buehler
 
 \param[out] abs_n2    WS Output
 \param[in]     abs_species WS Input
 \param[in]     abs_vmrs WS Input
 */
void abs_n2Set(Vector&            abs_n2,
               const ArrayOfArrayOfSpeciesTag& abs_species,
               const Matrix&    abs_vmrs,
               const Verbosity&)
{
  const Index n2_index 
    = find_first_species_tg( abs_species,
                             species_index_from_species_name("N2") );

  abs_n2.resize( abs_vmrs.ncols() );
  if ( n2_index < 0 )
    abs_n2 = -99;
  else
    abs_n2 = abs_vmrs(n2_index,Range(joker));   
}

//!  abs_o2Set.
/*!
 Sets abs_o2 to the profile of the first tag group containing
 molecular oxygen. See *abs_h2oSet* for more details.
 
 \author Mayuri Tatiya
 
 \param[out] abs_o2    WS Output
 \param[in]     abs_species WS Input
 \param[in]     abs_vmrs WS Input
 */
void abs_o2Set(Vector&            abs_o2,
               const ArrayOfArrayOfSpeciesTag& abs_species,
               const Matrix&    abs_vmrs,
               const Verbosity&)
{
  const Index o2_index 
    = find_first_species_tg( abs_species,
                             species_index_from_species_name("O2") );

  abs_o2.resize( abs_vmrs.ncols() );
  if ( o2_index < 0 )
    abs_o2 = -99;
  else
    abs_o2 = abs_vmrs(o2_index,Range(joker));   
}


/* Workspace method: Doxygen documentation will be auto-generated */
void AbsInputFromAtmFields (// WS Output:
                            Vector& abs_p,
                            Vector& abs_t,
                            Matrix& abs_vmrs,
                            // WS Input:
                            const Index& atmosphere_dim,
                            const Vector&  p_grid,
                            const Tensor3& t_field,
                            const Tensor4& vmr_field,
                            const Verbosity&
                           )
{
  // First, make sure that we really have a 1D atmosphere:
  if ( 1 != atmosphere_dim )
    {
      ostringstream os;
      os << "Atmospheric dimension must be 1D, but atmosphere_dim is "
         << atmosphere_dim << ".";
      throw runtime_error(os.str());
    }

  abs_p = p_grid;
  abs_t = t_field (joker, 0, 0);
  abs_vmrs = vmr_field (joker, joker, 0, 0);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_coefCalcFromXsec(// WS Output:
                          Matrix&              abs_coef,
                          Matrix&              src_coef,
                          ArrayOfMatrix&       dabs_coef_dx,
                          ArrayOfMatrix&       dsrc_coef_dx,
                          ArrayOfMatrix&       abs_coef_per_species,
                          ArrayOfMatrix&       src_coef_per_species,
                          // WS Input:         
                          const ArrayOfMatrix& abs_xsec_per_species,
                          const ArrayOfMatrix& src_xsec_per_species,
                          const ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx,
                          const ArrayOfArrayOfMatrix& dsrc_xsec_per_species_dx,
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          const ArrayOfRetrievalQuantity& jacobian_quantities,
                          const Matrix&        abs_vmrs,
                          const Vector&        abs_p,
                          const Vector&        abs_t,
                          const Verbosity&     verbosity)
{
  CREATE_OUT3;
  
  // Check that abs_vmrs and abs_xsec_per_species really have compatible
  // dimensions. In abs_vmrs there should be one row for each tg:
  if ( abs_vmrs.nrows() != abs_xsec_per_species.nelem() )
    {
      ostringstream os;
      os << "Variable abs_vmrs must have compatible dimension to abs_xsec_per_species.\n"
         << "abs_vmrs.nrows() = " << abs_vmrs.nrows() << "\n"
         << "abs_xsec_per_species.nelem() = " << abs_xsec_per_species.nelem();
      throw runtime_error(os.str());
    }

  // Check that number of altitudes are compatible. We only check the
  // first element, this is possilble because within arts all elements
  // are on the same altitude grid.
  if ( abs_vmrs.ncols() != abs_xsec_per_species[0].ncols() )
    {
      ostringstream os;
      os << "Variable abs_vmrs must have same numbers of altitudes as abs_xsec_per_species.\n"
         << "abs_vmrs.ncols() = " << abs_vmrs.ncols() << "\n"
         << "abs_xsec_per_species[0].ncols() = " << abs_xsec_per_species[0].ncols();
      throw runtime_error(os.str());
    }  

  // Check dimensions of abs_p and abs_t:
  chk_size("abs_p", abs_p, abs_vmrs.ncols());
  chk_size("abs_t", abs_t, abs_vmrs.ncols());

  
  // Initialize abs_coef and abs_coef_per_species. The array dimension of abs_coef_per_species
  // is the same as that of abs_xsec_per_species. The dimension of abs_coef should
  // be equal to one of the abs_xsec_per_species enries.
  abs_coef.resize( abs_xsec_per_species[0].nrows(), abs_xsec_per_species[0].ncols() );
  abs_coef = 0;  
  
  
  bool do_src = false;
  if( src_xsec_per_species.nelem() == abs_xsec_per_species.nelem() )
  {
    do_src = src_xsec_per_species[0].nrows() && src_xsec_per_species[0].ncols();
    src_coef.resize( src_xsec_per_species[0].nrows(), src_xsec_per_species[0].ncols() );
  }                    // Matpack can set all elements like this.
  
  const ArrayOfIndex jacobian_quantities_position = equivlent_propmattype_indexes(jacobian_quantities);
  dabs_coef_dx.resize(jacobian_quantities_position.nelem());
  for(Index ii=0;ii<jacobian_quantities_position.nelem();ii++)
      if(jacobian_quantities[jacobian_quantities_position[ii]] not_eq JacPropMatType::NotPropagationMatrixType)
      {
          dabs_coef_dx[ii].resize(abs_xsec_per_species[0].nrows(), abs_xsec_per_species[0].ncols());
          dabs_coef_dx[ii] = 0.0;
          if(do_src)
          {
              dsrc_coef_dx[ii].resize(src_xsec_per_species[0].nrows(), src_xsec_per_species[0].ncols());
              dsrc_coef_dx[ii] = 0.0;
          }
      }
  
  if(do_src)
    src_coef = 0;
  
  abs_coef_per_species.resize( abs_xsec_per_species.nelem() );
  if(do_src)
    src_coef_per_species.resize( src_xsec_per_species.nelem() );

  out3 << "  Computing abs_coef and abs_coef_per_species from abs_xsec_per_species.\n";
  // Loop through all tag groups
  for ( Index i=0; i<abs_xsec_per_species.nelem(); ++i )
    {
      out3 << "  Tag group " << i << "\n";

      // Make this element of abs_xsec_per_species the right size:
      abs_coef_per_species[i].resize( abs_xsec_per_species[i].nrows(), abs_xsec_per_species[i].ncols() );
      abs_coef_per_species[i] = 0;        // Initialize all elements to 0.
      if(do_src)
      {
        src_coef_per_species[i].resize( src_xsec_per_species[i].nrows(), src_xsec_per_species[i].ncols() );
        src_coef_per_species[i] = 0;        // Initialize all elements to 0.
      }

      // Loop through all altitudes
      for ( Index j=0; j<abs_xsec_per_species[i].ncols(); j++)
        {
          // Calculate total number density from pressure and temperature.
          const Numeric n = number_density(abs_p[j],abs_t[j]);
          const Numeric dn_dT = dnumber_density_dt(abs_p[j],abs_t[j]); 
          // Wasted calculations when Jacobians are not calculated... 
          // Though this is called seldom enough that it this fine?  value is -1/t*n

          // Loop through all frequencies
          for ( Index k=0; k<abs_xsec_per_species[i].nrows(); k++)
            {
              abs_coef_per_species[i](k,j) = abs_xsec_per_species[i](k,j) * n * abs_vmrs(i,j);
              if(do_src)
                src_coef_per_species[i](k,j) = src_xsec_per_species[i](k,j) * n * abs_vmrs(i,j);
              
              for(Index iq=0;iq<jacobian_quantities_position.nelem();iq++)
              {
                  if(jacobian_quantities[jacobian_quantities_position[iq]]==JacPropMatType::Temperature)
                  {
                      dabs_coef_dx[iq](k,j) += (dabs_xsec_per_species_dx[i][iq](k,j) * n + 
                      abs_xsec_per_species[i](k,j) * dn_dT) * abs_vmrs(i,j);
                      if(do_src)
                          dsrc_coef_dx[iq](k,j) += (dsrc_xsec_per_species_dx[i][iq](k,j) * n + 
                          src_xsec_per_species[i](k,j) * dn_dT) * abs_vmrs(i,j);
                  }
                  else if(jacobian_quantities[jacobian_quantities_position[iq]] == JacPropMatType::VMR)  // FIXME: Test that this works as expected using perturbations...
                  {
                    if(species_match(jacobian_quantities[jacobian_quantities_position[iq]], abs_species[i])) {
                      dabs_coef_dx[iq](k,j) += ( dabs_xsec_per_species_dx[i][iq](k,j) * abs_vmrs(i,j) +
                      abs_xsec_per_species[i](k,j) ) * n;
                      if(do_src)
                          dsrc_coef_dx[iq](k,j) += ( dsrc_xsec_per_species_dx[i][iq](k,j) * abs_vmrs(i,j) +
                          src_xsec_per_species[i](k,j) ) * n;
                    }
                  }
                  else if(jacobian_quantities[jacobian_quantities_position[iq]] not_eq JacPropMatType::NotPropagationMatrixType)
                  {
                      dabs_coef_dx[iq](k,j) += dabs_xsec_per_species_dx[i][iq](k,j) * n * abs_vmrs(i,j);
                      if(do_src)
                          dsrc_coef_dx[iq](k,j) += dsrc_xsec_per_species_dx[i][iq](k,j) * n * abs_vmrs(i,j);
                  }
              }
              
            }
        }

      // Add up to the total absorption:
      abs_coef += abs_coef_per_species[i];     // In Matpack you can use the +=
                                               // operator to do elementwise addition.
      if(do_src)
        src_coef += src_coef_per_species[i];
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesInit(// WS Output:
                              ArrayOfMatrix&   abs_xsec_per_species,
                              ArrayOfMatrix&   src_xsec_per_species,
                              ArrayOfArrayOfMatrix&   dabs_xsec_per_species_dx,
                              ArrayOfArrayOfMatrix&   dsrc_xsec_per_species_dx,
                              // WS Input:
                              const ArrayOfArrayOfSpeciesTag& tgs,
                              const ArrayOfRetrievalQuantity& jacobian_quantities,
                              const ArrayOfIndex& abs_species_active,
                              const Vector&    f_grid,
                              const Vector&    abs_p,
                              const Index&     abs_xsec_agenda_checked,
                              const Index&     nlte_do,
                              const Verbosity& verbosity
                              )
{
  CREATE_OUT3;

  if (!abs_xsec_agenda_checked)
    throw runtime_error("You must call *abs_xsec_agenda_checkedCalc* before calling this method.");

  // We need to check that abs_species_active doesn't have more elements than
  // abs_species (abs_xsec_agenda_checkedCalc doesn't know abs_species_active.
  // Usually we come here through an agenda call, where abs_species_active has
  // been properly created somewhere internally. But we might get here by
  // direct call, and then need to be safe!).
  if ( tgs.nelem() < abs_species_active.nelem() )
    {
      ostringstream os;
      os << "abs_species_active (n=" << abs_species_active.nelem()
         << ") not allowed to have more elements than abs_species (n="
         << tgs.nelem() << ")!\n";
      throw runtime_error(os.str());
    }

  // Initialize abs_xsec_per_species. The array dimension of abs_xsec_per_species
  // is the same as that of abs_lines_per_species.
  abs_xsec_per_species.resize( tgs.nelem() );
  src_xsec_per_species.resize( tgs.nelem() ); 
  
  const bool do_jac = supports_propmat_clearsky(jacobian_quantities);  //minus one is a flag for a non-species...
  const ArrayOfIndex jacobian_quantities_position = equivlent_propmattype_indexes(jacobian_quantities);
  
  dabs_xsec_per_species_dx.resize(do_jac?tgs.nelem():0);
  dsrc_xsec_per_species_dx.resize(do_jac?tgs.nelem():0);

  // Loop abs_xsec_per_species and make each matrix the right size,
  // initializing to zero.
  // But skip inactive species, loop only over the active ones.
  for ( Index ii=0; ii<abs_species_active.nelem(); ++ii )
    {
      const Index i = abs_species_active[ii];
      // Check that abs_species_active index is not higher than the number
      // of species
      if (i >= tgs.nelem())
      {
          ostringstream os;
          os << "*abs_species_active* contains an invalid species index.\n"
          << "Species index must be between 0 and " << tgs.nelem()-1;
          throw std::runtime_error(os.str());
      }
      // Make this element of abs_xsec_per_species the right size:
      abs_xsec_per_species[i].resize( f_grid.nelem(), abs_p.nelem() );
      abs_xsec_per_species[i] = 0;       // Matpack can set all elements like this.
      if (nlte_do)
      {
          src_xsec_per_species[i].resize( f_grid.nelem(), abs_p.nelem() );
          src_xsec_per_species[i] = 0;
      }
      else
      {
          src_xsec_per_species[i].resize(0, 0);
      }
     
     if(do_jac) {
       dabs_xsec_per_species_dx[ii] = ArrayOfMatrix(jacobian_quantities_position.nelem(), Matrix(f_grid.nelem(), abs_p.nelem(), 0.0));
        if(nlte_do)
          dsrc_xsec_per_species_dx[ii] = ArrayOfMatrix(jacobian_quantities_position.nelem(), Matrix(f_grid.nelem(), abs_p.nelem(), 0.0));
     }
     
    }

  ostringstream os;
  os << "  Initialized abs_xsec_per_species.\n"
     << "  Number of frequencies        : " << f_grid.nelem() << "\n"
     << "  Number of pressure levels    : " << abs_p.nelem() << "\n";
  out3 << os.str();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddLines(// WS Output:
                                  ArrayOfMatrix&                   abs_xsec_per_species,
				  ArrayOfMatrix&                   src_xsec_per_species,
                                  ArrayOfArrayOfMatrix&            dabs_xsec_per_species_dx,
                                  ArrayOfArrayOfMatrix&            dsrc_xsec_per_species_dx,
                                  // WS Input:             
                                  const ArrayOfArrayOfSpeciesTag&  tgs,
                                  const ArrayOfRetrievalQuantity& jacobian_quantities,
                                  const ArrayOfIndex& abs_species_active,
                                  const Vector&                    f_grid,
                                  const Vector&                    abs_p,
                                  const Vector&                    abs_t,
                                  const Matrix&                    abs_nlte,
                                  const Numeric&                   lm_p_lim,
                                  const Matrix&                    abs_vmrs,
                                  const ArrayOfArrayOfLineRecord&  abs_lines_per_species,
                                  const ArrayOfLineshapeSpec&      abs_lineshape,
                                  const SpeciesAuxData&            isotopologue_ratios,
                                  const SpeciesAuxData&            partition_functions,
                                  const Verbosity&                 verbosity)
{
  CREATE_OUT3;
  
  // Check that correct isotopologue ratios are defined for the species
  // we want to calculate
  checkIsotopologueRatios(tgs, isotopologue_ratios);
  
  // Check that all temperatures are at least 0 K. (Negative Kelvin
  // temperatures are unphysical.)  
  if ( min(abs_t) < 0 )
    {
      ostringstream os;
      os << "Temperature must be at least 0 K. But you request an absorption\n"
         << "calculation at " << min(abs_t) << " K!"; 
      throw std::runtime_error(os.str());
    }

  // Check that all parameters that should have the number of tag
  // groups as a dimension are consistent.
  {
    const Index n_tgs    = tgs.nelem();
    const Index n_xsec   = abs_xsec_per_species.nelem();
    const Index n_vmrs   = abs_vmrs.nrows();
    const Index n_lines  = abs_lines_per_species.nelem();
    const Index n_shapes = abs_lineshape.nelem();

    if ( n_tgs != n_xsec  ||
         n_tgs != n_vmrs  ||
         n_tgs != n_lines ||
         ( n_tgs != n_shapes &&
           1     != n_shapes ) )
      {
        ostringstream os;
        os << "The following variables must all have the same dimension:\n"
           << "tgs:          " << tgs.nelem() << "\n"
           << "abs_xsec_per_species:  " << abs_xsec_per_species.nelem() << "\n"
           << "abs_vmrs:         " << abs_vmrs.nrows() << "\n"
           << "abs_lines_per_species: " << abs_lines_per_species.nelem() << "\n"
           << "abs_lineshape:    " << abs_lineshape.nelem() << "\n"
           << "(As a special case, abs_lineshape is allowed to have only one element.)";
        throw std::runtime_error(os.str());
      }
  }  
  
  // Print information:
  //
  out3 << "  Calculating line spectra.\n";
  //
  // Uncomment the part below if you temporarily want detailed info about 
  // transitions to be done

      // The variables defined here (in particular the frequency
      // conversion) are just to make the output nice. They are not used
      // in subsequent calculations.
  //    cout << "  Transitions to do: \n";
  //    Index nlines = 0;
  //    String funit;
  //    Numeric ffac;
  //    if ( f_grid[0] < 3e12 )
  //      {
  //        funit = "GHz"; ffac = 1e9;
  //      }
  //    else
  //      {
  //        extern const Numeric SPEED_OF_LIGHT;
  //        funit = "cm-1"; ffac = SPEED_OF_LIGHT*100;
  //      }
  //    for ( Index i=0; i<abs_lines_per_species.nelem(); ++i )
  //      {
  //        for ( Index l=0; l<abs_lines_per_species[i].nelem(); ++l )
  //          {
  //          cout << "    " << abs_lines_per_species[i][l].Name() << " @ " 
  //               << abs_lines_per_species[i][l].F()/ffac  << " " << funit << " ("
  //               << abs_lines_per_species[i][l].I0() << "  "
  //               << abs_lines_per_species[i][l].Agam() << ")\n"; 
  //          nlines++;
  //        }
  //    }
  //  out2 << "  Total number of transistions : " << nlines << "\n";

  // Call xsec_species for each tag group.
  for ( Index ii=0; ii<abs_species_active.nelem(); ++ii )
    {
      const Index i = abs_species_active[ii];

      // Get a pointer to the line list for the current species. This
      // is just so that we don't have to type abs_lines_per_species[i] over
      // and over again.
      const ArrayOfLineRecord& ll = abs_lines_per_species[i];

      // Also get a pointer to the lineshape specification. This
      // requires special treatment: If there is only 1 lineshape
      // given, the same line shape should be used for all species.
      LineshapeSpec ls;
      if (1==abs_lineshape.nelem())
        ls = abs_lineshape[0];
      else
        ls = abs_lineshape[i];
      
      
      // We do the LBL calculation only if:
      // - The line list is not empty, and
      // - The species is not a Zeeman species.
      if ( 0 < ll.nelem() && tgs[i].nelem() && !is_zeeman(tgs[i]) )
        {
          // As a safety check, check that the species of the first
          // line matches the species we should have according to
          // tgs. (This in case the order in tgs has been changed and
          // abs_lines_per_species has not been changed consistently.)
          if (ll[0].Species() != tgs[i][0].Species() )
            {
              ostringstream os;
              os << "The species in the line list does not match the species\n"
                 << "for which you want to calculate absorption:\n"
                 << "abs_species:           " << get_tag_group_name(tgs[i]) << "\n"
                 << "abs_lines_per_species: " << ll[0].Name();
              throw std::runtime_error(os.str());
            }

          // Get the name of the species. The member function name of a
          // LineRecord returns the full name (species + isotopologue). So
          // it is for us more convenient to get the species index
          // from the LineRecord member function Species(), and then
          // use this to look up the species name in species_data.
          using global_data::species_data;
          String species_name = species_data[ll[0].Species()].Name();
          ArrayOfIndex jacobian_quantities_position = equivlent_propmattype_indexes(jacobian_quantities);
          const bool do_jac = jacobian_quantities_position.nelem();
          
          // Get the name of the lineshape. For that we use the member
          // function Ind_ls() to the lineshape data ls, which returns
          // an index. With that index we can go into lineshape_data
          // to get the name.
          // We will need this for safety checks later on.
          using global_data::lineshape_data;
          String lineshape_name = lineshape_data[ls.Ind_ls()].Name();


          // The use of overlap parameters for oxygen lines only makes
          // sense if the special Rosenkranz lineshape is used
          // (Rosenkranz_Voigt_Drayson or Rosenkranz_Voigt_Kuntz6). 
          // Likewise, the use of these lineshapes only makes sense if
          // overlap parameters are available. We will test both these
          // conditions.

          // First of all, let's find out if the species we are
          // dealing with is oxygen. 
          if ( "O2" == species_name )
            {
              // Do we have overlap parameters in the aux fields of
              // the LineRecord?
              if ( 0 != ll[0].Naux() )
                {
                  // Yes. So let's make sure that the lineshape is one
                  // that can use these parameters. 
                  if ( "Rosenkranz_Voigt_Drayson" != lineshape_name &&
                       "Rosenkranz_Voigt_Kuntz6"  != lineshape_name    )
                    {
                      ostringstream os;
                      os 
                        << "You are using a line catalogue that contains auxiliary parameters to\n"
                        << "take care of overlap for oxygen lines. But you are not using a\n"
                        << "lineshape that uses these parameters. Use Rosenkranz_Voigt_Drayson or\n"
                        << "Rosenkranz_Voigt_Kuntz6.";
                      throw runtime_error(os.str());                  
                    }
                }               
            }

          // Now we go the opposite way. Let's see if the Rosenkranz
          // lineshapes are used.
          if ( "Rosenkranz_Voigt_Drayson" == lineshape_name ||
               "Rosenkranz_Voigt_Kuntz6"  == lineshape_name    )
            {
              // Is the species oxygen, as it should be?
              if ( "O2" != species_name )
                {
                  ostringstream os;
                  os 
                    << "You are trying to use one of the special Rosenkranz lineshapes with\n"
                    << "overlap correction for a species other than oxygen. Your species is\n"
                    << species_name << ".\n"
                    << "Select another lineshape for this species.";
                    throw runtime_error(os.str());                    
                }
              else
                {
                  // Do we have overlap parameters in the aux fields of
                  // the LineRecord?
                  if ( 0 == ll[0].Naux() )
                    {
                      ostringstream os;
                      os 
                        << "You are trying to use one of the special Rosenkranz lineshapes with\n"
                        << "overlap correction. But your line file does not contain aux\n"
                        << "parameters. (I've checked only the first LineRecord). Use a line file\n"
                        << "with overlap parameters.";
                        throw runtime_error(os.str());                
                    }
                }
            }
            if (tgs[i][0].LineMixing() == SpeciesTag::LINE_MIXING_OFF && !do_jac )
            {
                Matrix dummy_phase;
                xsec_species(abs_xsec_per_species[i],
                             src_xsec_per_species[i],
                             dummy_phase,
                             f_grid,
                             abs_p,
                             abs_t,
                             abs_nlte,
                             abs_vmrs,
                             tgs,
                             i,
                             ll,
                             ls.Ind_ls(),
                             ls.Ind_lsn(),
                             ls.Cutoff(),
                             isotopologue_ratios,
                             partition_functions,
                             verbosity );
            }
          else
            {
                Matrix dummy_phase(abs_xsec_per_species[i].nrows(),
                                   abs_xsec_per_species[i].ncols(),
                                   0.);
                // This is a bad idea?
                ArrayOfMatrix dummy;
                
                xsec_species_line_mixing_wrapper(abs_xsec_per_species[i],
                                                 src_xsec_per_species[i],
                                                 dummy_phase,
                                                 do_jac?dabs_xsec_per_species_dx[i]:dummy,
                                                 abs_nlte.empty()?dummy:dsrc_xsec_per_species_dx[i],
                                                 dummy,
                                                 jacobian_quantities,
                                                 jacobian_quantities_position,
                                                 f_grid,
                                                 abs_p,
                                                 abs_t,
                                                 abs_nlte,
                                                 abs_vmrs,
                                                 tgs,
                                                 i,
                                                 ll,
                                                 0.0,
                                                 ls.Ind_ls(),
                                                 ls.Ind_lsn(),
                                                 lm_p_lim,
                                                 ls.Cutoff(),
                                                 isotopologue_ratios,
                                                 partition_functions,
                                                 verbosity);
            }
          // Note that we call xsec_species with a row of abs_vmrs,
          // selected by the above Matpack expression. This is
          // possible, because xsec_species is using Views.
          
        }

      if (out3.sufficient_priority())
        {
          ostringstream os;
          os << "  Tag group " << i
             << " (" << get_tag_group_name(tgs[i]) << "): "
             << ll.nelem() << " transitions\n";
          out3 << os.str();
        }

    } // End of species for loop.
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddConts(// WS Output:
                                  ArrayOfMatrix&                   abs_xsec_per_species,
                                  ArrayOfArrayOfMatrix&            dabs_xsec_per_species_dx,
                                  // WS Input:             
                                  const ArrayOfArrayOfSpeciesTag&  tgs,
                                  const ArrayOfRetrievalQuantity& jacobian_quantities,
                                  const ArrayOfIndex& abs_species_active,
                                  const Vector&                    f_grid,
                                  const Vector&                    abs_p,
                                  const Vector&                    abs_t,
                                  const Matrix&                    abs_vmrs,
                                  const ArrayOfString&             abs_cont_names,
                                  const ArrayOfVector&             abs_cont_parameters,
                                  const ArrayOfString&             abs_cont_models,
                                  const Verbosity&                 verbosity)
{
  CREATE_OUT3;
  
  // Needed for some continua, and set here from abs_vmrs:
  Vector abs_h2o, abs_n2, abs_o2;
  
  // Check that all paramters that should have the number of tag
  // groups as a dimension are consistent.
  {
    const Index n_tgs    = tgs.nelem();
    const Index n_xsec   = abs_xsec_per_species.nelem();
    const Index n_vmrs   = abs_vmrs.nrows();

    if ( n_tgs != n_xsec || n_tgs != n_vmrs )
      {
        ostringstream os;
        os << "The following variables must all have the same dimension:\n"
           << "tgs:          " << tgs.nelem() << "\n"
           << "abs_xsec_per_species:  " << abs_xsec_per_species.nelem() << "\n"
           << "abs_vmrs.nrows():      " << abs_vmrs.nrows();
        throw runtime_error(os.str());
      }
  }
  
  // Jacobian overhead START
  /* NOTE: if any of the functions inside continuum tags could 
           be made to give partial derivatives, then that would 
           speed things up.  Also be aware that line specific
           parameters cannot be retrieved while using these 
           models. */
  const bool do_jac = supports_continuum(jacobian_quantities);  // Throws runtime error if line parameters are wanted since we cannot know if the line is in the Continuum...
  const bool do_freq_jac = do_frequency_jacobian(jacobian_quantities);
  const bool do_temp_jac = do_temperature_jacobian(jacobian_quantities);
  Vector dfreq, dabs_t;
  const Numeric df = frequency_perturbation(jacobian_quantities);
  const Numeric dt = temperature_perturbation(jacobian_quantities);
  const ArrayOfIndex jacobian_quantities_position = equivlent_propmattype_indexes(jacobian_quantities);
  
  if(do_freq_jac)
  {
      dfreq.resize(f_grid.nelem());
      for(Index iv=0; iv<f_grid.nelem(); iv++)
          dfreq[iv] = f_grid[iv]+df;
  }
  if(do_temp_jac)
  {
      dabs_t.resize(abs_t.nelem());
      for(Index it=0; it<abs_t.nelem(); it++)
          dabs_t[it] = abs_t[it]+dt;
  }
  
  Matrix jacs_df, jacs_dt, normal;
  if(do_jac)
  {
      if(do_freq_jac) jacs_df.resize(f_grid.nelem(), abs_p.nelem());
      if(do_temp_jac) jacs_dt.resize(f_grid.nelem(), abs_p.nelem());
      normal.resize(f_grid.nelem(), abs_p.nelem());
  }
  // Jacobian overhead END

  // Check, that dimensions of abs_cont_names and
  // abs_cont_parameters are consistent...
  if ( abs_cont_names.nelem() !=
       abs_cont_parameters.nelem() )
    {
      ostringstream os;

      for (Index i=0; i < abs_cont_names.nelem(); ++i) 
        os << "abs_xsec_per_speciesAddConts: " << i << " name : " << abs_cont_names[i] << "\n";

      for (Index i=0; i < abs_cont_parameters.nelem(); ++i) 
        os << "abs_xsec_per_speciesAddConts: " << i << " param: " << abs_cont_parameters[i] << "\n";

      for (Index i=0; i < abs_cont_models.nelem(); ++i) 
        os << "abs_xsec_per_speciesAddConts: " << i << " option: " << abs_cont_models[i] << "\n";

      os << "The following variables must have the same dimension:\n"
         << "abs_cont_names:      " << abs_cont_names.nelem() << "\n"
         << "abs_cont_parameters: " << abs_cont_parameters.nelem();

      throw runtime_error(os.str());
    }

  // Check that abs_p, abs_t, and abs_vmrs have the same
  // dimension. This could be a user error, so we throw a
  // runtime_error. 

  if ( abs_t.nelem() != abs_p.nelem() )
    {
      ostringstream os;
      os << "Variable abs_t must have the same dimension as abs_p.\n"
         << "abs_t.nelem() = " << abs_t.nelem() << '\n'
         << "abs_p.nelem() = " << abs_p.nelem();
      throw runtime_error(os.str());
    }

  if ( abs_vmrs.ncols() != abs_p.nelem() )
    {
      ostringstream os;
      os << "Variable dimension abs_vmrs.ncols() must\n"
         << "be the same as abs_p.nelem().\n"
         << "abs_vmrs.ncols() = " << abs_vmrs.ncols() << '\n'
         << "abs_p.nelem() = " << abs_p.nelem();
      throw runtime_error(os.str());
    }

  // We set abs_h2o, abs_n2, and abs_o2 later, because we only want to
  // do it if the parameters are really needed.


  out3 << "  Calculating continuum spectra.\n";

  // Loop tag groups:
  for ( Index ii=0; ii<abs_species_active.nelem(); ++ii )
    {
      const Index i = abs_species_active[ii];

      using global_data::species_data;

      // Go through the tags in the current tag group to see if they
      // are continuum tags:  
      for ( Index s=0; s<tgs[i].nelem(); ++s )
        {
          // Continuum tags in the sense that we talk about here
          // (including complete absorption models) are marked by a special type.
          if (tgs[i][s].Type() == SpeciesTag::TYPE_PREDEF)
            {
                  // We have identified a continuum tag!

                  // Get only the continuum name. The full tag name is something like:
                  // H2O-HITRAN96Self-*-*. We want only the `H2O-HITRAN96Self' part:
                  const String name =
                    species_data[tgs[i][s].Species()].Name() + "-"
                    + species_data[tgs[i][s].Species()].Isotopologue()[tgs[i][s].Isotopologue()].Name();
  
                  // Check, if we have parameters for this model. For
                  // this, the model name must be listed in
                  // abs_cont_names.
                  const Index n =
                    find( abs_cont_names.begin(),
                          abs_cont_names.end(),
                          name ) - abs_cont_names.begin();

                  // n==abs_cont_names.nelem() indicates that
                  // the name was not found.
                  if ( n==abs_cont_names.nelem() )
                    {
                      ostringstream os;
                      os << "Cannot find model " << name
                         << " in abs_cont_names.";
                      throw runtime_error(os.str());                  
                    }

                  // Ok, the tag specifies a valid continuum model and
                  // we have continuum parameters.

                  if (out3.sufficient_priority())
                    {
                      ostringstream os;
                      os << "  Adding " << name
                         << " to tag group " << i << ".\n";
                      out3 << os.str();
                    }

                  // find the options for this continuum tag from the input array
                  // of options. The actual field of the array is n:
                  const String ContOption = abs_cont_models[n];

                  // Set abs_h2o, abs_n2, and abs_o2 from the first matching species.
                  abs_h2oSet(abs_h2o, tgs, abs_vmrs, verbosity);
                  abs_n2Set(abs_n2, tgs, abs_vmrs, verbosity);
                  abs_o2Set(abs_o2, tgs, abs_vmrs, verbosity);
 
                  // Add the continuum for this tag. The parameters in
                  // this call should be clear. The vmr is in
                  // abs_vmrs(i,Range(joker)). The other vmr variables,
                  // abs_h2o, abs_n2, and abs_o2 contains the real vmr of H2O,
                  // N2, nad O2, which are needed as additional information for
                  // certain continua:
                  // abs_h2o for
                  //   O2-PWR88, O2-PWR93, O2-PWR98,
                  //   O2-MPM85, O2-MPM87, O2-MPM89, O2-MPM92, O2-MPM93,
                  //   O2-TRE05,
                  //   O2-SelfContStandardType, O2-SelfContMPM93, O2-SelfContPWR93,
                  //   N2-SelfContMPM93, N2-DryContATM01,
                  //   N2-CIArotCKDMT252, N2-CIAfunCKDMT252
                  // abs_n2 for
                  //   H2O-SelfContCKD24, H2O-ForeignContCKD24,
                  //   O2-v0v0CKDMT100,
                  //   CO2-ForeignContPWR93, CO2-ForeignContHo66
                  // abs_o2 for
                  //   N2-CIArotCKDMT252, N2-CIAfunCKDMT252
                  if(!do_jac)
                    xsec_continuum_tag( abs_xsec_per_species[i],
                                        name,
                                        abs_cont_parameters[n],
                                        abs_cont_models[n], 
                                        f_grid,
                                        abs_p,
                                        abs_t,
                                        abs_n2,
                                        abs_h2o,
                                        abs_o2,
                                        abs_vmrs(i,Range(joker)),
                                        verbosity );
                  else // The Jacobian block
                  {
                      // Needs a reseted block here...
                      for(Index iv = 0; iv<f_grid.nelem(); iv++)
                      {
                          for(Index ip = 0; ip<abs_p.nelem(); ip++)
                          {
                              if(do_freq_jac) jacs_df(iv,ip) = 0.0;
                              if(do_temp_jac) jacs_dt(iv,ip) = 0.0;
                              normal(iv,ip)                  = 0.0;
                          }
                      }
                      
                      // Normal calculations
                      xsec_continuum_tag( normal,name,abs_cont_parameters[n],
                                          abs_cont_models[n], f_grid,abs_p,abs_t,
                                          abs_n2,abs_h2o,abs_o2,abs_vmrs(i,Range(joker)),verbosity );
                      
                      // Frequency calculations
                      if(do_freq_jac)
                          xsec_continuum_tag( jacs_df,name,abs_cont_parameters[n],
                                              abs_cont_models[n], dfreq, abs_p,abs_t,
                                              abs_n2,abs_h2o,abs_o2,abs_vmrs(i,Range(joker)),verbosity );
                      
                      //Temperature calculations
                      if(do_temp_jac)
                          xsec_continuum_tag( jacs_dt,name,abs_cont_parameters[n],
                                              abs_cont_models[n], f_grid, abs_p,dabs_t,
                                              abs_n2,abs_h2o,abs_o2,abs_vmrs(i,Range(joker)),verbosity );
                      for(Index iv = 0; iv<f_grid.nelem(); iv++)
                      {
                          for(Index ip = 0; ip<abs_p.nelem(); ip++)
                          {
                              abs_xsec_per_species[i](iv,ip) += normal(iv,ip);
                              for(Index iq=0; iq<jacobian_quantities_position.nelem(); iq++)
                              {
                                  if(is_frequency_parameter(jacobian_quantities[jacobian_quantities_position[iq]]))
                                      dabs_xsec_per_species_dx[i][iq](iv,ip) += (jacs_df(iv,ip)-normal(iv,ip))*(1./df);
                                  else if(jacobian_quantities[jacobian_quantities_position[iq]] == JacPropMatType::Temperature)
                                      dabs_xsec_per_species_dx[i][iq](iv,ip) += (jacs_dt(iv,ip)-normal(iv,ip))*(1./dt);
                              }
                          }
                      }
                  }
                  // Calling this function with a row of Matrix abs_vmrs
                  // is possible because it uses Views.
                }
            }
        }
    

}



//======================================================================
//             Methods related to continua
//======================================================================

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_cont_descriptionInit(// WS Output:
                              ArrayOfString& abs_cont_names,
                              ArrayOfString& abs_cont_options,
                              ArrayOfVector& abs_cont_parameters,
                              const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  abs_cont_names.resize(0);
  abs_cont_options.resize(0);
  abs_cont_parameters.resize(0);
  out2 << "  Initialized abs_cont_names \n"
          "              abs_cont_models\n"
          "              abs_cont_parameters.\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_cont_descriptionAppend(// WS Output:
                                ArrayOfString& abs_cont_names,
                                ArrayOfString& abs_cont_models,
                                ArrayOfVector& abs_cont_parameters,
                                // Control Parameters:
                                const String& tagname,
                                const String& model,
                                const Vector& userparameters,
                                const Verbosity&)
{
  // First we have to check that name matches a continuum species tag.
  check_continuum_model(tagname);

  //cout << "   + tagname:    " << tagname << "\n";
  //cout << "   + model:      " << model << "\n";
  //cout << "   + parameters: " << userparameters << "\n";

  // Add name and parameters to the apropriate variables:
  abs_cont_names.push_back(tagname);
  abs_cont_models.push_back(model);
  abs_cont_parameters.push_back(userparameters);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void nlte_sourceFromTemperatureAndSrcCoefPerSpecies(// WS Output:
                                                    ArrayOfStokesVector& nlte_source,
                                                    ArrayOfStokesVector& dnlte_dx_source,
                                                    ArrayOfStokesVector& nlte_dsource_dx,
                                                    // WS Input:
                                                    const ArrayOfMatrix& src_coef_per_species,
                                                    const ArrayOfMatrix& dsrc_coef_dx,
                                                    const ArrayOfRetrievalQuantity& jacobian_quantities,
                                                    const Vector&        f_grid,
                                                    const Numeric&       rtp_temperature,
                                                    const Verbosity&)
{
  // nlte_source has format
  // [ abs_species, f_grid, stokes_dim ].
  // src_coef_per_species has format ArrayOfMatrix (over species),
  // where for each species the matrix has format [f_grid, abs_p].
  
  Index n_species = src_coef_per_species.nelem(); // # species

  if (0==n_species)
    {
      ostringstream os;
      os << "Must have at least one species.";
      throw runtime_error(os.str());
    }

  Index n_f       = src_coef_per_species[0].nrows(); // # frequencies

  // # pressures must be 1:
  if (1!=src_coef_per_species[0].ncols())
    {
      ostringstream os;
      os << "Must have exactly one pressure.";
      throw runtime_error(os.str());
    }
  
  // Check species dimension of propmat_clearsky
  if ( nlte_source.nelem()!=n_species )
  {
    ostringstream os;
    os << "Species dimension of propmat_clearsky does not\n"
       << "match src_coef_per_species.";
    throw std::runtime_error( os.str() );
  }
  
  // Check frequency dimension of propmat_clearsky
  if ( nlte_source[0].NumberOfFrequencies()!=n_f )
  {
    ostringstream os;
    os << "Frequency dimension of propmat_clearsky does not\n"
       << "match abs_coef_per_species.";
    throw runtime_error( os.str() );
  }
  
  const ArrayOfIndex jacobian_quantities_position = equivlent_propmattype_indexes(jacobian_quantities);
  Vector B(n_f);
  
  for(Index iv=0; iv<n_f; iv++)
    B[iv] = planck(f_grid[iv], rtp_temperature);
  
  StokesVector sv(n_f, nlte_source[0].StokesDimensions());
  for ( Index si=0; si<n_species; ++si )
  {
    sv.Kjj() = src_coef_per_species[si](joker, 0);
    sv *= B;
    nlte_source[si].Kjj() += sv.Kjj();
  }
  
  // Jacobian
  for(Index ii = 0; ii<jacobian_quantities_position.nelem(); ii++)
  {
    if(jacobian_quantities[jacobian_quantities_position[ii]]==JacPropMatType::Temperature)
    {
      Vector dB(n_f);
      for(Index iv=0; iv<n_f; iv++)
        dB[iv] = dplanck_dt(f_grid[iv], rtp_temperature);
      
      for( Index si=0; si<n_species; ++si )
      {
        sv.Kjj() = src_coef_per_species[si](joker, 0);
        sv *= dB;
        nlte_dsource_dx[ii].Kjj() += sv.Kjj();
      }
      
      sv.Kjj() = dsrc_coef_dx[ii](joker,0);
      sv *= B;
      dnlte_dx_source[ii].Kjj() += sv.Kjj();
    }
    else if(is_frequency_parameter(jacobian_quantities[jacobian_quantities_position[ii]]))
    {
      Vector dB(n_f);
      for(Index iv=0; iv<n_f; iv++)
        dB[iv] = dplanck_df(f_grid[iv],rtp_temperature);
      
      for( Index si=0; si<n_species; ++si )
      {
        sv.Kjj() = src_coef_per_species[si](joker, 0);
        sv *= dB;
        nlte_dsource_dx[ii].Kjj() += sv.Kjj();
      }
      
      sv.Kjj() = dsrc_coef_dx[ii](joker,0);
      sv *= B;
      dnlte_dx_source[ii].Kjj() += sv.Kjj();
    }
    else if(jacobian_quantities[jacobian_quantities_position[ii]] not_eq JacPropMatType::NotPropagationMatrixType)
    {
      sv.Kjj() = dsrc_coef_dx[ii](joker, 0);
      sv *= B;
      dnlte_dx_source[ii].Kjj() += sv.Kjj();
    }
    else{/* All is fine! */}
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddFromAbsCoefPerSpecies(// WS Output:
                               ArrayOfPropagationMatrix&  propmat_clearsky,
                               ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                               // WS Input:
                               const ArrayOfMatrix& abs_coef_per_species,
                               const ArrayOfMatrix& dabs_coef_dx,
                               const Verbosity&)
{
  // propmat_clearsky has format
  // [ abs_species, f_grid, stokes_dim, stokes_dim ].
  // abs_coef_per_species has format ArrayOfMatrix (over species),
  // where for each species the matrix has format [f_grid, abs_p].
  
  Index n_species = abs_coef_per_species.nelem(); // # species

  if (0==n_species)
  {
    ostringstream os;
    os << "Must have at least one species.";
    throw runtime_error(os.str());
  }

  Index n_f = abs_coef_per_species[0].nrows(); // # frequencies

  // # pressures must be 1:
  if (1 not_eq abs_coef_per_species[0].ncols())
    {
      ostringstream os;
      os << "Must have exactly one pressure.";
      throw runtime_error(os.str());
    }
  
  // Check species dimension of propmat_clearsky
  if ( propmat_clearsky.nelem() not_eq n_species )
  {
    ostringstream os;
    os << "Species dimension of propmat_clearsky does not\n"
       << "match abs_coef_per_species.";
    throw runtime_error( os.str() );
  }
  
  // Check frequency dimension of propmat_clearsky
  if ( propmat_clearsky[0].NumberOfFrequencies() not_eq n_f )
  {
    ostringstream os;
    os << "Frequency dimension of propmat_clearsky does not\n"
       << "match abs_coef_per_species.";
    throw runtime_error( os.str() );
  }
  
  // Loop species and stokes dimensions, and add to propmat_clearsky:
  for ( Index si=0; si<n_species; ++si )
      propmat_clearsky[si].Kjj() += abs_coef_per_species[si](joker,0);
  
  for(Index iqn=0; iqn<dabs_coef_dx.nelem();iqn++)
  {
    if(dabs_coef_dx[iqn].nrows()==n_f)
    {
      if(dabs_coef_dx[iqn].ncols()==1)
      {
        dpropmat_clearsky_dx[iqn].Kjj() += dabs_coef_dx[iqn](joker,0);
      }
      else throw std::runtime_error( "Must have exactly one pressure.");
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyInit(   //WS Output
                             ArrayOfPropagationMatrix& propmat_clearsky,
                             ArrayOfStokesVector& nlte_source,
                             ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                             ArrayOfStokesVector& dnlte_dx_source,
                             ArrayOfStokesVector& nlte_dsource_dx,
                             //WS Input
                             const ArrayOfArrayOfSpeciesTag& abs_species,
                             const ArrayOfRetrievalQuantity& jacobian_quantities,
                             const Vector&                   f_grid,
                             const Index&                    stokes_dim,
                             const Index&                    propmat_clearsky_agenda_checked,
                             const Index&                    nlte_do,
                             const Verbosity&
                            )
{
  if (!propmat_clearsky_agenda_checked)
    throw runtime_error("You must call *propmat_clearsky_agenda_checkedCalc* before calling this method.");
  
  Index nf = f_grid.nelem();
  
  if(abs_species.nelem() > 0 )
  {
    if(nf > 0)
    {
      if(stokes_dim > 0)
      {
        propmat_clearsky = ArrayOfPropagationMatrix(abs_species.nelem(), PropagationMatrix(nf, stokes_dim)); 
        
        if (nlte_do)
          nlte_source = ArrayOfStokesVector(abs_species.nelem(), StokesVector(nf, stokes_dim));
        else
          nlte_source.resize(0);
      }
      else throw  runtime_error("stokes_dim = 0");
    }
    else throw runtime_error("nf = 0");
  }
  else throw runtime_error("abs_species.nelem() = 0");
  
  if(supports_propmat_clearsky(jacobian_quantities))
  {
    const Index nq = equivlent_propmattype_indexes(jacobian_quantities).nelem();
    dpropmat_clearsky_dx = ArrayOfPropagationMatrix(nq, PropagationMatrix(nf, stokes_dim));
    if(nlte_do)
    {
      dnlte_dx_source = ArrayOfStokesVector(nq, StokesVector(nf, stokes_dim));
      nlte_dsource_dx = ArrayOfStokesVector(nq, StokesVector(nf, stokes_dim));
    }
    else
    {
      dnlte_dx_source.resize(0);
      nlte_dsource_dx.resize(0);
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddFaraday(
         ArrayOfPropagationMatrix& propmat_clearsky,
         ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
   const Index&                    stokes_dim,
   const Index&                    atmosphere_dim,
   const Vector&                   f_grid,
   const ArrayOfArrayOfSpeciesTag& abs_species,
   const ArrayOfRetrievalQuantity& jacobian_quantities,
   const Vector&                   rtp_vmr,
   const Vector&                   rtp_los,
   const Vector&                   rtp_mag,
   const Verbosity& )
{
  // All the physical constants joined into one static constant:
  // (abs as e defined as negative)
  static const Numeric FRconst = abs( 
                        ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE / 
                        ( 8 * PI * PI * SPEED_OF_LIGHT * VACUUM_PERMITTIVITY * 
                          ELECTRON_MASS * ELECTRON_MASS ) );
  
  if( stokes_dim < 3 )
    throw runtime_error( 
                 "To include Faraday rotation, stokes_dim >= 3 is required." );

  if( atmosphere_dim==1 && rtp_los.nelem() < 1 )
    {
       ostringstream os; 
       os << "For applying propmat_clearskyAddFaraday, los needs to be specified\n"
          << "(at least zenith angle component for atmosphere_dim==1),\n"
          << "but it is not.\n";
       throw runtime_error( os.str() );
    }
  else if(  atmosphere_dim>1 && rtp_los.nelem() < 2 )
    {
       ostringstream os; 
       os << "For applying propmat_clearskyAddFaraday, los needs to be specified\n"
          << "(both zenith and azimuth angle components for atmosphere_dim>1),\n"
          << "but it is not.\n";
       throw runtime_error( os.str() );
    }
    
  const bool do_jac = supports_faraday(jacobian_quantities);
  const bool do_magn_jac = do_magnetic_jacobian(jacobian_quantities);
  const Numeric dmag = magnetic_field_perturbation(jacobian_quantities);
  const ArrayOfIndex jacobian_quantities_position = equivlent_propmattype_indexes(jacobian_quantities);
  
  Index ife = -1;  
  for( Index sp = 0; sp < abs_species.nelem() && ife < 0; sp++ )
    {
      if (abs_species[sp][0].Type() == SpeciesTag::TYPE_FREE_ELECTRONS)
        { ife = sp; }
    }

  if( ife < 0 )
    {
      throw runtime_error( "Free electrons not found in *abs_species* and "
                           "Faraday rotation can not be calculated." );
    }
  else
    {
      const Numeric ne = rtp_vmr[ife];

      if( ne!=0  &&  ( rtp_mag[0]!=0 || rtp_mag[1]!=0 || rtp_mag[2]!=0 ) )
        {
          // Include remaining terms, beside /f^2
          const Numeric c1 = 2 * FRconst * ne * dotprod_with_los( 
                 rtp_los, rtp_mag[0], rtp_mag[1], rtp_mag[2], atmosphere_dim );
          
          Numeric dc1_u=0.0, dc1_v=0.0, dc1_w=0.0;
          if(do_magn_jac)
          {
              dc1_u = (2 * FRconst * ne * dotprod_with_los( 
              rtp_los, rtp_mag[0]+dmag, rtp_mag[1], rtp_mag[2], atmosphere_dim ) - c1)/dmag;
              dc1_v = (2 * FRconst * ne * dotprod_with_los( 
              rtp_los, rtp_mag[0], rtp_mag[1]+dmag, rtp_mag[2], atmosphere_dim ) - c1)/dmag;
              dc1_w = (2 * FRconst * ne * dotprod_with_los( 
              rtp_los, rtp_mag[0], rtp_mag[1], rtp_mag[2]+dmag, atmosphere_dim ) - c1)/dmag;
          }

          if(not do_jac)
          {
              for( Index iv=0; iv<f_grid.nelem(); iv++ )
              {
                  const Numeric r = c1 / ( f_grid[iv] * f_grid[iv] );
                  propmat_clearsky[ife].SetFaraday(r, iv);
              }
          }
          else
          {
            for( Index iv=0; iv<f_grid.nelem(); iv++ )
            {
                const Numeric f2 = f_grid[iv] * f_grid[iv];
                const Numeric r = c1 / f2;
                propmat_clearsky[ife].SetFaraday(r, iv);
                
                // The Jacobian loop
                for(Index iq=0;iq<jacobian_quantities_position.nelem();iq++)
                {
                    if(is_frequency_parameter(jacobian_quantities[jacobian_quantities_position[iq]]))
                    {
                        dpropmat_clearsky_dx[iq].AddFaraday(-2.0 * r / f_grid[iv], iv);
                    }
                    else if(jacobian_quantities[jacobian_quantities_position[iq]]==JacPropMatType::MagneticU)
                    { 
                      dpropmat_clearsky_dx[iq].AddFaraday(dc1_u / f2, iv);
                    }
                    else if(jacobian_quantities[jacobian_quantities_position[iq]]==JacPropMatType::MagneticV)
                    { 
                      dpropmat_clearsky_dx[iq].AddFaraday(dc1_v / f2, iv);
                    }
                    else if(jacobian_quantities[jacobian_quantities_position[iq]]==JacPropMatType::MagneticW)
                    { 
                      dpropmat_clearsky_dx[iq].AddFaraday(dc1_w / f2, iv);
                    }
                    else if(jacobian_quantities[jacobian_quantities_position[iq]]==JacPropMatType::Electrons)
                    {
                      dpropmat_clearsky_dx[iq].AddFaraday(r/ne, iv);
                    }
                }
            }
          }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddParticles(
                                    // WS Output:
                                    ArrayOfPropagationMatrix& propmat_clearsky,
                                    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                                    // WS Input:
                                    const Index& stokes_dim,
                                    const Index& atmosphere_dim,
                                    const Vector& f_grid,
                                    const ArrayOfArrayOfSpeciesTag& abs_species,
                                    const ArrayOfRetrievalQuantity& jacobian_quantities,
                                    const Vector& rtp_vmr,
                                    const Vector& rtp_los,
                                    const Numeric& rtp_temperature,
                                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                                    const Index& scat_data_checked,
                                    const Index& use_abs_as_ext,
                                    // Verbosity object:
                                    const Verbosity& verbosity)
{
  CREATE_OUT2;

  // (i)yCalc only checks scat_data_checked if cloudbox is on. It is off here,
  // though, i.e. we need to check it here explicitly. (Also, cloudboxOff sets
  // scat_data_checked=0 as it does not check it and as we ususally don't need
  // scat_data for clearsky cases, hence don't want to check them by
  // scat_data_checkedCalc in that case. This approach seems to be the more
  // handy compared to cloudboxOff setting scat_data_checked=1 without checking
  // it assuming we won't use it anyways.)
  if( scat_data_checked != 1 )
    throw runtime_error( "The scat_data must be flagged to have "
                         "passed a consistency check (scat_data_checked=1)." );

  const Index ns = TotalNumberOfElements(scat_data);
  Index np = 0;
  for( Index sp = 0; sp < abs_species.nelem(); sp++ )
    {
      if (abs_species[sp][0].Type() == SpeciesTag::TYPE_PARTICLES)
        {
          np++;
        }
    }

  if( np == 0 )
    {
       ostringstream os; 
       os << "For applying propmat_clearskyAddParticles, *abs_species* needs to"
          << "contain species 'particles', but it does not.\n";
       throw runtime_error( os.str() );
    }

  if ( ns != np )
    {
      ostringstream os; 
      os << "Number of 'particles' entries in abs_species and of elements in\n"
         << "*scat_data* needs to be identical. But you have " << np
         << " 'particles' entries\n"
         << "and " << ns << " *scat_data* elements.\n";
      throw runtime_error( os.str() );
    }

  if( atmosphere_dim==1 && rtp_los.nelem() < 1 )
    {
       ostringstream os; 
       os << "For applying *propmat_clearskyAddParticles*, *rtp_los* needs to be specified\n"
          << "(at least zenith angle component for atmosphere_dim==1),\n"
          << "but it is not.\n";
       throw runtime_error( os.str() );
    }
  else if(  atmosphere_dim>1 && rtp_los.nelem() < 2 )
    {
       ostringstream os; 
       os << "For applying *propmat_clearskyAddParticles*, *rtp_los* needs to be specified\n"
          << "(both zenith and azimuth angle components for atmosphere_dim>1),\n"
          << "but it is not.\n";
       throw runtime_error( os.str() );
    }
    
    
  // Use for rescaling vmr of particulates
  Numeric rtp_vmr_sum = 0.0;
  
  // Tests and setup partial derivatives
  const bool do_jac = supports_particles(jacobian_quantities);
  const bool do_jac_temperature = do_temperature_jacobian(jacobian_quantities);
  const bool do_jac_frequencies = do_frequency_jacobian(jacobian_quantities);
  const ArrayOfIndex jacobian_quantities_position = equivlent_propmattype_indexes(jacobian_quantities);
  const Numeric dT = temperature_perturbation(jacobian_quantities);
  const Numeric df = frequency_perturbation(jacobian_quantities);
  bool first_loop = true;
      
  const Index nf = f_grid.nelem();
  const Index na = abs_species.nelem();
  Vector rtp_los_back;
  mirror_los( rtp_los_back, rtp_los, atmosphere_dim );
  ArrayOfArrayOfSingleScatteringData scat_data_mono, scat_data_mono_df;
  Matrix pnd_ext_mat(stokes_dim,stokes_dim), pnd_ext_mat_df(stokes_dim,stokes_dim), pnd_ext_mat_dt(stokes_dim,stokes_dim),
         tmp(stokes_dim, stokes_dim), tmp2(stokes_dim, stokes_dim);
  Vector pnd_abs_vec(stokes_dim), pnd_abs_vec_df(stokes_dim), pnd_abs_vec_dt(stokes_dim);
  Vector f_grid_df;
  
  // 170918 JM: along with transition to use of new-type (aka
  // pre-f_grid-interpolated) scat_data, freq perturbation switched off. Typical
  // clear-sky freq perturbations yield insignificant effects in particle
  // properties. Hence, the effort to make scat_data_monoCalc safely usable with
  // new-type scat_data is too high, and therefore this feature is neglected
  // here.
  if(do_jac_frequencies)
  {
    out2 << "WARNING:\n"
         << "Frequency perturbation not available for absorbing particles.\n";
    //f_grid_df = f_grid;
    //f_grid_df += ppd.Frequency_Perturbation();
  }
  
  for( Index iv=0; iv<nf; iv++ )
    { 
      // first, get the scat_data_single at the required frequency. we can do that for
      // all scattering elements at once.
      scat_data_monoExtract( scat_data_mono, scat_data, iv, verbosity );
      // 170918 JM: switched off (see above)
      //if(ppd.do_frequency())
      //  scat_data_monoCalc( scat_data_mono_df, scat_data, f_grid_df, iv, 
      //                      verbosity );

      // then we loop over the scat_data and link them with correct vmr_field
      // entry according to the position of the particle type entries in
      // abs_species.
      Index sp = 0;
      for( Index i_ss=0; i_ss<scat_data.nelem(); i_ss++ )
      {
        for( Index i_se=0; i_se<scat_data[i_ss].nelem(); i_se++ )
        {
          // forward to next particle entry in abs_species
          while ( sp < na && 
            abs_species[sp][0].Type() != SpeciesTag::TYPE_PARTICLES )
            sp++;
          
          // running beyond number of abs_species entries when looking for
          // next particle entry. shouldn't happen, though.
          assert ( sp < na );
          
          if ( rtp_vmr[sp] > 0. )
          {
            // get extinction matrix and absorption vector at
            // required temperature and direction for the individual
            // scattering element and multiply with their occurence.
            opt_propExtract(pnd_ext_mat, pnd_abs_vec,
                            scat_data_mono[i_ss][i_se],
                            rtp_los_back[0], rtp_los_back[1],
                            rtp_temperature, stokes_dim, verbosity);
            //cout << "absvec[0] = " << pnd_abs_vec[0] << "\n";
            pnd_ext_mat *= rtp_vmr[sp];
            pnd_abs_vec *= rtp_vmr[sp];
            
            // For wind derivatives
            /*
            if(ppd.do_frequency())
            {
              opt_propExtract(pnd_ext_mat_df, pnd_abs_vec_df,
                              scat_data_mono_df[i_ss][i_se],
                              rtp_los_back[0], rtp_los_back[1],
                              rtp_temperature, stokes_dim, verbosity);
              pnd_ext_mat_df *= rtp_vmr[sp];
              pnd_abs_vec_df *= rtp_vmr[sp];
            }
            */
            
            // For temperature derivatives
            if(do_jac_temperature)
            {
              opt_propExtract(pnd_ext_mat_dt, pnd_abs_vec_dt,
                              scat_data_mono[i_ss][i_se],
                              rtp_los_back[0], rtp_los_back[1],
                              rtp_temperature+dT, 
                              stokes_dim, verbosity);
              pnd_ext_mat_dt *= rtp_vmr[sp];
              pnd_abs_vec_dt *= rtp_vmr[sp];
            }
            
            // For number density derivatives
            if(do_jac && first_loop)
              rtp_vmr_sum += rtp_vmr[sp];
            
            if (use_abs_as_ext)
            {
              propmat_clearsky[sp].AddAbsorptionVectorAtPosition(pnd_abs_vec, iv);
            }
            else
            {
              propmat_clearsky[sp].SetAtPosition(pnd_ext_mat, iv);
              
            }
            
            for(Index iq=0; iq<jacobian_quantities_position.nelem(); iq++)
            {     
              if(is_frequency_parameter(jacobian_quantities[jacobian_quantities_position[iq]]))
              {
                if(use_abs_as_ext)
                {
                  ext_matFromabs_vec(tmp, pnd_abs_vec_df, stokes_dim);
                }
                else
                  tmp = pnd_ext_mat_df;
                
                propmat_clearsky[sp].MatrixAtPosition(tmp2, iv);
                tmp -= tmp2;
                tmp /= df;
                dpropmat_clearsky_dx[iq].AddAtPosition(tmp, iv);
              }
              else if(jacobian_quantities[jacobian_quantities_position[iq]] == JacPropMatType::Temperature)
              {
                if(use_abs_as_ext)
                {
                  ext_matFromabs_vec(tmp, pnd_abs_vec_dt, stokes_dim);
                }
                else
                  tmp = pnd_ext_mat_dt;
                propmat_clearsky[sp].MatrixAtPosition(tmp2, iv);
                tmp -= tmp2;
                tmp /= dT;
                dpropmat_clearsky_dx[iq].AddAtPosition(tmp, iv);
              }
              else if(jacobian_quantities[jacobian_quantities_position[iq]] == JacPropMatType::Particulates)
              {
                dpropmat_clearsky_dx[iq].AddAtPosition(propmat_clearsky[sp], iv);
              }
            }
          }
          sp++;
        }
      }
      //checking that no further 'particle' entry left after all scat_data
      //entries are processes. this is basically not necessary. but checking it
      //anyway to really be safe. remove later, when more extensively tested.
      while (sp < na)
      {
        assert ( abs_species[sp][0].Type() != SpeciesTag::TYPE_PARTICLES );
        sp++;
      }
        
      // Not first loop anymore
      first_loop = false;
    }
    
    if(rtp_vmr_sum != 0.0)
    {
      for(Index iq=0; iq<jacobian_quantities_position.nelem(); iq++)
      {
        if(jacobian_quantities[jacobian_quantities_position[iq]] == JacPropMatType::Particulates)
        {
          dpropmat_clearsky_dx[iq] /= rtp_vmr_sum;
        }
      }
    };
}


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddParticles2(
                                    // WS Output:
                                    ArrayOfPropagationMatrix& propmat_clearsky,
                                    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                                    // WS Input:
                                    const Index& stokes_dim,
                                    const Index& atmosphere_dim,
                                    const Vector& f_grid,
                                    const ArrayOfArrayOfSpeciesTag& abs_species,
                                    const ArrayOfRetrievalQuantity& jacobian_quantities,
                                    const Vector& rtp_vmr,
                                    const Vector& rtp_los,
                                    const Numeric& rtp_temperature,
                                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                                    const Index& scat_data_checked,
                                    const Index& use_abs_as_ext,
                                    // Verbosity object:
                                    const Verbosity& verbosity)
{
  CREATE_OUT1;

  // (i)yCalc only checks scat_data_checked if cloudbox is on. It is off here,
  // though, i.e. we need to check it here explicitly. (Also, cloudboxOff sets
  // scat_data_checked=0 as it does not check it and as we ususally don't need
  // scat_data for clearsky cases, hence don't want to check them by
  // scat_data_checkedCalc in that case. This approach seems to be the more
  // handy compared to cloudboxOff setting scat_data_checked=1 without checking
  // it assuming we won't use it anyways.)
  if( scat_data_checked != 1 )
    throw runtime_error( "The scat_data must be flagged to have "
                         "passed a consistency check (scat_data_checked=1)." );

  const Index ns = TotalNumberOfElements(scat_data);
  Index np = 0;
  for( Index sp = 0; sp < abs_species.nelem(); sp++ )
    {
      if (abs_species[sp][0].Type() == SpeciesTag::TYPE_PARTICLES)
        {
          np++;
        }
    }

  if( np == 0 )
    {
       ostringstream os; 
       os << "For applying propmat_clearskyAddParticles, *abs_species* needs to"
          << "contain species 'particles', but it does not.\n";
       throw runtime_error( os.str() );
    }

  if ( ns != np )
    {
      ostringstream os; 
      os << "Number of 'particles' entries in abs_species and of elements in\n"
         << "*scat_data* needs to be identical. But you have " << np
         << " 'particles' entries\n"
         << "and " << ns << " *scat_data* elements.\n";
      throw runtime_error( os.str() );
    }

  if( atmosphere_dim==1 && rtp_los.nelem() < 1 )
    {
       ostringstream os; 
       os << "For applying *propmat_clearskyAddParticles*, *rtp_los* needs to be specified\n"
          << "(at least zenith angle component for atmosphere_dim==1),\n"
          << "but it is not.\n";
       throw runtime_error( os.str() );
    }
  else if(  atmosphere_dim>1 && rtp_los.nelem() < 2 )
    {
       ostringstream os; 
       os << "For applying *propmat_clearskyAddParticles*, *rtp_los* needs to be specified\n"
          << "(both zenith and azimuth angle components for atmosphere_dim>1),\n"
          << "but it is not.\n";
       throw runtime_error( os.str() );
    }
    
    
  // Use for rescaling vmr of particulates
  Numeric rtp_vmr_sum = 0.0;
  
  // Tests and setup partial derivatives
  const bool do_jac = supports_particles(jacobian_quantities);
  const bool do_jac_temperature = do_temperature_jacobian(jacobian_quantities);
  const bool do_jac_frequencies = do_frequency_jacobian(jacobian_quantities);
  const ArrayOfIndex jacobian_quantities_position = equivlent_propmattype_indexes(jacobian_quantities);
  const Numeric dT = temperature_perturbation(jacobian_quantities);
      
  const Index na = abs_species.nelem();
  Vector rtp_los_back;
  mirror_los( rtp_los_back, rtp_los, atmosphere_dim );
  
  // 170918 JM: along with transition to use of new-type (aka
  // pre-f_grid-interpolated) scat_data, freq perturbation switched off. Typical
  // clear-sky freq perturbations yield insignificant effects in particle
  // properties. Hence, this feature is neglected here.
  if(do_jac_frequencies)
  {
    out1 << "WARNING:\n"
         << "Frequency perturbation not available for absorbing particles.\n";
  }
  
  // creating temporary output containers
  ArrayOfArrayOfTensor5 ext_mat_Nse;
  ArrayOfArrayOfTensor4 abs_vec_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;

  // preparing input in format needed
  Vector T_array;
  if(do_jac_temperature)
  {
    T_array.resize(2);
    T_array = rtp_temperature;
    T_array[1] += dT;
  }
  else
  {
    T_array.resize(1);
    T_array = rtp_temperature;
  }
  Matrix dir_array(1,2);
  dir_array(0,joker) = rtp_los_back;

  // ext/abs per scat element for all freqs at once
  opt_prop_NScatElems( ext_mat_Nse, abs_vec_Nse, ptypes_Nse, t_ok,
                       scat_data, stokes_dim, T_array, dir_array, -1 );

  const Index nf = abs_vec_Nse[0][0].nbooks();
  Tensor3 tmp(nf,stokes_dim,stokes_dim);

  // loop over the scat_data and link them with correct vmr_field entry according
 // to the position of the particle type entries in abs_species.
 Index sp = 0;
 Index i_se_flat = 0;
 for( Index i_ss=0; i_ss<scat_data.nelem(); i_ss++ )
 {
  for( Index i_se=0; i_se<scat_data[i_ss].nelem(); i_se++ )
  {
    // forward to next particle entry in abs_species
    while ( sp < na && 
            abs_species[sp][0].Type() != SpeciesTag::TYPE_PARTICLES )
      sp++;
          
    // running beyond number of abs_species entries when looking for next
    // particle entry. shouldn't happen, though.
    assert ( sp < na );
    if ( rtp_vmr[sp] < 0. )
    {
      ostringstream os;
      os << "Negative absorbing particle 'vmr' (aka number density)"
         << " encountered:\n"
         << "scat species #" << i_ss << ", scat elem #" << i_se
         << " (vmr_field entry #" << sp << ")\n";
      throw runtime_error( os.str() );
    }
    else if ( rtp_vmr[sp] > 0. )
    {
      if( t_ok(i_se_flat,0) < 0. )
      {
        ostringstream os;
        os << "Temperature interpolation error:\n"
           << "scat species #" << i_ss << ", scat elem #" << i_se << "\n";
        throw runtime_error( os.str() );
      }
      else
      {        
        if( use_abs_as_ext )
        {
          if( nf > 1 )
            for( Index iv=0; iv<f_grid.nelem(); iv++ )
              propmat_clearsky[sp].AddAbsorptionVectorAtPosition(
                abs_vec_Nse[i_ss][i_se](iv,0,0,joker), iv);
          else
            for( Index iv=0; iv<f_grid.nelem(); iv++ )
              propmat_clearsky[sp].AddAbsorptionVectorAtPosition(
                abs_vec_Nse[i_ss][i_se](0,0,0,joker), iv);
        }
        else
        {
          if( nf > 1 )
            for( Index iv=0; iv<f_grid.nelem(); iv++ )
              propmat_clearsky[sp].SetAtPosition(
                ext_mat_Nse[i_ss][i_se](iv,0,0,joker,joker), iv);
          else
            for( Index iv=0; iv<f_grid.nelem(); iv++ )
              propmat_clearsky[sp].SetAtPosition(
                ext_mat_Nse[i_ss][i_se](0,0,0,joker,joker), iv);
        }
        propmat_clearsky[sp] *= rtp_vmr[sp];
      }

      // For temperature derivatives (so we don't need to check it in jac loop)
      if( do_jac_temperature )
      {
        if( t_ok(i_se_flat,1) < 0. )
        {
          ostringstream os;
          os << "Temperature interpolation error (in perturbation):\n"
             << "scat species #" << i_ss << ", scat elem #" << i_se << "\n";
          throw runtime_error( os.str() );
        }
      }

      // For number density derivatives
      if( do_jac )
        rtp_vmr_sum += rtp_vmr[sp];
            
      for( Index iq=0; iq<jacobian_quantities_position.nelem(); iq++ )
      {
        // we don't do freq perturbations here, i.e. nothing to do.
        /*
        if( ppd(iq)==JQT_frequency || ppd(iq)==JQT_wind_magnitude || 
            ppd(iq)==JQT_wind_u || ppd(iq)==JQT_wind_v || ppd(iq)==JQT_wind_w )
        {
          dpropmat_clearsky_dx[iq] = 0.;
        }

        else if( ppd(iq) == JQT_temperature )
        */
        if( jacobian_quantities[jacobian_quantities_position[iq]] == JacPropMatType::Temperature )
        {
          if( use_abs_as_ext )
          {
            tmp(joker,joker,0) = abs_vec_Nse[i_ss][i_se](joker,1,0,joker);
            tmp(joker,joker,0) -= abs_vec_Nse[i_ss][i_se](joker,0,0,joker);
          }
          else
          {
            tmp = ext_mat_Nse[i_ss][i_se](joker,1,0,joker,joker);
            tmp -= ext_mat_Nse[i_ss][i_se](joker,0,0,joker,joker);
          }

          tmp *= rtp_vmr[sp];
          tmp /= dT;

          if( nf > 1 )
            for( Index iv=0; iv<f_grid.nelem(); iv++ )
              if ( use_abs_as_ext )
                dpropmat_clearsky_dx[iq].AddAbsorptionVectorAtPosition(
                  tmp(iv,joker,0), iv);
              else
                dpropmat_clearsky_dx[iq].AddAtPosition(
                  tmp(iv,joker,joker), iv);
          else
            for( Index iv=0; iv<f_grid.nelem(); iv++ )
              if ( use_abs_as_ext )
                dpropmat_clearsky_dx[iq].AddAbsorptionVectorAtPosition(
                  tmp(0,joker,0), iv);
              else
                dpropmat_clearsky_dx[iq].AddAtPosition(
                  tmp(0,joker,joker), iv);
        }
        /* // alternative version
        else if( ppd(iq) == JQT_temperature )
        {
          if( nf > 1 )
            for( Index iv=0; iv<f_grid.nelem(); iv++ )
              if( use_abs_as_ext )
              {
                atmp = abs_vec_Nse[i_ss][i_se](iv,1,0,joker);
                atmp -= abs_vec_Nse[i_ss][i_se](iv,0,0,joker);
                atmp *= (rtp_vmr[sp] / ppd.Temperature_Perturbation());
                dpropmat_clearsky_dx[iq].AddAbsorptionVectorAtPosition(atmp, iv);
              }
              else
              {
                etmp = ext_mat_Nse[i_ss][i_se](iv,1,0,joker,joker) -
                      ext_mat_Nse[i_ss][i_se](iv,0,0,joker,joker);
                etmp *= (rtp_vmr[sp] / ppd.Temperature_Perturbation());
                dpropmat_clearsky_dx[iq].AddAtPosition(etmp, iv);
              }
          else
            if( use_abs_as_ext )
            {
              atmp = (abs_vec_Nse[i_ss][i_se](0,1,0,joker) -
                      abs_vec_Nse[i_ss][i_se](0,0,0,joker));
              atmp *= (rtp_vmr[sp] / ppd.Temperature_Perturbation());
              for( Index iv=0; iv<f_grid.nelem(); iv++ )
                dpropmat_clearsky_dx[iq].AddAbsorptionVectorAtPosition(atmp, iv);
            }
            else
            {
              etmp = ext_mat_Nse[i_ss][i_se](0,1,0,joker,joker) -
                    ext_mat_Nse[i_ss][i_se](0,0,0,joker,joker);
              etmp *= (rtp_vmr[sp] / ppd.Temperature_Perturbation());
              for( Index iv=0; iv<f_grid.nelem(); iv++ )
                dpropmat_clearsky_dx[iq].AddAtPosition(etmp, iv);
            }
        }
        */

        else if(jacobian_quantities[jacobian_quantities_position[iq]] == JacPropMatType::Particulates )
        {
          for( Index iv=0; iv<f_grid.nelem(); iv++ )
            dpropmat_clearsky_dx[iq].AddAtPosition(propmat_clearsky[sp], iv);
        }
      }
    }
    sp++;
    i_se_flat++;
  }
 }

 //checking that no further 'particle' entry left after all scat_data entries
 //are processes. this is basically not necessary. but checking it anyway to
 //really be safe. remove later, when more extensively tested.
 while (sp < na)
 {
    assert ( abs_species[sp][0].Type() != SpeciesTag::TYPE_PARTICLES );
    sp++;
 }

 if(rtp_vmr_sum != 0.0)
 {
   for(Index iq=0; iq<jacobian_quantities_position.nelem(); iq++)
    {
      if(jacobian_quantities[jacobian_quantities_position[iq]] == JacPropMatType::Particulates)
      {
        dpropmat_clearsky_dx[iq] /= rtp_vmr_sum;
      }
    }
 }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddOnTheFly(// Workspace reference:
                                    Workspace& ws,
                                    // WS Output:
                                    ArrayOfPropagationMatrix& propmat_clearsky,
                                    ArrayOfStokesVector& nlte_source,
                                    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                                    ArrayOfStokesVector& dnlte_dx_source,
                                    ArrayOfStokesVector& nlte_dsource_dx,
                                    // WS Input:
                                    const Vector& f_grid,
                                    const ArrayOfArrayOfSpeciesTag& abs_species,
                                    const ArrayOfRetrievalQuantity& jacobian_quantities,
                                    const Numeric& rtp_pressure,
                                    const Numeric& rtp_temperature,
                                    const Vector& rtp_nlte,
                                    const Vector& rtp_vmr,
                                    const Agenda& abs_xsec_agenda,
                                    // Verbosity object:
                                    const Verbosity& verbosity)
{
  CREATE_OUT3;

  // Define communication variables for the actual absorption calculation:
  
  // Output of AbsInputFromRteScalars:
  Vector        abs_p;
  Vector        abs_t;
  Matrix        abs_nlte;
  Matrix        abs_vmrs;
  // Output of abs_h2oSet:
  Vector          abs_h2o;
  // Output of abs_coefCalc:
  Matrix                         abs_coef, src_coef;
  ArrayOfMatrix                  abs_coef_per_species, src_coef_per_species, dabs_coef_dx, dsrc_coef_dx;
    
  AbsInputFromRteScalars(abs_p,
                         abs_t,
                         abs_nlte,
                         abs_vmrs,
                         rtp_pressure,
                         rtp_temperature,
                         rtp_nlte,
                         rtp_vmr,
                         verbosity);
  
  // Absorption cross sections per tag group.
  ArrayOfMatrix abs_xsec_per_species;
  ArrayOfMatrix src_xsec_per_species;
  ArrayOfArrayOfMatrix dabs_xsec_per_species_dx, dsrc_xsec_per_species_dx;

  // Make all species active:
  ArrayOfIndex abs_species_active(abs_species.nelem());
  for (Index i=0; i<abs_species.nelem(); ++i)
    abs_species_active[i] = i;
  
  // Call agenda to calculate absorption:
  abs_xsec_agendaExecute(ws,
                         abs_xsec_per_species,
                         src_xsec_per_species,
                         dabs_xsec_per_species_dx,
                         dsrc_xsec_per_species_dx,
                         abs_species,
                         jacobian_quantities,
                         abs_species_active,
                         f_grid,
                         abs_p,
                         abs_t,
                         abs_nlte,
                         abs_vmrs,
                         abs_xsec_agenda);
  
  // Calculate absorption coefficients from cross sections:
  abs_coefCalcFromXsec(abs_coef, src_coef, dabs_coef_dx, dsrc_coef_dx, 
                       abs_coef_per_species,  src_coef_per_species, 
                       abs_xsec_per_species, src_xsec_per_species, 
                       dabs_xsec_per_species_dx, dsrc_xsec_per_species_dx, 
                       abs_species, jacobian_quantities,
                       abs_vmrs, abs_p, abs_t, verbosity);
  
  // Now add abs_coef_per_species to propmat_clearsky:
  propmat_clearskyAddFromAbsCoefPerSpecies(propmat_clearsky, dpropmat_clearsky_dx, abs_coef_per_species, dabs_coef_dx, verbosity);

  // Now turn nlte_source from absorption into a proper source function
  if(not nlte_source.empty())
      nlte_sourceFromTemperatureAndSrcCoefPerSpecies(nlte_source, dnlte_dx_source, nlte_dsource_dx, 
                                                     src_coef_per_species, dsrc_coef_dx, jacobian_quantities,
                                                     f_grid, rtp_temperature, verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyZero(
         ArrayOfPropagationMatrix&    propmat_clearsky,
   const Vector&     f_grid,
   const Index&      stokes_dim,
   const Verbosity& )
{
  propmat_clearsky.resize(1);
  propmat_clearsky[0] = PropagationMatrix(f_grid.nelem(), stokes_dim);
  propmat_clearsky[0].SetZero();
}


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyForceNegativeToZero(
    ArrayOfPropagationMatrix&    propmat_clearsky,
    const Verbosity& )
{
    for(auto& pm : propmat_clearsky)
      for(Index i = 0; i < pm.NumberOfFrequencies(); i++)
        if(pm.Kjj()[i] < 0.0)
          pm.SetAtPosition(0.0, i);
          
}


/* Workspace method: Doxygen documentation will be auto-generated */
void isotopologue_ratiosInitFromBuiltin(SpeciesAuxData& isotopologue_ratios,
                                        const Verbosity&)
{
    fillSpeciesAuxDataWithIsotopologueRatiosFromSpeciesData(isotopologue_ratios);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void partition_functionsInitFromBuiltin(SpeciesAuxData& partition_functions,
                                        const Verbosity&)
{
    fillSpeciesAuxDataWithPartitionFunctionsFromSpeciesData(partition_functions);
}


#ifdef ENABLE_NETCDF
/* Workspace method: Doxygen documentation will be auto-generated */
/* Included by Claudia Emde, 20100707 */
void WriteMolTau(//WS Input
                 const Vector& f_grid, 
                 const Tensor3& z_field,
                 const Tensor7& propmat_clearsky_field,
                 const Index& atmosphere_dim,
                 //Keyword
                 const String& filename,
                 const Verbosity&)
{
  
  int retval, ncid;
  int nlev_dimid, nlyr_dimid, nwvl_dimid, stokes_dimid, none_dimid;
  int dimids[4];
  int wvlmin_varid, wvlmax_varid, z_varid, wvl_varid, tau_varid;
  
  if (atmosphere_dim != 1)
    throw runtime_error("WriteMolTau can only be used for atmsophere_dim=1");

#pragma omp critical(netcdf__critical_region)
    {
  // Open file
  if ((retval = nc_create(filename.c_str(), NC_CLOBBER, &ncid)))
    nca_error (retval, "nc_create");
  
  // Define dimensions
  if ((retval = nc_def_dim(ncid, "nlev", (int) z_field.npages(), &nlev_dimid)))
    nca_error (retval, "nc_def_dim");
  
  if ((retval = nc_def_dim(ncid, "nlyr", (int) z_field.npages() - 1, &nlyr_dimid)))
    nca_error (retval, "nc_def_dim");
  
  if ((retval = nc_def_dim(ncid, "nwvl", (int) f_grid.nelem(), &nwvl_dimid)))
    nca_error (retval, "nc_def_dim");
  
  if ((retval = nc_def_dim(ncid, "none", 1, &none_dimid)))
    nca_error (retval, "nc_def_dim");

  if ((retval = nc_def_dim(ncid, "nstk", (int) propmat_clearsky_field.nbooks(), &stokes_dimid)))
    nca_error (retval, "nc_def_dim");

  // Define variables
  if ((retval = nc_def_var(ncid, "wvlmin", NC_DOUBLE, 1,  &none_dimid, &wvlmin_varid)))
    nca_error (retval, "nc_def_var wvlmin");
  
  if ((retval = nc_def_var(ncid, "wvlmax", NC_DOUBLE, 1,  &none_dimid, &wvlmax_varid)))
    nca_error (retval, "nc_def_var wvlmax");
  
  if ((retval = nc_def_var(ncid, "z", NC_DOUBLE, 1,  &nlev_dimid, &z_varid)))
    nca_error (retval, "nc_def_var z");
  
  if ((retval = nc_def_var(ncid, "wvl", NC_DOUBLE, 1,  &nwvl_dimid, &wvl_varid)))
    nca_error (retval, "nc_def_var wvl");
  
  dimids[0]=nlyr_dimid;
  dimids[1]=nwvl_dimid; 
  dimids[2]=stokes_dimid;
  dimids[3]=stokes_dimid;

  if ((retval = nc_def_var(ncid, "tau", NC_DOUBLE, 4, &dimids[0], &tau_varid)))
    nca_error (retval, "nc_def_var tau");
  
  // Units
  if ((retval = nc_put_att_text(ncid, wvlmin_varid, "units", 2, "nm")))
    nca_error (retval, "nc_put_att_text");

  if ((retval = nc_put_att_text(ncid, wvlmax_varid, "units", 2, "nm")))
    nca_error (retval, "nc_put_att_text");
  
  if ((retval = nc_put_att_text(ncid, z_varid, "units", 2, "km")))
    nca_error (retval, "nc_put_att_text");

  if ((retval = nc_put_att_text(ncid, wvl_varid, "units", 2, "nm")))
     nca_error (retval, "nc_put_att_text");
  
  if ((retval = nc_put_att_text(ncid, tau_varid, "units", 1, "-")))
     nca_error (retval, "nc_put_att_text");
  
  // End define mode. This tells netCDF we are done defining
  // metadata. 
  if ((retval = nc_enddef(ncid)))
    nca_error (retval, "nc_enddef");
  
  // Assign data
  double wvlmin[1];
  wvlmin[0]= SPEED_OF_LIGHT/f_grid[f_grid.nelem()-1]*1e9;
  if ((retval = nc_put_var_double (ncid, wvlmin_varid, &wvlmin[0])))
    nca_error (retval, "nc_put_var");

  double wvlmax[1];
  wvlmax[0]= SPEED_OF_LIGHT/f_grid[0]*1e9;
  if ((retval = nc_put_var_double (ncid, wvlmax_varid, &wvlmax[0])))
    nca_error (retval, "nc_put_var");
  
  double z[z_field.npages()];
  for (int iz=0; iz<z_field.npages(); iz++)
    z[iz]=z_field(z_field.npages()-1-iz, 0, 0)*1e-3;
  
  if ((retval = nc_put_var_double (ncid, z_varid, &z[0])))
    nca_error (retval, "nc_put_var");
  
  double wvl[f_grid.nelem()];
  for (int iv=0; iv<f_grid.nelem(); iv++)
    wvl[iv]=SPEED_OF_LIGHT/f_grid[f_grid.nelem()-1-iv]*1e9;
  
  if ((retval = nc_put_var_double (ncid, wvl_varid, &wvl[0])))
    nca_error (retval, "nc_put_var");

  const Index zfnp = z_field.npages()-1;
  const Index fgne = f_grid.nelem();
  const Index amfnb = propmat_clearsky_field.nbooks();

  Tensor4 tau(zfnp, fgne, amfnb, amfnb, 0.);

  // Calculate average tau for layers
  for (int is=0; is<propmat_clearsky_field.nlibraries(); is++)
    for (int iz=0; iz<zfnp; iz++)
      for (int iv=0; iv<fgne; iv++)
          for (int is1=0; is1<amfnb; is1++)
              for (int is2=0; is2<amfnb; is2++)
                // sum up all species
                tau(iz, iv, is1, is2) += 0.5 * (propmat_clearsky_field(is,f_grid.nelem()-1-iv,is1,is2,z_field.npages()-1-iz,0,0)+
                                    propmat_clearsky_field(is,f_grid.nelem()-1-iv,is1,is2,z_field.npages()-2-iz,0,0))
                *(z_field(z_field.npages()-1-iz,0,0)-z_field(z_field.npages()-2-iz,0,0));

  
  if ((retval = nc_put_var_double (ncid, tau_varid, tau.get_c_array())))
    nca_error (retval, "nc_put_var");
  
  // Close the file
  if ((retval = nc_close(ncid)))
    nca_error (retval, "nc_close");
    }
}

#else

void WriteMolTau(//WS Input
                 const Vector& f_grid _U_,
                 const Tensor3& z_field _U_,
                 const Tensor7& propmat_clearsky_field _U_,
                 const Index& atmosphere_dim _U_,
                 //Keyword
                 const String& filename _U_,
                 const Verbosity&)
{
  throw runtime_error("The workspace method WriteMolTau is not available"
                      "because ARTS was compiled without NetCDF support.");
}

#endif /* ENABLE_NETCDF */

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddLines2(// WS Output:
                                   ArrayOfMatrix& abs_xsec_per_species,
                                   ArrayOfMatrix& src_xsec_per_species,
                                   ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx,
                                   ArrayOfArrayOfMatrix& dsrc_xsec_per_species_dx,
                                   // WS Input:             
                                   const ArrayOfArrayOfSpeciesTag& tgs,
                                   const ArrayOfRetrievalQuantity& jacobian_quantities,
                                   const ArrayOfIndex& abs_species_active,
                                   const Vector& f_grid,
                                   const Vector& abs_p,
                                   const Vector& abs_t,
                                   const Matrix& abs_nlte,
                                   const Numeric& lm_p_lim,
                                   const Index& xsec_speedup_switch,
                                   const Matrix& abs_vmrs,
                                   const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                   const SpeciesAuxData& isotopologue_ratios,
                                   const SpeciesAuxData& partition_functions,
                                   const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Check that correct isotopologue ratios are defined for the species
  // we want to calculate
  checkIsotopologueRatios(tgs, isotopologue_ratios);
  
  // Check that all temperatures are at least 0 K. (Negative Kelvin
  // temperatures are unphysical.)  
  if ( min(abs_t) < 0 )
  {
    ostringstream os;
    os << "Temperature must be at least 0 K. But you request an absorption\n"
    << "calculation at " << min(abs_t) << " K!"; 
    throw std::runtime_error(os.str());
  }
  
  // Check that all parameters that should have the number of tag
  // groups as a dimension are consistent.
  {
    const Index n_tgs    = tgs.nelem();
    const Index n_xsec   = abs_xsec_per_species.nelem();
    const Index n_vmrs   = abs_vmrs.nrows();
    const Index n_lines  = abs_lines_per_species.nelem();
    
    if ( n_tgs not_eq n_xsec or n_tgs not_eq n_vmrs or n_tgs not_eq n_lines)
    {
      ostringstream os;
      os << "The following variables must all have the same dimension:\n"
      << "tgs:          " << tgs.nelem() << "\n"
      << "abs_xsec_per_species:  " << abs_xsec_per_species.nelem() << "\n"
      << "abs_vmrs:         " << abs_vmrs.nrows() << "\n"
      << "abs_lines_per_species: " << abs_lines_per_species.nelem() << "\n"
      << "(As a special case, abs_lineshape is allowed to have only one element.)";
      throw std::runtime_error(os.str());
    }
  }
  
  // If xsec_speedup_switch is true the length of f_grid is 2^N+1, find N:
  Index binary_speedup = 0;
  if(xsec_speedup_switch)
  {
    Index n = f_grid.nelem() - 1;
    if(n < 4)
      throw std::runtime_error("Requesting speedup requires atleast 5 frequency points in total.");
    do
    {
      if(n & 0x01)
        throw std::runtime_error("Requesting speedup requires 2^N+1 frequency points");
      n = n >> 1;
      binary_speedup += 1;
    } while(n not_eq 1);
  }
  
  const bool do_jac = supports_propmat_clearsky(jacobian_quantities);
  const bool do_lte = abs_nlte.empty();
  const ArrayOfIndex jac_pos = equivlent_propmattype_indexes(jacobian_quantities);
  
  // Skipping uninteresting data
  static ArrayOfMatrix dummy(0);
  
  // Print information:
  //
  out3 << "  Calculating line spectra.\n";
  
  // Call xsec_species for each tag group.
  for ( Index ii=0; ii<abs_species_active.nelem(); ++ii )
  {
    const Index i = abs_species_active[ii];
    
    // Get a pointer to the line list for the current species. This
    // is just so that we don't have to type abs_lines_per_species[i] over
    // and over again.
    const ArrayOfLineRecord& ll = abs_lines_per_species[i];
    
    // We do the LBL calculation only if:
    // - The line list is not empty, and
    // - The species is not a Zeeman species.
    if ( 0 < ll.nelem() and tgs[i].nelem() and not is_zeeman(tgs[i]) )
    {
      // As a safety check, check that the species of the first
      // line matches the species we should have according to
      // tgs. (This in case the order in tgs has been changed and
      // abs_lines_per_species has not been changed consistently.)
      if (ll[0].Species() != tgs[i][0].Species() )
      {
        ostringstream os;
        os << "The species in the line list does not match the species\n"
        << "for which you want to calculate absorption:\n"
        << "abs_species:           " << get_tag_group_name(tgs[i]) << "\n"
        << "abs_lines_per_species: " << ll[0].Name();
        throw std::runtime_error(os.str());
      }
      
      xsec_species2(abs_xsec_per_species[i],
                    src_xsec_per_species[i],
                    Matrix(0, 0),
                    do_jac?dabs_xsec_per_species_dx[i]:dummy,
                    (do_jac and not do_lte)?dsrc_xsec_per_species_dx[i]:dummy,
                    dummy,
                    jacobian_quantities,
                    jac_pos,
                    f_grid,
                    abs_p,
                    abs_t,
                    abs_nlte,
                    abs_vmrs,
                    tgs,
                    i,
                    ll,
                    0.0,
                    lm_p_lim,
                    isotopologue_ratios,
                    partition_functions,
                    binary_speedup, 
                    verbosity);
    }
    
    if (out3.sufficient_priority())
    {
      ostringstream os;
      os << "  Tag group " << i
      << " (" << get_tag_group_name(tgs[i]) << "): "
      << ll.nelem() << " transitions\n";
      out3 << os.str();
    }
    
  } // End of species for loop.
}

