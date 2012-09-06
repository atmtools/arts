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
#include "absorption.h"
#include "auto_md.h"
#include "math_funcs.h"
#include "make_array.h"
#include "physics_funcs.h"
#include "continua.h"
#include "make_vector.h"
#include "check_input.h"
#include "xml_io.h"
#include "parameters.h"

#ifdef ENABLE_NETCDF
#include <netcdf.h>
#include "nc_io.h"
#endif

extern const Numeric SPEED_OF_LIGHT;

/* Workspace method: Doxygen documentation will be auto-generated */
void AbsInputFromRteScalars(// WS Output:
                            Vector&        abs_p,
                            Vector&        abs_t,
                            Matrix&        abs_vmrs,
                            // WS Input:
                            const Numeric& rte_pressure,
                            const Numeric& rte_temperature,
                            const Vector&  rte_vmr_list,
                            const Verbosity&)
{
  // Prepare abs_p:
  abs_p.resize(1);
  abs_p = rte_pressure;

  // Prepare abs_t:
  abs_t.resize(1);
  abs_t = rte_temperature;

  // Prepare abs_vmrs:
  abs_vmrs.resize(rte_vmr_list.nelem(),1);
  abs_vmrs = rte_vmr_list;
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
    // Loop over all lines, use member function to do conversion.
#pragma omp parallel for    \
if(!arts_omp_in_parallel())
    for ( Index i=0; i<abs_lines.nelem(); ++i )
        abs_lines[i].ARTSCAT4FromARTSCAT3();
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
      if ( lr.ReadFromHitranStream(is, verbosity) )
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
                abs_lines.push_back(lr);
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

  out2 << "  Reading file: " << filename << "\n";
  open_input_file(is, filename);

  bool go_on = true;
  while ( go_on )
    {
      LineRecord lr;
      if ( lr.ReadFromHitran2004Stream(is, verbosity) )
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
                abs_lines.push_back(lr);
              else
                go_on = false;
            }
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
              abs_lines.push_back(lr);
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
                abs_lines.push_back(lr);
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
  extern const Array<SpeciesRecord> species_data;
  
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
      WriteXML(output_file_format, *it, species_filename, "", "", verbosity);
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
  extern const Array<SpeciesRecord> species_data;
 
  // Build a set of species indices. Duplicates are ignored.
  set<Index> unique_species;
  for (ArrayOfArrayOfSpeciesTag::const_iterator asp = abs_species.begin();
       asp != abs_species.end(); asp++)
    for (ArrayOfSpeciesTag::const_iterator sp = asp->begin();
         sp != asp->end(); sp++)
    {
      // If the Isotopologue number is equal to the number of Isotopologues for that species,
      // it means 'all isotopologues', in which case we have to read a catalog.
      // Continua don't need a catalog.
      if (sp->Isotopologue() == species_data[sp->Species()].Isotopologue().nelem()
          || !species_data[sp->Species()].Isotopologue()[sp->Isotopologue()].isContinuum())
          {
            unique_species.insert(sp->Species());
          }
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
              abs_lines.push_back(lr);
            }
        }
    }
  out2 << "  Read " << abs_lines.nelem() << " lines.\n";
}


/* FIXME: OLE: Do we need this function? */
void linesElowToJoule(// WS Output:
                      ArrayOfLineRecord& abs_lines )
{
  for ( Index i=0; i<abs_lines.nelem(); ++i )
    abs_lines[i].melow = wavenumber_to_joule(abs_lines[i].melow); 
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

  MakeArray<String>      real_filenames ( filenames[0] );
  MakeArray<String>      real_formats   ( formats[0]   );
  MakeArray<Numeric>     real_fmin      ( fmin[0]      );
  MakeArray<Numeric>     real_fmax      ( fmax[0]      );

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
              real_tgs.push_back( MakeArray<Index>(i) );

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
          abs_lines_per_species[real_tgs[i][s]] = these_abs_lines_per_species[s];
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
  extern const Array<SpeciesRecord> species_data;

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
                   << "why do you not include line "
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


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesAddMirrorLines(// WS Output:
                                         ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                         const Verbosity&)
{
  // We will simply append the mirror lines after the original
  // lines. This way we don't have to make a backup copy of the
  // original lines. 

  for ( Index i=0; i<abs_lines_per_species.nelem(); ++i )
    {
      // Get a reference to the current list of lines to save typing:
      ArrayOfLineRecord& ll = abs_lines_per_species[i];
      
      // For increased efficiency, reserve the necessary space:
      ll.reserve(2*ll.nelem());

      // Loop through all lines of this tag group:
      {
        // It is important that we determine the size of ll *before*
        // we start the loop. After all, we are adding elements. And
        // we cerainly don't want to continue looping the newly added
        // elements, we want to loop only the original elements.
        Index n=ll.nelem();
        for ( Index j=0; j<n; ++j )
          {
            LineRecord dummy = ll[j];
            dummy.setF( -dummy.F() );
            //      cout << "Adding ML at f = " << dummy.F() << "\n";
            ll.push_back(dummy);
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
                  nll.push_back (**j);
                }
              // Overwrite the old array with the new one
              ll.resize (nll.nelem ());
              ll = nll;
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesDefineAllInScenario(// WS Output:
                                    ArrayOfArrayOfSpeciesTag& tgs,
                                    // Control Parameters:
                                    const String& basename,
                                    const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  // Species lookup data:
  extern const Array<SpeciesRecord> species_data;

  // We want to make lists of included and excluded species:
  ArrayOfString included(0), excluded(0);
  bool found_file;

  // Command line parameters which give us the include search path.
  extern const Parameters parameters;
  ArrayOfString allpaths = parameters.includepath;
  allpaths.insert(allpaths.end(),
                  parameters.datapath.begin(),
                  parameters.datapath.end());

  tgs.resize(0);

  for ( Index i=0; i<species_data.nelem(); ++i )
    {
      const String specname = species_data[i].Name();
      String filename = basename + "." + specname;

      found_file = find_file(filename, ".xml", allpaths);
      if (!found_file) found_file = find_file(filename, ".xml.gz", allpaths);
      if (!found_file) found_file = find_file(filename, ".gz", allpaths);

      if (found_file)
        {
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
      else
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
  extern const Array<LineshapeRecord> lineshape_data;
  extern const Array<LineshapeNormRecord> lineshape_norm_data;


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
  extern const Array<LineshapeRecord> lineshape_data;
  extern const Array<LineshapeNormRecord> lineshape_norm_data;

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


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_h2oSet(Vector&          abs_h2o,
                const ArrayOfArrayOfSpeciesTag& abs_species,
                const Matrix&    abs_vmrs,
                const Verbosity&)
{
  const Index h2o_index 
    = find_first_species_tg( abs_species,
                             species_index_from_species_name("H2O") );

  if ( h2o_index < 0 )
    throw runtime_error("No tag group contains water!");
  
  abs_h2o.resize( abs_vmrs.ncols() );
  abs_h2o = abs_vmrs(h2o_index,Range(joker));   
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_n2Set(Vector&            abs_n2,
               const ArrayOfArrayOfSpeciesTag& abs_species,
               const Matrix&    abs_vmrs,
               const Verbosity&)
{
  const Index n2_index 
    = find_first_species_tg( abs_species,
                             species_index_from_species_name("N2") );

  if ( n2_index < 0 )
    throw runtime_error("No tag group contains nitrogen!");

  abs_n2.resize( abs_vmrs.ncols() );
  abs_n2 = abs_vmrs(n2_index,Range(joker));   
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
void abs_coefCalc(// WS Output:
                  Matrix&                         abs_coef,
                  ArrayOfMatrix&                  abs_coef_per_species,
                  // WS Input:                 
                  const ArrayOfArrayOfSpeciesTag& tgs,
                  const Vector&                   f_grid,
                  const Vector&                   abs_p,
                  const Vector&                   abs_t,
                  const Vector&                   abs_n2,
                  const Vector&                   abs_h2o,
                  const Matrix&                   abs_vmrs,
                  const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                  const ArrayOfLineshapeSpec&     abs_lineshape,
                  const ArrayOfString&            abs_cont_names,
                  const ArrayOfString&            abs_cont_models,
                  const ArrayOfVector&            abs_cont_parameters,
                  const ArrayOfVector&            isotopologue_ratios,
                  const Verbosity&                verbosity)
{
  // Dimension checks are performed in the executed functions

  // allocate local variable to hold the cross sections per tag group
  ArrayOfMatrix abs_xsec_per_species;
  
  abs_xsec_per_speciesInit(abs_xsec_per_species, tgs, f_grid, abs_p, verbosity);

  abs_xsec_per_speciesAddLines(abs_xsec_per_species,
                               tgs,
                               f_grid,
                               abs_p,
                               abs_t,
                               abs_vmrs,
                               abs_lines_per_species,
                               abs_lineshape,
                               isotopologue_ratios,
                               verbosity);

  abs_xsec_per_speciesAddConts(abs_xsec_per_species,
                               tgs,
                               f_grid,
                               abs_p,
                               abs_t,
                               abs_n2,
                               abs_h2o,
                               abs_vmrs,
                               abs_cont_names,
                               abs_cont_parameters,
                               abs_cont_models,
                               verbosity);

  abs_coefCalcFromXsec(abs_coef,
                       abs_coef_per_species,
                       abs_xsec_per_species,
                       abs_vmrs,
                       abs_p,
                       abs_t,
                       verbosity);

}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_coefCalcSaveMemory(// WS Output:
                            Matrix&                         abs_coef,
                            // WS Input:                 
                            const ArrayOfArrayOfSpeciesTag& tgs,
                            const Vector&                   f_grid,
                            const Vector&                   abs_p,
                            const Vector&                   abs_t,
                            const Vector&                   abs_n2,
                            const Vector&                   abs_h2o,
                            const Matrix&                   abs_vmrs,
                            const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                            const ArrayOfLineshapeSpec&     abs_lineshape,
                            const ArrayOfString&            abs_cont_names,
                            const ArrayOfString&            abs_cont_models,
                            const ArrayOfVector&            abs_cont_parameters,
                            const ArrayOfVector&            isotopologue_ratios,
                            const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Dimension checks are performed in the executed functions

  // Allocate local variable to hold the cross sections per tag group:
  ArrayOfMatrix abs_xsec_per_species;

  // Allocate local variable to hold the absorption for each tag group:
  Matrix this_abs;

  // Allocate local variable to hold abs_coef_per_species for each tag
  // group. This is just a dummy, we need it, since it is a formal
  // argument of abs_coefCalcFromXsec.
  ArrayOfMatrix this_abs_coef_per_species;

  // Define variable to hold a dummy list of tag groups with only 1 element:
  ArrayOfArrayOfSpeciesTag this_tg(1);

  // Local list of VMRs, only 1 element:
  Matrix this_vmr(1,abs_p.nelem());

  // Local abs_lines_per_species, only 1 element:
  ArrayOfArrayOfLineRecord these_lines(1);

  // Local lineshape list, only 1 element:
  ArrayOfLineshapeSpec     this_lineshape(1);

  // Initialize the output variable abs_coef:
  abs_coef.resize( f_grid.nelem(), abs_p.nelem() );
  abs_coef = 0;                      // Matpack can set all elements like this.

  out3 << "  Number of tag groups to do: " << tgs.nelem() << "\n";

  // Loop over all species:
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
      out3 << "  Doing tag group " << i << ".\n";

      // Get a dummy list of tag groups with only the current element:
      this_tg[0].resize(tgs[i].nelem());
      this_tg[0] = tgs[i];

      // VMR for this species:
      this_vmr(0,joker) = abs_vmrs(i,joker);

      // List of lines:
      these_lines[0].resize(abs_lines_per_species[i].nelem());
      these_lines[0] = abs_lines_per_species[i];

      // List of lineshapes:
      this_lineshape[0] = abs_lineshape[i];

      abs_xsec_per_speciesInit(abs_xsec_per_species, this_tg, f_grid, abs_p, verbosity);

      abs_xsec_per_speciesAddLines(abs_xsec_per_species,
                                   this_tg,
                                   f_grid,
                                   abs_p,
                                   abs_t,
                                   this_vmr,
                                   these_lines,
                                   this_lineshape,
                                   isotopologue_ratios,
                                   verbosity);

      abs_xsec_per_speciesAddConts(abs_xsec_per_species,
                                   this_tg,
                                   f_grid,
                                   abs_p,
                                   abs_t,
                                   abs_n2,
                                   abs_h2o,
                                   this_vmr,
                                   abs_cont_names,
                                   abs_cont_parameters,
                                   abs_cont_models,
                                   verbosity);

      abs_coefCalcFromXsec(this_abs,
                           this_abs_coef_per_species,
                           abs_xsec_per_species,
                           this_vmr,
                           abs_p,
                           abs_t,
                           verbosity);

      // Add absorption of this species to total absorption:
      assert(abs_coef.nrows()==this_abs.nrows());
      assert(abs_coef.ncols()==this_abs.ncols());
      abs_coef += this_abs;
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_coefCalcFromXsec(// WS Output:
                          Matrix&              abs_coef,
                          ArrayOfMatrix&       abs_coef_per_species,
                          // WS Input:         
                          const ArrayOfMatrix& abs_xsec_per_species,
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
  abs_coef = 0;                      // Matpack can set all elements like this.

  abs_coef_per_species.resize( abs_xsec_per_species.nelem() );

  out3 << "  Computing abs_coef and abs_coef_per_species from abs_xsec_per_species.\n";
  // Loop through all tag groups
  for ( Index i=0; i<abs_xsec_per_species.nelem(); ++i )
    {
      out3 << "  Tag group " << i << "\n";

      // Make this element of abs_xsec_per_species the right size:
      abs_coef_per_species[i].resize( abs_xsec_per_species[i].nrows(), abs_xsec_per_species[i].ncols() );
      abs_coef_per_species[i] = 0;        // Initialize all elements to 0.

      // Loop through all altitudes
      for ( Index j=0; j<abs_xsec_per_species[i].ncols(); j++)
        {
          // Calculate total number density from pressure and temperature.
          const Numeric n = number_density(abs_p[j],abs_t[j]);

          // Loop through all frequencies
          for ( Index k=0; k<abs_xsec_per_species[i].nrows(); k++)
            {
              abs_coef_per_species[i](k,j) = abs_xsec_per_species[i](k,j) * n * abs_vmrs(i,j);
            }
        }

      // Add up to the total absorption:
      abs_coef += abs_coef_per_species[i];     // In Matpack you can use the +=
                                               // operator to do elementwise addition.
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesInit(// WS Output:
                              ArrayOfMatrix&   abs_xsec_per_species,
                              // WS Input:
                              const ArrayOfArrayOfSpeciesTag& tgs,
                              const Vector&    f_grid,
                              const Vector&    abs_p,
                              const Verbosity& verbosity
                              )
{
  CREATE_OUT3;
  
  // Initialize abs_xsec_per_species. The array dimension of abs_xsec_per_species
  // is the same as that of abs_lines_per_species.
  abs_xsec_per_species.resize( tgs.nelem() );

  // Loop abs_xsec_per_species and make each matrix the right size,
  // initializing to zero:
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
      // Make this element of abs_coef_per_species the right size:
      abs_xsec_per_species[i].resize( f_grid.nelem(), abs_p.nelem() );
      abs_xsec_per_species[i] = 0;       // Matpack can set all elements like this.
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
                                  // WS Input:             
                                  const ArrayOfArrayOfSpeciesTag&  tgs,
                                  const Vector&                    f_grid,
                                  const Vector&                    abs_p,
                                  const Vector&                    abs_t,
                                  const Matrix&                    abs_vmrs,
                                  const ArrayOfArrayOfLineRecord&  abs_lines_per_species,
                                  const ArrayOfLineshapeSpec&      abs_lineshape,
                                  const ArrayOfVector&             isotopologue_ratios,
                                  const Verbosity&                 verbosity)
{
  CREATE_OUT3;
  
  // Check that all temperatures are at least 0 K. (Negative Kelvin
  // temperatures are unphysical.)  
  if ( min(abs_t) < 0 )
    {
      ostringstream os;
      os << "Temperature must be at least 0 K. But you request an absorption\n"
         << "calculation at " << min(abs_t) << " K!"; 
      throw runtime_error(os.str());
    }

  // Check that all paramters that should have the number of tag
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
        throw runtime_error(os.str());
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
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
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
      if ( 0 < ll.nelem() && tgs[i].nelem() && !tgs[i][0].Zeeman() )
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
              throw runtime_error(os.str());
            }

          // Get the name of the species. The member function name of a
          // LineRecord returns the full name (species + isotopologue). So
          // it is for us more convenient to get the species index
          // from the LineRecord member function Species(), and then
          // use this to look up the species name in species_data.
          extern const Array<SpeciesRecord> species_data;
          String species_name = species_data[ll[0].Species()].Name();

          // Get the name of the lineshape. For that we use the member
          // function Ind_ls() to the lineshape data ls, which returns
          // an index. With that index we can go into lineshape_data
          // to get the name.
          // We will need this for safety checks later on.
          extern const Array<LineshapeRecord> lineshape_data;
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

          xsec_species( abs_xsec_per_species[i],
                        f_grid,
                        abs_p,
                        abs_t,
                        abs_vmrs,
                        tgs,
                        i,
                        ll,
                        ls.Ind_ls(),
                        ls.Ind_lsn(),
                        ls.Cutoff(),
                        isotopologue_ratios,
                        verbosity );
          // Note that we call xsec_species with a row of abs_vmrs,
          // selected by the above Matpack expression. This is
          // possible, because xsec_species is using Views.
        }

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
                                  // WS Input:             
                                  const ArrayOfArrayOfSpeciesTag&  tgs,
                                  const Vector&                    f_grid,
                                  const Vector&                    abs_p,
                                  const Vector&                    abs_t,
                                  const Vector&                    abs_n2,
                                  const Vector&                    abs_h2o,
                                  const Matrix&                    abs_vmrs,
                                  const ArrayOfString&             abs_cont_names,
                                  const ArrayOfVector&             abs_cont_parameters,
                                  const ArrayOfString&             abs_cont_models,
                                  const Verbosity&                 verbosity)
{
  CREATE_OUT3;
  
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

  // ...and that indeed the names match valid continuum models:
  for ( Index i=0; i<abs_cont_names.nelem(); ++i )
    {
      check_continuum_model(abs_cont_names[i]);
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

  // We do checks on abs_h2o and abs_n2 later, because we only want to
  // do the check if the parameter are really needed.


  out3 << "  Calculating continuum spectra.\n";

  // Loop tag groups:
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
      extern const Array<SpeciesRecord> species_data; 

      // Go through the tags in the current tag group to see if they
      // are continuum tags:  
      for ( Index s=0; s<tgs[i].nelem(); ++s )
        {
          // First of all, we have to make sure that this is not a
          // tag that means `all isotopologues', because this should not
          // include continuum. For such tags, tag.Isotopologue() will
          // return the number of isotopologues (i.e., one more than the
          // allowed index range).
          if ( tgs[i][s].Isotopologue() <
               species_data[tgs[i][s].Species()].Isotopologue().nelem() )
            {
              // If we get here, it means that the tag describes a
              // specific isotopologue. Could be a continuum tag!
                
              // The if clause below checks for continuum tags.
              // It does the following:
              //
              // 1. species_data contains the lookup table of species
              //          specific data. We need the right entry in this
              //          table. The index of this is obtained by calling member function
              //          Species() on the tag. Thus we have:
              //          species_data[tgs[i][s].Species()].
              //
              // 2. Member function Isotopologue() on the above biest gives
              //    us the array of isotopologue specific data. This we have
              //    to subscribe with the right isotopologue index, which we
              //    get by calling member function Isotopologue on the
              //    tag. Thus we have:
              //    Isotopologue()[tgs[i][s].Isotopologue()]
              //
              // 3. Finally, from the isotopologue data we need to get the flag
              //    whether this is a continuum.
              if ( species_data[tgs[i][s].Species()].Isotopologue()[tgs[i][s].Isotopologue()].isContinuum() )
                {
                  // We have identified a continuum tag!

                  // Get only the continuum name. The full tag name is something like:
                  // H2O-HITRAN96Self-*-*. We want only the `H2O-HITRAN96Self' part:
                  const String name =
                    species_data[tgs[i][s].Species()].Name() + "-"
                    + species_data[tgs[i][s].Species()].Isotopologue()[tgs[i][s].Isotopologue()].Name();
  
                  // Check that this corresponds to a valid continuum
                  // model:
                  check_continuum_model(name);

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
                  
                  {
                    ostringstream os;
                    os << "  Adding " << name
                       << " to tag group " << i << ".\n";
                    out3 << os.str();
                  }

                  // find the options for this continuum tag from the input array
                  // of options. The actual field of the array is n:
                  const String ContOption = abs_cont_models[n];


                  // ------------------------------------------------------------------
                  // Now is the time to check whether abs_h2o and
                  // abs_n2 are ok!

                  // abs_h2o has a global scalar default value of -1. We throw an
                  // appropriate error message here if we find this, since most
                  // continuum models require it. (abs_h2o is the H2O VMR to be used
                  // for the continua of other species, such as O2.)

                  if ( -.99 > abs_h2o[0] )
                    {
                      ostringstream os;
                      os << "The variable abs_h2o seems to be set to its global default\n"
                         << "value of -1. You have to set this to a H2O VMR profile if\n"
                         << "you want to use absorption contiua. If you are calling\n"
                         << "absorption routines directly, or on the fly, you could\n"
                         << "use for example the method *abs_h2oSet*.\n"
                         << "If you are generating an absorption lookup table with\n"
                         << "abs_lookupCreate, it should be enough to add a H2O species\n"
                         << "to your calculation to fix this problem.";
                      throw runtime_error(os.str());
                    }

                  // If h2o_abs is not set to the default value, it
                  // must have the same size as the pressure grid:
                  if ( abs_h2o.nelem() != abs_p.nelem() )
                    {
                      ostringstream os;
                      os << "Variable abs_h2o must have the same dimension as abs_p.\n"
                         << "abs_h2o.nelem() = " << abs_h2o.nelem() << '\n'
                         << "abs_p.nelem() = " << abs_p.nelem();
                      throw runtime_error(os.str());
                    }

                  // For abs_n2 the situation is slightly
                  // different. The global scalar default value is a
                  // reasonable estimate for the N2 profile, so we
                  // just have to expand it to a vector here. Because
                  // we cannot modify abs_n2, we have to make a local
                  // copy in any case.

                    Vector n2_prof(abs_p.nelem());
                    if ( abs_n2.nelem() == abs_p.nelem() )
                      {
                        n2_prof = abs_n2;
                      }
                    else
                      {
                        if (1==abs_n2.nelem())
                          {
                            // We seem to have found the global
                            // default value. Expand this to a vector
                            // with the right length, by copying it to
                            // all elements of n2_prof.
                            n2_prof = abs_n2[0];         
                          }
                        else
                          {
                            ostringstream os;
                            os << "Variable abs_n2 must have dimension 1, or\n"
                               << "the same dimension as abs_p.\n"
                               << "abs_n2.nelem() = " << abs_n2.nelem() << '\n'
                               << "abs_p.nelem() = " << abs_p.nelem();
                            throw runtime_error(os.str());
                          }
                      }

                  // ------------------------------------------------------------------


                  // Add the continuum for this tag. The parameters in
                  // this call should be clear. The vmr is in
                  // abs_vmrs(i,Range(joker)). The other vmr variable, `abs_h2o'
                  // contains the real H2O vmr, which is needed for
                  // the oxygen continuum.
                  xsec_continuum_tag( abs_xsec_per_species[i],
                                      name,
                                      abs_cont_parameters[n],
                                      abs_cont_models[n], 
                                      f_grid,
                                      abs_p,
                                      abs_t,
                                      n2_prof,
                                      abs_h2o,
                                      abs_vmrs(i,Range(joker)),
                                      verbosity );
                  // Calling this function with a row of Matrix abs_vmrs
                  // is possible because it uses Views.
                }
            }
        }
    }

}



//======================================================================
//             Methods related to continua
//======================================================================

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_cont_descriptionInit(// WS Output:
                              ArrayOfString& names,
                              ArrayOfString& options, 
                              ArrayOfVector& parameters,
                              const Verbosity& verbosity)
{
  CREATE_OUT2;
  
  names.resize(0);
  options.resize(0);
  parameters.resize(0);
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
void abs_mat_per_speciesAddFromAbsCoefPerSpecies(// WS Output:
                               Tensor4&       abs_mat_per_species,
                               // WS Input:
                               const ArrayOfMatrix& abs_coef_per_species,
                               const Verbosity&)
{
  // abs_mat_per_species has format
  // [ abs_species, f_grid, stokes_dim, stokes_dim ].
  // abs_coef_per_species has format ArrayOfMatrix (over species),
  // where for each species the matrix has format [f_grid, abs_p].

  
  // Set stokes_dim (and check that the last two dimensions of
  // abs_mat_per_species really are equal).
  Index nr, nc, stokes_dim;
  // In the two stokes dimensions, abs_mat_per_species should be a
  // square matrix of stokes_dim*stokes_dim. Check this, and determine
  // stokes_dim:
  nr = abs_mat_per_species.nrows();
  nc = abs_mat_per_species.ncols();
  if ( nr!=nc )
  {
    ostringstream os;
    os << "The last two dimensions of abs_mat_per_species must be equal (stokes_dim).\n"
    << "But here they are " << nr << " and " << nc << ".";
    throw runtime_error( os.str() );
  }
  stokes_dim = nr;       // Could be nc here too, since they are the same.

  
  Index n_species = abs_coef_per_species.nelem(); // # species

  if (0==n_species)
    {
      ostringstream os;
      os << "Must have at least one species.";
      throw runtime_error(os.str());
    }

  Index n_f       = abs_coef_per_species[0].nrows(); // # frequencies

  // # pressures must be 1:
  if (1!=abs_coef_per_species[0].ncols())
    {
      ostringstream os;
      os << "Must have exactly one pressure.";
      throw runtime_error(os.str());
    }
  
  // Check species dimension of abs_mat_per_species
  if ( abs_mat_per_species.nbooks()!=n_species )
  {
    ostringstream os;
    os << "Species dimension of abs_mat_per_species does not\n"
       << "match abs_coef_per_species.";
    throw runtime_error( os.str() );
  }
  
  // Check frequency dimension of abs_mat_per_species
  if ( abs_mat_per_species.npages()!=n_f )
  {
    ostringstream os;
    os << "Frequency dimension of abs_mat_per_species does not\n"
       << "match abs_coef_per_species.";
    throw runtime_error( os.str() );
  }
  
  // Loop species and stokes dimensions, and add to abs_mat_per_species:
  for ( Index si=0; si<n_species; ++si )
    for ( Index ii=0; ii<stokes_dim; ++ii )
      abs_mat_per_species(si,joker,ii, ii) += abs_coef_per_species[si](joker,0);

}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_mat_per_speciesInit(//WS Output
                             Tensor4&                        abs_mat_per_species,
                             //WS Input
                             const ArrayOfArrayOfSpeciesTag& abs_species,
                             const Vector&                   f_grid,
                             const Index&                    stokes_dim,
                             const Verbosity&                
                            )
{
    
    Index nf = f_grid.nelem();
    
    if(abs_species.nelem() > 0 )
    {
        if(nf > 0)
        {
            if(stokes_dim > 0)
            {
                abs_mat_per_species.resize(abs_species.nelem(),nf, stokes_dim, stokes_dim);
                abs_mat_per_species = 0;
            }
            else throw  runtime_error("stokes_dim = 0");
        }
        else throw runtime_error("nf = 0");
    }
    else throw runtime_error("abs_species.nelem() = 0");

}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_mat_per_speciesAddLBL(// WS Output:
                           Tensor4& abs_mat_per_species,
                           // WS Input:
                           const Vector& f_grid,
                           const ArrayOfArrayOfSpeciesTag& abs_species,
                           const Vector& abs_n2,
                           const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                           const ArrayOfLineshapeSpec& abs_lineshape,
                           const ArrayOfString& abs_cont_names,
                           const ArrayOfString& abs_cont_models,
                           const ArrayOfVector& abs_cont_parameters,
                           const ArrayOfVector& isotopologue_ratios,
                           const Numeric& rte_pressure,
                           const Numeric& rte_temperature,
                           const Vector& rte_vmr_list,
                           const Numeric& rte_doppler,
                           const Verbosity& verbosity)
{
  CREATE_OUT3;

  
  // Define communication variables for the actual absorption calculation:
  
  // Output of AbsInputFromRteScalars:
  Vector        abs_p;
  Vector        abs_t;
  Matrix        abs_vmrs;
  // Output of abs_h2oSet:
  Vector          abs_h2o;
  // Output of abs_coefCalc:
  Matrix                         abs_coef;
  ArrayOfMatrix                  abs_coef_per_species;
  
  
  // This construct is needed for the Doppler treatment,
  // since that also modifies the local frequency grid.)
  Vector local_f_grid;
  const Vector* f_grid_pointer;
  
  // Make pointer point to original.
  f_grid_pointer = &f_grid;
  
  // Doppler treatment, do this only if there is a non-zero Doppler
  // shift. We do this after the frequency selection, so in the case
  // that we have only a single frequency, we have to shift only that!
  
  // Unfortunately, we need yet another local copy of f_grid. In
  // constrast to the frequency selection, we here want to modify the
  // actual frequency values inside!
  Vector local_doppler_f_grid;
  if (0==rte_doppler)
  {
    out3 << "  Doppler shift: None\n";
  }
  else
  {
    ostringstream os;
    os << "  Doppler shift: " << rte_doppler << " Hz\n";
    out3 << os.str();
    
    Numeric local_doppler;
    NumericScale( local_doppler, rte_doppler, -1, verbosity );
    // I could just have multiplied by -1 directly, but I like using
    // the WSM here.
    
    VectorAddScalar( local_doppler_f_grid,  *f_grid_pointer, local_doppler, verbosity );
    
    // Make pointer point to the doppler shifted frequency grid.
    f_grid_pointer = &local_doppler_f_grid;
  }
  
  AbsInputFromRteScalars(abs_p,
                         abs_t,
                         abs_vmrs,
                         rte_pressure,
                         rte_temperature,
                         rte_vmr_list,
                         verbosity);
  
  abs_h2oSet(abs_h2o, abs_species, abs_vmrs, verbosity);
  
  abs_coefCalc(abs_coef,
               abs_coef_per_species,
               abs_species,
               *f_grid_pointer,
               abs_p,
               abs_t,
               abs_n2,
               abs_h2o,
               abs_vmrs,
               abs_lines_per_species,
               abs_lineshape,
               abs_cont_names,
               abs_cont_models,
               abs_cont_parameters,
               isotopologue_ratios,
               verbosity);
  
  // Now add abs_coef_per_species to abs_mat_per_species:
  abs_mat_per_speciesAddFromAbsCoefPerSpecies(abs_mat_per_species,
                                              abs_coef_per_species,
                                              verbosity);

}



/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridSelectFIndex(// WS Output:
                        Vector& f_grid,
                        // WS Input:
                        const Index& f_index,
                        const Verbosity&)
{
  // Prepare f_grid. f_index < 0 means retain all frequencies, but
  // f_index >= 0 means to retain only that frequency. 
  if ( f_index >= 0 )
    {
      // Check that f_index is inside f_grid:
      if ( f_index >= f_grid.nelem() )
      {
        ostringstream os;
        os << "The frequency index you want is outside f_grid.\n"
           << "You have " << f_index
           << ", the largest allowed value is " << f_grid.nelem()-1 << ".";
        throw runtime_error( os.str() );
      }

      Numeric this_f = f_grid[f_index];
      f_grid.resize(1);
      f_grid = this_f;
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void isotopologue_ratiosInitFromBuiltin(ArrayOfVector& isotopologue_ratios,
                                        const Verbosity&)
{
    extern Array<SpeciesRecord> species_data;

    isotopologue_ratios.resize(species_data.nelem());

    for (Index isd = 0; isd < species_data.nelem(); isd++)
    {
        const Array<IsotopologueRecord>& air = species_data[isd].Isotopologue();
        isotopologue_ratios[isd].resize(air.nelem());
        for (Index iir = 0; iir < air.nelem(); iir++)
        {
            isotopologue_ratios[isd][iir] = air[iir].Abundance();
        }
    }
}


#ifdef ENABLE_NETCDF
/* Workspace method: Doxygen documentation will be auto-generated */
/* Included by Claudia Emde, 20100707 */
void WriteMolTau(//WS Input
                 const Vector& f_grid, 
                 const Tensor3& z_field,
                 const Tensor7& abs_mat_field,
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

  if ((retval = nc_def_dim(ncid, "nstk", (int) abs_mat_field.nbooks(), &stokes_dimid)))
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

  double tau[z_field.npages()-1][f_grid.nelem()][abs_mat_field.nbooks()][abs_mat_field.nbooks()];

  // Initialize tau
    for (int iz=0; iz<z_field.npages()-1; iz++)
        for (int iv=0; iv<f_grid.nelem(); iv++)
            for (int is1=0; iv<abs_mat_field.nbooks(); iv++)
                for (int is2=0; iv<abs_mat_field.nbooks(); iv++)
                    tau[iz][iv][is1][is2] = 0.0;

  // Calculate average tau for layers
  for (int is=0; is<abs_mat_field.nlibraries(); is++)
    for (int iz=0; iz<z_field.npages()-1; iz++)
      for (int iv=0; iv<f_grid.nelem(); iv++)
          for (int is1=0; iv<abs_mat_field.nbooks(); iv++)
              for (int is2=0; iv<abs_mat_field.nbooks(); iv++)
                // sum up all species
                tau[iz][iv][is1][is2] += 0.5 * (abs_mat_field(is,f_grid.nelem()-1-iv,is1,is2,z_field.npages()-1-iz,0,0)+
                                    abs_mat_field(is,f_grid.nelem()-1-iv,is1,is2,z_field.npages()-2-iz,0,0))
                *(z_field(z_field.npages()-1-iz,0,0)-z_field(z_field.npages()-2-iz,0,0));
  
  
  if ((retval = nc_put_var_double (ncid, tau_varid, &tau[0][0][0][0])))
    nca_error (retval, "nc_put_var");
  
  // Close the file
  if ((retval = nc_close(ncid)))
    nca_error (retval, "nc_close");

}

#else

void WriteMolTau(//WS Input
                 const Vector& f_grid _U_,
                 const Tensor3& z_field _U_,
                 const Tensor7& abs_mat_field _U_,
                 const Index& atmosphere_dim _U_,
                 //Keyword
                 const String& filename _U_,
                 const Verbosity&)
{
  throw runtime_error("The workspace method WriteMolTau is not available"
                      "because ARTS was compiled without NetCDF support.");
}
                 
#endif /* ENABLE_NETCDF */

