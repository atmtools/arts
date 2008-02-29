/* Copyright (C) 2000-2007
   Stefan Buehler   <sbuehler@ltu.se>
   Patrick Eriksson <patrick@rss.chalmers.se>
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

   This is the file from arts-1-0, back-ported to arts-1-1.

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
#include "auto_wsv.h"
#include "auto_md.h"
#include "math_funcs.h"
#include "make_array.h"
#include "physics_funcs.h"
#include "continua.h"
#include "make_vector.h"
#include "check_input.h"


/* Workspace method: Doxygen documentation will be auto-generated */
void AbsInputFromRteScalars(// WS Output:
                            Vector&        f_grid,
                            Vector&        abs_p,
                            Vector&        abs_t,
                            Matrix&        abs_vmrs,
                            // WS Input:
                            const Index&   f_index,
                            const Numeric& rte_pressure,
                            const Numeric& rte_temperature,
                            const Vector&  rte_vmr_list)
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

      Vector dummy = f_grid;
      f_grid.resize(1);
      f_grid = dummy[f_index];
    }

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
                          const ArrayOfArrayOfSpeciesTag& tgs)
{
  // Make abs_lines_per_species the right size:
  abs_lines_per_species.resize( tgs.nelem() );
  
  for (Index i=0; i<tgs.nelem(); ++i)
    {
      abs_lines_per_species[i].resize(0);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReadFromHitran(// WS Output:
                         ArrayOfLineRecord& abs_lines,
                          // Control Parameters:
                         const String& filename,
                         const Numeric& fmin,
                         const Numeric& fmax)
{
  ifstream is;

  // Reset lines in case it already existed:
  abs_lines.resize(0);

  out2 << "  Reading file: " << filename << "\n";
  open_input_file(is, filename);

  bool go_on = true;
  while ( go_on )
    {
      LineRecord lr;
      if ( lr.ReadFromHitranStream(is) )
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
void abs_linesReadFromHitran2004(// WS Output:
                             ArrayOfLineRecord& abs_lines,
                             // Control Parameters:
                             const String& filename,
                             const Numeric& fmin,
                             const Numeric& fmax)
{
  ifstream is;

  // Reset lines in case it already existed:
  abs_lines.resize(0);

  out2 << "  Reading file: " << filename << "\n";
  open_input_file(is, filename);

  bool go_on = true;
  while ( go_on )
    {
      LineRecord lr;
      if ( lr.ReadFromHitran2004Stream(is) )
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
                          const Numeric& fmax)
{
  ifstream is;

  // Reset lines in case it already existed:
  abs_lines.resize(0);

  out2 << "  Reading file: " << filename << "\n";
  open_input_file(is, filename);

  bool go_on = true;
  while ( go_on )
    {
      LineRecord lr;
      if ( lr.ReadFromMytran2Stream(is) )
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
                      const Numeric& fmax)
{
  ifstream is;

  // Reset lines in case it already existed:
  abs_lines.resize(0);

  out2 << "  Reading file: " << filename << "\n";
  open_input_file(is, filename);

  bool go_on = true;
  while ( go_on )
    {
      LineRecord lr;
      if ( lr.ReadFromJplStream(is) )
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
                       const Numeric& fmax)
{
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
    if ( v!=lr.Version() )
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
      if ( lr.ReadFromArtsStream(is) )
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


/* FIXME OLE: Do we need this function? */
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
                                    const Vector& fmax)
{
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
          abs_linesReadFromHitran( abs_lines, real_filenames[i], real_fmin[i], real_fmax[i] );
        }
      else if ( "HITRAN04"==real_formats[i] )
        {
          abs_linesReadFromHitran2004( abs_lines, real_filenames[i], real_fmin[i], real_fmax[i] );
        }
      else if ( "MYTRAN2"==real_formats[i] )
        {
          abs_linesReadFromMytran2( abs_lines, real_filenames[i], real_fmin[i], real_fmax[i] );
        }
      else if ( "JPL"==real_formats[i] )
        {
          abs_linesReadFromJpl( abs_lines, real_filenames[i], real_fmin[i], real_fmax[i] );
        }
      else if ( "ARTS"==real_formats[i] )
        {
          abs_linesReadFromArts( abs_lines, real_filenames[i], real_fmin[i], real_fmax[i] );
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
      abs_lines_per_speciesCreateFromLines( these_abs_lines_per_species, abs_lines, these_tgs );

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
                                  const ArrayOfLineRecord&   abs_lines,
                                  const ArrayOfArrayOfSpeciesTag&           tgs)
{
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

              // Test isotope. The isotope can either match directly, or
              // the Isotope of the tag can be one larger than the
              // number of isotopes, which means `all'. Test the second
              // condition first, since this will probably be more often
              // used.
              if ( this_tag.Isotope() != this_line.SpeciesData().Isotope().nelem() )
                if ( this_tag.Isotope() != this_line.Isotope() )
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
                                ArrayOfArrayOfLineRecord& abs_lines_per_species)
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
                         const Vector& f_grid)
{

  // Make sure abs_lines_per_species and abs_lineshape have the same dimension:
  if ( abs_lines_per_species.nelem() != abs_lineshape.nelem() ) 
    {
      ostringstream os;
      os << "Dimension of abs_lines_per_species does\n"
         << "not match that of abs_lineshape.";
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
      // Get cutoff frequency of this tag group:
      Numeric cutoff = abs_lineshape[i].Cutoff();

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
                                    const String& basename)
{
  // Species lookup data:
  extern const Array<SpeciesRecord> species_data;

  // We want to make lists of included and excluded species:
  ArrayOfString included(0), excluded(0);

  tgs.resize(0);

  for ( Index i=0; i<species_data.nelem(); ++i )
    {
      const String specname = species_data[i].Name();
      const String filename = basename + "." + specname + ".xml";

      // Try to open VMR file:
      try
        {
          ifstream file;
          open_input_file(file, filename);

          // Ok, if we get here the file was found.

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
      catch (runtime_error x)
        {
          // Ok, the file for the species could not be found.
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
                         const Numeric&           cutoff)
{
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
                            const Vector&         cutoff )
{
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

#if 0

void raw_vmrsReadFromScenario(// WS Output:
                              ArrayOfMatrix&   raw_vmrs,
                              // WS Input:
                              const ArrayOfArrayOfSpeciesTag&     tgs,
                              // Control Parameters:
                               const String&    basename)
{
  // The species lookup data:
  extern const Array<SpeciesRecord> species_data;

  // We need to read one profile for each tag group.
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
      // Determine the name.
      String name =
        basename + "." +
        species_data[tgs[i][0].Species()].Name() + ".aa";
      
      // Add an element for this tag group to the vmr profiles:
      raw_vmrs.push_back(Matrix());
      
      // Read the VMR:
      // (We use the workspace method MatrixReadAscii for this.)
      MatrixReadAscii(raw_vmrs[i],"",name);
      
      // state the source of profile.
      out3 << "  " << species_data[tgs[i][0].Species()].Name()
           << " profile read from file: " << name << "\n";
    }
}

/** Reads in the profiles from the specified files in filenames for the tag list seltags
    and for all the other tags the atmospheric profile from the general scenario stated in 
    basename.

   \param    raw_vmrs       volume mixing ratios per tag group    (output)
   \param    tgs            the list of tag groups                (input)
   \param    seltags        selected tags for special input files (input)
   \param    filenames      specific files for list of seltags    (input)
   \param    basename       general scenario base name            (input)

   \author Thomas Kuhn / Stefan Buehler
   \date 2001-08-02 / 2001-09-19
 */ 
void raw_vmrsReadFromFiles(// WS Output:
                           ArrayOfMatrix&   raw_vmrs,
                           // WS Input:
                           const ArrayOfArrayOfSpeciesTag& tgs,
                           // Control Parameters:
                           const ArrayOfString&  seltags,
                           const ArrayOfString&  filenames,
                           const String&         basename)
{
  // The species lookup data:
  extern const Array<SpeciesRecord> species_data;
  ArrayOfIndex i_seltags;

  // Get indices of seltags in tgs. The function will throw an error
  // unless each tg in seltags corresponds exactly to a tg
  // in tgs. 
  get_tagindex_for_Strings( i_seltags, tgs, seltags );

  // Now we have to build an Array of the filenames to use for all the
  // tag groups. Initialize it with basename...
  ArrayOfString true_filenames(tgs.nelem());
  true_filenames = basename + '.';

  // Initialize to the standard scenario (see function
  // raw_vmrsReadFromScenario): 
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
     true_filenames[i] +=
        species_data[tgs[i][0].Species()].Name() + ".aa";
     // Should be identical to how the filenames are constructed in
     // raw_vmrsReadFromScenario! 
    }

  // Replace the names of the tgs given in seltags with the names in
  // filenames: 
  for ( Index i=0; i<seltags.nelem(); ++i )
    {
      true_filenames[i_seltags[i]] = filenames[i];
    }
    
  // Read the files:
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
      // Add an element for this tag group to the vmr profiles:
      raw_vmrs.push_back(Matrix());
      
      // Read the VMR:
      // (We use the workspace method MatrixReadAscii for this.)
      MatrixReadAscii(raw_vmrs[i],"",true_filenames[i]);
      
      // state the source of profile.
      out3 << "  " << species_data[tgs[i][0].Species()].Name()
           << " profile read from file: " << true_filenames[i] << "\n";
    }
}

/** Calculates the water vapor saturation volume mixing ratio (VMR) in the 
    vertical range where liquid or ice clouds are in the atmosphere.
    At the pressure/altitude grid points where the liquid water content (LWC) 
    or ice water content (IWC) is larger than zero the H2O-VMR is set to 
    liquid water/ice saturation VMR (=e_s/p_tot).<br>
    The saturation pressure (e_s) is calculated according to the Goff-Gratch equations 
    (see functions WVSatPressureLiquidWater and WVSatPressureIce in file continua.cc).

    Before adjustment:<br>
    <table border=0 frame=void rules=none>
    <tr>
      <td>VMR<sub>H2O</sub></td> <td>=</td> <td>H2O volume mixing ratio</td>
    </tr>
    <tr>
      <td>P<sub>tot</sub></td>   <td>=</td> <td>total air pressure</td>
    </tr>
    <tr>
      <td>P<sub>H2O</sub></td>  <td>=</td> <td>H2O partial pressure</td>
    </tr>
    <tr>
      <td>P<sub>dry</sub></td>   <td>=</td> <td>dry air partial pressure</td>
    </tr>
    </table>

    after adjustment:<br>
    <table border=0 frame=void rules=none>
    <tr>
       <td>VMR'<sub>H2O</sub></td>   <td>=</td>   <td>H2O volume mixing ratio</td>  
    </tr>
    <tr>
       <td>P'<sub>tot</sub></td>     <td>=</td>   <td>total air pressure</td>  
    </tr>
    <tr>
       <td>P'<sub>H2O</sub></td>     <td>=</td>   <td>H2O partial pressure</td>  
    </tr>
    </table>

    H2O saturation pressure over water / ice:<br>
    e<sub>s</sub><br>

    calculation of the adjustment:<br>
    <table border=0 frame=void rules=none>
    <tr>
       <td>P<sub>H2O</sub></td>    <td>=</td>  <td>VMR<sub>H2O</sub> * P<sub>tot</sub></td>
    </tr>
    <tr>
       <td>P<sub>dry</sub></td>    <td>=</td>  <td>P<sub>tot</sub> * (1.0 - VMR<sub>H2O</sub>)</td>
    </tr>
    <tr>
       <td>P'<sub>H2O</sub></td>   <td>=</td>  <td>e<sub>s</sub></td>
    </tr>
    <tr>
       <td>P'<sub>tot</sub></td>   <td>=</td>  <td>e<sub>s</sub> + P<sub>dry</sub></td>
    </tr>
    <tr>
       <td>VMR'<sub>H2O</sub></td> <td>=</td>  <td>P'<sub>H2O</sub> / P'<sub>tot</sub></td>
    </tr>
    </table>

   \param    abs_vmrs      [1]     volume mixing ratios per tag group (input/output)
   \param    abs_t     [K]     temperature at pressure level      (input)
   \param    abs_p     [Pa]    pressure levels                    (input)
   \param    tgs       [1]     the list of tag groups             (input)

   \author Thomas Kuhn
   \date 2001-08-02
 */ 
void WaterVaporSaturationInClouds( // WS Input/Output
                                  Matrix&           abs_vmrs,  // manipulates this WS
                                  Vector&           abs_p, // manipulates this WS
                                  // WS Input
                                  const Vector&     abs_t, // constant
                                  const ArrayOfArrayOfSpeciesTag&  tgs  ) // constant
{

  // make sure that the VMR and pressure grid are the same
  assert ( abs_vmrs.ncols() == abs_p.nelem() );

  // The species lookup data
  extern const Array<SpeciesRecord> species_data;
  // cloud tag numbers:
  Index liquid_index = 1+tgs.nelem();
  Index ice_index    = 1+tgs.nelem();
  Index h2o_index[tgs.nelem()];

  // check size of input vectors.
  assert ( abs_t.nelem() == abs_p.nelem() ); 
  assert ( abs_vmrs.ncols()  == abs_p.nelem() ); 

  // find tags for clouds and for water vapor
  Index u = 0;
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
      String tag_name = species_data[tgs[i][0].Species()].Name();
      if (tag_name == "liquidcloud") // <=== liquid water clouds tag name
        {
          liquid_index = i;
        }
      if (tag_name == "icecloud") // <====== ice water clouds tag name
        {
          ice_index = i;
        }
      if (tag_name == "H2O") // <=========== water vapor tags to change VMR
        {
          h2o_index[u++] = i;
          //cout << "tag_name=" << tag_name << ",  tag=" << i << ",  u=" << u 
          //     << ",  h2o_index[u]=" << h2o_index[u] << "\n";
        }
    }


  // if no water vapor profile is in use do not go further
  if (u < 1)
    {
      out2 << "  WaterVaporSaturationInClouds: no H2O profile found to adjust for clouds.\n"
           << "  Therefore no saturation calculation is performed\n";
      return;
    }

  // ------------------------< saturation over liquid water >------------------------
  // modify the water vapor profiles for saturation over liquid water in a liquid water cloud
  if ( (liquid_index >= 0) && (liquid_index < tgs.nelem()) )
    {
      // sauration over liquid water 
      for (Index uu=0; uu<u; ++uu)              // --- loop over all H2O tags
        {
          for (Index i=0; i<abs_vmrs.ncols() ; ++i) // --- loop over altitude grid
            {
              if (abs_vmrs(liquid_index,i) > 0.000) // --- cloud present or not?
                {
                  // dry air part of the total pressure
                  Numeric p_dry         =  abs_p[i] * ( 1.000e0 - abs_vmrs(h2o_index[uu],i) );
                  // water vapor saturation over liquid water 
                  Numeric e_s           =  WVSatPressureLiquidWater( abs_t[i] );
                  // new total pressure = dry air pressure +  water vapor saturation pressure
                  abs_p[i]              =  e_s + p_dry;
                  // new water vapor volume mixing ratio
                  abs_vmrs(h2o_index[uu],i) =  e_s / abs_p[i];
                  // check if VMR has strange value
                  if ( (abs_vmrs(h2o_index[uu],i) < 0.000e0) || (abs_vmrs(h2o_index[uu],i) > 1.000e0) ) 
                    {
                      ostringstream os;
                      os << "WaterVaporSaturationInClouds: The water vapor VMR value " 
                         << abs_vmrs(h2o_index[uu],i) << "\n"
                         << " looks strange after setting it to the saturation pressure over liquid water.";
                      throw runtime_error(os.str());
                      return;
                    }
                }
            }
        }
    }


  // ------------------------< saturation over ice water >------------------------
  // modify the water vapor profiles for saturation over ice water in a ice water cloud
  if ( (ice_index >= 0) && (ice_index < tgs.nelem()) )
    {
      for (Index uu=0; uu<u; ++uu)              // --- loop over all H2O tags
        {
          for (Index i=0; i<abs_vmrs.ncols() ; ++i) // --- loop over altitude grid
            {
              if (abs_vmrs(ice_index,i) > 0.000)    // --- cloud present or not?
                {
                  // dry air part of the total pressure
                  Numeric p_dry         =  abs_p[i] * ( 1.000e0 - abs_vmrs(h2o_index[uu],i) );
                  // water vapor saturation over ice
                  Numeric e_s           =  WVSatPressureIce( abs_t[i] );
                  // new total pressure = dry air pressure +  water vapor saturation pressure
                  abs_p[i]              =  e_s + p_dry;
                  // new water vapor volume mixing ratio
                  abs_vmrs(h2o_index[uu],i) =  e_s / abs_p[i];
                  // check if VMR has strange value
                  if ( (abs_vmrs(h2o_index[uu],i) < 0.000e0)  || (abs_vmrs(h2o_index[uu],i) > 1.000e0) )
                    {
                      ostringstream os;
                      os << "WaterVaporSaturationInClouds: The water vapor VMR value " 
                         << abs_vmrs(h2o_index[uu],i) << "\n"
                         << " looks strange after setting it to the saturation pressure over ice.";
                      throw runtime_error(os.str());
                      return;
                    }
                }
            }
        }
    }
  
  return;
}



/** Interpolate atmospheric quantities from their individual grids to
    the common abs_p grid. 

    See also arts -d online documentation.

    This function does the following:
    1. Interpolation of temperature and altitude
    2. Interpolation of VMR profiles
    3. Saturation adjustment VMR profiles of H2O tags in clouds

    Step 3 is only carried out if keyword CloudSatWV is set to "yes".
 */
void AtmFromRaw(// WS Output:
                  Vector&        abs_t,
                  Vector&        z_abs,
                  Matrix&        abs_vmrs,
                  // WS Input:      
                  const ArrayOfArrayOfSpeciesTag&       tgs,
                  const Vector&          abs_p,
                  const Matrix&          raw_ptz,
                  const ArrayOfMatrix&   raw_vmrs)
{
  
  //---------------< 1. Interpolation of temperature and altitude >---------------
  {  
    // Safety check: Make sure that raw_ptz really is a [x,3] matrix:
    if ( 3 != raw_ptz.ncols() )
      {
        ostringstream os;
        os << "The variable raw_ptz does not have the right dimensions,\n"
           << "ncols() should be 3, but is actually "<< raw_ptz.ncols();
        throw runtime_error(os.str());
      }


    // Contents of raw_ptz: p_raw is column 1,  tz_raw is column 2-3.

    // Interpolate tz_raw to abs_p grid.
    // The reason why we take tz_raw as a matrix is that the
    // interpolation can then be done simultaneously, hence slightly
    // more efficient.

    // For the interpolated profiles:
    Matrix tz_intp( 2, abs_p.nelem() );

    interpp( tz_intp,
             raw_ptz(Range(joker),0),
             transpose(raw_ptz(Range(joker),Range(1,joker))),
             abs_p );

    // The first Matpack expression selects the first column of
    // raw_ptz as a vector. The second Matpack expression gives the
    // transpose of the last two columns of raw_ptz. The function
    // interpp can be called with these selections directly.

    // Extract abs_t:
    abs_t.resize( tz_intp.ncols() );
    abs_t = tz_intp(0,Range(joker));    // Matpack can copy the first row of
                                        // tz_intp to abs_t like this. But
                                        // abs_t has to have the right size!

    // Extract z_abs:
    z_abs.resize( tz_intp.ncols() );
    z_abs = tz_intp(1,Range(joker));    // Matpack can copy the second row of
                                        // tz_intp to abs_t like this. But
                                        // abs_t has to have the right size!
  }

  //---------------< 2. Interpolation of VMR profiles >---------------
  {
    // check size of input String vectors.
    assert ( tgs.nelem() == raw_vmrs.nelem() ); 

    // Make abs_vmrs the right size:
    abs_vmrs.resize( raw_vmrs.nelem(), abs_p.nelem() );
    
    // For sure, we need to loop through all VMR profiles:
    for ( Index j=0; j<raw_vmrs.nelem(); ++j )
      {

        // Get a reference to the profile we are concerned with:
        const Matrix& raw = raw_vmrs[j];

        // Raw should be a matrix with dimension [x,2], the first column
        // is the raw pressure grid, the second column the VMR values.
      
        // Safety check to ensure this:
        if ( 2 != raw.ncols() )
          {
            ostringstream os;
            os << "The variable raw_vmrs("
               << j
               << ") does not have the right dimensions,\n"
               << "ncols() should be 2, but is actually "<< raw.ncols();
            throw runtime_error(os.str());
          }

        // Interpolate the profile to the predefined pressure grid:
        //String tag_name = species_data[tgs[j][0].Species()].Name(); // name of the tag
        // if ( (tag_name == "liquidcloud") || (tag_name == "icecloud") )
        //   {
        //     // Interpolate linearly the cloud profiles
        //     interpp_cloud( abs_vmrs(j,Range(joker)),
        //                 raw(Range(joker),0),
        //                 raw(Range(joker),1),
        //                 abs_p );
            //  out3 << "This VMR: " << abs_vmrs(j,Range(joker)) << "\n";
        //   }
        // else
        //   {
            // Interpolate VMRs:
            interpp( abs_vmrs(j,Range(joker)),
                     raw(Range(joker),0),
                     raw(Range(joker),1),
                     abs_p );
            // out3 << "This VMR: " << abs_vmrs(j,Range(joker)) << "\n";
        //   }
        // The calls to interpp_cloud and inerpp contain some nice
        // Matpack  features:
        // 1. abs_vmrs(j,Range(joker)) selects the jth row of abs_vmrs.
        // 2. raw(Range(joker),0) and raw(Range(joker),1) select the
        // first and second column of raw. We don't need transpose
        // here, since the selected objects are vectors. 
        //
        // Note that you can call the interpolation functions directly
        // with these selections. 
      }
  }
}


/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-18
*/
void hseSet(
          Vector&    hse,
    const Numeric&   pref,
    const Numeric&   zref,
    const Numeric&   g0,
    const Index&     niter )
{
  hse.resize( 5 );
  
  hse[0] = 1;
  hse[1] = pref;
  hse[2] = zref;
  hse[3] = g0;
  hse[4] = Numeric( niter );
}

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Axel von Engeln
   \date   2003-07-23
*/
void hseSetFromLatitude(
                        Vector&    hse,
                const Numeric&   pref,
                const Numeric&   zref,
                const Numeric&   latitude,
                const Index&     niter )
{

  hse.resize( 5 );
  
  hse[0] = 1;
  hse[1] = pref;
  hse[2] = zref;
  hse[3] = g_of_lat(latitude);
  hse[4] = Numeric( niter );
}

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Axel von Engeln
   \date   2003-07-24
*/
void hseSetFromLatitudeIndex(
                              Vector&    hse,
                      const Vector&    abs_p,
                      const Vector&    z_abs,
                      const Numeric&   latitude,
                      const Index&     index,
                      const Index&     niter )
{

  /* check index range */
  check_if_in_range( 0, abs_p.nelem()-1, index, "index" );

  hse.resize( 5 );
  
  hse[0] = 1;
  hse[1] = abs_p[index];
  hse[2] = z_abs[index];
  hse[3] = g_of_lat(latitude);
  hse[4] = Numeric( niter );
}


/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-18
*/
void hseFromBottom(
          Vector&    hse,
    const Vector&    abs_p,
    const Vector&    z_abs,
    const Numeric&   g0,
    const Index&       niter )
{
  hseSet( hse, abs_p[0], z_abs[0], g0, niter );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-18
*/
void hseOff(
          Vector&    hse )
{
  hse.resize( 1 );
  hse[0] = 0;
}



// Algorithm based on equations from my (PE) thesis (page 274) and the book
// Meteorology today for scientists and engineers by R.B. Stull (pages 9-10). 
//
/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-19
*/
void hseCalc(
          Vector&    z_abs,
    const Vector&    abs_p,
    const Vector&    abs_t,
    const Vector&    abs_h2o,
    const Numeric&   r_geoid,   
    const Vector&    hse )
{
  check_if_bool( static_cast<Index>(hse[0]), 
                                        "the HSE flag (first element of hse)");
  
  if ( hse[0] )
  {
    if ( hse.nelem() != 5 )
      throw runtime_error("The length of the *hse* vector must be 5.");

    const Index   np = abs_p.nelem();
          Index   i;                     // altitude index
          Numeric  g;                     // gravitational acceleration
          Numeric  r;                     // water mixing ratio in gram/gram
          Numeric  tv;                    // virtual temperature
          Numeric  dz;                    // step geometrical altitude
          Vector   ztmp(np);              // temporary storage for z_abs
  
    // Pick out values from hse
    const Numeric   pref  = hse[1];
    const Numeric   zref  = hse[2];
    const Numeric   g0    = hse[3];
    const Index     niter = Index( hse[4] );
  
    check_lengths( z_abs, "z_abs", abs_t, "abs_t" );  
    check_lengths( z_abs, "z_abs", abs_h2o, "abs_h2o" );  
    if ( niter < 1 )
      throw runtime_error("The number of iterations must be > 0.");
  
    for ( Index iter=0; iter<niter; iter++ )
    {
      // Init ztmp
      ztmp[0] = z_abs[0];
  
      // Calculate new altitudes (relative z_abs(1)) and store in ztmp
      for ( i=0; i<np-1; i++ )
      {
        // Calculate g 
        g  = ( g_of_z(r_geoid,g0,z_abs[i]) + 
               g_of_z(r_geoid,g0,z_abs[i+1]) ) / 2.0;
  
        // Calculate weight mixing ratio for water assuming constant average
        // molecular weight of the air
        r  = 18/28.96 * (abs_h2o[i]+abs_h2o[i+1])/2;
  
        // The virtual temperature (no liquid water)
        tv = (1+0.61*r) * (abs_t[i]+abs_t[i+1])/2;
  
        // The change in vertical altitude from i to i+1 
        dz = 287.053*tv/g * log( abs_p[i]/abs_p[i+1] );
        ztmp[i+1] = ztmp[i] + dz;
      }
  
      // Match the altitude of the reference point
      dz = interpp( abs_p, ztmp, pref ) - zref;

      //  z_abs = ztmp - dz;
      z_abs = ztmp;
      z_abs -= dz;              
    }
  }
}


/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Carlos Jimenez 
   \date   2001-08-14
*/
void vmrsScale(
               Matrix&                abs_vmrs,
               const ArrayOfArrayOfSpeciesTag&       tgs,
               const ArrayOfString&   scaltgs,
               const Vector&          scalfac)
{
  Index                            itag;
  ArrayOfIndex                     tagindex;      

  if ( scalfac.nelem() != scaltgs.nelem()  )
    throw runtime_error("vmrScale: Number of tgs and fac are different!");
  
  get_tagindex_for_Strings( tagindex, tgs, scaltgs );

  const Index   n = tagindex.nelem();

  for ( itag=0; itag<n; itag++ )
    {
      //out2 << scalfac[itag] << ".\n";
      abs_vmrs(tagindex[itag],Range(joker)) *= scalfac[itag]; 
      // Matpack can multiply all elements of a vector with a constant
      // factor like this. In this case the vector is the selected row
      // of Matrix abs_vmrs.
  
      //out2 << abs_vmrs(tagindex[itag],Range(joker)) << ".\n";
    }
}

#endif


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_h2oSet(Vector&          abs_h2o,
                const ArrayOfArrayOfSpeciesTag& abs_species,
                const Matrix&    abs_vmrs)
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
               const   ArrayOfArrayOfSpeciesTag& abs_species,
               const   Matrix&    abs_vmrs)
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
                            const Tensor4& vmr_field
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
                  const ArrayOfVector&            abs_cont_parameters)
{
  // Dimension checks are performed in the executed functions

  // allocate local variable to hold the cross sections per tag group
  ArrayOfMatrix abs_xsec_per_species;
  
  abs_xsec_per_speciesInit( abs_xsec_per_species, tgs, f_grid, abs_p );

  abs_xsec_per_speciesAddLines( abs_xsec_per_species,
                       tgs,
                       f_grid,
                       abs_p,
                       abs_t,
                       abs_h2o,
                       abs_vmrs,
                       abs_lines_per_species,
                       abs_lineshape );

  abs_xsec_per_speciesAddConts( abs_xsec_per_species,
                       tgs,
                       f_grid,
                       abs_p,
                       abs_t,
                       abs_n2,
                       abs_h2o,
                       abs_vmrs,
                       abs_cont_names,
                       abs_cont_parameters,
                       abs_cont_models);

  abs_coefCalcFromXsec(abs_coef,
                       abs_coef_per_species,
                       abs_xsec_per_species,
                       abs_vmrs,
                       abs_p,
                       abs_t);

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
                            const ArrayOfString&   abs_cont_names,
                            const ArrayOfString&   abs_cont_models,
                            const ArrayOfVector&   abs_cont_parameters)
{
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

  out2 << "  Number of tag groups to do: " << tgs.nelem() << "\n";

  // Loop over all species:
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
      out2 << "  Doing tag group " << i << ".\n";

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

      abs_xsec_per_speciesInit( abs_xsec_per_species, this_tg, f_grid, abs_p );

      abs_xsec_per_speciesAddLines( abs_xsec_per_species,
                           this_tg,
                           f_grid,
                           abs_p,
                           abs_t,
                           abs_h2o,
                           this_vmr,
                           these_lines,
                           this_lineshape );

      abs_xsec_per_speciesAddConts( abs_xsec_per_species,
                           this_tg,
                           f_grid,
                           abs_p,
                           abs_t,
                           abs_n2,
                           abs_h2o,
                           this_vmr,
                           abs_cont_names,
                           abs_cont_parameters,
                           abs_cont_models);

      abs_coefCalcFromXsec(this_abs,
                           this_abs_coef_per_species,
                           abs_xsec_per_species,
                           this_vmr,
                           abs_p,
                           abs_t);

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
                          const Vector&        abs_t)
{
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

  out2 << "  Computing abs_coef and abs_coef_per_species from abs_xsec_per_species.\n";

  // Loop through all tag groups
  for ( Index i=0; i<abs_xsec_per_species.nelem(); ++i )
    {
      out2 << "  Tag group " << i << "\n";

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
                     const Vector&    abs_p
                     )
{
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

  out3 << "  Initialized abs_xsec_per_species.\n"
       << "  Number of frequencies        : " << f_grid.nelem() << "\n"
       << "  Number of pressure levels    : " << abs_p.nelem() << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddLines(// WS Output:
                         ArrayOfMatrix&                   abs_xsec_per_species,
                         // WS Input:             
                         const ArrayOfArrayOfSpeciesTag&                 tgs,
                         const Vector&                    f_grid,
                         const Vector&                    abs_p,
                         const Vector&                    abs_t,
                         const Vector&                    abs_h2o,
                         const Matrix&                    abs_vmrs,
                         const ArrayOfArrayOfLineRecord&  abs_lines_per_species,
                         const ArrayOfLineshapeSpec&      abs_lineshape)
{
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
      out3 << "  Tag group " << i
           << " (" << get_tag_group_name(tgs[i]) << "): ";
      
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
      
      // Skip the call to abs_xsec_per_species if the line list is empty.
      if ( 0 < ll.nelem() )
        {

          // Get the name of the species. The member function name of a
          // LineRecord returns the full name (species + isotope). So
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

          out3 << ll.nelem() << " transitions\n";
          xsec_species( abs_xsec_per_species[i],
                        f_grid,
                        abs_p,
                        abs_t,
                        abs_h2o,
                        abs_vmrs(i,Range(joker)),
                        ll,
                        ls.Ind_ls(),
                        ls.Ind_lsn(),
                        ls.Cutoff());
          // Note that we call xsec_species with a row of abs_vmrs,
          // selected by the above Matpack expression. This is
          // possible, because xsec_species is using Views.
        }
      else
        {
          out3 << ll.nelem() << " transitions, skipping\n";
        }
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddConts(// WS Output:
                         ArrayOfMatrix&                   abs_xsec_per_species,
                         // WS Input:             
                         const ArrayOfArrayOfSpeciesTag&                 tgs,
                         const Vector&                    f_grid,
                         const Vector&                    abs_p,
                         const Vector&                    abs_t,
                         const Vector&                    abs_n2,
                         const Vector&                    abs_h2o,
                         const Matrix&                    abs_vmrs,
                         const ArrayOfString&             abs_cont_names,
                         const ArrayOfVector&             abs_cont_parameters,
                         const ArrayOfString&             abs_cont_models )
{
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
      for (Index i=0; i < abs_cont_names.nelem(); ++i) 
        {
          cout << "abs_xsec_per_speciesAddConts: " << i << " name : " << abs_cont_names[i] << "\n";
        }
      for (Index i=0; i < abs_cont_parameters.nelem(); ++i) 
        {
          cout << "abs_xsec_per_speciesAddConts: " << i << " param: " << abs_cont_parameters[i] << "\n";
        }
      for (Index i=0; i < abs_cont_models.nelem(); ++i) 
        {
          cout << "abs_xsec_per_speciesAddConts: " << i << " option: " << abs_cont_models[i] << "\n";
        }
      ostringstream os;
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
          // tag that means `all isotopes', because this should not
          // include continuum. For such tags, tag.Isotope() will
          // return the number of isotopes (i.e., one more than the
          // allowed index range).
          if ( tgs[i][s].Isotope() <
               species_data[tgs[i][s].Species()].Isotope().nelem() )
            {
              // If we get here, it means that the tag describes a
              // specific isotope. Could be a continuum tag!
                
              // The if clause below checks whether the abundance of this tag
              // is negative. (Negative abundance marks continuum tags.)
              // It does the following:
              //
              // 1. species_data contains the lookup table of species
              //          specific data. We need the right entry in this
              //          table. The index of this is obtained by calling member function
              //          Species() on the tag. Thus we have:
              //          species_data[tgs[i][s].Species()].
              //
              // 2. Member function Isotope() on the above biest gives
              //    us the array of isotope specific data. This we have
              //    to subscribe with the right isotope index, which we
              //    get by calling member function Isotope on the
              //    tag. Thus we have:
              //    Isotope()[tgs[i][s].Isotope()]
              //
              // 3. Finally, from the isotope data we need to get the mixing
              //    ratio by calling the member function Abundance().
              if ( 0 >
                   species_data[tgs[i][s].Species()].Isotope()[tgs[i][s].Isotope()].Abundance() )
                {
                  // We have identified a continuum tag!

                  // Get only the continuum name. The full tag name is something like:
                  // H2O-HITRAN96Self-*-*. We want only the `H2O-HITRAN96Self' part:
                  const String name =
                    species_data[tgs[i][s].Species()].Name() + "-"
                    + species_data[tgs[i][s].Species()].Isotope()[tgs[i][s].Isotope()].Name();
  
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
                  
                  out3 << "  Adding " << name
                       << " to tag group " << i << ".\n";

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
                                      abs_vmrs(i,Range(joker)) );
                  // Calling this function with a row of Matrix abs_vmrs
                  // is possible because it uses Views.
                }
            }
        }
    }

}

#if 0

/** Reduces the size of abs_coef_per_species.  Only absorption coefficients for
    which weighting functions are calculated are kept in memory.
    
    \retval abs_coef_per_species absorption coefficients
    \param  tgs        all selected tag groups
    \param  wfs_tgs    the tag groups for which we want weighting functions.

    \author Axel von Engeln and Stefan Buehler */
void abs_coef_per_speciesReduce(// WS Output:
                      ArrayOfMatrix&         abs_coef_per_species,
                      // WS Input:
                      const ArrayOfArrayOfSpeciesTag&       tgs,
                      const ArrayOfArrayOfSpeciesTag&       wfs_tgs)
{

  // Make a safety check that the dimensions of tgs and
  // abs_coef_per_species are the same (could happen that we call this workspace
  // method twice by accident).
  if ( abs_coef_per_species.nelem()!=tgs.nelem() )
    throw(runtime_error("The variables abs_coef_per_species and tgs must\n"
                        "have the same dimension."));

  // Erase could be very inefficient in this case, since elements
  // behind the erased one are copied in order to fill the
  // gap. Therefore, we will construct a new abs_coef_per_species, and finally
  // use it to replace the old one.
  ArrayOfMatrix abs_coef_per_species_out( wfs_tgs.nelem() );

  // Go through the weighting function tag groups:
  for ( Index i=0; i<wfs_tgs.nelem(); ++i )
    {
      // Index to the elements of wfs_tgs in tgs:
      Index n;
      get_tag_group_index_for_tag_group( n, tgs, wfs_tgs[i] );

      abs_coef_per_species_out[i].resize( abs_coef_per_species[n].nrows(), abs_coef_per_species[n].ncols() );
      abs_coef_per_species_out[i] = abs_coef_per_species[n]; // Matpack can copy the contents of
                                         // matrices like this. The dimensions
                                         // must be the same! 
    }  

  // Copy the generated matrices back to abs_coef_per_species
  abs_coef_per_species.resize( wfs_tgs.nelem() );
  abs_coef_per_species = abs_coef_per_species_out;  // FIXME: It should be checked whether
                                // this works correctly. Does my Array
                                // implementation work as it should? 
}



//======================================================================
//             Methods related to refraction
//======================================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-19
*/
void refrSet( 
              Index&    refr,
              Index&    refr_lfac,
              String&   refr_model,
        const Index&    on,
        const String&   model,
        const Index&    lfac )
{
  check_if_bool( on, "on" );
  if ( lfac < 1 )
    throw runtime_error("The length factor cannot be smaller than 1.");      

  refr       = on;
  refr_lfac  = lfac;
  refr_model = model;
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-22
*/
void refrOff( 
              Index&    refr,
              Index&    refr_lfac,
              String&   refr_model )
{
  refrSet( refr, refr_lfac, refr_model, 0, "", 1 );
  refr_lfac = 0;
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-19
*/
void refrCalc (
                    Vector&   refr_index,
              const Vector&   abs_p,
              const Vector&   abs_t,
              const Vector&   abs_h2o,
              const Index&    refr,
              const String&   refr_model )
{
  check_if_bool( refr, "refr" );

  if ( refr == 0 )
    refr_index.resize( 0 );

  else
  {
    if ( refr_model == "Unity" )
    {
      const Index n = abs_p.nelem();
      refr_index.resize( n );
      refr_index = 1.0;
    }
    
    else if ( refr_model == "Boudouris" )
      refr_index_Boudouris( refr_index, abs_p, abs_t, abs_h2o );

    else if ( refr_model == "BoudourisDryAir" )
      refr_index_BoudourisDryAir( refr_index, abs_p, abs_t );

    else
    {
      ostringstream os;
      os << "Unknown refraction parameterization: " << refr_model << "\n"
         << "Existing parameterizations are: \n"
         << "Unity\n"
         << "Boudouris\n"
         << "BoudourisDryAir";
      throw runtime_error(os.str());
    }
  }
}

#endif

//======================================================================
//             Methods related to continua
//======================================================================

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_cont_descriptionInit(// WS Output:
                          ArrayOfString& names,
                          ArrayOfString& options, 
                          ArrayOfVector& parameters)
{
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
                       const Vector& userparameters)
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
void abs_scalar_gasFromAbsCoef(// WS Output:
                               Matrix&       abs_scalar_gas,
                               // WS Input:
                               const ArrayOfMatrix& abs_coef_per_species)
{
  // abs_scalar_gas has format [f_grid, abs_species].
  // abs_coef_per_species has format ArrayOfMatrix (over species),
  // where for each species the matrix has format [f_grid, abs_p].

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
  
  abs_scalar_gas.resize(n_f,n_species);

  // Loop species:
  for ( Index si=0; si<n_species; ++si )
    abs_scalar_gas(Range(joker),si) = abs_coef_per_species[si](Range(joker),0);

}

