/* Copyright (C) 2000, 2001 Stefan Buehler   <sbuehler@uni-bremen.de>
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

   \author Stefan Buehler
   \date   2001-03-12
*/

#include <math.h>
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
#include "atm_funcs.h"
#include "continua.h"
#include "make_vector.h"

void linesReadFromHitran(// WS Output:
                         ArrayOfLineRecord& lines,
                          // Control Parameters:
                         const String& filename,
                         const Numeric& fmin,
                         const Numeric& fmax)
{
  ifstream is;

  out2 << "  Reading file: " << filename << '\n';
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
		lines.push_back(lr);
	      else
		go_on = false;
	    }
	}
    }
  out2 << "  Read " << lines.nelem() << " lines.\n";
}


void linesReadFromMytran2(// WS Output:
			  ArrayOfLineRecord& lines,
                          // Control Parameters:
			  const String& filename,
			  const Numeric& fmin,
			  const Numeric& fmax)
{
  ifstream is;

  out2 << "  Reading file: " << filename << '\n';
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
	      lines.push_back(lr);
	}
    }
  out2 << "  Read " << lines.nelem() << " lines.\n";
}

void linesReadFromJpl(// WS Output:
		      ArrayOfLineRecord& lines,
		      // Control Parameters:
		      const String& filename,
		      const Numeric& fmin,
		      const Numeric& fmax)
{
  ifstream is;

  out2 << "  Reading file: " << filename << '\n';
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
		lines.push_back(lr);
	      else
		go_on = false;
	    }
	}
    }
  out2 << "  Read " << lines.nelem() << " lines.\n";
}


void linesReadFromArts(// WS Output:
		       ArrayOfLineRecord& lines,
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

  out2 << "  Reading file: " << filename << '\n';
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
	  if ( fmin <= lr.F() )
	    {
	      // lines are not necessarily frequency sorted 
	      if ( fmin <= lr.F() )
		if ( lr.F() <= fmax )
		  {
		    lines.push_back(lr);
		    //		    out3 << lr << "\n";
		  }
	    }
	}
    }
  out2 << "  Read " << lines.nelem() << " lines.\n";
}

void linesElowToJoule(// WS Output:
		      ArrayOfLineRecord& lines )
{
  for ( Index i=0; i<lines.nelem(); ++i )
    lines[i].melow = wavenumber_to_joule(lines[i].melow); 
}

/**
  This method can read lines from different line catalogues. For each
  tag group, you can specify which catalogue to use. Because the
  method creates lines_per_tg directly, it replaces for example the
  following two method calls:
  - linesReadFromHitran
  - lines_per_tgCreateFromLines

  This method needs as input WSVs the list of tag groups. Keyword
  parameters must specify the names of the catalogue files to use and
  the matching formats. Names can be anything, formats can currently
  be HITRAN96, MYTRAN2, JPL, or ARTS. Furthermore, keyword parameters
  have to specify minimum and maximum frequency for each tag group. To
  safe typing, if there are less elements in the keyword parameters
  than there are tag groups, the last parameters are applied to all
  following tag groups.

  Example usage:

  lines_per_tgReadFromCatalogues{
  	filenames = [ "../data/cat1.dat", "../data/cat2.dat" ]
	formats   = [ "MYTRAN2",          "HITRAN96"         ]
	fmin      = [ 0,                  0                  ]
	fmax      = [ 2000e9,             100e9              ]
  }

  In this example, lines for the first tag group will be taken from
  cat1, lines for all other tag groups will be taken from cat2.

  This methods allows you for example to use a special line file just
  for water vapor lines. This could be the improved water vapor line
  file generated by Thomas Kuhn.

  Catalogues are only read once, even if several tag groups have the
  same catalogue. However, in that case the frequency ranges MUST be
  the same. (If you want to do fine-tuning of the frequency ranges,
  you can do this inside the tag definitions, e.g., "H2O-*-0-2000e9".)

  This function uses the various reading routines
  (linesReadFromHitran, etc.), as well as
  lines_per_tgCreateFromLines. 
 
  \author Stefan Buehler
  \date 2000-01-19 */
void lines_per_tgReadFromCatalogues(// WS Output:
				    ArrayOfArrayOfLineRecord& lines_per_tg,
				    // WS Input:
				    const TagGroups& tgs,
                                    // Control Parameters:
                                    const ArrayOfString& filenames,
                                    const ArrayOfString& formats,
                                    const Vector& fmin,
                                    const Vector& fmax)
{
  const Index n_tg   = tgs.nelem();	// # tag groups
  const Index n_cat = filenames.nelem();	// # tag Catalogues

  // Check that dimensions of the keyword parameters are consistent
  // (must all be the same). 

  if ( n_cat != formats.nelem() ||
       n_cat != fmin.nelem() ||
       n_cat != fmax.nelem() )
    {
      ostringstream os;
      os << "lines_per_tgReadFromCatalogues: All keyword\n"
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
      os << "lines_per_tgReadFromCatalogues: You specified more\n"
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
      os << "lines_per_tgReadFromCatalogues: You must have at\n"
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

  MakeArray<String> 	 real_filenames ( filenames[0] );
  MakeArray<String> 	 real_formats   ( formats[0]   );
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
		   fmax[i]    != real_fmax[that_cat] 	   )
		{
		  ostringstream os;
		  os << "lines_per_tgReadFromCatalogues: If you specify the\n"
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

	      last_cat = i;	// assign remainder of lines to this
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

  // Make lines_per_tg the right size:
  lines_per_tg.resize( tgs.nelem() );

  // Loop through the catalogues to read:
  for ( Index i=0; i<n_real_cat; ++i )
    {
      ArrayOfLineRecord   lines;

      // Read catalogue:

      if ( "HITRAN96"==real_formats[i] )
	{
	  linesReadFromHitran( lines, real_filenames[i], real_fmin[i], real_fmax[i] );
	}
      else if ( "MYTRAN2"==real_formats[i] )
	{
	  linesReadFromMytran2( lines, real_filenames[i], real_fmin[i], real_fmax[i] );
	}
      else if ( "JPL"==real_formats[i] )
	{
	  linesReadFromJpl( lines, real_filenames[i], real_fmin[i], real_fmax[i] );
	}
      else if ( "ARTS"==real_formats[i] )
	{
	  linesReadFromArts( lines, real_filenames[i], real_fmin[i], real_fmax[i] );
	}
      else
	{
	  ostringstream os;
	  os << "lines_per_tgReadFromCatalogues: You specified the\n"
             << "format `" << real_formats[i] << "', which is unknown.\n"
	     << "Allowd formats are: HITRAN96, MYTRAN2, JPL, ARTS.";
	  throw runtime_error(os.str());
	}

      // We need to make subset tgs for the groups that should
      // be read from this catalogue.
      TagGroups  these_tgs(real_tgs[i].nelem());
      for ( Index s=0; s<real_tgs[i].nelem(); ++s )
	{
	  these_tgs[s] = tgs[real_tgs[i][s]];
	}

      // Create these_lines_per_tg:
      ArrayOfArrayOfLineRecord these_lines_per_tg;
      lines_per_tgCreateFromLines( these_lines_per_tg, lines, these_tgs );

      // Put these lines in the right place in lines_per_tg:
      for ( Index s=0; s<real_tgs[i].nelem(); ++s )
	{
	  lines_per_tg[real_tgs[i][s]] = these_lines_per_tg[s];
	}
    }
}


void lines_per_tgCreateFromLines(// WS Output:
                                  ArrayOfArrayOfLineRecord& lines_per_tg,
                                  // WS Input:
                                  const ArrayOfLineRecord&   lines,
                                  const TagGroups&           tgs)
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
      
  // Make lines_per_tg the right size:
  lines_per_tg = ArrayOfArrayOfLineRecord(tgs.nelem());

  // Unfortunately, MTL conatains a bug that leads to all elements of
  // the outer Array of an Array<Array>> pointing to the same data
  // after creation. So we need to fix this explicitly:
  for ( Index i=0; i<lines_per_tg.nelem(); ++i )
    lines_per_tg[i] = ArrayOfLineRecord();

  // Loop all lines in the input line list:
  for ( Index i=0; i<lines.nelem(); ++i )
    {
      // Get a convenient reference to the current line:
      const LineRecord& this_line = lines[i];

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
	      const OneTag& this_tag = tgs[j][k];

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
	  lines_per_tg[j-1].push_back(this_line);

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
	      out0 << "Your tags include other lines of species "
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

	out3 << ": " << lines_per_tg[i].nelem() << " lines\n";
    }

}

/** Adds mirror lines at negative frequencies to the lines_per_tg.
    For each line at frequency +f in lines_per_tg a corresponding entry at
    frequency -f is added to lines_per_tg.

    \retval lines_per_tg The array of arrays of lines for each tag group.
    
    \author Axel von Engeln and Stefan Buehler */
void lines_per_tgAddMirrorLines(// WS Output:
                                ArrayOfArrayOfLineRecord& lines_per_tg)
{
  // We will simply append the mirror lines after the original
  // lines. This way we don't have to make a backup copy of the
  // original lines. 

  for ( Index i=0; i<lines_per_tg.nelem(); ++i )
    {
      // Get a reference to the current list of lines to save typing:
      ArrayOfLineRecord& ll = lines_per_tg[i];
      
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
	    //	    cout << "Adding ML at f = " << dummy.F() << "\n";
	    ll.push_back(dummy);
	  }
      }
    }

}

/** Removes all lines outside the defined lineshape cutoff frequency
    from the lines_per_tg, in order to save computation time.
    It should be particularly useful to call this method after
    lines_per_tgAddMirrorLines.

    \retval lines_per_tg the old and newly compacted line list
    \param  lineshape the lineshape spceifications
    \param  f_mono the frequency grid

    \author Axel von Engeln and Stefan Buehler */
void lines_per_tgCompact(// WS Output:
			 ArrayOfArrayOfLineRecord& lines_per_tg,
			 // WS Input:
			 const ArrayOfLineshapeSpec& lineshape,
			 const Vector& f_mono)
{

  // Make sure lines_per_tg and lineshape have the same dimension:
  if ( lines_per_tg.nelem() != lineshape.nelem() ) 
    {
      ostringstream os;
      os << "Dimension of lines_per_tg does\n"
	 << "not match that of lineshape.";
      throw runtime_error(os.str());
    }
  
  // Make sure that the frequency grid is properly sorted:
  for ( Index s=0; s<f_mono.nelem()-1; ++s )
    {
      if ( f_mono[s+1] <= f_mono[s] )
	{
	  ostringstream os;
	  os << "The frequency grid f_mono is not properly sorted.\n"
	     << "f_mono[" << s << "] = " << f_mono[s] << "\n"
	     << "f_mono[" << s+1 << "] = " << f_mono[s+1];
	  throw runtime_error(os.str());
	}
    }

  // Cycle through all tag groups:
  for ( Index i=0; i<lines_per_tg.nelem(); ++i )
    {
      // Get cutoff frequency of this tag group:
      Numeric cutoff = lineshape[i].Cutoff();

      // Check whether cutoff is defined:
      if ( cutoff != -1)
	{
	  // Get a reference to the current list of lines to save typing:
	  ArrayOfLineRecord& ll = lines_per_tg[i];

	  // Calculate the borders:
	  Numeric upp = f_mono[f_mono.nelem()-1] + cutoff;
	  Numeric low = f_mono[0] - cutoff;

	  // Cycle through all lines within this tag group. 
	  for ( ArrayOfLineRecord::iterator j=ll.begin(); j<ll.end(); ++j )
	    {
	      // Center frequency:
	      const Numeric F0 = j->F();

	      if ( ( F0 < low) || ( F0 > upp) )
		{
		  j = ll.erase(j) - 1;
		  // We need the -1 here, otherwise due to the
		  // following increment we would miss the element
		  // behind the erased one, which is now at the
		  // position of the erased one.
		}
	    }
	}
    }
}


void linesWriteAscii(// WS Input:
		      const ArrayOfLineRecord& lines,
		      // Control Parameters:
		      const String& f)
{
  String filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const String basename;                       
      filename = basename+".lines.aa";
    }

  ofstream os;

  out2 << "  Writing file: " << filename << '\n';
  open_output_file(os, filename);

  write_lines_to_stream(os,lines);
}


void lines_per_tgWriteAscii(// WS Input:
			      const ArrayOfArrayOfLineRecord& lines_per_tg,
			      // Control Parameters:
			      const String& f)
{
  String filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const String basename;                       
      filename = basename+".lines_per_tg.aa";
    }

  ofstream os;

  out2 << "  Writing file: " << filename << '\n';
  open_output_file(os, filename);

  os << lines_per_tg.nelem() << '\n';

  for ( Index i=0; i<lines_per_tg.nelem(); ++i )
    {
      const ArrayOfLineRecord& lines = lines_per_tg[i];
      write_lines_to_stream( os, lines );
    }
}


void tgsDefine(// WS Output:
	       TagGroups& tgs,
	       // Control Parameters:
	       const ArrayOfString& tags)
{
  tgs.resize(tags.nelem());

  //cout << "Tags: " << tags << "\n";

  // Each element of the array of Strings tags defines one tag
  // group. Let's work through them one by one.
  for ( Index i=0; i<tags.nelem(); ++i )
    {
      // There can be a comma separated list of tag definitions, so we
      // need to break the String apart at the commas.
      ArrayOfString tag_def;

      bool go_on = true;
      String these_tags = tags[i];
      while (go_on)
	{
	  Index n = these_tags.find(',');
	  if ( n == these_tags.npos ) // npos indicates `not found'
	    {
	      // There are no more commas.
	      //	      cout << "these_tags: (" << these_tags << ")\n";
	      tag_def.push_back(these_tags);
	      go_on = false;
	    }
	  else
	    {
	      tag_def.push_back(these_tags.substr(0,n));
	      these_tags.erase(0,n+1);
	    }
	}

      // tag_def now holds the different tag Strings for this group.
      //      cout << "tag_def =\n" << tag_def << endl;

      // Set size to zero, in case the method has been called before.
      tgs[i].resize(0);

      for ( Index s=0; s<tag_def.nelem(); ++s )
	{
	  // Remove leading whitespace, if there is any:
	  while ( ' '  == tag_def[s][0] ||
		  '\t' == tag_def[s][0]    )	tag_def[s].erase(0,1);

	  OneTag this_tag(tag_def[s]);

	  // Safety check: For s>0 check that the tags belong to the same species.
	  if (s>0)
	    if ( tgs[i][0].Species() != this_tag.Species() )
	      throw runtime_error("Tags in a tag group must belong to the same species.");

	  tgs[i].push_back(this_tag);
	}
    }

  // Print list of tag groups to the most verbose output stream:
  out3 << "  Defined tag groups:";
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
      out3 << "\n  " << i << ":";
      for ( Index s=0; s<tgs[i].nelem(); ++s )
	{
	  out3 << " " << tgs[i][s].Name();
	}
    }
  out3 << '\n';

//  cout << endl << endl << tgs << endl;
}



void lineshapeDefine(// WS Output:
		     ArrayOfLineshapeSpec&    lineshape,
		     // WS Input:
		     const TagGroups&         tgs,
		     const String&            shape,
		     const String&            normalizationfactor,
		     const Numeric&           cutoff)
{
  // Make lineshape and normalization factor data visible:
  extern const Array<LineshapeRecord> lineshape_data;
  extern const Array<LineshapeNormRecord> lineshape_norm_data;


  // generate the right number of elements
  Index tag_sz = tgs.nelem();
  lineshape.resize(tag_sz);

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
      lineshape[i].SetInd_ls( found0 );
      lineshape[i].SetInd_lsn( found1 );
      lineshape[i].SetCutoff( cutoff );
    }
}

void lineshape_per_tgDefine(// WS Output:
			    ArrayOfLineshapeSpec& lineshape,
			    // WS Input:
			    const TagGroups&      tgs,
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
      os << "lineshape_per_tgDefine: number of elements does\n"
	 << "not match the number of tag groups defined.";
      throw runtime_error(os.str());
    }
      

  // generate the right number of elements
  lineshape.resize(tg_sz);

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
      lineshape[k].SetInd_ls( found0 );
      lineshape[k].SetInd_lsn( found1 );
      lineshape[k].SetCutoff( cutoff[k] );
    }
}



void raw_vmrsReadFromScenario(// WS Output:
                              ArrayOfMatrix&   raw_vmrs,
                              // WS Input:
                              const TagGroups&     tgs,
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
			   const TagGroups& tgs,
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
    liquid water/ice saturation VMR (=p_saturation/p_tot)

   \param    vmrs           volume mixing ratios per tag group (input/output)
   \param    t_abs          temperature at pressure level      (input)
   \param    p_abs          pressure levels                    (input)
   \param    tgs            the list of tag groups             (input)

   \author Thomas Kuhn
   \date 2001-08-02
 */ 
Index WVsatinClouds( Matrix&    vmrs,  // manipulates this array
		     const Vector&     t_abs, // constant
		     const Vector&     p_abs, // constant
		     const TagGroups&  tgs  ) // constant
{


  // The species lookup data
  extern const Array<SpeciesRecord> species_data;
  // cloud tag numbers:
  Index liquid_index = 1+tgs.nelem();
  Index ice_index    = 1+tgs.nelem();
  Index h2o_index[tgs.nelem()];

  // check size of input vectors.
  assert ( t_abs.nelem() == p_abs.nelem() ); 
  assert ( vmrs.ncols()  == p_abs.nelem() ); 

  // find tags for clouds and for water vapor
  Index u = 0;
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
      String tag_name = species_data[tgs[i][0].Species()].Name();
      if (tag_name == "liquidcloud")
	{
	  liquid_index = i;
	}
      if (tag_name == "icecloud")
	{
	  ice_index = i;
	}
      if (tag_name == "H2O")
	{
	  h2o_index[u++] = i;
	  //cout << "tag_name=" << tag_name << ",  tag=" << i << ",  u=" << u 
          //     << ",  h2o_index[u]=" << h2o_index[u] << "\n";
	}
    }

  // if no water vapor profile is in use do not go further
  if (u < 1)
    {
      cout << "WVsatinClouds: no H2O profile found to adjust for clouds..\n";
      return 0;
    }

  // modify the water vapor profiles for saturation over liquid water in a liquid water cloud
  if ( (liquid_index >= 0) && (liquid_index < tgs.nelem()) )
    {
      // sauration over liquid water 
      for (Index i=0; i<vmrs.ncols() ; ++i)
	{
	  if (vmrs(liquid_index,i) > 0.000)
	    {
	      for (Index uu=0; uu<u; ++uu)
		{
		  //cout << "uu=" << uu  << "  (" << p_abs.nelem() << ")" << "\n";
		  //cout << " h2o_index[uu]=" << h2o_index[uu] << ",  i=" << i << "\n";
		  //cout << "tag name=" << species_data[tgs[h2o_index[uu]][0].Species()].Name() << "\n";
		  //cout << "0 vmrs(h2o_index[uu],i)=" << vmrs(h2o_index[uu],i); 
		  vmrs(h2o_index[uu],i) = ( WVSatPressureLiquidWater( t_abs[i] ) / p_abs[i] );
		  if (vmrs(h2o_index[uu],i) < 0.000) return 1;
		  //cout << ",  1 vmrs(h2o_index[uu],i)=" << vmrs(h2o_index[uu],i) << "\n";
		  //cout << "T=" << t_abs[i] << "K,  ptot=" << p_abs[i] << "Pa,  psat="
                  //     << WVSatPressureLiquidWater( t_abs[i] ) << "\n";
		}
	    }
	}
    }

  // modify the water vapor profiles for saturation over ice water in a ice water cloud
  if ( (ice_index >= 0) && (ice_index < tgs.nelem()) )
    {
      for (Index i=0; i<vmrs.ncols() ; ++i)
	{
	  if (vmrs(ice_index,i) > 0.000)
	    {
	      // sauration over ice water 
	      for (Index uu=0; uu<u; ++uu)
		{
		  //cout << "uu=" << uu  << "  (" << p_abs.nelem() << ")" << "\n";
		  //cout << " h2o_index[uu]=" << h2o_index[uu] << ",  i=" << i << "\n";
		  //cout << "tag name=" << species_data[tgs[h2o_index[uu]][0].Species()].Name() << "\n";
		  vmrs(h2o_index[uu],i) = ( WVSatPressureIce( t_abs[i] ) / p_abs[i] );
		  if (vmrs(h2o_index[uu],i) < 0.000) return 1;
		}
	    }
	}
    }

  return 0;
} // end of WVsatinClouds --------------------------------------------------------------------

/** Interpolate atmospheric quantities from their individual grids to
    the common p_abs grid. 

    See also arts -d online documentation.

    This function does the following:
    1. Interpolation of temperature and altitude
    2. Interpolation of VMR profiles
    3. Saturation adjustment VMR profiles of H2O tags in clouds

    Step 3 is only carried out if keyword CloudSatWV is set to "yes".
 */
void AtmFromRaw(// WS Output:
		  Vector& 	 t_abs,
		  Vector& 	 z_abs,
		  Matrix&        vmrs,
		  // WS Input:      
		  const TagGroups&       tgs,
		  const Vector&  	 p_abs,
		  const Matrix&  	 raw_ptz,
		  const ArrayOfMatrix&   raw_vmrs,
		  // Control Parameters:
		  const String&          CloudSatWV)
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

    // Interpolate tz_raw to p_abs grid.
    // The reason why we take tz_raw as a matrix is that the
    // interpolation can then be done simultaneously, hence slightly
    // more efficient.

    // For the interpolated profiles:
    Matrix tz_intp( 2, p_abs.nelem() );

    interpp( tz_intp,
	     raw_ptz(Range(joker),0),
	     transpose(raw_ptz(Range(joker),Range(1,joker))),
	     p_abs );
    // The first Matpack expression selects the first column of
    // raw_ptz as a vector. The second Matpack expression gives the
    // transpose of the last two columns of raw_ptz. The function
    // interpp can be called with these selections directly.

    // Extract t_abs:
    t_abs.resize( tz_intp.ncols() );
    t_abs = tz_intp(0,Range(joker));	// Matpack can copy the first row of
					// tz_intp to t_abs like this. But
					// t_abs has to have the right size!

    // Extract z_abs:
    z_abs.resize( tz_intp.ncols() );
    z_abs = tz_intp(1,Range(joker));	// Matpack can copy the second row of
					// tz_intp to t_abs like this. But
					// t_abs has to have the right size!
  }

  //---------------< 2. Interpolation of VMR profiles >---------------
  {
    // The species lookup data
    extern const Array<SpeciesRecord> species_data;

    // check size of input String vectors.
    assert ( tgs.nelem() == raw_vmrs.nelem() ); 

    // Make vmrs the right size:
    vmrs.resize( raw_vmrs.nelem(), p_abs.nelem() );
    
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
	String tag_name = species_data[tgs[j][0].Species()].Name(); // name of the tag
	if ( (tag_name == "liquidcloud") || (tag_name == "icecloud") )
	  {
	    // Interpolate linearly the cloud profiles
	    interpp_cloud( vmrs(j,Range(joker)),
			   raw(Range(joker),0),
			   raw(Range(joker),1),
			   p_abs );
	    //	out3 << "This VMR: " << vmrs(j,Range(joker)) << "\n";
	  }
	else
	  {
	    // Interpolate VMRs:
	    interpp( vmrs(j,Range(joker)),
		     raw(Range(joker),0),
		     raw(Range(joker),1),
		     p_abs );
	    // out3 << "This VMR: " << vmrs(j,Range(joker)) << "\n";
	  }
	// The calls to interpp_cloud and inerpp contain some nice
	// Matpack  features:
	// 1. vmrs(j,Range(joker)) selects the jth row of vmrs.
	// 2. raw(Range(joker),0) and raw(Range(joker),1) select the
	// first and second column of raw. We don't need transpose
	// here, since the selected objects are vectors. 
	//
	// Note that you can call the interpolation functions directly
	// with these selections. 
      }
  }

  //---------< 3. Saturation adjustment VMR profiles of H2O tags in clouds >-----------
  {
    if( "yes" == CloudSatWV )
      {
	assert ( vmrs.ncols() == p_abs.nelem() );

	out2 << "Performing water vapor saturation adjustment for clouds.\n";

	Index satflag = WVsatinClouds( vmrs,  // manipulates this array for H2O tags
				       t_abs, // constant
				       p_abs, // constant
				       tgs);  // constant
	if (satflag != 0)
	  {
	    ostringstream os;
	    os << "WRONG WATER VAPOR SATURATION CALCULATION IN CLOUD REGION.";
	    throw runtime_error(os.str());
	  }
      }
    else if ( "no" == CloudSatWV )
      {
	out2 << "No water vapor saturation adjustment for clouds.\n";
      }
    else
      {
	out2 << "Normally it should be yes or no for CloudSatWV.\n"
	     << "But this will be moved out of this method anyway.\n";
// 	ostringstream os;
// 	os << "The keyword CloudSatWV must be yes or no.";
// 	throw runtime_error(os.str());
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
    const Index&       niter )
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

   \author Patrick Eriksson
   \date   2001-04-18
*/
void hseFromBottom(
          Vector&    hse,
    const Vector&    p_abs,
    const Vector&    z_abs,
    const Numeric&   g0,
    const Index&       niter )
{
  hseSet( hse, p_abs[0], z_abs[0], g0, niter );
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



// Algorithm based on equation from my (PE) thesis (page 274) and the book
// Meteorology today for scientists and engineers by R.B. Stull (pages 9-10). 
//
/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-19
*/
void hseCalc(
          Vector&    z_abs,
    const Vector&    p_abs,
    const Vector&    t_abs,
    const Vector&    h2o_abs,
    const Numeric&   r_geoid,   
    const Vector&    hse )
{
  if ( !isbool( static_cast<Index>(hse[0]) ) )  
    throw runtime_error(
        "The HSE flag (first element of hse) must either be 0 or 1.");
  
  if ( hse[0] )
  {
    if ( hse.nelem() != 5 )
    throw runtime_error(
        "The length of the hse vector must be 5.");

    const Index   np = p_abs.nelem();
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
  
    if ( (z_abs.nelem()!=np) || (t_abs.nelem()!=np) || (h2o_abs.nelem()!=np) )
      throw runtime_error(
                         "The input vectors do not all have the same length.");
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
	r  = 18/28.96 * (h2o_abs[i]+h2o_abs[i+1])/2;
  
	// The virtual temperature (no liquid water)
	tv = (1+0.61*r) * (t_abs[i]+t_abs[i+1])/2;
  
	// The change in vertical altitude from i to i+1 
	dz = 287.053*tv/g * log( p_abs[i]/p_abs[i+1] );
	ztmp[i+1] = ztmp[i] + dz;
      }
  
      // Match the altitude of the reference point
      dz = interpp( p_abs, ztmp, pref ) - zref;

      //  z_abs = ztmp - dz;
      z_abs = ztmp;
      z_abs -= dz;		// Note the new Matpack operations =
				// and -=
    }
  }
}


/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-18
*/
void h2o_absSet(
		Vector&          h2o_abs,
		const TagGroups& tgs,
		const Matrix&    vmrs )
{
  const Index   n = tgs.nelem();
  Index   found = -1;
  String  s;

  for( Index i=0; i<n && found<0; i++ ) 
  {
    s = tgs[i][0].Name();
    
    if ( s.substr(0,3) == "H2O" )
      found = i;
  }

  if ( found < 0 )
    throw runtime_error("h2o_absSet: No tag group contains water!");
  
  h2o_abs.resize( vmrs.ncols() );
  h2o_abs = vmrs(found,Range(joker));	
  // Matpack can copy the contents of vectors like this. The
  // dimensions must be the same! The expression
  // vmrs(found,Range(joker)) selects the row with index corresponding
  // to found.
}




/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Carlos Jimenez 
   \date   2001-08-14
*/
void vmrsScale(
	       Matrix&                vmrs,
	       const TagGroups&       tgs,
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
      vmrs(tagindex[itag],Range(joker)) *= scalfac[itag]; 
      // Matpack can multiply all elements of a vector with a constant
      // factor like this. In this case the vector is the selected row
      // of Matrix vmrs.
  
      //out2 << vmrs(tagindex[itag],Range(joker)) << ".\n";
    }
}







/**
   See the the online help (arts -d FUNCTION_NAME)
   Just a copy of the function 'h2o_absSet' 
   but now for nitrogen.

   \author Patrick Eriksson
   \date   2001-01-18
*/
void n2_absSet(
	       Vector&            n2_abs,
	       const   TagGroups& tgs,
	       const   Matrix&    vmrs )
{
  const Index   n = tgs.nelem();
  Index     found = -1;
  String  s;

  for( Index i=0; i<n && found<0; i++ ) 
  {
    s = tgs[i][0].Name();
    
    if ( s.substr(0,2) == "N2" )
      found = i;
  }

  if ( found < 0 )
    throw runtime_error("n2_absSet: No tag group contains nitrogen!");
  
  n2_abs.resize( vmrs.ncols() );
  n2_abs = vmrs(found,Range(joker));
  // Matpack can copy the contents of vectors like this. The
  // dimensions must be the same! The expression
  // vmrs(found,Range(joker)) selects the row with index corresponding
  // to found.
}


/** Calculates the absorption coefficients by first calculating the
   cross sections per tag group and then the absorption from the cross
   sections.

   \retval   abs            absorption coefficients
   \retval   abs_per_tg     absorption coefficients per tag group

   \param    tgs     the list of tag groups 
   \param    f_mono         monochromatic frequency grid
   \param    p_abs          pressure levels 
   \param    t_abs          temperature at pressure level
   \param    h2o_abs        total volume mixing ratio of water vapor
   \param    vmrs           volume mixing ratios per tag group
   \param    lines_per_tg   transition lines per tag group
   \param    lineshape      lineshape specifications to use per tag group
   \param    cont_description_names names of different continuum
                                    models
   \param    cont_description_parameters continuum parameters for the
                                         models listed in
					 cont_description_names 

   \author Axel von Engeln
   \date 2001-01-11

   \author Stefan Buehler
   \date 2001-03-13
 */
void absCalc(// WS Output:
             Matrix&        		     abs,
             ArrayOfMatrix& 		     abs_per_tg,
             // WS Input:		  
	     const TagGroups&                tgs,
             const Vector&  		     f_mono,
             const Vector&  		     p_abs,
             const Vector&  		     t_abs,
	     const Vector&  		     n2_abs,
	     const Vector&  		     h2o_abs,
             const Matrix&                   vmrs,
             const ArrayOfArrayOfLineRecord& lines_per_tg,
	     const ArrayOfLineshapeSpec&     lineshape,
	     const ArrayOfString&            cont_description_names,
	     const ArrayOfVector& 	     cont_description_parameters)
{
  // Dimension checks are performed in the executed functions

  // allocate local variable to hold the cross sections per tag group
  ArrayOfMatrix xsec_per_tg;

  xsec_per_tgInit( xsec_per_tg, tgs, f_mono, p_abs );

  xsec_per_tgAddLines( xsec_per_tg,
		       tgs,
		       f_mono,
		       p_abs,
		       t_abs,
		       h2o_abs,
		       vmrs,
		       lines_per_tg,
		       lineshape );

  xsec_per_tgAddConts( xsec_per_tg,
		       tgs,
		       f_mono,
		       p_abs,
		       t_abs,
		       n2_abs,
		       h2o_abs,
		       vmrs,
		       cont_description_names,
		       cont_description_parameters );

  absCalcFromXsec(abs,
		  abs_per_tg,
		  xsec_per_tg,
		  vmrs );

}


/**
   Calculates the absorption from a given cross section.

   Only the cross section and the vmrs are required, it is assumed
   that the vmrs are in the order of the cross sections, only the
   dimension is checked.

   \retval   abs            absorption coefficients
   \retval   abs_per_tg     absorption coefficients per tag group
   \param    xsec_per_tg    cross sections per tag group
   \param    vmrs           volume mixing ratios per tag group

   \author Stefan Bühler and Axel von Engeln
   \date   2001-01-11
*/
void absCalcFromXsec(// WS Output:
		     Matrix&        		     abs,
		     ArrayOfMatrix& 		     abs_per_tg,
		     // WS Input:		  
		     const ArrayOfMatrix&            xsec_per_tg,
		     const Matrix&                   vmrs)
{
  // Check that vmrs and xsec_per_tg really have compatible
  // dimensions. In vmrs there should be one row for each tg:
  if ( vmrs.nrows() != xsec_per_tg.nelem() )
    {
      ostringstream os;
      os << "Variable vmrs must have compatible dimension to xsec_per_tg.\n"
	 << "vmrs.nrows() = " << vmrs.nrows() << '\n'
	 << "xsec_per_tg.nelem() = " << xsec_per_tg.nelem();
      throw runtime_error(os.str());
    }

  // Check that number of altitudes are compatible. We only check the
  // first element, this is possilble because within arts all elements
  // are on the same altitude grid.
  if ( vmrs.ncols() != xsec_per_tg[0].ncols() )
    {
      ostringstream os;
      os << "Variable vmrs must have same numbers of altitudes as xsec_per_tg.\n"
	 << "vmrs.ncols() = " << vmrs.ncols() << '\n'
	 << "xsec_per_tg[0].ncols() = " << xsec_per_tg[0].ncols();
      throw runtime_error(os.str());
    }  
  
  // Initialize abs and abs_per_tg. The array dimension of abs_per_tg
  // is the same as that of xsec_per_tg. The dimension of abs should
  // be equal to one of the xsec_per_tg enries.
  abs.resize( xsec_per_tg[0].nrows(), xsec_per_tg[0].ncols() );
  abs = 0;			// Matpack can set all elements like this.

  abs_per_tg.resize( xsec_per_tg.nelem() );

  out2 << "  Computing abs and abs_per_tg from xsec_per_tg.\n";

  // Loop through all tag groups
  for ( Index i=0; i<xsec_per_tg.nelem(); ++i )
    {
      out2 << "  Tag group " << i << '\n';

      // Make this element of xsec_per_tg the right size:
      abs_per_tg[i].resize( xsec_per_tg[i].nrows(), xsec_per_tg[i].ncols() );
      abs_per_tg[i] = 0;	// Initialize all elements to 0.

      // Loop through all altitudes
      for ( Index j=0; j<xsec_per_tg[i].ncols(); j++)
	{
	  // Loop through all frequencies
	  for ( Index k=0; k<xsec_per_tg[i].nrows(); k++)
	    {
	      abs_per_tg[i](k,j) = xsec_per_tg[i](k,j) * vmrs(i,j);
	    }
	}

      // Add up to the total absorption:
      abs += abs_per_tg[i];	// In Matpack you can use the +=
				// operator to do elementwise addition.
    }
}


/**
   Initialize xsec_per_tg. The initialization is
   necessary, because methods `xsec_per_tgAddLines'
   and `xsec_per_tgAddConts' just add to xsec_per_tg.

   \retval   xsec_per_tg    cross section per tag group
   \param    tgs            the list of tag groups
   \param    f_mono         monochromatic frequency grid
   \param    p_abs          pressure levels 

   \author Stefan Buehler
   \date   2001-03-12
*/
void xsec_per_tgInit(// WS Output:
		     ArrayOfMatrix&   xsec_per_tg,
		     // WS Input:
		     const TagGroups& tgs,
		     const Vector&    f_mono,
		     const Vector&    p_abs
		     )
{
  // Initialize xsec_per_tg. The array dimension of xsec_per_tg
  // is the same as that of lines_per_tg.
  xsec_per_tg.resize( tgs.nelem() );

  // Loop xsec_per_tg and make each matrix the right size,
  // initializing to zero:
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
      // Make this element of abs_per_tg the right size:
      xsec_per_tg[i].resize( f_mono.nelem(), p_abs.nelem() );
      xsec_per_tg[i] = 0;	// Matpack can set all elements like this.
    }

  out2 << "  Initialized xsec_per_tg.\n"
       << "  Number of frequencies        : " << f_mono.nelem() << "\n"
       << "  Number of pressure levels    : " << p_abs.nelem() << "\n";
}

/**
   Calculates the line spectrum for each tag group and adds it to
   xsec_per_tg. 

   \retval   xsec_per_tg    cross section per tag group
   \param    tgs            the list of tag groups
   \param    f_mono         monochromatic frequency grid
   \param    p_abs          pressure levels 
   \param    t_abs          temperature at pressure level
   \param    h2o_abs        total volume mixing ratio of water vapor
   \param    vmrs           volume mixing ratios per tag group
   \param    lines_per_tg   transition lines per tag group
   \param    lineshape      lineshape specifications to use per tag group

   \author Stefan Bühler and Axel von Engeln
   \date   2001-01-11
*/
void xsec_per_tgAddLines(// WS Output:
			 ArrayOfMatrix& 		  xsec_per_tg,
			 // WS Input:		  
			 const TagGroups&                 tgs,
			 const Vector&  		  f_mono,
			 const Vector&  		  p_abs,
			 const Vector&  		  t_abs,
			 const Vector&  		  h2o_abs,
			 const Matrix&                    vmrs,
			 const ArrayOfArrayOfLineRecord&  lines_per_tg,
			 const ArrayOfLineshapeSpec&      lineshape)
{
  // Check that all paramters that should have the number of tag
  // groups as a dimension are consistent.
  {
    const Index n_tgs    = tgs.nelem();
    const Index n_xsec   = xsec_per_tg.nelem();
    const Index n_vmrs   = vmrs.nrows();
    const Index n_lines  = lines_per_tg.nelem();
    const Index n_shapes = lineshape.nelem();

    if ( n_tgs != n_xsec  ||
	 n_tgs != n_vmrs  ||
	 n_tgs != n_lines ||
	 n_tgs != n_shapes   )
      {
	ostringstream os;
	os << "The following variables must all have the same dimension:\n"
	   << "tgs:          " << tgs.nelem() << '\n'
	   << "xsec_per_tg:  " << xsec_per_tg.nelem() << '\n'
	   << "vmrs:         " << vmrs.nrows() << '\n'
	   << "lines_per_tg: " << lines_per_tg.nelem() << '\n'
	   << "lineshape:    " << lineshape.nelem();
	throw runtime_error(os.str());
      }
  }  

  // Print information:
  {
    // The variables defined here (in particular the frequency
    // conversion) are just to make the output nice. They are not used
    // in subsequent calculations.
    out2 << "  Calculating line spectra.\n";
    out3 << "  Transitions to do: \n";
    Index nlines = 0;
    String funit;
    Numeric ffac;
    if ( f_mono[0] < 3e12 )
      {
	funit = "GHz"; ffac = 1e9;
      }
    else
      {
	extern const Numeric SPEED_OF_LIGHT;
	funit = "cm-1"; ffac = SPEED_OF_LIGHT*100;
      }
    for ( Index i=0; i<lines_per_tg.nelem(); ++i )
      {
	for ( Index l=0; l<lines_per_tg[i].nelem(); ++l )
	  {
	    out3 << "    " << lines_per_tg[i][l].Name() << " @ " 
		 << lines_per_tg[i][l].F()/ffac  << " " << funit << " ("
		 << lines_per_tg[i][l].I0() << ")\n"; 
	    nlines++;
	  }
      }
    out2 << "  Total number of transistions : " << nlines << "\n";
  }

  // Call xsec_species for each tag group.
  for ( Index i=0; i<tgs.nelem(); ++i )
    {
      out2 << "  Tag group " << i
	   << " (" << get_tag_group_name(tgs[i]) << "): ";
      
      // Get a pointer to the line list for the current species. This
      // is just so that we don't have to type lines_per_tg[i] over
      // and over again.
      const ArrayOfLineRecord& ll = lines_per_tg[i];

      // Also get a pointer to the lineshape specification:
      const LineshapeSpec& ls = lineshape[i];
      
      // Skip the call to xsec_per_tg if the line list is empty.
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

	  out2 << ll.nelem() << " transitions\n";
	  xsec_species( xsec_per_tg[i],
			f_mono,
			p_abs,
			t_abs,
			h2o_abs,
			vmrs(i,Range(joker)),
			ll,
			ls.Ind_ls(),
			ls.Ind_lsn(),
			ls.Cutoff());
	  // Note that we call xsec_species with a row of vmrs,
	  // selected by the above Matpack expression. This is
	  // possible, because xsec_species is using Views.
	}
      else
	{
	  out2 << ll.nelem() << " transitions, skipping\n";
	}
    }
}

/**
   Calculates the continuum for each tag group and adds it to
   xsec_per_tg. 

   \retval   xsec_per_tg    cross section per tag group
   \param    tgs            the list of tag groups
   \param    f_mono         monochromatic frequency grid
   \param    p_abs          pressure levels 
   \param    t_abs          temperature at pressure level
   \param    n2_abs         total volume mixing ratio of nitrogen
   \param    h2o_abs        total volume mixing ratio of water vapor
   \param    vmrs           volume mixing ratios per tag group
   \param    cont_description_names names of different continuum
                                    models
   \param    cont_description_parameters continuum parameters for the
                                         models listed in
					 cont_description_names 

   \throw runtime_error something went wrong

   \author Stefan Bühler
   \date   2001-03-12
*/
void xsec_per_tgAddConts(// WS Output:
			 ArrayOfMatrix& 		  xsec_per_tg,
			 // WS Input:		  
			 const TagGroups&                 tgs,
			 const Vector&  		  f_mono,
			 const Vector&  		  p_abs,
			 const Vector&  		  t_abs,
			 const Vector&  		  n2_abs,
			 const Vector&  		  h2o_abs,
			 const Matrix&                    vmrs,
                         const ArrayOfString&             cont_description_names,
                         const ArrayOfVector& 		  cont_description_parameters )
{
  // Check that all paramters that should have the number of tag
  // groups as a dimension are consistent.
  {
    const Index n_tgs    = tgs.nelem();
    const Index n_xsec   = xsec_per_tg.nelem();
    const Index n_vmrs   = vmrs.nrows();

    if ( n_tgs != n_xsec || n_tgs != n_vmrs )
      {
	ostringstream os;
	os << "The following variables must all have the same dimension:\n"
	   << "tgs:          " << tgs.nelem() << '\n'
	   << "xsec_per_tg:  " << xsec_per_tg.nelem() << '\n'
	   << "vmrs:         " << vmrs.nrows();
	throw runtime_error(os.str());
      }
  }

  // Check, that dimensions of cont_description_names and
  // cont_description_parameters are consistent...
  if ( cont_description_names.nelem() !=
       cont_description_parameters.nelem() )
    {
	ostringstream os;
	os << "The following variables must have the same dimension:\n"
	   << "cont_description_names:      " << cont_description_names.nelem() << '\n'
	   << "cont_description_parameters: " << cont_description_parameters.nelem();
	throw runtime_error(os.str());
    }

  // ...and that indeed the names match valid continuum models:
  for ( Index i=0; i<cont_description_names.nelem(); ++i )
    {
      check_continuum_model(cont_description_names[i]);
    }

  out2 << "  Calculating continuum spectra.\n";

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
	      // 	  specific data. We need the right entry in this
	      // 	  table. The index of this is obtained by calling member function
	      // 	  Species() on the tag. Thus we have:
	      // 	  species_data[tgs[i][s].Species()].
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
		  // cont_description_names.
		  const Index n =
		    find( cont_description_names.begin(),
			  cont_description_names.end(),
			  name ) - cont_description_names.begin();

		  // n==cont_description_names.nelem() indicates that
		  // the name was not found.
		  if ( n==cont_description_names.nelem() )
		    {
		      ostringstream os;
		      os << "Cannot find model " << name
			 << " in cont_description_names.";
		      throw runtime_error(os.str());		      
		    }

		  // Ok, the tag specifies a valid continuum model and
		  // we have continuum parameters.
		  
		  out2 << "  Adding " << name
		       << " to tag group " << i << ".\n";

		  // Add the continuum for this tag. The parameters in
		  // this call should be clear. The vmr is in
		  // vmrs(i,Range(joker)). The other vmr variable, `h2o_abs'
		  // contains the real H2O vmr, which is needed for
		  // the oxygen continuum.
		  xsec_continuum_tag( xsec_per_tg[i],
				      name,
				      cont_description_parameters[n],
				      f_mono,
				      p_abs,
				      t_abs,
				      n2_abs,
				      h2o_abs,
				      vmrs(i,Range(joker)) );
		  // Calling this function with a row of Matrix vmrs
		  // is possible because it uses Views.
		}
	    }
	}
    }

}

/** Reduces the size of abs_per_tg.  Only absorption coefficients for
    which weighting functions are calculated are kept in memory.
    
    \retval abs_per_tg absorption coefficients
    \param  tgs        all selected tag groups
    \param  wfs_tgs    the tag groups for which we want weighting functions.

    \author Axel von Engeln and Stefan Buehler */
void abs_per_tgReduce(// WS Output:
                      ArrayOfMatrix&         abs_per_tg,
                      // WS Input:
                      const TagGroups&       tgs,
                      const TagGroups&       wfs_tgs)
{

  // Make a safety check that the dimensions of tgs and
  // abs_per_tg are the same (could happen that we call this workspace
  // method twice by accident).
  if ( abs_per_tg.nelem()!=tgs.nelem() )
    throw(runtime_error("The variables abs_per_tg and tgs must\n"
			"have the same dimension."));

  // Erase could be very inefficient in this case, since elements
  // behind the erased one are copied in order to fill the
  // gap. Therefore, we will construct a new abs_per_tg, and finally
  // use it to replace the old one.
  ArrayOfMatrix abs_per_tg_out( wfs_tgs.nelem() );

  // Go through the weighting function tag groups:
  for ( Index i=0; i<wfs_tgs.nelem(); ++i )
    {
      // Index to the elements of wfs_tgs in tgs:
      Index n;
      get_tag_group_index_for_tag_group( n, tgs, wfs_tgs[i] );

      abs_per_tg_out[i].resize( abs_per_tg[n].nrows(), abs_per_tg[n].ncols() );
      abs_per_tg_out[i] = abs_per_tg[n]; // Matpack can copy the contents of
					 // matrices like this. The dimensions
					 // must be the same! 
    }  

  // Copy the generated matrices back to abs_per_tg
  abs_per_tg.resize( wfs_tgs.nelem() );
  abs_per_tg = abs_per_tg_out;	// FIXME: It should be checked whether
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
              Index&      refr,
              Index&      refr_lfac,
              String&   refr_model,
        const Index&      on,
        const String&   model,
        const Index&      lfac )
{
  if ( !isbool( on ) )  
    throw runtime_error("The on/off flag must either be 0 or 1.");
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
              Index&      refr,
              Index&      refr_lfac,
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
              const Vector&   p_abs,
              const Vector&   t_abs,
              const Vector&   h2o_abs,
              const Index&      refr,
              const String&   refr_model )
{
  if ( !isbool( refr ) )  
    throw runtime_error("The refraction flag must either be 0 or 1.");

  if ( refr == 0 )
    refr_index.resize( 0 );

  else
  {
    if ( refr_model == "Unity" )
    {
      cout << "DOING Unity \n";
      const Index n = p_abs.nelem();
      refr_index.resize( n );
      refr_index = 1.0;		// Matpack can set all elements like this.
    }
    
    else if ( refr_model == "Boudouris" )
      refr_indexBoudouris( refr_index, p_abs, t_abs, h2o_abs );

    else if ( refr_model == "BoudourisDryAir" )
      refr_indexBoudourisDryAir( refr_index, p_abs, t_abs );

    else
    {
      ostringstream os;
      os << "Unknown refraction model: " << refr_model;
      throw runtime_error(os.str());
    }
  }
}




//======================================================================
//             Methods related to continua
//======================================================================

/**
   Initializes the two continuum description WSVs,
   `cont_description_names' and `cont_description_parameters'.  

   This method does not really do anything, except setting the two
   variables to empty Arrays. It is just necessary because the method
   `cont_descriptionAppend' wants to append to the variables.

   Formally, the continuum description WSVs are required by the
   absorption calculation methods (e.g., `absCalc'). Therefore you
   always have to call `cont_descriptionInit'.

   Usage example: cont_descriptionInit{}

   \author Stefan Buehler
   \date 2001-03-12 */
void cont_descriptionInit(// WS Output:
                          ArrayOfString& names,
                          ArrayOfVector& parameters)
{
  names.resize(0);
  parameters.resize(0);
  out2 << "  Initialized cont_description_names and cont_description_parameters.\n";
}


/**
   Append a continuum description to `cont_description_names' and
   `cont_description_parameters'. 

   It is checked that the name given is indeed the name of an allowed
   continuum model. This is done by looking in the species_data
   lookup table.

   \author Stefan Buehler
   \date 2001-03-12 */
void cont_descriptionAppend(// WS Output:
		       ArrayOfString& cont_description_names,
		       ArrayOfVector& cont_description_parameters,
		       // Control Parameters:
		       const String& name,
		       const Vector& parameters)
{

  // First we have to check that name matches a continuum species tag.
  check_continuum_model(name);

  // Add name and parameters to the apropriate variables:
  cont_description_names.push_back(name);
  cont_description_parameters.push_back(parameters);
}
