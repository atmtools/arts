/* Copyright (C) 2000 Stefan Buehler   <sbuehler@uni-bremen.de>
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


#include "arts.h"
#include "vecmat.h"
#include "messages.h"
#include "file.h"
#include "absorption.h"
#include "wsv.h"
#include "md.h"
#include "math_funcs.h"
#include "make_array.h"
#include "atm_funcs.h"
#include "continua.h"

void linesReadFromHitran(// WS Output:
                         ARRAYofLineRecord& lines,
                          // Control Parameters:
                         const string& filename,
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
  out2 << "  Read " << lines.size() << " lines.\n";
}


void linesReadFromMytran2(// WS Output:
			  ARRAYofLineRecord& lines,
                          // Control Parameters:
			  const string& filename,
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
  out2 << "  Read " << lines.size() << " lines.\n";
}

void linesReadFromJpl(// WS Output:
		      ARRAYofLineRecord& lines,
		      // Control Parameters:
		      const string& filename,
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
  out2 << "  Read " << lines.size() << " lines.\n";
}


void linesReadFromArts(// WS Output:
		       ARRAYofLineRecord& lines,
		       // Control Parameters:
		       const string& filename,
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
    string v;
    is >> v;
    if ( v!=lr.Version() )
      {
	ostringstream os;
	
	// If what we read is the version string, it should have at elast 9 characters.
	if ( 9 <= v.size() )
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
  out2 << "  Read " << lines.size() << " lines.\n";
}

void linesElowToJoule(// WS Output:
		      ARRAYofLineRecord& lines )
{
  for ( size_t i=0; i<lines.size(); ++i )
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
				    ARRAYofARRAYofLineRecord& lines_per_tg,
				    // WS Input:
				    const TagGroups& tgs,
                                    // Control Parameters:
                                    const ARRAY<string>& filenames,
                                    const ARRAY<string>& formats,
                                    const VECTOR& fmin,
                                    const VECTOR& fmax)
{
  const size_t n_tg   = tgs.size();	// # tag groups
  const size_t n_cat = filenames.size();	// # tag Catalogues

  // Check that dimensions of the keyword parameters are consistent
  // (must all be the same). 

  if ( n_cat != formats.size() ||
       n_cat != fmin.size() ||
       n_cat != fmax.size() )
    {
      ostringstream os;
      os << "lines_per_tgReadFromCatalogues: All keyword\n"
	 << "parameters must get the same number of arguments.\n"
	 << "You specified:\n"
	 << "filenames: " << n_cat         << "\n"
	 << "formats:   " << formats.size() << "\n"
	 << "fmin:      " << fmin.size()    << "\n"
	 << "fmax:      " << fmax.size();
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
  // real catalogues, along with an data structure that tells us which
  // tag groups should use this catalogue.

  ARRAY< string > real_filenames ( 1, filenames[0]    );
  ARRAY< string > real_formats   ( 1, formats[0]      );
  VECTOR real_fmin               ( 1, fmin[0]         );
  VECTOR real_fmax               ( 1, fmax[0]         );

  ARRAY< ARRAY <size_t> > real_tgs( 1, make_array<size_t>(0) );

  // The last specified catalogue, to which we should assign all
  // remaining lines. Index of this one in real_ arrays.
  size_t last_cat = 0;		

  for ( size_t i=1; i<n_tg; ++i )
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
	  const size_t that_cat = find( real_filenames.begin(),
					real_filenames.end(),
					filenames[i] ) - real_filenames.begin();
	  if ( that_cat < real_filenames.size() )
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
	      real_tgs.push_back( make_array<size_t>(i) );

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

  size_t n_real_cat = real_filenames.size(); // # real catalogues to read

  // Some output to low priority stream:
  out3 << "  Catalogues to read and tgs for which these will be used:\n";
  for ( size_t i=0; i<n_real_cat; ++i )
    {
      out3 << "  " << real_filenames[i] << ":";
      for ( size_t s=0; s<real_tgs[i].size(); ++s )
	out3 << " " << real_tgs[i][s];
      out3 << "\n";
    }

  // Make lines_per_tg the right size:
  resize( lines_per_tg, tgs.size() );

  // Loop through the catalogues to read:
  for ( size_t i=0; i<n_real_cat; ++i )
    {
      ARRAYofLineRecord   lines;

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
      TagGroups  these_tgs(real_tgs[i].size());
      for ( size_t s=0; s<real_tgs[i].size(); ++s )
	{
	  these_tgs[s] = tgs[real_tgs[i][s]];
	}

      // Create these_lines_per_tg:
      ARRAYofARRAYofLineRecord these_lines_per_tg;
      lines_per_tgCreateFromLines( these_lines_per_tg, lines, these_tgs );

      // Put these lines in the right place in lines_per_tg:
      for ( size_t s=0; s<real_tgs[i].size(); ++s )
	{
	  lines_per_tg[real_tgs[i][s]] = these_lines_per_tg[s];
	}
    }
}


void lines_per_tgCreateFromLines(// WS Output:
                                  ARRAYofARRAYofLineRecord& lines_per_tg,
                                  // WS Input:
                                  const ARRAYofLineRecord&   lines,
                                  const TagGroups&           tgs)
{
  // The species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;

  // As a safety feature, we will watch out for the case that a
  // species is included in the calculation, but not all lines are
  // used. For this we need an array to flag the used species:
  
  // For some weird reason, ARRAYs of bool do not work, although all
  // other types seem to work fine. So in this single case, I'll use
  // the stl vector directly. The other place where this is done is in
  // the function executor in main.cc.
  // FIXME: Fix this when ARRAY<bool> works.
  std::vector<bool> species_used (species_data.size(),false);
      
  // Make lines_per_tg the right size:
  lines_per_tg = ARRAYofARRAYofLineRecord(tgs.size());

  // Unfortunately, MTL conatains a bug that leads to all elements of
  // the outer ARRAY of an ARRAY<ARRAY>> pointing to the same data
  // after creation. So we need to fix this explicitly:
  for ( size_t i=0; i<lines_per_tg.size(); ++i )
    lines_per_tg[i] = ARRAYofLineRecord();

  // Loop all lines in the input line list:
  for ( size_t i=0; i<lines.size(); ++i )
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
      size_t j;

      // Loop the tag groups:
      for ( j=0; j<tgs.size() && !found ; ++j ) 
	{
	  // A tag group can contain several tags:
	  for ( size_t k=0; k<tgs[j].size() && !found; ++k )
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
	      if ( this_tag.Isotope() != this_line.SpeciesData().Isotope().size() )
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
  for (size_t i=0; i<tgs.size(); ++i)
    {
	out3 << "  " << i << ":";

	for (size_t s=0; s<tgs[i].size(); ++s)
	  out3 << " " << tgs[i][s].Name();

	out3 << ": " << lines_per_tg[i].size() << " lines\n";
    }

}

/** Adds mirror lines at negative frequencies to the lines_per_tg.
    For each line at frequency +f in lines_per_tg a corresponding entry at
    frequency -f is added to lines_per_tg.

    \retval lines_per_tg The array of arrays of lines for each tag group.
    
    \author Axel von Engeln and Stefan Buehler */
void lines_per_tgAddMirrorLines(// WS Output:
                                ARRAYofARRAYofLineRecord& lines_per_tg)
{
  // We will simply append the mirror lines after the original
  // lines. This way we don't have to make a backup copy of the
  // original lines. 

  for ( size_t i=0; i<lines_per_tg.size(); ++i )
    {
      // Get a reference to the current list of lines to save typing:
      ARRAYofLineRecord& ll = lines_per_tg[i];
      
      // For increased efficiency, reserve the necessary space:
      ll.reserve(2*ll.size());

      // Loop through all lines of this tag group:
      {
	// It is important that we determine the size of ll *before*
	// we start the loop. After all, we are adding elements. And
	// we cerainly don't want to continue looping the newly added
	// elements, we want to loop only the original elements.
	size_t n=ll.size();
	for ( size_t j=0; j<n; ++j )
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
			 ARRAYofARRAYofLineRecord& lines_per_tg,
			 // WS Input:
			 const ARRAYofLineshapeSpec& lineshape,
			 const VECTOR& f_mono)
{

  // Make sure lines_per_tg and lineshape have the same dimension:
  if ( lines_per_tg.size() != lineshape.size() ) 
    {
      ostringstream os;
      os << "Dimension of lines_per_tg does\n"
	 << "not match that of lineshape.";
      throw runtime_error(os.str());
    }
  
  // Make sure that the frequency grid is properly sorted:
  for ( INDEX s=0; s<f_mono.size()-1; ++s )
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
  for ( INDEX i=0; i<lines_per_tg.size(); ++i )
    {
      // Get cutoff frequency of this tag group:
      Numeric cutoff = lineshape[i].Cutoff();

      // Check whether cutoff is defined:
      if ( cutoff != -1)
	{
	  // Get a reference to the current list of lines to save typing:
	  ARRAYofLineRecord& ll = lines_per_tg[i];

	  // Calculate the borders:
	  Numeric upp = f_mono[f_mono.size()-1] + cutoff;
	  Numeric low = f_mono[0] - cutoff;

	  // Cycle through all lines within this tag group. 
	  for ( ARRAYofLineRecord::iterator j=ll.begin(); j<ll.end(); ++j )
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
		      const ARRAYofLineRecord& lines,
		      // Control Parameters:
		      const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+".lines.aa";
    }

  ofstream os;

  out2 << "  Writing file: " << filename << '\n';
  open_output_file(os, filename);

  write_lines_to_stream(os,lines);
}


void lines_per_tgWriteAscii(// WS Input:
			      const ARRAYofARRAYofLineRecord& lines_per_tg,
			      // Control Parameters:
			      const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+".lines_per_tg.aa";
    }

  ofstream os;

  out2 << "  Writing file: " << filename << '\n';
  open_output_file(os, filename);

  os << lines_per_tg.size() << '\n';

  for ( size_t i=0; i<lines_per_tg.size(); ++i )
    {
      const ARRAYofLineRecord& lines = lines_per_tg[i];
      write_lines_to_stream( os, lines );
    }
}


void tgsDefine(// WS Output:
	       TagGroups& tgs,
	       // Control Parameters:
	       const ARRAYofstring& tags)
{
  tgs = TagGroups(tags.size());

  //  cout << "Tags: " << tags << "\n";

  // Each element of the array of strings tags defines one tag
  // group. Let's work through them one by one.
  for ( size_t i=0; i<tags.size(); ++i )
    {
      // There can be a comma separated list of tag definitions, so we
      // need to break the string apart at the commas.
      ARRAY<string> tag_def;

      bool go_on = true;
      string these_tags = tags[i];
      while (go_on)
	{
	  size_t n = these_tags.find(',');
	  if ( n >= these_tags.size() )
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

      // tag_def now holds the different tag strings for this group.
      //      cout << "tag_def =\n" << tag_def << endl;


      // Unfortunately, MTL conatains a bug that leads to all elements of
      // the outer ARRAY of an ARRAY<ARRAY>> pointing to the same data
      // after creation. So we need to fix this explicitly:
      tgs[i] = ARRAY<OneTag>();

      for ( size_t s=0; s<tag_def.size(); ++s )
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
  for ( size_t i=0; i<tgs.size(); ++i )
    {
      out3 << "\n  " << i << ":";
      for ( size_t s=0; s<tgs[i].size(); ++s )
	{
	  out3 << " " << tgs[i][s].Name();
	}
    }
  out3 << '\n';

//  cout << endl << endl << tgs << endl;
}



void lineshapeDefine(// WS Output:
		     ARRAYofLineshapeSpec&    lineshape,
		     // WS Input:
		     const TagGroups&         tgs,
		     const string&            shape,
		     const string&            normalizationfactor,
		     const Numeric&           cutoff)
{
  // Make lineshape and normalization factor data visible:
  extern const ARRAY<LineshapeRecord> lineshape_data;
  extern const ARRAY<LineshapeNormRecord> lineshape_norm_data;


  // generate the right number of elements
  size_t tag_sz = tgs.size();
  resize(lineshape,tag_sz);

  // Is this lineshape available?
  int found0=-1;
  for ( size_t i=0; i<lineshape_data.size() && (found0 == -1) ; ++i )
    {
      const string& str = lineshape_data[i].Name();
      if (str == shape) 
	{
	  out2 << "  Selected lineshape: " << str << "\n";
	  found0=i;
	}
    }

  // Is this normalization to the lineshape available?
  int found1=-1;
  for ( size_t i=0; i<lineshape_norm_data.size() && (found1 == -1); ++i )
    {
      const string& str = lineshape_norm_data[i].Name();
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
  for (size_t i=0; i<tag_sz; i++)
    {
      lineshape[i].SetInd_ls( found0 );
      lineshape[i].SetInd_lsn( found1 );
      lineshape[i].SetCutoff( cutoff );
    }
}

void lineshape_per_tgDefine(// WS Output:
			    ARRAYofLineshapeSpec& lineshape,
			    // WS Input:
			    const TagGroups&      tgs,
			    const ARRAYofstring&  shape,
			    const ARRAYofstring&  normalizationfactor,
			    const VECTOR&         cutoff )
{
  // Make lineshape and normalization factor data visible:
  extern const ARRAY<LineshapeRecord> lineshape_data;
  extern const ARRAY<LineshapeNormRecord> lineshape_norm_data;

  // check that the number of elements are equal
  size_t tg_sz = tgs.size();
  if ( (tg_sz != shape.size()) ||
       (tg_sz != normalizationfactor.size()) || 
       (tg_sz != cutoff.size()) )
    {
      ostringstream os;
      os << "lineshape_per_tgDefine: number of elements does\n"
	 << "not match the number of tag groups defined.";
      throw runtime_error(os.str());
    }
      

  // generate the right number of elements
  resize(lineshape,tg_sz);

  // Is this lineshape available?
  for (size_t k=0; k<tg_sz; ++k)
    {
      int found0=-1;
      for ( size_t i=0; i<lineshape_data.size() && (found0 == -1); ++i )
	{
	  const string& str = lineshape_data[i].Name();
	  if (str == shape[k]) 
	    {
	      out2 << "  Tag Group: [";
	      for (size_t s=0; s<tgs[k].size()-1; ++s)
		out2 << tgs[k][s].Name() << ", "; 
	      out2 << tgs[k][tgs[k].size()-1].Name() << "]\n";
	      out2 << "  Selected lineshape: " << str << "\n";
	      found0=i;
	    }
	}

      // Is this normalization to the lineshape available?
      int found1=-1;
      for ( size_t i=0; i<lineshape_norm_data.size() && (found1 == -1); ++i )
	{
	  const string& str = lineshape_norm_data[i].Name();
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



void raw_vmrs_1dReadFromScenario(// WS Output:
                                 ARRAYofMATRIX&   raw_vmrs_1d,
                                 // WS Input:
                                 const TagGroups&     tgs,
                                 // Control Parameters:
                                 const string&    basename)
{
  // The species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;

  // We need to read one profile for each tag group.
  for ( size_t i=0; i<tgs.size(); ++i )
    {
      // Determine the name.
      string name =
	basename + "." +
	species_data[tgs[i][0].Species()].Name() + ".am";
      
      // Add an element for this tag group to the vmr profiles:
      raw_vmrs_1d.push_back(MATRIX());
      
      // Read the VMR:
      // (We use the workspace method MatrixReadAscii for this.)
      MatrixReadAscii(raw_vmrs_1d[i],"",name);
      
      // state the source of profile.
      //cout << species_data[tgs[i][0].Species()].Name() << " profile read from file: " << name << "\n";
    }
}

void raw_vmrs_1dReadFromFiles(// WS Output:
			      ARRAYofMATRIX&   raw_vmrs_1d,
			      // WS Input:
			      const TagGroups& tgs,
			      // Control Parameters:
			      const ARRAYofstring&  seltags,
			      const ARRAYofstring&  filenames,
			      const string&         basename)
{
  // The species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;
  //int isoflag = 0;  // wilde card set for isotope -> all isotopes
  const string wildcard = "*";

  // check size of input string vectors.
  assert ( filenames.size() == seltags.size() ); 

  for ( size_t i=0; i<tgs.size(); ++i )
    {
      string tname = species_data[tgs[i][0].Species()].Name();

      // loop over tags.
      int flagtag = 1;

      // loop over seltags.
      for ( size_t s=0; s<seltags.size(); ++s )
	{
	  if (flagtag == 1)
	    {
	      // check if more than one isotope wild card "*".
	      //assert ( seltags[s].find("-*") == seltags[s].rfind("-*"));
	  
	      // species name of actual seltags.
	      string sname = seltags[s];
	      string iname =  wildcard;
	      if ((sname.find("-") > 0) && (sname.find("-") < sname.length()))
		{
		  sname.erase(sname.find("-"),sname.length()-1);
		  //cout << "sname = " << sname << "\n";
		  
		  // species isotope name of actual seltags.
		  string iname = seltags[s];
		  iname.erase(0,iname.find("-")+1);
		  //cout << "iname = " << iname << "\n";
		}
	      // specify if specific isotope is given or not.
	      //if (iname == wildcard) {
	      //  isoflag = 0;  // wilde card set for isotope -> all isotopes
	      //} else {
	      //  isoflag = 1;
	      //}
	      
	      // seltag equal tag ?
	      if (tname == sname) 
		{
		  flagtag = 0; // seltags and tag mached
		  // cout << "no iso tag="<< tgs[i][0].Isotope() << "\n";
		  // for ( size_t k=0; k<tgs[i][0].Isotope(); ++k ) 
		  //{
		  //  cout << "  iso  tag["<< i << "] = " << 
		  //    species_data[tgs[i][0].Species()].Isotope()[k].Name() << "\n";
		  //	}
		  
		  string fname = filenames[s];
		  // Add an element for this tag group to the vmr profiles:
		  raw_vmrs_1d.push_back(MATRIX());
		  
		  // Read the VMR:
		  // (We use the workspace method MatrixReadAscii for this.)
		  MatrixReadAscii(raw_vmrs_1d[i],"",fname);
		  
		  // state the source of profile.
		  //cout << "  " << tname << " profile read from file: " << fname << "\n";
		}
	    }
	}
      if (flagtag == 1)
	{
	  // read from base name if seltags is not find 
	  // Determine the name.
	  string aname = basename + "." + tname + ".am";
	  // Add an element for this tag group to the vmr profiles:
	  raw_vmrs_1d.push_back(MATRIX());
	  
	  // Read the VMR:
	  // (We use the workspace method MatrixReadAscii for this.)
	  MatrixReadAscii(raw_vmrs_1d[i],"",aname);
	  
	  // state the source of profile.
	  //cout << "  " << tname << " profile read from file: " << aname << "\n";
	}
    }
}

// void Atm2dFromRaw1D(// WS Output:
//                     ARRAYofVECTOR& 	 t_abs_2d,
//                     ARRAYofVECTOR& 	 z_abs_2d,
//                     ARRAYofMATRIX& 	 vmrs_2d,
//                     // WS Input:      
//                     const VECTOR&  	 p_abs,
//                     const MATRIX&  	 raw_ptz_1d,
//                     const ARRAYofMATRIX& raw_vmrs_1d)
// {
//   // FIXME: This function is terrible! Make this better using MTL functionality.
//   // Does this function work at all? Has it ever been used? I think this is garbage. 

//   // This function uses a lot of copying. Rather inefficient. The
//   // problem is that the raw matrices do not directly fit the
//   // interpolation routines. If this turns out to be too slow, it
//   // should be replace by a routine using element-wise
//   // interpolation. This can be also efficient, if the search for the
//   // right interpolation point is done efficiently (see function
//   // `hunt' in numerical recipies).

//   // Also, I'm sure the copying could be done more elegantly, but I
//   // just wasted a lot of time trying to do this with matrix / vector
//   // notation. 


//   //---------------< 1. Interpolation of temperature and altitude >---------------
//   {  
//     // Safety check: Make sure that raw_ptz_1d really is a [x,3] matrix:
//     if ( 3 != raw_ptz_1d.ncols() )
//       {
// 	ostringstream os;
// 	os << "The variable raw_ptz_1d does not have the right dimensions,\n"
// 	   << "ncols() should be 3, but is actually "<< raw_ptz_1d.ncols();
// 	throw runtime_error(os.str());
//       }

//     // Break up raw_ptz_1d in p_raw, tz_raw.
//     // The reason why we take tz_raw as a matrix is that the
//     // interpolation can then be done simultaneously, hence slightly
//     // more efficient.

//     // p_raw is column 1:
//     VECTOR p_raw;
//     col( p_raw, 1, raw_ptz_1d );

//     // tz_raw is column 2-3:
//     MATRIX tz_raw;
//     col( tz_raw, 2, 3, raw_ptz_1d );

//     // Now interpolate tz_raw to p_abs grid:
//     MATRIX tz_intp( p_abs.size(), tz_raw.ncols() );
//     interp_lin_col( tz_intp,
// 		    p_raw, tz_raw, p_abs );

//     // Extract t_abs_2d:
//     resize(t_abs_2d,0);
//     t_abs_2d.push_back(VECTOR());
//     col( t_abs_2d[0], 1, tz_intp );

//     // Extract z_abs_2d:
//     resize(z_abs_2d,0);
//     z_abs_2d.push_back(VECTOR());
//     col( z_abs_2d[0], 2, tz_intp );
//   }

//   //---------------< 2. Interpolation of VMR profiles >---------------
//   {
//     // We will write everything to the first array element of vmrs_2d
//     // (more array elements would only be used in a 2D calculation).

//     // Get room for our results:
//     resize(vmrs_2d,0);
//     vmrs_2d.push_back(MATRIX());
  
//     // Get a convenient reference:
//     MATRIX& intp = vmrs_2d[0];

//     // Set dimensions.
//     // The first dimension is the number of profiles (= the number of
//     // tag groups). The second dimension is the dimension of the new
//     // pressure grid.
//     resize( intp, raw_vmrs_1d.size() , p_abs.size() );
  
//     // For sure, we need to loop through all VMR profiles:
//     for ( size_t j=0; j<raw_vmrs_1d.size(); ++j )
//       {
// 	// Get a reference to the profile we are concerned with
// 	const MATRIX& raw = raw_vmrs_1d[j];

// 	// Raw should be a matrix with dimension [x,2], the first column
// 	// is the raw pressure grid, the second column the VMR values.
      
// 	// Safety check to ensure this:
// 	if ( 2 != raw.ncols() )
// 	  {
// 	    ostringstream os;
// 	    os << "The variable raw_vmrs_1d("
// 	       << j
// 	       << ") does not have the right dimensions,\n"
// 	       << "ncols() should be 2, but is actually "<< raw.ncols();
// 	    throw runtime_error(os.str());
// 	  }

// 	// Interpolate:
// 	interp_lin( intp[j],
// 		    raw( 0, raw.nrows(), 0, 1 ),
// 		    raw( 0, raw.nrows(), 1, 2 ),
// 		    p_abs );
//       }
//   }
// }

void AtmFromRaw1D(// WS Output:
		  VECTOR& 	 t_abs,
		  VECTOR& 	 z_abs,
		  ARRAYofVECTOR& vmrs,
		  // WS Input:      
		  const TagGroups&       tgs,
		  const VECTOR&  	 p_abs,
		  const MATRIX&  	 raw_ptz_1d,
		  const ARRAYofMATRIX&   raw_vmrs_1d)
{
  
  //---------------< 1. Interpolation of temperature and altitude >---------------
  {  
    // Safety check: Make sure that raw_ptz_1d really is a [x,3] matrix:
    if ( 3 != raw_ptz_1d.ncols() )
      {
	ostringstream os;
	os << "The variable raw_ptz_1d does not have the right dimensions,\n"
	   << "ncols() should be 3, but is actually "<< raw_ptz_1d.ncols();
	throw runtime_error(os.str());
      }


    // Contents of raw_ptz_1d: p_raw is column 1,  tz_raw is column 2-3.

    // Interpolate tz_raw to p_abs grid.
    // The reason why we take tz_raw as a matrix is that the
    // interpolation can then be done simultaneously, hence slightly
    // more efficient.

    // We need to make temporary copies of things. Unnefficient!
    // FIXME: Improve this when we have better MATRIX functionality.

    // For the interpolated profiles:
    MATRIX tz_intp( 2, p_abs.size() );

    // The original pressure grid:
    VECTOR p0( raw_ptz_1d.nrows() );
    copy(columns(raw_ptz_1d)[0],p0);

    // The matrix to interpolate:
    MATRIX tz( 2, raw_ptz_1d.nrows() );
    copy(trans(raw_ptz_1d.sub_matrix(0,
				     raw_ptz_1d.nrows(),
				     1,
				     raw_ptz_1d.ncols())),
	 tz);

    interpp( tz_intp,
	     p0,
	     tz,
	     p_abs );

    // Extract t_abs:
    //    col( t_abs, 1, tz_intp );
    resize( t_abs, tz_intp.ncols() );
    copy( tz_intp[0], t_abs );

    // Extract z_abs:
    //    col( z_abs, 2, tz_intp );
    resize( z_abs, tz_intp.ncols() );
    copy( tz_intp[1], z_abs );
  }

  //---------------< 2. Interpolation of VMR profiles >---------------
  {
    // check size of input string vectors.
    assert ( tgs.size() == raw_vmrs_1d.size() ); 

    // Make vmrs the right size:
    resize( vmrs, raw_vmrs_1d.size() );
    
    // For sure, we need to loop through all VMR profiles:
    for ( size_t j=0; j<raw_vmrs_1d.size(); ++j )
      {

	// Get a reference to the profile we are concerned with:
	const MATRIX& raw = raw_vmrs_1d[j];

	// Raw should be a matrix with dimension [x,2], the first column
	// is the raw pressure grid, the second column the VMR values.
      
	// Safety check to ensure this:
	if ( 2 != raw.ncols() )
	  {
	    ostringstream os;
	    os << "The variable raw_vmrs_1d("
	       << j
	       << ") does not have the right dimensions,\n"
	       << "ncols() should be 2, but is actually "<< raw.ncols();
	    throw runtime_error(os.str());
	  }

	// Split the matrix into pressure and vmr vector:
	VECTOR p_raw(raw.nrows()), vmr_raw(raw.nrows());
	copy(columns(raw)[0],p_raw);
	copy(columns(raw)[1],vmr_raw);

	// Make vmrs[j] the right size:
	resize( vmrs[j], p_abs.size() );
	
	// Interpolate:
	interpp( vmrs[j],
		 p_raw,
		 vmr_raw,
		 p_abs );
	//	out3 << "This VMR: " << vmrs[j] << "\n";
      }
  }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-18
*/
void hseSet(
          VECTOR&    hse,
    const Numeric&   pref,
    const Numeric&   zref,
    const Numeric&   g0,
    const int&       niter )
{
  resize( hse, 5 );
  
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
          VECTOR&    hse,
    const VECTOR&    p_abs,
    const VECTOR&    z_abs,
    const Numeric&   g0,
    const int&       niter )
{
  hseSet( hse, p_abs[0], z_abs[0], g0, niter );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-18
*/
void hseOff(
	  VECTOR&    hse )
{
  resize( hse, 1 );
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
          VECTOR&    z_abs,
    const VECTOR&    p_abs,
    const VECTOR&    t_abs,
    const VECTOR&    h2o_abs,
    const Numeric&   r_geoid,   
    const VECTOR&    hse )
{
  if ( !isbool( int(hse[0]) ) )  
    throw runtime_error(
        "The HSE flag (first element of hse) must either be 0 or 1.");
  
  if ( hse[0] )
  {
    if ( hse.size() != 5 )
    throw runtime_error(
        "The length of the hse vector must be 5.");

    const size_t   np = p_abs.size();
	  size_t   i;                     // altitude index
	  Numeric  g;                     // gravitational acceleration
	  Numeric  r;                     // water mixing ratio in gram/gram
	  Numeric  tv;                    // virtual temperature
	  Numeric  dz;                    // step geometrical altitude
	  VECTOR   ztmp(np);              // temporary storage for z_abs
  
    // Pick out values from hse
    const Numeric   pref  = hse[1];
    const Numeric   zref  = hse[2];
    const Numeric   g0    = hse[3];
    const INDEX     niter = INDEX( hse[4] );
  
    if ( (z_abs.size()!=np) || (t_abs.size()!=np) || (h2o_abs.size()!=np) )
      throw runtime_error(
                         "The input vectors do not all have the same length.");
    if ( niter < 1 )
      throw runtime_error("The number of iterations must be > 0.");
  
    for ( INDEX iter=0; iter<niter; iter++ )
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
      setto(z_abs,-dz);
      add(ztmp,z_abs);		//  z_abs = ztmp - dz;
    }
  }
}



/**
   This method sets the H2o VMR profile to the saturated H2O VMR profile
   according to the saturation pressure calculated with WVSatPressureLiquidWater 
   for vapor-water and WVSatPressureIce for vapor-ice phase mixtures.

   \author Thomas Kuhn
   \date   2001-06-22
*/
void h2o_SatVMR( const TagGroups&       tgs,
                 const ARRAYofVECTOR&   vmrs )
{
  const INDEX   n = tgs.size();
  int   found = -1;
  string  s;
  VECTOR  h2o_sat;
  for( INDEX i=0; i<n && found<0; i++ ) 
  {
    s = tgs[i][0].Name();
    
    if ( s.substr(0,3) == "H2O" )
      found = int(i);
  }

  if ( found < 0 )
    throw runtime_error("h2o_SatVMR: No tag group contains water vapor!");
  
  // set the appropriate size of the h2o_sat vector
  resize( h2o_sat, vmrs[found].size() );

  //


  // copy the saturation water vapor VMR to the VMR matrix
  copy( h2o_sat, vmrs[found] );
}




/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-18
*/
void h2o_absSet(
              VECTOR&          h2o_abs,
        const TagGroups&       tgs,
        const ARRAYofVECTOR&   vmrs )
{
  const INDEX   n = tgs.size();
  int   found = -1;
  string  s;

  for( INDEX i=0; i<n && found<0; i++ ) 
  {
    s = tgs[i][0].Name();
    
    if ( s.substr(0,3) == "H2O" )
      found = int(i);
  }

  if ( found < 0 )
    throw runtime_error("h2o_absSet: No tag group contains water!");
  
  resize( h2o_abs, vmrs[found].size() );
  copy( vmrs[found], h2o_abs );
}




/**
   See the the online help (arts -d FUNCTION_NAME)
   Just a copy of the function 'h2o_absSet' 
   but now for nitrogen.

   \author Patrick Eriksson
   \date   2001-01-18
*/
void n2_absSet(
	       VECTOR&          n2_abs,
         const TagGroups&       tgs,
         const ARRAYofVECTOR&   vmrs )
{
  const INDEX   n = tgs.size();
  int     found = -1;
  string  s;

  for( INDEX i=0; i<n && found<0; i++ ) 
  {
    s = tgs[i][0].Name();
    
    if ( s.substr(0,2) == "N2" )
      found = int(i);
  }

  if ( found < 0 )
    throw runtime_error("n2_absSet: No tag group contains nitrogen!");
  
  resize( n2_abs, vmrs[found].size() );
  copy( vmrs[found], n2_abs );
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
             MATRIX&        		     abs,
             ARRAYofMATRIX& 		     abs_per_tg,
             // WS Input:		  
	     const TagGroups&                tgs,
             const VECTOR&  		     f_mono,
             const VECTOR&  		     p_abs,
             const VECTOR&  		     t_abs,
	     const VECTOR&  		     n2_abs,
	     const VECTOR&  		     h2o_abs,
             const ARRAYofVECTOR&            vmrs,
             const ARRAYofARRAYofLineRecord& lines_per_tg,
	     const ARRAYofLineshapeSpec&     lineshape,
	     const ARRAYofstring&            cont_description_names,
	     const ARRAYofVECTOR& 	     cont_description_parameters)
{
  // Dimension checks are performed in the executed functions

  // allocate local variable to hold the cross sections per tag group
  ARRAYofMATRIX xsec_per_tg;

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
		     MATRIX&        		     abs,
		     ARRAYofMATRIX& 		     abs_per_tg,
		     // WS Input:		  
		     const ARRAYofMATRIX&            xsec_per_tg,
		     const ARRAYofVECTOR&            vmrs)
{
  // Check that vmrs and xsec_per_tg really have compatible
  // dimensions. In vmrs there should be one VECTOR for each tg:
  if ( vmrs.size() != xsec_per_tg.size() )
    {
      ostringstream os;
      os << "Variable vmrs must have compatible dimension to xsec_per_tg.\n"
	 << "vmrs.size() = " << vmrs.size() << '\n'
	 << "xsec_per_tg.size() = " << xsec_per_tg.size();
      throw runtime_error(os.str());
    }

  // Check that number of altitudes are compatible. We only check the
  // first element, this is possilble because within arts all elements
  // are on the same altitude grid.
  if ( vmrs[0].size() != xsec_per_tg[0].ncols() )
    {
      ostringstream os;
      os << "Variable vmrs must have same numbers of altitudes as xsec_per_tg.\n"
	 << "vmrs[0].size() = " << vmrs[0].size() << '\n'
	 << "xsec_per_tg[0].ncols() = " << xsec_per_tg[0].ncols();
      throw runtime_error(os.str());
    }  
  
  // Initialize abs and abs_per_tg. The array dimension of abs_per_tg
  // is the same as that of xsec_per_tg. The dimension of abs should
  // be equal to one of the xsec_per_tg enries.
  resize( abs, xsec_per_tg[0].nrows(), xsec_per_tg[0].ncols() );
  setto( abs, 0);
  resize( abs_per_tg, xsec_per_tg.size() );

  out2 << "  Computing abs and abs_per_tg from xsec_per_tg.\n";

  // Loop through all tag groups
  for ( INDEX i=0; i<xsec_per_tg.size(); ++i )
    {
      out2 << "  Tag group " << i << '\n';

      // Make this element of xsec_per_tg the right size:
      resize( abs_per_tg[i], xsec_per_tg[i].nrows(), xsec_per_tg[i].ncols() );
      setto( abs_per_tg[i], 0 );

      // Loop through all altitudes
      for ( INDEX j=0; j<xsec_per_tg[i].ncols(); j++)
	{

	  // Loop through all frequencies
	  for ( INDEX k=0; k<xsec_per_tg[i].nrows(); k++)
	    {
	      abs_per_tg[i][k][j] = xsec_per_tg[i][k][j] * vmrs[i][j];
	    }
	}

	  // Add up to the total absorption:
	  add( abs_per_tg[i], abs );
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
		     ARRAYofMATRIX&   xsec_per_tg,
		     // WS Input:
		     const TagGroups& tgs,
		     const VECTOR&    f_mono,
		     const VECTOR&    p_abs
		     )
{
  // Initialize xsec_per_tg. The array dimension of xsec_per_tg
  // is the same as that of lines_per_tg.
  resize( xsec_per_tg, tgs.size() );

  // Loop xsec_per_tg and make each matrix the right size,
  // initializing to zero:
  for ( size_t i=0; i<tgs.size(); ++i )
    {
      // Make this element of abs_per_tg the right size:
      resize( xsec_per_tg[i], f_mono.size(), p_abs.size() );
      setto(  xsec_per_tg[i], 0 );
    }

  out2 << "  Initialized xsec_per_tg.\n"
       << "  Number of frequencies        : " << f_mono.size() << "\n"
       << "  Number of pressure levels    : " << p_abs.size() << "\n";
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
			 ARRAYofMATRIX& 		  xsec_per_tg,
			 // WS Input:		  
			 const TagGroups&                 tgs,
			 const VECTOR&  		  f_mono,
			 const VECTOR&  		  p_abs,
			 const VECTOR&  		  t_abs,
			 const VECTOR&  		  h2o_abs,
			 const ARRAYofVECTOR&             vmrs,
			 const ARRAYofARRAYofLineRecord&  lines_per_tg,
			 const ARRAYofLineshapeSpec&      lineshape)
{
  // Check that all paramters that should have the number of tag
  // groups as a dimension are consistent.
  {
    const size_t n_tgs    = tgs.size();
    const size_t n_xsec   = xsec_per_tg.size();
    const size_t n_vmrs   = vmrs.size();
    const size_t n_lines  = lines_per_tg.size();
    const size_t n_shapes = lineshape.size();

    if ( n_tgs != n_xsec  ||
	 n_tgs != n_vmrs  ||
	 n_tgs != n_lines ||
	 n_tgs != n_shapes   )
      {
	ostringstream os;
	os << "The following variables must all have the same dimension:\n"
	   << "tgs:          " << tgs.size() << '\n'
	   << "xsec_per_tg:  " << xsec_per_tg.size() << '\n'
	   << "vmrs:         " << vmrs.size() << '\n'
	   << "lines_per_tg: " << lines_per_tg.size() << '\n'
	   << "lineshape:    " << lineshape.size();
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
    size_t nlines = 0;
    string funit;
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
    for ( size_t i=0; i<lines_per_tg.size(); ++i )
      {
	for ( size_t l=0; l<lines_per_tg[i].size(); ++l )
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
  for ( size_t i=0; i<tgs.size(); ++i )
    {
      out2 << "  Tag group " << i
	   << " (" << get_tag_group_name(tgs[i]) << "): ";
      
      // Get a pointer to the line list for the current species. This
      // is just so that we don't have to type lines_per_tg[i] over
      // and over again.
      const ARRAYofLineRecord& ll = lines_per_tg[i];

      // Also get a pointer to the lineshape specification:
      const LineshapeSpec& ls = lineshape[i];
      
      // Skip the call to xsec_per_tg if the line list is empty.
      if ( 0 < ll.size() )
	{

	  // Get the name of the species. The member function name of a
	  // LineRecord returns the full name (species + isotope). So
	  // it is for us more convenient to get the species index
	  // from the LineRecord member function Species(), and then
	  // use this to look up the species name in species_data.
	  extern const ARRAY<SpeciesRecord> species_data;
	  string species_name = species_data[ll[0].Species()].Name();

	  // Get the name of the lineshape. For that we use the member
	  // function Ind_ls() to the lineshape data ls, which returns
	  // an index. With that index we can go into lineshape_data
	  // to get the name.
	  // We will need this for safety checks later on.
	  extern const ARRAY<LineshapeRecord> lineshape_data;
	  string lineshape_name = lineshape_data[ls.Ind_ls()].Name();


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

	  out2 << ll.size() << " transitions\n";
	  xsec_species( xsec_per_tg[i],
			f_mono,
			p_abs,
			t_abs,
			h2o_abs,
			vmrs[i],
			ll,
			ls.Ind_ls(),
			ls.Ind_lsn(),
			ls.Cutoff());
	}
      else
	{
	  out2 << ll.size() << " transitions, skipping\n";
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
			 ARRAYofMATRIX& 		  xsec_per_tg,
			 // WS Input:		  
			 const TagGroups&                 tgs,
			 const VECTOR&  		  f_mono,
			 const VECTOR&  		  p_abs,
			 const VECTOR&  		  t_abs,
			 const VECTOR&  		  n2_abs,
			 const VECTOR&  		  h2o_abs,
			 const ARRAYofVECTOR&             vmrs,
                         const ARRAYofstring&             cont_description_names,
                         const ARRAYofVECTOR& 		  cont_description_parameters )
{
  // Check that all paramters that should have the number of tag
  // groups as a dimension are consistent.
  {
    const size_t n_tgs    = tgs.size();
    const size_t n_xsec   = xsec_per_tg.size();
    const size_t n_vmrs   = vmrs.size();

    if ( n_tgs != n_xsec || n_tgs != n_vmrs )
      {
	ostringstream os;
	os << "The following variables must all have the same dimension:\n"
	   << "tgs:          " << tgs.size() << '\n'
	   << "xsec_per_tg:  " << xsec_per_tg.size() << '\n'
	   << "vmrs:         " << vmrs.size();
	throw runtime_error(os.str());
      }
  }

  // Check, that dimensions of cont_description_names and
  // cont_description_parameters are consistent...
  if ( cont_description_names.size() !=
       cont_description_parameters.size() )
    {
	ostringstream os;
	os << "The following variables must have the same dimension:\n"
	   << "cont_description_names:      " << cont_description_names.size() << '\n'
	   << "cont_description_parameters: " << cont_description_parameters.size();
	throw runtime_error(os.str());
    }

  // ...and that indeed the names match valid continuum models:
  for ( size_t i=0; i<cont_description_names.size(); ++i )
    {
      check_continuum_model(cont_description_names[i]);
    }

  out2 << "  Calculating continuum spectra.\n";

  // Loop tag groups:
  for ( size_t i=0; i<tgs.size(); ++i )
    {
      extern const ARRAY<SpeciesRecord> species_data; 

      // Go through the tags in the current tag group to see if they
      // are continuum tags:  
      for ( size_t s=0; s<tgs[i].size(); ++s )
	{
	  // First of all, we have to make sure that this is not a
	  // tag that means `all isotopes', because this should not
	  // include continuum. For such tags, tag.Isotope() will
	  // return the number of isotopes (i.e., one more than the
	  // allowed index range).
	  if ( tgs[i][s].Isotope() <
	       species_data[tgs[i][s].Species()].Isotope().size() )
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
		  const string name =
		    species_data[tgs[i][s].Species()].Name() + "-"
		    + species_data[tgs[i][s].Species()].Isotope()[tgs[i][s].Isotope()].Name();
  
		  // Check that this corresponds to a valid continuum
		  // model:
		  check_continuum_model(name);

		  // Check, if we have parameters for this model. For
		  // this, the model name must be listed in
		  // cont_description_names.
		  const size_t n =
		    find( cont_description_names.begin(),
			  cont_description_names.end(),
			  name ) - cont_description_names.begin();

		  // n==cont_description_names.size() indicates that
		  // the name was not found.
		  if ( n==cont_description_names.size() )
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
		  // vmrs[i]. The other vmr variable, `h2o_abs'
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
				      vmrs[i] );
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
                      ARRAYofMATRIX&         abs_per_tg,
                      // WS Input:
                      const TagGroups&       tgs,
                      const TagGroups&       wfs_tgs)
{

  // Make a safety check that the dimensions of tgs and
  // abs_per_tg are the same (could happen that we call this workspace
  // method twice by accident).
  if ( abs_per_tg.size()!=tgs.size() )
    throw(runtime_error("The variables abs_per_tg and tgs must\n"
			"have the same dimension."));

  // Erase could be very inefficient in this case, since elements
  // behind the erased one are copied in order to fill the
  // gap. Therefore, we will construct a new abs_per_tg, and finally
  // use it to replace the old one.
  ARRAYofMATRIX abs_per_tg_out( wfs_tgs.size() );

  // Go through the weighting function tag groups:
  for ( size_t i=0; i<wfs_tgs.size(); ++i )
    {
      // Index to the elements of wfs_tgs in tgs:
      size_t n;
      get_tag_group_index_for_tag_group( n, tgs, wfs_tgs[i] );

      resize( abs_per_tg_out[i], abs_per_tg[n].nrows(), abs_per_tg[n].ncols() );
      copy( abs_per_tg[n], abs_per_tg_out[i] );
    }  

  // Copy the generated matrices back to abs_per_tg
  //  resize( abs_per_tg, wfs_tgs.size() );
  //  copy( abs_per_tg_out, abs_per_tg );
  // FIXME: Replace by our own swap function? Or better make std::swap
  // work also for our own Vector and Matrix.
  std::swap( abs_per_tg_out, abs_per_tg );
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
              int&      refr,
              int&      refr_lfac,
              string&   refr_model,
        const int&      on,
        const string&   model,
        const int&      lfac )
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
              int&      refr,
              int&      refr_lfac,
              string&   refr_model )
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
                    VECTOR&   refr_index,
              const VECTOR&   p_abs,
              const VECTOR&   t_abs,
              const VECTOR&   h2o_abs,
              const int&      refr,
              const string&   refr_model )
{
  if ( !isbool( refr ) )  
    throw runtime_error("The refraction flag must either be 0 or 1.");

  if ( refr == 0 )
    resize( refr_index, 0 );

  else
  {
    if ( refr_model == "Unity" )
    {
      cout << "DOING Unity \n";
      const INDEX n = p_abs.size();
      resize( refr_index, n );
      setto( refr_index, 1.0 );
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
   variables to empty ARRAYs. It is just necessary because the method
   `cont_descriptionAppend' wants to append to the variables.

   Formally, the continuum description WSVs are required by the
   absorption calculation methods (e.g., `absCalc'). Therefore you
   always have to call `cont_descriptionInit'.

   Usage example: cont_descriptionInit{}

   \author Stefan Buehler
   \date 2001-03-12 */
void cont_descriptionInit(// WS Output:
                          ARRAYofstring& names,
                          ARRAYofVECTOR& parameters)
{
  resize(names,0);
  resize(parameters,0);
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
		       ARRAYofstring& cont_description_names,
		       ARRAYofVECTOR& cont_description_parameters,
		       // Control Parameters:
		       const string& name,
		       const VECTOR& parameters)
{

  // First we have to check that name matches a continuum species tag.
  check_continuum_model(name);

  // Add name and parameters to the apropriate variables:
  cont_description_names.push_back(name);
  cont_description_parameters.push_back(parameters);
}
