// Stuff related to the calculation of absorption coefficients.

#include "arts.h"
#include "vecmat.h"
#include "messages.h"
#include "file.h"
#include "absorption.h"
#include "workspace.h"
#include "md.h"
#include "math_funcs.h"

void linesReadFromHitran(// WS Output:
                         ARRAYofLineRecord& lines,
                        // WS Generic Output Names:
                         const string& lines_name,
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
}


void linesWriteToFile(// WS Input:
		      const ARRAYofLineRecord& lines,
		      // Control Parameters:
		      const string& f)
{
  string filename = f;
  
  // If filename is empty then set the default filename:
  if ( "" == filename )
    {
      extern const string basename;                       
      filename = basename+".lines.al";
    }

  ofstream os;

  out2 << "  Writing file: " << filename << '\n';
  open_output_file(os, filename);

  for ( size_t i=0; i<lines.size(); ++i )
    {
      os << lines[i] << '\n';
    }
}


void tag_groupsDefine(// WS Output:
                      TagGroups& tag_groups,
                      // Control Parameters:
                      const ARRAY<string>& tags)
{
  tag_groups.resize(tags.size());

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
//    cout << "tag_def =\n" << tag_def << endl;

      for ( size_t s=0; s<tag_def.size(); ++s )
	{
	  // Remove leading whitespace, if there is any:
	  while ( ' '  == tag_def[s][0] ||
		  '\t' == tag_def[s][0]    )	tag_def[s].erase(0,1);

	  OneTag this_tag(tag_def[s]);

	  // Safety check: For s>0 check that the tags belong to the same species.
	  if (s>0)
	    if ( tag_groups[i][0].Species() != this_tag.Species() )
	      throw runtime_error("Tags in a tag group must belong to the same species.");

	  tag_groups[i].push_back(this_tag);
	}
    }

  // Print list of tag groups to the most verbose output stream:
  out3 << "  Defined tag groups:";
  for ( size_t i=0; i<tag_groups.size(); ++i )
    {
      out3 << "\n  " << i+1 << ":";
      for ( size_t s=0; s<tag_groups[i].size(); ++s )
	{
	  out3 << " " << tag_groups[i][s].Name();
	}
    }
  out3 << '\n';

//  cout << endl << endl << tag_groups << endl;
}

void raw_vmrs_1dReadFromScenario(// WS Output:
                                 ARRAYofMATRIX&   raw_vmrs_1d,
                                 // WS Input:
                                 const TagGroups& tag_groups,
                                 // Control Parameters:
                                 const string&    basename)
{
  // The species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;

  // We need to read one profile for each tag group.
  for ( size_t i=0; i<tag_groups.size(); ++i )
    {
      // Determine the name.
      string name =
	basename + "." +
	species_data[tag_groups[i][0].Species()].Name() + ".am";

      // Add an element for this tag group to the vmr profiles:
      raw_vmrs_1d.push_back(MATRIX());

      // Read the VMR:
      // (We use the workspace method MatrixReadFromFile for this.)
      MatrixReadFromFile(raw_vmrs_1d[i],"",name);
    }
}

void Atm2dFromRaw1D(// WS Output:
                    ARRAYofVECTOR& 	 t_abs_2d,
                    ARRAYofVECTOR& 	 z_abs_2d,
                    ARRAYofMATRIX& 	 vmrs_2d,
                    // WS Input:      
                    const VECTOR&  	 p_abs,
                    const MATRIX&  	 raw_ptz_1d,
                    const ARRAYofMATRIX& raw_vmrs_1d)
{
  // This function uses a lot of copying. Rather inefficient. The
  // problem is that the raw matrices do not directly fit the
  // interpolation routines. If this turns out to be too slow, it
  // should be replace by a routine using element-wise
  // interpolation. This can be also efficient, if the search for the
  // right interpolation point is done efficiently (see function
  // `hunt' in numerical recipies).

  // Also, I'm sure the copying could be done more elegantly, but I
  // just wasted a lot of time trying to do this with matrix / vector
  // notation. 


  //---------------< 1. Interpolation of temperature and altitude >---------------
  {  
    // Safety check: Make sure that raw_ptz_1d really is a [x,3] matrix:
    if ( 3 != raw_ptz_1d.dim(2) )
      {
	ostringstream os;
	os << "The variable raw_ptz_1d does not have the right dimensions,\n"
	   << "dim(2) should be 3, but is actually "<< raw_ptz_1d.dim(2);
	throw runtime_error(os.str());
      }

    // Break up raw_ptz_1d in p_raw, tz_raw.
    // The reason why we take tz_raw as a matrix is that the
    // interpolation can then be done simultaneously, hence slightly
    // more efficient.

    // p_raw is column 1:
    VECTOR p_raw;
    col( p_raw, 1, raw_ptz_1d );

    // tz_raw is column 2-3:
    MATRIX tz_raw;
    col( tz_raw, 2, 3, raw_ptz_1d );

    // Now interpolate tz_raw to p_abs grid:
    MATRIX tz_intp;
    interp_lin_col( tz_intp,
		    p_raw, tz_raw, p_abs );

    // Extract t_abs_2d:
    t_abs_2d.clear();
    t_abs_2d.push_back(VECTOR());
    col( t_abs_2d(1), 1, tz_intp );

    // Extract z_abs_2d:
    z_abs_2d.clear();
    z_abs_2d.push_back(VECTOR());
    col( z_abs_2d(1), 2, tz_intp );
  }

  //---------------< 2. Interpolation of VMR profiles >---------------
  {
    // We will write everything to the first array element of vmrs_2d
    // (more array elements would only be used in a 2D calculation).

    // Get room for our results:
    vmrs_2d.clear();
    vmrs_2d.push_back(MATRIX());
  
    // Get a convenient reference:
    MATRIX& intp = vmrs_2d(1);

    // Set dimensions.
    // The first dimension is the number of profiles (= the number of
    // tag groups). The second dimension is the dimension of the new
    // pressure grid.
    intp.newsize( raw_vmrs_1d.dim() , p_abs.dim() );
  
    // We need this for each profile, therefore we define it here:
    VECTOR target;

    // For sure, we need to loop through all VMR profiles:
    for ( size_t j=1; j<=raw_vmrs_1d.dim(); ++j )
      {
	// Get a reference to the profile we are concerned with
	const MATRIX& raw = raw_vmrs_1d(j);

	// Raw should be a matrix with dimension [x,2], the first column
	// is the raw pressure grid, the second column the VMR values.
      
	// Safety check to ensure this:
	if ( 2 != raw.dim(2) )
	  {
	    ostringstream os;
	    os << "The variable raw_vmrs_1d("
	       << j
	       << ") does not have the right dimensions,\n"
	       << "dim(2) should be 2, but is actually "<< raw.dim(2);
	    throw runtime_error(os.str());
	  }

	// Extract p_raw and vmr_raw:
	VECTOR p_raw;
	col( p_raw, 1, raw );

	VECTOR vmr_raw;
	col( vmr_raw, 2, raw );

	// Interpolate:
	interp_lin( target,
		    p_raw, vmr_raw, p_abs );

	// Put the result in the apropriate row of intp:
	for ( size_t i=1; i<=p_abs.dim(); ++i )
	  {
	    intp(j,i) = target(i);
	  }
      }
  }
}

void AllAbsExample(// WS Output:
                   VECTOR& f_abs,
                   VECTOR& p_abs,
                   VECTOR& t_abs,
                   MATRIX& abs)
{
  // Patrick, you can set these variables to reasonable
  // values here to have a test case.
  f_abs = VECTOR(3,"500e9 501e9 502e9");
  p_abs = VECTOR(4,"1000 100 10 1");
  t_abs = VECTOR(4,"300 250 260 290");
  abs   = MATRIX(f_abs.size(),p_abs.size(),
		 "1000 800 600 400 "
		 "1001 801 601 401 "
		 "1002 802 602 402 ");
  
  // Safety check:
  if ( p_abs.size() != t_abs.size() )
    {
      ostringstream os;
      os << "Dimensions of p_abs and t_abs must be the same!\n"
	 << "p_abs.size() = " << p_abs.size() << '\n'
	 << "t_abs.size() = " << t_abs.size() << '\n';
      throw runtime_error(os.str());
    }

//   out3 << "f_abs:\n" << f_abs << '\n';
//   out3 << "p_abs:\n" << p_abs << '\n';
//   out3 << "t_abs:\n" << t_abs << '\n';
//   out3 << "abs:\n"   << abs << '\n';
}


