// Stuff related to the calculation of absorption coefficients.

#include "arts.h"
#include "vecmat.h"
#include "messages.h"
#include "file.h"
#include "absorption.h"
#include "workspace.h"

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


void linesWriteToNamedFile(// WS Input:
                           const ARRAYofLineRecord& lines,
                           // Control Parameters:
                           const string& filename)
{
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

void raw_vmr_profilesReadFromScenario(// WS Output:
                                      ARRAYofMATRIX&   raw_vmr_profiles,
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

      cout << name << endl;





    }
}
