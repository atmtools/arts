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
      std::ostrstream os;
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
