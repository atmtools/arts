// Stuff related to the calculation of absorption coefficients.

#include "arts.h"
#include "vecmat.h"
#include "messages.h"

void AllAbsExample(// WS Output:
                   VECTOR& f_abs,
                   VECTOR& p_abs,
                   VECTOR& t_abs,
                   MATRIX& abs)
{
  out2 << "AllAbsExample!\n"
       << "Patrick, you can set these variables to reasonable\n"
       << "values here to have a test case.\n";
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

  out3 << "f_abs:\n" << f_abs << '\n';
  out3 << "p_abs:\n" << p_abs << '\n';
  out3 << "t_abs:\n" << t_abs << '\n';
  out3 << "abs:\n"   << abs << '\n';
}
