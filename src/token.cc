#include "arts.h"
#include "messages.h"
#include "token.h"

/** The name of the type associated with the different tokens. This
    has to be the name exactly as it appears in declarations of these
    variables in the program, because it is used by make_md_h.cc to
    automatically generate declarations for method functions. */
string TokValTypeName[6] = {"string", "int", "Numeric",
			    "ARRAY<String>", "ARRAY<Int>", "VECTOR"};


// Conversion functions to read TokVal for the 6 different types: 
  
TokVal::operator string() const {
  assert (mtype == str_);
  return ms;
}

TokVal::operator int() const {
  assert (mtype == int_);
  return mn;
}
  
TokVal::operator Numeric() const {
  assert (mtype == num_);
  return mx;
}


TokVal::operator ARRAY<string>() const {
  assert (mtype == strvec_);
  return msv;
}

TokVal::operator ARRAY<int>() const {
  assert (mtype == intvec_);
  return mnv;
}
  
TokVal::operator ARRAY<Numeric>() const {
  assert (mtype == numvec_);
  return mxv;
}


ostream& operator<<(ostream& os, const TokVal& a)
{
  switch (a.mtype)
    {
    case str_:
      os << a.ms;
      break;
    case int_:
      os << a.mn;
      break;
    case num_:
      os << a.mx;
      break;
    case strvec_:
      os << a.msv;
      break;
    case intvec_:
      os << a.mnv;
      break;
    case numvec_:
      os << a.mxv;
      break;
    default:
      out0 << "Undefined token type.\n";
      exit(1);
    }
  return os;
}


// main()
// {
//   string a("Test");
//   TokVal tv(a);

//   string b=tv;
//   cout << b << '\n';
//   Numeric c = 3.8;
//   TokVal tvtv(c);
//   tv = tvtv;
//   Numeric d = tv;
//   cout << d << '\n';
//   b = tv;			// should cause an error because of wrong type
// }
