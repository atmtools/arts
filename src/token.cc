#include "arts.h"
#include "messages.h"
#include "token.h"

string TokValTypeName[6] = {"string", "int", "Numeric",
			    "StringVector", "IntVector", "NumVector"};


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
