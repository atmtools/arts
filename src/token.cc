/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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

#include "arts.h"
#include "messages.h"
#include "token.h"

/** The name of the type associated with the different tokens. This
    has to be the name exactly as it appears in declarations of these
    variables in the program, because it is used by make_md_h.cc to
    automatically generate declarations for method functions. */
String TokValTypeName[7] = {"String", "int", "Numeric",
			    "Array<String>", "Array<int>", "Vector",
                            "undefined"};


// Conversion functions to read TokVal for the 6 different types: 
  
TokVal::operator String() const {
  assert (mtype == String_t);
  return ms;
}

TokVal::operator int() const {
  assert (mtype == int_t);
  return mn;
}
  
TokVal::operator Numeric() const {
  assert (mtype == Numeric_t);
  return mx;
}


TokVal::operator Array<String>() const {
  assert (mtype == Array_String_t);
  return msv;
}

TokVal::operator Array<int>() const {
  assert (mtype == Array_int_t);
  return mnv;
}
  
TokVal::operator Vector() const {
  assert (mtype == Vector_t);
  return mxv;
}


ostream& operator<<(ostream& os, const TokVal& a)
{
  switch (a.mtype)
    {
    case String_t:
      os << a.ms;
      break;
    case int_t:
      os << a.mn;
      break;
    case Numeric_t:
      os << a.mx;
      break;
    case Array_String_t:
      print_vector(os,a.msv);
      break;
    case Array_int_t:
      print_vector(os,a.mnv);
      break;
    case Vector_t:
      print_vector(os,a.mxv);
      break;
    default:
      out0 << "Undefined token type.\n";
      exit(1);
    }
  return os;
}


// main()
// {
//   String a("Test");
//   TokVal tv(a);

//   String b=tv;
//   cout << b << '\n';
//   Numeric c = 3.8;
//   TokVal tvtv(c);
//   tv = tvtv;
//   Numeric d = tv;
//   cout << d << '\n';
//   b = tv;			// should cause an error because of wrong type
// }
