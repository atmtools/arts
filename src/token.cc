/* Copyright (C) 2000-2012 Stefan Buehler <sbuehler@ltu.se>

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

#include "token.h"
#include "arts.h"
#include "messages.h"

/** The name of the type associated with the different tokens. This
    has to be the name exactly as it appears in declarations of these
    variables in the program, because it is used by make_md_h.cc to
    automatically generate declarations for method functions. */
String TokValTypeName[8] = {"String",
                            "Index",
                            "Numeric",
                            "ArrayOfString",
                            "ArrayOfIndex",
                            "Vector",
                            "Matrix",
                            "undefined"};

// Conversion functions to read TokVal for the 6 different types:

TokVal::operator String() const {
  ARTS_ASSERT(mtype == String_t);
  return ms;
}

TokVal::operator Index() const {
  ARTS_ASSERT(mtype == Index_t);
  return mn;
}

TokVal::operator Numeric() const {
  ARTS_ASSERT(mtype == Numeric_t);
  return mx;
}

TokVal::operator ArrayOfString() const {
  ARTS_ASSERT(mtype == Array_String_t);
  return msv;
}

TokVal::operator ArrayOfIndex() const {
  ARTS_ASSERT(mtype == Array_Index_t);
  return mnv;
}

TokVal::operator ArrayOfSpeciesTag() const {
  ARTS_ASSERT(mtype == Array_SpeciesTa_t);
  return mnst;
}

TokVal::operator Vector() const {
  ARTS_ASSERT(mtype == Vector_t);
  return mxv;
}

TokVal::operator Matrix() const {
  ARTS_ASSERT(mtype == Matrix_t);
  return mm;
}

ostream& operator<<(ostream& os, const TokVal& a) {
  // This is needed for nice formating:
  bool first = true;

  switch (a.mtype) {
    case String_t:
      os << "\"" << a.ms << "\"";
      break;
    case Index_t:
      os << a.mn;
      break;
    case Numeric_t:
      os << a.mx;
      break;
    case Array_String_t:
      os << "[";
      for (Index i = 0; i < a.msv.nelem(); ++i) {
        if (first)
          first = false;
        else
          os << ",";
        os << "\"" << a.msv[i] << "\"";
      }
      os << "]";
      break;
    case Array_Index_t:
      os << "[";
      for (Index i = 0; i < a.mnv.nelem(); ++i) {
        if (first)
          first = false;
        else
          os << ",";

        os << a.mnv[i];
      }
      os << "]";
      break;
    case Vector_t:
      os << "[";
      for (Index i = 0; i < a.mxv.nelem(); ++i) {
        if (first)
          first = false;
        else
          os << ",";

        os << a.mxv[i];
      }
      os << "]";
      break;
    case Matrix_t:
      os << "[";
      for (Index i = 0; i < a.mm.nrows(); ++i) {
        for (Index j = 0; i < a.mm.ncols(); ++j) {
          if (first)
            first = false;
          else {
            if (j == 0)
              os << ";";
            else
              os << ",";
          }

          os << a.mm(i, j);
        }
      }
      os << "]";
      break;
    default:
      cerr << "Undefined token type.\n";
      arts_exit();
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
//   b = tv;                    // should cause an error because of wrong type
// }
