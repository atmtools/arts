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

#ifndef token_h
#define token_h

#include "array.h"
#include "abs_species_tags.h"
#include "matpackI.h"
#include "mystring.h"

/** The different token value types. These are the types that keyword
    parameters in the controlfile can have. */
enum TokValType {
  String_t,
  Index_t,
  Numeric_t,
  Array_String_t,
  Array_Index_t,
  Array_SpeciesTag_t,
  Vector_t,
  Matrix_t,
  undefined_t
};

/** This stores arbitrary token values and remembers the type. Only
    the correct type can be extracted again. */
class TokVal {
 public:
  /** Default Constructor. (Sets type to undefined_t) */
  TokVal()
      : mtype(undefined_t), ms(), mn(-1), mx(0.), msv(), mnv(), mnst(), mxv(), mm() {}

  /** To set TokVal from String (C - style). */
  TokVal(const char s[])
      : mtype(String_t), ms(s), mn(-1), mx(0.), msv(), mnv(), mnst(), mxv(), mm() {}

  /** To set TokVal from String (C++ - style). */
  TokVal(const String& s)
      : mtype(String_t), ms(s), mn(-1), mx(0.), msv(), mnv(), mnst(), mxv(), mm() {}

  /** To set TokVal from an integer. */
  TokVal(Index n)
      : mtype(Index_t), ms(), mn(n), mx(0.), msv(), mnv(), mnst(), mxv(), mm() {}

  /** To set TokVal from a Numeric. */
  TokVal(Numeric x)
      : mtype(Numeric_t), ms(), mn(-1), mx(x), msv(), mnv(), mnst(), mxv(), mm() {}

  /** To set TokVal from an array of Strings. */
  TokVal(ArrayOfString sv)
      : mtype(Array_String_t),
        ms(),
        mn(-1),
        mx(0.),
        msv(sv),
        mnv(),
        mnst(),
        mxv(),
        mm() {}

  /** To set TokVal from an array of integers. */
  TokVal(ArrayOfIndex nv)
      : mtype(Array_Index_t),
        ms(),
        mn(-1),
        mx(0.),
        msv(),
        mnv(nv),
        mnst(),
        mxv(),
        mm() {}

  /** To set TokVal from an array of species tags. */
  TokVal(ArrayOfSpeciesTag nst)
      : mtype(Array_SpeciesTag_t),
        ms(),
        mn(-1),
        mx(0.),
        msv(),
        mnv(),
        mnst(nst),
        mxv(),
        mm() {}

  /** To set TokVal from a Vector. */
  TokVal(Vector xv)
      : mtype(Vector_t), ms(), mn(-1), mx(0.), msv(), mnv(), mnst(), mxv(xv), mm() {}

  /** To set TokVal from a Matrix. */
  TokVal(Matrix m)
      : mtype(Matrix_t), ms(), mn(-1), mx(0.), msv(), mnv(), mnst(), mxv(), mm(m) {}

  // Conversion functions to return TokVal for the 6 different types:

  /** Return String. */
  operator String() const;
  /** Return Index. */
  operator Index() const;
  /** Return Numeric. */
  operator Numeric() const;

  /** Return array of Strings. */
  operator ArrayOfString() const;
  /** Return array of integers. */
  operator ArrayOfIndex() const;
  /** Return array of integers. */
  operator ArrayOfSpeciesTag() const;
  /** Return Vector. */
  operator Vector() const;
  /** Return Matrix. */
  operator Matrix() const;

  /** Output operator. */
  friend std::ostream& operator<<(std::ostream& os, const TokVal& a);

 private:
  TokValType mtype;
  String ms;
  Index mn;
  Numeric mx;
  ArrayOfString msv;
  ArrayOfIndex mnv;
  ArrayOfSpeciesTag mnst;
  Vector mxv;
  Matrix mm;
};

// typedef Array<TokValType> TokValTypeVector;
// typedef Array<TokVal>     TokValVector;

#endif  // token_h
