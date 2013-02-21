/* Copyright (C) 2001-2012 Stefan Buehler <sbuehler@ltu.se>

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

/*!
  \file   make_array.h
  \brief  Implements the class MakeArray, which is a derived class of
  Array, allowing explicit initialization.

  \author Stefan Buehler
  \date   2001-09-13
*/

#ifndef make_array_h
#define make_array_h

#include "array.h"

/**
   Explicit construction of Arrays.

   The only purpose of this class is to provide constructors with which
   Arrays can be initialized explicitly. Example:

   Array<Index> b = MakeArray<Index>(1,2,3);

   will create an Array of Index with elements 1, 2, and 3. It is not
   possible to have such constructors for the class Array itself, due to
   the clash with the constructor setting the size. (For Index Arrays it
   the constructor setting the size could be interpreted as an explicit
   constructor for an Array with one element.) 
   
   Just use this class instead of Array whenever you want explicit
   initialization. The method information lookup table (see file
   methods.cc) is for example built that way.
*/
template<class base>
class MakeArray : public Array<base>
{
public:
  MakeArray();
  MakeArray(Array<base>& narray);
  MakeArray(
        const base& a0
        );
  MakeArray(
        const base& a0,
        const base& a1
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24,
        const base& a25
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24,
        const base& a25,
        const base& a26
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24,
        const base& a25,
        const base& a26,
        const base& a27
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24,
        const base& a25,
        const base& a26,
        const base& a27,
        const base& a28
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24,
        const base& a25,
        const base& a26,
        const base& a27,
        const base& a28,
        const base& a29
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24,
        const base& a25,
        const base& a26,
        const base& a27,
        const base& a28,
        const base& a29,
        const base& a30
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24,
        const base& a25,
        const base& a26,
        const base& a27,
        const base& a28,
        const base& a29,
        const base& a30,
        const base& a31
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24,
        const base& a25,
        const base& a26,
        const base& a27,
        const base& a28,
        const base& a29,
        const base& a30,
        const base& a31,
        const base& a32
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24,
        const base& a25,
        const base& a26,
        const base& a27,
        const base& a28,
        const base& a29,
        const base& a30,
        const base& a31,
        const base& a32,
        const base& a33
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24,
        const base& a25,
        const base& a26,
        const base& a27,
        const base& a28,
        const base& a29,
        const base& a30,
        const base& a31,
        const base& a32,
        const base& a33,
        const base& a34
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24,
        const base& a25,
        const base& a26,
        const base& a27,
        const base& a28,
        const base& a29,
        const base& a30,
        const base& a31,
        const base& a32,
        const base& a33,
        const base& a34,
        const base& a35
        );
  MakeArray(
        const base& a0,
        const base& a1,
        const base& a2,
        const base& a3,
        const base& a4,
        const base& a5,
        const base& a6,
        const base& a7,
        const base& a8,
        const base& a9,
        const base& a10,
        const base& a11,
        const base& a12,
        const base& a13,
        const base& a14,
        const base& a15,
        const base& a16,
        const base& a17,
        const base& a18,
        const base& a19,
        const base& a20,
        const base& a21,
        const base& a22,
        const base& a23,
        const base& a24,
        const base& a25,
        const base& a26,
        const base& a27,
        const base& a28,
        const base& a29,
        const base& a30,
        const base& a31,
        const base& a32,
        const base& a33,
        const base& a34,
        const base& a35,
        const base& a36
        );
};


// Define the functions here, to avoid inlining:

template<class base>
MakeArray<base>::MakeArray() : Array<base>(0)
{
  // Just an empty array.
}
template<class base>
MakeArray<base>::MakeArray(Array<base>& narray) : Array<base> (narray)
{
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0
                     ) : Array<base>(1)
{
  vector<base>::operator[](0) = a0;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1
                     ) : Array<base>(2)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2
                     ) : Array<base>(3)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3
                     ) : Array<base>(4)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4
                     ) : Array<base>(5)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5
                     ) : Array<base>(6)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6
                     ) : Array<base>(7)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7
                     ) : Array<base>(8)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8
                     ) : Array<base>(9)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9
                     ) : Array<base>(10)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10
                     ) : Array<base>(11)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11
                     ) : Array<base>(12)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12
                     ) : Array<base>(13)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13
                     ) : Array<base>(14)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14
                     ) : Array<base>(15)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15
                     ) : Array<base>(16)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16
                     ) : Array<base>(17)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17
                     ) : Array<base>(18)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18
                     ) : Array<base>(19)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19
                     ) : Array<base>(20)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20
                     ) : Array<base>(21)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21
                     ) : Array<base>(22)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22
                     ) : Array<base>(23)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23
                     ) : Array<base>(24)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24
                     ) : Array<base>(25)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24,
                     const base& a25
                     ) : Array<base>(26)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
  vector<base>::operator[](25) = a25;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24,
                     const base& a25,
                     const base& a26
                     ) : Array<base>(27)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
  vector<base>::operator[](25) = a25;
  vector<base>::operator[](26) = a26;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24,
                     const base& a25,
                     const base& a26,
                     const base& a27
                     ) : Array<base>(28)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
  vector<base>::operator[](25) = a25;
  vector<base>::operator[](26) = a26;
  vector<base>::operator[](27) = a27;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24,
                     const base& a25,
                     const base& a26,
                     const base& a27,
                     const base& a28
                     ) : Array<base>(29)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
  vector<base>::operator[](25) = a25;
  vector<base>::operator[](26) = a26;
  vector<base>::operator[](27) = a27;
  vector<base>::operator[](28) = a28;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24,
                     const base& a25,
                     const base& a26,
                     const base& a27,
                     const base& a28,
                     const base& a29
                     ) : Array<base>(30)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
  vector<base>::operator[](25) = a25;
  vector<base>::operator[](26) = a26;
  vector<base>::operator[](27) = a27;
  vector<base>::operator[](28) = a28;
  vector<base>::operator[](29) = a29;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24,
                     const base& a25,
                     const base& a26,
                     const base& a27,
                     const base& a28,
                     const base& a29,
                     const base& a30
                     ) : Array<base>(31)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
  vector<base>::operator[](25) = a25;
  vector<base>::operator[](26) = a26;
  vector<base>::operator[](27) = a27;
  vector<base>::operator[](28) = a28;
  vector<base>::operator[](29) = a29;
  vector<base>::operator[](30) = a30;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24,
                     const base& a25,
                     const base& a26,
                     const base& a27,
                     const base& a28,
                     const base& a29,
                     const base& a30,
                     const base& a31
                     ) : Array<base>(32)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
  vector<base>::operator[](25) = a25;
  vector<base>::operator[](26) = a26;
  vector<base>::operator[](27) = a27;
  vector<base>::operator[](28) = a28;
  vector<base>::operator[](29) = a29;
  vector<base>::operator[](30) = a30;
  vector<base>::operator[](31) = a31;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24,
                     const base& a25,
                     const base& a26,
                     const base& a27,
                     const base& a28,
                     const base& a29,
                     const base& a30,
                     const base& a31,
                     const base& a32
                      ) : Array<base>(33)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
  vector<base>::operator[](25) = a25;
  vector<base>::operator[](26) = a26;
  vector<base>::operator[](27) = a27;
  vector<base>::operator[](28) = a28;
  vector<base>::operator[](29) = a29;
  vector<base>::operator[](30) = a30;
  vector<base>::operator[](31) = a31;
  vector<base>::operator[](32) = a32;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24,
                     const base& a25,
                     const base& a26,
                     const base& a27,
                     const base& a28,
                     const base& a29,
                     const base& a30,
                     const base& a31,
                     const base& a32,
                     const base& a33
                      ) : Array<base>(34)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
  vector<base>::operator[](25) = a25;
  vector<base>::operator[](26) = a26;
  vector<base>::operator[](27) = a27;
  vector<base>::operator[](28) = a28;
  vector<base>::operator[](29) = a29;
  vector<base>::operator[](30) = a30;
  vector<base>::operator[](31) = a31;
  vector<base>::operator[](32) = a32;
  vector<base>::operator[](33) = a33;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24,
                     const base& a25,
                     const base& a26,
                     const base& a27,
                     const base& a28,
                     const base& a29,
                     const base& a30,
                     const base& a31,
                     const base& a32,
                     const base& a33,
                     const base& a34
                      ) : Array<base>(35)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
  vector<base>::operator[](25) = a25;
  vector<base>::operator[](26) = a26;
  vector<base>::operator[](27) = a27;
  vector<base>::operator[](28) = a28;
  vector<base>::operator[](29) = a29;
  vector<base>::operator[](30) = a30;
  vector<base>::operator[](31) = a31;
  vector<base>::operator[](32) = a32;
  vector<base>::operator[](33) = a33;
  vector<base>::operator[](34) = a34;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24,
                     const base& a25,
                     const base& a26,
                     const base& a27,
                     const base& a28,
                     const base& a29,
                     const base& a30,
                     const base& a31,
                     const base& a32,
                     const base& a33,
                     const base& a34,
                     const base& a35
                      ) : Array<base>(36)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
  vector<base>::operator[](25) = a25;
  vector<base>::operator[](26) = a26;
  vector<base>::operator[](27) = a27;
  vector<base>::operator[](28) = a28;
  vector<base>::operator[](29) = a29;
  vector<base>::operator[](30) = a30;
  vector<base>::operator[](31) = a31;
  vector<base>::operator[](32) = a32;
  vector<base>::operator[](33) = a33;
  vector<base>::operator[](34) = a34;
  vector<base>::operator[](35) = a35;
}
template<class base>
MakeArray<base>::MakeArray(
                     const base& a0,
                     const base& a1,
                     const base& a2,
                     const base& a3,
                     const base& a4,
                     const base& a5,
                     const base& a6,
                     const base& a7,
                     const base& a8,
                     const base& a9,
                     const base& a10,
                     const base& a11,
                     const base& a12,
                     const base& a13,
                     const base& a14,
                     const base& a15,
                     const base& a16,
                     const base& a17,
                     const base& a18,
                     const base& a19,
                     const base& a20,
                     const base& a21,
                     const base& a22,
                     const base& a23,
                     const base& a24,
                     const base& a25,
                     const base& a26,
                     const base& a27,
                     const base& a28,
                     const base& a29,
                     const base& a30,
                     const base& a31,
                     const base& a32,
                     const base& a33,
                     const base& a34,
                     const base& a35,
                     const base& a36
                      ) : Array<base>(37)
{
  vector<base>::operator[](0) = a0;
  vector<base>::operator[](1) = a1;
  vector<base>::operator[](2) = a2;
  vector<base>::operator[](3) = a3;
  vector<base>::operator[](4) = a4;
  vector<base>::operator[](5) = a5;
  vector<base>::operator[](6) = a6;
  vector<base>::operator[](7) = a7;
  vector<base>::operator[](8) = a8;
  vector<base>::operator[](9) = a9;
  vector<base>::operator[](10) = a10;
  vector<base>::operator[](11) = a11;
  vector<base>::operator[](12) = a12;
  vector<base>::operator[](13) = a13;
  vector<base>::operator[](14) = a14;
  vector<base>::operator[](15) = a15;
  vector<base>::operator[](16) = a16;
  vector<base>::operator[](17) = a17;
  vector<base>::operator[](18) = a18;
  vector<base>::operator[](19) = a19;
  vector<base>::operator[](20) = a20;
  vector<base>::operator[](21) = a21;
  vector<base>::operator[](22) = a22;
  vector<base>::operator[](23) = a23;
  vector<base>::operator[](24) = a24;
  vector<base>::operator[](25) = a25;
  vector<base>::operator[](26) = a26;
  vector<base>::operator[](27) = a27;
  vector<base>::operator[](28) = a28;
  vector<base>::operator[](29) = a29;
  vector<base>::operator[](30) = a30;
  vector<base>::operator[](31) = a31;
  vector<base>::operator[](32) = a32;
  vector<base>::operator[](33) = a33;
  vector<base>::operator[](34) = a34;
  vector<base>::operator[](35) = a35;
  vector<base>::operator[](36) = a36;
}



#endif  // make_array_h
