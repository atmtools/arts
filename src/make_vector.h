/* Copyright (C) 2000-2007 Stefan Buehler <sbuehler@ltu.se>

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
  \file   make_vector.h
  \brief  The class MakeVector is a special kind of Vector that can be
  initialized explicitly from one or more arguments of type Numeric.

  Usage is as simple as:

  \code
  MakeVector x(33.2, 17.3, 3.8);
  \endcode

  which will generate a Vector with the desired arguments. You can use
  a MakeVector like a Vector.

  \author Stefan Buehler
  \date   2001-09-15
*/

#ifndef make_vector_h
#define make_vector_h

#include "matpackI.h"

class MakeVector : public Vector
{
public:
  MakeVector();
  MakeVector(
        Numeric a0
        );
  MakeVector(
        Numeric a0,
        Numeric a1
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17,
        Numeric a18
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17,
        Numeric a18,
        Numeric a19
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17,
        Numeric a18,
        Numeric a19,
        Numeric a20
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17,
        Numeric a18,
        Numeric a19,
        Numeric a20,
        Numeric a21
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17,
        Numeric a18,
        Numeric a19,
        Numeric a20,
        Numeric a21,
        Numeric a22
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17,
        Numeric a18,
        Numeric a19,
        Numeric a20,
        Numeric a21,
        Numeric a22,
        Numeric a23
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17,
        Numeric a18,
        Numeric a19,
        Numeric a20,
        Numeric a21,
        Numeric a22,
        Numeric a23,
        Numeric a24
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17,
        Numeric a18,
        Numeric a19,
        Numeric a20,
        Numeric a21,
        Numeric a22,
        Numeric a23,
        Numeric a24,
        Numeric a25
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17,
        Numeric a18,
        Numeric a19,
        Numeric a20,
        Numeric a21,
        Numeric a22,
        Numeric a23,
        Numeric a24,
        Numeric a25,
        Numeric a26
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17,
        Numeric a18,
        Numeric a19,
        Numeric a20,
        Numeric a21,
        Numeric a22,
        Numeric a23,
        Numeric a24,
        Numeric a25,
        Numeric a26,
        Numeric a27
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17,
        Numeric a18,
        Numeric a19,
        Numeric a20,
        Numeric a21,
        Numeric a22,
        Numeric a23,
        Numeric a24,
        Numeric a25,
        Numeric a26,
        Numeric a27,
        Numeric a28
        );
  MakeVector(
        Numeric a0,
        Numeric a1,
        Numeric a2,
        Numeric a3,
        Numeric a4,
        Numeric a5,
        Numeric a6,
        Numeric a7,
        Numeric a8,
        Numeric a9,
        Numeric a10,
        Numeric a11,
        Numeric a12,
        Numeric a13,
        Numeric a14,
        Numeric a15,
        Numeric a16,
        Numeric a17,
        Numeric a18,
        Numeric a19,
        Numeric a20,
        Numeric a21,
        Numeric a22,
        Numeric a23,
        Numeric a24,
        Numeric a25,
        Numeric a26,
        Numeric a27,
        Numeric a28,
        Numeric a29
        );
};


#endif  // make_vector_h
