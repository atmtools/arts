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

#include "make_vector.h"


MakeVector::MakeVector() : Vector(0)
{
  // Just an empty Vector.
}

MakeVector::MakeVector(
		     Numeric a0
		     ) : Vector(1)
{
  Vector::operator[](0) = a0;
}

MakeVector::MakeVector(
		     Numeric a0,
		     Numeric a1
		     ) : Vector(2)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
}

MakeVector::MakeVector(
		     Numeric a0,
		     Numeric a1,
		     Numeric a2
		     ) : Vector(3)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
}

MakeVector::MakeVector(
		     Numeric a0,
		     Numeric a1,
		     Numeric a2,
		     Numeric a3
		     ) : Vector(4)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
}

MakeVector::MakeVector(
		     Numeric a0,
		     Numeric a1,
		     Numeric a2,
		     Numeric a3,
		     Numeric a4
		     ) : Vector(5)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
}

MakeVector::MakeVector(
		     Numeric a0,
		     Numeric a1,
		     Numeric a2,
		     Numeric a3,
		     Numeric a4,
		     Numeric a5
		     ) : Vector(6)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
}

MakeVector::MakeVector(
		     Numeric a0,
		     Numeric a1,
		     Numeric a2,
		     Numeric a3,
		     Numeric a4,
		     Numeric a5,
		     Numeric a6
		     ) : Vector(7)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
}

MakeVector::MakeVector(
		     Numeric a0,
		     Numeric a1,
		     Numeric a2,
		     Numeric a3,
		     Numeric a4,
		     Numeric a5,
		     Numeric a6,
		     Numeric a7
		     ) : Vector(8)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
}

MakeVector::MakeVector(
		     Numeric a0,
		     Numeric a1,
		     Numeric a2,
		     Numeric a3,
		     Numeric a4,
		     Numeric a5,
		     Numeric a6,
		     Numeric a7,
		     Numeric a8
		     ) : Vector(9)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
}

MakeVector::MakeVector(
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
		     ) : Vector(10)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
}

MakeVector::MakeVector(
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
		     ) : Vector(11)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
}

MakeVector::MakeVector(
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
		     ) : Vector(12)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
}

MakeVector::MakeVector(
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
		     ) : Vector(13)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
}

MakeVector::MakeVector(
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
		     ) : Vector(14)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
}

MakeVector::MakeVector(
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
		     ) : Vector(15)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
}

MakeVector::MakeVector(
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
		     ) : Vector(16)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
}

MakeVector::MakeVector(
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
		     ) : Vector(17)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
}

MakeVector::MakeVector(
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
		     ) : Vector(18)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
}

MakeVector::MakeVector(
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
		     ) : Vector(19)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
  Vector::operator[](18) = a18;
}

MakeVector::MakeVector(
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
		     ) : Vector(20)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
  Vector::operator[](18) = a18;
  Vector::operator[](19) = a19;
}

MakeVector::MakeVector(
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
		     ) : Vector(21)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
  Vector::operator[](18) = a18;
  Vector::operator[](19) = a19;
  Vector::operator[](20) = a20;
}

MakeVector::MakeVector(
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
		     ) : Vector(22)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
  Vector::operator[](18) = a18;
  Vector::operator[](19) = a19;
  Vector::operator[](20) = a20;
  Vector::operator[](21) = a21;
}

MakeVector::MakeVector(
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
		     ) : Vector(23)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
  Vector::operator[](18) = a18;
  Vector::operator[](19) = a19;
  Vector::operator[](20) = a20;
  Vector::operator[](21) = a21;
  Vector::operator[](22) = a22;
}

MakeVector::MakeVector(
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
		     ) : Vector(24)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
  Vector::operator[](18) = a18;
  Vector::operator[](19) = a19;
  Vector::operator[](20) = a20;
  Vector::operator[](21) = a21;
  Vector::operator[](22) = a22;
  Vector::operator[](23) = a23;
}

MakeVector::MakeVector(
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
		     ) : Vector(25)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
  Vector::operator[](18) = a18;
  Vector::operator[](19) = a19;
  Vector::operator[](20) = a20;
  Vector::operator[](21) = a21;
  Vector::operator[](22) = a22;
  Vector::operator[](23) = a23;
  Vector::operator[](24) = a24;
}

MakeVector::MakeVector(
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
		     ) : Vector(26)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
  Vector::operator[](18) = a18;
  Vector::operator[](19) = a19;
  Vector::operator[](20) = a20;
  Vector::operator[](21) = a21;
  Vector::operator[](22) = a22;
  Vector::operator[](23) = a23;
  Vector::operator[](24) = a24;
  Vector::operator[](25) = a25;
}

MakeVector::MakeVector(
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
		     ) : Vector(27)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
  Vector::operator[](18) = a18;
  Vector::operator[](19) = a19;
  Vector::operator[](20) = a20;
  Vector::operator[](21) = a21;
  Vector::operator[](22) = a22;
  Vector::operator[](23) = a23;
  Vector::operator[](24) = a24;
  Vector::operator[](25) = a25;
  Vector::operator[](26) = a26;
}

MakeVector::MakeVector(
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
		     ) : Vector(28)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
  Vector::operator[](18) = a18;
  Vector::operator[](19) = a19;
  Vector::operator[](20) = a20;
  Vector::operator[](21) = a21;
  Vector::operator[](22) = a22;
  Vector::operator[](23) = a23;
  Vector::operator[](24) = a24;
  Vector::operator[](25) = a25;
  Vector::operator[](26) = a26;
  Vector::operator[](27) = a27;
}

MakeVector::MakeVector(
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
		     ) : Vector(29)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
  Vector::operator[](18) = a18;
  Vector::operator[](19) = a19;
  Vector::operator[](20) = a20;
  Vector::operator[](21) = a21;
  Vector::operator[](22) = a22;
  Vector::operator[](23) = a23;
  Vector::operator[](24) = a24;
  Vector::operator[](25) = a25;
  Vector::operator[](26) = a26;
  Vector::operator[](27) = a27;
  Vector::operator[](28) = a28;
}

MakeVector::MakeVector(
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
		     ) : Vector(30)
{
  Vector::operator[](0) = a0;
  Vector::operator[](1) = a1;
  Vector::operator[](2) = a2;
  Vector::operator[](3) = a3;
  Vector::operator[](4) = a4;
  Vector::operator[](5) = a5;
  Vector::operator[](6) = a6;
  Vector::operator[](7) = a7;
  Vector::operator[](8) = a8;
  Vector::operator[](9) = a9;
  Vector::operator[](10) = a10;
  Vector::operator[](11) = a11;
  Vector::operator[](12) = a12;
  Vector::operator[](13) = a13;
  Vector::operator[](14) = a14;
  Vector::operator[](15) = a15;
  Vector::operator[](16) = a16;
  Vector::operator[](17) = a17;
  Vector::operator[](18) = a18;
  Vector::operator[](19) = a19;
  Vector::operator[](20) = a20;
  Vector::operator[](21) = a21;
  Vector::operator[](22) = a22;
  Vector::operator[](23) = a23;
  Vector::operator[](24) = a24;
  Vector::operator[](25) = a25;
  Vector::operator[](26) = a26;
  Vector::operator[](27) = a27;
  Vector::operator[](28) = a28;
  Vector::operator[](29) = a29;
}

