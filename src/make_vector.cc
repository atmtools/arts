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

// For VECTORs:
// ------------

VECTOR make_vector()
{
  VECTOR v(0);

  return v;
}

VECTOR make_vector(const Numeric& a1)
{
  VECTOR v(1);

  v(1) = a1;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2)
{
  VECTOR v(2);

  v(1) = a1;
  v(2) = a2;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3)
{
  VECTOR v(3);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4)
{
  VECTOR v(4);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5)
{
  VECTOR v(5);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6)
{
  VECTOR v(6);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7)
{
  VECTOR v(7);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8)
{
  VECTOR v(8);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9)
{
  VECTOR v(9);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;
  v(9) = a9;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10)
{
  VECTOR v(10);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;
  v(9) = a9;
  v(10) = a10;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11)
{
  VECTOR v(11);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;
  v(9) = a9;
  v(10) = a10;
  v(11) = a11;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12)
{
  VECTOR v(12);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;
  v(9) = a9;
  v(10) = a10;
  v(11) = a11;
  v(12) = a12;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13)
{
  VECTOR v(13);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;
  v(9) = a9;
  v(10) = a10;
  v(11) = a11;
  v(12) = a12;
  v(13) = a13;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14)
{
  VECTOR v(14);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;
  v(9) = a9;
  v(10) = a10;
  v(11) = a11;
  v(12) = a12;
  v(13) = a13;
  v(14) = a14;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14,
		   const Numeric& a15)
{
  VECTOR v(15);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;
  v(9) = a9;
  v(10) = a10;
  v(11) = a11;
  v(12) = a12;
  v(13) = a13;
  v(14) = a14;
  v(15) = a15;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14,
		   const Numeric& a15,
		   const Numeric& a16)
{
  VECTOR v(16);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;
  v(9) = a9;
  v(10) = a10;
  v(11) = a11;
  v(12) = a12;
  v(13) = a13;
  v(14) = a14;
  v(15) = a15;
  v(16) = a16;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14,
		   const Numeric& a15,
		   const Numeric& a16,
		   const Numeric& a17)
{
  VECTOR v(17);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;
  v(9) = a9;
  v(10) = a10;
  v(11) = a11;
  v(12) = a12;
  v(13) = a13;
  v(14) = a14;
  v(15) = a15;
  v(16) = a16;
  v(17) = a17;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14,
		   const Numeric& a15,
		   const Numeric& a16,
		   const Numeric& a17,
		   const Numeric& a18)
{
  VECTOR v(18);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;
  v(9) = a9;
  v(10) = a10;
  v(11) = a11;
  v(12) = a12;
  v(13) = a13;
  v(14) = a14;
  v(15) = a15;
  v(16) = a16;
  v(17) = a17;
  v(18) = a18;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14,
		   const Numeric& a15,
		   const Numeric& a16,
		   const Numeric& a17,
		   const Numeric& a18,
		   const Numeric& a19)
{
  VECTOR v(19);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;
  v(9) = a9;
  v(10) = a10;
  v(11) = a11;
  v(12) = a12;
  v(13) = a13;
  v(14) = a14;
  v(15) = a15;
  v(16) = a16;
  v(17) = a17;
  v(18) = a18;
  v(19) = a19;

  return v;
}

VECTOR make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11,
		   const Numeric& a12,
		   const Numeric& a13,
		   const Numeric& a14,
		   const Numeric& a15,
		   const Numeric& a16,
		   const Numeric& a17,
		   const Numeric& a18,
		   const Numeric& a19,
		   const Numeric& a20)
{
  VECTOR v(20);

  v(1) = a1;
  v(2) = a2;
  v(3) = a3;
  v(4) = a4;
  v(5) = a5;
  v(6) = a6;
  v(7) = a7;
  v(8) = a8;
  v(9) = a9;
  v(10) = a10;
  v(11) = a11;
  v(12) = a12;
  v(13) = a13;
  v(14) = a14;
  v(15) = a15;
  v(16) = a16;
  v(17) = a17;
  v(18) = a18;
  v(19) = a19;
  v(20) = a20;

  return v;
}



// For ARRAYs:
// -----------

// template<class T>
// ARRAY<T> make_array()
// {
//   ARRAY<T> v(0,T());

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1)
// {
//   ARRAY<T> v(1,T());

//   v[0] = a1;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2)
// {
//   ARRAY<T> v(2,T());

//   v[0] = a1;
//   v[1] = a2;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3)
// {
//   ARRAY<T> v(3,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4)
// {
//   ARRAY<T> v(4,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5)
// {
//   ARRAY<T> v(5,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6)
// {
//   ARRAY<T> v(6,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7)
// {
//   ARRAY<T> v(7,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8)
// {
//   ARRAY<T> v(8,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9)
// {
//   ARRAY<T> v(9,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;
//   v[8] = a9;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10)
// {
//   ARRAY<T> v(10,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;
//   v[8] = a9;
//   v[9] = a10;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11)
// {
//   ARRAY<T> v(11,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;
//   v[8] = a9;
//   v[9] = a10;
//   v[10] = a11;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12)
// {
//   ARRAY<T> v(12,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;
//   v[8] = a9;
//   v[9] = a10;
//   v[10] = a11;
//   v[11] = a12;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13)
// {
//   ARRAY<T> v(13,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;
//   v[8] = a9;
//   v[9] = a10;
//   v[10] = a11;
//   v[11] = a12;
//   v[12] = a13;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14)
// {
//   ARRAY<T> v(14,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;
//   v[8] = a9;
//   v[9] = a10;
//   v[10] = a11;
//   v[11] = a12;
//   v[12] = a13;
//   v[13] = a14;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14,
// 		    const T& a15)
// {
//   ARRAY<T> v(15,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;
//   v[8] = a9;
//   v[9] = a10;
//   v[10] = a11;
//   v[11] = a12;
//   v[12] = a13;
//   v[13] = a14;
//   v[14] = a15;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14,
// 		    const T& a15,
// 		    const T& a16)
// {
//   ARRAY<T> v(16,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;
//   v[8] = a9;
//   v[9] = a10;
//   v[10] = a11;
//   v[11] = a12;
//   v[12] = a13;
//   v[13] = a14;
//   v[14] = a15;
//   v[15] = a16;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14,
// 		    const T& a15,
// 		    const T& a16,
// 		    const T& a17)
// {
//   ARRAY<T> v(17,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;
//   v[8] = a9;
//   v[9] = a10;
//   v[10] = a11;
//   v[11] = a12;
//   v[12] = a13;
//   v[13] = a14;
//   v[14] = a15;
//   v[15] = a16;
//   v[16] = a17;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14,
// 		    const T& a15,
// 		    const T& a16,
// 		    const T& a17,
// 		    const T& a18)
// {
//   ARRAY<T> v(18,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;
//   v[8] = a9;
//   v[9] = a10;
//   v[10] = a11;
//   v[11] = a12;
//   v[12] = a13;
//   v[13] = a14;
//   v[14] = a15;
//   v[15] = a16;
//   v[16] = a17;
//   v[17] = a18;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14,
// 		    const T& a15,
// 		    const T& a16,
// 		    const T& a17,
// 		    const T& a18,
// 		    const T& a19)
// {
//   ARRAY<T> v(19,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;
//   v[8] = a9;
//   v[9] = a10;
//   v[10] = a11;
//   v[11] = a12;
//   v[12] = a13;
//   v[13] = a14;
//   v[14] = a15;
//   v[15] = a16;
//   v[16] = a17;
//   v[17] = a18;
//   v[18] = a19;

//   return v;
// }

// template<class T>
// ARRAY<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11,
// 		    const T& a12,
// 		    const T& a13,
// 		    const T& a14,
// 		    const T& a15,
// 		    const T& a16,
// 		    const T& a17,
// 		    const T& a18,
// 		    const T& a19,
// 		    const T& a20)
// {
//   ARRAY<T> v(20,T());

//   v[0] = a1;
//   v[1] = a2;
//   v[2] = a3;
//   v[3] = a4;
//   v[4] = a5;
//   v[5] = a6;
//   v[6] = a7;
//   v[7] = a8;
//   v[8] = a9;
//   v[9] = a10;
//   v[10] = a11;
//   v[11] = a12;
//   v[12] = a13;
//   v[13] = a14;
//   v[14] = a15;
//   v[15] = a16;
//   v[16] = a17;
//   v[17] = a18;
//   v[18] = a19;
//   v[19] = a20;

//   return v;
// }
