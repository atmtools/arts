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

/*  */

/*!
  \file   make_array.h
  \brief  The make_array objects can be used to initialize
          ARRAYS in the source code. 

  Usage is as simple as:

  \code
     ARRAY<string> y = make_array<string>("ab", "cde");
  \endcode

  which will make an ARRAY with the desired arguments. At the
  moment, up to 20 arguments are possible. If more are needed this
  can be easily extended.

  This has to be solved with a class, since template *functions*
  would have to be explicitly instantiated (because we want to use
  make_array in different source files). Putting functions
  definitions in a .h file does not work (leads to linker errors).
   
  \author Stefan Buehler
  \date   2000-03-28
*/

#ifndef make_array_h
#define make_array_h

#include "vecmat.h"

template<class T>
class make_array {
 public:
  // All the constructors for the different dimensions:
  make_array() : v(0,T()) {}
  make_array(const T& a1) : v(1,T()) {
    v[0] = a1;
  }
  make_array(const T& a1,
	     const T& a2) : v(2,T()) {
    v[0] = a1;
    v[1] = a2;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3) : v(3,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4) : v(4,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5) : v(5,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6) : v(6,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7) : v(7,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8) : v(8,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9) : v(9,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10) : v(10,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11) : v(11,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12) : v(12,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13) : v(13,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14) : v(14,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15) : v(15,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16) : v(16,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17) : v(17,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18) : v(18,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18,
	     const T& a19) : v(19,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
    v[18] = a19;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18,
	     const T& a19,
	     const T& a20) : v(20,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
    v[18] = a19;
    v[19] = a20;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18,
	     const T& a19,
	     const T& a20,
	     const T& a21) : v(21,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
    v[18] = a19;
    v[19] = a20;
    v[20] = a21;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18,
	     const T& a19,
	     const T& a20,
	     const T& a21,
	     const T& a22) : v(22,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
    v[18] = a19;
    v[19] = a20;
    v[20] = a21;
    v[21] = a22;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18,
	     const T& a19,
	     const T& a20,
	     const T& a21,
	     const T& a22,
	     const T& a23) : v(23,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
    v[18] = a19;
    v[19] = a20;
    v[20] = a21;
    v[21] = a22;
    v[22] = a23;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18,
	     const T& a19,
	     const T& a20,
	     const T& a21,
	     const T& a22,
	     const T& a23,
             const T& a24) : v(24,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
    v[18] = a19;
    v[19] = a20;
    v[20] = a21;
    v[21] = a22;
    v[22] = a23;
    v[23] = a24;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18,
	     const T& a19,
	     const T& a20,
	     const T& a21,
	     const T& a22,
	     const T& a23,
             const T& a24,
	     const T& a25) : v(25,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
    v[18] = a19;
    v[19] = a20;
    v[20] = a21;
    v[21] = a22;
    v[22] = a23;
    v[23] = a24;
    v[24] = a25;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18,
	     const T& a19,
	     const T& a20,
	     const T& a21,
	     const T& a22,
	     const T& a23,
             const T& a24,
             const T& a25,
	     const T& a26) : v(26,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
    v[18] = a19;
    v[19] = a20;
    v[20] = a21;
    v[21] = a22;
    v[22] = a23;
    v[23] = a24;
    v[24] = a25;
    v[25] = a26;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18,
	     const T& a19,
	     const T& a20,
	     const T& a21,
	     const T& a22,
	     const T& a23,
             const T& a24,
             const T& a25,
             const T& a26,
	     const T& a27) : v(27,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
    v[18] = a19;
    v[19] = a20;
    v[20] = a21;
    v[21] = a22;
    v[22] = a23;
    v[23] = a24;
    v[24] = a25;
    v[25] = a26;
    v[26] = a27;
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18,
	     const T& a19,
	     const T& a20,
	     const T& a21,
	     const T& a22,
	     const T& a23,
             const T& a24,
             const T& a25,
             const T& a26,
             const T& a27,
	     const T& a28) : v(28,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
    v[18] = a19;
    v[19] = a20;
    v[20] = a21;
    v[21] = a22;
    v[22] = a23;
    v[23] = a24;
    v[24] = a25;
    v[25] = a26;
    v[26] = a27;
    v[27] = a28;  
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18,
	     const T& a19,
	     const T& a20,
	     const T& a21,
	     const T& a22,
	     const T& a23,
             const T& a24,
             const T& a25,
             const T& a26,
             const T& a27,
             const T& a28,
	     const T& a29) : v(29,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
    v[18] = a19;
    v[19] = a20;
    v[20] = a21;
    v[21] = a22;
    v[22] = a23;
    v[23] = a24;
    v[24] = a25;
    v[25] = a26;
    v[26] = a27;
    v[27] = a28;  
    v[28] = a29;  
  }
  make_array(const T& a1,
	     const T& a2,
	     const T& a3,
	     const T& a4,
	     const T& a5,
	     const T& a6,
	     const T& a7,
	     const T& a8,
	     const T& a9,
	     const T& a10,
	     const T& a11,
	     const T& a12,
	     const T& a13,
	     const T& a14,
	     const T& a15,
	     const T& a16,
	     const T& a17,
	     const T& a18,
	     const T& a19,
	     const T& a20,
	     const T& a21,
	     const T& a22,
	     const T& a23,
             const T& a24,
             const T& a25,
             const T& a26,
             const T& a27,
             const T& a28,
             const T& a29,
	     const T& a30) : v(30,T()) {
    v[0] = a1;
    v[1] = a2;
    v[2] = a3;
    v[3] = a4;
    v[4] = a5;
    v[5] = a6;
    v[6] = a7;
    v[7] = a8;
    v[8] = a9;
    v[9] = a10;
    v[10] = a11;
    v[11] = a12;
    v[12] = a13;
    v[13] = a14;
    v[14] = a15;
    v[15] = a16;
    v[16] = a17;
    v[17] = a18;
    v[18] = a19;
    v[19] = a20;
    v[20] = a21;
    v[21] = a22;
    v[22] = a23;
    v[23] = a24;
    v[24] = a25;
    v[25] = a26;
    v[26] = a27;
    v[27] = a28;  
    v[28] = a29;  
    v[29] = a30;  
  }

  // Conversion operator to ARRAY:
  operator ARRAY<T>() const {
    return v;
  }
 private:
  ARRAY<T> v;
};


#endif  // make_array_h
