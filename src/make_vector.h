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

/*!
  \file   make_vector.h
  \brief  The make_vector functions can be used to initialize
          VectorS in the source code. 

  Usage is as simple as:

  \code
  Vector x        = make_vector(33.2, 17.3, 3.8);
  \endcode

  which will generate a Vector with the desired arguments. At the
  moment, up to 20 arguments are possible. If more are needed this
  can be easily extended.

  \author Stefan Buehler
  \date   1999-08-02

  Moved make_array to separate file.
  \author Stefan Buehler
  \date   2000-03-28
*/

#ifndef make_vector_h
#define make_vector_h

#include "vecmat.h"

// For Vectors:
// ------------

Vector make_vector();

Vector make_vector(const Numeric& a1);

Vector make_vector(const Numeric& a1,
		   const Numeric& a2);

Vector make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3);

Vector make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4);

Vector make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5);

Vector make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6);

Vector make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7);

Vector make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8);

Vector make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9);

Vector make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10);

Vector make_vector(const Numeric& a1,
		   const Numeric& a2,
		   const Numeric& a3,
		   const Numeric& a4,
		   const Numeric& a5,
		   const Numeric& a6,
		   const Numeric& a7,
		   const Numeric& a8,
		   const Numeric& a9,
		   const Numeric& a10,
		   const Numeric& a11);

Vector make_vector(const Numeric& a1,
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
		   const Numeric& a12);

Vector make_vector(const Numeric& a1,
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
		   const Numeric& a13);

Vector make_vector(const Numeric& a1,
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
		   const Numeric& a14);

Vector make_vector(const Numeric& a1,
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
		   const Numeric& a15);

Vector make_vector(const Numeric& a1,
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
		   const Numeric& a16);

Vector make_vector(const Numeric& a1,
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
		   const Numeric& a17);

Vector make_vector(const Numeric& a1,
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
		   const Numeric& a18);

Vector make_vector(const Numeric& a1,
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
		   const Numeric& a19);

Vector make_vector(const Numeric& a1,
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
		   const Numeric& a20);



// For Arrays:
// -----------

// template<class T>
// Array<T> make_array();

// template<class T>
// Array<T> make_array(const T& a1);

// template<class T>
// Array<T> make_array(const T& a1,
// 		    const T& a2);

// template<class T>
// Array<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3);

// template<class T>
// Array<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4);

// template<class T>
// Array<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5);

// template<class T>
// Array<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6);

// template<class T>
// Array<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7);

// template<class T>
// Array<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8);

// template<class T>
// Array<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9);

// template<class T>
// Array<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10);

// template<class T>
// Array<T> make_array(const T& a1,
// 		    const T& a2,
// 		    const T& a3,
// 		    const T& a4,
// 		    const T& a5,
// 		    const T& a6,
// 		    const T& a7,
// 		    const T& a8,
// 		    const T& a9,
// 		    const T& a10,
// 		    const T& a11);

// template<class T>
// Array<T> make_array(const T& a1,
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
// 		    const T& a12);

// template<class T>
// Array<T> make_array(const T& a1,
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
// 		    const T& a13);

// template<class T>
// Array<T> make_array(const T& a1,
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
// 		    const T& a14);

// template<class T>
// Array<T> make_array(const T& a1,
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
// 		    const T& a15);

// template<class T>
// Array<T> make_array(const T& a1,
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
// 		    const T& a16);

// template<class T>
// Array<T> make_array(const T& a1,
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
// 		    const T& a17);

// template<class T>
// Array<T> make_array(const T& a1,
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
// 		    const T& a18);

// template<class T>
// Array<T> make_array(const T& a1,
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
// 		    const T& a19);

// template<class T>
// Array<T> make_array(const T& a1,
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
// 		    const T& a20);


#endif  // make_vector_h
