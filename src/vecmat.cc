/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>
                      Patrick Eriksson <patrick@rss.chalmers.se>

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
  \file   vecmat.cc
  \brief  Implementation of some MATRIX/VECTOR functions.

  \author Stefan Buehler
  \date   2000-08-16
*/

#include "vecmat.h"

MATRIX emult(const MATRIX &A, const MATRIX &B)
  {
    SUBSCRIPT M = A.dim(1);
    SUBSCRIPT N = A.dim(2);

    assert(M==B.dim(1));
    assert(N==B.dim(2));

    MATRIX tmp(M,N);
    SUBSCRIPT i,j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
	  tmp(i,j) = A(i,j) * B(i,j);

    return tmp;
  }        

MATRIX ediv(const MATRIX &A, const MATRIX &B)
  {
    SUBSCRIPT M = A.dim(1);
    SUBSCRIPT N = A.dim(2);

    assert(M==B.dim(1));
    assert(N==B.dim(2));

    MATRIX tmp(M,N);
    SUBSCRIPT i,j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            tmp(i,j) = A(i,j) / B(i,j);

    return tmp;
  }        

  VECTOR emult( const VECTOR &A, const VECTOR &B)
  {
    SUBSCRIPT N = A.dim();

    assert(N==B.dim());

    VECTOR tmp(N);
    SUBSCRIPT i;

    for (i=1; i<=N; i++)
            tmp(i) = A(i) * B(i);

    return tmp;
  }       

  VECTOR ediv( const VECTOR &A, const VECTOR &B)
  {
    SUBSCRIPT N = A.dim();

    assert(N==B.dim());

    VECTOR tmp(N);
    SUBSCRIPT i;

    for (i=1; i<=N; i++)
            tmp(i) = A(i) / B(i);

    return tmp;
  }       

  MATRIX operator+(const MATRIX &A, const Numeric scalar)
  {
    SUBSCRIPT M = A.dim(1);
    SUBSCRIPT N = A.dim(2);

    MATRIX tmp(M,N);
    SUBSCRIPT i,j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            tmp(i,j) = A(i,j) + scalar;

    return tmp;
  }

  MATRIX operator+(const Numeric scalar, const MATRIX &A)
  {
    return ( A + scalar );
  }                 

// - 
//
  MATRIX operator-(const MATRIX &A, const Numeric scalar)
  {
    SUBSCRIPT M = A.dim(1);
    SUBSCRIPT N = A.dim(2);

    MATRIX tmp(M,N);
    SUBSCRIPT i,j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            tmp(i,j) = A(i,j) - scalar;

    return tmp;
  }

  MATRIX operator-(const Numeric scalar, const MATRIX &A)
  {
    SUBSCRIPT M = A.dim(1);
    SUBSCRIPT N = A.dim(2);

    MATRIX tmp(M,N);
    SUBSCRIPT i,j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            tmp(i,j) =  scalar - A(i,j);

    return tmp;
  }


// *
//
  MATRIX operator*(const MATRIX &A, const Numeric scalar)
  {
    SUBSCRIPT M = A.dim(1);
    SUBSCRIPT N = A.dim(2);

    MATRIX tmp(M,N);
    SUBSCRIPT i,j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            tmp(i,j) = A(i,j) * scalar;

    return tmp;
  }

  MATRIX operator*(const Numeric scalar, const MATRIX &A)
  {
    return ( A * scalar );
  }                 

//
// /
  MATRIX operator/(const MATRIX &A, const Numeric scalar)
  {
    SUBSCRIPT M = A.dim(1);
    SUBSCRIPT N = A.dim(2);

    MATRIX tmp(M,N);
    SUBSCRIPT i,j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            tmp(i,j) = A(i,j) / scalar;

    return tmp;
  }                       

  MATRIX operator/(const Numeric scalar, const MATRIX &A)
  {
    SUBSCRIPT M = A.dim(1);
    SUBSCRIPT N = A.dim(2);

    MATRIX tmp(M,N);
    SUBSCRIPT i,j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            tmp(i,j) = scalar / A(i,j);

    return tmp; 
  }        


// VECTOR - SCALAR ---------------------------------------------------------

// + 
//
  VECTOR operator+(const VECTOR &A, const Numeric scalar)
  {
    SUBSCRIPT N = A.dim();

    VECTOR tmp(N);
    SUBSCRIPT i;

    for (i=1; i<=N; i++)
      tmp(i) = A(i) + scalar;

    return tmp;
  }

  VECTOR operator+(const Numeric scalar, const VECTOR &A)
  {
    return ( A + scalar );
  }

// -
//
  VECTOR operator-(const VECTOR &A, const Numeric scalar)
  {
    SUBSCRIPT N = A.dim();

    VECTOR tmp(N);
    SUBSCRIPT i;

    for (i=1; i<=N; i++)
      tmp(i) = A(i) - scalar;

    return tmp;
  }

  VECTOR operator-(const Numeric scalar, const VECTOR &A)
  {
    SUBSCRIPT N = A.dim();

    VECTOR tmp(N);
    SUBSCRIPT i;

    for (i=1; i<=N; i++)
      tmp(i) = scalar - A(i);

    return tmp;
  }

// *
//
  VECTOR operator*(const VECTOR &A, const Numeric scalar)
  {
    SUBSCRIPT N = A.dim();

    VECTOR tmp(N);
    SUBSCRIPT i;

    for (i=1; i<=N; i++)
      tmp(i) = A(i) * scalar;

    return tmp;
  }

  VECTOR operator*(const Numeric scalar, const VECTOR &A)
  {
    return ( A * scalar );
  }

// /
//
  VECTOR operator/(const VECTOR &A, const Numeric scalar)
  {
    SUBSCRIPT N = A.dim();

    VECTOR tmp(N);
    SUBSCRIPT i;

    for (i=1; i<=N; i++)
      tmp(i) = A(i) / scalar;

    return tmp;
  }

  VECTOR operator/(const Numeric scalar, const VECTOR &A)
  {
    SUBSCRIPT N = A.dim();

    VECTOR tmp(N);
    SUBSCRIPT i;

    for (i=1; i<=N; i++)
      tmp(i) = scalar / A(i);

    return tmp;
  }
