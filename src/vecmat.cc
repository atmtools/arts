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



////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
  \file   vecmat.cc
  
   Implementation of some MATRIX/VECTOR functions.

  \date   2001-01-08
  \author Stefan Buehler
*/



#include "vecmat.h"



////////////////////////////////////////////////////////////////////////////
//   The functions
////////////////////////////////////////////////////////////////////////////


// Output operator for VECTOR:
std::ostream& operator<<(std::ostream& s, const VECTOR& A)
{
  print_vector(s, A);
  return s;
}

// Output operator for VECTOR subrange:
std::ostream& operator<<(std::ostream& s, const VECTOR::subrange_type& A)
{
  print_vector(s, A);
  return s;
}

// Output operator for MATRIX:
std::ostream& operator<<(std::ostream& s, const MATRIX& A)
{
  print_all_matrix(s, A);
  return s;
}

// Output operator for sub MATRIX:
std::ostream& operator<<(std::ostream& s, const MATRIX::submatrix_type& A)
{
  print_all_matrix(s, A);
  return s;
}

// Output operator for SPARSE:
std::ostream& operator<<(std::ostream& s, const SPARSE& A)
{
  print_all_matrix(s, A);
  return s;
}

// Output operator for SYMMETRIC:
std::ostream& operator<<(std::ostream& s, const SYMMETRIC& A)
{
  print_all_matrix(s, A);
  return s;
}



// Input operator for VECTOR:
std::istream& operator>>(std::istream& s, VECTOR& A)
{
  read_vector_from_stream(A,s);
  return s;
}

// Input operator for MATRIX:
std::istream& operator>>(std::istream& s, MATRIX& A)
{
  read_matrix_from_stream(A,s);
  return s;
}

// Resize functions for VECTOR

void resize(VECTOR& x, INDEX n)
{
  if ( n!=x.size() )
    x = VECTOR(n);
}

// Resize function for string. This is just for consistent notation within ARTS.

void resize(string& x, INDEX n)
{
  if ( n!=x.size() )
    x.resize(n);
}

// Resize functions for MATRIX

void resize(MATRIX& x, INDEX r, INDEX c)
{
  if ( r!=x.nrows() || c!=x.ncols() )
    x = MATRIX(r,c);
}

// Resize functions for SPARSE

void resize(SPARSE& x, int r, int c)
{
  // For some reason, nrows() and ncols() for this type return int!
  if ( r!=x.nrows() || c!=x.ncols() )
    x = SPARSE(r,c);
}

// Resize functions for SYMMETRIC

void resize(SYMMETRIC& x, INDEX r, INDEX c)
{
  if ( r!=x.nrows() || c!=x.ncols() )
    x = SYMMETRIC(r,c);
}


//----------------------------------------------------------------------
// 	Transform functions
//----------------------------------------------------------------------

/** A generic transform function for vectors, which can be used to
    implement mathematical functions operating on all
    elements. Because we have this, we don't need explicit functions
    like sqrt for vectors! The type of the mathematical function is
    double (&my_func)(double). Numeric would not work here, since
    mathematical functions for float do not exist!

   \retval   y   the results of the function acting on each element of x
   \param    my_func a function (e.g., sqrt)
   \param    x   a vector
    
    \date   2000-12-27
    \author Stefan Buehler */
void transf( const VECTOR& x,
		       double (&my_func)(double),
		       VECTOR& y )
{
  assert( x.size()==y.size() );
  mtl_algo::transform(x.begin(),x.end(),y.begin(),my_func);
}

/** Return version of \ref transf.

   \return   the results of the function acting on each element of x
   \param    my_func a function (e.g., sqrt)
   \param    x   a vector
    
    \date   2000-12-27
    \author Stefan Buehler */
VECTOR transf( const VECTOR& x,
			 double (&my_func)(double) )
{
  VECTOR y(x.size());
  transf( x, my_func, y );
  return y; 
}


// Matrices:

/** A generic transform function for matrices, which can be used to
    implement mathematical functions operating on all
    elements. Because we have this, we don't need explicit functions
    like sqrt for matrices! The type of the mathematical function is
    double (&my_func)(double). Numeric would not work here, since
    mathematical functions for float do not exist!

   \retval   y   the results of the function acting on each element of x
   \param    my_func a function (e.g., sqrt)
   \param    x   a matrix
    
    \date   2000-12-27
    \author Stefan Buehler */
void transf( const MATRIX& x,
	     double (&my_func)(double),
	     MATRIX& y )
{
  // This code is adapted from the function twod_copy_default in file
  // mtl.h. The algorithm should also work for sparse
  // matrices. However, in that case the function would be applied
  // only to those elements which are occupied.

  assert( y.nrows()==x.nrows() );
  assert( y.ncols()==x.ncols() );

  MATRIX::const_iterator i;
  MATRIX::OneD::const_iterator j, jend;

  for (i = x.begin(); i != x.end(); ++i)
    {
      j = (*i).begin(); jend = (*i).end();
      for (; j != jend; ++j)
	y(j.row(),j.column()) = my_func(*j);
    }
}

/** Return version of transform for matrices.

   \return   the results of the function acting on each element of x
   \param    my_func a function (e.g., sqrt)
   \param    x   a matrix
    
    \date   2000-12-27
    \author Stefan Buehler */
MATRIX transf( const MATRIX& x,
	       double (&my_func)(double) )
{
  MATRIX y( x.nrows(), x.ncols() );
  transf( x, my_func, y );
  return y;
}
