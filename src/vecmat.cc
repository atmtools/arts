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
  
   Implementation of some Matrix/Vector functions.

  \date   2001-01-08
  \author Stefan Buehler
*/



#include "vecmat.h"



////////////////////////////////////////////////////////////////////////////
//   The functions
////////////////////////////////////////////////////////////////////////////


// Output operator for Vector:
std::ostream& operator<<(std::ostream& s, const Vector& A)
{
  print_vector(s, A);
  return s;
}

// Output operator for Vector subrange:
std::ostream& operator<<(std::ostream& s, const Vector::subrange_type& A)
{
  print_vector(s, A);
  return s;
}

// Output operator for Matrix:
std::ostream& operator<<(std::ostream& s, const Matrix& A)
{
  print_all_matrix(s, A);
  return s;
}

// Output operator for sub Matrix:
std::ostream& operator<<(std::ostream& s, const Matrix::submatrix_type& A)
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



// Input operator for Vector:
std::istream& operator>>(std::istream& s, Vector& A)
{
  read_vector_from_stream(A,s);
  return s;
}

// Input operator for Matrix:
std::istream& operator>>(std::istream& s, Matrix& A)
{
  read_matrix_from_stream(A,s);
  return s;
}

// Resize functions for Vector

void resize(Vector& x, Index n)
{
  if ( n!=x.size() )
    x = Vector(n);
}

// Resize function for string. This is just for consistent notation within ARTS.

void resize(string& x, Index n)
{
  if ( n!=x.size() )
    x.resize(n);
}

// Resize functions for Matrix

void resize(Matrix& x, Index r, Index c)
{
  if ( r!=x.nrows() || c!=x.ncols() )
    x = Matrix(r,c);
}

// Resize functions for SPARSE

void resize(SPARSE& x, int r, int c)
{
  // For some reason, nrows() and ncols() for this type return int!
  if ( r!=x.nrows() || c!=x.ncols() )
    x = SPARSE(r,c);
}

// Resize functions for SYMMETRIC

void resize(SYMMETRIC& x, Index r, Index c)
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
void transf( const Vector& x,
		       double (&my_func)(double),
		       Vector& y )
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
Vector transf( const Vector& x,
			 double (&my_func)(double) )
{
  Vector y(x.size());
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
void transf( const Matrix& x,
	     double (&my_func)(double),
	     Matrix& y )
{
  // This code is adapted from the function twod_copy_default in file
  // mtl.h. The algorithm should also work for sparse
  // matrices. However, in that case the function would be applied
  // only to those elements which are occupied.

  assert( y.nrows()==x.nrows() );
  assert( y.ncols()==x.ncols() );

  Matrix::const_iterator i;
  Matrix::OneD::const_iterator j, jend;

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
Matrix transf( const Matrix& x,
	       double (&my_func)(double) )
{
  Matrix y( x.nrows(), x.ncols() );
  transf( x, my_func, y );
  return y;
}
