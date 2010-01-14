/* Copyright (C) 2000, 2001 Patrick Eriksson <patrick@rss.chalmers.se>
                            Stefan Buehler   <sbuehler@uni-bremen.de>

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
   \file   math_funcs.cc

   Contains the code of basic mathematical and vector/matrix functions.

   Example on types of functions:
   \begin{enumerate}
    \item First and last element of a vector
    \item Boolean functions                         
    \item Creation of common vectors                
    \item Interpolation routines            
    \item Check of function input
   \end{enumerate}

   \author Patrick Eriksson
   \date 2000-09-18 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include <time.h>
#include <math.h>
#include <stdexcept>
#include "arts.h"
#include "math_funcs.h"
#include "make_vector.h"
#include "array.h"



////////////////////////////////////////////////////////////////////////////
//// first and last /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//
/** Gives the first value of a vector.

    \return       the first value of x
    \param    x   a vector

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric first( ConstVectorView x )
{
  return x[0]; 
}

/** Gives the last value of a vector.

    \return      the last value of x
    \param   x   a vector

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric last( ConstVectorView x )
{
  return x[x.nelem()-1]; 
}



////////////////////////////////////////////////////////////////////////////
//   Logical functions
////////////////////////////////////////////////////////////////////////////

//// any ///////////////////////////////////////////////////////////////////
//
/** True if any element of a boolean vector, b, is not 0.

    \return       a boolean, true if any x is != 0
    \param    x   a vector

    \author Patrick Eriksson
    \date   2000-06-27
*/
bool any( const ArrayOfIndex& x ) 
{
  for ( Index i=0; i<x.nelem(); i++ ) {
    if ( x[i] )
      return true;
  }
  return false;
}




////////////////////////////////////////////////////////////////////////////
//   Functions to generate vectors
////////////////////////////////////////////////////////////////////////////

//// linspace //////////////////////////////////////////////////////////////
//
/** Linearly spaced vector with specified spacing. 

    The first element of x is always start. The next value is start+step etc.
    Note that the last value can deviate from stop.
    The step can be both positive and negative. 
    (in Matlab notation: start:step:stop)

    Size of result is adjusted within this function!

    \retval   x       linearly spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    step    distance between values in x

    \author Patrick Eriksson
    \date   2000-06-27
*/
void linspace(                      
              Vector&     x,           
              const Numeric     start,    
              const Numeric     stop,        
              const Numeric     step )
{
  Index n = (Index) floor( (stop-start)/step ) + 1;
  if ( n<1 )
    n=1;
  x.resize(n);
  for ( Index i=0; i<n; i++ )
    x[i] = start + (Numeric)i*step;
}

/** Linearly spaced vector with specified length. 

    Returns a vector equally and linearly spaced between start and stop 
    of length n. (equals the Matlab function linspace)

    The length must be > 1.

    \retval   x       linearly spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    n       length of x

    \author Patrick Eriksson
    \date   2000-06-27
*/
void nlinspace(         
               Vector&     x, 
               const Numeric     start,     
               const Numeric     stop,        
               const Index       n )
{
  assert( 1<n );                // Number of points must be greatere 1.
  x.resize(n);
  Numeric step = (stop-start)/(Numeric)(n-1) ;
  for ( Index i=0; i<n; i++ )
    x[i] = start + (Numeric)i*step;
}

//// nlogspace /////////////////////////////////////////////////////////////
//
/** Logarithmically spaced vector with specified length. 

    Returns a vector logarithmically spaced vector between start and 
    stop of length n (equals the Matlab function logspace)

    The length must be > 1.

    \retval   x       logarithmically spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    n       length of x

    \author Patrick Eriksson
    \date   2000-06-27
*/
void nlogspace(         
               Vector&     x, 
               const Numeric     start,     
               const Numeric     stop,        
               const Index         n )
{
  // Number of points must be greatere 1:
  assert( 1<n );        
  // Only positive numbers are allowed for start and stop:
  assert( 0<start );
  assert( 0<stop );

  x.resize(n);
  Numeric a = log(start);
  Numeric step = (log(stop)-a)/(Numeric)(n-1);
  x[0] = start;
  for ( Index i=1; i<n-1; i++ )
    x[i] = exp(a + (Numeric)i*step);
  x[n-1] = stop;
}

/** Logarithmically spaced vector with specified length (return version). 

    Returns a vector logarithmically spaced vector between start and 
    stop of length n (equals the Matlab function logspace)

    The length must be > 1.

    \return   x       logarithmically spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    n       length of x

    \author Patrick Eriksson
    \date   2000-06-27
*/
Vector nlogspace(  
                 const Numeric start, 
                 const Numeric stop,  
                 const Index     n )
{
  Vector x;
  nlogspace( x, start, stop, n );
  return x; 
}                     




////////////////////////////////////////////////////////////////////////////
//   Interpolation routines
////////////////////////////////////////////////////////////////////////////

/** Local help function to check input grids.

    Tolerates the new grid to extend half a grid distance outside the
    original grid.

    \date   2000-06-29, 2003-08-06
    \author Patrick Eriksson, Stefan Buehler 
*/
Index interp_check( ConstVectorView  x,        
                    ConstVectorView  xi,
                    const Index   n_y )
{
  const Index  n  = x.nelem();
  const Index  ni = xi.nelem();
  Index  order=1;          // flag for pos. or neg. order, assume pos.

  // Determine the order, -1=decreasing and 1=increasing
  if ( x[0] > x[n-1] )
    order = -1;

  if ( n < 2 )
    throw runtime_error("Vector length for interpolation must be >= 2.");

  if ( n != n_y ) 
    throw runtime_error("Sizes of input data to interpolation do not match.");

  // Make lower and upper bounds tolerate the new point to be half a
  // grid distance outside the original grid.
  const Numeric lower_bound = x[0]   - 0.5*(x[1]-x[0]);
  const Numeric upper_bound = x[n-1] + 0.5*(x[n-1]-x[n-2]);

  if ( ((Numeric)order*xi[0]<(Numeric)order*lower_bound)
       || ((Numeric)order*xi[ni-1]>(Numeric)order*upper_bound) ) 
    {
      ostringstream os;
      os << "Interpolation points must be not more than\n"
         << "half a grid spacing outside the original range.\n"
         << "Int.:  xi[0] = " << xi[0] << ", xi[ni-1] = " << xi[ni-1] << '\n'
         << "Orig.: x[0]  = " << x[0]  << ", x[n-1]   = " << x[n-1];
      throw runtime_error(os.str());
    }

  for ( Index i=0; i<n-1; i++ )
  {
    if ( (Numeric)order*x[i+1] < (Numeric)order*x[i] ) 
      throw runtime_error("Original interpolation grid must be ordered");
  }

  for ( Index i=0; i<ni-1; i++ )
  {
    if ( (Numeric)order*xi[i+1] < (Numeric)order*xi[i] ) 
      throw runtime_error("Interpolation points must be ordered");
  }

  return order;
}

/** Multiple linear interpolation of a vector. Because Views are used,
    you can use matrix rows or columns, or vector sub-ranges inside
    the argument of this function.

    The vector x specifies the points at which the data y is given. It
    must be strictly monotonically increasing. (No two x values must
    be the same.)

    The size of yi has to be the same as for xi.

    Interpolation works also for new grid points just outside the
    original grid range.

    \retval  yi      interpolated values 
    \param   x       the x grid
    \param   y       the function to interpolate
    \param   xi      interpolation points

    \date   2001-01-05, 2003-08-06
    \author Stefan Buehler, Stefan Buehler
*/
void interp_lin_vector( VectorView       yi,
                        ConstVectorView  x, 
                        ConstVectorView  y, 
                        ConstVectorView  xi )
{
  // Check grids and get order of grids
  Index order = interp_check( x, xi, y.nelem() ); 

  Index        i, j=0;
  const Index  n=xi.nelem();
  const Index  nx=x.nelem();
  Numeric      w;

  // Check that output vector has the right size:
  assert( n==yi.nelem() ); 

  for ( i=0; i<n; i++ )
  {
    // Find right place:
    while ( j<nx-2 && (Numeric)order*x[j+1]<(Numeric)order*xi[i] )
      {
        ++j;
      }
    // j should now point to the place in the old grid below the
    // interpolation point. If the interpolation point is below the
    // original grid, j=0. If the interpolation point is above the
    // original grid, j=nx-2.
    
    assert( x[j+1]!=x[j] );
    w = (xi[i]-x[j]) / (x[j+1]-x[j]);
    // This expression should also be correct for points just outside
    // the original grid. In that case w can be negative or larger
    // than 1.
    yi[i] = y[j] + w * (y[j+1]-y[j]); 
  }
}      

/** Multiple linear interpolation of matrix rows. Works with subranges
    inside the argument. Use the transpose(A) function if you want to
    interpolate columns of A instead of rows.

    The vector x specifies the points at which the data y is given. 

    Interpolation works also for new grid points just outside the
    original grid range.

    \retval  Yi      interpolated values (matrix)
    \param   x       the x grid (vector)
    \param   Y       the function to interpolate (matrix)
    \param   xi      interpolation points (vector)

    \date   2001-01-06, 2003-08-06
    \author Stefan Buehler, Stefan Buehler
*/
void interp_lin_matrix(    
                       MatrixView        Yi,
                       ConstVectorView   x, 
                       ConstMatrixView   Y, 
                       ConstVectorView   xi )
{
  // Check grids and get order of grids
  Index order = interp_check( x, xi, Y.ncols() ); 

  Index j=0;
  const Index n=xi.nelem(), nrow=Y.nrows(), nx=x.nelem();
  Numeric w;

  assert( nrow == Yi.nrows() );
  assert( n    == Yi.ncols() );

  for (Index i=0; i<n; i++ )
  {
    // Find right place:
    while ( j<nx-2 && (Numeric)order*x[j+1]<(Numeric)order*xi[i] )
      {
        ++j;
      }
    // j should now point to the place in the old grid below the
    // interpolation point. If the interpolation point is below the
    // original grid, j=0. If the interpolation point is above the
    // original grid, j=nx-2.

    assert( x[j+1]!=x[j] );
    w = (xi[i]-x[j]) / (x[j+1]-x[j]);
    // This expression should also be correct for points just outside
    // the original grid. In that case w can be negative or larger
    // than 1.
    for( Index k=0; k<nrow; k++ )
    {
      Yi(k,i) = Y(k,j) + w * (Y(k,j+1)-Y(k,j));
    }
  }
}        

/** Single linear interpolation of a vector (return version).

    The vector x specifies the points at which the data y is given. 

    \return          interpolated value
    \param   x       the x grid
    \param   y       the function to interpolate
    \param   xi      interpolation point

    \author Patrick Eriksson
    \date   2000-06-29
*/
Numeric interp_lin(
                   ConstVectorView  x, 
                   ConstVectorView  y, 
                   const Numeric  xi )
{
  Vector Yi(1);
  MakeVector Xi(xi);   // MakeVector is a special kind of vector that
                       // can be initialized explicitly with one or
                       // more arguments of type Numeric. 

  interp_lin_vector( Yi, x, y, Xi );
  return Yi[0];
}        




/////////////////////////////////////////////////////////////////////////////
//   Check of function input
/////////////////////////////////////////////////////////////////////////////

//// check_if_bool ////////////////////////////////////////////////////////////
//
/** Checks if an integer is 0 or 1.

    A runtime error is thrown if the integer is not a boolean.

    \param    x        an integer
    \param    x_name   the name of the integer

    \author Patrick Eriksson
    \date   2001-09-19
*/
void check_if_bool( const Index& x, const String& x_name ) 
{
  if ( !(x==0 || x==1) )
  {
    ostringstream os;
    os << "The boolean *" << x_name <<  "* must either be 0 or 1.\n" 
       << "The present value of *"<< x_name <<  "* is " << x << ".";
    throw runtime_error( os.str() );
  }
}



//// check_if_in_range ////////////////////////////////////////////////////////
//
/** Checks if a numeric variable is inside a specified range.

    A runtime error is thrown if the variable is outside the range.

    \param    x_low    lower limit
    \param    x_high   upper limit
    \param    x        a numeric
    \param    x_name   the name of the numeric

    \author Patrick Eriksson
    \date   2001-09-26
*/
void check_if_in_range( 
   const Numeric& x_low, 
   const Numeric& x_high, 
   const Numeric& x, 
   const String&  x_name ) 
{
  if ( (x<x_low) || (x>x_high) )
  {
    ostringstream os;
    os << "The variable *" << x_name <<  "* must fulfill:\n"
       << "   " << x_low << " <= " << x_name << " <= " << x_high << "\n" 
       << "The present value of *"<< x_name <<  "* is " << x << ".";
    throw runtime_error( os.str() );
  }
}



//// check_lengths (vector-vector) ///////////////////////////////////////////
//
/** Checks that two vectors have the same length

    A runtime error is thrown if the lengths of the vectors differ.

    \param    x1        vector 1
    \param    x1_name   the name of vector1
    \param    x2        vector 2
    \param    x2_name   the name of vector2

    \author Patrick Eriksson
    \date   2001-09-19
*/
void check_lengths( const Vector& x1, const String& x1_name,
                    const Vector& x2, const String& x2_name ) 
{
  if ( x1.nelem() != x2.nelem() )
  {
    ostringstream os;
    os << "The vectors *" << x1_name <<  "* and *" << x2_name << "*\n"
       << "must have the same lengths. \nThe lengths are: \n"
       << x1_name << ": " << x1.nelem() << "\n"
       << x2_name << ": " << x2.nelem();
    throw runtime_error( os.str() );
  }
}



//// check_length_nrow  /////////////////////////////////////////////////////
//
/** Checks that the length of a vector and the number of rows of a matrix
    match.

    A runtime error is thrown if the length of the vectors differs from
    the number of rows.

    \param    x        a vector
    \param    x_name   the name of the vector
    \param    A        a matrix
    \param    A_name   the name of the matrix

    \author Patrick Eriksson
    \date   2001-09-19
*/
void check_length_nrow( const Vector& x, const String& x_name,
                        const Matrix& A, const String& A_name ) 
{
  if ( x.nelem() != A.nrows() )
  {
    ostringstream os;
    os << "The length of vector *" << x_name <<  "* must be the\n"
       << "same as the number of rows of *" << A_name << "*.\n"
       << "The length of *" << x_name <<  "* is " << x.nelem() << ".\n"
       << "The number of rows of *" << A_name <<  "* is " << A.nrows() 
       << ".";
    throw runtime_error( os.str() );
  }
}



//// check_length_ncol  /////////////////////////////////////////////////////
//
/** Checkss that the length of a vector and the number of columns of a matrix
    match.

    A runtime error is thrown if the length of the vectors differs from
    the number of columns.

    \param    x        a vector
    \param    x_name   the name of the vector
    \param    A        a matrix
    \param    A_name   the name of the matrix

    \author Patrick Eriksson
    \date   2001-09-19
*/
void check_length_ncol( const Vector& x, const String& x_name,
                        const Matrix& A, const String& A_name ) 
{
  if ( x.nelem() != A.ncols() )
  {
    ostringstream os;
    os << "The length of vector *" << x_name <<  "* must be the\n"
       << "same as the number of columns of *" << A_name << "*.\n"
       << "The length of *" << x_name <<  "* is " << x.nelem() << ".\n"
       << "The number of columns of *" << A_name <<  "* is " << A.ncols()
       << ".";
    throw runtime_error( os.str() );
  }
}

//// check_ncol_nrow  /////////////////////////////////////////////////////
//
/** Checks that the number of columns of the first matrix is the same as the
    number of rows of the second matrix.

    A runtime error is thrown otherwise.

    \param    A        first Matrix
    \param    A_name   the name of the first Matrix
    \param    B        second matrix
    \param    B_name   the name of the second matrix

    \author Stefan Buehler
    \date   2001-10-02
*/
void check_ncol_nrow( const Matrix& A, const String& A_name,
                      const Matrix& B, const String& B_name ) 
{
  if ( A.ncols() != B.nrows() )
  {
    ostringstream os;
    os << "The number of columns of *" << A_name << "* must be the\n"
       << "same as the number of rows of *" << B_name << "*."
       << "The number of columns of *" << A_name <<  "* is " << A.ncols()
       << ".\n"
       << "The number of rows of *" << B_name <<  "* is " << B.nrows()
       << ".";
    throw runtime_error( os.str() );
  }
}
