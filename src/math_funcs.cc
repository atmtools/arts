/* Copyright (C) 2000 Patrick Eriksson <patrick@rss.chalmers.se>

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
    \item Element-wise application of common scalar functions
    \item Boolean functions                         
    \item Creation of common vectors                
    \item Interpolation routines                    
    \item Integration routines                      
    \item Conversion between vector and matrix types
   \end{enumerate}

   \author Patrick Eriksson
   \date 2000-09-18 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include <time.h>
#include "arts.h"
#include "math_funcs.h"


////////////////////////////////////////////////////////////////////////////
//   Basic mathematical vector and vector functions
////////////////////////////////////////////////////////////////////////////

//// sqrt //////////////////////////////////////////////////////////////////
//
/** Gives the elementwise square root of a vector.
   
   \retval   y   the square root of x
   \param    x   a vector

   \author Patrick Eriksson
   \date   2000-06-27
*/
void sqrt( VECTOR& y, const VECTOR& x ) 
{
  int n = x.size();
  y.resize(n);
  for ( int i=1; i<=n; i++ )
    y(i) = sqrt(x(i));
}

/** Gives the elementwise square root of a vector (return version).
   
   \return       the square root of x
   \param    x   a vector

   \author Patrick Eriksson
   \date   2000-06-27
*/
VECTOR sqrt( const VECTOR& x )     
{
  VECTOR y;
  sqrt( y, x );
  return y; 
}



//// exp ///////////////////////////////////////////////////////////////////
//
/** Gives the elementwise exponential of a vector.
   
   \retval   y   the exponential of x
   \param    x   a vector

   \author Patrick Eriksson
   \date   2000-06-27
*/
void exp( VECTOR& y, const VECTOR& x ) 
{
  int n = x.size();
  y.resize(n);
  for ( int i=1; i<=n; i++ )
    y(i) = exp(x(i));
}

/** Gives the elementwise exponential of a vector (return version).
   
   \return       the exponential of x
   \param    x   a vector

   \author Patrick Eriksson
   \date   2000-06-27
*/
VECTOR exp( const VECTOR& x )    
{
  VECTOR y;
  exp( y, x );
  return y; 
}

/** Gives the elementwise exponential of a matrix.
   
   \retval   y   the exponential of x
   \param    x   a matrix

   \author Patrick Eriksson
   \date   2000-06-27
*/
void exp( MATRIX& Y, const MATRIX& X )   
{
  int m = X.dim(1);
  int n = X.dim(2);
  Y.resize(m,n);
  for ( int i=1; i<=m; i++ ) {
    for ( int j=1; j<=n; j++ )
      Y(i,j) = exp(X(i,j)); }
}

/** Gives the elementwise exponential of a matrix (return version).
   
   \return       the exponential of x
   \param    x   a matrix

   \author Patrick Eriksson
   \date   2000-06-27
*/
MATRIX exp( const MATRIX& X )  
{
  MATRIX Y;
  exp( Y, X );
  return Y; 
}



//// log ////////////////////////////////////////////////////////////////////
//
/** Gives the elementwise natural logarithm of a vector.
   
   \retval   y   the natural logarithm of x
   \param    x   a vector

   \author Patrick Eriksson
   \date   2000-06-27
*/
void log( VECTOR& y, const VECTOR& x )   // vector parameter version
{
  int n = x.size();
  y.resize(n);
  for ( int i=1; i<=n; i++ )
    y(i) = log(x(i));
}

/** Gives the elementwise natural logarithm of a vector (return version).
   
   \return       the natural logarithm of x
   \param    x   a vector

   \author Patrick Eriksson
   \date   2000-06-27
*/
VECTOR log( const VECTOR& x )            // vector return version
{
  VECTOR y;
  log( y, x );
  return y; 
}

/** Gives the elementwise natural logarithm of a matrix.
   
   \retval   y   the natural logarithm of x
   \param    x   a matrix

   \author Patrick Eriksson
   \date   2000-06-27
*/
void log( MATRIX& Y, const MATRIX& X )   // matrix parameter version
{
  int m = X.dim(1);
  int n = X.dim(2);
  Y.resize(m,n);
  for ( int i=1; i<=m; i++ ) {
    for ( int j=1; j<=n; j++ )
      Y(i,j) = log(X(i,j)); }
}

/** Gives the elementwise natural logarithm of a matrix (return version).
   
   \return       the natural logarithm of x
   \param    x   a matrix

   \author Patrick Eriksson
   \date   2000-06-27
*/
MATRIX log( const MATRIX& X )            // matrix return version
{
  MATRIX Y;
  log( Y, X );
  return Y; 
}



//// mean and standard deviation ////////////////////////////////////////////
//
/** Calculates the mean of the rows of a matrix.

    \retval  m   row means
    \param   x   a matrix

    \author Patrick Eriksson 
    \date   2000-12-06
*/
void mean_row( VECTOR& m, const MATRIX& x )
{
  const size_t  nrows = x.dim(1);
  const size_t  ncols = x.dim(2);
        size_t  col,row;
        Numeric rowsum;

  m.resize(nrows);
  for ( row=1; row<=nrows; row++ )
  {
    rowsum = 0;
    for ( col=1; col<=ncols; col++ )
      rowsum += x(row,col);
    m(row) = rowsum/ncols;
  }
}

/** Calculates the standard deviation for the rows of a matrix.

    \retval  s   row standard deviations
    \param   x   a matrix
    \param   m   row means

    \author Patrick Eriksson 
    \date   2000-12-06
*/
void std_row( VECTOR& s, const MATRIX& x, const VECTOR& m  )
{
  const size_t  nrows = x.dim(1);
  const size_t  ncols = x.dim(2);
        size_t  col,row;
        Numeric rowsum;

  if ( nrows != m.size() )
    throw runtime_error("The size of the given mean profile does not match the data."); 

  s.resize(nrows);
  for ( row=1; row<=nrows; row++ )
  {
    rowsum = 0;
    for ( col=1; col<=ncols; col++ )
      rowsum += pow( x(row,col)-m(row), 2 );
    s(row) = rowsum/(ncols-1);
  }
}



//// min and max ////////////////////////////////////////////////////////////
//
/** Gives the minimum value of a vector.

    \return      the minimum value of x
    \param   x   a vector

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric min( const VECTOR& x )
{
  int n = x.size();
  Numeric y=x(1);
  for ( int i=2; i<=n; i++ )
  {
    if ( x(i) < y )
      y = x(i);
  }
  return y; 
}

/** Gives the minimum value of a matrix.

    \return      the minimum value of A
    \param   A   a matrix

    \autho  Stefan Buehler
    \date   2000-06-27
*/
Numeric min( const MATRIX& A )
{
  int r = A.dim(1);
  int c = A.dim(2);

  Numeric y=A(1,1);

  for ( int i=2; i<=r; i++ )
    for ( int s=2; s<=c; s++ )
      {
	if ( A(i,s) < y )
	  y = A(i,s);
      }
  return y; 
}

/** Gives the maximum value of a vector.

    \return      the maximum value of x
    \param   x   a vector

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric max( const VECTOR& x )
{
  int n = x.size();
  Numeric y=x(1);
  for ( int i=2; i<=n; i++ )
  {
    if ( x(i) > y )
      y = x(i);
  }
  return y; 
}

/** Gives the maximum value of a matrix.

    \return       the maximum value of A
    \param    A   a matrix

    \autho  Stefan Buehler
    \date   2000-06-27
*/
Numeric max( const MATRIX& A )
{
  int r = A.dim(1);
  int c = A.dim(2);

  Numeric y=A(1,1);

  for ( int i=2; i<=r; i++ )
    for ( int s=2; s<=c; s++ )
      {
	if ( A(i,s) > y )
	  y = A(i,s);
      }
  return y; 
}



//// first and last /////////////////////////////////////////////////////////
//
/** Gives the first value of a vector.

    \return       the first value of x
    \param    x   a vector

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric first( const VECTOR& x )
{
  return x(1); 
}

/** Gives the last value of a vector.

    \return      the last value of x
    \param   x   a vector

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric last( const VECTOR& x )
{
  return x(x.size()); 
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
bool any( const ARRAYofsizet& x ) 
{
  for ( size_t i=0; i<x.size(); i++ ) {
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

    \retval   x       linearly spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    step    distance between values in x

    \author Patrick Eriksson
    \date   2000-06-27
*/
void linspace(                      
              VECTOR&     x,           
        const Numeric     start,    
        const Numeric     stop,        
        const Numeric     step )
{
  int n = (int) floor( (stop-start)/step ) + 1;
  if ( n<1 )
    n=1;
  x.resize(n);
  for ( int i=1; i<=n; i++ )
    x(i) = start + (i-1)*step;
}

/** Linearly spaced vector with specified spacing (return version). 

    The first element of the returned vector is always start. 
    The next value is start+step etc.
    Note that the last value can deviate from stop.
    The step can be both positive and negative. 
    (in Matlab notation: start:step:stop)

    \return           linearly spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    step    distance between values in x

    \author Patrick Eriksson
    \date   2000-06-27
*/
VECTOR linspace(  
        const Numeric start, 
        const Numeric stop,  
        const Numeric step )
{
  VECTOR x;
  linspace( x, start, stop, step );
  return x; 
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
              VECTOR&     x, 
        const Numeric     start,     
        const Numeric     stop,        
        const int         n )
{
  if ( n<2 )
    throw runtime_error("NLINSPACE: The number of points must be > 1"); 
  x.resize(n);
  Numeric step = (stop-start)/(n-1) ;
  for ( int i=1; i<=n; i++ )
    x(i) = start + (i-1)*step;
}

/** Linearly spaced vector with specified length (return version). 

    Returns a vector equally and linearly spaced between start and stop 
    of length n. (equals the Matlab function linspace)

    The length must be > 1.

    \return   x       linearly spaced vector
    \param    start   first value in x
    \param    stop    last value of x <= stop
    \param    n       length of x

    \author Patrick Eriksson
    \date   2000-06-27
*/
VECTOR nlinspace(                 // As above but return version
        const Numeric start, 
        const Numeric stop,  
        const int     n )
{
  VECTOR x;
  nlinspace( x, start, stop, n );
  return x; 
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
              VECTOR&     x, 
        const Numeric     start,     
        const Numeric     stop,        
        const int         n )
{
  if ( n<2 )
    throw runtime_error("NLOGSPACE: The number of points must be > 1"); 
  if ( (start<=0) || (stop<=0) )
    throw runtime_error("NLOGSPACE: Only positive numbers are allowed"); 
  x.resize(n);
  Numeric a = log(start);
  Numeric step = (log(stop)-a)/(n-1) ;
  x(1) = start;
  for ( int i=2; i<n; i++ )
    x(i) = exp(a + (i-1)*step);
  x(n) = stop;
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
VECTOR nlogspace(  
        const Numeric start, 
        const Numeric stop,  
        const int     n )
{
  VECTOR x;
  nlogspace( x, start, stop, n );
  return x; 
}                     



////////////////////////////////////////////////////////////////////////////
//   Interpolation routines
//
//     All the functions assume that the interpolation points, XI, are ordered.
//     However, the values can be in increasing or decreasing order.
//
////////////////////////////////////////////////////////////////////////////

// Local help function to check input grids
//
// Patrick Eriksson 2000-06-29
//
int interp_check(                  
	const VECTOR&  x,        
        const VECTOR&  xi,
        const size_t   n_y )
{
  size_t  n  = x.size();
  size_t  ni = xi.size();
  int     order=1;          // flag for pos. or neg. order, assume pos.

  // Determine the order, -1=decreasing and 1=increasing
  if ( x(1) > x(n) )
    order = -1;

  if ( n < 2 )
    throw runtime_error("Vector length for interpolation must be >= 2");

  if ( n != n_y ) 
    throw runtime_error("Sizes of input data to interpolation do not match");

  if ( (order*xi(1)<order*x(1)) || (order*xi(ni)>order*x(n)) ) 
    {
      ostringstream os;
      os << "Interpolation points must be inside the original range.\n"
	 << "Int.:  xi(1) = " << xi(1) << ", xi(ni) = " << xi(ni) << '\n'
	 << "Orig.: x(1)  = " << x(1)  << ", x(n)   = " << x(n);
      throw runtime_error(os.str());
    }

  for (size_t i=1; i<n; i++ )
  {
    if ( order*x(i+1) < order*x(i) ) 
      throw runtime_error("Original interpolation grid must be ordered");
  }

  for (size_t i=1; i<ni; i++ )
  {
    if ( order*xi(i+1) < order*xi(i) ) 
      throw runtime_error("Interpolation points must be ordered");
  }

  return order;
}



//// interp_lin ////////////////////////////////////////////////////////////
//
/** Multiple linear interpolation of a vector.

    The vector x specifies the points at which the data y is given. 

    The size of yi is the same as for xi.

    \retval  yi      interpolated values 
    \param   x       the x grid
    \param   y       the function to interpolate
    \param   xi      interpolation points

    \author Patrick Eriksson
    \date   2000-06-29
*/
void interp_lin(          
              VECTOR&  yi,
        const VECTOR&  x, 
        const VECTOR&  y, 
        const VECTOR&  xi )
{
  // Check grids and get order of grids
  int order = interp_check( x, xi, y.size() ); 

  size_t        i, j=1, n=xi.size();
  Numeric       w;
  yi.resize(n); 

  for ( i=1; i<=n; i++ )
  {
    for( ;  order*x(j+1) < order*xi(i); j++ ) {}
    w = (xi(i)-x(j)) / (x(j+1)-x(j));
    yi(i) = y(j) + w * (y(j+1)-y(j)); 
  }
}      

/** Multiple linear interpolation of a vector (return version).

    The vector x specifies the points at which the data y is given. 

    The size of yi is the same as for xi.

    \return          interpolated values 
    \param   x       the x grid
    \param   y       the function to interpolate
    \param   xi      interpolation points

    \author Patrick Eriksson
    \date   2000-06-29
*/
VECTOR interp_lin( 
        const VECTOR&  x, 
        const VECTOR&  y, 
        const VECTOR&  xi )
{
  VECTOR yi; 
  interp_lin( yi, x, y, xi );
  return yi;
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
        const VECTOR&  x, 
        const VECTOR&  y, 
        const Numeric  xi )
{
  VECTOR yi; 
  interp_lin( yi, x, y, VECTOR(1,xi) );
  return yi(1);
}        

/** Multiple linear interpolation of matrix rows.

    The vector x specifies the points at which the data y is given. 

    \retval  Yi      interpolated values 
    \param   x       the x grid
    \param   Y       the function to interpolate
    \param   xi      interpolation points

    \author Patrick Eriksson
    \date   2000-06-29
*/
void interp_lin_row(    
              MATRIX&  Yi,
        const VECTOR&  x, 
        const MATRIX&  Y, 
        const VECTOR&  xi )
{
  // Check grids and get order of grids
  int order = interp_check( x, xi, Y.dim(2) ); 

  int        k, j=1, n=xi.size(), nrow=Y.dim(1);
  Numeric    w;
  Yi.resize( nrow, n ); 

  for (int i=1; i<=n; i++ )
  {
    for( ;  order*x(j+1) < order*xi(i); j++ ) {}
    w = (xi(i)-x(j)) / (x(j+1)-x(j));
    for( k=1; k<=nrow; k++ )
      Yi(k,i) = Y(k,j) + w * (Y(k,j+1)-Y(k,j)); 
  }
}        

/** Multiple linear interpolation of matrix rows (return version).

    The vector x specifies the points at which the data y is given. 

    \retval  Yi      interpolated values 
    \param   x       the x grid
    \param   Y       the function to interpolate
    \param   xi      interpolation points

    \author Patrick Eriksson
    \date   2000-06-29
*/
MATRIX interp_lin_row(       // As above but return version
        const VECTOR&  x, 
        const MATRIX&  Y, 
        const VECTOR&  xi )
{
  MATRIX Yi; 
  interp_lin_row( Yi, x, Y, xi );
  return Yi;
}        

/** Multiple linear interpolation of matrix columns.

    The vector x specifies the points at which the data y is given. 

    \retval  Yi      interpolated values 
    \param   x       the x grid
    \param   Y       the function to interpolate
    \param   xi      interpolation points

    \author Stefan Buehler
    \date   2000-06-29
*/
void interp_lin_col(    
              MATRIX&  Yi,
        const VECTOR&  x, 
        const MATRIX&  Y, 
        const VECTOR&  xi )
{
  // Check grids and get order of grids
  int order = interp_check( x, xi, Y.dim(1) ); 

  int        k, j=1, n=xi.size(), ncol=Y.dim(2);
  Numeric    w;
  Yi.resize( n, ncol ); 

  for (int i=1; i<=n; i++ )
  {
    for( ;  order*x(j+1) < order*xi(i); j++ ) {}
    w = (xi(i)-x(j)) / (x(j+1)-x(j));
    for( k=1; k<=ncol; k++ )
      Yi(i,k) = Y(j,k) + w * (Y(j+1,k)-Y(j,k)); 
  }
}        

/** Multiple linear interpolation of matrix columns (return version).

    The vector x specifies the points at which the data y is given. 

    \return          interpolated values 
    \param   x       the x grid
    \param   Y       the function to interpolate
    \param   xi      interpolation points

    \author Stefan Buehler
    \date   2000-06-29
*/
MATRIX interp_lin_col(       // As above but return version
        const VECTOR&  x, 
        const MATRIX&  Y, 
        const VECTOR&  xi )
{
  MATRIX Yi; 
  interp_lin_col( Yi, x, Y, xi );
  return Yi;
}        



////////////////////////////////////////////////////////////////////////////
//   Intergration routines
//
//     For the moment, not used.
//
////////////////////////////////////////////////////////////////////////////
/*
//
// Patrick Eriksson 2000-06-29
// 
Numeric integr_lin(                // Integrates Y over X assuming that Y is
        const VECTOR&  x,          // linear between the given points
        const VECTOR&  y )
{
  size_t i, n = x.size();
  Numeric w=0.0; 

  if ( n < 2 )
    throw runtime_error("INTEGR_LIN: Vector length must be >= 2");
  if ( n != y.size() )
    throw runtime_error("INTEGR_LIN: Sizes of input data do not match");

  for( i=1; i<n; i++ )
    w += (x(i+1)-x(i)) * (y(i)+y(i+1))/2.0;
  return w;
}

void integr_lin(                    // As above but parameter version
              Numeric&  w,
        const VECTOR&   x,  
        const VECTOR&   y ) 
{
  w = integr_lin( x, y );
}

void integr_lin(                // Integrates the rows of M assuming that
              MATRIX&  W,       // the functions are linear between the
        const VECTOR&  x,       // given points
        const MATRIX&  M )   
{
  size_t i,j,rows= M.dim(1), cols= M.dim(2);
  Numeric w; 
  W.resize(rows,1);
  W = 0.0;

  if ( cols < 2 )
    throw runtime_error("Vector length for integration must be >= 2");
  if ( cols != x.size() )
    throw runtime_error("Sizes of input data for integration do not match");

  for ( i=1; i<cols; i++ ) 
  {
    w = ( x(i+1) - x(i) ) / 2.0;
    for ( j=1; j<=rows; j++ )
      W(j,1) += w * ( M(j,i) + M(j,i+1) );
  }
}

MATRIX integr_lin(              // As above but return version
        const VECTOR&  x,  
        const MATRIX&  M ) 
{
  MATRIX W;
  integr_lin( W, x, M );
  return W;
}
*/



/////////////////////////////////////////////////////////////////////////////
//   Conversions between VECTOR and MATRIX types
/////////////////////////////////////////////////////////////////////////////

//// to_matrix //////////////////////////////////////////////////////////////
//
/** Converts a vector to a matrix. 

    For a vector of length n, the dimension of the matrix is [n,1], in other 
    words, the vector is interpreted as a column vector.

    \retval   W       the matrix, size n x 1
    \param    x       a vector of length n

    \author Stefan Buehler 
    \date 2000-09-01
*/
void to_matrix(MATRIX& W, const VECTOR& x)
{
  // FIXME: I'm sure this can be made more efficient when TNT has more
  // functionality.  
  W.resize(x.size(),1);
  for (size_t i=1; i<=x.size() ; ++i)
    {
      W(i,1) = x(i);
    }
}

/** Converts a vector to a matrix (return version). 

    For a vector of length n, the dimension of the matrix is [n,1], in other 
    words, the vector is interpreted as a column vector.

    \return           the matrix, size n x 1
    \param    x       a vector of length n

    \author Stefan Buehler 
    \date 2000-09-01
*/
MATRIX to_matrix(const VECTOR& x)
{
  MATRIX W;
  to_matrix(W,x);
  return W;
}



//// to_vector //////////////////////////////////////////////////////////////
//
/** Conversion of a matrix to a vector.

    The matrix can either be a column (n x 1) or row (1 x n) vector.

    Uses the functions row and col.

    \retval  x       the vector, length n
    \param   W       matrix, size (n x 1) or (1 x n)

    @exception runtime_error None of the dimensions of W was 1.

    @see row col

    \author Stefan Buehler 
    \date 2000-09-01
*/
void to_vector(VECTOR& x, const MATRIX& W)
{
  // FIXME: I'm sure this can be made more efficient when TNT has more
  // functionality.

  // Check if one of the dimensions of W is 1:
  if ( 1 == W.dim(2) )
    {
      col(x,1,W);
    }
  else if ( 1 == W.dim(1) )
    {
      row(x,1,W);
    }
  else
    throw runtime_error("You tried to convert a matrix to a vector,\n"
			"but none of the dimensions is 1.");
}

/** Conversion of a matrix to a vector (return version).

    The matrix can either be a column (n x 1) or row (1 x n) vector.

    Uses the functions row and col.

    \return  x       the vector, length n
    \param   W       matrix, size (n x 1) or (1 x n)

    @exception runtime_error None of the dimensions of W was 1.

    @see row col

    \author Stefan Buehler 
    \date 2000-09-01
*/
VECTOR to_vector(const MATRIX& W)
{
  VECTOR x;
  to_vector(x,W);
  return x;
}



/////////////////////////////////////////////////////////////////////////////
//   Extraction of matrix columns and rows
/////////////////////////////////////////////////////////////////////////////

//// row ////////////////////////////////////////////////////////////////////
//
/** Extracts row i of MATRIX A.

    \retval   x   row i of A as a vector
    \param    i   row index
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
void row(VECTOR& x,
	 size_t i,
	 const MATRIX& A)
{
  // Make sure that i is legal:
  assert ( i >  0        );
  assert ( i <= A.dim(1) );

  const size_t n = A.dim(2);
  x.resize(n);
  for ( size_t j=1; j<=n; ++j )
    {
      x(j) = A(i,j);
    }
}

/** Extracts row i of MATRIX A (return version).

    \return       row i of A as a vector
    \param    i   row index
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
VECTOR row(size_t i,
	   const MATRIX& A)
{
  VECTOR x;
  row(x,i,A);
  return x;
}

/** Extracts row i of MATRIX A.

    \retval   x   column i of A as a vector
    \param    i   column index
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
void row(MATRIX& X,
	 size_t i,
	 size_t k,
	 const MATRIX& A)
{
  // Make sure that i and k are legal:
  assert ( i >  0        );
  assert ( i <= A.dim(1) );
  assert ( k >  0        );
  assert ( k <= A.dim(1) );
  assert ( k >= i        );

  const size_t n = A.dim(2);
  const size_t m = k-i + 1;
  X.resize(m,n);
  for ( size_t j=1; j<=n; ++j )
    for ( size_t l=1; l<=m; ++l )
      {
	X(l,j) = A(i+l-1,j);
      }
}

/** Extracts row i of MATRIX A (return version).

    \return       column i of A as a vector
    \param    i   column index
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
MATRIX row(size_t i,
	    size_t k,
	    const MATRIX& A)
{
  MATRIX X;
  row(X,i,k,A);
  return X;
}



//// column /////////////////////////////////////////////////////////////////
//
/** Extracts column i of MATRIX A.

    \retval   X   the extracted part of A
    \param    i   first row to extract
    \param    k   last row to extract
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
void col(VECTOR& x,
	 size_t i,
	 const MATRIX& A)
{
  // Make sure that i is legal:
  assert ( i >  0        );
  assert ( i <= A.dim(2) );

  const size_t n = A.dim(1);
  x.resize(n);
  for ( size_t j=1; j<=n; ++j )
    {
      x(j) = A(j,i);
    }
}

/** Extracts column i of MATRIX A (return version).

    \return       the extracted part of A
    \param    i   first row to extract
    \param    k   last row to extract
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
VECTOR col(size_t i,
	   const MATRIX& A)
{
  VECTOR x;
  col(x,i,A);
  return x;
}

/** Extracts coulmns i to k of MATRIX A.

    \retval   X   the extracted part of A
    \param    i   first column to extract
    \param    k   last column to extract
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
void col(MATRIX& X,
	  size_t i,
	  size_t k,
	  const MATRIX& A)
{
  // Make sure that i and k are legal:
  assert ( i >  0        );
  assert ( i <= A.dim(2) );
  assert ( k >  0        );
  assert ( k <= A.dim(2) );
  assert ( k >= i        );

  const size_t n = A.dim(1);
  const size_t m = k-i + 1;
  X.resize(n,m);
  for ( size_t j=1; j<=n; ++j )
    for ( size_t l=1; l<=m; ++l )
      {
	X(j,l) = A(j,i+l-1);
      }
}

/** Extracts coulmns i to k of MATRIX A (return version).

    \return       the extracted part of A
    \param    i   first column to extract
    \param    k   last column to extract
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
MATRIX col(size_t i,
	    size_t k,
	    const MATRIX& A)
{
  MATRIX X;
  col(X,i,k,A);
  return X;
}



/////////////////////////////////////////////////////////////////////////////
//   Putting data in a matrix column or and row
/////////////////////////////////////////////////////////////////////////////

//// put_in_col /////////////////////////////////////////////////////////////
//
/** Puts a vector in column i of a matrix

    \retval   A   a matrix
    \param    i   column index
    \param    x   a vector

    \author Patrick Eriksson 
    \date   2000-12-07
*/
void put_in_col(
              MATRIX& A,
	      size_t  i, 
        const VECTOR& x )
{
  const size_t n = A.dim(1);

  // Make sure that i is legal:
  assert ( i >  0        );
  assert ( i <= A.dim(2) );

  // Check length of x
  assert ( n == x.size()  );

  for ( size_t j=1; j<=n; ++j )
    A(j,i) =  x(j);
}



/////////////////////////////////////////////////////////////////////////////
//   Random data
/////////////////////////////////////////////////////////////////////////////

//// rand_uniform ///////////////////////////////////////////////////////////
/**
   Creates a vector with random data uniformerly distributed between
   the lower and higher limit given.

   The random data is uncorrelated.

   \retval   y          random vector
   \param    n          length of the random vector to generate
   \param    x_low      lower limit for the random values
   \param    x_high     upper limit for the random data

   \author Patrick Eriksson
   \date   2000-11-26
*/
void rand_uniform(
              VECTOR&    r,
        const size_t     n,
        const Numeric&   x_low,
        const Numeric&   x_high )
{
  const Numeric dx = x_high-x_low;

  // Init seed in a "random" way
  srand( (unsigned int) time( NULL );

  r.resize(n);
  for ( size_t i=1; i<=n; i ++)
    r(i) = x_low + dx* (double(rand())/double(RAND_MAX));
}



//// rand_gaussian ///////////////////////////////////////////////////////////
/**
   Creates a gaussian random vector with zero mean and 
   the standard deviation given.

   The random data is uncorrelated.

   The algorith is taken from Numerical Recipies, Section 7.2. 
   See www.nr.com.

   \retval   y          random vector
   \param    n          length of the random vector to generate
   \param    s          standard deviation

   \author Patrick Eriksson
   \date   2000-11-27
*/
void rand_gaussian(
              VECTOR&    r,
        const size_t     n,
        const Numeric&   s )
{
  VECTOR  z;    // A vector of length 2 with uniform PDF between -1 and 1
  Numeric rad;  // The radius cooresponding to z
  Numeric fac;  // Normalisation factor
 
  r.resize(n);
  for ( size_t i=0; i<n; )
  {
    rand_uniform( z, 2, -1, 1 );
    rad = z(1)*z(1) + z(2)*z(2);
    if ( (rad<1) && (rad>0) )
    {
      fac = sqrt( -2*log(rad)/rad );
      i++;
      r(i) = s*fac*z(1);
      i++;
      if ( i <=n )
        r(i) = s*fac*z(2);        
    }
  }
}



//// rand_matrix_uniform ////////////////////////////////////////////////////
/**
   Creates a matrix with random data uniformerly distributed between
   the lower and higher limit given.

   The random data is uncorrelated.

   \retval   m          random matrix
   \param    nrows      number of rows of m
   \param    ncols      number of cols of m
   \param    x_low      lower limit for the random values
   \param    x_high     upper limit for the random data

   \author Patrick Eriksson
   \date   2000-12-07
*/
void rand_matrix_uniform(
              MATRIX&    m,
        const size_t&    nrows,
        const size_t&    ncols,
        const Numeric&   x_low,
        const Numeric&   x_high )
{
  VECTOR r;
  size_t row,col;
  m.resize(nrows,ncols);
  for ( col=1; col<=ncols; col++)
  {
    rand_uniform( r, nrows, x_low, x_high );
    for ( row=1; row<=nrows; row++)
      m(row,col) = r(row);
  }
}



//// rand_matrix_gaussian ////////////////////////////////////////////////////
/**
   Creates a gaussian random matrix with zero mean and 
   the standard deviation given.

   The random data is uncorrelated.

   See further rand_gaussian

   \retval   y          random vector
   \param    nrows      number of rows of m
   \param    ncols      number of cols of m
   \param    s          standard deviation

   \author Patrick Eriksson
   \date   2000-12-07
*/
void rand_matrix_gaussian(
              MATRIX&    m,
        const size_t&    nrows,
        const size_t&    ncols,
        const Numeric&   s )
{
  VECTOR r;
  size_t row,col;
  m.resize(nrows,ncols);
  for ( col=1; col<=ncols; col++)
  {
    rand_gaussian( r, nrows, s );
    for ( row=1; row<=nrows; row++)
      m(row,col) = r(row);
  }
}



//// rand_data_gaussian ////////////////////////////////////////////////////
/**
   Creates a matrix with random data fulfilling the given criteria.

   The mean of the data is z0. Standard deviations and correlation 
   follow the given covariance matrix.

   Each realisation is stored as a column in z.

   \retval   z          random matrix
   \param    n          number of random vectors
   \param    z0         mean vector
   \param    s          covariance matrix

   \author Patrick Eriksson
   \date   2000-12-07
*/
void rand_data_gaussian(
              MATRIX&    z,
        const size_t&    n,
        const VECTOR&    z0,
        const MATRIX&    s )
{
  const size_t   nrows = z0.size();
        size_t   row,col;

  if ( nrows != s.dim(1) )
    throw runtime_error("The length of the mean vector and the size of the covariance matrix do not match."); 

  // Make Cholesky decomposition of s, l'*l=s
  // Until we have MTL, only the diagonal elements are considered.
  MATRIX   l(nrows,nrows,0.0);
  for ( row=1; row<=nrows; row++ )
    l(row,row) = sqrt( s(row,row) );

  // Create matrix with gaussian data having zero mean and standard deviation 1
  MATRIX   r;
  rand_matrix_gaussian( r, nrows, n, 1 );

  // Multiplicate l and r to get z
  z = l*r;

  // Add mean vector
  for ( row=1; row<=nrows; row++ )
  {
    for ( col=1; col<=n; col++ )
      z(row,col) += z0(row);
  }  
}



