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
   \end{enumerate}

   \author Patrick Eriksson
   \date 2000-09-18 

   Adapted to MTL. Gone from 1-based to 0-based. Made non-return
   versions do no resize. Used more generic algorithms in some cases.
   \date 2000-12-25
   \author Stefan Buehler
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include <time.h>
#include "arts.h"
#include "math_funcs.h"
#include "make_vector.h"


//// mean and standard deviation ////////////////////////////////////////////
//
/** Calculates the mean of the rows of a matrix.

    \retval  m   row means
    \param   x   a matrix
    
    \author Patrick Eriksson 
    \date   2000-12-06
    
    Made implementation `generic' using MTL sum algorithm.
    \author Stefan Buehler
    \date   2000-12-24
*/
void mean_row( VECTOR& m, const MATRIX& x )
{
  m = VECTOR(x.nrows());

  for ( INDEX i=0; i<x.nrows(); ++i ) 
    {
      m[i] = mtl::sum(x[i])/x.ncols();
    }
}

/** Calculates the standard deviation for the rows of a matrix.

    \retval  s   row standard deviations
    \param   x   a matrix
    \param   m   row means

    \author Patrick Eriksson 
    \date   2000-12-06

    Made implementation `generic' using MTL sum algorithm.
    \author Stefan Buehler
    \date   2000-12-24
*/
void std_row( VECTOR& s, const MATRIX& x, const VECTOR& m  )
{
  VECTOR d(x.ncols());		// We need this to store the deviation
				// from the mean.

  if ( x.nrows() != m.size() )
    throw runtime_error("std_row: The size of the given mean profile "
			"does not match the data.");
  
  s = VECTOR(x.nrows());
  
  for ( INDEX i=0; i<x.nrows(); ++i ) 
    {
      // Store -mean in d:
      setto(d,-m[i]);
      add(x[i],d);
      ele_mult(d,d,d);
      s[i] = mtl::sum(d) / (x.ncols()-1);
    }  
}



//// min and max ////////////////////////////////////////////////////////////
//

// Min of a vector is already part of MTL.

/** Gives the minimum value of a matrix.

    \return      the minimum value of A
    \param   A   a matrix

    \author  Stefan Buehler
    \date   2000-06-27
*/
// FIXME: Maybe we do need this for matrix.
// Numeric min( const MATRIX& A )
// {
//   int r = A.nrows();
//   int c = A.ncols();

//   Numeric y=A(1,1);

//   for ( int i=2; i<=r; i++ )
//     for ( int s=2; s<=c; s++ )
//       {
// 	if ( A(i,s) < y )
// 	  y = A(i,s);
//       }
//   return y; 
// }

// Max of a vector is part of MTL.


/** Gives the maximum value of a matrix.

    \return       the maximum value of A
    \param    A   a matrix

    \autho  Stefan Buehler
    \date   2000-06-27
*/
// FIXME: Maybe we do need this for matrix.
// Numeric max( const MATRIX& A )
// {
//   int r = A.nrows();
//   int c = A.ncols();

//   Numeric y=A(1,1);

//   for ( int i=2; i<=r; i++ )
//     for ( int s=2; s<=c; s++ )
//       {
// 	if ( A(i,s) > y )
// 	  y = A(i,s);
//       }
//   return y; 
// }



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
  return x[0]; 
}

/** Gives the last value of a vector.

    \return      the last value of x
    \param   x   a vector

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric last( const VECTOR& x )
{
  return x[x.size()-1]; 
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

   Adapted to MTL. Gone from 1-based to 0-based. 
   \date 2000-12-25
   \author Stefan Buehler
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

    Adapted to MTL. Gone from 1-based to 0-based. 
    \date 2000-12-25
    \author Stefan Buehler
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
  x = VECTOR(n);
  for ( int i=0; i<n; i++ )
    x[i] = start + i*step;
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

    Adapted to MTL. Gone from 1-based to 0-based. 
    \date 2000-12-25
    \author Stefan Buehler
*/
void nlinspace(         
	       VECTOR&     x, 
	       const Numeric     start,     
	       const Numeric     stop,        
	       const int         n )
{
  if ( n<2 )
    throw runtime_error("nlinspace: The number of points must be > 1"); 
  x = VECTOR(n);
  Numeric step = (stop-start)/(n-1) ;
  for ( int i=0; i<n; i++ )
    x[i] = start + i*step;
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

    Adapted to MTL. Gone from 1-based to 0-based. 
    \date 2000-12-25
    \author Stefan Buehler
*/
void nlogspace(         
	       VECTOR&     x, 
	       const Numeric     start,     
	       const Numeric     stop,        
	       const int         n )
{
  if ( n<2 )
    throw runtime_error("nlogspace: The number of points must be > 1"); 
  if ( (start<=0) || (stop<=0) )
    throw runtime_error("nlogspace: Only positive numbers are allowed"); 
  x = VECTOR(n);
  Numeric a = log(start);
  Numeric step = (log(stop)-a)/(n-1) ;
  x[0] = start;
  for ( int i=1; i<n-1; i++ )
    x[i] = exp(a + i*step);
  x[n-1] = stop;
  //  cout << "x: " << x << "\n";
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




//// interp_lin ////////////////////////////////////////////////////////////
//
/** Multiple linear interpolation of a vector.

    The vector x specifies the points at which the data y is given. 

    The size of yi has to be the same as for xi.

    \retval  yi      interpolated values 
    \param   x       the x grid
    \param   y       the function to interpolate
    \param   xi      interpolation points

    \author Patrick Eriksson
    \date   2000-06-29

    Adapted to MTL. Gone from 1-based to 0-based. No resize!
    \date 2000-12-25
    \author Stefan Buehler
*/
void interp_lin(          
		VECTOR&  yi,
		const VECTOR&  x, 
		const VECTOR&  y, 
		const VECTOR&  xi )
{
  // Check grids and get order of grids
  int order = interp_check( x, xi, y.size() ); 

  size_t        i, j=0, n=xi.size();
  Numeric       w;

  assert( n==yi.size() ); 

  for ( i=0; i<n; i++ )
  {
    for( ;  order*x[j+1] < order*xi[i]; j++ ) {}
    w = (xi[i]-x[j]) / (x[j+1]-x[j]);
    yi[i] = y[j] + w * (y[j+1]-y[j]); 
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
  VECTOR yi(xi.size()); 
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
  VECTOR yi(1); 
  interp_lin( yi, x, y, make_vector(xi) );
  return yi[0];
}        

/** Multiple linear interpolation of matrix rows.

    The vector x specifies the points at which the data y is given. 

    \retval  Yi      interpolated values 
    \param   x       the x grid
    \param   Y       the function to interpolate
    \param   xi      interpolation points

    \author Patrick Eriksson
    \date   2000-06-29

    Adapted to MTL. Gone from 1-based to 0-based. No resize!
    \date 2000-12-25
    \author Stefan Buehler
*/
void interp_lin_row(    
		    MATRIX&  Yi,
		    const VECTOR&  x, 
		    const MATRIX&  Y, 
		    const VECTOR&  xi )
{
  // Check grids and get order of grids
  int order = interp_check( x, xi, Y.ncols() ); 

  INDEX      j=0, n=xi.size(), nrow=Y.nrows();
  Numeric    w;

  assert( nrow == Yi.nrows() );
  assert( n    == Yi.ncols() );

  for (INDEX i=0; i<n; i++ )
  {
    for( ;  order*x[j+1] < order*xi[i]; j++ ) {}
    w = (xi[i]-x[j]) / (x[j+1]-x[j]);
    for( INDEX k=0; k<nrow; k++ )
      Yi[k][i] = Y[k][j] + w * (Y[k][j+1]-Y[k][j]); 
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
  MATRIX Yi( Y.nrows(), xi.size() ); 
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

    Adapted to MTL. Gone from 1-based to 0-based. No resize!
    \date 2000-12-25
    \author Stefan Buehler
*/
void interp_lin_col(    
		    MATRIX&  Yi,
		    const VECTOR&  x, 
		    const MATRIX&  Y, 
		    const VECTOR&  xi )
{
  // Check grids and get order of grids
  int order = interp_check( x, xi, Y.nrows() ); 

  INDEX        j=0, n=xi.size(), ncol=Y.ncols();
  Numeric    w;

  assert( n    == Yi.nrows() );
  assert( ncol == Yi.ncols() );

  for (INDEX i=0; i<n; i++ )
  {
    for( ;  order*x[j+1] < order*xi[i]; j++ ) {}
    w = (xi[i]-x[j]) / (x[j+1]-x[j]);
    for( INDEX k=0; k<ncol; k++ )
      Yi[i][k] = Y[j][k] + w * (Y[j+1][k]-Y[j][k]); 
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
  MATRIX Yi( xi.size(), Y.ncols() ); 
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
  size_t i,j,rows= M.nrows(), cols= M.ncols();
  Numeric w; 
  resize(W,rows,1);
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
//   Random data
/////////////////////////////////////////////////////////////////////////////

//// rand_uniform ///////////////////////////////////////////////////////////
// Is now a template function, hence has to reside in math_funcs.h

//// rand_gaussian ///////////////////////////////////////////////////////////
// Now a template function and hence in math_funcs.h

//// rand_matrix_uniform ////////////////////////////////////////////////////
/**
   Creates a matrix with random data uniformerly distributed between
   the lower and higher limit given.

   The random data is uncorrelated.

   \retval   m          random matrix
   \param    x_low      lower limit for the random values
   \param    x_high     upper limit for the random data

   \author Patrick Eriksson
   \date   2000-12-07
*/
void rand_matrix_uniform(
		  MATRIX&    m,
		  const Numeric&   x_low,
		  const Numeric&   x_high )
{
  for ( INDEX i=0; i<m.nrows(); ++i )
  {
    rand_uniform( m[i], x_low, x_high );
  }
}



//// rand_matrix_gaussian ////////////////////////////////////////////////////
/**
   Creates a gaussian random matrix with zero mean and 
   the standard deviation given.

   The random data is uncorrelated.

   See further rand_gaussian

   \retval   y          random vector
   \param    s          standard deviation

   \author Patrick Eriksson
   \date   2000-12-07
*/
void rand_matrix_gaussian(
		   MATRIX&    m,
		   const Numeric&   s )
{
  for ( INDEX i=0; i<m.nrows(); ++i )
  {
    rand_gaussian( m[i], s );
  }
}



//// rand_data_gaussian ////////////////////////////////////////////////////
/**
   Creates a matrix with random data fulfilling the given criteria.

   The mean of the data is z0. Standard deviations and correlation 
   follow the given covariance matrix.

   Each realisation is stored as a column in z. The number of columns
   in z hence determines the number of realizations generated.

   \retval   z          random matrix
   \param    z0         mean vector
   \param    s          covariance matrix

   \author Patrick Eriksson
   \date   2000-12-07

   Adapted to MTL. Gone from 1-based to 0-based. No resize for
   non-return versions.
   \date 2000-12-25
   \author Stefan Buehler
*/
void rand_data_gaussian(
			MATRIX&    z,
			const VECTOR&    z0,
			const MATRIX&    s )
{
  INDEX n = z.ncols();

  const size_t   nrows = z0.size();
        size_t   row,col;

  if ( nrows != s.nrows() )
    throw runtime_error("The length of the mean vector and the size of the covariance matrix do not match."); 

  // Make Cholesky decomposition of s, l'*l=s
  // Until we have MTL, only the diagonal elements are considered.
  MATRIX   l(nrows,nrows);
  setto(l,0.0);
  for ( row=0; row<nrows; row++ )
    l[row][row] = sqrt( s[row][row] );

  // Create matrix with gaussian data having zero mean and standard deviation 1
  MATRIX   r(nrows,n);
  rand_matrix_gaussian( r, 1 );

  // Multiply l and r to get z
  mult(l,r,z);

  // Add mean vector
  for ( col=0; col<n; col++ )
    add(z0,columns(z)[col]);
}


//// to_vector //////////////////////////////////////////////////////////////
//
/** Conversion of a matrix to a vector.

    The matrix can either be a column (n x 1) or row (1 x n) vector.

    Uses MTL functionality. Adapted from the old to_vector.

    \retval  x       the vector, length n
    \param   W       matrix, size (n x 1) or (1 x n)

    @exception runtime_error None of the dimensions of W was 1.

    \author Stefan Buehler 
    \date   2001-01-06
*/
void to_vector(VECTOR& x, const MATRIX& W)
{
  // Check if one of the dimensions of W is 1:
  if ( 1 == W.ncols() )
    {
      //      col(x,1,W);
      resize( x, W.nrows() );
      copy( columns(W)[0], x );
    }
  else if ( 1 == W.nrows() )
    {
      //      row(x,1,W);
      resize( x, W.ncols() );
      copy( W[0], x );
    }
  else
    throw runtime_error("You tried to convert a matrix to a vector,\n"
			"but none of the dimensions is 1.");
}
