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



//// isbool ///////////////////////////////////////////////////////////////////
//
/** Checks if an integer is either 0 or 1.

    \return       a boolean, true if x is either 0 or 1
    \param    x   an integer

    \author Patrick Eriksson
    \date   2001-04-19
*/
bool isbool( const int& x ) 
{
  if ( x==0 || x==1 )
    return true;
  else
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
		              MATRIX&       z,
			const VECTOR&       z0,
			const SYMMETRIC&    s )
{
  INDEX n = z.ncols();

  const size_t   nrows = z0.size();
        size_t   col;

  if ( nrows != s.nrows() )
    throw runtime_error("The length of the mean vector and the size of the covariance matrix do not match."); 

  // Make Cholesky decomposition of s, l'*l=s
  MATRIX   l(nrows,nrows);
  setto(l,0.0);
  chol(l,s);

  // Create matrix with gaussian data having zero mean and standard deviation 1
  MATRIX   r(nrows,n);
  rand_matrix_gaussian( r, 1 );

  // Multiply l and r to get z
  mult(l,r,z);

  // Add mean vector
  for ( col=0; col<n; col++ )
    add(z0,columns(z)[col]);
}



/////////////////////////////////////////////////////////////////////////////
//   Conversions between VECTOR and MATRIX types
/////////////////////////////////////////////////////////////////////////////

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



////////////////////////////////////////////////////////////////////////////
//   Interpolation routines
////////////////////////////////////////////////////////////////////////////

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
  interp_lin_vector( yi, x, y, make_vector(xi) );
  return yi[0];
}        



/////////////////////////////////////////////////////////////////////////////
//   Factorization of matrices
/////////////////////////////////////////////////////////////////////////////


//// chol //// //////////////////////////////////////////////////////////////

/** Choleski factorization (columnwise version). 

    Given s positive definite, the upper triangular matrix z with 
    positive diagonal elements such as c=r'*r is calculated.
    Algorithm used from Numerical Methods, Åke Björck, 1990, p46.

    \retval   r       Choleski factor of c
    \param    c       matrix to be factorized
 
    \author Carlos Jimenez
    \date   2001-02-14

    Adapted to MTL.
*/

void chol(
                MATRIX&       r, 
          const SYMMETRIC&    c )
{
  const INDEX nrows = c.nrows(), ncols = c.ncols();
  INDEX j, i, k;
  Numeric a = 0;

  assert( nrows == r.nrows());
  assert( ncols == r.ncols());

  for (j=0; j<nrows; ++j)
  {
    if (j>0)
    {
      for (i=0; i<j; ++i)
      {
        a = 0;
        if (i>0)
	{
          for (k=0; k<i; ++k)
	      a = a + r[k][i] * r[k][j];
	}  
        r[i][j] = ( c[i][j] - a ) / r[i][i];       
      }
    
      a = 0;
      for (k=0; k<j; ++k)
	 a = a + r[k][j] * r[k][j];
      
    }
    r[j][j] = sqrt( c[j][j] - a);
    
  }


  // Checking that it works, if r is positive definite it does not
  for (i=0; i<nrows; ++i)
    for (j=0; j<nrows; ++j)    
      if ( isnan(r[i][j]) & isinf(r[i][j]) )
      {
        ostringstream os;
        os << "Choleski decomposition does not work, s positive definite? \n";
        throw runtime_error(os.str());
      }

}
