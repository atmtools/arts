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
    \item Element-wise application of common scalar functions
    \item Boolean functions                         
    \item Creation of common vectors                
    \item Interpolation routines                    
    \item Integration routines                      
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


//// mean and standard deviation ////////////////////////////////////////////
//
/** Calculates the mean of the rows of a matrix.

    Dimensions of m and x must match!

    \retval  m   row means
    \param   x   a matrix
    
    \author Patrick Eriksson 
    \date   2000-12-06
*/
void mean_row( VectorView m, ConstMatrixView x )
{
  assert( m.nelem()==x.nrows() );

  for ( Index i=0; i<x.nrows(); ++i ) 
    {
      m[i] = x(i,Range(joker)).sum() / x.ncols();
      // x(i,Range(joker)) picks out the ith row of x. The member
      // function .sum computes the sum of all elements. This we just
      // have to divide by the number of columns to get the mean.
    }
}

/** Calculates the standard deviation for the rows of a matrix.

    Dimensions of s, x, and m must match! 

    \retval  s   row standard deviations
    \param   x   a matrix
    \param   m   row means

    \author Patrick Eriksson 
    \date   2000-12-06
*/
void std_row( VectorView s, ConstMatrixView x, ConstVectorView m  )
{
  Vector d(x.ncols());		// We need this to store the deviation
				// from the mean.

  // Does the given mean provile match the number of rows of the matrix?
  assert( m.nelem()==x.nrows() );
  
  // Does the result vector match the number of rows of the matrix?
  assert( s.nelem()==x.nrows() );
  
  for ( Index i=0; i<x.nrows(); ++i ) 
    {
      // Put the deviation form the mean into d
      // d[j] = x[i,j]-m[j]
      d = x(i,Range(joker));
      d -= m[i];

      // Now we have to compute the square of d element-vise:
      d *= d;
	
      // Finally, take the sum and divide by the number of columns - 1.
      // FIXME: Patrick, I know this is not new, but: Should there
      // really be the -1 here?
      s[i] = d.sum() / (x.ncols()-1);
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



//// isbool ///////////////////////////////////////////////////////////////////
//
/** Checks if an integer is either 0 or 1.

    \return       a boolean, true if x is either 0 or 1
    \param    x   an integer

    \author Patrick Eriksson
    \date   2001-04-19
*/
bool isbool( const Index x ) 
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
    x[i] = start + i*step;
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
  assert( 1<n );		// Number of points must be greatere 1.
  x.resize(n);
  Numeric step = (stop-start)/(n-1) ;
  for ( Index i=0; i<n; i++ )
    x[i] = start + i*step;
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
  Numeric step = (log(stop)-a)/(n-1);
  x[0] = start;
  for ( Index i=1; i<n-1; i++ )
    x[i] = exp(a + i*step);
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




/////////////////////////////////////////////////////////////////////////////
//   Random data
/////////////////////////////////////////////////////////////////////////////

//// rand_uniform ////////////////////////////////////////////////////
/**
   Creates a vector with random data uniformerly distributed between
   the lower and higher limit given.

   The random data is uncorrelated. The length of the random vector is
   taken from r.nelem().

   \retval   r          random vector
   \param    x_low      lower limit for the random values
   \param    x_high     upper limit for the random data

   \author Patrick Eriksson
   \date   2000-11-26
*/
void rand_uniform(
		  VectorView  r,
		  const Numeric     x_low,
		  const Numeric     x_high )
{
  Numeric dx = x_high-x_low;

  for ( Index i=0; i<r.nelem(); i++ )
    r[i] = x_low + dx * (Numeric(rand())/Numeric(RAND_MAX));
}


//// rand_gaussian ////////////////////////////////////////////////////
/**
   Creates a gaussian random vector with zero mean and 
   the standard deviation given.

   The random data is uncorrelated. The length of the random vector to
   generate is taken from r.nelem().

   The algorith is taken from Numerical Recipies, Section 7.2. 
   See www.nr.com.

   \retval   r          random vector
   \param    s          standard deviation

   \author Patrick Eriksson
   \date   2000-11-27
*/
void rand_gaussian(
		   VectorView       r,
		   const Numeric   s )
{
  Vector  z(2);    // A vector of length 2 with uniform PDF between -1 and 1
  Numeric rad;     // The radius cooresponding to z
  Numeric fac;     // Normalisation factor
 
  const Index n = r.nelem();

  for ( Index i=0; i<n; )
  {
    rand_uniform( z, -1, 1 );
    rad = z[0]*z[0] + z[1]*z[1];

    if ( (rad<1) && (rad>0) )
    {
      fac = sqrt( -2*log(rad)/rad );
      r[i] = s*fac*z[0];
      i++;
      if ( i < n )
      {
        r[i] = s*fac*z[1];        
        i++;
      }
    }
  }
}

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
			 MatrixView       m,
			 const Numeric   x_low,
			 const Numeric   x_high )
{
  for ( Index i=0; i<m.nrows(); ++i )
  {
    rand_uniform( m(i,Range(joker)), x_low, x_high );
    // Matpack: m(i,Range(joker)) picks out the ith row of m. Because
    // rand_uniform takes an argument of type VectorView, it is
    // perfectly ok to call it with the row of a matrix. No
    // temporaries needed!
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
			  MatrixView    m,
			  const Numeric      s )
{
  for ( Index i=0; i<m.nrows(); ++i )
  {
    rand_gaussian( m(i,Range(joker)), s );
    // Matpack: m(i,Range(joker)) picks out the ith row of m. Because
    // rand_uniform takes an argument of type VectorView, it is
    // perfectly ok to call it with the row of a matrix. No
    // temporaries needed!
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
*/
void rand_data_gaussian(
			MatrixView       z,
			ConstVectorView  z0,
			ConstMatrixView   s )
{
  const Index n = z.ncols();

  const Index   nrows = z0.nelem();
  Index         col;

  // Check that s really is a square matrix:
  assert( s.ncols()==s.nrows() );
  
  // Check that the length of the mean vector is consistent with s: 
  assert ( nrows==s.nrows() );

  // Make Cholesky decomposition of s, l'*l=s
  Matrix   l(nrows,nrows);
  l = 0;			// Matpack can assign a scalar to all
				// elements of a matrix like this.
  chol(l,s);

  // Create matrix with gaussian data having zero mean and standard deviation 1
  Matrix   r(nrows,n);
  rand_matrix_gaussian( r, 1 );

  // Multiply l and r to get z. Note that the order in Matpack is
  // different from how it used to be with MTL.
  mult(z,l,r);

  // Add mean vector
  for ( col=0; col<n; col++ )
    z(Range(joker),col) += z0;
  // z(Range(joker),col) picks out a column of z. The += operator adds
  // the vector z0 to this element-vise.
}


////////////////////////////////////////////////////////////////////////////
//   Interpolation routines
////////////////////////////////////////////////////////////////////////////

/** Local help function to check input grids.

    \date   2000-06-29
    \author Patrick Eriksson 
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

  if ( (order*xi[0]<order*x[0]) || (order*xi[ni-1]>order*x[n-1]) ) 
    {
      ostringstream os;
      os << "Interpolation points must be inside the original range.\n"
	 << "Int.:  xi[0] = " << xi[0] << ", xi[ni-1] = " << xi[ni-1] << '\n'
	 << "Orig.: x[0]  = " << x[0]  << ", x[n-1]   = " << x[n-1];
      throw runtime_error(os.str());
    }

  for ( Index i=0; i<n-1; i++ )
  {
    if ( order*x[i+1] < order*x[i] ) 
      throw runtime_error("Original interpolation grid must be ordered");
  }

  for ( Index i=0; i<ni-1; i++ )
  {
    if ( order*xi[i+1] < order*xi[i] ) 
      throw runtime_error("Interpolation points must be ordered");
  }

  return order;
}

/** Multiple linear interpolation of a vector. Because Views are used,
    you can use matrix rows or columns, or vector sub-ranges inside
    the argument of this function.

    The vector x specifies the points at which the data y is given. 

    The size of yi has to be the same as for xi.

    \retval  yi      interpolated values 
    \param   x       the x grid
    \param   y       the function to interpolate
    \param   xi      interpolation points

    \date   2001-01-05
    \author Stefan Buehler
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
  Numeric      w;

  // Check that output vector has the right size:
  assert( n==yi.nelem() ); 

  for ( i=0; i<n; i++ )
  {
    for( ;  order*x[j+1] < order*xi[i]; j++ ) {}
    w = (xi[i]-x[j]) / (x[j+1]-x[j]);
    yi[i] = y[j] + w * (y[j+1]-y[j]); 
  }
}      

/** Multiple linear interpolation of matrix rows. Works with subranges
    inside the argument. Use the transpose(A) function if you want to
    interpolate columns of A instead of rows.

    The vector x specifies the points at which the data y is given. 

    \retval  Yi      interpolated values (matrix)
    \param   x       the x grid (vector)
    \param   Y       the function to interpolate (matrix)
    \param   xi      interpolation points (vector)

    \date   2001-01-06
    \author Stefan Buehler
*/
void interp_lin_matrix(    
		       MatrixView        Yi,
		       ConstVectorView   x, 
		       ConstMatrixView   Y, 
		       ConstVectorView   xi )
{
  // Check grids and get order of grids
  Index order = interp_check( x, xi, Y.ncols() ); 

  Index      j=0, n=xi.nelem(), nrow=Y.nrows();
  Numeric    w;

  assert( nrow == Yi.nrows() );
  assert( n    == Yi.ncols() );

  for (Index i=0; i<n; i++ )
  {
    for( ;  order*x[j+1] < order*xi[i]; j++ ) {}
    w = (xi[i]-x[j]) / (x[j+1]-x[j]);
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
//   Factorization of matrices
/////////////////////////////////////////////////////////////////////////////


//// chol //// //////////////////////////////////////////////////////////////

/** Choleski factorization (columnwise version). 

    Given c positive definite, the upper triangular matrix r with 
    positive diagonal elements such that c=r'*r is calculated.
    Algorithm used from Numerical Methods, Åke Björck, 1990, p46.

    \retval   r       Choleski factor of c
    \param    c       matrix to be factorized
 
    \author Carlos Jimenez
    \date   2001-02-14
*/

void chol(
	  MatrixView         r, 
          ConstMatrixView    c )
{
  const Index nrows = c.nrows(), ncols = c.ncols();
  Index j, i, k;
  Numeric a = 0;

  assert( nrows == ncols    );	// This makes only sense for square
				// matrices, doesn't it?
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
	      a = a + r(k,i) * r(k,j);
	}  
        r(i,j) = ( c(i,j) - a ) / r(i,i);       
      }
    
      a = 0;
      for (k=0; k<j; ++k)
	 a = a + r(k,j) * r(k,j);
      
    }
    r(j,j) = sqrt( c(j,j) - a);
    
  }

  // Checking that it works, if c is not positive definite it does not.
  for (i=0; i<nrows; ++i)
    for (j=0; j<nrows; ++j)    
      if ( isnan(r(i,j)) || isinf(r(i,j)) )
      {
	// 2001-09-15 FIXME: Patrick, I changed the condition above to
	// or. Before it was just a single &, I believe this was a
	// bug. Please verify. - Stefan
        ostringstream os;
        os << "Choleski decomposition does not work, c positive definite? \n";
        throw runtime_error(os.str());
      }
}

/////////////////////////////////////////////////////////////////////////////
//   Assert functions
/////////////////////////////////////////////////////////////////////////////

//// assert_bool ////////////////////////////////////////////////////////////
//
/** Asserts that an integer is 0 or 1.

    A runtime error is thrown if the integer is not a boolean.

    \param    x        an integer
    \param    x_name   the name of the integer (in upper case letters)

    \author Patrick Eriksson
    \date   2001-09-19
*/
void assert_bool( const Index& x, const String& x_name ) 
{
  if ( !(x==0 || x==1) )
  {
    ostringstream os;
    os << "The boolean " << x_name <<  " must either be 0 or 1.\n" 
       << "The present value of "<< x_name <<  " is " << x << ".";
    throw runtime_error( os.str() );
  }
}



//// assert_lengths (vector-vector) ///////////////////////////////////////////
//
/** Asserts that two vectors have the same length

    A runtime error is thrown if the lengths of the vector differ.

    \param    x1        vector 1
    \param    x1_name   the name of vector1 (in upper case letters)
    \param    x2        vector 2
    \param    x2_name   the name of vector2 (in upper case letters)

    \author Patrick Eriksson
    \date   2001-09-19
*/
void assert_lengths( const Vector& x1, const String& x1_name,
                     const Vector& x2, const String& x2_name ) 
{
  if ( x1.nelem() != x2.nelem() )
  {
    ostringstream os;
    os << "The vectors " << x1_name <<  " and " << x2_name << "\n"
       << "must have the same lengths. \nThe lengths are: \n"
       << x1_name << ": " << x1.nelem() << "\n"
       << x2_name << ": " << x2.nelem() << "\n";
    throw runtime_error( os.str() );
  }
}



//// assert_length_nrow  /////////////////////////////////////////////////////
//
/** Asserts that the length of a vector and the number of rows of a matrix
    match.

    A runtime error is thrown if the length of the vector differs from
    the number of rows.

    \param    x        a vector
    \param    x_name   the name of the vector (in upper case letters)
    \param    A        a matrix
    \param    A_name   the name of the matrix (in upper case letters)

    \author Patrick Eriksson
    \date   2001-09-19
*/
void assert_length_nrow( const Vector& x, const String& x_name,
                         const Matrix& A, const String& A_name ) 
{
  if ( x.nelem() != A.nrows() )
  {
    ostringstream os;
    os << "The length of vector " << x_name <<  " must be the same as \n"
       << "the number of rows of " << A_name << ".\n"
       << "The length of " << x_name <<  " is " << x.nelem() << ".\n"
       << "The number of rows of " << A_name <<  " is " << A.nrows() << ".\n";
    throw runtime_error( os.str() );
  }
}



//// assert_length_ncol  /////////////////////////////////////////////////////
//
/** Asserts that the length of a vector and the number of columns of a matrix
    match.

    A runtime error is thrown if the length of the vector differs from
    the number of columns.

    \param    x        a vector
    \param    x_name   the name of the vector (in upper case letters)
    \param    A        a matrix
    \param    A_name   the name of the matrix (in upper case letters)

    \author Patrick Eriksson
    \date   2001-09-19
*/
void assert_length_ncol( const Vector& x, const String& x_name,
                         const Matrix& A, const String& A_name ) 
{
  if ( x.nelem() != A.ncols() )
  {
    ostringstream os;
    os << "The length of vector " << x_name <<  " must be the same as \n"
       << "the number of columns of " << A_name << ".\n"
       << "The length of " << x_name <<  " is " << x.nelem() << ".\n"
       << "The number of columns of " << A_name <<  " is " << A.ncols()
       << ".\n";
    throw runtime_error( os.str() );
  }
}
