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
   \file   math_funcs.h

   Contains declerations of basic mathematical and vector/matrix functions.

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



#ifndef math_funcs_h
#define math_funcs_h

////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "vecmat.h"



////////////////////////////////////////////////////////////////////////////
//   Basic mathematical vector and vector functions
////////////////////////////////////////////////////////////////////////////

/** Gives the elementwise square root of a vector.
   
   \retval   y   the square root of x
   \param    x   a vector

   \author Patrick Eriksson
   \date   2000-06-27
*/
void sqrt( VECTOR& y, const VECTOR& x );

/** Gives the elementwise square root of a vector (return version).
   
   \return       the square root of x
   \param    x   a vector

   \author Patrick Eriksson
   \date   2000-06-27
*/
VECTOR sqrt( const VECTOR& x );



/** Gives the elementwise exponential of a vector.
   
   \retval   y   the exponential of x
   \param    x   a vector

   \author Patrick Eriksson
   \date   2000-06-27
*/
void exp( VECTOR& y, const VECTOR& x );

/** Gives the elementwise exponential of a vector (return version).
   
   \return       the exponential of x
   \param    x   a vector

   \author Patrick Eriksson
   \date   2000-06-27
*/
VECTOR exp( const VECTOR& x );

/** Gives the elementwise exponential of a matrix.
   
   \retval   y   the exponential of x
   \param    x   a matrix

   \author Patrick Eriksson
   \date   2000-06-27
*/
void exp( MATRIX& Y, const MATRIX& X );

/** Gives the elementwise exponential of a matrix (return version).
   
   \return       the exponential of x
   \param    x   a matrix

   \author Patrick Eriksson
   \date   2000-06-27
*/
MATRIX exp( const MATRIX& X );



/** Gives the elementwise natural logarithm of a vector.
   
   \retval   y   the natural logarithm of x
   \param    x   a vector

   \author Patrick Eriksson
   \date   2000-06-27
*/
void log( VECTOR& y, const VECTOR& x );

/** Gives the elementwise natural logarithm of a vector (return version).
   
   \return       the natural logarithm of x
   \param    x   a vector

   \author Patrick Eriksson
   \date   2000-06-27
*/
VECTOR log( const VECTOR& x );

/** Gives the elementwise natural logarithm of a matrix.
   
   \retval   y   the natural logarithm of x
   \param    x   a matrix

   \author Patrick Eriksson
   \date   2000-06-27
*/
void log( MATRIX& Y, const MATRIX& X );

/** Gives the elementwise natural logarithm of a matrix (return version).
   
   \return       the natural logarithm of x
   \param    x   a matrix

   \author Patrick Eriksson
   \date   2000-06-27
*/
MATRIX log( const MATRIX& X );


/** Gives the minimum value of a vector.

    \return      the minimum value of x
    \param   x   a vector

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric min( const VECTOR& x );

/** Gives the minimum value of a matrix.

    \return      the minimum value of A
    \param   A   a matrix

    \autho  Stefan Buehler
    \date   2000-06-27
*/
Numeric min( const MATRIX& A );

/** Gives the maximum value of a vector.

    \return      the maximum value of x
    \param   x   a vector

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric max( const VECTOR& x );

/** Gives the maximum value of a matrix.

    \return       the maximum value of A
    \param    A   a matrix

    \autho  Stefan Buehler
    \date   2000-06-27
*/
Numeric max( const MATRIX& A );


/** Gives the maximum value of an array.

    Because this is a template function, the definition has to be also 
    in the header file, and not in file math_func.cc.

    \return      the maximum value of x
    \param   x   an array

    \author Stefan Buehler
    \date   2000-06-27
*/
template<class T>
T max( const ARRAY<T>& x )
{
  size_t n = x.dim();
  T y=x(1);
  for ( size_t i=2; i<=n; i++ )
  {
    if ( x(i) > y )
      y = x(i);
  }
  return y; 
}

/** Gives the minimum value of an array.

    Because this is a template function, the definition has to be also 
    in the header file, and not in file math_func.cc.

    \return      the minimum value of x
    \param   x   an array

    \author Stefan Buehler
    \date   2000-06-27
*/
template<class T>
T min( const ARRAY<T>& x )
{
  size_t n = x.dim();
  T y=x(1);
  for ( size_t i=2; i<=n; i++ )
  {
    if ( x(i) < y )
      y = x(i);
  }
  return y; 
}


/** Gives the first value of a vector.

    \return       the first value of x
    \param    x   a vector

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric first( const VECTOR& x );

/** Gives the last value of a vector.

    \return      the last value of x
    \param   x   a vector

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric last( const VECTOR& x );



////////////////////////////////////////////////////////////////////////////
//   Logical functions
////////////////////////////////////////////////////////////////////////////

/** True if any element of a boolean vector, b, is not 0.

    \return       a boolean, true if any x is != 0
    \param    x   a vector

    \author Patrick Eriksson
    \date   2000-06-27
*/
bool any( const ARRAY<int>& x ); 



////////////////////////////////////////////////////////////////////////////
// Functions to generate vectors
////////////////////////////////////////////////////////////////////////////

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
        const Numeric  start,    
        const Numeric  stop,        
        const Numeric  step );

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
        const Numeric  start, 
        const Numeric  stop,  
        const Numeric  step );



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
        const int         n );

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
VECTOR nlinspace(         
        const Numeric     start, 
        const Numeric     stop,  
        const int         n );



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
        const int         n );

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
        const Numeric     start, 
        const Numeric     stop,  
        const int         n );



////////////////////////////////////////////////////////////////////////////
//   Interpolation routines
////////////////////////////////////////////////////////////////////////////

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
        const VECTOR&  xi );

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
        const VECTOR&  xi );

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
        const Numeric  xi );


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
        const VECTOR&  xi );

/** Multiple linear interpolation of matrix rows (return version).

    The vector x specifies the points at which the data y is given. 

    \retval  Yi      interpolated values 
    \param   x       the x grid
    \param   Y       the function to interpolate
    \param   xi      interpolation points

    \author Patrick Eriksson
    \date   2000-06-29
*/
MATRIX interp_lin_row(      
        const VECTOR&  x, 
        const MATRIX&  Y, 
        const VECTOR&  xi );

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
        const VECTOR&  xi );

/** Multiple linear interpolation of matrix columns (return version).

    The vector x specifies the points at which the data y is given. 

    \return          interpolated values 
    \param   x       the x grid
    \param   Y       the function to interpolate
    \param   xi      interpolation points

    \author Stefan Buehler
    \date   2000-06-29
*/
MATRIX interp_lin_col(      
        const VECTOR&  x, 
        const MATRIX&  Y, 
        const VECTOR&  xi );



/////////////////////////////////////////////////////////////////////////////
//   Integration functions for vectors and matrices
//     These functions are not used for the moment
/////////////////////////////////////////////////////////////////////////////
/*
Numeric integr_lin(        
        const VECTOR&  x,  
        const VECTOR&  y );
void integr_lin(            
              Numeric&  w,
        const VECTOR&   x,  
        const VECTOR&   y ); 
void integr_lin(           
              MATRIX&  W,  
        const VECTOR&  x,  
        const MATRIX&  M );   
MATRIX integr_lin(         
        const VECTOR&  x,  
        const MATRIX&  M ); */



/////////////////////////////////////////////////////////////////////////////
//   Conversions between VECTOR and MATRIX types
/////////////////////////////////////////////////////////////////////////////

/** Converts a vector to a matrix. 

    For a vector of length n, the dimension of the matrix is [n,1], in other 
    words, the vector is interpreted as a column vector.

    \retval   W       the matrix, size n x 1
    \param    x       a vector of length n

    \author Stefan Buehler 
    \date 2000-09-01
*/
void to_matrix(MATRIX& W, const VECTOR& x);

/** Converts a vector to a matrix (return version). 

    For a vector of length n, the dimension of the matrix is [n,1], in other 
    words, the vector is interpreted as a column vector.

    \return           the matrix, size n x 1
    \param    x       a vector of length n

    \author Stefan Buehler 
    \date 2000-09-01
*/
MATRIX to_matrix(const VECTOR& x);



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
void to_vector(VECTOR& x, const MATRIX& W);

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
VECTOR to_vector(const MATRIX& W);



/////////////////////////////////////////////////////////////////////////////
//   Extraction of matrix columns and rows
//
//     x = row(i,A), row(x,i,A)
//             Generates a vector which contains row i of A.
//
//     x = col(i,A), col(x,i,A)
//             Generates a vector which contains column i of A.
//
/////////////////////////////////////////////////////////////////////////////

/** Extracts row i of MATRIX A.

    \retval   x   row i of A as a vector
    \param    i   row index
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
void row(VECTOR& x,
	 size_t i,
	 const MATRIX& A);

/** Extracts row i of MATRIX A (return version).

    \return       row i of A as a vector
    \param    i   row index
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
VECTOR row(size_t i,
	   const MATRIX& A);

/** Extracts column i of MATRIX A.

    \retval   x   column i of A as a vector
    \param    i   column index
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
void col(VECTOR& x,
	 size_t i,
	 const MATRIX& A);

/** Extracts column i of MATRIX A (return version).

    \return       column i of A as a vector
    \param    i   column index
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
VECTOR col(size_t i,
	   const MATRIX& A);



/** Extracts rows i to k of MATRIX A.

    \retval   X   the extracted part of A
    \param    i   first row to extract
    \param    k   last row to extract
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
void row(MATRIX& X,
	 size_t i,
	 size_t k,
	 const MATRIX& A);

/** Extracts rows i to k of MATRIX A (return version).

    \return       the extracted part of A
    \param    i   first row to extract
    \param    k   last row to extract
    \param    A   a matrix.

    \author Stefan Buehler 
    \date   2000-09-01
*/
MATRIX row(size_t i,
	   size_t k,
	   const MATRIX& A);

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
	 const MATRIX& A);

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
	   const MATRIX& A);



#endif  // math_funcs_h
