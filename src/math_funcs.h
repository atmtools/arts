#ifndef math_funcs_h
#define math_funcs_h

#include "vecmat.h"



//
// Basic mathematical vector and vector functions:
//   Vectors: SQRT, EXP, LOG, MIN, MAX, FIRST and LAST
//   Matrices: EXP and LOG
//

/** Gives the elementwise square root of a vector.
    Both return and parameter versions exist, i.e. you can type 

    sqrt(y,x)

    or

    y = sqrt(x)

    @param   y   Output: the square root of x
    @param   x   The input vector 

    @author Patrick Eriksson 27.06.99 */
void sqrt( VECTOR& y, const VECTOR& x );
VECTOR sqrt( const VECTOR& x );


/** Gives the elementwise exponential of a vector or matrix.
    Both return and parameter versions exist, i.e. you can type 

    exp(y,x)

    or

    y = exp(x)

    Note that x can be both a vector or a matrix.

    @param   y   Output: y, the exponential of x
    @param   x   The input vector or matrix

    @author Patrick Eriksson 27.06.99 */
void exp( VECTOR& y, const VECTOR& x );
VECTOR exp( const VECTOR& x );
void exp( MATRIX& Y, const MATRIX& X );
MATRIX exp( const MATRIX& X );


/** Gives the elementwise natural logarith of a vector or matrix.
    Both return and parameter versions exist, i.e. you can type 

    log(y,x)

    or

    y = log(x)

    Note that x can be both a vector or a matrix.

    @param   y   Output: y, the natural logarith of x
    @param   x   The input vector or matrix

    @author Patrick Eriksson 27.06.99 */
void log( VECTOR& y, const VECTOR& x );
VECTOR log( const VECTOR& x );
void log( MATRIX& Y, const MATRIX& X );
MATRIX log( const MATRIX& X );


/** Gives the minimum value of a vector, v = min(x).

    @param   v   Returns: the minimum value of x
    @param   x   The input vector

    @author Patrick Eriksson 27.06.99 */
Numeric min( const VECTOR& x );

/** Gives the minimum value of a matrix, v = min(A).

    @param   v   Returns: the minimum value of x
    @param   A   The input matrix

    @author Stefan Buehler 08.05.2000 */
Numeric min( const MATRIX& A );

/** Gives the minimum value of an array, v = min(x).

    @param   v   Returns: the minimum value of x
    @param   x   The input array

    @author Stefan Buehler 08.05.2000 */
template<class T>
T min( const ARRAY<T>& x );


/** Gives the maximum value of a vector, v = min(x).

    @param   v   Returns: the maximum value of x
    @param   x   The input vector

    @author Patrick Eriksson 27.06.99 */
Numeric max( const VECTOR& x );

/** Gives the maximum value of a matrix, v = max(A).

    @param   v   Returns: the maximum value of x
    @param   A   The input matrix

    @author Stefan Buehler 08.05.2000 */
Numeric max( const MATRIX& A );


/** Gives the maximum value of an array, v = max(x). Because this is a
    template function, the definition has to be also in the header
    file, and not in file math_func.cc.

    @param   v   Returns: the maximum value of x
    @param   x   The input array

    @author Stefan Buehler 08.05.2000 */
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


/** Gives the first value of a vector, v = first(x).

    @param   v   Returns: the first value of x
    @param   x   The input vector

    @author Patrick Eriksson 27.06.99 */
Numeric first( const VECTOR& x );


/** Gives the last value of a vector, v = last(x).

    @param   v   Returns: the last value of x
    @param   x   The input vector

    @author Patrick Eriksson 27.06.99 */
Numeric last( const VECTOR& x );



//
// Logical vector functions
//   ANY
//

/** True if any element of a boolean vector, b, is not 0, b = any(x).

    @param   v   Returns: a boolean, true if any x is != 0
    @param   x   The input vector

    @author Patrick Eriksson 27.06.99 */
bool any( const ARRAY<int>& x ); 



//
// Functions to generate vectors
//   LINSPACE, NLINSPACE and NLOGSPACE
//   

/** Linearly spaced vector with specified spacing. 
    Both return and parameter versions exist:

    linspace(x,start,stop,step)

    x = linspace(start,stop,step)

    The first element of x is always start. The next value is start+step etc.
    Note that the last value can deviate from stop.
    The step can be both positive and negative. 
    (in Matlab notation: start:step:stop)

    @param   x       Output: a linearly spaced vector 
    @param   start   The start value
    @param   stop    The stop value
    @param   step    The spacing of the vector

    @author Patrick Eriksson 27.06.99 */
void linspace(                      
              VECTOR&     x,           
        const Numeric     start,    
        const Numeric     stop,        
        const Numeric     step );
VECTOR linspace(             
        const Numeric     start, 
        const Numeric     stop,  
        const Numeric     step );


/** Linearly spaced vector with specified length. 
    Both return and parameter versions exist:

    nlinspace(x,start,stop,n)

    x = nlinspace(start,stop,n)

    Returns a vector equally and linearly spaced between start and stop 
    of length n. (equals the Matlab function linspace)
    The length must be > 1.

    @param   x       Output: a linearly spaced vector 
    @param   start   The start value
    @param   stop    The stop value
    @param   n       The vector length (>1)

    @author Patrick Eriksson 27.06.99 */
void nlinspace(         
              VECTOR&     x, 
        const Numeric     start,     
        const Numeric     stop,        
        const int         n );
VECTOR nlinspace(         
        const Numeric     start, 
        const Numeric     stop,  
        const int         n );


/** Logarithmically spaced vector with specified length. 
    Both return and parameter versions exist:

    nlogspace(x,start,stop,n)

    x = nlogspace(start,stop,n)

    Returns a vector logarithmically spaced vector between start and 
    stop of length n (equals the Matlab function logspace)
    The length must be > 1.

    @param   x       Output: a logarithmically spaced vector 
    @param   start   The start value
    @param   stop    The stop value
    @param   n       The vector length (>1)

    @author Patrick Eriksson 07.04.00 */
void nlogspace(         
              VECTOR&     x, 
        const Numeric     start,     
        const Numeric     stop,        
        const int         n );
VECTOR nlogspace(  
        const Numeric     start, 
        const Numeric     stop,  
        const int         n );



//
// Interpolation routines.
//   Vectors:  INTERP_LIN
//   Matrices: INTERP_LIN_ROW
//

/** Linear interpolation of a vector.
    Both return and parameter versions exist:

    interp_lin(yi,x,y,xi)

    yi = interp_lin(x,y,xi)

    The vector x specifies the points at which the data y is given. 
    The interpolation points, xi, can either be a scalar (Numeric) or a vector.
    The type of yi is the same as for xi.

    @param   yi      Output: interpolated values 
    @param   x       The x grid
    @param   y       The function to interpolate
    @param   xi      Interpolation point(s)

    @author Patrick Eriksson 29.06.99 */
void interp_lin(            
              VECTOR&  yi,
        const VECTOR&  x, 
        const VECTOR&  y, 
        const VECTOR&  xi );
VECTOR interp_lin(          
        const VECTOR&  x, 
        const VECTOR&  y, 
        const VECTOR&  xi );
Numeric interp_lin(         
        const VECTOR&  x, 
        const VECTOR&  y, 
        const Numeric  xi );


/** Linear interpolation of the rows of a matrix.
    Both return and parameter versions exist:

    interp_lin_row(Yi,x,Y,xi)

    Yi = interp_lin_row(x,Y,xi)

    The vector x specifies the points at which the data y is given. 

    @param   Yi      Output: interpolated values 
    @param   x       The x grid
    @param   Y       The functions to interpolate
    @param   xi      Interpolation points

    @author Patrick Eriksson 29.06.99 */
void interp_lin_row(    
              MATRIX&  Yi,
        const VECTOR&  x, 
        const MATRIX&  Y, 
        const VECTOR&  xi );
MATRIX interp_lin_row(      
        const VECTOR&  x, 
        const MATRIX&  Y, 
        const VECTOR&  xi );



//
// Integration functions for vectors and matrices
//    INTEGR_LIN
//

/** Integration of vectors and matrices assuming piecewise linear functions.
    Both return and parameter versions exist:

    integrs_lin(w,x,y)

    w = interp_lin(x,y)

    The vector x specifies the points at which the data y is given. 
    The function y is treated to be linear between the points of x and zero outside 
    the defined range.

    @param   w       Output: the integrated value 
    @param   x       The x grid
    @param   y       The function to integrate

    @author Patrick Eriksson 12.04.00 */
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
        const MATRIX&  M );



//
// Vector to matrix conversion:
//   to_matrix
//

/** Converts a vector to a matrix. 
    For a vector of length n, the dimension of the matrix is [n,1], in other 
    words, the vector is interpreted as a column vector.
    Both return and parameter versions exist:

    to_matrix(W,x)

    W = to_matrix(x)

    @param   W       Output: the matrix, size n x 1
    @param   x       The vector, length n

    @author Stefan Buehler */
void to_matrix(MATRIX& W, const VECTOR& x);
MATRIX to_matrix(const VECTOR& x);



//
// Matrix to vector conversion:
//   to_vector
//

/** Conversion of a matrix to a vector.
    The matrix can either be a column (n x 1) or row (1 x n) vector.
    Both return and parameter versions exist:

    to_vector(x,W)

    x = to_vector(W)

    @param   x       Output: the vector, length n
    @param   W       The matrix, size (n x 1) or (1 x n)

    @exception runtime_error None of the dimensions of W was 1.

    @author Stefan Buehler */
void to_vector(VECTOR& x, const MATRIX& W);
VECTOR to_vector(const MATRIX& W);


#endif  // math_funcs_h
