/*-----------------------------------------------------------------------
FILE:      math_funcs.cc

INCLUDES:  This file contains basic mathematical functions.

           The functions are of different types:
	    1. Element-wise application of common scalar functions 
            2. Boolean functions
	    3. Creation of common vectors
	    4. Interpolation routines
	    5. Integration routines
            6. Conversion between vector and matrix types

FUNCTIONS: sqrt, exp, log, min, max, first, last 
           any
           linspace, nlinspace and nlogspace
           interp_lin, interp_lin_row
           integr_lin
           to_matrix, to_vector

HISTORY:   27.06.1999 Created by Patrick Eriksson.
	   23.03.2000 Stefan Buehler: Adapted to new ARRAY<>, VECTOR and
	              MATRIX types.
           07.04.2000 Patrick Eriksson: Added text for doc++ and nlogspace.
-----------------------------------------------------------------------*/

#include "arts.h"
#include "math_funcs.h"


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
void sqrt( VECTOR& y, const VECTOR& x )   // vector parameter version
{
  int n = x.size();
  y.newsize(n);
  for ( int i=1; i<=n; i++ )
    y(i) = sqrt(x(i));
}

VECTOR sqrt( const VECTOR& x )            // vector return version
{
  VECTOR y;
  sqrt( y, x );
  return y; 
}


/** Gives the elementwise exponential of a vector or matrix.
    Both return and parameter versions exist, i.e. you can type 

    exp(y,x)

    or

    y = exp(x)

    Note that x can be both a vector or a matrix.

    @param   y   Output: y, the exponential of x
    @param   x   The input vector or matrix

    @author Patrick Eriksson 27.06.99 */
void exp( VECTOR& y, const VECTOR& x )   // vector parameter version
{
  int n = x.size();
  y.newsize(n);
  for ( int i=1; i<=n; i++ )
    y(i) = exp(x(i));
}

VECTOR exp( const VECTOR& x )            // vector return version
{
  VECTOR y;
  exp( y, x );
  return y; 
}

void exp( MATRIX& Y, const MATRIX& X )   // matrix parameter version
{
  int m = X.dim(1);
  int n = X.dim(2);
  Y.newsize(m,n);
  for ( int i=1; i<=m; i++ ) {
    for ( int j=1; j<=n; j++ )
      Y(i,j) = exp(X(i,j)); }
}

MATRIX exp( const MATRIX& X )            // matrix return version
{
  MATRIX Y;
  exp( Y, X );
  return Y; 
}


/** Gives the elementwise natural logarith of a vector or matrix.
    Both return and parameter versions exist, i.e. you can type 

    log(y,x)

    or

    y = log(x)

    Note that x can be both a vector or a matrix.

    @param   y   Output: y, the natural logarith of x
    @param   x   The input vector or matrix

    @author Patrick Eriksson 27.06.99 */
void log( VECTOR& y, const VECTOR& x )   // vector parameter version
{
  int n = x.size();
  y.newsize(n);
  for ( int i=1; i<=n; i++ )
    y(i) = log(x(i));
}

VECTOR log( const VECTOR& x )            // vector return version
{
  VECTOR y;
  log( y, x );
  return y; 
}

void log( MATRIX& Y, const MATRIX& X )   // matrix parameter version
{
  int m = X.dim(1);
  int n = X.dim(2);
  Y.newsize(m,n);
  for ( int i=1; i<=m; i++ ) {
    for ( int j=1; j<=n; j++ )
      Y(i,j) = log(X(i,j)); }
}

MATRIX log( const MATRIX& X )            // matrix return version
{
  MATRIX Y;
  log( Y, X );
  return Y; 
}


/** Gives the minimum value of a vector, v = min(x).

    @param   v   Returns: the minimum value of x
    @param   x   The input vector

    @author Patrick Eriksson 27.06.99 */
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


/** Gives the maximum value of a vector, v = min(x).

    @param   v   Returns: the maximum value of x
    @param   x   The input vector

    @author Patrick Eriksson 27.06.99 */
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


/** Gives the first value of a vector, v = first(x).

    @param   v   Returns: the first value of x
    @param   x   The input vector

    @author Patrick Eriksson 27.06.99 */
Numeric first( const VECTOR& x )
{
  return x(1); 
}


/** Gives the last value of a vector, v = last(x).

    @param   v   Returns: the last value of x
    @param   x   The input vector

    @author Patrick Eriksson 27.06.99 */
Numeric last( const VECTOR& x )
{
  return x(x.size()); 
}




//
// Logical vector functions
//   ANY
//

/** True if any element of a boolean vector, b, is not 0, b = any(x).

    @param   v   Returns: a boolean, true if any x is != 0
    @param   x   The input vector

    @author Patrick Eriksson 27.06.99 */
bool any( const ARRAY<int>& x )         // True if any element of x != 0
{
  for ( size_t i=1; i<=x.size(); i++ ) {
    if ( x(i) )
      return true;
  }
  return false;
}




//
// Functions to generate vectors
//   LINSPACE, NLINSPACE, and LOGNSPACE
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
        const Numeric     step )
{
  int n = (int) floor( (stop-start)/step ) + 1;
  if ( n<1 )
    n=1;
  x.newsize(n);
  for ( int i=1; i<=n; i++ )
    x(i) = start + (i-1)*step;
}

VECTOR linspace(                 // As above but return version
        const Numeric start, 
        const Numeric stop,  
        const Numeric step )
{
  VECTOR x;
  linspace( x, start, stop, step );
  return x; 
}
      

/** Linearly spaced vector with specified length. 
    Both return and parameter versions exist:

    nlinspace(x,start,stop,n)

    x = nlinspace(start,stop,n)

    Returns a vector equally and linearly spaced between start and stop of length n. 
    (equals the Matlab function linspace)
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
        const int         n )
{
  if ( n<2 )
    throw runtime_error("NLINSPACE: The number of points must be > 1"); 
  x.newsize(n);
  Numeric step = (stop-start)/(n-1) ;
  for ( int i=1; i<=n; i++ )
    x(i) = start + (i-1)*step;
}

VECTOR nlinspace(                 // As above but return version
        const Numeric start, 
        const Numeric stop,  
        const int     n )
{
  VECTOR x;
  nlinspace( x, start, stop, n );
  return x; 
}                     


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
        const int         n )
{
  if ( n<2 )
    throw runtime_error("NLOGSPACE: The number of points must be > 1"); 
  if ( (start<=0) || (stop<=0) )
    throw runtime_error("NLOGSPACE: Only positive numbers are allowed"); 
  x.newsize(n);
  Numeric a = log(start);
  Numeric step = (log(stop)-a)/(n-1) ;
  x(1) = start;
  for ( int i=2; i<n; i++ )
    x(i) = exp(a + (i-1)*step);
  x(n) = stop;
}

VECTOR nlogspace(                 // As above but return version
        const Numeric start, 
        const Numeric stop,  
        const int     n )
{
  VECTOR x;
  nlogspace( x, start, stop, n );
  return x; 
}                     




//
// Interpolation routines.
//   Vectors:  INTERP_LIN
//   Matrices: INTERP_LIN_ROW
//
//   All thefunctions assume that the interpolation points, XI, are ordered.
//

// Local help function to check input grids (PE 29.06.99)
void interp_check(                  
	const VECTOR&  x,        
        const VECTOR&  xi,
        const int      n_y )
{
  int n  = x.size();
  int ni = xi.size();

  if ( (xi(1)<x(1)) || (xi(ni)>x(n)) ) 
   throw runtime_error("INTERPOLATION: Interpolation points must be inside the range of X");

  for (int i=1; i<ni; i++ )
  {
    if ( xi(i+1) < xi(i) ) 
      throw runtime_error("INTERPOLATION: Interpolation points must be ordered");
  }

  if ( n != n_y ) 
    throw runtime_error("INTERPOLATION: Sizes of input data do not match");
}

    
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
void interp_lin(                  // Linear interpolation of a vector
              VECTOR&  yi,
        const VECTOR&  x, 
        const VECTOR&  y, 
        const VECTOR&  xi )
{
  interp_check( x, xi, y.size() ); 

  int        j=1, n=xi.size();
  Numeric    w;
  yi.newsize(n); 
  
  for (int i=1; i<=n; i++ )
  {
    for( ;  x(j+1) < xi(i); j++ ) {}
    w = (xi(i)-x(j)) / (x(j+1)-x(j));
    yi(i) = y(j) + w * (y(j+1)-y(j)); 
  }
}      

VECTOR interp_lin(               // As above but return version
        const VECTOR&  x, 
        const VECTOR&  y, 
        const VECTOR&  xi )
{
  VECTOR yi; 
  interp_lin( yi, x, y, xi );
  return yi;
}        

Numeric interp_lin(              // As above but for only one point
        const VECTOR&  x, 
        const VECTOR&  y, 
        const Numeric  xi )
{
  VECTOR yi(1,xi); 
  interp_lin( yi, x, y, xi );
  return yi(1);
}        

/** Linear interpolation of the rows of a matrix.
    Both return and parameter versions exist:

    interp_lin_row(Yi,x,Y,xi)

    Yi = interp_lin(x,Y,xi)

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
        const VECTOR&  xi )
{
  interp_check( x, xi, Y.dim(2) ); 

  int        k, j=1, n=xi.size(), nrow=Y.dim(1);
  Numeric    w;
  Yi.newsize( nrow, n ); 

  for (int i=1; i<=n; i++ )
  {
    for( ;  x(j+1) < xi(i); j++ ) {}
    w = (xi(i)-x(j)) / (x(j+1)-x(j));
    for( k=1; k<=nrow; k++ )
      Yi(k,i) = Y(k,j) + w * (Y(k,j+1)-Y(k,j)); 
  }
}        

MATRIX interp_lin_row(       // As above but return version
        const VECTOR&  x, 
        const MATRIX&  Y, 
        const VECTOR&  xi )
{
  MATRIX Yi; 
  interp_lin_row( Yi, x, Y, xi );
  return Yi;
}        




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

    @author Patrick Eriksson 29.06.99 */
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
  W.newsize(rows,1);
  W = 0.0;

  if ( cols < 2 )
    throw runtime_error("INTEGR_LIN: Vector length must be >= 2");
  if ( cols != x.size() )
    throw runtime_error("INTEGR_LIN: Sizes of input data do not match");

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
void to_matrix(MATRIX& W, const VECTOR& x)
{
  // FIXME: I'm sure this can be made more efficient when TNT has more
  // functionality.  
  W.newsize(x.dim(),1);
  //  cout << "W size = " << W.num_rows() << " " << W.num_cols() << endl;
  for (size_t i=1; i<=x.dim() ; ++i)
    {
      //      cout << "i = " << i << endl;
      W(i,1) = x(i);
    }
}

MATRIX to_matrix(const VECTOR& x)
{
  // FIXME: I'm sure this can be made more efficient when TNT has more
  // functionality.
  MATRIX W(x.dim(),0);
  for (size_t i=1; i<=x.dim() ; ++i)
    W(i,0) = x(i);
  return W;
}




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
//
//  Matrix to vector conversion:
//
void to_vector(VECTOR& x, const MATRIX& W)
{
  // FIXME: I'm sure this can be made more efficient when TNT has more
  // functionality.

  // Check if one of the dimensions of W is 1:
  if ( 1 == W.num_cols() )
    {
      x.newsize(W.num_rows());
      for (size_t i=1; i<=x.dim() ; ++i)
	x(i) = W(i,1);
    }
  else if ( 1 == W.num_rows() )
    {
      x.newsize(W.num_cols());
      for (size_t i=1; i<=x.dim() ; ++i)
	x(i) = W(1,i);
    }
  else
    throw runtime_error("You tried to convert a matrix to a vector,\n"
			"but none of the dimensions is 1.");
}

VECTOR to_vector(const MATRIX& W)
{
  // FIXME: I'm sure this can be made more efficient when TNT has more
  // functionality.

  VECTOR x;

  // Check if one of the dimensions of W is 1:
  if ( 1 == W.num_cols() )
    {
      x.newsize(W.num_rows());
      for (size_t i=1; i<=x.dim() ; ++i)
	x(i) = W(i,1);
    }
  else if ( 1 == W.num_rows() )
    {
      x.newsize(W.num_cols());
      for (size_t i=1; i<=x.dim() ; ++i)
	x(i) = W(1,i);
    }
  else
    throw runtime_error("You tried to convert a matrix to a vector,\n"
			"but none of the dimensions is 1.");

  return x;
}


