#include "arts.h"
#include "math_funcs.h"

//
// Basic mathematical vector functions
//
// 27.06.99 Patrick Eriksson   Included:
//          SQRT, EXP, MIN, MAX, FIRST and LAST
//
void sqrt( VECTOR& y, const VECTOR& x )
{
  int n = x.size();
  y.newsize(n);
  for ( int i=1; i<=n; i++ )
    y(i) = sqrt(x(i));
}

VECTOR sqrt( const VECTOR& x )
{
  VECTOR y;
  sqrt( y, x );
  return y; 
}


void exp( VECTOR& y, const VECTOR& x )
{
  int n = x.size();
  y.newsize(n);
  for ( int i=1; i<=n; i++ )
    y(i) = exp(x(i));
}

VECTOR exp( const VECTOR& x )
{
  VECTOR y;
  exp( y, x );
  return y; 
}

Numeric min( const VECTOR& x )    // Gives the minimum element of a vector
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


Numeric max( const VECTOR& x )     // Gives the maximum element of a vector
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


Numeric first( const VECTOR& x )    // Gives first element of a vector
{
  return x(1); 
}


Numeric last( const VECTOR& x )     // Gives last element of a vector
{
  return x(x.size()); 
}




//
// Logical vector functions
// 27.06.99 Patrick Eriksson   Included: ANY
//
bool any( const ARRAY<bool>& x )         // True if any element of x != 0
{
  for ( size_t i=1; i<=x.size(); i++ ) {
    if ( x(i) )
      return true;
  }
  return false;
}




//
// Functions to generate vectors
//
// 28.06.99 Patrick Eriksson   Included: LINSPACE and NLINSPACE
//       
void linspace(                      // Linearly spaced vector with spacing STEP
              VECTOR&  x,        // i.e. START:STEP:STOP. Note that last
        const Numeric     start,    // element can deviate from STOP.
        const Numeric     stop,        
        const Numeric     step )
{
  int n = (int) floor( (stop-start)/step ) + 1;
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
      

void nlinspace(                      // Linearly spaced vector of length N with
              VECTOR&  x,         // equally spaced points between START 
        const Numeric     start,     // and STOP.
        const Numeric     stop,        
        const int         n )
{
  #ifdef DO_CHECKS //--------  Safety check  ---------------------------------
    if ( n<2 )
      throw runtime_error("NLINSPACE: The number of points must be > 1"); 
  #endif //------------------  End safety check -------------------------------
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



//
// Basic mathematical matrix functions
//
// 27.06.99 Patrick Eriksson   Included: EXP and LOG
//

void exp( MATRIX& Y, const MATRIX& X )        // Elementwise exp
{
  int m = X.dim(1);
  int n = X.dim(2);
  Y.newsize(m,n);
  for ( int i=1; i<=m; i++ ) {
    for ( int j=1; j<=n; j++ )
      Y(i,j) = exp(X(i,j)); }
}

MATRIX exp( const MATRIX& X )                // Elementwise exp
{
  MATRIX Y;
  exp( Y, X );
  return Y; 
}


void log( MATRIX& Y, const MATRIX& X )       // Elementwise natural log.
{
  int m = X.dim(1);
  int n = X.dim(2);
  Y.newsize(m,n);
  for ( int i=1; i<=m; i++ ) {
    for ( int j=1; j<=n; j++ )
      Y(i,j) = log(X(i,j)); }
}

MATRIX log( const MATRIX& X )                // Elementwise natural log.
{
  MATRIX Y;
  log( Y, X );
  return Y; 
}




//
// Interpolation routines.
// All functions assume that the interpolation points, XI, are ordered.
//
// 27.06.99 Patrick Eriksson   Included:
//          INTERP_CHECK, INTERP_LIN and INTERP_LIN_ROW
//
void interp_check(                  // Local help function to check input
	const VECTOR&  x,        // grids.
        const VECTOR&  xi,
        const int         n_y )
{
  #ifdef DO_CHECKS //--------  Safety check  ---------------------------------

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

  #endif //------------------  End safety check  -----------------------------
}

    
void interp_lin(            // Linear interpolation of a vector
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

VECTOR interp_lin(            // As above but return version
        const VECTOR&  x, 
        const VECTOR&  y, 
        const VECTOR&  xi )
{
  VECTOR yi; 
  interp_lin( yi, x, y, xi );
  return yi;
}        

Numeric interp_lin(            // As above but for only one point
        const VECTOR&  x, 
        const VECTOR&  y, 
        const Numeric     xi )
{
  VECTOR yi(1,xi); 
  interp_lin( yi, x, y, xi );
  return yi(1);
}        

void interp_lin_row(       // Row-wise linear interpolation of a matrix
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
// Integration functions
//
// 27.06.99 Patrick Eriksson   Included: INTEGR_LIN
//
Numeric integr_lin(                // Integrates Y over X assuming that Y is
        const VECTOR&  x,       // linear between the given points
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


void integr_lin(                   // Integrates the rows of M assuming that
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
//

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


