/*-----------------------------------------------------------------------
FILE:       math_funcs.h


INCLUDES:   See MATH_FUNCS.H for inventory


         This file contains basic mathematical functions.

         The functions are of different types, e.g.:
	 1. Element-wise application of common scalar functions 
	 2. Creation of common vector and matrices
	 3. Interpolation routines
	 4. Integration routines

HISTORY: 27.06.1999 Created by Patrick Eriksson.
         See further below.
	 23.03.2000 Stefan Buehler: Adapted to new ARRAY<>, VECTOR and
	 MATRIX types.
  
-----------------------------------------------------------------------*/

#ifndef math_funcs_h
#define math_funcs_h

#include "vecmat.h"

//
// Basic mathematical vector functions
//
void   	sqrt( VECTOR& y, const VECTOR& x ); // Elementwise sqrt
VECTOR 	sqrt( const VECTOR& x );            // Elementwise sqrt
void   	exp( VECTOR& y, const VECTOR& x );  // Elementwise exp
VECTOR 	exp( const VECTOR& x );             // Elementwise exp
Numeric min( const VECTOR& x );       	    // Gives min. element of a vector
Numeric max( const VECTOR& x );       	    // Gives max. element of a vector
Numeric first( const VECTOR& x );     	    // Gives first element of a vector
Numeric last( const VECTOR& x );      	    // Gives last element of a vector


//
// Logical vector functions
//
bool      any( const ARRAY<bool>& x );       // True if any element of x != 0


//
// Functions to generate vectors
//
void linspace(                      // Linearly spaced vector with spacing STEP
              VECTOR&  x,        // i.e. START:STEP:STOP. Note that last
        const Numeric     start,    // element can deviate from STOP.
        const Numeric     stop,        
        const Numeric     step );
VECTOR linspace(                 // As above but return version
      const Numeric  start, 
      const Numeric  stop, 
      const Numeric  step );

void nlinspace(                      // Linearly spaced vector of length N with
              VECTOR&  x,         // equally spaced points between START 
        const Numeric     start,     // and STOP.
        const Numeric     stop,        
        const int         n );
VECTOR nlinspace(                 // As above but return version
        const Numeric start, 
        const Numeric stop,  
        const int     n );


//
// Basic mathematical matrix functions
//
void      exp( MATRIX& Y, const MATRIX& X );    // Elementwise exp
MATRIX exp( const MATRIX& x );                  // Elementwise exp
void      log( MATRIX& Y, const MATRIX& X );    // Elementwise natural log.
MATRIX log( const MATRIX& x );                  // Elementwise natural log.


//
// Interpolation routines.
// All functions assume that the interpolation points, XI, are ordered.
//
void interp_lin(                // Linear interpolation of a vector
		VECTOR& yi,  // Length of x and y must be equal.
       const VECTOR& x, 
       const VECTOR& y, 
       const VECTOR& xi );        
VECTOR interp_lin(           // As above but return version
       const VECTOR& x,  
       const VECTOR& y, 
       const VECTOR& xi );        
Numeric interp_lin(            // As above but for only one point
        const VECTOR&  x, 
        const VECTOR&  y, 
        const Numeric     xi );

void interp_lin_row(            // Row-wise linear interpolation of a matrix,
             MATRIX& YI,     // i.e. the length of x shall equal the
       const VECTOR& x,      // number of columns of Y.
       const MATRIX& Y, 
       const VECTOR& xi );        
MATRIX interp_lin_row(       // As above but return version
       const VECTOR& x,  
       const MATRIX& Y,  
       const VECTOR& xi );        


//
// Integration functions
//
Numeric integr_lin( const VECTOR& x, const VECTOR& y );
                           // Integrates Y over X assuming that Y is
                           // linear between the given points
void integr_lin( MATRIX&, const VECTOR& x, const MATRIX& M );
                           // As above but works with a matrix instead
                           // of a vector
MATRIX integr_lin( const VECTOR& x, const MATRIX& M );
                           // As above but return version


#endif  // math_funcs_h
