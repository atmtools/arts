#ifndef math_funcs_h
#define math_funcs_h

#include "vecmat.h"

//
// Basic mathematical vector and vector functions:
//   Vectors: SQRT, EXP, LOG, MIN, MAX, FIRST and LAST
//   Matrices: EXP and LOG
//
void sqrt( VECTOR& y, const VECTOR& x );
VECTOR sqrt( const VECTOR& x );
//
void exp( VECTOR& y, const VECTOR& x );
VECTOR exp( const VECTOR& x );
void exp( MATRIX& Y, const MATRIX& X );
MATRIX exp( const MATRIX& X );
//
void log( VECTOR& y, const VECTOR& x );
VECTOR log( const VECTOR& x );
void log( MATRIX& Y, const MATRIX& X );
MATRIX log( const MATRIX& X );
//
Numeric min( const VECTOR& x );
Numeric max( const VECTOR& x );
Numeric first( const VECTOR& x );
Numeric last( const VECTOR& x );


//
// Logical vector functions
//   ANY
//
bool any( const ARRAY<int>& x ); 


//
// Functions to generate vectors
//   LINSPACE, NLINSPACE and NLOGSPACE
//   
void linspace(                      
              VECTOR&     x,           
        const Numeric     start,    
        const Numeric     stop,        
        const Numeric     step );
VECTOR linspace(             
        const Numeric     start, 
        const Numeric     stop,  
        const Numeric     step );
void nlinspace(         
              VECTOR&     x, 
        const Numeric     start,     
        const Numeric     stop,        
        const int         n );
VECTOR nlinspace(         
        const Numeric     start, 
        const Numeric     stop,  
        const int         n );
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
void to_matrix(MATRIX& W, const VECTOR& x);
MATRIX to_matrix(const VECTOR& x);



//
// Matrix to vector conversion:
//   to_vector
//
void to_vector(VECTOR& x, const MATRIX& W);
VECTOR to_vector(const MATRIX& W);


#endif  // math_funcs_h
