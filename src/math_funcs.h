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

//// sqrt //////////////////////////////////////////////////////////////////

void sqrt( VECTOR& y, const VECTOR& x );

VECTOR sqrt( const VECTOR& x );



//// exp ///////////////////////////////////////////////////////////////////

void exp( VECTOR& y, const VECTOR& x );

VECTOR exp( const VECTOR& x );

void exp( MATRIX& Y, const MATRIX& X );

MATRIX exp( const MATRIX& X );



//// log ////////////////////////////////////////////////////////////////////

void log( VECTOR& y, const VECTOR& x );

VECTOR log( const VECTOR& x );

void log( MATRIX& Y, const MATRIX& X );

MATRIX log( const MATRIX& X );



//// mean and standard deviation ////////////////////////////////////////////

void mean_row( VECTOR& m, const MATRIX& x );

void std_row( VECTOR& s, const MATRIX& x, const VECTOR& m  );



//// min and max ////////////////////////////////////////////////////////////

Numeric min( const VECTOR& x );

Numeric min( const MATRIX& A );

Numeric max( const VECTOR& x );

Numeric max( const MATRIX& A );

// Max and min are not needed for array, due to MTL builtin functions.

/* Gives the maximum value of an array.

    Because this is a template function, the definition has to be also 
    in the header file, and not in file math_func.cc.

    \return      the maximum value of x
    \param   x   an array

    \author Stefan Buehler
    \date   2000-06-27
*/
// template<class T>
// T max( const ARRAY<T>& x )
// {
//   size_t n = x.size();
//   T y=x[0];
//   for ( size_t i=1; i<n; i++ )
//     {
//       if ( x[i] > y )
// 	y = x[i];
//     }
//   return y; 
// }

/* Gives the minimum value of an array.

    Because this is a template function, the definition has to be also 
    in the header file, and not in file math_func.cc.

    \return      the minimum value of x
    \param   x   an array

    \author Stefan Buehler
    \date   2000-06-27
*/
// template<class T>
// T min( const ARRAY<T>& x )
// {
//   size_t n = x.size();
//   T y=x(1);
//   for ( size_t i=2; i<=n; i++ )
//   {
//     if ( x(i) < y )
//       y = x(i);
//   }
//   return y; 
// }



//// first and last /////////////////////////////////////////////////////////

Numeric first( const VECTOR& x );

Numeric last( const VECTOR& x );



////////////////////////////////////////////////////////////////////////////
//   Logical functions
////////////////////////////////////////////////////////////////////////////

bool any( const ARRAYofsizet& x ); 



////////////////////////////////////////////////////////////////////////////
// Functions to generate vectors
////////////////////////////////////////////////////////////////////////////

void linspace(                      
              VECTOR&     x,           
        const Numeric  start,    
        const Numeric  stop,        
        const Numeric  step );

VECTOR linspace(             
        const Numeric  start, 
        const Numeric  stop,  
        const Numeric  step );

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



////////////////////////////////////////////////////////////////////////////
//   Interpolation routines
////////////////////////////////////////////////////////////////////////////

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

void interp_lin_col(    
              MATRIX&  Yi,
        const VECTOR&  x, 
        const MATRIX&  Y, 
        const VECTOR&  xi );

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

void to_matrix(MATRIX& W, const VECTOR& x);

MATRIX to_matrix(const VECTOR& x);



void to_vector(VECTOR& x, const MATRIX& W);

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

void row(VECTOR& x,
	 size_t i,
	 const MATRIX& A);

VECTOR row(size_t i,
	   const MATRIX& A);

void col(VECTOR& x,
	 size_t i,
	 const MATRIX& A);

VECTOR col(size_t i,
	   const MATRIX& A);

void row(MATRIX& X,
	 size_t i,
	 size_t k,
	 const MATRIX& A);

MATRIX row(size_t i,
	   size_t k,
	   const MATRIX& A);

void col(MATRIX& X,
	 size_t i,
	 size_t k,
	 const MATRIX& A);

MATRIX col(size_t i,
	   size_t k,
	   const MATRIX& A);



/////////////////////////////////////////////////////////////////////////////
//   Putting data in a matrix column or and row
/////////////////////////////////////////////////////////////////////////////

void put_in_col(
              MATRIX& A,
	      size_t  i, 
        const VECTOR& x );



/////////////////////////////////////////////////////////////////////////////
//   Random data
/////////////////////////////////////////////////////////////////////////////

void rand_uniform(
              VECTOR&    r,
        const size_t     n,
        const Numeric&   x_low,
        const Numeric&   x_high );

void rand_gaussian(
              VECTOR&    r,
        const size_t     n,
        const Numeric&   s );

void rand_matrix_uniform(
              MATRIX&    m,
        const size_t&    nrows,
        const size_t&    ncols,
        const Numeric&   x_low,
        const Numeric&   x_high );

void rand_matrix_gaussian(
              MATRIX&    r,
        const size_t&    nrows,
        const size_t&    ncols,
        const Numeric&   s );

void rand_data_gaussian(
              MATRIX&    z,
        const size_t&    n,
        const VECTOR&    z0,
        const MATRIX&    s );

#endif  // math_funcs_h
