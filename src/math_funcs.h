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

#include "matpackI.h"
#include "mystring.h"

////////////////////////////////////////////////////////////////////////////
//// mean and standard deviation ////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void mean_row( VectorView m, ConstMatrixView x );

void std_row( VectorView s, ConstMatrixView x, ConstVectorView m );



////////////////////////////////////////////////////////////////////////////
//// first and last /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

Numeric first( ConstVectorView x );

Numeric last( ConstVectorView x );



////////////////////////////////////////////////////////////////////////////
//// Logical functions
////////////////////////////////////////////////////////////////////////////

bool any( const ArrayOfIndex& x ); 

bool isbool( const Index x );


////////////////////////////////////////////////////////////////////////////
// Functions to generate vectors
////////////////////////////////////////////////////////////////////////////

void linspace(                      
              Vector&     x,           
	      const Numeric  start,    
	      const Numeric  stop,        
	      const Numeric  step );

Vector linspace(             
		const Numeric  start, 
		const Numeric  stop,  
		const Numeric  step );

void nlinspace(         
              Vector&     x, 
	      const Numeric     start,     
	      const Numeric     stop,        
	      const Index       n );

Vector nlinspace(         
		 const Numeric     start, 
		 const Numeric     stop,  
		 const Index         n );

void nlogspace(         
	       Vector&     x, 
	       const Numeric     start,     
	       const Numeric     stop,        
	       const Index         n );

Vector nlogspace(  
		 const Numeric     start, 
		 const Numeric     stop,  
		 const Index         n );



/////////////////////////////////////////////////////////////////////////////
//   Random data
/////////////////////////////////////////////////////////////////////////////


void rand_uniform(
		  VectorView r,
		  const Numeric    x_low,
		  const Numeric    x_high );

void rand_gaussian(
		   VectorView r,
		   const Numeric    s );

void rand_matrix_uniform(
			 MatrixView m,
			 const Numeric   x_low,
			 const Numeric   x_high );

void rand_matrix_gaussian(
			  MatrixView r,
			  const Numeric    s );

void rand_data_gaussian(
			MatrixView         z,
			ConstVectorView    z0,
			ConstMatrixView    s );


////////////////////////////////////////////////////////////////////////////
//   Interpolation routines
////////////////////////////////////////////////////////////////////////////


void interp_lin_vector( VectorView       yi,
			ConstVectorView  x, 
			ConstVectorView  y, 
			ConstVectorView  xi );

void interp_lin_matrix(    
		       MatrixView        Yi,
		       ConstVectorView   x, 
		       ConstMatrixView   Y, 
		       ConstVectorView   xi );

Numeric interp_lin(         
		   ConstVectorView  x, 
		   ConstVectorView  y, 
		   const Numeric  xi );


/////////////////////////////////////////////////////////////////////////////
//   Factorization of matrices
/////////////////////////////////////////////////////////////////////////////

void chol(
	  MatrixView      r, 
	  ConstMatrixView c );

/////////////////////////////////////////////////////////////////////////////
//   Assert functions
/////////////////////////////////////////////////////////////////////////////

void assert_bool( const Index& x, const String& x_name );

void assert_lengths( const Vector& x1, const String& x1_name,
                     const Vector& x2, const String& x2_name );

void assert_length_nrow( const Vector& x, const String& x_name,
                         const Matrix& A, const String& A_name );

void assert_length_ncol( const Vector& x, const String& x_name,
                         const Matrix& A, const String& A_name );

#endif  // math_funcs_h
