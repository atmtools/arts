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
   \file   math_funcs.h

   Contains the decleration of the functions in math_funcs.cc.

   \author Patrick Eriksson
   \date 2000-09-18 
*/



#ifndef math_funcs_h
#define math_funcs_h


#include <time.h>
#include <math.h>
#include <stdexcept>
#include "array.h"
#include "arts.h"
#include "matpackI.h"
#include "mystring.h"



Index is_bool( 
       const Index&    x );

Numeric last( ConstVectorView x );

Index is_sorted( 
       ConstVectorView&    x );

Index is_increasing( 
        ConstVectorView&    x );

Index is_decreasing( 
        ConstVectorView&    x );

Index last( const ArrayOfIndex& x );




//
// Old functions:
//



////////////////////////////////////////////////////////////////////////////
//// Logical functions
////////////////////////////////////////////////////////////////////////////

bool any( const ArrayOfIndex& x ); 



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
//   Check of function input
/////////////////////////////////////////////////////////////////////////////

void check_if_bool( const Index& x, const String& x_name );

void check_if_in_range( 
   const Numeric& x_low, 
   const Numeric& x_high, 
   const Numeric& x, 
   const String&  x_name );

void check_lengths( const Vector& x1, const String& x1_name,
                     const Vector& x2, const String& x2_name );

void check_length_nrow( const Vector& x, const String& x_name,
                         const Matrix& A, const String& A_name );

void check_length_ncol( const Vector& x, const String& x_name,
                         const Matrix& A, const String& A_name );

void check_ncol_nrow( const Matrix& A, const String& A_name,
		      const Matrix& B, const String& B_name );

#endif  // math_funcs_h
