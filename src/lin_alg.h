/* Copyright (C) 2002 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                      
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
   USA. 
*/

   /*!
     \file   lin_alg.h
     \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
     \date   Thu May  2 14:34:05 2002
     
     \brief  Linear algebra functions.
     
   */

   
#ifndef linalg_h
#define linalg_h

#include <stdexcept>
#include "arts.h"
#include "matpackI.h"
#include "make_vector.h"
#include "array.h"


//! Linear equation sytem solver.
/*! 
  Linear equation systems of the type Kx=y are solved using the LU-decomposition method. 
  (This function should also be used for calculating K^(-1)y, as it is more effective
  than calculating the inverse of K and then multipying with y.)

  \param x Output: Solution vector x.
  \param K Input: Matrix K.
  \param y Input: "Right-hand-side-vector" y
*/
void lusolve(VectorView x, ConstMatrixView K, ConstVectorView y);



//! LU decomposition.
/*!
  Given a matrix A this routine replaces it by the LU decompostion of a rowwise permutation of 
  itself. (Compare Numerical Recipies in C, pages 36-48.)
  
  \param LU Output: returns L and U in one matrix
  \param indx Output: Vector that records the row permutation.
  \param A Input: Matrix for which the LU decomposition is performed
*/
void ludcmp(MatrixView& LU, ArrayOfIndex& indx, MatrixView& A); 


//! LU backsubstitution
/*! 
  Solves a set of linear equations Ax=b. It is neccessairy to do a LU decomposition using 
  the function ludcp before using this function.
  
  \param x Output: Solution vector of the equation system. 
  \param LU Input: LU decomposition of the matrix (output of function ludcp).
  \param b  Input: Right-hand-side vector of equation system.
  \param indx Input: Pivoting information (output of function ludcp).
*/
void lubacksub(VectorView& x, MatrixView& LU, ConstVectorView& b, ArrayOfIndex& indx);


#endif    // linalg_h
