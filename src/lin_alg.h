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


#include "matpackI.h"

// LU decomposition
void
ludcmp(MatrixView LU, 
       ArrayOfIndex& indx,
       ConstMatrixView A); 


// LU backsubstitution
void 
lubacksub(VectorView x, 
          ConstMatrixView LU,
          ConstVectorView b,
          const ArrayOfIndex& indx);


// Exponential of a Matrix
void 
matrix_exp(MatrixView F,
           ConstMatrixView A, 
           const Index& q);


// Maximum absolute row sum norm 
Numeric 
norm_inf(ConstMatrixView A);


// Identity Matrix
void
id_mat(MatrixView I);

#endif    // linalg_h
