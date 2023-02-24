/* Copyright (C) 2002-2012 Claudia Emde <claudia.emde@dlr.de>

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
     \author Claudia Emde <claudia.emde@dlr.de>
     \date   Thu May  2 14:34:05 2002

     \brief  Linear algebra functions.

   */

#ifndef linalg_h
#define linalg_h

#include "array.h"
#include "matpack_complex.h"
#include "matpack_data.h"

// LU decomposition
void ludcmp(Matrix& LU, ArrayOfIndex& indx, ConstMatrixView A);

// LU backsubstitution
void lubacksub(VectorView x,
               ConstMatrixView LU,
               ConstVectorView b,
               const ArrayOfIndex& indx);

// Solve linear system
void solve(VectorView x, ConstMatrixView A, ConstVectorView b);

// Matrix inverse
void inv(MatrixView Ainv, ConstMatrixView A);

// Matrix inverse
void inv(ComplexMatrixView Ainv, const ConstComplexMatrixView A);

// Matrix diagonalization with lapack
void diagonalize(MatrixView P, VectorView WR, VectorView WI, ConstMatrixView A);

// Matrix diagonalization with lapack
void diagonalize(ComplexMatrixView P,
                 ComplexVectorView W,
                 const ConstComplexMatrixView A);

// Exponential of a Matrix
void matrix_exp(MatrixView F, ConstMatrixView A, const Index& q = 10);

// 2-norm of a vector
Numeric norm2(ConstVectorView v);

// Maximum absolute row sum norm
Numeric norm_inf(ConstMatrixView A);

// Identity Matrix
void id_mat(MatrixView I);

Numeric det(ConstMatrixView A);

void linreg(Vector& p, ConstVectorView x, ConstVectorView y);


/** Least squares fitting by solving x for known A and y
 * 
 * (A^T A)x = A^T y
 * 
 * Returns the squared residual, i.e., <(A^T A)x-A^T y|(A^T A)x-A^T y>.
 * 
 * @param[in]  x   As equation
 * @param[in]  A   As equation
 * @param[in]  y   As equation
 * @param[in]  residual (optional) Returns the residual if true
 * @return Squared residual or 0
 */
Numeric lsf(VectorView x, ConstMatrixView A, ConstVectorView y, bool residual=true) noexcept;

#endif  // linalg_h
