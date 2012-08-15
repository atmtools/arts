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
  \file   lin_alg.cc
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   Thu May  2 10:59:55 2002
  
  \brief  Linear algebra functions.
  
  This file contains mathematical tools to solve the vector radiative transfer 
  equation. 
*/

#include "lin_alg.h"

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include <cmath>
#include "arts.h"
#include "make_vector.h"
#include "array.h"
#include "logic.h"

//! LU decomposition.
/*!
  This function performes a LU Decomposition of the matrix A.
  (Compare Numerical Recipies in C, pages 36-48.)
  
  \param LU Output: returns L and U in one matrix
  \param indx Output: Vector that records the row permutation.
  \param A Input: Matrix for which the LU decomposition is performed
*/
void 
ludcmp(MatrixView LU,
       ArrayOfIndex& indx, 
       ConstMatrixView A) 
{
  Index imax = 0;
  const Numeric TINY=1.0e-20;
  Numeric big, dum, sum, temp, d;
 
  LU = A;
  const Index dim = A.nrows();
  assert(is_size(A,dim,dim)); //check, if A is quadratic

  Vector vv(dim);  //stores implicit scaling of each row
  d = 1.0;
  for (Index i=0; i<dim; i++)
    {
      indx[i]= i;
      big = 0.0;
      for (Index j=0; j<dim; j++)
        {
          if ((temp = abs(LU(i,j))) > big)
            big = temp;
        }
      if (big == 0.)
        throw runtime_error("ludcmp: Matrix is singular");
      vv[i] = 1.0/big; // save scaling
      }
  
  for (Index j=0; j<dim; j++)
    {
      for (Index i=0; i<j; i++)
        {
          sum = LU(i,j);
          for (Index k=0; k<i; k++) 
            sum -= LU(i,k)*LU(k,j);
          LU(i,j) = sum;
        }
      big = 0.0;
      for( Index i=j; i<dim; i++)
        {
          sum = LU(i,j);
          for (Index k=0; k<j; k++)
            sum -= LU(i,k)*LU(k,j);
          LU(i,j) = sum;
          if( (dum = vv[i]*fabs(sum)) >= big) 
            {
              big = dum;
              imax = i;
            }
        }
      if (j!=imax)
        {
          for(Index k=0; k<dim; k++)
            {
              dum = LU(imax,k);
              LU(imax,k) = LU(j,k);
              LU(j,k) = dum;
            }
          d = -d;
          vv[imax] = vv[j];
          indx[j] = imax;
          indx[imax] =j;
         }
     
      if(LU(j,j) == 0.0) LU(j,j) = TINY;
      
      if (j != dim) 
        {
          dum=1.0/LU(j,j);
          for (Index i=j+1; i<dim; i++)
            LU(i,j) *=dum;
        }
    }
}

 
 


//! LU backsubstitution
/*! 
  Solves a set of linear equations Ax=b. It is neccessairy to do a LU          
  decomposition using 
  the function ludcp before using this function.
  
  \param x Output: Solution vector of the equation system. 
  \param LU Input: LU decomposition of the matrix (output of function ludcp).
  \param b  Input: Right-hand-side vector of equation system.
  \param indx Input: Pivoting information (output of function ludcp).
*/
void 
lubacksub(VectorView x, 
          ConstMatrixView LU, 
          ConstVectorView b, 
          const ArrayOfIndex& indx)
{
  Index dim = LU.nrows(); 

  /* Check if the dimensions of the input matrix and vectors agree and if LU is a quadratic matrix.*/
  assert(is_size(LU, dim, dim));
  assert(is_size(b, dim));
  assert(is_size(indx, dim));

  
  for(Index i=0; i<dim; i++)
     {
       x[indx[i]] = b[i];
     }
 
  for (Index i=0; i<dim; i++)
    {
     Numeric sum = x[i];
     for (Index j=0; j<=i-1; j++)
       sum -= LU(i,j)*x[j]; 
      x[i] = sum;
     }

  for(Index i=dim-1; i>=0; i--)
    {
      Numeric sum = x[i];
      for (Index j=i+1; j<dim; j++)
        sum -= LU(i,j)*x[j];
      x[i] = sum/LU(i,i);
    }
}


//! Exponential of a Matrix
/*! 
  The exponential of a matrix is computed using the Pade-Approximation. 
  The method is decribed in:
  Golub, G. H. and C. F. Van Loan, Matrix Computation, p. 384, 
  Johns Hopkins University Press, 1983.
  
  \param F Output: The matrix exponential of A (Has to be initialized before
  calling the function.
  \param A Input: arbitrary square matrix
  \param q Input: Parameter for the accuracy of the computation
*/
void 
matrix_exp(MatrixView F,
           ConstMatrixView A,
           const Index& q)
{
  const Index n = A.ncols();

  /* Check if A and F are a quadratic and of the same dimension. */
  assert(is_size(A,n,n));
  assert(is_size(F,n,n));

  if ( is_diagonal(A) )
  {
      F = 0.;
      for(Index ii = 0; ii < n; ii++)
      {
          F(ii, ii) = exp(A(ii, ii));
      }
      return;
  }

  Numeric A_norm_inf, c;
  Numeric j;
  Matrix D(n,n), N(n,n), X(n,n), cX(n,n,0.0), B(n,n,0.0);
  Vector N_col_vec(n,0.), F_col_vec(n,0.);

  A_norm_inf = norm_inf(A);

  // This formular is derived in the book by Golub and Van Loan.
  j = 1 +  floor(1./log(2.)*log(A_norm_inf));
  
  if(j<0) j=0.;
  Index j_index = (Index)(j);

  // Scale matrix
  F = A;
  F /= pow(2,j);
 
  /* The higher q the more accurate is the computation, 
     see user guide for accuracy */
  //  Index q = 8;
  Numeric q_n = (Numeric)(q);
  id_mat(D);
  id_mat(N);
  id_mat(X);
  c = 1.;

  for(Index k=0; k<q; k++)
    {
      Numeric k_n = (Numeric)(k+1);
      c *= (q_n-k_n+1)/((2*q_n-k_n+1)*k_n);
      mult(B,F,X);              // X = F * X
      X = B;
      cX = X;                   
      cX *= c;                  // cX = X*c
      N += cX;                  // N = N + X*c
      cX *= pow(-1,k_n);        // cX = (-1)^k*c*X
      D += cX;                  // D = D + (-1)^k*c*X
    }

  /*Solve the equation system DF=N for F using LU decomposition,
   use the backsubstitution routine for columns of N*/ 

   /* Now use X  for the LU decomposition matrix of D.*/
  ArrayOfIndex indx(n);

  ludcmp(X, indx, D);

  for(Index i=0; i<n; i++)
    {
      N_col_vec = N(joker,i);  // extract column vectors of N
      lubacksub(F_col_vec, X, N_col_vec, indx);
      F(joker,i) = F_col_vec;  // construct F matrix  from column vectors 
    }
  
  /* The square of F gives the result. */
    for(Index k=0; k<j_index; k++)
    {
      mult(B,F,F);              // F = F^2
      F = B;
    }
 
}


//! Maximum absolute row sum norm 
/*! 
  This function returns the maximum absolute row sum norm of a 
  matrix A (see user guide for the definition).

  \param A Input: arbitrary matrix
  
  \return Maximum absolute row sum norm 
*/
Numeric norm_inf(ConstMatrixView A)
{
  Numeric norm_inf = 0;
  
  for(Index j=0; j<A.nrows(); j++)
    {
      Numeric row_sum = 0;
      //Calculate the row sum for all rows
      for(Index i=0; i<A.ncols(); i++)
        row_sum += abs(A(i,j));
      //Pick out the row with the highest row sum
      if( norm_inf < row_sum)
        norm_inf = row_sum;
    }
  return norm_inf;
}


//! Identity Matrix
/*! 
  \param I Output: identity matrix
*/
void 
id_mat(MatrixView I)
{

  const Index n = I.ncols();
  assert(n == I.nrows());
  
  I = 0;
  for(Index i=0; i<n; i++)
    I(i,i) = 1.;
}

/*!
    Determinant of N by N matrix. Simple recursive method.

    \param  A   In:    Matrix of size NxN.

    \author Richard Larsson
    \date   2012-08-03
*/
Numeric det(ConstMatrixView A)
{
    const Index dim = A.nrows();
    assert(dim == A.ncols());

    if(dim == 3)
        return A(0,0)*A(1,1)*A(2,2) + A(0,1)*A(1,2)*A(2,0) +
               A(0,2)*A(1,0)*A(2,1) - A(0,2)*A(1,1)*A(2,0) -
               A(0,1)*A(1,0)*A(2,2) - A(0,0)*A(1,2)*A(2,1);
    else if(dim == 2)
        return A(0,0) * A(1,1) - A(0,1) * A(1,0);
    else if(dim == 1)
        return A(0,0);

    Numeric ret_val = 0.;

    for(Index j = 0; j < dim; j++)
    {
        Matrix temp(dim-1,dim-1);
        for(Index I = 1; I < dim; I++)
            for(Index J = 0; J < dim; J++)
            {
                if(J < j)
                    temp(I-1,J) = A(I,J);
                else if(J > j)
                    temp(I-1,J-1) = A(I,J);
            }

        Numeric tempNum = det(temp);

        ret_val += ((j%2 == 0)?-1.:1.) * tempNum * A(0,j);
    }
    return ret_val;
}
