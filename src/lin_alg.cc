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
  \file   lin_alg.cc
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
  \date   Thu May  2 10:59:55 2002
  
  \brief  Linear algebra functions.
  
  This file contains mathematical tools to solve the vector radiative transfer 
  equation. 
*/

#include "lin_alg.h"


//! Linear equation sytem solver.
/*! 
  Linear equations of the type Kx=y are solved using the LU-decomposition method. 
  (This function should also be used for calculating K^(-1)y, as it is more effective
  than calculating the inverse of K and then multipying with y.)

  \param x Output: Solution vector x.
  \param K Input: Matrix K.
  \param y Input: "Right-hand-side-vector" y
*/
void lusolve(VectorView x, ConstMatrixView K, ConstVectorView y)
{

  /* For the LU decomposition K has to be quadratic and the "right hand side vector" x   
     has to be of the same dimension. */
  assert(K.nrows() == K.ncols());
  assert(K.nrows() == y.nelem());

  static Matrix A;
  if (A.nrows() != K.nrows() || A.ncols() != K.ncols())
    A.resize(K.nrows(), K.ncols());

  static Matrix LU;
  if (LU.nrows() != K.nrows() || LU.ncols() != K.ncols())
    LU.resize(K.nrows(), K.ncols());

  static ArrayOfIndex indx;
  if (K.nrows() != indx.nelem())
    indx.resize(K.nrows());

  /*Copy K to the variable used in the functions ludcmp.*/
  A = K;
  

  /* Perform a LU-decomposition.*/
  ludcmp(LU, indx, A);
  
  /* Perform the backsubstitution to solve the equation system. */
  lubacksub(x, LU, y, indx);

}






//! LU decomposition.
/*!
  Given a matrix A this routine replaces it by the LU decompostion of a rowwise permutation of 
  itself. (Compare Numerical Recipies in C, pages 36-48.)
  
  \param LU Output: returns L and U in one matrix
  \param indx Output: Vector that records the row permutation.
  \param A Input: Matrix for which the LU decomposition is performed
*/
void ludcmp(MatrixView& LU, ArrayOfIndex& indx, MatrixView& A) 
{
  int imax, dim;
  const Numeric TINY=1.0e-20;
  Numeric big, dum, sum, temp, d;
 
  LU = A;
  dim = A.nrows();
  Vector vv(dim);  //stores implicit scaling of each row
  d = 1.0;
  for (Index i=0; i<dim; i++)
    {
      indx[i]= i;
      big = 0.0;
      for (Index j=0; j<dim; j++)
        if ((temp = fabs(LU(i,j))) > big) big = temp;
      if (big == 0)
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
      //  indx[j] = j;

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

  vv[Range(0,dim)] = 0.0;
}

 
 


//! LU backsubstitution
/*! 
  Solves a set of linear equations Ax=b. It is neccessairy to do a LU decomposition using 
  the function ludcp before using this function.
  
  \param x Output: Solution vector of the equation system. 
  \param LU Input: LU decomposition of the matrix (output of function ludcp).
  \param b  Input: Right-hand-side vector of equation system.
  \param indx Input: Pivoting information (output of function ludcp).
*/
void lubacksub(VectorView& x, MatrixView& LU, ConstVectorView& b, ArrayOfIndex& indx)
{
  int ip,dim;
  float sum;
 
  dim = LU.nrows(); 
  
  for(Index i=0; i<dim; i++)
     {
       ip = indx[i];
       x[ip] = b[i];
     }
 
  for (Index i=0; i<dim; i++)
    {
      sum = x[i];
     for (Index j=0; j<=i-1; j++)
       sum -= LU(i,j)*x[j]; 
      x[i] = sum;
     }

  for(Index i=dim-1; i>=0; i--)
    {
      sum = x[i];
      for (Index j=i+1; j<dim; j++)
	sum -= LU(i,j)*x[j];
      x[i] = sum/LU(i,i);
    }
}

