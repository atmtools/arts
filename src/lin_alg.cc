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
#include "lapack.h"
#include "matpackIII.h"
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
void ludcmp( Matrix &LU,
             ArrayOfIndex& indx,
             ConstMatrixView A )
{

  // Assert that A is quadratic.
  const Index n = A.nrows();
  assert( is_size( LU, n, n) );

  int n_int, info;
  int *ipiv = new int[n];

  LU = transpose(A);

  // Standard case: The arts matrix is not transposed, the leading
  // dimension is the row stride of the matrix.
  n_int = (int) n;

  // Compute LU decomposition using LAPACK dgetrf_.
  lapack::dgetrf_( &n_int, &n_int, LU.mdata, &n_int, ipiv, &info );

  // Copy pivot array to pivot vector.
  for ( Index i = 0; i < n; i++ )
  {
      indx[i] = ipiv[i];
  }
  delete [] ipiv;

}

//! LU backsubstitution
/*!
  Solves a set of linear equations Ax=b. It is neccessairy to do a L
  decomposition using the function ludcp before using this function. The
  backsubstitution is in-place, i.e. x and b may be the same vector.

  \param x Output: Solution vector of the equation system.
  \param LU Input: LU decomposition of the matrix (output of function ludcp).
  \param b  Input: Right-hand-side vector of equation system.
  \param indx Input: Pivoting information (output of function ludcp).
*/
void lubacksub( VectorView x,
                ConstMatrixView LU,
                ConstVectorView b,
                const ArrayOfIndex& indx )
{
  Index n = LU.nrows();

  /* Check if the dimensions of the input matrix and vectors agree and if LU
     is a quadratic matrix.*/
  Index column_stride = LU.mcr.get_stride();
  Index vec_stride = b.mrange.get_stride();

  assert(is_size(LU, n, n));
  assert(is_size(b, n));
  assert(is_size(indx, n));
  assert(column_stride == 1);
  assert(vec_stride == 1);

  char trans = 'N';
  int n_int = (int) n;
  int one = (int) 1;
  int ipiv[n];
  double rhs[n];
  int info;

  for ( Index i = 0; i < n; i++ )
  {
      ipiv[i] = (int) indx[i];
      rhs[i] = b[i];
  }

  lapack::dgetrs_( &trans,
                   &n_int,
                   &one,
                   LU.mdata,
                   &n_int,
                   ipiv,
                   rhs,
                   &n_int,
                   &info );

  for ( Index i = 0; i < n; i++ )
  {
      x[i] = rhs[i];
  }
}

//! Solve a linear system.
/*!
  Solves the linear system A*x = b for a general matrix A. For the solution of
  the system an additional n-times-n matrix and a size-n index vector are
  allocated.

  \param x The solution vector x.
  \param A The matrix A defining the system.
  \param b The vector b.
*/
void solve( VectorView x,
            ConstMatrixView A,
            ConstVectorView b )
{
    Index n = A.ncols();

    // Check dimensions of the system.
    assert( n == A.nrows() );
    assert( n == x.nelem() );
    assert( n == b.nelem() );

    // Allocate matrix and index vector for the LU decomposition.
    Matrix LU = Matrix(n,n);
    ArrayOfIndex indx(n);

    // Perform LU decomposition.
    ludcmp(LU, indx, A);

    // Solve the system using backsubstitution.
    lubacksub( x, LU, b, indx );
}

//! Matrix Inverse
/*!
  Compute the inverse of a matrix such that I = Ainv*A = A*Ainv. Both
  MatrixViews must be square and have the same size n. During the inversion one
  additional n times n Matrix is allocated and work space memory for faster
  inversion is allocated and freed.

  \param[out] Ainv The MatrixView to contain the inverse of A.
  \param[in] A The matrix to be inverted.
*/
void inv( MatrixView Ainv,
          ConstMatrixView A)
{
    // A must be a square matrix.
    assert(A.ncols() == A.nrows());

    Index n = A.ncols();
    Matrix LU( A );

    int n_int, info;
    int *ipiv = new int[n];

    n_int = (int) n;

    // Compute LU decomposition using LAPACK dgetrf_.
    lapack::dgetrf_( &n_int, &n_int, LU.mdata, &n_int, ipiv, &info );

    // Invert matrix.
    int lwork = n_int;
    double *work = new double[lwork];
    if (!work)
    {
        throw runtime_error( "Error inverting matrix: Could not allocate workspace memory." );
    }

    lapack::dgetri_( &n_int, LU.mdata, &n_int, ipiv, work, &lwork, &info );
    delete[] work;

    // Check for success.
    if (info == 0)
    {
        Ainv = LU;
    }
    else
    {
        throw runtime_error( "Error inverting matrix: Matrix not of full rank." );
    }
}

//! General exponential of a Matrix
/*!

  The exponential of a matrix is computed using the Pade-Approximation. The
  method is decribed in: Golub, G. H. and C. F. Van Loan, Matrix Computation,
  p. 384, Johns Hopkins University Press, 1983.

  The Pade-approximation is applied on all cases. If a faster option can be
  applied has to be checked before calling the function.

  \param F Output: The matrix exponential of A (Has to be initialized before
  calling the function.
  \param A Input: arbitrary square matrix
  \param q Input: Parameter for the accuracy of the computation
*/
void matrix_exp(
         MatrixView      F,
         ConstMatrixView A,
   const Index&          q )
{
  const Index n = A.ncols();

  /* Check if A and F are a quadratic and of the same dimension. */
  assert( is_size(A,n,n) );
  assert( is_size(F,n,n) );

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


//! General exponential of a Matrix with their derivatives
/*!
 *
 * The exponential of a matrix is computed using the Pade-Approximation. The
 * method is decribed in: Golub, G. H. and C. F. Van Loan, Matrix Computation,
 * p. 384, Johns Hopkins University Press, 1983.
 *
 * The extension of Pade-Approximation to the derivative is explained by
 * Lubomír Brančík, MATLAB PROGRAMS FOR MATRIX EXPONENTIAL FUNCTION
 * DERIVATIVE EVALUATION, 2008
 *
 * The Pade-approximation is applied on all cases. If a faster option can be
 * applied has to be checked before calling the function.
 *
 * \param F Output: The matrix exponential of A (Has to be initialized before
 * calling the function.
 * \param dF Output: The derivative  of the matrix exponential of A (Has to be initialized before
 * calling the function.  Page dimension is for each derivative.
 * \param A Input:  arbitrary square matrix.
 * \param dA Input: derivative of A.   Page dimension is for each derivative.
 * \param q Input: Parameter for the accuracy of the computation,  Matlab default is 6.
 */
void special_matrix_exp_and_dmatrix_exp_dx_for_rt(
    MatrixView           F,
    Tensor3View         dF_upp,
    Tensor3View         dF_low,
    ConstMatrixView      A,
    ConstTensor3View    dA_upp,
    ConstTensor3View    dA_low,
    const Index&         q )
{
    const Index n_partials = dA_upp.npages();

    const Index n = A.ncols();

    /* Check if A and F are a quadratic and of the same dimension. */
    assert( is_size(A,n,n) );
    assert( is_size(F,n,n) );
    assert(n_partials==dF_upp.npages());
    assert(n_partials==dF_low.npages());
    assert(n_partials==dA_low.npages());
    for(Index ii=0;ii<n_partials;ii++)
    {
        assert( is_size(dA_upp(ii,joker,joker),n,n) );
        assert( is_size(dA_low(ii,joker,joker),n,n) );
        assert( is_size(dF_upp(ii,joker,joker),n,n) );
        assert( is_size(dF_low(ii,joker,joker),n,n) );
    }

    // This sets up some cnstants
    Numeric A_norm_inf, e;
    A_norm_inf = norm_inf(A);
    e = 1. +  floor(1./log(2.)*log(A_norm_inf));
    Numeric r = (e+1.)>0.?(e+1.):0., pow2rm1=1./pow(2,r);
    Numeric c = 0.5;

    // For non-derivatives
    Matrix M=A, X(n,n), cX(n,n), D(n,n);
    M  *= pow2rm1;
    X   = M;
    cX  = X;
    cX*=c;
    id_mat(F); F+=cX;
    id_mat(D); D-=cX;

    // For derivatives
    Tensor3 dM_upp(n_partials,n,n),dM_low(n_partials,n,n),
    dD_upp(n_partials,n,n),dD_low(n_partials,n,n),
    Y_upp(n_partials,n,n),Y_low(n_partials,n,n),
    cY_upp(n_partials,n,n),cY_low(n_partials,n,n);
    for(Index ii=0;ii<n_partials;ii++)
    {
        for(Index jj=0;jj<n;jj++)
        {
            for(Index kk=0;kk<n;kk++)
            {
                dM_upp(ii,jj,kk)  = dA_upp(ii,jj,kk) * pow2rm1;
                dM_low(ii,jj,kk)  = dA_low(ii,jj,kk) * pow2rm1;

                Y_upp(ii,jj,kk)  = dM_upp(ii,jj,kk);
                Y_low(ii,jj,kk)  = dM_low(ii,jj,kk);

                cY_upp(ii,jj,kk)  = c*Y_upp(ii,jj,kk);
                cY_low(ii,jj,kk)  = c*Y_low(ii,jj,kk);

                dF_upp(ii,jj,kk)  = cY_upp(ii,jj,kk);
                dF_low(ii,jj,kk)  = cY_low(ii,jj,kk);

                dD_upp(ii,jj,kk)  =-cY_upp(ii,jj,kk);
                dD_low(ii,jj,kk)  =-cY_low(ii,jj,kk);
            }
        }
    }

    // NOTE: MATLAB paper sets q = 6 but we allow other numbers
    Matrix tmp1(n,n), tmp2(n,n),
    tmp_low1(n,n), tmp_upp1(n,n),
    tmp_low2(n,n), tmp_upp2(n,n);

    for(Index k=2; k<=q; k++)
    {
        c *= (Numeric)(q-k+1)/(Numeric)((k)*(2*q-k+1));

        // For partials
        for(Index ii=0;ii<n_partials;ii++)
        {
            // Y = dM*X + M*Y
            mult(tmp_upp1, dM_upp(ii,joker,joker), X);
            mult(tmp_upp2,  M, Y_upp(ii,joker,joker));
            mult(tmp_low1, dM_low(ii,joker,joker), X);
            mult(tmp_low2,  M, Y_low(ii,joker,joker));

            for(Index jj=0;jj<n;jj++)
            {
                for(Index kk=0;kk<n;kk++)
                {
                    Y_upp(ii,jj,kk) = tmp_upp1(jj,kk)+tmp_upp2(jj,kk);
                    Y_low(ii,jj,kk) = tmp_low1(jj,kk)+tmp_low2(jj,kk);

                    // cY = c*Y
                    cY_upp(ii,jj,kk) = c * Y_upp(ii,jj,kk);
                    cY_low(ii,jj,kk) = c * Y_low(ii,jj,kk);

                    // dF = dF + cY
                    dF_upp(ii,jj,kk) += cY_upp(ii,jj,kk);
                    dF_low(ii,jj,kk) += cY_low(ii,jj,kk);

                    if(k%2==0)//For even numbers, add.
                    {
                        // dD = dD + cY
                        dD_upp(ii,jj,kk)+=cY_upp(ii,jj,kk);
                        dD_low(ii,jj,kk)+=cY_low(ii,jj,kk);
                    }
                    else
                    {
                        // dD = dD - cY
                        dD_upp(ii,jj,kk)-=cY_upp(ii,jj,kk);
                        dD_low(ii,jj,kk)-=cY_low(ii,jj,kk);
                    }
                }
            }
        }

        //For full derivative (Note X in partials above)
        // X=M*X
        mult(tmp1, M, X);
        for(Index jj=0;jj<n;jj++)
        {
            for(Index kk=0;kk<n;kk++)
            {
                X(jj,kk) = tmp1(jj,kk);
                cX(jj,kk) = tmp1(jj,kk)*c;
                F(jj,kk) += cX(jj,kk);

                if(k%2==0)//For even numbers, D = D + cX
                    D(jj,kk)+=cX(jj,kk);
                else//For odd numbers, D = D - cX
                    D(jj,kk)-=cX(jj,kk);
            }
        }
    }

    // D^-1
    inv(tmp1,D);

    // F = D\F, or D^-1*F
    mult(tmp2,tmp1,F);
    F = tmp2;

    // For partials
    for(Index ii=0;ii<n_partials;ii++)
    {
        //dF = D \ (dF - dF*F), or D^-1 * (dF - dF*F)
        mult(tmp_upp1, dD_upp(ii,joker,joker), F);
        mult(tmp_low1, dD_low(ii,joker,joker), F);

        for(Index jj=0;jj<n;jj++)
        {
            for(Index kk=0;kk<n;kk++)
            {
                dF_upp(ii,jj,kk)-=tmp_upp1(jj,kk);// dF - dF * F
                dF_low(ii,jj,kk)-=tmp_low1(jj,kk);// dF - dF * F
            }
        }

        mult(tmp_upp2,tmp1,dF_upp(ii,joker,joker));
        mult(tmp_low2,tmp1,dF_low(ii,joker,joker));

        for(Index jj=0;jj<n;jj++)
        {
            for(Index kk=0;kk<n;kk++)
            {
                dF_upp(ii,jj,kk)=tmp_upp2(jj,kk);
                dF_low(ii,jj,kk)=tmp_low2(jj,kk);
            }
        }

    }

    for(Index k=1; k<=r; k++)
    {
        for(Index ii=0;ii<n_partials;ii++)
        {
            // dF=F*dF+dF*F
            mult(tmp_upp1,F,dF_upp(ii,joker,joker));//F*dF
            mult(tmp_upp2,dF_upp(ii,joker,joker),F);//dF*F
            mult(tmp_low1,F,dF_low(ii,joker,joker));//F*dF
            mult(tmp_low2,dF_low(ii,joker,joker),F);//dF*F

            for(Index jj=0;jj<n;jj++)
            {
                for(Index kk=0;kk<n;kk++)
                {
                    dF_upp(ii,jj,kk)=tmp_upp1(jj,kk)+tmp_upp2(jj,kk);
                    dF_low(ii,jj,kk)=tmp_low1(jj,kk)+tmp_low2(jj,kk);
                }
            }
        }

        // F=F*F
        mult(tmp1,F,F);
        F=tmp1;
    }
}


//! General exponential of a Matrix with their derivatives
/*!
 *
 * The exponential of a matrix is computed using the Pade-Approximation. The
 * method is decribed in: Golub, G. H. and C. F. Van Loan, Matrix Computation,
 * p. 384, Johns Hopkins University Press, 1983.
 *
 * The extension of Pade-Approximation to the derivative is explained by
 * Lubomír Brančík, MATLAB PROGRAMS FOR MATRIX EXPONENTIAL FUNCTION
 * DERIVATIVE EVALUATION, 2008
 *
 * The Pade-approximation is applied on all cases. If a faster option can be
 * applied has to be checked before calling the function.
 *
 * \param F Output: The matrix exponential of A (Has to be initialized before
 * calling the function.
 * \param dF Output: The derivative  of the matrix exponential of A (Has to be initialized before
 * calling the function.  Page dimension is for each derivative.
 * \param A Input:  arbitrary square matrix.
 * \param dA Input: derivative of A.   Page dimension is for each derivative.
 * \param q Input: Parameter for the accuracy of the computation,  Matlab default is 6.
 */
void matrix_exp_dmatrix_exp(
    MatrixView           F,
    Tensor3View         dF,
    ConstMatrixView      A,
    ConstTensor3View    dA,
    const Index&         q )
{
    const Index n_partials = dA.npages();

    const Index n = A.ncols();

    /* Check if A and F are a quadratic and of the same dimension. */
    assert( is_size(A,n,n) );
    assert( is_size(F,n,n) );
    assert(n_partials==dF.npages());
    for(Index ii=0;ii<n_partials;ii++)
    {
        assert( is_size(dA(ii,joker,joker),n,n) );
        assert( is_size(dF(ii,joker,joker),n,n) );
    }

    // This sets up some cnstants
    Numeric A_norm_inf, e;
    A_norm_inf = norm_inf(A);
    e = 1. +  floor(1./log(2.)*log(A_norm_inf));
    Numeric r = (e+1.)>0.?(e+1.):0., pow2rm1=1./pow(2,r);
    Numeric c = 0.5;

    // For non-derivatives
    Matrix M=A, X(n,n), cX(n,n), D(n,n);
    M  *= pow2rm1;
    X   = M;
    cX  = X;
    cX*=c;
    id_mat(F); F+=cX;
    id_mat(D); D-=cX;

    // For derivatives
    Tensor3 dM(n_partials,n,n),Y(n_partials,n,n),cY(n_partials,n,n), dD(n_partials,n,n);
    for(Index ii=0;ii<n_partials;ii++)
    {
        for(Index jj=0;jj<n;jj++)
        {
            for(Index kk=0;kk<n;kk++)
            {
                dM(ii,jj,kk)  = dA(ii,jj,kk) * pow2rm1;

                 Y(ii,jj,kk)  = dM(ii,jj,kk);

                cY(ii,jj,kk)  = c*Y(ii,jj,kk);

                dF(ii,jj,kk)  = cY(ii,jj,kk);

                dD(ii,jj,kk)  =-cY(ii,jj,kk);
            }
        }
    }

    // NOTE: MATLAB paper sets q = 6 but we allow other numbers

    Matrix tmp1(n,n), tmp2(n,n);

    for(Index k=2; k<=q; k++)
    {
        c *= (Numeric)(q-k+1)/(Numeric)((k)*(2*q-k+1));

        // For partials
        for(Index ii=0;ii<n_partials;ii++)
        {
            // Y = dM*X + M*Y
            mult(tmp1, dM(ii,joker,joker), X);
            mult(tmp2,  M,                 Y(ii,joker,joker));


            for(Index jj=0;jj<n;jj++)
                for(Index kk=0;kk<n;kk++)
                {
                    Y(ii,jj,kk) = tmp1(jj,kk)+tmp2(jj,kk);

                    // cY = c*Y
                    cY(ii,jj,kk) = c * Y(ii,jj,kk);

                    // dF = dF + cY
                    dF(ii,jj,kk) += cY(ii,jj,kk);

                    if(k%2==0)//For even numbers, add.
                    {
                        // dD = dD + cY
                        dD(ii,jj,kk)+=cY(ii,jj,kk);
                    }
                    else
                    {
                        // dD = dD - cY
                        dD(ii,jj,kk)-=cY(ii,jj,kk);
                    }
                }
        }

        //For full derivative (Note X in partials above)
        // X=M*X
        mult(tmp1, M, X);
        X = tmp1;

        //cX = c*X
        cX=X; cX*=c;

        // F = F + cX
        F  += cX;

        if(k%2==0)//For even numbers, D = D + cX
            D+=cX;
        else//For odd numbers, D = D - cX
            D-=cX;
    }

    // D^-1
    inv(tmp1,D);

    // F = D\F, or D^-1*F
    mult(tmp2,tmp1,F);
    F = tmp2;

    // For partials
    for(Index ii=0;ii<n_partials;ii++)
    {
        //dF = D \ (dF - dF*F), or D^-1 * (dF - dF*F)
        mult(tmp2, dD(ii,joker,joker), F);// dF * F
        dF(ii,joker,joker)-=tmp2;// dF - dF * F
        mult(tmp2,tmp1,dF(ii,joker,joker));
        dF(ii,joker,joker)=tmp2;
    }

    for(Index k=1; k<=r; k++)
    {
        for(Index ii=0;ii<n_partials;ii++)
        {
            // dF=F*dF+dF*F
            mult(tmp1,F,dF(ii,joker,joker));//F*dF
            mult(tmp2,dF(ii,joker,joker),F);//dF*F
            dF(ii,joker,joker)=tmp1;dF(ii,joker,joker)+=tmp2;
        }

        // F=F*F
        mult(tmp1,F,F);
        F=tmp1;
    }
}


//! General exponential of a Matrix with their derivatives
/*!
 *
 * The exponential of a matrix is computed using the Pade-Approximation. The
 * method is decribed in: Golub, G. H. and C. F. Van Loan, Matrix Computation,
 * p. 384, Johns Hopkins University Press, 1983.
 *
 * The extension of Pade-Approximation to the derivative is explained by
 * Lubomír Brančík, MATLAB PROGRAMS FOR MATRIX EXPONENTIAL FUNCTION
 * DERIVATIVE EVALUATION, 2008
 *
 * The Pade-approximation is applied on all cases. If a faster option can be
 * applied has to be checked before calling the function.
 *
 * \param F Output: The matrix exponential of A (Has to be initialized before
 * calling the function.
 * \param dF Output: The derivative  of the matrix exponential of A (Has to be initialized before
 * calling the function.
 * \param A Input:  arbitrary square matrix
 * \param dA Input: derivative of arbitrary square matrix
 * \param q Input: Parameter for the accuracy of the computation
 */
void matrix_exp_dmatrix_exp(
    MatrixView      F,
    MatrixView      dF,
    ConstMatrixView A,
    ConstMatrixView dA,
    const Index&          q )
{
    const Index n = A.ncols();

    /* Check if A and F are a quadratic and of the same dimension. */
    assert( is_size(A,n,n) );
    assert( is_size(F,n,n) );
    assert( is_size(dA,n,n) );
    assert( is_size(dF,n,n) );

    // This is the definition of how to scale
    Numeric A_norm_inf, e;
    A_norm_inf = norm_inf(A);
    e = 1. +  floor(1./log(2.)*log(A_norm_inf));
    Numeric r = (e+1.)>0.?(e+1.):0., pow2rm1=1./pow(2,r);

    // A and dA are scaled
    Matrix M=A, dM=dA;
    M  *= pow2rm1;
    dM *= pow2rm1;

    // These variables will hold the multiplication calculations
    Matrix X(n,n), Y(n,n);
    X =  M;
    Y = dM;

    // cX = c*M
    // cY = c*dM
    Matrix cX=X, cY=Y;
    Numeric c = 0.5;
    cX*=c;
    cY*=c;

    Matrix D(n,n), dD(n,n);
    // F = I + c*M
    id_mat(F); F+=cX;

    // dF = c*dM;
    dF =  cY;

    //D = I -c*M
    id_mat(D); D-=cX;

    // dD = -c*dM
    dD =  cY; dD*=-1.;

    // NOTE: MATLAB paper sets q = 6 but we allow other numbers

    Matrix tmp1(n,n), tmp2(n,n);

    for(Index k=2; k<=q; k++)
    {
        c *= (Numeric)(q-k+1)/(Numeric)((k)*(2*q-k+1));

        // Y = dM*X + M*Y
        mult(tmp1, dM, X);
        mult(tmp2, M, Y);
        Y = tmp1; Y+= tmp2;

        // X=M*X
        mult(tmp1, M, X);
        X = tmp1;

        //cX = c*X
        cX=X; cX*=c;

        // cY = c*Y
        cY=Y; cY*=c;

        // F = F + cX
        F  += cX;

        // dF = dF + cY
        dF += cY;

        if(k%2==0)//For even numbers, add.
        {
            // D = D + cX
            D+=cX;

            // dD = dD + cY
            dD+=cY;
        }
        else//For odd numbers, subtract
        {
            // D = D - cX
            D-=cX;

            // dD = dD - cY
            dD-=cY;
        }
    }

    // D^-1
    inv(tmp1,D);

    // F = D\F, or D^-1*F
    mult(tmp2,tmp1,F);
    F = tmp2;

    //dF = D \ (dF - dF*F), or D^-1 * (dF - dF*F)
    mult(tmp2, dD, F);// dF * F
    dF-=tmp2;// dF - dF * F
    mult(tmp2,tmp1,dF);
    dF=tmp2;

    for(Index k=1; k<=r; k++)
    {
        // dF=F*dF+dF*F
        mult(tmp1,F,dF);//F*dF
        mult(tmp2,dF,F);//dF*F
        dF=tmp1;dF+=tmp2;

        // F=F*F
        mult(tmp1,F,F);
        F=tmp1;
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
void id_mat(MatrixView I)
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



/*!
    Determines coefficients for linear regression

    Performs a least squares estimation of the model

       y = p[0] + p[1] * x

    \param  p   Out: Fitted coefficients.
    \param  x   In: x-value of data points
    \param  y   In: y-value of data points

    \author Patrick Eriksson
    \date   2013-01-25
*/
void linreg(
       Vector&    p,
  ConstVectorView x,
  ConstVectorView y )
{
  const Index n = x.nelem();

  assert( y.nelem() == n );

  p.resize(2);

  // Basic algorithm found at e.g.
  // http://en.wikipedia.org/wiki/Simple_linear_regression
  // The basic algorithm is as follows:
  /*
  Numeric s1=0, s2=0, s3=0, s4=0;
  for( Index i=0; i<n; i++ )
    {
      s1 += x[i] * y[i];
      s2 += x[i];
      s3 += y[i];
      s4 += x[i] * x[i];
    }

  p[1] = ( s1 - (s2*s3)/n ) / ( s4 - s2*s2/n );
  p[0] = s3/n - p[1]*s2/n;
  */

  // A version abit more numerical stable:
  // Mean value of x is removed before the fit: x' = (x-mean(x))
  // This corresponds to that s2 in version above becomes 0
  // y = a + b*x'
  // p[1] = b
  // p[0] = a - p[1]*mean(x)

  Numeric s1=0, xm=0, s3=0, s4=0;

  for( Index i=0; i<n; i++ )
    { xm += x[i]/Numeric(n); }

  for( Index i=0; i<n; i++ )
    {
      const Numeric xv = x[i] - xm;
      s1 += xv * y[i];
      s3 += y[i];
      s4 += xv * xv;
    }

  p[1] = s1 / s4;
  p[0] = s3/Numeric(n) - p[1]*xm;
}
