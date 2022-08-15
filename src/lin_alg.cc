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

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <new>
#include <stdexcept>

#include "array.h"
#include "arts.h"
#include "arts_omp.h"
#include "lapack.h"
#include "logic.h"
#include "matpackIII.h"

//! LU decomposition.
/*!
  This function performes a LU Decomposition of the matrix A.
  (Compare Numerical Recipies in C, pages 36-48.)

  \param LU Output: returns L and U in one matrix
  \param indx Output: Vector that records the row permutation.
  \param A Input: Matrix for which the LU decomposition is performed
*/
void ludcmp(Matrix& LU, ArrayOfIndex& indx, ConstMatrixView A) {
  // Assert that A is quadratic.
  const Index n = A.nrows();
  ARTS_ASSERT(is_size(LU, n, n));

  int n_int, info;
  int* ipiv = new int[n];

  LU = transpose(A);

  // Standard case: The arts matrix is not transposed, the leading
  // dimension is the row stride of the matrix.
  n_int = (int)n;

  // Compute LU decomposition using LAPACK dgetrf_.
  lapack::dgetrf_(&n_int, &n_int, LU.mdata, &n_int, ipiv, &info);

  // Copy pivot array to pivot vector.
  for (Index i = 0; i < n; i++) {
    indx[i] = ipiv[i];
  }
  delete[] ipiv;
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
void lubacksub(VectorView x,
               ConstMatrixView LU,
               ConstVectorView b,
               const ArrayOfIndex& indx) {
  Index n = LU.nrows();

  /* Check if the dimensions of the input matrix and vectors agree and if LU
     is a quadratic matrix.*/
  DEBUG_ONLY(Index column_stride = LU.mcr.get_stride());
  DEBUG_ONLY(Index vec_stride = b.mrange.get_stride());

  ARTS_ASSERT(is_size(LU, n, n));
  ARTS_ASSERT(is_size(b, n));
  ARTS_ASSERT(is_size(indx, n));
  ARTS_ASSERT(column_stride == 1);
  ARTS_ASSERT(vec_stride == 1);

  char trans = 'N';
  int n_int = (int)n;
  int one = (int)1;
  std::vector<int> ipiv(n);
  std::vector<double> rhs(n);
  int info;

  for (Index i = 0; i < n; i++) {
    ipiv[i] = (int)indx[i];
    rhs[i] = b[i];
  }

  lapack::dgetrs_(
      &trans, &n_int, &one, LU.mdata, &n_int, ipiv.data(), rhs.data(), &n_int, &info);

  for (Index i = 0; i < n; i++) {
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
void solve(VectorView x, ConstMatrixView A, ConstVectorView b) {
  Index n = A.ncols();

  // Check dimensions of the system.
  ARTS_ASSERT(n == A.nrows());
  ARTS_ASSERT(n == x.nelem());
  ARTS_ASSERT(n == b.nelem());

  // Allocate matrix and index vector for the LU decomposition.
  Matrix LU = Matrix(n, n);
  ArrayOfIndex indx(n);

  // Perform LU decomposition.
  ludcmp(LU, indx, A);

  // Solve the system using backsubstitution.
  lubacksub(x, LU, b, indx);
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
void inv(MatrixView Ainv, ConstMatrixView A) {
  // A must be a square matrix.
  ARTS_ASSERT(A.ncols() == A.nrows());

  Index n = A.ncols();
  Matrix LU(A);

  int n_int, info;
  int* ipiv = new int[n];

  n_int = (int)n;

  // Compute LU decomposition using LAPACK dgetrf_.
  lapack::dgetrf_(&n_int, &n_int, LU.mdata, &n_int, ipiv, &info);

  // Invert matrix.
  int lwork = n_int;
  double* work;

  try {
    work = new double[lwork];
  } catch (std::bad_alloc& ba) {
    ARTS_USER_ERROR(
        "Error inverting matrix: Could not allocate workspace memory.");
  }

  lapack::dgetri_(&n_int, LU.mdata, &n_int, ipiv, work, &lwork, &info);
  delete[] work;
  delete[] ipiv;

  // Check for success.
  ARTS_USER_ERROR_IF(info not_eq 0,
                     "Error inverting matrix: Matrix not of full rank.");
  Ainv = LU;
}

void inv(ComplexMatrixView Ainv, const ConstComplexMatrixView A) {
  // A must be a square matrix.
  ARTS_ASSERT(A.ncols() == A.nrows());

  Index n = A.ncols();

  // Workdata
  Ainv = A;
  int n_int = int(n), lwork = n_int, info;
  std::vector<int> ipiv(n);
  ComplexVector work(lwork);

  // Compute LU decomposition using LAPACK dgetrf_.
  lapack::zgetrf_(
      &n_int, &n_int, Ainv.get_c_array(), &n_int, ipiv.data(), &info);
  lapack::zgetri_(&n_int,
                  Ainv.get_c_array(),
                  &n_int,
                  ipiv.data(),
                  work.get_c_array(),
                  &lwork,
                  &info);

  // Check for success.
  ARTS_USER_ERROR_IF(info not_eq 0,
                     "Error inverting matrix: Matrix not of full rank.");
}

//! Matrix Diagonalization
/*!
 * Return P and W from A in the statement diag(P^-1*A*P)-W == 0.
 * The real function will require some manipulation if 
 * the eigenvalues are imaginary.
 * 
 * The real version returns WR and WI as returned by dgeev. 
 * The complex version just returns W.
 * 
 * The function makes many copies and is thereby not fast.
 * There are no tests on the condition of the returned matrix,
 * so nan and inf can occur.
 * 
 * \param[out] P The right eigenvectors.
 * \param[out] WR The real eigenvalues.
 * \param[out] WI The imaginary eigenvalues.
 * \param[in]  A The matrix to diagonalize.
 */
void diagonalize(MatrixView P,
                 VectorView WR,
                 VectorView WI,
                 ConstMatrixView A) {
  Index n = A.ncols();

  // A must be a square matrix.
  ARTS_ASSERT(n == A.nrows());
  ARTS_ASSERT(n == WR.nelem());
  ARTS_ASSERT(n == WI.nelem());
  ARTS_ASSERT(n == P.nrows());
  ARTS_ASSERT(n == P.ncols());

  Matrix A_tmp = A;
  Matrix P2 = P;
  Vector WR2 = WR;
  Vector WI2 = WI;

  // Integers
  int LDA, LDA_L, LDA_R, n_int, info;
  n_int = (int)n;
  LDA = LDA_L = LDA_R = (int)A.mcr.get_extent();

  // We want to calculate RP not LP
  char l_eig = 'N', r_eig = 'V';

  // Work matrix
  int lwork = 2 * n_int + n_int * n_int;
  double *work, *rwork;
  try {
    rwork = new double[2 * n_int];
    work = new double[lwork];
  } catch (std::bad_alloc& ba) {
    ARTS_USER_ERROR(
        "Error diagonalizing: Could not allocate workspace memory.");
  }

  // Memory references
  double* adata = A_tmp.mdata;
  double* rpdata = P2.mdata;
  auto* lpdata = new double[0];  //To not confuse the compiler
  double* wrdata = WR2.mdata;
  double* widata = WI2.mdata;

  // Main calculations.  Note that errors in the output is ignored
  lapack::dgeev_(&l_eig,
                 &r_eig,
                 &n_int,
                 adata,
                 &LDA,
                 wrdata,
                 widata,
                 lpdata,
                 &LDA_L,
                 rpdata,
                 &LDA_R,
                 work,
                 &lwork,
                 rwork,
                 &info);

  // Free memory.  Can these be sent in to speed things up?
  delete[] work;
  delete[] rwork;
  delete[] lpdata;

  // Re-order.  This can be done better?
  for (Index i = 0; i < n; i++)
    for (Index j = 0; j < n; j++) P(j, i) = P2(i, j);
  WI = WI2;
  WR = WR2;
}

//! Matrix Diagonalization
/*!
 * Return P and W from A in the statement diag(P^-1*A*P)-W == 0.
 * 
 * The function makes many copies and is thereby not fast.
 * There are no tests on the condition of the returned matrix,
 * so nan and inf can occur.
 * 
 * \param[out] P The right eigenvectors.
 * \param[out] W The eigenvalues.
 * \param[in]  A The matrix to diagonalize.
 */
void diagonalize(ComplexMatrixView P,
                 ComplexVectorView W,
                 const ConstComplexMatrixView A) {
  Index n = A.ncols();

  // A must be a square matrix.
  ARTS_ASSERT(n == A.nrows());
  ARTS_ASSERT(n == W.nelem());
  ARTS_ASSERT(n == P.nrows());
  ARTS_ASSERT(n == P.ncols());

  ComplexMatrix A_tmp = A;

  // Integers
  int LDA = int(A.get_column_extent()), LDA_L = int(A.get_column_extent()),
      LDA_R = int(A.get_column_extent()), n_int = int(n), info;

  // We want to calculate RP not LP
  char l_eig = 'N', r_eig = 'V';

  // Work matrix
  int lwork = 2 * n_int + n_int * n_int;
  ComplexVector work(lwork);
  ComplexVector lpdata(0);
  Vector rwork(2 * n_int);

  // Main calculations.  Note that errors in the output is ignored
  lapack::zgeev_(&l_eig,
                 &r_eig,
                 &n_int,
                 A_tmp.get_c_array(),
                 &LDA,
                 W.get_c_array(),
                 lpdata.get_c_array(),
                 &LDA_L,
                 P.get_c_array(),
                 &LDA_R,
                 work.get_c_array(),
                 &lwork,
                 rwork.get_c_array(),
                 &info);

  for (Index i = 0; i < n; i++)
    for (Index j = 0; j <= i; j++) std::swap(P(j, i), P(i, j));
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
void matrix_exp(MatrixView F, ConstMatrixView A, const Index& q) {
  const Index n = A.ncols();

  /* Check if A and F are a quadratic and of the same dimension. */
  ARTS_ASSERT(is_size(A, n, n));
  ARTS_ASSERT(is_size(F, n, n));

  Numeric A_norm_inf, c;
  Numeric j;
  Matrix D(n, n), N(n, n), X(n, n), cX(n, n, 0.0), B(n, n, 0.0);
  Vector N_col_vec(n, 0.), F_col_vec(n, 0.);

  A_norm_inf = norm_inf(A);

  // This formular is derived in the book by Golub and Van Loan.
  j = 1 + floor(1. / log(2.) * log(A_norm_inf));

  if (j < 0) j = 0.;
  auto j_index = (Index)(j);

  // Scale matrix
  F = A;
  F /= pow(2, j);

  /* The higher q the more accurate is the computation,
     see user guide for accuracy */
  //  Index q = 8;
  auto q_n = (Numeric)(q);
  id_mat(D);
  id_mat(N);
  id_mat(X);
  c = 1.;

  for (Index k = 0; k < q; k++) {
    auto k_n = (Numeric)(k + 1);
    c *= (q_n - k_n + 1) / ((2 * q_n - k_n + 1) * k_n);
    mult(B, F, X);  // X = F * X
    X = B;
    cX = X;
    cX *= c;             // cX = X*c
    N += cX;             // N = N + X*c
    cX *= pow(-1, k_n);  // cX = (-1)^k*c*X
    D += cX;             // D = D + (-1)^k*c*X
  }

  /*Solve the equation system DF=N for F using LU decomposition,
   use the backsubstitution routine for columns of N*/

  /* Now use X  for the LU decomposition matrix of D.*/
  ArrayOfIndex indx(n);

  ludcmp(X, indx, D);

  for (Index i = 0; i < n; i++) {
    N_col_vec = N(joker, i);  // extract column vectors of N
    lubacksub(F_col_vec, X, N_col_vec, indx);
    F(joker, i) = F_col_vec;  // construct F matrix  from column vectors
  }

  /* The square of F gives the result. */
  for (Index k = 0; k < j_index; k++) {
    mult(B, F, F);  // F = F^2
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
Numeric norm_inf(ConstMatrixView A) {
  Numeric norm_inf = 0;

  for (Index j = 0; j < A.nrows(); j++) {
    Numeric row_sum = 0;
    //Calculate the row sum for all rows
    for (Index i = 0; i < A.ncols(); i++) row_sum += abs(A(i, j));
    //Pick out the row with the highest row sum
    if (norm_inf < row_sum) norm_inf = row_sum;
  }
  return norm_inf;
}

//! Identity Matrix
/*!
  \param I Output: identity matrix
*/
void id_mat(MatrixView I) {
  const Index n = I.ncols();
  ARTS_ASSERT(n == I.nrows());

  I = 0;
  for (Index i = 0; i < n; i++) I(i, i) = 1.;
}

/*!
    Determinant of N by N matrix. Simple recursive method.

    \param  A   In:    Matrix of size NxN.

    \author Richard Larsson
    \date   2012-08-03
*/
Numeric det(ConstMatrixView A) {
  const Index dim = A.nrows();
  ARTS_ASSERT(dim == A.ncols());

  if (dim == 3)
    return A(0, 0) * A(1, 1) * A(2, 2) + A(0, 1) * A(1, 2) * A(2, 0) +
           A(0, 2) * A(1, 0) * A(2, 1) - A(0, 2) * A(1, 1) * A(2, 0) -
           A(0, 1) * A(1, 0) * A(2, 2) - A(0, 0) * A(1, 2) * A(2, 1);
  if (dim == 2) return A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
  if (dim == 1) return A(0, 0);

  Numeric ret_val = 0.;

  for (Index j = 0; j < dim; j++) {
    Matrix temp(dim - 1, dim - 1);
    for (Index I = 1; I < dim; I++)
      for (Index J = 0; J < dim; J++) {
        if (J < j)
          temp(I - 1, J) = A(I, J);
        else if (J > j)
          temp(I - 1, J - 1) = A(I, J);
      }

    Numeric tempNum = det(temp);

    ret_val += ((j % 2 == 0) ? -1. : 1.) * tempNum * A(0, j);
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
void linreg(Vector& p, ConstVectorView x, ConstVectorView y) {
  const Index n = x.nelem();

  ARTS_ASSERT(y.nelem() == n);

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

  Numeric s1 = 0, xm = 0, s3 = 0, s4 = 0;

  for (Index i = 0; i < n; i++) {
    xm += x[i] / Numeric(n);
  }

  for (Index i = 0; i < n; i++) {
    const Numeric xv = x[i] - xm;
    s1 += xv * y[i];
    s3 += y[i];
    s4 += xv * xv;
  }

  p[1] = s1 / s4;
  p[0] = s3 / Numeric(n) - p[1] * xm;
}

Numeric lsf(VectorView x,
            ConstMatrixView A,
            ConstVectorView y,
            bool residual) noexcept {
  // Size of the problem
  const Index n = x.nelem();
  Matrix AT, ATA(n, n);
  Vector ATy(n);

  // Solver
  AT = transpose(A);
  mult(ATA, AT, A);
  mult(ATy, AT, y);
  solve(x, ATA, ATy);

  // Residual
  if (residual) {
    Vector r(n);
    mult(r, ATA, x);
    r -= ATy;
    return r * r;
  }

  return 0;
}
