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

#include <array.h>
#include <debug.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string_view>
#include <vector>

#include "lapack.h"
#include "matpack_mdspan_helpers_eigen.h"
#include "matpack_mdspan_helpers_matrix.h"
#include "matpack_mdspan_helpers_reduce.h"

namespace {

/** Check the status returned by a LAPACK LU factorization or inversion.
 *
 * A negative status identifies an invalid LAPACK argument.  A positive
 * status identifies the one-based diagonal element of U that is exactly
 * zero, so the matrix is singular.
 */
void check_lu_info(const int info, const std::string_view routine) {
  ARTS_USER_ERROR_IF(info < 0, "{} received an illegal value for argument {}.", routine, -info);
  ARTS_USER_ERROR_IF(info > 0, "{} found a singular matrix: U({}, {}) is exactly zero.", routine, info, info);
}

/** Check the status returned by a LAPACK LU back-substitution.
 *
 * DGETRS reports invalid input through a negative status.  A positive status
 * is not specified and is treated as an unexpected LAPACK failure.
 */
void check_lu_solve_info(const int info, const std::string_view routine) {
  ARTS_USER_ERROR_IF(info < 0, "{} received an illegal value for argument {}.", routine, -info);
  ARTS_USER_ERROR_IF(info > 0, "{} returned an unexpected positive INFO value {}.", routine, info);
}

int complex_lapack_size(const Index n) {
  ARTS_USER_ERROR_IF(n < 0 or n > std::numeric_limits<int>::max() / 2,
                     "Complex LAPACK dimension {} is outside the supported integer range.",
                     n);
  return static_cast<int>(n);
}

bool finite_complex(const Complex value) { return std::isfinite(value.real()) and std::isfinite(value.imag()); }

void check_eigen_info(const int info) {
  ARTS_USER_ERROR_IF(info < 0, "ZGEEV received an illegal value for argument {}.", -info);
  ARTS_USER_ERROR_IF(info > 0, "ZGEEV failed to converge (INFO={}); no right eigenvectors were computed.", info);
}

}  // namespace

void svd(Matrix& U, Vector& s, Matrix& V, ConstMatrixView A, bool full_matrices) {
  ARTS_USER_ERROR_IF(A.nrows() <= 0 or A.ncols() <= 0 or A.nrows() > std::numeric_limits<int>::max() or
                         A.ncols() > std::numeric_limits<int>::max(),
                     "SVD requires nonempty dimensions within the LAPACK integer range.")
  int       m = static_cast<int>(A.nrows()), n = static_cast<int>(A.ncols());
  const int p = std::min(m, n);
  // ARTS uses row-major storage; LAPACK sees these transposed buffers as columns.
  const int left_columns = full_matrices ? m : p;
  Matrix    a(n, m), ut(left_columns, m);
  a                 = transpose(A);
  int right_columns = full_matrices ? n : p;
  V.resize(n, right_columns);  // The column-major VT output is row-major V.
  s.resize(p);
  char    jobu = full_matrices ? 'A' : 'S', jobvt = jobu;
  int     info = 0, lwork = -1;
  Numeric query = 0;
  lapack::dgesvd_(&jobu,
                  &jobvt,
                  &m,
                  &n,
                  a.data_handle(),
                  &m,
                  s.data_handle(),
                  ut.data_handle(),
                  &m,
                  V.data_handle(),
                  &right_columns,
                  &query,
                  &lwork,
                  &info);
  ARTS_USER_ERROR_IF(info != 0 or not std::isfinite(query) or query < 1 or query > std::numeric_limits<int>::max(),
                     "LAPACK DGESVD workspace query failed (INFO={}).",
                     info)
  lwork = static_cast<int>(query);
  std::vector<Numeric> work(lwork);
  lapack::dgesvd_(&jobu,
                  &jobvt,
                  &m,
                  &n,
                  a.data_handle(),
                  &m,
                  s.data_handle(),
                  ut.data_handle(),
                  &m,
                  V.data_handle(),
                  &right_columns,
                  work.data(),
                  &lwork,
                  &info);
  ARTS_USER_ERROR_IF(info != 0, "LAPACK DGESVD failed (INFO={}).", info)
  U.resize(m, left_columns);
  U = transpose(ut);
}

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
  assert((LU.shape() == std::array{n, n}));

  int  n_int, info;
  auto ipiv = std::vector<int>(n);

  LU = transpose(A);

  // Standard case: The arts matrix is not transposed, the leading
  // dimension is the row stride of the matrix.
  n_int = (int)n;

  // Compute LU decomposition using LAPACK dgetrf_.
  lapack::dgetrf_(&n_int, &n_int, LU.data_handle(), &n_int, ipiv.data(), &info);
  check_lu_info(info, "DGETRF");

  // Copy pivot array to pivot vector.
  for (Index i = 0; i < n; i++) { indx[i] = ipiv[i]; }
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
void lubacksub(VectorView x, ConstMatrixView LU, ConstVectorView b, const ArrayOfIndex& indx) {
  Index n = LU.nrows();

  /* Check if the dimensions of the input matrix and vectors agree and if LU
     is a quadratic matrix.*/

  assert((LU.shape() == std::array{n, n}));
  assert(b.size() == static_cast<Size>(n));
  assert(indx.size() == static_cast<Size>(n));
  assert(LU.stride(1) == 1);
  assert(b.stride(0) == 1);

  char                trans = 'N';
  int                 n_int = (int)n;
  int                 one   = (int)1;
  std::vector<int>    ipiv(n);
  std::vector<double> rhs(n);
  int                 info;

  for (Index i = 0; i < n; i++) {
    ipiv[i] = (int)indx[i];
    rhs[i]  = b[i];
  }

  lapack::dgetrs_(
      &trans, &n_int, &one, const_cast<Numeric*>(LU.data_handle()), &n_int, ipiv.data(), rhs.data(), &n_int, &info);
  check_lu_solve_info(info, "DGETRS");

  for (Index i = 0; i < n; i++) { x[i] = rhs[i]; }
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
  assert(n == A.nrows());
  assert(n == static_cast<Index>(x.size()));
  assert(n == static_cast<Index>(b.size()));

  // Allocate matrix and index vector for the LU decomposition.
  Matrix       LU = Matrix(n, n);
  ArrayOfIndex indx(n);

  // Perform LU decomposition.
  ludcmp(LU, indx, A);

  // Solve the system using backsubstitution.
  lubacksub(x, LU, b, indx);
}

Numeric solve(StridedComplexVectorView      x,
              StridedConstComplexMatrixView A,
              StridedConstComplexVectorView b,
              const Numeric                 min_rcond) {
  const Index n = A.ncols();
  complex_lapack_size(n);
  ARTS_USER_ERROR_IF(A.nrows() != n or x.size() != static_cast<Size>(n) or b.size() != static_cast<Size>(n),
                     "Complex solve requires a square matrix and matching input/output vector dimensions.");
  ComplexMatrix rhs(n, 1), solution(n, 1);
  for (Index i = 0; i < n; ++i) rhs[i, 0] = b[i];
  const Numeric rcond = solve(solution, A, rhs, min_rcond);
  for (Index i = 0; i < n; ++i) x[i] = solution[i, 0];
  return rcond;
}

Numeric solve(StridedComplexMatrixView      X,
              StridedConstComplexMatrixView A,
              StridedConstComplexMatrixView B,
              const Numeric                 min_rcond) {
  const Index n     = A.ncols();
  int         n_int = complex_lapack_size(n);
  ARTS_USER_ERROR_IF(A.nrows() != n or B.nrows() != n or X.shape() != B.shape(),
                     "Complex solve requires a square matrix and matching input/output matrix dimensions.");
  ARTS_USER_ERROR_IF(B.ncols() > std::numeric_limits<int>::max(),
                     "Complex solve has too many right-hand sides for LAPACK.");
  ARTS_USER_ERROR_IF(not std::isfinite(min_rcond) or min_rcond < 0 or min_rcond > 1,
                     "Complex solve requires a finite minimum reciprocal condition number in [0, 1].");
  if (n == 0) return 1;

  // Transpose into contiguous storage so LAPACK sees the original matrix.
  ComplexMatrix lu(n, n);
  lu = transpose(A);
  ComplexMatrix    rhs{transpose(B)};
  std::vector<int> ipiv(n);
  Numeric          anorm = 0;
  for (Index j = 0; j < n; ++j) {
    Numeric column_sum = 0;
    for (Index i = 0; i < n; ++i) column_sum += std::abs(A[i, j]);
    anorm = std::max(anorm, column_sum);
  }
  ARTS_USER_ERROR_IF(not std::isfinite(anorm), "Complex solve matrix norm overflowed.");

  int info = 0;
  lapack::zgetrf_(&n_int, &n_int, lu.data_handle(), &n_int, ipiv.data(), &info);
  check_lu_info(info, "ZGETRF");

  char          norm  = '1';
  Numeric       rcond = 0;
  ComplexVector work(2 * n);
  Vector        rwork(2 * n);
  lapack::zgecon_(
      &norm, &n_int, lu.data_handle(), &n_int, &anorm, &rcond, work.data_handle(), rwork.data_handle(), &info);
  check_lu_solve_info(info, "ZGECON");
  ARTS_USER_ERROR_IF(not std::isfinite(rcond) or rcond <= 0 or rcond < min_rcond,
                     "Complex solve matrix is singular or ill-conditioned (rcond={}, minimum={}).",
                     rcond,
                     min_rcond);

  char trans = 'N';
  int  nrhs  = static_cast<int>(B.ncols());
  if (nrhs > 0)
    lapack::zgetrs_(&trans, &n_int, &nrhs, lu.data_handle(), &n_int, ipiv.data(), rhs.data_handle(), &n_int, &info);
  check_lu_solve_info(info, "ZGETRS");
  X = transpose(rhs);
  return rcond;
}

void inv_inplace(MatrixView A, inv_workdata& wo) {
  Index n = A.ncols();

  // A must be a square matrix.
  assert(n == A.nrows());
  assert(n == static_cast<Index>(wo.N));

  int info;
  int n_int = (int)n;

  // Compute LU decomposition using LAPACK dgetrf_.
  lapack::dgetrf_(&n_int, &n_int, A.data_handle(), &n_int, wo.ipiv.data(), &info);
  check_lu_info(info, "DGETRF");

  // Invert matrix.
  int lwork = n_int;

  lapack::dgetri_(&n_int, A.data_handle(), &n_int, wo.ipiv.data(), wo.work.data(), &lwork, &info);
  check_lu_info(info, "DGETRI");
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
void inv_inplace(MatrixView A) {
  inv_workdata wo(A.ncols());
  inv_inplace(A, wo);
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
  Matrix LU(A);
  inv_inplace(LU);
  Ainv = LU;
}

void inv(ComplexMatrixView Ainv, const ConstComplexMatrixView A) {
  // A must be a square matrix.
  assert(A.ncols() == A.nrows());

  Index n = A.ncols();

  // Workdata
  Ainv                       = A;
  int                  n_int = int(n), lwork = n_int, info;
  std::vector<int>     ipiv(n);
  std::vector<Complex> work(lwork);

  // Compute LU decomposition using LAPACK dgetrf_.
  lapack::zgetrf_(&n_int, &n_int, const_cast<Complex*>(Ainv.data_handle()), &n_int, ipiv.data(), &info);
  check_lu_info(info, "ZGETRF");
  lapack::zgetri_(&n_int, const_cast<Complex*>(Ainv.data_handle()), &n_int, ipiv.data(), work.data(), &lwork, &info);
  check_lu_info(info, "ZGETRI");
}

void diagonalize_inplace(MatrixView P, VectorView WR, VectorView WI, MatrixView A, diagonalize_workdata& wo) {
  Index n = A.ncols();

  // A must be a square matrix.
  assert(n == A.nrows());
  assert(n == static_cast<Index>(WR.size()));
  assert(n == static_cast<Index>(WI.size()));
  assert(n == P.nrows());
  assert(n == P.ncols());
  assert(n == static_cast<Index>(wo.N));

  inplace_transpose(A);

  // Integers
  int LDA, LDA_L, LDA_R, n_int, info = 0;
  n_int = (int)n;
  LDA = LDA_L = LDA_R = (int)A.extent(0);

  // We want to calculate RP not LP
  char l_eig = 'N', r_eig = 'V';

  // Work matrix
  int lwork = std::max(4 * n_int, 2 * n_int + n_int * n_int);

  // Memory references
  double* adata  = A.data_handle();
  double* rpdata = P.data_handle();
  double* wrdata = WR.data_handle();
  double* widata = WI.data_handle();

  // Main calculations.  Note that errors in the output is ignored
  lapack::dgeev_(
      &l_eig, &r_eig, &n_int, adata, &LDA, wrdata, widata, nullptr, &LDA_L, rpdata, &LDA_R, wo.work(), &lwork, &info);

  ARTS_USER_ERROR_IF(info != 0, "DGEEV failed while diagonalizing a {}x{} matrix (INFO={})", n, n, info);

  inplace_transpose(P);
}

void diagonalize_inplace(MatrixView P, VectorView WR, VectorView WI, MatrixView A) {
  diagonalize_workdata wo(A.ncols());
  diagonalize_inplace(P, WR, WI, A, wo);
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
void diagonalize(MatrixView P, VectorView WR, VectorView WI, ConstMatrixView A) {
  Matrix A_tmp{A};
  Matrix P2{P};
  Vector WR2{WR};
  Vector WI2{WI};

  diagonalize_inplace(P2, WR2, WI2, A_tmp);

  // Re-order.  This can be done better?
  P  = P2;
  WI = WI2;
  WR = WR2;
}

void diagonalize(MatrixView P, VectorView WR, VectorView WI, ConstMatrixView A, diagonalize_workdata& wo) {
  Matrix A_tmp{A};
  Matrix P2{P};
  Vector WR2{WR};
  Vector WI2{WI};

  diagonalize_inplace(P2, WR2, WI2, A_tmp, wo);

  // Re-order.  This can be done better?
  P  = P2;
  WI = WI2;
  WR = WR2;
}

complex_diagonalize_workdata::complex_diagonalize_workdata(const Index n) {
  complex_lapack_size(n);
  N = n;
  matrix.resize(n, n);
  eigenvectors.resize(n, n);
  eigenvalues.resize(n);
  rwork.resize(2 * n);
}

void diagonalize(StridedComplexMatrixView P, StridedComplexVectorView W, StridedConstComplexMatrixView A) {
  ARTS_USER_ERROR_IF(A.nrows() != A.ncols() or P.shape() != A.shape() or W.size() != static_cast<Size>(A.ncols()),
                     "Complex diagonalization requires a square matrix and matching output dimensions.");
  complex_diagonalize_workdata workdata(A.ncols());
  diagonalize(P, W, A, workdata);
}

void diagonalize(StridedComplexMatrixView      P,
                 StridedComplexVectorView      W,
                 StridedConstComplexMatrixView A,
                 complex_diagonalize_workdata& workdata) {
  const Index n     = A.ncols();
  int         n_int = complex_lapack_size(n);
  ARTS_USER_ERROR_IF(A.nrows() != n or P.shape() != A.shape() or W.size() != static_cast<Size>(n),
                     "Complex diagonalization requires a square matrix and matching output dimensions.");
  ARTS_USER_ERROR_IF(
      workdata.N != n or workdata.matrix.shape() != A.shape() or workdata.eigenvectors.shape() != A.shape() or
          workdata.eigenvalues.size() != static_cast<Size>(n) or workdata.rwork.size() < static_cast<Size>(2 * n),
      "Complex diagonalization workspace has incompatible dimensions.");
  if (n == 0) return;

  // All LAPACK buffers are contiguous. Copy outputs only after successful
  // completion, so strided views and input/output aliasing are supported.
  workdata.matrix = transpose(A);
  int     info = 0, lda_l = 1;
  char    l_eig = 'N', r_eig = 'V';
  Complex left_eigenvector_dummy{};

  if (workdata.work.empty()) {
    int     lwork = -1;
    Complex query{};
    lapack::zgeev_(&l_eig,
                   &r_eig,
                   &n_int,
                   workdata.matrix.data_handle(),
                   &n_int,
                   workdata.eigenvalues.data_handle(),
                   &left_eigenvector_dummy,
                   &lda_l,
                   workdata.eigenvectors.data_handle(),
                   &n_int,
                   &query,
                   &lwork,
                   workdata.rwork.data_handle(),
                   &info);
    check_eigen_info(info);
    ARTS_USER_ERROR_IF(
        not finite_complex(query) or query.real() < 2 * n_int or query.real() > std::numeric_limits<int>::max(),
        "ZGEEV returned an invalid workspace size {}.",
        query);
    workdata.work.resize(static_cast<Index>(query.real()));
  }

  ARTS_USER_ERROR_IF(workdata.work.size() < static_cast<Size>(2 * n) or
                         workdata.work.size() > static_cast<Size>(std::numeric_limits<int>::max()),
                     "Complex diagonalization workspace size is outside the LAPACK range.");
  int lwork = static_cast<int>(workdata.work.size());
  lapack::zgeev_(&l_eig,
                 &r_eig,
                 &n_int,
                 workdata.matrix.data_handle(),
                 &n_int,
                 workdata.eigenvalues.data_handle(),
                 &left_eigenvector_dummy,
                 &lda_l,
                 workdata.eigenvectors.data_handle(),
                 &n_int,
                 workdata.work.data_handle(),
                 &lwork,
                 workdata.rwork.data_handle(),
                 &info);
  check_eigen_info(info);
  P = transpose(workdata.eigenvectors);
  W = workdata.eigenvalues;
}

void diagonalize(StridedComplexMatrixView      P,
                 StridedComplexVectorView      W,
                 StridedComplexMatrixView      dP,
                 StridedComplexVectorView      dW,
                 StridedConstComplexMatrixView A,
                 StridedConstComplexMatrixView dA) {
  complex_diagonalize_workdata workdata(A.ncols());
  diagonalize(P, W, dP, dW, A, dA, workdata);
}

void diagonalize(StridedComplexMatrixView      P,
                 StridedComplexVectorView      W,
                 StridedComplexMatrixView      dP,
                 StridedComplexVectorView      dW,
                 StridedConstComplexMatrixView A,
                 StridedConstComplexMatrixView dA,
                 complex_diagonalize_workdata& workdata) {
  const Index n = A.ncols();
  ARTS_USER_ERROR_IF(A.nrows() != n or P.shape() != A.shape() or dP.shape() != A.shape() or dA.shape() != A.shape() or
                         W.size() != static_cast<Size>(n) or dW.size() != static_cast<Size>(n),
                     "Complex eigendecomposition derivative requires matching square matrix and vector dimensions.");
  ComplexTensor3 directions(1, n, n), derivatives(1, n, n);
  ComplexMatrix  value_derivatives(1, n);
  directions[0] = dA;
  diagonalize(P, W, derivatives, value_derivatives, A, directions, workdata);
  dP = derivatives[0];
  dW = value_derivatives[0];
}

void diagonalize(StridedComplexMatrixView       P,
                 StridedComplexVectorView       W,
                 StridedComplexTensor3View      dP,
                 StridedComplexMatrixView       dW,
                 StridedConstComplexMatrixView  A,
                 StridedConstComplexTensor3View dA) {
  complex_diagonalize_workdata workdata(A.ncols());
  diagonalize(P, W, dP, dW, A, dA, workdata);
}

void diagonalize(StridedComplexMatrixView       P,
                 StridedComplexVectorView       W,
                 StridedComplexTensor3View      dP,
                 StridedComplexMatrixView       dW,
                 StridedConstComplexMatrixView  A,
                 StridedConstComplexTensor3View dA,
                 complex_diagonalize_workdata&  workdata) {
  const Index n = A.ncols(), nq = dA.npages();
  ARTS_USER_ERROR_IF(A.nrows() != n or P.shape() != A.shape() or dP.shape() != dA.shape() or dA.nrows() != n or
                         dA.ncols() != n or W.size() != static_cast<Size>(n) or dW.nrows() != nq or dW.ncols() != n,
                     "Batched complex eigendecomposition derivatives require matching matrix and Jacobian dimensions.");
  if (nq == 0) {
    diagonalize(P, W, A, workdata);
    return;
  }
  ARTS_USER_ERROR_IF(n > 0 and nq > std::numeric_limits<int>::max() / n,
                     "Batched eigendecomposition has too many derivative right-hand sides.");
  ComplexMatrix  vectors(n, n), value_derivatives(nq, n);
  ComplexTensor3 derivatives(nq, n, n, 0);
  ComplexVector  values(n);
  diagonalize(vectors, values, A, workdata);
  if (n == 0) return;
  // Subtract a scalar carrier when measuring the gap scale. The derivative
  // should not become ill-conditioned merely by adding a large multiple of I.
  Numeric scale = 0;
  for (Index i = 0; i < n; ++i) {
    Numeric row_sum = 0;
    for (Index j = 0; j < n; ++j) { row_sum += std::abs(A[i, j] - (i == j ? A[0, 0] : Complex{})); }
    scale = std::max(scale, row_sum);
  }
  ARTS_USER_ERROR_IF(not std::isfinite(scale), "Complex eigendecomposition derivative norm overflowed.");
  const Numeric gap_tolerance = 64 * std::numeric_limits<Numeric>::epsilon() * scale;
  for (Index i = 0; i < n; ++i)
    for (Index j = 0; j < i; ++j)
      ARTS_USER_ERROR_IF(std::abs(values[i] - values[j]) <= gap_tolerance,
                         "Complex eigendecomposition derivative requires separated eigenvalues; modes {} and {} "
                         "have gap {} (tolerance {}).",
                         i,
                         j,
                         std::abs(values[i] - values[j]),
                         gap_tolerance);

  // Solve P*B=dA*P with one factorization and all nq*n right-hand sides, without
  // forming P^-1. The local buffer stores the transpose of all RHS columns.
  ComplexMatrix transformed(nq * n, n, 0);
  for (Index q = 0; q < nq; ++q)
    for (Index col = 0; col < n; ++col)
      for (Index row = 0; row < n; ++row)
        for (Index k = 0; k < n; ++k) transformed[q * n + col, row] += dA[q, row, k] * vectors[k, col];
  solve(transpose(transformed), vectors, transpose(transformed), 1e-12);

  for (Index q = 0; q < nq; ++q) {
    for (Index i = 0; i < n; ++i) {
      value_derivatives[q, i] = transformed[q * n + i, i];
      for (Index j = 0; j < n; ++j) {
        if (i == j) continue;
        const Complex coefficient = transformed[q * n + i, j] / (values[i] - values[j]);
        for (Index row = 0; row < n; ++row) derivatives[q, row, i] += vectors[row, j] * coefficient;
      }
      Complex projection   = 0;
      Numeric squared_norm = 0;
      for (Index row = 0; row < n; ++row) {
        projection   += std::conj(vectors[row, i]) * derivatives[q, row, i];
        squared_norm += std::norm(vectors[row, i]);
      }
      for (Index row = 0; row < n; ++row) derivatives[q, row, i] -= vectors[row, i] * projection / squared_norm;
    }
  }
  P  = vectors;
  W  = values;
  dP = derivatives;
  dW = value_derivatives;
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
  assert((A.shape() == std::array{n, n}));
  assert((F.shape() == std::array{n, n}));

  Numeric A_norm_inf, c;
  Numeric j;
  Matrix  D(n, n), N(n, n), X(n, n), cX(n, n, 0.0), B(n, n, 0.0);
  Vector  N_col_vec(n, 0.), F_col_vec(n, 0.);

  A_norm_inf = norm_inf(A);

  // This formular is derived in the book by Golub and Van Loan.
  j = 1 + floor(1. / log(2.) * log(A_norm_inf));

  if (j < 0) j = 0.;
  auto j_index = (Index)(j);

  // Scale matrix
  F  = A;
  F /= nonstd::pow(2.0, j);

  /* The higher q the more accurate is the computation,
     see user guide for accuracy */
  //  Index q = 8;
  auto q_n = (Numeric)(q);
  id_mat(D);
  id_mat(N);
  id_mat(X);
  c = 1.;

  for (Index k = 0; k < q; k++) {
    auto k_n  = (Numeric)(k + 1);
    c        *= (q_n - k_n + 1) / ((2 * q_n - k_n + 1) * k_n);
    mult(B, F, X);  // X = F * X
    X   = B;
    cX  = X;
    cX *= c;                     // cX = X*c
    N  += cX;                    // N = N + X*c
    cX *= nonstd::pow(-1, k_n);  // cX = (-1)^k*c*X
    D  += cX;                    // D = D + (-1)^k*c*X
  }

  /*Solve the equation system DF=N for F using LU decomposition,
   use the backsubstitution routine for columns of N*/

  /* Now use X  for the LU decomposition matrix of D.*/
  ArrayOfIndex indx(n);

  ludcmp(X, indx, D);

  for (Index i = 0; i < n; i++) {
    N_col_vec = N[joker, i];  // extract column vectors of N
    lubacksub(F_col_vec, X, N_col_vec, indx);
    F[joker, i] = F_col_vec;  // construct F matrix  from column vectors
  }

  /* The square of F gives the result. */
  for (Index k = 0; k < j_index; k++) {
    mult(B, F, F);  // F = F^2
    F = B;
  }
}

//! 2-norm of a vector
/*!

  \param v Input: vector, with length>=1

  \return Norm
*/
Numeric norm2(ConstVectorView v) {
  assert(v.size());
  return sqrt(dot(v, v));
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
    for (Index i = 0; i < A.ncols(); i++) row_sum += std::abs(A[i, j]);
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
  assert(n == I.nrows());

  I = 0;
  for (Index i = 0; i < n; i++) I[i, i] = 1.;
}

/*!
    Determinant of N by N matrix. Simple recursive method.

    \param  A   In:    Matrix of size NxN.

    \author Richard Larsson
    \date   2012-08-03
*/
Numeric det(ConstMatrixView A) { return matpack::eigen::as_eigen(A).determinant(); }

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
  const Size n = x.size();

  assert(y.size() == n);

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

  for (Size i = 0; i < n; i++) { xm += x[i] / Numeric(n); }

  for (Size i = 0; i < n; i++) {
    const Numeric xv  = x[i] - xm;
    s1               += xv * y[i];
    s3               += y[i];
    s4               += xv * xv;
  }

  p[1] = s1 / s4;
  p[0] = s3 / Numeric(n) - p[1] * xm;
}

Numeric lsf(VectorView x, ConstMatrixView A, ConstVectorView y, bool residual) noexcept {
  // Size of the problem
  const Index n = x.size();
  Matrix      AT, ATA(n, n);
  Vector      ATy(n);

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
    return dot(r, r);
  }

  return 0;
}

void solve_inplace(VectorView X, MatrixView A, solve_workdata& wo) {
  // Assert that A is quadratic.
  const Index n = A.nrows();
  assert((A.shape() == std::array{n, n}));

  char trans = 'N';
  int  n_int = static_cast<int>(n);
  int  info{};
  int  one = 1;

  // Standard case: The arts matrix is not transposed, the leading
  // dimension is the row stride of the matrix.
  n_int = (int)n;

  // Compute LU decomposition using LAPACK dgetrf_.
  lapack::dgetrf_(
      &n_int, &n_int, const_cast<Numeric*>(inplace_transpose(A).data_handle()), &n_int, wo.ipiv.data(), &info);
  check_lu_info(info, "DGETRF");

  lapack::dgetrs_(&trans, &n_int, &one, A.data_handle(), &n_int, wo.ipiv.data(), X.data_handle(), &n_int, &info);
  check_lu_solve_info(info, "DGETRS");
}

void solve_inplace(VectorView X, MatrixView A) {
  solve_workdata wo(A.ncols());
  solve_inplace(X, A, wo);
}
