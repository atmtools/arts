/*!
     \file   lin_alg.h
     \author Claudia Emde <claudia.emde@dlr.de>
     \date   Thu May  2 14:34:05 2002

     \brief  Linear algebra functions.

   */

#ifndef linalg_h
#define linalg_h

#include "matpack_arrays.h"

// LU decomposition
void ludcmp(Matrix& LU, ArrayOfIndex& indx, ConstMatrixView A);

// LU backsubstitution
void lubacksub(VectorView x, ConstMatrixView LU, ConstVectorView b, const ArrayOfIndex& indx);

// Solve linear system
void solve(VectorView x, ConstMatrixView A, ConstVectorView b);

/** Solve A*x=b and return an estimate of the reciprocal 1-norm condition
 * number. Inputs are preserved, except that x may alias b. Strided views are
 * supported. Throws for a singular A or rcond < min_rcond, and leaves x
 * unchanged on failure. An empty system returns 1.
 *
 * A nonfinite A is rejected through its 1-norm. A nonfinite b is not searched
 * for and propagates into x: callers detect it on their own result rather than
 * paying for a scan of every right-hand side here.
 */
Numeric solve(StridedComplexVectorView      x,
              StridedConstComplexMatrixView A,
              StridedConstComplexVectorView b,
              Numeric                       min_rcond = 0);

/** Solve A*X=B for all right-hand-side columns with one LU factorization and
 * return the reciprocal 1-norm condition estimate. Supports strided views and
 * X aliasing B. Throws for a singular A or rcond < min_rcond and leaves X
 * unchanged on failure. An empty matrix A returns 1.
 *
 * A nonfinite A is rejected through its 1-norm. A nonfinite B is not searched
 * for and propagates into X.
 */
Numeric solve(StridedComplexMatrixView      X,
              StridedConstComplexMatrixView A,
              StridedConstComplexMatrixView B,
              Numeric                       min_rcond = 0);

/** A = U Sigma V^T, via LAPACK. s contains the min(m,n) singular values.
 * U and V are square when full_matrices is true; otherwise both have min(m,n)
 * columns. Input is preserved.
 */
void svd(Matrix& U, Vector& s, Matrix& V, ConstMatrixView A, bool full_matrices = true);

struct solve_workdata {
  std::size_t      N{};
  std::vector<int> ipiv{};

  constexpr solve_workdata()                                 = default;
  constexpr solve_workdata(const solve_workdata&)            = default;
  constexpr solve_workdata(solve_workdata&&)                 = default;
  constexpr solve_workdata& operator=(const solve_workdata&) = default;
  constexpr solve_workdata& operator=(solve_workdata&&)      = default;

  constexpr solve_workdata(std::size_t N_) : N(N_), ipiv(N) {}
  constexpr void resize(std::size_t N_) {
    N = N_;
    ipiv.resize(N);
  }
};

/*! Solves A X = B inplace using dgesv.
  * 
  * Returns the Lapack ipiv array.
  *
  * @param[in,out] X   As equation, on input it is B on output is is X
  * @param[in]     A   As equation, it is destroyed on output (LU decomposition)
  * @throws If the system cannot be solved according to Lapack info
  */
void solve_inplace(VectorView X, MatrixView A, solve_workdata& wo);

//! As above but allocates WO
void solve_inplace(VectorView X, MatrixView A);

struct inv_workdata {
  std::size_t          N{};
  std::vector<int>     ipiv{};
  std::vector<Numeric> work{};

  constexpr inv_workdata()                               = default;
  constexpr inv_workdata(const inv_workdata&)            = default;
  constexpr inv_workdata(inv_workdata&&)                 = default;
  constexpr inv_workdata& operator=(const inv_workdata&) = default;
  constexpr inv_workdata& operator=(inv_workdata&&)      = default;

  constexpr inv_workdata(std::size_t N_) : N(N_), ipiv(N), work(N) {}
  constexpr void resize(std::size_t N_) {
    N = N_;
    ipiv.resize(N);
    work.resize(N);
  }
};

// Matrix inverse
void inv(MatrixView Ainv, ConstMatrixView A);

// Matrix inverse in place with destructive consequences
void inv_inplace(MatrixView A);

// Matrix inverse in place with destructive consequences
void inv_inplace(MatrixView A, inv_workdata& wo);

// Matrix inverse
void inv(ComplexMatrixView Ainv, const ConstComplexMatrixView A);

struct diagonalize_workdata {
  std::size_t          N{};
  std::vector<Numeric> w{};

  constexpr diagonalize_workdata()                                       = default;
  constexpr diagonalize_workdata(const diagonalize_workdata&)            = default;
  constexpr diagonalize_workdata(diagonalize_workdata&&)                 = default;
  constexpr diagonalize_workdata& operator=(const diagonalize_workdata&) = default;
  constexpr diagonalize_workdata& operator=(diagonalize_workdata&&)      = default;

  constexpr diagonalize_workdata(std::size_t N_) : N(N_), w(4 * N + N * N) {}
  constexpr Numeric* work() { return w.data(); }
  constexpr Numeric* rwork() { return w.data() + 2 * N; }
};

struct complex_diagonalize_workdata {
  Index         N{};
  ComplexMatrix matrix;
  ComplexMatrix eigenvectors;
  ComplexVector eigenvalues;
  ComplexVector work;
  Vector        rwork;

  complex_diagonalize_workdata() = default;
  explicit complex_diagonalize_workdata(Index n);
};

// Matrix diagonalization with lapack
void diagonalize(MatrixView P, VectorView WR, VectorView WI, ConstMatrixView A);

// Matrix diagonalization with lapack
void diagonalize(MatrixView P, VectorView WR, VectorView WI, ConstMatrixView A, diagonalize_workdata& wo);

// Same as diagonalize but inplace manilpulation of input with destructive consqeuences
void diagonalize_inplace(MatrixView P, VectorView WR, VectorView WI, MatrixView A);

// Same as diagonalize but inplace manilpulation of input with destructive consqeuences
void diagonalize_inplace(MatrixView P, VectorView WR, VectorView WI, MatrixView A, diagonalize_workdata& wo);

/** Complex eigendecomposition A*P=P*diag(W), with right eigenvectors in the
 * columns of P. Supports strided views and preserves A unless it aliases an
 * output. Throws on invalid dimensions or LAPACK failure; outputs are
 * unchanged on failure. Empty matrices are accepted. A nonfinite A is not
 * searched for: LAPACK reports what it detects, and anything else propagates
 * into P and W.
 */
void diagonalize(StridedComplexMatrixView P, StridedComplexVectorView W, StridedConstComplexMatrixView A);
void diagonalize(StridedComplexMatrixView      P,
                 StridedComplexVectorView      W,
                 StridedConstComplexMatrixView A,
                 complex_diagonalize_workdata& workdata);

/** Complex eigendecomposition and its directional derivative for dA.
 * Requires distinct eigenvalues and a well-conditioned eigenvector basis;
 * throws when the eigenvalue gaps are numerically unresolved or rcond(P)<1e-12.
 * The derivative uses the unit-norm, parallel-transport gauge P_i^H*dP_i=0.
 * Thus dP need not differentiate LAPACK's phase convention; phase-invariant
 * quantities and dW are independent of that convention. Supports strided views,
 * preserves inputs unless they alias outputs, and leaves outputs unchanged on
 * failure. The four outputs must not overlap each other. A nonfinite A is
 * rejected through the gap scale; a nonfinite dA propagates into dP and dW.
 */
void diagonalize(StridedComplexMatrixView      P,
                 StridedComplexVectorView      W,
                 StridedComplexMatrixView      dP,
                 StridedComplexVectorView      dW,
                 StridedConstComplexMatrixView A,
                 StridedConstComplexMatrixView dA);
void diagonalize(StridedComplexMatrixView      P,
                 StridedComplexVectorView      W,
                 StridedComplexMatrixView      dP,
                 StridedComplexVectorView      dW,
                 StridedConstComplexMatrixView A,
                 StridedConstComplexMatrixView dA,
                 complex_diagonalize_workdata& workdata);

/** Batched version: dA[q,:,:] and dP[q,:,:] are matrix directions and
 * eigenvector derivatives, while dW[q,:] contains eigenvalue derivatives.
 * Computes the eigendecomposition and eigenvector LU factorization once for
 * all directions. With zero directions, this is ordinary diagonalization.
 * The scalar overload above delegates to a one-direction batch.
 */
void diagonalize(StridedComplexMatrixView       P,
                 StridedComplexVectorView       W,
                 StridedComplexTensor3View      dP,
                 StridedComplexMatrixView       dW,
                 StridedConstComplexMatrixView  A,
                 StridedConstComplexTensor3View dA);
void diagonalize(StridedComplexMatrixView       P,
                 StridedComplexVectorView       W,
                 StridedComplexTensor3View      dP,
                 StridedComplexMatrixView       dW,
                 StridedConstComplexMatrixView  A,
                 StridedConstComplexTensor3View dA,
                 complex_diagonalize_workdata&  workdata);

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
Numeric lsf(VectorView x, ConstMatrixView A, ConstVectorView y, bool residual = true) noexcept;

#endif  // linalg_h
