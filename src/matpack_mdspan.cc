#include <matpack_mdspan.h>

#include "blas.h"

void mult_fast(matpack::md::simple_view<Numeric, 1, false> y,
               const matpack::md::simple_view<Numeric, 2, true> &A,
               const matpack::md::simple_view<Numeric, 1, true> &x,
               Numeric alpha, Numeric beta) {
  int M = static_cast<int>(x.nelem());
  int N = static_cast<int>(y.nelem());

  ARTS_ASSERT(A.shape() == (std::array<Index, 2>{N, M}), A.nrows(), ", ",
              A.ncols(), " != ", N, ", ", M)

  if (N == 0 or M == 0)
    return;

  char trans = 'T';
  int LDA = M;
  int incx = 1;
  int incy = 1;

  dgemv_(&trans, &M, &N, &alpha, A.data_handle(), &LDA, x.data_handle(), &incx,
         &beta, y.data_handle(), &incy);
}

void mult(matpack::md::strided_view<Numeric, 1, false> y,
          const matpack::md::strided_view<Numeric, 2, true> &A,
          const matpack::md::strided_view<Numeric, 1, true> &x) {
  ARTS_ASSERT(A.nrows() == y.nelem() and A.ncols() == x.nelem())

  const Index N = y.nelem();
  for (Index i = 0; i < N; i++)
    y[i] = A[i] * x;
}

void mult_fast(matpack::md::simple_view<Numeric, 2, false> C,
               const matpack::md::simple_view<Numeric, 2, true> &A,
               const matpack::md::simple_view<Numeric, 2, true> &B,
               Numeric alpha, Numeric beta) {
  int M = static_cast<int>(B.ncols());
  int N = static_cast<int>(A.nrows());
  int K = static_cast<int>(A.ncols());

  ARTS_ASSERT(C.shape() == (std::array<Index, 2>{N, M}), C.nrows(), ", ",
              C.ncols(), " != ", N, ", ", M);
  ARTS_ASSERT(A.shape() == (std::array<Index, 2>{N, K}), A.nrows(), ", ",
              A.ncols(), " != ", N, ", ", K);
  ARTS_ASSERT(B.shape() == (std::array<Index, 2>{K, M}), B.nrows(), ", ",
              B.ncols(), " != ", K, ", ", M);

  if (N == 0 or M == 0 or K == 0)
    return;

  char transa = 'N';
  char transb = 'N';
  int LDA = M;
  int LDB = K;
  int LDC = M;

  //! Note that ARTS is row-major, so while the expression is C = A * B, we have
  //! to write C = (A*B)^T, or C = B^T * A^T
  dgemm_(&transa, &transb, &M, &N, &K, &alpha, B.data_handle(), &LDA,
         A.data_handle(), &LDB, &beta, C.data_handle(), &LDC);
}

void mult(matpack::md::strided_view<Numeric, 2, false> C,
          const matpack::md::strided_view<Numeric, 2, true> &A,
          const matpack::md::strided_view<Numeric, 2, true> &B) {
  const auto [N, M] = C.shape();

  for (Index i = 0; i < N; i++) {
    for (Index j = 0; j < M; j++) {
      C(i, j) = A[i] * B(joker, j);
    }
  }
}

void mult_fast(matpack::md::simple_view<Complex, 1, false> y,
               const matpack::md::simple_view<Complex, 2, true> &A,
               const matpack::md::simple_view<Complex, 1, true> &x,
               Complex alpha, Complex beta) {
  int M = static_cast<int>(x.nelem());
  int N = static_cast<int>(y.nelem());

  ARTS_ASSERT(A.shape() == (std::array<Index, 2>{N, M}), A.nrows(), ", ",
              A.ncols(), " != ", N, ", ", M)

  if (N == 0 or M == 0)
    return;

  char trans = 'T';
  int LDA = M;
  int incx = 1;
  int incy = 1;

  zgemv_(&trans, &M, &N, &alpha, A.data_handle(), &LDA, x.data_handle(), &incx,
         &beta, y.data_handle(), &incy);
}

void mult(matpack::md::strided_view<Complex, 1, false> y,
          const matpack::md::strided_view<Complex, 2, true> &A,
          const matpack::md::strided_view<Complex, 1, true> &x) {
  ARTS_ASSERT(A.nrows() == y.nelem() and A.ncols() == x.nelem())

  const Index N = y.nelem();
  for (Index i = 0; i < N; i++)
    y[i] = A[i] * x;
}

void mult_fast(matpack::md::simple_view<Complex, 2, false> C,
               const matpack::md::simple_view<Complex, 2, true> &A,
               const matpack::md::simple_view<Complex, 2, true> &B,
               Complex alpha, Complex beta) {
  int M = static_cast<int>(B.ncols());
  int N = static_cast<int>(A.nrows());
  int K = static_cast<int>(A.ncols());

  ARTS_ASSERT(C.shape() == (std::array<Index, 2>{N, M}), C.nrows(), ", ",
              C.ncols(), " != ", N, ", ", M);
  ARTS_ASSERT(A.shape() == (std::array<Index, 2>{N, K}), A.nrows(), ", ",
              A.ncols(), " != ", N, ", ", K);
  ARTS_ASSERT(B.shape() == (std::array<Index, 2>{K, M}), B.nrows(), ", ",
              B.ncols(), " != ", K, ", ", M);

  if (N == 0 or M == 0 or K == 0)
    return;

  char transa = 'N';
  char transb = 'N';
  int LDA = M;
  int LDB = K;
  int LDC = M;

  //! Note that ARTS is row-major, so while the expression is C = A * B, we have
  //! to write C = (A*B)^T, or C = B^T * A^T
  zgemm_(&transa, &transb, &M, &N, &K, &alpha, B.data_handle(), &LDA,
         A.data_handle(), &LDB, &beta, C.data_handle(), &LDC);
}

void mult(matpack::md::strided_view<Complex, 2, false> C,
          const matpack::md::strided_view<Complex, 2, true> &A,
          const matpack::md::strided_view<Complex, 2, true> &B) {
  const auto [N, M] = C.shape();

  for (Index i = 0; i < N; i++) {
    for (Index j = 0; j < M; j++) {
      C(i, j) = A[i] * B(joker, j);
    }
  }
}
