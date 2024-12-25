#include "matpack_mdspan_math.h"

#include "blas.h"
#include "matpack_mdspan_helpers_eigen.h"

void mult(StridedMatrixView A,
          const StridedConstMatrixView &B,
          const StridedConstMatrixView &C,
          Numeric alpha,
          Numeric beta) {
  // Check dimensions:
  assert(A.nrows() == B.nrows());
  assert(A.ncols() == C.ncols());
  assert(B.ncols() == C.nrows());

  // Catch trivial case if one of the matrices is empty.
  if ((B.nrows() == 0) || (B.ncols() == 0) || (C.ncols() == 0)) return;

  // Matrices B and C must be continuous in at least on dimension,  C
  // must be continuous along the second dimension.
  if (((B.stride(0) == 1) || (B.stride(1) == 1)) &&
      ((C.stride(0) == 1) || (C.stride(1) == 1)) && (A.stride(1) == 1)) {
    // BLAS uses column-major order while arts uses row-major order.
    // Hence instead of C = A * B we compute C^T = A^T * B^T!

    int k, m, n;

    k = (int)B.ncols();
    m = (int)C.ncols();
    n = (int)B.nrows();

    // Note also the clash in nomenclature: BLAS uses C = A * B while
    // arts uses A = B * C. Taking into accout this and the difference in
    // memory layouts, we need to map the MatrixViews A, B and C to the BLAS
    // arguments as follows:
    // A (arts) -> C (BLAS)
    // B (arts) -> B (BLAS)
    // C (arts) -> A (BLAS)

    // Char indicating whether A (BLAS) or B (BLAS) should be transposed.
    char transa, transb;
    // Sizes of the matrices along the direction in which they are
    // traversed.
    int lda, ldb, ldc;

    // Check if C (arts) is transposed.
    if (C.stride(0) == 1) {
      transa = 'T';
      lda    = (int)C.stride(1);
    } else {
      transa = 'N';
      lda    = (int)C.stride(0);
    }

    // Check if B (arts) is transposed.
    if (B.stride(0) == 1) {
      transb = 'T';
      ldb    = (int)B.stride(1);
    } else {
      transb = 'N';
      ldb    = (int)B.stride(0);
    }

    // In case B (arts) has only one column, column and row stride are 1.
    // We therefore need to set ldb to k, since dgemm_ requires lda to be at
    // least k / m if A is non-transposed / transposed.
    if ((B.stride(1) == 1) && (B.stride(0) == 1)) {
      transb = 'N';
      ldb    = k;
    }

    // The same holds for C (arts).
    if ((C.stride(1) == 1) && (C.stride(0) == 1)) {
      transa = 'N';
      lda    = m;
    }

    ldc = (int)A.stride(0);
    // The same holds for A (arts).
    if ((A.stride(1) == 1) && (A.stride(0) == 1)) {
      ldc = m;
    }

    dgemm_(&transa,
           &transb,
           &m,
           &n,
           &k,
           &alpha,
           const_cast<Numeric *>(C.data_handle()),
           &lda,
           const_cast<Numeric *>(B.data_handle()),
           &ldb,
           &beta,
           A.data_handle(),
           &ldc);

  } else {
    if (beta == 0.0) {
      matpack::eigen::as_eigen(A).noalias() = alpha * B * C;
    } else {
      A                                     *= beta;
      matpack::eigen::as_eigen(A).noalias() += alpha * B * C;
    }
  }
}

void mult(StridedComplexMatrixView A,
          const StridedConstComplexMatrixView &B,
          const StridedConstComplexMatrixView &C,
          Complex alpha,
          Complex beta) {
  // Check dimensions:
  assert(A.nrows() == B.nrows());
  assert(A.ncols() == C.ncols());
  assert(B.ncols() == C.nrows());

  // Catch trivial case if one of the matrices is empty.
  if ((B.nrows() == 0) || (B.ncols() == 0) || (C.ncols() == 0)) return;

  // Matrices B and C must be continuous in at least on dimension,  C
  // must be continuous along the second dimension.
  if (((B.stride(0) == 1) || (B.stride(1) == 1)) &&
      ((C.stride(0) == 1) || (C.stride(1) == 1)) && (A.stride(1) == 1)) {
    // BLAS uses column-major order while arts uses row-major order.
    // Hence instead of C = A * B we compute C^T = A^T * B^T!

    int k, m, n;

    k = (int)B.ncols();
    m = (int)C.ncols();
    n = (int)B.nrows();

    // Note also the clash in nomenclature: BLAS uses C = A * B while
    // arts uses A = B * C. Taking into accout this and the difference in
    // memory layouts, we need to map the MatrixViews A, B and C to the BLAS
    // arguments as follows:
    // A (arts) -> C (BLAS)
    // B (arts) -> B (BLAS)
    // C (arts) -> A (BLAS)

    // Char indicating whether A (BLAS) or B (BLAS) should be transposed.
    char transa, transb;
    // Sizes of the matrices along the direction in which they are
    // traversed.
    int lda, ldb, ldc;

    // Check if C (arts) is transposed.
    if (C.stride(0) == 1) {
      transa = 'T';
      lda    = (int)C.stride(1);
    } else {
      transa = 'N';
      lda    = (int)C.stride(0);
    }

    // Check if B (arts) is transposed.
    if (B.stride(0) == 1) {
      transb = 'T';
      ldb    = (int)B.stride(1);
    } else {
      transb = 'N';
      ldb    = (int)B.stride(0);
    }

    // In case B (arts) has only one column, column and row stride are 1.
    // We therefore need to set ldb to k, since dgemm_ requires lda to be at
    // least k / m if A is non-transposed / transposed.
    if ((B.stride(1) == 1) && (B.stride(0) == 1)) {
      transb = 'N';
      ldb    = k;
    }

    // The same holds for C (arts).
    if ((C.stride(1) == 1) && (C.stride(0) == 1)) {
      transa = 'N';
      lda    = m;
    }

    ldc = (int)A.stride(0);
    // The same holds for A (arts).
    if ((A.stride(1) == 1) && (A.stride(0) == 1)) {
      ldc = m;
    }

    zgemm_(&transa,
           &transb,
           &m,
           &n,
           &k,
           &alpha,
           const_cast<Complex *>(C.data_handle()),
           &lda,
           const_cast<Complex *>(B.data_handle()),
           &ldb,
           &beta,
           A.data_handle(),
           &ldc);

  } else {
    if (beta == 0.0) {
      matpack::eigen::as_eigen(A).noalias() = alpha * B * C;
    } else {
      A                                     *= beta;
      matpack::eigen::as_eigen(A).noalias() += alpha * B * C;
    }
  }
}

void mult(StridedVectorView y,
          const StridedConstMatrixView &M,
          const StridedConstVectorView &x,
          Numeric alpha,
          Numeric beta) {
  assert(y.size() == static_cast<Size>(M.nrows()));
  assert(static_cast<Size>(M.ncols()) == x.size());
  assert(not M.empty());

  if ((M.stride(0) == 1) || (M.stride(1) == 1)) {
    char trans;
    int m, n;
    int LDA, incx, incy;

    if (M.stride(1) != 1) {
      trans = 'n';
      m     = (int)M.nrows();
      n     = (int)M.ncols();
      LDA   = (int)M.stride(1);
    } else {
      trans = 't';
      m     = (int)M.ncols();
      n     = (int)M.nrows();
      LDA   = (int)M.stride(0);
      if (M.stride(0) == 1) LDA = m;
    }

    incx = (int)x.stride(0);
    incy = (int)y.stride(0);

    double *mstart = const_cast<Numeric *>(M.data_handle());
    double *ystart = y.data_handle();
    double *xstart = const_cast<Numeric *>(x.data_handle());

    dgemv_(&trans,
           &m,
           &n,
           &alpha,
           mstart,
           &LDA,
           xstart,
           &incx,
           &beta,
           ystart,
           &incy);
  } else {
    if (beta == 0.0) {
      matpack::eigen::as_eigen(y).noalias() = alpha * M * x;
    } else {
      y                                     *= beta;
      matpack::eigen::as_eigen(y).noalias() += alpha * M * x;
    }
  }
}

void mult(StridedComplexVectorView y,
          const StridedConstComplexMatrixView &M,
          const StridedConstComplexVectorView &x,
          Complex alpha,
          Complex beta) {
  assert(y.size() == static_cast<Size>(M.nrows()));
  assert(static_cast<Size>(M.ncols()) == x.size());
  assert(not M.empty());

  if ((M.stride(0) == 1) || (M.stride(1) == 1)) {
    char trans;
    int m, n;
    int LDA, incx, incy;

    if (M.stride(1) != 1) {
      trans = 'n';
      m     = (int)M.nrows();
      n     = (int)M.ncols();
      LDA   = (int)M.stride(1);
    } else {
      trans = 't';
      m     = (int)M.ncols();
      n     = (int)M.nrows();
      LDA   = (int)M.stride(0);
      if (M.stride(0) == 1) LDA = m;
    }

    incx = (int)x.stride(0);
    incy = (int)y.stride(0);

    auto *mstart = const_cast<Complex *>(M.data_handle());
    auto *ystart = y.data_handle();
    auto *xstart = const_cast<Complex *>(x.data_handle());

    zgemv_(&trans,
           &m,
           &n,
           &alpha,
           mstart,
           &LDA,
           xstart,
           &incx,
           &beta,
           ystart,
           &incy);
  } else {
    if (beta == 0.0) {
      matpack::eigen::as_eigen(y).noalias() = alpha * M * x;
    } else {
      y                                     *= beta;
      matpack::eigen::as_eigen(y).noalias() += alpha * M * x;
    }
  }
}
