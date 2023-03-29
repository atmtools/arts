#include "matpack_eigen.h"
#include "matpack_math.h"

#include "blas.h"
#include "lapack.h"

#include <algorithm>

void mult(MatrixView A, const ConstMatrixView& B, const ConstMatrixView& C) {
  // Check dimensions:
  ARTS_ASSERT(A.nrows() == B.nrows());
  ARTS_ASSERT(A.ncols() == C.ncols());
  ARTS_ASSERT(B.ncols() == C.nrows());

  // Catch trivial case if one of the matrices is empty.
  if ((B.nrows() == 0) || (B.ncols() == 0) || (C.ncols() == 0)) return;

  // Matrices B and C must be continuous in at least on dimension,  C
  // must be continuous along the second dimension.
  if (((B.stride(0) == 1) || (B.stride(1) == 1)) &&
      ((C.stride(0) == 1) || (C.stride(1) == 1)) &&
      (A.stride(1) == 1)) {
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
      lda = (int)C.stride(1);
    } else {
      transa = 'N';
      lda = (int)C.stride(0);
    }

    // Check if B (arts) is transposed.
    if (B.stride(0) == 1) {
      transb = 'T';
      ldb = (int)B.stride(1);
    } else {
      transb = 'N';
      ldb = (int)B.stride(0);
    }

    // In case B (arts) has only one column, column and row stride are 1.
    // We therefore need to set ldb to k, since dgemm_ requires lda to be at
    // least k / m if A is non-transposed / transposed.
    if ((B.stride(1) == 1) && (B.stride(0) == 1)) {
      transb = 'N';
      ldb = k;
    }

    // The same holds for C (arts).
    if ((C.stride(1) == 1) && (C.stride(0) == 1)) {
      transa = 'N';
      lda = m;
    }

    ldc = (int)A.stride(0);
    // The same holds for A (arts).
    if ((A.stride(1) == 1) && (A.stride(0) == 1)) {
      ldc = m;
    }
    double alpha = 1.0, beta = 0.0;

    dgemm_(&transa,
           &transb,
           &m,
           &n,
           &k,
           &alpha,
           C.unsafe_data_handle(),
           &lda,
           B.unsafe_data_handle(),
           &ldb,
           &beta,
           A.unsafe_data_handle(),
           &ldc);

  } else {
    matpack::eigen::as_eigen(A).noalias() = B * C;
  }
}

void mult(ComplexMatrixView A, const ConstComplexMatrixView& B, const ConstComplexMatrixView& C) {
  // Check dimensions:
  ARTS_ASSERT(A.nrows() == B.nrows());
  ARTS_ASSERT(A.ncols() == C.ncols());
  ARTS_ASSERT(B.ncols() == C.nrows());

  // Catch trivial case if one of the matrices is empty.
  if ((B.nrows() == 0) || (B.ncols() == 0) || (C.ncols() == 0)) return;

  // Matrices B and C must be continuous in at least on dimension,  C
  // must be continuous along the second dimension.
  if (((B.stride(0) == 1) || (B.stride(1) == 1)) &&
      ((C.stride(0) == 1) || (C.stride(1) == 1)) &&
      (A.stride(1) == 1)) {
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
      lda = (int)C.stride(1);
    } else {
      transa = 'N';
      lda = (int)C.stride(0);
    }

    // Check if B (arts) is transposed.
    if (B.stride(0) == 1) {
      transb = 'T';
      ldb = (int)B.stride(1);
    } else {
      transb = 'N';
      ldb = (int)B.stride(0);
    }

    // In case B (arts) has only one column, column and row stride are 1.
    // We therefore need to set ldb to k, since dgemm_ requires lda to be at
    // least k / m if A is non-transposed / transposed.
    if ((B.stride(1) == 1) && (B.stride(0) == 1)) {
      transb = 'N';
      ldb = k;
    }

    // The same holds for C (arts).
    if ((C.stride(1) == 1) && (C.stride(0) == 1)) {
      transa = 'N';
      lda = m;
    }

    ldc = (int)A.stride(0);
    // The same holds for A (arts).
    if ((A.stride(1) == 1) && (A.stride(0) == 1)) {
      ldc = m;
    }
    std::complex<double> alpha = 1.0, beta = 0.0;

    zgemm_(&transa,
           &transb,
           &m,
           &n,
           &k,
           &alpha,
           C.unsafe_data_handle(),
           &lda,
           B.unsafe_data_handle(),
           &ldb,
           &beta,
           A.unsafe_data_handle(),
           &ldc);

  } else {
    matpack::eigen::as_eigen(A).noalias() = B * C;
  }
}

void mult(VectorView y, const ConstMatrixView &M, const ConstVectorView &x) {
  ARTS_ASSERT(y.nelem() == M.nrows());
  ARTS_ASSERT(M.ncols() == x.nelem());
  ARTS_ASSERT(not M.empty());

  if ((M.stride(0) == 1) || (M.stride(1) == 1)) {
    char trans;
    int m, n;
    double zero = 0.0;
    double one = 1.0;
    int LDA, incx, incy;

    if (M.stride(1) != 1) {
      trans = 'n';
      m = (int)M.nrows();
      n = (int)M.ncols();
      LDA = (int)M.stride(1);
    } else {
      trans = 't';
      m = (int)M.ncols();
      n = (int)M.nrows();
      LDA = (int)M.stride(0);
      if (M.stride(0) == 1)
        LDA = m;
    }

    incx = (int)x.stride(0);
    incy = (int)y.stride(0);

    double *mstart = M.unsafe_data_handle();
    double *ystart = y.unsafe_data_handle();
    double *xstart = x.unsafe_data_handle();

    dgemv_(&trans, &m, &n, &one, mstart, &LDA, xstart, &incx, &zero, ystart,
           &incy);
  } else {
    matpack::eigen::as_eigen(y).noalias() = M * x;
  }
}

void mult(ComplexVectorView y, const ConstComplexMatrixView &M, const ConstComplexVectorView &x) {
  ARTS_ASSERT(y.nelem() == M.nrows());
  ARTS_ASSERT(M.ncols() == x.nelem());
  ARTS_ASSERT(not M.empty());

  if ((M.stride(0) == 1) || (M.stride(1) == 1)) {
    char trans;
    int m, n;
    std::complex<double> zero = 0.0;
    std::complex<double> one = 1.0;
    int LDA, incx, incy;

    if (M.stride(1) != 1) {
      trans = 'n';
      m = (int)M.nrows();
      n = (int)M.ncols();
      LDA = (int)M.stride(1);
    } else {
      trans = 't';
      m = (int)M.ncols();
      n = (int)M.nrows();
      LDA = (int)M.stride(0);
      if (M.stride(0) == 1)
        LDA = m;
    }

    incx = (int)x.stride(0);
    incy = (int)y.stride(0);

    auto *mstart = M.unsafe_data_handle();
    auto *ystart = y.unsafe_data_handle();
    auto *xstart = x.unsafe_data_handle();

    zgemv_(&trans, &m, &n, &one, mstart, &LDA, xstart, &incx, &zero, ystart,
           &incy);
  } else {
    matpack::eigen::as_eigen(y).noalias() = M * x;
  }
}

Vector uniform_grid(Numeric x0, Index N, Numeric dx) {
  Vector out(N);
  std::generate(out.begin(), out.end(), [x=x0, dx]() mutable {auto xd=x; x+=dx; return xd;});
  return out;
}

ComplexVector uniform_grid(Complex x0, Index N, Complex dx) {
  ComplexVector out(N);
  std::generate(out.begin(), out.end(), [x=x0, dx]() mutable {auto xd=x; x+=dx; return xd;});
  return out;
}

Vector diagonal(const ConstMatrixView& A) {
  using namespace matpack::eigen;
  return Vector{as_eigen(A).diagonal()};
}

void cross3(VectorView c, const ConstVectorView& a, const ConstVectorView& b) {
  ARTS_ASSERT(a.nelem() == 3, a.nelem(), " vs 3");
  ARTS_ASSERT(b.nelem() == 3, b.nelem(), " vs 3");
  ARTS_ASSERT(c.nelem() == 3, c.nelem(), " vs 3");

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

void cross3(ComplexVectorView c, const ConstComplexVectorView& a, const ConstComplexVectorView& b) {
  ARTS_ASSERT(a.nelem() == 3, a.nelem(), " vs 3");
  ARTS_ASSERT(b.nelem() == 3, b.nelem(), " vs 3");
  ARTS_ASSERT(c.nelem() == 3, c.nelem(), " vs 3");

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}
