#include "matpack_band_matrix.h"

extern "C" void dgbsv_(int* N,
                       int* KL,
                       int* KU,
                       int* NRHS,
                       double* AB,
                       int* LDAB,
                       int* IPIV,
                       double* B,
                       int* LDB,
                       int* INFO);

namespace matpack {
Matrix band_matrix::mat(Index KL, Index KU, Index N) {
  return Matrix(N, 2 * KL + KU + 1);
}

Index band_matrix::end_row(Index j) const {
  return std::min<Index>(M, j + KL + 1);
}

Index band_matrix::start_row(Index j) const {
  return std::max<Index>(0, j - KU);
}

band_matrix::band_matrix(Index ku, Index kl, Index m, Index n)
    : KU(ku), KL(kl), M(m), N(n), AB(mat(KL, KU, N)), ipiv(N) {}

band_matrix::band_matrix(const Matrix& ab)
    : KU(0), KL(0), M(ab.nrows()), N(ab.ncols()), AB(0, 0), ipiv(N) {
  for (Index i = 0; i < M; ++i) {
    for (Index j = 0; j < N; ++j) {
      if (i == j) continue;

      if (ab(i, j) != 0) {
        if (j > i)
          KU = std::max(j - i, KU);
        else if (j < i)
          KL = std::max(i - j, KL);
      }
    }
  }

  AB = mat(KL, KU, N);

  for (Index j = 0; j < N; j++) {
    const Index maxi = end_row(j);
    for (Index i = start_row(j); i < maxi; i++) {
      this->operator()(i, j) = ab(i, j);
    }
  }
}

Numeric& band_matrix::operator()(Index i, Index j) {
  ARTS_ASSERT(
      i >= start_row(j), "Out of lower bound. limit: {} vs {}", start_row(j), i)
  ARTS_ASSERT(
      i < end_row(j), "Out of upper bound. limit: {} vs {}", end_row(j), i)
  return AB(j, KU + KL + i - j);
}

int band_matrix::solve(Vector& bx) {
  int n = static_cast<int>(N);
  int kl = static_cast<int>(KL);
  int ku = static_cast<int>(KU);
  int nrhs = 1;
  int ldab = 2 * kl + ku + 1;
  int info = 0;

  dgbsv_(&n,
         &kl,
         &ku,
         &nrhs,
         AB.data_handle(),
         &ldab,
         ipiv.data(),
         bx.data_handle(),
         &n,
         &info);

  return info;
}
}  // namespace matpack