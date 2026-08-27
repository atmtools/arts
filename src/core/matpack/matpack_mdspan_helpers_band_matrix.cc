#include "matpack_mdspan_helpers_band_matrix.h"

#include <time_report.h>

extern "C" void dgbsv_(
    int* N, int* KL, int* KU, int* NRHS, double* AB, int* LDAB, int* IPIV, double* B, int* LDB, int* INFO);
extern "C" void zgbsv_(
    int* N, int* KL, int* KU, int* NRHS, Complex* AB, int* LDAB, int* IPIV, Complex* B, int* LDB, int* INFO);

namespace matpack {
template <> int basic_band_matrix<Numeric>::solve(Vector& bx) {
  ARTS_TIME_REPORT

  int n    = static_cast<int>(N);
  int kl   = static_cast<int>(KL);
  int ku   = static_cast<int>(KU);
  int nrhs = 1;
  int ldab = 2 * kl + ku + 1;
  int info = 0;

  dgbsv_(&n, &kl, &ku, &nrhs, AB.data_handle(), &ldab, ipiv.data(), bx.data_handle(), &n, &info);

  return info;
}

template <> int basic_band_matrix<Complex>::solve(ComplexVector& bx) {
  ARTS_TIME_REPORT

  int n    = static_cast<int>(N);
  int kl   = static_cast<int>(KL);
  int ku   = static_cast<int>(KU);
  int nrhs = 1;
  int ldab = 2 * kl + ku + 1;
  int info = 0;

  zgbsv_(&n, &kl, &ku, &nrhs, AB.data_handle(), &ldab, ipiv.data(), bx.data_handle(), &n, &info);
  return info;
}
}  // namespace matpack
