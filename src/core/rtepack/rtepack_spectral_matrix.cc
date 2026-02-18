#include "rtepack_spectral_matrix.h"

#include <Faddeeva.hh>

namespace rtepack {
specmat dawson(const specmat &A) {
  using Faddeeva::Dawson;

  return specmat{Dawson(A[0, 0]),
                 Dawson(A[0, 1]),
                 Dawson(A[0, 2]),
                 Dawson(A[0, 3]),
                 Dawson(A[1, 0]),
                 Dawson(A[1, 1]),
                 Dawson(A[1, 2]),
                 Dawson(A[1, 3]),
                 Dawson(A[2, 0]),
                 Dawson(A[2, 1]),
                 Dawson(A[2, 2]),
                 Dawson(A[2, 3]),
                 Dawson(A[3, 0]),
                 Dawson(A[3, 1]),
                 Dawson(A[3, 2]),
                 Dawson(A[3, 3])};
}

specmat erfcx(const specmat &m) {
  using Faddeeva::erfcx;
  return specmat{erfcx(m[0, 0]),
                 erfcx(m[0, 1]),
                 erfcx(m[0, 2]),
                 erfcx(m[0, 3]),
                 erfcx(m[1, 0]),
                 erfcx(m[1, 1]),
                 erfcx(m[1, 2]),
                 erfcx(m[1, 3]),
                 erfcx(m[2, 0]),
                 erfcx(m[2, 1]),
                 erfcx(m[2, 2]),
                 erfcx(m[2, 3]),
                 erfcx(m[3, 0]),
                 erfcx(m[3, 1]),
                 erfcx(m[3, 2]),
                 erfcx(m[3, 3])};
}
}  // namespace rtepack
