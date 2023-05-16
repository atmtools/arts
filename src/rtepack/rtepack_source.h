#pragma once

#include "debug.h"
#include "matpack_view.h"

#include "rtepack/rtepack_mueller.h"
#include "rtepack_multitype.h"
#include <algorithm>

namespace rtepack::source {
constexpr stokes level_lte(Numeric B) { return B; }

constexpr stokes level_lte(stokes_view &dj, Numeric B,
                           const ExhaustiveConstVectorView &dB) {
  ARTS_ASSERT(dj.nelem() == dB.nelem())
  std::transform(dB.elem_begin(), dB.elem_end(), dj.elem_begin(),
                 [](auto &db) -> stokes {
                   return db;
                 });
  return level_lte(B);
}

constexpr stokes level_nlte(Numeric B, const propmat &k, const stokes &n) {
  return inv(k) * (absvec(k) * B + n);
}

constexpr stokes level_nlte(stokes_view &dj, Numeric B,
                            const ExhaustiveConstVectorView &dB,
                            const propmat &k, const propmat_view &dk,
                            const stokes &n, const stokes_view &dn) {
  const Index N = dj.nelem();
  ARTS_ASSERT(N == dB.nelem())
  ARTS_ASSERT(N == dk.nelem())
  ARTS_ASSERT(N == dn.nelem())

  const auto inv_k = inv(k);
  const auto a = absvec(k);
  stokes out = inv_k * (a * B + n);

  for (Index i = 0; i < N; i++) {
    dj[i] = inv_k * (dk[i] * out + absvec(dk[i]) * B + a * dB[i] + dn[i]);
  }

  return out;
}

constexpr stokes unpolarized_basis_vector() { return 1; }

constexpr stokes unpolarized_basis_vector(stokes_view &dj) {
  std::fill(dj.elem_begin(), dj.elem_end(), stokes{});
  return unpolarized_basis_vector();
}
} // namespace rtepack::source
