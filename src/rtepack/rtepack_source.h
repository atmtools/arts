#pragma once

#include "debug.h"
#include "matpack_view.h"

#include "rtepack_multitype.h"
#include <algorithm>

namespace rtepack::source {
constexpr stokvec level_lte(Numeric B) { return B; }

constexpr stokvec level_lte(stokvec_vector_view &dj, Numeric B,
                            const ExhaustiveConstVectorView &dB) {
  ARTS_ASSERT(dj.nelem() == dB.nelem())
  std::transform(dB.elem_begin(), dB.elem_end(), dj.elem_begin(),
                 [](auto &db) -> stokvec { return db; });
  return level_lte(B);
}

constexpr stokvec level_nlte(Numeric B, const propmat &k, const stokvec &n) {
  return inv(k) * (absvec(k) * B + n);
}

constexpr stokvec level_nlte(stokvec_vector_view &dj, Numeric B,
                             const ExhaustiveConstVectorView &dB,
                             const propmat &k, const propmat_vector_view &dk,
                             const stokvec &n, const stokvec_vector_view &dn) {
  const Index N = dj.nelem();
  ARTS_ASSERT(N == dB.nelem())
  ARTS_ASSERT(N == dk.nelem())
  ARTS_ASSERT(N == dn.nelem())

  const auto inv_k = inv(k);
  const auto a = absvec(k);
  stokvec out = inv_k * (a * B + n);

  for (Index i = 0; i < N; i++) {
    dj[i] = inv_k * (dk[i] * out + absvec(dk[i]) * B + a * dB[i] + dn[i]);
  }

  return out;
}

constexpr stokvec unpolarized_basis_vector() { return 1; }

constexpr stokvec unpolarized_basis_vector(stokvec_vector_view &dj) {
  std::fill(dj.elem_begin(), dj.elem_end(), stokvec{});
  return unpolarized_basis_vector();
}
} // namespace rtepack::source
