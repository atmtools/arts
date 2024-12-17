#pragma once

#include <algorithm>

#include "rtepack_multitype.h"

namespace rtepack::source {
constexpr stokvec level_lte(Numeric B) { return stokvec{B, 0, 0, 0}; }

constexpr stokvec level_lte(stokvec_vector_view dj,
                            Numeric B,
                            const ExhaustiveConstVectorView &dB) {
  ARTS_ASSERT(dj.size() == dB.size())
  std::transform(dB.elem_begin(), dB.elem_end(), dj.elem_begin(), [](auto &db) {
    return stokvec{db};
  });
  return level_lte(B);
}

constexpr stokvec level_nlte(Numeric B, const propmat &k, const stokvec &n) {
  return inv(k) * n + stokvec{B, 0, 0, 0};
}

constexpr stokvec level_nlte(stokvec_vector_view dj,
                             Numeric B,
                             const ExhaustiveConstVectorView &dB,
                             const propmat &k,
                             const propmat_vector_view &dk,
                             const stokvec &n,
                             const stokvec_vector_view &dn) {
  const Index N = dj.size();
  ARTS_ASSERT(N == dB.size())
  ARTS_ASSERT(N == dk.size())
  ARTS_ASSERT(N == dn.size())

  const auto inv_k = inv(k);
  const auto a = absvec(k);
  stokvec out = inv_k * (a * B + n);

  for (Index i = 0; i < N; i++) {
    dj[i] = inv_k * (dk[i] * out + absvec(dk[i]) * B + a * dB[i] + dn[i]);
  }

  return out;
}

constexpr stokvec unpolarized_basis_vector() { return stokvec{1}; }

constexpr stokvec unpolarized_basis_vector(stokvec_vector_view dj) {
  std::fill(dj.elem_begin(), dj.elem_end(), unpolarized_basis_vector());
  return unpolarized_basis_vector();
}

void level_nlte_and_scattering_and_sun(stokvec_vector_view J,
                                       stokvec_matrix_view dJ,
                                       const stokvec_vector_const_view &J_add,
                                       const propmat_vector_const_view &K,
                                       const stokvec_vector_const_view &a,
                                       const stokvec_vector_const_view &S,
                                       const propmat_matrix_const_view &dK,
                                       const stokvec_matrix_const_view &da,
                                       const stokvec_matrix_const_view &dS,
                                       const ExhaustiveConstVectorView &B,
                                       const ExhaustiveConstMatrixView &dB);

void level_nlte_and_scattering(stokvec_vector_view J,
                               stokvec_matrix_view dJ,
                               const propmat_vector_const_view &K,
                               const stokvec_vector_const_view &a,
                               const stokvec_vector_const_view &S,
                               const propmat_matrix_const_view &dK,
                               const stokvec_matrix_const_view &da,
                               const stokvec_matrix_const_view &dS,
                               const ExhaustiveConstVectorView &B,
                               const ExhaustiveConstMatrixView &dB);

void level_nlte_and_scattering(stokvec_vector_view J,
                               const propmat_vector_const_view &K,
                               const stokvec_vector_const_view &a,
                               const stokvec_vector_const_view &S,
                               const ExhaustiveConstVectorView &B);

void level_nlte(stokvec_vector_view J,
                stokvec_matrix_view dJ,
                const propmat_vector_const_view &K,
                const stokvec_vector_const_view &S,
                const propmat_matrix_const_view &dK,
                const stokvec_matrix_const_view &dS,
                const ExhaustiveConstVectorView &f,
                const Numeric& t,
                const Index& it);
}  // namespace rtepack::source
