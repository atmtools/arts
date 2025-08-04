#include "rtepack_source.h"

#include <debug.h>
#include <physics_funcs.h>

#include "rtepack_multitype.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack::source {

stokvec level_lte(Numeric B) { return stokvec{B, 0, 0, 0}; }

stokvec level_lte(stokvec_vector_view dj,
                  Numeric B,
                  const ConstVectorView &dB) {
  assert(dj.size() == dB.size());

  std::transform(dB.elem_begin(), dB.elem_end(), dj.elem_begin(), [](auto &db) {
    return stokvec{db};
  });
  return level_lte(B);
}

stokvec level_nlte(Numeric B, const propmat &k, const stokvec &n) {
  return inv(k) * n + stokvec{B, 0, 0, 0};
}

stokvec level_nlte(stokvec_vector_view dj,
                   Numeric B,
                   const ConstVectorView &dB,
                   const propmat &k,
                   const propmat_vector_view &dk,
                   const stokvec &n,
                   const stokvec_vector_view &dn) {
  const Size N = dj.size();
  assert(N == dB.size());
  assert(N == dk.size());
  assert(N == dn.size());

  const auto inv_k = inv(k);
  const auto a     = absvec(k);
  stokvec out      = inv_k * (a * B + n);

  for (Size i = 0; i < N; i++) {
    dj[i] = inv_k * (dk[i] * out + absvec(dk[i]) * B + a * dB[i] + dn[i]);
  }

  return out;
}

stokvec unpolarized_basis_vector() { return stokvec{1}; }

stokvec unpolarized_basis_vector(stokvec_vector_view dj) {
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
                                       const ConstVectorView &B,
                                       const ConstMatrixView &dB) {
  const Size N = J.size();
  assert(N == static_cast<Size>(dJ.ncols()));
  assert(N == J_add.size());
  assert(N == K.size());
  assert(N == a.size());
  assert(N == S.size());
  assert(N == static_cast<Size>(dK.ncols()));
  assert(N == static_cast<Size>(da.ncols()));
  assert(N == static_cast<Size>(dS.ncols()));
  assert(N == B.size());
  assert(N == static_cast<Size>(dB.ncols()));

  const Index M = dJ.nrows();
  assert(M == dK.nrows());
  assert(M == da.nrows());
  assert(M == dS.nrows());
  assert(M == dB.nrows());

  for (Size i = 0; i < N; i++) {
    if (K[i].is_rotational()) {
      J[i]         = 0.0;
      dJ[joker, i] = 0.0;
    } else {
      const auto invK = inv(K[i]);
      const auto src  = a[i] * B[i] + S[i] + J_add[i];
      J[i]            = invK * src;

      // TODO: Add jacobians dJ_add of additional source
      for (Index j = 0; j < M; j++) {
        const auto dsrc = dS[j, i] + da[j, i] * B[i] + a[i] * dB[j, i];
        dJ[j, i]        = 0.5 * (invK * (dsrc - dK[j, i] * J[i]));
      }
    }
  }
}

void level_nlte_and_scattering(stokvec_vector_view J,
                               stokvec_matrix_view dJ,
                               const propmat_vector_const_view &K,
                               const stokvec_vector_const_view &a,
                               const stokvec_vector_const_view &S,
                               const propmat_matrix_const_view &dK,
                               const stokvec_matrix_const_view &da,
                               const stokvec_matrix_const_view &dS,
                               const ConstVectorView &B,
                               const ConstMatrixView &dB) {
  const Size N = J.size();
  assert(N == static_cast<Size>(dJ.ncols()));
  assert(N == K.size());
  assert(N == a.size());
  assert(N == S.size());
  assert(N == static_cast<Size>(dK.ncols()));
  assert(N == static_cast<Size>(da.ncols()));
  assert(N == static_cast<Size>(dS.ncols()));
  assert(N == B.size());
  assert(N == static_cast<Size>(dB.ncols()));

  const Index M = dJ.nrows();
  assert(M == dK.nrows());
  assert(M == da.nrows());
  assert(M == dS.nrows());
  assert(M == dB.nrows());

  for (Size i = 0; i < N; i++) {
    if (K[i].is_rotational()) {
      J[i]         = 0.0;
      dJ[joker, i] = 0.0;
    } else {
      const auto invK = inv(K[i]);
      const auto src  = a[i] * B[i] + S[i];
      J[i]            = invK * src;

      for (Index j = 0; j < M; j++) {
        const auto dsrc = dS[j, i] + da[j, i] * B[i] + a[i] * dB[j, i];
        dJ[j, i]        = 0.5 * (invK * (dsrc - dK[j, i] * J[i]));
      }
    }
  }
}

void level_nlte_and_scattering(stokvec_vector_view J,
                               const propmat_vector_const_view &K,
                               const stokvec_vector_const_view &a,
                               const stokvec_vector_const_view &S,
                               const ConstVectorView &B) {
  const Size N = J.size();
  assert(N == K.size());
  assert(N == a.size());
  assert(N == S.size());
  assert(N == B.size());

  for (Size i = 0; i < N; i++) {
    J[i] = K[i].is_rotational() ? stokvec{0, 0, 0, 0}
                                : inv(K[i]) * (a[i] * B[i] + S[i]);
  }
}

void level_nlte(stokvec_vector_view J,
                stokvec_matrix_view dJ,
                const propmat_vector_const_view &K,
                const stokvec_vector_const_view &S,
                const propmat_matrix_const_view &dK,
                const stokvec_matrix_const_view &dS,
                const ConstVectorView &f,
                const Numeric &t,
                const Index &it) {
  const Size N = J.size();
  assert(N == static_cast<Size>(dJ.ncols()));
  assert(N == K.size());
  assert(N == S.size());
  assert(N == static_cast<Size>(dK.ncols()));
  assert(N == static_cast<Size>(dS.ncols()));
  assert(N == f.size());

  const Index M = dJ.nrows();
  assert(M == dK.nrows());
  assert(M == dS.nrows());
  assert(M > it);

  for (Size i = 0; i < N; i++) {
    if (K[i].is_rotational()) {
      J[i]         = 0.0;
      dJ[joker, i] = 0.0;
    } else {
      const auto b    = stokvec{planck(f[i], t), 0, 0, 0};
      const auto invK = inv(K[i]);
      J[i]            = b + invK * S[i];

      for (Index j = 0; j < M; j++) {
        dJ[j, i]  = {it == j ? dplanck_dt(f[i], t) : 0, 0, 0, 0};
        dJ[j, i] -= invK * (dK[j, i] * invK * S[i] - dS[j, i]);
      }
    }
  }
}
}  // namespace rtepack::source
