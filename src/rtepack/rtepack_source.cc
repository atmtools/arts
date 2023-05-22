#include "rtepack_source.h"
#include "rtepack_multitype.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack::source  {
void level_nlte_and_scattering_and_sun(
    stokvec_vector_view J, stokvec_matrix_view dJ,
    const stokvec_vector_const_view &J_add, const propmat_vector_const_view &K,
    const stokvec_vector_const_view &a, const stokvec_vector_const_view &S,
    const propmat_matrix_const_view &dK, const stokvec_matrix_const_view &da,
    const stokvec_matrix_const_view &dS, const ExhaustiveConstVectorView &B,
    const ExhaustiveConstMatrixView &dB) {
  const Index N = J.nelem();
  ARTS_ASSERT(N == dJ.ncols())
  ARTS_ASSERT(N == J_add.nelem())
  ARTS_ASSERT(N == K.nelem())
  ARTS_ASSERT(N == a.nelem())
  ARTS_ASSERT(N == S.nelem())
  ARTS_ASSERT(N == dK.ncols())
  ARTS_ASSERT(N == da.ncols())
  ARTS_ASSERT(N == dS.ncols())
  ARTS_ASSERT(N == B.nelem())
  ARTS_ASSERT(N == dB.ncols())

  const Index M = dJ.nrows();
  ARTS_ASSERT(M == dK.nrows())
  ARTS_ASSERT(M == da.nrows())
  ARTS_ASSERT(M == dS.nrows())
  ARTS_ASSERT(M == dB.nrows())

  for (Index i = 0; i < N; i++) {
    if (K[i].is_rotational()) {
      J[i] = 0.0;
      dJ(joker, i) = 0.0;
    } else {
      const auto invK = inv(K[i]);
      const auto src = a[i] * B[i] + S[i] + J_add[i];
      J[i] = invK * src;

      // TODO: Add jacobians dJ_add of additional source
      for (Index j = 0; j < M; j++) {
        const auto dsrc = dS(j, i) + da(j, i) * B[i] + a[i] * dB(j, i);
        dJ(j, i) = 0.5 * (invK * (dsrc - dK(j, i) * J[i]));
      }
    }
  }
}

void level_nlte_and_scattering(
    stokvec_vector_view J, stokvec_matrix_view dJ,
    const propmat_vector_const_view &K,
    const stokvec_vector_const_view &a, const stokvec_vector_const_view &S,
    const propmat_matrix_const_view &dK, const stokvec_matrix_const_view &da,
    const stokvec_matrix_const_view &dS, const ExhaustiveConstVectorView &B,
    const ExhaustiveConstMatrixView &dB) {
  const Index N = J.nelem();
  ARTS_ASSERT(N == dJ.ncols())
  ARTS_ASSERT(N == K.nelem())
  ARTS_ASSERT(N == a.nelem())
  ARTS_ASSERT(N == S.nelem())
  ARTS_ASSERT(N == dK.ncols())
  ARTS_ASSERT(N == da.ncols())
  ARTS_ASSERT(N == dS.ncols())
  ARTS_ASSERT(N == B.nelem())
  ARTS_ASSERT(N == dB.ncols())

  const Index M = dJ.nrows();
  ARTS_ASSERT(M == dK.nrows())
  ARTS_ASSERT(M == da.nrows())
  ARTS_ASSERT(M == dS.nrows())
  ARTS_ASSERT(M == dB.nrows())

  for (Index i = 0; i < N; i++) {
    if (K[i].is_rotational()) {
      J[i] = 0.0;
      dJ(joker, i) = 0.0;
    } else {
      const auto invK = inv(K[i]);
      const auto src = a[i] * B[i] + S[i];
      J[i] = invK * src;

      for (Index j = 0; j < M; j++) {
        const auto dsrc = dS(j, i) + da(j, i) * B[i] + a[i] * dB(j, i);
        dJ(j, i) = 0.5 * (invK * (dsrc - dK(j, i) * J[i]));
      }
    }
  }
}

void level_nlte_and_scattering(stokvec_vector_view J,
                               const propmat_vector_const_view &K,
                               const stokvec_vector_const_view &a,
                               const stokvec_vector_const_view &S,
                               const ExhaustiveConstVectorView &B) {
  const Index N = J.nelem();
  ARTS_ASSERT(N == K.nelem())
  ARTS_ASSERT(N == a.nelem())
  ARTS_ASSERT(N == S.nelem())
  ARTS_ASSERT(N == B.nelem())

  for (Index i = 0; i < N; i++) {
    J[i] = K[i].is_rotational() ? stokvec{0, 0, 0, 0}
                                : inv(K[i]) * (a[i] * B[i] + S[i]);
  }
}

void level_nlte(stokvec_vector_view J, stokvec_matrix_view dJ,
                const propmat_vector_const_view &K,
                const stokvec_vector_const_view &S,
                const propmat_matrix_const_view &dK,
                const stokvec_matrix_const_view &dS,
                const ExhaustiveConstVectorView &B,
                const ExhaustiveConstMatrixView &dB) {
  const Index N = J.nelem();
  ARTS_ASSERT(N == dJ.ncols())
  ARTS_ASSERT(N == K.nelem())
  ARTS_ASSERT(N == S.nelem())
  ARTS_ASSERT(N == dK.ncols())
  ARTS_ASSERT(N == dS.ncols())
  ARTS_ASSERT(N == B.nelem())
  ARTS_ASSERT(N == dB.ncols())

  const Index M = dJ.nrows();
  ARTS_ASSERT(M == dK.nrows())
  ARTS_ASSERT(M == dS.nrows())
  ARTS_ASSERT(M == dB.nrows())

  for (Index i = 0; i < N; i++) {
    if (K[i].is_rotational()) {
      J[i] = 0.0;
      dJ(joker, i) = 0.0;
    } else {
      const auto invK = inv(K[i]);
      const auto b = stokvec{B[i], 0, 0, 0};
      J[i] = b + invK * S[i];

      for (Index j = 0; j < M; j++) {
        const auto db = stokvec{dB(j, i), 0, 0, 0};
        const auto JB = stokvec{J[i].I() + B[i], J[i].Q(), J[i].U(), J[i].V()};
        dJ(j, i) = 0.5 * (db + invK * (dS(j, i) - dK(j, i) * JB));
      }
    }
  }
}
} // namespace rtepack::source
