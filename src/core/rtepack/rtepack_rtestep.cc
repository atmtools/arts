#include "rtepack_rtestep.h"

#include <array_algo.h>
#include <arts_omp.h>
#include <physics_funcs.h>

#include "rtepack_multitype.h"
#include "rtepack_transmission.h"

namespace rtepack {
void two_level_linear_emission_step_by_step_full(
    stokvec_vector &I,
    std::vector<stokvec_matrix> &dI,
    const std::vector<muelmat_vector> &Ts,
    const std::vector<muelmat_vector> &Pi,
    const std::vector<muelmat_tensor3> &dTs,
    const std::vector<stokvec_vector> &Js,
    const std::vector<stokvec_matrix> &dJs,
    const stokvec_vector &I0) {
  const Size nv = I0.size();
  const Size N  = Ts.size();

  I = I0;
  dI.resize(N);

  if (N == 0) return;

  const Size nq = dJs.front().nrows();

  ARTS_USER_ERROR_IF(not arr::same_size(Ts, Pi, dTs, Js, dJs),
                     R"(
  Mismatched number of levels ({}) in Ts:s, Pi:s, dTs:s, Js:s, and dJs:s:
  Ts  size: {}
  Pi  size: {}
  dTs size: {}
  Js  size: {}
  dJs size: {}
)",
                     N,
                     Ts.size(),
                     Pi.size(),
                     dTs.size(),
                     Js.size(),
                     dJs.size());

  ARTS_USER_ERROR_IF(
      not all_same_shape(I0, Ts, Pi, Js),
      R"(Must have same number of frequency elements ({}) in all Ts:s, Pi:s, and Js:s)",
      nv);

  ARTS_USER_ERROR_IF(
      not all_same_shape({nq, nv}, dJs),
      R"(Must have same number of elements ({}, {}) in all dJs:s)",
      nq,
      nv);

  ARTS_USER_ERROR_IF(
      not all_same_shape({2, nq, nv}, dTs),
      R"(Must have same number of elements ({}, {}, {}) in all dTs:s)",
      2,
      nq,
      nv);

  for (auto &x : dI) {
    x.resize(nq, nv);
    x = 0.0;
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    stokvec &Iv = I[iv];
    for (Size i = N - 2; i < N; i--) {
      const muelmat &T = Ts[i + 1][iv];
      const stokvec Jv = avg(Js[i][iv], Js[i + 1][iv]);

      Iv -= Jv;

      if (nq) {
        auto dI0 = dI[i][joker, iv];
        auto dI1 = dI[i + 1][joker, iv];

        const auto dT0 = dTs[i][0, joker, iv];
        const auto dJ0 = dJs[i][joker, iv];

        const auto dT1 = dTs[i + 1][1, joker, iv];
        const auto dJ1 = dJs[i + 1][joker, iv];

        const auto &P = Pi[i][iv];

        for (Size iq = 0; iq < nq; iq++) {
          dI0[iq] += P * (dT0[iq] * Iv + avg(dJ0[iq], -(T * dJ0[iq])));
          dI1[iq] += P * (dT1[iq] * Iv + avg(dJ1[iq], -(T * dJ1[iq])));
        }
      }

      Iv = T * Iv + Jv;
    }
  }
}

void two_level_linear_in_J_step_by_step_full(
    stokvec_vector &I,
    std::vector<stokvec_matrix> &dI,
    const std::vector<muelmat_vector> &Ts,
    const std::vector<muelmat_vector> &L0s,
    const std::vector<muelmat_vector> &Pi,
    const std::vector<muelmat_tensor3> &dTs,
    const std::vector<muelmat_tensor3> &dL0s,
    const std::vector<stokvec_vector> &Js,
    const std::vector<stokvec_matrix> &dJs,
    const stokvec_vector &I0) {
  const Size nv = I0.size();
  const Size N  = Ts.size();

  I = I0;
  dI.resize(N);

  if (N == 0) return;

  const Size nq = dJs.front().nrows();

  ARTS_USER_ERROR_IF(not arr::same_size(Ts, L0s, Pi, dTs, dL0s, Js, dJs),
                     R"(
  Mismatched number of levels ({}) in Ts:s, L0s:s, Pi:s, dTs:s, dL0s:s, Js:s, and dJs:s:
  Ts   size: {}
  L0s  size: {}
  Pi   size: {}
  dTs  size: {}
  dL0s size: {}
  Js   size: {}
  dJs  size: {}
)",
                     N,
                     Ts.size(),
                     L0s.size(),
                     Pi.size(),
                     dTs.size(),
                     dL0s.size(),
                     Js.size(),
                     dJs.size());

  ARTS_USER_ERROR_IF(
      not all_same_shape(I0, Ts, L0s, Pi, Js),
      R"(Must have same number of frequency elements ({}) in all Ts:s, L0s:s, Pi:s, and Js:s)",
      nv);

  ARTS_USER_ERROR_IF(
      not all_same_shape({nq, nv}, dJs),
      R"(Must have same number of elements ({}, {}) in all dJs:s)",
      nq,
      nv);

  ARTS_USER_ERROR_IF(
      not all_same_shape({2, nq, nv}, dTs, dL0s),
      R"(Must have same number of elements ({}, {}, {}) in all dTs:s and dL0s:s)",
      2,
      nq,
      nv);

  for (auto &x : dI) {
    x.resize(nq, nv);
    x = 0.0;
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    stokvec &IM = I[iv];
    for (Size i = N - 2; i < N; i--) {
      const stokvec &JN = Js[i + 1][iv];
      const stokvec &JM = Js[i][iv];
      const muelmat &T  = Ts[i + 1][iv];
      const muelmat &L0 = L0s[i + 1][iv];

      const stokvec IMmJM = IM - JM;
      const stokvec JMmJN = JM - JN;

      if (nq) {
        auto dI0 = dI[i][joker, iv];
        auto dI1 = dI[i + 1][joker, iv];

        const auto dT0  = dTs[i][0, joker, iv];
        const auto dL00 = dL0s[i][0, joker, iv];
        const auto dJ0  = dJs[i][joker, iv];

        const auto dT1  = dTs[i + 1][1, joker, iv];
        const auto dL01 = dL0s[i + 1][1, joker, iv];
        const auto dJ1  = dJs[i + 1][joker, iv];

        const auto &P = Pi[i][iv];

        for (Size iq = 0; iq < nq; iq++) {
          dI0[iq] +=
              P * (dJ0[iq] - L0 * dJ0[iq] + dT0[iq] * IMmJM + dL00[iq] * JMmJN);
          dI1[iq] += P * (dT1[iq] * IMmJM + dL01[iq] * JMmJN + L0 * dJ1[iq] -
                          T * dJ1[iq]);
        }
      }

      IM = T * IMmJM + L0 * JMmJN + JN;
    }
  }
}

void two_level_linear_in_J_and_K_step_by_step_full(
    stokvec_vector &I,
    std::vector<stokvec_matrix> &dI,
    const std::vector<muelmat_vector> &Ts,
    const std::vector<muelmat_vector> &L0s,
    const std::vector<muelmat_vector> &L1s,
    const std::vector<muelmat_vector> &Pi,
    const std::vector<muelmat_tensor3> &dTs,
    const std::vector<muelmat_tensor3> &dL0s,
    const std::vector<muelmat_tensor3> &dL1s,
    const std::vector<stokvec_vector> &Js,
    const std::vector<stokvec_matrix> &dJs,
    const stokvec_vector &I0) {
  const Size nv = I0.size();
  const Size N  = Ts.size();

  I = I0;
  dI.resize(N);

  if (N == 0) return;

  const Size nq = dJs.front().nrows();

  ARTS_USER_ERROR_IF(
      not arr::same_size(Ts, L0s, L1s, Pi, dTs, dL0s, dL1s, Js, dJs),
      R"(
  Mismatched number of levels ({}) in Ts:s, L0s:s, L1s:s, Pi:s, dTs:s, dL0s:s, dL1s:s, Js:s, and dJs:s:
  Ts   size: {}
  L0s  size: {}
  L1s  size: {}
  Pi   size: {}
  dTs  size: {}
  dL0s size: {}
  dL1s size: {}
  Js   size: {}
  dJs  size: {}
)",
      N,
      Ts.size(),
      L0s.size(),
      L1s.size(),
      Pi.size(),
      dTs.size(),
      dL0s.size(),
      dL1s.size(),
      Js.size(),
      dJs.size());

  ARTS_USER_ERROR_IF(
      not all_same_shape(I0, Ts, L0s, L1s, Pi, Js),
      R"(Must have same number of frequency elements ({}) in all Ts:s, L0s:s, L1s:s, Pi:s, and Js:s)",
      nv);

  ARTS_USER_ERROR_IF(
      not all_same_shape({nq, nv}, dJs),
      R"(Must have same number of elements ({}, {}) in all dJs:s)",
      nq,
      nv);

  ARTS_USER_ERROR_IF(
      not all_same_shape({2, nq, nv}, dTs, dL0s, dL1s),
      R"(Must have same number of elements ({}, {}, {}) in all dTs:s, dL0s:s, and dL1s:s)",
      2,
      nq,
      nv);

  for (auto &x : dI) {
    x.resize(nq, nv);
    x = 0.0;
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    stokvec &IM = I[iv];
    for (Size i = N - 2; i < N; i--) {
      const stokvec &JN = Js[i + 1][iv];
      const stokvec &JM = Js[i][iv];
      const muelmat &T  = Ts[i + 1][iv];
      const muelmat &L0 = L0s[i + 1][iv];
      const muelmat &L1 = L1s[i + 1][iv];

      const stokvec IMmJM = IM - JM;
      const stokvec JMmJN = JM - JN;

      if (nq) {
        auto dI0 = dI[i][joker, iv];
        auto dI1 = dI[i + 1][joker, iv];

        const auto dT0  = dTs[i][0, joker, iv];
        const auto dL00 = dL0s[i][0, joker, iv];
        const auto dL10 = dL1s[i][0, joker, iv];
        const auto dJ0  = dJs[i][joker, iv];

        const auto dT1  = dTs[i + 1][1, joker, iv];
        const auto dL01 = dL0s[i + 1][1, joker, iv];
        const auto dL11 = dL1s[i + 1][1, joker, iv];
        const auto dJ1  = dJs[i + 1][joker, iv];

        const auto &P = Pi[i][iv];

        for (Size iq = 0; iq < nq; iq++) {
          dI0[iq] +=
              P * (dJ0[iq] - L0 * dJ0[iq] - L1 * dJ0[iq] + dT0[iq] * IMmJM +
                   dL00[iq] * JMmJN + dL10[iq] * JMmJN);

          dI1[iq] +=
              P * (dT1[iq] * IMmJM + dL01[iq] * JMmJN + dL11[iq] * JMmJN +
                   L0 * dJ1[iq] + L1 * dJ1[iq] - T * dJ1[iq]);
        }
      }

      IM = T * IMmJM + L0 * JMmJN + L1 * JMmJN + JN;
    }
  }
}

void two_level_linear_emission_step_by_step_full(
    std::vector<stokvec_vector> &Is,
    const std::vector<muelmat_vector> &Ts,
    const std::vector<stokvec_vector> &Js) {
  const Size N = Ts.size();

  if (N == 0) return;

  const Size nv = Is.front().size();

  ARTS_USER_ERROR_IF(not arr::same_size(Is, Ts, Js),
                     R"(Not same sizes:

Is.size() = {},
Ts.size() = {},
Js.size() = {}
)",
                     Is.size(),
                     Ts.size(),
                     Js.size())

  ARTS_USER_ERROR_IF(
      not arr::elemwise_same_size(Is, Ts, Js),
      "Not all elements have the same number of frequency elements")

  for (Size i = N - 2; i < N; i--) {
    stokvec_vector &I0       = Is[i + 1];
    stokvec_vector &I1       = Is[i];
    const stokvec_vector &J0 = Js[i + 1];
    const stokvec_vector &J1 = Js[i];
    const muelmat_vector &T  = Ts[i + 1];

#pragma omp parallel for if (not arts_omp_in_parallel())
    for (Size iv = 0; iv < nv; iv++) {
      const stokvec Jv = avg(J0[iv], J1[iv]);
      I1[iv]           = T[iv] * (I0[iv] - Jv) + Jv;
    }
  }
}

void two_level_linear_transmission_step(stokvec_vector &I,
                                        std::vector<stokvec_matrix> &dI,
                                        const std::vector<muelmat_vector> &Ts,
                                        const std::vector<muelmat_vector> &Pi,
                                        const std::vector<muelmat_tensor3> &dTs,
                                        const stokvec_vector &I0) {
  const Size nv = I0.size();
  const Size N  = Ts.size();

  I = I0;
  dI.resize(N);

  if (N == 0) return;

  const Size nq = dTs.front().nrows();

  ARTS_USER_ERROR_IF(not arr::same_size(Ts, Pi, dTs),
                     R"(
  Mismatched number of levels ({}) in Ts:s, Pi:s, and dTs:s:
  Ts size:  {}
  Pi size:  {}
  dTs size: {}
)",
                     N,
                     Ts.size(),
                     Pi.size(),
                     dTs.size());

  ARTS_USER_ERROR_IF(
      not all_same_shape(I0, Ts, Pi),
      R"(Must have same number of frequency elements ({}) in all Ts:s, and Pi:s)",
      nv);

  ARTS_USER_ERROR_IF(
      not all_same_shape({2, nq, nv}, dTs),
      R"(Must have same number of elements ({}, {}, {}) in all dTs:s)",
      2,
      nq,
      nv);

  for (auto &x : dI) {
    x.resize(nq, nv);
    x = 0.0;
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    I[iv] = Pi.back()[iv] * I0[iv];
  }

  if (nq == 0 or N == 1) return;

  // Add transmittance background
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    const auto src = I0[iv];

    muelmat P = 1.0;
    for (Size i = N - 2; i > 0; i--) {
      const muelmat R = Ts[i + 1][iv] * P;
      for (Size iq = 0; iq < nq; iq++) {
        dI[i][iq, iv] += Pi[i][iv] * dTs[i][0, iq, iv] * P * src;
        dI[i][iq, iv] += Pi[i - 1][iv] * dTs[i][1, iq, iv] * R * src;
      }
      P = R;
    }

    for (Size iq = 0; iq < nq; iq++) {
      dI[0][iq, iv]     += Pi[0][iv] * dTs[0][0, iq, iv] * P * src;
      dI[N - 1][iq, iv] += Pi[N - 2][iv] * dTs[N - 1][1, iq, iv] * src;
    }
  }
}

void nlte_step(stokvec_vector_view I,
               const Vector &f,
               const propmat_vector_const_view &K0,
               const propmat_vector_const_view &K1,
               const stokvec_vector_const_view &J0,
               const stokvec_vector_const_view &J1,
               const Numeric &T0,
               const Numeric &T1,
               const Numeric &r) {
  const Size NV = I.size();
  assert(arr::same_size(I, f, K0, K1, J0, J1));

  for (Size iv = 0; iv < NV; iv++) {
    const Numeric B0 = planck(f[iv], T0);
    const Numeric B1 = planck(f[iv], T1);
    const stokvec J =
        std::midpoint(B0, B1) + avg(inv(K0[iv]) * J0[iv], inv(K1[iv]) * J1[iv]);

    I[iv] = tran(K0[iv], K1[iv], r)() * (I[iv] - J) + J;
  }
}
}  // namespace rtepack