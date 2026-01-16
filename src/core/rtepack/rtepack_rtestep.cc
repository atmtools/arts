#include "rtepack_rtestep.h"

#include <array_algo.h>
#include <arts_omp.h>
#include <debug.h>
#include <physics_funcs.h>

#include "rtepack_mueller_matrix.h"
#include "rtepack_multitype.h"
#include "rtepack_stokes_vector.h"
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
        auto dI0       = dI[i][joker, iv];
        auto dI1       = dI[i + 1][joker, iv];
        const auto dT0 = dTs[i][0, joker, iv];
        const auto dJ0 = dJs[i][joker, iv];
        const auto dT1 = dTs[i + 1][1, joker, iv];
        const auto dJ1 = dJs[i + 1][joker, iv];
        const auto &P  = Pi[i][iv];

        for (Size iq = 0; iq < nq; iq++) {
          dI0[iq] += P * (dT0[iq] * Iv + avg(dJ0[iq], -(T * dJ0[iq])));
          dI1[iq] += P * (dT1[iq] * Iv + avg(dJ1[iq], -(T * dJ1[iq])));
        }
      }

      Iv = T * Iv + Jv;
    }
  }
}

void two_level_linear_evolution_step_by_step_full(
    stokvec_vector &I,
    std::vector<stokvec_matrix> &dI,
    const std::vector<muelmat_vector> &Ts,
    const std::vector<muelmat_vector> &Ls,
    const std::vector<muelmat_vector> &Pi,
    const std::vector<muelmat_tensor3> &dTs,
    const std::vector<muelmat_tensor3> &dLs,
    const std::vector<stokvec_vector> &Js,
    const std::vector<stokvec_matrix> &dJs,
    const stokvec_vector &I0) {
  const Size nv = I0.size();
  const Size N  = Ts.size();

  I = I0;
  dI.resize(N);

  if (N == 0) return;

  const Size nq = dJs.front().nrows();

  ARTS_USER_ERROR_IF(not arr::same_size(Ts, Ls, Pi, dTs, dLs, Js, dJs),
                     R"(
  Mismatched number of levels ({}) in Ts:s, Ls:s, Pi:s, dTs:s, dLs:s, Js:s, and dJs:s:
  Ts  size: {}
  Ls  size: {}
  Pi  size: {}
  dTs size: {}
  dLs size: {}
  Js  size: {}
  dJs size: {}
)",
                     N,
                     Ts.size(),
                     Ls.size(),
                     Pi.size(),
                     dTs.size(),
                     dLs.size(),
                     Js.size(),
                     dJs.size());

  ARTS_USER_ERROR_IF(
      not all_same_shape(I0, Ts, Ls, Pi, Js),
      R"(Must have same number of frequency elements ({}) in all Ts:s, Ls:s, Pi:s, and Js:s)",
      nv);

  ARTS_USER_ERROR_IF(
      not all_same_shape({nq, nv}, dJs),
      R"(Must have same number of elements ({}, {}) in all dJs:s)",
      nq,
      nv);

  ARTS_USER_ERROR_IF(
      not all_same_shape({2, nq, nv}, dTs, dLs),
      R"(Must have same number of elements ({}, {}, {}) in all dTs:s and dLs:s)",
      2,
      nq,
      nv);

  for (auto &x : dI) {
    x.resize(nq, nv);
    x = 0.0;
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    stokvec &Ii = I[iv];

    for (Size i = N - 2; i < N; i--) {
      const stokvec &J0 = Js[i + 1][iv];
      const stokvec &J1 = Js[i][iv];
      const muelmat &T  = Ts[i + 1][iv];
      const muelmat &L  = Ls[i + 1][iv];

      const stokvec ImJ0  = Ii - J0;
      const stokvec J0mJ1 = J0 - J1;

      if (nq) {
        auto dI0       = dI[i][joker, iv];
        auto dI1       = dI[i + 1][joker, iv];
        const auto dT0 = dTs[i][0, joker, iv];
        const auto dL0 = dLs[i][0, joker, iv];
        const auto dJ0 = dJs[i][joker, iv];
        const auto dT1 = dTs[i + 1][1, joker, iv];
        const auto dL1 = dLs[i + 1][1, joker, iv];
        const auto dJ1 = dJs[i + 1][joker, iv];
        const auto &P  = Pi[i][iv];

        for (Size iq = 0; iq < nq; iq++) {
          dI0[iq] +=
              P * (dJ1[iq] - L * dJ0[iq] + dT0[iq] * ImJ0 + dL0[iq] * J0mJ1);
          dI1[iq] += P * (dT1[iq] * ImJ0 + dL1[iq] * J0mJ1 + L * dJ1[iq] -
                          T * dJ0[iq]);
        }
      }

      Ii = T * ImJ0 + L * J0mJ1 + J1;
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

namespace {
void constant(stokvec_vector_view &Is,
              stokvec_tensor3_view &dIs,
              const muelmat_matrix_const_view &Ts,
              const muelmat_matrix_const_view &Pi,
              const muelmat_tensor3_const_view &dTs0,
              const muelmat_tensor3_const_view &dTs1,
              const stokvec_matrix_const_view &Js,
              const stokvec_tensor3_const_view &dJs) {
  const Size nv = dIs.npages();
  const Size np = dIs.nrows();
  const Size nq = dIs.ncols();

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    stokvec &Iv                          = Is[iv];
    stokvec_matrix_view dIv              = dIs[iv];
    const muelmat_vector_const_view Tv   = Ts[iv];
    const muelmat_matrix_const_view dT0v = dTs0[iv];
    const muelmat_matrix_const_view dT1v = dTs1[iv];
    const stokvec_vector_const_view Jv   = Js[iv];
    const stokvec_matrix_const_view dJv  = dJs[iv];

    for (Size i = np - 2; i < np; i--) {
      const muelmat &T = Tv[i + 1];
      const stokvec J  = avg(Jv[i], Jv[i + 1]);

      Iv -= J;

      if (nq) {
        auto &&dI0       = dIv[i];
        auto &&dI1       = dIv[i + 1];
        const auto &&dT0 = dT0v[i];
        const auto &&dT1 = dT1v[i + 1];
        const auto &&dJ0 = dJv[i];
        const auto &&dJ1 = dJv[i + 1];
        const auto &P    = Pi[iv, i];

        for (Size iq = 0; iq < nq; iq++) {
          dI0[iq] += P * (dT0[iq] * Iv + avg(dJ0[iq], -(T * dJ0[iq])));
          dI1[iq] += P * (dT1[iq] * Iv + avg(dJ1[iq], -(T * dJ1[iq])));
        }
      }

      Iv = T * Iv + J;
    }
  }
}

void linevo(stokvec_vector_view &Is,
            stokvec_tensor3_view &dIs,
            const muelmat_matrix_const_view &Ts,
            const muelmat_matrix_const_view &Ls,
            const muelmat_matrix_const_view &Pi,
            const muelmat_tensor3_const_view &dTs0,
            const muelmat_tensor3_const_view &dTs1,
            const muelmat_tensor3_const_view &dLs0,
            const muelmat_tensor3_const_view &dLs1,
            const stokvec_matrix_const_view &Js,
            const stokvec_tensor3_const_view &dJs) {
  const Size nv = dIs.npages();
  const Size np = dIs.nrows();
  const Size nq = dIs.ncols();

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    stokvec &Iv                          = Is[iv];
    stokvec_matrix_view dIv              = dIs[iv];
    const muelmat_vector_const_view Tv   = Ts[iv];
    const muelmat_vector_const_view Lv   = Ls[iv];
    const muelmat_matrix_const_view dT0v = dTs0[iv];
    const muelmat_matrix_const_view dT1v = dTs1[iv];
    const muelmat_matrix_const_view dL0v = dLs0[iv];
    const muelmat_matrix_const_view dL1v = dLs1[iv];
    const stokvec_vector_const_view Jv   = Js[iv];
    const stokvec_matrix_const_view dJv  = dJs[iv];

    for (Size i = np - 2; i < np; i--) {
      const muelmat &T    = Tv[i + 1];
      const muelmat &L    = Lv[i + 1];
      const stokvec &J0   = Jv[i + 1];
      const stokvec &J1   = Jv[i];
      const stokvec ImJ0  = Iv - J0;
      const stokvec J0mJ1 = J0 - J1;

      if (nq) {
        auto &&dI0       = dIv[i];
        auto &&dI1       = dIv[i + 1];
        const auto &&dT0 = dT0v[i];
        const auto &&dT1 = dT1v[i + 1];
        const auto &&dL0 = dL0v[i];
        const auto &&dL1 = dL1v[i + 1];
        const auto &&dJ0 = dJv[i];
        const auto &&dJ1 = dJv[i + 1];
        const auto &P    = Pi[iv, i];

        for (Size iq = 0; iq < nq; iq++) {
          dI0[iq] +=
              P * (dJ1[iq] - L * dJ0[iq] + dT0[iq] * ImJ0 + dL0[iq] * J0mJ1);
          dI1[iq] += P * (dT1[iq] * ImJ0 + dL1[iq] * J0mJ1 + L * dJ1[iq] -
                          T * dJ0[iq]);
        }
      }

      Iv = T * ImJ0 + L * J0mJ1 + J1;
    }
  }
}
}  // namespace

void rte_emission(stokvec_vector_view I,
                  stokvec_tensor3_view dI,
                  const TransmittanceMatrix &tramat,
                  const SourceVector &srcvec) {
  switch (tramat.option) {
    case TransmittanceOption::constant:
      constant(I,
               dI,
               tramat.T,
               tramat.P,
               tramat.dT[0],
               tramat.dT[1],
               srcvec.J,
               srcvec.dJ);
      break;
    case TransmittanceOption::linsrc:
    case TransmittanceOption::linprop:
      linevo(I,
             dI,
             tramat.T,
             tramat.L,
             tramat.P,
             tramat.dT[0],
             tramat.dT[1],
             tramat.dL[0],
             tramat.dL[1],
             srcvec.J,
             srcvec.dJ);
      break;
  }
}

namespace {
void constant(stokvec_matrix_view Is,
              const muelmat_matrix_const_view &Ts,
              const stokvec_matrix_const_view &Js) {
  const Size nv = Is.nrows();
  const Size np = Is.ncols();

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    for (Size i = np - 2; i < np; i--) {
      const muelmat &T = Ts[iv, i + 1];
      const stokvec J  = avg(Js[iv, i], Js[iv, i + 1]);
      Is[iv, i]        = T * (Is[iv, i + 1] - J) + J;
    }
  }
}

void linevo(stokvec_matrix_view Is,
            const muelmat_matrix_const_view &Ts,
            const muelmat_matrix_const_view &Ls,
            const stokvec_matrix_const_view &Js) {
  const Size nv = Is.nrows();
  const Size np = Is.ncols();

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    for (Size i = np - 2; i < np; i--) {
      const muelmat &T  = Ts[iv, i + 1];
      const muelmat &L  = Ls[iv, i + 1];
      const stokvec &J0 = Js[iv, i + 1];
      const stokvec &J1 = Js[iv, i];
      Is[iv, i]         = T * (Is[iv, i + 1] - J0) + L * (J0 - J1) + J1;
    }
  }
}
}  // namespace

void rte_emission_path(stokvec_matrix_view Is,
                       const TransmittanceMatrix &Ts,
                       const SourceVector &Js) {
  if (Is.size() == 0) return;

  switch (Ts.option) {
    case TransmittanceOption::constant: constant(Is, Ts.T, Js.J); break;
    case TransmittanceOption::linsrc:
    case TransmittanceOption::linprop:  linevo(Is, Ts.T, Ts.L, Js.J); break;
  }
}

namespace {
void all(stokvec_vector_view I,
         stokvec_tensor3_view dI,
         const muelmat_matrix_const_view &Ts,
         const muelmat_matrix_const_view &Pi,
         const muelmat_tensor4_const_view &dTs,
         const stokvec_vector_const_view &I0) {
  const Size nv = I.size();
  const Size N  = Ts.ncols();
  const Size nq = dTs.ncols();

  if (N == 0) return;

  dI = 0;

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) I[iv] = Pi[iv, N - 1] * I0[iv];

  if (nq == 0 or N == 1) return;

  // Add transmittance background
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    const auto src = I0[iv];

    muelmat P = 1.0;
    for (Size i = N - 2; i > 0; i--) {
      const muelmat R = Ts[i + 1][iv] * P;
      for (Size iq = 0; iq < nq; iq++) {
        dI[iv, i, iq]     += Pi[iv, i] * dTs[0, iv, i, iq] * P * src;
        dI[iv, i - 1, iq] += Pi[iv, i - 1] * dTs[1, iv, i, iq] * R * src;
      }
      P = R;
    }

    for (Size iq = 0; iq < nq; iq++) {
      dI[iv, 0, iq]     += Pi[iv, 0] * dTs[0, iv, 0, iq] * P * src;
      dI[iv, N - 1, iq] += Pi[iv, N - 2] * dTs[1, iv, N - 1, iq] * src;
    }
  }
}
}  // namespace

void rte_transmission(stokvec_vector_view I,
                      stokvec_tensor3_view dI,
                      const TransmittanceMatrix &tramat,
                      const stokvec_vector_const_view &I0) {
  all(I, dI, tramat.T, tramat.P, tramat.dT, I0);
}
}  // namespace rtepack