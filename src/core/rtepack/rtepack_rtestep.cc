#include "rtepack_rtestep.h"

#include <array_algo.h>
#include <arts_omp.h>
#include <physics_funcs.h>

#include "rtepack_transmission.h"

namespace rtepack {
stokvec linear_step(const muelmat &T, const stokvec &I, const stokvec &J) {
  return T * (I - J) + J;
}

stokvec two_level_linear_step(stokvec_vector_view &dI1,
                              stokvec_vector_view &dI2,
                              const muelmat &T,
                              const muelmat &PiT,
                              const stokvec &I,
                              const stokvec &J1,
                              const stokvec &J2,
                              const muelmat_vector_const_view &dT1,
                              const muelmat_vector_const_view &dT2,
                              const stokvec_vector_const_view &dJ1,
                              const stokvec_vector_const_view &dJ2) {
  assert(dI1.size() == dI2.size() and dT1.size() == dT2.size() and
         dJ1.size() == dJ2.size() and dI1.size() == dT1.size());

  const auto J = avg(J1, J2);

  for (Size k = 0; k < dT1.size(); k++) {
    dI1[k] += PiT * (dT1[k] * (I - J) + dJ1[k] - T * dJ1[k]);
    dI2[k] += PiT * (dT2[k] * (I - J) + dJ2[k] - T * dJ2[k]);
  }

  return linear_step(T, I, J);
}

void two_level_linear_emission_step(stokvec_vector_view I,
                                    stokvec_matrix_view dI1,
                                    stokvec_matrix_view dI2,
                                    const stokvec_vector_const_view &J1,
                                    const stokvec_vector_const_view &J2,
                                    const stokvec_matrix_const_view &dJ1,
                                    const stokvec_matrix_const_view &dJ2,
                                    const muelmat_vector_const_view &T,
                                    const muelmat_vector_const_view &PiT,
                                    const muelmat_matrix_const_view &dT1,
                                    const muelmat_matrix_const_view &dT2) {
  const Size N = I.size();
  const Size M = dI1.nrows();

  assert(N == static_cast<Size>(dI1.ncols()) and
         N == static_cast<Size>(dI2.ncols()) and N == J1.size() and
         N == J2.size() and N == static_cast<Size>(dJ1.ncols()) and
         N == static_cast<Size>(dJ2.ncols()) and N == T.size() and
         N == PiT.size() and N == static_cast<Size>(dT1.ncols()) and
         N == static_cast<Size>(dT2.ncols()));
  assert(M == static_cast<Size>(dI2.nrows()) and
         M == static_cast<Size>(dJ1.nrows()) and
         M == static_cast<Size>(dJ2.nrows()) and
         M == static_cast<Size>(dT1.nrows()) and
         M == static_cast<Size>(dT2.nrows()));

  for (Size i = 0; i < N; i++) {
    const auto J = avg(J1[i], J2[i]);
    I[i]         = I[i] - J;

    for (Size j = 0; j < M; j++) {
      dI1[j, i] += PiT[i] * dT1[j, i] * I[i] + dJ1[j, i] - T[i] * dJ1[j, i];
      dI2[j, i] += PiT[i] * dT2[j, i] * I[i] + dJ2[j, i] - T[i] * dJ2[j, i];
    }

    I[i] = T[i] * I[i] + J;
  }
}

void two_level_linear_emission_step(stokvec_vector_view I,
                                    const stokvec_vector_const_view &J1,
                                    const stokvec_vector_const_view &J2,
                                    const muelmat_vector_const_view &T) {
  const Size N = I.size();
  assert(N == J1.size() and N == J2.size() and N == T.size());

  for (Size i = 0; i < N; i++) {
    const auto J = avg(J1[i], J2[i]);
    I[i]         = T[i] * (I[i] - J) + J;
  }
}

void two_level_linear_transmission_step(stokvec_vector_view I,
                                        stokvec_matrix_view dI1,
                                        stokvec_matrix_view dI2,
                                        const muelmat_vector_const_view &T,
                                        const muelmat_vector_const_view &PiT,
                                        const muelmat_matrix_const_view &dT1,
                                        const muelmat_matrix_const_view &dT2) {
  const Size N = I.size();
  const Size M = dI1.nrows();
  assert(N == static_cast<Size>(dI1.ncols()) and
         N == static_cast<Size>(dI2.ncols()) and N == T.size() and
         N == PiT.size() and N == static_cast<Size>(dT1.ncols()) and
         N == static_cast<Size>(dT2.ncols()));
  assert(M == static_cast<Size>(dI2.nrows()) and
         M == static_cast<Size>(dT1.nrows()) and
         M == static_cast<Size>(dT2.nrows()));

  for (Size j = 0; j < M; j++) {
    for (Size i = 0; i < N; i++) {
      dI1[j, i] += PiT[i] * dT1[j, i] * I[i];
      dI2[j, i] += PiT[i] * dT2[j, i] * I[i];
    }
  }

  for (Size i = 0; i < N; i++) {
    I[i] = T[i] * I[i];
  }
}

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

  ARTS_USER_ERROR_IF(
      N != dTs.size(), "Must have same number of levels ({}) in Ts and dTs", N);

  ARTS_USER_ERROR_IF(
      N != Js.size(), "Must have same number of levels ({}) in Ts and Js", N);

  ARTS_USER_ERROR_IF(
      N != dJs.size(), "Must have same number of levels ({}) in Ts and dJs", N);

  I = I0;
  dI.resize(N);

  const Index nq = dJs[0].nrows();

  for (auto &x : dI) {
    x.resize(nq, nv);
    x = 0.0;
  }

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(elemwise_range(Ts),
                          Cmp::ne(nv),
                          [](const auto &x) { return x.size(); }),
      "Must have same number of frequency elements ({}) in all Ts:s",
      nv);

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(elemwise_range(Js),
                          Cmp::ne(nv),
                          [](const auto &x) { return x.size(); }),
      "Must have same number of frequency elements ({}) in all Js:s",
      nv);

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(
          dTs, Cmp::ne(static_cast<Index>(nv)), &muelmat_tensor3::ncols) or
          std::ranges::any_of(dTs, Cmp::ne(nq), &muelmat_tensor3::nrows) or
          std::ranges::any_of(dTs, Cmp::ne(2), &muelmat_tensor3::npages),
      "Must have same number of derivative elements (2, {}, {}) in all dTs:s",
      nq,
      nv);

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(
          dJs, Cmp::ne(static_cast<Index>(nv)), &stokvec_matrix::ncols) or
          std::ranges::any_of(dJs, Cmp::ne(nq), &stokvec_matrix::nrows),
      "Must have same number of derivative elements ({}, {}) in all dJs:s",
      nq,
      nv);

  if (N == 0) return;

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    stokvec &Iv = I[iv];
    for (Size i = N - 2; i < N; i--) {
      const stokvec Jv = avg(Js[i][iv], Js[i + 1][iv]);

      Iv -= Jv;

      for (Index iq = 0; iq < nq; iq++) {
        dI[i][iq, iv] +=
            Pi[i][iv] *
            (dTs[i][0, iq, iv] * Iv +
             0.5 * (dJs[i][iq, iv] - Ts[i + 1][iv] * dJs[i][iq, iv]));
        dI[i + 1][iq, iv] +=
            Pi[i][iv] *
            (dTs[i + 1][1, iq, iv] * Iv +
             0.5 * (dJs[i + 1][iq, iv] - Ts[i + 1][iv] * dJs[i + 1][iq, iv]));
      }

      Iv = Ts[i + 1][iv] * Iv + Jv;
    }
  }
}

void two_level_linear_emission_step_by_step_full(
    std::vector<stokvec_vector> &Is,
    const std::vector<muelmat_vector> &Ts,
    const std::vector<stokvec_vector> &Js) {
  const Size N = Ts.size();

  ARTS_USER_ERROR_IF(not arr::same_size(Is, Ts, Js),
                     R"(Not same sizes:

Is.size() = {},
Ts.size() = {},
Js.size() = {}
)",
                     Is.size(),
                     Ts.size(),
                     Js.size())

  if (N == 0) return;

  ARTS_USER_ERROR_IF(
      not arr::elemwise_same_size(Is, Ts, Js),
      "Not all elements have the same number of frequency elements")

  const Size nv = Is.front().size();

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

  ARTS_USER_ERROR_IF(
      N != dTs.size(), "Must have same number of levels ({}) in Ts and dTs", N);

  I.resize(nv);
  dI.resize(N);

  const Index nq = dTs[0].nrows();

  for (auto &x : dI) {
    x.resize(nq, nv);
    x = 0.0;
  }

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(elemwise_range(Ts),
                          Cmp::ne(nv),
                          [](const auto &x) { return x.size(); }),
      "Must have same number of frequency elements ({}) in all Ts:s",
      nv);

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(
          dTs, Cmp::ne(static_cast<Index>(nv)), &muelmat_tensor3::ncols) or
          std::ranges::any_of(dTs, Cmp::ne(nq), &muelmat_tensor3::nrows) or
          std::ranges::any_of(dTs, Cmp::ne(2), &muelmat_tensor3::npages),
      "Must have same number of derivative elements (2, {}, {}) in all dTs:s",
      nq,
      nv);

  if (N == 0) {
    I = I0;
    return;
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    I[iv] = Pi.back()[iv] * I0[iv];
  }

  if (nq == 0) return;

  // Add transmittance background
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    const auto src = I0[iv];

    muelmat P = 1.0;
    for (Size i = N - 2; i > 0; i--) {
      const muelmat R = Ts[i + 1][iv] * P;
      for (Index iq = 0; iq < nq; iq++) {
        dI[i][iq, iv] += Pi[i][iv] * dTs[i][0, iq, iv] * P * src;
        dI[i][iq, iv] += Pi[i - 1][iv] * dTs[i][1, iq, iv] * R * src;
      }
      P = R;
    }

    for (Index iq = 0; iq < nq; iq++) {
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