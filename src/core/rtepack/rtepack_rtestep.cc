#include "rtepack_rtestep.h"

#include <arts_omp.h>

namespace rtepack {
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
  const Index N = I.nelem();
  const Index M = dI1.nrows();

  ARTS_ASSERT(N == dI1.ncols() and N == dI2.ncols() and N == J1.nelem() and
              N == J2.nelem() and N == dJ1.ncols() and N == dJ2.ncols() and
              N == T.nelem() and N == PiT.nelem() and N == dT1.ncols() and
              N == dT2.ncols())
  ARTS_ASSERT(M == dI2.nrows() and M == dJ1.nrows() and M == dJ2.nrows() and
              M == dT1.nrows() and M == dT2.nrows())

  for (Index i = 0; i < N; i++) {
    const auto J = avg(J1[i], J2[i]);
    I[i]         = I[i] - J;

    for (Index j = 0; j < M; j++) {
      dI1(j, i) += PiT[i] * dT1(j, i) * I[i] + dJ1(j, i) - T[i] * dJ1(j, i);
      dI2(j, i) += PiT[i] * dT2(j, i) * I[i] + dJ2(j, i) - T[i] * dJ2(j, i);
    }

    I[i] = T[i] * I[i] + J;
  }
}

void two_level_linear_emission_step(stokvec_vector_view I,
                                    const stokvec_vector_const_view &J1,
                                    const stokvec_vector_const_view &J2,
                                    const muelmat_vector_const_view &T) {
  const Index N = I.nelem();
  ARTS_ASSERT(N == J1.nelem() and N == J2.nelem() and N == T.nelem())

  for (Index i = 0; i < N; i++) {
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
  const Index N = I.nelem();
  const Index M = dI1.nrows();
  ARTS_ASSERT(N == dI1.ncols() and N == dI2.ncols() and N == T.nelem() and
              N == PiT.nelem() and N == dT1.ncols() and N == dT2.ncols())
  ARTS_ASSERT(M == dI2.nrows() and M == dT1.nrows() and M == dT2.nrows())

  for (Index j = 0; j < M; j++) {
    for (Index i = 0; i < N; i++) {
      dI1(j, i) += PiT[i] * dT1(j, i) * I[i];
      dI2(j, i) += PiT[i] * dT2(j, i) * I[i];
    }
  }

  for (Index i = 0; i < N; i++) {
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
  const Index nv = I0.size();
  const Size N   = Ts.size();

  ARTS_USER_ERROR_IF(N != dTs.size(),
                     "Must have same number of levels (",
                     N,
                     ") in Ts and dTs");

  ARTS_USER_ERROR_IF(
      N != Js.size(), "Must have same number of levels (", N, ") in Ts and Js");

  ARTS_USER_ERROR_IF(N != dJs.size(),
                     "Must have same number of levels (",
                     N,
                     ") in Ts and dJs");

  I = I0;
  dI.resize(N);

  const Index nq = dJs[0].nrows();

  for (auto &x : dI) {
    x.resize(nq, nv);
    x = 0.0;
  }

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(Ts, Cmp::ne(nv), &muelmat_vector::size),
      "Must have same number of frequency elements (",
      nv,
      ") in all Ts:s");

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(Js, Cmp::ne(nv), &stokvec_vector::size),
      "Must have same number of frequency elements (",
      nv,
      ") in all Js:s");

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(dTs, Cmp::ne(nv), &muelmat_tensor3::ncols) or
          std::ranges::any_of(dTs, Cmp::ne(nq), &muelmat_tensor3::nrows) or
          std::ranges::any_of(dTs, Cmp::ne(2), &muelmat_tensor3::npages),
      "Must have same number of derivative elements (",
      2,
      ", ",
      nq,
      ", ",
      nv,
      ") in all dTs:s");

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(dJs, Cmp::ne(nv), &stokvec_matrix::ncols) or
          std::ranges::any_of(dJs, Cmp::ne(nq), &stokvec_matrix::nrows),
      "Must have same number of derivative elements (",
      nq,
      ", ",
      nv,
      ") in all dJs:s");

  if (N == 0) return;

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Index iv = 0; iv < nv; iv++) {
    stokvec &Iv = I[iv];
    for (Size i = N - 2; i < N; i--) {
      const stokvec Jv = avg(Js[i][iv], Js[i + 1][iv]);

      Iv -= Jv;

      for (Index iq = 0; iq < nq; iq++) {
        dI[i](iq, iv) +=
            Pi[i][iv] * (dTs[i](0, iq, iv) * Iv +
                         (dJs[i](iq, iv) - Ts[i + 1][iv] * dJs[i](iq, iv)));
        dI[i + 1](iq, iv) +=
            Pi[i][iv] *
            (dTs[i + 1](1, iq, iv) * Iv +
             (dJs[i + 1](iq, iv) - Ts[i + 1][iv] * dJs[i + 1](iq, iv)));
      }

      Iv = Ts[i + 1][iv] * Iv + Jv;
    }
  }
}

void two_level_linear_emission_cumulative_full(
    stokvec_vector &I,
    std::vector<stokvec_matrix> &dI,
    const std::vector<muelmat_vector> &Ts,
    const std::vector<muelmat_vector> &Pi,
    const std::vector<muelmat_tensor3> &dTs,
    const std::vector<stokvec_vector> &Js,
    const std::vector<stokvec_matrix> &dJs,
    const stokvec_vector &I0) {
  const Index nv = I0.size();
  const Size N   = Ts.size();

  ARTS_USER_ERROR_IF(N != dTs.size(),
                     "Must have same number of levels (",
                     N,
                     ") in Ts and dTs");

  ARTS_USER_ERROR_IF(
      N != Js.size(), "Must have same number of levels (", N, ") in Ts and Js");

  ARTS_USER_ERROR_IF(N != dJs.size(),
                     "Must have same number of levels (",
                     N,
                     ") in Ts and dJs");

  I.resize(nv);
  dI.resize(N);

  const Index nq = dJs[0].nrows();

  for (auto &x : dI) {
    x.resize(nq, nv);
    x = 0.0;
  }

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(Ts, Cmp::ne(nv), &muelmat_vector::size),
      "Must have same number of frequency elements (",
      nv,
      ") in all Ts:s");

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(Js, Cmp::ne(nv), &stokvec_vector::size),
      "Must have same number of frequency elements (",
      nv,
      ") in all Js:s");

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(dTs, Cmp::ne(nv), &muelmat_tensor3::ncols) or
          std::ranges::any_of(dTs, Cmp::ne(nq), &muelmat_tensor3::nrows) or
          std::ranges::any_of(dTs, Cmp::ne(2), &muelmat_tensor3::npages),
      "Must have same number of derivative elements (",
      2,
      ", ",
      nq,
      ", ",
      nv,
      ") in all dTs:s");

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(dJs, Cmp::ne(nv), &stokvec_matrix::ncols) or
          std::ranges::any_of(dJs, Cmp::ne(nq), &stokvec_matrix::nrows),
      "Must have same number of derivative elements (",
      nq,
      ", ",
      nv,
      ") in all dJs:s");

  if (N == 0) {
    I = I0;
    return;
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Index iv = 0; iv < nv; iv++) {
    I[iv] = Pi.back()[iv] * (I0[iv] - 0.5 * (Js[N - 2][iv] + Js[N - 1][iv])) +
            0.5 * (Js[0][iv] + Js[1][iv]);

    for (Size i = 1; i < N - 1; i++) {
      I[iv] += 0.5 * Pi[i][iv] * (Js[i + 1][iv] - Js[i - 1][iv]);
    }
  }

  if (nq == 0) return;

    // Add non-transmittance
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Index iv = 0; iv < nv; iv++) {
    for (Index iq = 0; iq < nq; iq++) {
      dI[0](iq, iv)     += 0.5 * dJs[0](iq, iv);
      dI[1](iq, iv)     += 0.5 * dJs[1](iq, iv);
      dI[N - 2](iq, iv) -= 0.5 * Pi.back()[iv] * dJs[N - 2](iq, iv);
      dI[N - 1](iq, iv) -= 0.5 * Pi.back()[iv] * dJs[N - 1](iq, iv);

      for (Size j = 1; j < N - 1; j++) {
        dI[j + 1](iq, iv) += 0.5 * Pi[j][iv] * dJs[j + 1](iq, iv);
        dI[j - 1](iq, iv) -= 0.5 * Pi[j][iv] * dJs[j - 1](iq, iv);
      }
    }
  }

  // Add transmittance background
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Index iv = 0; iv < nv; iv++) {
    const auto src = (I0[iv] - 0.5 * (Js[N - 2][iv] + Js[N - 1][iv]));

    muelmat P = 1.0;
    for (Size i = N - 2; i > 0; i--) {
      const muelmat R = Ts[i + 1][iv] * P;
      for (Index iq = 0; iq < nq; iq++) {
        dI[i](iq, iv) += Pi[i][iv] * dTs[i](0, iq, iv) * P * src;
        dI[i](iq, iv) += Pi[i - 1][iv] * dTs[i](1, iq, iv) * R * src;
      }
      P = R;
    }

    for (Index iq = 0; iq < nq; iq++) {
      dI[0](iq, iv)     += Pi[0][iv] * dTs[0](0, iq, iv) * P * src;
      dI[N - 1](iq, iv) += Pi[N - 2][iv] * dTs[N - 1](1, iq, iv) * src;
    }
  }

  // Add transmittance layers
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Index iv = 0; iv < nv; iv++) {
    for (Size j = 1; j < N - 1; j++) {
      const auto jsrc = 0.5 * (Js[j + 1][iv] - Js[j - 1][iv]);

      muelmat P = 1.0;
      for (Size i = j - 1; i > 0; i--) {
        const muelmat R = Ts[i + 1][iv] * P;
        for (Index iq = 0; iq < nq; iq++) {
          dI[i](iq, iv) += Pi[i][iv] * dTs[i](0, iq, iv) * P * jsrc;
          dI[i](iq, iv) += Pi[i - 1][iv] * dTs[i](1, iq, iv) * R * jsrc;
        }
        P = R;
      }

      for (Index iq = 0; iq < nq; iq++) {
        dI[j](iq, iv) += Pi[j - 1][iv] * dTs[j](1, iq, iv) * jsrc;
        dI[0](iq, iv) += (dTs[0](0, iq, iv) + dTs[0](1, iq, iv)) * P * jsrc;
      }
    }
  }
}

void two_level_linear_transmission_step(stokvec_vector &I,
                                        std::vector<stokvec_matrix> &dI,
                                        const std::vector<muelmat_vector> &Ts,
                                        const std::vector<muelmat_vector> &Pi,
                                        const std::vector<muelmat_tensor3> &dTs,
                                        const stokvec_vector &I0) {
  const Index nv = I0.size();
  const Size N   = Ts.size();

  ARTS_USER_ERROR_IF(N != dTs.size(),
                     "Must have same number of levels (",
                     N,
                     ") in Ts and dTs");

  I.resize(nv);
  dI.resize(N);

  const Index nq = dTs[0].nrows();

  for (auto &x : dI) {
    x.resize(nq, nv);
    x = 0.0;
  }

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(Ts, Cmp::ne(nv), &muelmat_vector::size),
      "Must have same number of frequency elements (",
      nv,
      ") in all Ts:s");

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(dTs, Cmp::ne(nv), &muelmat_tensor3::ncols) or
          std::ranges::any_of(dTs, Cmp::ne(nq), &muelmat_tensor3::nrows) or
          std::ranges::any_of(dTs, Cmp::ne(2), &muelmat_tensor3::npages),
      "Must have same number of derivative elements (",
      2,
      ", ",
      nq,
      ", ",
      nv,
      ") in all dTs:s");

  if (N == 0) {
    I = I0;
    return;
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Index iv = 0; iv < nv; iv++) {
    I[iv] = Pi.back()[iv] * I0[iv];
  }

  if (nq == 0) return;

    // Add transmittance background
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Index iv = 0; iv < nv; iv++) {
    const auto src = I0[iv];

    muelmat P = 1.0;
    for (Size i = N - 2; i > 0; i--) {
      const muelmat R = Ts[i + 1][iv] * P;
      for (Index iq = 0; iq < nq; iq++) {
        dI[i](iq, iv) += Pi[i][iv] * dTs[i](0, iq, iv) * P * src;
        dI[i](iq, iv) += Pi[i - 1][iv] * dTs[i](1, iq, iv) * R * src;
      }
      P = R;
    }

    for (Index iq = 0; iq < nq; iq++) {
      dI[0](iq, iv)     += Pi[0][iv] * dTs[0](0, iq, iv) * P * src;
      dI[N - 1](iq, iv) += Pi[N - 2][iv] * dTs[N - 1](1, iq, iv) * src;
    }
  }
}
}  // namespace rtepack