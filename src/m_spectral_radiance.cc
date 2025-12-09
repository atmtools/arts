#include <array_algo.h>
#include <arts_omp.h>
#include <atm_field.h>
#include <auto_wsm.h>
#include <configtypes.h>
#include <enumsAtmKey.h>
#include <jacobian.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <rtepack_transmission.h>
#include <surf_field.h>
#include <time_report.h>
#include <workspace.h>

void spectral_tramat_pathFromPath(
    ArrayOfMuelmatVector& spectral_tramat_path,
    ArrayOfMuelmatTensor3& spectral_tramat_jac_path,
    const ArrayOfPropmatVector& spectral_propmat_path,
    const ArrayOfPropmatMatrix& spectral_propmat_jac_path,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfAtmPoint& atm_path,
    const SurfaceField& surf_field,
    const JacobianTargets& jac_targets,
    const Index& hse_derivative) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  ARTS_USER_ERROR_IF(ray_path.size() == 0, "Empty path.");

  ARTS_USER_ERROR_IF(
      not arr::same_size(
          ray_path, spectral_propmat_path, spectral_propmat_jac_path, atm_path),
      R"(Not same sizes:

ray_path.size()                  = {},
spectral_propmat_path.size()     = {},
spectral_propmat_jac_path.size() = {},
atm_path.size()                  = {}
)",
      ray_path.size(),
      spectral_propmat_path.size(),
      spectral_propmat_jac_path.size(),
      atm_path.size());

  // HSE variables
  const Index temperature_derivative_position =
      jac_targets.target_position(AtmKey::t);

  const Size N = ray_path.size();

  spectral_tramat_path.resize(N);
  spectral_tramat_jac_path.resize(N);

  if (N == 0) return;

  const Index nq = jac_targets.target_count();

  const Vector ray_path_distance = distance(ray_path, surf_field.ellipsoid);
  Tensor3 ray_path_distance_jacobian(2, N - 1, nq, 0.0);

  // FIXME: UNKNOWN ORDER OF "ip" AND "ip - 1" FOR TEMPERATURE
  if (hse_derivative and temperature_derivative_position >= 0) {
    for (Size ip = 0; ip < N - 1; ip++) {
      ray_path_distance_jacobian[0, ip, temperature_derivative_position] =
          ray_path_distance[ip] / (2.0 * atm_path[ip].temperature);
      ray_path_distance_jacobian[1, ip, temperature_derivative_position] =
          ray_path_distance[ip] / (2.0 * atm_path[ip + 1].temperature);
    }
  }

  rtepack::two_level_exp(spectral_tramat_path,
                         spectral_tramat_jac_path,
                         spectral_propmat_path,
                         spectral_propmat_jac_path,
                         ray_path_distance,
                         ray_path_distance_jacobian);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radStepByStepEmission(
    StokvecVector& spectral_rad,
    ArrayOfStokvecMatrix& spectral_rad_jac_path,
    const ArrayOfMuelmatVector& spectral_tramat_path,
    const ArrayOfMuelmatVector& spectral_tramat_cumulative_path,
    const ArrayOfMuelmatTensor3& spectral_tramat_jac_path,
    const ArrayOfStokvecVector& spectral_rad_srcvec_path,
    const ArrayOfStokvecMatrix& spectral_rad_srcvec_jac_path,
    const StokvecVector& spectral_rad_bkg) try {
  ARTS_TIME_REPORT

  rtepack::two_level_linear_emission_step_by_step_full(
      spectral_rad,
      spectral_rad_jac_path,
      spectral_tramat_path,
      spectral_tramat_cumulative_path,
      spectral_tramat_jac_path,
      spectral_rad_srcvec_path,
      spectral_rad_srcvec_jac_path,
      spectral_rad_bkg);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radCumulativeTransmission(
    StokvecVector& spectral_rad,
    ArrayOfStokvecMatrix& spectral_rad_jac_path,
    const ArrayOfMuelmatVector& spectral_tramat_path,
    const ArrayOfMuelmatVector& spectral_tramat_cumulative_path,
    const ArrayOfMuelmatTensor3& spectral_tramat_jac_path,
    const StokvecVector& spectral_rad_bkg) try {
  ARTS_TIME_REPORT

  rtepack::two_level_linear_transmission_step(spectral_rad,
                                              spectral_rad_jac_path,
                                              spectral_tramat_path,
                                              spectral_tramat_cumulative_path,
                                              spectral_tramat_jac_path,
                                              spectral_rad_bkg);
}
ARTS_METHOD_ERROR_CATCH

namespace {
void single_tramat_pathFromPath(MuelmatVector& single_tramat_path,
                                MuelmatTensor3& single_tramat_jac_path,
                                const PropmatVector& single_propmat_path,
                                const PropmatMatrix& single_propmat_jac_path,
                                const ArrayOfPropagationPathPoint& ray_path,
                                const ArrayOfAtmPoint& atm_path,
                                const SurfaceField& surf_field,
                                const JacobianTargets& jac_targets,
                                const Index& hse_derivative) try {
  ARTS_TIME_REPORT

  const Size N  = ray_path.size();
  const Size nq = jac_targets.target_count();

  ARTS_USER_ERROR_IF(
      single_propmat_path.size() != N or
          single_propmat_jac_path.nrows() != static_cast<Index>(N) or
          single_propmat_jac_path.ncols() != static_cast<Index>(nq) or
          atm_path.size() != N,
      R"(Not same sizes:

ray_path:                (N)     = [{}]
jac_targets:             (nq)    = [{}]
atm_path:                (N)     = [{}]
single_propmat_path:     (N)     = {:B,}
single_propmat_jac_path: (N, nq) = {:B,}
)",
      N,
      nq,
      atm_path.size(),
      single_propmat_path.shape(),
      single_propmat_jac_path.shape());

  single_tramat_path.resize(N);
  single_tramat_jac_path.resize(N, 2, nq);
  single_tramat_path     = 1.0;
  single_tramat_jac_path = 0.0;

  const Vector ray_path_distance = distance(ray_path, surf_field.ellipsoid);
  Tensor3 ray_path_distance_jacobian(2, N - 1, nq, 0.0);

  const Index temperature_derivative_position =
      jac_targets.target_position(AtmKey::t);
  if (hse_derivative and temperature_derivative_position >= 0) {
    for (Size ip = 0; ip < N - 1; ip++) {
      ray_path_distance_jacobian[0, ip, temperature_derivative_position] =
          ray_path_distance[ip] / (2.0 * atm_path[ip].temperature);
      ray_path_distance_jacobian[1, ip, temperature_derivative_position] =
          ray_path_distance[ip] / (2.0 * atm_path[ip + 1].temperature);
    }
  }

  for (Size i = 0; i < N - 1; i++) {
    auto& T        = single_tramat_path[i];
    auto dT1       = single_tramat_jac_path[i, 0];
    auto dT2       = single_tramat_jac_path[i + 1, 1];
    const auto& K1 = single_propmat_path[i];
    const auto& K2 = single_propmat_path[i + 1];
    const auto dK1 = single_propmat_jac_path[i];
    const auto dK2 = single_propmat_jac_path[i + 1];
    const auto r   = ray_path_distance[i];
    const auto dr1 = ray_path_distance_jacobian[0, i];
    const auto dr2 = ray_path_distance_jacobian[1, i];

    const rtepack::tran t(K1, K2, r);
    T = t();

    for (Size j = 0; j < nq; j++) {
      dT1[j] += t.deriv(T, K1, K2, dK1[j], r, dr1[j]);
      dT2[j] += t.deriv(T, K1, K2, dK2[j], r, dr2[j]);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void single_tramat_cumulative_pathFromPath(
    MuelmatVector& single_tramat_cumulative_path,
    const MuelmatVector& single_tramat_path) try {
  ARTS_TIME_REPORT

  const Size N = single_tramat_path.size();

  single_tramat_cumulative_path.resize(N);

  if (N == 0) return;

  auto& Pi = single_tramat_cumulative_path;
  auto& T  = single_tramat_path;

  Pi.front() = T.front();
  for (Size i = 1; i < N; i++) Pi[i] = Pi[i - 1] * T[i];
}
ARTS_METHOD_ERROR_CATCH

void single_rad_srcvec_jacFromPropmat(
    StokvecVector& single_rad_srcvec_jac,
    StokvecMatrix& single_rad_srcvec_jac_path,
    const PropmatVector& single_propmat_path,
    const StokvecVector& single_nlte_srcvec_path,
    const PropmatMatrix& single_propmat_jac_path,
    const StokvecMatrix& single_nlte_srcvec_jac_path,
    const Vector& single_freq_path,
    const ArrayOfAtmPoint& atm_path,
    const JacobianTargets& jac_targets) try {
  ARTS_TIME_REPORT

  const Index np = atm_path.size();
  const Index nq = jac_targets.target_count();

  ARTS_USER_ERROR_IF(np != single_freq_path.ncols() or
                         np != single_propmat_path.ncols() or
                         np != single_nlte_srcvec_path.ncols() or
                         np != single_propmat_jac_path.nrows() or
                         nq != single_propmat_jac_path.ncols() or
                         np != single_nlte_srcvec_jac_path.nrows() or
                         nq != single_nlte_srcvec_jac_path.ncols(),
                     R"(Not same sizes:

jac_targets                  (nq)     = [{}]
single_freq_path             (np)     = {:B,}
atm_path                     (np)     = [{}]
single_propmat_path:         (np)     = {:B,}
single_nlte_srcvec_path:     (np)     = {:B,}
single_propmat_jac_path:     (np, nq) = {:B,}
single_nlte_srcvec_jac_path: (np, nq) = {:B,}
)",
                     nq,
                     single_freq_path.shape(),
                     np,
                     single_propmat_path.shape(),
                     single_nlte_srcvec_path.shape(),
                     single_propmat_jac_path.shape(),
                     single_nlte_srcvec_jac_path.shape())

  single_rad_srcvec_jac.resize(np);
  single_rad_srcvec_jac_path.resize(np, nq);
  single_rad_srcvec_jac      = 0.0;
  single_rad_srcvec_jac_path = 0.0;

  if (np == 0) return;

  const Index it = jac_targets.target_position(AtmKey::t);

  for (Index i = 0; i < np; i++) {
    auto& J       = single_rad_srcvec_jac[i];
    auto dJ       = single_rad_srcvec_jac_path[i];
    const auto& K = single_propmat_path[i];
    const auto& S = single_nlte_srcvec_path[i];
    const auto& f = single_freq_path[i];
    const auto t  = atm_path[i].temperature;
    const auto dK = single_propmat_jac_path[i];
    const auto dS = single_nlte_srcvec_jac_path[i];

    if (K.is_rotational()) {
      J  = 0.0;
      dJ = 0.0;
    } else {
      const auto b    = Stokvec{planck(f, t), 0, 0, 0};
      const auto invK = inv(K);
      J               = b + invK * S;

      for (Index j = 0; j < nq; j++) {
        dJ[j]  = {it == j ? dplanck_dt(f, t) : 0, 0, 0, 0};
        dJ[j] -= invK * (dK[j] * invK * S - dS[j]);
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void single_radStepByStepEmission(
    Stokvec& single_rad,
    StokvecMatrix& single_rad_jac_path,
    const MuelmatVector& single_tramat_path,
    const MuelmatVector& single_tramat_cumulative_path,
    const MuelmatTensor3& single_tramat_jac_path,
    const StokvecVector& single_rad_srcvec_jac,
    const StokvecMatrix& single_rad_srcvec_jac_path) try {
  ARTS_TIME_REPORT

  const Size N  = single_tramat_path.size();
  const Size nq = single_rad_srcvec_jac_path.ncols();

  ARTS_USER_ERROR_IF(
      single_tramat_path.size() != N or
          single_tramat_cumulative_path.size() != N or
          single_rad_srcvec_jac.size() != N or
          single_rad_srcvec_jac_path.nrows() != static_cast<Index>(N) or
          single_rad_srcvec_jac_path.ncols() != static_cast<Index>(nq) or
          single_tramat_jac_path.npages() != static_cast<Index>(N) or
          single_tramat_jac_path.nrows() != 2,
      R"(Not same sizes:

single_rad_srcvec_jac         (N)        = {:B,}
single_tramat_path            (N)        = {:B,}
single_rad_srcvec_jac_path    (N, nq)    = {:B,}
single_tramat_jac_path        (N, 2, nq) = {:B,}
single_tramat_cumulative_path (N)        = {:B,}
)",
      single_rad_srcvec_jac.shape(),
      single_tramat_path.shape(),
      single_rad_srcvec_jac_path.shape(),
      single_tramat_jac_path.shape(),
      single_tramat_cumulative_path.shape());

  single_rad_jac_path.resize(N, nq);
  single_rad_jac_path = 0.0;

  if (N == 0) return;

  auto& I         = single_rad;
  auto& dI        = single_rad_jac_path;
  const auto& Js  = single_rad_srcvec_jac;
  const auto& dJs = single_rad_srcvec_jac_path;
  const auto& Pi  = single_tramat_cumulative_path;
  const auto& Ts  = single_tramat_path;
  const auto& dTs = single_tramat_jac_path;

  for (Size i = N - 2; i < N; i--) {
    const Stokvec J = avg(Js[i], Js[i + 1]);

    I -= J;

    for (Size iq = 0; iq < nq; iq++) {
      dI[i, iq] += Pi[i] * (dTs[i, 0, iq] * I +
                            0.5 * (dJs[i, iq] - Ts[i + 1] * dJs[i, iq]));
      dI[i + 1, iq] +=
          Pi[i] * (dTs[i + 1, 1, iq] * I +
                   0.5 * (dJs[i + 1, iq] - Ts[i + 1] * dJs[i + 1, iq]));
    }

    I = Ts[i + 1] * I + J;
  }
}
ARTS_METHOD_ERROR_CATCH

void single_rad_jacFromBackground(StokvecVector& single_rad_jac,
                                  const Muelmat& background_transmittance) try {
  ARTS_TIME_REPORT
  for (auto& dI : single_rad_jac) dI = background_transmittance * dI;
}
ARTS_METHOD_ERROR_CATCH

void single_rad_jacAddPathPropagation(
    StokvecVector& single_rad_jac,
    const StokvecMatrix& single_rad_jac_path,
    const JacobianTargets& jac_targets,
    const AtmField& atm_field,
    const ArrayOfPropagationPathPoint& ray_path) try {
  ARTS_TIME_REPORT

  const Size np  = ray_path.size();
  const Index nt = jac_targets.target_count();
  const Index nx = jac_targets.x_size();

  ARTS_USER_ERROR_IF(
      single_rad_jac.ncols() != nx or
          single_rad_jac_path.nrows() != static_cast<Index>(np) or
          single_rad_jac_path.ncols() != nt,
      R"(Not same sizes:

ray_path:            (np)     = [{}]
jac_targets:         (nq)     = [{}]
single_rad_jac:      (nx)     = {:B,}
single_rad_jac_path: (np, nq) = {:B,}
)",
      np,
      jac_targets.target_count(),
      single_rad_jac.shape(),
      single_rad_jac_path.shape());

  for (auto& atm_block : jac_targets.atm) {
    const auto& data = atm_field[atm_block.type];
    for (Size ip = 0; ip < np; ip++) {
      const auto weights = data.flat_weight(ray_path[ip].pos);
      const auto local   = single_rad_jac_path[ip];

      for (auto& w : weights) {
        if (w.second != 0.0) {
          const auto i       = w.first + atm_block.x_start;
          single_rad_jac[i] += w.second * local[atm_block.target_pos];
        }
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void single_radFromPropagation(Stokvec& single_rad,
                               StokvecVector& single_rad_jac,
                               StokvecVector& single_rad_srcvec_jac,
                               StokvecMatrix& single_rad_srcvec_jac_path,
                               StokvecMatrix& single_rad_jac_path,
                               MuelmatVector& single_tramat_path,
                               MuelmatVector& single_tramat_cumulative_path,
                               MuelmatTensor3& single_tramat_jac_path,
                               const JacobianTargets& jac_targets,
                               const ArrayOfPropagationPathPoint& ray_path,
                               const Vector& single_freq_path,
                               const ArrayOfAtmPoint& atm_path,
                               const PropmatVector& single_propmat_path,
                               const PropmatMatrix& single_propmat_jac_path,
                               const StokvecVector& single_nlte_srcvec_path,
                               const StokvecMatrix& single_nlte_srcvec_jac_path,
                               const SurfaceField& surf_field,
                               const AtmField& atm_field,
                               const Index& hse_derivative) try {
  ARTS_TIME_REPORT

  single_tramat_pathFromPath(single_tramat_path,
                             single_tramat_jac_path,
                             single_propmat_path,
                             single_propmat_jac_path,
                             ray_path,
                             atm_path,
                             surf_field,
                             jac_targets,
                             hse_derivative);
  single_tramat_cumulative_pathFromPath(single_tramat_cumulative_path,
                                        single_tramat_path);
  const auto& transmission_single_matrix_background =
      single_tramat_cumulative_path.back();
  single_rad_srcvec_jacFromPropmat(single_rad_srcvec_jac,
                                   single_rad_srcvec_jac_path,
                                   single_propmat_path,
                                   single_nlte_srcvec_path,
                                   single_propmat_jac_path,
                                   single_nlte_srcvec_jac_path,
                                   single_freq_path,
                                   atm_path,
                                   jac_targets);
  single_radStepByStepEmission(single_rad,
                               single_rad_jac_path,
                               single_tramat_path,
                               single_tramat_cumulative_path,
                               single_tramat_jac_path,
                               single_rad_srcvec_jac,
                               single_rad_srcvec_jac_path);
  single_rad_jacFromBackground(single_rad_jac,
                               transmission_single_matrix_background);
  single_rad_jacAddPathPropagation(
      single_rad_jac, single_rad_jac_path, jac_targets, atm_field, ray_path);
}
ARTS_METHOD_ERROR_CATCH
}  // namespace

void spectral_radSetToBackground(
    StokvecVector& spectral_rad,
    StokvecMatrix& spectral_rad_jac,
    const StokvecVector& spectral_rad_bkg,
    const StokvecMatrix& spectral_rad_bkg_jac) try {
  ARTS_TIME_REPORT

  spectral_rad     = spectral_rad_bkg;
  spectral_rad_jac = spectral_rad_bkg_jac;
}
ARTS_METHOD_ERROR_CATCH

void spectral_radSinglePathEmissionFrequencyLoop(
    StokvecVector& spectral_rad,
    StokvecMatrix& spectral_rad_jac,
    const JacobianTargets& jac_targets,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfAscendingGrid& freq_grid_path,
    const ArrayOfAtmPoint& atm_path,
    const ArrayOfPropmatVector& spectral_propmat_path,
    const ArrayOfStokvecVector& spectral_nlte_srcvec_path,
    const ArrayOfPropmatMatrix& spectral_propmat_jac_path,
    const ArrayOfStokvecMatrix& spectral_nlte_srcvec_jac_path,
    const SurfaceField& surf_field,
    const AtmField& atm_field,
    const Index& hse_derivative) try {
  ARTS_TIME_REPORT

  const Size nf = spectral_rad.size();
  const Size N  = ray_path.size();
  const Size nq = jac_targets.target_count();
  const Size nx = jac_targets.x_size();

  if (N == 0) return;

  ARTS_USER_ERROR_IF(not matpack::same_shape<2>({nx, nf}, spectral_rad_jac),
                     R"(spectral_rad_jac has wrong shape:
Expected: ({}, {})
Actual:   ({}, {})
)",
                     nx,
                     nf,
                     spectral_rad_jac.nrows(),
                     spectral_rad_jac.ncols());

  ARTS_USER_ERROR_IF(not arr::same_size(ray_path,
                                        freq_grid_path,
                                        atm_path,
                                        spectral_propmat_path,
                                        spectral_nlte_srcvec_path,
                                        spectral_propmat_jac_path,
                                        spectral_nlte_srcvec_jac_path),
                     R"(Not same sizes:

ray_path.size()                      = {},
freq_grid_path.size()                = {},
atm_path.size()                      = {},
spectral_propmat_path.size()         = {},
spectral_nlte_srcvec_path.size()     = {},
spectral_propmat_jac_path.size()     = {},
spectral_nlte_srcvec_jac_path.size() = {}
)",
                     ray_path.size(),
                     freq_grid_path.size(),
                     atm_path.size(),
                     spectral_propmat_path.size(),
                     spectral_nlte_srcvec_path.size(),
                     spectral_propmat_jac_path.size(),
                     spectral_nlte_srcvec_jac_path.size());

  ARTS_USER_ERROR_IF(
      not arr::elemwise_same_size(
          freq_grid_path, spectral_propmat_path, spectral_nlte_srcvec_path),
      R"(Not same sizes elemwise:
freq_grid_path.size()             = {}
spectral_propmat_path.shape()     = {:B,}
spectral_nlte_srcvec_path.shape() = {:B,}
)",
      freq_grid_path.size(),
      spectral_propmat_path[0].shape(),
      spectral_nlte_srcvec_path[0].shape());

  ARTS_USER_ERROR_IF(not arr::elemwise_same_size(spectral_propmat_jac_path,
                                                 spectral_nlte_srcvec_jac_path),
                     R"(Not same sizes elemwise:
spectral_propmat_jac_path.shape()     = {:B,}
spectral_nlte_srcvec_jac_path.shape() = {:B,}
)",
                     spectral_propmat_jac_path[0].shape(),
                     spectral_nlte_srcvec_jac_path[0].shape());

  Vector single_freq_path(N);
  StokvecVector single_rad_jac(nx);
  StokvecVector single_rad_srcvec_jac(N);
  StokvecMatrix single_rad_srcvec_jac_path(N, nq);
  StokvecMatrix single_rad_jac_path(N, nq);
  MuelmatVector single_tramat_path(N);
  MuelmatVector single_tramat_cumulative_path(N);
  MuelmatTensor3 single_tramat_jac_path(N, 2, nq);
  PropmatVector single_propmat_path(N);
  PropmatMatrix single_propmat_jac_path(N, nq);
  StokvecVector single_nlte_srcvec_path(N);
  StokvecMatrix single_nlte_srcvec_jac_path(N, nq);

  std::string error{};

#pragma omp parallel for firstprivate(     \
        single_freq_path,                  \
            single_rad_jac,                \
            single_rad_srcvec_jac,         \
            single_rad_srcvec_jac_path,    \
            single_rad_jac_path,           \
            single_tramat_path,            \
            single_tramat_cumulative_path, \
            single_tramat_jac_path,        \
            single_propmat_path,           \
            single_propmat_jac_path,       \
            single_nlte_srcvec_path,       \
            single_nlte_srcvec_jac_path) if (arts_omp_parallel())
  for (Size f = 0; f < nf; f++) {
    try {
      Stokvec& single_rad = spectral_rad[f];
      for (Size i = 0; i < N; i++) {
        single_freq_path[i]        = freq_grid_path[i][f];
        single_propmat_path[i]     = spectral_propmat_path[i][f];
        single_nlte_srcvec_path[i] = spectral_nlte_srcvec_path[i][f];

        single_propmat_jac_path[i] = spectral_propmat_jac_path[i][joker, f];
        single_nlte_srcvec_jac_path[i] =
            spectral_nlte_srcvec_jac_path[i][joker, f];
      }

      single_rad_jac = spectral_rad_jac[joker, f];

      single_radFromPropagation(single_rad,
                                single_rad_jac,
                                single_rad_srcvec_jac,
                                single_rad_srcvec_jac_path,
                                single_rad_jac_path,
                                single_tramat_path,
                                single_tramat_cumulative_path,
                                single_tramat_jac_path,
                                jac_targets,
                                ray_path,
                                single_freq_path,
                                atm_path,
                                single_propmat_path,
                                single_propmat_jac_path,
                                single_nlte_srcvec_path,
                                single_nlte_srcvec_jac_path,
                                surf_field,
                                atm_field,
                                hse_derivative);

      spectral_rad_jac[joker, f] = single_rad_jac;
    } catch (std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(
      not error.empty(),
      "Error in spectral_radSinglePathEmissionFrequencyLoop: {:s}",
      error);
}
ARTS_METHOD_ERROR_CATCH

void single_radClearskyEmissionPropagation(
    const Workspace& ws,
    Stokvec& single_rad,
    StokvecVector& single_rad_jac,
    ArrayOfPropagationPathPoint& ray_path,
    const AtmField& atm_field,
    const Numeric& frequency_,
    const JacobianTargets& jac_targets,
    const Agenda& single_rad_space_agenda,
    const Agenda& single_rad_surface_agenda,
    const Agenda& single_propmat_agenda,
    const Agenda& ray_point_back_propagation_agenda,
    const SubsurfaceField& subsurf_field,
    const SurfaceField& surf_field,
    const Vector3& obs_pos,
    const Vector2& obs_los,
    const Numeric& max_stepsize,
    const Propmat& polarization,
    const Numeric& max_tau,
    const Numeric& cutoff_tau,
    const Index& hse_derivative,
    const Index& N) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      N <= 0 or max_stepsize <= 0.0 or max_tau <= 0.0 or cutoff_tau <= 0.0,
      "Must be strictly positive: N = {}, max_stepsize = {}, max_tau = {}, cutoff_tau = {}",
      N,
      max_stepsize,
      max_tau,
      cutoff_tau);

  const Size nt = jac_targets.x_size();
  single_rad_jac.resize(nt);
  single_rad_jac = 0.0;

  Numeric single_dispersion;
  Vector single_dispersion_jac;

  Vector single_freq_path;
  ArrayOfAtmPoint atm_path;
  PropmatVector single_propmat_path;
  ArrayOfPropmatVector single_propmat_jac_path;
  StokvecVector single_nlte_srcvec_path;
  ArrayOfStokvecVector single_nlte_srcvec_jac_path;

  single_freq_path.reserve(N);
  atm_path.reserve(N);
  single_propmat_path.reserve(N);
  single_propmat_jac_path.reserve(N);
  single_nlte_srcvec_path.reserve(N);
  single_nlte_srcvec_jac_path.reserve(N);

  ray_path.clear();
  ray_path.reserve(N);
  ray_path.resize(1);
  ray_path[0] =
      path::init_with_lostype(obs_pos, obs_los, atm_field, surf_field, true);

  if (ray_path.back().los_type == PathPositionType::space) {
    single_rad_space_agendaExecute(ws,
                                   single_rad,
                                   single_rad_jac,
                                   frequency_,
                                   jac_targets,
                                   ray_path.back(),
                                   single_rad_space_agenda);
    return;
  }

  if (ray_path.back().los_type == PathPositionType::surface) {
    single_rad_surface_agendaExecute(ws,
                                     single_rad,
                                     single_rad_jac,
                                     frequency_,
                                     jac_targets,
                                     ray_path.back(),
                                     surf_field,
                                     subsurf_field,
                                     single_rad_surface_agenda);
    return;
  }

  Numeric sumAr = 0.0;
  while (true) {
    const auto& ray_point = ray_path.back();
    atm_path.emplace_back(atm_field.at(ray_point.pos));

    Numeric frequency = frequency_;
    Vector3 freq_wind_shift_jac;
    freqWindShift(frequency, freq_wind_shift_jac, atm_path.back(), ray_point);
    single_freq_path.push_back(frequency);

    single_propmat_agendaExecute(ws,
                                 single_propmat_path.emplace_back(),
                                 single_nlte_srcvec_path.emplace_back(),
                                 single_dispersion,
                                 single_propmat_jac_path.emplace_back(),
                                 single_nlte_srcvec_jac_path.emplace_back(),
                                 single_dispersion_jac,
                                 frequency,
                                 freq_wind_shift_jac,
                                 jac_targets,
                                 "AIR"_spec,
                                 ray_point,
                                 atm_path.back(),
                                 single_propmat_agenda);
    single_dispersion += dot(polarization, single_propmat_path.back());

    const Numeric A = single_propmat_path.back().A();
    if (A == 0.0) {
      PropagationPathPoint tmp_point;
      ray_point_back_propagation_agendaExecute(
          ws,
          tmp_point,
          ray_path,
          single_dispersion,
          single_propmat_path.back(),
          max_stepsize,
          ray_point_back_propagation_agenda);
      ray_path.push_back(tmp_point);
    } else {
      PropagationPathPoint tmp_point;
      ray_point_back_propagation_agendaExecute(
          ws,
          tmp_point,
          ray_path,
          single_dispersion,
          single_propmat_path.back(),
          std::min(max_tau / A, max_stepsize),
          ray_point_back_propagation_agenda);
      ray_path.push_back(tmp_point);
    }

    if (ray_path.back().los_type == PathPositionType::space) {
      single_rad_space_agendaExecute(ws,
                                     single_rad,
                                     single_rad_jac,
                                     frequency,
                                     jac_targets,
                                     ray_path.back(),
                                     single_rad_space_agenda);
      break;
    }

    if (ray_path.back().los_type == PathPositionType::surface) {
      single_rad_surface_agendaExecute(ws,
                                       single_rad,
                                       single_rad_jac,
                                       frequency,
                                       jac_targets,
                                       ray_path.back(),
                                       surf_field,
                                       subsurf_field,
                                       single_rad_surface_agenda);
      break;
    }

    sumAr += A * path::distance(ray_path.back().pos,
                                ray_path[ray_path.size() - 2].pos,
                                surf_field.ellipsoid);

    if (sumAr >= cutoff_tau) {
      single_rad = planck(frequency_, atm_path.back().temperature);
      break;
    }
  }

  ray_path.pop_back();
  StokvecVector single_rad_srcvec_jac;
  StokvecMatrix single_rad_srcvec_jac_path;
  StokvecMatrix single_rad_jac_path;
  MuelmatVector single_tramat_path;
  MuelmatVector single_tramat_cumulative_path;
  MuelmatTensor3 single_tramat_jac_path;
  single_radFromPropagation(single_rad,
                            single_rad_jac,
                            single_rad_srcvec_jac,
                            single_rad_srcvec_jac_path,
                            single_rad_jac_path,
                            single_tramat_path,
                            single_tramat_cumulative_path,
                            single_tramat_jac_path,
                            jac_targets,
                            ray_path,
                            single_freq_path,
                            atm_path,
                            single_propmat_path,
                            matpack::create(single_propmat_jac_path),
                            single_nlte_srcvec_path,
                            matpack::create(single_nlte_srcvec_jac_path),
                            surf_field,
                            atm_field,
                            hse_derivative);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radClearskyEmissionFrequencyDependentPropagation(
    const Workspace& ws,
    StokvecVector& spectral_rad,
    StokvecMatrix& spectral_rad_jac,
    ArrayOfArrayOfPropagationPathPoint& spectral_ray_path,
    const AtmField& atm_field,
    const AscendingGrid& freq_grid,
    const JacobianTargets& jac_targets,
    const Agenda& single_rad_space_agenda,
    const Agenda& single_rad_surface_agenda,
    const Agenda& single_propmat_agenda,
    const Agenda& ray_point_back_propagation_agenda,
    const SubsurfaceField& subsurf_field,
    const SurfaceField& surf_field,
    const Vector3& obs_pos,
    const Vector2& obs_los,
    const Numeric& max_stepsize,
    const Propmat& polarization,
    const Numeric& max_tau,
    const Numeric& cutoff_tau,
    const Index& hse_derivative,
    const Index& N) try {
  ARTS_TIME_REPORT

  const Size nf = freq_grid.size();
  spectral_rad.resize(nf);
  spectral_rad_jac.resize(jac_targets.x_size(), nf);
  spectral_rad = 0.0;
  spectral_ray_path.clear();
  spectral_ray_path.resize(nf);

  std::string error{};
  StokvecVector single_rad_jac;

#pragma omp parallel for firstprivate(single_rad_jac) if (arts_omp_parallel())
  for (Size f = 0; f < nf; f++) {
    try {
      single_radClearskyEmissionPropagation(ws,
                                            spectral_rad[f],
                                            single_rad_jac,
                                            spectral_ray_path[f],
                                            atm_field,
                                            freq_grid[f],
                                            jac_targets,
                                            single_rad_space_agenda,
                                            single_rad_surface_agenda,
                                            single_propmat_agenda,
                                            ray_point_back_propagation_agenda,
                                            subsurf_field,
                                            surf_field,
                                            obs_pos,
                                            obs_los,
                                            max_stepsize,
                                            polarization,
                                            max_tau,
                                            cutoff_tau,
                                            hse_derivative,
                                            N);

      spectral_rad_jac[joker, f] = single_rad_jac;
    } catch (std::exception& e) {
#pragma omp critical
      if (error.empty()) {
        error = std::format("Fail at index {}:\n{}", f, e.what());
      }
    }
  }

  ARTS_USER_ERROR_IF(not error.empty(), "{}", error)
}
ARTS_METHOD_ERROR_CATCH

void freq_gridFromSingleFrequency(AscendingGrid& freq_grid,
                                  const Numeric& frequency) try {
  ARTS_TIME_REPORT

  freq_grid = Vector{frequency};
}
ARTS_METHOD_ERROR_CATCH

void single_radFromVector(Stokvec& single_rad,
                          StokvecVector& single_rad_jac,
                          const StokvecVector& spectral_rad,
                          const StokvecMatrix& spectral_rad_jac,
                          const Index& f) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      spectral_rad.ncols() < f or spectral_rad_jac.ncols() < f,
      R"(Index f={} is out of bounds for spectral_rad with shape {:B,} and spectral_rad_jac with shape {:B,})",
      f,
      spectral_rad.shape(),
      spectral_rad_jac.shape());

  single_rad = spectral_rad[f];
  single_rad_jac.resize(spectral_rad_jac.nrows());
  single_rad_jac = spectral_rad_jac[joker, f];
}
ARTS_METHOD_ERROR_CATCH
