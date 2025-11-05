#include <array_algo.h>
#include <arts_omp.h>
#include <atm_field.h>
#include <auto_wsm.h>
#include <configtypes.h>
#include <jacobian.h>
#include <path_point.h>
#include <rtepack.h>
#include <rtepack_transmission.h>
#include <surface_field.h>
#include <time_report.h>
#include <workspace.h>

#include <ranges>

#include "auto_wsa.h"
#include "enumsAtmKey.h"
#include "matpack_mdspan_helpers_grid_t.h"
#include "physics_funcs.h"
#include "rtepack_stokes_vector.h"

void ray_path_transmission_matrixFromPath(
    ArrayOfMuelmatVector& ray_path_transmission_matrix,
    ArrayOfMuelmatTensor3& ray_path_transmission_matrix_jacobian,
    const ArrayOfPropmatVector& ray_path_propagation_matrix,
    const ArrayOfPropmatMatrix& ray_path_propagation_matrix_jacobian,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const SurfaceField& surface_field,
    const JacobianTargets& jacobian_targets,
    const Index& hse_derivative) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surface_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surface_field.ellipsoid)

  ARTS_USER_ERROR_IF(ray_path.size() == 0, "Empty path.");

  ARTS_USER_ERROR_IF(not arr::same_size(ray_path,
                                        ray_path_propagation_matrix,
                                        ray_path_propagation_matrix_jacobian,
                                        ray_path_atmospheric_point),
                     R"(Not same sizes:

ray_path.size()                             = {},
ray_path_propagation_matrix.size()          = {},
ray_path_propagation_matrix_jacobian.size() = {},
ray_path_atmospheric_point.size()           = {})",
                     ray_path.size(),
                     ray_path_propagation_matrix.size(),
                     ray_path_propagation_matrix_jacobian.size(),
                     ray_path_atmospheric_point.size());

  // HSE variables
  const Index temperature_derivative_position =
      jacobian_targets.target_position(AtmKey::t);

  const Size N = ray_path.size();

  ray_path_transmission_matrix.resize(N);
  ray_path_transmission_matrix_jacobian.resize(N);

  if (N == 0) return;

  const Index nq = jacobian_targets.target_count();

  const Vector ray_path_distance = distance(ray_path, surface_field.ellipsoid);
  Tensor3 ray_path_distance_jacobian(2, N - 1, nq, 0.0);

  // FIXME: UNKNOWN ORDER OF "ip" AND "ip - 1" FOR TEMPERATURE
  if (hse_derivative and temperature_derivative_position >= 0) {
    for (Size ip = 0; ip < N - 1; ip++) {
      ray_path_distance_jacobian[0, ip, temperature_derivative_position] =
          ray_path_distance[ip] /
          (2.0 * ray_path_atmospheric_point[ip].temperature);
      ray_path_distance_jacobian[1, ip, temperature_derivative_position] =
          ray_path_distance[ip] /
          (2.0 * ray_path_atmospheric_point[ip + 1].temperature);
    }
  }

  rtepack::two_level_exp(ray_path_transmission_matrix,
                         ray_path_transmission_matrix_jacobian,
                         ray_path_propagation_matrix,
                         ray_path_propagation_matrix_jacobian,
                         ray_path_distance,
                         ray_path_distance_jacobian);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceStepByStepEmission(
    StokvecVector& spectral_radiance,
    ArrayOfStokvecMatrix& ray_path_spectral_radiance_jacobian,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix_cumulative,
    const ArrayOfMuelmatTensor3& ray_path_transmission_matrix_jacobian,
    const ArrayOfStokvecVector& ray_path_spectral_radiance_source,
    const ArrayOfStokvecMatrix& ray_path_spectral_radiance_source_jacobian,
    const StokvecVector& spectral_radiance_background) try {
  ARTS_TIME_REPORT

  rtepack::two_level_linear_emission_step_by_step_full(
      spectral_radiance,
      ray_path_spectral_radiance_jacobian,
      ray_path_transmission_matrix,
      ray_path_transmission_matrix_cumulative,
      ray_path_transmission_matrix_jacobian,
      ray_path_spectral_radiance_source,
      ray_path_spectral_radiance_source_jacobian,
      spectral_radiance_background);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceCumulativeTransmission(
    StokvecVector& spectral_radiance,
    ArrayOfStokvecMatrix& ray_path_spectral_radiance_jacobian,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix_cumulative,
    const ArrayOfMuelmatTensor3& ray_path_transmission_matrix_jacobian,
    const StokvecVector& spectral_radiance_background) try {
  ARTS_TIME_REPORT

  rtepack::two_level_linear_transmission_step(
      spectral_radiance,
      ray_path_spectral_radiance_jacobian,
      ray_path_transmission_matrix,
      ray_path_transmission_matrix_cumulative,
      ray_path_transmission_matrix_jacobian,
      spectral_radiance_background);
}
ARTS_METHOD_ERROR_CATCH

namespace {
void ray_path_single_transmission_matrixFromPath(
    MuelmatVector& ray_path_single_transmission_matrix,
    MuelmatTensor3& ray_path_single_transmission_matrix_jacobian,
    const PropmatVector& ray_path_single_propagation_matrix,
    const PropmatMatrix& ray_path_single_propagation_matrix_jacobian,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const SurfaceField& surface_field,
    const JacobianTargets& jacobian_targets,
    const Index& hse_derivative) try {
  ARTS_TIME_REPORT

  const Size N  = ray_path.size();
  const Size nq = jacobian_targets.target_count();

  ARTS_USER_ERROR_IF(ray_path_single_propagation_matrix.size() != N or
                         ray_path_single_propagation_matrix_jacobian.nrows() !=
                             static_cast<Index>(N) or
                         ray_path_single_propagation_matrix_jacobian.ncols() !=
                             static_cast<Index>(nq) or
                         ray_path_atmospheric_point.size() != N,
                     R"(Not same sizes:

ray_path:                                    (N)     = [{}]
jacobian_targets:                            (nq)    = [{}]
ray_path_atmospheric_point:                  (N)     = [{}]
ray_path_single_propagation_matrix:          (N)     = {:B,}
ray_path_single_propagation_matrix_jacobian: (N, nq) = {:B,}
)",
                     N,
                     nq,
                     ray_path_atmospheric_point.size(),
                     ray_path_single_propagation_matrix.shape(),
                     ray_path_single_propagation_matrix_jacobian.shape());

  ray_path_single_transmission_matrix.resize(N);
  ray_path_single_transmission_matrix_jacobian.resize(N, 2, nq);
  ray_path_single_transmission_matrix          = 1.0;
  ray_path_single_transmission_matrix_jacobian = 0.0;

  const Vector ray_path_distance = distance(ray_path, surface_field.ellipsoid);
  Tensor3 ray_path_distance_jacobian(2, N - 1, nq, 0.0);

  const Index temperature_derivative_position =
      jacobian_targets.target_position(AtmKey::t);
  if (hse_derivative and temperature_derivative_position >= 0) {
    for (Size ip = 0; ip < N - 1; ip++) {
      ray_path_distance_jacobian[0, ip, temperature_derivative_position] =
          ray_path_distance[ip] /
          (2.0 * ray_path_atmospheric_point[ip].temperature);
      ray_path_distance_jacobian[1, ip, temperature_derivative_position] =
          ray_path_distance[ip] /
          (2.0 * ray_path_atmospheric_point[ip + 1].temperature);
    }
  }

  for (Size i = 0; i < N - 1; i++) {
    auto& T        = ray_path_single_transmission_matrix[i];
    auto dT1       = ray_path_single_transmission_matrix_jacobian[i, 0];
    auto dT2       = ray_path_single_transmission_matrix_jacobian[i + 1, 1];
    const auto& K1 = ray_path_single_propagation_matrix[i];
    const auto& K2 = ray_path_single_propagation_matrix[i + 1];
    const auto dK1 = ray_path_single_propagation_matrix_jacobian[i];
    const auto dK2 = ray_path_single_propagation_matrix_jacobian[i + 1];
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

void ray_path_single_transmission_matrix_cumulativeFromPath(
    MuelmatVector& ray_path_single_transmission_matrix_cumulative,
    const MuelmatVector& ray_path_single_transmission_matrix) try {
  ARTS_TIME_REPORT

  const Size N = ray_path_single_transmission_matrix.size();

  ray_path_single_transmission_matrix_cumulative.resize(N);

  if (N == 0) return;

  auto& Pi = ray_path_single_transmission_matrix_cumulative;
  auto& T  = ray_path_single_transmission_matrix;

  Pi.front() = T.front();
  for (Size i = 1; i < N; i++) Pi[i] = Pi[i - 1] * T[i];
}
ARTS_METHOD_ERROR_CATCH

void ray_path_single_radiance_sourceFromPropmat(
    StokvecVector& ray_path_single_radiance_source,
    StokvecMatrix& ray_path_single_radiance_source_jacobian,
    const PropmatVector& ray_path_single_propagation_matrix,
    const StokvecVector& ray_path_single_propagation_matrix_nonlte,
    const PropmatMatrix& ray_path_single_propagation_matrix_jacobian,
    const StokvecMatrix& ray_path_single_propagation_matrix_nonlte_jacobian,
    const Vector& ray_path_single_frequency,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const JacobianTargets& jacobian_targets) try {
  ARTS_TIME_REPORT

  const Index np = ray_path_atmospheric_point.size();
  const Index nq = jacobian_targets.target_count();

  ARTS_USER_ERROR_IF(
      np != ray_path_single_frequency.ncols() or
          np != ray_path_single_propagation_matrix.ncols() or
          np != ray_path_single_propagation_matrix_nonlte.ncols() or
          np != ray_path_single_propagation_matrix_jacobian.nrows() or
          nq != ray_path_single_propagation_matrix_jacobian.ncols() or
          np != ray_path_single_propagation_matrix_nonlte_jacobian.nrows() or
          nq != ray_path_single_propagation_matrix_nonlte_jacobian.ncols(),
      R"(Not same sizes:

jacobian_targets                                    (nq)     = [{}]
ray_path_single_frequency                           (np)     = {:B,}
ray_path_atmospheric_point                          (np)     = [{}]
ray_path_single_propagation_matrix:                 (np)     = {:B,}
ray_path_single_propagation_matrix_nonlte:          (np)     = {:B,}
ray_path_single_propagation_matrix_jacobian:        (np, nq) = {:B,}
ray_path_single_propagation_matrix_nonlte_jacobian: (np, nq) = {:B,}
)",
      nq,
      ray_path_single_frequency.shape(),
      np,
      ray_path_single_propagation_matrix.shape(),
      ray_path_single_propagation_matrix_nonlte.shape(),
      ray_path_single_propagation_matrix_jacobian.shape(),
      ray_path_single_propagation_matrix_nonlte_jacobian.shape())

  ray_path_single_radiance_source.resize(np);
  ray_path_single_radiance_source_jacobian.resize(np, nq);
  ray_path_single_radiance_source          = 0.0;
  ray_path_single_radiance_source_jacobian = 0.0;

  if (np == 0) return;

  const Index it = jacobian_targets.target_position(AtmKey::t);

  for (Index i = 0; i < np; i++) {
    auto& J       = ray_path_single_radiance_source[i];
    auto dJ       = ray_path_single_radiance_source_jacobian[i];
    const auto& K = ray_path_single_propagation_matrix[i];
    const auto& S = ray_path_single_propagation_matrix_nonlte[i];
    const auto& f = ray_path_single_frequency[i];
    const auto t  = ray_path_atmospheric_point[i].temperature;
    const auto dK = ray_path_single_propagation_matrix_jacobian[i];
    const auto dS = ray_path_single_propagation_matrix_nonlte_jacobian[i];

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

void single_radianceStepByStepEmission(
    Stokvec& single_radiance,
    StokvecMatrix& ray_path_single_radiance_jacobian,
    const MuelmatVector& ray_path_single_transmission_matrix,
    const MuelmatVector& ray_path_single_transmission_matrix_cumulative,
    const MuelmatTensor3& ray_path_single_transmission_matrix_jacobian,
    const StokvecVector& ray_path_single_radiance_source,
    const StokvecMatrix& ray_path_single_radiance_source_jacobian) try {
  ARTS_TIME_REPORT

  const Size N  = ray_path_single_transmission_matrix.size();
  const Size nq = ray_path_single_radiance_source_jacobian.ncols();

  ARTS_USER_ERROR_IF(
      ray_path_single_transmission_matrix.size() != N or
          ray_path_single_transmission_matrix_cumulative.size() != N or
          ray_path_single_radiance_source.size() != N or
          ray_path_single_radiance_source_jacobian.nrows() !=
              static_cast<Index>(N) or
          ray_path_single_radiance_source_jacobian.ncols() !=
              static_cast<Index>(nq) or
          ray_path_single_transmission_matrix_jacobian.npages() !=
              static_cast<Index>(N) or
          ray_path_single_transmission_matrix_jacobian.nrows() != 2,
      R"(Not same sizes:

ray_path_single_radiance_source                (N)        = {:B,}
ray_path_single_transmission_matrix            (N)        = {:B,}
ray_path_single_radiance_source_jacobian       (N, nq)    = {:B,}
ray_path_single_transmission_matrix_jacobian   (N, 2, nq) = {:B,}
ray_path_single_transmission_matrix_cumulative (N)        = {:B,}
)",
      ray_path_single_radiance_source.shape(),
      ray_path_single_transmission_matrix.shape(),
      ray_path_single_radiance_source_jacobian.shape(),
      ray_path_single_transmission_matrix_jacobian.shape(),
      ray_path_single_transmission_matrix_cumulative.shape());

  ray_path_single_radiance_jacobian.resize(N, nq);
  ray_path_single_radiance_jacobian = 0.0;

  if (N == 0) return;

  auto& I         = single_radiance;
  auto& dI        = ray_path_single_radiance_jacobian;
  const auto& Js  = ray_path_single_radiance_source;
  const auto& dJs = ray_path_single_radiance_source_jacobian;
  const auto& Pi  = ray_path_single_transmission_matrix_cumulative;
  const auto& Ts  = ray_path_single_transmission_matrix;
  const auto& dTs = ray_path_single_transmission_matrix_jacobian;

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

void single_radiance_jacobianFromBackground(
    StokvecVector& single_radiance_jacobian,
    const Muelmat& background_transmittance) try {
  ARTS_TIME_REPORT
  for (auto& dI : single_radiance_jacobian) dI = background_transmittance * dI;
}
ARTS_METHOD_ERROR_CATCH

void single_radiance_jacobianAddPathPropagation(
    StokvecVector& single_radiance_jacobian,
    const StokvecMatrix& ray_path_single_radiance_jacobian,
    const JacobianTargets& jacobian_targets,
    const AtmField& atmospheric_field,
    const ArrayOfPropagationPathPoint& ray_path) try {
  ARTS_TIME_REPORT

  const Size np  = ray_path.size();
  const Index nt = jacobian_targets.target_count();

  ARTS_USER_ERROR_IF(
      single_radiance_jacobian.ncols() != nt or
          ray_path_single_radiance_jacobian.nrows() != static_cast<Index>(np) or
          ray_path_single_radiance_jacobian.ncols() != nt,
      R"(Not same sizes:

ray_path:                            (np)     = [{}]
jacobian_targets:                    (nq)     = [{}]
single_radiance_jacobian:            (nq)     = {:B,}
ray_path_single_radiance_jacobian:   (np, nq) = {:B,}
)",
      np,
      jacobian_targets.target_count(),
      single_radiance_jacobian.shape(),
      ray_path_single_radiance_jacobian.shape());

  for (auto& atm_block : jacobian_targets.atm) {
    const auto& data = atmospheric_field[atm_block.type];
    for (Size ip = 0; ip < np; ip++) {
      const auto weights = data.flat_weight(ray_path[ip].pos);
      const auto local   = ray_path_single_radiance_jacobian[ip];

      for (auto& w : weights) {
        if (w.second != 0.0) {
          const auto i                 = w.first + atm_block.x_start;
          single_radiance_jacobian[i] += w.second * local[atm_block.target_pos];
        }
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void single_radianceFromPropagation(
    Stokvec& single_radiance,
    StokvecVector& single_radiance_jacobian,
    StokvecVector& ray_path_single_radiance_source,
    StokvecMatrix& ray_path_single_radiance_source_jacobian,
    StokvecMatrix& ray_path_single_radiance_jacobian,
    MuelmatVector& ray_path_single_transmission_matrix,
    MuelmatVector& ray_path_single_transmission_matrix_cumulative,
    MuelmatTensor3& ray_path_single_transmission_matrix_jacobian,
    const JacobianTargets& jacobian_targets,
    const ArrayOfPropagationPathPoint& ray_path,
    const Vector& ray_path_single_frequency,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const PropmatVector& ray_path_single_propagation_matrix,
    const PropmatMatrix& ray_path_single_propagation_matrix_jacobian,
    const StokvecVector& ray_path_single_propagation_matrix_nonlte,
    const StokvecMatrix& ray_path_single_propagation_matrix_nonlte_jacobian,
    const SurfaceField& surface_field,
    const AtmField& atmospheric_field,
    const Index& hse_derivative) try {
  ARTS_TIME_REPORT

  ray_path_single_transmission_matrixFromPath(
      ray_path_single_transmission_matrix,
      ray_path_single_transmission_matrix_jacobian,
      ray_path_single_propagation_matrix,
      ray_path_single_propagation_matrix_jacobian,
      ray_path,
      ray_path_atmospheric_point,
      surface_field,
      jacobian_targets,
      hse_derivative);
  ray_path_single_transmission_matrix_cumulativeFromPath(
      ray_path_single_transmission_matrix_cumulative,
      ray_path_single_transmission_matrix);
  const auto& transmission_single_matrix_background =
      ray_path_single_transmission_matrix_cumulative.back();
  ray_path_single_radiance_sourceFromPropmat(
      ray_path_single_radiance_source,
      ray_path_single_radiance_source_jacobian,
      ray_path_single_propagation_matrix,
      ray_path_single_propagation_matrix_nonlte,
      ray_path_single_propagation_matrix_jacobian,
      ray_path_single_propagation_matrix_nonlte_jacobian,
      ray_path_single_frequency,
      ray_path_atmospheric_point,
      jacobian_targets);
  single_radianceStepByStepEmission(
      single_radiance,
      ray_path_single_radiance_jacobian,
      ray_path_single_transmission_matrix,
      ray_path_single_transmission_matrix_cumulative,
      ray_path_single_transmission_matrix_jacobian,
      ray_path_single_radiance_source,
      ray_path_single_radiance_source_jacobian);
  single_radiance_jacobianFromBackground(single_radiance_jacobian,
                                         transmission_single_matrix_background);
  single_radiance_jacobianAddPathPropagation(single_radiance_jacobian,
                                             ray_path_single_radiance_jacobian,
                                             jacobian_targets,
                                             atmospheric_field,
                                             ray_path);
}
ARTS_METHOD_ERROR_CATCH
}  // namespace

void spectral_radianceSetToBackground(
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_jacobian,
    const StokvecVector& spectral_radiance_background,
    const StokvecMatrix& spectral_radiance_background_jacobian) try {
  ARTS_TIME_REPORT

  spectral_radiance          = spectral_radiance_background;
  spectral_radiance_jacobian = spectral_radiance_background_jacobian;
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceSinglePathEmissionFrequencyLoop(
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_jacobian,
    const JacobianTargets& jacobian_targets,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfAscendingGrid& ray_path_frequency_grid,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const ArrayOfPropmatVector& ray_path_propagation_matrix,
    const ArrayOfStokvecVector&
        ray_path_propagation_matrix_source_vector_nonlte,
    const ArrayOfPropmatMatrix& ray_path_propagation_matrix_jacobian,
    const ArrayOfStokvecMatrix&
        ray_path_propagation_matrix_source_vector_nonlte_jacobian,
    const SurfaceField& surface_field,
    const AtmField& atmospheric_field,
    const Index& hse_derivative) try {
  ARTS_TIME_REPORT

  const Size nf = spectral_radiance.size();
  const Size N  = ray_path.size();
  const Size nq = jacobian_targets.target_count();
  const Size nx = jacobian_targets.x_size();

  if (N == 0) return;

  ARTS_USER_ERROR_IF(
      not matpack::same_shape<2>({nx, nf}, spectral_radiance_jacobian),
      R"(spectral_radiance_jacobian has wrong shape:
Expected: ({}, {})
Actual:   ({}, {})
)",
      nx,
      nf,
      spectral_radiance_jacobian.nrows(),
      spectral_radiance_jacobian.ncols());

  ARTS_USER_ERROR_IF(
      not arr::same_size(
          ray_path,
          ray_path_frequency_grid,
          ray_path_atmospheric_point,
          ray_path_propagation_matrix,
          ray_path_propagation_matrix_source_vector_nonlte,
          ray_path_propagation_matrix_jacobian,
          ray_path_propagation_matrix_source_vector_nonlte_jacobian),
      R"(Not same sizes:

ray_path.size()                                                  = {},
ray_path_frequency_grid.size()                                   = {},
ray_path_atmospheric_point.size()                                = {},
ray_path_propagation_matrix.size()                               = {},
ray_path_propagation_matrix_source_vector_nonlte.size()          = {},
ray_path_propagation_matrix_jacobian.size()                      = {},
ray_path_propagation_matrix_source_vector_nonlte_jacobian.size() = {}
)",
      ray_path.size(),
      ray_path_frequency_grid.size(),
      ray_path_atmospheric_point.size(),
      ray_path_propagation_matrix.size(),
      ray_path_propagation_matrix_source_vector_nonlte.size(),
      ray_path_propagation_matrix_jacobian.size(),
      ray_path_propagation_matrix_source_vector_nonlte_jacobian.size());

  ARTS_USER_ERROR_IF(
      not arr::elemwise_same_size(
          ray_path_frequency_grid,
          ray_path_propagation_matrix,
          ray_path_propagation_matrix_source_vector_nonlte),
      R"(Not same sizes elemwise:
ray_path_frequency_grid.size()                           = {}
ray_path_propagation_matrix.shape()                      = {:B,}
ray_path_propagation_matrix_source_vector_nonlte.shape() = {:B,}
)",
      ray_path_frequency_grid.size(),
      ray_path_propagation_matrix[0].shape(),
      ray_path_propagation_matrix_source_vector_nonlte[0].shape());

  ARTS_USER_ERROR_IF(
      not arr::elemwise_same_size(
          ray_path_propagation_matrix_jacobian,
          ray_path_propagation_matrix_source_vector_nonlte_jacobian),
      R"(Not same sizes elemwise:
ray_path_propagation_matrix_jacobian.shape()                      = {:B,}
ray_path_propagation_matrix_source_vector_nonlte_jacobian.shape() = {:B,}
)",
      ray_path_propagation_matrix_jacobian[0].shape(),
      ray_path_propagation_matrix_source_vector_nonlte_jacobian[0].shape());

  Vector ray_path_single_frequency(N);
  StokvecVector single_radiance_jacobian(nx);
  StokvecVector ray_path_single_radiance_source(N);
  StokvecMatrix ray_path_single_radiance_source_jacobian(N, nq);
  StokvecMatrix ray_path_single_radiance_jacobian(N, nq);
  MuelmatVector ray_path_single_transmission_matrix(N);
  MuelmatVector ray_path_single_transmission_matrix_cumulative(N);
  MuelmatTensor3 ray_path_single_transmission_matrix_jacobian(N, 2, nq);
  PropmatVector ray_path_single_propagation_matrix(N);
  PropmatMatrix ray_path_single_propagation_matrix_jacobian(N, nq);
  StokvecVector ray_path_single_propagation_matrix_nonlte(N);
  StokvecMatrix ray_path_single_propagation_matrix_nonlte_jacobian(N, nq);

  std::string error{};

#pragma omp parallel for firstprivate(                      \
        ray_path_single_frequency,                          \
            single_radiance_jacobian,                       \
            ray_path_single_radiance_source,                \
            ray_path_single_radiance_source_jacobian,       \
            ray_path_single_radiance_jacobian,              \
            ray_path_single_transmission_matrix,            \
            ray_path_single_transmission_matrix_cumulative, \
            ray_path_single_transmission_matrix_jacobian,   \
            ray_path_single_propagation_matrix,             \
            ray_path_single_propagation_matrix_jacobian,    \
            ray_path_single_propagation_matrix_nonlte,      \
            ray_path_single_propagation_matrix_nonlte_jacobian) if (arts_omp_parallel())
  for (Size f = 0; f < nf; f++) {
    try {
      Stokvec& single_radiance = spectral_radiance[f];
      for (Size i = 0; i < N; i++) {
        ray_path_single_frequency[i] = ray_path_frequency_grid[i][f];
        ray_path_single_propagation_matrix[i] =
            ray_path_propagation_matrix[i][f];
        ray_path_single_propagation_matrix_nonlte[i] =
            ray_path_propagation_matrix_source_vector_nonlte[i][f];

        ray_path_single_propagation_matrix_jacobian[i] =
            ray_path_propagation_matrix_jacobian[i][joker, f];
        ray_path_single_propagation_matrix_nonlte_jacobian[i] =
            ray_path_propagation_matrix_source_vector_nonlte_jacobian[i]
                                                                     [joker, f];
      }

      single_radiance_jacobian = spectral_radiance_jacobian[joker, f];

      single_radianceFromPropagation(
          single_radiance,
          single_radiance_jacobian,
          ray_path_single_radiance_source,
          ray_path_single_radiance_source_jacobian,
          ray_path_single_radiance_jacobian,
          ray_path_single_transmission_matrix,
          ray_path_single_transmission_matrix_cumulative,
          ray_path_single_transmission_matrix_jacobian,
          jacobian_targets,
          ray_path,
          ray_path_single_frequency,
          ray_path_atmospheric_point,
          ray_path_single_propagation_matrix,
          ray_path_single_propagation_matrix_jacobian,
          ray_path_single_propagation_matrix_nonlte,
          ray_path_single_propagation_matrix_nonlte_jacobian,
          surface_field,
          atmospheric_field,
          hse_derivative);

      spectral_radiance_jacobian[joker, f] = single_radiance_jacobian;
    } catch (std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(
      not error.empty(),
      "Error in spectral_radianceSinglePathEmissionFrequencyLoop: {:s}",
      error);
}
ARTS_METHOD_ERROR_CATCH

void single_spectral_radianceClearskyEmissionPropagation(
    const Workspace& ws,
    Stokvec& single_spectral_radiance,
    StokvecVector& single_spectral_radiance_jacobian,
    ArrayOfPropagationPathPoint& ray_path,
    const AtmField& atmospheric_field,
    const Numeric& frequency_,
    const JacobianTargets& jacobian_targets,
    const Agenda& single_spectral_radiance_space_agenda,
    const Agenda& single_spectral_radiance_surface_agenda,
    const Agenda& propagation_matrix_single_agenda,
    const Agenda& ray_path_point_back_propagation_agenda,
    const SubsurfaceField& subsurface_field,
    const SurfaceField& surface_field,
    const Vector3& spectral_radiance_observer_position,
    const Vector2& spectral_radiance_observer_line_of_sight,
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

  const Size nt = jacobian_targets.x_size();
  single_spectral_radiance_jacobian.resize(nt);
  single_spectral_radiance_jacobian = 0.0;

  Numeric dispersion_single;
  Vector dispersion_single_jacobian;

  Vector ray_path_single_frequency;
  ArrayOfAtmPoint ray_path_atmospheric_point;
  PropmatVector ray_path_single_propagation_matrix;
  ArrayOfPropmatVector ray_path_single_propagation_matrix_jacobian;
  StokvecVector ray_path_single_propagation_matrix_nonlte;
  ArrayOfStokvecVector ray_path_single_propagation_matrix_nonlte_jacobian;

  ray_path_single_frequency.reserve(N);
  ray_path_atmospheric_point.reserve(N);
  ray_path_single_propagation_matrix.reserve(N);
  ray_path_single_propagation_matrix_jacobian.reserve(N);
  ray_path_single_propagation_matrix_nonlte.reserve(N);
  ray_path_single_propagation_matrix_nonlte_jacobian.reserve(N);

  ray_path.clear();
  ray_path.reserve(N);
  ray_path.resize(1);
  ray_path[0] =
      path::init_with_lostype(spectral_radiance_observer_position,
                              spectral_radiance_observer_line_of_sight,
                              atmospheric_field,
                              surface_field,
                              true);

  if (ray_path.back().los_type == PathPositionType::space) {
    single_spectral_radiance_space_agendaExecute(
        ws,
        single_spectral_radiance,
        single_spectral_radiance_jacobian,
        frequency_,
        jacobian_targets,
        ray_path.back(),
        single_spectral_radiance_space_agenda);
    return;
  }

  if (ray_path.back().los_type == PathPositionType::surface) {
    single_spectral_radiance_surface_agendaExecute(
        ws,
        single_spectral_radiance,
        single_spectral_radiance_jacobian,
        frequency_,
        jacobian_targets,
        ray_path.back(),
        surface_field,
        subsurface_field,
        single_spectral_radiance_surface_agenda);
    return;
  }

  Numeric sumAr = 0.0;
  while (true) {
    const auto& ray_path_point = ray_path.back();
    ray_path_atmospheric_point.emplace_back(
        atmospheric_field.at(ray_path_point.pos));

    Numeric frequency = frequency_;
    Vector3 frequency_wind_shift_jacobian;
    frequencyWindShift(frequency,
                       frequency_wind_shift_jacobian,
                       ray_path_atmospheric_point.back(),
                       ray_path_point);
    ray_path_single_frequency.push_back(frequency);

    propagation_matrix_single_agendaExecute(
        ws,
        ray_path_single_propagation_matrix.emplace_back(),
        ray_path_single_propagation_matrix_nonlte.emplace_back(),
        dispersion_single,
        ray_path_single_propagation_matrix_jacobian.emplace_back(),
        ray_path_single_propagation_matrix_nonlte_jacobian.emplace_back(),
        dispersion_single_jacobian,
        frequency,
        frequency_wind_shift_jacobian,
        jacobian_targets,
        "AIR"_spec,
        ray_path_point,
        ray_path_atmospheric_point.back(),
        propagation_matrix_single_agenda);
    dispersion_single +=
        dot(polarization, ray_path_single_propagation_matrix.back());

    const Numeric A = ray_path_single_propagation_matrix.back().A();
    if (A == 0.0) {
      PropagationPathPoint tmp_point;
      ray_path_point_back_propagation_agendaExecute(
          ws,
          tmp_point,
          ray_path,
          dispersion_single,
          ray_path_single_propagation_matrix.back(),
          max_stepsize,
          ray_path_point_back_propagation_agenda);
      ray_path.push_back(tmp_point);
    } else {
      PropagationPathPoint tmp_point;
      ray_path_point_back_propagation_agendaExecute(
          ws,
          tmp_point,
          ray_path,
          dispersion_single,
          ray_path_single_propagation_matrix.back(),
          std::min(max_tau / A, max_stepsize),
          ray_path_point_back_propagation_agenda);
      ray_path.push_back(tmp_point);
    }

    if (ray_path.back().los_type == PathPositionType::space) {
      single_spectral_radiance_space_agendaExecute(
          ws,
          single_spectral_radiance,
          single_spectral_radiance_jacobian,
          frequency,
          jacobian_targets,
          ray_path.back(),
          single_spectral_radiance_space_agenda);
      break;
    }

    if (ray_path.back().los_type == PathPositionType::surface) {
      single_spectral_radiance_surface_agendaExecute(
          ws,
          single_spectral_radiance,
          single_spectral_radiance_jacobian,
          frequency,
          jacobian_targets,
          ray_path.back(),
          surface_field,
          subsurface_field,
          single_spectral_radiance_surface_agenda);
      break;
    }

    sumAr += A * path::distance(ray_path.back().pos,
                                ray_path[ray_path.size() - 2].pos,
                                surface_field.ellipsoid);

    if (sumAr >= cutoff_tau) {
      single_spectral_radiance =
          planck(frequency_, ray_path_atmospheric_point.back().temperature);
      break;
    }
  }

  ray_path.pop_back();
  StokvecVector ray_path_single_radiance_source;
  StokvecMatrix ray_path_single_radiance_source_jacobian;
  StokvecMatrix ray_path_single_radiance_jacobian;
  MuelmatVector ray_path_single_transmission_matrix;
  MuelmatVector ray_path_single_transmission_matrix_cumulative;
  MuelmatTensor3 ray_path_single_transmission_matrix_jacobian;
  single_radianceFromPropagation(
      single_spectral_radiance,
      single_spectral_radiance_jacobian,
      ray_path_single_radiance_source,
      ray_path_single_radiance_source_jacobian,
      ray_path_single_radiance_jacobian,
      ray_path_single_transmission_matrix,
      ray_path_single_transmission_matrix_cumulative,
      ray_path_single_transmission_matrix_jacobian,
      jacobian_targets,
      ray_path,
      ray_path_single_frequency,
      ray_path_atmospheric_point,
      ray_path_single_propagation_matrix,
      matpack::create(ray_path_single_propagation_matrix_jacobian),
      ray_path_single_propagation_matrix_nonlte,
      matpack::create(ray_path_single_propagation_matrix_nonlte_jacobian),
      surface_field,
      atmospheric_field,
      hse_derivative);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceClearskyEmissionFrequencyDependentPropagation(
    const Workspace& ws,
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_jacobian,
    ArrayOfArrayOfPropagationPathPoint& ray_paths,
    const AtmField& atmospheric_field,
    const AscendingGrid& frequency_grid,
    const JacobianTargets& jacobian_targets,
    const Agenda& single_spectral_radiance_space_agenda,
    const Agenda& single_spectral_radiance_surface_agenda,
    const Agenda& propagation_matrix_single_agenda,
    const Agenda& ray_path_point_back_propagation_agenda,
    const SubsurfaceField& subsurface_field,
    const SurfaceField& surface_field,
    const Vector3& spectral_radiance_observer_position,
    const Vector2& spectral_radiance_observer_line_of_sight,
    const Numeric& max_stepsize,
    const Propmat& polarization,
    const Numeric& max_tau,
    const Numeric& cutoff_tau,
    const Index& hse_derivative,
    const Index& N) try {
  ARTS_TIME_REPORT

  const Size nf = frequency_grid.size();
  spectral_radiance.resize(nf);
  spectral_radiance_jacobian.resize(jacobian_targets.x_size(), nf);
  spectral_radiance = 0.0;
  ray_paths.clear();
  ray_paths.resize(nf);

  std::string error{};
  StokvecVector single_spectral_radiance_jacobian;

#pragma omp parallel for firstprivate( \
        single_spectral_radiance_jacobian) if (arts_omp_parallel())
  for (Size f = 0; f < nf; f++) {
    try {
      single_spectral_radianceClearskyEmissionPropagation(
          ws,
          spectral_radiance[f],
          single_spectral_radiance_jacobian,
          ray_paths[f],
          atmospheric_field,
          frequency_grid[f],
          jacobian_targets,
          single_spectral_radiance_space_agenda,
          single_spectral_radiance_surface_agenda,
          propagation_matrix_single_agenda,
          ray_path_point_back_propagation_agenda,
          subsurface_field,
          surface_field,
          spectral_radiance_observer_position,
          spectral_radiance_observer_line_of_sight,
          max_stepsize,
          polarization,
          max_tau,
          cutoff_tau,
          hse_derivative,
          N);

      spectral_radiance_jacobian[joker, f] = single_spectral_radiance_jacobian;
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

void frequency_gridFromSingleFrequency(AscendingGrid& frequency_grid,
                                       const Numeric& frequency) try {
  ARTS_TIME_REPORT

  frequency_grid = Vector{frequency};
}
ARTS_METHOD_ERROR_CATCH

void single_spectral_radianceFromVector(
    Stokvec& single_spectral_radiance,
    StokvecVector& single_spectral_radiance_jacobian,
    const StokvecVector& spectral_radiance,
    const StokvecMatrix& spectral_radiance_jacobian,
    const Index& f) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      spectral_radiance.ncols() < f or spectral_radiance_jacobian.ncols() < f,
      R"(Index f={} is out of bounds for spectral_radiance with shape {:B,} and spectral_radiance_jacobian with shape {:B,})",
      f,
      spectral_radiance.shape(),
      spectral_radiance_jacobian.shape());

  single_spectral_radiance = spectral_radiance[f];
  single_spectral_radiance_jacobian.resize(spectral_radiance_jacobian.nrows());
  single_spectral_radiance_jacobian = spectral_radiance_jacobian[joker, f];
}
ARTS_METHOD_ERROR_CATCH
