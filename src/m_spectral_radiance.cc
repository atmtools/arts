#include <workspace.h>

#include <algorithm>

#include "configtypes.h"
#include "debug.h"
#include "rtepack.h"

#include "mh_checks.h"

void spectral_radianceStepByStep(
    StokvecVector& spectral_radiance,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix,
    const ArrayOfStokvecVector& ray_path_spectral_radiance_source,
    const StokvecVector& spectral_radiance_background) try {
  ARTS_USER_ERROR_IF(not all_same_shape(spectral_radiance_background.shape(),
                                        ray_path_spectral_radiance_source,
                                        ray_path_transmission_matrix),
                     std::format(
                         R"(Not same inner shape: {:B,}

Variables: spectral_radiance_background
           ray_path_spectral_radiance_source
           ray_path_transmission_matrix_cumulative)",
                         spectral_radiance_background.shape()));

  const Size N   = ray_path_transmission_matrix.size();
  const Index nv = spectral_radiance_background.size();

  ARTS_USER_ERROR_IF(ray_path_spectral_radiance_source.size() != N,
                     std::format(
                         R"(Size mismatch.

ray_path_spectral_radiance_source.size() = {},
ray_path_transmission_matrix.size()      = {})",
                         ray_path_spectral_radiance_source.size(),
                         ray_path_transmission_matrix.size()));

  constexpr auto step = [](Stokvec I0, Stokvec J, Muelmat T) {
    return J + T * (I0 - J);
  };

  spectral_radiance = spectral_radiance_background;
  for (Size i = N - 2; i < N; i--) {
    for (Index iv = 0; iv < nv; iv++) {
      spectral_radiance[iv] =
          step(spectral_radiance[iv],
               0.5 * (ray_path_spectral_radiance_source[i][iv] +
                      ray_path_spectral_radiance_source[i + 1][iv]),
               ray_path_transmission_matrix[i + 1][iv]);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

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
  ARTS_USER_ERROR_IF(ray_path.size() == 0, "Empty path.");

  ARTS_USER_ERROR_IF(not all_same_size(ray_path,
                                       ray_path_propagation_matrix,
                                       ray_path_propagation_matrix_jacobian,
                                       ray_path_atmospheric_point),
                     std::format(
                         R"(Not same sizes:

ray_path.size()                             = {},
ray_path_propagation_matrix.size()          = {},
ray_path_propagation_matrix_jacobian.size() = {},
ray_path_atmospheric_point.size()           = {})",
                         ray_path.size(),
                         ray_path_propagation_matrix.size(),
                         ray_path_propagation_matrix_jacobian.size(),
                         ray_path_atmospheric_point.size()));

  // HSE variables
  const Index temperature_derivative_position =
      jacobian_targets.target_position<Jacobian::AtmTarget>(AtmKey::t);

  const Size N = ray_path.size();

  ray_path_transmission_matrix.resize(N);
  ray_path_transmission_matrix_jacobian.resize(N);

  if (N == 0) return;

  const Index nq = jacobian_targets.target_count();

  Vector ray_path_distance(N, 0.0);
  Tensor3 ray_path_distance_jacobian(2, N, nq, 0.0);

  for (Size ip = 1; ip < N; ip++) {
    ray_path_distance = path::distance(
        ray_path[ip - 1].pos, ray_path[ip].pos, surface_field.ellipsoid);
  }

  // FIXME: UNKNOWN ORDER OF "ip" AND "ip - 1" FOR TEMPERATURE
  if (hse_derivative and temperature_derivative_position >= 0) {
    for (Size ip = 1; ip < N; ip++) {
      ray_path_distance_jacobian(0, ip, temperature_derivative_position) =
          ray_path_distance[ip] /
          (2.0 * ray_path_atmospheric_point[ip - 1].temperature);
      ray_path_distance_jacobian(1, ip, temperature_derivative_position) =
          ray_path_distance[ip] /
          (2.0 * ray_path_atmospheric_point[ip].temperature);
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

void spectral_radianceCumulativeEmission(
    StokvecVector& spectral_radiance,
    ArrayOfStokvecMatrix& ray_path_spectral_radiance_jacobian,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix_cumulative,
    const ArrayOfMuelmatTensor3& ray_path_transmission_matrix_jacobian,
    const ArrayOfStokvecVector& ray_path_spectral_radiance_source,
    const ArrayOfStokvecMatrix& ray_path_spectral_radiance_source_jacobian,
    const StokvecVector& spectral_radiance_background) try {
  rtepack::two_level_linear_emission_step(
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
  rtepack::two_level_linear_transmission_step(
      spectral_radiance,
      ray_path_spectral_radiance_jacobian,
      ray_path_transmission_matrix,
      ray_path_transmission_matrix_cumulative,
      ray_path_transmission_matrix_jacobian,
      spectral_radiance_background);
}
ARTS_METHOD_ERROR_CATCH
