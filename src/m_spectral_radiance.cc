#include <array_algo.h>
#include <workspace.h>

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

  Vector ray_path_distance(N, 0.0);
  Tensor3 ray_path_distance_jacobian(2, N, nq, 0.0);

  for (Size ip = 1; ip < N; ip++) {
    ray_path_distance = path::distance(
        ray_path[ip - 1].pos, ray_path[ip].pos, surface_field.ellipsoid);
  }

  // FIXME: UNKNOWN ORDER OF "ip" AND "ip - 1" FOR TEMPERATURE
  if (hse_derivative and temperature_derivative_position >= 0) {
    for (Size ip = 1; ip < N; ip++) {
      ray_path_distance_jacobian[0, ip, temperature_derivative_position] =
          ray_path_distance[ip] /
          (2.0 * ray_path_atmospheric_point[ip - 1].temperature);
      ray_path_distance_jacobian[1, ip, temperature_derivative_position] =
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
