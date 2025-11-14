#include <arts_conversions.h>
#include <arts_omp.h>
#include <atm_path.h>
#include <debug.h>
#include <jacobian.h>
#include <path_point.h>
#include <rtepack.h>
#include <surf.h>
#include <workspace.h>

#include <stdexcept>

void ray_path_zeeman_magnetic_fieldFromPath(
    ArrayOfVector3 &ray_path_zeeman_magnetic_field,
    const ArrayOfPropagationPathPoint &ray_path,
    const ArrayOfAtmPoint &ray_path_atm_point) try {
  ARTS_TIME_REPORT

  const Size np = ray_path_atm_point.size();
  ARTS_USER_ERROR_IF(np != ray_path.size(),
                     "ray_path and ray_path_atm_point must have the same size")

  ray_path_zeeman_magnetic_field.resize(np);

  std::transform(ray_path.begin(),
                 ray_path.end(),
                 ray_path_atm_point.begin(),
                 ray_path_zeeman_magnetic_field.begin(),
                 [](const auto &p, const auto &a) -> Vector3 {
                   using Conversion::rad2deg;
                   const auto zz = lbl::zeeman::magnetic_angles(a.mag, p.los);
                   return {zz.H, rad2deg(zz.theta()), rad2deg(zz.eta())};
                 });
}
ARTS_METHOD_ERROR_CATCH

void ray_path_spectral_radiance_sourceFromPropmat(
    ArrayOfStokvecVector &ray_path_spectral_radiance_source,
    ArrayOfStokvecMatrix &ray_path_spectral_radiance_source_jacobian,
    const ArrayOfPropmatVector &ray_path_propagation_matrix,
    const ArrayOfStokvecVector &ray_path_source_vector_nonlte,
    const ArrayOfPropmatMatrix &ray_path_propagation_matrix_jacobian,
    const ArrayOfStokvecMatrix &ray_path_source_vector_nonlte_jacobian,
    const ArrayOfAscendingGrid &freq_grid_path,
    const ArrayOfAtmPoint &ray_path_atm_point,
    const JacobianTargets &jacobian_targets) try {
  ARTS_TIME_REPORT

  String error{};

  const Index np = ray_path_atm_point.size();
  if (np == 0) {
    ray_path_spectral_radiance_source.resize(0);
    ray_path_spectral_radiance_source_jacobian.resize(0);
    return;
  }

  const Index nf = ray_path_propagation_matrix.front().size();
  const Index nq = jacobian_targets.target_count();

  const Index it = jacobian_targets.target_position(AtmKey::t);

  ray_path_spectral_radiance_source.resize(np);
  for (auto &t : ray_path_spectral_radiance_source) {
    t.resize(nf);
    t = 0;
  }
  ray_path_spectral_radiance_source_jacobian.resize(np);
  for (auto &t : ray_path_spectral_radiance_source_jacobian) {
    t.resize(nq, nf);
    t = 0;
  }

  // Loop ppath points and determine radiative properties
#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Index ip = 0; ip < np; ip++) {
    try {
      rtepack::source::level_nlte(
          ray_path_spectral_radiance_source[ip],
          ray_path_spectral_radiance_source_jacobian[ip],
          ray_path_propagation_matrix[ip],
          ray_path_source_vector_nonlte[ip],
          ray_path_propagation_matrix_jacobian[ip],
          ray_path_source_vector_nonlte_jacobian[ip],
          freq_grid_path[ip],
          ray_path_atm_point[ip].temperature,
          it);
    } catch (const std::runtime_error &e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH

void ray_path_atm_pointFromPath(ArrayOfAtmPoint &ray_path_atm_point,
                                const ArrayOfPropagationPathPoint &ray_path,
                                const AtmField &atm_field) try {
  ARTS_TIME_REPORT

  forward_atm_path(
      atm_path_resize(ray_path_atm_point, ray_path), ray_path, atm_field);
}
ARTS_METHOD_ERROR_CATCH

void freq_grid_pathFromPath(ArrayOfAscendingGrid &freq_grid_path,
                            ArrayOfVector3 &freq_wind_shift_jac_path,
                            const AscendingGrid &freq_grid,
                            const ArrayOfPropagationPathPoint &ray_path,
                            const ArrayOfAtmPoint &ray_path_atm_point) try {
  ARTS_TIME_REPORT

  std::string error;

  freq_grid_path.resize(ray_path.size());
  freq_wind_shift_jac_path.resize(ray_path.size());

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size ip = 0; ip < ray_path.size(); ip++) {
    try {
      freq_gridWindShift(freq_grid_path[ip] = freq_grid,
                         freq_wind_shift_jac_path[ip],
                         ray_path_atm_point[ip],
                         ray_path[ip]);
    } catch (std::exception &e) {
#pragma omp critical
      error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH

void ray_path_transmission_matrix_cumulativeFromPath(
    ArrayOfMuelmatVector &ray_path_transmission_matrix_cumulative,
    const ArrayOfMuelmatVector &ray_path_transmission_matrix) try {
  ARTS_TIME_REPORT

  forward_cumulative_transmission(ray_path_transmission_matrix_cumulative,
                                  ray_path_transmission_matrix);
}
ARTS_METHOD_ERROR_CATCH
