#include <arts_omp.h>
#include <atm_path.h>
#include <jacobian.h>
#include <path_point.h>
#include <surf.h>
#include <workspace.h>

#include <ranges>
#include <stdexcept>

#include "arts_conversions.h"
#include "debug.h"
#include "rtepack.h"
#include "rtepack_mueller_matrix.h"

void ray_path_spectral_radianceCalcTransmission(
    ArrayOfStokvecVector &ray_path_spectral_radiance,
    ArrayOfStokvecMatrix &ray_path_spectral_radiance_jacobian,
    const StokvecVector &spectral_radiance_background,
    const ArrayOfMuelmatVector &ray_path_transmission_matrix,
    const ArrayOfMuelmatVector &ray_path_transmission_matrix_cumulative,
    const ArrayOfArrayOfMuelmatMatrix
        &ray_path_transmission_matrix_jacobian) try {
  ARTS_TIME_REPORT

  const Size np = ray_path_transmission_matrix.size();
  const Index nq =
      ray_path_transmission_matrix_jacobian.front().front().nrows();
  const Index nf =
      ray_path_transmission_matrix_jacobian.front().front().ncols();

  ARTS_USER_ERROR_IF(np != ray_path_transmission_matrix.size(),
                     "ray_path_transmission_matrix must have (np) elements")
  ARTS_USER_ERROR_IF(
      np != ray_path_transmission_matrix_cumulative.size(),
      "ray_path_transmission_matrix_cumulative must have (np) elements")
  ARTS_USER_ERROR_IF(
      2 != ray_path_transmission_matrix_jacobian.size() or
          ray_path_transmission_matrix_jacobian.front().size() !=
              ray_path_transmission_matrix_jacobian.back().size() or
          ray_path_transmission_matrix_jacobian.front().size() != np,
      "ray_path_transmission_matrix_jacobian must (2 x np) elements")

  ARTS_USER_ERROR_IF(
      static_cast<Index>(spectral_radiance_background.size()) != nf,
      "spectral_radiance_background must have (nf) elements. Should have ({}"
      ") vs have ({})",
      nf,
      spectral_radiance_background.size())

  if (np == 0) {
    ray_path_spectral_radiance.resize(0);
    ray_path_spectral_radiance_jacobian.resize(0);
    return;
  }

  const auto test_nf = [nf](auto &v) {
    return static_cast<Index>(v.size()) != nf;
  };
  ARTS_USER_ERROR_IF(
      std::any_of(ray_path_transmission_matrix.begin(),
                  ray_path_transmission_matrix.end(),
                  test_nf),
      "ray_path_transmission_matrix must have (nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ray_path_transmission_matrix_cumulative.begin(),
                  ray_path_transmission_matrix_cumulative.end(),
                  test_nf),
      "ray_path_spectral_radiance_source must have (nf) inner elements")

  const auto test_nfnq = [nf, nq](auto &v) {
    return v.ncols() != nf or v.nrows() != nq;
  };
  ARTS_USER_ERROR_IF(
      std::any_of(ray_path_transmission_matrix_jacobian.front().begin(),
                  ray_path_transmission_matrix_jacobian.front().end(),
                  test_nfnq),
      "ray_path_transmission_matrix_jacobian must have (nq x nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ray_path_transmission_matrix_jacobian.back().begin(),
                  ray_path_transmission_matrix_jacobian.back().end(),
                  test_nfnq),
      "ray_path_transmission_matrix_jacobian must have (nq x nf) inner elements")

  ray_path_spectral_radiance.resize(np);
  ray_path_spectral_radiance.back() = spectral_radiance_background;
  ray_path_spectral_radiance_jacobian.resize(np);
  for (auto &t : ray_path_spectral_radiance_jacobian) {
    t.resize(nq, nf);
    t = 0;
  }
  for (Index ip = np - 2; ip >= 0; ip--) {
    ray_path_spectral_radiance[ip] = ray_path_spectral_radiance[ip + 1];
    two_level_linear_transmission_step(
        ray_path_spectral_radiance[ip],
        ray_path_spectral_radiance_jacobian[ip],
        ray_path_spectral_radiance_jacobian[ip + 1],
        ray_path_transmission_matrix[ip + 1],
        ray_path_transmission_matrix_cumulative[ip],
        ray_path_transmission_matrix_jacobian[0][ip + 1],
        ray_path_transmission_matrix_jacobian[1][ip + 1]);
  }
}
ARTS_METHOD_ERROR_CATCH

void ray_path_zeeman_magnetic_fieldFromPath(
    ArrayOfVector3 &ray_path_zeeman_magnetic_field,
    const ArrayOfPropagationPathPoint &ray_path,
    const ArrayOfAtmPoint &ray_path_atmospheric_point) try {
  ARTS_TIME_REPORT

  const Size np = ray_path_atmospheric_point.size();
  ARTS_USER_ERROR_IF(
      np != ray_path.size(),
      "ray_path and ray_path_atmospheric_point must have the same size")

  ray_path_zeeman_magnetic_field.resize(np);

  std::transform(ray_path.begin(),
                 ray_path.end(),
                 ray_path_atmospheric_point.begin(),
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
    const ArrayOfAscendingGrid &ray_path_frequency_grid,
    const ArrayOfAtmPoint &ray_path_atmospheric_point,
    const JacobianTargets &jacobian_targets) try {
  ARTS_TIME_REPORT

  String error{};

  const Index np = ray_path_atmospheric_point.size();
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
          ray_path_frequency_grid[ip],
          ray_path_atmospheric_point[ip].temperature,
          it);
    } catch (const std::runtime_error &e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH

void ray_path_atmospheric_pointFromPath(
    ArrayOfAtmPoint &ray_path_atmospheric_point,
    const ArrayOfPropagationPathPoint &ray_path,
    const AtmField &atm_field) try {
  ARTS_TIME_REPORT

  forward_atm_path(atm_path_resize(ray_path_atmospheric_point, ray_path),
                   ray_path,
                   atm_field);
}
ARTS_METHOD_ERROR_CATCH

void ray_path_frequency_gridFromPath(
    ArrayOfAscendingGrid &ray_path_frequency_grid,
    ArrayOfVector3 &ray_path_frequency_grid_wind_shift_jacobian,
    const AscendingGrid &frequency_grid,
    const ArrayOfPropagationPathPoint &ray_path,
    const ArrayOfAtmPoint &ray_path_atmospheric_point) try {
  ARTS_TIME_REPORT

  std::string error;

  ray_path_frequency_grid.resize(ray_path.size());
  ray_path_frequency_grid_wind_shift_jacobian.resize(ray_path.size());

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size ip = 0; ip < ray_path.size(); ip++) {
    try {
      frequency_gridWindShift(ray_path_frequency_grid[ip] = frequency_grid,
                              ray_path_frequency_grid_wind_shift_jacobian[ip],
                              ray_path_atmospheric_point[ip],
                              ray_path[ip]);
    } catch (std::exception &e) {
#pragma omp critical
      error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH

void ray_path_spectral_radianceCalcClearsky(
    ArrayOfStokvecVector &ray_path_spectral_radiance,
    ArrayOfStokvecMatrix &ray_path_spectral_radiance_jacobian,
    const StokvecVector &background_rad,
    const ArrayOfStokvecVector &ray_path_spectral_radiance_source,
    const ArrayOfStokvecMatrix &ray_path_spectral_radiance_source_jacobian,
    const ArrayOfMuelmatVector &ray_path_transmission_matrix,
    const ArrayOfMuelmatVector &ray_path_transmission_matrix_cumulative,
    const ArrayOfArrayOfMuelmatMatrix
        &ray_path_transmission_matrix_jacobian) try {
  ARTS_TIME_REPORT

  const Size np = ray_path_spectral_radiance_source.size();

  ARTS_USER_ERROR_IF(
      np != ray_path_spectral_radiance_source_jacobian.size(),
      "ray_path_spectral_radiance_source_jacobian must have (np) elements")
  ARTS_USER_ERROR_IF(np != ray_path_transmission_matrix.size(),
                     "ray_path_transmission_matrix must have (np) elements")
  ARTS_USER_ERROR_IF(
      np != ray_path_transmission_matrix_cumulative.size(),
      "ray_path_transmission_matrix_cumulative must have (np) elements")
  ARTS_USER_ERROR_IF(
      2 != ray_path_transmission_matrix_jacobian.size() or
          ray_path_transmission_matrix_jacobian.front().size() !=
              ray_path_transmission_matrix_jacobian.back().size() or
          ray_path_transmission_matrix_jacobian.front().size() != np,
      "ray_path_transmission_matrix_jacobian must (2 x np) elements")

  if (np == 0) {
    ray_path_spectral_radiance.resize(0);
    ray_path_spectral_radiance_jacobian.resize(0);
    return;
  }

  const Index nq = ray_path_spectral_radiance_source_jacobian.front().nrows();
  const Index nf = ray_path_spectral_radiance_source_jacobian.front().ncols();

  const auto test_nf = [nf](auto &v) {
    return static_cast<Index>(v.size()) != nf;
  };
  ARTS_USER_ERROR_IF(static_cast<Size>(nf) != background_rad.size(),
                     "background_rad must have nf elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ray_path_spectral_radiance_source.begin(),
                  ray_path_spectral_radiance_source.end(),
                  test_nf),
      "ray_path_spectral_radiance_source must have (nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ray_path_transmission_matrix.begin(),
                  ray_path_transmission_matrix.end(),
                  test_nf),
      "ray_path_transmission_matrix must have (nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ray_path_transmission_matrix_cumulative.begin(),
                  ray_path_transmission_matrix_cumulative.end(),
                  test_nf),
      "ray_path_spectral_radiance_source must have (nf) inner elements")

  const auto test_nfnq = [nf, nq](auto &v) {
    return v.ncols() != nf or v.nrows() != nq;
  };
  ARTS_USER_ERROR_IF(
      std::any_of(ray_path_spectral_radiance_source_jacobian.begin(),
                  ray_path_spectral_radiance_source_jacobian.end(),
                  test_nfnq),
      "ray_path_spectral_radiance_source_jacobian must have (nq x nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ray_path_transmission_matrix_jacobian.front().begin(),
                  ray_path_transmission_matrix_jacobian.front().end(),
                  test_nfnq),
      "ray_path_transmission_matrix_jacobian must have (nq x nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ray_path_transmission_matrix_jacobian.back().begin(),
                  ray_path_transmission_matrix_jacobian.back().end(),
                  test_nfnq),
      "ray_path_transmission_matrix_jacobian must have (nq x nf) inner elements")

  ray_path_spectral_radiance.resize(np);
  ray_path_spectral_radiance.back() = background_rad;
  ray_path_spectral_radiance_jacobian.resize(np);
  for (auto &t : ray_path_spectral_radiance_jacobian) {
    t.resize(nq, nf);
    t = 0;
  }
  for (Index ip = np - 2; ip >= 0; ip--) {
    ray_path_spectral_radiance[ip] = ray_path_spectral_radiance[ip + 1];
    two_level_linear_emission_step(
        ray_path_spectral_radiance[ip],
        ray_path_spectral_radiance_jacobian[ip],
        ray_path_spectral_radiance_jacobian[ip + 1],
        ray_path_spectral_radiance_source[ip],
        ray_path_spectral_radiance_source[ip + 1],
        ray_path_spectral_radiance_source_jacobian[ip],
        ray_path_spectral_radiance_source_jacobian[ip + 1],
        ray_path_transmission_matrix[ip + 1],
        ray_path_transmission_matrix_cumulative[ip],
        ray_path_transmission_matrix_jacobian[0][ip + 1],
        ray_path_transmission_matrix_jacobian[1][ip + 1]);
  }
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
