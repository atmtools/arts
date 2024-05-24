#include <arts_omp.h>
#include <atm_path.h>
#include <jacobian.h>
#include <path_point.h>
#include <rte.h>
#include <surf.h>
#include <workspace.h>

#include "debug.h"
#include "sorted_grid.h"

void ray_path_spectral_radianceCalcTransmission(
    ArrayOfStokvecVector &ray_path_spectral_radiance,
    ArrayOfStokvecMatrix &ray_path_spectral_radiance_jacobian,
    const StokvecVector &spectral_radiance_background,
    const ArrayOfMuelmatVector &ray_path_transmission_matrix,
    const ArrayOfMuelmatVector &ray_path_transmission_matrix_cumulative,
    const ArrayOfArrayOfMuelmatMatrix
        &ray_path_transmission_matrix_jacobian) try {
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
      spectral_radiance_background.size() != nf,
      "spectral_radiance_background must have (nf) elements. Should have (",
      nf,
      ") vs have (",
      spectral_radiance_background.size(),
      ")")

  if (np == 0) {
    ray_path_spectral_radiance.resize(0);
    ray_path_spectral_radiance_jacobian.resize(0);
    return;
  }

  const auto test_nf = [nf](auto &v) { return v.size() != nf; };
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

void ray_path_propagation_matrixFromPath(
    const Workspace &ws,
    ArrayOfPropmatVector &ray_path_propagation_matrix,
    ArrayOfStokvecVector &ray_path_source_vector_nonlte,
    ArrayOfPropmatMatrix &ray_path_propagation_matrix_jacobian,
    ArrayOfStokvecMatrix &ray_path_source_vector_nonlte_jacobian,
    const Agenda &propagation_matrix_agenda,
    const JacobianTargets &jacobian_targets,
    const ArrayOfAscendingGrid &ray_path_frequency_grid,
    const ArrayOfPropagationPathPoint &ray_path,
    const ArrayOfAtmPoint &ray_path_atmospheric_point) try {
  const Size np = ray_path.size();
  if (np == 0) {
    ray_path_propagation_matrix.resize(0);
    ray_path_source_vector_nonlte.resize(0);
    ray_path_propagation_matrix_jacobian.resize(0);
    ray_path_source_vector_nonlte_jacobian.resize(0);
    return;
  }

  ray_path_propagation_matrix.resize(np);
  ray_path_source_vector_nonlte.resize(np);
  ray_path_propagation_matrix_jacobian.resize(np);
  ray_path_source_vector_nonlte_jacobian.resize(np);

  if (arts_omp_in_parallel()) {
    for (Size ip = 0; ip < np; ip++) {
      get_stepwise_clearsky_propmat(ws,
                                    ray_path_propagation_matrix[ip],
                                    ray_path_source_vector_nonlte[ip],
                                    ray_path_propagation_matrix_jacobian[ip],
                                    ray_path_source_vector_nonlte_jacobian[ip],
                                    propagation_matrix_agenda,
                                    jacobian_targets,
                                    ray_path_frequency_grid[ip],
                                    ray_path[ip],
                                    ray_path_atmospheric_point[ip]);
    }
  } else {
    ArrayOfString fail_msg;
    bool do_abort = false;
#pragma omp parallel for if (!arts_omp_in_parallel())
    for (Size ip = 0; ip < np; ip++) {
      if (do_abort) continue;
      try {
        get_stepwise_clearsky_propmat(
            ws,
            ray_path_propagation_matrix[ip],
            ray_path_source_vector_nonlte[ip],
            ray_path_propagation_matrix_jacobian[ip],
            ray_path_source_vector_nonlte_jacobian[ip],
            propagation_matrix_agenda,
            jacobian_targets,
            ray_path_frequency_grid[ip],
            ray_path[ip],
            ray_path_atmospheric_point[ip]);
      } catch (const std::runtime_error &e) {
#pragma omp critical(iyEmissionStandard_source)
        {
          do_abort = true;
          fail_msg.push_back(
              var_string("Runtime-error in propagation radiative "
                         "properties calculation at index ",
                         ip,
                         ": \n",
                         e.what()));
        }
      }
    }

    ARTS_USER_ERROR_IF(
        do_abort, "Error messages from failed cases:\n", fail_msg)
  }
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
  ArrayOfString fail_msg;
  bool do_abort = false;

  const Index np = ray_path_atmospheric_point.size();
  if (np == 0) {
    ray_path_spectral_radiance_source.resize(0);
    ray_path_spectral_radiance_source_jacobian.resize(0);
    return;
  }

  const Index nf = ray_path_propagation_matrix.front().size();
  const Index nq = jacobian_targets.target_count();

  const Index it =
      jacobian_targets.target_position<Jacobian::AtmTarget>(AtmKey::t);

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
    if (do_abort) continue;
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
#pragma omp critical(iyEmissionStandard_source)
      {
        do_abort = true;
        fail_msg.push_back(
            var_string("Runtime-error in source calculation at index ",
                       ip,
                       ": \n",
                       e.what()));
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void ray_path_transmission_matrixFromPath(
    ArrayOfMuelmatVector &ray_path_transmission_matrix,
    ArrayOfArrayOfMuelmatMatrix &ray_path_transmission_matrix_jacobian,
    const ArrayOfPropmatVector &ray_path_propagation_matrix,
    const ArrayOfPropmatMatrix &ray_path_propagation_matrix_jacobian,
    const ArrayOfPropagationPathPoint &ray_path,
    const ArrayOfAtmPoint &ray_path_atmospheric_point,
    const SurfaceField &surface_field,
    const JacobianTargets &jacobian_targets,
    const Index &hse_derivative) try {
  // HSE variables
  const Index temperature_derivative_position =
      jacobian_targets.target_position<Jacobian::AtmTarget>(AtmKey::t);

  const Size np = ray_path.size();

  if (np == 0) {
    ray_path_transmission_matrix.resize(0);
    ray_path_transmission_matrix_jacobian.resize(2);
    ray_path_transmission_matrix_jacobian[0].resize(0);
    ray_path_transmission_matrix_jacobian[1].resize(0);
    return;
  }

  const Index nf = ray_path_propagation_matrix.front().size();
  const Index nq = jacobian_targets.target_count();

  Vector ray_path_distance_jacobian1(nq, 0.0);
  Vector ray_path_distance_jacobian2(nq, 0.0);

  ray_path_transmission_matrix.resize(np);
  for (auto &t : ray_path_transmission_matrix) {
    t.resize(nf);
    t = 0.0;
  }
  ray_path_transmission_matrix_jacobian.resize(2);
  for (auto &t : ray_path_transmission_matrix_jacobian) {
    t.resize(np);
    for (auto &tt : t) {
      tt.resize(nq, nf);
      tt = 0.0;
    }
  }

  if (arts_omp_in_parallel()) {
    for (Size ip = 1; ip < np; ip++) {
      const Numeric ray_path_distance = path::distance(
          ray_path[ip - 1].pos, ray_path[ip].pos, surface_field.ellipsoid);
      if (hse_derivative and temperature_derivative_position >= 0) {
        ray_path_distance_jacobian1[temperature_derivative_position] =
            ray_path_distance /
            (2.0 * ray_path_atmospheric_point[ip - 1].temperature);
        ray_path_distance_jacobian2[temperature_derivative_position] =
            ray_path_distance /
            (2.0 * ray_path_atmospheric_point[ip].temperature);
      }

      two_level_exp(ray_path_transmission_matrix[ip],
                    ray_path_transmission_matrix_jacobian[0][ip],
                    ray_path_transmission_matrix_jacobian[1][ip],
                    ray_path_propagation_matrix[ip - 1],
                    ray_path_propagation_matrix[ip],
                    ray_path_propagation_matrix_jacobian[ip - 1],
                    ray_path_propagation_matrix_jacobian[ip],
                    ray_path_distance,
                    ray_path_distance_jacobian1,
                    ray_path_distance_jacobian2);
    }
  } else {
    ArrayOfString fail_msg;
    bool do_abort = false;

#pragma omp parallel for if (!arts_omp_in_parallel()) \
    firstprivate(ray_path_distance_jacobian1, ray_path_distance_jacobian2)
    for (Size ip = 1; ip < np; ip++) {
      if (do_abort) continue;
      try {
        const Numeric ray_path_distance = path::distance(
            ray_path[ip - 1].pos, ray_path[ip].pos, surface_field.ellipsoid);
        if (hse_derivative and temperature_derivative_position >= 0) {
          ray_path_distance_jacobian1[temperature_derivative_position] =
              ray_path_distance /
              (2.0 * ray_path_atmospheric_point[ip - 1].temperature);
          ray_path_distance_jacobian2[temperature_derivative_position] =
              ray_path_distance /
              (2.0 * ray_path_atmospheric_point[ip].temperature);
        }

        two_level_exp(ray_path_transmission_matrix[ip],
                      ray_path_transmission_matrix_jacobian[0][ip],
                      ray_path_transmission_matrix_jacobian[1][ip],
                      ray_path_propagation_matrix[ip - 1],
                      ray_path_propagation_matrix[ip],
                      ray_path_propagation_matrix_jacobian[ip - 1],
                      ray_path_propagation_matrix_jacobian[ip],
                      ray_path_distance,
                      ray_path_distance_jacobian1,
                      ray_path_distance_jacobian2);
      } catch (const std::runtime_error &e) {
#pragma omp critical(iyEmissionStandard_transmission)
        {
          do_abort = true;
          fail_msg.push_back(
              var_string("Runtime-error in transmission calculation at index ",
                         ip,
                         ": \n",
                         e.what()));
        }
      }
    }

    ARTS_USER_ERROR_IF(
        do_abort, "Error messages from failed cases:\n", fail_msg)
  }
}
ARTS_METHOD_ERROR_CATCH

void ray_path_atmospheric_pointFromPath(
    ArrayOfAtmPoint &ray_path_atmospheric_point,
    const ArrayOfPropagationPathPoint &ray_path,
    const AtmField &atm_field) try {
  forward_atm_path(atm_path_resize(ray_path_atmospheric_point, ray_path),
                   ray_path,
                   atm_field);
}
ARTS_METHOD_ERROR_CATCH

void ray_path_frequency_gridFromPath(
    ArrayOfAscendingGrid &ray_path_frequency_grid,
    const AscendingGrid &frequency_grid,
    const ArrayOfPropagationPathPoint &ray_path,
    const ArrayOfAtmPoint &ray_path_atmospheric_point,
    const Numeric &rte_alonglos_v) try {
  forward_path_freq(
      path_freq_resize(
          ray_path_frequency_grid, frequency_grid, ray_path_atmospheric_point),
      frequency_grid,
      ray_path,
      ray_path_atmospheric_point,
      rte_alonglos_v);
}
ARTS_METHOD_ERROR_CATCH

void ray_path_spectral_radianceCalcEmission(
    ArrayOfStokvecVector &ray_path_spectral_radiance,
    ArrayOfStokvecMatrix &ray_path_spectral_radiance_jacobian,
    const StokvecVector &background_rad,
    const ArrayOfStokvecVector &ray_path_spectral_radiance_source,
    const ArrayOfStokvecMatrix &ray_path_spectral_radiance_source_jacobian,
    const ArrayOfMuelmatVector &ray_path_transmission_matrix,
    const ArrayOfMuelmatVector &ray_path_transmission_matrix_cumulative,
    const ArrayOfArrayOfMuelmatMatrix
        &ray_path_transmission_matrix_jacobian) try {
  const Size np = ray_path_spectral_radiance_source.size();

  ARTS_USER_ERROR_IF(np != ray_path_spectral_radiance_source_jacobian.size(),
                     "ray_path_spectral_radiance_source_jacobian must have (",
                     np,
                     ") elements, has ",
                     ray_path_spectral_radiance_source_jacobian.size())
  ARTS_USER_ERROR_IF(np != ray_path_transmission_matrix.size(),
                     "ray_path_transmission_matrix must have (",
                     np,
                     ") elements, has ",
                     ray_path_transmission_matrix.size())
  ARTS_USER_ERROR_IF(np != ray_path_transmission_matrix_cumulative.size(),
                     "ray_path_transmission_matrix_cumulative must have (",
                     np,
                     ") elements, has ",
                     ray_path_transmission_matrix_cumulative.size())
  ARTS_USER_ERROR_IF(
      2 != ray_path_transmission_matrix_jacobian.size(),
      "Outer size of ray_path_transmission_matrix_jacobian must be 2, is ",
      ray_path_transmission_matrix_jacobian.size())
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(ray_path_transmission_matrix_jacobian,
                          Cmp::ne(np),
                          &ArrayOfMuelmatMatrix::size),
      "Inner size of ray_path_transmission_matrix_jacobian must be (",
      np,
      "), isn't for some element.  Good luck!")

  if (np == 0) {
    ray_path_spectral_radiance.resize(0);
    ray_path_spectral_radiance_jacobian.resize(0);
    return;
  }

  const Index nq = ray_path_spectral_radiance_source_jacobian.front().nrows();
  const Index nf = ray_path_spectral_radiance_source_jacobian.front().ncols();

  const auto test_nf = [nf](auto &v) { return v.size() != nf; };
  ARTS_USER_ERROR_IF(nf != background_rad.size(),
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

  const auto test_nf = [nf](auto &v) { return v.size() != nf; };
  ARTS_USER_ERROR_IF(nf != background_rad.size(),
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

void ray_path_transmission_matrix_cumulativeForward(
    ArrayOfMuelmatVector &ray_path_transmission_matrix_cumulative,
    const ArrayOfMuelmatVector &ray_path_transmission_matrix) try {
  ray_path_transmission_matrix_cumulative =
      forward_cumulative_transmission(ray_path_transmission_matrix);
}
ARTS_METHOD_ERROR_CATCH

void ray_path_transmission_matrix_cumulativeReverse(
    ArrayOfMuelmatVector &ray_path_transmission_matrix_cumulative,
    const ArrayOfMuelmatVector &ray_path_transmission_matrix) try {
  ray_path_transmission_matrix_cumulative =
      reverse_cumulative_transmission(ray_path_transmission_matrix);
}
ARTS_METHOD_ERROR_CATCH
