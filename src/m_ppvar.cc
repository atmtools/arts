#include <arts_omp.h>
#include <atm_path.h>
#include <jacobian.h>
#include <path_point.h>
#include <rte.h>
#include <surf.h>
#include <workspace.h>
#include "sorted_grid.h"

void propagation_path_spectral_radianceCalcTransmission(
    ArrayOfStokvecVector &propagation_path_spectral_radiance,
    ArrayOfStokvecMatrix &propagation_path_spectral_radiance_jacobian,
    const ArrayOfMuelmatVector &propagation_path_transmission_matrix,
    const ArrayOfMuelmatVector &propagation_path_transmission_matrix_cumulative,
    const ArrayOfArrayOfMuelmatMatrix
        &propagation_path_transmission_matrix_jacobian) try {
  const Size np = propagation_path_transmission_matrix.size();

  ARTS_USER_ERROR_IF(
      np not_eq propagation_path_transmission_matrix.size(),
      "propagation_path_transmission_matrix must have (np) elements")
  ARTS_USER_ERROR_IF(
      np not_eq propagation_path_transmission_matrix_cumulative.size(),
      "propagation_path_transmission_matrix_cumulative must have (np) elements")
  ARTS_USER_ERROR_IF(
      2 not_eq propagation_path_transmission_matrix_jacobian.size() or
          propagation_path_transmission_matrix_jacobian.front().size() not_eq
              propagation_path_transmission_matrix_jacobian.back().size() or
          propagation_path_transmission_matrix_jacobian.front().size() not_eq
              np,
      "propagation_path_transmission_matrix_jacobian must (2 x np) elements")

  if (np == 0) {
    propagation_path_spectral_radiance.resize(0);
    propagation_path_spectral_radiance_jacobian.resize(0);
    return;
  }

  const Index nq =
      propagation_path_transmission_matrix_jacobian.front().front().nrows();
  const Index nf =
      propagation_path_transmission_matrix_jacobian.front().front().ncols();

  const auto test_nf = [nf](auto &v) { return v.size() not_eq nf; };
  ARTS_USER_ERROR_IF(
      std::any_of(propagation_path_transmission_matrix.begin(),
                  propagation_path_transmission_matrix.end(),
                  test_nf),
      "propagation_path_transmission_matrix must have (nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(propagation_path_transmission_matrix_cumulative.begin(),
                  propagation_path_transmission_matrix_cumulative.end(),
                  test_nf),
      "propagation_path_spectral_radiance_source must have (nf) inner elements")

  const auto test_nfnq = [nf, nq](auto &v) {
    return v.ncols() not_eq nf or v.nrows() not_eq nq;
  };
  ARTS_USER_ERROR_IF(
      std::any_of(propagation_path_transmission_matrix_jacobian.front().begin(),
                  propagation_path_transmission_matrix_jacobian.front().end(),
                  test_nfnq),
      "propagation_path_transmission_matrix_jacobian must have (nq x nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(propagation_path_transmission_matrix_jacobian.back().begin(),
                  propagation_path_transmission_matrix_jacobian.back().end(),
                  test_nfnq),
      "propagation_path_transmission_matrix_jacobian must have (nq x nf) inner elements")

  propagation_path_spectral_radiance.resize(
      np, StokvecVector(nf, Stokvec{1, 0, 0, 0}));
  propagation_path_spectral_radiance_jacobian.resize(np, StokvecMatrix(nq, nf));
  for (Index ip = np - 2; ip >= 0; ip--) {
    propagation_path_spectral_radiance[ip] =
        propagation_path_spectral_radiance[ip + 1];
    two_level_linear_transmission_step(
        propagation_path_spectral_radiance[ip],
        propagation_path_spectral_radiance_jacobian[ip],
        propagation_path_spectral_radiance_jacobian[ip + 1],
        propagation_path_transmission_matrix[ip + 1],
        propagation_path_transmission_matrix_cumulative[ip],
        propagation_path_transmission_matrix_jacobian[0][ip + 1],
        propagation_path_transmission_matrix_jacobian[1][ip + 1]);
  }
}
ARTS_METHOD_ERROR_CATCH

void propagation_path_propagation_matrixCalc(
    const Workspace &ws,
    ArrayOfPropmatVector &propagation_path_propagation_matrix,
    ArrayOfStokvecVector &propagation_path_source_vector_nonlte,
    ArrayOfPropmatMatrix &propagation_path_propagation_matrix_jacobian,
    ArrayOfStokvecMatrix &propagation_path_source_vector_nonlte_jacobian,
    const Agenda &propagation_matrix_agenda,
    const JacobianTargets &jacobian_targets,
    const ArrayOfAscendingGrid &propagation_path_frequency_grid,
    const ArrayOfPropagationPathPoint &rad_path,
    const ArrayOfAtmPoint &propagation_path_atmospheric_point) try {
  ArrayOfString fail_msg;
  bool do_abort = false;

  const Size np = rad_path.size();
  if (np == 0) {
    propagation_path_propagation_matrix.resize(0);
    propagation_path_source_vector_nonlte.resize(0);
    propagation_path_propagation_matrix_jacobian.resize(0);
    propagation_path_source_vector_nonlte_jacobian.resize(0);
    return;
  }

  propagation_path_propagation_matrix.resize(np);
  propagation_path_source_vector_nonlte.resize(np);
  propagation_path_propagation_matrix_jacobian.resize(np);
  propagation_path_source_vector_nonlte_jacobian.resize(np);

  // Loop ppath points and determine radiative properties
#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size ip = 0; ip < np; ip++) {
    if (do_abort) continue;
    try {
      //! FIXME: Send in full path point instead of this
      get_stepwise_clearsky_propmat(
          ws,
          propagation_path_propagation_matrix[ip],
          propagation_path_source_vector_nonlte[ip],
          propagation_path_propagation_matrix_jacobian[ip],
          propagation_path_source_vector_nonlte_jacobian[ip],
          propagation_matrix_agenda,
          jacobian_targets,
          propagation_path_frequency_grid[ip],
          rad_path[ip],
          propagation_path_atmospheric_point[ip]);
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

  ARTS_USER_ERROR_IF(do_abort, "Error messages from failed cases:\n", fail_msg)
}
ARTS_METHOD_ERROR_CATCH

void propagation_path_spectral_radiance_sourceFromPropmat(
    ArrayOfStokvecVector &propagation_path_spectral_radiance_source,
    ArrayOfStokvecMatrix &propagation_path_spectral_radiance_source_jacobian,
    const ArrayOfPropmatVector &propagation_path_propagation_matrix,
    const ArrayOfStokvecVector &propagation_path_source_vector_nonlte,
    const ArrayOfPropmatMatrix &propagation_path_propagation_matrix_jacobian,
    const ArrayOfStokvecMatrix &propagation_path_source_vector_nonlte_jacobian,
    const ArrayOfAscendingGrid &propagation_path_frequency_grid,
    const ArrayOfAtmPoint &propagation_path_atmospheric_point,
    const JacobianTargets &jacobian_targets) try {
  ArrayOfString fail_msg;
  bool do_abort = false;

  const Index np = propagation_path_atmospheric_point.size();
  if (np == 0) {
    propagation_path_spectral_radiance_source.resize(0);
    propagation_path_spectral_radiance_source_jacobian.resize(0);
    return;
  }

  const Index nf = propagation_path_propagation_matrix.front().size();
  const Index nq = jacobian_targets.target_count();

  const Index it =
      jacobian_targets.target_position<Jacobian::AtmTarget>(Atm::Key::t);

  propagation_path_spectral_radiance_source.resize(np, StokvecVector(nf));
  propagation_path_spectral_radiance_source_jacobian.resize(
      np, StokvecMatrix(nq, nf));

  // Loop ppath points and determine radiative properties
#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Index ip = 0; ip < np; ip++) {
    if (do_abort) continue;
    try {
      rtepack::source::level_nlte(
          propagation_path_spectral_radiance_source[ip],
          propagation_path_spectral_radiance_source_jacobian[ip],
          propagation_path_propagation_matrix[ip],
          propagation_path_source_vector_nonlte[ip],
          propagation_path_propagation_matrix_jacobian[ip],
          propagation_path_source_vector_nonlte_jacobian[ip],
          propagation_path_frequency_grid[ip],
          propagation_path_atmospheric_point[ip].temperature,
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

void propagation_path_transmission_matrixCalc(
    ArrayOfMuelmatVector &propagation_path_transmission_matrix,
    ArrayOfArrayOfMuelmatMatrix &propagation_path_transmission_matrix_jacobian,
    Vector &propagation_path_distance,
    ArrayOfArrayOfVector &propagation_path_distance_jacobian,
    const ArrayOfPropmatVector &propagation_path_propagation_matrix,
    const ArrayOfPropmatMatrix &propagation_path_propagation_matrix_jacobian,
    const ArrayOfPropagationPathPoint &rad_path,
    const ArrayOfAtmPoint &propagation_path_atmospheric_point,
    const SurfaceField &surface_field,
    const JacobianTargets &jacobian_targets,
    const Index &hse_derivative) try {
  ArrayOfString fail_msg;
  bool do_abort = false;

  // HSE variables
  const Index temperature_derivative_position =
      jacobian_targets.target_position<Jacobian::AtmTarget>(Atm::Key::t);

  const Size np = rad_path.size();

  if (np == 0) {
    propagation_path_transmission_matrix.resize(0);
    propagation_path_transmission_matrix_jacobian.resize(
        2, ArrayOfMuelmatMatrix{});
    propagation_path_distance.resize(0);
    propagation_path_distance_jacobian.resize(2, ArrayOfVector{});
    return;
  }

  const Index nf = propagation_path_propagation_matrix.front().size();
  const Index nq = jacobian_targets.target_count();

  propagation_path_transmission_matrix.resize(np, MuelmatVector(nf));
  propagation_path_transmission_matrix_jacobian.resize(
      2, ArrayOfMuelmatMatrix(np, MuelmatMatrix(nq, nf)));
  propagation_path_distance.resize(np);
  propagation_path_distance_jacobian.resize(2,
                                            ArrayOfVector(np, Vector(nq, 0)));

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size ip = 1; ip < np; ip++) {
    if (do_abort) continue;
    try {
      propagation_path_distance[ip - 1] =
          path::distance(rad_path[ip - 1].pos, rad_path[ip].pos, surface_field);
      if (hse_derivative and temperature_derivative_position >= 0) {
        propagation_path_distance_jacobian
            [0][ip][temperature_derivative_position] =
                propagation_path_distance[ip - 1] /
                (2.0 * propagation_path_atmospheric_point[ip - 1].temperature);
        propagation_path_distance_jacobian
            [1][ip][temperature_derivative_position] =
                propagation_path_distance[ip - 1] /
                (2.0 * propagation_path_atmospheric_point[ip].temperature);
      }

      two_level_exp(propagation_path_transmission_matrix[ip],
                    propagation_path_transmission_matrix_jacobian[0][ip],
                    propagation_path_transmission_matrix_jacobian[1][ip],
                    propagation_path_propagation_matrix[ip - 1],
                    propagation_path_propagation_matrix[ip],
                    propagation_path_propagation_matrix_jacobian[ip - 1],
                    propagation_path_propagation_matrix_jacobian[ip],
                    propagation_path_distance[ip - 1],
                    propagation_path_distance_jacobian[0][ip],
                    propagation_path_distance_jacobian[1][ip]);
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

  ARTS_USER_ERROR_IF(do_abort, "Error messages from failed cases:\n", fail_msg)
}
ARTS_METHOD_ERROR_CATCH

void propagation_path_atmospheric_pointFromPath(
    ArrayOfAtmPoint &propagation_path_atmospheric_point,
    const ArrayOfPropagationPathPoint &rad_path,
    const AtmField &atm_field) try {
  forward_atm_path(
      atm_path_resize(propagation_path_atmospheric_point, rad_path),
      rad_path,
      atm_field);
}
ARTS_METHOD_ERROR_CATCH

void propagation_path_frequency_gridFromPath(
    ArrayOfAscendingGrid &propagation_path_frequency_grid,
    const AscendingGrid &frequency_grid,
    const ArrayOfPropagationPathPoint &rad_path,
    const ArrayOfAtmPoint &propagation_path_atmospheric_point,
    const Numeric &rte_alonglos_v) try {
  forward_path_freq(path_freq_resize(propagation_path_frequency_grid,
                                     frequency_grid,
                                     propagation_path_atmospheric_point),
                    frequency_grid,
                    rad_path,
                    propagation_path_atmospheric_point,
                    rte_alonglos_v);
}
ARTS_METHOD_ERROR_CATCH

void propagation_path_spectral_radianceCalcEmission(
    ArrayOfStokvecVector &propagation_path_spectral_radiance,
    ArrayOfStokvecMatrix &propagation_path_spectral_radiance_jacobian,
    const StokvecVector &background_rad,
    const ArrayOfStokvecVector &propagation_path_spectral_radiance_source,
    const ArrayOfStokvecMatrix
        &propagation_path_spectral_radiance_source_jacobian,
    const ArrayOfMuelmatVector &propagation_path_transmission_matrix,
    const ArrayOfMuelmatVector &propagation_path_transmission_matrix_cumulative,
    const ArrayOfArrayOfMuelmatMatrix
        &propagation_path_transmission_matrix_jacobian) try {
  const Size np = propagation_path_spectral_radiance_source.size();

  ARTS_USER_ERROR_IF(
      np not_eq propagation_path_spectral_radiance_source_jacobian.size(),
      "propagation_path_spectral_radiance_source_jacobian must have (np) elements")
  ARTS_USER_ERROR_IF(
      np not_eq propagation_path_transmission_matrix.size(),
      "propagation_path_transmission_matrix must have (np) elements")
  ARTS_USER_ERROR_IF(
      np not_eq propagation_path_transmission_matrix_cumulative.size(),
      "propagation_path_transmission_matrix_cumulative must have (np) elements")
  ARTS_USER_ERROR_IF(
      2 not_eq propagation_path_transmission_matrix_jacobian.size() or
          propagation_path_transmission_matrix_jacobian.front().size() not_eq
              propagation_path_transmission_matrix_jacobian.back().size() or
          propagation_path_transmission_matrix_jacobian.front().size() not_eq
              np,
      "propagation_path_transmission_matrix_jacobian must (2 x np) elements")

  if (np == 0) {
    propagation_path_spectral_radiance.resize(0);
    propagation_path_spectral_radiance_jacobian.resize(0);
    return;
  }

  const Index nq =
      propagation_path_spectral_radiance_source_jacobian.front().nrows();
  const Index nf =
      propagation_path_spectral_radiance_source_jacobian.front().ncols();

  const auto test_nf = [nf](auto &v) { return v.size() not_eq nf; };
  ARTS_USER_ERROR_IF(nf not_eq background_rad.size(),
                     "background_rad must have nf elements")
  ARTS_USER_ERROR_IF(
      std::any_of(propagation_path_spectral_radiance_source.begin(),
                  propagation_path_spectral_radiance_source.end(),
                  test_nf),
      "propagation_path_spectral_radiance_source must have (nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(propagation_path_transmission_matrix.begin(),
                  propagation_path_transmission_matrix.end(),
                  test_nf),
      "propagation_path_transmission_matrix must have (nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(propagation_path_transmission_matrix_cumulative.begin(),
                  propagation_path_transmission_matrix_cumulative.end(),
                  test_nf),
      "propagation_path_spectral_radiance_source must have (nf) inner elements")

  const auto test_nfnq = [nf, nq](auto &v) {
    return v.ncols() not_eq nf or v.nrows() not_eq nq;
  };
  ARTS_USER_ERROR_IF(
      std::any_of(propagation_path_spectral_radiance_source_jacobian.begin(),
                  propagation_path_spectral_radiance_source_jacobian.end(),
                  test_nfnq),
      "propagation_path_spectral_radiance_source_jacobian must have (nq x nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(propagation_path_transmission_matrix_jacobian.front().begin(),
                  propagation_path_transmission_matrix_jacobian.front().end(),
                  test_nfnq),
      "propagation_path_transmission_matrix_jacobian must have (nq x nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(propagation_path_transmission_matrix_jacobian.back().begin(),
                  propagation_path_transmission_matrix_jacobian.back().end(),
                  test_nfnq),
      "propagation_path_transmission_matrix_jacobian must have (nq x nf) inner elements")

  propagation_path_spectral_radiance.resize(np, background_rad);
  propagation_path_spectral_radiance_jacobian.resize(np, StokvecMatrix(nq, nf));
  for (Index ip = np - 2; ip >= 0; ip--) {
    propagation_path_spectral_radiance[ip] =
        propagation_path_spectral_radiance[ip + 1];
    two_level_linear_emission_step(
        propagation_path_spectral_radiance[ip],
        propagation_path_spectral_radiance_jacobian[ip],
        propagation_path_spectral_radiance_jacobian[ip + 1],
        propagation_path_spectral_radiance_source[ip],
        propagation_path_spectral_radiance_source[ip + 1],
        propagation_path_spectral_radiance_source_jacobian[ip],
        propagation_path_spectral_radiance_source_jacobian[ip + 1],
        propagation_path_transmission_matrix[ip + 1],
        propagation_path_transmission_matrix_cumulative[ip],
        propagation_path_transmission_matrix_jacobian[0][ip + 1],
        propagation_path_transmission_matrix_jacobian[1][ip + 1]);
  }
}
ARTS_METHOD_ERROR_CATCH

void propagation_path_transmission_matrix_cumulativeForward(
    ArrayOfMuelmatVector &propagation_path_transmission_matrix_cumulative,
    const ArrayOfMuelmatVector &propagation_path_transmission_matrix) {
  propagation_path_transmission_matrix_cumulative =
      forward_cumulative_transmission(propagation_path_transmission_matrix);
}

void propagation_path_transmission_matrix_cumulativeReverse(
    ArrayOfMuelmatVector &propagation_path_transmission_matrix_cumulative,
    const ArrayOfMuelmatVector &propagation_path_transmission_matrix) {
  propagation_path_transmission_matrix_cumulative =
      reverse_cumulative_transmission(propagation_path_transmission_matrix);
}
