#include <arts_omp.h>
#include <atm_path.h>
#include <jacobian.h>
#include <path_point.h>
#include <rte.h>
#include <surf.h>
#include <workspace.h>

void spectral_radiance_pathCalcTransmission(
    ArrayOfStokvecVector &spectral_radiance_path,
    ArrayOfStokvecMatrix &spectral_radiance_path_jacobian,
    const ArrayOfMuelmatVector &ppvar_tramat,
    const ArrayOfMuelmatVector &ppvar_cumtramat,
    const ArrayOfArrayOfMuelmatMatrix &ppvar_dtramat) try {
  const Size np = ppvar_tramat.size();

  ARTS_USER_ERROR_IF(np not_eq ppvar_tramat.size(),
                     "ppvar_tramat must have (np) elements")
  ARTS_USER_ERROR_IF(np not_eq ppvar_cumtramat.size(),
                     "ppvar_cumtramat must have (np) elements")
  ARTS_USER_ERROR_IF(
      2 not_eq ppvar_dtramat.size() or
          ppvar_dtramat.front().size() not_eq ppvar_dtramat.back().size() or
          ppvar_dtramat.front().size() not_eq np,
      "ppvar_dtramat must (2 x np) elements")

  if (np == 0) {
    spectral_radiance_path.resize(0);
    spectral_radiance_path_jacobian.resize(0);
    return;
  }

  const Index nq = ppvar_dtramat.front().front().nrows();
  const Index nf = ppvar_dtramat.front().front().ncols();

  const auto test_nf = [nf](auto &v) { return v.size() not_eq nf; };
  ARTS_USER_ERROR_IF(
      std::any_of(ppvar_tramat.begin(), ppvar_tramat.end(), test_nf),
      "ppvar_tramat must have (nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ppvar_cumtramat.begin(), ppvar_cumtramat.end(), test_nf),
      "spectral_radiance_path_source must have (nf) inner elements")

  const auto test_nfnq = [nf, nq](auto &v) {
    return v.ncols() not_eq nf or v.nrows() not_eq nq;
  };
  ARTS_USER_ERROR_IF(std::any_of(ppvar_dtramat.front().begin(),
                                 ppvar_dtramat.front().end(),
                                 test_nfnq),
                     "ppvar_dtramat must have (nq x nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(
          ppvar_dtramat.back().begin(), ppvar_dtramat.back().end(), test_nfnq),
      "ppvar_dtramat must have (nq x nf) inner elements")

  spectral_radiance_path.resize(np, StokvecVector(nf, Stokvec{1, 0, 0, 0}));
  spectral_radiance_path_jacobian.resize(np, StokvecMatrix(nq, nf));
  for (Index ip = np - 2; ip >= 0; ip--) {
    spectral_radiance_path[ip] = spectral_radiance_path[ip + 1];
    two_level_linear_transmission_step(spectral_radiance_path[ip],
                                       spectral_radiance_path_jacobian[ip],
                                       spectral_radiance_path_jacobian[ip + 1],
                                       ppvar_tramat[ip + 1],
                                       ppvar_cumtramat[ip],
                                       ppvar_dtramat[0][ip + 1],
                                       ppvar_dtramat[1][ip + 1]);
  }
}
ARTS_METHOD_ERROR_CATCH

void ppvar_propmatCalc(const Workspace &ws,
                       ArrayOfPropmatVector &ppvar_propmat,
                       ArrayOfStokvecVector &ppvar_nlte,
                       ArrayOfPropmatMatrix &ppvar_dpropmat,
                       ArrayOfStokvecMatrix &ppvar_dnlte,
                       const Agenda &propmat_clearsky_agenda,
                       const JacobianTargets &jacobian_targets,
                       const ArrayOfVector &ppvar_f,
                       const ArrayOfPropagationPathPoint &rad_path,
                       const ArrayOfAtmPoint &ppvar_atm) try {
  ArrayOfString fail_msg;
  bool do_abort = false;

  const Size np = rad_path.size();
  if (np == 0) {
    ppvar_propmat.resize(0);
    ppvar_nlte.resize(0);
    ppvar_dpropmat.resize(0);
    ppvar_dnlte.resize(0);
    return;
  }

  ppvar_propmat.resize(np);
  ppvar_nlte.resize(np);
  ppvar_dpropmat.resize(np);
  ppvar_dnlte.resize(np);

  // Loop ppath points and determine radiative properties
#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size ip = 0; ip < np; ip++) {
    if (do_abort) continue;
    try {
      //! FIXME: Send in full path point instead of this
      get_stepwise_clearsky_propmat(ws,
                                    ppvar_propmat[ip],
                                    ppvar_nlte[ip],
                                    ppvar_dpropmat[ip],
                                    ppvar_dnlte[ip],
                                    propmat_clearsky_agenda,
                                    jacobian_targets,
                                    ppvar_f[ip],
                                    rad_path[ip],
                                    ppvar_atm[ip]);
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

void spectral_radiance_path_sourceFromPropmat(
    ArrayOfStokvecVector &spectral_radiance_path_source,
    ArrayOfStokvecMatrix &spectral_radiance_path_source_jacobian,
    const ArrayOfPropmatVector &ppvar_propmat,
    const ArrayOfStokvecVector &ppvar_nlte,
    const ArrayOfPropmatMatrix &ppvar_dpropmat,
    const ArrayOfStokvecMatrix &ppvar_dnlte,
    const ArrayOfVector &ppvar_f,
    const ArrayOfAtmPoint &ppvar_atm,
    const JacobianTargets &jacobian_targets) try {
  ArrayOfString fail_msg;
  bool do_abort = false;

  const Index np = ppvar_atm.size();
  if (np == 0) {
    spectral_radiance_path_source.resize(0);
    spectral_radiance_path_source_jacobian.resize(0);
    return;
  }

  const Index nf = ppvar_propmat.front().size();
  const Index nq = jacobian_targets.target_count();

  const Index it =
      jacobian_targets.target_position<Jacobian::AtmTarget>(Atm::Key::t);

  spectral_radiance_path_source.resize(np, StokvecVector(nf));
  spectral_radiance_path_source_jacobian.resize(np, StokvecMatrix(nq, nf));

  // Loop ppath points and determine radiative properties
#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Index ip = 0; ip < np; ip++) {
    if (do_abort) continue;
    try {
      rtepack::source::level_nlte(spectral_radiance_path_source[ip],
                                  spectral_radiance_path_source_jacobian[ip],
                                  ppvar_propmat[ip],
                                  ppvar_nlte[ip],
                                  ppvar_dpropmat[ip],
                                  ppvar_dnlte[ip],
                                  ppvar_f[ip],
                                  ppvar_atm[ip].temperature,
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

void ppvar_tramatCalc(ArrayOfMuelmatVector &ppvar_tramat,
                      ArrayOfArrayOfMuelmatMatrix &ppvar_dtramat,
                      Vector &ppvar_distance,
                      ArrayOfArrayOfVector &ppvar_ddistance,
                      const ArrayOfPropmatVector &ppvar_propmat,
                      const ArrayOfPropmatMatrix &ppvar_dpropmat,
                      const ArrayOfPropagationPathPoint &rad_path,
                      const ArrayOfAtmPoint &ppvar_atm,
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
    ppvar_tramat.resize(0);
    ppvar_dtramat.resize(2, ArrayOfMuelmatMatrix{});
    ppvar_distance.resize(0);
    ppvar_ddistance.resize(2, ArrayOfVector{});
    return;
  }

  const Index nf = ppvar_propmat.front().size();
  const Index nq = jacobian_targets.target_count();

  ppvar_tramat.resize(np, MuelmatVector(nf));
  ppvar_dtramat.resize(2, ArrayOfMuelmatMatrix(np, MuelmatMatrix(nq, nf)));
  ppvar_distance.resize(np);
  ppvar_ddistance.resize(2, ArrayOfVector(np, Vector(nq, 0)));

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size ip = 1; ip < np; ip++) {
    if (do_abort) continue;
    try {
      ppvar_distance[ip - 1] =
          path::distance(rad_path[ip - 1].pos, rad_path[ip].pos, surface_field);
      if (hse_derivative and temperature_derivative_position >= 0) {
        ppvar_ddistance[0][ip][temperature_derivative_position] =
            ppvar_distance[ip - 1] / (2.0 * ppvar_atm[ip - 1].temperature);
        ppvar_ddistance[1][ip][temperature_derivative_position] =
            ppvar_distance[ip - 1] / (2.0 * ppvar_atm[ip].temperature);
      }

      two_level_exp(ppvar_tramat[ip],
                    ppvar_dtramat[0][ip],
                    ppvar_dtramat[1][ip],
                    ppvar_propmat[ip - 1],
                    ppvar_propmat[ip],
                    ppvar_dpropmat[ip - 1],
                    ppvar_dpropmat[ip],
                    ppvar_distance[ip - 1],
                    ppvar_ddistance[0][ip],
                    ppvar_ddistance[1][ip]);
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

void ppvar_atmFromPath(ArrayOfAtmPoint &ppvar_atm,
                       const ArrayOfPropagationPathPoint &rad_path,
                       const AtmField &atm_field) try {
  forward_atm_path(atm_path_resize(ppvar_atm, rad_path), rad_path, atm_field);
}
ARTS_METHOD_ERROR_CATCH

void ppvar_fFromPath(ArrayOfVector &ppvar_f,
                     const Vector &f_grid,
                     const ArrayOfPropagationPathPoint &rad_path,
                     const ArrayOfAtmPoint &ppvar_atm,
                     const Numeric &rte_alonglos_v) try {
  forward_path_freq(path_freq_resize(ppvar_f, f_grid, ppvar_atm),
                    f_grid,
                    rad_path,
                    ppvar_atm,
                    rte_alonglos_v);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radiance_pathCalcEmission(
    ArrayOfStokvecVector &spectral_radiance_path,
    ArrayOfStokvecMatrix &spectral_radiance_path_jacobian,
    const StokvecVector &background_rad,
    const ArrayOfStokvecVector &spectral_radiance_path_source,
    const ArrayOfStokvecMatrix &spectral_radiance_path_source_jacobian,
    const ArrayOfMuelmatVector &ppvar_tramat,
    const ArrayOfMuelmatVector &ppvar_cumtramat,
    const ArrayOfArrayOfMuelmatMatrix &ppvar_dtramat) try {
  const Size np = spectral_radiance_path_source.size();

  ARTS_USER_ERROR_IF(
      np not_eq spectral_radiance_path_source_jacobian.size(),
      "spectral_radiance_path_source_jacobian must have (np) elements")
  ARTS_USER_ERROR_IF(np not_eq ppvar_tramat.size(),
                     "ppvar_tramat must have (np) elements")
  ARTS_USER_ERROR_IF(np not_eq ppvar_cumtramat.size(),
                     "ppvar_cumtramat must have (np) elements")
  ARTS_USER_ERROR_IF(
      2 not_eq ppvar_dtramat.size() or
          ppvar_dtramat.front().size() not_eq ppvar_dtramat.back().size() or
          ppvar_dtramat.front().size() not_eq np,
      "ppvar_dtramat must (2 x np) elements")

  if (np == 0) {
    spectral_radiance_path.resize(0);
    spectral_radiance_path_jacobian.resize(0);
    return;
  }

  const Index nq = spectral_radiance_path_source_jacobian.front().nrows();
  const Index nf = spectral_radiance_path_source_jacobian.front().ncols();

  const auto test_nf = [nf](auto &v) { return v.size() not_eq nf; };
  ARTS_USER_ERROR_IF(nf not_eq background_rad.size(),
                     "background_rad must have nf elements")
  ARTS_USER_ERROR_IF(
      std::any_of(spectral_radiance_path_source.begin(),
                  spectral_radiance_path_source.end(),
                  test_nf),
      "spectral_radiance_path_source must have (nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ppvar_tramat.begin(), ppvar_tramat.end(), test_nf),
      "ppvar_tramat must have (nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ppvar_cumtramat.begin(), ppvar_cumtramat.end(), test_nf),
      "spectral_radiance_path_source must have (nf) inner elements")

  const auto test_nfnq = [nf, nq](auto &v) {
    return v.ncols() not_eq nf or v.nrows() not_eq nq;
  };
  ARTS_USER_ERROR_IF(
      std::any_of(spectral_radiance_path_source_jacobian.begin(),
                  spectral_radiance_path_source_jacobian.end(),
                  test_nfnq),
      "spectral_radiance_path_source_jacobian must have (nq x nf) inner elements")
  ARTS_USER_ERROR_IF(std::any_of(ppvar_dtramat.front().begin(),
                                 ppvar_dtramat.front().end(),
                                 test_nfnq),
                     "ppvar_dtramat must have (nq x nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(
          ppvar_dtramat.back().begin(), ppvar_dtramat.back().end(), test_nfnq),
      "ppvar_dtramat must have (nq x nf) inner elements")

  spectral_radiance_path.resize(np, background_rad);
  spectral_radiance_path_jacobian.resize(np, StokvecMatrix(nq, nf));
  for (Index ip = np - 2; ip >= 0; ip--) {
    spectral_radiance_path[ip] = spectral_radiance_path[ip + 1];
    two_level_linear_emission_step(
        spectral_radiance_path[ip],
        spectral_radiance_path_jacobian[ip],
        spectral_radiance_path_jacobian[ip + 1],
        spectral_radiance_path_source[ip],
        spectral_radiance_path_source[ip + 1],
        spectral_radiance_path_source_jacobian[ip],
        spectral_radiance_path_source_jacobian[ip + 1],
        ppvar_tramat[ip + 1],
        ppvar_cumtramat[ip],
        ppvar_dtramat[0][ip + 1],
        ppvar_dtramat[1][ip + 1]);
  }
}
ARTS_METHOD_ERROR_CATCH

void ppvar_cumtramatForward(ArrayOfMuelmatVector &ppvar_cumtramat,
                            const ArrayOfMuelmatVector &ppvar_tramat) {
  ppvar_cumtramat = forward_cumulative_transmission(ppvar_tramat);
}

void ppvar_cumtramatReverse(ArrayOfMuelmatVector &ppvar_cumtramat,
                            const ArrayOfMuelmatVector &ppvar_tramat) {
  ppvar_cumtramat = reverse_cumulative_transmission(ppvar_tramat);
}
