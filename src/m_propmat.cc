#include <array_algo.h>
#include <arts_omp.h>
#include <workspace.h>

void spectral_propmat_pathFromPath(
    const Workspace &ws,
    ArrayOfPropmatVector &spectral_propmat_path,
    ArrayOfStokvecVector &spectral_nlte_srcvec_path,
    ArrayOfPropmatMatrix &spectral_propmat_jac_path,
    ArrayOfStokvecMatrix &spectral_nlte_srcvec_jac_path,
    const Agenda &spectral_propmat_agenda,
    const ArrayOfAscendingGrid &freq_grid_path,
    const ArrayOfVector3 &freq_wind_shift_jac_path,
    const JacobianTargets &jac_targets,
    const ArrayOfPropagationPathPoint &ray_path,
    const ArrayOfAtmPoint &atm_path) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      not arr::same_size(
          ray_path, atm_path, freq_grid_path, freq_wind_shift_jac_path),
      R"(Not same size:

ray_path                 size: {} element(s)
atm_path                 size: {} element(s)
freq_grid_path           size: {} element(s)
freq_wind_shift_jac_path size: {} element(s)
)",
      ray_path.size(),
      atm_path.size(),
      freq_grid_path.size(),
      freq_wind_shift_jac_path.size())

  const Size np = ray_path.size();
  spectral_propmat_path.resize(np);
  spectral_nlte_srcvec_path.resize(np);
  spectral_propmat_jac_path.resize(np);
  spectral_nlte_srcvec_jac_path.resize(np);

  String error{};

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size ip = 0; ip < np; ip++) {
    try {
      spectral_propmat_agendaExecute(ws,
                                     spectral_propmat_path[ip],
                                     spectral_nlte_srcvec_path[ip],
                                     spectral_propmat_jac_path[ip],
                                     spectral_nlte_srcvec_jac_path[ip],
                                     freq_grid_path[ip],
                                     freq_wind_shift_jac_path[ip],
                                     jac_targets,
                                     {},
                                     ray_path[ip],
                                     atm_path[ip],
                                     spectral_propmat_agenda);
    } catch (const std::runtime_error &e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH

void spectral_propmat_pathAdaptiveHalfPath(
    const Workspace &ws,
    ArrayOfPropmatVector &spectral_propmat_path,
    ArrayOfStokvecVector &spectral_nlte_srcvec_path,
    ArrayOfPropmatMatrix &spectral_propmat_jac_path,
    ArrayOfStokvecMatrix &spectral_nlte_srcvec_jac_path,
    ArrayOfAscendingGrid &freq_grid_path,
    ArrayOfVector3 &freq_wind_shift_jac_path,
    ArrayOfPropagationPathPoint &ray_path,
    ArrayOfAtmPoint &atm_path,
    const Agenda &spectral_propmat_agenda,
    const JacobianTargets &jac_targets,
    const AscendingGrid &freq_grid,
    const AtmField &atm_field,
    const SurfaceField &surf_field,
    const Numeric &max_stepsize,
    const Numeric &max_tau,
    const Numeric &cutoff_tau) try {
  ARTS_TIME_REPORT

  spectral_propmat_pathFromPath(ws,
                                spectral_propmat_path,
                                spectral_nlte_srcvec_path,
                                spectral_propmat_jac_path,
                                spectral_nlte_srcvec_jac_path,
                                spectral_propmat_agenda,
                                freq_grid_path,
                                freq_wind_shift_jac_path,
                                jac_targets,
                                ray_path,
                                atm_path);

  const Size nf = freq_grid.size();

  if (ray_path.size() < 2) return;

  const auto &ellipsoid = surf_field.ellipsoid;

  Vector sum_tau(nf, 0.0);
  ArrayOfPropagationPathPoint nrp;
  ArrayOfAtmPoint nap;
  ArrayOfAscendingGrid nfgp;
  ArrayOfVector3 nfwsjp;
  ArrayOfPropmatVector nspp;
  ArrayOfStokvecVector nsnsp;
  ArrayOfPropmatMatrix nspjp;
  ArrayOfStokvecMatrix nsnsjp;

  for (Size ip = 0; ip < ray_path.size() - 1; ip++) {
    const auto &p0  = spectral_propmat_path[ip];
    const auto &p1  = spectral_propmat_path[ip + 1];
    const auto &rp0 = ray_path[ip];
    const auto &rp1 = ray_path[ip + 1];
    const Numeric r = path::distance(rp0.pos, rp1.pos, ellipsoid);

    Numeric max_dtau = -1.0;
    for (Size f = 0; f < nf; f++) {
      const Numeric dtau0 = p0[f].A() * r;
      const Numeric dtau1 = p1[f].A() * r;
      const Numeric dtau  = std::max(dtau0, dtau1);
      if (sum_tau[f] < cutoff_tau and dtau > max_tau and dtau > max_dtau) {
        max_dtau = dtau;
      }
      sum_tau[f] += dtau;
    }

    if (max_dtau > 0.0) {
      nrp.resize(2);
      nrp.front() = rp0;
      nrp.back()  = rp1;
      const Numeric new_r =
          std::min(std::nextafter(0.5 * r, 0.0),  // at least 3 points
                   std::max(max_stepsize, max_tau / max_dtau * r));
      path::fill_geometric_by_half_steps(nrp, surf_field, new_r);
      nrp.erase(nrp.begin());
      nrp.pop_back();
      if (not nrp.empty()) {
        atm_pathFromPath(nap, nrp, atm_field);
        freq_grid_pathFromPath(nfgp, nfwsjp, freq_grid, nrp, nap);
        spectral_propmat_pathFromPath(ws,
                                      nspp,
                                      nsnsp,
                                      nspjp,
                                      nsnsjp,
                                      spectral_propmat_agenda,
                                      nfgp,
                                      nfwsjp,
                                      jac_targets,
                                      nrp,
                                      nap);
        ray_path.insert_range(ray_path.begin() + ip+1, nrp);
        atm_path.insert_range(atm_path.begin() + ip+1 , nap);
        freq_grid_path.insert_range(freq_grid_path.begin() + ip+1 , nfgp);
        freq_wind_shift_jac_path.insert_range(
            freq_wind_shift_jac_path.begin() + ip+1 , nfwsjp);
        spectral_propmat_path.insert_range(
            spectral_propmat_path.begin() + ip+1 , nspp);
        spectral_nlte_srcvec_path.insert_range(
            spectral_nlte_srcvec_path.begin() + ip+1 , nsnsp);
        spectral_propmat_jac_path.insert_range(
            spectral_propmat_jac_path.begin() + ip+1 , nspjp);
        spectral_nlte_srcvec_jac_path.insert_range(
            spectral_nlte_srcvec_jac_path.begin() + ip+1 , nsnsjp);
        ip += nrp.size();
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_propmat_path_species_splitFromPath(
    const Workspace &ws,
    ArrayOfArrayOfPropmatVector &spectral_propmat_path_species_split,
    ArrayOfArrayOfStokvecVector &spectral_nlte_srcvec_path_species_split,
    ArrayOfArrayOfPropmatMatrix &spectral_propmat_jac_path_species_split,
    ArrayOfArrayOfStokvecMatrix &spectral_nlte_srcvec_jac_path_species_split,
    const Agenda &spectral_propmat_agenda,
    const ArrayOfAscendingGrid &freq_grid_path,
    const ArrayOfVector3 &freq_wind_shift_jac_path,
    const JacobianTargets &jac_targets,
    const ArrayOfPropagationPathPoint &ray_path,
    const ArrayOfAtmPoint &atm_path,
    const ArrayOfSpeciesEnum &select_species_list) try {
  ARTS_TIME_REPORT

  const Size ns = select_species_list.size();
  const Size np = ray_path.size();

  spectral_propmat_path_species_split.resize(ns);
  spectral_nlte_srcvec_path_species_split.resize(ns);
  spectral_propmat_jac_path_species_split.resize(ns);
  spectral_nlte_srcvec_jac_path_species_split.resize(ns);
  for (auto &s : spectral_propmat_path_species_split) s.resize(np);
  for (auto &s : spectral_nlte_srcvec_path_species_split) s.resize(np);
  for (auto &s : spectral_propmat_jac_path_species_split) s.resize(np);
  for (auto &s : spectral_nlte_srcvec_jac_path_species_split) s.resize(np);

  String error{};

#pragma omp parallel for if (!arts_omp_in_parallel()) collapse(2)
  for (Size is = 0; is < ns; is++) {
    for (Size ip = 0; ip < np; ip++) {
      try {
        spectral_propmat_agendaExecute(
            ws,
            spectral_propmat_path_species_split[is][ip],
            spectral_nlte_srcvec_path_species_split[is][ip],
            spectral_propmat_jac_path_species_split[is][ip],
            spectral_nlte_srcvec_jac_path_species_split[is][ip],
            freq_grid_path[ip],
            freq_wind_shift_jac_path[ip],
            jac_targets,
            select_species_list[is],
            ray_path[ip],
            atm_path[ip],
            spectral_propmat_agenda);
      } catch (const std::runtime_error &e) {
#pragma omp critical
        if (error.empty()) error = e.what();
      }
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH
