#include <array_algo.h>
#include <arts_omp.h>
#include <workspace.h>

void spectral_propmat_pathFromPath(
    const Workspace &ws,
    ArrayOfPropmatVector &spectral_propmat_path,
    ArrayOfStokvecVector &ray_path_source_vector_nonlte,
    ArrayOfPropmatMatrix &spectral_propmat_jac_path,
    ArrayOfStokvecMatrix &ray_path_source_vector_nonlte_jacobian,
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

ray_path                               size: {} element(s)
atm_path             size: {} element(s)
freq_grid_path                size: {} element(s)
freq_wind_shift_jac_path size: {} element(s)
)",
      ray_path.size(),
      atm_path.size(),
      freq_grid_path.size(),
      freq_wind_shift_jac_path.size())

  const Size np = ray_path.size();
  spectral_propmat_path.resize(np);
  ray_path_source_vector_nonlte.resize(np);
  spectral_propmat_jac_path.resize(np);
  ray_path_source_vector_nonlte_jacobian.resize(np);

  String error{};

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size ip = 0; ip < np; ip++) {
    try {
      spectral_propmat_agendaExecute(ws,
                                     spectral_propmat_path[ip],
                                     ray_path_source_vector_nonlte[ip],
                                     spectral_propmat_jac_path[ip],
                                     ray_path_source_vector_nonlte_jacobian[ip],
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
