#include <array_algo.h>
#include <arts_omp.h>
#include <workspace.h>

void ray_path_propagation_matrixFromPath(
    const Workspace &ws,
    ArrayOfPropmatVector &ray_path_propagation_matrix,
    ArrayOfStokvecVector &ray_path_source_vector_nonlte,
    ArrayOfPropmatMatrix &ray_path_propagation_matrix_jacobian,
    ArrayOfStokvecMatrix &ray_path_source_vector_nonlte_jacobian,
    const Agenda &propagation_matrix_agenda,
    const ArrayOfAscendingGrid &ray_path_frequency_grid,
    const ArrayOfVector3 &ray_path_frequency_wind_shift_jacobian,
    const JacobianTargets &jacobian_targets,
    const ArrayOfPropagationPathPoint &ray_path,
    const ArrayOfAtmPoint &ray_path_atm_point) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not arr::same_size(ray_path,
                                        ray_path_atm_point,
                                        ray_path_frequency_grid,
                                        ray_path_frequency_wind_shift_jacobian),
                     R"(Not same size:

ray_path                               size: {} element(s)
ray_path_atm_point             size: {} element(s)
ray_path_frequency_grid                size: {} element(s)
ray_path_frequency_wind_shift_jacobian size: {} element(s)
)",
                     ray_path.size(),
                     ray_path_atm_point.size(),
                     ray_path_frequency_grid.size(),
                     ray_path_frequency_wind_shift_jacobian.size())

  const Size np = ray_path.size();
  ray_path_propagation_matrix.resize(np);
  ray_path_source_vector_nonlte.resize(np);
  ray_path_propagation_matrix_jacobian.resize(np);
  ray_path_source_vector_nonlte_jacobian.resize(np);

  String error{};

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Size ip = 0; ip < np; ip++) {
    try {
      propagation_matrix_agendaExecute(
          ws,
          ray_path_propagation_matrix[ip],
          ray_path_source_vector_nonlte[ip],
          ray_path_propagation_matrix_jacobian[ip],
          ray_path_source_vector_nonlte_jacobian[ip],
          ray_path_frequency_grid[ip],
          ray_path_frequency_wind_shift_jacobian[ip],
          jacobian_targets,
          {},
          ray_path[ip],
          ray_path_atm_point[ip],
          propagation_matrix_agenda);
    } catch (const std::runtime_error &e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH

void ray_path_propagation_matrix_species_splitFromPath(
    const Workspace &ws,
    ArrayOfArrayOfPropmatVector &ray_path_propagation_matrix_species_split,
    ArrayOfArrayOfStokvecVector
        &ray_path_propagation_matrix_source_vector_nonlte_species_split,
    ArrayOfArrayOfPropmatMatrix
        &ray_path_propagation_matrix_jacobian_species_split,
    ArrayOfArrayOfStokvecMatrix &
        ray_path_propagation_matrix_source_vector_nonlte_jacobian_species_split,
    const Agenda &propagation_matrix_agenda,
    const ArrayOfAscendingGrid &ray_path_frequency_grid,
    const ArrayOfVector3 &ray_path_frequency_wind_shift_jacobian,
    const JacobianTargets &jacobian_targets,
    const ArrayOfPropagationPathPoint &ray_path,
    const ArrayOfAtmPoint &ray_path_atm_point,
    const ArrayOfSpeciesEnum &select_species_list) try {
  ARTS_TIME_REPORT

  const Size ns = select_species_list.size();
  const Size np = ray_path.size();

  ray_path_propagation_matrix_species_split.resize(ns);
  ray_path_propagation_matrix_source_vector_nonlte_species_split.resize(ns);
  ray_path_propagation_matrix_jacobian_species_split.resize(ns);
  ray_path_propagation_matrix_source_vector_nonlte_jacobian_species_split
      .resize(ns);
  for (auto &s : ray_path_propagation_matrix_species_split) s.resize(np);
  for (auto &s : ray_path_propagation_matrix_source_vector_nonlte_species_split)
    s.resize(np);
  for (auto &s : ray_path_propagation_matrix_jacobian_species_split)
    s.resize(np);
  for (auto &s :
       ray_path_propagation_matrix_source_vector_nonlte_jacobian_species_split)
    s.resize(np);

  String error{};

#pragma omp parallel for if (!arts_omp_in_parallel()) collapse(2)
  for (Size is = 0; is < ns; is++) {
    for (Size ip = 0; ip < np; ip++) {
      try {
        propagation_matrix_agendaExecute(
            ws,
            ray_path_propagation_matrix_species_split[is][ip],
            ray_path_propagation_matrix_source_vector_nonlte_species_split[is]
                                                                          [ip],
            ray_path_propagation_matrix_jacobian_species_split[is][ip],
            ray_path_propagation_matrix_source_vector_nonlte_jacobian_species_split
                [is][ip],
            ray_path_frequency_grid[ip],
            ray_path_frequency_wind_shift_jacobian[ip],
            jacobian_targets,
            select_species_list[is],
            ray_path[ip],
            ray_path_atm_point[ip],
            propagation_matrix_agenda);
      } catch (const std::runtime_error &e) {
#pragma omp critical
        if (error.empty()) error = e.what();
      }
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH
