#include <array_algo.h>
#include <arts_omp.h>
#include <debug.h>
#include <workspace.h>

void legendre_degreeFromDisortSettings(Index& legendre_degree,
                                       const DisortSettings& disort_settings) {
  legendre_degree = disort_settings.legendre_polynomial_dimension - 1;
  ARTS_USER_ERROR_IF(legendre_degree < 0,
                     "The legendre_degree must be non-negative, is {}",
                     legendre_degree)
}

void scattering_speciesInit(ArrayOfScatteringSpecies& scattering_species) {
  scattering_species = ArrayOfScatteringSpecies{};
}

void propagation_matrix_scatteringSpectralInit(
    PropmatVector& propagation_matrix_scattering,
    StokvecVector& absorption_vector_scattering,
    SpecmatMatrix& phase_matrix_scattering_spectral,
    const AscendingGrid& frequency_grid,
    const Index& legendre_degree) {
  ARTS_USER_ERROR_IF(legendre_degree < 0,
                     "The legendre_degree must be non-negative, is {}",
                     legendre_degree)

  propagation_matrix_scattering.resize(frequency_grid.size());
  absorption_vector_scattering.resize(frequency_grid.size());
  phase_matrix_scattering_spectral.resize(frequency_grid.size(),
                                          legendre_degree + 1);

  propagation_matrix_scattering    = 0.0;
  absorption_vector_scattering     = 0.0;
  phase_matrix_scattering_spectral = Complex{0.0};
}

void propagation_matrix_scatteringAddSpectralScatteringSpeciesTRO(
    PropmatVector& propagation_matrix_scattering,
    StokvecVector& absorption_vector_scattering,
    SpecmatMatrix& phase_matrix_scattering_spectral,
    const AscendingGrid& frequency_grid,
    const AtmPoint& atm_point,
    const ArrayOfScatteringSpecies& scattering_species) try {
  const Index L = phase_matrix_scattering_spectral.ncols();

  ARTS_USER_ERROR_IF(L < 1, "Need at least one Legendre coefficient")
  ARTS_USER_ERROR_IF(scattering_species.species.empty(),
                     "No scattering species")

  const auto [phase_matrix_opt, extinction_matrix, absorption_vector] =
      scattering_species.get_bulk_scattering_properties_tro_spectral(
          atm_point, frequency_grid, L - 1);

  ARTS_USER_ERROR_IF(not phase_matrix_opt.has_value(), "No phase matrix")

  const auto& phase_matrix = phase_matrix_opt.value();

  ARTS_USER_ERROR_IF(
      not same_shape<2>(phase_matrix, phase_matrix_scattering_spectral) or
          not same_shape<1>(extinction_matrix,
                            absorption_vector,
                            propagation_matrix_scattering,
                            absorption_vector_scattering),
      R"(Shape mismatch in return from scattering_species.get_bulk_scattering_properties_tro_spectral():

Phase matrix shape (must match):
  phase_matrix.shape(): [OUTPUT]            {0:B,}
  phase_matrix_scattering_spectral.shape(): {1:B,}

Extinction matrix shape (must match):
  extinction_matrix.shape(): [OUTPUT]       {2:B,}
  propagation_matrix_scattering.shape():    {3:B,}

Absorption vector shape (must match):
  absorption_vector.shape(): [OUTPUT]       {4:B,}
  absorption_vector_scattering.shape():     {5:B,}

Supporting variable sizes:
  frequency_grid.size():                    {6}
  scattering_species.size():                {7}
)",
      phase_matrix.shape(),
      phase_matrix_scattering_spectral.shape(),
      extinction_matrix.shape(),
      propagation_matrix_scattering.shape(),
      absorption_vector.shape(),
      absorption_vector_scattering.shape(),
      frequency_grid.size(),
      scattering_species.species.size())

  propagation_matrix_scattering    += extinction_matrix;
  phase_matrix_scattering_spectral += phase_matrix;
  absorption_vector_scattering     += absorption_vector;
}
ARTS_METHOD_ERROR_CATCH

void ray_path_propagation_matrixAddScattering(
    ArrayOfPropmatVector& ray_path_propagation_matrix,
    const ArrayOfPropmatVector& ray_path_propagation_matrix_scattering) try {
  const Size N = ray_path_propagation_matrix.size();

  ARTS_USER_ERROR_IF(
      N != ray_path_propagation_matrix_scattering.size(),
      R"(The size of ray_path_propagation_matrix and ray_path_propagation_matrix_scattering must be the same.
  ray_path_propagation_matrix.size():            {}
  ray_path_propagation_matrix_scattering.size(): {}
)",
      ray_path_propagation_matrix.size(),
      ray_path_propagation_matrix_scattering.size());

  if (N == 0) return;

  ARTS_USER_ERROR_IF(
      not all_same_shape<1>(ray_path_propagation_matrix.front().shape(),
                            ray_path_propagation_matrix,
                            ray_path_propagation_matrix_scattering),
      "The inner shapes of ray_path_propagation_matrix and ray_path_propagation_matrix_scattering must be the same (first elem shape: {:B,}).",
      ray_path_propagation_matrix.front().shape())

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < ray_path_propagation_matrix.size(); i++) {
    ray_path_propagation_matrix[i] += ray_path_propagation_matrix_scattering[i];
  }
}
ARTS_METHOD_ERROR_CATCH

void ray_path_propagation_matrix_scatteringFromSpectralAgenda(
    const Workspace& ws,
    ArrayOfPropmatVector& ray_path_propagation_matrix_scattering,
    ArrayOfStokvecVector& ray_path_absorption_vector_scattering,
    ArrayOfSpecmatMatrix& ray_path_phase_matrix_scattering_spectral,
    const ArrayOfAscendingGrid& ray_path_frequency_grid,
    const ArrayOfAtmPoint& ray_path_atm_point,
    const Index& legendre_degree,
    const Agenda& propagation_matrix_scattering_spectral_agenda) try {
  const Size N = ray_path_frequency_grid.size();

  ARTS_USER_ERROR_IF(
      not arr::same_size(ray_path_frequency_grid, ray_path_atm_point),
      R"(The size of ray_path_frequency_grid and ray_path_atm_point must be the same.
  ray_path_frequency_grid.size():    {}
  ray_path_atm_point.size(): {}
)",
      ray_path_frequency_grid.size(),
      ray_path_atm_point.size());

  ray_path_propagation_matrix_scattering.resize(N);
  ray_path_absorption_vector_scattering.resize(N);
  ray_path_phase_matrix_scattering_spectral.resize(N);

  std::string error{};
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < N; i++) {
    try {
      propagation_matrix_scattering_spectral_agendaExecute(
          ws,
          ray_path_propagation_matrix_scattering[i],
          ray_path_absorption_vector_scattering[i],
          ray_path_phase_matrix_scattering_spectral[i],
          ray_path_frequency_grid[i],
          ray_path_atm_point[i],
          legendre_degree,
          propagation_matrix_scattering_spectral_agenda);
    } catch (const std::exception& e) {
#pragma omp critical
      {
        if (error.empty()) {
          error = std::format(
              "Error in parallel-for. Index i={} (range: 0<=i<{}):\n{}\n",
              i,
              N,
              e.what());
        }
      }
    }
  }

  if (not error.empty()) throw std::runtime_error(error);
}
ARTS_METHOD_ERROR_CATCH
