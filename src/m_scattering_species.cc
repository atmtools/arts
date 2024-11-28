#include <arts_omp.h>
#include <workspace.h>

#include "debug.h"
#include "mh_checks.h"

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
    ComplexMuelmatMatrix& phase_matrix_scattering_spectral,
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
    ComplexMuelmatMatrix& phase_matrix_scattering_spectral,
    const AscendingGrid& frequency_grid,
    const AtmPoint& atmospheric_point,
    const ArrayOfScatteringSpecies& scattering_species) try {
  const Index F = frequency_grid.size();
  const Index L = phase_matrix_scattering_spectral.ncols();

  ARTS_USER_ERROR_IF(L < 1, "Must have at least one legendre coefficient")
  ARTS_USER_ERROR_IF(
      not same_shape(std::array{F}, propagation_matrix_scattering),
      "The shape of propagation_matrix_scattering must be {:B,}, is {:B,}",
      std::array{F},
      propagation_matrix_scattering.shape())
  ARTS_USER_ERROR_IF(
      not same_shape(std::array{F}, absorption_vector_scattering),
      "The shape of absorption_vector_scattering must be {:B,}, is {:B,}",
      std::array{F},
      absorption_vector_scattering.shape())
  ARTS_USER_ERROR_IF(
      not same_shape(std::array{F, L}, phase_matrix_scattering_spectral),
      "The shape of phase_matrix_scattering_spectral must be {:B,}, is {:B,}",
      std::array{F, L},
      phase_matrix_scattering_spectral.shape())

  const auto [phase_matrix_opt, extinction_matrix, absorption_vector] =
      scattering_species.get_bulk_scattering_properties_tro_spectral(
          atmospheric_point, frequency_grid, L - 1);

  ARTS_USER_ERROR_IF(not phase_matrix_opt.has_value(),
                     "Expects phase_matrix to have a value")

  const auto& phase_matrix = phase_matrix_opt.value();

  ARTS_USER_ERROR_IF(
      not same_shape(std::array{Index{1}, F, L, Index{6}}, phase_matrix),
      "The shape of phase_matrix from scattering_species must be {:B,}, is {:B,}",
      std::array{Index{1}, F, L, Index{6}},
      phase_matrix.shape())
  for (Index i = 0; i < F; i++) {
    for (Index j = 0; j < L; j++) {
      phase_matrix_scattering_spectral(i, j)(0, 0) += phase_matrix(0, i, j, 0);
      phase_matrix_scattering_spectral(i, j)(1, 1) += phase_matrix(0, i, j, 0);
      phase_matrix_scattering_spectral(i, j)(2, 2) += phase_matrix(0, i, j, 0);
      phase_matrix_scattering_spectral(i, j)(3, 3) += phase_matrix(0, i, j, 0);
    }
  }

  ARTS_USER_ERROR_IF(
      not same_shape(std::array{Index{1}, F, Index{1}}, extinction_matrix),
      "The shape of extinction_matrix must be {:B,}, is {:B,}",
      std::array{Index{1}, F, Index{1}},
      extinction_matrix.shape())
  for (Index iv = 0; iv < F; iv++) {
    propagation_matrix_scattering[iv].A() += extinction_matrix(0, iv, 0);
  }

  ARTS_USER_ERROR_IF(
      not same_shape(std::array{Index{1}, F, Index{1}}, absorption_vector),
      "The shape of absorption_vector must be {:B,}, is {:B,}",
      std::array{Index{1}, F, Index{1}},
      absorption_vector.shape())
  for (Index iv = 0; iv < F; iv++) {
    absorption_vector_scattering[iv].I() += absorption_vector(0, iv, 0);
  }
}
ARTS_METHOD_ERROR_CATCH

void ray_path_propagation_matrixAddTotallyRandomOrientationSpectral(
    ArrayOfPropmatVector& ray_path_propagation_matrix,
    const ArrayOfPropmatVector& ray_path_propagation_matrix_scattering) try {
  const Size N = ray_path_propagation_matrix.size();

  ARTS_USER_ERROR_IF(
      N != ray_path_propagation_matrix_scattering.size(),
      R"(The size of ray_path_propagation_matrix and ray_path_propagation_matrix_scattering must be the same.
  ray_path_propagation_matrix.size():                                                {}
  ray_path_propagation_matrix_scattering.size(): {}
)",
      ray_path_propagation_matrix.size(),
      ray_path_propagation_matrix_scattering.size());

  if (N == 0) return;

  ARTS_USER_ERROR_IF(
      not all_same_shape(ray_path_propagation_matrix.front().shape(),
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
    ArrayOfComplexMuelmatMatrix& ray_path_phase_matrix_scattering_spectral,
    const ArrayOfAscendingGrid& ray_path_frequency_grid,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const Index& legendre_degree,
    const Agenda& propagation_matrix_scattering_spectral_agenda) try {
  const Size N = ray_path_frequency_grid.size();

  ARTS_USER_ERROR_IF(
      not all_same_size(ray_path_frequency_grid, ray_path_atmospheric_point),
      R"(The size of ray_path_frequency_grid and ray_path_atmospheric_point must be the same.
  ray_path_frequency_grid.size():    {}
  ray_path_atmospheric_point.size(): {}
)",
      ray_path_frequency_grid.size(),
      ray_path_atmospheric_point.size());

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
          ray_path_atmospheric_point[i],
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

  ARTS_USER_ERROR_IF(not error.empty(), "{}", error);
}
ARTS_METHOD_ERROR_CATCH
