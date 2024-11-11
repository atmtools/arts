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

void propagation_matrix_scattering_totally_random_orientation_spectralInit(
    PropmatVector&
        propagation_matrix_scattering_totally_random_orientation_spectral,
    StokvecVector&
        absorption_vector_scattering_totally_random_orientation_spectral,
    Matrix& phase_matrix_scattering_totally_random_orientation_spectral,
    const AscendingGrid& frequency_grid,
    const Index& legendre_degree) {
  ARTS_USER_ERROR_IF(legendre_degree < 0,
                     "The legendre_degree must be non-negative, is {}",
                     legendre_degree)

  propagation_matrix_scattering_totally_random_orientation_spectral.resize(
      frequency_grid.size());
  absorption_vector_scattering_totally_random_orientation_spectral.resize(
      frequency_grid.size());
  phase_matrix_scattering_totally_random_orientation_spectral.resize(
      frequency_grid.size(), legendre_degree + 1);

  propagation_matrix_scattering_totally_random_orientation_spectral = 0.0;
  absorption_vector_scattering_totally_random_orientation_spectral  = 0.0;
  phase_matrix_scattering_totally_random_orientation_spectral       = 0.0;
}

void propagation_matrix_scattering_totally_random_orientation_spectralAddScatteringSpecies(
    PropmatVector&
        propagation_matrix_scattering_totally_random_orientation_spectral,
    StokvecVector&
        absorption_vector_scattering_totally_random_orientation_spectral,
    Matrix& phase_matrix_scattering_totally_random_orientation_spectral,
    const AscendingGrid& frequency_grid,
    const AtmPoint& atmospheric_point,
    const ArrayOfScatteringSpecies& scattering_species) try {
  const Index F = frequency_grid.size();
  const Index L =
      phase_matrix_scattering_totally_random_orientation_spectral.ncols();

  ARTS_USER_ERROR_IF(
      not same_shape(
          std::array{F},
          propagation_matrix_scattering_totally_random_orientation_spectral),
      "The shape of propagation_matrix_scattering_totally_random_orientation_spectral must be {:B,}, is {:B,}",
      std::array{F},
      propagation_matrix_scattering_totally_random_orientation_spectral.shape())
  ARTS_USER_ERROR_IF(
      not same_shape(
          std::array{F},
          absorption_vector_scattering_totally_random_orientation_spectral),
      "The shape of absorption_vector_scattering_totally_random_orientation_spectral must be {:B,}, is {:B,}",
      std::array{F},
      absorption_vector_scattering_totally_random_orientation_spectral.shape())
  ARTS_USER_ERROR_IF(
      F != phase_matrix_scattering_totally_random_orientation_spectral.nrows(),
      "The shape of phase_matrix_scattering_totally_random_orientation_spectral must be {:B,}, is {:B,}",
      std::array{F, L},
      phase_matrix_scattering_totally_random_orientation_spectral.shape())

  //! FIXME: The following is a dummy implementation until the stokes-dim thing is fixed to stokes_dim=4
  const auto [phase_matrix, extinction_matrix, absorption_vector] =
      scattering_species.get_bulk_scattering_properties_tro_spectral<1>(
          atmospheric_point, frequency_grid, L - 1);

  if (phase_matrix.has_value()) {
    ARTS_USER_ERROR_IF(
        not same_shape(std::array{Index{1}, F, L, Index{1}},
                       phase_matrix.value()),
        "The shape of phase_matrix from scattering_species must be {:B,}, is {:B,}",
        std::array{Index{1}, F, L, Index{1}},
        phase_matrix->shape())
    phase_matrix_scattering_totally_random_orientation_spectral +=
        phase_matrix.value()(0, joker, joker, 0).real();
  }

  ARTS_USER_ERROR_IF(
      not same_shape(std::array{Index{1}, F, Index{1}}, extinction_matrix),
      "The shape of extinction_matrix must be {:B,}, is {:B,}",
      std::array{Index{1}, F, Index{1}},
      extinction_matrix.shape())
  for (Index iv = 0; iv < F; iv++) {
    propagation_matrix_scattering_totally_random_orientation_spectral[iv].A() +=
        extinction_matrix(0, iv, 0);
  }

  ARTS_USER_ERROR_IF(
      not same_shape(std::array{Index{1}, F, Index{1}}, absorption_vector),
      "The shape of absorption_vector must be {:B,}, is {:B,}",
      std::array{Index{1}, F, Index{1}},
      absorption_vector.shape())
  for (Index iv = 0; iv < F; iv++) {
    absorption_vector_scattering_totally_random_orientation_spectral[iv].I() +=
        absorption_vector(0, iv, 0);
  }
}
ARTS_METHOD_ERROR_CATCH

void ray_path_propagation_matrixAddTotallyRandomOrientationSpectral(
    ArrayOfPropmatVector& ray_path_propagation_matrix,
    const ArrayOfPropmatVector&
        ray_path_propagation_matrix_scattering_totally_random_orientation_spectral) try {
  const Size N = ray_path_propagation_matrix.size();

  ARTS_USER_ERROR_IF(
      N !=
          ray_path_propagation_matrix_scattering_totally_random_orientation_spectral
              .size(),
      R"(The size of ray_path_propagation_matrix and ray_path_propagation_matrix_scattering_totally_random_orientation_spectral must be the same.
  ray_path_propagation_matrix.size():                                                {}
  ray_path_propagation_matrix_scattering_totally_random_orientation_spectral.size(): {}
)",
      ray_path_propagation_matrix.size(),
      ray_path_propagation_matrix_scattering_totally_random_orientation_spectral
          .size());

  if (N == 0) return;

  ARTS_USER_ERROR_IF(
      not all_same_shape(
          ray_path_propagation_matrix.front().shape(),
          ray_path_propagation_matrix,
          ray_path_propagation_matrix_scattering_totally_random_orientation_spectral),
      "The inner shapes of ray_path_propagation_matrix and ray_path_propagation_matrix_scattering_totally_random_orientation_spectral must be the same (first elem shape: {:B,}).",
      ray_path_propagation_matrix.front().shape())

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < ray_path_propagation_matrix.size(); i++) {
    ray_path_propagation_matrix[i] +=
        ray_path_propagation_matrix_scattering_totally_random_orientation_spectral
            [i];
  }
}
ARTS_METHOD_ERROR_CATCH

void ray_path_propagation_matrix_scattering_totally_random_orientation_spectralFromAgenda(
    const Workspace& ws,
    ArrayOfPropmatVector&
        ray_path_propagation_matrix_scattering_totally_random_orientation_spectral,
    ArrayOfStokvecVector&
        ray_path_absorption_vector_scattering_totally_random_orientation_spectral,
    ArrayOfMatrix&
        ray_path_phase_matrix_scattering_totally_random_orientation_spectral,
    const ArrayOfAscendingGrid& ray_path_frequency_grid,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const Index& legendre_degree,
    const Agenda&
        propagation_matrix_scattering_totally_random_orientation_spectral_agenda) try {
  const Size N = ray_path_frequency_grid.size();

  ARTS_USER_ERROR_IF(
      not all_same_size(ray_path_frequency_grid, ray_path_atmospheric_point),
      R"(The size of ray_path_frequency_grid and ray_path_atmospheric_point must be the same.
  ray_path_frequency_grid.size():    {}
  ray_path_atmospheric_point.size(): {}
)",
      ray_path_frequency_grid.size(),
      ray_path_atmospheric_point.size());

  ray_path_propagation_matrix_scattering_totally_random_orientation_spectral
      .resize(N);
  ray_path_absorption_vector_scattering_totally_random_orientation_spectral
      .resize(N);
  ray_path_phase_matrix_scattering_totally_random_orientation_spectral.resize(
      N);

  std::string error{};
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < N; i++) {
    try {
      propagation_matrix_scattering_totally_random_orientation_spectral_agendaExecute(
          ws,
          ray_path_propagation_matrix_scattering_totally_random_orientation_spectral
              [i],
          ray_path_absorption_vector_scattering_totally_random_orientation_spectral
              [i],
          ray_path_phase_matrix_scattering_totally_random_orientation_spectral
              [i],
          ray_path_frequency_grid[i],
          ray_path_atmospheric_point[i],
          legendre_degree,
          propagation_matrix_scattering_totally_random_orientation_spectral_agenda);
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
