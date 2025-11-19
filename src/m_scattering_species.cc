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

void scat_speciesInit(ArrayOfScatteringSpecies& scattering_species) {
  scattering_species = ArrayOfScatteringSpecies{};
}

void spectral_propmat_scatSpectralInit(PropmatVector& spectral_propmat_scat,
                                       StokvecVector& spectral_absvec_scat,
                                       SpecmatMatrix& spectral_phamat_spectral,
                                       const AscendingGrid& freq_grid,
                                       const Index& legendre_degree) {
  ARTS_USER_ERROR_IF(legendre_degree < 0,
                     "The legendre_degree must be non-negative, is {}",
                     legendre_degree)

  spectral_propmat_scat.resize(freq_grid.size());
  spectral_absvec_scat.resize(freq_grid.size());
  spectral_phamat_spectral.resize(freq_grid.size(), legendre_degree + 1);

  spectral_propmat_scat    = 0.0;
  spectral_absvec_scat     = 0.0;
  spectral_phamat_spectral = Complex{0.0};
}

void spectral_propmat_scatAddSpectralScatteringSpeciesTRO(
    PropmatVector& spectral_propmat_scat,
    StokvecVector& spectral_absvec_scat,
    SpecmatMatrix& spectral_phamat_spectral,
    const AscendingGrid& freq_grid,
    const AtmPoint& atm_point,
    const ArrayOfScatteringSpecies& scattering_species) try {
  const Index L = spectral_phamat_spectral.ncols();

  ARTS_USER_ERROR_IF(L < 1, "Need at least one Legendre coefficient")
  ARTS_USER_ERROR_IF(scattering_species.species.empty(),
                     "No scattering species")

  const auto [phase_matrix_opt, extinction_matrix, absorption_vector] =
      scattering_species.get_bulk_scattering_properties_tro_spectral(
          atm_point, freq_grid, L - 1);

  ARTS_USER_ERROR_IF(not phase_matrix_opt.has_value(), "No phase matrix")

  const auto& phase_matrix = phase_matrix_opt.value();

  ARTS_USER_ERROR_IF(
      not same_shape<2>(phase_matrix, spectral_phamat_spectral) or
          not same_shape<1>(extinction_matrix,
                            absorption_vector,
                            spectral_propmat_scat,
                            spectral_absvec_scat),
      R"(Shape mismatch in return from scattering_species.get_bulk_scattering_properties_tro_spectral():

Phase matrix shape (must match):
  phase_matrix.shape(): [OUTPUT]    {0:B,}
  spectral_phamat_spectral.shape(): {1:B,}

Extinction matrix shape (must match):
  extinction_matrix.shape(): [OUTPUT] {2:B,}
  spectral_propmat_scat.shape():      {3:B,}

Absorption vector shape (must match):
  absorption_vector.shape(): [OUTPUT] {4:B,}
  spectral_absvec_scat.shape():       {5:B,}

Supporting variable sizes:
  freq_grid.size():          {6}
  scattering_species.size(): {7}
)",
      phase_matrix.shape(),
      spectral_phamat_spectral.shape(),
      extinction_matrix.shape(),
      spectral_propmat_scat.shape(),
      absorption_vector.shape(),
      spectral_absvec_scat.shape(),
      freq_grid.size(),
      scattering_species.species.size())

  spectral_propmat_scat    += extinction_matrix;
  spectral_phamat_spectral += phase_matrix;
  spectral_absvec_scat     += absorption_vector;
}
ARTS_METHOD_ERROR_CATCH

void spectral_propmat_pathAddScattering(
    ArrayOfPropmatVector& spectral_propmat_path,
    const ArrayOfPropmatVector& spectral_propmat_scat_path) try {
  const Size N = spectral_propmat_path.size();

  ARTS_USER_ERROR_IF(
      N != spectral_propmat_scat_path.size(),
      R"(The size of spectral_propmat_path and spectral_propmat_scat_path must be the same.
  spectral_propmat_path.size():      {}
  spectral_propmat_scat_path.size(): {}
)",
      spectral_propmat_path.size(),
      spectral_propmat_scat_path.size());

  if (N == 0) return;

  ARTS_USER_ERROR_IF(
      not all_same_shape<1>(spectral_propmat_path.front().shape(),
                            spectral_propmat_path,
                            spectral_propmat_scat_path),
      "The inner shapes of spectral_propmat_path and spectral_propmat_scat_path must be the same (first elem shape: {:B,}).",
      spectral_propmat_path.front().shape())

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < spectral_propmat_path.size(); i++) {
    spectral_propmat_path[i] += spectral_propmat_scat_path[i];
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_propmat_scat_pathFromSpectralAgenda(
    const Workspace& ws,
    ArrayOfPropmatVector& spectral_propmat_scat_path,
    ArrayOfStokvecVector& spectral_absvec_scat_path,
    ArrayOfSpecmatMatrix& spectral_phamat_spectral_path,
    const ArrayOfAscendingGrid& freq_grid_path,
    const ArrayOfAtmPoint& atm_path,
    const Index& legendre_degree,
    const Agenda& spectral_propmat_scat_spectral_agenda) try {
  const Size N = freq_grid_path.size();

  ARTS_USER_ERROR_IF(
      not arr::same_size(freq_grid_path, atm_path),
      R"(The size of freq_grid_path and atm_path must be the same.
  freq_grid_path.size(): {}
  atm_path.size():       {}
)",
      freq_grid_path.size(),
      atm_path.size());

  spectral_propmat_scat_path.resize(N);
  spectral_absvec_scat_path.resize(N);
  spectral_phamat_spectral_path.resize(N);

  std::string error{};
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < N; i++) {
    try {
      spectral_propmat_scat_spectral_agendaExecute(
          ws,
          spectral_propmat_scat_path[i],
          spectral_absvec_scat_path[i],
          spectral_phamat_spectral_path[i],
          freq_grid_path[i],
          atm_path[i],
          legendre_degree,
          spectral_propmat_scat_spectral_agenda);
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
