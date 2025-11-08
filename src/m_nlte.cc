#include <array_algo.h>
#include <arts_constants.h>
#include <arts_omp.h>
#include <auto_wsm.h>
#include <jacobian.h>
#include <workspace.h>

void abs_bandsSetNonLTE(AbsorptionBands& abs_bands) {
  ARTS_TIME_REPORT

  for (auto& [_, band] : abs_bands) {
    band.lineshape = LineByLineLineshape::VP_LINE_NLTE;
  }
}

void atm_fieldInitializeNonLTE(AtmField& atm_field,
                               const AbsorptionBands& abs_bands,
                               const Numeric& normalizing_factor) try {
  ARTS_TIME_REPORT

  atm_field.nlte =
      lbl::nlte::from_lte(atm_field, abs_bands, normalizing_factor);
}
ARTS_METHOD_ERROR_CATCH

void frequency_gridFitNonLTE(AscendingGrid& frequency_grid,
                             const AbsorptionBands& abs_bands,
                             const Numeric& df,
                             const Index& nf) try {
  ARTS_TIME_REPORT

  Vector freq;
  freq.reserve(count_lines(abs_bands) * nf);
  for (auto f0 : abs_bands | stdv::values |
                     stdv::transform([](auto& x) { return x.lines; }) |
                     stdv::join |
                     stdv::transform([](auto& x) { return x.f0; })) {
    const Size n = freq.size();
    freq.resize(n + nf);
    nlinspace(freq[Range(n, nf)], f0 * (1 - df), f0 * (1 + df), nf);
  }

  stdr::sort(freq);

  frequency_grid = AscendingGrid{std::move(freq)};
}
ARTS_METHOD_ERROR_CATCH

void atm_profileFitNonLTE(
    const Workspace& ws,
    ArrayOfAtmPoint& atm_profile,
    const AbsorptionBands& abs_bands,
    const Agenda& propagation_matrix_agenda,
    const SurfaceField& surface_field,
    const AscendingGrid& frequency_grid,
    const AscendingGrid& altitude_grid,
    const Numeric& latitude,
    const Numeric& longitude,
    const QuantumIdentifierGriddedField1Map& collision_data,
    const ArrayOfQuantumLevelIdentifier& levels,
    const Stokvec& pol,
    const Numeric& azimuth,
    const Numeric& dza,
    const Numeric& convergence_limit,
    const Index& iteration_limit,
    const Index& consider_limb) try {
  ARTS_TIME_REPORT

  using namespace lbl::nlte;

  ARTS_USER_ERROR_IF(
      not arr::same_size(altitude_grid, atm_profile),
      "Altitude grid and atmospheric point grid must have the same size")
  ARTS_USER_ERROR_IF(convergence_limit <= 0 or iteration_limit <= 0,
                     "Convergence limit and iteration limit must be positive")
  ARTS_USER_ERROR_IF(levels.empty(), "Need energy levels")

  const auto Aij = createAij(abs_bands);
  const auto Bij = createBij(abs_bands);
  const auto Bji = createBji(Bij, abs_bands);
  const auto Cij = createCij(abs_bands, collision_data, atm_profile);
  const auto Cji = createCji(Cij, abs_bands, atm_profile);

  const auto band_level_map = band_level_mapFromLevelKeys(abs_bands, levels);
  const Size nlevels        = level_count(band_level_map);
  const Vector r_sum        = nlte_ratio_sum(atm_profile, levels);
  const Size unique_level   = band_level_mapUniquestIndex(band_level_map);

  Matrix A;
  Matrix spectral_flux_profile;
  QuantumIdentifierVectorMap nlte_line_flux_profile;
  Vector r(nlevels, 0.0), x(nlevels, 0.0);

  int i              = 0;
  Numeric max_change = 1e99;
  while (max_change > convergence_limit and i < iteration_limit) {
    i++;
    max_change = 0.0;

    spectral_flux_profilePseudo2D(ws,
                                  spectral_flux_profile,
                                  altitude_grid,
                                  atm_profile,
                                  frequency_grid,
                                  latitude,
                                  longitude,
                                  propagation_matrix_agenda,
                                  surface_field,
                                  pol,
                                  dza,
                                  consider_limb,
                                  azimuth);

    nlte_line_flux_profileIntegrate(nlte_line_flux_profile,
                                    spectral_flux_profile,
                                    abs_bands,
                                    atm_profile,
                                    frequency_grid);

    for (Size atmi = 0; atmi < altitude_grid.size(); ++atmi) {
      A = statistical_equilibrium_equation(Aij,
                                           Bij,
                                           Bji,
                                           Cij,
                                           Cji,
                                           nlte_line_flux_profile,
                                           band_level_map,
                                           atmi,
                                           nlevels);

      A[unique_level] = 1.0;
      r[unique_level] = r_sum[atmi];

      solve(x, A, r);

      const Numeric max_change_local = set_nlte(atm_profile[atmi], levels, x);

      max_change = std::max(max_change, max_change_local);
    }
  }
}
ARTS_METHOD_ERROR_CATCH
