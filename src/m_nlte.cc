#include <workspace.h>

#include "auto_wsm.h"
#include "jacobian.h"

void absorption_bandsSetNonLTE(AbsorptionBands& absorption_bands) {
  for (auto& [_, band] : absorption_bands) {
    band.lineshape = LineByLineLineshape::VP_LINE_NLTE;
  }
}

void atmospheric_fieldInitializeNonLTE(AtmField& atmospheric_field,
                                       const AbsorptionBands& absorption_bands,
                                       const Numeric& normalizing_factor) try {
  atmospheric_field.nlte() = lbl::nlte::from_lte(
      atmospheric_field, absorption_bands, normalizing_factor);
}
ARTS_METHOD_ERROR_CATCH

QuantumIdentifierNumericMap createAij(
    const AbsorptionBands& absorption_bands) try {
  QuantumIdentifierNumericMap Aij;
  Aij.reserve(absorption_bands.size());

  for (const auto& [key, data] : absorption_bands) {
    assert(data.size() == 1);
    Aij[key] = data.lines.front().a;
  }

  return Aij;
}
ARTS_METHOD_ERROR_CATCH

QuantumIdentifierNumericMap createBij(
    const AbsorptionBands& absorption_bands) try {
  constexpr Numeric c0 = 2.0 * Constant::h / Math::pow2(Constant::c);

  // Size of problem
  QuantumIdentifierNumericMap Bij;
  Bij.reserve(absorption_bands.size());

  // Base equation for single state:  B21 = A21 c^2 / 2 h f^3  (nb. SI, don't use this without checking your need)
  for (const auto& [key, data] : absorption_bands) {
    assert(data.size() == 1);
    const auto& line = data.lines.front();
    Bij[key]         = line.a / (c0 * Math::pow3(line.f0));
  }

  return Bij;
}
ARTS_METHOD_ERROR_CATCH

QuantumIdentifierNumericMap createBji(
    const QuantumIdentifierNumericMap& Bij,
    const AbsorptionBands& absorption_bands) try {
  QuantumIdentifierNumericMap Bji(Bij);

  // Base equation for single state:  B12 = B21 g2 / g1
  for (const auto& [key, data] : absorption_bands) {
    assert(data.size() == 1);
    const auto& line  = data.lines.front();
    Bji[key]         *= line.gu / line.gl;
  }

  return Bji;
}
ARTS_METHOD_ERROR_CATCH

QuantumIdentifierVectorMap createCij(
    const AbsorptionBands& absorption_bands,
    const QuantumIdentifierGriddedField1Map& collision_data,
    const ArrayOfAtmPoint& ray_path_atmospheric_point) try {
  using lag = FixedLagrangeInterpolation<1>;

  QuantumIdentifierVectorMap Cij(absorption_bands.size());
  for (const auto& [key, data] : absorption_bands) {
    assert(data.size() == 1);
    const auto& coll = collision_data.at(key);
    Vector& x        = Cij[key];

    for (auto& atmospheric_point : ray_path_atmospheric_point) {
      const auto numden = atmospheric_point.number_density(key.Isotopologue());
      x.push_back(coll.interp<lag>(atmospheric_point.temperature) * numden);
    }
  }

  return Cij;
}
ARTS_METHOD_ERROR_CATCH

QuantumIdentifierVectorMap createCji(
    const QuantumIdentifierVectorMap& Cij,
    const AbsorptionBands& absorption_bands,
    const ArrayOfAtmPoint& ray_path_atmospheric_point) try {
  using Constant::h, Constant::k;

  QuantumIdentifierVectorMap Cji(Cij);
  for (const auto& [key, data] : absorption_bands) {
    assert(data.size() == 1);
    const auto& line = data.lines.front();

    auto& cji = Cji.at(key);
    for (Size i = 0; i < ray_path_atmospheric_point.size(); i++) {
      const Numeric dkT = 1.0 / (k * ray_path_atmospheric_point[i].temperature);
      cji[i] *= std::exp(-h * line.f0 * dkT) * line.gu / line.gl;
    }
  }

  return Cji;
}
ARTS_METHOD_ERROR_CATCH

const GriddedField3::grids_t& nlte_grid(const AtmField& atmospheric_field) try {
  ARTS_USER_ERROR_IF(atmospheric_field.nlte().empty(), "No level data provided")

  const auto beg         = atmospheric_field.nlte().begin();
  const auto& ld0        = beg->first;
  const auto& front_data = beg->second.data;
  ARTS_USER_ERROR_IF(not std::holds_alternative<GriddedField3>(front_data),
                     "The first \"{}\" level data is not a GriddedField3",
                     ld0)
  const auto& gf3   = std::get<GriddedField3>(front_data);
  const Vector& alt = gf3.grid<0>();
  const Vector& lat = gf3.grid<1>();
  const Vector& lon = gf3.grid<2>();
  ARTS_USER_ERROR_IF(lat.size() != 1 or lon.size() != 1,
                     "Latitude and longitude data must be 1D")

  for (auto& [ld, nlte_data] : atmospheric_field.nlte() | stdv::drop(1)) {
    ARTS_USER_ERROR_IF(not atmospheric_field.has(ld),
                       "Level data \"{}\" not found in atmospheric field",
                       ld)
    const auto& data = nlte_data.data;
    ARTS_USER_ERROR_IF(not std::holds_alternative<GriddedField3>(data),
                       "Level data \"{}\" not a GriddedField3 in atmosphere",
                       ld)
    const auto& gf3i   = std::get<GriddedField3>(data);
    const Vector& alti = gf3i.template grid<0>();
    const Vector& lati = gf3i.template grid<1>();
    const Vector& loni = gf3i.template grid<2>();
    ARTS_USER_ERROR_IF(alt != alti or lat != lati or lon != loni,
                       R"(Level data grids do not match

  Altitude[~0]:  {:B,}
  Altitude[~i]:  {:B,}

  Latitude[~0]:  {:B,}
  Latitude[~i]:  {:B,}

  Longitude[~0]: {:B,}
  Longitude[~i]: {:B,}

  level[~0]:     {}
  level[~i]:     {}
)",
                       alt,
                       alti,
                       lat,
                       lati,
                       lon,
                       loni,
                       ld0,
                       ld)
  }

  return gf3.grids;
}
ARTS_METHOD_ERROR_CATCH

ArrayOfAtmPoint extract(const AtmField& atmospheric_field,
                        const Vector& altitude_grid,
                        const Vector& latitude_grid,
                        const Vector& longitude_grid) try {
  ArrayOfAtmPoint atm_point;
  atm_point.reserve(altitude_grid.size() * latitude_grid.size() *
                    longitude_grid.size());

  for (const auto& alt : altitude_grid) {
    for (const auto& lat : latitude_grid) {
      for (const auto& lon : longitude_grid) {
        atm_point.push_back(atmospheric_field.at(alt, lat, lon));
      }
    }
  }

  return atm_point;
}
ARTS_METHOD_ERROR_CATCH

struct UppLow {
  Size upp, low;
};

std::unordered_map<QuantumIdentifier, UppLow> band_level_mapFromLevelKeys(
    const AbsorptionBands& absorption_bands,
    const ArrayOfQuantumIdentifier& level_keys) try {
  std::unordered_map<QuantumIdentifier, UppLow> band_level_map;
  band_level_map.reserve(absorption_bands.size());

  for (const auto& key : absorption_bands | stdv::keys) {
    UppLow& ul = band_level_map[key];

    const auto lower      = key.LowerLevel();
    const auto lower_find = stdr::find(level_keys, lower);
    ul.low                = std::distance(stdr::begin(level_keys), lower_find);
    ARTS_USER_ERROR_IF(lower_find == stdr::end(level_keys),
                       "Lower level {} not found in level keys",
                       lower)

    const auto upper      = key.UpperLevel();
    const auto upper_find = stdr::find(level_keys, upper);
    ul.upp                = std::distance(stdr::begin(level_keys), upper_find);
    ARTS_USER_ERROR_IF(upper_find == stdr::end(level_keys),
                       "Upper level {} not found in level keys",
                       upper)
  }

  return band_level_map;
}
ARTS_METHOD_ERROR_CATCH

Size level_count(
    const std::unordered_map<QuantumIdentifier, UppLow>& band_level_map) try {
  Size nlevels = 0;

  for (auto [i, j] : band_level_map | stdv::values) {
    nlevels = std::max({nlevels, i, j});
  }

  return nlevels + 1;
}
ARTS_METHOD_ERROR_CATCH

Matrix statistical_equilibrium_equation(
    const QuantumIdentifierNumericMap& Aij,
    const QuantumIdentifierNumericMap& Bij,
    const QuantumIdentifierNumericMap& Bji,
    const QuantumIdentifierVectorMap& Cij,
    const QuantumIdentifierVectorMap& Cji,
    const QuantumIdentifierVectorMap& Jij,
    const std::unordered_map<QuantumIdentifier, UppLow>& level_map,
    const Size atmi,
    const Size nlevels) try {
  assert(arr::same_size(Aij, Bij, Bji, Cij, Jij, level_map));

  Matrix A(nlevels, nlevels, 0.0);
  for (const auto& [key, ul] : level_map) {
    const auto i = ul.upp;
    const auto j = ul.low;

    assert(i < nlevels and j < nlevels);

    const auto& aij = Aij.at(key);
    const auto& bij = Bij.at(key);
    const auto& bji = Bji.at(key);
    const auto& cij = Cij.at(key)[atmi];
    const auto& cji = Cji.at(key)[atmi];
    const auto jij  = Jij.at(key)[atmi] * Constant::inv_two_pi;

    A[j, j] -= bji * jij + cji;
    A[i, i] -= aij + bij * jij + cij;

    A[j, i] += aij + bij * jij + cij;
    A[i, j] += bji * jij + cji;
  }

  return A;
}
ARTS_METHOD_ERROR_CATCH

Vector nlte_ratio_sum(const ArrayOfAtmPoint& ray_path_atmospheric_point) try {
  return Vector(
      std::from_range,
      ray_path_atmospheric_point | stdv::transform([](const AtmPoint& atm) {
        Numeric s{0.0};
        for (Numeric x : atm.nlte | stdv::values) s += x;
        return s;
      }));
}
ARTS_METHOD_ERROR_CATCH

Numeric set_atmospheric_field(AtmField& atmospheric_field,
                              const ArrayOfQuantumIdentifier& level_keys,
                              const Vector& x,
                              Size atmi) try {
  assert(x.size() == level_keys.size());
  Numeric max_change = 0.0;

  for (Size i = 0; i < level_keys.size(); i++) {
    auto& key = level_keys[i];
    auto& v = std::get<GriddedField3>(atmospheric_field[key].data)[atmi, 0, 0];
    max_change = std::max(max_change, std::abs(v - x[i]));
    v          = x[i];
  }

  return max_change;
}
ARTS_METHOD_ERROR_CATCH

Numeric set_atmospheric_point(AtmPoint& atmospheric_point,
                              const ArrayOfQuantumIdentifier& level_keys,
                              const Vector& x) try {
  assert(x.size() == level_keys.size());
  Numeric max_change = 0.0;

  for (Size i = 0; i < level_keys.size(); i++) {
    auto& key  = level_keys[i];
    auto& v    = atmospheric_point[key];
    max_change = std::max(max_change, std::abs(v - x[i]));
    v          = x[i];
  }

  return max_change;
}
ARTS_METHOD_ERROR_CATCH

Size band_level_mapUniquestIndex(
    const std::unordered_map<QuantumIdentifier, UppLow>& band_level_map) try {
  std::unordered_map<Size, Size> unique_level_map;
  for (auto& [upp, low] : band_level_map | stdv::values) {
    unique_level_map[upp]++;
    unique_level_map[low]++;
  }

  return stdr::min_element(
             unique_level_map, {}, [](auto& x) { return x.second; })
      ->first;
}
ARTS_METHOD_ERROR_CATCH

void frequency_gridFitNonLTE(AscendingGrid& frequency_grid,
                             const AbsorptionBands& absorption_bands,
                             const Numeric& df,
                             const Index& nf) try {
  Vector freq;
  freq.reserve(count_lines(absorption_bands) * nf);
  for (auto f0 : absorption_bands | stdv::values |
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

void atmospheric_fieldFitNonLTE(
    const Workspace& ws,
    AtmField& atmospheric_field,
    const AbsorptionBands& absorption_bands,
    const ArrayOfArrayOfPropagationPathPoint& ray_path_field,
    const Agenda& propagation_matrix_agenda,
    const Agenda& spectral_radiance_space_agenda,
    const Agenda& spectral_radiance_surface_agenda,
    const SurfaceField& surface_field,
    const AscendingGrid& frequency_grid,
    const QuantumIdentifierGriddedField1Map& collision_data,
    const ArrayOfQuantumIdentifier& levels,
    const Numeric& convergence_limit,
    const Index& iteration_limit) try {
  ARTS_USER_ERROR_IF(convergence_limit <= 0 or iteration_limit <= 0,
                     "Convergence limit and iteration limit must be positive")
  ARTS_USER_ERROR_IF(levels.empty(), "Need energy levels")

  const auto& grid{nlte_grid(atmospheric_field)};
  const AscendingGrid altitude_grid{std::get<0>(grid)};
  const ArrayOfAtmPoint ray_path_atmospheric_point = extract(
      atmospheric_field, altitude_grid, std::get<1>(grid), std::get<2>(grid));

  const auto Aij = createAij(absorption_bands);
  const auto Bij = createBij(absorption_bands);
  const auto Bji = createBji(Bij, absorption_bands);
  const auto Cij =
      createCij(absorption_bands, collision_data, ray_path_atmospheric_point);
  const auto Cji = createCji(Cij, absorption_bands, ray_path_atmospheric_point);

  const auto band_level_map =
      band_level_mapFromLevelKeys(absorption_bands, levels);
  const Size nlevels      = level_count(band_level_map);
  const Vector r_sum      = nlte_ratio_sum(ray_path_atmospheric_point);
  const Size unique_level = band_level_mapUniquestIndex(band_level_map);

  Matrix A;
  Matrix spectral_flux_profile;
  QuantumIdentifierVectorMap nlte_line_flux_profile;
  Vector r(nlevels, 0.0), x(nlevels, 0.0);

  int i              = 0;
  Numeric max_change = 1e99;
  while (max_change > convergence_limit and i < iteration_limit) {
    i++;
    max_change = 0.0;

    spectral_flux_profileFromPathField(ws,
                                       spectral_flux_profile,
                                       ray_path_field,
                                       atmospheric_field,
                                       propagation_matrix_agenda,
                                       spectral_radiance_space_agenda,
                                       spectral_radiance_surface_agenda,
                                       surface_field,
                                       frequency_grid,
                                       altitude_grid);

    nlte_line_flux_profileIntegrate(nlte_line_flux_profile,
                                    spectral_flux_profile,
                                    absorption_bands,
                                    ray_path_atmospheric_point,
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

      const Numeric max_change_local =
          set_atmospheric_field(atmospheric_field, levels, x, atmi);
      max_change = std::max(max_change, max_change_local);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void atmospheric_fieldDisortFitNonLTE(
    const Workspace& ws,
    AtmField& atmospheric_field,
    DisortSettings& disort_settings,
    const AbsorptionBands& absorption_bands,
    const ArrayOfPropagationPathPoint& ray_path,
    const Agenda& propagation_matrix_agenda,
    const AscendingGrid& frequency_grid,
    const QuantumIdentifierGriddedField1Map& collision_data,
    const ArrayOfQuantumIdentifier& levels,
    const Numeric& convergence_limit,
    const Index& iteration_limit) try {
  ARTS_USER_ERROR_IF(convergence_limit <= 0 or iteration_limit <= 0,
                     "Convergence limit and iteration limit must be positive")
  ARTS_USER_ERROR_IF(levels.empty(), "Need energy levels")

  const auto& grid{nlte_grid(atmospheric_field)};
  const AscendingGrid altitude_grid{std::get<0>(grid)};
  ArrayOfAtmPoint ray_path_atmospheric_point = extract(
      atmospheric_field, altitude_grid, std::get<1>(grid), std::get<2>(grid));

  const auto Aij = createAij(absorption_bands);
  const auto Bij = createBij(absorption_bands);
  const auto Bji = createBji(Bij, absorption_bands);
  const auto Cij =
      createCij(absorption_bands, collision_data, ray_path_atmospheric_point);
  const auto Cji = createCji(Cij, absorption_bands, ray_path_atmospheric_point);

  const auto band_level_map =
      band_level_mapFromLevelKeys(absorption_bands, levels);
  const Size nlevels      = level_count(band_level_map);
  const Vector r_sum      = nlte_ratio_sum(ray_path_atmospheric_point);
  const Size unique_level = band_level_mapUniquestIndex(band_level_map);

  ArrayOfPropmatVector ray_path_propagation_matrix{};
  ArrayOfStokvecVector ray_path_propagation_matrix_source_vector_nonlte{};
  ArrayOfPropmatMatrix d0{};
  ArrayOfStokvecMatrix d1{};
  ArrayOfAscendingGrid ray_path_frequency_grid(ray_path.size(), frequency_grid);
  ArrayOfVector3 ray_path_frequency_grid_wind_shift_jacobian(ray_path.size());
  Tensor3 disort_spectral_flux_field{};
  Matrix spectral_flux_profile(ray_path_atmospheric_point.size(),
                               frequency_grid.size());
  QuantumIdentifierVectorMap nlte_line_flux_profile{};
  JacobianTargets jacobian_targets{};
  Matrix A{};
  Vector r(nlevels, 0.0), x(nlevels, 0.0);

  int i              = 0;
  Numeric max_change = 1e99;
  while (max_change > convergence_limit and i < iteration_limit) {
    max_change = 0.0;

    if (i++ != 0) {
      ray_path_propagation_matrixFromPath(
          ws,
          ray_path_propagation_matrix,
          ray_path_propagation_matrix_source_vector_nonlte,
          d0,
          d1,
          propagation_matrix_agenda,
          ray_path_frequency_grid,
          ray_path_frequency_grid_wind_shift_jacobian,
          jacobian_targets,
          ray_path,
          ray_path_atmospheric_point);

      disort_settingsOpticalThicknessFromPath(
          disort_settings, ray_path, ray_path_propagation_matrix, 1);

      disort_settingsLayerNonThermalEmissionLinearInTau(
          disort_settings,
          ray_path_atmospheric_point,
          ray_path_propagation_matrix,
          ray_path_propagation_matrix_source_vector_nonlte,
          frequency_grid);
    }

    disort_spectral_flux_fieldCalc(disort_spectral_flux_field, disort_settings);

std::println("{:B,}", disort_spectral_flux_field[200, joker, joker]);
    StridedMatrixView v = transpose(spectral_flux_profile);
    v  = transpose(disort_spectral_flux_field[joker, 0, joker]);
    v += transpose(disort_spectral_flux_field[joker, 1, joker]);
    v += transpose(disort_spectral_flux_field[joker, 2, joker]);

    nlte_line_flux_profileIntegrate(nlte_line_flux_profile,
                                    spectral_flux_profile,
                                    absorption_bands,
                                    ray_path_atmospheric_point,
                                    frequency_grid);

                                    for (auto x:nlte_line_flux_profile|stdv::values){
                                      std::println("a.append({:B,})", x);
                                    }std::println("");

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

      const Numeric max_change_local =
          set_atmospheric_point(ray_path_atmospheric_point[atmi], levels, x);
      max_change = std::max(max_change, max_change_local);
    }
  }

  const AtmField new_atm = atm_from_profile(ray_path_atmospheric_point,
                                            altitude_grid,
                                            InterpolationExtrapolation::Linear);
  for (auto& [key, data] : new_atm.nlte()) {
    atmospheric_field[key].data = data.data;
  }
}
ARTS_METHOD_ERROR_CATCH