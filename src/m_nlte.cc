#include <arts_omp.h>
#include <workspace.h>

#include "arts_constants.h"
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

std::unordered_map<QuantumIdentifier, UppLow> band_level_mapFromKeys(
    const AbsorptionBands& absorption_bands) try {
  std::unordered_map<QuantumIdentifier, UppLow> band_level_map;
  band_level_map.reserve(absorption_bands.size());

  std::vector<QuantumIdentifier> level_keys;
  level_keys.reserve(absorption_bands.size());

  for (const auto& key : absorption_bands | stdv::keys) {
    UppLow& ul = band_level_map[key];

    const auto lower      = key.LowerLevel();
    const auto lower_find = stdr::find(level_keys, lower);
    ul.low                = std::distance(stdr::begin(level_keys), lower_find);
    if (lower_find == stdr::end(level_keys)) level_keys.push_back(lower);

    const auto upper      = key.UpperLevel();
    const auto upper_find = stdr::find(level_keys, upper);
    ul.upp                = std::distance(stdr::begin(level_keys), upper_find);
    if (upper_find == stdr::end(level_keys)) level_keys.push_back(upper);
  }

  return band_level_map;
}
ARTS_METHOD_ERROR_CATCH

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

Vector nlte_ratio_sum(const AtmField& atmospheric_field, const Size nalt) try {
  Vector r(nalt);

  for (Size i = 0; i < nalt; i++) {
    r[i] = std::transform_reduce(
        atmospheric_field.nlte().begin(),
        atmospheric_field.nlte().end(),
        0.0,
        std::plus{},
        [i](const auto& x) {
          return std::get<GriddedField3>(x.second.data)[i, 0, 0];
        });
  }

  return r;
}
ARTS_METHOD_ERROR_CATCH

Numeric set_atmospheric_field(
    AtmField& atmospheric_field,
    const std::unordered_map<QuantumIdentifier, UppLow>& band_level_map,
    const Vector& x,
    Size atmi) try {
  std::vector<bool> changed(atmospheric_field.nlte().size(), false);

  Numeric max_change = 0.0;

  for (auto& [key, ul] : band_level_map) {
    if (not changed[ul.low]) {
      changed[ul.low] = true;
      const auto low  = key.LowerLevel();
      auto& data = std::get<GriddedField3>(atmospheric_field.nlte()[low].data);
      Numeric& a = data[atmi, 0, 0];
      const Numeric& b = x[ul.low];
      max_change       = std::max(max_change, std::abs(a - b));
      a                = b;
    }

    if (not changed[ul.upp]) {
      changed[ul.upp] = true;
      const auto upp  = key.UpperLevel();
      auto& data = std::get<GriddedField3>(atmospheric_field.nlte()[upp].data);
      Numeric& a = data[atmi, 0, 0];
      const Numeric& b = x[ul.upp];
      max_change       = std::max(max_change, std::abs(a - b));
      a                = b;
    }
  }

  return max_change;
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

void atmospheric_profileFitNonLTE(
    const Workspace& ws,
    ArrayOfAtmPoint& atmospheric_profile,
    const AbsorptionBands& absorption_bands,
    const Agenda& propagation_matrix_agenda,
    const SurfaceField& surface_field,
    const AscendingGrid& frequency_grid,
    const AscendingGrid& altitude_grid,
    const Numeric& latitude,
    const Numeric& longitude,
    const QuantumIdentifierGriddedField1Map& collision_data,
    const ArrayOfQuantumIdentifier& levels,
    const Stokvec& pol,
    const Numeric& azimuth,
    const Numeric& dza,
    const Numeric& convergence_limit,
    const Index& iteration_limit,
    const Index& consider_limb) try {
  ARTS_USER_ERROR_IF(
      not arr::same_size(altitude_grid, atmospheric_profile),
      "Altitude grid and atmospheric point grid must have the same size")
  ARTS_USER_ERROR_IF(convergence_limit <= 0 or iteration_limit <= 0,
                     "Convergence limit and iteration limit must be positive")
  ARTS_USER_ERROR_IF(levels.empty(), "Need energy levels")

  const auto Aij = createAij(absorption_bands);
  const auto Bij = createBij(absorption_bands);
  const auto Bji = createBji(Bij, absorption_bands);
  const auto Cij =
      createCij(absorption_bands, collision_data, atmospheric_profile);
  const auto Cji = createCji(Cij, absorption_bands, atmospheric_profile);

  const auto band_level_map =
      band_level_mapFromLevelKeys(absorption_bands, levels);
  const Size nlevels      = level_count(band_level_map);
  const Vector r_sum      = nlte_ratio_sum(atmospheric_profile);
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
    spectral_flux_profilePseudo2D(ws,
                                  spectral_flux_profile,
                                  altitude_grid,
                                  atmospheric_profile,
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
                                    absorption_bands,
                                    atmospheric_profile,
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
          set_atmospheric_point(atmospheric_profile[atmi], levels, x);

      max_change = std::max(max_change, max_change_local);
    }
  }
}
ARTS_METHOD_ERROR_CATCH
