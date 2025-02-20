#include "lbl_nlte.h"

#include <partfun.h>

#include "matpack_mdspan_elemwise_mditer.h"

namespace {
Numeric level_density(Numeric T,
                      Numeric g,
                      Numeric E,
                      const SpeciesIsotope& spec) {
  using Constant::k;
  return g * std::exp(-E / (k * T)) / PartitionFunctions::Q(T, spec);
}

AtmFunctionalData level_density(const AtmFunctionalData& T,
                                Numeric g,
                                Numeric E,
                                const SpeciesIsotope& spec) {
  return {
      .f = [g, E, spec, T](Numeric alt, Numeric lat, Numeric lon) -> Numeric {
        return level_density(T(alt, lat, lon), g, E, spec);
      }};
}

GriddedField3 level_density(const GriddedField3& T,
                            Numeric g,
                            Numeric E,
                            const SpeciesIsotope& spec) {
  GriddedField3 out(T);
  out.data_name = "NLTE";

  stdr::transform(matpack::elemwise_range(T.data),
                  out.data.elem_begin(),
                  [g, E, spec](Numeric T) -> Numeric {
                    return level_density(T, g, E, spec);
                  });

  return out;
}
}  // namespace

namespace lbl::nlte {
std::unordered_map<QuantumIdentifier, AtmData> from_lte(
    const AtmField& atmospheric_field,
    const AbsorptionBands& absorption_bands,
    const Numeric& normalizing_factor) {
  std::unordered_map<QuantumIdentifier, AtmData> out;

  ARTS_USER_ERROR_IF(not atmospheric_field.contains(AtmKey::t),
                     "Atmospheric field does not contain temperature data");

  const AtmData& t_data = atmospheric_field[AtmKey::t];

  for (const auto& [key, band] : absorption_bands) {
    ARTS_USER_ERROR_IF(band.size() != 1, "Only for single-line bands");

    const auto upp = key.UpperLevel();
    const auto low = key.LowerLevel();

    const auto& line = band.lines.front();
    const auto spec  = key.Isotopologue();

    if (not out.contains(upp)) {
      const Numeric g = line.gu;
      const Numeric E = line.e0 + Constant::h * line.f0;
      out[upp]        = std::visit(
          [g, E, spec](auto& T) {
            return AtmData{::level_density(T, g, E, spec)};
          },
          t_data.data);
    }

    if (not out.contains(low)) {
      const Numeric g = line.gl;
      const Numeric E = line.e0;
      out[low]        = std::visit(
          [g, E, spec](auto& T) {
            return AtmData{::level_density(T, g, E, spec)};
          },
          t_data.data);
    }
  }

  if (normalizing_factor > 0.0) {
    std::unordered_map<SpeciesIsotope, AtmData> normalization;

    for (const auto& [key, data] : out) {
      if (not normalization.contains(key.Isotopologue())) {
        normalization[key.Isotopologue()] = out[key];
      } else {
        std::visit(
            []<typename T, typename U>(T& a, const U& b) {
              if constexpr (std::same_as<T, U>) {
                if constexpr (std::is_same_v<T, GriddedField3>) {
                  a.data += b.data;
                } else if constexpr (std::is_same_v<T, Numeric>) {
                  a += b;
                } else {
                  throw std::runtime_error(
                      "Unsupported type, please use GriddedField3 or Numeric for NonLTE");
                }
              } else {
                throw std::logic_error(
                    "Mismatched types - this should never happen");
              }
            },
            normalization[key.Isotopologue()].data,
            data.data);
      }
    }

    for (auto& [_, data] : normalization) {
      std::visit(
          [normalizing_factor]<typename T>(T& a) {
            if constexpr (std::is_same_v<T, GriddedField3>) {
              a.data /= normalizing_factor;
            } else if constexpr (std::is_same_v<T, Numeric>) {
              a /= normalizing_factor;
            } else {
              throw std::runtime_error(
                  "Unsupported type, please use GriddedField3 or Numeric for NonLTE");
            }
          },
          data.data);
    }

    for (auto& [key, data] : out) {
      std::visit(
          []<typename T, typename U>(T& a, const U& b) {
            if constexpr (std::same_as<T, U>) {
              if constexpr (std::is_same_v<T, GriddedField3>) {
                a.data /= b.data;
              } else if constexpr (std::is_same_v<T, Numeric>) {
                a /= b;
              } else {
                throw std::runtime_error(
                    "Unsupported type, please use GriddedField3 or Numeric for NonLTE");
              }
            } else {
              throw std::logic_error(
                  "Mismatched types - this should never happen");
            }
          },
          data.data,
          normalization[key.Isotopologue()].data);
    }
  }

  return out;
}

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

Vector nlte_ratio_sum(const ArrayOfAtmPoint& ray_path_atmospheric_point,
                      const ArrayOfQuantumIdentifier& levels) try {
  return Vector(std::from_range,
                ray_path_atmospheric_point |
                    stdv::transform([&levels](const AtmPoint& atm) {
                      Numeric s{0.0};
                      for (auto& x : levels) s += atm.nlte.at(x);
                      return s;
                    }));
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

Numeric set_nlte(AtmPoint& atmospheric_point,
                 const ArrayOfQuantumIdentifier& level_keys,
                 const Vector& x) try {
  assert(x.size() == level_keys.size());
  Numeric max_change = 0.0;

  for (Size i = 0; i < level_keys.size(); i++) {
    auto& key  = level_keys[i];
    auto& v    = atmospheric_point.nlte[key];
    max_change = std::max(max_change, std::abs(v - x[i]));
    v          = x[i];
  }

  return max_change;
}
ARTS_METHOD_ERROR_CATCH
}  // namespace lbl::nlte
