#include "lbl_nlte.h"

#include <array_algo.h>
#include <lagrange_interp.h>
#include <matpack_mdspan_elemwise_mditer.h>
#include <partfun.h>

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

GeodeticField3 level_density(const GeodeticField3& T,
                             Numeric g,
                             Numeric E,
                             const SpeciesIsotope& spec) {
  GeodeticField3 out(T);
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
std::unordered_map<QuantumLevelIdentifier, AtmData> from_lte(
    const AtmField& atm_field,
    const AbsorptionBands& abs_bands,
    const Numeric& normalizing_factor) {
  std::unordered_map<QuantumLevelIdentifier, AtmData> out;

  ARTS_USER_ERROR_IF(not atm_field.contains(AtmKey::t),
                     "Atmospheric field does not contain temperature data");

  const AtmData& t_data = atm_field[AtmKey::t];

  for (const auto& [key, band] : abs_bands) {
    ARTS_USER_ERROR_IF(band.size() != 1, "Only for single-line bands");

    const auto upp = key.upper();
    const auto low = key.lower();

    const auto& line = band.lines.front();
    const auto spec  = key.isot;

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
      if (not normalization.contains(key.isot)) {
        normalization[key.isot] = out[key];
      } else {
        auto& x = normalization[key.isot];

        auto* srhs = std::get_if<GeodeticField3>(&data.data);
        if (auto* lhs = std::get_if<GeodeticField3>(&x.data); lhs and srhs) {
          lhs->data += srhs->data;
          continue;
        }

        auto* nrhs = std::get_if<Numeric>(&data.data);
        if (auto* lhs = std::get_if<Numeric>(&x.data); lhs and nrhs) {
          *lhs += *nrhs;
          continue;
        }

        throw std::runtime_error(
            std::format("Unsupported types or type mismatch {} and {}, "
                        "please use GeodeticField3 or Numeric for NonLTE",
                        data.data_type(),
                        x.data_type()));
      }
    }

    for (auto& [_, data] : normalization) {
      if (auto* ptr = std::get_if<GeodeticField3>(&data.data); ptr) {
        ptr->data /= normalizing_factor;
        continue;
      }

      if (auto* ptr = std::get_if<Numeric>(&data.data); ptr) {
        *ptr /= normalizing_factor;
        continue;
      }

      throw std::runtime_error(
          std::format("Unsupported type: {}\n"
                      "Please use GeodeticField3 or Numeric for NonLTE",
                      data.data_type()));
    }

    for (auto& [key, data] : out) {
      auto& x = normalization[key.isot];

      if (auto lhs = std::get_if<GeodeticField3>(&data.data),
          rhs      = std::get_if<GeodeticField3>(&x.data);
          lhs and rhs) {
        lhs->data /= rhs->data;
        continue;
      }

      if (auto lhs = std::get_if<Numeric>(&data.data),
          rhs      = std::get_if<Numeric>(&x.data);
          lhs and rhs) {
        *lhs /= *rhs;
        continue;
      }

      throw std::runtime_error(
          std::format("Unsupported types or type mismatch {} and {}, "
                      "please use GeodeticField3 or Numeric for NonLTE",
                      data.data_type(),
                      x.data_type()));
    }
  }

  return out;
}

QuantumIdentifierNumericMap createAij(const AbsorptionBands& abs_bands) try {
  QuantumIdentifierNumericMap Aij;
  Aij.reserve(abs_bands.size());

  for (const auto& [key, data] : abs_bands) {
    assert(data.size() == 1);
    Aij[key] = data.lines.front().a;
  }

  return Aij;
}
ARTS_METHOD_ERROR_CATCH

QuantumIdentifierNumericMap createBij(const AbsorptionBands& abs_bands) try {
  constexpr Numeric c0 = 2.0 * Constant::h / Math::pow2(Constant::c);

  // Size of problem
  QuantumIdentifierNumericMap Bij;
  Bij.reserve(abs_bands.size());

  // Base equation for single state:  B21 = A21 c^2 / 2 h f^3  (nb. SI, don't use this without checking your need)
  for (const auto& [key, data] : abs_bands) {
    assert(data.size() == 1);
    const auto& line = data.lines.front();
    Bij[key]         = line.a / (c0 * Math::pow3(line.f0));
  }

  return Bij;
}
ARTS_METHOD_ERROR_CATCH

QuantumIdentifierNumericMap createBji(const QuantumIdentifierNumericMap& Bij,
                                      const AbsorptionBands& abs_bands) try {
  QuantumIdentifierNumericMap Bji(Bij);

  // Base equation for single state:  B12 = B21 g2 / g1
  for (const auto& [key, data] : abs_bands) {
    assert(data.size() == 1);
    const auto& line  = data.lines.front();
    Bji[key]         *= line.gu / line.gl;
  }

  return Bji;
}
ARTS_METHOD_ERROR_CATCH

QuantumIdentifierVectorMap createCij(
    const AbsorptionBands& abs_bands,
    const QuantumIdentifierGriddedField1Map& collision_data,
    const ArrayOfAtmPoint& atm_path) try {
  QuantumIdentifierVectorMap Cij(abs_bands.size());
  for (const auto& [key, data] : abs_bands) {
    assert(data.size() == 1);
    const auto& coll = collision_data.at(key);
    Vector& x        = Cij[key];

    for (auto& atm_point : atm_path) {
      const auto numden = atm_point.number_density(key.isot);

      x.push_back(interp(coll.data, coll.lag<0, 1>(atm_point.temperature)) *
                  numden);
    }
  }

  return Cij;
}
ARTS_METHOD_ERROR_CATCH

QuantumIdentifierVectorMap createCji(const QuantumIdentifierVectorMap& Cij,
                                     const AbsorptionBands& abs_bands,
                                     const ArrayOfAtmPoint& atm_path) try {
  using Constant::h, Constant::k;

  QuantumIdentifierVectorMap Cji(Cij);
  for (const auto& [key, data] : abs_bands) {
    assert(data.size() == 1);
    const auto& line = data.lines.front();

    auto& cji = Cji.at(key);
    for (Size i = 0; i < atm_path.size(); i++) {
      const Numeric dkT  = 1.0 / (k * atm_path[i].temperature);
      cji[i]            *= std::exp(-h * line.f0 * dkT) * line.gu / line.gl;
    }
  }

  return Cji;
}
ARTS_METHOD_ERROR_CATCH

Vector nlte_ratio_sum(const ArrayOfAtmPoint& atm_path,
                      const ArrayOfQuantumLevelIdentifier& levels) try {
  return Vector(std::from_range,
                atm_path | stdv::transform([&levels](const AtmPoint& atm) {
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
    const AbsorptionBands& abs_bands,
    const ArrayOfQuantumLevelIdentifier& level_keys) try {
  std::unordered_map<QuantumIdentifier, UppLow> band_level_map;
  band_level_map.reserve(abs_bands.size());

  for (const auto& key : abs_bands | stdv::keys) {
    UppLow& ul = band_level_map[key];

    const auto lower      = key.lower();
    const auto lower_find = stdr::find(level_keys, lower);
    ul.low                = std::distance(stdr::begin(level_keys), lower_find);
    ARTS_USER_ERROR_IF(lower_find == stdr::end(level_keys),
                       "Lower level {} not found in level keys",
                       lower)

    const auto upper      = key.upper();
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

Numeric set_nlte(AtmPoint& atm_point,
                 const ArrayOfQuantumLevelIdentifier& level_keys,
                 const Vector& x) try {
  assert(x.size() == level_keys.size());
  Numeric max_change = 0.0;

  for (Size i = 0; i < level_keys.size(); i++) {
    auto& key  = level_keys[i];
    auto& v    = atm_point.nlte[key];
    max_change = std::max(max_change, std::abs(v - x[i]));
    v          = x[i];
  }

  return max_change;
}
ARTS_METHOD_ERROR_CATCH
}  // namespace lbl::nlte
