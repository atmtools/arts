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
}  // namespace lbl::nlte
