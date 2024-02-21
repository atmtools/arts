#ifndef ARTS_CORE_PROPERTIES_H_
#define ARTS_CORE_PROPERTIES_H_

#include <array.h>

#include <algorithm>
#include <functional>
#include <string_view>

#include "enums.h"
#include "mystring.h"

/** Particulate properties
 *
 * Represents the physical quantities used to represent properties of scattering
 * particles.
 */
ENUMCLASS(ParticulateProperty,
          unsigned char,
          MassDensity,
          NumberDensity,
          DMax,
          DVeq,
          ShapeParameter,
          InterceptParameter)

/*** ScatteringSpeciesProperty
 *
 * Used to uniquely identify an atmospheric field that holds properties
 * of a given scattering species.
 *
 * */
struct ScatteringSpeciesProperty {
  std::string species_name;
  ParticulateProperty pproperty;

  constexpr auto operator<=>(const ScatteringSpeciesProperty& other) const =
      default;
};

inline std::ostream& operator<<(std::ostream& os,
                                const ScatteringSpeciesProperty& ssp) {
  return os << ssp.species_name << "_" << ssp.pproperty;
}

namespace std {
template <>
struct hash<ScatteringSpeciesProperty> {
  std::size_t operator()(const ScatteringSpeciesProperty& ssp) const {
    return std::hash<String>{}(ssp.species_name +
                               std::string(toString(ssp.pproperty)));
  }
};
}  // namespace std

#endif
