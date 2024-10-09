#pragma once

#include <arts_constants.h>
#include <configtypes.h>
#include <debug.h>
#include <matpack_math.h>

#include <cmath>
#include <format>
#include <iostream>
#include <variant>

#include "properties.h"
#include "henyey_greenstein.h"
#include "psd.h"
#include "particle_habit.h"

namespace scattering {

using PSD = std::variant<MGDSingleMoment>;


/*** A scattering habit
 *
 * A scattering habit combines a particle habit with an additional PSD
 * and thus defines a mapping between atmospheric scattering species properties
 * and corresponding bulk skattering properties.
 */
class ScatteringHabit {
 public:
  ScatteringHabit(){};

 private:
  ParticleHabit particle_habit;
  PSD psd;
};


using Species = std::variant<HenyeyGreensteinScatterer, ScatteringHabit>;

}  // namespace scattering

using ScatteringSpecies = scattering::Species;

class ArrayOfScatteringSpecies : std::vector<scattering::Species> {
 public:
  void add(const scattering::Species& species) { push_back(species); }
};

inline std::ostream& operator<<(std::ostream& os,
                                const ArrayOfScatteringSpecies& /*species*/) {
  os << "An array of scattering species." << std::endl;
  return os;
}

template<>
struct std::formatter<ArrayOfScatteringSpecies> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const ArrayOfScatteringSpecies&,
                              FmtContext& ctx) const {

    return ctx.out();
  }
};

using HenyeyGreensteinScatterer = scattering::HenyeyGreensteinScatterer;
using ParticleHabit = scattering::ParticleHabit;
using ScatteringHabit = scattering::ScatteringHabit;
using PSD = scattering::PSD;

std::ostream& operator<<(
    std::ostream& os,
    const std::variant<HenyeyGreensteinScatterer,
                       scattering::ScatteringHabit>& /*species*/);
