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
#include "psd.h"

namespace Scattering {

using PSD = std::variant<MGDSingleMoment>;

/*** A particle habit
 *
 * A particle habit represents a collection of scattering particles each with
 * a given size, mass, and specific scattering properties.
 */
class ParticleHabit {
 public:
  ParticleHabit(){};
  ParticleHabit(std::string scattering_data){};
};

/*** A scattering habit
 *
 * A scattering habit combines a particle habit with an additional PSD
 * and thus defines a mapping between atmospheric scattering species properties
 * and corresponding bulk skattering properties.
 */
class ScatteringHabit {
 public:
  ScatteringHabit(){};
  ScatteringHabit(ParticleHabit particle_habit_, PSD psd_)
      : particle_habit(particle_habit_), psd(psd_){};

 private:
  ParticleHabit particle_habit;
  PSD psd;
};

class HenyeyGreenstein {
  Numeric g = 0.0;

 public:
  HenyeyGreenstein(){};
  HenyeyGreenstein(const Numeric& g_);
  Numeric evaluate_phase_function(const Numeric& theta);
  Vector evaluate_phase_function(const Vector& theta);

  Numeric get_g() const { return g; };
  void set_g(const Numeric& g_) { g = g_; };

  friend std::ostream& operator<<(std::ostream& os,
                                  const HenyeyGreenstein& scatterer);
};

using Species = std::variant<HenyeyGreenstein, ScatteringHabit>;

}  // namespace Scattering

using ScatteringSpecies = Scattering::Species;

class ArrayOfScatteringSpecies : std::vector<Scattering::Species> {
 public:
  void add(const Scattering::Species& species) { push_back(species); }
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
  FmtContext::iterator format(const ArrayOfScatteringSpecies& v,
                              FmtContext& ctx) const {

    return ctx.out();
  }
};

using HenyeyGreenstein = Scattering::HenyeyGreenstein;
using ParticleHabit = Scattering::ParticleHabit;
using ScatteringHabit = Scattering::ScatteringHabit;
using PSD = Scattering::PSD;

std::ostream& operator<<(
    std::ostream& os,
    const std::variant<HenyeyGreenstein,
                       Scattering::ScatteringHabit>& /*species*/);
