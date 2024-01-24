#pragma once

#include <arts_constants.h>
#include <configtypes.h>
#include <debug.h>
#include <matpack_math.h>

#include <cmath>
#include <format>
#include <iostream>
#include <variant>

namespace Scattering {

class ScatteringHabit {};
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

using HenyeyGreenstein = Scattering::HenyeyGreenstein;

std::ostream& operator<<(
    std::ostream& os,
    const std::variant<HenyeyGreenstein,
                       Scattering::ScatteringHabit>& /*species*/);
