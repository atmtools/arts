#include "arts_conversions.h"

#include <array>

#include "nonstd.h"

/** Namespace containing several practical unit conversions, physical and mathematical **/
namespace Conversion {
std::pair<char, Numeric> metric_prefix(Numeric x) {
  if (x == 0.0) return {' ', 0.0};

  constexpr std::array<std::pair<char, Numeric>, 11> high_units = {
      {{'Q', 1e30},
       {'R', 1e27},
       {'Y', 1e24},
       {'Z', 1e21},
       {'E', 1e18},
       {'P', 1e15},
       {'T', 1e12},
       {'G', 1e9},
       {'M', 1e6},
       {'k', 1e3},
       {'h', 1e2}}};

  for (const auto& [unit, factor] : high_units) {
    if (nonstd::abs(x) > factor) return {unit, x / factor};
  }

  constexpr std::array<std::pair<char, Numeric>, 12> low_units = {
      {{'q', 1e-30 * 1000},
       {'r', 1e-27 * 1000},
       {'y', 1e-24 * 1000},
       {'z', 1e-21 * 1000},
       {'a', 1e-18 * 1000},
       {'f', 1e-15 * 1000},
       {'p', 1e-12 * 1000},
       {'n', 1e-9 * 1000},
       {'u', 1e-6 * 1000},
       {'m', 1e-3 * 1000},
       {'c', 1e-2 * 1000},
       {'d', 1e-1 * 1000}}};

  for (const auto& [unit, factor] : low_units) {
    if (nonstd::abs(x) < factor) return {unit, 1000.0 * x / factor};
  }

  return {' ', 0.0};
}
};  // namespace Conversion
