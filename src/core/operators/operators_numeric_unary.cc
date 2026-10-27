#include "operators_numeric_unary.h"

#include <arts_constants.h>

#include <cmath>

Numeric WaterEquivalentPressureMK05::operator()(Numeric t) const {
  if (only_liquid or t > Constant::temperature_at_0c) {
    return std::exp(54.842763 - 6763.22 / t - 4.21 * std::log(t) + 0.000367 * t +
                    std::tanh(0.0415 * (t - 218.8)) * (53.878 - 1331.22 / t - 9.44523 * std::log(t) + 0.014025 * t));
  }

  return std::exp(9.550426 - 5723.265 / t + 3.53068 * std::log(t) - 0.00728332 * t);
}
