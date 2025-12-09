#include <workspace.h>

void water_equivalent_pressure_operatorMK05(
    NumericUnaryOperator& water_equivalent_pressure_operator,
    const Index& only_liquid) {
  ARTS_TIME_REPORT

  if (only_liquid) {
    water_equivalent_pressure_operator = NumericUnaryOperator{[](Numeric t) {
      return std::exp(
          54.842763 - 6763.22 / t - 4.21 * std::log(t) + 0.000367 * t +
          std::tanh(0.0415 * (t - 218.8)) *
              (53.878 - 1331.22 / t - 9.44523 * std::log(t) + 0.014025 * t));
    }};
  } else {
    water_equivalent_pressure_operator = NumericUnaryOperator{[](Numeric t) {
      static constexpr Numeric TEMP_0_C = Constant::temperature_at_0c;
      if (t > TEMP_0_C) {
        return std::exp(
            54.842763 - 6763.22 / t - 4.21 * std::log(t) + 0.000367 * t +
            std::tanh(0.0415 * (t - 218.8)) *
                (53.878 - 1331.22 / t - 9.44523 * std::log(t) + 0.014025 * t));
      }
      return std::exp(9.550426 - 5723.265 / t + 3.53068 * std::log(t) -
                      0.00728332 * t);
    }};
  }
}
