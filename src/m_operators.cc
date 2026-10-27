#include <workspace.h>

void water_equivalent_pressure_operatorMK05(NumericUnaryOperator& water_equivalent_pressure_operator,
                                            const Index&          only_liquid) {
  ARTS_TIME_REPORT

  water_equivalent_pressure_operator = NumericUnaryOperator{WaterEquivalentPressureMK05{.only_liquid = only_liquid != 0}};
}
