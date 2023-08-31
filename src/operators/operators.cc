#include "operators.h"

std::ostream& operator<<(std::ostream& os, const NumericUnaryOperator&) {
  return os << "Numeric := op(Numeric)";
}
