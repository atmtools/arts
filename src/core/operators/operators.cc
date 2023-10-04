#include "operators.h"

std::ostream& operator<<(std::ostream& os, const NumericUnaryOperator&) {
  return os << "Numeric := op(Numeric)";
}

std::ostream& operator<<(std::ostream& os, const NumericTernaryOperator&) {
  return os << "Numeric := op(Numeric, Numeric, Numeric)";
}
