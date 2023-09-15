#include "array.h"

namespace std {
std::ostream& operator<<(std::ostream& os, const ArrayOfIndex& x) {
  for (auto& a : x) os << a << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os,
                                const ArrayOfArrayOfIndex& x) {
  for (auto& a : x) os << a << '\n';
  return os;
}
}  // namespace std

std::ostream& operator<<(std::ostream& os, const ArrayOfNumeric& x) {
  for (auto& a : x) os << a << '\n';
  return os;
}