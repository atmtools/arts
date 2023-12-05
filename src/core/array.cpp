#include "array.h"

#include <string_view>

namespace std {
std::ostream& operator<<(std::ostream& os, const ArrayOfIndex& x) {
  std::string_view sp = "";
  for (auto& a : x) {
    os << sp << a;
    std::exchange(sp, " ");
  }
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfIndex& x) {
  for (auto& a : x) os << a << '\n';
  return os;
}
}  // namespace std

std::ostream& operator<<(std::ostream& os, const ArrayOfNumeric& x) {
  std::string_view sp = "";
  for (auto& a : x) {
    os << sp << a;
    std::exchange(sp, " ");
  }
  return os;
}