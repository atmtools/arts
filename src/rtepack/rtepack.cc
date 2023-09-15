#include "rtepack.h"

namespace rtepack {
std::ostream& operator<<(std::ostream& os, const Array<propmat>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const Array<Array<propmat>>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const Array<muelmat>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const Array<Array<muelmat>>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const Array<stokvec>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const Array<Array<stokvec>>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}
}  // namespace rtepack