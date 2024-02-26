#include "species.h"

std::ostream& operator<<(std::ostream& os, const ArrayOfSpeciesEnum& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfSpeciesEnum& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}
