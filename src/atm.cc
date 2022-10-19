#include "atm.h"

#include <algorithm>
#include <ostream>
#include <variant>

#include "gridded_fields.h"
#include "matpackIV.h"

namespace Atm {
std::ostream& operator<<(std::ostream& os, const Point& pnt) {
  os << "Temperature: " << pnt.temperature << " K,\n";
  os << "Pressure: " << pnt.pressure << " Pa,\n";
  os << "Wind Field: [u: " << pnt.wind[0] << ", v: " << pnt.wind[1] << ", w: "
     << pnt.wind[2] << "] m/s,\n";
  os << "Magnetic Field: [u: " << pnt.mag[0] << ", v: " << pnt.mag[1] << ", w: "
     << pnt.mag[2] << "] T";
  for (auto& spec : pnt.species_content) {
    os << ",\n" << spec.first << ": " << spec.second;
  }
  return os;
}
}  // namespace Atm
