#include "functional_atm.h"

namespace Atm {
Numeric External::operator()(Numeric alt, Numeric lat, Numeric lon) const {
  return f(alt, lat, lon);
}

VectorView External::x() {
  auto s = mx();
  return VectorView{matpack::mdview_t<Numeric, 1>{s.data(), s.size()}};
}

ConstVectorView External::x() const {
  auto s = cx();
  return ConstVectorView{
      matpack::mdview_t<const Numeric, 1>{s.data(), s.size()}};
}

std::vector<std::pair<Index, Numeric>> External::w(Numeric alt,
                                                   Numeric lat,
                                                   Numeric lon) const {
  return cw(alt, lat, lon);
}
}  // namespace Atm
