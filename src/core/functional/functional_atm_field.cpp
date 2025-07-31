#include "functional_atm_field.h"

#include <geodetic.h>

#include <utility>

#include "functional_atm_field_interp.h"

namespace Atm {
Numeric MagnitudeField::operator()(Numeric alt,
                                   Numeric lat,
                                   Numeric lon) const {
  const Numeric M = Atm::interp::get(magnitude, alt, lat, lon);
  const Numeric T = Atm::interp::get(theta, alt, lat, lon);
  const Numeric P = Atm::interp::get(phi, alt, lat, lon);

  const Vector3 x = geocentric2ecef({M, T, P});

  switch (component) {
    case FieldComponent::u: return x[0];
    case FieldComponent::v: return x[1];
    case FieldComponent::w: return x[2];
  }
  std::unreachable();
}

ConstVectorView MagnitudeField::x() const {
  return magnitude.data.view_as(magnitude.data.size());
}

VectorView MagnitudeField::x() {
  return magnitude.data.view_as(magnitude.data.size());
}

std::vector<std::pair<Index, Numeric>> MagnitudeField::w(Numeric alt,
                                                         Numeric lat,
                                                         Numeric lon) const {
  ARTS_USER_ERROR_IF(
      magnitude.grids != theta.grids or magnitude.grids != phi.grids,
      "MagnitudeField grids do not match: {:Bs,} != {:Bs,} != {:Bs,}",
      magnitude.grids,
      theta.grids,
      phi.grids);

  std::vector<std::pair<Index, Numeric>> out =
      interp::flat_weight(magnitude, alt, lat, lon);

  switch (component) {
    case FieldComponent::u:
      for (auto& [i, v] : out)
        v *= geocentric2ecef(
            {1.0, theta.data.elem_at(i), phi.data.elem_at(i)})[0];
      break;
    case FieldComponent::v:
      for (auto& [i, v] : out)
        v *= geocentric2ecef(
            {1.0, theta.data.elem_at(i), phi.data.elem_at(i)})[1];
      break;
    case FieldComponent::w:
      for (auto& [i, v] : out)
        v *= geocentric2ecef(
            {1.0, theta.data.elem_at(i), phi.data.elem_at(i)})[2];
      break;
  }

  return out;
}
}  // namespace Atm
