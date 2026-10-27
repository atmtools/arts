#include "functional_atm_field.h"

#include <geodetic.h>

#include <utility>

#include "functional_field_interp.h"

namespace Atm {
Numeric MagnitudeField::operator()(Numeric alt, Numeric lat, Numeric lon) const {
  const Numeric M = field_interp::interp(magnitude, {alt, lat, lon});
  const Numeric T = field_interp::interp(theta, {alt, lat, lon});
  const Numeric P = field_interp::interp(phi, {alt, lat, lon});

  const Vector3 x = geocentric2ecef({M, T, P});

  switch (component) {
    case FieldComponent::u: return x[0];
    case FieldComponent::v: return x[1];
    case FieldComponent::w: return x[2];
  }
  std::unreachable();
}

ConstVectorView MagnitudeField::x() const { return magnitude.data.view_as(magnitude.data.size()); }

VectorView MagnitudeField::x() { return magnitude.data.view_as(magnitude.data.size()); }

std::vector<std::pair<Index, Numeric>> MagnitudeField::w(Numeric alt, Numeric lat, Numeric lon) const {
  if (magnitude.grids != theta.grids or magnitude.grids != phi.grids)
    arts::user_error(
        "MagnitudeField grids do not match: {:Bs,} != {:Bs,} != {:Bs,}", magnitude.grids, theta.grids, phi.grids);

  std::vector<std::pair<Index, Numeric>> out = field_interp::flat_weight(magnitude, {alt, lat, lon});

  switch (component) {
    case FieldComponent::u:
      for (auto& [i, v] : out) v *= geocentric2ecef({1.0, theta.data.elem_at(i), phi.data.elem_at(i)})[0];
      break;
    case FieldComponent::v:
      for (auto& [i, v] : out) v *= geocentric2ecef({1.0, theta.data.elem_at(i), phi.data.elem_at(i)})[1];
      break;
    case FieldComponent::w:
      for (auto& [i, v] : out) v *= geocentric2ecef({1.0, theta.data.elem_at(i), phi.data.elem_at(i)})[2];
      break;
  }

  return out;
}
}  // namespace Atm
