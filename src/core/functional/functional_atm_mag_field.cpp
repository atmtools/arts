#include "functional_atm_mag_field.h"

#include <igrf13.h>
#include <planet_data.h>

namespace Atm {
Numeric IGRF13::operator()(Numeric a, Numeric la, Numeric lo) const {
  constexpr Vector2 ell{Body::Earth::a, Body::Earth::b};

  const Vector3 mag = IGRF::igrf({a, la, lo}, ell, time);

  switch (component) {
    case FieldComponent::u: return mag[0];
    case FieldComponent::v: return mag[1];
    case FieldComponent::w: return mag[2];
  }

  return NAN;
}
}  // namespace Atm
