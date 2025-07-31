#include "functional_atm_mag_field.h"

#include <configtypes.h>
#include <debug.h>
#include <geodetic.h>
#include <igrf13.h>
#include <legendre.h>
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

SchmidthLegendre from(const IGRF13& x) {
  SchmidthLegendre result;
  const auto [g, h] = IGRF::igrf_coefficients(x.time);
  ARTS_USER_ERROR_IF(g.shape() != h.shape(),
                     "IGRF13 coefficients g and h must have the same shape, "
                     "got {} and {}",
                     g.shape(),
                     h.shape());
  result.gh.resize(2, g.nrows(), g.ncols());
  result.gh[0]      = g;
  result.gh[1]      = h;
  result.r0         = 6371.2e3;
  result.ell        = {Body::Earth::a, Body::Earth::b};
  result.component  = x.component;
  result.gh        *= 1e-9;  // Convert to nT
  return result;
}

Numeric SchmidthLegendre::operator()(Numeric a, Numeric la, Numeric lo) const {
  ARTS_USER_ERROR_IF(gh.nrows() != gh.ncols(),
                     "Bad shape, inner size must be square: hg: {:B,}",
                     gh.shape());

  using Conversion::cosd, Conversion::sind;

  const Vector3 geoc = geodetic2geocentric({a, la, lo}, ell);
  const Vector3 mag  = Legendre::schmidt_fieldcalc(gh[0], gh[1], r0, geoc);

  const Numeric ang =
      sind(la) * sind(90.0 - geoc[1]) - cosd(la) * cosd(90.0 - geoc[1]);
  const Numeric ca = std::cos(ang);
  const Numeric sa = std::sin(ang);

  switch (component) {
    case FieldComponent::u: return mag[2];
    case FieldComponent::v: return -ca * mag[1] - sa * mag[0];
    case FieldComponent::w: return -sa * mag[1] + ca * mag[0];
  }
}

ConstVectorView SchmidthLegendre::x() const { return gh.view_as(gh.size()); }

VectorView SchmidthLegendre::x() { return gh.view_as(gh.size()); }

std::vector<std::pair<Index, Numeric>> SchmidthLegendre::w(Numeric alt,
                                                           Numeric lat,
                                                           Numeric lon) const {
  using Legendre::dschmidt_fieldcalc;

  const Index c = component == FieldComponent::u   ? 0
                  : component == FieldComponent::v ? 1
                                                   : 2;

  const Size N = gh.ncols();
  const Size n = N * N;

  const Tensor3 v{dschmidt_fieldcalc(N, r0, {alt, lat, lon})[joker, c]};

  std::vector<std::pair<Index, Numeric>> out(2 * n);
  for (Size i = 0; i < n; i++) {
    out[i]     = {i, v[0, c].elem_at(i)};
    out[i + n] = {i + n, v[1, c].elem_at(i)};
  }

  return out;
}
}  // namespace Atm
