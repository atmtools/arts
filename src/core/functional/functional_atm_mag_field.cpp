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

  result.gh.resize(2, g.x.size());
  result.N          = g.N;
  result.gh[0]      = g.x;
  result.gh[1]      = h.x;
  result.r0         = 6371.2e3;
  result.ell        = {Body::Earth::a, Body::Earth::b};
  result.component  = x.component;
  result.gh        *= 1e-9;  // Convert to nT
  std::println("{:B}", result.gh);
  return result;
}

Numeric SchmidthLegendre::operator()(Numeric a, Numeric la, Numeric lo) const {
  using Conversion::cosd, Conversion::sind;

  const Vector3 geoc = geodetic2geocentric({a, la, lo}, ell);
  const Vector3 mag =
      Legendre::schmidt_fieldcalc({N, gh[0]}, {N, gh[1]}, r0, geoc);

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
                                                           Numeric lon) const
    try {
  using Legendre::dschmidt_fieldcalc;

  using Conversion::cosd, Conversion::sind;

  const Vector3 geoc = geodetic2geocentric({alt, lat, lon}, ell);

  const Numeric ang =
      sind(lat) * sind(90.0 - geoc[1]) - cosd(lat) * cosd(90.0 - geoc[1]);
  const Numeric ca = std::cos(ang);
  const Numeric sa = std::sin(ang);

  const Size n = gh.ncols();

  const auto v{dschmidt_fieldcalc(N, r0, {alt, lat, lon})};

  assert(v.size() == 3);
  assert(v[0].first.x.size() == n);
  assert(v[1].first.x.size() == n);
  assert(v[2].first.x.size() == n);
  assert(v[0].second.x.size() == n);
  assert(v[1].second.x.size() == n);
  assert(v[2].second.x.size() == n);

  std::vector<std::pair<Index, Numeric>> out(gh.size());
  switch (component) {
    case FieldComponent::u:
      for (Size i = 0; i < n; i++) {
        const Numeric dBdg = v[2].first.x[i];
        const Numeric dBdh = v[2].second.x[i];
        out[i]             = {i, dBdg};
        out[i + n]         = {i + n, dBdh};
      }
      break;
    case FieldComponent::v:
      for (Size i = 0; i < n; i++) {
        const Numeric dBdg = -ca * v[1].first.x[i] - sa * v[0].first.x[i];
        const Numeric dBdh = -ca * v[1].second.x[i] - sa * v[0].second.x[i];
        out[i]             = {i, dBdg};
        out[i + n]         = {i + n, dBdh};
      }
      break;
    case FieldComponent::w:
      for (Size i = 0; i < n; i++) {
        const Numeric dBdg = -sa * v[1].first.x[i] + ca * v[0].first.x[i];
        const Numeric dBdh = -sa * v[1].second.x[i] + ca * v[0].second.x[i];
        out[i]             = {i, dBdg};
        out[i + n]         = {i + n, dBdh};
      }
  }

  return out;
}
ARTS_METHOD_ERROR_CATCH
}  // namespace Atm
