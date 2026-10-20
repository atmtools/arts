#include "functional_atm_mag_field.h"

#include <arts_conversions.h>
#include <configtypes.h>
#include <debug.h>
#include <geodetic.h>
#include <igrf.h>
#include <legendre.h>
#include <planet_data.h>

#include <tuple>

namespace Atm {
Numeric IGRF::operator()(Numeric a, Numeric la, Numeric lo) const {
  constexpr Vector2 ell{Body::Earth::a, Body::Earth::b};

  const Vector3 mag = ::IGRF::igrf({a, la, lo}, ell, time);

  switch (component) {
    case FieldComponent::u: return mag[0];
    case FieldComponent::v: return mag[1];
    case FieldComponent::w: return mag[2];
  }

  return NAN;
}

SchmidthLegendre from(const IGRF& x) {
  const auto [g, h] = ::IGRF::coefficients(x.time);

  SchmidthLegendre result{.gh        = Matrix(2, g.flat().size()),
                          .N         = static_cast<Size>(g.nrows()),
                          .r0        = ::IGRF::r0,
                          .ell       = {Body::Earth::a, Body::Earth::b},
                          .component = x.component};
  result.gh[0]  = g.flat();
  result.gh[1]  = h.flat();
  result.gh    *= 1e-9;  // Convert from nT to T
  return result;
}

std::pair<ConstTriangularMatrixView, ConstTriangularMatrixView> SchmidthLegendre::coefficients() const {
  ARTS_USER_ERROR_IF(
      gh.nrows() != 2 or static_cast<Size>(gh.ncols()) != matpack::triangular_size(N),
      "The coefficients must be 2 x {}, g and h of {} x {} lower triangular matrices, but they are {:B,}",
      matpack::triangular_size(N),
      N,
      N,
      gh.shape())

  const Index n = static_cast<Index>(N);
  return {ConstTriangularMatrixView{gh[0].data_handle(), n, n}, ConstTriangularMatrixView{gh[1].data_handle(), n, n}};
}

Numeric SchmidthLegendre::operator()(Numeric a, Numeric la, Numeric lo) const {
  const auto [g, h] = coefficients();

  const Vector3 geoc = geodetic2geocentric({a, la, lo}, ell);
  const Vector3 mag  = Legendre::schmidt_fieldcalc(g, h, r0, geoc);

  // Rotate from the geocentric to the geodetic vertical
  const Numeric ang = Conversion::deg2rad(la - geoc[1]);
  const Numeric ca  = std::cos(ang);
  const Numeric sa  = std::sin(ang);

  switch (component) {
    case FieldComponent::u: return mag[2];
    case FieldComponent::v: return -ca * mag[1] - sa * mag[0];
    case FieldComponent::w: return -sa * mag[1] + ca * mag[0];
  }

  return NAN;  // Should never happen
}

ConstVectorView SchmidthLegendre::x() const { return gh.view_as(gh.size()); }

VectorView SchmidthLegendre::x() { return gh.view_as(gh.size()); }

std::vector<std::pair<Index, Numeric>> SchmidthLegendre::w(Numeric alt, Numeric lat, Numeric lon) const try {
  using Legendre::dschmidt_fieldcalc;

  std::ignore = coefficients();

  const Vector3 geoc = geodetic2geocentric({alt, lat, lon}, ell);

  // Rotate from the geocentric to the geodetic vertical
  const Numeric ang = Conversion::deg2rad(lat - geoc[1]);
  const Numeric ca  = std::cos(ang);
  const Numeric sa  = std::sin(ang);

  const Size n = gh.ncols();

  const auto v{dschmidt_fieldcalc(N, r0, geoc)};

  assert(v.size() == 3);
  assert(v[0].first.flat().size() == n);
  assert(v[1].first.flat().size() == n);
  assert(v[2].first.flat().size() == n);
  assert(v[0].second.flat().size() == n);
  assert(v[1].second.flat().size() == n);
  assert(v[2].second.flat().size() == n);

  const auto dr  = std::pair{v[0].first.flat(), v[0].second.flat()};
  const auto dth = std::pair{v[1].first.flat(), v[1].second.flat()};
  const auto dph = std::pair{v[2].first.flat(), v[2].second.flat()};

  std::vector<std::pair<Index, Numeric>> out(gh.size());
  switch (component) {
    case FieldComponent::u:
      for (Size i = 0; i < n; i++) {
        const Numeric dBdg = dph.first[i];
        const Numeric dBdh = dph.second[i];
        out[i]             = {i, dBdg};
        out[i + n]         = {i + n, dBdh};
      }
      break;
    case FieldComponent::v:
      for (Size i = 0; i < n; i++) {
        const Numeric dBdg = -ca * dth.first[i] - sa * dr.first[i];
        const Numeric dBdh = -ca * dth.second[i] - sa * dr.second[i];
        out[i]             = {i, dBdg};
        out[i + n]         = {i + n, dBdh};
      }
      break;
    case FieldComponent::w:
      for (Size i = 0; i < n; i++) {
        const Numeric dBdg = -sa * dth.first[i] + ca * dr.first[i];
        const Numeric dBdh = -sa * dth.second[i] + ca * dr.second[i];
        out[i]             = {i, dBdg};
        out[i + n]         = {i + n, dBdh};
      }
  }

  return out;
}
ARTS_METHOD_ERROR_CATCH
}  // namespace Atm
