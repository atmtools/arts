#include "functional_atm_mag_field.h"

#include <arts_conversions.h>
#include <configtypes.h>
#include <debug.h>
#include <geodetic.h>
#include <igrf.h>
#include <legendre.h>
#include <planet_data.h>

#include <limits>

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

namespace {
//! The number of coefficients that affect the field component, see SchmidthLegendre::gh
Size ncoefficients(Size N, FieldComponent component) {
  return component == FieldComponent::u ? N * (N - 1) : N * N - 1;
}

/** Calls f(i, n, m, is_h) for each coefficient that affects the field component
 *
 * i is the position of the coefficient in SchmidthLegendre::gh.  The east component
 * does not depend on the m = 0 coefficients.
 */
template <typename F> void for_each_coefficient(Size N, FieldComponent component, F&& f) {
  const Index m0 = component == FieldComponent::u ? 1 : 0;

  Index i = 0;
  for (Index n = 1; n < static_cast<Index>(N); n++) {
    for (Index m = m0; m <= n; m++) {
      f(i++, n, m, false);
      if (m > 0) f(i++, n, m, true);
    }
  }
}

void check_size(const SchmidthLegendre& x) {
  ARTS_USER_ERROR_IF(x.N == 0 or static_cast<Size>(x.gh.size()) != ncoefficients(x.N, x.component),
                     "There must be {} coefficients of the {} component for N = {}, but there are {}",
                     x.N == 0 ? 0 : ncoefficients(x.N, x.component),
                     x.component,
                     x.N,
                     x.gh.size())
}
}  // namespace

SchmidthLegendre from(const IGRF& x) {
  const auto [g, h] = ::IGRF::coefficients(x.time);
  const Size N      = static_cast<Size>(g.nrows());

  SchmidthLegendre result{.gh        = Vector(ncoefficients(N, x.component)),
                          .N         = N,
                          .r0        = ::IGRF::r0,
                          .ell       = {Body::Earth::a, Body::Earth::b},
                          .component = x.component};

  // Convert from nT to T
  for_each_coefficient(
      N, x.component, [&](Index i, Index n, Index m, bool is_h) { result.gh[i] = 1e-9 * (is_h ? h : g)[n, m]; });

  return result;
}

std::pair<TriangularMatrix, TriangularMatrix> SchmidthLegendre::coefficients() const {
  check_size(*this);

  // The coefficients that do not affect the field component are NaN
  constexpr Numeric nan = std::numeric_limits<Numeric>::quiet_NaN();
  std::pair         out{TriangularMatrix(static_cast<Index>(N), nan), TriangularMatrix(static_cast<Index>(N), nan)};
  for_each_coefficient(
      N, component, [&](Index i, Index n, Index m, bool is_h) { (is_h ? out.second : out.first)[n, m] = gh[i]; });
  return out;
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

ConstVectorView SchmidthLegendre::x() const { return gh; }

VectorView SchmidthLegendre::x() { return gh; }

std::vector<std::pair<Index, Numeric>> SchmidthLegendre::w(Numeric alt, Numeric lat, Numeric lon) const try {
  check_size(*this);

  const Vector3 geoc = geodetic2geocentric({alt, lat, lon}, ell);

  // Rotate from the geocentric to the geodetic vertical
  const Numeric ang = Conversion::deg2rad(lat - geoc[1]);
  const Numeric ca  = std::cos(ang);
  const Numeric sa  = std::sin(ang);

  const auto [dr, dth, dph] = Legendre::dschmidt_fieldcalc(N, r0, geoc);

  const auto derivative = [&](Index n, Index m, bool is_h) -> Numeric {
    const TriangularMatrix& r  = is_h ? dr.second : dr.first;
    const TriangularMatrix& th = is_h ? dth.second : dth.first;
    const TriangularMatrix& ph = is_h ? dph.second : dph.first;

    switch (component) {
      case FieldComponent::u: return ph[n, m];
      case FieldComponent::v: return -ca * th[n, m] - sa * r[n, m];
      case FieldComponent::w: return -sa * th[n, m] + ca * r[n, m];
    }

    return NAN;  // Should never happen
  };

  std::vector<std::pair<Index, Numeric>> out(gh.size());
  for_each_coefficient(
      N, component, [&](Index i, Index n, Index m, bool is_h) { out[i] = {i, derivative(n, m, is_h)}; });
  return out;
}
ARTS_METHOD_ERROR_CATCH
}  // namespace Atm
