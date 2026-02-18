#include "rtepack_stokes_vector.h"

#include <physics_funcs.h>

#include <utility>

namespace rtepack {
std::function<stokvec(const stokvec, const Numeric)> unit_converter(
    const SpectralRadianceUnitType type, const Numeric n) {
  using enum SpectralRadianceUnitType;
  switch (type) {
    case unit:
      return [n2 = n * n](const stokvec v, const Numeric) { return n2 * v; };
    case RJBT:
      return [](const stokvec v, const Numeric f) {
        return invrayjean(1.0, f) * v;
      };
    case PlanckBT:
      return [](const stokvec v, const Numeric f) {
        return stokvec{invplanck(v.I(), f),
                       invplanck(0.5 * (v.I() + v.Q()), f) -
                           invplanck(0.5 * (v.I() - v.Q()), f),
                       invplanck(0.5 * (v.I() + v.U()), f) -
                           invplanck(0.5 * (v.I() - v.U()), f),
                       invplanck(0.5 * (v.I() + v.V()), f) -
                           invplanck(0.5 * (v.I() - v.V()), f)};
      };
    case SpectralRadianceUnitType::W_m2_m_sr:
      return [n2 = n * n](const stokvec v, const Numeric f) {
        return n2 * v * (f * (f / Constant::c));
      };
    case SpectralRadianceUnitType::W_m2_m1_sr:
      return [n2 = n * n](const stokvec v, const Numeric) {
        return n2 * v * Constant::c;
      };
  }

  std::unreachable();
}

std::function<stokvec(const stokvec, const stokvec, const Numeric)>
dunit_converter(const SpectralRadianceUnitType type, const Numeric n) {
  using enum SpectralRadianceUnitType;
  switch (type) {
    case unit:
      return [n2 = n * n](const stokvec dv, const stokvec, const Numeric) {
        return n2 * dv;
      };
    case RJBT:
      return [](const stokvec dv, const stokvec, const Numeric f) {
        return invrayjean(1.0, f) * dv;
      };
    case PlanckBT:
      return [](const stokvec dv, const stokvec v, const Numeric f) {
        return stokvec{dinvplanckdI(v.I(), f) * dv.I(),
                       (dinvplanckdI(0.5 * (v.I() + v.Q()), f) -
                        dinvplanckdI(0.5 * (v.I() - v.Q()), f)) *
                           dv.Q(),
                       (dinvplanckdI(0.5 * (v.I() + v.U()), f) -
                        dinvplanckdI(0.5 * (v.I() - v.U()), f)) *
                           dv.U(),
                       (dinvplanckdI(0.5 * (v.I() + v.V()), f) -
                        dinvplanckdI(0.5 * (v.I() - v.V()), f)) *
                           dv.V()};
      };
    case SpectralRadianceUnitType::W_m2_m_sr:
      return [n2 = n * n](const stokvec dv, const stokvec, const Numeric f) {
        return n2 * dv * (f * (f / Constant::c));
      };
    case SpectralRadianceUnitType::W_m2_m1_sr:
      return [n2 = n * n](const stokvec dv, const stokvec, const Numeric) {
        return n2 * dv * Constant::c;
      };
  }

  std::unreachable();
}
}  // namespace rtepack
