#include "lbl_voigt.h"

#include <Faddeeva.hh>

bool is_voigt(LineByLineLineshape lsm) {
  using enum LineByLineLineshape;

  // If this gives you a warning, you need to update this function
  switch (lsm) {
    case VP_ECS_HARTMANN:
    case VP_ECS_MAKAROV:
    case VP_LTE:
    case VP_LINE_NLTE:
    case VP_LTE_MIRROR:   return true;
  }

  return false;
}

namespace lbl {
void compute_voigt(VectorView y,
           const line& l,
           const AscendingGrid& f_grid,
           const AtmPoint& atm,
           const Numeric mass) {
  assert(y.size() == f_grid.size());

  constexpr Numeric dop = Constant::doppler_broadening_const_squared;

  const Numeric G0    = l.ls.G0(atm);
  const Numeric D0    = l.ls.D0(atm);
  const Numeric DV    = l.ls.DV(atm);
  const Numeric Y     = l.ls.Y(atm);
  const Numeric G     = l.ls.G(atm);
  const Numeric f0    = l.f0 + D0 + DV;
  const Numeric invGD = 1.0 / (std::sqrt(dop * atm.temperature / mass) * f0);

  const Complex z0{-f0 * invGD, G0 * invGD};
  const Complex s{Constant::inv_sqrt_pi * invGD * Complex{1 + G, -Y}};

  stdr::transform(f_grid, y.begin(), [s, z0, invGD](Numeric f) {
    return std::real(s * Faddeeva::w(f * invGD + z0));
  });
}
}  // namespace lbl
