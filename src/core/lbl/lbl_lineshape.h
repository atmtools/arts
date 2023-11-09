#include "lbl_lineshape_voigt.h"

namespace lbl {
  void calculate(PropmatVector& pm,
               PropmatMatrix& dpm,
               const Vector& f_grid,
               const Jacobian::Targets& jacobian_targets,
               const bands& bnds,
               const AtmPoint& atm,
               const Vector2 los={},
               const bool zeeman=false);
}
