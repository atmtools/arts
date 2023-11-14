#include "lbl_lineshape_voigt.h"

namespace lbl {
void calculate(PropmatVectorView pm,
               matpack::matpack_view<Propmat, 2, false, true> dpm,
               const ExhaustiveConstVectorView& f_grid,
               const Jacobian::Targets& jacobian_targets,
               const std::span<const lbl::band>& bnds,
               const AtmPoint& atm,
               const Vector2 los = {},
               const bool zeeman = false);
}
