#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

//! FIXME: These functions should be elsewhere?
namespace Jacobian {
struct Targets;
}  // namespace Jacobian

namespace lbl {
void calculate(PropmatVectorView pm,
               StokvecVectorView sv,
               PropmatMatrixView dpm,
               StokvecMatrixView dsv,
               const ConstVectorView f_grid,
               const Range& f_range,
               const Jacobian::Targets& jac_targets,
               const SpeciesEnum species,
               const AbsorptionBands& bnds,
               const LinemixingEcsData& ecs_data,
               const AtmPoint& atm,
               const Vector2 los,
               const bool no_negative_absorption);
}  // namespace lbl