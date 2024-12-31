#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

//! FIXME: These functions should be elsewhere?
namespace Jacobian {
struct Targets;
}  // namespace Jacobian

namespace lbl {
//! NOTE: dpm and dsv are strided as input because the outer dimension is jacobian targets, however, the inner frequency dimension must be contiguous, or the code will terminate.
void calculate(PropmatVectorView pm,
               StokvecVectorView sv,
               PropmatMatrixView dpm,
               StokvecMatrixView dsv,
               const ConstVectorView f_grid,
               const Range& f_range,
               const Jacobian::Targets& jacobian_targets,
               const SpeciesEnum species,
               const AbsorptionBands& bnds,
               const linemixing::isot_map& ecs_data,
               const AtmPoint& atm,
               const Vector2 los,
               const bool no_negative_absorption);
}  // namespace lbl