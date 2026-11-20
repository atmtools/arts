#pragma once

#include "lbl_data.h"
#include "lbl_lineshape_linemixing.h"

//! FIXME: These functions should be elsewhere?
namespace Jacobian {
struct Targets;
}  // namespace Jacobian

namespace lbl {
void calculate(PropmatVectorView        pm,
               StokvecVectorView        sv,
               PropmatMatrixView        dpm,
               StokvecMatrixView        dsv,
               const ConstVectorView    f_grid,
               const Range&             f_range,
               const Jacobian::Targets& jac_targets,
               const SpeciesEnum        species,
               const AbsorptionBands&   bnds,
               const LinemixingEcsData& abs_ecs_data,
               const AtmPoint&          atm,
               const Vector2            los,
               const bool               no_negative_absorption);

//! Whether calculate adds to the non-LTE source for species, which it does for a band of a non-LTE line shape
[[nodiscard]] bool adds_nlte(const AbsorptionBands& bnds, const SpeciesEnum species);
}  // namespace lbl
