#pragma once

#include <atm.h>
#include <jacobian.h>
#include <predef_data.h>
#include <rtepack.h>

namespace Absorption::PredefinedModel {

//! Returns true if the model can be computed
bool can_compute(const SpeciesIsotope& model);

/** Compute the predefined model
 *
 * The tag is checked, so this should just be looped over by all available species
 * 
 * @param[inout] propmat_clearsky As WSV
 * @param[inout] dpropmat_clearsky_dx As WSV
 * @param[in] tag An isotope record
 * @param[in] f_grid As WSV
 * @param[in] atm_point An atmospheric point object
 * @param[in] jacobian_targets As WSV
 * @param[in] predefined_model_data As WSV
 */
void compute(
    PropmatVector& propmat_clearsky,
    PropmatMatrix& dpropmat_clearsky_dx,
    const SpeciesIsotope& tag,
    const Vector& f_grid,
    const AtmPoint& atm_point,
    const JacobianTargets& jacobian_targets,
    const Absorption::PredefinedModel::ModelVariant& predefined_model_data);
}  // namespace Absorption::PredefinedModel
