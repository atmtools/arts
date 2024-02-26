/*!
 * @file   fullmodel.h
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */

#ifndef fullmodel_h
#define fullmodel_h

#include <rtepack.h>

#include "atm.h"
#include "jacobian.h"
#include "predefined/predef_data.h"
#include "species.h"

namespace Absorption::PredefinedModel {
/** Contains known required VMR values
 *
 *  Note: If you add a species here, add it to the Jacobian
 *  wrapper in the compute function.
 */
struct VMRS {
  static constexpr std::array species = {SpeciesEnum::CarbonDioxide,
                                         SpeciesEnum::Oxygen,
                                         SpeciesEnum::Nitrogen,
                                         SpeciesEnum::Water,
                                         SpeciesEnum::liquidcloud};

  Numeric CO2{0};
  Numeric O2{0};
  Numeric N2{0};
  Numeric H2O{0};
  Numeric LWC{0};

  /**  Construct a new VMRS object
   *
   * @param atm An Atmosphere
   */
  VMRS(const AtmPoint& atm)
      : CO2(atm[species[0]]),
        O2(atm[species[1]]),
        N2(atm[species[2]]),
        H2O(atm[species[3]]),
        LWC(atm[species[4]]) {}

  constexpr VMRS() = default;

  friend std::ostream& operator<<(std::ostream& os, const VMRS& vmrs);
};

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
 * @param[in] rtp_pressure As WSV
 * @param[in] rtp_temperature As WSV
 * @param[in] vmr The VMRS defined from WSVs abs_species and rtp_vmr
 * @param[in] jacobian_targets As WSV
 * @param[in] predefined_model_data As WSV
 */
void compute(
    PropmatVector& propmat_clearsky,
    PropmatMatrix& dpropmat_clearsky_dx,
    const SpeciesIsotope& tag,
    const Vector& f_grid,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const VMRS& vmr,
    const JacobianTargets& jacobian_targets,
    const Absorption::PredefinedModel::ModelVariant& predefined_model_data);
}  // namespace Absorption::PredefinedModel

#endif  // fullmodel_h
