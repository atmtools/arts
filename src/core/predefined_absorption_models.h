/*!
 * @file   fullmodel.h
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */


#ifndef fullmodel_h
#define fullmodel_h

#include "atm.h"
#include "jacobian.h"
#include "predefined/predef_data.h"
#include "species.h"
#include <algorithm>
#include <random>
#include <rtepack.h>

namespace Absorption::PredefinedModel {
/** Contains known required VMR values
 *
 *  Note: If you add a species here, add it to the Jacobian
 *  wrapper in the compute function.
 */
struct VMRS {
  Numeric CO2{0};
  Numeric O2{0};
  Numeric N2{0};
  Numeric H2O{0};
  Numeric LWC{0};

  /**  Construct a new VMRS object
   * 
   * @param abs_species As WSV
   * @param rtp_vmr As WSV
   */
  VMRS(const ArrayOfArrayOfSpeciesTag& abs_species, const Vector& rtp_vmr) :
  CO2(Species::first_vmr(abs_species, rtp_vmr, Species::Species::CarbonDioxide)),
  O2(Species::first_vmr(abs_species, rtp_vmr, Species::Species::Oxygen)),
  N2(Species::first_vmr(abs_species, rtp_vmr, Species::Species::Nitrogen)),
  H2O(Species::first_vmr(abs_species, rtp_vmr, Species::Species::Water)),
  LWC(Species::first_vmr(abs_species, rtp_vmr, Species::Species::liquidcloud))
  {}

  /**  Construct a new VMRS object
   *
   * @param atm An Atmosphere
   */
  VMRS(const AtmPoint &atm)
      : CO2(atm[Species::Species::CarbonDioxide]),
        O2(atm[Species::Species::Oxygen]), N2(atm[Species::Species::Nitrogen]),
        H2O(atm[Species::Species::Water]),
        LWC(atm[Species::Species::liquidcloud]) {}

  constexpr VMRS() = default;

  friend std::ostream& operator<<(std::ostream& os, const VMRS& vmrs) {
    return os << "O2: " << vmrs.O2 << '\n' <<
                 "N2: " << vmrs.N2 << '\n' <<
                 "H2O: " << vmrs.H2O << '\n';
  }
};

//! Returns true if the model can be computed
bool can_compute(const SpeciesIsotopeRecord& model);

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
 * @param[in] jacobian_quantities As WSV
 * @param[in] predefined_model_data As WSV
 */
void compute(
    PropmatVector& propmat_clearsky,
    PropmatMatrix& dpropmat_clearsky_dx,
    const SpeciesIsotopeRecord& tag,
    const Vector& f_grid,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const VMRS& vmr,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const PredefinedModelData& predefined_model_data);
} // namespace Absorption::PredefinedModel

#endif  // fullmodel_h
