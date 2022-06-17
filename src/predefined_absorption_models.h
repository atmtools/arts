/* Copyright (C) 2020
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/*!
 * @file   fullmodel.h
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */


#ifndef fullmodel_h
#define fullmodel_h

#include "jacobian.h"
#include "predefined/predef_data.h"
#include "propagationmatrix.h"
#include "species.h"
#include <algorithm>
#include <random>

namespace Absorption::PredefinedModel {
/** Contains known required VMR values
 *
 *  Note: If you add a species here, add it to the Jacobian
 *  wrapper in the compute function.
 */
struct VMRS {
  Numeric O2{0};
  Numeric H2O{0};

  /**  Construct a new VMRS object
   * 
   * @param abs_species 
   * @param rtp_vmr 
   */
  VMRS(const ArrayOfArrayOfSpeciesTag& abs_species, const Vector& rtp_vmr) :
  O2(Species::first_vmr(abs_species, rtp_vmr, Species::Species::Oxygen)),
  H2O(Species::first_vmr(abs_species, rtp_vmr, Species::Species::Water))
  {}

  constexpr VMRS() = default;

  friend std::ostream& operator<<(std::ostream& os, const VMRS& vmrs) {
    return os << "O2: " << vmrs.O2 << '\n' <<
                 "H2O: " << vmrs.H2O << '\n';
  }
};

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
    PropagationMatrix& propmat_clearsky,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    const SpeciesIsotopeRecord& tag,
    const Vector& f_grid,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const VMRS& vmr,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const PredefinedModelData& predefined_model_data);
} // namespace Absorption::PredefinedModel

#endif  // fullmodel_h
