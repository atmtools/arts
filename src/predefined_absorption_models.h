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
#include "propagationmatrix.h"
#include "species.h"
#include <algorithm>
#include <random>

namespace Absorption::PredefinedModel {
/** Contains known required VMR values */
struct VMRS {
  Numeric O2{0};
  Numeric H2O{0};

  VMRS(const ArrayOfArrayOfSpeciesTag& specs, const Vector& rtp_vmr) :
  O2(Species::first_vmr(specs, rtp_vmr, Species::Species::Oxygen)),
  H2O(Species::first_vmr(specs, rtp_vmr, Species::Species::Water))
  {}
  
  constexpr VMRS() = default;
};

void compute(
    PropagationMatrix& propmat_clearsky,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    const SpeciesIsotopeRecord& tag,
    const Vector& f_grid,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const VMRS& vmr,
    const ArrayOfRetrievalQuantity& jacobian_quantities);
} // namespace Absorption::PredefinedModel

#endif  // fullmodel_h
