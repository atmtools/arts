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

namespace Absorption::PredefinedModel {

/** Adds Makarov MPM2020 O2 absorption lines to the propagation matrix
 * 
 * Adds negative values outside of bounds.  Works for Earth only
 * 
 * Water is 10% more effective at broadening than dry air so must
 * be included.  Does not deal with Zeeman effect and ignores negative
 * absorption far from the line center (this can be fixed but requires
 * clear-cut motivation)
 * 
 * @param[in,out] propmat_clearsky As WSV
 * @param[in,out] dpropmat_clearsky_dx As WSV
 * @param[in]     f Frequency grid of computations (size: [f])
 * @param[in]     p Pressure of computations
 * @param[in]     t Temperature of computations
 * @param[in]     oxygen_vmr O2 volume mixing ratio (only for strength scaling)
 * @param[in]     water_vmr Water volume mixing ratio
 * @param[in]     self_vmr Molecular oxygen volume mixing ratio
 * @param[in]     jacs The Jacobian descriptions (size: [greater than max(jacs_pos)])
 */
void makarov2020_o2_lines_mpm(PropagationMatrix& propmat_clearsky,
                              ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                              const Vector& f,
                              const Numeric& p,
                              const Numeric& t,
                              const Numeric& oxygen_vmr,
                              const Numeric& water_vmr,
                              const ArrayOfRetrievalQuantity& jacs);
} // namespace Absorption::PredefinedModel

#endif  // fullmodel_h
