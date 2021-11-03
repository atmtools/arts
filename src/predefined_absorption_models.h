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

namespace Absorption::PredefinedModel::Makarov2020etal {

/** Adds Makarov MPM2020 O2 absorption lines to the propagation matrix
 *
 * Adds the 60 GHz O2-66 band, including the 118.75 GHz line.  No Zeeman
 * effect is considered.
 * 
 * Adds no negative values outside of bounds.
 *
 * Jacobian values are computed by perturbations, only
 * frequency, temperature, and the two styles of VMR are supported
 * 
 * @param[in,out] propmat_clearsky As WSV
 * @param[in,out] dpropmat_clearsky_dx As WSV
 * @param[in]     f_grid As WSV
 * @param[in]     rtp_pressure As WSV
 * @param[in]     rtp_temperature As WSV
 * @param[in]     oxygen_vmr O2 volume mixing ratio (only for strength scaling)
 * @param[in]     jacobian_quantities As WSV
 */
void compute(PropagationMatrix& propmat_clearsky,
             ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
             const Vector& f_grid,
             const Numeric& rtp_pressure,
             const Numeric& rtp_temperature,
             const Numeric& oxygen_vmr,
             const ArrayOfRetrievalQuantity& jacobian_quantities) ARTS_NOEXCEPT;
} // namespace Absorption::PredefinedModel::Makarov2020etal

#endif  // fullmodel_h
