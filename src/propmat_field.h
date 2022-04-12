/* Copyright (C) 2019 Richard Larsson <ric.larsson@gmail.com>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/**
  * @file   propmat_field.h
  * @author Richard Larsson
  * @date   2019-02-26
  * 
  * @brief Implements a propagation matrix field
  *
  * This file contains internal code to speed up propagation field
  * calculations by sacrificing memory and by interpolating from
  * potentially coarser resolution.
*/

#ifndef PROPAGATION_FIELD_HEADER
#define PROPAGATION_FIELD_HEADER

#include "energylevelmap.h"
#include "field.h"
#include "transmissionmatrix.h"

using FieldOfTransmissionMatrix = Field3D<TransmissionMatrix>;
using FieldOfPropagationMatrix = Field3D<PropagationMatrix>;
using FieldOfStokesVector = Field3D<StokesVector>;

/** Creates a field of propagation matrices, absorption vectors, and source vectors
 * 
 * @param[in] ws A workspace
 * @param[out] propmat_field A 3D field of propagation matrices
 * @param[out] absorption_field A 3D field of absorption vectors
 * @param[out] additional_source_field A 3D field of source vectors
 * @param[in] stokes_dim As WSV
 * @param[in] f_grid As WSV
 * @param[in] p_grid As WSV
 * @param[in] z_field As WSV
 * @param[in] t_field As WSV
 * @param[in] nlte_field As WSV
 * @param[in] vmr_field As WSV
 * @param[in] jacobian_quantities As WSV
 * @param[in] propmat_clearsky_agenda As WSA
 */
void field_of_propagation(Workspace& ws,
                          FieldOfPropagationMatrix& propmat_field,
                          FieldOfStokesVector& absorption_field,
                          FieldOfStokesVector& additional_source_field,
                          const Index& stokes_dim,
                          const Vector& f_grid,
                          const Vector& p_grid,
                          const Tensor3& z_field,
                          const Tensor3& t_field,
                          const EnergyLevelMap& nlte_field,
                          const Tensor4& vmr_field,
                          const ArrayOfRetrievalQuantity& jacobian_quantities,
                          const Agenda& propmat_clearsky_agenda);

/** Get a field of transmission matrices from the propagation matrix field
 * 
 * @param[in] propmat_field A 3D field of propagation matrices
 * @param[in] r The distance
 * @return FieldOfTransmissionMatrix 
 */
FieldOfTransmissionMatrix transmat_field_calc_from_propmat_field(
    const FieldOfPropagationMatrix& propmat_field, const Numeric& r = 1.0);

/** Computes the radiation and transmission from fields of atmospheric propagation
 * 
 * Only for 1D atmospheres for now.
 * 
 * Computes The forward simulations by interpolating the fields of
 * radiative properties to the selected propagation path.
 * 
 * Not well-tested.
 * 
 * @param[in] ws A workspace
 * @param[out] lvl_rad Level by level radiation
 * @param[out] src_rad Level by level source function
 * @param[out] lyr_tra Layered transmission
 * @param[out] tot_tra Total transmission from layer to background
 * @param[in] propmat_field 3D field of propagation matrices
 * @param[in] absorption_field A 3D field of absorption vectors
 * @param[in] additional_source_field A 3D field of source vectors
 * @param[in] f_grid As WSV
 * @param[in] t_field As WSV
 * @param[in] nlte_field As WSV
 * @param[in] ppath As WSV
 * @param[in] iy_main_agenda As WSA
 * @param[in] iy_space_agenda As WSA
 * @param[in] iy_surface_agenda As WSA
 * @param[in] iy_cloudbox_agenda As WSA
 * @param[in] surface_props_data As WSV
 * @param[in] verbosity Level of verbosity in underlying calls
 */
void emission_from_propmat_field(
    Workspace& ws,
    ArrayOfRadiationVector& lvl_rad,
    ArrayOfRadiationVector& src_rad,
    ArrayOfTransmissionMatrix& lyr_tra,
    ArrayOfTransmissionMatrix& tot_tra,
    const FieldOfPropagationMatrix& propmat_field,
    const FieldOfStokesVector& absorption_field,
    const FieldOfStokesVector& additional_source_field,
    const Vector& f_grid,
    const Tensor3& t_field,
    const EnergyLevelMap& nlte_field,
    const Ppath& ppath,
    const Agenda& iy_main_agenda,
    const Agenda& iy_space_agenda,
    const Agenda& iy_surface_agenda,
    const Agenda& iy_cloudbox_agenda,
    const Tensor3& surface_props_data,
    const Verbosity& verbosity);
void emission_from_propmat_field(
  Workspace& ws,
  ArrayOfRadiationVector& lvl_rad,
  ArrayOfRadiationVector& src_rad,
  ArrayOfTransmissionMatrix& lyr_tra,
  ArrayOfTransmissionMatrix& tot_tra,
  const FieldOfPropagationMatrix& propmat_field,
  const FieldOfStokesVector& absorption_field,
  const FieldOfStokesVector& additional_source_field,
  const Vector& f_grid,
  const Tensor3& t_field,
  const EnergyLevelMap& nlte_field,
  const Ppath& ppath,
  const Agenda& iy_main_agenda,
  const Agenda& iy_space_agenda,
  const Agenda& iy_surface_agenda,
  const Agenda& iy_cloudbox_agenda,
  const Tensor3& surface_props_data,
  const Verbosity& verbosity);

#endif  // PROPAGATION_FIELD_HEADER
