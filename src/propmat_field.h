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
   \file   propmat_field.h

   This file contains internal code to speed up propagation field
   calculations by sacrificing memory and by interpolating from
   potentially coarser resolution.

   \author Richard Larsson
   \date   2019-02-26
*/

#ifndef PROPAGATION_FIELD_HEADER
#define PROPAGATION_FIELD_HEADER

#include "propagationmatrix.h"
#include "field.h"
#include "rte.h"

typedef Field3D<PropagationMatrix> FieldOfPropagationMatrix;
typedef Field3D<StokesVector> FieldOfStokesVector;

void field_of_propagation(Workspace&                        ws,
                          FieldOfPropagationMatrix&         propmat_field,
                          FieldOfStokesVector&              absorption_field,
                          FieldOfStokesVector&              additional_source_field,
                          const Index&                      stokes_dim,
                          const Vector&                     f_grid,
                          const Vector&                     p_grid,
                          const Tensor3&                    z_field,
                          const Tensor3&                    t_field,
                          const Tensor4&                    nlte_field,
                          const Tensor4&                    vmr_field,
                          const ArrayOfRetrievalQuantity&   jacobian_quantities,
                          const Agenda&                     propmat_clearsky_agenda);

#endif  // PROPAGATION_FIELD_HEADER
