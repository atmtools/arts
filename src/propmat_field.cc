/* Copyright (C) 2019 Richard Larsson <ric.larsson@gmail.com>
 * 
 T his pr*ogram is free software; you can redistribute it and/or modify it
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


#include "propmat_field.h"

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
                          const Agenda&                     propmat_clearsky_agenda)
{
  const Index nalt = z_field.npages();
  const Index nlat = z_field.nrows();
  const Index nlon = z_field.ncols();
  const Index nq = jacobian_quantities.nelem();
  const Index nf = f_grid.nelem();
  
  if(nq)
    throw std::runtime_error("Does not support Jacobian calculations at this time");
  if(stokes_dim not_eq 1)
    throw std::runtime_error("Only for stokes_dim 1 at this time.");

  // Compute variables  
  const Vector mag_field = Vector(3, 0);
  const Vector los = Vector(2, 0);
  const ArrayOfIndex spec_jac(nq);
  const Vector tmp(0);
  ArrayOfIndex lte(nalt, 0);
  ArrayOfStokesVector dS_dx(nq);
  ArrayOfPropagationMatrix dK_dx(nq);
  
  propmat_field           = FieldOfPropagationMatrix(nalt, nlat, nlon, PropagationMatrix(nf, stokes_dim));
  absorption_field        = FieldOfStokesVector(     nalt, nlat, nlon, StokesVector(     nf, stokes_dim));
  additional_source_field = FieldOfStokesVector(     nalt, nlat, nlon, StokesVector(     nf, stokes_dim));
  for(Index i=0; i<nalt; i++) {
    for(Index j=0; j<nlat; j++) {
      for(Index k=0; k<nlon; k++) {
        get_stepwise_clearsky_propmat(ws, propmat_field(i, j, k),
                                      additional_source_field(i, j, k),
                                      lte[i], dK_dx, dS_dx,
                                      propmat_clearsky_agenda,
                                      jacobian_quantities,
                                      f_grid, mag_field, los,
                                      nlte_field.empty()?tmp[joker]:nlte_field(joker, i, j, k),
                                      vmr_field(joker, i, j, k),
                                      t_field(i, j, k), p_grid[i],
                                      spec_jac, 0);
        absorption_field(i, j, k) = propmat_field(i, j, k);
      }
    }
  }
}
