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

/**
  * @file   propmat_field.c
  * @author Richard Larsson
  * @date   2019-02-26
  * 
  * @brief Implements a propagation matrix field
  *
  * This file contains internal code to speed up propagation field
  * calculations by sacrificing memory and by interpolating from
  * potentially coarser resolution.
*/

#include "propmat_field.h"
#include "rte.h"
#include "special_interp.h"
#include "transmissionmatrix.h"

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
                          const Agenda& propmat_clearsky_agenda)
{
  const Index nalt = z_field.npages();
  const Index nlat = z_field.nrows();
  const Index nlon = z_field.ncols();
  const Index nq = jacobian_quantities.nelem();
  const Index nf = f_grid.nelem();

  ARTS_USER_ERROR_IF (nq,
        "Does not support Jacobian calculations at this time");
  ARTS_USER_ERROR_IF (stokes_dim not_eq 1,
        "Only for stokes_dim 1 at this time.");

  // Compute variables
  const Vector mag_field = Vector(3, 0);
  const Vector los = Vector(2, 0);
  const Vector tmp(0);
  ArrayOfStokesVector dS_dx(nq);
  ArrayOfPropagationMatrix dK_dx(nq);

  propmat_field = FieldOfPropagationMatrix(
      nalt, nlat, nlon, PropagationMatrix(nf, stokes_dim));
  absorption_field =
      FieldOfStokesVector(nalt, nlat, nlon, StokesVector(nf, stokes_dim));
  additional_source_field =
      FieldOfStokesVector(nalt, nlat, nlon, StokesVector(nf, stokes_dim));

  Workspace l_ws(ws);
  Agenda l_propmat_clearsky_agenda(propmat_clearsky_agenda);

#pragma omp parallel for if (not arts_omp_in_parallel()) schedule(guided) \
    firstprivate(l_ws, l_propmat_clearsky_agenda)
  for (Index i = 0; i < nalt; i++) {
    for (Index j = 0; j < nlat; j++) {
      for (Index k = 0; k < nlon; k++) {
        thread_local Index itmp;
        get_stepwise_clearsky_propmat(
            l_ws,
            propmat_field(i, j, k),
            additional_source_field(i, j, k),
            itmp,
            dK_dx,
            dS_dx,
            l_propmat_clearsky_agenda,
            jacobian_quantities,
            f_grid,
            mag_field,
            los,
            nlte_field(i, j, k),
            vmr_field(joker, i, j, k),
            t_field(i, j, k),
            p_grid[i],
            0);
        absorption_field(i, j, k) = propmat_field(i, j, k);
      }
    }
  }
}

FieldOfTransmissionMatrix transmat_field_calc_from_propmat_field(
    const FieldOfPropagationMatrix& propmat_field, const Numeric& r)
{
  FieldOfTransmissionMatrix transmat_field(
      propmat_field.npages(), propmat_field.nrows(), propmat_field.ncols());
  for (size_t ip = 0; ip < propmat_field.npages(); ip++)
    for (size_t ir = 0; ir < propmat_field.nrows(); ir++)
      for (size_t ic = 0; ic < propmat_field.ncols(); ic++)
        transmat_field(ip, ir, ic) =
            TransmissionMatrix(propmat_field(ip, ir, ic), r);
  return transmat_field;
}

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
    const Verbosity& verbosity)
{
  // Size of problem
  const Index nf = f_grid.nelem();
  const Index ns = propmat_field(0, 0, 0).StokesDimensions();
  const Index np = ppath.np;

  // Current limitations
  ARTS_USER_ERROR_IF (ns not_eq 1, "Only for stokes_dim 1");
  ARTS_USER_ERROR_IF (ppath.dim not_eq 1, "Only for atmosphere_dim 1");

  // Size of compute variables
  lvl_rad = ArrayOfRadiationVector(np, RadiationVector(nf, ns));
  src_rad = ArrayOfRadiationVector(np, RadiationVector(nf, ns));
  lyr_tra = ArrayOfTransmissionMatrix(np, TransmissionMatrix(nf, ns));

  // Size radiative variables always used
  Vector B(nf);
  PropagationMatrix K_this(nf, ns), K_past(nf, ns);

  // Temporary empty variables to fit available function handles
  Vector vtmp(0);
  ArrayOfTensor3 t3tmp;
  ArrayOfTransmissionMatrix tmtmp(0);
  ArrayOfRadiationVector rvtmp(0);
  ArrayOfPropagationMatrix pmtmp(0);
  ArrayOfStokesVector svtmp(0);
  ArrayOfRetrievalQuantity rqtmp(0);

  // Loop ppath points and determine radiative properties
  for (Index ip = 0; ip < np; ip++) {
    get_stepwise_blackbody_radiation(
        B, vtmp, f_grid, interp_atmfield_by_gp(1, t_field, ppath.gp_p[ip]), 0);
    K_this = propmat_field(ppath.gp_p[ip]);
    const StokesVector S(additional_source_field(ppath.gp_p[ip]));
    const StokesVector a(absorption_field(ppath.gp_p[ip]));

    if (ip)
      stepwise_transmission(lyr_tra[ip],
                            tmtmp,
                            tmtmp,
                            K_past,
                            K_this,
                            pmtmp,
                            pmtmp,
                            ppath.lstep[ip - 1],
                            0,
                            0,
                            -1);

    stepwise_source(src_rad[ip],
                    rvtmp,
                    K_this,
                    a,
                    S,
                    pmtmp,
                    svtmp,
                    svtmp,
                    B,
                    vtmp,
                    rqtmp,
                    0);

    swap(K_past, K_this);
  }

  // In case of backwards RT necessary
  tot_tra = cumulative_transmission(lyr_tra, CumulativeTransmission::Forward);
  const Tensor3 iy_trans_new = tot_tra[np - 1];

  // Radiative background
  Matrix iy;
  get_iy_of_background(ws,
                       iy,
                       t3tmp,
                       iy_trans_new,
                       0,
                       0,
                       rqtmp,
                       ppath,
                       {0},
                       1,
                       nlte_field,
                       0,
                       1,
                       f_grid,
                       "1",
                       surface_props_data,
                       iy_main_agenda,
                       iy_space_agenda,
                       iy_surface_agenda,
                       iy_cloudbox_agenda,
                       1,
                       verbosity);
  lvl_rad[np - 1] = iy;

  // Radiative transfer calculations
  for (Index ip = np - 2; ip >= 0; ip--)
    update_radiation_vector(lvl_rad[ip] = lvl_rad[ip + 1],
                            rvtmp,
                            rvtmp,
                            src_rad[ip],
                            src_rad[ip + 1],
                            rvtmp,
                            rvtmp,
                            lyr_tra[ip + 1],
                            tot_tra[ip],
                            tmtmp,
                            tmtmp,
                            PropagationMatrix(),
                            PropagationMatrix(),
                            ArrayOfPropagationMatrix(),
                            ArrayOfPropagationMatrix(),
                            Numeric(),
                            Vector(),
                            Vector(),
                            0,
                            0,
                            RadiativeTransferSolver::Emission);
}
