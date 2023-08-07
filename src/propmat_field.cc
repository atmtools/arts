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
#include "matpack_data.h"
#include "physics_funcs.h"
#include "rte.h"
#include "special_interp.h"
#include "arts_omp.h"

void field_of_propagation(const Workspace& ws,
                          FieldOfPropmatVector& propmat_field,
                          FieldOfStokvecVector& absorption_field,
                          FieldOfStokvecVector& additional_source_field,
                          const Vector& f_grid,
                          const AtmField& atm_field,
                          const ArrayOfRetrievalQuantity& jacobian_quantities,
                          const Agenda& propmat_clearsky_agenda)
{
//  const Index nalt = atm_field.regularized_shape()[0];
//  const Index nlat = atm_field.regularized_shape()[1];
//  const Index nlon = atm_field.regularized_shape()[2];
ARTS_USER_ERROR("ERROR")
Index nalt, nlat, nlon;
  const Index nq = jacobian_quantities.nelem();
  const Index nf = f_grid.nelem();

  ARTS_USER_ERROR_IF (nq,
        "Does not support Jacobian calculations at this time");

  // Compute variables
  const Vector mag_field = Vector(3, 0);
  const Vector los = Vector(2, 0);
  const Vector tmp(0);
  StokvecMatrix dS_dx(nq, nf);
  PropmatMatrix dK_dx(nq, nf);

  propmat_field = FieldOfPropmatVector(
      nalt, nlat, nlon, PropmatVector(nf));
  absorption_field =
      FieldOfStokvecVector(nalt, nlat, nlon, StokvecVector(nf));
  additional_source_field =
      FieldOfStokvecVector(nalt, nlat, nlon, StokvecVector(nf));

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(3)
  for (Index i = 0; i < nalt; i++) {
    for (Index j = 0; j < nlat; j++) {
      for (Index k = 0; k < nlon; k++) {
        ARTS_USER_ERROR("ERROR")
        get_stepwise_clearsky_propmat(
            ws,
            propmat_field(i, j, k),
            additional_source_field(i, j, k),
            dK_dx,
            dS_dx,
            propmat_clearsky_agenda,
            jacobian_quantities,
            f_grid,
            los,
       AtmPoint{},  //   atm_field.at({atm_field.grid[0][i]}, {atm_field.grid[1][j]}, {atm_field.grid[2][k]})[0],
            false);
        absorption_field(i, j, k) = absvec(propmat_field(i, j, k));
      }
    }
  }
}

FieldOfMuelmatVector transmat_field_calc_from_propmat_field(
    const FieldOfPropmatVector& propmat_field, const Numeric& r)
{
  FieldOfMuelmatVector transmat_field(
      propmat_field.npages(), propmat_field.nrows(), propmat_field.ncols());
  for (size_t ip = 0; ip < propmat_field.npages(); ip++)
    for (size_t ir = 0; ir < propmat_field.nrows(); ir++)
      for (size_t ic = 0; ic < propmat_field.ncols(); ic++)
        for (Index iv = 0; iv < propmat_field(0, 0, 0).nelem(); iv++)
          transmat_field(ip, ir, ic)[iv] = rtepack::exp(-r * propmat_field(ip, ir, ic)[iv]);
  return transmat_field;
}

void emission_from_propmat_field(
    const Workspace& ws,
    ArrayOfStokvecVector& lvl_rad,
    ArrayOfStokvecVector& src_rad,
    ArrayOfMuelmatVector& lyr_tra,
    ArrayOfMuelmatVector& tot_tra,
    const FieldOfPropmatVector& propmat_field,
    const FieldOfStokvecVector& absorption_field,
    const FieldOfStokvecVector& additional_source_field,
    const Vector& f_grid,
    const AtmField& atm_field,
    const Ppath& ppath,
    const Agenda& iy_main_agenda,
    const Agenda& iy_space_agenda,
    const Agenda& iy_surface_agenda,
    const Agenda& iy_cloudbox_agenda,
    const SurfaceField& surface_field)
{
  // Size of problem
  const Index nf = f_grid.nelem();
  const Index np = ppath.np;

  // Current limitations
  ARTS_USER_ERROR ("Only for 1D atmospheres at this time");

  // Size of compute variables
  lvl_rad.resize(np, nf);
  src_rad.resize(np, nf);
  lyr_tra.resize(np, nf);

  // Size radiative variables always used
  Vector B(nf);
  PropmatVector K_this(nf), K_past(nf);

  // Temporary empty variables to fit available function handles
  Vector vtmp(0);
  ArrayOfTensor3 t3tmp;
  const ArrayOfRetrievalQuantity rqtmp(0);

  // Loop ppath points and determine radiative properties
  for (Index ip = 0; ip < np; ip++) {
    std::transform(f_grid.elem_begin(), f_grid.elem_end(), B.elem_begin(),
                   [T = atm_field[Atm::Key::t].at(
                        ppath.pos[ip][0], ppath.pos[ip][1], ppath.pos[ip][2])](
                       const Numeric f) { return planck(f, T); });

    K_this = propmat_field(ppath.gp_p[ip]);
    const StokvecVector S(additional_source_field(ppath.gp_p[ip]));
    const StokvecVector a(absorption_field(ppath.gp_p[ip]));

    if (ip)
      two_level_exp(lyr_tra[ip], K_past, K_this, ppath.lstep[ip - 1]);
    rtepack::source::level_nlte_and_scattering(src_rad[ip], K_this, a, S, B);

    swap(K_past, K_this);
  }

  // In case of backwards RT necessary
  tot_tra = forward_cumulative_transmission(lyr_tra);
  const Tensor3 iy_trans_new = to_tensor3(tot_tra[np - 1]);

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
                       Vector{0},
                       atm_field,
                       0,
                       f_grid,
                       "1",
                       surface_field,
                       iy_main_agenda,
                       iy_space_agenda,
                       iy_surface_agenda,
                       iy_cloudbox_agenda,
                       1);
  lvl_rad[np - 1] = rtepack::to_stokvec_vector(iy);

  // Radiative transfer calculations
  for (Index ip = np - 2; ip >= 0; ip--)
  two_level_linear_emission_step(lvl_rad[ip] = lvl_rad[ip + 1],
                            src_rad[ip],
                            src_rad[ip + 1], lyr_tra[ip + 1]);
}
