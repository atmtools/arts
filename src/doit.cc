/* Copyright (C) 2002-2012 Claudia Emde <claudia.emde@dlr.de>
                      
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
   USA. 
*/

/*!
  \file   doit.cc
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   Wed Jun 04 11:03:57 2003
  
  \brief  This file contains functions to calculate the radiative transfer
  inside the cloudbox using the DOIT method.
  
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "doit.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "agenda_class.h"
#include "array.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include "atm.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "debug.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "propagationmatrix.h"
#include "rte.h"
#include "sorting.h"
#include "special_interp.h"
#include "species_tags.h"
#include "xml_io.h"

inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);

//FIXME function name of 'rte_step_doit_replacement' should be replaced by
//proper name
void rte_step_doit_replacement(  //Output and Input:
    VectorView stokes_vec,
    MatrixView trans_mat,
    //Input
    const PropagationMatrix& ext_mat_av,
    const StokesVector& abs_vec_av,
    ConstVectorView sca_vec_av,
    const Numeric& lstep,
    const Numeric& rtp_planck_value,
    const bool& trans_is_precalc) {
  //Stokes dimension:
  Index stokes_dim = stokes_vec.nelem();

  //Test sizes
  ARTS_ASSERT(ext_mat_av.NumberOfFrequencies() == 1 and
         abs_vec_av.NumberOfFrequencies() == 1);

  //Check inputs:
  ARTS_ASSERT(is_size(trans_mat, stokes_dim, stokes_dim));
  ARTS_ASSERT(stokes_dim == ext_mat_av.StokesDimensions() and
         stokes_dim == abs_vec_av.StokesDimensions());
  ARTS_ASSERT(is_size(sca_vec_av, stokes_dim));
  ARTS_ASSERT(rtp_planck_value >= 0);
  ARTS_ASSERT(lstep >= 0);
  //ARTS_ASSERT (not ext_mat_av.AnySingular());  This is ARTS_ASSERTed at a later time in this version...

  // Check, if only the first component of abs_vec is non-zero:
//   const bool unpol_abs_vec = abs_vec_av.IsUnpolarized(0);

//   bool unpol_sca_vec = true;

//   for (Index i = 1; i < stokes_dim; i++)
//     if (sca_vec_av[i] != 0) unpol_sca_vec = false;

  // Calculate transmission by general function, if not precalculated
//   Index extmat_case = 0;
  if (!trans_is_precalc) {
    compute_transmission_matrix_from_averaged_matrix_at_frequency(
        trans_mat, lstep, ext_mat_av, 0);
  }

  //--- Scalar case: ---------------------------------------------------------
  if (stokes_dim == 1) {
    stokes_vec[0] = stokes_vec[0] * trans_mat(0, 0) +
                    (abs_vec_av.Kjj()[0] * rtp_planck_value + sca_vec_av[0]) /
                        ext_mat_av.Kjj()[0] * (1 - trans_mat(0, 0));
  }

  //--- Vector case: ---------------------------------------------------------

  // We have here two cases, diagonal or non-diagonal ext_mat_gas
  // For diagonal ext_mat_gas, we expect abs_vec_gas to only have a
  // non-zero value in position 1.

  //- Unpolarised
//   else if (extmat_case == 1 && unpol_abs_vec && unpol_sca_vec) {
//     const Numeric invK = 1.0 / ext_mat_av.Kjj()[0];
//     // Stokes dim 1
//     stokes_vec[0] = stokes_vec[0] * trans_mat(0, 0) +
//                     (abs_vec_av.Kjj()[0] * rtp_planck_value + sca_vec_av[0]) *
//                         invK * (1 - trans_mat(0, 0));
// 
//     // Stokes dims > 1
//     for (Index i = 1; i < stokes_dim; i++) {
//       stokes_vec[i] = stokes_vec[i] * trans_mat(i, i) +
//                       sca_vec_av[i] * invK * (1 - trans_mat(i, i));
//     }
//   }

  //- General case
  else {
    Matrix invK(stokes_dim, stokes_dim);
    ext_mat_av.MatrixInverseAtPosition(invK);

    Vector source{abs_vec_av.VectorAtPosition()};
    source *= rtp_planck_value;

    for (Index i = 0; i < stokes_dim; i++)
      source[i] += sca_vec_av[i];  // b = abs_vec * B + sca_vec

    // solve K^(-1)*b = x
    Vector x(stokes_dim);
    mult(x, invK, source);

    Vector term1(stokes_dim);
    Vector term2(stokes_dim);

    Matrix ImT(stokes_dim, stokes_dim);
    id_mat(ImT);
    ImT -= trans_mat;
    mult(term2, ImT, x);  // term2: second term of the solution of the RTE with
                          //fixed scattered field

    // term1: first term of solution of the RTE with fixed scattered field
    mult(term1, trans_mat, stokes_vec);

    for (Index i = 0; i < stokes_dim; i++)
      stokes_vec[i] = term1[i] + term2[i];  // Compute the new Stokes Vector
  }
}

void cloud_fieldsCalc(Workspace& ws,
                      // Output and Input:
                      Tensor5View ext_mat_field,
                      Tensor4View abs_vec_field,
                      // Input:
                      const Agenda& spt_calc_agenda,
                      const Index& za_index,
                      const Index& aa_index,
                      const ArrayOfIndex& cloudbox_limits,
                      ConstTensor3View t_field,
                      ConstTensor4View pnd_field) {
  // Input variables are checked in the WSMs i_fieldUpdateSeqXXX, from
  // where this function is called.

  const Index atmosphere_dim = cloudbox_limits.nelem() / 2;
  const Index N_se = pnd_field.nbooks();
  const Index stokes_dim = ext_mat_field.ncols();

  ARTS_ASSERT(atmosphere_dim == 1 || atmosphere_dim == 3);
  ARTS_ASSERT(ext_mat_field.ncols() == ext_mat_field.nrows() &&
         ext_mat_field.ncols() == abs_vec_field.ncols());

  const Index Np_cloud = cloudbox_limits[1] - cloudbox_limits[0] + 1;

  // If atmosohere_dim == 1
  Index Nlat_cloud = 1;
  Index Nlon_cloud = 1;

  if (atmosphere_dim == 3) {
    Nlat_cloud = cloudbox_limits[3] - cloudbox_limits[2] + 1;
    Nlon_cloud = cloudbox_limits[5] - cloudbox_limits[4] + 1;
  }

  // Initialize ext_mat(_spt), abs_vec(_spt)
  // Resize and initialize variables for storing optical properties
  // of all scattering elements.
  ArrayOfStokesVector abs_vec_spt_local(N_se);
  for (auto& av : abs_vec_spt_local) {
    av = StokesVector(1, stokes_dim);
    av.SetZero();
  }
  ArrayOfPropagationMatrix ext_mat_spt_local(N_se);
  for (auto& pm : ext_mat_spt_local) {
    pm = PropagationMatrix(1, stokes_dim);
    pm.SetZero();
  }

  StokesVector abs_vec_local;
  PropagationMatrix ext_mat_local;
  Numeric rtp_temperature_local;

  // Calculate ext_mat, abs_vec for all points inside the cloudbox.
  // sca_vec can be obtained from the workspace variable doit_scat_field.
  // As we need the average for each layer, it makes sense to calculte
  // the coefficients once and store them in an array instead of
  // calculating at each point the coefficient of the point above and
  // the point below.
  // To use special interpolation functions for atmospheric fields we
  // use ext_mat_field and abs_vec_field:

  // Loop over all positions inside the cloudbox defined by the
  // cloudbox_limits.
  for (Index scat_p_index_local = 0; scat_p_index_local < Np_cloud;
       scat_p_index_local++) {
    for (Index scat_lat_index_local = 0; scat_lat_index_local < Nlat_cloud;
         scat_lat_index_local++) {
      for (Index scat_lon_index_local = 0; scat_lon_index_local < Nlon_cloud;
           scat_lon_index_local++) {
        if (atmosphere_dim == 1)
          rtp_temperature_local =
              t_field(scat_p_index_local + cloudbox_limits[0], 0, 0);
        else
          rtp_temperature_local =
              t_field(scat_p_index_local + cloudbox_limits[0],
                      scat_lat_index_local + cloudbox_limits[2],
                      scat_lon_index_local + cloudbox_limits[4]);

        //Calculate optical properties for individual scattering elements:
        //( Execute agendas silently. )
        spt_calc_agendaExecute(ws,
                               ext_mat_spt_local,
                               abs_vec_spt_local,
                               scat_p_index_local,
                               scat_lat_index_local,
                               scat_lon_index_local,
                               rtp_temperature_local,
                               za_index,
                               aa_index,
                               spt_calc_agenda);
        /*
// so far missing here (accessed through workspace within agenda):
// - scat_data
// - za_grid, aa_grid
// - f_index
              opt_prop_sptFromScat_data(ext_mat_spt_local, abs_vec_spt_local,
                                        scat_data, 1,
                                        za_grid, aa_grid,
                                        za_index, aa_index,
                                        f_index,
                                        rtp_temperature_local,
                                        pnd_field, 
                                        scat_p_index_local,
                                        scat_lat_index_local,
                                        scat_lon_index_local,
                                        verbosity);
*/

        opt_prop_bulkCalc(ext_mat_local,
                          abs_vec_local,
                          ext_mat_spt_local,
                          abs_vec_spt_local,
                          Tensor4{pnd_field},
                          scat_p_index_local,
                          scat_lat_index_local,
                          scat_lon_index_local);

        // Store coefficients in arrays for the whole cloudbox.
        abs_vec_field(scat_p_index_local,
                      scat_lat_index_local,
                      scat_lon_index_local,
                      joker) = abs_vec_local.VectorAtPosition();

        ext_mat_local.MatrixAtPosition(ext_mat_field(scat_p_index_local,
                                                     scat_lat_index_local,
                                                     scat_lon_index_local,
                                                     joker,
                                                     joker));
      }
    }
  }
}

void cloud_ppath_update1D(Workspace& ws,
                          // Input and output
                          Tensor6View cloudbox_field_mono,
                          // ppath_step_agenda:
                          const Index& p_index,
                          const Index& za_index,
                          ConstVectorView za_grid,
                          const ArrayOfIndex& cloudbox_limits,
                          ConstTensor6View doit_scat_field,
                          // Calculate scalar gas absorption:
                          const Agenda& propmat_clearsky_agenda,
                          const AtmField& atm_field,
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          // Propagation path calculation:
                          const Agenda& ppath_step_agenda,
                          const Numeric& ppath_lmax,
                          const Numeric& ppath_lraytrace,
                          ConstVectorView refellipsoid,
                          // Calculate thermal emission:
                          ConstVectorView f_grid,
                          const Index& f_index,
                          //particle optical properties
                          ConstTensor5View ext_mat_field,
                          ConstTensor4View abs_vec_field,
                          const Agenda& surface_rtprop_agenda,
                          //const Agenda& surface_rtprop_agenda,
                          const Index& scat_za_interp) {
  Ppath ppath_step;
  // Input variables are checked in the WSMs i_fieldUpdateSeqXXX, from
  // where this function is called.

  //Initialize ppath for 1D.
  //ppath_init_structure(ppath_step, 1, 1);
  // See documentation of ppath_init_structure for understanding
  // the parameters.

//FIXME: MUST HAVE REGULAR FIELD
  Vector z_grid;
  Tensor3 t_field;
  Tensor4 vmr_field;
  ARTS_USER_ERROR("ERROR")
  //ARTS_USER_ERROR_IF(not atm_field.regularized, "Must have regular grid atmospheric field")
  //const auto& z_grid = atm_field.grid[0];
  //const auto& t_field = atm_field[Atm::Key::t].get<const Tensor3&>();
  //const auto vmr_field = Atm::extract_specs_content(atm_field, abs_species);

  // Assign value to ppath.pos:
  ppath_step.pos(0, 0) = z_grid[p_index];
  ppath_step.r[0] = refellipsoid[0] + z_grid[p_index];

  // Define the direction:
  ppath_step.los(0, 0) = za_grid[za_index];

  // Define the grid positions:
  ppath_step.gp_p[0].idx = p_index;
  ppath_step.gp_p[0].fd[0] = 0;
  ppath_step.gp_p[0].fd[1] = 1;

  // Call ppath_step_agenda:
  ppath_step_agendaExecute(ws,
                           ppath_step,
                           ppath_lmax,
                           ppath_lraytrace,
                           Vector(1, f_grid[f_index]),
                           ppath_step_agenda);

  // Check whether the next point is inside or outside the
  // cloudbox. Only if the next point lies inside the
  // cloudbox a radiative transfer step caclulation has to
  // be performed.

  if ((cloudbox_limits[0] <= ppath_step.gp_p[1].idx &&
       cloudbox_limits[1] > ppath_step.gp_p[1].idx) ||
      (cloudbox_limits[1] == ppath_step.gp_p[1].idx &&
       abs(ppath_step.gp_p[1].fd[0]) < 1e-6)) {
    // Stokes dimension
    const Index stokes_dim = cloudbox_field_mono.ncols();
    // Number of species
    const Index N_species = abs_species.nelem();

    // Ppath_step normally has 2 points, the starting
    // point and the intersection point.
    // But there can be points in between, when a maximum
    // lstep is given. We have to interpolate on all the
    // points in the ppath_step.

    // Initialize variables for interpolated values
    Tensor3 ext_mat_int(stokes_dim, stokes_dim, ppath_step.np, 0.);
    Matrix abs_vec_int(stokes_dim, ppath_step.np, 0.);
    Matrix sca_vec_int(stokes_dim, ppath_step.np, 0.);
    Matrix cloudbox_field_mono_int(stokes_dim, ppath_step.np, 0.);
    Vector t_int(ppath_step.np, 0.);
    Matrix vmr_list_int(N_species, ppath_step.np, 0.);
    Vector p_int(ppath_step.np, 0.);

    interp_cloud_coeff1D(ext_mat_int,
                         abs_vec_int,
                         sca_vec_int,
                         cloudbox_field_mono_int,
                         t_int,
                         vmr_list_int,
                         p_int,
                         ext_mat_field,
                         abs_vec_field,
                         doit_scat_field,
                         cloudbox_field_mono,
                         t_field,
                         vmr_field,
                         z_grid,
                         ppath_step,
                         cloudbox_limits,
                         za_grid,
                         scat_za_interp);

    // ppath_what_background(ppath_step) tells the
    // radiative background.  More information in the
    // function get_iy_of_background.
    // if there is no background we proceed the RT
    Index bkgr;// = ppath_what_background(ppath_step);
  ARTS_USER_ERROR("ERROR")
    // Radiative transfer from one layer to the next, starting
    // at the intersection with the next layer and propagating
    // to the considered point.
    cloud_RT_no_background(ws,
                           cloudbox_field_mono,
                           propmat_clearsky_agenda,
                           ppath_step,
                           t_int,
                           vmr_list_int,
                           ext_mat_int,
                           abs_vec_int,
                           sca_vec_int,
                           cloudbox_field_mono_int,
                           p_int,
                           cloudbox_limits,
                           f_grid,
                           f_index,
                           p_index,
                           0,
                           0,
                           za_index,
                           0);

    // bkgr=2 indicates that the background is the surface
    if (bkgr == 2) {
      // cout << "hit surface "<< ppath_step.gp_p << endl;
      cloud_RT_surface(ws,
                       cloudbox_field_mono,
                       surface_rtprop_agenda,
                       f_grid,
                       f_index,
                       stokes_dim,
                       ppath_step,
                       cloudbox_limits,
                       za_grid,
                       za_index);
    }

  }  //end if inside cloudbox
}

void cloud_ppath_update1D_noseq(Workspace& ws,
                                // Output
                                Tensor6View cloudbox_field_mono,
                                // ppath_step_agenda:
                                const Index& p_index,
                                const Index& za_index,
                                ConstVectorView za_grid,
                                const ArrayOfIndex& cloudbox_limits,
                                ConstTensor6View cloudbox_field_mono_old,
                                ConstTensor6View doit_scat_field,
                                // Calculate scalar gas absorption:
                                const Agenda& propmat_clearsky_agenda,
                                ConstTensor4View vmr_field,
                                // Propagation path calculation:
                                const Agenda& ppath_step_agenda,
                                const Numeric& ppath_lmax,
                                const Numeric& ppath_lraytrace,
                                ConstVectorView p_grid,
                                ConstTensor3View z_field,
                                ConstVectorView refellipsoid,
                                // Calculate thermal emission:
                                ConstTensor3View t_field,
                                ConstVectorView f_grid,
                                // used for surface ?
                                const Index& f_index,
                                //particle optical properties
                                ConstTensor5View ext_mat_field,
                                ConstTensor4View abs_vec_field,
                                const Agenda& surface_rtprop_agenda,
                                const Index& scat_za_interp) {
  Ppath ppath_step;
  // Input variables are checked in the WSMs i_fieldUpdateSeqXXX, from
  // where this function is called.

  //Initialize ppath for 1D.
  //ppath_init_structure(ppath_step, 1, 1);
  ARTS_USER_ERROR("ERROR")
  // See documentation of ppath_init_structure for understanding
  // the parameters.

  // Assign value to ppath.pos:
  ppath_step.pos(0, 0) = z_field(p_index, 0, 0);
  ppath_step.r[0] = refellipsoid[0] + z_field(p_index, 0, 0);

  // Define the direction:
  ppath_step.los(0, 0) = za_grid[za_index];

  // Define the grid positions:
  ppath_step.gp_p[0].idx = p_index;
  ppath_step.gp_p[0].fd[0] = 0;
  ppath_step.gp_p[0].fd[1] = 1;

  // Call ppath_step_agenda:
  ppath_step_agendaExecute(ws,
                           ppath_step,
                           ppath_lmax,
                           ppath_lraytrace,
                           Vector(1, f_grid[f_index]),
                           ppath_step_agenda);

  // Check whether the next point is inside or outside the
  // cloudbox. Only if the next point lies inside the
  // cloudbox a radiative transfer step caclulation has to
  // be performed.

  if ((cloudbox_limits[0] <= ppath_step.gp_p[1].idx &&
       cloudbox_limits[1] > ppath_step.gp_p[1].idx) ||
      (cloudbox_limits[1] == ppath_step.gp_p[1].idx &&
       abs(ppath_step.gp_p[1].fd[0]) < 1e-6)) {
    // Stokes dimension
    const Index stokes_dim = cloudbox_field_mono.ncols();
    // Number of species
    const Index N_species = vmr_field.nbooks();

    // Ppath_step normally has 2 points, the starting
    // point and the intersection point.
    // But there can be points in between, when a maximum
    // lstep is given. We have to interpolate on all the
    // points in the ppath_step.

    // Initialize variables for interpolated values
    Tensor3 ext_mat_int(stokes_dim, stokes_dim, ppath_step.np, 0.);
    Matrix abs_vec_int(stokes_dim, ppath_step.np, 0.);
    Matrix sca_vec_int(stokes_dim, ppath_step.np, 0.);
    Matrix cloudbox_field_mono_int(stokes_dim, ppath_step.np, 0.);
    Vector t_int(ppath_step.np, 0.);
    Matrix vmr_list_int(N_species, ppath_step.np, 0.);
    Vector p_int(ppath_step.np, 0.);

    interp_cloud_coeff1D(ext_mat_int,
                         abs_vec_int,
                         sca_vec_int,
                         cloudbox_field_mono_int,
                         t_int,
                         vmr_list_int,
                         p_int,
                         ext_mat_field,
                         abs_vec_field,
                         doit_scat_field,
                         cloudbox_field_mono_old,
                         t_field,
                         vmr_field,
                         p_grid,
                         ppath_step,
                         cloudbox_limits,
                         za_grid,
                         scat_za_interp);

    // ppath_what_background(ppath_step) tells the
    // radiative background.  More information in the
    // function get_iy_of_background.
    // if there is no background we proceed the RT
    Index bkgr;// = ppath_what_background(ppath_step);
ARTS_USER_ERROR("ERROR")

    // if 0, there is no background
    // do this in any case. cause we need downwelling cloudbox_field_mono
    // at the surface for calculating surface scattering

    // Radiative transfer from one layer to the next, starting
    // at the intersection with the next layer and propagating
    // to the considered point.
    cloud_RT_no_background(ws,
                           cloudbox_field_mono,
                           propmat_clearsky_agenda,
                           ppath_step,
                           t_int,
                           vmr_list_int,
                           ext_mat_int,
                           abs_vec_int,
                           sca_vec_int,
                           cloudbox_field_mono_int,
                           p_int,
                           cloudbox_limits,
                           f_grid,
                           f_index,
                           p_index,
                           0,
                           0,
                           za_index,
                           0);

    if (bkgr == 2) {
      cloud_RT_surface(ws,
                       cloudbox_field_mono,
                       surface_rtprop_agenda,
                       f_grid,
                       f_index,
                       stokes_dim,
                       ppath_step,
                       cloudbox_limits,
                       za_grid,
                       za_index);
    }
  }  //end if inside cloudbox
}

void cloud_ppath_update1D_planeparallel(Workspace& ws,
                                        Tensor6View cloudbox_field_mono,
                                        const Index& p_index,
                                        const Index& za_index,
                                        ConstVectorView za_grid,
                                        const ArrayOfIndex& cloudbox_limits,
                                        ConstTensor6View doit_scat_field,
                                        // Calculate scalar gas absorption:
                                        const Agenda& propmat_clearsky_agenda,
                                        ConstTensor4View vmr_field,
                                        // Propagation path calculation:
                                        ConstVectorView z_grid,
                                        ConstTensor3View p_field,
                                        // Calculate thermal emission:
                                        ConstTensor3View t_field,
                                        ConstVectorView f_grid,
                                        const Index& f_index,
                                        //particle opticla properties
                                        ConstTensor5View ext_mat_field,
                                        ConstTensor4View abs_vec_field) {
  const Index N_species = vmr_field.nbooks();
  const Index stokes_dim = cloudbox_field_mono.ncols();
  const Index atmosphere_dim = 1;
  PropagationMatrix propmat_clearsky;
  PropagationMatrix ext_mat;
  StokesVector abs_vec;
  Matrix matrix_tmp(stokes_dim, stokes_dim);
  Vector vector_tmp(stokes_dim);
  Vector rtp_vmr(N_species, 0.);
  Vector sca_vec_av(stokes_dim, 0);

  // Radiative transfer from one layer to the next, starting
  // at the intersection with the next layer and propagating
  // to the considered point.
  Vector stokes_vec(stokes_dim);
  Index bkgr;
  if ((p_index == 0) && (za_grid[za_index] > 90)) {
    bkgr = 2;
  } else {
    bkgr = 0;
  }

  // if 0, there is no background
  if (bkgr == 0) {
    if (za_grid[za_index] <= 90.0) {
      stokes_vec = cloudbox_field_mono(
          p_index - cloudbox_limits[0] + 1, 0, 0, za_index, 0, joker);
      Numeric z_field_above = z_grid[p_index + 1];
      Numeric z_field_0 = z_grid[p_index];

      Numeric cos_theta, lstep;
      if (za_grid[za_index] == 90.0) {
        //cos_theta = cos((za_grid[za_index] -1)*PI/180.);
        //                  cos_theta = cos(89.999999999999*PI/180.);
        cos_theta = 1e-20;

      } else {
        cos_theta = cos(za_grid[za_index] * PI / 180.0);
      }
      Numeric dz = (z_field_above - z_field_0);

      lstep = abs(dz / cos_theta);

      // Average temperature
      Numeric rtp_temperature =
          0.5 * (t_field(p_index, 0, 0) + t_field(p_index + 1, 0, 0));
      //
      // Average pressure
      Numeric rtp_pressure = 0.5 * (p_field[p_index][0][0] + p_field[p_index + 1][0][0]);

      // Average vmrs
      for (Index j = 0; j < N_species; j++)
        rtp_vmr[j] = 0.5 * (vmr_field(j, p_index, 0, 0) +
                            vmr_field(j, p_index + 1, 0, 0));
      //
      // Calculate scalar gas absorption and add it to abs_vec
      // and ext_mat.
      //

      const Vector rtp_mag_dummy(3, 0);
      const Vector ppath_los_dummy;

      StokesVector nlte_dummy;  //FIXME: do this right?
      ArrayOfPropagationMatrix partial_dummy;  // This is right since there should be only clearsky partials
      ArrayOfStokesVector partial_nlte_dummy;  // This is right since there should be only clearsky partials
      propmat_clearsky_agendaExecute(ws,
                                     propmat_clearsky,
                                     nlte_dummy,
                                     partial_dummy,
                                     partial_nlte_dummy,
                                     ArrayOfRetrievalQuantity(0),
                                     {},
                                     Vector{f_grid[Range(f_index, 1)]},
                                     ppath_los_dummy,
                                     AtmPoint{},  // FIXME: DUMMY VALUE
                                     propmat_clearsky_agenda);

      opt_prop_sum_propmat_clearsky(ext_mat, abs_vec, propmat_clearsky);

      for (Index k = 0; k < stokes_dim; k++) {
        //
        // Averaging of sca_vec:
        //
        sca_vec_av[k] +=
            0.5 *
            (doit_scat_field(
                 p_index - cloudbox_limits[0], 0, 0, za_index, 0, k) +
             doit_scat_field(
                 p_index - cloudbox_limits[0] + 1, 0, 0, za_index, 0, k));
      }

      //
      // Add average particle absorption to abs_vec.
      //
      abs_vec.AddAverageAtPosition(
          abs_vec_field(p_index - cloudbox_limits[0], 0, 0, joker),
          abs_vec_field(p_index - cloudbox_limits[0] + 1, 0, 0, joker));

      //
      // Add average particle extinction to ext_mat.
      ext_mat.AddAverageAtPosition(
          ext_mat_field(p_index - cloudbox_limits[0], 0, 0, joker, joker),
          ext_mat_field(p_index - cloudbox_limits[0] + 1, 0, 0, joker, joker));

      // Frequency
      Numeric f = f_grid[f_index];
      //
      // Calculate Planck function
      //
      Numeric rte_planck_value = planck(f, rtp_temperature);

      // Radiative transfer step calculation. The Stokes vector
      // is updated until the considered point is reached.
      rte_step_doit_replacement(stokes_vec,
                                Matrix(stokes_dim, stokes_dim),
                                ext_mat,
                                abs_vec,
                                sca_vec_av,
                                lstep,
                                rte_planck_value);

      // Assign calculated Stokes Vector to cloudbox_field_mono.
      cloudbox_field_mono(
          p_index - cloudbox_limits[0], 0, 0, za_index, 0, joker) =
          stokes_vec;
    }
    if (za_grid[za_index] > 90) {
      stokes_vec = cloudbox_field_mono(
          p_index - cloudbox_limits[0] - 1, 0, 0, za_index, 0, joker);
      Numeric z_field_below = z_grid[p_index - 1];
      Numeric z_field_0 = z_grid[p_index];

      Numeric cos_theta, lstep;
      if (za_grid[za_index] == 90.0) {
        cos_theta = cos((za_grid[za_index] - 1) * PI / 180.);
        //cos_theta = cos(90.00000001*PI/180.);
        //cout<<"cos_theta"<<cos_theta<<endl;
      } else {
        cos_theta = cos(za_grid[za_index] * PI / 180.0);
      }
      Numeric dz = (z_field_0 - z_field_below);
      lstep = abs(dz / cos_theta);

      // Average temperature
      Numeric rtp_temperature =
          0.5 * (t_field(p_index, 0, 0) + t_field(p_index - 1, 0, 0));
      //
      // Average pressure
      Numeric rtp_pressure = 0.5 * (p_field[p_index][0][0] + p_field[p_index - 1][0][0]);

      //
      // Average vmrs
      for (Index k = 0; k < N_species; k++)
        rtp_vmr[k] = 0.5 * (vmr_field(k, p_index, 0, 0) +
                            vmr_field(k, p_index - 1, 0, 0));
      //
      // Calculate scalar gas absorption and add it to abs_vec
      // and ext_mat.
      //

      const Vector rtp_mag_dummy(3, 0);
      const Vector ppath_los_dummy;
      StokesVector nlte_dummy;  //FIXME: do this right?
      ArrayOfPropagationMatrix partial_dummy;  // This is right since there should be only clearsky partials
      ArrayOfStokesVector partial_nlte_dummy;  // This is right since there should be only clearsky partials
      propmat_clearsky_agendaExecute(ws,
                                     propmat_clearsky,
                                     nlte_dummy,
                                     partial_dummy,
                                     partial_nlte_dummy,
                                     ArrayOfRetrievalQuantity(0),
                                     {},
                                     Vector{f_grid[Range(f_index, 1)]},
                                     ppath_los_dummy,
                                     AtmPoint{},  // FIXME: DUMMY VALUE
                                     propmat_clearsky_agenda);

      opt_prop_sum_propmat_clearsky(ext_mat, abs_vec, propmat_clearsky);

      //
      // Add average particle extinction to ext_mat.
      //

      // cout<<"cloudbox_limits[1]"<<cloudbox_limits[1]<<endl;
      //               cout<<"p_index - cloudbox_limits[0]"<<p_index - cloudbox_limits[0]<<endl;
      for (Index k = 0; k < stokes_dim; k++) {
        //
        // Averaging of sca_vec:
        //
        sca_vec_av[k] +=
            0.5 *
            (doit_scat_field(
                 p_index - cloudbox_limits[0], 0, 0, za_index, 0, k) +
             doit_scat_field(
                 p_index - cloudbox_limits[0] - 1, 0, 0, za_index, 0, k));
      }
      //
      // Add average particle absorption to abs_vec.
      //
      abs_vec.AddAverageAtPosition(
          abs_vec_field(p_index - cloudbox_limits[0], 0, 0, joker),
          abs_vec_field(p_index - cloudbox_limits[0] - 1, 0, 0, joker));

      //
      // Add average particle extinction to ext_mat.
      ext_mat.AddAverageAtPosition(
          ext_mat_field(p_index - cloudbox_limits[0], 0, 0, joker, joker),
          ext_mat_field(p_index - cloudbox_limits[0] + 1, 0, 0, joker, joker));

      // Frequency
      Numeric f = f_grid[f_index];
      //
      // Calculate Planck function
      //
      Numeric rte_planck_value = planck(f, rtp_temperature);

      // Radiative transfer step calculation. The Stokes vector
      // is updated until the considered point is reached.
      rte_step_doit_replacement(stokes_vec,
                                Matrix(stokes_dim, stokes_dim),
                                ext_mat,
                                abs_vec,
                                sca_vec_av,
                                lstep,
                                rte_planck_value);

      // Assign calculated Stokes Vector to cloudbox_field_mono.
      cloudbox_field_mono(
          p_index - cloudbox_limits[0], 0, 0, za_index, 0, joker) =
          stokes_vec;
    }

  }  // if loop end - for non_ground background

  // bkgr=2 indicates that the background is the surface
  else if (bkgr == 2) {
    //Set rte_pos, rte_gp_p and rte_los to match the last point
    //in ppath.
    //pos
    Vector rte_pos(atmosphere_dim);
    Numeric z_field_0 = z_grid[0];
    rte_pos = z_field_0;  //ppath_step.pos(np-1,Range(0,atmosphere_dim));
    //los
    Vector rte_los(1);
    rte_los = za_grid[za_index];  //ppath_step.los(np-1,joker);
    //gp_p
    //GridPos rte_gp_p;
    //rte_gp_p.idx   = p_index;
    //rte_gp_p.fd[0] = 0;
    //rte_gp_p.fd[1] = 1;
    //gridpos_copy( rte_gp_p, ppath_step.gp_p[np-1] );
    // Executes the surface agenda
    // FIXME: Convert to new agenda scheme before using
    // surface_rtprop_agenda.execute();

    ARTS_USER_ERROR (
        "Surface reflections inside cloud box not yet handled.");
    /*
        See comment in function above
          // Check returned variables
          if( surface_emission.nrows() != f_grid.nelem()  ||
              surface_emission.ncols() != stokes_dim )
            throw runtime_error(
                  "The size of the created *surface_emission* is not correct.");

          Index nlos = surface_los.nrows();

          // Define a local vector cloudbox_field_mono_sum which adds the
          // products of groudnd_refl_coeffs with the downwelling
          // radiation for each elements of surface_los
          Vector cloudbox_field_mono_sum(stokes_dim,0);
          // Loop over the surface_los elements
          for( Index ilos=0; ilos < nlos; ilos++ )
            {
              if( stokes_dim == 1 )
                {
                  cloudbox_field_mono_sum[0] += surface_refl_coeffs(ilos,f_index,0,0) *
                    cloudbox_field_mono(cloudbox_limits[0],
                            0, 0,
                            (za_grid.nelem() -1 - za_index), 0,
                            0);
                }
              else
                {
                  Vector stokes_vec2(stokes_dim);
                  mult( stokes_vec2,
                        surface_refl_coeffs(ilos,0,joker,joker),
                        cloudbox_field_mono(cloudbox_limits[0],
                                0, 0,
                                (za_grid.nelem() -1 - za_index), 0,
                                joker));
                  for( Index is=0; is < stokes_dim; is++ )
                    {
                      cloudbox_field_mono_sum[is] += stokes_vec2[is];
                    }

                }
            }
          // Copy from *cloudbox_field_mono_sum* to *cloudbox_field_mono*, and add the surface emission
          for( Index is=0; is < stokes_dim; is++ )
            {
              cloudbox_field_mono (cloudbox_limits[0],
                       0, 0,
                       za_index, 0,
                       is) = cloudbox_field_mono_sum[is] + surface_emission(f_index,is);
            }

          //cout<<"za_grid"<<za_grid[za_index]<<endl;
          //cout<<"p_index"<<p_index<<endl;
          //cout<<"cloudboxlimit[0]"<<cloudbox_limits[0]<<endl;
          // now the RT is done to the next point in the path.
          //
          Vector stokes_vec_local;
          stokes_vec_local = cloudbox_field_mono (0,
                                      0, 0,
                                      za_index, 0,
                                      joker);

          Numeric z_field_above = z_field(p_index +1, 0, 0);
          //Numeric z_field_0 = z_field(p_index, 0, 0);
          Numeric cos_theta;
          if (za_grid[za_index] == 90.0)
            {
              //cos_theta = cos((za_grid[za_index] -1)*PI/180.);
              cos_theta = cos(90.00000001*PI/180.);
            }
          else
            {
              cos_theta = cos(za_grid[za_index]* PI/180.0);
            }
          Numeric dz = (z_field_above -  z_field_0);

          Numeric lstep =  abs(dz/cos_theta) ;

          // Average temperature
          Numeric rtp_temperature =   0.5 * (t_field(p_index,0,0)
                                             + t_field(p_index + 1,0,0));

          //
          // Average pressure
          Numeric rtp_pressure = 0.5 * (p_grid[p_index]
                                        + p_grid[p_index + 1]);

          //
          const Index N_species = vmr_field.nbooks();
          // Average vmrs
          for (Index k = 0; k < N_species; k++)
            {
              rtp_vmr[k] = 0.5 * (vmr_field(k,p_index,0,0) +
                                       vmr_field(k,p_index + 1,0,0));
            }
          //
          // Calculate scalar gas absorption and add it to abs_vec
          // and ext_mat.
          //

          // FIXME: Convert to new agenda scheme before using
          //abs_scalar_gas_agenda.execute(p_index);

          opt_prop_gas_agendaExecute(ext_mat, abs_vec, abs_scalar_gas,
                                     opt_prop_gas_agenda);

          //
          // Add average particle extinction to ext_mat.
          //
          //cout<<"Reached Here ????????????????????????????????????????????????";
          for (Index k = 0; k < stokes_dim; k++)
            {
              for (Index j = 0; j < stokes_dim; j++)
                {
                  ext_mat(0,k,j) += 0.5 *
                    (ext_mat_field(p_index - cloudbox_limits[0],0,0,k,j) +
                     ext_mat_field(p_index - cloudbox_limits[0]+ 1,0,0,k,j));
                }


              //
              //
              // Add average particle absorption to abs_vec.
              //
              abs_vec(0,k) += 0.5 *
                (abs_vec_field(p_index - cloudbox_limits[0],0,0,k) +
                 abs_vec_field(p_index - cloudbox_limits[0]+1,0,0,k));
              //
              // Averaging of sca_vec:
              //
              sca_vec_av[k] = 0.5*( doit_scat_field(p_index- cloudbox_limits[0],
                                               0, 0, za_index, 0, k)
                                    + doit_scat_field(p_index- cloudbox_limits[0]+1,
                                                 0, 0, za_index, 0, k)) ;

            }
          // Frequency
          Numeric f = f_grid[f_index];
          //
          // Calculate Planck function
          //
          Numeric rte_planck_value = planck(f, rtp_temperature);

          ARTS_ASSERT (!is_singular( ext_mat(0,joker,joker)));

          // Radiative transfer step calculation. The Stokes vector
          // is updated until the considered point is reached.
          rte_step_doit(stokes_vec_local, ext_mat(0,joker,joker),
                   abs_vec(0,joker),
                   sca_vec_av, lstep, rte_planck_value);
          // Assign calculated Stokes Vector to cloudbox_field_mono.
          cloudbox_field_mono(p_index - cloudbox_limits[0],
                  0, 0,
                  za_index, 0,
                  joker) = stokes_vec_local;
      */
  }  //end else loop over surface
}

void cloud_ppath_update3D(Workspace& ws,
                          Tensor6View cloudbox_field_mono,
                          // ppath_step_agenda:
                          const Index& p_index,
                          const Index& lat_index,
                          const Index& lon_index,
                          const Index& za_index,
                          const Index& aa_index,
                          ConstVectorView za_grid,
                          ConstVectorView aa_grid,
                          const ArrayOfIndex& cloudbox_limits,
                          ConstTensor6View doit_scat_field,
                          // Calculate scalar gas absorption:
                          const Agenda& propmat_clearsky_agenda,
                          const AtmField& atm_field,
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          // Propagation path calculation:
                          const Agenda& ppath_step_agenda,
                          const Numeric& ppath_lmax,
                          const Numeric& ppath_lraytrace,
                          ConstVectorView refellipsoid,
                          // Calculate thermal emission:
                          ConstVectorView f_grid,
                          const Index& f_index,
                          //particle optical properties
                          ConstTensor5View ext_mat_field,
                          ConstTensor4View abs_vec_field,
                          const Index&) {
//FIXME: MUST HAVE REGULAR FIELD
 Vector z_grid, lat_grid, lon_grid;
 Tensor3 p_field, t_field;
 Tensor4 vmr_field;
//  ARTS_USER_ERROR_IF(not atm_field.regularized, "Must have regular grid field")
  //const Vector& z_grid = atm_field.grid[0];
 // const Vector& lat_grid = atm_field.grid[1];
 // const Vector& lon_grid = atm_field.grid[2];
 // const auto& t_field = atm_field[Atm::Key::t].get<const Tensor3&>();
 // const auto& p_field = atm_field[Atm::Key::p].get<const Tensor3&>();
 // const Tensor4 vmr_field = Atm::extract_specs_content(atm_field, abs_species);

  Ppath ppath_step;
  const Index stokes_dim = cloudbox_field_mono.ncols();

  Vector sca_vec_av(stokes_dim, 0);
  Vector aa_g(aa_grid.nelem());

  for (Index i = 0; i < aa_grid.nelem(); i++)
    aa_g[i] = aa_grid[i] - 180.;

  //Initialize ppath for 3D.
  //ppath_init_structure(ppath_step, 3, 1);
  ARTS_USER_ERROR("ERROR")
  // See documentation of ppath_init_structure for
  // understanding the parameters.

  // The first dimension of pos are the points in
  // the propagation path.
  // Here we initialize the first point.
  // The second is: radius, latitude, longitude

  ppath_step.pos(0, 2) = lon_grid[lon_index];
  ppath_step.pos(0, 1) = lat_grid[lat_index];
  ppath_step.pos(0, 0) = z_grid[p_index];
  // As always on top of the lat. grid positions, OK to call refell2r:
 // ppath_step.r[0] =
   //   refell2r(refellipsoid, ppath_step.pos(0, 1)) + ppath_step.pos(0, 0);
ARTS_USER_ERROR("ERROR")

  // Define the direction:
  ppath_step.los(0, 0) = za_grid[za_index];
  ppath_step.los(0, 1) = aa_g[aa_index];

  // Define the grid positions:
  ppath_step.gp_p[0].idx = p_index;
  ppath_step.gp_p[0].fd[0] = 0.;
  ppath_step.gp_p[0].fd[1] = 1.;

  ppath_step.gp_lat[0].idx = lat_index;
  ppath_step.gp_lat[0].fd[0] = 0.;
  ppath_step.gp_lat[0].fd[1] = 1.;

  ppath_step.gp_lon[0].idx = lon_index;
  ppath_step.gp_lon[0].fd[0] = 0.;
  ppath_step.gp_lon[0].fd[1] = 1.;

  // Call ppath_step_agenda:
  ppath_step_agendaExecute(ws,
                           ppath_step,
                           ppath_lmax,
                           ppath_lraytrace,
                           Vector(1, f_grid[f_index]),
                           ppath_step_agenda);

  // Check whether the next point is inside or outside the
  // cloudbox. Only if the next point lies inside the
  // cloudbox a radiative transfer step caclulation has to
  // be performed.
  if (is_inside_cloudbox(ppath_step, cloudbox_limits, true)) {
    // Gridpositions inside the cloudbox.
    // The optical properties are stored only inside the
    // cloudbox. For interpolation we use grids
    // inside the cloudbox.

    ArrayOfGridPos cloud_gp_p = ppath_step.gp_p;
    ArrayOfGridPos cloud_gp_lat = ppath_step.gp_lat;
    ArrayOfGridPos cloud_gp_lon = ppath_step.gp_lon;

    for (Index i = 0; i < ppath_step.np; i++) {
      cloud_gp_p[i].idx -= cloudbox_limits[0];
      cloud_gp_lat[i].idx -= cloudbox_limits[2];
      cloud_gp_lon[i].idx -= cloudbox_limits[4];
    }
    const Index n1 = cloudbox_limits[1] - cloudbox_limits[0];
    const Index n2 = cloudbox_limits[3] - cloudbox_limits[2];
    const Index n3 = cloudbox_limits[5] - cloudbox_limits[4];
    gridpos_upperend_check(cloud_gp_p[0], n1);
    gridpos_upperend_check(cloud_gp_p[ppath_step.np - 1], n1);
    gridpos_upperend_check(cloud_gp_lat[0], n2);
    gridpos_upperend_check(cloud_gp_lat[ppath_step.np - 1], n2);
    gridpos_upperend_check(cloud_gp_lon[0], n3);
    gridpos_upperend_check(cloud_gp_lon[ppath_step.np - 1], n3);

    Matrix itw(ppath_step.np, 8);
    interpweights(itw, cloud_gp_p, cloud_gp_lat, cloud_gp_lon);

    Matrix itw_p(ppath_step.np, 2);
    interpweights(itw_p, cloud_gp_p);

    // The zenith angles and azimuth of the propagation path are
    // needed as we have to
    // interpolate the intensity field and the scattered field on the
    // right angles.
    VectorView los_grid_za = ppath_step.los(joker, 0);
    VectorView los_grid_aa = ppath_step.los(joker, 1);

    for (Index i = 0; i < los_grid_aa.nelem(); i++)
      los_grid_aa[i] = los_grid_aa[i] + 180.;

    ArrayOfGridPos gp_za(los_grid_za.nelem());
    gridpos(gp_za, za_grid, los_grid_za);

    ArrayOfGridPos gp_aa(los_grid_aa.nelem());
    gridpos(gp_aa, aa_grid, los_grid_aa);

    Matrix itw_p_za(ppath_step.np, 32);
    interpweights(
        itw_p_za, cloud_gp_p, cloud_gp_lat, cloud_gp_lon, gp_za, gp_aa);

    // Ppath_step normally has 2 points, the starting
    // point and the intersection point.
    // But there can be points in between, when a maximum
    // lstep is given. We have to interpolate on all the
    // points in the ppath_step.

    Tensor3 ext_mat_int(stokes_dim, stokes_dim, ppath_step.np);
    Matrix abs_vec_int(stokes_dim, ppath_step.np);
    Matrix sca_vec_int(stokes_dim, ppath_step.np, 0.);
    Matrix cloudbox_field_mono_int(stokes_dim, ppath_step.np, 0.);
    Vector t_int(ppath_step.np);
    Vector vmr_int(ppath_step.np);
    Vector p_int(ppath_step.np);
    Vector stokes_vec(stokes_dim);
    //Tensor3 ext_mat_gas(stokes_dim, stokes_dim, ppath_step.np);
    //Matrix abs_vec_gas(stokes_dim, ppath_step.np);

    // Calculate the average of the coefficients for the layers
    // to be considered in the
    // radiative transfer calculation.

    for (Index i = 0; i < stokes_dim; i++) {
      // Extinction matrix requires a second loop
      // over stokes_dim
      for (Index j = 0; j < stokes_dim; j++) {
        //
        // Interpolation of ext_mat
        //
        interp(ext_mat_int(i, j, joker),
               itw,
               ext_mat_field(joker, joker, joker, i, j),
               cloud_gp_p,
               cloud_gp_lat,
               cloud_gp_lon);
      }
      // Absorption vector:
      //
      // Interpolation of abs_vec
      //
      interp(abs_vec_int(i, joker),
             itw,
             abs_vec_field(joker, joker, joker, i),
             cloud_gp_p,
             cloud_gp_lat,
             cloud_gp_lon);
      //
      // Scattered field:
      //
      // Interpolation of sca_vec:
      //
      interp(sca_vec_int(i, joker),
             itw_p_za,
             doit_scat_field(joker, joker, joker, joker, joker, i),
             cloud_gp_p,
             cloud_gp_lat,
             cloud_gp_lon,
             gp_za,
             gp_aa);
      interp(cloudbox_field_mono_int(i, joker),
             itw_p_za,
             cloudbox_field_mono(joker, joker, joker, joker, joker, i),
             cloud_gp_p,
             cloud_gp_lat,
             cloud_gp_lon,
             gp_za,
             gp_aa);
    }
    //
    // Planck function
    //
    // Interpolate temperature field
    //
    interp(t_int,
           itw,
           t_field(joker, joker, joker),
           ppath_step.gp_p,
           ppath_step.gp_lat,
           ppath_step.gp_lon);

    //
    // The vmr_field is needed for the gaseous absorption
    // calculation.
    //
    const Index N_species = vmr_field.nbooks();
    //
    // Interpolated vmr_list, holds a vmr_list for each point in
    // ppath_step.
    //
    Matrix vmr_list_int(N_species, ppath_step.np);

    for (Index i = 0; i < N_species; i++) {
      interp(vmr_int,
             itw,
             vmr_field(i, joker, joker, joker),
             ppath_step.gp_p,
             ppath_step.gp_lat,
             ppath_step.gp_lon);

      vmr_list_int(i, joker) = vmr_int;
    }

    // Presssure (needed for the calculation of gas absorption)
    itw2p(p_int, VectorView{p_field(p_index, lat_index, lon_index)}, ppath_step.gp_p, itw_p);

    cloud_RT_no_background(ws,
                           cloudbox_field_mono,
                           propmat_clearsky_agenda,
                           ppath_step,
                           t_int,
                           vmr_list_int,
                           ext_mat_int,
                           abs_vec_int,
                           sca_vec_int,
                           cloudbox_field_mono_int,
                           p_int,
                           cloudbox_limits,
                           f_grid,
                           f_index,
                           p_index,
                           lat_index,
                           lon_index,
                           za_index,
                           aa_index);
  }  //end if inside cloudbox
}

void cloud_RT_no_background(Workspace& ws,
                            //Output
                            Tensor6View cloudbox_field_mono,
                            // Input
                            const Agenda& propmat_clearsky_agenda,
                            const Ppath& ppath_step,
                            ConstVectorView t_int,
                            ConstMatrixView vmr_list_int,
                            ConstTensor3View ext_mat_int,
                            ConstMatrixView abs_vec_int,
                            ConstMatrixView sca_vec_int,
                            ConstMatrixView cloudbox_field_mono_int,
                            ConstVectorView p_int,
                            const ArrayOfIndex& cloudbox_limits,
                            ConstVectorView f_grid,
                            const Index& f_index,
                            const Index& p_index,
                            const Index& lat_index,
                            const Index& lon_index,
                            const Index& za_index,
                            const Index& aa_index) {
  const Index N_species = vmr_list_int.nrows();
  const Index stokes_dim = cloudbox_field_mono.ncols();
  const Index atmosphere_dim = cloudbox_limits.nelem() / 2;

  Vector sca_vec_av(stokes_dim, 0);
  Vector stokes_vec(stokes_dim, 0.);
  Vector rtp_vmr_local(N_species, 0.);

  // Two propmat_clearsky to average between
  PropagationMatrix cur_propmat_clearsky;
  PropagationMatrix prev_propmat_clearsky;

  PropagationMatrix ext_mat_local;
  StokesVector abs_vec_local;
  Matrix matrix_tmp(stokes_dim, stokes_dim);
  Vector vector_tmp(stokes_dim);

  // Incoming stokes vector
  stokes_vec = cloudbox_field_mono_int(joker, ppath_step.np - 1);

  for (Index k = ppath_step.np - 1; k >= 0; k--) {
    // Save propmat_clearsky from previous level by
    // swapping it with current level
    swap(cur_propmat_clearsky, prev_propmat_clearsky);

    //
    // Calculate scalar gas absorption
    //
    const Vector rtp_mag_dummy(3, 0);
    const Vector ppath_los_dummy;

    StokesVector nlte_dummy;  //FIXME: do this right?
    ArrayOfPropagationMatrix partial_dummy;  // This is right since there should be only clearsky partials
    ArrayOfStokesVector partial_nlte_dummy;  // This is right since there should be only clearsky partials
    propmat_clearsky_agendaExecute(ws,
                                   cur_propmat_clearsky,
                                   nlte_dummy,
                                   partial_dummy,
                                   partial_nlte_dummy,
                                   ArrayOfRetrievalQuantity(0),
                                   {},
                                   Vector{f_grid[Range(f_index, 1)]},
                                   ppath_los_dummy,
                                   AtmPoint{},  // FIXME: DUMMY VALUE
                                   propmat_clearsky_agenda);

    // Skip any further calculations for the first point.
    // We need values at two ppath points before we can average.
    if (k == ppath_step.np - 1) continue;

    // Average prev_propmat_clearsky with cur_propmat_clearsky
    prev_propmat_clearsky += cur_propmat_clearsky;
    prev_propmat_clearsky *= 0.5;

    opt_prop_sum_propmat_clearsky(
        ext_mat_local, abs_vec_local, prev_propmat_clearsky);

    for (Index i = 0; i < stokes_dim; i++) {
      //
      // Averaging of sca_vec:
      //
      sca_vec_av[i] = 0.5 * (sca_vec_int(i, k) + sca_vec_int(i, k + 1));
    }

    //
    // Add average particle absorption to abs_vec.
    //
    abs_vec_local.AddAverageAtPosition(abs_vec_int(joker, k),
                                       abs_vec_int(joker, k + 1));

    //
    // Add average particle extinction to ext_mat.
    //
    ext_mat_local.AddAverageAtPosition(ext_mat_int(joker, joker, k),
                                       ext_mat_int(joker, joker, k + 1));

    // Frequency
    Numeric f = f_grid[f_index];
    //
    // Calculate Planck function
    //
    Numeric rte_planck_value = planck(f, 0.5 * (t_int[k] + t_int[k + 1]));

    // Length of the path between the two layers.
    Numeric lstep = ppath_step.lstep[k];

    // Radiative transfer step calculation. The Stokes vector
    // is updated until the considered point is reached.
    rte_step_doit_replacement(stokes_vec,
                              Matrix(stokes_dim, stokes_dim),
                              ext_mat_local,
                              abs_vec_local,
                              sca_vec_av,
                              lstep,
                              rte_planck_value);

  }  // End of loop over ppath_step.
  // Assign calculated Stokes Vector to cloudbox_field_mono.
  if (atmosphere_dim == 1)
    cloudbox_field_mono(
        p_index - cloudbox_limits[0], 0, 0, za_index, 0, joker) =
        stokes_vec;
  else if (atmosphere_dim == 3)
    cloudbox_field_mono(p_index - cloudbox_limits[0],
                      lat_index - cloudbox_limits[2],
                      lon_index - cloudbox_limits[4],
                      za_index,
                      aa_index,
                      joker) = stokes_vec;
}

void cloud_RT_surface(Workspace& ws,
                      //Output
                      Tensor6View cloudbox_field_mono,
                      //Input
                      const Agenda& surface_rtprop_agenda,
                      ConstVectorView f_grid,
                      const Index& f_index,
                      const Index& stokes_dim,
                      const Ppath& ppath_step,
                      const ArrayOfIndex& cloudbox_limits,
                      ConstVectorView za_grid,
                      const Index& za_index) {
  chk_not_empty("surface_rtprop_agenda", surface_rtprop_agenda);

  Matrix iy;

  // Local output of surface_rtprop_agenda.
  SurfacePoint surface_point;
  Matrix surface_emission;
  Matrix surface_los;
  Tensor4 surface_rmatrix;

  //Set rte_pos and rte_los to match the last point in ppath.

  Index np = ppath_step.np;

  Vector rte_pos;  // ppath_step.pos contains two columns for 1D
  rte_pos.resize(ppath_step.dim);
  rte_pos = ppath_step.pos(np - 1, Range(0, ppath_step.dim));

  Vector rte_los;
  rte_los.resize(ppath_step.los.ncols());
  rte_los = ppath_step.los(np - 1, joker);

  //Execute the surface_rtprop_agenda which gives the surface
  //parameters.
  surface_rtprop_agendaExecute(ws,
                               surface_point,
                               surface_emission,
                               surface_los,
                               surface_rmatrix,
                               Vector(1, f_grid[f_index]),
                               rte_pos,
                               rte_los,
                               surface_rtprop_agenda);

  iy = surface_emission;

  Index nlos = surface_los.nrows();

  if (nlos > 0) {
    Vector rtmp(stokes_dim);  // Reflected Stokes vector for 1 frequency

    for (Index ilos = 0; ilos < nlos; ilos++) {
      // Several things needs to be fixed here. As far as I understand it,
      // this works only for specular cases and if the lower cloudbox limit
      // is exactly at the surface (PE, 120828)

      mult(rtmp,
           surface_rmatrix(ilos, 0, joker, joker),
           cloudbox_field_mono(cloudbox_limits[0],
                             0,
                             0,
                             (za_grid.nelem() - 1 - za_index),
                             0,
                             joker));
      iy(0, joker) += rtmp;
    }
  }
  cloudbox_field_mono(cloudbox_limits[0], 0, 0, za_index, 0, joker) =
      iy(0, joker);
}

void cloudbox_field_ngAcceleration(Tensor6& cloudbox_field_mono,
                                 const ArrayOfTensor6& acceleration_input,
                                 const Index& accelerated) {
  const Index N_p = cloudbox_field_mono.nvitrines();
  const Index N_za = cloudbox_field_mono.npages();

  // Loop over 4 components of Stokes Vector
  for (Index i = 0; i < accelerated; ++i) {
    ConstMatrixView S1 = acceleration_input[0](joker, 0, 0, joker, 0, i);
    ConstMatrixView S2 = acceleration_input[1](joker, 0, 0, joker, 0, i);
    ConstMatrixView S3 = acceleration_input[2](joker, 0, 0, joker, 0, i);
    ConstMatrixView S4 = acceleration_input[3](joker, 0, 0, joker, 0, i);

    ConstMatrixView J = S4;
    Matrix Q1;
    Matrix Q2;
    Matrix Q3;
    Numeric A1 = 0;
    Numeric A2B1 = 0;
    Numeric B2 = 0;
    Numeric C1 = 0;
    Numeric C2 = 0;
    Numeric NGA = 0;
    Numeric NGB = 0;

    // Q1 = -2*S3 + S4 + S2

    Q1 = S3;
    Q1 *= -2;
    Q1 += S4;
    Q1 += S2;

    // Q2 = S4 - S3 - S2 + S1
    Q2 = S4;
    Q2 -= S3;
    Q2 -= S2;
    Q2 += S1;

    //Q3 = S4 - S3
    Q3 = S4;
    Q3 -= S3;

    for (Index p_index = 0; p_index < N_p; ++p_index) {
      for (Index za_index = 0; za_index < N_za; ++za_index) {
        A1 += Q1(p_index, za_index) * Q1(p_index, za_index) *
              J(p_index, za_index);
        A2B1 += Q2(p_index, za_index) * Q1(p_index, za_index) *
                J(p_index, za_index);
        B2 += Q2(p_index, za_index) * Q2(p_index, za_index) *
              J(p_index, za_index);
        C1 += Q1(p_index, za_index) * Q3(p_index, za_index) *
              J(p_index, za_index);
        C2 += Q2(p_index, za_index) * Q3(p_index, za_index) *
              J(p_index, za_index);
      }
    }

    NGA = (C1 * B2 - C2 * A2B1) / (A1 * B2 - A2B1 * A2B1);
    NGB = (C2 * A1 - C1 * A2B1) / (A1 * B2 - A2B1 * A2B1);

    if (!std::isnan(NGB) && !std::isnan(NGA)) {
      // Calculating the accelerated field
      for (Index p_index = 0; p_index < N_p; ++p_index) {
        for (Index za_index = 0; za_index < N_za; ++za_index) {
          Q1(p_index, za_index) = (1 - NGA - NGB) * S4(p_index, za_index) +
                                  NGA * S3(p_index, za_index) +
                                  NGB * S2(p_index, za_index);
        }
      }
      cloudbox_field_mono(joker, 0, 0, joker, 0, i) = Q1;
    }
  }
}

void interp_cloud_coeff1D(  //Output
    Tensor3View ext_mat_int,
    MatrixView abs_vec_int,
    MatrixView sca_vec_int,
    MatrixView cloudbox_field_mono_int,
    VectorView t_int,
    MatrixView vmr_list_int,
    VectorView p_int,
    //Input
    ConstTensor5View ext_mat_field,
    ConstTensor4View abs_vec_field,
    ConstTensor6View doit_scat_field,
    ConstTensor6View cloudbox_field_mono,
    ConstTensor3View t_field,
    ConstTensor4View vmr_field,
    ConstVectorView p_grid,
    const Ppath& ppath_step,
    const ArrayOfIndex& cloudbox_limits,
    ConstVectorView za_grid,
    const Index& scat_za_interp) {
  // Stokes dimension
  const Index stokes_dim = cloudbox_field_mono.ncols();

  // Gridpositions inside the cloudbox.
  // The optical properties are stored only inside the
  // cloudbox. For interpolation we use grids
  // inside the cloudbox.
  ArrayOfGridPos cloud_gp_p = ppath_step.gp_p;

  for (Index i = 0; i < ppath_step.np; i++)
    cloud_gp_p[i].idx -= cloudbox_limits[0];

  // Grid index for points at upper limit of cloudbox must be shifted
  const Index n1 = cloudbox_limits[1] - cloudbox_limits[0];
  gridpos_upperend_check(cloud_gp_p[0], n1);
  gridpos_upperend_check(cloud_gp_p[ppath_step.np - 1], n1);

  Matrix itw(cloud_gp_p.nelem(), 2);
  interpweights(itw, cloud_gp_p);

  // The zenith angles of the propagation path are needed as we have to
  // interpolate the intensity field and the scattered field on the
  // right angles.
  Vector los_grid{ppath_step.los(joker, 0)};

  ArrayOfGridPos gp_za(los_grid.nelem());
  gridpos(gp_za, za_grid, los_grid);

  Matrix itw_p_za(cloud_gp_p.nelem(), 4);
  interpweights(itw_p_za, cloud_gp_p, gp_za);

  // Calculate the average of the coefficients for the layers
  // to be considered in the
  // radiative transfer calculation.

  for (Index i = 0; i < stokes_dim; i++) {
    // Extinction matrix requires a second loop
    // over stokes_dim
    for (Index j = 0; j < stokes_dim; j++) {
      //
      // Interpolation of ext_mat
      //
      interp(ext_mat_int(i, j, joker),
             itw,
             ext_mat_field(joker, 0, 0, i, j),
             cloud_gp_p);
    }
    // Particle absorption vector:
    //
    // Interpolation of abs_vec
    //  //
    interp(
        abs_vec_int(i, joker), itw, abs_vec_field(joker, 0, 0, i), cloud_gp_p);
    //
    // Scattered field:
    //
    //

    if (scat_za_interp == 0)  // linear interpolation
    {
      interp(sca_vec_int(i, joker),
             itw_p_za,
             doit_scat_field(joker, 0, 0, joker, 0, i),
             cloud_gp_p,
             gp_za);
      interp(cloudbox_field_mono_int(i, joker),
             itw_p_za,
             cloudbox_field_mono(joker, 0, 0, joker, 0, i),
             cloud_gp_p,
             gp_za);
    } else if (scat_za_interp == 1)  //polynomial interpolation
    {
      // These intermediate variables are needed because polynomial
      // interpolation is not implemented as multidimensional
      // interpolation.
      Tensor3 sca_vec_int_za(
          stokes_dim, ppath_step.np, za_grid.nelem(), 0.);
      Tensor3 cloudbox_field_mono_int_za(
          stokes_dim, ppath_step.np, za_grid.nelem(), 0.);
      for (Index za = 0; za < za_grid.nelem(); za++) {
        interp(sca_vec_int_za(i, joker, za),
               itw,
               doit_scat_field(joker, 0, 0, za, 0, i),
               cloud_gp_p);
        interp(cloudbox_field_mono_int_za(i, joker, za),
               itw,
               cloudbox_field_mono(joker, 0, 0, za, 0, i),
               cloud_gp_p);
      }
      for (Index ip = 0; ip < ppath_step.np; ip++) {
        sca_vec_int(i, ip) = interp_poly(za_grid,
                                         sca_vec_int_za(i, ip, joker),
                                         los_grid[ip],
                                         gp_za[ip]);
        cloudbox_field_mono_int(i, ip) =
            interp_poly(za_grid,
                        cloudbox_field_mono_int_za(i, ip, joker),
                        los_grid[ip],
                        gp_za[ip]);
      }
    }
  }
  //
  // Planck function
  //
  // Interpolate temperature field
  //
  interp(t_int, itw, t_field(joker, 0, 0), ppath_step.gp_p);
  //
  // The vmr_field is needed for the gaseous absorption
  // calculation.
  //
  const Index N_species = vmr_field.nbooks();
  //
  // Interpolated vmr_list, holds a vmr_list for each point in
  // ppath_step.
  //
  Vector vmr_int(ppath_step.np);

  for (Index i_sp = 0; i_sp < N_species; i_sp++) {
    interp(vmr_int, itw, vmr_field(i_sp, joker, 0, 0), ppath_step.gp_p);
    vmr_list_int(i_sp, joker) = vmr_int;
  }
  //
  // Interpolate pressure
  //
  itw2p(p_int, p_grid, ppath_step.gp_p, itw);
}

void za_gridOpt(  //Output:
    Vector& za_grid_opt,
    Matrix& cloudbox_field_opt,
    // Input
    ConstVectorView za_grid_fine,
    ConstTensor6View cloudbox_field_mono,
    const Numeric& acc,
    const Index& scat_za_interp) {
  Index N_za = za_grid_fine.nelem();

  ARTS_ASSERT(cloudbox_field_mono.npages() == N_za);

  Index N_p = cloudbox_field_mono.nvitrines();

  Vector i_approx_interp(N_za);
  Vector za_reduced(2);

  ArrayOfIndex idx;
  idx.push_back(0);
  idx.push_back(N_za - 1);
  ArrayOfIndex idx_unsorted;

  Numeric max_diff = 100;

  ArrayOfGridPos gp_za(N_za);
  Matrix itw(za_grid_fine.nelem(), 2);

  ArrayOfIndex i_sort;
  Vector diff_vec(N_za);
  Vector max_diff_za(N_p);
  ArrayOfIndex ind_za(N_p);
  Numeric max_diff_p;
  Index ind_p = 0;

  while (max_diff > acc) {
    za_reduced.resize(idx.nelem());
    cloudbox_field_opt.resize(N_p, idx.nelem());
    max_diff_za = 0.;
    max_diff_p = 0.;

    // Interpolate reduced intensity field on fine za_grid for
    // all pressure levels
    for (Index i_p = 0; i_p < N_p; i_p++) {
      for (Index i_za_red = 0; i_za_red < idx.nelem(); i_za_red++) {
        za_reduced[i_za_red] = za_grid_fine[idx[i_za_red]];
        cloudbox_field_opt(i_p, i_za_red) =
            cloudbox_field_mono(i_p, 0, 0, idx[i_za_red], 0, 0);
      }
      // Calculate grid positions
      gridpos(gp_za, za_reduced, za_grid_fine);
      //linear interpolation
      if (scat_za_interp == 0 || idx.nelem() < 3) {
        interpweights(itw, gp_za);
        interp(i_approx_interp, itw, cloudbox_field_opt(i_p, joker), gp_za);
      } else if (scat_za_interp == 1) {
        for (Index i_za = 0; i_za < N_za; i_za++) {
          i_approx_interp[i_za] = interp_poly(za_reduced,
                                              cloudbox_field_opt(i_p, joker),
                                              za_grid_fine[i_za],
                                              gp_za[i_za]);
        }
      } else
        // Interpolation method not defined
        ARTS_ASSERT(false);

      // Calculate differences between approximated i-vector and
      // exact i_vector for the i_p pressure level
      for (Index i_za = 0; i_za < N_za; i_za++) {
        diff_vec[i_za] = abs(cloudbox_field_mono(i_p, 0, 0, i_za, 0, 0) -
                             i_approx_interp[i_za]);
        if (diff_vec[i_za] > max_diff_za[i_p]) {
          max_diff_za[i_p] = diff_vec[i_za];
          ind_za[i_p] = i_za;
        }
      }
      // Take maximum value of max_diff_za
      if (max_diff_za[i_p] > max_diff_p) {
        max_diff_p = max_diff_za[i_p];
        ind_p = i_p;
      }
    }

    //Transform in %
    max_diff =
        max_diff_p / cloudbox_field_mono(ind_p, 0, 0, ind_za[ind_p], 0, 0) * 100.;

    idx.push_back(ind_za[ind_p]);
    idx_unsorted = idx;

    i_sort.resize(idx_unsorted.nelem());
    get_sorted_indexes(i_sort, idx_unsorted);

    for (Index i = 0; i < idx_unsorted.nelem(); i++)
      idx[i] = idx_unsorted[i_sort[i]];

    za_reduced.resize(idx.nelem());
  }

  za_grid_opt.resize(idx.nelem());
  cloudbox_field_opt.resize(N_p, idx.nelem());
  for (Index i = 0; i < idx.nelem(); i++) {
    za_grid_opt[i] = za_grid_fine[idx[i]];
    cloudbox_field_opt(joker, i) = cloudbox_field_mono(joker, 0, 0, idx[i], 0, 0);
  }
}

void doit_scat_fieldNormalize(Workspace& ws,
                              Tensor6& doit_scat_field,
                              const Tensor6& cloudbox_field_mono,
                              const ArrayOfIndex& cloudbox_limits,
                              const Agenda& spt_calc_agenda,
                              const Index& atmosphere_dim,
                              const Vector& za_grid,
                              const Vector& aa_grid,
                              const Tensor4& pnd_field,
                              const Tensor3& t_field,
                              const Numeric& norm_error_threshold,
                              const Index& norm_debug) {
  ARTS_USER_ERROR_IF (atmosphere_dim != 1,
                      "Only 1D is supported here for now");

  // Number of zenith angles.
  const Index Nza = za_grid.nelem();

  ARTS_USER_ERROR_IF (za_grid[0] != 0. || za_grid[Nza - 1] != 180.,
                      "The range of *za_grid* must [0 180].");

  // Number of azimuth angles.
  const Index Naa = aa_grid.nelem();

  ARTS_USER_ERROR_IF (Naa > 1 && (aa_grid[0] != 0. || aa_grid[Naa - 1] != 360.),
                      "The range of *aa_grid* must [0 360].");

  // Get stokes dimension from *doit_scat_field*:
  const Index stokes_dim = doit_scat_field.ncols();
  ARTS_ASSERT(stokes_dim > 0 || stokes_dim < 5);

  // To use special interpolation functions for atmospheric fields we
  // use ext_mat_field and abs_vec_field:
  Tensor5 ext_mat_field(cloudbox_limits[1] - cloudbox_limits[0] + 1,
                        1,
                        1,
                        stokes_dim,
                        stokes_dim,
                        0.);
  Tensor4 abs_vec_field(
      cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1, stokes_dim, 0.);

  const Index Np = doit_scat_field.nvitrines();

  Tensor5 doit_scat_ext_field(doit_scat_field.nvitrines(),
                              doit_scat_field.nshelves(),
                              doit_scat_field.nbooks(),
                              doit_scat_field.npages(),
                              doit_scat_field.nrows(),
                              0.);

  Index aa_index_local = 0;

  // Calculate scattering extinction field
  for (Index za_index_local = 0; za_index_local < Nza;
       za_index_local++) {
    // This function has to be called inside the angular loop, as
    // spt_calc_agenda takes *za_index* and *aa_index*
    // from the workspace.
    cloud_fieldsCalc(ws,
                     ext_mat_field,
                     abs_vec_field,
                     spt_calc_agenda,
                     za_index_local,
                     aa_index_local,
                     cloudbox_limits,
                     t_field,
                     pnd_field);

    for (Index p_index = 0;
         p_index <= (cloudbox_limits[1] - cloudbox_limits[0]);
         p_index++) {
      // For all in p_grid (in cloudbox):
      // I_ext = (ext_mat_field - abs_vec_field) * cloudbox_field_mono
      // equivalent to:
      // I_ext = I * (K11-a1) + Q * (K12 - a2) + U * (K13 - a3) + V * (K14 - a4)
      for (Index i = 0; i < stokes_dim; i++) {
        doit_scat_ext_field(p_index, 0, 0, za_index_local, 0) +=
            cloudbox_field_mono(p_index, 0, 0, za_index_local, 0, i) *
            (ext_mat_field(p_index, 0, 0, 0, i) -
             abs_vec_field(p_index, 0, 0, i));
      }
    }
  }

  Numeric corr_max = .0;
  Index corr_max_p_index = -1;

  for (Index p_index = 0; p_index < Np; p_index++) {
    // Calculate scattering integrals
    const Numeric scat_int = AngIntegrate_trapezoid(
        doit_scat_field(p_index, 0, 0, joker, 0, 0), za_grid);

    const Numeric scat_ext_int = AngIntegrate_trapezoid(
        doit_scat_ext_field(p_index, 0, 0, joker, 0), za_grid);

    // Calculate factor between scattered extinction field integral
    // and scattered field integral
    const Numeric corr_factor = scat_ext_int / scat_int;

    // If no scattering is present, the correction factor can become
    // inf or nan. We just don't apply it for those cases.
    if (!std::isnan(corr_factor) && !std::isinf(corr_factor)) {
      if (abs(corr_factor) > abs(corr_max)) {
        corr_max = corr_factor;
        corr_max_p_index = p_index;
      }
      if (norm_debug) {
      }
      ARTS_USER_ERROR_IF (abs(1. - corr_factor) > norm_error_threshold,
          "ERROR: DOIT correction factor exceeds threshold (=",
          norm_error_threshold, "): ", setprecision(4),
          1. - corr_factor, " at p_index ", p_index, "\n")
      if (abs(1. - corr_factor) > norm_error_threshold / 2.) {
      }

      // Scale scattered field with correction factor
      doit_scat_field(p_index, 0, 0, joker, 0, joker) *= corr_factor;
    } else if (norm_debug) {
    }
  }

  ostringstream os;
  if (corr_max_p_index != -1) {
    os << "  Max. DOIT correction factor in this iteration: " << 1. - corr_max
       << " at p_index " << corr_max_p_index << "\n";
  } else {
    os << "  No DOIT correction performed in this iteration.\n";
  }
}
