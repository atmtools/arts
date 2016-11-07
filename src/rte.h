/* Copyright (C) 2002-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
                            
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



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   rte.h
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2002-05-29

  \brief  Declaration of functions in rte.cc.
*/



#ifndef rte_h
#define rte_h

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "agenda_class.h"
#include "arts.h"
#include "auto_md.h"
#include "complex.h"          
#include "jacobian.h"
#include "ppath.h"
#include "matpackI.h"
#include "matpackII.h"
#include "matpackIII.h"
#include "optproperties.h"



/*===========================================================================
  === Functions in rte.cc
  ===========================================================================*/

void adjust_los( 
         VectorView   los, 
   const Index &      atmosphere_dim );

void apply_iy_unit( 
            MatrixView   iy, 
         const String&   y_unit, 
       ConstVectorView   f_grid,
   const Numeric&        n,
   const ArrayOfIndex&   i_pol );

void apply_iy_unit2( 
   Tensor3View           J,
   ConstMatrixView       iy, 
   const String&         y_unit, 
   ConstVectorView       f_grid,
   const Numeric&        n,
   const ArrayOfIndex&   i_pol );

void bending_angle1d( 
        Numeric&   alpha,
  const Ppath&     ppath );

void defocusing_general( 
        Workspace&   ws,
        Numeric&     dlf,
  const Agenda&      ppath_step_agenda,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstVectorView    lat_grid,
  ConstVectorView    lon_grid,
  ConstTensor3View   t_field,
  ConstTensor3View   z_field,
  ConstTensor4View   vmr_field,
  ConstVectorView    f_grid,
  ConstVectorView    refellipsoid,
  ConstMatrixView    z_surface,
  const Ppath&       ppath,
  const Numeric&     ppath_lmax,
  const Numeric&     ppath_lraytrace,
  const Numeric&     dza,
  const Verbosity&   verbosity );

void defocusing_sat2sat( 
        Workspace&   ws,
        Numeric&     dlf,
  const Agenda&      ppath_step_agenda,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstVectorView    lat_grid,
  ConstVectorView    lon_grid,
  ConstTensor3View   t_field,
  ConstTensor3View   z_field,
  ConstTensor4View   vmr_field,
  ConstVectorView    f_grid,
  ConstVectorView    refellipsoid,
  ConstMatrixView    z_surface,
  const Ppath&       ppath,
  const Numeric&     ppath_lmax,
  const Numeric&     ppath_lraytrace,
  const Numeric&     dza,
  const Verbosity&   verbosity );

Numeric dotprod_with_los(
  ConstVectorView   los, 
  const Numeric&    u,
  const Numeric&    v,
  const Numeric&    w,
  const Index&      atmosphere_dim );

void emission_rtstep(
          Matrix&         iy,
    const Index&          stokes_dim,
    ConstVectorView       bbar,
    const ArrayOfIndex&   extmat_case,
    ConstTensor3View      t,
    const bool&           nonlte,
    ConstTensor3View      extbar,
    ConstMatrixView       sourcebar );

void ext2trans(
         MatrixView   trans_mat,
         Index&       icase,
   ConstMatrixView    ext_mat_av,
   const Numeric&     l_step );

void ext2trans_and_ext2dtrans_dx(
    MatrixView   trans_mat,
    Tensor3View  dtrans_mat_dx_upp,
    Tensor3View  dtrans_mat_dx_low,
    Index&       icase,
    ConstMatrixView    ext_mat,
    ConstTensor3View   dext_mat_dx_upp,
    ConstTensor3View   dext_mat_dx_low,
    const Numeric&     lstep );

void get_iy(
         Workspace&   ws,
         Matrix&      iy,
   ConstTensor3View   t_field,
   ConstTensor3View   z_field,
   ConstTensor4View   vmr_field,
   const Index&       cloudbox_on,
   ConstVectorView    f_grid,
   ConstVectorView    rte_pos,
   ConstVectorView    rte_los,
   ConstVectorView    rte_pos2,
   const String&      iy_unit,
   const Agenda&      iy_main_agenda );

void get_iy_of_background(
        Workspace&        ws,
        Matrix&           iy,
        ArrayOfTensor3&   diy_dx,
  ConstTensor3View        iy_transmission,
  const Index&            jacobian_do,
  const Ppath&            ppath,
  ConstVectorView         rte_pos2,
  const Index&            atmosphere_dim,
  ConstTensor3View        t_field,
  ConstTensor3View        z_field,
  ConstTensor4View        vmr_field,
  const Index&            cloudbox_on,
  const Index&            stokes_dim,
  ConstVectorView         f_grid,
  const String&           iy_unit,  
  const Agenda&           iy_main_agenda,
  const Agenda&           iy_space_agenda,
  const Agenda&           iy_surface_agenda,
  const Agenda&           iy_cloudbox_agenda,
  const Verbosity&        verbosity);

void get_ppath_atmvars( 
        Vector&      ppath_p, 
        Vector&      ppath_t, 
        Matrix&      ppath_t_nlte,
        Matrix&      ppath_vmr, 
        Matrix&      ppath_wind, 
        Matrix&      ppath_mag,
  const Ppath&       ppath,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstTensor3View   t_field,
  ConstTensor4View   t_nlte_field,
  ConstTensor4View   vmr_field,
  ConstTensor3View   wind_u_field,
  ConstTensor3View   wind_v_field,
  ConstTensor3View   wind_w_field,
  ConstTensor3View   mag_u_field,
  ConstTensor3View   mag_v_field,
  ConstTensor3View   mag_w_field );

void get_ppath_pmat( 
        Workspace&      ws,
        Tensor4&        ppath_ext,
        Tensor3&        ppath_nlte_source,
        ArrayOfIndex&   lte,
        Tensor5&        abs_per_species,
        Tensor5&        dppath_ext_dx,
        Tensor4&        dppath_nlte_source_dx,
  const Agenda&         propmat_clearsky_agenda,
  const ArrayOfRetrievalQuantity& jacobian_quantities,
  const Ppath&          ppath,
  ConstVectorView       ppath_p, 
  ConstVectorView       ppath_t, 
  ConstMatrixView       ppath_t_nlte, 
  ConstMatrixView       ppath_vmr, 
  ConstMatrixView       ppath_f, 
  ConstMatrixView       ppath_mag,
  ConstVectorView       f_grid, 
  const Index&          stokes_dim,
  const ArrayOfIndex&   ispecies );

void get_ppath_blackrad( 
        Matrix&      ppath_blackrad,
  const Ppath&       ppath,
  ConstVectorView    ppath_t, 
  ConstMatrixView    ppath_f );

void get_dppath_blackrad_dt( 
        Matrix&             dppath_blackrad_dt,
        ConstVectorView     ppath_t, 
        ConstMatrixView     ppath_f,
        const ArrayOfIndex& jac_is_t,
        const bool&         j_analytical_do);

void get_ppath_ext( 
        ArrayOfIndex&                  clear2cloudbox,
        Tensor3&                       pnd_abs_vec, 
        Tensor4&                       pnd_ext_mat, 
  Array<ArrayOfArrayOfSingleScatteringData>&  scat_data_single,
        Matrix&                        ppath_pnd,
  const Ppath&                         ppath,
  ConstVectorView                      ppath_t, 
  const Index&                         stokes_dim,
  ConstMatrixView                      ppath_f, 
  const Index&                         atmosphere_dim,
  const ArrayOfIndex&                  cloudbox_limits,
  const Tensor4&                       pnd_field,
  const Index&                         use_mean_scat_data,
  const ArrayOfArrayOfSingleScatteringData&   scat_data,
  const Verbosity&                     verbosity );

void get_ppath_f( 
        Matrix&    ppath_f,
  const Ppath&     ppath,
  ConstVectorView  f_grid, 
  const Index&     atmosphere_dim,
  const Numeric&   rte_alonglos_v,
  ConstMatrixView  ppath_wind );

void get_ppath_f_partials( 
Matrix&    ppath_f_partials,
const Index& component,
const Ppath&     ppath,
ConstVectorView  f_grid, 
const Index&     atmosphere_dim);

void get_ppath_trans( 
        Tensor4&               trans_partial,
        ArrayOfArrayOfIndex&   extmat_case,
        Tensor4&               trans_cumulat,
        Vector&                scalar_tau,
  const Ppath&                 ppath,
  ConstTensor4View&            ppath_ext,
  ConstVectorView              f_grid, 
  const Index&                 stokes_dim );

void get_ppath_trans_and_dppath_trans_dx( 
        Tensor4&               trans_partial,
        Tensor5&               dtrans_partial_dx_from_above,
        Tensor5&               dtrans_partial_dx_from_below,
        ArrayOfArrayOfIndex&   extmat_case,
        Tensor4&               trans_cumulat,
        Vector&                scalar_tau,
  const Ppath&                 ppath,
  ConstTensor4View&            ppath_ext,
  ConstTensor5View&            dppath_ext_dx,
  const ArrayOfRetrievalQuantity& jacobian_quantities,
  ConstVectorView              f_grid, 
  const ArrayOfIndex&          for_flux,
  const Index&                 stokes_dim );

void get_ppath_trans2( 
        Tensor4&               trans_partial,
        ArrayOfArrayOfIndex&   extmat_case,
        Tensor4&               trans_cumulat,
        Vector&                scalar_tau,
  const Ppath&                 ppath,
  ConstTensor4View&            ppath_ext,
  ConstVectorView              f_grid, 
  const Index&                 stokes_dim,
  const ArrayOfIndex&          clear2cloudbox,
  ConstTensor4View             pnd_ext_mat );

void get_ppath_trans2_and_dppath_trans_dx(  Tensor4&               trans_partial,
                                            Tensor5&               dtrans_partial_dx_from_above,
                                            Tensor5&               dtrans_partial_dx_from_below,
                                            ArrayOfArrayOfIndex&   extmat_case,
                                            Tensor4&               trans_cumulat,
                                            Vector&                scalar_tau,
                                            const Ppath&                 ppath,
                                            ConstTensor4View&            ppath_ext,
                                            ConstTensor5View&            dppath_ext_dx,
                                            const ArrayOfRetrievalQuantity& jacobian_quantities,
                                            ConstVectorView              f_grid, 
                                            const Index&                 stokes_dim,
                                            const ArrayOfIndex&          clear2cloudbox,
                                            ConstTensor4View             pnd_ext_mat );

Range get_rowindex_for_mblock( 
  const Sparse&   sensor_response, 
  const Index&    imblock );

void iy_transmission_mult( 
       Tensor3&      iy_trans_total,
  ConstTensor3View   iy_trans_old,
  ConstTensor3View   iy_trans_new );

void get_ppath_pmat_and_tmat( 
                            Workspace&      ws,
                            Tensor4&        ppath_ext,
                            Tensor3&        ppath_nlte_source,
                            ArrayOfIndex&   lte,
                            Tensor5&        abs_per_species,
                            Tensor5&        dppath_ext_dx,
                            Tensor4&        dppath_nlte_source_dx,
                            Tensor4&               trans_partial,
                            Tensor5&               dtrans_partial_dx_above,
                            Tensor5&               dtrans_partial_dx_below,
                            ArrayOfArrayOfIndex&   extmat_case,
                            ArrayOfIndex&   clear2cloudbox,
                            Tensor4&               trans_cumulat,
                            Vector&                scalar_tau,
                            Tensor4&               pnd_ext_mat,
                            Matrix&                ppath_pnd,
                            const Agenda&         propmat_clearsky_agenda,
                            const ArrayOfRetrievalQuantity& jacobian_quantities,
                            const PropmatPartialsData&      ppd,
                            const Ppath&          ppath,
                            ConstVectorView       ppath_p, 
                            ConstVectorView       ppath_t, 
                            ConstMatrixView       ppath_t_nlte, 
                            ConstMatrixView       ppath_vmr, 
                            ConstMatrixView       ppath_mag,
                            ConstMatrixView       ppath_wind,
                            ConstMatrixView       ppath_f, 
                            ConstVectorView       f_grid, 
                            const ArrayOfIndex&   jac_species_i,
                            const ArrayOfIndex&   jac_is_t,
                            const ArrayOfIndex&   jac_wind_i,
                            const ArrayOfIndex&   jac_mag_i,
                            const ArrayOfIndex&   for_flux,
                            const ArrayOfIndex&   jac_other,
                            const ArrayOfIndex&   ispecies,
                            const ArrayOfArrayOfSingleScatteringData scat_data,
                            const Tensor4&        pnd_field,
                            const ArrayOfIndex&   cloudbox_limits,
                            const Index&          use_mean_scat_data,
                            const Numeric&        rte_alonglos_v,
                            const Index&          atmosphere_dim,
                            const Index&          stokes_dim,
                            const bool&           jacobian_do,
                            const bool&           cloudbox_on,
                            const Verbosity&      verbosity);

void iyb_calc(
        Workspace&                  ws,
        Vector&                     iyb,
        ArrayOfVector&              iyb_aux,
        ArrayOfMatrix&              diyb_dx,
        Matrix&                     geo_pos_matrix,
  const Index&                      imblock,
  const Index&                      atmosphere_dim,
  ConstTensor3View                  t_field,
  ConstTensor3View                  z_field,
  ConstTensor4View                  vmr_field,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  ConstVectorView                   f_grid,
  ConstMatrixView                   sensor_pos,
  ConstMatrixView                   sensor_los,
  ConstMatrixView                   transmitter_pos,
  ConstMatrixView                   mblock_dlos_grid,
  const String&                     iy_unit,  
  const Agenda&                     iy_main_agenda,
  const Agenda&                     geo_pos_agenda,
  const Index&                      j_analytical_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices,
  const ArrayOfString&              iy_aux_vars,
  const Verbosity&                  verbosity );

void mirror_los(
        Vector&     los_mirrored,
  ConstVectorView   los, 
  const Index&      atmosphere_dim );

void pos2true_latlon( 
          Numeric&     lat,
          Numeric&     lon,
    const Index&       atmosphere_dim,
    ConstVectorView    lat_grid,
    ConstVectorView    lat_true,
    ConstVectorView    lon_true,
    ConstVectorView    pos );

#endif  // rte_h
