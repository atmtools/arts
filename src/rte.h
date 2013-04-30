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
          Matrix&      iy,
    const Index&       stokes_dim,
    ConstVectorView    bbar,
    ConstTensor3View   t );

void ext2trans(
         MatrixView   trans_mat,
   ConstMatrixView    ext_mat_av,
   const Numeric&     l_step );

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
  const Agenda&           iy_main_agenda,
  const Agenda&           iy_space_agenda,
  const Agenda&           iy_surface_agenda,
  const Agenda&           iy_cloudbox_agenda,
  const Verbosity&        verbosity);

void get_ppath_atmvars( 
        Vector&      ppath_p, 
        Vector&      ppath_t, 
        Matrix&      ppath_vmr, 
        Matrix&      ppath_wind, 
        Matrix&      ppath_mag,
  const Ppath&       ppath,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstTensor3View   t_field,
  ConstTensor4View   vmr_field,
  ConstTensor3View   wind_u_field,
  ConstTensor3View   wind_v_field,
  ConstTensor3View   wind_w_field,
  ConstTensor3View   mag_u_field,
  ConstTensor3View   mag_v_field,
  ConstTensor3View   mag_w_field );

void get_ppath_abs( 
        Workspace&      ws,
        Tensor5&        ppath_abs,
  const Agenda&         propmat_clearsky_agenda,
  const Ppath&          ppath,
  ConstVectorView       ppath_p, 
  ConstVectorView       ppath_t, 
  ConstMatrixView       ppath_vmr, 
  ConstMatrixView       ppath_f, 
  ConstMatrixView       ppath_mag,
  ConstVectorView       f_grid, 
  const Index&          stokes_dim,
  const bool&           only_sum_abs );

void get_ppath_blackrad( 
        Workspace&   ws,
        Matrix&      ppath_blackrad,
  const Agenda&      blackbody_radiation_agenda,
  const Ppath&       ppath,
  ConstVectorView    ppath_t, 
  ConstMatrixView    ppath_f );

void get_ppath_ext( 
        ArrayOfIndex&                  clear2cloudbox,
        Tensor3&                       pnd_abs_vec, 
        Tensor4&                       pnd_ext_mat, 
  Array<ArrayOfSingleScatteringData>&  scat_data,
        Matrix&                        ppath_pnd,
  const Ppath&                         ppath,
  ConstVectorView                      ppath_t, 
  const Index&                         stokes_dim,
  ConstMatrixView                      ppath_f, 
  const Index&                         atmosphere_dim,
  const ArrayOfIndex&                  cloudbox_limits,
  const Tensor4&                       pnd_field,
  const Index&                         use_mean_scat_data,
  const ArrayOfSingleScatteringData&   scat_data_raw,
  const Verbosity&                     verbosity );

void get_ppath_f( 
        Matrix&    ppath_f,
  const Ppath&     ppath,
  ConstVectorView  f_grid, 
  const Index&     atmosphere_dim,
  const Numeric&   rte_alonglos_v,
  ConstMatrixView  ppath_wind );

void get_ppath_trans( 
        Tensor4&        trans_partial,
        Tensor4&        trans_cumulat,
        Vector&         scalar_tau,
  const Ppath&          ppath,
  ConstTensor5View&     ppath_abs,
  ConstVectorView       f_grid, 
  const Index&          stokes_dim );

void get_ppath_trans2( 
        Tensor4&        trans_partial,
        Tensor4&        trans_cumulat,
        Vector&         scalar_tau,
  const Ppath&          ppath,
  ConstTensor5View&     ppath_abs,
  ConstVectorView       f_grid, 
  const Index&          stokes_dim,
  const ArrayOfIndex&   clear2cloudbox,
  ConstTensor4View      pnd_ext_mat );

Range get_rowindex_for_mblock( 
  const Sparse&   sensor_response, 
  const Index&    imblock );

void iy_transmission_mult( 
       Tensor3&      iy_trans_total,
  ConstTensor3View   iy_trans_old,
  ConstTensor3View   iy_trans_new );

void iyb_calc(
        Workspace&                  ws,
        Vector&                     iyb,
        ArrayOfVector&              iyb_aux,
        ArrayOfMatrix&              diyb_dx,
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
  ConstVectorView                   mblock_za_grid,
  ConstVectorView                   mblock_aa_grid,
  const Index&                      antenna_dim,
  const Agenda&                     iy_main_agenda,
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

void surface_calc(
              Matrix&         iy,
        ConstTensor3View      I,
        ConstMatrixView       surface_los,
        ConstTensor4View      surface_rmatrix,
        ConstMatrixView       surface_emission );

#endif  // rte_h
