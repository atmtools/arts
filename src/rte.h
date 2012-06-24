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
#include "complex.h"          
#include "jacobian.h"
#include "ppath.h"
#include "matpackI.h"
#include "matpackII.h"
#include "matpackIII.h"
#include "optproperties.h"



// A macro to loop analytical jacobian quantities
#define FOR_ANALYTICAL_JACOBIANS_DO(what_to_do) \
  for( Index iq=0; iq<jacobian_quantities.nelem(); iq++ ) \
    { \
      if( jacobian_quantities[iq].Analytical() ) \
        { what_to_do } \
    } 



/*===========================================================================
  === Functions in rte.cc
  ===========================================================================*/

void adjust_los( 
         VectorView   los, 
   const Index &      atmosphere_dim );

void apply_y_unit( 
            MatrixView   iy, 
         const String&   y_unit, 
       ConstVectorView   f_grid,
   const ArrayOfIndex&   i_pol );

void apply_y_unit2( 
   Tensor3View           J,
   ConstMatrixView       iy, 
   const String&         y_unit, 
   ConstVectorView       f_grid,
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
  ConstTensor3View   edensity_field,
  const Index&       f_index,
  ConstVectorView    refellipsoid,
  ConstMatrixView    z_surface,
  const Ppath&       ppath,
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
  ConstTensor3View   edensity_field,
  const Index&       f_index,
  ConstVectorView    refellipsoid,
  ConstMatrixView    z_surface,
  const Ppath&       ppath,
  const Verbosity&   verbosity );

void ext2trans(
         MatrixView   trans_mat,
   ConstMatrixView    ext_mat_av,
   const Numeric&     l_step );

void get_iy_of_background(
        Workspace&        ws,
        Matrix&           iy,
        Matrix&           iy_error,
        Index&            iy_error_type,
        Matrix&           iy_aux,
        ArrayOfTensor3&   diy_dx,
  ConstTensor3View        iy_transmission,
  const Index&            jacobian_do,
  const Ppath&            ppath,
  const Index&            atmosphere_dim,
  ConstTensor3View        t_field,
  ConstTensor3View        z_field,
  ConstTensor4View        vmr_field,
  const Index&            cloudbox_on,
  const Index&            stokes_dim,
  ConstVectorView         f_grid,
  const Agenda&           iy_clearsky_agenda,
  const Agenda&           iy_space_agenda,
  const Agenda&           iy_surface_agenda,
  const Agenda&           iy_cloudbox_agenda,
  const Verbosity&        verbosity);

void get_ppath_atmvars( 
        Vector&      ppath_p, 
        Vector&      ppath_t, 
        Matrix&      ppath_vmr, 
        Vector&      ppath_wind_u, 
        Vector&      ppath_wind_v, 
        Vector&      ppath_wind_w, 
  const Ppath&       ppath,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstTensor3View   t_field,
  ConstTensor4View   vmr_field,
  ConstTensor3View   wind_u_field,
  ConstTensor3View   wind_v_field,
  ConstTensor3View   wind_w_field );

void get_ppath_cloudrtvars(Workspace&       ws,
                           Tensor3&         ppath_asp_abs_vec, 
                           Tensor4&         ppath_asp_ext_mat, 
                           Tensor3&         ppath_pnd_abs_vec, 
                           Tensor4&         ppath_pnd_ext_mat, 
                           Tensor4&         ppath_transmission,
                           Tensor3&         total_transmission,
                           Matrix&          ppath_emission, 
                           Array<ArrayOfSingleScatteringData>&  scat_data,
                           const Agenda&    abs_scalar_gas_agenda,
                           const Agenda&    emission_agenda,
                           const Agenda&    opt_prop_gas_agenda,
                           const Ppath&     ppath,
                           ConstVectorView  ppath_p, 
                           ConstVectorView  ppath_t, 
                           ConstMatrixView  ppath_vmr, 
                           ConstVectorView  ppath_wind_u, 
                           ConstVectorView  ppath_wind_v, 
                           ConstVectorView  ppath_wind_w, 
                           ConstMatrixView  ppath_pnd, 
                           const Index&     use_mean_scat_data,
                           const ArrayOfSingleScatteringData&   scat_data_raw,
                           const Index&     stokes_dim,
                           ConstVectorView  f_grid, 
                           const Index&     atmosphere_dim,
                           const Index&     emission_do,
                           const Verbosity& verbosity);

void get_ppath_pnd( 
        Matrix&         ppath_pnd, 
  const Ppath&          ppath,
  const Index&          atmosphere_dim,
  const ArrayOfIndex&   cloudbox_limits,
  ConstTensor4View      pnd_field );

void get_ppath_rtvars( 
        Workspace&   ws,
        Tensor3&     ppath_abs_scalar, 
        Matrix&      ppath_tau,
        Vector&      total_tau,
        Matrix&      ppath_emission, 
  const Agenda&      abs_scalar_agenda,
  const Agenda&      emission_agenda,
  const Ppath&       ppath,
  ConstVectorView    ppath_p, 
  ConstVectorView    ppath_t, 
  ConstMatrixView    ppath_vmr, 
  ConstVectorView    ppath_wind_u, 
  ConstVectorView    ppath_wind_v, 
  ConstVectorView    ppath_wind_w, 
  ConstVectorView    f_grid, 
  const Index&       f_index,
  const Index&       atmosphere_dim,
  const Index&       emission_do );

Range get_rowindex_for_mblock( 
  const Sparse&   sensor_response, 
  const Index&    imblock );

void iy_transmission_for_scalar_tau( 
       Tensor3&     iy_transmission,
  const Index&      stokes_dim,
  ConstVectorView   tau );

void iy_transmission_mult( 
       Tensor3&      iy_trans_new,
  ConstTensor3View   iy_transmission,
  ConstTensor3View   trans );

void iy_transmission_mult_scalar_tau( 
       Tensor3&      iy_trans_new,
  ConstTensor3View   iy_transmission,
  ConstVectorView    tau );

void iyb_calc(
        Workspace&                  ws,
        Vector&                     iyb,
        Vector&                     iyb_error,
        Index&                      iy_error_type,
        Vector&                     iyb_aux,
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
  ConstVectorView                   mblock_za_grid,
  ConstVectorView                   mblock_aa_grid,
  const Index&                      antenna_dim,
  const Agenda&                     iy_clearsky_agenda,
  const String&                     y_unit,
  const Index&                      j_analytical_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices,
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

void rte_step_std(
         //Output and Input:
         VectorView stokes_vec,
         MatrixView trans_mat,
         //Input
         ConstMatrixView ext_mat_av,
         ConstVectorView abs_vec_av,
         ConstVectorView sca_vec_av, 
         const Numeric& l_step,
         const Numeric& rte_planck_value );

void surface_calc(
              Matrix&         iy,
        ConstTensor3View      I,
        ConstMatrixView       surface_los,
        ConstTensor4View      surface_rmatrix,
        ConstMatrixView       surface_emission );

void trans_step_std(//Output and Input:
              VectorView stokes_vec,
              MatrixView trans_mat,
              //Input
              ConstMatrixView ext_mat_av,
              const Numeric& l_step );

#endif  // rte_h
