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
         const String&   iy_unit, 
       ConstVectorView   f_grid,
   const Numeric&        n,
   const ArrayOfIndex&   i_pol );

void apply_iy_unit2( 
   Tensor3View           J,
   ConstMatrixView       iy, 
   const String&         iy_unit, 
   ConstVectorView       f_grid,
   const Numeric&        n,
   const ArrayOfIndex&   i_pol );

void ze_cfac(
         Vector&    fac,
   const Vector&    f_grid,
   const Numeric&   ze_tref,
   const Numeric&   k2 );

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

void ext2trans(
         MatrixView   trans_mat,
         Index&       icase,
   ConstMatrixView    ext_mat_av,
   const Numeric&     l_step );

void get_iy(
         Workspace&   ws,
         Matrix&      iy,
   ConstTensor3View   t_field,
   ConstTensor3View   z_field,
   ConstTensor4View   vmr_field,
   ConstTensor4View   nlte_field,
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
  const Index&            iy_id,
  const Index&            jacobian_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
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
  const Tensor3&          surface_props_data,
  const Agenda&           iy_main_agenda,
  const Agenda&           iy_space_agenda,
  const Agenda&           iy_surface_agenda,
  const Agenda&           iy_cloudbox_agenda,
  const Index&            iy_agenda_call1,
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

void get_ppath_cloudvars( 
        ArrayOfIndex&                  clear2cloudy,
        Matrix&                        ppath_pnd,
        ArrayOfMatrix&                 ppath_dpnd_dx,
  const Ppath&                         ppath,
  const Index&                         atmosphere_dim,
  const ArrayOfIndex&                  cloudbox_limits,
  const Tensor4&                       pnd_field,
  const ArrayOfTensor4&                dpnd_field_dx );

void get_ppath_f( 
        Matrix&    ppath_f,
  const Ppath&     ppath,
  ConstVectorView  f_grid, 
  const Index&     atmosphere_dim,
  const Numeric&   rte_alonglos_v,
  ConstMatrixView  ppath_wind );

Range get_rowindex_for_mblock( 
  const Sparse&   sensor_response, 
  const Index&    imblock );

void iy_transmission_mult( 
       Tensor3&      iy_trans_total,
  ConstTensor3View   iy_trans_old,
  ConstTensor3View   iy_trans_new );

void iy_transmission_mult( 
       Matrix&       iy_new,
  ConstTensor3View   iy_trans,
  ConstMatrixView    iy_old );

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
  ConstTensor4View                  nlte_field,
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

void ext_mat_case(Index& icase, ConstMatrixView ext_mat, const Index stokes_dim);

void get_stepwise_frequency_grid(VectorView ppath_f_grid,
                                 ConstVectorView  f_grid,
                                 ConstVectorView ppath_wind,
                                 ConstVectorView ppath_line_of_sight,
                                 const Numeric& rte_alonglos_v,
                                 const Index& atmosphere_dim);

void get_stepwise_f_partials(Vector& f_partials,
                             const Index& component,
                             ConstVectorView& line_of_sight,
                             ConstVectorView f_grid, 
                             const Index& atmosphere_dim);

void get_stepwise_blackbody_radiation(VectorView B,
                                      VectorView dB_dT,
                                      ConstVectorView ppath_f_grid,
                                      const Numeric& ppath_temperature,
                                      const bool& do_temperature_derivative);

void get_stepwise_clearsky_propmat(Workspace& ws,
                                   PropagationMatrix& K,
                                   StokesVector& S,
                                   Index& lte,
                                   ArrayOfPropagationMatrix& dK_dx,
                                   ArrayOfStokesVector& dS_dx,
                                   const Agenda& propmat_clearsky_agenda,
                                   const ArrayOfRetrievalQuantity& jacobian_quantities,
                                   const PropmatPartialsData& partial_derivatives,
                                   ConstVectorView ppath_f_grid,
                                   ConstVectorView ppath_magnetic_field,
                                   ConstVectorView ppath_line_of_sight,
                                   ConstVectorView ppath_nlte_temperatures,
                                   ConstVectorView ppath_vmrs,
                                   const Numeric& ppath_temperature,
                                   const Numeric& ppath_pressure,
                                   const ArrayOfIndex& jacobian_species,
                                   const bool& jacobian_do);

void adapt_stepwise_partial_derivatives(ArrayOfPropagationMatrix& dK_dx,
                                        ArrayOfStokesVector& dS_dx,
                                        const ArrayOfRetrievalQuantity& jacobian_quantities,
                                        ConstVectorView ppath_f_grid,
                                        ConstVectorView ppath_line_of_sight,
                                        ConstVectorView ppath_vmrs,
                                        const Numeric& ppath_temperature,
                                        const Numeric& ppath_pressure,
                                        const ArrayOfIndex& jacobian_species,
                                        const ArrayOfIndex& jacobian_wind,
                                        const Index& lte,
                                        const Index& atmosphere_dim,
                                        const bool& jacobian_do);

void get_stepwise_transmission_matrix(Tensor3View cumulative_transmission,
                                      Tensor3View T,
                                      Tensor4View dT_dx_close,
                                      Tensor4View dT_dx_far,
                                      ConstTensor3View cumulative_transmission_close,
                                      const PropagationMatrix& K_close,
                                      const PropagationMatrix& K_far,
                                      const ArrayOfPropagationMatrix& dK_close_dx,
                                      const ArrayOfPropagationMatrix& dK_far_dx,
                                      const Numeric& ppath_distance,
                                      const bool& first_level);

void sum_stepwise_scalar_tau_and_extmat_case(VectorView scalar_tau,
                                             ArrayOfIndex& extmat_case,
                                             const PropagationMatrix& upper_level,
                                             const PropagationMatrix& lower_level,
                                             const Numeric& distance);

void get_stepwise_scattersky_propmat(StokesVector& ap,
                                     PropagationMatrix& Kp,
                                     ArrayOfStokesVector& dap_dx,
                                     ArrayOfPropagationMatrix& dKp_dx,
                                     const ArrayOfRetrievalQuantity& jacobian_quantities,
                                     ConstMatrixView ppath_1p_pnd,       // the ppath_pnd at this ppath point
                                     const ArrayOfMatrix& ppath_dpnd_dx, // the full ppath_dpnd_dx, ie all ppath points
                                     const Index ppath_1p_id,
                                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                                     ConstVectorView ppath_line_of_sight,
                                     ConstVectorView ppath_temperature,
                                     const Index& atmosphere_dim,
                                     const bool& jacobian_do);

void get_stepwise_scattersky_source(StokesVector& Sp,
                                    ArrayOfStokesVector& dSp_dx,
                                    const ArrayOfRetrievalQuantity& jacobian_quantities,
                                    ConstVectorView ppath_1p_pnd,
                                    const ArrayOfMatrix& ppath_dpnd_dx,
                                    const Index ppath_1p_id,
                                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                                    ConstTensor7View doit_i_field,
                                    ConstVectorView scat_za_grid,
                                    ConstVectorView scat_aa_grid,
                                    ConstMatrixView ppath_line_of_sight,
                                    const GridPos& ppath_pressure,
                                    const Vector& temperature,
                                    const Index& atmosphere_dim,
                                    const bool& jacobian_do,
                                    const Index& t_interp_order=1);

void get_stepwise_effective_source(MatrixView J,
                                   Tensor3View dJ_dx,
                                   const PropagationMatrix& K,
                                   const StokesVector& a,
                                   const StokesVector& S,
                                   const ArrayOfPropagationMatrix& dK_dx,
                                   const ArrayOfStokesVector& da_dx,
                                   const ArrayOfStokesVector& dS_dx,
                                   ConstVectorView B,
                                   ConstVectorView dB_dT,
                                   const ArrayOfRetrievalQuantity& jacobian_quantities,
                                   const bool& jacobian_do);

void rtmethods_jacobian_init(
         ArrayOfIndex&               jac_species_i,
         ArrayOfIndex&               jac_scat_i,
         ArrayOfIndex&               jac_is_t,
         ArrayOfIndex&               jac_wind_i,
         ArrayOfIndex&               jac_mag_i,
         ArrayOfIndex&               jac_other,
         ArrayOfTensor3&             diy_dx,
         ArrayOfTensor3&             diy_dpath,         
   const Index&                      ns,
   const Index&                      nf,
   const Index&                      np,
   const Index&                      nq,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const ArrayOfString&              scat_species,         
   const ArrayOfTensor4&             dpnd_field_dx,
   const PropmatPartialsData&        ppd,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const Index&                      iy_agenda_call1,
   const bool                        is_active = false );

void rtmethods_jacobian_finalisation(
         Workspace&                  ws,
         ArrayOfTensor3&             diy_dx,
         ArrayOfTensor3&             diy_dpath,  
   const Index&                      ns,
   const Index&                      nf,
   const Index&                      np,
   const Index&                      atmosphere_dim,
   const Ppath&                      ppath,
   const Vector&                     ppvar_p,
   const Vector&                     ppvar_t,
   const Matrix&                     ppvar_vmr,
   const Index&                      iy_agenda_call1,         
   const Tensor3&                    iy_transmission,
   const Agenda&                     water_psat_agenda,   
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfIndex                jac_species_i,
   const ArrayOfIndex                jac_is_t);

void rtmethods_unit_conversion(
         Matrix&                     iy,
         ArrayOfTensor3&             diy_dx,
         Tensor3&                    ppvar_iy,  
   const Index&                      ns,
   const Index&                      np,
   const Vector&                     f_grid,         
   const Ppath&                      ppath,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const Index&                      j_analytical_do,
   const String&                     iy_unit );

#endif  // rte_h
