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
  const Index&            iy_id,
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
        ArrayOfPropagationMatrix&        ppath_ext,
        ArrayOfStokesVector&        ppath_nlte_source,
        ArrayOfIndex&   lte,
        ArrayOfArrayOfPropagationMatrix&        abs_per_species,
        ArrayOfArrayOfPropagationMatrix&        dppath_ext_dx,
        ArrayOfArrayOfStokesVector&        dppath_nlte_source_dx,
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

void get_ppath_partopt( 
        Tensor3&                       pnd_abs_vec, 
        ArrayOfPropagationMatrix&      pnd_ext_mat, 
  Array<ArrayOfArrayOfSingleScatteringData>&  scat_data_single,
  const ArrayOfIndex&                  clear2cloudy,
  const Matrix&                        ppath_pnd,
  const Ppath&                         ppath,
  ConstVectorView                      ppath_t, 
  const Index&                         stokes_dim,
  ConstMatrixView                      ppath_f, 
  const Index&                         atmosphere_dim,
  const Index&                         use_mean_scat_data,
  const ArrayOfArrayOfSingleScatteringData&   scat_data,
  const Index&                         scat_data_checked,
  const Verbosity&                     verbosity );

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
  const ArrayOfPropagationMatrix& ppath_ext,
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
  const ArrayOfPropagationMatrix& ppath_ext,
  const ArrayOfArrayOfPropagationMatrix& dppath_ext_dx,
  const ArrayOfRetrievalQuantity& jacobian_quantities,
  ConstVectorView              f_grid, 
  const ArrayOfIndex&          for_distance_integration,
  const Index&                 stokes_dim );

void get_ppath_trans2( 
        Tensor4&               trans_partial,
        ArrayOfArrayOfIndex&   extmat_case,
        Tensor4&               trans_cumulat,
        Vector&                scalar_tau,
  const Ppath&                 ppath,
  const ArrayOfPropagationMatrix& ppath_ext,
  ConstVectorView              f_grid, 
  const Index&                 stokes_dim,
  const ArrayOfIndex&          clear2cloudy,
  const ArrayOfPropagationMatrix& pnd_ext_mat );

void get_ppath_trans2_and_dppath_trans_dx(  Tensor4&               trans_partial,
                                            Tensor5&               dtrans_partial_dx_from_above,
                                            Tensor5&               dtrans_partial_dx_from_below,
                                            ArrayOfArrayOfIndex&   extmat_case,
                                            Tensor4&               trans_cumulat,
                                            Vector&                scalar_tau,
                                            const Ppath&                 ppath,
                                            const ArrayOfPropagationMatrix& ppath_ext,
                                            const ArrayOfArrayOfPropagationMatrix& dppath_ext_dx,
                                            const ArrayOfRetrievalQuantity& jacobian_quantities,
                                            ConstVectorView              f_grid, 
                                            const Index&                 stokes_dim,
                                            const ArrayOfIndex&          clear2cloudy,
                                            const ArrayOfPropagationMatrix& pnd_ext_mat );

Range get_rowindex_for_mblock( 
  const Sparse&   sensor_response, 
  const Index&    imblock );

void iy_transmission_mult( 
       Tensor3&      iy_trans_total,
  ConstTensor3View   iy_trans_old,
  ConstTensor3View   iy_trans_new );

void get_ppath_pmat_and_tmat( 
                            Workspace&      ws,
                            ArrayOfPropagationMatrix&        ppath_ext,
                            ArrayOfStokesVector&        ppath_nlte_source,
                            ArrayOfIndex&   lte,
                            ArrayOfArrayOfPropagationMatrix&        abs_per_species,
                            ArrayOfArrayOfPropagationMatrix&        dppath_ext_dx,
                            ArrayOfArrayOfStokesVector&        dppath_nlte_source_dx,
                            Tensor4&               trans_partial,
                            Tensor5&               dtrans_partial_dx_above,
                            Tensor5&               dtrans_partial_dx_below,
                            ArrayOfArrayOfIndex&   extmat_case,
                            ArrayOfIndex&   clear2cloudy,
                            Tensor4&               trans_cumulat,
                            Vector&                scalar_tau,
                            ArrayOfPropagationMatrix&               pnd_ext_mat,
                            Tensor3&               pnd_abs_vec,
                            Matrix&                ppath_pnd,
                            ArrayOfMatrix&         ppath_dpnd_dx,
                            Array<ArrayOfArrayOfSingleScatteringData>& scat_data_single,
                            const Agenda&         propmat_clearsky_agenda,
                            const ArrayOfRetrievalQuantity& jacobian_quantities,
                            const PropmatPartialsData&      ppd,
                            const Ppath&          ppath,
                            ConstVectorView       ppath_p, 
                            ConstVectorView       ppath_t, 
                            ConstMatrixView       ppath_t_nlte, 
                            ConstMatrixView       ppath_vmr, 
                            ConstMatrixView       ppath_mag,
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
                            const Index&          scat_data_checked,
                            const Tensor4&        pnd_field,
                            const ArrayOfTensor4& dpnd_field_dx,
                            const ArrayOfIndex&   cloudbox_limits,
                            const Index&          use_mean_scat_data,
                            const Index&          atmosphere_dim,
                            const Index&          stokes_dim,
                            const bool&           jacobian_do,
                            const bool&           cloudbox_on,
                            const Verbosity&      verbosity);

void get_ppath_scat_source(
                           Tensor4&         ppath_scat_source,
                           const ArrayOfArrayOfSingleScatteringData scat_data,
                           const Index&     scat_data_checked,
                           ConstTensor7View doit_i_field,
                           ConstVectorView  scat_za_grid,
                           ConstVectorView  f_grid, 
                           const Index&     stokes_dim,
                           const Ppath&     ppath,
                           ConstVectorView  ppath_t, 
                           ConstMatrixView  ppath_pnd,
                           const Index&     j_analytical_do,
                           const Index&     Naa,
                           const Verbosity& verbosity );

void get_ppath_scat_source_fixT(
                           Tensor4&         ppath_scat_source,
                           const ArrayOfArrayOfSingleScatteringData scat_data,
                           const Index&     scat_data_checked,
                           ConstTensor7View doit_i_field,
                           ConstVectorView  scat_za_grid,
                           ConstVectorView  f_grid, 
                           const Index&     stokes_dim,
                           const Ppath&     ppath,
                           ConstMatrixView  ppath_pnd,
                           const Index&     j_analytical_do,
                           const Index&     Naa,
                           const Numeric&   rtp_temp,
                           const Verbosity& verbosity );

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

void ext_mat_case(Index& icase, ConstMatrixView ext_mat, const Index stokes_dim);

void emission_rtstep_replacement( MatrixView iy,
                                  const Index stokes_dim,
                                  ConstVectorView planck_emission,
                                  const ArrayOfIndex&   extmat_case,
                                  ConstTensor3View transmission,
                                  const bool nonlte,
                                  const PropagationMatrix& propagation_matrix,
                                  const StokesVector& source_vector);

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
                                     ConstVectorView ppath_1p_pnd,       // the ppath_pnd at this ppath point
                                     const ArrayOfMatrix& ppath_dpnd_dx, // the full ppath_dpnd_dx, ie all ppath points
                                     const Index ppath_1p_id,
                                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                                     ConstVectorView ppath_line_of_sight,
                                     const Numeric& ppath_temperature,
                                     const Index& atmosphere_dim,
                                     const bool& jacobian_do,
                                     const Verbosity& verbosity);

void get_stepwise_scattersky_source(StokesVector& Sp,
                                    ArrayOfStokesVector& dSp_dx,
                                    const ArrayOfRetrievalQuantity& jacobian_quantities,
                                    ConstVectorView ppath_1p_pnd,       // the ppath_pnd at this ppath point
                                    const ArrayOfMatrix& ppath_dpnd_dx, // the full ppath_dpnd_dx, ie all ppath points
                                    const Index ppath_1p_id,
                                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                                    const Index& scat_data_checked,
                                    ConstTensor7View doit_i_field,
                                    ConstVectorView scat_za_grid,
                                    ConstVectorView scat_aa_grid,
                                    ConstVectorView ppath_line_of_sight,
                                    const GridPos& ppath_pressure,
                                    const Numeric& ppath_temperature,
                                    const Index& atmosphere_dim,
                                    const bool& jacobian_do,
                                    const Verbosity& verbosity);

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
         ArrayOfIndex&               jac_to_integrate,
         ArrayOfTensor3&             diy_dx,
         ArrayOfTensor3&             diy_dpath,         
   const Index&                      ns,
   const Index&                      nf,
   const Index&                      np,
   const Index&                      nq,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const ArrayOfString&              scat_species,         
   const PropmatPartialsData&        ppd,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,   
   const ArrayOfArrayOfIndex&        jacobian_indices, 
   const Index&                      iy_agenda_call1 );

void rtmethods_jacobian_finalisation(
         ArrayOfTensor3&             diy_dx,
   const Index&                      ns,
   const Index&                      nf,
   const Index&                      atmosphere_dim,
   const Ppath&                      ppath,
   const Vector&                     ppath_p,
   const Index&                      iy_agenda_call1,         
   const Tensor3&                    iy_transmission,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfIndex&               jac_to_integrate, 
   const ArrayOfTensor3&             diy_dpath );

void rtmethods_unit_conversion(
         Matrix&                     iy,
         ArrayOfTensor3&             diy_dx,
   const Index&                      ns,
   const Index&                      np,
   const Vector&                     f_grid,         
   const Ppath&                      ppath,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const Index&                      j_analytical_do,
   const String&                     iy_unit );

#endif  // rte_h
