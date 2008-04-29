/* Copyright (C) 2003-2008 Cory Davis <cory@met.ed.ac.uk>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

#ifndef montecarlo_h
#define montecarlo_h

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <stdexcept>
#include <cmath>
#include "messages.h"
#include "arts.h"
#include "ppath.h"
#include "matpackI.h"
#include "special_interp.h"
#include "check_input.h"
#include "rte.h"
#include "lin_alg.h"
#include "logic.h"
#include "optproperties.h"
#include "physics_funcs.h"
#include "xml_io.h"
#include "rng.h"
#include "scatrte.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;

void clear_rt_vars_at_gp(
                         MatrixView&             ext_mat_mono,
                         VectorView&             abs_vec_mono,
                         Numeric&                temperature,
                         const Agenda&           opt_prop_gas_agenda,
                         const Agenda&           abs_scalar_gas_agenda,
                         const Index&            f_index,
                         const GridPos&          gp_p,
                         const GridPos&          gp_lat,
                         const GridPos&          gp_lon,
                         const ConstVectorView   p_grid,
                         const ConstVectorView   lat_grid,
                         const ConstVectorView   lon_grid,
                         const ConstTensor3View  t_field,
                         const ConstTensor4View  vmr_field);

void cloud_atm_vars_by_gp(
                          VectorView pressure,
                          VectorView temperature,
                          MatrixView vmr,
                          MatrixView pnd,
                          const ArrayOfGridPos& gp_p,
                          const ArrayOfGridPos& gp_lat,
                          const ArrayOfGridPos& gp_lon,
                          const ArrayOfIndex& cloudbox_limits,
                          const ConstVectorView p_grid_cloud,
                          const ConstVectorView lat_grid_cloud,
                          const ConstVectorView lon_grid_cloud,
                          const ConstTensor3View   t_field_cloud,
                          const ConstTensor4View   vmr_field_cloud,
                          const ConstTensor4View   pnd_field
);

void Cloudbox_ppathCalc(
                        //  Output:
                        Ppath&          ppath,
                        Ppath&          ppath_step,
                        //  Input:
                        const Agenda&         ppath_step_agenda,
                        const Index&          atmosphere_dim,
                        const Vector&         p_grid,
                        const Vector&         lat_grid,
                        const Vector&         lon_grid,
                        const Tensor3&        z_field,
                        const Matrix&         r_geoid,
                        const Matrix&         z_surface,
                        const ArrayOfIndex&   cloudbox_limits,
                        const Vector&         rte_pos,
                        const Vector&         rte_los,
                        const Index& z_field_is_1D);


void Cloudbox_ppath_rteCalc(
                             Ppath&                ppathcloud,
                             Ppath&                ppath,
                             Ppath&                ppath_step,
                             //Vector&               ppath_p,
                             //Vector&               ppath_t,
                             //Matrix&               ppath_vmr,
                             Vector&               rte_pos,
                             Vector&               rte_los,
                             Vector&               cum_l_step,
                             ArrayOfMatrix&        TArray,
                             ArrayOfMatrix&        ext_matArray,
                             ArrayOfVector&        abs_vecArray,
                             Vector&               t_ppath,
                             //Vector&               scat_za_grid,
                             //Vector&               scat_aa_grid,
                             Tensor3&              ext_mat,
                             Matrix&               abs_vec,
                             Numeric&              rte_pressure,
                             Numeric&              rte_temperature,
                             Vector&               rte_vmr_list,
                             Matrix&               iy,
                             //Matrix&               i_space,
                             //Matrix&               ground_emission,
                             //Matrix&               ground_los, 
                             //Tensor4&              ground_refl_coeffs,
                             Matrix&               pnd_ppath,
                             const Agenda&         ppath_step_agenda,
                             const Index&          atmosphere_dim,
                             const Vector&         p_grid,
                             const Vector&         lat_grid,
                             const Vector&         lon_grid,
                             const Tensor3&        z_field,
                             const Matrix&         r_geoid,
                             const Matrix&         z_surface,
                             const ArrayOfIndex&   cloudbox_limits,
                             const Index&          record_ppathcloud,
                             const Index&          record_ppath,
                             const Agenda&         opt_prop_gas_agenda,
                             const Agenda&         abs_scalar_gas_agenda,
                             const Index&          f_index,
                             const Index&          stokes_dim,
                             const Tensor3&        t_field,
                             const Tensor4&        vmr_field,
                             const Agenda&         rte_agenda,
                             const Agenda&         iy_space_agenda,
                             const Agenda&         surface_prop_agenda,
                             const Agenda&         iy_cloudbox_agenda,
                             const Vector&         f_grid,
                             const Index&          photon_number,
                             const Index&          scattering_order,
                             const Tensor4&        pnd_field,
                             const ArrayOfSingleScatteringData& scat_data_mono,
                             const Index& z_field_is_1D
);


void cloudbox_ppath_start_stepping(
                                   Ppath&          ppath,
                                   const Index&          atmosphere_dim,
                                   ConstVectorView       p_grid,
                                   ConstVectorView       lat_grid,
                                   ConstVectorView       lon_grid,
                                   ConstTensor3View      z_field,
                                   ConstMatrixView       r_geoid,
                                   ConstMatrixView       z_surface,
                                   ConstVectorView       rte_pos,
                                   ConstVectorView       rte_los,
                                   const Index& z_field_is_1D );
          

void cloudy_rt_vars_at_gp(
                       MatrixView&           ext_mat_mono,
                       VectorView&           abs_vec_mono,
                       VectorView&           pnd_vec,
                       Numeric&              temperature,
                       const Agenda&         opt_prop_gas_agenda,
                       const Agenda&         abs_scalar_gas_agenda,
                       const Index&          stokes_dim,
                       const Index&          f_index,
                       const GridPos         gp_p,
                       const GridPos         gp_lat,
                       const GridPos         gp_lon,
                       const ConstVectorView    p_grid_cloud,
                       const ConstVectorView    lat_grid_cloud,
                       const ConstVectorView    lon_grid_cloud,
                       const ConstTensor3View   t_field_cloud,
                       const ConstTensor4View   vmr_field_cloud,
                       const Tensor4&           pnd_field,
                       const ArrayOfSingleScatteringData& scat_data_mono,
                       const ArrayOfIndex&                cloudbox_limits,
                       const Vector&            rte_los
                       );

void cum_l_stepCalc(
                      Vector& cum_l_step,
                      const Ppath& ppath
                      );

void findZ11max(Vector& Z11maxvector,
        const  ArrayOfSingleScatteringData& scat_data_mono);

bool is_anyptype30(const ArrayOfSingleScatteringData& scat_data_mono);

void iwp_cloud_opt_pathCalc(
                            Numeric& iwp,
                            Numeric& cloud_opt_path,
                            //input
                            const Vector&         rte_pos,
                            const Vector&         rte_los,
                            const Agenda&         ppath_step_agenda,
                            const Vector&         p_grid,
                            const Vector&         lat_grid, 
                            const Vector&         lon_grid, 
                            const Matrix&         r_geoid, 
                            const Matrix&         z_surface,
                            const Tensor3&        z_field, 
                            const Tensor3&        t_field, 
                            const Tensor4&        vmr_field, 
                            const ArrayOfIndex&   cloudbox_limits, 
                            const Tensor4&        pnd_field,
                            const ArrayOfSingleScatteringData& scat_data_mono,
                            const Vector&          particle_masses
                            );


void matrix_exp_p30(MatrixView& M,
                    ConstMatrixView& A);

void mcPathTrace(MatrixView&           evol_op,
                 VectorView& abs_vec_mono,
                 Numeric& temperature,
                 MatrixView& ext_mat_mono,
                 Rng&                  rng,
                 Vector&               rte_pos,
                 Vector&               rte_los,
                 Vector& pnd_vec,
                 Numeric&    g,
                 bool&                 left_cloudbox,
                 const Agenda&         opt_prop_gas_agenda,
                 const Agenda&         abs_scalar_gas_agenda,
                 const Index&          stokes_dim,
                 const Index&          f_index,
                 const Vector&         p_grid,
                 const Vector&         lat_grid,
                 const Vector&         lon_grid,
                 const Tensor3&        z_field,
                 const Matrix&         r_geoid,
                 const Matrix&         z_surface,
                 const Tensor3&   t_field,
                 const Tensor4&   vmr_field,
                 const ArrayOfIndex&   cloudbox_limits,
                 const Tensor4&   pnd_field,
                 const ArrayOfSingleScatteringData& scat_data_mono,
                 const Index& z_field_is_1D);

void mcPathTraceGeneral(MatrixView&           evol_op,
                        Vector&               abs_vec_mono,
                        Numeric&              temperature,
                        MatrixView&           ext_mat_mono,
                        Rng&                  rng,
                        Vector&               rte_pos,
                        Vector&               rte_los,
                        Vector&               pnd_vec,
                        Numeric&              g,
                        Ppath&                ppath_step,
                        Index&                termination_flag,
                        bool&                 inside_cloud,
                        //Numeric&              rte_pressure,
                        //Vector&               rte_vmr_list,
                        const Agenda&         opt_prop_gas_agenda,
                        const Agenda&         abs_scalar_gas_agenda,
                        const Index&          stokes_dim,
                        const Index&          f_index,
                        const Vector&         p_grid,
                        const Vector&         lat_grid,
                        const Vector&         lon_grid,
                        const Tensor3&        z_field,
                        const Matrix&         r_geoid,
                        const Matrix&         z_surface,
                        const Tensor3&        t_field,
                        const Tensor4&        vmr_field,
                        const ArrayOfIndex&   cloudbox_limits,
                        const Tensor4&        pnd_field,
                        const ArrayOfSingleScatteringData& scat_data_mono,
                        const Index&          z_field_is_1D);

void mcPathTraceIPA(MatrixView&           evol_op,
                    Vector&               abs_vec_mono,
                    Numeric&              temperature,
                    MatrixView&           ext_mat_mono,
                    Rng&                  rng,
                    Vector&               rte_pos,
                    Vector&               rte_los,
                    Vector&               pnd_vec,
                    Numeric&              g,
                    Index&                termination_flag,
                    bool&                 inside_cloud,
                    //Numeric&              rte_pressure,
                    //Vector&               rte_vmr_list,
                    const Agenda&         opt_prop_gas_agenda,
                    const Agenda&         abs_scalar_gas_agenda,
                    const Index&          stokes_dim,
                    const Index&          f_index,
                    const Vector&         p_grid,
                    const Vector&         lat_grid,
                    const Vector&         lon_grid,
                    const Tensor3&        z_field,
                    const Matrix&         r_geoid,
                    const Matrix&         z_surface,
                    const Tensor3&        t_field,
                    const Tensor4&        vmr_field,
                    const ArrayOfIndex&   cloudbox_limits,
                    const Tensor4&        pnd_field,
                    const ArrayOfSingleScatteringData& scat_data_mono,
                    const Index&          z_field_is_1D,
                    const Ppath&          ppath);

void montecarloGetIncoming(
                           Matrix&               iy,
                           Vector&               rte_pos,
                           Vector&               rte_los,
                           Ppath&                ppath,
                           //Vector&               ppath_p,
                           //Vector&               ppath_t,
                           //Matrix&               ppath_vmr,
                           const Agenda&         ppath_step_agenda,
                           const Agenda&         rte_agenda,
                           const Agenda&         iy_space_agenda,
                           const Agenda&         surface_prop_agenda,
                           const Agenda&         iy_cloudbox_agenda,
                           const Vector&         p_grid,
                           const Vector&         lat_grid,
                           const Vector&         lon_grid,
                           const Tensor3&        z_field,
                           const Tensor3&        t_field,
                           const Tensor4&        vmr_field,
                           const Matrix&         r_geoid,
                           const Matrix&         z_surface,
                           const ArrayOfIndex&   cloudbox_limits,
                           const Index&          atmosphere_dim,
                           const Vector&         f_grid,
                           const Index&          stokes_dim
                           );

Numeric opt_depth_calc(
                       Tensor3& ext_mat,
                       Matrix&  abs_vec,
                       Numeric&   rte_pressure,
                       Numeric&   rte_temperature,
                       Vector&    rte_vmr_list,
                       const Ppath&     ppath,
                       const Agenda& opt_prop_gas_agenda,
                       const Agenda& abs_scalar_gas_agenda,
                       const Index&     f_index,
                       const Vector&    p_grid,
                       const Vector&    lat_grid,
                       const Vector&    lon_grid,
                       const Tensor3&   t_field,
                       const Tensor4&   vmr_field,
                       const Index&     atmosphere_dim);

void opt_propCalc(
                  MatrixView& K,
                  VectorView& K_abs,
                  const Numeric za,
                  const Numeric aa,
                  const ArrayOfSingleScatteringData& scat_data_mono,
                  const Index&          stokes_dim,
                  const VectorView& pnd_vec,
                  const Numeric& rte_temperature
                  );

void opt_propExtract(
                     MatrixView& K_spt,
                     VectorView& K_abs_spt,
                     const SingleScatteringData& scat_data,
                     const Numeric& za,
                     const Numeric& aa,
                     const Numeric& rte_temperature,
                     const Index& stokes_dim
                     );

void pha_mat_singleCalc(
                        MatrixView& Z,                  
                        Numeric za_sca, 
                        Numeric aa_sca, 
                        Numeric za_inc, 
                        Numeric aa_inc,
                        const ArrayOfSingleScatteringData& scat_data_mono,
                        const Index&          stokes_dim,
                        const VectorView& pnd_vec,
                        const Numeric& rte_temperature
                        );

void pha_mat_singleExtract(
                           MatrixView& Z_spt,
                           const SingleScatteringData& scat_data,
                           const Numeric& za_sca,
                           const Numeric& aa_sca,
                           const Numeric& za_inc,
                           const Numeric& aa_inc,
                           const Numeric& rte_temperature,
                           const Index& stokes_dim
                           );


void Sample_los (
                   VectorView& new_rte_los,
                   Numeric& g_los_csc_theta,
                   MatrixView& Z,
                   Rng& rng,
                   const VectorView& rte_los,
                   const ArrayOfSingleScatteringData& scat_data_mono,
                   const Index&          stokes_dim,
                   const VectorView& pnd_vec,
                   const bool& anyptype30,
                   const VectorView& Z11maxvector,
                   Numeric Csca,
                   const Numeric& rte_temperature
                   );

void Sample_ppathlengthLOS (
                         Numeric& pathlength, 
                         Numeric& g,
                         Rng& rng,
                         const ArrayOfMatrix& TArray,
                         const ConstVectorView& cum_l_step
                         );


void TArrayCalc(
                //output
                ArrayOfMatrix& TArray,
                ArrayOfMatrix& ext_matArray,
                ArrayOfVector& abs_vecArray,
                Vector& t_ppath,
                Tensor3& ext_mat,
                Matrix& abs_vec,
                Numeric&   rte_pressure,
                Numeric&   rte_temperature,
                Vector&    rte_vmr_list,
                Matrix&  pnd_ppath,
                //input
                const Ppath& ppath,
                const Agenda& opt_prop_gas_agenda,
                const Agenda& abs_scalar_gas_agenda,
                const Index& f_index,
                const Index& stokes_dim,
                const ConstVectorView&    p_grid_cloud,
                const ConstVectorView&    lat_grid_cloud,
                const ConstVectorView&    lon_grid_cloud,
                const ConstTensor3View&   t_field_cloud,
                const ConstTensor4View&   vmr_field_cloud,
                const Tensor4&   pnd_field,
                const ArrayOfSingleScatteringData& scat_data_mono,
                const ArrayOfIndex& cloudbox_limits
                );

#endif  // montecarlo_h
