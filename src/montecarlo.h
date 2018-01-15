/* Copyright (C) 2003-2012 Cory Davis <cory@met.ed.ac.uk>
  
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
#include "cloudbox.h"
#include "optproperties.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;

void clear_rt_vars_at_gp (Workspace&          ws,
                          MatrixView          ext_mat_mono,
                          VectorView          abs_vec_mono,
                          Numeric&            temperature,
                          const Agenda&       propmat_clearsky_agenda,
                          const Numeric&      f_mono,
                          const GridPos&      gp_p,
                          const GridPos&      gp_lat,
                          const GridPos&      gp_lon,
                          ConstVectorView     p_grid,
                          ConstTensor3View    t_field,
                          ConstTensor4View    vmr_field);

void cloudy_rt_vars_at_gp (Workspace&            ws,
                           MatrixView            ext_mat_mono,
                           VectorView            abs_vec_mono,
                           VectorView            pnd_vec,
                           Numeric&              temperature,
                           const Agenda&         propmat_clearsky_agenda,
                           const Index           stokes_dim,
                           const Numeric&        f_mono,
                           const GridPos&        gp_p,
                           const GridPos&        gp_lat,
                           const GridPos&        gp_lon,
                           ConstVectorView       p_grid_cloud,
                           ConstTensor3View      t_field_cloud,
                           ConstTensor4View      vmr_field_cloud,
                           const Tensor4&        pnd_field,
                           const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                           const ArrayOfIndex&   cloudbox_limits,
                           const Vector&         rte_los,
                           const Verbosity&      verbosity);

void cloud_atm_vars_by_gp (VectorView pressure,
                           VectorView temperature,
                           MatrixView vmr,
                           MatrixView pnd,
                           const ArrayOfGridPos& gp_p,
                           const ArrayOfGridPos& gp_lat,
                           const ArrayOfGridPos& gp_lon,
                           const ArrayOfIndex&   cloudbox_limits,
                           ConstVectorView    p_grid_cloud,
                           ConstTensor3View   t_field_cloud,
                           ConstTensor4View   vmr_field_cloud,
                           ConstTensor4View   pnd_field);

void get_ppath_transmat (Workspace&      ws,
                         MatrixView&     trans_mat,
                   const Ppath&          ppath,
                   const Agenda&         propmat_clearsky_agenda,
                   const Index           stokes_dim,
                   const Numeric&        f_mono,
                   const Vector&         p_grid,
                   const Tensor3&        t_field,
                   const Tensor4&        vmr_field,
                   const ArrayOfIndex&   cloudbox_limits,
                   const Tensor4&        pnd_field,
                   const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                   const Verbosity&      verbosity );

bool is_anyptype30 (const ArrayOfArrayOfSingleScatteringData& scat_data_mono);


void mcPathTraceGeneral (Workspace&            ws,
                         MatrixView            evol_op,
                         Vector&               abs_vec_mono,
                         Numeric&              temperature,
                         MatrixView            ext_mat_mono,
                         Rng&                  rng,
                         Vector&               rte_pos,
                         Vector&               rte_los,
                         Vector&               pnd_vec,
                         Numeric&              g,
                         Ppath&                ppath_step,
                         Index&                termination_flag,
                         bool&                 inside_cloud,
                         const Agenda&         ppath_step_agenda,
                         const Numeric&        ppath_lmax,
                         const Numeric&        ppath_lraytrace,
                         const Numeric&        taustep_limit,
                         const Agenda&         propmat_clearsky_agenda,
                         const Index           stokes_dim,
                         const Numeric&        f_mono,
                         const Vector&         p_grid,
                         const Vector&         lat_grid,
                         const Vector&         lon_grid,
                         const Tensor3&        z_field,
                         const Vector&         refellipsoid,
                         const Matrix&         z_surface,
                         const Tensor3&        t_field,
                         const Tensor4&        vmr_field,
                         const ArrayOfIndex&   cloudbox_limits,
                         const Tensor4&        pnd_field,
                         const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                         const Verbosity&      verbosity);

void mcPathTraceRadar (Workspace&            ws,
                       MatrixView            evol_op,
                       Vector&               abs_vec_mono,
                       Numeric&              temperature,
                       MatrixView            ext_mat_mono,
                       Rng&                  rng,
                       Vector&               rte_pos,
                       Vector&               rte_los,
                       Vector&               pnd_vec,
                       Numeric&              stot,
                       Numeric&              ttot,
                       Ppath&                ppath_step,
                       Index&                termination_flag,
                       bool&                 inside_cloud,
                       const Agenda&         ppath_step_agenda,
                       const Numeric&        ppath_lmax,
                       const Numeric&        ppath_lraytrace,
                       const Agenda&         propmat_clearsky_agenda,
                       const bool&           anyptype30,
                       const Index           stokes_dim,
                       const Numeric&        f_mono,
                       const Vector&         I_inc,
                       const Vector&         p_grid,
                       const Vector&         lat_grid,
                       const Vector&         lon_grid,
                       const Tensor3&        z_field,
                       const Vector&         refellipsoid,
                       const Matrix&         z_surface,
                       const Tensor3&        t_field,
                       const Tensor4&        vmr_field,
                       const ArrayOfIndex&   cloudbox_limits,
                       const Tensor4&        pnd_field,
                       const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                       const Verbosity&      verbosity);

void opt_propCalc (MatrixView      K,
                   VectorView      K_abs,
                   const Numeric   za,
                   const Numeric   aa,
                   const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                   const Index     stokes_dim,
                   ConstVectorView pnd_vec,
                   const Numeric   rtp_temperature,
                   const Verbosity& verbosity);

void opt_propExtract (MatrixView     K_spt,
                      VectorView     K_abs_spt,
                      const SingleScatteringData& scat_data_single,
                      const Numeric  za,
                      const Numeric  aa,
                      const Numeric  rtp_temperature,
                      const Index    stokes_dim,
                      const Verbosity& verbosity);

void pha_mat_singleCalc (MatrixView Z,                  
                         const Numeric    za_sca, 
                         const Numeric    aa_sca, 
                         const Numeric    za_inc, 
                         const Numeric    aa_inc,
                         const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                         const Index      stokes_dim,
                         ConstVectorView  pnd_vec,
                         const Numeric    rtp_temperature,
                         const Verbosity& verbosity);

void pha_mat_singleCalcScatElement(
                        Tensor3View Z,                  
                        const Numeric    za_sca, 
                        const Numeric    aa_sca, 
                        const Numeric    za_inc, 
                        const Numeric    aa_inc,
                        const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                        const Index      stokes_dim,
                        ConstVectorView  pnd_vec,
                        const Numeric    rtp_temperature,
                        const Verbosity& verbosity);

void pha_mat_singleExtract (MatrixView Z_spt,
                            const SingleScatteringData& scat_data_single,
                            const Numeric za_sca,
                            const Numeric aa_sca,
                            const Numeric za_inc,
                            const Numeric aa_inc,
                            const Numeric rtp_temperature,
                            const Index   stokes_dim,
                            const Verbosity& verbosity);

void Sample_los (VectorView       new_rte_los,
                 Numeric&         g_los_csc_theta,
                 MatrixView       Z,
                 Rng&             rng,
                 ConstVectorView  rte_los,
                 const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                 const Index      stokes_dim,
                 ConstVectorView  pnd_vec,
                 ConstVectorView  Z11maxvector,
                 const Numeric    Csca,
                 const Numeric    rtp_temperature,
                 const Verbosity& verbosity);

void Sample_los_uniform (VectorView    new_rte_los,
                         Rng&          rng);

#endif  // montecarlo_h
