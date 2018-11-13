/* Copyright (C) 2006-2012 Claudia Emde <claudia.emde@dlr.de>

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

/**
 \file   disort.h
 \author Claudia Emde <claudia.emde@dlr.de>
 \date   Tue Feb  7 11:48:17 2006
  
 \brief  Functions for disort interface
 * 
 * 
 */

#ifndef disort_h
#define disort_h

#include "agenda_class.h"
#include "matpackIV.h"
#include "mystring.h"
#include "optproperties.h"


void check_disort_input( // Input
                         const Index& cloudbox_on,
                         const Index& atmfields_checked,
                         const Index& atmgeom_checked,
                         const Index& cloudbox_checked,
                         const Index& scat_data_checked,
                         const Index& atmosphere_dim,
                         const Index& stokes_dim,
                         const ArrayOfIndex& cloudbox_limits, 
                         const ArrayOfArrayOfSingleScatteringData& scat_data,
                         ConstVectorView scat_za_grid,
                         const Index& nstreams,
                         const String& pfct_method,
                         const Index& pnd_ncols );

void init_ifield( // Output
                  Tensor7& doit_i_field,
                  // Input
                  const Vector& f_grid,
                  const ArrayOfIndex& cloudbox_limits, 
                  const Index& nang,
                  const Index& stokes_dim );

void get_disortsurf_props( // Output
                        Vector& albedo,
                        Numeric& btemp,
                        // Input
                        ConstVectorView f_grid,
                        const Numeric& surface_skin_t,
                        ConstVectorView surface_scalar_reflectivity );

void run_disort( Workspace& ws,
              // Output
              Tensor7& doit_i_field,
              // Input
              ConstVectorView f_grid,
              ConstVectorView p_grid,
              ConstTensor3View z_field,
              ConstTensor3View t_field,
              ConstTensor4View vmr_field,
              ConstTensor4View pnd_field,
              const ArrayOfArrayOfSingleScatteringData& scat_data,
              const Agenda& propmat_clearsky_agenda, 
              const ArrayOfIndex& cloudbox_limits,
              Numeric& surface_skin_t,
              Vector& surface_scalar_reflectivity,
              ConstVectorView scat_za_grid,
              const Index& nstreams,
              const Index& do_deltam,
              const String& pfct_method,
              const Verbosity& verbosity );

void run_disort2( Workspace& ws,
              // Output
              Tensor7& doit_i_field,
              // Input
              ConstVectorView f_grid,
              ConstVectorView p_grid,
              ConstTensor3View z_field,
              ConstTensor3View t_field,
              ConstTensor4View vmr_field,
              ConstTensor4View pnd_field,
              const ArrayOfArrayOfSingleScatteringData& scat_data,
              const Agenda& propmat_clearsky_agenda, 
              const ArrayOfIndex& cloudbox_limits,
              Numeric& surface_skin_t,
              Vector& surface_scalar_reflectivity,
              ConstVectorView scat_za_grid,
              const Index& nstreams,
              const Index& do_deltam,                  
              const Index& Npfct,
              const Verbosity& verbosity );

void dtauc_ssalbCalc( Workspace &ws,
                      VectorView dtauc,
                      VectorView ssalb,
                      const ArrayOfArrayOfSingleScatteringData& scat_data,
                      const Index& f_index,
                      const Agenda& propmat_clearsky_agenda,
                      ConstTensor4View pnd_field,
                      ConstTensor3View t_field,
                      ConstTensor3View z_field, 
                      ConstTensor4View vmr_field,
                      ConstVectorView p_grid,
                      const ArrayOfIndex& cloudbox_limits,
                      ConstVectorView f_mono,
                      const Verbosity& verbosity );

void get_gasoptprop( Workspace &ws,
                     MatrixView ext_bulk_gas,
                     const Agenda& propmat_clearsky_agenda,
                     ConstVectorView t_field,
                     ConstMatrixView vmr_field,
                     ConstVectorView p_grid,
                     ConstVectorView f_grid );

void get_paroptprop( MatrixView ext_bulk_par,
                     MatrixView abs_bulk_par,
                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                     ConstMatrixView pnd_field,
                     ConstVectorView t_field,
                     ConstVectorView p_grid,
                     const ArrayOfIndex& cloudbox_limits,
                     ConstVectorView f_grid );

void get_dtauc_ssalb( MatrixView dtauc,
                      MatrixView ssalb,
                      ConstMatrixView ext_bulk_gas,
                      ConstMatrixView ext_bulk_par,
                      ConstMatrixView abs_bulk_par,
                      ConstVectorView z_field );

void get_angs( Vector& pfct_angs,
               const ArrayOfArrayOfSingleScatteringData& scat_data,
               const Index& Npfct );

void get_parZ( Tensor3& pha_bulk_par,
               const ArrayOfArrayOfSingleScatteringData& scat_data,
               ConstMatrixView pnd_field,
               ConstVectorView t_field,
               ConstVectorView pfct_angs,
               const ArrayOfIndex& cloudbox_limits );

void get_pfct( Tensor3& pfct_bulk_par,
               ConstTensor3View& pha_bulk_par,
               ConstMatrixView ext_bulk_par,
               ConstMatrixView abs_bulk_par,
               const ArrayOfIndex& cloudbox_limits );

void get_pmom( Tensor3View pmom,
               ConstTensor3View pfct_bulk_par,
               ConstVectorView pfct_angs,
               const Index& Nlegendre );

void dtauc_ssalbCalc2( Workspace &ws,
                      MatrixView dtauc,
                      MatrixView ssalb,
                      const Agenda& propmat_clearsky_agenda,
                      const ArrayOfArrayOfSingleScatteringData& scat_data,
                      ConstMatrixView pnd_field,
                      ConstVectorView t_field,
                      ConstVectorView z_field, 
                      ConstMatrixView vmr_field,
                      ConstVectorView p_grid,
                      const ArrayOfIndex& cloudbox_limits,
                      ConstVectorView f_grid );

void phase_functionCalc2( //Output
                          MatrixView phase_function,
                          //Input
                          const ArrayOfArrayOfSingleScatteringData& scat_data,
                          const Index& f_index,
                          ConstTensor4View pnd_field,
                          ConstTensor3View t_field,
                          const ArrayOfIndex& cloudbox_limits,
                          const Index& pfct_za_grid_size,
                          const Verbosity& verbosity );

void phase_functionCalc( //Output
                         MatrixView phase_function,
                         //Input
                         const ArrayOfArrayOfSingleScatteringData& scat_data,
                         const Index& f_index,
                         ConstTensor4View pnd_field,
                         const ArrayOfIndex& cloudbox_limits,
                         const String pfct_method );

void pmomCalc2( //Output
                MatrixView pmom,
                //Input
                ConstMatrixView phase_function, 
                ConstVectorView scat_angle_grid,
                const Index n_legendre,
                const Verbosity& verbosity );

void pmomCalc( //Output
               MatrixView pmom,
               //Input
               ConstMatrixView phase_function, 
               ConstVectorView scat_angle_grid,
               const Index n_legendre,
               const Verbosity& verbosity );

void surf_albedoCalc( Workspace& ws, 
                      //Output
                      VectorView albedo,
                      Numeric& btemp,
                      //Input
                      const Agenda& surface_rtprop_agenda,
                      ConstVectorView f_grid,
                      ConstVectorView scat_za_grid,
                      const Numeric& surf_alt,
                      const Verbosity& verbosity );

#ifdef ENABLE_DISORT
void get_cb_inc_field( Workspace&      ws,
                       Matrix&         cb_inc_field,
                       const Agenda&   iy_main_agenda,
                       const Tensor3&  z_field,
                       const Tensor3&  t_field,
                       const Tensor4&  vmr_field,
                       const Tensor4&  nlte_field,
                       const ArrayOfIndex&   cloudbox_limits,
                       const Vector&   f_grid,
                       const Vector&   scat_za_grid,
                       const Index&    nstreams );
#else /* ENABLE_DISORT */
void get_cb_inc_field( Workspace&,
                       Matrix&,
                       const Agenda&,
                       const Tensor3&,
                       const Tensor3&,
                       const Tensor4&,
                       const Tensor4&,
                       const Index&,
                       const ArrayOfIndex&,
                       const Vector&,
                       const Vector&,
                       const Index& );

#endif /* ENABLE_DISORT */

#endif /* disort_h */

