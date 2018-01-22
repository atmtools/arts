/* Copyright (C) 2016 Jana Mendrok <jana.mendrok@gmail.com>

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
  \file   rt4.h
  \author Jana Mendrok <jana.mendrok@gmail.com>
  \date   2016-05-24
  
  \brief  Contains functions related to application of scattering solver RT4.
  
*/

#ifndef rt4_h
#define rt4_h

#ifdef ENABLE_RT4

#include "messages.h"

void check_rt4_input( // Output
                      Index& nhstreams,
                      Index& nhza,
                      Index& nummu,
                      // Input
                      const Index& cloudbox_on,
                      const Index& atmfields_checked,
                      const Index& atmgeom_checked,
                      const Index& cloudbox_checked,
                      const Index& scat_data_checked,
                      const ArrayOfIndex& cloudbox_limits, 
                      const ArrayOfArrayOfSingleScatteringData& scat_data,
                      const Index& atmosphere_dim,
                      const Index& stokes_dim,
                      const Index& nstreams,
                      const String& quad_type,
                      const Index& add_straight_angles,
                      const Index& pnd_ncols );

void get_quad_angles( // Output
                      VectorView mu_values,
                      VectorView quad_weights,
                      Vector& scat_za_grid,
                      Vector& scat_aa_grid,
                      // Input
                      const String& quad_type,
                      const Index& nhstreams,
                      const Index& nhza,
                      const Index& nummu );

void get_rt4surf_props( // Output
                        Vector& ground_albedo,
                        Tensor3& ground_reflec,
                        ComplexVector& ground_index,
                        // Input
                        ConstVectorView f_grid,
                        const String& ground_type,
                        const Numeric& surface_skin_t,
                        ConstVectorView surface_scalar_reflectivity,
                        ConstTensor3View surface_reflectivity,
                        const GriddedField3& surface_complex_refr_index,
                        const Index& stokes_dim );

void run_rt4( Workspace& ws,
              // Output
              Tensor7& doit_i_field,
              Vector& scat_za_grid,
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
              const Index& stokes_dim,
              const Index& nummu,
              const Index& nhza,
              const String& ground_type,
              const Numeric& surface_skin_t,
              ConstVectorView ground_albedo,
              ConstTensor3View ground_reflec,
              ConstComplexVectorView ground_index,
              ConstTensor5View surf_refl_mat,
              ConstTensor3View surf_emis_vec,
              const Agenda& surface_rtprop_agenda,
              const Numeric& surf_altitude,
              const String& quad_type,
              Vector& mu_values,
              ConstVectorView quad_weights,
              const Index& auto_inc_nstreams,
              const Index& robust,
              const Index& za_interp_order,
              const Index& cos_za_interp,
              const String& pfct_method,
              const Index& pfct_aa_grid_size,
              const Numeric& pfct_threshold,
              const Numeric& max_delta_tau,
              const Index& new_optprop,
              const Verbosity& verbosity );

void scat_za_grid_adjust( // Output
                          Vector& scat_za_grid,
                          // Input
                          ConstVectorView mu_values,
                          const Index& nummu );

void gas_optpropCalc( Workspace& ws,
                      //Output
                      VectorView gas_extinct,
                      //Input
                      const Agenda& propmat_clearsky_agenda,
                      ConstTensor3View t_field, 
                      ConstTensor4View vmr_field,
                      ConstVectorView p_grid,
                      ConstVectorView f_mono );

void par_optpropCalc( //Output
                      Tensor4View emis_vector,
                      Tensor5View extinct_matrix,
                      //VectorView scatlayers,
                      //Input
                      const ArrayOfArrayOfSingleScatteringData& scat_data,
                      const Vector& scat_za_grid,
                      const Index& f_index,
                      ConstTensor4View pnd_field,
                      ConstTensor3View t_field,
                      const ArrayOfIndex& cloudbox_limits,
                      const Index& stokes_dim,
                      const Index& nummu,
                      const Verbosity& verbosity );

void par_optpropCalc2( //Output
                      Tensor5View emis_vector,
                      Tensor6View extinct_matrix,
                      //VectorView scatlayers,
                      //Input
                      const ArrayOfArrayOfSingleScatteringData& scat_data,
                      const Vector& scat_za_grid,
                      const Index& f_index,
                      ConstTensor4View pnd_field,
                      ConstTensor3View t_field,
                      const ArrayOfIndex& cloudbox_limits,
                      const Index& stokes_dim );

void sca_optpropCalc( //Output
                      Tensor6View scatter_matrix,
                      Index& pfct_failed,
                      //Input
                      ConstTensor4View emis_vector,
                      ConstTensor5View extinct_matrix,
                      const Index& f_index,
                      const ArrayOfArrayOfSingleScatteringData& scat_data,
                      ConstTensor4View pnd_field,
                      const Index& stokes_dim,
                      const Vector& scat_za_grid,
                      ConstVectorView quad_weights,
                      const String& pfct_method,
                      const Index& pfct_aa_grid_size,
                      const Numeric& pfct_threshold,
                      const Index& auto_inc_nstreams,
                      const Verbosity& verbosity );

void surf_optpropCalc( Workspace& ws,
                       //Output
                       Tensor5View surf_refl_mat,
                       Tensor3View surf_emis_vec,
                       //Input
                       const Agenda& surface_rtprop_agenda,
                       ConstVectorView f_grid,
                       ConstVectorView scat_za_grid, 
                       ConstVectorView mu_values,
                       ConstVectorView quad_weights,
                       const Index& stokes_dim,
                       const Numeric& surf_alt );

void rt4_test( Tensor4& out_rad,
               const String& datapath,
               const Verbosity& verbosity );


extern "C" {

    void radtrano_( const Index&   nstokes,
                    const Index&   nummu,
                    const Index&   nuummu,
                    const Numeric& max_delta_tau,
                    const char*    quad_type,
                    const Numeric& ground_temp,
                    const char*    ground_type,
                    const Numeric& ground_albedo,
                    const Complex& ground_index,
                    const Numeric* ground_reflec,
                    const Numeric* sre_data,
                    const Numeric* sem_data,
                    const Numeric& sky_temp,
                    const Numeric& wavelength,
                    const Index&   num_layers,
                    const Numeric* height,
                    const Numeric* temperatures,
                    const Numeric* gas_extinct,
                    const Index&   num_scatlayers,
                    const Numeric* scatlayers,
                    const Numeric* ext_data,
                    const Numeric* abs_data,
                    const Numeric* sca_data,
                    //const Index&   noutlevels,
                    //const Index*   outlevels,
                    Numeric* mu_values,
                    Numeric* up_rad,
                    Numeric* down_rad
                    );

    void double_gauss_quadrature_( const Index& nummu,
                                   Numeric* mu_values,
                                   Numeric* quad_weights
                                   );

    void lobatto_quadrature_( const Index& nummu,
                              Numeric* mu_values,
                              Numeric* quad_weights
                              );

    void gauss_legendre_quadrature_( const Index& nummu,
                                     Numeric* mu_values,
                                     Numeric* quad_weights
                                     );

    void planck_function_( const Numeric& temp,
                           const char* units,
                           const Numeric& wavelength,
                           Numeric& planck
                           );
}

#endif /* ENABLE_RT4 */

#endif /* rt4_h */

