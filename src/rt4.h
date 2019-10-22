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

//! Checks that input of RT4Calc* is sane.
/*!

  \param[out]  nhstreams Number of single hemisphere streams (quadrature angles)
  \param[out]  nhza Number of single hemisphere additional angles with RT output
  \param[out]  nummu Total number of single hemisphere angles with RT output
  \param[in]   cloudbox_on Flag to activate the cloud box
  \param[in]   atmfields_checked OK-flag for atmospheric grids and (physical)
               fields
  \param[in]   atmgeom_checked OK-flag for the geometry of the model atmosphere
  \param[in]   cloudbox_checked OK-flag for variables associated with the
               cloudbox
  \param[in]   scat_data_checked OK-flag for scat_data
  \param[in]   cloudbox_limits Cloudbox limits
  \param[in]   scat_data Array of single scattering data
  \param[in]   atmosphere_dim Dimension of atmosphere
  \param[in]   stokes_dim Dimension of Stokes vector
  \param[in]   nstreams Total number of quadrature angles (both hemispheres).
  \param[in]   quad_type Quadrature method.
  \param[in]   add_straight_angles Flag whether to include nadir and zenith as
               explicit directions
  \param[in]   pnd_ncols Number of columns (latitude points) in *pnd_field*.

  \author Jana Mendrok
  \date   2017-02-22
*/
void check_rt4_input(  // Output
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
    const Index& pnd_ncols);

//! Derive the quadrature angles
/*!
  Derive the quadrature angles related to selected RT4 quadrature type and set
  scat_za_grid accordingly.

  \param[out]  mu_values Quadrature angle cosines.
  \param[out]  quad_weights Quadrature weights associated with mu_values.
  \param[out]  scat_za_grid Zenith angle grid
  \param[out]  scat_aa_grid Azimuth angle grid
  \param[in]   quad_type Quadrature method.
  \param[in]   nhstreams Number of single hemisphere streams (quadrature angles)
  \param[in]   nhza Number of single hemisphere additional angles with RT output
  \param[in]   nummu Total number of single hemisphere angles with RT output

  \author Jana Mendrok
  \date   2017-02-22
*/
void get_quad_angles(  // Output
    VectorView mu_values,
    VectorView quad_weights,
    Vector& scat_za_grid,
    Vector& scat_aa_grid,
    // Input
    const String& quad_type,
    const Index& nhstreams,
    const Index& nhza,
    const Index& nummu);


//! Derive surface property input for RT4's proprietary surface handling
/*!
  Derive surface property input for RT4's proprietary surface handling depending
  on surface reflection type.

  \param[out]  ground_albedo Scalar surface albedo (for ground_type=L).
  \param[out]  ground_reflec Vector surface relfectivity (for ground_type=S).
  \param[out]  ground_index Surface complex refractive index (for ground_type=F).
  \param[in]   f_grid Frequency grid
  \param[in]   ground_type Surface reflection type flag.
  \param[in]   surface_skin_t Surface skin temperature
  \param[in]   surface_scalar_reflectivity Surface scalar reflectivity
               (used with ground_type=L).
  \param[in]   surface_reflectivity  Surface reflectivity
               (used with ground_type=S).
  \param[in]   surface_complex_refr_index  Surface complex refractive index
               (used with ground_type=F).
  \param[in]   stokes_dim Dimension of Stokes vector

  \author Jana Mendrok
  \date   2017-02-22
*/
void get_rt4surf_props(  // Output
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
    const Index& stokes_dim);

//! Runs RT4
/*!
  Prepares actual input variables for RT4, runs it, and sorts the output into
  doit_i_field.

  \param[in,out] ws Current workspace
  \param[out]    doit_i_field Radiation field
  \param[out]    scat_za_grid Zenith angle grid
  \param[in]     f_grid Frequency grid
  \param[in]     p_grid Pressure rid
  \param[in]     z_profile Profile of geometrical altitudes
  \param[in]     t_profile Temperature profile
  \param[in]     vmr_profiles VMR profiles
  \param[in]     pnd_profiles PND profiles
  \param[in]     scat_data Array of single scattering data
  \param[in]     propmat_clearsky_agenda calculates the absorption coefficient
                 matrix
  \param[in]     cloudbox_limits Cloudbox limits
  \param[in]     stokes_dim Dimension of Stokes vector
  \param[in]     nummu Total number of single hemisphere angles with RT output.
  \param[in]     nhza Number of single hemisphere additional angles with RT output.
  \param[in]     ground_type Surface reflection type flag.
  \param[in]     surface_skin_t Surface scalar temperature
  \param[in]     ground_albedo Scalar surface albedo (for ground_type=L).
  \param[in]     ground_reflec Vector surface relfectivity (for ground_type=S).
  \param[in]     ground_index Surface complex refractive index (for ground_type=F).
  \param[in]     surf_refl_mat Surface reflection matrix (for ground_type=A).
  \param[in]     surf_emis_vec Surface emission vector (for ground_type=A).
  \param[in]     surface_rtprop_agenda Provides radiative properties of the
                 surface
  \param[in]     surf_altitude Surface altitude (lowest level of z_field)
  \param[in]     quad_type Quadrature method.
  \param[in]     mu_values Quadrature angle cosines.
  \param[in]     quad_weights Quadrature weights associated with mu_values.
  \param[in]     auto_inc_nstreams Flag whether to internally increase nstreams
  \param[in]     robust if auto_inc_nstreams>0,flag whether to not fail even if
                 scattering matrix norm is not preserved when maximum stream
                 number is reached
  \param[in]     za_interp_order if auto_inc_nstreams>0, polar angle
                 interpolation order
  \param[in]     cos_za_interp if auto_inc_nstreams>0, flag whether to do polar
                 angle interpolation in cosine (='mu') space
  \param[in]     pfct_method Flag which method to apply to derive phase function
  \param[in]     pfct_aa_grid_size Number of azimuthal angle grid points to
                 consider in Fourier series decomposition of scattering matrix
  \param[in]     pfct_threshold Requested scatter_matrix norm accuracy
                 (in terms of single scat albedo).
  \param[in]     max_delta_tau Maximum optical depth of infinitesimal layer
  \param[in]     verbosity Verbosity setting

  \author Jana Mendrok
  \date   2017-02-22
*/
void run_rt4(Workspace& ws,
             // Output
             Tensor7& doit_i_field,
             Vector& scat_za_grid,
             // Input
             ConstVectorView f_grid,
             ConstVectorView p_grid,
             ConstVectorView z_profile,
             ConstVectorView t_profile,
             ConstMatrixView vmr_profiles,
             ConstMatrixView pnd_profiles,
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
             const Verbosity& verbosity);

//! Reset scat_za_grid such that it is consistent with ARTS
/*!
  Reset scat_za_grid such that it is consistent with ARTS' scat_za_grid
  requirements (instead of with RT4 as in input state).

  \param[out]  scat_za_grid Zenith angle grid
  \param[in]   mu_values Quadrature angle cosines.
  \param[in]   nummu Number of single hemisphere polar angles.

  \author Jana Mendrok
  \date   2017-02-22
*/
void scat_za_grid_adjust(  // Output
    Vector& scat_za_grid,
    // Input
    ConstVectorView mu_values,
    const Index& nummu);

//! Calculates layer averaged gaseous extinction
/*!
  Calculates layer averaged gaseous extinction (gas_extinct). This variable is
  required as input for the RT4 subroutine.

  \param[in,out]  ws Current workspace
  \param[out]     gas_extinct Layer averaged gas extinction for all layers
  \param[in]      propmat_clearsky_agenda propmat_clearsky_agenda calculates
                  the absorption coefficient matrix
  \param[in]      t_profile Temperature profile
  \param[in]      vmr_profiles VMR profiles
  \param[in]      p_grid Pressure grid
  \param[in]      f_mono Frequency (single entry vector)

  \author Jana Mendrok
  \date   2016-08-08
*/
void gas_optpropCalc(Workspace& ws,
                     //Output
                     VectorView gas_extinct,
                     //Input
                     const Agenda& propmat_clearsky_agenda,
                     ConstVectorView t_field,
                     ConstMatrixView vmr_field,
                     ConstVectorView p_grid,
                     ConstVectorView f_mono);


//! Calculates layer averaged particle extinction and absorption
/*!
  Calculates layer averaged particle extinction and absorption (extinct_matrix
  and emis_vector)). These variables are required as input for the RT4 subroutine.

  \param[out] emis_vector Layer averaged particle absorption for all particle
              layers
  \param[out] extinct_matrix  Layer averaged particle extinction for all
              particle layers
  \param[in]  scat_data Array of single scattering data
  \param[in]  scat_za_grid Zenith angle grid
  \param[in]  f_index Index of frequency grid point handeled
  \param[in]  pnd_profiles PND profiles
  \param[in]  t_profile Temperature profile
  \param[in]  cloudbox_limits Cloudbox limits
  \param[in]  stokes_dim Dimension of Stokes vector

  \author Jana Mendrok
  \date   2016-08-08
*/
void par_optpropCalc(  //Output
                     Tensor5View emis_vector,
                     Tensor6View extinct_matrix,
                     //VectorView scatlayers,
                     //Input
                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                     const Vector& scat_za_grid,
                     const Index& f_index,
                     ConstMatrixView pnd_profiles,
                     ConstVectorView t_profile,
                     const ArrayOfIndex& cloudbox_limits,
                     const Index& stokes_dim);

//! Calculates layer (and azimuthal) averaged phase matrix
/*!
  Calculates layer (and azimuthal) averaged phase matrix (scatter_matrix). This
  variable is required as input for the RT4 subroutine.

  \param[out] scatter_matrix  Layer averaged scattering matrix (azimuth mode 0)
              for all particle layers
  \param[out] pfct_failed Flag whether norm of scatter_matrix is sufficiently
              accurate
  \param[in]  emis_vector Layer averaged particle absorption for all particle
              layers
  \param[in]  extinct_matrix Layer averaged particle extinction for all particle
              layers
  \param[in]  f_index Frequency index
  \param[in]  scat_data Array of single scattering data
              (new-type, f_grid prepared)
  \param[in]  pnd_profiles PND profiles
  \param[in]  stokes_dim Dimension of Stokes vector
  \param[in]  scat_za_grid Zenith angle grid
  \param[in]  quad_weights Quadrature weights associated with scat_za_grid
  \param[in]  pfct_method Method for scattering matrix temperature dependance
              handling
  \param[in]  pfct_aa_grid_size Number of azimuthal grid points in Fourier
              series decomposition of randomly oriented particles
  \param[in]  pfct_threshold Requested scatter_matrix norm accuracy
              (in terms of single scat albedo)
  \param[in]  auto_inc_nstreams Flag whether to internally increase nstreams
  \param[in]  verbosity Verbosity setting

  \author Jana Mendrok
  \date   2016-08-08
*/
void sca_optpropCalc(  //Output
                     Tensor6View scatter_matrix,
                     Index& pfct_failed,
                     //Input
                     ConstTensor4View emis_vector,
                     ConstTensor5View extinct_matrix,
                     const Index& f_index,
                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                     ConstMatrixView pnd_profiles,
                     const Index& stokes_dim,
                     const Vector& scat_za_grid,
                     ConstVectorView quad_weights,
                     const String& pfct_method,
                     const Index& pfct_aa_grid_size,
                     const Numeric& pfct_threshold,
                     const Index& auto_inc_nstreams,
                     const Verbosity& verbosity);

//! Calculates bidirectional surface reflection matrices and emission direction
/*!
  Calculates bidirectional surface reflection matrices and emission direction
  dependent surface emission terms as required as input for the RT4 subroutine.

  \param[in,out]  ws Current workspace
  \param[out]     surf_refl_mat Bidirectional surface reflection matrices on
                  RT4 stream directions.
  \param[out]     surf_emis_vec Directional surface emission vector on RT4
                  stream directions.
  \param[in]      surface_rtprop_agenda surface_rtprop_agenda Provides radiative
                  properties of the surface
  \param[in]      f_grid Frequency grid
  \param[in]      scat_za_grid Zenith angle grid
  \param[in]      mu_values Cosines of scat_za_grid angles.
  \param[in]      quad_weights Quadrature weights associated with mu_values.
  \param[in]      stokes_dim Dimension of Stokes vector
  \param[in]      surf_alt Surface altitude.

  \author Jana Mendrok
  \date   2017-02-09
*/
void surf_optpropCalc(Workspace& ws,
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
                      const Numeric& surf_alt);

//! Calculate radiation field using RT4
/*!
  Calculate radiation field using Evans' RT4 model (part of PolRadTran).

  This is a direct interface to the (almost orignal) RT4 FORTRAN code. No checks
  of input are made. Function is only to be called through other
  functions/methods, which have to ensure input consistency.

  \param[out] out_rad FIXMEDOC
  \param[in]  datapath FIXMEDOC
  \param[in]  verbosity Verbosity setting

  \author Jana Mendrok
  \date 2016-05-24
*/
void rt4_test(Tensor4& out_rad,
              const String& datapath,
              const Verbosity& verbosity);

extern "C" {

void radtrano_(const Index& nstokes,
               const Index& nummu,
               const Index& nuummu,
               const Numeric& max_delta_tau,
               const char* quad_type,
               const Numeric& ground_temp,
               const char* ground_type,
               const Numeric& ground_albedo,
               const Complex& ground_index,
               const Numeric* ground_reflec,
               const Numeric* sre_data,
               const Numeric* sem_data,
               const Numeric& sky_temp,
               const Numeric& wavelength,
               const Index& num_layers,
               const Numeric* height,
               const Numeric* temperatures,
               const Numeric* gas_extinct,
               const Index& num_scatlayers,
               const Numeric* scatlayers,
               const Numeric* ext_data,
               const Numeric* abs_data,
               const Numeric* sca_data,
               //const Index&   noutlevels,
               //const Index*   outlevels,
               Numeric* mu_values,
               Numeric* up_rad,
               Numeric* down_rad);

void double_gauss_quadrature_(const Index& nummu,
                              Numeric* mu_values,
                              Numeric* quad_weights);

void lobatto_quadrature_(const Index& nummu,
                         Numeric* mu_values,
                         Numeric* quad_weights);

void gauss_legendre_quadrature_(const Index& nummu,
                                Numeric* mu_values,
                                Numeric* quad_weights);

void planck_function_(const Numeric& temp,
                      const char* units,
                      const Numeric& wavelength,
                      Numeric& planck);
}

#endif /* ENABLE_RT4 */

#endif /* rt4_h */
