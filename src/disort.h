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
 * @file   disort.h
 * @author Claudia Emde <claudia.emde@dlr.de>
 * @date   Tue Feb  7 11:48:17 2006
  
 * @brief  Functions for disort interface.
 * 
 */

#ifndef disort_h
#define disort_h

#include "agenda_class.h"
#include "matpackIV.h"
#include "mystring.h"
#include "optproperties.h"

/** check_disort_input. *** FIXMEDOC *** in disort.cc, line 197
 *
 * Checks that input of DisortCalc* is sane.
 *
 * @param[in]  cloudbox_on           as the WSV.
 * @param[in]  disort_is_initialized as the WSV.
 * @param[in]  atmfields_checked     as the WSV.
 * @param[in]  atmgeom_checked       as the WSV.
 * @param[in]  cloudbox_checked      as the WSV.
 * @param[in]  scat_data             as the WSV.
 * @param[in]  scat_za_grid          as the WSV.
 * @param[in]  nstreams              Number of quadrature angles (both hemispheres).
 * @param[in]  pfct_method           see DisortCalc doc.
 *
 * @author     Jana Mendrok
 * @date       2017-02-23
 */
void check_disort_input(  // Input
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
    const String& pfct_method);

/** init_ifield.
 *
 * Initialize doit_i_field with the right size and NaN values.
 *
 * @param[out] doit_i_field       As the WSV.
 * @param[in]  f_grid             As the WSV.
 * @param[in]  cloudbox_limits    As the WSV.
 * @param[in]  nang               Total number of angles with RT output.
 * @param[in]  stokes_dim         As the WSV.
 *
 * @author     Jana Mendrok
 * @date       2017-03-06
 */
void init_ifield(  // Output
    Tensor7& doit_i_field,
    // Input
    const Vector& f_grid,
    const ArrayOfIndex& cloudbox_limits,
    const Index& nang,
    const Index& stokes_dim);

/** get_disortsurf_props. *** FIXMEDOC *** input/output
 *
 * Derive surface property input for RT4's proprietary surface handling depending
 * on surface reflection type.
 *
 * @param[out] ground_albedo               Scalar surface albedo (for ground_type=L).
 * @param[out] ground_reflec               Vector surface relfectivity (for ground_type=S).
 * @param[out] ground_index                Surface complex refractive index (for ground_type=F).
 * @param[out] btemp                       Surface temperature
 * @param[in]  f_grid                      As the WSV.
 * @param[in]  ground_type                 Surface reflection type flag.
 * @param[in]  surface_skin_t              As the WSV.
 * @param[in]  surface_scalar_reflectivity As the WSV (used with ground_type=L).
 * @param[in]  surface_reflectivity        As the WSV (used with ground_type=S).
 * @param[in]  surface_complex_refr_index  As the WSV (used with ground_type=F).
 * @param[in]  stokes_dim                  As the WSV.
 *
 * @author Jana Mendrok
 * @date   2017-02-23
 */
void get_disortsurf_props(  // Output
    Vector& albedo,
    Numeric& btemp,
    // Input
    ConstVectorView f_grid,
    const Numeric& surface_skin_t,
    ConstVectorView surface_scalar_reflectivity);

/** Calculate doit_i_feild with Disort.
 *
 * Prepares actual input variables for Disort, runs it, and sorts the output
 * into doit_i_field.
 *
 * This version uses the C implementation of Disort based on ::run_disort.
 *
 * Altitudes, temperatures, VMRs and PNDs shall be provided with lat and lon
 * dimensions removed
 *
 * @param[in,out] ws Current workspace
 * @param[out]    doit_i_field Radiation field
 * @param[in]     f_grid Frequency grid
 * @param[in]     p_grid Pressure grid
 * @param[in]     z_profile Profile of geometric altitudes.
 * @param[in]     z_surface Surface altitude.
 * @param[in]     t_profile Temperature profile.
 * @param[in]     vmr_profiles VMR profiles.
 * @param[in]     pnd_profiles PND profiles.
 * @param[in]     scat_data Array of single scattering data
 * @param[in]     propmat_clearsky_agenda calculates the absorption coefficient
                  matrix
 * @param[in]     cloudbox_limits Cloudbox limits
 * @param[in]     surface_skin_t Surface skin temperature
 * @param[in]     surface_scalar_reflectivity Surface scalar reflectivity
 * @param[in]     scat_za_grid Zenith angle grid
 * @param[in]     nstreams Number of quadrature angles (both hemispheres).
 * @param[in]     Npfct Number of angular grid points to calculate bulk phase
 *                function
 * @param[in]     quiet Silence warnings
 * @param[in]     verbosity Verbosity setting
 *
 * @author        Oliver Lemke
 * @date          2019-09-19
 */
void run_cdisort(Workspace& ws,
                 // Output
                 Tensor7& doit_i_field,
                 // Input
                 ConstVectorView f_grid,
                 ConstVectorView p_grid,
                 ConstVectorView z_profile,
                 const Numeric& z_surface,
                 ConstVectorView t_profile,
                 ConstMatrixView vmr_profiles,
                 ConstMatrixView pnd_profiles,
                 const ArrayOfArrayOfSingleScatteringData& scat_data,
                 const Agenda& propmat_clearsky_agenda,
                 const ArrayOfIndex& cloudbox_limits,
                 const Numeric& surface_skin_t,
                 const Vector& surface_scalar_reflectivity,
                 ConstVectorView scat_za_grid,
                 const Index& nstreams,
                 const Index& Npfct,
                 const Index& quiet,
                 const Verbosity& verbosity);

/** get_gasoptprop.
 *
 * Derives level-based gas bulk optical properties (extinction).
 *
 * @param[in,out] ws                      Current workspace.
 * @param[out]    ext_bulk_gas            Gas bulk extinction (all levels & freqs).
 * @param[in]     propmat_clearsky_agenda As the WSV.
 * @param[in]     t_profile               Temperature profile
 * @param[in]     vmr_profiles            VMR profiles
 * @param[in]     p_grid                  As the WSV.
 * @param[in]     f_grid                  As the WSV.
 *
 * @author        Jana Mendrok
 * @date          2018-04-04
 */
void get_gasoptprop(Workspace& ws,
                    MatrixView ext_bulk_gas,
                    const Agenda& propmat_clearsky_agenda,
                    ConstVectorView t_profile,
                    ConstMatrixView vmr_profiles,
                    ConstVectorView p_grid,
                    ConstVectorView f_grid);

/** get_paroptprop.
 *
 * Derives level-based particle bulk optical properties (extinction and
 * absorption).
 *
 * @param[out] ext_bulk_par     Particle bulk extinction (all levels & freqs).
 * @param[out] abs_bulk_par     Particle bulk absorption (all levels & freqs).
 * @param[in]  scat_data        As the WSV.
 * @param[in]  pnd_profiles     PND profiles.
 * @param[in]  t_profile        Temperature profile
 * @param[in]  p_grid           As the WSV.
 * @param[in]  cloudbox_limits  As the WSV.
 * @param[in]  f_grid           As the WSV.
 *
 * @author     Jana Mendrok
 * @date       2018-04-04
 */
void get_paroptprop(MatrixView ext_bulk_par,
                    MatrixView abs_bulk_par,
                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                    ConstMatrixView pnd_profiles,
                    ConstVectorView t_profile,
                    ConstVectorView p_grid,
                    const ArrayOfIndex& cloudbox_limits,
                    ConstVectorView f_grid);

/** get_dtauc_ssalb
 *
 * Calculates layer averaged cloud optical depth (dtauc) and
 * single scattering albedo (ssalb) as required as DISORT subroutine input from
 * level-based gas extinction and particle extinction and absorption.
 *
 * @param[out] dtauc         Optical depths for all layers.
 * @param[out] ssalb         Single scattering albedos for all layers.
 * @param[in]  ext_bulk_gas  See get_gasoptprop.
 * @param[in]  ext_bulk_par  See get_paroptprop.
 * @param[in]  abs_bulk_par  See get_paroptprop.
 * @param[in]  z_profile     Profile of geometrical altitudes.
 *
 * @author     Jana Mendrok
 * @date       2018-04-04
 */
void get_dtauc_ssalb(MatrixView dtauc,
                     MatrixView ssalb,
                     ConstMatrixView ext_bulk_gas,
                     ConstMatrixView ext_bulk_par,
                     ConstMatrixView abs_bulk_par,
                     ConstVectorView z_profile);

/** get_angs.
 *
 * Derives angular grid to derive bulk phase matrix/function data on for further
 * Legendre decomposition.
 *
 * @param[out] pfct_angs  Angular grid of pfct_bulk_par.
 * @param[in]  scat_data  As the WSV.
 * @param[in]  nang       Number of angular grid points in pfct_angs. If<0,
 *
 * pfct_angs is taken from scat_data (the finest za_grid used over the scat elems),
 * else an equidistant grid with nang grid points is used.
 *
 * @author Jana Mendrok
 * @date   2018-04-04
 */
void get_angs(Vector& pfct_angs,
              const ArrayOfArrayOfSingleScatteringData& scat_data,
              const Index& Npfct);

/** get_parZ.
 *
 * Derives level-based particle bulk phase matrix Z (Csca scaled).
 * NOTE: Provided on ssd's freq grid (i.e. for nf=1 only of ssd.f_grid.nelem==1)
 * in order to avoid duplicate calculations in get_pmom (instead we duplicate the
 * results there to the RT calc's f_grid).
 *
 * @param[out] pha_bulk_par     Particle bulk phase function (all levels & ssd freqs).
 * @param[out] pfct_angs        Angular grid of pfct_bulk_par.
 * @param[in]  scat_data        As the WSV.
 * @param[in]  pnd_profiles     PND profiles.
 * @param[in]  t_profile        Temperature profile.
 * @param[in]  p_grid           As the WSV.
 * @param[in]  cloudbox_limits  As the WSV.
 * @param[in]  ext_bulk_par     See get_paroptprop.
 * @param[in]  abs_bulk_par     See get_paroptprop.
 * @param[in]  nang             Number of angular grid points in pfct_angs. If<0,
 *
 * pfct_angs is taken from scat_data (the finest za_grid used over the scat elems),
 * else an equidistant grid with nang grid points is used.
 *
 * @author   Jana Mendrok
 * @date     2018-04-04
 */
void get_parZ(Tensor3& pha_bulk_par,
              const ArrayOfArrayOfSingleScatteringData& scat_data,
              ConstMatrixView pnd_profiles,
              ConstVectorView t_profile,
              ConstVectorView pfct_angs,
              const ArrayOfIndex& cloudbox_limits);

/** get_pfct.
 *
 * Derives layer averaged particle bulk phase function P (4Pi scaled)
 * NOTE: Provided on ssd's freq grid (i.e. for nf=1 only if ssd.f_grid.nelem==1)
 * in order to avoid duplicate calculations in get_pmom (instead we duplicate the
 * results there to the RT calc's f_grid).
 *
 * @param[out] pfct_bulk_par  Particle bulk phase function (all levels & ssd freqs).
 * @param[in]  pha_bulk_par   See get_parZ.
 * @param[in]  ext_bulk_par   See get_paroptprop.
 * @param[in]  abs_bulk_par   See get_paroptprop.
 *
 * @author Jana Mendrok
 * @date   2018-04-04
 */
void get_pfct(Tensor3& pfct_bulk_par,
              ConstTensor3View& pha_bulk_par,
              ConstMatrixView ext_bulk_par,
              ConstMatrixView abs_bulk_par,
              const ArrayOfIndex& cloudbox_limits);

/** get_pmom
 *
 * Calculates Legendre moments of the layer averaged phase functionss (pmom) as
 * required as DISORT subroutine input from level-based bulk particle phase
 * function (4-Pi normalized scalar phase matrix).
 *
 * @param[out] pmom           Legendre moments for all layers.
 * @param[in]  pfct_bulk_par  See get_pfct.
 * @param[in]  pfct_angs      See get_parZ.
 * @param[in]  Nlegendre      Number of Legendre moments to derive.
 *
 * @author     Jana Mendrok
 * @date       2018-04-04
 */
void get_pmom(Tensor3View pmom,
              ConstTensor3View pfct_bulk_par,
              ConstVectorView pfct_angs,
              const Index& Nlegendre);

#endif /* disort_h */
