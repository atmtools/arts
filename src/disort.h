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
 *** FIXMEDOC ***: get_cb_inc_field (x2);get_disortsurf_props;run_disort2
 *** FIXMEDOC ***: functions usually x2; probably version2 should be employed?
 *** FIXMEDOC ***: within disort.cc plenty of text is commented out; delete?
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
 * @param[in]  pnd_ncols             Number of columns (latitude points) in *pnd_field*.
 * @param[in]  ifield_npages         Number of pages (polar angle points) in *doit_i_field*.
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
    const String& pfct_method,
    const Index& pnd_ncols);

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

/** run_disort.
 *
 * Prepares actual input variables for Disort, runs it, and sorts the output into
 * doit_i_field.
 *
 * @param[in,out] ws                          Current workspace.
 * @param[out]    doit_i_field                As the WSV.
 * @param[in]     f_grid                      As the WSV.
 * @param[in]     p_grid                      As the WSV.
 * @param[in]     z_field                     As the WSV.
 * @param[in]     t_field                     As the WSV.
 * @param[in]     vmr_field                   As the WSV.
 * @param[in]     pnd_field                   As the WSV.
 * @param[in]     scat_data                   As the WSV.
 * @param[in]     propmat_clearsky_agenda     As the WSA.
 * @param[in]     cloudbox_limits             As the WSV.
 * @param[in]     surface_skin_t              As the WSV.
 * @param[in]     surface_scalar_reflectivity As the WSM.
 * @param[in]     scat_za_grid                As the WSV.
 * @param[in]     nstreams                    Number of quadrature angles (both hemispheres).
 * @param[in]     do_deltam                   See DisortCalc doc.
 * @param[in]     non_iso_inc                 See DisortCalc doc.
 * @param[in]     pfct_method                 See DisortCalc doc.
 *
 * @author        Jana Mendrok
 * @date          2017-02-23
 */
void run_disort(Workspace& ws,
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
                const Verbosity& verbosity);

/** run_disort2 *** FIXMEDOC ***
 *
 * Prepares actual input variables for Disort, runs it, and sorts the output into
 * doit_i_field.
 * This version using unified optprop extraction scheme. supposed to replace
 * run_disort
 *
 * @param[in,out] ws                          Current workspace.
 * @param[out]    doit_i_field                As the WSV.
 * @param[in]     f_grid                      As the WSV.
 * @param[in]     p_grid                      As the WSV.
 * @param[in]     z_field                     As the WSV.
 * @param[in]     t_field                     As the WSV.
 * @param[in]     vmr_field                   As the WSV.
 * @param[in]     pnd_field                   As the WSV.
 * @param[in]     scat_data                   As the WSV.
 * @param[in]     propmat_clearsky_agenda     As the WSA.
 * @param[in]     cloudbox_limits             As the WSV.
 * @param[in]     surface_skin_t              As the WSV.
 * @param[in]     surface_scalar_reflectivity As the WSM.
 * @param[in]     scat_za_grid                As the WSV.
 * @param[in]     nstreams                    Number of quadrature angles (both hemispheres).
 * @param[in]     do_deltam                   See DisortCalc doc.
 * @param[in]     non_iso_inc                 See DisortCalc doc.
 *
 * @author        Jana Mendrok
 * @date          2017-02-23
 */
void run_disort2(Workspace& ws,
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
                 const Verbosity& verbosity);

/** dtauc_ssalbCalc.
 *
 * Calculates layer averaged cloud optical depth (dtauc) and single
 * scattering albedo (ssalb). These variables are required as input
 * for the DISORT subroutine.
 *
 * @param[in,out] ws                  Current workspace.
 * @param[out] dtauc                  Optical depths for all layers.
 * @param[out] ssalb                  Single scattering albedos for all layers.
 * @param[in] scat_data               As the WSV.
 * @param[in] f_index                 Index of frequency grid point handeled.
 * @param[in] propmat_clearsky_agenda As the WSA.
 * @param[in] pnd_field               As the WSV.
 * @param[in] t_field                 As the WSV.
 * @param[in] z_field                 As the WSV.
 * @param[in] vmr_field               As the WSV.
 * @param[in] p_grid                  As the WSV.
 * @param[in] cloudbox_limits         As the WSV.
 * @param[in] f_mono                  Frequency (single entry vector).
 *
 * @author Claudia Emde, Jana Mendrok
 * @date   2006-02-10
 */
void dtauc_ssalbCalc(Workspace& ws,
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
                     const Verbosity& verbosity);

/** get_gasoptprop.
 *
 * Derives level-based gas bulk optical properties (extinction).
 *
 * @param[in,out] ws                      Current workspace.
 * @param[out]    ext_bulk_gas            Gas bulk extinction (all levels & freqs).
 * @param[in]     propmat_clearsky_agenda As the WSA.
 * @param[in]     t_field                 As the WSV.
 * @param[in]     vmr_field               As the WSV.
 * @param[in]     p_grid                  As the WSV.
 * @param[in]     f_grid                  As the WSV.
 *
 * @author        Jana Mendrok
 * @date          2018-04-04
 */
void get_gasoptprop(Workspace& ws,
                    MatrixView ext_bulk_gas,
                    const Agenda& propmat_clearsky_agenda,
                    ConstVectorView t_field,
                    ConstMatrixView vmr_field,
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
 * @param[in]  pnd_field        As the WSV.
 * @param[in]  t_field          As the WSV.
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
                    ConstMatrixView pnd_field,
                    ConstVectorView t_field,
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
 * @param[in]  z_field       Ss the WSV.
 *
 * @author     Jana Mendrok
 * @date       2018-04-04
 */
void get_dtauc_ssalb(MatrixView dtauc,
                     MatrixView ssalb,
                     ConstMatrixView ext_bulk_gas,
                     ConstMatrixView ext_bulk_par,
                     ConstMatrixView abs_bulk_par,
                     ConstVectorView z_field);

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
 * @param[in]  pnd_field        As the WSV.
 * @param[in]  t_field          As the WSV.
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
              ConstMatrixView pnd_field,
              ConstVectorView t_field,
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

/** dtauc_ssalbCalc2.
 *
 * Calculates layer averaged cloud optical depth (dtauc) and
 * single scattering albedo (ssalb). These variables are required as
 * input for the DISORT subroutine
 *
 * @param[in,out] ws                       Current workspace.
 * @param[out]    dtauc                    Optical depths for all layers.
 * @param[out]    ssalb                    Single scattering albedos for all layers.
 * @param[in]     scat_data                As the WSV.
 * @param[in]     propmat_clearsky_agenda  As the WSA.
 * @param[in]     pnd_field                As the WSV.
 * @param[in]     t_field                  As the WSV.
 * @param[in]     z_field                  As the WSV.
 * @param[in]     vmr_field                As the WSV.
 * @param[in]     p_grid                   As the WSV.
 * @param[in]     cloudbox_limits          As the WSV.
 * @param[in]     f_grid                   As the WSV.
 *
 * @author Jana Mendrok
 * @date   2018-02-18
 */
void dtauc_ssalbCalc2(Workspace& ws,
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
                      ConstVectorView f_grid);

/** phase_functionCalc2
 *
 * Calculates layer averaged normalized phase functions from the phase matrix
 * in SingleScatteringData. Temperature and angle grid interpolations are applied.
 *
 * @param[out] phase_function     Normalized layer-averaged bulk phase function.
 * @param[in]  scat_data          As the WSV.
 * @param[in]  f_index            Index of frequency grid point handeled.
 * @param[in]  pnd_field          As the WSV.
 * @param[in]  t_field            As the WSV.
 * @param[in]  cloudbox_limits    As the WSV.
 * @param[in]  pfct_za_grid_size  Number of equidistant scatt. angles in 0-180deg.
 *
 * @author Claudia Emde, Jana Mendrok
 * @date   2006-02-10
 */
void phase_functionCalc2(  //Output
    MatrixView phase_function,
    //Input
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Index& f_index,
    ConstTensor4View pnd_field,
    ConstTensor3View t_field,
    const ArrayOfIndex& cloudbox_limits,
    const Index& pfct_za_grid_size,
    const Verbosity& verbosity);

/** phase_functionCalc.
 *
 * Calculates layer averaged normalized phase functions from the phase matrix
 * in SingleScatteringData. The scattering angle grid is taken from the data.
 * It is required that all scattering elements are given on the same scattering
 * angle grid. No temperature interpolation done.
 *
 * @param[out] phase_function   Normalized phase function.
 * @param[in]  scat_data        As the WSV.
 * @param[in]  f_index          Index of frequency grid point handeled.
 * @param[in]  pnd_field        As the WSV.
 * @param[in]  cloudbox_limits  As the WSV.
 *
 * @author     Claudia Emde, Jana Mendrok
 * @date       2006-02-10
 */
void phase_functionCalc(  //Output
    MatrixView phase_function,
    //Input
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Index& f_index,
    ConstTensor4View pnd_field,
    const ArrayOfIndex& cloudbox_limits,
    const String pfct_method);

/** pmomCalc2.
 *
 * Calculates Legendre polynomials of phase functions for each layer.
 * The Legendre polynomial are required as input for DISORT.
 *
 * @param[out] pmom            Legendre polynomial of phase functions
 * @param[in]  phase_function  Normalized phase function
 * @param[in]  scat_angle_grid Scattering angle grid corresponding to phase
 functions
 * @param[in]  Nlegendre       Number of Legendre polynomials to be calculated
 *
 * @author Claudia Emde, Jana Mendrok
 * @date   2006-02-10
 */
void pmomCalc2(  //Output
    MatrixView pmom,
    //Input
    ConstMatrixView phase_function,
    ConstVectorView scat_angle_grid,
    const Index n_legendre,
    const Verbosity& verbosity);

/** pmomCalc.
 *
 * Calculates Legendre polynomials of phase functions for each layer.
 * The Legendre polynomial are required as input for DISORT.
 *
 * @param[out] pmom            Legendre polynomial of phase functions
 * @param[in]  phase_function  Normalized phase function
 * @param[in]  scat_angle_grid Scattering angle grid corresponding to phase
 functions
 * @param[in]  Nlegendre       Number of Legendre polynomials to be calculated
 *
 * @author Claudia Emde, Jana Mendrok
 * @date   2006-02-10
 */
void pmomCalc(  //Output
    MatrixView pmom,
    //Input
    ConstMatrixView phase_function,
    ConstVectorView scat_angle_grid,
    const Index n_legendre,
    const Verbosity& verbosity);

/** surf_albedoCalc
 *
 * Calculates the diffuse power reflection coefficient (an estimate of the total
 * surface albedo, equivalent to ARTS' surface_scalar_reflectivity) from
 * reflection matrices according to *surface_rt_prop_agenda* settings for use as
 * input parameter albedo to a Disort calculation (internally applying a
 * Lambertian surface).
 *
 * @param[in,out] ws                     Current workspace.
 * @param[out]    albedo                 Diffuse power reflection coefficient.
 * @param[out]    btemp                  Surface temperature
 * @param[in]     surface_rtprop_agenda  As the WSA.
 * @param[in]     f_grid                 As the WSV.
 * @param[in]     scat_za_grid           As the WSV.
 * @param[in]     surf_alt               Surface altitude.
 *
 * @author        Jana Mendrok
 * @date          2017-02-16
 */
void surf_albedoCalc(Workspace& ws,
                     //Output
                     VectorView albedo,
                     Numeric& btemp,
                     //Input
                     const Agenda& surface_rtprop_agenda,
                     ConstVectorView f_grid,
                     ConstVectorView scat_za_grid,
                     const Numeric& surf_alt,
                     const Verbosity& verbosity);

#ifdef ENABLE_DISORT
/** get_cb_inc_field
 *
 *** FIXMEDOC ***
 *
 *
 * @param[in,out] ws               Current workspace.
 * @param[out]    cb_inc_field     *** FIXMEDOC ***
 * @param[in]     iy_main_agenda   As the WSA.
 * @param[in]     z_field          As the WSV.
 * @param[in]     t_field          As the WSV.
 * @param[in]     vmr_field        As the WSV.
 * @param[in]     nlte_field       *** FIXMEDOC ***
 * @param[in]     cloudbox_limits  As the WSV.
 * @param[in]     f_grid           As the WSV.
 * @param[in]     scat_za_grid     As the WSV.
 * @param[in]     nstreams         Number of quadrature angles (both hemispheres).
 *
 * @author        Jana Mendrok
 * @date          2017-02-16
 */
void get_cb_inc_field(Workspace& ws,
                      Matrix& cb_inc_field,
                      const Agenda& iy_main_agenda,
                      const Tensor3& z_field,
                      const Tensor3& t_field,
                      const Tensor4& vmr_field,
                      const Tensor4& nlte_field,
                      const ArrayOfIndex& cloudbox_limits,
                      const Vector& f_grid,
                      const Vector& scat_za_grid,
                      const Index& nstreams);

#else /* ENABLE_DISORT */

/** get_cb_inc_field
 *
 *** FIXMEDOC ***
 *
 *
 * @param[in,out] ws               Current workspace.
 * @param[out]                      *** FIXMEDOC ***
 * @param[in]                       *** FIXMEDOC ***
 * @param[in]                       *** FIXMEDOC ***
 * @param[in]                       *** FIXMEDOC ***
 * @param[in]                       *** FIXMEDOC ***
 * @param[in]                       *** FIXMEDOC ***
 * @param[in]                       *** FIXMEDOC ***
 * @param[in]                       *** FIXMEDOC ***
 * @param[in]                       *** FIXMEDOC ***
 * @param[in]                       *** FIXMEDOC ***
 *
 * @author        Jana Mendrok
 * @date          2017-02-16
 */
void get_cb_inc_field(Workspace&,
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
                      const Index&);

#endif /* ENABLE_DISORT */

#endif /* disort_h */
