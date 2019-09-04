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

/**
 * @file   montecarlo.cc
 * @author Cory Davis <cory@met.ed.ac.uk>
 * @date   2003-06-19
 *
 * @brief  functions used by MCGeneral
 *
 *** FIXMEDOC ***: mcPathTraceRadar, opt_propExtract, Sample_los, Sample_los_uniform
 *                 cloudy_rt_vars_at_gp (JUST DATE CHECK), clear_rt_vars_at_gp(JUST
 *                 DATE CHECK).
 *** FIXMEDOC ***: In .cc some root-finding helper functions (for MCRadar) that don't
 *                 need visibility outside .cc, meaning that are not in .h
 *                 JUST CLARIFICATION => brent_zero, ext_I. Further docomentation ???
 *
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "check_input.h"
#include "cloudbox.h"
#include "lin_alg.h"
#include "logic.h"
#include "matpackI.h"
#include "messages.h"
#include "optproperties.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rng.h"
#include "rte.h"
#include "special_interp.h"
#include "xml_io.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;

/** clear_rt_vars_at_gp.
 *
 * Calculates a bunch of atmospheric variables at the end of a ppath.
 *
 * @param[in,out] ws                      Current workspace.
 * @param[out]    ext_mat_mono            Total monochromatic extinction matrix.
 * @param[out]    abs_vec_mono            Total monochromatic absorption vector.
 * @param[out]    temperature             a vector of temperatures
 * @param[in]     propmat_clearsky_agenda As the WSA.
 * @param[in]     f_mono                  Frequency (single entry vector).
 * @param[in]     gp_p                    An array of pressure gridpoints.
 * @param[in]     gp_lat                  An array of latitude gridpoints.
 * @param[in]     gp_lon                  An array of longitude gridpoints.
 * @param[in]     p_grid                  The pressure grid.
 * @param[in]     t_field                 The temperature grid.
 * @param[in]     vmr_field               VMR field.
 *
 * @author Cory Davis
 * @date 2005-02-19? *** FIXMEDOC ***
 */
void clear_rt_vars_at_gp(Workspace& ws,
                         MatrixView ext_mat_mono,
                         VectorView abs_vec_mono,
                         Numeric& temperature,
                         const Agenda& propmat_clearsky_agenda,
                         const Numeric& f_mono,
                         const GridPos& gp_p,
                         const GridPos& gp_lat,
                         const GridPos& gp_lon,
                         ConstVectorView p_grid,
                         ConstTensor3View t_field,
                         ConstTensor4View vmr_field);

/** cloudy_rt_vars_at_gp.
 *
 * Calculates a bunch of atmospheric variables at the end of a ppath.
 *
 * @param[in,out] ws                      Current workspace.
 * @param[out]    ext_mat_mono            Total monochromatic extinction matrix.
 * @param[out]    abs_vec_mono            Total monochromatic absorption vector.
 * @param[out]    pnd_vec                 Vector of particle number densities (one element per scattering element).
 * @param[out]    temperature             A vector of temperatures
 * @param[in]     propmat_clearsky_agenda Agenda calculating the absorption coefficient matrices.
 * @param[in]     stokes_dim              The dimensionality of the Stokes vector (1-4).
 * @param[in]     f_index                 Index of frequency grid point handeled.
 * @param[in]     f_grid                  Frequency grid for monochromatic pencil beam calculations.
 * @param[in]     gp_p                    An array of pressure gridpoints.
 * @param[in]     gp_lat                  An array of latitude gridpoints.
 * @param[in]     gp_lon                  An array of longitude gridpoints.
 * @param[in]     p_grid_cloud            The subset of p_grid within the cloudbox.
 * @param[in]     t_field_cloud           The subset of t_field within the cloudbox.
 * @param[in]     vmr_field_cloud         The subset of vmr_field within the cloudbox.
 * @param[in]     pnd_field               Particle number density field.
 * @param[in]     scat_data               Array of single scattering data.
 * @param[in]     cloudbox_limits         The limits of the cloud box.
 * @param[in]     rte_los
 *
 * @author Cory Davis
 * @date 2005-02-19? *** FIXMEDOC ***
 */
void cloudy_rt_vars_at_gp(Workspace& ws,
                          MatrixView ext_mat_mono,
                          VectorView abs_vec_mono,
                          VectorView pnd_vec,
                          Numeric& temperature,
                          const Agenda& propmat_clearsky_agenda,
                          const Index stokes_dim,
                          const Index f_index,
                          const Vector& f_grid,
                          const GridPos& gp_p,
                          const GridPos& gp_lat,
                          const GridPos& gp_lon,
                          ConstVectorView p_grid_cloud,
                          ConstTensor3View t_field_cloud,
                          ConstTensor4View vmr_field_cloud,
                          const Tensor4& pnd_field,
                          const ArrayOfArrayOfSingleScatteringData& scat_data,
                          const ArrayOfIndex& cloudbox_limits,
                          const Vector& rte_los);

/** cloud_atm_vars_by_gp.
 *
 *
 * Returns pressure, temperature, VMRs and PNDs, at points corresponding
 * to arrays of gridpositions gp_p, gp_lat, and gp_lon.  The field and grid
 * input variables all span only the cloudbox
 *
 * @param[out] pressure         A vector of pressures.
 * @param[out] temperature      A vector of temperatures.
 * @param[out] vmr              A n_species by n_p matrix of VMRs.
 * @param[out] pnd              A n_scatelem by n_p matrix of PNDs.
 * @param[in]  gp_p             An array of pressure gridpoints.
 * @param[in]  gp_lat           An array of latitude gridpoints.
 * @param[in]  gp_lon           An array of longitude gridpoints.
 * @param[in]  cloudbox_limits  The limits of the cloud box.
 * @param[in]  p_grid_cloud     The subset of p_grid within the cloudbox.
 * @param[in]  t_field_cloud    The subset of t_field within the cloudbox.
 * @param[in]  vmr_field_cloud  The subset of vmr_field within the cloudbox.
 * @param[in]  pnd_field        Particle number density field.
 *
 * @author Cory Davis
 * @date 2005-06-07
 */
void cloud_atm_vars_by_gp(VectorView pressure,
                          VectorView temperature,
                          MatrixView vmr,
                          MatrixView pnd,
                          const ArrayOfGridPos& gp_p,
                          const ArrayOfGridPos& gp_lat,
                          const ArrayOfGridPos& gp_lon,
                          const ArrayOfIndex& cloudbox_limits,
                          ConstVectorView p_grid_cloud,
                          ConstTensor3View t_field_cloud,
                          ConstTensor4View vmr_field_cloud,
                          ConstTensor4View pnd_field);

/** get_ppath_transmat.
 *
 * Routine to get the transmission matrix along a pre-defined propagation path.
 * This is based on mcPathTraceGeneral using the routines from this source file.
 * Routines from rte.cc require wind and magnetic field data that has not been
 * typically passed to the Monte Carlo routines.
 *
 * @param[in,out] ws                      Current workspace.
 * @param[out]    trans_mat               Matrix defining transmission over the ppath
 *                                        direction multiplied by sin(za)
 * @param[in]     ppath                   Propagation path over which transmission matrix is desired
 * @param[in]     propmat_clearsky_agenda Agenda calculating the absorption coefficient matrices.
 * @param[in]     stokes_dim              The dimensionality of the Stokes vector (1-4).
 * @param[in]     f_index                 Index of frequency grid point handeled.
 * @param[in]     f_grid                  Frequency grid for monochromatic pencil beam calculations.
 * @param[in]     p_grid                  The pressure grid.
 * @param[in]     t_field                 The temperature grid.
 * @param[in]     vmr_field               VMR field.
 * @param[in]     cloudbox_limits         The limits of the cloud box.
 * @param[in]     pnd_field               Particle number density field.
 * @param[in]     scat_data               Array of single scattering data.
 * @param[in]     verbosity               Verbosity variable to dynamically control the reporting
 *                                        level during runtime.
 *
 *
 * @author        Ian S. Adams
 * @date          2015-09-15
 */
void get_ppath_transmat(
    Workspace& ws,
    MatrixView& trans_mat,
    const Ppath& ppath,
    const Agenda& propmat_clearsky_agenda,
    const Index stokes_dim,
    const Index f_index,
    const Vector& f_grid,
    const Vector& p_grid,
    const Tensor3& t_field,
    const Tensor4& vmr_field,
    const ArrayOfIndex& cloudbox_limits,
    const Tensor4& pnd_field,
    const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
    const Verbosity& verbosity);

/** is_anyptype_nonTotRan.
 *
 * Some operations in Monte Carlo simulations are different depending on the
 * ptype of the scattering elements. This function searches scat_data
 * to determine if any of the scattering elements have ptype=30.
 *
 * @param[in]  scat_data_mono  Monochromatic single scattering data.
 *
 * @author Cory Davis
 * @date 2004-1-31
 */
bool is_anyptype_nonTotRan(
    const ArrayOfArrayOfSingleScatteringData& scat_data_mono);

/** mcPathTraceGeneral.
 *
 * Performs the tasks of pathlength sampling.
 *
 * Ray tracing done (but now only as far as determined by pathlength
 * sampling) and calculation of the evolution operator and several
 * atmospheric variables at the new point.
 *
 * The end point of the ray tracing is returned by ppath_step, where the
 * point of concern has index ppath_step.np-1. However, a somehwat dirty trick
 * is used here to avoid copying of data. Only ppath.np is adjusted, and
 * ppath_step can contain additional points (that should not be used).
 *
 * @param[in,out] ws                      Current workspace.
 * @param[out]    evol_op                 Evolution operator (Stokes attenuation operator; exp(-tau)).
 * @param[out]    abs_vec_mono            Total monochromatic absorption vector.
 * @param[out]    temperature             A vector of temperatures.
 * @param[out]    ext_mat_mono            Total monochromatic extinction matrix.
 * @param[out]    rng                     Random number generator instance.
 * @param[out]    rte_pos                 Position for starting radiative transfer simulations.
 * @param[out]    rte_los                 Incident line of sight for subsequent ray-tracing.
 * @param[out]    g                       randomly chosen extinction path.
 * @param[out]    ppath_step              A propagation path step.
 * @param[out]    termination_flag        Flag defining whether the path of the photon is terminated
 * @param[in]     inside_cloud            Flag defining inside or not the cloud box
 * @param[in]     ppath_step_agenda       Calculation of a propagation path step.
 * @param[in]     ppath_lmax              Maximum length between points describing propagation paths.
 * @param[in]     ppath_lraytrace         Maximum length of ray tracing steps when determining
 *                                        propagation paths.
 * @param[in]     taustep_limit           Defines an upper step length in terms of optical thickness
 *                                        for Monte Carlo calculations.
 * @param[in]     propmat_clearsky_agenda Agenda calculating the absorption coefficient matrices.
 * @param[in]     stokes_dim              The dimensionality of the Stokes vector (1-4).
 * @param[in]     f_index                 Index of frequency grid point handeled.
 * @param[in]     f_grid                  Frequency grid for monochromatic pencil beam calculations.
 * @param[in]     p_grid                  Pressure grid.
 * @param[in]     lat_grid                Latitude grid.
 * @param[in]     lon_grid                Longitude grid.
 * @param[in]     z_field                 The field of geometrical altitudes.
 * @param[in]     refellipsoid            Reference ellipsoid.
 * @param[in]     z_surface               The surface altitude.
 * @param[in]     t_field                 The temperature grid.
 * @param[in]     vmr_field               VMR field.
 * @param[in]     cloudbox_limits         The limits of the cloud box.
 * @param[in]     pnd_field               Particle number density field.
 * @param[in]     scat_data               Array of single scattering data.
 * @param[in]     verbosity               Verbosity variable to dynamically control the reporting
 *                                        level during runtime.
 *
 *
 * @author        Cory Davis
 * @date          2005-2-21
 */
void mcPathTraceGeneral(Workspace& ws,
                        MatrixView evol_op,
                        Vector& abs_vec_mono,
                        Numeric& temperature,
                        MatrixView ext_mat_mono,
                        Rng& rng,
                        Vector& rte_pos,
                        Vector& rte_los,
                        Vector& pnd_vec,
                        Numeric& g,
                        Ppath& ppath_step,
                        Index& termination_flag,
                        bool& inside_cloud,
                        const Agenda& ppath_step_agenda,
                        const Numeric& ppath_lmax,
                        const Numeric& ppath_lraytrace,
                        const Numeric& taustep_limit,
                        const Agenda& propmat_clearsky_agenda,
                        const Index stokes_dim,
                        const Index f_index,
                        const Vector& f_grid,
                        const Vector& p_grid,
                        const Vector& lat_grid,
                        const Vector& lon_grid,
                        const Tensor3& z_field,
                        const Vector& refellipsoid,
                        const Matrix& z_surface,
                        const Tensor3& t_field,
                        const Tensor4& vmr_field,
                        const ArrayOfIndex& cloudbox_limits,
                        const Tensor4& pnd_field,
                        const ArrayOfArrayOfSingleScatteringData& scat_data,
                        const Verbosity& verbosity);

/** mcPathTraceRadar.
 *
 * Performs the tasks of pathlength sampling.
 *
 * Ray tracing done (but now only as far as determined by pathlength
 * sampling) and calculation of the evolution operator and several
 * atmospheric variables at the new point.
 *
 * The end point of the ray tracing is returned by ppath_step, where the
 * point of concern has index ppath_step.np-1. However, a somehwat dirty trick
 * is used here to avoid copying of data. Only ppath.np is adjusted, and
 * ppath_step can contain additional points (that should not be used).
 *
 * @param[in,out] ws                      Current workspace.
 * @param[out]    evol_op                 Evolution operator (Stokes attenuation operator; exp(-tau)).
 * @param[out]    abs_vec_mono            Total monochromatic absorption vector.
 * @param[out]    temperature             A vector of temperatures.
 * @param[out]    ext_mat_mono            Total monochromatic extinction matrix.
 * @param[out]    rng                     Random number generator instance.
 * @param[out]    rte_pos                 Position for starting radiative transfer simulations.
 * @param[out]    rte_los                 Incident line of sight for subsequent ray-tracing.
 * @param[out]    stot                    Starting point of total path lenght
 * @param[out]    ttot                    Total path length *** FIXMEDOC ***
 * @param[out]    ppath_step              A propagation path step.
 * @param[out]    termination_flag        Flag defining whether the path of the photon is terminated
 * @param[in]     inside_cloud            Flag defining inside or not the cloud box
 * @param[in]     ppath_step_agenda       Calculation of a propagation path step.
 * @param[in]     ppath_lmax              Maximum length between points describing propagation paths.
 * @param[in]     ppath_lraytrace         Maximum length of ray tracing steps when determining
 *                                        propagation paths.
 * @param[in]     taustep_limit           Defines an upper step length in terms of optical thickness
 *                                        for Monte Carlo calculations.
 * @param[in]     propmat_clearsky_agenda Agenda calculating the absorption coefficient matrices.
 * @param[in]     anyptype_nonTotRan      Flag definining any particle type, but for totally random oriented.
 * @param[in]     stokes_dim              The dimensionality of the Stokes vector (1-4).
 * @param[in]     f_index                 Index of frequency grid point handeled.
 * @param[in]     f_grid                  Frequency grid for monochromatic pencil beam calculations.
 * @param[in]     Iprop                   Incident I component of the Stokes vector *** FIXMEDOC ***
 * @param[in]     p_grid                  Pressure grid.
 * @param[in]     lat_grid                Latitude grid.
 * @param[in]     lon_grid                Longitude grid.
 * @param[in]     z_field                 The field of geometrical altitudes.
 * @param[in]     refellipsoid            Reference ellipsoid.
 * @param[in]     z_surface               The surface altitude.
 * @param[in]     t_field                 The temperature grid.
 * @param[in]     vmr_field               VMR field.
 * @param[in]     cloudbox_limits         The limits of the cloud box.
 * @param[in]     pnd_field               Particle number density field.
 * @param[in]     scat_data               Array of single scattering data.
 * @param[in]     verbosity               Verbosity variable to dynamically control the reporting
 *                                        level during runtime.
 *
 * @author        Cory Davis (mcPathTraceGeneral), Ian S. Adams
 * @date          2015-09-08
 */
void mcPathTraceRadar(Workspace& ws,
                      MatrixView evol_op,
                      Vector& abs_vec_mono,
                      Numeric& temperature,
                      MatrixView ext_mat_mono,
                      Rng& rng,
                      Vector& rte_pos,
                      Vector& rte_los,
                      Vector& pnd_vec,
                      Numeric& stot,
                      Numeric& ttot,
                      Ppath& ppath_step,
                      Index& termination_flag,
                      bool& inside_cloud,
                      const Agenda& ppath_step_agenda,
                      const Numeric& ppath_lmax,
                      const Numeric& ppath_lraytrace,
                      const Agenda& propmat_clearsky_agenda,
                      const bool& anyptype_nonTotRan,
                      const Index stokes_dim,
                      const Index f_index,
                      const Vector& f_grid,
                      const Vector& Iprop,
                      const Vector& p_grid,
                      const Vector& lat_grid,
                      const Vector& lon_grid,
                      const Tensor3& z_field,
                      const Vector& refellipsoid,
                      const Matrix& z_surface,
                      const Tensor3& t_field,
                      const Tensor4& vmr_field,
                      const ArrayOfIndex& cloudbox_limits,
                      const Tensor4& pnd_field,
                      const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                      const Verbosity& verbosity);
/* opt_propCalc.
 *
 * Returns the total monochromatic extinction matrix and absorption vector over all
 * scattering elements at a specific atmospheric location.
 *
 * @return    ext_mat_mono    Total monochromatic extinction matrix.
 * @return    abs_vec_mono    Total monochromatic absorption vector.
 * @param[in] za              Zenith angle of propagation direction.
 * @param[in] aa              Azimuthal angle of propagation.
 * @param[in] scat_data_mono  As the WSV.
 * @param[in] stokes_dim      As the WSV.
 * @param[in] pnd_vec         Vector of particle number densities (one element per
 *                            scattering element).
 * @param[in] rtp_temperature As the WSV
 *
 * @author    Cory Davis
 * @date      2004-7-16
 */
void opt_propCalc(MatrixView K,
                  VectorView K_abs,
                  const Numeric za,
                  const Numeric aa,
                  const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                  const Index stokes_dim,
                  ConstVectorView pnd_vec,
                  const Numeric rtp_temperature,
                  const Verbosity& verbosity);

/* opt_propExtract.
 *
 * Extracts the total monochromatic extinction matrix and absorption vector over all
 * scattering elements at a specific atmospheric location.
 *
 * @return    ext_mat_mono_spt *** FIXMEDOC ***    Total monochromatic extinction matrix.
 * @return    abs_vec_mono_spt *** FIXMEDOC ***    Total monochromatic absorption vector.
 * @param[in]  scat_data_single A monochromatic SingleScatteringData object.
 * @param[in] za               Zenith angle of propagation direction.
 * @param[in] aa (_U_,)        *** FIXMEDOC ***       Azimuthal angle of propagation.
 * @param[in] rtp_temperature  As the WSV
 * @param[in] stokes_dim       As the WSV.
 * @param[in] verbosity        Verbosity variable to dynamically control the reporting
 *                             level during runtime.
 *
 * @author    *** FIXMEDOC ***
 * @date      *** FIXMEDOC ***
 */
void opt_propExtract(MatrixView ext_mat_mono_spt,
                     VectorView abs_vec_mono_spt,
                     const SingleScatteringData& scat_data_single,
                     const Numeric za,
                     const Numeric aa,
                     const Numeric rtp_temperature,
                     const Index stokes_dim,
                     const Verbosity& verbosity);

/** pha_mat_singleCalc.
 *
 * Returns the total phase matrix for given incident and scattered directions.
 * It requires a vector of particle number densities to be precalculated.
 *
 * @param[out] Z               Phase matrix.
 * @param[in]  za_sca          Zenith angle of scattering direction.
 * @param[in]  aa_sca          Azimuth angle of scattering direction.
 * @param[in]  za_inc          Zenith angle of incident direction.
 * @param[in]  aa_inc          Azimuth angle of incident direction.
 * @param[in]  scat_data_mono  As the WSV.
 * @param[in]  stokes_dim      As the WSV.
 * @param[in]  pnd_vec         Vector of particle number densities at the point
 *                             in question.
 * @param[in]  rtp_temperature As the WSV.
 *
 * @author     Cory Davis
 * @date       2003-11-27
 */
void pha_mat_singleCalc(
    MatrixView Z,
    const Numeric za_sca,
    const Numeric aa_sca,
    const Numeric za_inc,
    const Numeric aa_inc,
    const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
    const Index stokes_dim,
    ConstVectorView pnd_vec,
    const Numeric rtp_temperature,
    const Verbosity& verbosity);

/** pha_mat_singleCalcScatElement.
 *
 * Returns the phase matrix for each scattering element, given incident and
 * scattered directions. It requires a vector of particle number densities
 * to be precalculated.
 *
 * @param[out] Z               Phase matrix.
 * @param[in]  za_sca          Zenith angle of scattering direction.
 * @param[in]  aa_sca          Azimuth angle of scattering direction.
 * @param[in]  za_inc          Zenith angle of incident direction.
 * @param[in]  aa_inc          Azimuth angle of incident direction.
 * @param[in]  scat_data_mono  As the WSV.
 * @param[in]  stokes_dim      As the WSV.
 * @param[in]  pnd_vec         Vector of particle number densities at the point
 *                             in question.
 * @param[in]  rtp_temperature As the WSV.
 *
 * @author     Cory Davis
 * @date       2003-11-27
 */
void pha_mat_singleCalcScatElement(
    Tensor3View Z,
    const Numeric za_sca,
    const Numeric aa_sca,
    const Numeric za_inc,
    const Numeric aa_inc,
    const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
    const Index stokes_dim,
    ConstVectorView pnd_vec,
    const Numeric rtp_temperature,
    const Verbosity& verbosity);

/** pha_mat_singleExtract.
 *
 * Extract the phase matrix from a monochromatic SingleScatteringData object
 *
 * Given a monochromatic SingleScatteringData object, incident and scattered
 * directions, and the temperature, this function returns the phase matrix in
 * the laboratory frame.
 *
 * @param[out] Z_spt            The phase matrix.
 * @param[in]  scat_data_single A monochromatic SingleScatteringData object.
 * @param[in]  za_sca           Zenith angle of scattering direction.
 * @param[in]  aa_sca           Azimuth angle of scattering direction.
 * @param[in]  za_inc           Zenith angle of incident direction.
 * @param[in]  aa_inc           Azimuth angle of incident direction.
 * @param[in]  rtp_temperature  As the WSV.
 * @param[in]  stokes_dim       As the WSV.
 *
 * @author     Cory Davis
 * @date       2004-07-16
 */
void pha_mat_singleExtract(MatrixView Z_spt,
                           const SingleScatteringData& scat_data_single,
                           const Numeric za_sca,
                           const Numeric aa_sca,
                           const Numeric za_inc,
                           const Numeric aa_inc,
                           const Numeric rtp_temperature,
                           const Index stokes_dim,
                           const Verbosity& verbosity);

/** Sample_los.
 *
 * *** FIXMEDOC ***: 2011-06-17 Documentation removed by Gerrit (severely out of date)
 * *** FIXMEDOC ***: Sampling the new direction???
 *
 *
 * @param[out]    new_rte_los      Incident line of sight for subsequent.
 * @param[out]    g_los_csc_theta  Probability density for the chosen
 *                                 direction multiplied by sin(za)
 * @param[out]    Z                Bulk phase matrix in Stokes notation.
 * @param[in,out] rng              Rng random number generator instance.
 * @param[in]     rte_los          Incident line of sight for subsequent
 *                                 ray-tracing.
 * @param[in]     scat_data        As the WSV.
 * @param[in]     stokes_dim       As the WSV.
 * @param[out]    pnd_vec          Vector of particle number densities (one element per scattering element).
 * @param[in]     Z11maxvector     *** FIXMEDOC *** Strong Forward peak???
 * @param[in]     Csca             *** FIXMEDOC *** Scattering cross section???
 * @param[in]     rtp_temperature  As the WSV.
 *
 * @author Cory Davis
 * @date   2003-06-19
 */
void Sample_los(VectorView new_rte_los,
                Numeric& g_los_csc_theta,
                MatrixView Z,
                Rng& rng,
                ConstVectorView rte_los,
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const Index stokes_dim,
                const Index f_index,
                ConstVectorView pnd_vec,
                ConstVectorView Z11maxvector,
                const Numeric Csca,
                const Numeric rtp_temperature,
                const Index t_interp_order = 1);

/** Sample_los_uniform.
 *
 * *** FIXMEDOC ***
 * *** FIXMEDOC ***: Sampling the new direction uniformly???
 *
 * @param[out] new_rte_los Incident line of sight for subsequent.
 * @param[out] rng         Rng random number generator instance.
 *
 * @author     *** FIXMEDOC ***
 * @date       *** FIXMEDOC ***
 */
void Sample_los_uniform(VectorView new_rte_los, Rng& rng);

#endif  // montecarlo_h
