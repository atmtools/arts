#ifndef montecarlo_h
#define montecarlo_h

/**
 * @file   montecarlo.cc
 * @author Cory Davis <cory@met.ed.ac.uk>
 * @date   2003-06-19
 *
 * @brief  functions used by MCGeneral
 *
 *** FIXMEDOC ***: mcPathTraceRadar, opt_propExtract, Sample_los_uniform (JUST author)
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

#include <workspace.h>

#include "check_input.h"
#include "optproperties.h"

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
void clear_rt_vars_at_gp(const Workspace& ws,
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
void cloudy_rt_vars_at_gp(const Workspace& ws,
                          MatrixView ext_mat_mono,
                          VectorView abs_vec_mono,
                          VectorView pnd_vec,
                          Numeric& temperature,
                          const Agenda& propmat_clearsky_agenda,
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


/** Sample_los.
 *
 *  Sampling the new incident direction according to a rejection method.
 *  Calculation of the bulk scattering phase matrix.
 *
 * @param[out]    new_rte_los      Incident line of sight for subsequent.
 * @param[out]    g_los_csc_theta  Probability density for the chosen
 *                                 direction multiplied by sin(za)
 * @param[out]    Z                Bulk phase matrix in Stokes notation.
 * @param[in,out] rng              Rng random number generator instance.
 * @param[in]     rte_los          Incident line of sight for subsequent
 *                                 ray-tracing.
 * @param[in]     scat_data        As the WSV.
 * @param[out]    pnd_vec          Vector of particle number densities (one element per scattering element).
 * @param[in]     Z11maxvector     Vector holding the maximum phase function for each scattering element.
 * @param[in]     Csca             Scattering cross section
 * @param[in]     rtp_temperature  As the WSV.
 *
 * @author Cory Davis
 * @date   2003-06-19
 */
void Sample_los(VectorView new_rte_los,
                Numeric& g_los_csc_theta,
                MatrixView Z,
                RandomNumberGenerator<>& rng,
                ConstVectorView rte_los,
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const Index f_index,
                ConstVectorView pnd_vec,
                ConstVectorView Z11maxvector,
                const Numeric Csca,
                const Numeric rtp_temperature,
                const Index t_interp_order = 1);

/** Sample_los_uniform.
 *
 * Sampling the new direction uniformly
 *
 * @param[out] new_rte_los Incident line of sight for subsequent.
 * @param[out] rng         Rng random number generator instance.
 *
 * @author     *** FIXMEDOC ***
 * @date       *** FIXMEDOC ***
 */
void Sample_los_uniform(VectorView new_rte_los, RandomNumberGenerator<>& rng);

#endif  // montecarlo_h
