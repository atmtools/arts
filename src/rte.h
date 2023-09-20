/*===========================================================================
  === File description 
  ===========================================================================*/

/**
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

#include <workspace.h>

#include "ppath.h"

/*===========================================================================
  === Functions in rte.cc
  ===========================================================================*/

/** Adapts clearsky partial derivatives.
 *
 * The following fields:
 * 
 *   Wind
 *   VMR
 * 
 *  Adaptation means changing unit by user input
 * 
 * @param[in,out] dK_dx Propagation matrix derivatives at propagation path point, adapted for wind and VMR units
 * @param[in,out] dS_dx NLTE source adjustment derivatives at propagation path point, adapted for wind and VMR units
 * @param[in] jacobian_quantities as WSV
 * @param[in] ppath_f_grid Wind-adjusted frequency grid at propagation path point
 * @param[in] ppath_line_of_sight Line of sight at propagation path point
 * @param[in] lte Boolean index for whether or not the atmosphere is in LTE at propagation path point
 * @param[in] jacobian_do As WSV
 * 
 * @author Richard Larsson 
 * @date   2017-09-21
 */
void adapt_stepwise_partial_derivatives(
    PropmatMatrix& dK_dx,
    StokvecMatrix& dS_dx,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ConstVectorView& ppath_f_grid,
    const ConstVectorView& ppath_line_of_sight);

/** Performs conversion from radiance to other units, as well as applies
    refractive index to fulfill the n2-law of radiance.

    Use *apply_iy_unit2* for conversion of jacobian data.

    @param[in,out]   iy   Tensor3 with data to be converted, where 
                          column dimension corresponds to Stokes dimensionality
                          and row dimension corresponds to frequency.
    @param[in]   iy_unit  As the WSV.
    @param[in]   f_grid   As the WSV.
    @param[in]   n        Refractive index at the observation position.
    @param[in]   i_pol    Polarisation indexes. See documentation of y_pol.

    @author Patrick Eriksson 
    @date   2010-04-07
 */
void apply_iy_unit(MatrixView iy,
                   const String& iy_unit,
                   const ConstVectorView& f_grid,
                   const Numeric& n,
                   const ArrayOfIndex& i_pol);

/** Largely as *apply_iy_unit* but operates on jacobian data.

    The associated spectrum data *iy* must be in radiance. That is, the
    spectrum can only be converted to Tb after the jacobian data. 

    @param[in,out]   J    Tensor3 with data to be converted, where 
                          column dimension corresponds to Stokes dimensionality
                          and row dimension corresponds to frequency.
    @param[in]   iy       Associated radiance data.
    @param[in]   iy_unit  As the WSV.
    @param[in]   f_grid   As the WSV.
    @param[in]   n        Refractive index at the observation position.
    @param[in]   i_pol    Polarisation indexes. See documentation of y_pol.

    @author Patrick Eriksson 
    @date   2010-04-10
 */
void apply_iy_unit2(Tensor3View J,
                    const ConstMatrixView& iy,
                    const String& iy_unit,
                    const ConstVectorView& f_grid,
                    const Numeric& n,
                    const ArrayOfIndex& i_pol);

/** Calculates the dot product between a field and a LOS

    The line-of-sight shall be given as in the ppath structure (i.e. the
    viewing direction), but the dot product is calculated for the photon
    direction. The field is specified by its three components.

    The returned value can be written as |f|*cos(theta), where |f| is the field
    strength, and theta the angle between the field and photon vectors.

    @param[in]   los               Pppath line-of-sight.
    @param[in]   u                 U-component of field.
    @param[in]   v                 V-component of field.
    @param[in]   w                 W-component of field.

    @return   The result of the dot product

    @author Patrick Eriksson 
    @date   2012-12-12
*/
Numeric dotprod_with_los(const ConstVectorView& los,
                         const Numeric& u,
                         const Numeric& v,
                         const Numeric& w);

/** Basic call of *iy_main_agenda*.

    This function is an interface to *iy_main_agenda* that can be used when
    only *iy* is of interest. That is, jacobian and auxilary parts are
    deactivated/ignored.

    @param[in,out]   ws                The workspace
    @param[out]   iy                   As the WSV.
    @param[in]   cloudbox_on           As the WSV.
    @param[in]   f_grid                As the WSV.
    @param[in]   nlte_field            As the WSV.
    @param[in]   rte_pos               As the WSV.
    @param[in]   rte_los               As the WSV.
    @param[in]   iy_unit               As the WSV.
    @param[in]   iy_main_agenda        As the WSV.

    @author Patrick Eriksson 
    @date   2012-08-08
 */
void get_iy(const Workspace& ws,
            Matrix& iy,
            const Index& cloudbox_on,
            const Vector& f_grid,
            const AtmField& atm_field,
            const Vector& rte_pos,
            const Vector& rte_los,
            const Vector& rte_pos2,
            const String& iy_unit,
            const Agenda& iy_main_agenda);

/** Determines iy of the "background" of a propgation path.

    The task is to determine *iy* and related variables for the
    background, or to continue the raditiave calculations
    "backwards". The details here depends on the method selected for
    the agendas.

    Each background is handled by an agenda. Several of these agandes
    can involve recursive calls of *iy_main_agenda*. 

    @param[in,out] ws                  The workspace
    @param[out]  iy                    As the WSV.
    @param[out]  diy_dx                As the WSV.
    @param[in]   iy_transmittance       As the WSV.
    @param[in]   jacobian_do           As the WSV.
    @param[in]   ppath                 As the WSV.
    @param[in]   cloudbox_on           As the WSV.
    @param[in]   f_grid                As the WSV.
    @param[in]   iy_unit               As the WSV.    
    @param[in]   surface_props_data    As the WSV.    
    @param[in]   iy_main_agenda        As the WSV.
    @param[in]   iy_space_agenda       As the WSV.
    @param[in]   iy_surface_agenda     As the WSV.
    @param[in]   iy_cloudbox_agenda    As the WSV.

    @author Patrick Eriksson 
    @date   2009-10-08
 */
void get_iy_of_background(const Workspace& ws,
                          Matrix& iy,
                          ArrayOfTensor3& diy_dx,
                          const Tensor3& iy_transmittance,
                          const Index& iy_id,
                          const Index& jacobian_do,
                          const ArrayOfRetrievalQuantity& jacobian_quantities,
                          const Ppath& ppath,
                          const Vector& rte_pos2,
                          const AtmField& atm_field,
                          const Index& cloudbox_on,
                          const Vector& f_grid,
                          const String& iy_unit,
                          const SurfaceField& surface_field,
                          const Agenda& iy_main_agenda,
                          const Agenda& iy_space_agenda,
                          const Agenda& iy_surface_agenda,
                          const Agenda& iy_cloudbox_agenda,
                          const Index& iy_agenda_call1y);

/** Determines the particle fields along a propagation path.

    @param[out]  clear2cloudy        Mapping of index. See code for details. 
    @param[out]  ppath_pnd           The particle number density for each
                                     path point (also outside cloudbox).
    @param[out]  ppath_dpnd_dx       dpnd_field_dx for each path point
                                     (also outside cloudbox).
    @param[in]   ppath               As the WSV.    
    @param[in]   cloubox_limits      As the WSV.    
    @param[in]   pnd_field           As the WSV.    
    @param[in]   dpnd_field_dx       As the WSV.    

    @author Jana Mendrok, Patrick Eriksson 
    @date   2017-09-18
 */
void get_ppath_cloudvars(ArrayOfIndex& clear2cloudy,
                         Matrix& ppath_pnd,
                         ArrayOfMatrix& ppath_dpnd_dx,
                         const Ppath& ppath,
                         const ArrayOfIndex& cloudbox_limits,
                         const Tensor4& pnd_field,
                         const ArrayOfTensor4& dpnd_field_dx);

/** Determines the Doppler shifted frequencies along the propagation path.

    ppath_doppler[ nf,np]

    @param[out]  ppath_f          Doppler shifted f_grid
    @param[in]   ppath            Propagation path.
    @param[in]   f_grid           Original f_grid.
    @param[in]   rte_alonglos_v   As the WSV.
    @param[in]   ppath_wind       See get_ppath_atmvars.

    @author Patrick Eriksson 
    @date   2013-02-21
 */
void get_ppath_f(Matrix& ppath_f,
                 const Ppath& ppath,
                 const ConstVectorView& f_grid,
                 const Numeric& rte_alonglos_v,
                 const ConstMatrixView& ppath_wind);

/** Get the blackbody radiation at propagation path point
 * 
 * @param[in,out] B Blackbody radiation at propagation path point
 * @param[in,out] dB_dT Blackbody radiation temperature derivative at propagation path point
 * @param[in] ppath_f_grid Wind-adjusted frequency grid at propagation path point
 * @param[in] ppath_temperature Temperature of atmosphere at propagation path point
 * @param[in] do_temperature_derivative Fill dB_dT?
 * 
 * @author Richard Larsson 
 * @date   2017-09-21
 */
void get_stepwise_blackbody_radiation(VectorView B,
                                      VectorView dB_dT,
                                      const ConstVectorView& ppath_f_grid,
                                      const Numeric& ppath_temperature,
                                      const bool& do_temperature_derivative);

/** Gets the clearsky propgation matrix and NLTE contributions
 * 
 * Basically a wrapper for calls to the propagation clearsky agenda
 * 
 * @param[in] ws The workspace
 * @param[in,out] K Propagation matrix at propagation path point
 * @param[in,out] S NLTE source adjustment at propagation path point
 * @param[in,out] dK_dx Propagation matrix derivatives at propagation path point
 * @param[in,out] dS_dx NLTE source adjustment derivatives at propagation path point
 * @param[in] propmat_clearsky_agenda As WSA
 * @param[in] jacobian_quantities As WSV
 * @param[in] ppath_f_grid Wind-adjusted frequency grid at propagation path point
 * @param[in] ppath_line_of_sight Line of sight at propagation path point
 * @param[in] atm_point As WSV
 * @param[in] jacobian_do As WSV
 * 
 * @author Richard Larsson 
 * @date   2017-09-21
 */
void get_stepwise_clearsky_propmat(
  const Workspace& ws,
  PropmatVector& K,
  StokvecVector& S,
  PropmatMatrix& dK_dx,
  StokvecMatrix& dS_dx,
  const Agenda& propmat_clearsky_agenda,
  const ArrayOfRetrievalQuantity& jacobian_quantities,
  const Vector& ppath_f_grid,
  const Vector& ppath_line_of_sight,
  const AtmPoint& atm_point,
  const bool jacobian_do);

/** Computes the ratio that a partial derivative with regards to frequency
 *  relates to the wind of come component
 * 
 * @param[in,out] ppath_f_grid Wind-adjusted frequency grid wind derivative at propagation path point
 * @param[in] ppath_line_of_sight Line of sight at propagation path point
 * @param[in] f_grid As WSV
 * @param[in] wind_type The wind component
 * 
 * @author Richard Larsson 
 * @date   2017-09-21
 */
Vector get_stepwise_f_partials(const ConstVectorView& ppath_line_of_sight,
                               const ConstVectorView& f_grid,
                               const Jacobian::Atm wind_type);

/** Computes the contribution by scattering at propagation path point
 * 
 * @param[in,out] ap Absorption vector scattersky at propagation path point
 * @param[in,out] Kp Propagation matrix scattersky at propagation path point
 * @param[in,out] dap_dx Absorption vector scattersky derivatives at propagation path point
 * @param[in,out] dKp_dx Propagation matrix scattersky derivatives at propagation path point
 * @param[in] jacobian_quantities As WSV
 * @param[in] ppath_1p_pnd Particulate number density at propagation path point
 * @param[in] scat_data As WSV
 * @param[in] ppath_line_of_sight Line of sight at propagation path point
 * @param[in] ppath_temperature Temperature at propagation path point
 * @param[in] jacobian_do As WSV
 * 
 *  @author Jana Mendrok, Richard Larsson 
 *  @date   2017-09-21
 */
void get_stepwise_scattersky_propmat(
    StokvecVector& ap,
    PropmatVector& Kp,
    StokvecMatrix& dap_dx,
    PropmatMatrix& dKp_dx,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ConstMatrixView& ppath_1p_pnd,  // the ppath_pnd at this ppath point
    const ArrayOfMatrix&
        ppath_dpnd_dx,  // the full ppath_dpnd_dx, ie all ppath points
    const Index ppath_1p_id,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const ConstVectorView& ppath_line_of_sight,
    const ConstVectorView& ppath_temperature,
    const bool& jacobian_do);

/**
 *  Calculates the stepwise scattering source terms.
 *  Uses new, unified phase matrix extraction scheme.
 * 
 *  @param[in]   Sp                    Out: The scattering source term
 *  @param[in]   dSp_dx                Out: The derivative of the scattering source term
 ...
 * 
 *  @author Jana Mendrok 
 *  \adapted from non-stepwise function
 *  @date   2018-03-29
 */
void get_stepwise_scattersky_source(
    StokvecVector& Sp,
    StokvecMatrix& dSp_dx,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ConstVectorView& ppath_1p_pnd,
    const ArrayOfMatrix& ppath_dpnd_dx,
    const Index ppath_1p_id,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const ConstTensor7View& cloudbox_field,
    const ConstVectorView& za_grid,
    const ConstVectorView& aa_grid,
    const ConstMatrixView& ppath_line_of_sight,
    const GridPos& ppath_pressure,
    const Vector& temperature,
    const bool& jacobian_do,
    const Index& t_interp_order = 1);

/** Computes layer transmission matrix and cumulative transmission
 * 
 * FIXME: This function should be removed
 * 
 * @param[in,out] cumulative_transmission Present accumulation of transmission for Jacobian computations
 * @param[in,out] T Layer transmission
 * @param[in,out] dT_close_dx Layer transmission derivative due to the close propagation path point
 * @param[in,out] dT_far_dx Layer transmission derivative due to the far propagation path point
 * @param[in] cumulative_transmission_close Past accumulation of transmission for Jacobian computations
 * @param[in] K_close Level propagation matrix due to the close propagation path point
 * @param[in] K_far Level propagation matrix due to the far propagation path point
 * @param[in] dK_close_dx Level propagation matrix derivatives due to the close propagation path point
 * @param[in] dK_far_dx Level propagation matrix derivatives due to the far propagation path point
 * @param[in] ppath_distance Thickness of the layer
 * @param[in] first_level Boolean for if this is the first level, i.e., there is no cumulative_transmission_close
 * @param[in] dr_dT_close Thickness of the layer derivative due to temperature of the close propagation path point
 * @param[in] dr_dT_far Thickness of the layer derivative due to temperature of the far propagation path point
 * @param[in] it Index to temperature derivatives
 * 
 * @author Richard Larsson 
 * @date   2017-09-21
 */
void get_stepwise_transmission_matrix(
    Tensor3View cumulative_transmission,
    Tensor3View T,
    Tensor4View dT_dx_close,
    Tensor4View dT_dx_far,
    const ConstTensor3View& cumulative_transmission_close,
    const PropmatVector& K_close,
    const PropmatVector& K_far,
    const PropmatMatrix& dK_close_dx,
    const PropmatMatrix& dK_far_dx,
    const Numeric& ppath_distance,
    const bool& first_level,
    const Numeric& dr_dT_close = 0,
    const Numeric& dr_dT_far = 0,
    const Index& it = -1);

/** Multiplicates iy_transmittance with transmissions.

    That is, a multiplication of *iy_transmittance* with another
    variable having same structure and holding transmission values.

    The "new path" is assumed to be further away from the sensor than 
    the propagtion path already included in iy_transmittance. That is,
    the operation can be written as:
    
       Ttotal = Told * Tnew

    where Told is the transmission corresponding to *iy_transmittance*
    and Tnew corresponds to *tau*.

    *iy_trans_new* is sized by the function.

    @param[out]   iy_trans_total    Updated version of *iy_transmittance*
    @param[in]   iy_trans_old      A variable matching *iy_transmittance*.
    @param[in]   iy_trans_new      A variable matching *iy_transmittance*.

    @author Patrick Eriksson 
    @date   2009-10-06
*/
void iy_transmittance_mult(Tensor3& iy_trans_total,
                          const ConstTensor3View& iy_trans_old,
                          const ConstTensor3View& iy_trans_new);

/** Multiplicates iy_transmittance with iy-variable.

    The operation can be written as:
    
       iy_new = T * iy_old

    where T is the transmission corresponding to *iy_transmittance*
    and iy_old is a variable matching iy.

    *iy_new* is sized by the function.

    @param[out]   iy_new        Updated version of iy 
    @param[in]   iy_trans      A variable matching *iy_transmittance*.
    @param[in]   iy_old        A variable matching *iy*.

    @author Patrick Eriksson 
    @date   2018-04-10
*/
void iy_transmittance_mult(Matrix& iy_new,
                           const ConstTensor3View& iy_trans,
                           const ConstMatrixView& iy_old);

/** Determines the backward direction for a given line-of-sight.

    ZZZ Use los_reverse in geodetic.h ZZZ

    This function can be used to get the LOS to apply for extracting single
    scattering properties, if the propagation path LOS is given.

    A viewing direction of aa=0 is assumed for 1D. This corresponds to 
    positive za for 2D.

    @param[out]  los_mirrored      The line-of-sight for reversed direction.
    @param[in]   los               A line-of-sight

    @author Patrick Eriksson 
    @date   2011-07-15
*/
void mirror_los(Vector& los_mirrored,
                const ConstVectorView& los);

//! mueller_modif2stokes
/*!
   Returns the Mueller matrix for transformation of a modified Stokes vector to
   its standard counterpart.

   See ARTS Theory document, section "Change of the Stokes basis" for details.

   \param   Cs          Mueller matrix

   \author Patrick Eriksson
   \date   2021-12-22
*/
void mueller_modif2stokes(Matrix& Cs);

//! mueller_rotation
/*!
   Returns the Mueller matrix for a rotation of the coordinate system defining
   H and V directions.

   As muellersparse_rotation, besides that this function returns a matrix (not
   Sparse) and the matrix is sized by the function.

   \param   L           Mueller matrix
   \param   rotangle    Rotation angle.

   \author Patrick Eriksson
   \date   2021-12-22
*/
void mueller_rotation(Matrix& L,
                      const Numeric& rotangle);

//! mueller_stokes2modif
/*!
   Returns the Mueller matrix for transformation from Stokes to its modified
   counterpart.

   See ARTS Theory document, section "Change of the Stokes basis" for details.

   \param   Cm          Mueller matrix

   \author Patrick Eriksson
   \date   2021-12-22
*/
void mueller_stokes2modif(Matrix& Cm);

/** This function fixes the last steps to made on the Jacobian in some
    radiative transfer WSMs. The method applies iy_transmittance, maps from
    ppath to the retrieval grids and applies non-standard Jacobian units.

    See iyEmissonStandard for usage example.

    @author Patrick Eriksson 
    @date   2017-11-19
*/
void rtmethods_jacobian_finalisation(
    const Workspace& ws,
    ArrayOfTensor3& diy_dx,
    ArrayOfTensor3& diy_dpath,
    const Index& nf,
    const Index& np,
    const Ppath& ppath,
    const ArrayOfAtmPoint& ppvar_atm,
    const Index& iy_agenda_call1,
    const Tensor3& iy_transmittance,
    const Agenda& water_p_eq_agenda,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfIndex& jac_species_i);

/** This function handles the unit conversion to be done at the end of some
    radiative transfer WSMs. 

    The method hanldes both *iy* and analytical parts of the Jacobian.

    See iyEmissonStandard for usage example.

    @author Patrick Eriksson 
    @date   2017-11-19
 */
void rtmethods_unit_conversion(
    Matrix& iy,
    ArrayOfTensor3& diy_dx,
    Tensor3& ppvar_iy,
    const Vector& f_grid,
    const Ppath& ppath,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Index& j_analytical_do,
    const String& iy_unit);

/** Calculates factor to convert back-scattering to Ze

   The vector *fac* shall be sized to match f_grid, before calling the
   function.

   If k2 <= 0, the K" factor is calculated. Otherwise the input k2 is applied
   as "hard-coded".

   @param[out]   fac     Vector with factors.
   @param[in]   f_grid   As the WSV.
   @param[in]   z_tref   Reference temperature for conversion to Ze.
   @param[in]   k2       Reference dielectric factor.

   @author Patrick Eriksson
   @date   2002-05-20
*/
void ze_cfac(Vector& fac,
             const Vector& f_grid,
             const Numeric& ze_tref,
             const Numeric& k2);

/** Get the stepwise blackbody radiation object
 * 
 * @param B Plank function at propagation path point [size is (ppath_f_grid.nelem())]
 * @param dB Derivative of Plank function at propagation path point [size is (jacobian_quantities.nelem() x ppath_f_grid.nelem())]
 * @param ppath_f_grid Frequency grid at propagation path point
 * @param ppath_temperature Temperature of atmosphere at propagation path point
 * @param jacobian_quantities As WSV
 * @param j_analytical_do Flag for analytical Jacobian
 */
void get_stepwise_blackbody_radiation(
    Vector &B, Matrix &dB, const Vector &ppath_f_grid,
    const Numeric &ppath_temperature,
    const ArrayOfRetrievalQuantity &jacobian_quantities,
                            const bool j_analytical_do);

#endif  // rte_h
