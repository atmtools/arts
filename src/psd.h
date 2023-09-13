/**
 * @file   psd.h
 * @author Jana Mendrok Patrick Eriksson
 * @date   2017-11-05 
 * 
 * @brief  Internal functions associated with size distributions
 */

#ifndef psd_h
#define psd_h

#include "array.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "matpack_data.h"
#include "optproperties.h"
#include "ppath.h"

// ------------------------------------------------------
// Macros to avoid duplication of code
// ------------------------------------------------------

#define START_OF_PSD_METHODS()                                                 \
  const Index nin = pnd_agenda_input_names.size();                            \
  const Index ndx = dpnd_data_dx_names.size();                                \
  const Index np = pnd_agenda_input.nrows();                                   \
  const Index nsi = psd_size_grid.size();                                     \
  ArrayOfIndex dx2in(ndx);                                                     \
                                                                               \
  if (pnd_agenda_input.ncols() != nin)                                         \
    throw std::runtime_error(                                                       \
        "Length of *pnd_agenda_input_names* and number of "                    \
        "columns in *pnd_agenda_input* must be equal.");                       \
  if (ndx) {                                                                   \
    if (ndx > nin)                                                             \
      throw std::runtime_error(                                                     \
          "The length of *dpnd_data_dx_names* can not "                        \
          "exceed the one of *pnd_agenda_input_names*.");                      \
    for (Index i = 0; i < ndx; i++) {                                          \
      dx2in[i] = find_first(pnd_agenda_input_names, dpnd_data_dx_names[i]);    \
      if (dx2in[i] < 0) {                                                      \
        std::ostringstream os;                                                      \
        os << "dpnd_data_dx_names[" << i << "] is " << dpnd_data_dx_names[i]   \
           << "\nThis string could not be found in *pnd_agenda_input_names*."; \
        throw std::runtime_error(os.str());                                    \
      }                                                                        \
    }                                                                          \
  }                                                                            \
                                                                               \
  psd_data.resize(np, nsi);                                                    \
  psd_data = 0.0;                                                              \
  if (ndx) {                                                                   \
    dpsd_data_dx.resize(ndx, np, nsi);                                         \
    dpsd_data_dx = 0.0;                                                        \
  } else {                                                                     \
    dpsd_data_dx.resize(0, 0, 0);                                              \
  }

/** The MH97 cloud ice PSD
 *  
 * Handles a vector of sizes at a time. Implicitly assumes particles of water
 * ice. Strictly requires IWC and T to be positive, i.e. calling method needs
 * to ensure this.
 *  
 * @param[out] psd      particle number density per size interval [#/m3*m]
 * @param[in] diameter  size of the scattering elements (supposed to be mass (aka
 *                      volume) equivalent diameter of pure ice particle) [m]
 * @param[in] iwc       atmospheric ice water content [kg/m3]
 * @param t             atmospheric temperature [K]
 * @param noisy         flag whether to add noise onto PSD parameters according to
 *                      their reported error statistics
 * 
 * @author Jana Mendrok, Daniel Kreyling
 * @date 2017-06-07
 */
void psd_cloudice_MH97(Vector& psd,
                       const Vector& diameter,
                       const Numeric& iwc,
                       const Numeric& t,
                       const bool noisy);

/** Code common to MGD PSD involving the integrated mass
 *
 * Valid choices for *something* is "mean size", "median size", "mean particle
 * mass" and "Ntot".
 *
 * @param[out] psd_data                 As the WSV wih same name
 * @param[out] dpsd_data_dx             As the WSV wih same name
 * @param[in]  something                Text string gving the second moment
 * @param[in]  psd_size_grid            As the WSV wih same name
 * @param[in]  pnd_agenda_input_t       As the WSV wih same name
 * @param[in]  pnd_agenda_input         As the WSV wih same name
 * @param[in]  pnd_agenda_input_names   As the WSV wih same name
 * @param[in]  dpnd_data_dx_names       As the WSV wih same name
 * @param[in]  scat_species_a           As the WSV wih same name
 * @param[in]  scat_species_b           As the WSV wih same name
 * @param[in]  n0                       Selection for N0 parameter
 * @param[in]  mu                       Selection for mu parameter
 * @param[in]  la                       Selection for lambda parameter
 * @param[in]  ga                       Selection for gamma parameter
 * @param[in]  t_min                    PSD set to zero below this temperature
 * @param[in]  t_max                    PSD set to zero above this temperature
 * @param[in]  picky                    Triggers more check of input
 *
 * @author Patrick Eriksson
 * @date 2017-11-05
 */
void psd_mgd_mass_and_something(Matrix& psd_data,
                                Tensor3& dpsd_data_dx,
                                const String& something,
                                const Vector& psd_size_grid,
                                const Vector& pnd_agenda_input_t,
                                const Matrix& pnd_agenda_input,
                                const ArrayOfString& pnd_agenda_input_names,
                                const ArrayOfString& dpnd_data_dx_names,
                                const Numeric& scat_species_a,
                                const Numeric& scat_species_b,
                                const Numeric& n0,
                                const Numeric& mu,
                                const Numeric& la,
                                const Numeric& ga,
                                const Numeric& t_min,
                                const Numeric& t_max,
                                const Index& picky);

/** Code common to PSDs of mono type
 *
 * @param[out] psd_data                 As the WSV wih same name
 * @param[out] dpsd_data_dx             As the WSV wih same name
 * @param[in]  pnd_agenda_input_t       As the WSV wih same name
 * @param[in]  pnd_agenda_input         As the WSV wih same name
 * @param[in]  pnd_agenda_input_names   As the WSV wih same name
 * @param[in]  dpnd_data_dx_names       As the WSV wih same name
 * @param[in]  scat_meta                As the WSV wih same name
 * @param[in]  species_index            Index of the scattering species of concern
 * @param[in]  t_min                    PSD set to zero below this temperature
 * @param[in]  t_max                    PSD set to zero above this temperature
 * @param[in]  picky                    Triggers more check of input
 *
 * @author Patrick Eriksson
 * @date 2017-11-05
 */
void psd_mono_common(Matrix& psd_data,
                     Tensor3& dpsd_data_dx,
                     const String& type,
                     const Vector& pnd_agenda_input_t,
                     const Matrix& pnd_agenda_input,
                     const ArrayOfString& pnd_agenda_input_names,
                     const ArrayOfString& dpnd_data_dx_names,
                     const ArrayOfArrayOfScatteringMetaData& scat_meta,
                     const Index& species_index,
                     const Numeric& t_min,
                     const Numeric& t_max,
                     const Index& picky);

/** The Wang16 rain DSD DEPRECATED BY NEW MGD_SMM_COMMON
 *  Only included for compatibility with "old" pnd_fieldCalcFromscat_speciesField
 *  
 *  Uses rain water water content. PSD follows an exponential distribution.
 *  Handles a vector of sizes at a time.
 *  
 *  Reference: Wang et al., 2016, "Investigation of liquid cloud microphysical
 *  properties of deep convective systems: 1. Parameterization raindrop size
 *  distribution and its application for stratiform rain estimation".
 *  Ported from CloudArts matlab implementation.
 *
 *  @param[out] psd       particle number density per size interval [#/m3*m]
 *  @param[in]  diameter  size of the scattering elements (volume equivalent
 *                        diameter) [m]
 *  @param[in]  rwc       atmospheric rain water content [kg/m3]
 * 
 * @author Jana Mendrok, Bengt Rydberg
 * @date 2017-06-07
 */
void psd_rain_W16(Vector& psd, const Vector& diameter, const Numeric& rwc);

/** Code common to a number of modified gamma PSDs used with single-moment
 *  mass schemes. All PSDs take the form n(D) = n_alpha lam^n_b D^mu exp(-lam*D^gamma).
 *
 * @param[out] psd_data                 As the WSV wih same name
 * @param[out] dpsd_data_dx             As the WSV wih same name
 * @param[in]  psd_name                 Name of selected PSD
 * @param[in]  psd_size_grid            As the WSV wih same name
 * @param[in]  pnd_agenda_input_t       As the WSV wih same name
 * @param[in]  pnd_agenda_input         As the WSV wih same name
 * @param[in]  pnd_agenda_input_names   As the WSV wih same name
 * @param[in]  dpnd_data_dx_names       As the WSV wih same name
 * @param[in]  scat_species_a           As the WSV wih same name
 * @param[in]  scat_species_b           As the WSV wih same name
 * @param[in]  t_min                    PSD set to zero below this temperature
 * @param[in]  t_max                    PSD set to zero above this temperature
 * @param[in]  n_alpha_in               Value of n_alpha for psd_name "generic"
 * @param[in]  n_b_in                   Value of n_b for psd_name "generic"
 * @param[in]  mu_in                    Value of mu for psd_name "generic"
 * @param[in]  gamma_in                 Value of gamma for psd_name "generic"
 * @param[in]  picky                    Triggers more check of input
 *
 * @author Stuart Fox
 * @date 2019-10-02
 */
void psd_mgd_smm_common(Matrix& psd_data,
                        Tensor3& dpsd_data_dx,
                        const String& psd_name,
                        const Vector& psd_size_grid,
                        const Vector& pnd_agenda_input_t,
                        const Matrix& pnd_agenda_input,
                        const ArrayOfString& pnd_agenda_input_names,
                        const ArrayOfString& dpnd_data_dx_names,
                        const Numeric& scat_species_a,
                        const Numeric& scat_species_b,
                        const Numeric& n_alpha_in,
                        const Numeric& n_b_in,
                        const Numeric& mu_in,
                        const Numeric& gamma_in,
                        const Numeric& t_min,
                        const Numeric& t_max,
                        const Index& picky);

/** The F07 snow PSD
 *
 *  Handles a vector of sizes at a time.
 *  Strictly requires SWC and T to be positive and regime to be either "TR" or
 *  "ML", i.e. calling methods need to ensure these.
 *  No further limitations on the allowed temperatures here. Strictly valid it's
 *  only within -60<=t<=0C, the measured t-range the parametrization is based
 *  on. However, this is left to be handled by the calling methods.
 *
 *  @param[out] psd      particle number density per size interval [#/m3*m]
 *  @param[in] diameter  size of the scattering elements (supposed to be maximum
 *                       diameter of the ice particles) [m]
 *  @param[in] swc       atmospheric snow water content [kg/m^3]
 *  @param[in] t         atmospheric temperature [K]
 *  @param[in] alpha     mass-dimension relationship scaling factor
 *                       (m=alpha*(Dmax/D0)^beta) [kg]
 *  @param[in] beta      mass-dimension relationship exponent [-]
 *  @param[in] regime    parametrization regime to apply (TR=tropical, ML=midlatitude)
 * 
 * @author Manfred Brath, Jana Mendrok
 * @date 2017-06-13
 */
void psd_snow_F07(Vector& psd,
                  const Vector& diameter,
                  const Numeric& swc,
                  const Numeric& t,
                  const Numeric alpha,
                  const Numeric beta,
                  const String& regime);

/*! Calculates the particle number density field according
 *  to the two moment scheme of Seifert and Beheng, 2006b,a that is used 
 *  in the ICON model.
 *  One call of this function calculates one particle number density.
 
 \return dN particle number density per diameter interval [#/m3/m]

 \param mass   Mass of scattering particle [kg]
 \param N_tot  Total number of particles (0th moment) [#/m3/m/kg^mu]
 \param M      Total mass concentration of Particles (1st moment) [kg/m^3]
 \param psd_type string with a tag defining the (hydrometeor) scheme

 \author Manfred Brath
 \date 2015-01-19
 */
void psd_SB06(Vector& psd,
              Matrix& dpsd,
              const Vector& mass,
              const Numeric& N_tot,
              const Numeric& WC,
              const String& hydrometeor_type);

/*! Calculates the particle number density field according
 *  to the Milbrandt and Yau two moment scheme, which is used in the GEM model.
 *  See also milbrandt and yau, 2005.
 *  One call of this function calculates one particle number density.

 \return dN particle number density per diameter interval [#/m3/m]

 \param mass   Mass of scattering particle [kg]
 \param N_tot  Total number of particles (0th moment) [#/m3/m/kg^mu]
 \param M      Total mass concentration of Particles (1st moment) [kg/m^3]
 \param psd_type string with a tag defining the (hydrometeor) scheme

 \author Manfred Brath
 \date 2017-08-01
 */
void psd_MY05(Vector& psd,
              Matrix& dpsd,
              const Vector& diameter_max,
              const Numeric N_tot,
              const Numeric WC,
              const String psd_type);

/** Derives Dm from IWC and N0star
 *
 * For definition of Dm and N0star follows the DARDAR PSD.
 *
 * @param[in]   iwc   IWC-value  
 * @param[in]   n0    N0star value
 * @param[in]   rho   Density (of water or ice) 
 *
 * @return   Derived Dm-value
 *
 * @author Simon Pfreundschuh 
 * @date 20??-??-??
 */
Numeric dm_from_iwc_n0(Numeric iwc, Numeric n0, Numeric rho);

/** Derives N0star from IWC and Dm
 *
 * For definition of Dm and N0star follows the DARDAR PSD.
 *
 * @param[in]   iwc   IWC-value  
 * @param[in]   dm    Dm-value
 * @param[in]   rho   Density (of water or ice) 
 *
 * @return   Derived N0star-value
 *
 * @author Simon Pfreundschuh 
 * @date 20??-??-??
 */
Numeric n0_from_iwc_dm(Numeric iwc, Numeric dm, Numeric rho);

/** Sets N0star based on temperature
 *
 * For definition of N0star follows the DARDAR PSD.
 *
 * This a priori parameterisation is taken from Table 5 in "Normalized 
 * particle size distribution for remote sensing application" by Delanoë 
 * et al., JGR, 2014.
 *
 * @param[in]   t     temperature
 *
 * @return   Derived N0star-value
 *
 * @author Simon Pfreundschuh 
 * @date 20??-??-??
 */
Numeric n0_from_t(Numeric t);


#endif  //psd_h
