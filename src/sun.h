/*!
  \file   sun.h
  \author Jon Petersen <jon.petersen@studium.uni-hamburg.de>
          Manfred Brath  <manfred.brath@uni-hamburg.de>
  \date   2021-02-22

  \brief  Declaration of functions in star.cc.
*/

#ifndef star_h
#define star_h

/*===========================================================================
  === External declarations
  ===========================================================================*/


#include "arts.h"
#include "atm.h"
#include "gridded_fields.h"
#include "matpack_concepts.h"
#include "surf.h"
#include <matpack.h>
#include <rtepack.h>
#include "optproperties.h"
#include "ppath_struct.h"
#include "jacobian.h"


class Agenda;
class Workspace;

/*===========================================================================
  === structs/classes  in sun.h
  ===========================================================================*/

/** The structure to describe a propagation path and releated quantities.
 *
 *  The fields of the structure are described more in detail inside the ARTS
 *  user guide (AUG).
 */
struct Sun {
  /** Sun description */
  String description;
  /** Sun spectrum, monochrmatic radiance spectrum at the surface of the sun*/
  Matrix spectrum;
  /** Sun radius */
  Numeric radius;
  /** distance from center of planet to center of sun*/
  Numeric distance;
  /** latitude of the sun in the sky of the planet */
  Numeric latitude;
  /** longitude of the sun in the sky of the planet */
  Numeric longitude;

  friend std::ostream& operator<<(std::ostream& os, const Sun& sun);
};



/** An array of sun. */
using ArrayOfSun = Array<Sun>;


/*===========================================================================
  === Functions in sun.h
  ===========================================================================*/

/** Calculates the radiance spectrum of sun which is scattered by the atmospheric gases.
 *
 * @param[in,out] ws ARTS workspace.
 * @param[out] scattered_sunlight StokvecVector scattered monochromatic radiance
 *             spectrum of sun.
 * @param[in] f_grid Vector frequency grid.
 * @param[in] p Numeric pressure at location of scattering.
 * @param[in] T Numeric temperature at location of scattering.
 * @param[in] vmr Vector volume mixing ratios of absorption species at location
 *            of scattering.
 * @param[in] transmitted_sunlight Matrix transmitted monochromatic irradiance
 *             spectrum of sun at location of scattering.
 * @param[in] gas_scattering_los_in Vector incoming direction of the transmitted sun irradiance
 *            spectrum.
 * @param[in] gas_scattering_los_out outgoing direction of the transmitted sun irradiance
 *            spectrum.
 * @param[in] gas_scattering_agenda Agenda agenda calculating the gas scattering
 *            cross sectionand matrix.
 */
void get_scattered_sunsource(const Workspace& ws,
                              StokvecVector& scattered_sunlight,
                              const Vector& f_grid,
                              const AtmPoint& atm_point,
                              const Matrix& transmitted_sunlight,
                              const Vector& gas_scattering_los_in,
                              const Vector& gas_scattering_los_out,
                              const Agenda& gas_scattering_agenda);

/** Gets the sun background for a given ppath.
 *
 * iy is zero if there is no sun in the line of sight at TOA.
 *
 * @param[out] iy Matrix radiance spectrum of suns.
 * @param[out] suns_visible Index indicating if suns are in los
 * @param[in] suns Array of sun (structures).
 * @param[in] ppath Propagation path as the WSV.
 * @param[in] f_grid Vector as the WSV.
 * @param[in] refellipsoid Reference ellipsoid
 */
void get_sun_background(Matrix& iy,
                         Index& suns_visible,
                         const ArrayOfSun& suns,
                         const Ppath& ppath,
                         const Vector& f_grid,
                         const Vector2 refellipsoid);

/** Checks and adds sun radiance if sun is in line of sight.
 *
 * @param[in, out] iy Matrix of sun.
 * @param[out] suns_visible Index indicating if sun are in los
 * @param[in] sun Sun-structure.
 * @param[in] rtp_pos The position of the ppath point.
 * @param[in] rtp_los The line of sight of the ppath.
 * @param[in] refellipsoid Reference ellipsoid
  */
void get_sun_radiation(Matrix& iy,
                        Index& suns_visible,
                         const Sun& sun,
                         const Vector& rtp_pos,
                         const Vector& rtp_los,
                         const Vector2 refellipsoid);

/** Calculates the transmitted sun radiation at the end position of the ppath
 *
 * @param[in,out] ws ARTS workspace.
 * @param[out] direct_radiation Matrix Transmitted monochromatic irradiance
 *             spectrum of sun.
 * @param[out] ddirect_radiation_dx Array of Tensor3 Jacobian of transmitted
 *              monochromatic irradiance spectrum of sun.
 * @param[in] f_grid As the WSV.
 * @param[in] p_grid As the WSV.
 * @param[in] lat_grid As the WSV.
 * @param[in] lon_grid As the WSV.
 * @param[in] t_field As the WSV.
 * @param[in] nlte_field As the WSV.
 * @param[in] vmr_field As the WSV.
 * @param[in] abs_species As the WSV.
 * @param[in] wind_u_field As the WSV.
 * @param[in] wind_v_field As the WSV.
 * @param[in] wind_w_field As the WSV.
 * @param[in] mag_u_field As the WSV.
 * @param[in] mag_v_field As the WSV.
 * @param[in] mag_w_field As the WSV.
 * @param[in] cloudbox_on As the WSV.
 * @param[in] cloudbox_limits As the WSV.
 * @param[in] gas_scattering_do As the WSV.
 * @param[in] irradiance_flag Index Flag indicating if the transmitted radiation
 *                            is spectral irradiance (1) or if it is spectral
 *                            radiance (0).
 * @param[in] sun_ppaths ArrayOfPpath Propagation path towards each sun.
 * @param[in] suns As the WSV.
 * @param[in] suns_visible ArrayOfIndex Flag indicating if eah sun is visible.
 * @param[in] pnd_field As the WSV.
 * @param[in] dpnd_field_dx As the WSV.
 * @param[in] scat_species As the WSV.
 * @param[in] scat_data As the WSV.
 * @param[in] jacobian_do As the WSV.
 * @param[in] jacobian_quantities As the WSV.
 * @param[in] propmat_clearsky_agenda As the WSV.
 * @param[in] water_p_eq_agenda As the WSV.
 * @param[in] gas_scattering_agenda As the WSV.
 * @param[in] ppath_step_agenda As the WSV.
 * @param[in] rte_alonglos_v As the WSV.
 */
void get_direct_radiation(const Workspace& ws,
                     ArrayOfMatrix& direct_radiation,
                     ArrayOfArrayOfTensor3& ddirect_radiation_dx,
                     const Vector& f_grid,
                     const ArrayOfArrayOfSpeciesTag& abs_species,
                     const AtmField& atm_field,
                     const Index& cloudbox_on,
                     const ArrayOfIndex& cloudbox_limits,
                     const Index& gas_scattering_do,
                     const Index& irradiance_flag,
                     const ArrayOfPpath& sun_ppaths,
                     const ArrayOfSun& suns,
                     const ArrayOfIndex& suns_visible,
                     const Vector2 refellipsoid,
                     const Tensor4& pnd_field,
                     const ArrayOfTensor4& dpnd_field_dx,
                     const ArrayOfString& scat_species,
                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                     const Index& jacobian_do,
                     const ArrayOfRetrievalQuantity& jacobian_quantities,
                     const Agenda& propmat_clearsky_agenda,
                     const Agenda& water_p_eq_agenda,
                     const Agenda& gas_scattering_agenda,
                     const Numeric& rte_alonglos_v);

/** Calculates the ppath towards the suns from a given position and indicates
 *  if sun is visible or not.
 *
 * @param[in,out] ws ARTS workspace.
 * @param[out] sun_ppaths ArrayOfPpath Propagation path towards each sun.
 * @param[out] suns_visible ArrayOfIndex Flag indicating if eah sun is visible.
 * @param[out] sun_rte_los ArrayOfVector Incoming direction of the each
 *             transmitted sun radiation.
 * @param[in] rte_pos As the WSV.
 * @param[in] suns As the WSV.
 * @param[in] f_grid As the WSV.
 * @param[in] p_grid As the WSV.
 * @param[in] lat_grid As the WSV.
 * @param[in] lon_grid As the WSV.
 * @param[in] z_field As the WSV.
 * @param[in] z_surface As the WSV.
 * @param[in] ppath_lmax As the WSV.
 * @param[in] ppath_lraytrace As the WSV.
 * @param[in] ppath_step_agenda As the WSV.
 */
void get_sun_ppaths(const Workspace& ws,
                     ArrayOfPpath& sun_ppaths,
                     ArrayOfIndex& suns_visible,
                     ArrayOfVector& sun_rte_los,
                     const Vector& rte_pos,
                     const ArrayOfSun& suns,
                     const Vector& f_grid,
                     const AtmField& atm_field,
                     const SurfaceField& surface_field,
                     const Numeric& ppath_lmax,
                     const Numeric& ppath_lraytrace,
                     const Agenda& ppath_step_agenda);
 
 /** regrid_sun_spectrum
 *
 * Regrids a given spectrum from a griddedfield2 to the f_grid.
 * if the f_grid covers a larger range as the given one, one
 * can choose between two padding options:
 * zeros: Intensities outside the given spectrum are set to zero 
 * planck: Intensities outside the given spectrum are initilizied 
 *        with the black body value at that frequency.
 *
 * @param[in]  sun_spectrum_raw  gf2 of the given spectrum.
 * @param[in]  f_grid  f_grid for the calculation.
 * @param[in]  temperature  Temperature for the planck padding.
 *
 * @return     interpolated spectrum
 *
 * @author Jon Petersen
 * @date   2022-01-19
 */
Matrix regrid_sun_spectrum(const GriddedField2& sun_spectrum_raw,
                          const Vector &f_grid,
                          const Numeric &temperature);

#endif /* star_h */
