/* Copyright (C) 2021
   Jon Petersen <jon.petersen@studium.uni-hamburg.de>
   Manfred Brath  <manfred.brath@uni-hamburg.de>

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
  \file   star.h
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
#include "gridded_fields.h"
#include "matpack.h"
#include "transmissionmatrix.h"
#include "energylevelmap.h"
#include "optproperties.h"


class Agenda;
class Workspace;

/*===========================================================================
  === structs/classes  in star.h
  ===========================================================================*/

/** The structure to describe a propagation path and releated quantities.
 *
 *  The fields of the structure are described more in detail inside the ARTS
 *  user guide (AUG).
 */
struct Star {
  /** star description */
  String description;
  /** Star spectrum, monochrmatic radiance spectrum at the surface of the star*/
  Matrix spectrum;
  /** Star radius */
  Numeric radius;
  /** star distance from center of planet to center of star*/
  Numeric distance;
  /** latitude of the star in the sky of the planet */
  Numeric latitude;
  /** longitude of the star in the sky of the planet */
  Numeric longitude;

  friend std::ostream& operator<<(std::ostream& os, const Star& star);
};



/** An array of star. */
using ArrayOfStar = Array<Star>;


/*===========================================================================
  === Functions in star.h
  ===========================================================================*/

/** Calculates the radiance spectrum of star which is scattered by the atmospheric gases.
 *
 * @param[in,out] ws ARTS workspace.
 * @param[out] scattered_starlight RadiationVector scattered monochromatic radiance
 *             spectrum of star.
 * @param[in] f_grid Vector frequency grid.
 * @param[in] p Numeric pressure at location of scattering.
 * @param[in] T Numeric temperature at location of scattering.
 * @param[in] vmr Vector volume mixing ratios of absorption species at location
 *            of scattering.
 * @param[in] transmitted_starlight Matrix transmitted monochromatic irradiance
 *             spectrum of star at location of scattering.
 * @param[in] gas_scattering_los_in Vector incoming direction of the transmitted star irradiance
 *            spectrum.
 * @param[in] gas_scattering_los_out outgoing direction of the transmitted star irradiance
 *            spectrum.
 * @param[in] gas_scattering_agenda Agenda agenda calculating the gas scattering
 *            cross sectionand matrix.
 */
void get_scattered_starsource(Workspace& ws,
                              RadiationVector& scattered_starlight,
                              const Vector& f_grid,
                              const Numeric& p,
                              const Numeric& T,
                              const Vector& vmr,
                              const Matrix& transmitted_starlight,
                              const Vector& gas_scattering_los_in,
                              const Vector& gas_scattering_los_out,
                              const Agenda& gas_scattering_agenda);

/** Gets the star background for a given ppath.
 *
 * iy is zero if there is no star in the line of sight at TOA.
 *
 * @param[out] iy Matrix radiance spectrum of stars.
 * @param[out] stars_visible Index indicating if stars are in los
 * @param[in] star Star-structure.
 * @param[in] ppath Propagation path as the WSV.
 * @param[in] f_grid Vector as the WSV.
 * @param[in] stokes_dim Index as the WSV.
 * @param[in] atmosphere_dim Index as the WSV.
 * @param[in] refellipsoid Vector as the WSV.
 */
void get_star_background(Matrix& iy,
                         Index& stars_visible,
                         const ArrayOfStar& stars,
                         const Ppath& ppath,
                         const Vector& f_grid,
                         const Index& stokes_dim,
                         const Index& atmosphere_dim,
                         const Vector& refellipsoid);

/** Checks and adds star radiance if star is in line of sight.
 *
 * @param[in, out] iy Matrix of star.
 * @param[out] stars_visible Index indicating if star are in los
 * @param[in] star Star-structure.
 * @param[in] rtp_pos The position of the ppath point.
 * @param[in] rtp_los The line of sight of the ppath.
 * @param[in] refellipsoid As the WSV with the same name.
  */
void get_star_radiation(Matrix& iy,
                        Index& stars_visible,
                         const Star& star,
                         const Vector& rtp_pos,
                         const Vector& rtp_los,
                         const Vector& refellipsoid);

/** Calculates the transmitted star radiation at the end position of the ppath
 *
 * @param[in,out] ws ARTS workspace.
 * @param[out] direct_radiation Matrix Transmitted monochromatic irradiance
 *             spectrum of star.
 * @param[out] ddirect_radiation_dx Array of Tensor3 Jacobian of transmitted
 *              monochromatic irradiance spectrum of star.
 * @param[in] stokes_dim As the WSV.
 * @param[in] f_grid As the WSV.
 * @param[in] atmosphere_dim As the WSV.
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
 * @param[in] star_ppaths ArrayOfPpath Propagation path towards each star.
 * @param[in] stars As the WSV.
 * @param[in] stars_visible ArrayOfIndex Flag indicating if eah star is visible.
 * @param[in] refellipsoid As the WSV.
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
 * @param[in] verbosity Verbosity.
 */
void get_direct_radiation(Workspace& ws,
                     ArrayOfMatrix& direct_radiation,
                     ArrayOfArrayOfTensor3& ddirect_radiation_dx,
                     const Index& stokes_dim,
                     const Vector& f_grid,
                     const Index& atmosphere_dim,
                     const Vector& p_grid,
                     const Vector& lat_grid,
                     const Vector& lon_grid,
                     const Tensor3& t_field,
                     const EnergyLevelMap& nlte_field,
                     const Tensor4& vmr_field,
                     const ArrayOfArrayOfSpeciesTag& abs_species,
                     const Tensor3& wind_u_field,
                     const Tensor3& wind_v_field,
                     const Tensor3& wind_w_field,
                     const Tensor3& mag_u_field,
                     const Tensor3& mag_v_field,
                     const Tensor3& mag_w_field,
                     const Index& cloudbox_on,
                     const ArrayOfIndex& cloudbox_limits,
                     const Index& gas_scattering_do,
                     const Index& irradiance_flag,
                     const ArrayOfPpath& star_ppaths,
                     const ArrayOfStar& stars,
                     const ArrayOfIndex& stars_visible,
                     const Vector& refellipsoid,
                     const Tensor4& pnd_field,
                     const ArrayOfTensor4& dpnd_field_dx,
                     const ArrayOfString& scat_species,
                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                     const Index& jacobian_do,
                     const ArrayOfRetrievalQuantity& jacobian_quantities,
                     const Agenda& propmat_clearsky_agenda,
                     const Agenda& water_p_eq_agenda,
                     const Agenda& gas_scattering_agenda,
                     const Numeric& rte_alonglos_v,
                     const Verbosity& verbosity);

/** Calculates the ppath towards the stars from a given position and indicates
 *  if star is visible or not.
 *
 * @param[in,out] ws ARTS workspace.
 * @param[out] star_ppaths ArrayOfPpath Propagation path towards each star.
 * @param[out] stars_visible ArrayOfIndex Flag indicating if eah star is visible.
 * @param[out] star_rte_los ArrayOfVector Incoming direction of the each
 *             transmitted star radiation.
 * @param[in] rte_pos As the WSV.
 * @param[in] stars As the WSV.
 * @param[in] f_grid As the WSV.
 * @param[in] atmosphere_dim As the WSV.
 * @param[in] p_grid As the WSV.
 * @param[in] lat_grid As the WSV.
 * @param[in] lon_grid As the WSV.
 * @param[in] z_field As the WSV.
 * @param[in] z_surface As the WSV.
 * @param[in] refellipsoid As the WSV.
 * @param[in] ppath_lmax As the WSV.
 * @param[in] ppath_lraytrace As the WSV.
 * @param[in] ppath_step_agenda As the WSV.
 * @param[in] verbosity Verbosity.
 */
void get_star_ppaths(Workspace& ws,
                     ArrayOfPpath& star_ppaths,
                     ArrayOfIndex& stars_visible,
                     ArrayOfVector& star_rte_los,
                     const Vector& rte_pos,
                     const ArrayOfStar& stars,
                     const Vector& f_grid,
                     const Index& atmosphere_dim,
                     const Vector& p_grid,
                     const Vector& lat_grid,
                     const Vector& lon_grid,
                     const Tensor3& z_field,
                     const Matrix& z_surface,
                     const Vector& refellipsoid,
                     const Numeric& ppath_lmax,
                     const Numeric& ppath_lraytrace,
                     const Agenda& ppath_step_agenda,
                     const Verbosity& verbosity);
 
 /** regrid_star_spectrum
 *
 * Regrids a given spectrum from a griddedfield2 to the f_grid.
 * if the f_grid covers a larger range as the given one, one
 * can choose between two padding options:
 * zeros: Intensities outside the given spectrum are set to zero 
 * planck: Intensities outside the given spectrum are initilizied 
 *        with the black body value at that frequency.
 *
 * @param[in]  star_spectrum_raw  gf2 of the given spectrum.
 * @param[in]  f_grid  f_grid for the calculation.
 * @param[in]  stokes_dim  stokes_dim for the calculation.
 * @param[in]  temperature  Temperature for the planck padding.
 *
 * @return     interpolated spectrum
 *
 * @author Jon Petersen
 * @date   2022-01-19
 */
Matrix regrid_star_spectrum(const GriddedField2& star_spectrum_raw,
                          const Vector &f_grid,
                          const Index &stokes_dim,
                          const Numeric &temperature,
                          const Verbosity &verbosity);

#endif /* star_h */
