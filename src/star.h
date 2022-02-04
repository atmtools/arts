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
#include "agenda_class.h"
#include "gridded_fields.h"
#include "matpack.h"
#include "transmissionmatrix.h"
#include "energylevelmap.h"

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
};

std::ostream& operator<<(std::ostream& os, const Star& star);

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
 * @param[out] dscattered_starlight ArrayOfRadiationVector derivatives of
 *              scattered monochromatic radiance.
 * @param[in] f_grid Vector frequency grid.
 * @param[in] p Numeric pressure at location of scattering.
 * @param[in] T Numeric temperature at location of scattering.
 * @param[in] vmr Vector volume mixing ratios of absorption species at location
 *            of scattering.
 * @param[in] transmitted_starlight Matrix transmitted monochromatic irradiance
 *             spectrum of star at location of scattering.
 * @param[in] in_los Vector incoming direction of the transmitted star irradiance
 *            spectrum.
 * @param[in] out_los outgoing direction of the transmitted star irradiance
 *            spectrum.
 * @param[in] gas_scattering_agenda Agenda agenda calculating the gas scattering
 *            cross sectionand matrix.
 */
void get_scattered_starsource(
    Workspace& ws,
    RadiationVector& scattered_starlight,
    ArrayOfRadiationVector& dscattered_starlight,
    const Vector& f_grid,
    const Numeric& p,
    const Numeric& T,
    const Vector& vmr,
    const Matrix& transmitted_starlight,
    const Vector& in_los,
    const Vector& out_los,
    const Agenda& gas_scattering_agenda
);

/** Gets the star background for a given ppath.
 *
 * iy is zero if there is no star in the line of sight at TOA.
 *
 * @param[out] iy Matrix radiance spectrum of stars.
 * @param[in] star Star-structure.
 * @param ppath Propagation path as the WSV.
 * @param atmosphere_dim Index as the WSV.
 * @param f_grid Vector as the WSV.
 * @param stokes_dim Index as the WSV.
 * @param refellipsoid Vector as the WSV.
 */
void get_star_background(Matrix& iy,
                         const ArrayOfStar& stars,
                         const Ppath& ppath,
                         const Index& atmosphere_dim,
                         const Vector& f_grid,
                         const Index& stokes_dim,
                         const Vector& refellipsoid);

/** Checks and adds star radiance if star is in line of sight.
 *
 * @param[in, out] iy Matrix of star.
 * @param[in] star Star-structure.
 * @param[in] rtp_pos The position of the ppath point.
 * @param[in] rtp_los The line of sight of the ppath.
 * @param[in] refellipsoid As the WSV with the same name.
  */
void get_star_radiation(Matrix& iy,
                         const Star& stars,
                         const Vector& rtp_pos,
                         const Vector& rtp_los,
                         const Vector& refellipsoid);

/** Calculates the transmitted radiation of one star for a given ppath position
 *
 * @param[in,out] ws ARTS workspace.
 * @param[out] transmitted_starlight Matrix Transmitted monochromatic irradiance
 *             spectrum of star.
 * @param[out] dtransmitted_starlight Array of Tensor3 Jacobian of transmitted
*              monochromatic irradiance spectrum of star.
 * @param[out] star_rte_los Vector Incoming direction of the transmitted star irradiance
 *            spectrum.
 * @param[in] rte_pos As the WSV.
 * @param[in] i_star Index of star.
 * @param[in] stokes_dim As the WSV.
 * @param[in] f_grid As the WSV.
 * @param[in] atmosphere_dim As the WSV.
 * @param[in] p_grid As the WSV.
 * @param[in] lat_grid As the WSV.
 * @param[in] lon_grid As the WSV.
 * @param[in] z_field As the WSV.
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
 * @param[in] z_surface As the WSV.
 * @param[in] refellipsoid As the WSV.
 * @param[in] ppath_lmax As the WSV.
 * @param[in] ppath_lraytrace As the WSV.
 * @param[in] gas_scattering_do As the WSV.
 * @param[in] jacobian_do As the WSV.
 * @param[in] jacobian_quantities As the WSV.
 * @param[in] stars As the WSV.
 * @param[in] rte_alonglos_v As the WSV.
 * @param[in] propmat_clearsky_agenda As the WSV.
 * @param[in] water_p_eq_agenda As the WSV.
 * @param[in] gas_scattering_agenda As the WSV.
 * @param[in] ppath_step_agenda As the WSV.
 * @param[in] verbosity Verbosity.
 */
void get_transmitted_starlight(
    Workspace& ws,
    Matrix& transmitted_starlight,
    ArrayOfTensor3 dtransmitted_starlight,
    Vector& star_rte_los,
    Index& star_path_ok,
    const Vector& rte_pos,
    const Index& i_star,
    const Index& stokes_dim,
    const Vector& f_grid,
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const Tensor3& z_field,
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
    const Matrix& z_surface,
    const Vector& refellipsoid,
    const Numeric& ppath_lmax,
    const Numeric& ppath_lraytrace,
    const Index& gas_scattering_do,
    const Index& jacobian_do,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfStar& stars,
    const Numeric& rte_alonglos_v,
    const Agenda& propmat_clearsky_agenda,
    const Agenda& water_p_eq_agenda,
    const Agenda& gas_scattering_agenda,
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
