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

#include "agenda_class.h"
#include "transmissionmatrix.h"

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

/** Checks and adds star radiance to background if star is in line of sight.
 *
 * @param[in, out] iy Matrix of star background.
 * @param[in] star Star-structure.
 * @param[in] rtp_pos The position of the ppath point.
 * @param[in] rtp_los The line of sight of the ppath.
 * @param[in] refellipsoid As the WSV with the same name.
  */
void get_star_background(Matrix& iy,
                         const Star& stars,
                         const Vector& rtp_pos,
                         const Vector& rtp_los,
                         const Vector& refellipsoid);



#endif /* star_h */
