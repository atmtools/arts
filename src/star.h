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
  === Functions in star.cc
  ===========================================================================*/

/** Calculates the radiance spectrum of star which is scattered by the atmospheric gases.
 *
 * @param[in,out] ws ARTS workspace.
 * @param[out] scattered_starlight StokesVector scattered monochromatic radiance
 *             spectrum of star.
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
    StokesVector& scattered_starlight,
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
 * @param[in, out] iy Matrix of star background
 * @param[in] star_pos Position of star in geographical coordinate
 * @param[in] star_radius Radius of star
 * @param[in] star_spectrum Monochromatic irradiance
 *             spectrum of star
 * @param[in] rte_pos The position of the ppath point.
 * @param[in] rte_los The line of sight of the ppath.
 * @param[in] refellipsoid As the WSV with the same name.
  */
void get_star_background(Matrix& iy,
                         const Vector& star_pos,
                         const Numeric& star_radius,
                         const Matrix& star_spectrum,
                         const Vector& rte_pos,
                         const Vector& rte_los,
                         const Vector& refellipsoid,
                         const Index& atmosphere_dim);



#endif /* star_h */
