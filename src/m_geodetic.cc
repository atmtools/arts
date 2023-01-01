/* Copyright (C) 2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
                            
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

/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_geodetic.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2012-02-06

  \brief  Workspace functions of geodetic character.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "arts_conversions.h"
#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "matpack_data.h"
#include "messages.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidForAzimuth(Vector& refellipsoid,
                            const Numeric& latitude,
                            const Numeric& azimuth,
                            const Verbosity&) {
  ARTS_USER_ERROR_IF (refellipsoid.nelem() != 2,
                      "Input *refellispoid must be a vector of length 2*.");

  if (refellipsoid[1] > 0) {
    const Numeric e2 = refellipsoid[1] * refellipsoid[1];
    const Numeric a = 1 - e2 * pow(sin(DEG2RAD * latitude), 2.0);

    const Numeric rn = 1 / sqrt(a);
    const Numeric rm = (1 - e2) * (rn / a);

    const Numeric v = DEG2RAD * azimuth;

    refellipsoid[0] =
        refellipsoid[0] / (pow(cos(v), 2.0) / rm + pow(sin(v), 2.0) / rn);
    refellipsoid[1] = 0;
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void refellipsoidOrbitPlane(Vector& refellipsoid,
                            const Numeric& orbitinc,
                            const Verbosity&) {
  ARTS_USER_ERROR_IF (refellipsoid.nelem() != 2,
                      "Input *refellispoid must be a vector of length 2*.");
  chk_if_in_range("orbitinc", orbitinc, 0, 180);

  // Radius at maximum latitude
  const Numeric rp = refell2r(refellipsoid, orbitinc);

  // New eccentricity
  refellipsoid[1] = sqrt(1 - pow(rp / refellipsoid[0], 2.0));
}
