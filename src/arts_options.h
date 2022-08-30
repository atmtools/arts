/* Copyright (C) 2019
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/*!
 * \file   arts_options.h
 * \brief  Options for ARTS from enumeration (including error handling)
 * 
 * \author Richard Larsson
 * \date   2019-04-01
 */

#ifndef OPTIONS_IN_ARTS_H
#define OPTIONS_IN_ARTS_H

#include "enums.h"

namespace Options {
/** Keep time options available to switch over them */
ENUMCLASS(
    TimeStep, char, hour, hours, h, minute, minutes, min, second, seconds, s)

/** Possible line shape coefficients */
ENUMCLASS(LineShapeCoeff, char, X0, X1, X2, X3)

/** Possible Jacobian for Wind and Magnetic field */
ENUMCLASS(WindMagJacobian, char, u, v, w, strength)

/** Possible Jacobian for basic line parameters */
ENUMCLASS(BasicCatParamJacobian, char, LineStrength, LineCenter)

/** Possible HITRAN types of catalogs */
ENUMCLASS(
    HitranType,
    char,
    Pre2004,   // 2004 version changed the .par-length
    Post2004,  // New par length
    Online  // Onine expects a modern .par line followed by Upper then Lower quantum numbers
)

/** Possible AddLines Speedups */
ENUMCLASS(LblSpeedup, char, None, QuadraticIndependent, LinearIndependent)

ENUMCLASS(SortingOption, char, ByFrequency, ByEinstein)

/** Options for setting iy_main_agenda */
ENUMCLASS(iy_main_agendaDefaultOptions,
          char,
          Emission,
          Clearsky,
          Transmission,
          TransmissionUnitUnpolIntensity,
          TransmissionUnitPolIntensity,
          Freqloop,
          ScattMC)
}  // namespace Options

#endif
