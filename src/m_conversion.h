/* Copyright (C) 2012 Claudia Emde <claudia.emde@lmu.de>
 
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

/*!
  \file   m_conversion.h
  \author Claudia Emde <claudia.emde@lmu.de>
  \date   2010-07-21
  
  \brief  Implementation of unit conversion functions
  
*/

#ifndef m_conversion_h
#define m_conversion_h

#include "arts_constants.h"
#include "matpack_data.h"

inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;
inline constexpr Numeric PI=Constant::pi;

/* Workspace method: Doxygen documentation will be auto-generated */
inline void FrequencyFromWavelength(  // WS Generic Output
    Numeric& frequency,
    // WS Generic Input
    const Numeric& wavelength) {
  // Convert from wavelength to frequency
  frequency = SPEED_OF_LIGHT / wavelength;
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void FrequencyFromWavelength(  // WS Generic Output
    Vector& frequency,
    // WS Generic Input
    const Vector& wavelength) {
  frequency.resize(wavelength.nelem());
  // Convert from wavelength to frequency
  for (Index i = 0; i < wavelength.nelem(); i++)
    frequency[i] = SPEED_OF_LIGHT / wavelength[i];
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void FrequencyFromCGSAngularWavenumber(  // WS Generic Output
    Numeric& frequency,
    // WS Generic Input
    const Numeric& angular_wavenumber) {
  frequency = SPEED_OF_LIGHT * angular_wavenumber / (2 * PI) * 100;
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void FrequencyFromCGSAngularWavenumber(  // WS Generic Output
    Vector& frequency,
    // WS Generic Input
    const Vector& angular_wavenumber) {
  frequency.resize(angular_wavenumber.nelem());
  // Convert from angular wavenumber to frequency
  for (Index i = 0; i < angular_wavenumber.nelem(); i++)
    frequency[i] = SPEED_OF_LIGHT * angular_wavenumber[i] / (2 * PI) * 100;
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void FrequencyFromCGSKayserWavenumber(  // WS Generic Output
    Numeric& frequency,
    // WS Generic Input
    const Numeric& kayser_wavenumber) {
  frequency = SPEED_OF_LIGHT * kayser_wavenumber * 100;
}

/* Workspace method: Doxygen documentation will be auto-generated */
inline void FrequencyFromCGSKayserWavenumber(  // WS Generic Output
    Vector& frequency,
    // WS Generic Input
    const Vector& kayser_wavenumber) {
  frequency.resize(kayser_wavenumber.nelem());
  // Convert from Kayser wavenumber to frequency
  for (Index i = 0; i < kayser_wavenumber.nelem(); i++)
    frequency[i] = SPEED_OF_LIGHT * kayser_wavenumber[i] * 100;
}

#endif /* m_conversion_h */
