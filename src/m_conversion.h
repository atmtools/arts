/* Copyright (C) 2010 Claudia Emde <claudia.emde@lmu.de>
 
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

extern const Numeric SPEED_OF_LIGHT;

/* Workspace method: Doxygen documentation will be auto-generated */
void FrequencyFromWavelength(// WS Generic Output
                             Numeric& frequency,
                             // WS Generic Input
                             const Numeric& wavelength,
                             const Verbosity&)
{
 
  // Convert from wavelength to frequency
  frequency=SPEED_OF_LIGHT/wavelength;

}


/* Workspace method: Doxygen documentation will be auto-generated */
void FrequencyFromWavelength(// WS Generic Output
                             Vector& frequency,
                             // WS Generic Input
                             const Vector& wavelength,
                             const Verbosity&)
{
  frequency.resize(wavelength.nelem());
  // Convert from wavelength to frequency
  for (Index i=0; i<wavelength.nelem(); i++)
    frequency[i]=SPEED_OF_LIGHT/wavelength[i];
}

#endif /* m_conversion_h */

