/* Copyright (C) 2002-2012
   Sreerekha Ravi <rekha@sat.physik.uni-bremen.de>
   Stefan Buehler <sbuehler@ltu.se>                  

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
  \file   fastem.h
  \author Sreerekha Ravi <rekha@sat.physik.uni-bremen.de>
  \date   Tue Aug 11 18:09:31 2004
  
  \brief  This file contains functions that are adapted from FASTEM 
  code which is used to calculate surface emissivity.
*/

#ifndef fastem_h
#define fastem_h

#include "matpack_data.h"

void fastem(  // Output:
    Vector &emissivity,
    Vector &reflectivity,
    // Input:
    const Numeric frequency,
    const Numeric za,
    const Numeric temperature,
    const Numeric salinity,
    const Numeric wind_speed,
    const Numeric transmittance,
    const Numeric rel_azimuth,
    const Index fastem_version);

#endif  //fastem_h
