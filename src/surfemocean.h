/*!
  \file   surfemocean_h
  \author Sreerekha Ravi <rekha@sat.physik.uni-bremen.de>
  \date   Tue Aug 11 18:09:31 2004
  
  \brief  This file contains functions that are adapted from FASTEM 
  code which is used to calculate surface emissivity.
*/

#ifndef surfemocean_h
#define surfemocean_h

#include "matpack_data.h"

void surfemocean(  // Output:
    Vector &emissivity,
    Vector &reflectivity,
    // Input:
    const Numeric frequency,
    const Numeric za,
    const Numeric temperature,
    const Numeric salinity,
    const Numeric wind_speed,
    const Numeric transmittance,
    const Numeric rel_azimuth);

#endif  //surfemocean_h
