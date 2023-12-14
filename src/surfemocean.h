/*!
  \file   surfemocean_h
  \author Shaofei Wang <wangshaofei224@163.com>
  \date   Tue Dec 12 16:48:31 2023
  
  \brief  This file contains functions that are adapted from SURFEM-Ocean 
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
