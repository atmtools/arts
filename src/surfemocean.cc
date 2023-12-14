/*!
  \file   surfemocean.cc
  \author Shaofei Wang <wangshaofei224@163.com>
  \date   Tue Dec 12 16:48:31 2023
  
  \brief  This file contains functions that are adapted from SURFEM-Ocean 
  code which is used to calculate surface emissivity.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts_constants.h"
#include "arts_conversions.h"
#include "matpack_complex.h"
#include "exceptions.h"
#include "matpack_data.h"

using std::ostringstream;
using std::runtime_error;

inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);

#ifdef ENABLE_SURFEMOCEAN
extern "C" {
#endif

void surfem_ocean_(const Numeric& frequency,
                   const Numeric& za,
                   const Numeric& wind_speed,
                   const Numeric& temperature,
                   const Numeric& salinity,
                   const Numeric& rel_azimuth,
                   const Numeric& transmittance,
                   Numeric* emissivity,
                   Numeric* reflectivity);

#ifdef ENABLE_SURFEMOCEAN
}
#endif

// Define dummy function that throws a runtime error if ARTS is
// compiled without SURFEMOCEAN support.
#ifndef ENABLE_SURFEMOCEAN
void surfem_ocean_(const Numeric&,
                   const Numeric&,
                   const Numeric&,
                   const Numeric&,
                   const Numeric&,
                   const Numeric&,
                   const Numeric&,
                   Numeric*,
                   Numeric*) {
  throw std::runtime_error(
      "This version of ARTS was compiled without SURFEM-Ocean support.");
}

#endif

//! Calculate the surface emissivity using SURFEM-Ocean
/*! 
  Calculate surface emissivity using the SURFEM-Ocean model from RTTOV.

  This is a direct interface to the code from RTTOV. No checkls of input is
  made, to obtain this feature use SurfemOceanStandAlone that is also handling
  multiple frequencies.

  \param[out] emissivity      Calculated surface emissivity
  \param[out] reflectivity    Calculated surface reflectivity
  \param[in]  frequency       Frequency [Hz]
  \param[in]  za              Zenith angle of line-of-sigh 
  \param[in]  temperature     Temperature
  \param[in]  salinity        Salinity [0-1]
  \param[in]  wind_speed      Wind speed
  \param[in]  transmittance   Transmittance along downwelling direction.
  \param[in]  rel_azimuth     Relative azimuth angle (may not be used)

  \author Shaofei Wang
  \date 2023-12-12
*/
void surfemocean(  // Output:
    Vector& emissivity,
    Vector& reflectivity,
    // Input:
    const Numeric frequency,
    const Numeric za,
    const Numeric temperature,
    const Numeric salinity,
    const Numeric wind_speed,
    const Numeric transmittance,
    const Numeric rel_azimuth) {

  emissivity.resize(4);
  reflectivity.resize(4);

  surfem_ocean_(frequency / 1e9,
                180 - za,
                wind_speed,
                temperature,
                salinity * 1e3,
                rel_azimuth,
                transmittance,
                emissivity.data_handle(),
                reflectivity.data_handle());
}