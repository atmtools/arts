/*!
  \file   fastem.cc
  \author Sreerekha Ravi <rekha@sat.physik.uni-bremen.de>
  \date   Tue Aug 10 15:16:31 2004
  
  \brief  This file contains functions that are adapted from FASTEM 
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

#ifdef ENABLE_FASTEM
extern "C" {
#endif

void rttov_fastem5_(const Index& fastem_version,
                    const Numeric& frequency,
                    const Numeric& za,
                    const Numeric& temperature,
                    const Numeric& salinity,
                    const Numeric& wind_speed,
                    Numeric* emissivity,
                    Numeric* reflectivity,
                    const Numeric& transmittance,
                    const Numeric& rel_azimuth);

#ifdef ENABLE_FASTEM
}
#endif

// Define dummy function that throws a runtime error if ARTS is
// compiled without FASTEM support.
#ifndef ENABLE_FASTEM
void rttov_fastem5_(const Index&,
                    const Numeric&,
                    const Numeric&,
                    const Numeric&,
                    const Numeric&,
                    const Numeric&,
                    Numeric*,
                    Numeric*,
                    const Numeric&,
                    const Numeric&) {
  throw std::runtime_error(
      "This version of ARTS was compiled without FASTEM support.");
}

#endif

//! Calculate the surface emissivity using FASTEM
/*! 
  Calculate surface emissivity using the FASTEM model from RTTOV.

  This is a direct interface to the code from RTTOV. No checkls of input is
  made, to obtain this feature use FastemStandAlone that is also handling
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
  \param[in]  fastem_version  FASTEM version

  \author Oliver Lemke
  \date 2014-12-09
*/
void fastem(  // Output:
    Vector& emissivity,
    Vector& reflectivity,
    // Input:
    const Numeric frequency,
    const Numeric za,
    const Numeric temperature,
    const Numeric salinity,
    const Numeric wind_speed,
    const Numeric transmittance,
    const Numeric rel_azimuth,
    const Index fastem_version) {
  emissivity.resize(4);
  reflectivity.resize(4);

  rttov_fastem5_(fastem_version,
                 frequency / 1e9,
                 180 - za,
                 temperature,
                 salinity * 1e3,
                 wind_speed,
                 emissivity.data_handle(),
                 reflectivity.data_handle(),
                 transmittance,
                 rel_azimuth);
}
