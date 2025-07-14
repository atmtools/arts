#pragma once

#include <artstime.h>
#include <matpack.h>

namespace IGRF {
/** Computes the magnetic field based on IGRF13 coefficients
 * 
 * Only 2000-2020 is implemented.  Any time before uses pure 2000 data, 
 * and any time after uses pure 2020 data.
 * 
 * WARNING:  No conversion of ENU to geodetic equivalents are performed
 * Instead the assumption is that the spherical ENU is good enough.
 * 
 * @param[in] pos The position in [alt, lat, lon] (geodetic)
 * @param[in] ell The ellipsoid (a, b)
 * @param[in] time A time stamp
 * @return The magnetic field in ENU as described by the MagneticField struct
 */
Vector3 igrf(const Vector3 pos, const Vector2 ell, const Time& time=Time{});
}  // namespace IGRF
