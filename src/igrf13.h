#ifndef igrf13_h
#define igrf13_h

#include "artstime.h"
#include "matpack_constexpr.h"
#include "matpack_data.h"

namespace IGRF {
//! Magnetic field for the east (u), north (v), and up (w) components as the ENU-coordinate system
struct MagneticField {
  Tensor3 u;
  Tensor3 v;
  Tensor3 w;
  
  //! Explicitly set the values to zero
  MagneticField(Index p, Index r, Index c) noexcept :
  u(p, r, c, 0), v(p, r, c, 0), w(p, r, c, 0) {}
};

/** Computes the magnetic field based on IGRF13 coefficients
 * 
 * Only 2000-2020 is implemented.  Any time before uses pure 2000 data, 
 * and any time after uses pure 2020 data.
 * 
 * WARNING:  No conversion of ENU to geodetic equivalents are performed
 * Instead the assumption is that the spherical ENU is good enough.
 * 
 * @param[in] z_field As WSV
 * @param[in] lat_grid As WSV
 * @param[in] lon_grid As WSV
 * @param[in] time A time
 * @param[in] ell The ellipsoid
 * @return The magnetic field in ENU as described by the MagneticField struct
 */
MagneticField compute(const Tensor3& z_field, const Vector& lat_grid, const Vector& lon_grid, const Time& time, const Vector2 ell);
} // namespace IGRF

#endif  // igrf13_h
