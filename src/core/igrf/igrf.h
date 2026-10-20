#pragma once

#include <artstime.h>
#include <auto_igrf_data.h>
#include <matpack.h>

namespace IGRF {
//! The reference radius of the IGRF model [m]
inline constexpr Numeric r0 = 6371.2e3;

//! The g or h coefficients [nT] of the IGRF model, the (n, m) coefficient is at [n, m]
using Coefficients = Data::Coefficients;

/** The decimal year of a time stamp, from its date and time in UTC
 *
 * This is the year plus the elapsed part of it, so "2001-07-02 12:00:00" is
 * 2001.5.  The IGRF model is linear in the decimal year.
 *
 * @param[in] time A time stamp
 * @return The decimal year
 */
Numeric decimal_year(const Time& time);

/** Get the IGRF coefficients at a time
 *
 * The coefficients are interpolated linearly in the decimal year between the
 * epochs of the model, and extrapolated from the last epoch with the secular
 * variation of the model.
 *
 * Throws if the decimal year of the time is outside of Data::valid_range, the
 * validity of the coefficients generated from mag-field/igrf.txt in arts-cat-data.
 *
 * @param[in] time A time stamp
 * @return The g and h coefficients [nT]
 */
std::pair<Coefficients, Coefficients> coefficients(const Time& time);

/** Computes the magnetic field of the IGRF model
 *
 * Throws if the decimal year of the time is outside of Data::valid_range.
 *
 * @param[in] pos The position in [alt, lat, lon] (geodetic)
 * @param[in] ell The ellipsoid (a, b)
 * @param[in] time A time stamp
 * @return The magnetic field [T] in geodetic east, north, and up
 */
Vector3 igrf(const Vector3 pos, const Vector2 ell, const Time& time = Time{});
}  // namespace IGRF
