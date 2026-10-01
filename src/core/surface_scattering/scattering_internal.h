#pragma once

#include <arts_constants.h>
#include <configtypes.h>
#include <debug.h>
#include <enumsInterpolationExtrapolation.h>
#include <limits>

#include <lagrange_interp.h>
#include <cmath>

namespace surface_scattering {

//! Longitude cycler for [-180, 180] range
using lon_cycler = lagrange_interp::loncross;  // cycler<-180.0, 180.0>

/** Validate that longitude grid is within [-180, 180] range.
 *
 * @param lon_grid The longitude grid to validate
 * @throw ARTS_USER_ERROR if any grid value is outside [-180, 180]
 */
inline void validate_longitude_grid(const Vector& lon_grid) {
  for (const auto lon : lon_grid) {
    ARTS_USER_ERROR_IF(
        lon < -180.0 || lon > 180.0,
        "Longitude grid value ",
        lon,
        " is outside the supported range [-180, 180]. "
        "Only longitude grids in the [-180, 180] convention are supported.");
  }
}

/** Map InterpolationExtrapolation to frequency extrapolation limit.
 *
 * Determines how far the interpolation can extrapolate beyond grid boundaries:
 * - Linear: unlimited extrapolation (max limit)
 * - Nearest/None/Zero: clamp at boundaries (0.0 limit)
 *
 * @param extrap The extrapolation mode
 * @return Numeric extrapolation limit for lagrange_interp::lag()
 */
inline Numeric frequency_extrap_limit(InterpolationExtrapolation extrap) {
  switch (extrap) {
    case InterpolationExtrapolation::Linear:
      return std::numeric_limits<Numeric>::max();
    case InterpolationExtrapolation::Nearest:
    case InterpolationExtrapolation::None:
    case InterpolationExtrapolation::Zero:
      return 0.0;
  }
  std::unreachable();
}

}  // namespace surface_scattering
