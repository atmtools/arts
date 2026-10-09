#pragma once

#include <arts_constants.h>
#include <configtypes.h>
#include <debug.h>
#include <enumsInterpolationExtrapolation.h>
#include <limits>

#include <lagrange_interp.h>
#include <algorithm>
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
        "Longitude grid value {} is outside the supported range [-180, 180]. "
        "Only longitude grids in the [-180, 180] convention are supported.",
        lon);
  }
}

/** Minimum and maximum of a stored (sorted) frequency grid. */
inline std::pair<Numeric, Numeric> grid_range(const Vector& grid) {
  const auto [minp, maxp] = std::minmax_element(grid.begin(), grid.end());
  return {*minp, *maxp};
}

/** Prepare the query frequency grid for interpolation according to the
 * extrapolation mode.  The interpolation itself always runs unlimited; the
 * mode is enforced here and in extrap_postprocess():
 *
 * - Linear:  unlimited linear extrapolation (query grid unchanged)
 * - Nearest: query frequencies are clamped to the stored grid range, so
 *            outside values evaluate to the grid-edge value
 * - None:    any query frequency outside the stored grid range is a user error
 * - Zero:    query grid unchanged, values outside are zeroed afterwards
 *
 * @param f_grid The simulation frequency grid (query points)
 * @param gmin   Minimum of the stored frequency grid
 * @param gmax   Maximum of the stored frequency grid
 * @param extrap The extrapolation mode
 * @param what   Name of the data being interpolated (for error messages)
 * @return The (possibly clamped) query frequency grid
 */
inline Vector extrap_query_grid(const Vector&                    f_grid,
                                Numeric                          gmin,
                                Numeric                          gmax,
                                InterpolationExtrapolation       extrap,
                                const char*                      what) {
  if (extrap == InterpolationExtrapolation::None) {
    for (const auto f : f_grid) {
      ARTS_USER_ERROR_IF(f < gmin or f > gmax,
                         "Frequency {} is outside the stored frequency grid [{}, {}] of {}, "
                         "and interp_extrapolation is None.",
                         f,
                         gmin,
                         gmax,
                         what);
    }
    return f_grid;
  }

  if (extrap == InterpolationExtrapolation::Nearest) {
    Vector query(f_grid.size());
    for (Index i = 0; i < f_grid.size(); ++i) query[i] = std::clamp(f_grid[i], gmin, gmax);
    return query;
  }

  // Linear and Zero interpolate freely; Zero is applied in extrap_postprocess
  return f_grid;
}

/** Zero interpolated data at query frequencies outside the stored grid range
 * (Zero extrapolation mode); no-op for every other mode. */
inline void extrap_postprocess(Vector& data,
                               const Vector&                    f_grid,
                               Numeric                          gmin,
                               Numeric                          gmax,
                               InterpolationExtrapolation       extrap) {
  if (extrap != InterpolationExtrapolation::Zero) return;

  for (Index i = 0; i < data.size(); ++i) {
    if (f_grid[i] < gmin or f_grid[i] > gmax) data[i] = 0.0;
  }
}

}  // namespace surface_scattering
