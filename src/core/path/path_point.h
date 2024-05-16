#pragma once

#include <atm.h>
#include <enums.h>
#include <matpack.h>
#include <surf.h>

#include "matpack_constexpr.h"

namespace path {

//! A simple path-point of a propagation path
struct PropagationPathPoint {
  /*! The position type and line-of-sight type
   *
   * Note that when in a list of points, the position types and line-of-sight
   * types should match so that the previous point's los_type is the same as
   * the current point's pos_type, and that the current point's los_type is the
   * same as the next point's pos_type.
   */
  PathPositionType pos_type{PathPositionType::unknown}, los_type{PathPositionType::unknown};

  //! Position of the point: alt [m], lat [deg], lon [deg]
  Vector3 pos;

  //! Line-of-sight of the point: zenith [deg], azimuth [deg]
  Vector2 los;

  //! The real part of the refractive index at each path position
  Numeric nreal{1};

  //! The group index of refraction
  Numeric ngroup{1};

  [[nodiscard]] constexpr bool has(PathPositionType x) const noexcept {
    return pos_type == x or los_type == x;
  }

  friend std::ostream& operator<<(std::ostream& os,
                                  const PropagationPathPoint& p);

  [[nodiscard]] constexpr Numeric altitude() const noexcept { return pos[0]; }
  [[nodiscard]] constexpr Numeric& altitude() noexcept { return pos[0]; }
  [[nodiscard]] constexpr Numeric latitude() const noexcept { return pos[1]; }
  [[nodiscard]] constexpr Numeric& latitude() noexcept { return pos[1]; }
  [[nodiscard]] constexpr Numeric longitude() const noexcept { return pos[2]; }
  [[nodiscard]] constexpr Numeric& longitude() noexcept { return pos[2]; }
  [[nodiscard]] constexpr Numeric zenith() const noexcept { return los[0]; }
  [[nodiscard]] constexpr Numeric& zenith() noexcept { return los[0]; }
  [[nodiscard]] constexpr Numeric azimuth() const noexcept { return los[1]; }
  [[nodiscard]] constexpr Numeric& azimuth() noexcept { return los[1]; }
};

using ArrayOfPropagationPathPoint = std::vector<PropagationPathPoint>;
using ArrayOfArrayOfPropagationPathPoint = std::vector<ArrayOfPropagationPathPoint>;

std::ostream& operator<<(std::ostream& os,
                         const ArrayOfPropagationPathPoint& p);

std::ostream& operator<<(std::ostream& os,
                         const ArrayOfArrayOfPropagationPathPoint& p);

/** Mirror the line-of-sight
 * 
 * @param los A line-of-sight [zenith, azimuth] in degrees
 * @return The mirrored line-of-sight [zenith, azimuth] in degrees
 */
Vector2 mirror(const Vector2 los);

/** Initializes a propagation path point
 *
 * The point is initialized at the sensor or start of radiation as determined 
 * by the as_sensor flag.  If the path is initialized as a sensor
 * the line-of-sight is mirrored to fit the definition of a PropagationPathPoint.
 * 
 * @param pos The sensor or radiation position of the path
 * @param los The sensor or radiation line-of-sight
 * @param atm_field The atmospheric field (as the WSV) 
 * @param surface_field The surface field (as the WSV)
 * @param as_sensor Treat as sensor flag
 * @return A path point that may initialize a path
 */
PropagationPathPoint init(const Vector3& pos,
                          const Vector2& los,
                          const AtmField& atm_field,
                          const SurfaceField& surface_field,
                          bool as_sensor = true);

/** Set the geometric extremes object
 * 
 * Geometric extremes are intersections with the surface and the top-of-the-atmosphere.
 *
 * @param path The propagation path
 * @param atm_field The atmospheric field (as the WSV)
 * @param surface_field The surface field (as the WSV)
 * @param surface_search_accuracy The surface search accuracy in meters
 * @param surface_search_safe Flag to search the surface safely or fast
 * @return The input for piping
 */
ArrayOfPropagationPathPoint& set_geometric_extremes(
    ArrayOfPropagationPathPoint& path,
    const AtmField& atm_field,
    const SurfaceField& surface_field,
    const Numeric surface_search_accuracy = 0.1,
    const bool surface_search_safe = false);

/** Fills a propagation path with geometrically spaced points
 *
 * Only looks at propagation path point pairs that are both in the atmosphere
 * according to PropagationPathPoint::has(PathPositionType::atm)
 *
 * From the first point, a new point is added to the path by moving a distance
 * of max_step along the line-of-sight untill the last point is less than
 * max_step away from the end of the path.
 * 
 * @param path The propagation path
 * @param surface_field The surface field (as the WSV)
 * @param max_step The maximum step size in meters in path after this method completes
 * @return The input for piping
 */
ArrayOfPropagationPathPoint& fill_geometric_stepwise(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surface_field,
    const Numeric max_step);

/** Adds all altitude grid crossings to a propagation path
 * 
 * @param path The propagation path
 * @param surface_field The surface field (as the WSV)
 * @param alt_grid The altitude grid
 * @return The input for piping
 */
ArrayOfPropagationPathPoint& fill_geometric_altitude_crossings(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surface_field,
    const Vector& alt_grid);

/** Adds all latitude grid crossings to a propagation path
 * 
 * @param path The propagation path
 * @param surface_field The surface field (as the WSV)
 * @param lat_grid The latitude grid
 * @return The input for piping
 */
ArrayOfPropagationPathPoint& fill_geometric_latitude_crossings(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surface_field,
    const Vector& lat_grid);

/** Adds all longitude grid crossings to a propagation path
 * 
 * @param path The propagation path
 * @param surface_field The surface field (as the WSV)
 * @param lon_grid The longitude grid
 * @return The input for piping
 */
ArrayOfPropagationPathPoint& fill_geometric_longitude_crossings(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surface_field,
    const Vector& lon_grid);

/** Finds the geometric limb of a propagation path
 * 
 * @param path The propagation path
 * @param surface_field The surface field (as the WSV)
 * @return The geometric limb as a path point 
 */
PropagationPathPoint find_geometric_limb(
    const ArrayOfPropagationPathPoint& path, const SurfaceField& surface_field);

/** Fills the geometric limb of a propagation path
 *
 * Uses find_geometric_limb to find the limb and then inserts this to the path
 * as the lowest altitude point.
 *
 * @param path The propagation path
 * @param surface_field The surface field (as the WSV)
 * @return The input for piping
 */
ArrayOfPropagationPathPoint& fill_geometric_limb(
    ArrayOfPropagationPathPoint& path, const SurfaceField& surface_field);

/** Erases all points that are closer than min_dist to the surface
 *
 * @param path The propagation path
 * @param surface_field The surface field (as the WSV)
 * @param min_dist The minimum distance between two path points
 * @param first If true, the first path point is erased, otherwise the last
 * @return The input for piping
 */
ArrayOfPropagationPathPoint& erase_closeby(ArrayOfPropagationPathPoint& path,
                                           const SurfaceField& surface_field,
                                           const Numeric min_dist,
                                           const bool first = true);

/** Sums up the geometric path length of a propagation path
 *
 * Only uses the first and last atmospheric path points as defined by
 * PropagationPathPoint::has(PathPositionType::atm)
 * 
 * @param path The propagation path
 * @param surface_field The surface field (as the WSV)
 * @return Numeric Distance in meters
 */
Numeric total_geometric_path_length(const ArrayOfPropagationPathPoint& path,
                                    const SurfaceField& surface_field);

/** Distance between two positions in in meters
 *
 * The positions are given as [alt, lat, lon] in meters and degrees
 * 
 * @param pos1 The first position
 * @param pos2 The second position
 * @param ellipsoid The ellipsoid of the body
 * @return Numeric Distance in meters
 */
Numeric distance(const Vector3 pos1,
                 const Vector3 pos2,
                 const Vector2 ellipsoid);

/** Finds the zenith angle of the tangent limb at a given altitude as viewed by
 * a sensor at a given position
 *
 * Note that the tangent altitude must be strictly below the sensor altitude.
 *
 * The algorithm uses a binary-reduction approach to hone into the closest 
 * solution that can be represented by a Numeric
 *
 * @param[in] pos [alt, lat, lon] of the sensor
 * @param[in] surface_field The surface field (as the WSV)
 * @param[in] alt The tangent altitude
 * @param[in] azimuth The azimuth of the sensor
 * @return The zenith angle
 */
Numeric geometric_tangent_zenith(const Vector3 pos,
                                 const SurfaceField& surface_field,
                                 const Numeric alt,
                                 const Numeric azimuth = 0);

/*! Remove all propagation path points that are not looking at or are in the atmosphere
 *
 * The PropagationPathPoint::has function is used to determine if a point is
 * looking at or is in the atmosphere.
 *
 * @param[in] path The propagation path
 * @return The input for piping
 */
ArrayOfPropagationPathPoint& keep_only_atm(ArrayOfPropagationPathPoint& path);

/** Fix the up-down azimuth of a propagation path to the given value
 *
 * This compensates for the atan2 rules that:
 *
 *  If y is ±0 and x is negative or -0, ±π is returned.
 *  If y is ±0 and x is positive or +0, ±0 is returned.
 *  If x is ±0 and y is negative, -π/2 is returned.
 *  If x is ±0 and y is positive, +π/2 is returned. 
 * 
 * Since the line-of-sight ECEF vector may randomly choose between +0 and -0
 * for the output of a los vector such as [180, 0] that has been propagated
 * forward some distance, the azimuth of an uncompensated path may randomly
 * switch between an azimuth angle of -90, 0, and 90 degrees.
 *
 * This method checks that if all of the angles are 0, 90, or -90,
 * then all azimuth angles are set to the first path point azimuth input.
 * If this condition is false, the path is left unchanged.
 * 
 * @param[in] path The propagation path
 * @return The input for piping
 */
ArrayOfPropagationPathPoint& fix_updown_azimuth_to_first(
    ArrayOfPropagationPathPoint& path);
}  // namespace path

using PropagationPathPoint = path::PropagationPathPoint;
using ArrayOfPropagationPathPoint = path::ArrayOfPropagationPathPoint;
using ArrayOfArrayOfPropagationPathPoint = path::ArrayOfArrayOfPropagationPathPoint;
