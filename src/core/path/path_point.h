#pragma once

#include <atm.h>
#include <enumsPathPositionType.h>
#include <matpack.h>
#include <surf.h>

#include <format>

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
  PathPositionType pos_type{PathPositionType::unknown};
  PathPositionType los_type{PathPositionType::unknown};

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
using ArrayOfArrayOfPropagationPathPoint =
    std::vector<ArrayOfPropagationPathPoint>;
using ArrayOfArrayOfArrayOfPropagationPathPoint =
    std::vector<ArrayOfArrayOfPropagationPathPoint>;

/** Mirror the line-of-sight
 * 
 * @param los A line-of-sight [zenith, azimuth] in degrees
 * @return The mirrored line-of-sight [zenith, azimuth] in degrees
 */
constexpr Vector2 mirror(const Vector2 los) {
  Vector2 los_mirrored;
  los_mirrored[0] = 180 - los[0];
  los_mirrored[1] = los[1] + 180;
  if (los_mirrored[1] > 180) los_mirrored[1] -= 360;
  return los_mirrored;
}

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
 * @param safe_search_accuracy The surface search accuracy in meters
 * @param search_safe Flag to search the surface safely or fast
 * @return The input for piping
 */
ArrayOfPropagationPathPoint& set_geometric_extremes(
    ArrayOfPropagationPathPoint& path,
    const AtmField& atm_field,
    const SurfaceField& surface_field,
    const Numeric surface_search_accuracy = 0.1,
    const bool surface_search_safe        = false);

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

/** Fills a propagation path with geometrically spaced points
 *
 * Only looks at propagation path point pairs that are both in the atmosphere
 * according to PropagationPathPoint::has(PathPositionType::atm)
 *
 * From the first point, a new point is added halfway between two existing points
 * if the distance between them is larger than max_step.
 * 
 * @param path The propagation path
 * @param surface_field The surface field (as the WSV)
 * @param max_step The maximum step size in meters in path after this method completes
 * @return The input for piping
 */
ArrayOfPropagationPathPoint& fill_geometric_by_half_steps(
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

/** Adds all grid crossings to a propagation path
 * 
 * @param path The propagation path
 * @param surface_field The surface field (as the WSV)
 * @param alt_grid The altitude grid
 * @param lat_grid The latitude grid
 * @param lon_grid The longitude grid
 * @return The input for piping
 */
ArrayOfPropagationPathPoint& fill_geometric_crossings(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surface_field,
    const Vector& alt_grid,
    const Vector& lat_grid,
    const Vector& lon_grid);

/** Finds the geometric limb of a propagation path
 * 
 * @param path The propagation path
 * @param surface_field The surface field (as the WSV)
 * @return The geometric limb as a path point 
 */
PropagationPathPoint find_geometric_limb(
    const ArrayOfPropagationPathPoint& path, const SurfaceField& surface_field);

/** Find the two intersections of a line with an ellipsoid at a given altitude
 * 
 * @param alt The altitude
 * @param ecef The position in ECEF coordinates
 * @param decef The line-of-sight in ECEF coordinates
 * @param ell The ellipsoid of the body [a, b]
 * @return std::pair<Numeric, Numeric> The two intersection points' distances
 */
std::pair<Numeric, Numeric> line_ellipsoid_altitude_intersect(
    const Numeric alt,
    const Vector3 ecef,
    const Vector3 decef,
    const Vector2 ell);

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

/** Return distance between neighboring points in a path
 * 
 * @param path A list of path points
 * @param ellipsoid The ellipsoid of the body
 * @return Vector 
 */
Vector distance(const ArrayOfPropagationPathPoint& path,
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
 * @param[in] ell The ellipsoid of the body [a, b]
 * @param[in] alt The tangent altitude
 * @param[in] azimuth The azimuth of the sensor
 * @return The zenith angle
 */
Numeric geometric_tangent_zenith(const Vector3 pos,
                                 const Vector2& ell,
                                 const Numeric alt,
                                 const Numeric azimuth = 0);

/** Find the distance to propagate decef from ecef to find the intersection with
 * the input altitude.
 * 
 * @param ecef Position in ECEF coordinates
 * @param decef Line-of-sight in ECEF coordinates
 * @param refellipsoid Reference ellipsoid [a, b]
 * @param altitude Intersection altitude
 * @param l_min Safety distance to move away from "here" (e.g., for limb-finding)
 * @return Numeric 
 */
Numeric intersection_altitude(const Vector3 ecef,
                              const Vector3 decef,
                              const Vector2 refellipsoid,
                              const Numeric altitude,
                              const Numeric l_min);

std::pair<Vector3, Vector3> geodetic_poslos2ecef(const Vector3 pos,
                                                 const Vector2 los,
                                                 const Vector2 ell) noexcept;
std::pair<Vector3, Vector2> ecef2geodetic_poslos(
    const Vector3 ecef,
    const Vector3 decef,
    const Vector2 refellipsoid) noexcept;

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

bool is_valid_old_pos(const Vector3& pos);
bool is_valid_old_pos(const StridedConstVectorView& pos);
}  // namespace path

using PropagationPathPoint        = path::PropagationPathPoint;
using ArrayOfPropagationPathPoint = path::ArrayOfPropagationPathPoint;
using ArrayOfArrayOfPropagationPathPoint =
    path::ArrayOfArrayOfPropagationPathPoint;
using ArrayOfArrayOfArrayOfPropagationPathPoint =
    path::ArrayOfArrayOfArrayOfPropagationPathPoint;

template <>
struct std::formatter<PropagationPathPoint> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const PropagationPathPoint& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();

    tags.add_if_bracket(ctx, '[');
    if (tags.io) {
      tags.format(ctx, v.pos, sep, v.los, sep, v.nreal, sep, v.ngroup);
    } else {
      tags.format(ctx,
                  v.pos_type,
                  sep,
                  v.los_type,
                  sep,
                  v.pos,
                  sep,
                  v.los,
                  sep,
                  v.nreal,
                  sep,
                  v.ngroup);
    }
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

template <>
struct xml_io_stream<PropagationPathPoint> {
  static constexpr std::string_view type_name = "PropagationPathPoint"sv;

  static void write(std::ostream& os,
                    const PropagationPathPoint& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   PropagationPathPoint& x,
                   bifstream* pbifs = nullptr);
};
