#pragma once

#include <atm.h>
#include <enums.h>
#include <matpack.h>
#include <surf.h>

namespace path {
ENUMCLASS(PositionType, char, atm, space, subsurface, surface, unknown)

struct RefractionData {
  //! The real part of the refractive index at each path position
  Numeric nreal{1};

  //! The group index of refraction
  Numeric ngroup{1};
};

//! A simple path-point
struct PropagationPathPoint {
  /*! Start and end position types.
  
  These must match with previous and next positions
  when a list of path points is created */
  PositionType pos_type{PositionType::unknown}, los_type{PositionType::unknown};

  //! Position of the point: alt [m], lat [deg], lon [deg]
  Vector3 pos;

  //! Line-of-sight of the point: zenith [deg], azimuth [deg]
  Vector2 los;

  //! Data about refraction
  RefractionData n{};

  [[nodiscard]] constexpr bool has(PositionType x) const noexcept {
    return pos_type == x or los_type == x;
  }

  friend std::ostream& operator<<(std::ostream& os,
                                  const PropagationPathPoint& p);
};

using ArrayOfPropagationPathPoint = std::vector<PropagationPathPoint>;

std::ostream& operator<<(std::ostream& os,
                         const ArrayOfPropagationPathPoint& p);

Vector2 mirror(const Vector2 los);

bool valid_los_pos_pairs(const ArrayOfPropagationPathPoint& path);

ArrayOfPropagationPathPoint init(const Vector3& pos,
                                 const Vector2& los,
                                 const AtmField& atm_field,
                                 const SurfaceField& surface_field,
                                 bool as_sensor = true);

ArrayOfPropagationPathPoint& set_geometric_extremes(
    ArrayOfPropagationPathPoint& path,
    const AtmField& atm_field,
    const SurfaceField& surface_field,
    const Numeric surface_search_accuracy = 0.1,
    const bool surface_search_safe = false);

ArrayOfPropagationPathPoint& fill_geometric_stepwise(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surface_field,
    const Numeric max_step);

ArrayOfPropagationPathPoint& fill_geometric_altitude_crossings(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surface_field,
    const Vector& alt_grid);

ArrayOfPropagationPathPoint& fill_geometric_latitude_crossings(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surface_field,
    const Vector& lat_grid);

ArrayOfPropagationPathPoint& fill_geometric_longitude_crossings(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surface_field,
    const Vector& lon_grid);

PropagationPathPoint find_geometric_limb(
    const ArrayOfPropagationPathPoint& path, const SurfaceField& surface_field);

ArrayOfPropagationPathPoint& fill_geometric_limb(
    ArrayOfPropagationPathPoint& x, const SurfaceField& surface_field);

ArrayOfPropagationPathPoint& erase_closeby(ArrayOfPropagationPathPoint& path,
                                           const SurfaceField& surface_field,
                                           const Numeric min_dist,
                                           const bool first = true);

Numeric total_geometric_path_length(const ArrayOfPropagationPathPoint& path,
                                    const SurfaceField& surface_field);

Numeric geometric_tangent_zenith(const Vector3 pos,
                                 const SurfaceField& surface_field,
                                 const Numeric alt,
                                 const Numeric azimuth = 0);

ArrayOfPropagationPathPoint& keep_only_atm(ArrayOfPropagationPathPoint& path);
}  // namespace path

using ArrayOfPropagationPathPoint = path::ArrayOfPropagationPathPoint;
