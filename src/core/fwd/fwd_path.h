#pragma once

#include <matpack.h>
#include <path_point.h>

namespace fwd {
struct path {
  PropagationPathPoint point;
  Size alt_index;
  Size lat_index;
  Size lon_index;
  Numeric alt_weight;
  Numeric lat_weight;
  Numeric lon_weight;
  Numeric distance;  //! From the past point to this point
};

std::vector<path> geometric_planar(const Vector3 pos,
                                   const Vector2 los,
                                   const AscendingGrid& alt,
                                   const LatGrid& lat,
                                   const LonGrid& lon);

void path_from_propagation_path(
    std::vector<path>& out,
    const ArrayOfPropagationPathPoint& propagation_path,
    const AscendingGrid& alt,
    const LatGrid& lat,
    const LonGrid& lon,
    const Vector2 ellipsoid);

std::vector<path> path_from_propagation_path(
    const ArrayOfPropagationPathPoint& propagation_path,
    const AscendingGrid& alt,
    const LatGrid& lat,
    const LonGrid& lon,
    const Vector2 ellipsoid);
}  // namespace fwd
