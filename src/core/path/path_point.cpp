#include "path_point.h"

#include <arts_constexpr_math.h>
#include <arts_conversions.h>
#include <atm.h>
#include <configtypes.h>
#include <debug.h>
#include <geodetic.h>
#include <nonstd.h>
#include <surf.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <utility>

namespace path {
Numeric PropagationPathPoint::altitude() const noexcept { return pos[0]; }
Numeric& PropagationPathPoint::altitude() noexcept { return pos[0]; }
Numeric PropagationPathPoint::latitude() const noexcept { return pos[1]; }
Numeric& PropagationPathPoint::latitude() noexcept { return pos[1]; }
Numeric PropagationPathPoint::longitude() const noexcept { return pos[2]; }
Numeric& PropagationPathPoint::longitude() noexcept { return pos[2]; }
Numeric PropagationPathPoint::zenith() const noexcept { return los[0]; }
Numeric& PropagationPathPoint::zenith() noexcept { return los[0]; }
Numeric PropagationPathPoint::azimuth() const noexcept { return los[1]; }
Numeric& PropagationPathPoint::azimuth() noexcept { return los[1]; }

Vector2 mirror(const Vector2 los) {
  Vector2 los_mirrored;
  los_mirrored[0] = 180 - los[0];
  los_mirrored[1] = los[1] + 180;
  if (los_mirrored[1] > 180) los_mirrored[1] -= 360;
  return los_mirrored;
}

namespace {
Numeric surf_altitude(const SurfaceField& surf_field,
                      const Numeric lat,
                      const Numeric lon) {
  return surf_field.contains(SurfaceKey::h)
             ? surf_field.single_value(SurfaceKey::h, lat, lon)
             : 0.0;
}

std::pair<Numeric, Numeric> minmax_surf_altitude(
    const SurfaceField& surf_field) {
  return surf_field.contains(SurfaceKey::h)
             ? surf_field.minmax_single_value(SurfaceKey::h)
             : std::pair<Numeric, Numeric>{0.0, 0.0};
}

Numeric find_crossing_with_surf_z(const Vector3 pos,
                                  const Vector2 los,
                                  const Vector3 ecef,
                                  const Vector3 decef,
                                  const SurfaceField& surf_field,
                                  const Numeric& safe_search_accuracy,
                                  const bool search_safe) {
  // Find min and max surface altitude
  const auto [z_min, z_max] = minmax_surf_altitude(surf_field);

  // Catch upward looking cases that can not have a surface intersection
  if (pos[0] >= z_max && los[0] <= 90) {
    return -1;
  }

  // Check that observation position is above ground
  if (pos[0] < z_max) {
    const Numeric z_surf = surf_altitude(surf_field, pos[1], pos[2]);
    if (pos[0] < z_surf - safe_search_accuracy)
      ARTS_USER_ERROR(
          "The sensor is below the surface. Not allowed!\n"
          "The sensor altitude is at {} m\n"
          "The surface altitude is {} m\n"
          "The position is (lat, lon): ({}, {})",
          pos[0],
          z_surf,
          pos[1],
          pos[2]);
  }

  // Constant surface altitude (in comparison to *safe_search_accuracy*)
  if (z_max - z_min < safe_search_accuracy / 100) {
    // Catch cases with position on the ground, as they can fail if
    // intersection_altitude is used
    if (pos[0] <= z_max) {
      return 0.0;
    }
    return intersection_altitude(ecef, decef, surf_field.ellipsoid, z_min, 0);
  }

  // The general case
  // Find a distance that is guaranteed above or at surface
  // If below z_max, this distance is 0. Otherwise given by z_max
  Numeric l_min;
  if (pos[0] <= z_max)
    l_min = 0;
  else {
    l_min = intersection_altitude(ecef, decef, surf_field.ellipsoid, z_max, 0);
    // No intersection if not even z_max is reached
    if (l_min < 0) return -1;
  }
  // Find max distance for search.
  // If below z_max and upward, given by z_max
  // Otherwise in general given by z_min. If z_min not reached, the distance
  // is instead given by tangent point
  Numeric l_max;
  bool l_max_could_be_above_surface = false;
  if (pos[0] <= z_max && los[0] <= 90) {
    l_max = intersection_altitude(ecef, decef, surf_field.ellipsoid, z_max, 0);
    l_max_could_be_above_surface = true;
  } else {
    l_max = intersection_altitude(ecef, decef, surf_field.ellipsoid, z_min, 0);
  }
  if (l_max < 0) {
    const Vector3 ecef_tan =
        approx_geometrical_tangent_point(ecef, decef, surf_field.ellipsoid);
    l_max = ecef_distance(ecef, ecef_tan);
    // To not miss intersections just after the tangent point, we add a
    // a distance that depends om planet radius (for Earth 111 km).
    l_max += surf_field.ellipsoid[0] * Conversion::sind(1);
    l_max_could_be_above_surface = true;
  }

  // Safe but slow approach
  // ----------------------
  if (search_safe) {
    Numeric l_test =
        l_min - safe_search_accuracy / 2;  // Remove l/2 to get exact result
    bool above_surface = true;             // if true l_test is 0
    while (above_surface && l_test < l_max) {
      l_test += safe_search_accuracy;
      Vector3 local_pos =
          pos_at_distance(ecef, decef, surf_field.ellipsoid, l_test);
      const Numeric z_surf =
          surf_altitude(surf_field, local_pos[1], local_pos[2]);
      if (local_pos[0] < z_surf) above_surface = false;
    }
    if (above_surface) {
      return -1;
    }
    return l_test - safe_search_accuracy / 2;
  }

  // Bisection search
  // ----------------------
  // If l_max matches a point above the surface, we have no intersection
  // according to this search algorithm. And the search fails. So we need
  // to check that point if status unclear
  if (l_max_could_be_above_surface) {
    Vector3 local_pos =
        pos_at_distance(ecef, decef, surf_field.ellipsoid, l_max);
    Numeric z_surf = surf_altitude(surf_field, local_pos[1], local_pos[2]);
    if (local_pos[0] > z_surf) return -1;
  }
  // Start bisection
  while (l_max - l_min > 2 * safe_search_accuracy) {
    const Numeric l_test = std::midpoint(l_min, l_max);
    Vector3 local_pos =
        pos_at_distance(ecef, decef, surf_field.ellipsoid, l_test);
    Numeric z_surf = surf_altitude(surf_field, local_pos[1], local_pos[2]);
    if (local_pos[0] >= z_surf)
      l_min = l_test;
    else
      l_max = l_test;
  }
  return std::midpoint(l_min, l_max);
}
}  // namespace

PropagationPathPoint init(const Vector3& pos,
                          const Vector2& los,
                          const AtmField& atm_field,
                          const SurfaceField& surf_field,
                          bool as_sensor) {
  using enum PathPositionType;
  ARTS_USER_ERROR_IF(pos[1] > 90 or pos[1] < -90, "Non-polar coordinate")

  if (pos[0] > atm_field.top_of_atmosphere) {
    return PropagationPathPoint{.pos_type = space,
                                .los_type = unknown,
                                .pos      = pos,
                                .los      = as_sensor ? mirror(los) : los};
  }

  const Numeric surf_alt = surf_altitude(surf_field, pos[1], pos[2]);

  if (pos[0] < surf_alt) {
    return PropagationPathPoint{.pos_type = unknown,
                                .los_type = unknown,
                                .pos      = pos,
                                .los      = as_sensor ? mirror(los) : los};
  }

  return PropagationPathPoint{.pos_type = atm,
                              .los_type = unknown,
                              .pos      = pos,
                              .los      = as_sensor ? mirror(los) : los};
}

PropagationPathPoint init_with_lostype(const Vector3& pos,
                                       const Vector2& los,
                                       const AtmField& atm_field,
                                       const SurfaceField& surf_field,
                                       bool as_sensor) {
  PropagationPathPoint p = init(pos, los, atm_field, surf_field, as_sensor);

  using enum PathPositionType;
  switch (p.pos_type) {
    case atm:     [[fallthrough]];
    case space:   [[fallthrough]];
    case surface: {
      ArrayOfPropagationPathPoint path{p};
      set_geometric_extremes(path, atm_field, surf_field);
      for (auto& pt : path) {
        if (pt.los_type == atm) return pt;
      }
      return path[0];
    } break;
    case unknown:
      throw std::runtime_error("Cannot set los_type for unknown position type");
  }
  std::unreachable();
}

namespace {
constexpr Numeric nan = std::numeric_limits<Numeric>::quiet_NaN();

const Numeric& min_geq0(const Numeric& a, const Numeric& b) noexcept {
  const bool at = a > 0;
  const bool bt = b > 0;
  return at and bt ? std::min(a, b) : at ? a : bt ? b : nan;
}
}  // namespace

std::pair<Numeric, Numeric> line_ellipsoid_altitude_intersect(
    const Numeric alt,
    const Vector3 ecef,
    const Vector3 decef,
    const Vector2 ell) {
  using Math::pow2;

  // Normalize coordinates to the altitude-adjusted spheroid to keep numbers small
  const Numeric a = ell[0] + alt;  // equatorial radius
  const Numeric b = ell[1] + alt;  // polar radius

  const Numeric x0 = ecef[0] / a;
  const Numeric y0 = ecef[1] / a;
  const Numeric z0 = ecef[2] / b;
  const Numeric dx = decef[0] / a;
  const Numeric dy = decef[1] / a;
  const Numeric dz = decef[2] / b;

  // Quadratic coefficients in normalized space:
  // A l^2 + B' l + C = 0, with
  //   A = |d|^2, B' = 2 dot(d, x0), C = |x0|^2 - 1
  const Numeric A   = pow2(dx) + pow2(dy) + pow2(dz);
  const Numeric dot = dx * x0 + dy * y0 + dz * z0;  // half of B'
  const Numeric C   = pow2(x0) + pow2(y0) + pow2(z0) - 1;

  assert(std::isnormal(A));

  // Discriminant in terms of 'dot' (since B' = 2*dot): D = B'^2 - 4AC = 4(dot^2 - A*C)
  const Numeric Dhalf = dot * dot - A * C;  // equals D/4
  if (Dhalf < 0) return {nan, nan};

  const Numeric sqrtDhalf = std::sqrt(Dhalf);

  // Robust quadratic solution (q-method) to avoid catastrophic cancellation
  // For A l^2 + B' l + C = 0, let q' = -0.5*(B' + sign(B')*sqrt(D)),
  // then roots are l1 = q'/A, l2 = C/q'. Using dot and sqrtDhalf avoids large values.
  const Numeric q = -(dot + std::copysign(sqrtDhalf, dot));

  Numeric l1, l2;
  if (q == 0) {
    // Tangent or near-tangent case: both roots equal -B'/(2A) = -dot/A
    l1 = -dot / A;
    l2 = l1;
  } else {
    l1 = q / A;
    l2 = C / q;
  }

  // Select smallest positive root first
  const Numeric& l = min_geq0(l1, l2);
  if (Dhalf == 0) return {l, nan};
  return {l, (&l == &l1) ? l2 : l1};
}

namespace {
std::pair<Numeric, Numeric> line_ellipsoid_latitude_intersect(
    const Numeric lat,
    const Vector3 ecef,
    const Vector3 decef,
    const Vector2 ell) {
  // Algorithm based on:
  // lousodrome.net/blog/light/2017/01/03/intersection-of-a-ray-and-a-cone/
  // C: Position of tip of cone
  Vector3 C{0, 0, 0};
  if (not is_ellipsoid_spherical(ell)) {
    // Calculate normal to ellipsoid at (0,lat,0) and use it to determine z
    // of cone tip. The distance from (0,lat,0) to z-axis is l=x/dx, so
    // z of cone tip is z-l*dz.
    const auto [ecefn, n] = geodetic_los2ecef({0, lat, 0}, {0, 0}, ell);
    const Numeric l2axis  = ecefn[0] / n[0];
    C[2]                  = ecefn[2] - l2axis * n[2];
  }

  // V: Vector describing centre of cone
  const Vector3 V{0, 0, (lat > 0) ? 1. : -1.};

  // Angle term (cos(lat)^2)
  const Numeric costerm = Math::pow2(Conversion::cosd(90 - std::abs(lat)));

  // Rename to follow nomenclature on web page
  const Vector3& D = decef;

  // Vector from C to O
  const Vector3 CO = ecef_vector_distance(C, ecef);
  // Dot products repeated
  const Numeric DVdot  = dot(D, V);
  const Numeric COVdot = dot(CO, V);
  // The a, b, c and delta terms
  const Numeric a = DVdot * DVdot - costerm;
  const Numeric b = 2 * (DVdot * COVdot - dot(D, CO) * costerm);
  const Numeric c = COVdot * COVdot - dot(CO, CO) * costerm;
  const Numeric d = b * b - 4 * a * c;

  if (d < 0) {
    return {nan, nan};
  }

  if (d == 0) {
    return {-b / (2 * a), nan};
  }

  const Numeric sqrtd = std::sqrt(d);
  const Numeric aa    = 2 * a;

  Numeric l1 = (-b - sqrtd) / aa;
  // Check that crossing is not with -lat
  if (l1 > 0) {
    const Vector3 P  = ecef_at_distance(ecef, decef, l1);
    const Vector3 PC = ecef_vector_distance(C, P);
    if (dot(PC, V) <= 0) l1 = nan;
  }

  Numeric l2 = (-b + sqrtd) / aa;
  if (l2 > 0) {
    const Vector3 P  = ecef_at_distance(ecef, decef, l2);
    const Vector3 PC = ecef_vector_distance(C, P);
    if (dot(PC, V) <= 0) l2 = nan;
  }

  const Numeric& l = min_geq0(l1, l2);
  return {l, (&l == &l1) ? l2 : l1};
}

Numeric line_ellipsoid_longitude_intersect(const Numeric lon,
                                           const Vector3 pos,
                                           const Vector2 los,
                                           const Vector3 ecef,
                                           const Vector3 decef) {
  if ((los[1] == 0) or (nonstd::abs(los[1]) == 180) or
      (pos[2] > lon and los[1] > 0) or (pos[2] < lon and los[1] < 0)) {
    return nan;
  }

  const Numeric tanlon = Conversion::tand(lon);
  return (ecef[1] - ecef[0] * tanlon) / (decef[0] * tanlon - decef[1]);
}

template <bool mirror_los = true>
PropagationPathPoint path_at_distance(const Vector3 ecef,
                                      const Vector3 decef,
                                      const Vector2 ell,
                                      const Numeric l,
                                      const PathPositionType start,
                                      const PathPositionType end) {
  const auto [pos, los] = poslos_at_distance(ecef, decef, ell, l);

  if constexpr (mirror_los)
    return PropagationPathPoint{
        .pos_type = start, .los_type = end, .pos = pos, .los = mirror(los)};
  else
    return PropagationPathPoint{
        .pos_type = start, .los_type = end, .pos = pos, .los = los};
}

struct Intersections {
  PropagationPathPoint first;
  PropagationPathPoint second;
  bool second_is_valid;
};

Intersections pair_line_ellipsoid_intersect(const PropagationPathPoint& path,
                                            const AtmField& atm_field,
                                            const SurfaceField& surf_field,
                                            const Numeric safe_search_accuracy,
                                            const bool search_safe) {
  using enum PathPositionType;

  const auto [surf_h_min, surf_h_max] = minmax_surf_altitude(surf_field);
  ARTS_USER_ERROR_IF(
      surf_h_max > atm_field.top_of_atmosphere,
      "The top of the atmosphere must be above the surface. "
      "This is not the case for the current surface and atmospheric fields. "
      "The maximum altitude of the surface is {}"
      " m, while the top of the atmosphere is at {}"
      " m. This is not allowed.",
      surf_h_max,
      atm_field.top_of_atmosphere);

  const auto looking_los = mirror(path.los);
  const auto x = geodetic_los2ecef(path.pos, looking_los, surf_field.ellipsoid);

  const bool looking_down = path.los[0] < 90;

  const auto get_r_surface = [&, ecef = x.first, decef = x.second]() {
    return find_crossing_with_surf_z(path.pos,
                                     looking_los,
                                     ecef,
                                     decef,
                                     surf_field,
                                     safe_search_accuracy,
                                     search_safe);
  };

  const auto get_r_atm = [&, ecef = x.first, decef = x.second]() {
    return line_ellipsoid_altitude_intersect(
        atm_field.top_of_atmosphere, ecef, decef, surf_field.ellipsoid);
  };

  const auto get_point =
      [&, ecef = x.first, decef = x.second](
          Numeric r, PathPositionType start, PathPositionType end) {
        PropagationPathPoint p =
            path_at_distance(ecef, decef, surf_field.ellipsoid, r, start, end);
        if (end == surface)
          p.pos[0] = surf_altitude(surf_field, p.pos[1], p.pos[2]);
        else
          p.pos[0] = atm_field.top_of_atmosphere;
        return p;
      };

  const auto error = [&path](PathPositionType end) -> Intersections {
    PropagationPathPoint p = path;
    p.los_type             = end;
    return {.first = p, .second = p, .second_is_valid = false};
  };

  switch (path.pos_type) {
    case unknown: ARTS_USER_ERROR("Unknown start position");
    case space:   {
      if (not looking_down) return error(space);

      const auto [r_start_atm, r_end_atm] = get_r_atm();
      if (nonstd::isnan(r_start_atm) or r_start_atm < 0) return error(space);

      const Numeric r_surface = get_r_surface();

      if (r_end_atm < r_surface or r_surface < 0) {
        return {.first           = get_point(r_start_atm, space, atm),
                .second          = get_point(r_end_atm, atm, space),
                .second_is_valid = true};
      }

      return {.first           = get_point(r_start_atm, space, atm),
              .second          = get_point(r_surface, atm, surface),
              .second_is_valid = true};
    }
    case surface:
      if (looking_down) return error(surface);
      [[fallthrough]];
    case atm: {
      const Numeric r_surface = get_r_surface();
      if (r_surface > 0 or (looking_down and r_surface == 0.0)) {
        return {.first           = get_point(r_surface, atm, surface),
                .second          = path,
                .second_is_valid = false};
      }

      const auto [r0, r1] = get_r_atm();

      if (atm_field.top_of_atmosphere == path.altitude()) {
        return {.first           = get_point(std::max(r0, r1), atm, space),
                .second          = path,
                .second_is_valid = false};
      }

      if (nonstd::isnan(r0) or r0 < 0) {
        return {.first           = get_point(r1, atm, space),
                .second          = path,
                .second_is_valid = false};
      }

      if (nonstd::isnan(r1) or r1 < 0) {
        return {.first           = get_point(r0, atm, space),
                .second          = path,
                .second_is_valid = false};
      }

      return {.first           = get_point(std::min(r0, r1), atm, space),
              .second          = get_point(std::max(r0, r1), atm, space),
              .second_is_valid = true};
    }
  }
  ARTS_USER_ERROR("Invalid start position type");
}

bool not_looking_down(const Vector2 los) { return los[0] >= 90; }
}  // namespace

ArrayOfPropagationPathPoint& set_geometric_extremes(
    ArrayOfPropagationPathPoint& path,
    const AtmField& atm_field,
    const SurfaceField& surf_field,
    const Numeric safe_search_accuracy,
    const bool search_safe) {
  ARTS_USER_ERROR_IF(
      path.size() == 0,
      "Must have at least one path point, please call some init() first")
  ARTS_USER_ERROR_IF(
      path.back().los_type != PathPositionType::unknown,
      "Cannot set extremes for path that knows where it is looking")

  if (not_looking_down(path.back().los) and
      atm_field.top_of_atmosphere <= path.back().pos[0]) {
    path.back().los_type = PathPositionType::space;
    return path;
  }

  const auto [first, second, second_is_valid] = pair_line_ellipsoid_intersect(
      path.back(), atm_field, surf_field, safe_search_accuracy, search_safe);
  path.back().los_type = first.pos_type;
  path.push_back(first);
  if (second_is_valid) path.push_back(second);

  return path;
}

ArrayOfPropagationPathPoint& fill_geometric_stepwise(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surf_field,
    const Numeric max_step) {
  ARTS_USER_ERROR_IF(max_step <= 0, "Must move forward")
  if (path.size() == 0) return path;

  //! NOTE: grows path as it loops, so we cannot use iterators
  for (Size i = 0; i < path.size() - 1; ++i) {
    const auto& p1 = path[i];
    const auto& p2 = path[i + 1];
    if (p1.has(PathPositionType::atm) and p2.has(PathPositionType::atm)) {
      const auto [rad_start, rad_dir] =
          geodetic_los2ecef(p2.pos, p2.los, surf_field.ellipsoid);
      const Numeric distance =
          ecef_distance(rad_start, geodetic2ecef(p1.pos, surf_field.ellipsoid));
      if (distance <= max_step) continue;
      path.reserve(path.size() + static_cast<Size>(distance / max_step));
      for (Numeric d = distance - max_step; d > 0; d -= max_step) {
        path.insert(path.begin() + 1 + i,
                    path_at_distance<false>(rad_start,
                                            rad_dir,
                                            surf_field.ellipsoid,
                                            d,
                                            PathPositionType::atm,
                                            PathPositionType::atm));
        ++i;
      }
    }
  }

  return path;
}

ArrayOfPropagationPathPoint& fill_geometric_by_half_steps(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surf_field,
    const Numeric max_step) {
  ARTS_USER_ERROR_IF(max_step <= 0, "Must move forward")
  if (path.size() == 0) return path;

  //! NOTE: grows path as it loops, so we cannot use iterators
  for (Size i = 0; i < path.size() - 1; ++i) {
    const auto& p1 = path[i];
    const auto& p2 = path[i + 1];

    if (p1.has(PathPositionType::atm) and p2.has(PathPositionType::atm)) {
      const auto [rad_start, rad_dir] =
          geodetic_los2ecef(p2.pos, p2.los, surf_field.ellipsoid);
      const Numeric distance =
          ecef_distance(rad_start, geodetic2ecef(p1.pos, surf_field.ellipsoid));

      if (distance <= max_step) continue;

      path.insert(path.begin() + 1 + i,
                  path_at_distance<false>(rad_start,
                                          rad_dir,
                                          surf_field.ellipsoid,
                                          0.5 * distance,
                                          PathPositionType::atm,
                                          PathPositionType::atm));
      --i;
    }
  }

  return path;
}

ArrayOfPropagationPathPoint& fill_geometric_altitude_crossings(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surf_field,
    const Vector& alt_grid) {
  std::vector<std::pair<Numeric, Numeric>> potential_crossings(2 *
                                                               alt_grid.size());

  const auto filter_lt = [](const Numeric& r) {
    return std::ranges::views::filter(
        [r](const std::pair<Numeric, Numeric>& x) {
          return x.first > 0 and x.first < r;
        });
  };

  //! NOTE: grows path as it loops, so we cannot use iterators
  for (Size i = 0; i < path.size() - 1; ++i) {
    const auto& p1 = path[i];
    const auto& p2 = path[i + 1];
    if (p1.has(PathPositionType::atm) and p2.has(PathPositionType::atm)) {
      const auto [ecef, decef] =
          geodetic_los2ecef(p2.pos, p2.los, surf_field.ellipsoid);
      const Numeric distance =
          ecef_distance(geodetic2ecef(p1.pos, surf_field.ellipsoid), ecef);

      potential_crossings.resize(0);
      for (auto alt : alt_grid) {
        const auto [l0, l1] = line_ellipsoid_altitude_intersect(
            alt, ecef, decef, surf_field.ellipsoid);
        potential_crossings.emplace_back(nonstd::isnan(l0) ? -1 : l0, alt);
        potential_crossings.emplace_back(nonstd::isnan(l1) ? -1 : l1, alt);
      }

      std::ranges::sort(potential_crossings,
                        std::greater{},
                        &std::pair<Numeric, Numeric>::first);
      for (auto r : potential_crossings | filter_lt(distance)) {
        path.insert(path.begin() + 1 + i,
                    path_at_distance<false>(ecef,
                                            decef,
                                            surf_field.ellipsoid,
                                            r.first,
                                            PathPositionType::atm,
                                            PathPositionType::atm))
            ->pos[0] = r.second;
        ++i;
      }
    }
  }

  return path;
}

ArrayOfPropagationPathPoint& fill_geometric_latitude_crossings(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surf_field,
    const Vector& lat_grid) {
  std::vector<std::pair<Numeric, Numeric>> potential_crossings(2 *
                                                               lat_grid.size());

  const auto filter_lt = [](const Numeric& r) {
    return std::ranges::views::filter(
        [r](const std::pair<Numeric, Numeric>& x) {
          return x.first > 0 and x.first < r;
        });
  };

  //! NOTE: grows path as it loops, so we cannot use iterators
  for (Size i = 0; i < path.size() - 1; ++i) {
    const auto& p1 = path[i];
    const auto& p2 = path[i + 1];
    if (p1.has(PathPositionType::atm) and p2.has(PathPositionType::atm)) {
      const auto [ecef, decef] =
          geodetic_los2ecef(p2.pos, p2.los, surf_field.ellipsoid);
      const Numeric distance =
          ecef_distance(geodetic2ecef(p1.pos, surf_field.ellipsoid), ecef);

      potential_crossings.resize(0);
      for (auto lat : lat_grid) {
        const auto [l0, l1] = line_ellipsoid_latitude_intersect(
            lat, ecef, decef, surf_field.ellipsoid);
        potential_crossings.emplace_back(nonstd::isnan(l0) ? -1 : l0, lat);
        potential_crossings.emplace_back(nonstd::isnan(l1) ? -1 : l1, lat);
      }

      std::ranges::sort(potential_crossings,
                        std::greater{},
                        &std::pair<Numeric, Numeric>::first);
      for (auto r : potential_crossings | filter_lt(distance)) {
        path.insert(path.begin() + 1 + i,
                    path_at_distance<false>(ecef,
                                            decef,
                                            surf_field.ellipsoid,
                                            r.first,
                                            PathPositionType::atm,
                                            PathPositionType::atm))
            ->pos[1] = r.second;
        ++i;
      }
    }
  }

  return path;
}

ArrayOfPropagationPathPoint& fill_geometric_longitude_crossings(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surf_field,
    const Vector& lon_grid) {
  std::vector<std::pair<Numeric, Numeric>> potential_crossings(lon_grid.size());

  const auto filter_lt = [](const Numeric& r) {
    return std::ranges::views::filter(
        [r](const std::pair<Numeric, Numeric>& x) {
          return x.first > 0 and x.first < r;
        });
  };

  //! NOTE: grows path as it loops, so we cannot use iterators
  for (Size i = 0; i < path.size() - 1; ++i) {
    const auto& p1 = path[i];
    const auto& p2 = path[i + 1];
    if (p1.has(PathPositionType::atm) and p2.has(PathPositionType::atm)) {
      const auto [ecef, decef] =
          geodetic_los2ecef(p2.pos, p2.los, surf_field.ellipsoid);
      const Numeric distance =
          ecef_distance(geodetic2ecef(p1.pos, surf_field.ellipsoid), ecef);

      potential_crossings.resize(0);
      for (auto lon : lon_grid) {
        const auto l0 = line_ellipsoid_longitude_intersect(
            lon, p2.pos, p2.los, ecef, decef);
        potential_crossings.emplace_back(nonstd::isnan(l0) ? -1 : l0, lon);
      }

      std::ranges::sort(potential_crossings,
                        std::greater{},
                        &std::pair<Numeric, Numeric>::first);
      for (auto r : potential_crossings | filter_lt(distance)) {
        path.insert(path.begin() + 1 + i,
                    path_at_distance<false>(ecef,
                                            decef,
                                            surf_field.ellipsoid,
                                            r.first,
                                            PathPositionType::atm,
                                            PathPositionType::atm))
            ->pos[2] = r.second;
        ++i;
      }
    }
  }

  return path;
}

ArrayOfPropagationPathPoint& fill_geometric_crossings(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surf_field,
    const Vector& alt_grid,
    const Vector& lat_grid,
    const Vector& lon_grid) {
  if (alt_grid.size() > 1)
    fill_geometric_altitude_crossings(path, surf_field, alt_grid);
  if (lat_grid.size() > 1)
    fill_geometric_latitude_crossings(path, surf_field, lat_grid);
  if (lon_grid.size() > 1)
    fill_geometric_longitude_crossings(path, surf_field, lon_grid);
  return path;
}

PropagationPathPoint find_geometric_limb(
    const ArrayOfPropagationPathPoint& path, const SurfaceField& surf_field) {
  const auto pre_limb_point = std::adjacent_find(
      path.begin(),
      path.end(),
      [](const PropagationPathPoint& p1, const PropagationPathPoint& p2) {
        return (p1.los[0] >= 90 and p2.los[0] < 90) or
               (p1.los[0] < 90 and p2.los[0] >= 90);
      });

  ARTS_USER_ERROR_IF(
      pre_limb_point == path.end() or pre_limb_point == path.end() - 1,
      "Needs at least one point before limb point");
  const auto post_limb_point = pre_limb_point + 1;

  const auto post_ecef =
      geodetic2ecef(post_limb_point->pos, surf_field.ellipsoid);
  const auto pre_ecef =
      geodetic2ecef(pre_limb_point->pos, surf_field.ellipsoid);
  const Numeric distance = ecef_distance(pre_ecef, post_ecef);
  ARTS_USER_ERROR_IF(distance == 0,
                     "Distance between points is zero: {:B,} and {:B,}.",
                     post_ecef,
                     pre_ecef);

  const Vector3 decef{(post_ecef[0] - pre_ecef[0]) / distance,
                      (post_ecef[1] - pre_ecef[1]) / distance,
                      (post_ecef[2] - pre_ecef[2]) / distance};

  const auto get_limb_point = [ell   = surf_field.ellipsoid,
                               ecef  = pre_ecef,
                               decef = decef](Numeric dist) {
    return path_at_distance<false>(
        ecef, decef, ell, dist, PathPositionType::atm, PathPositionType::atm);
  };

  Numeric x0 = 0, x1 = distance, x = std::midpoint(x0, x1);
  PropagationPathPoint cur = *pre_limb_point, next = get_limb_point(x);
  while (cur.los[0] != 90.0 and
         std::nextafter(cur.los[0], next.los[0]) != next.los[0]) {
    cur                            = next;
    (cur.los[0] >= 90.0 ? x0 : x1) = x;
    x                              = std::midpoint<Numeric>(x0, x1);
    next                           = get_limb_point(x);
  }

  cur.los = mirror(cur.los);
  return cur;
}

ArrayOfPropagationPathPoint& fill_geometric_limb(
    ArrayOfPropagationPathPoint& x, const SurfaceField& surf_field) {
  const auto limb    = find_geometric_limb(x, surf_field);
  const auto min_pos = std::ranges::min_element(
      x, [](const auto& a, const auto& b) { return a.pos[0] < b.pos[0]; });

  if (min_pos == x.end() - 1) {
    x.push_back(limb);
  } else if ((min_pos + 1)->pos[0] < (min_pos - 1)->pos[0]) {
    x.insert(min_pos + 1, limb);
  } else {
    x.insert(min_pos, limb);
  }

  return x;
}

ArrayOfPropagationPathPoint& erase_closeby(ArrayOfPropagationPathPoint& path,
                                           const SurfaceField& surf_field,
                                           const Numeric min_dist,
                                           const bool first) {
  const auto next = [&]() {
    return std::ranges::adjacent_find(
        path,
        [&](const Vector3& a, const Vector3& b) {
          return ecef_distance(geodetic2ecef(a, surf_field.ellipsoid),
                               geodetic2ecef(b, surf_field.ellipsoid)) <
                 min_dist;
        },
        &PropagationPathPoint::pos);
  };

  for (auto it = next(); it != path.end(); it = next()) {
    path.erase(it + (not first));
  }

  return path;
}

Numeric total_geometric_path_length(const ArrayOfPropagationPathPoint& path,
                                    const SurfaceField& surf_field) {
  const auto in_atm = [](const PropagationPathPoint& p) {
    return p.has(PathPositionType::atm);
  };

  const auto first = std::ranges::find_if(path, in_atm);
  const auto last  = std::ranges::find_if_not(first, path.end(), in_atm) - 1;

  ARTS_USER_ERROR_IF(first == path.end(), "No path points in atmosphere")

  return ecef_distance(geodetic2ecef(first->pos, surf_field.ellipsoid),
                       geodetic2ecef(last->pos, surf_field.ellipsoid));
}

ArrayOfPropagationPathPoint& keep_only_atm(ArrayOfPropagationPathPoint& path) {
  using enum PathPositionType;
  path.erase(std::remove_if(
                 path.begin(),
                 path.end(),
                 [](const PropagationPathPoint& p) { return not p.has(atm); }),
             path.end());
  return path;
}

Numeric geometric_tangent_zenith(const Vector3 pos,
                                 const Vector2& ell,
                                 const Numeric alt,
                                 const Numeric azimuth) {
  ARTS_USER_ERROR_IF(alt >= pos[0],
                     "Tangent altitude cannot be above the position")

  const auto intersects = [pos, alt, azimuth, ell](const Numeric za) {
    const auto [ecef, decef] =
        geodetic_los2ecef(pos, mirror({za, azimuth}), ell);
    const auto [l0, l1] =
        line_ellipsoid_altitude_intersect(alt, ecef, decef, ell);
    return l0 > 0 and l1 > 0;
  };

  Numeric za0 = 0.0, za1 = 90.0;
  while (std::nextafter(za0, za1) != za1) {
    const Numeric za             = std::midpoint(za0, za1);
    (intersects(za) ? za0 : za1) = za;
  }

  return za0;
}

Numeric distance(const Vector3 pos1,
                 const Vector3 pos2,
                 const Vector2 ellipsoid) {
  return ecef_distance(geodetic2ecef(pos1, ellipsoid),
                       geodetic2ecef(pos2, ellipsoid));
}

Vector distance(const ArrayOfPropagationPathPoint& path,
                const Vector2 ellipsoid) {
  const Size N = path.size();
  assert(N != 0);

  Vector dists(N - 1);
  for (Size i = 0; i < N - 1; ++i) {
    dists[i] = distance(path[i].pos, path[i + 1].pos, ellipsoid);
  }

  return dists;
}

ArrayOfPropagationPathPoint& fix_updown_azimuth_to_first(
    ArrayOfPropagationPathPoint& path) {
  if (path.size() == 0) return path;

  constexpr auto atan2_failstate = [](const Vector2 los) {
    return los[1] == 0.0 or los[1] == 90.0 or los[1] == -90.0;
  };
  constexpr auto updown_failstate = [](const Vector2 los) {
    return los[0] == 0.0 or los[0] == 180.0;
  };

  if (std::ranges::all_of(path, atan2_failstate, &PropagationPathPoint::los) or
      std::ranges::all_of(path, updown_failstate, &PropagationPathPoint::los)) {
    std::ranges::for_each(
        path,
        [azimuth = path.front().los[1]](Vector2& los) { los[1] = azimuth; },
        &PropagationPathPoint::los);
  }

  return path;
}

bool is_valid_old_pos(const StridedConstVectorView& pos) {
  return pos.size() != 3 or pos[1] < -90 or pos[1] > 90 or pos[2] < -360 or
         pos[2] > 360;
}

bool is_valid_old_pos(const Vector3& pos) {
  return is_valid_old_pos(pos.view());
}

PropagationPathPoint past_geometric(const PropagationPathPoint& this_geometric,
                                    const AtmField& atm_field,
                                    const SurfaceField& surf_field,
                                    const Numeric max_step,
                                    const Numeric safe_search_accuracy,
                                    const bool search_safe) {
  ARTS_USER_ERROR_IF(max_step <= 0, "Must move forward")

  const auto [ecef, decef] = geodetic_los2ecef(this_geometric.pos,
                                               path::mirror(this_geometric.los),
                                               surf_field.ellipsoid);

  const auto [pos, los] =
      poslos_at_distance(ecef, decef, surf_field.ellipsoid, max_step);

  PropagationPathPoint past_geometric{.pos_type = PathPositionType::atm,
                                      .los_type = PathPositionType::atm,
                                      .pos      = pos,
                                      .los      = path::mirror(los),
                                      .nreal    = this_geometric.nreal,
                                      .ngroup   = this_geometric.ngroup};

  if (past_geometric.altitude() >= atm_field.top_of_atmosphere or
      past_geometric.altitude() <= surf_altitude(surf_field,
                                                 past_geometric.latitude(),
                                                 past_geometric.longitude())) {
    ArrayOfPropagationPathPoint path;
    path.reserve(3);
    path.emplace_back(init(
        this_geometric.pos, this_geometric.los, atm_field, surf_field, false));
    set_geometric_extremes(
        path, atm_field, surf_field, safe_search_accuracy, search_safe);
    past_geometric        = path.back();
    past_geometric.nreal  = this_geometric.nreal;
    past_geometric.ngroup = this_geometric.ngroup;
  }

  return past_geometric;
}
}  // namespace path

void xml_io_stream<PropagationPathPoint>::write(std::ostream& os,
                                                const PropagationPathPoint& x,
                                                bofstream* pbofs,
                                                std::string_view name) {
  XMLTag tag(type_name,
             "name",
             name,
             "pos_type",
             toString(x.pos_type),
             "los_type",
             toString(x.los_type));
  tag.write_to_stream(os);

  if (pbofs) {
    xml_io_stream<Vector3>::put({&x.pos, 1}, pbofs);
    xml_io_stream<Vector2>::put({&x.los, 1}, pbofs);
    xml_io_stream<Numeric>::put({&x.nreal, 1}, pbofs);
    xml_io_stream<Numeric>::put({&x.ngroup, 1}, pbofs);
  } else {
    xml_io_stream<Vector3>::write(os, x.pos, pbofs, "Pos"sv);
    xml_io_stream<Vector2>::write(os, x.los, pbofs, "Los"sv);
    xml_io_stream<Numeric>::write(os, x.nreal, pbofs, "Nreal"sv);
    xml_io_stream<Numeric>::write(os, x.ngroup, pbofs, "Nimag"sv);
  }

  tag.write_to_end_stream(os);
}

void xml_io_stream<PropagationPathPoint>::read(std::istream& is,
                                               PropagationPathPoint& x,
                                               bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  String str;
  tag.get_attribute_value("pos_type", str);
  x.pos_type = to<PathPositionType>(str);
  tag.get_attribute_value("los_type", str);
  x.los_type = to<PathPositionType>(str);

  if (pbifs) {
    xml_io_stream<Vector3>::get({&x.pos, 1}, pbifs);
    xml_io_stream<Vector2>::get({&x.los, 1}, pbifs);
    xml_io_stream<Numeric>::get({&x.nreal, 1}, pbifs);
    xml_io_stream<Numeric>::get({&x.ngroup, 1}, pbifs);
  } else {
    xml_io_stream<Vector3>::read(is, x.pos);
    xml_io_stream<Vector2>::read(is, x.los);
    xml_io_stream<Numeric>::read(is, x.nreal);
    xml_io_stream<Numeric>::read(is, x.ngroup);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
