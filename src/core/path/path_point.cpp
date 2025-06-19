#include "path_point.h"

#include <arts_constexpr_math.h>
#include <arts_conversions.h>
#include <atm.h>
#include <configtypes.h>
#include <debug.h>
#include <nonstd.h>
#include <surf.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <ranges>
#include <utility>

namespace path {
namespace {
/** Size of north and south poles
  
    Latitudes with an absolute value > POLELAT are considered to be on
    the south or north pole. This is needed for definition of azimuth.
*/
constexpr Numeric POLELAT = 90.0 - 1e-8;

/** Threshold for non-spherical ellipsoid

    If the two radii of an ellipsoid differ with less than this value, it is
    treated as spherical for efficiency reasons.
*/
constexpr Numeric ellipsoid_radii_threshold = 1e-3;

constexpr bool is_polar_ecef(const Vector3 ecef) {
  constexpr Numeric test = 1e-8;
  return nonstd::abs(ecef[0]) < test and nonstd::abs(ecef[1]) < test;
}

Vector3 los2enu(const Vector2 los) {
  const Numeric zarad = Conversion::deg2rad(los[0]);
  const Numeric aarad = Conversion::deg2rad(los[1]);
  const Numeric st    = std::sin(zarad);
  return {st * std::sin(aarad), st * std::cos(aarad), std::cos(zarad)};
}

Vector3 geocentric2ecef(const Vector3 pos) {
  const Numeric latrad = Conversion::deg2rad(pos[1]);
  const Numeric lonrad = Conversion::deg2rad(pos[2]);
  Vector3 ecef;
  ecef[0] = pos[0] * std::cos(latrad);  // Common term for x and z
  ecef[1] = ecef[0] * std::sin(lonrad);
  ecef[0] = ecef[0] * std::cos(lonrad);
  ecef[2] = pos[0] * std::sin(latrad);
  return ecef;
}

constexpr bool is_ellipsoid_spherical(const Vector2 ellipsoid) {
  return nonstd::abs(ellipsoid[0] - ellipsoid[1]) < ellipsoid_radii_threshold;
}

constexpr Vector3 geodetic2ecef(const Vector3 pos, const Vector2 refellipsoid) {
  Vector3 ecef;

  // Use geocentric function if geoid is spherical
  if (is_ellipsoid_spherical(refellipsoid)) {
    ecef = geocentric2ecef({pos[0] + refellipsoid[0], pos[1], pos[2]});
  } else {
    // See https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates
    const Numeric latrad = Conversion::deg2rad(pos[1]);
    const Numeric lonrad = Conversion::deg2rad(pos[2]);
    const Numeric sinlat = std::sin(latrad);
    const Numeric coslat = std::cos(latrad);
    const Numeric a2     = refellipsoid[0] * refellipsoid[0];
    const Numeric b2     = refellipsoid[1] * refellipsoid[1];
    const Numeric N =
        a2 / std::sqrt(a2 * coslat * coslat + b2 * sinlat * sinlat);
    const Numeric nhcos = (N + pos[0]) * coslat;
    ecef[0]             = nhcos * std::cos(lonrad);
    ecef[1]             = nhcos * std::sin(lonrad);
    ecef[2]             = ((b2 / a2) * N + pos[0]) * sinlat;
  }

  return ecef;
}
}  // namespace

std::pair<Vector3, Vector3> geodetic_poslos2ecef(const Vector3 pos,
                                                 const Vector2 los,
                                                 const Vector2 ell) noexcept {
  // lat = +-90
  // For lat = +- 90 the azimuth angle gives the longitude along which the
  // LOS goes
  // At the poles, no difference between geocentric and geodetic zenith
  Vector3 ecef, decef;
  if (nonstd::abs(pos[1]) > POLELAT) {
    const Numeric s     = sign(pos[1]);
    const Numeric zarad = Conversion::deg2rad(los[0]);
    const Numeric aarad = Conversion::deg2rad(los[1]);
    ecef[0]             = 0;
    ecef[1]             = 0;
    ecef[2]             = s * (pos[0] + ell[1]);
    decef[2]            = s * std::cos(zarad);
    decef[0]            = std::sin(zarad);
    decef[1]            = decef[0] * std::sin(aarad);
    decef[0]            = decef[0] * std::cos(aarad);
  }

  else {
    const Numeric latrad = Conversion::deg2rad(pos[1]);
    const Numeric lonrad = Conversion::deg2rad(pos[2]);
    const Numeric coslat = std::cos(latrad);
    const Numeric sinlat = std::sin(latrad);
    const Numeric coslon = std::cos(lonrad);
    const Numeric sinlon = std::sin(lonrad);

    ecef = geodetic2ecef(pos, ell);

    const Vector3 enu = los2enu(los);

    // See
    // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
    decef[0] =
        -sinlon * enu[0] - sinlat * coslon * enu[1] + coslat * coslon * enu[2];
    decef[1] =
        coslon * enu[0] - sinlat * sinlon * enu[1] + coslat * sinlon * enu[2];
    decef[2] = coslat * enu[1] + sinlat * enu[2];
  }

  return {ecef, decef};
}

Numeric intersection_altitude(const Vector3 ecef,
                              const Vector3 decef,
                              const Vector2 refellipsoid,
                              const Numeric altitude,
                              const Numeric l_min) {
  Numeric l;
  Vector2 ellipsoid{refellipsoid};
  ellipsoid += altitude;

  // Code taken from Atmlab's ellipsoid_intersection

  // Spherical case
  if (is_ellipsoid_spherical(ellipsoid)) {
    const Numeric p =
        ecef[0] * decef[0] + ecef[1] * decef[1] + ecef[2] * decef[2];
    const Numeric pp = p * p;
    const Numeric q  = ecef[0] * ecef[0] + ecef[1] * ecef[1] +
                      ecef[2] * ecef[2] - ellipsoid[0] * ellipsoid[0];
    if (q > pp)
      l = l_min - 1.0;
    else {
      const Numeric sq = std::sqrt(pp - q);
      l                = min_geq(-p - sq, -p + sq, l_min);
    }
  }

  // Ellipsoid case
  else {
    // Based on https://medium.com/@stephenhartzell/
    // satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6
    const Numeric a   = ellipsoid[0];
    const Numeric b   = ellipsoid[0];
    const Numeric c   = ellipsoid[1];
    const Numeric a2  = a * a;
    const Numeric b2  = b * b;
    const Numeric c2  = c * c;
    const Numeric x2  = ecef[0] * ecef[0];
    const Numeric y2  = ecef[1] * ecef[1];
    const Numeric z2  = ecef[2] * ecef[2];
    const Numeric dx2 = decef[0] * decef[0];
    const Numeric dy2 = decef[1] * decef[1];
    const Numeric dz2 = decef[2] * decef[2];
    const Numeric rad = a2 * b2 * dz2 + a2 * c2 * dy2 - a2 * dy2 * z2 +
                        2 * a2 * decef[1] * decef[2] * ecef[1] * ecef[2] -
                        a2 * dz2 * y2 + b2 * c2 * dx2 - b2 * dx2 * z2 +
                        2 * b2 * decef[0] * decef[2] * ecef[0] * ecef[2] -
                        b2 * dz2 * x2 - c2 * dx2 * y2 +
                        2 * c2 * decef[0] * decef[1] * ecef[0] * ecef[1] -
                        c2 * dy2 * x2;
    if (rad < 0)
      l = -1.0;
    else {
      const Numeric val = -a2 * b2 * decef[2] * ecef[2] -
                          a2 * c2 * decef[1] * ecef[1] -
                          b2 * c2 * decef[0] * ecef[0];
      const Numeric mag = a2 * b2 * dz2 + a2 * c2 * dy2 + b2 * c2 * dx2;
      const Numeric abc = a * b * c * std::sqrt(rad);
      l                 = min_geq((val - abc) / mag, (val + abc) / mag, l_min);
    }
  }
  return l;
}

namespace {
Vector3 ecef2geocentric(const Vector3 ecef) {
  assert(dot(ecef, ecef) > 0);

  Vector3 pos;
  pos[0] = std::hypot(ecef[0], ecef[1], ecef[2]);
  pos[1] = Conversion::asind(ecef[2] / pos[0]);
  pos[2] = Conversion::atan2d(ecef[1], ecef[0]);
  return pos;
}

constexpr Vector3 ecef2geodetic(const Vector3 ecef,
                                const Vector2 refellipsoid) {
  using Math::pow2, Math::pow3;

  assert(not nonstd::isnan(dot(ecef, ecef)));
  assert(stdr::all_of(refellipsoid, Cmp::gt(0)));

  Vector3 pos;
  // Use geocentric function if geoid is spherical
  if (is_ellipsoid_spherical(refellipsoid)) {
    pos     = ecef2geocentric(ecef);
    pos[0] -= refellipsoid[0];

    // The general algorithm not stable for lat=+-90. Catch these cases
  } else if (is_polar_ecef(ecef)) {
    pos[0] = nonstd::abs(ecef[2]) - refellipsoid[1];
    pos[1] = ecef[2] >= 0 ? 90 : -90;
    pos[2] = 0;

    // General algorithm
  } else {
    // From https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#The_application_of_Ferrari's_solution
    const Numeric a = refellipsoid[0];
    const Numeric b = refellipsoid[1];
    const Numeric X = ecef[0];
    const Numeric Y = ecef[1];
    const Numeric Z = ecef[2];

    const Numeric a2  = pow2(a);
    const Numeric b2  = pow2(b);
    const Numeric e2  = (a2 - b2) / a2;
    const Numeric DZ2 = (1 - e2) * Z * Z;
    const Numeric r   = std::hypot(X, Y);
    const Numeric e2p = (a2 - b2) / b2;
    const Numeric F   = 54 * pow2(b * Z);
    const Numeric G   = pow2(r) + DZ2 - e2 * (a2 - b2);
    const Numeric c   = pow2(e2) * F * pow2(r) / pow3(G);
    const Numeric s   = std::cbrt(1 + c + std::sqrt(pow2(c) + 2 * c));
    const Numeric fP  = F / (3 * pow2(G * (s + 1 / s + 1)));
    const Numeric Q   = std::sqrt(1 + 2 * pow2(e2) * fP);
    const Numeric r0  = (-fP * e2 * r) / (1 + Q) +
                       std::sqrt(0.5 * a2 * (1 + 1 / Q) -
                                 fP * DZ2 / (Q * (1 + Q)) - 0.5 * fP * pow2(r));
    const Numeric U  = std::hypot(r - e2 * r0, Z);
    const Numeric V  = std::sqrt(pow2(r - e2 * r0) + DZ2);
    const Numeric z0 = b2 * Z / (a * V);

    pos[0] = U * (1 - b2 / (a * V));
    pos[1] = Conversion::atan2d(Z + e2p * z0, r);
    pos[2] = Conversion::atan2d(Y, X);
  }

  return pos;
}

constexpr Vector3 ecef_at_distance(const Vector3 ecef0,
                                   const Vector3 decef,
                                   const Numeric l) {
  return {ecef0[0] + l * decef[0],
          ecef0[1] + l * decef[1],
          ecef0[2] + l * decef[2]};
}

constexpr Vector3 pos_at_distance(const Vector3 ecef,
                                  const Vector3 decef,
                                  const Vector2 refellipsoid,
                                  const Numeric l) {
  return ecef2geodetic(ecef_at_distance(ecef, decef, l), refellipsoid);
}

Numeric ecef_distance(const Vector3 ecef1, const Vector3 ecef2) {
  return std::hypot(
      ecef2[0] - ecef1[0], ecef2[1] - ecef1[1], ecef2[2] - ecef1[2]);
}

constexpr Vector3 approx_geometrical_tangent_point(const Vector3 ecef,
                                                   const Vector3 decef,
                                                   const Vector2 refellipsoid) {
  // Spherical case (length simply obtained by dot product)
  if (is_ellipsoid_spherical(refellipsoid)) {
    return ecef_at_distance(ecef, decef, -(dot(decef, ecef)));
  }

  // General case
  // The algorithm used for non-spherical cases is derived by Nick Lloyd at
  // University of Saskatchewan, Canada (nick.lloyd@usask.ca), and is part of
  // the operational code for both OSIRIS and SMR on-board- the Odin
  // satellite.

  // It seems that there is some numerical inaccuracy if the observation is
  // done from above one of the poles (lat = +-90deg)

  const Numeric a2 = refellipsoid[0] * refellipsoid[0];
  const Numeric b2 = refellipsoid[1] * refellipsoid[1];
  Vector3 yunit, zunit;

  zunit  = cross(decef, ecef);
  zunit /= std::hypot(zunit[0], zunit[1], zunit[2]);
  yunit  = cross(zunit, decef);
  yunit /= std::hypot(yunit[0], yunit[1], yunit[2]);

  const Numeric yr = dot(ecef, yunit);
  const Numeric xr = dot(ecef, decef);
  const Numeric B  = 2.0 * ((decef[0] * yunit[0] + decef[1] * yunit[1]) / a2 +
                           (decef[2] * yunit[2]) / b2);
  Numeric xx;
  if (B == 0.0) {
    xx = 0.0;
  } else {
    const Numeric A = (decef[0] * decef[0] + decef[1] * decef[1]) / a2 +
                      decef[2] * decef[2] / b2;
    const Numeric C = (yunit[0] * yunit[0] + yunit[1] * yunit[1]) / a2 +
                      yunit[2] * yunit[2] / b2;
    const Numeric K      = -2.0 * A / B;
    const Numeric factor = 1.0 / (A + (B + C * K) * K);
    xx                   = std::sqrt(factor);
    const Numeric yy     = K * ecef[0];
    const Numeric dist1  = (xr - xx) * (xr - xx) + (yr - yy) * (yr - yy);
    const Numeric dist2  = (xr + xx) * (xr + xx) + (yr + yy) * (yr + yy);
    if (dist1 > dist2) xx = -xx;
  }

  return {decef[0] * xx + yunit[0] * yr,
          decef[1] * xx + yunit[1] * yr,
          decef[2] * xx + yunit[2] * yr};
}

constexpr Numeric surface_altitude(const SurfaceField& surface_field,
                                   const Numeric lat,
                                   const Numeric lon) {
  return surface_field.has(SurfaceKey::h)
             ? surface_field.single_value(SurfaceKey::h, lat, lon)
             : 0.0;
}

constexpr std::pair<Numeric, Numeric> minmax_surface_altitude(
    const SurfaceField& surface_field) {
  return surface_field.has(SurfaceKey::h)
             ? surface_field.minmax_single_value(SurfaceKey::h)
             : std::pair<Numeric, Numeric>{0.0, 0.0};
}

constexpr Numeric find_crossing_with_surface_z(
    const Vector3 pos,
    const Vector2 los,
    const Vector3 ecef,
    const Vector3 decef,
    const SurfaceField& surface_field,
    const Numeric& safe_search_accuracy,
    const bool search_safe) {
  // Find min and max surface altitude
  const auto [z_min, z_max] = minmax_surface_altitude(surface_field);

  // Catch upward looking cases that can not have a surface intersection
  if (pos[0] >= z_max && los[0] <= 90) {
    return -1;
  }

  // Check that observation position is above ground
  if (pos[0] < z_max) {
    const Numeric z_surf = surface_altitude(surface_field, pos[1], pos[2]);
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
    return intersection_altitude(
        ecef, decef, surface_field.ellipsoid, z_min, 0);
  }

  // The general case
  // Find a distance that is guaranteed above or at surface
  // If below z_max, this distance is 0. Otherwise given by z_max
  Numeric l_min;
  if (pos[0] <= z_max)
    l_min = 0;
  else {
    l_min =
        intersection_altitude(ecef, decef, surface_field.ellipsoid, z_max, 0);
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
    l_max =
        intersection_altitude(ecef, decef, surface_field.ellipsoid, z_max, 0);
    l_max_could_be_above_surface = true;
  } else {
    l_max =
        intersection_altitude(ecef, decef, surface_field.ellipsoid, z_min, 0);
  }
  if (l_max < 0) {
    const Vector3 ecef_tan =
        approx_geometrical_tangent_point(ecef, decef, surface_field.ellipsoid);
    l_max = ecef_distance(ecef, ecef_tan);
    // To not miss intersections just after the tangent point, we add a
    // a distance that depends om planet radius (for Earth 111 km).
    l_max += surface_field.ellipsoid[0] * Conversion::sind(1);
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
          pos_at_distance(ecef, decef, surface_field.ellipsoid, l_test);
      const Numeric z_surf =
          surface_altitude(surface_field, local_pos[1], local_pos[2]);
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
        pos_at_distance(ecef, decef, surface_field.ellipsoid, l_max);
    Numeric z_surf =
        surface_altitude(surface_field, local_pos[1], local_pos[2]);
    if (local_pos[0] > z_surf) return -1;
  }
  // Start bisection
  while (l_max - l_min > 2 * safe_search_accuracy) {
    const Numeric l_test = std::midpoint(l_min, l_max);
    Vector3 local_pos =
        pos_at_distance(ecef, decef, surface_field.ellipsoid, l_test);
    Numeric z_surf =
        surface_altitude(surface_field, local_pos[1], local_pos[2]);
    if (local_pos[0] >= z_surf)
      l_min = l_test;
    else
      l_max = l_test;
  }
  return std::midpoint(l_min, l_max);
}

Vector2 enu2los(const Vector3 enu) {
  // los[0] came out as Nan for a case as enu[2] was just below -1
  // So let's be safe and normalise enu[2], and get a cheap assert for free
  const Numeric twonorm = std::hypot(enu[0], enu[1], enu[2]);
  assert(nonstd::abs(twonorm - 1.0) < 1e-6);
  return {Conversion::acosd(enu[2] / twonorm),
          Conversion::atan2d(enu[0], enu[1])};
}
}  // namespace

std::pair<Vector3, Vector2> ecef2geodetic_poslos(
    const Vector3 ecef,
    const Vector3 decef,
    const Vector2 refellipsoid) noexcept {
  // Catch nans and zero
  assert(nonstd::abs(dot(decef, decef)) > 0);
  assert(not nonstd::isnan(dot(ecef, ecef)));

  const Vector3 pos = ecef2geodetic(ecef, refellipsoid);

  const Numeric latrad = Conversion::deg2rad(pos[1]);
  const Numeric lonrad = Conversion::deg2rad(pos[2]);
  const Numeric coslat = std::cos(latrad);
  const Numeric sinlat = std::sin(latrad);
  const Numeric coslon = std::cos(lonrad);
  const Numeric sinlon = std::sin(lonrad);

  // See
  // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
  Vector3 enu;
  enu[0] = -sinlon * decef[0] + coslon * decef[1];
  enu[1] = -sinlat * coslon * decef[0] - sinlat * sinlon * decef[1] +
           coslat * decef[2];
  enu[2] = coslat * coslon * decef[0] + coslat * sinlon * decef[1] +
           sinlat * decef[2];

  Vector2 los = enu2los(enu);

  // Azimuth at poles needs special treatment
  if (nonstd::abs(pos[1]) > POLELAT) {
    los[1] = Conversion::atan2d(decef[1], decef[0]);
  }

  return {pos, los};
}

namespace {
std::pair<Vector3, Vector2> poslos_at_distance(const Vector3 ecef,
                                               const Vector3 decef,
                                               const Vector2 ell,
                                               const Numeric distance) {
  assert(not nonstd::isnan(distance));
  const Vector3 new_ecef{ecef[0] + distance * decef[0],
                         ecef[1] + distance * decef[1],
                         ecef[2] + distance * decef[2]};
  return ecef2geodetic_poslos(new_ecef, decef, ell);
}
}  // namespace

PropagationPathPoint init(const Vector3& pos,
                          const Vector2& los,
                          const AtmField& atm_field,
                          const SurfaceField& surface_field,
                          bool as_sensor) {
  using enum PathPositionType;
  ARTS_USER_ERROR_IF(pos[1] > 90 or pos[1] < -90, "Non-polar coordinate")

  if (pos[0] > atm_field.top_of_atmosphere) {
    return PropagationPathPoint{.pos_type = space,
                                .los_type = unknown,
                                .pos      = pos,
                                .los      = as_sensor ? mirror(los) : los};
  }

  const Numeric surface_alt = surface_altitude(surface_field, pos[1], pos[2]);

  if (pos[0] < surface_alt) {
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

namespace {
constexpr Numeric nan = std::numeric_limits<Numeric>::quiet_NaN();

constexpr const Numeric& min_geq0(const Numeric& a, const Numeric& b) noexcept {
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

  const Numeric x0  = ecef[0];
  const Numeric y0  = ecef[1];
  const Numeric z0  = ecef[2];
  const Numeric dx  = decef[0];
  const Numeric dy  = decef[1];
  const Numeric dz  = decef[2];
  const Numeric a   = ell[0] + alt;
  const Numeric b   = ell[1] + alt;
  const Numeric a2  = pow2(a);
  const Numeric b2  = pow2(b);
  const Numeric den = a2 * pow2(dz) + b2 * (pow2(dx) + pow2(dy));
  const Numeric sqr_a =
      a2 * (den - pow2(dy * z0 - dz * y0) - pow2(dx * z0 - dz * x0));
  const Numeric sqr_b  = b2 * pow2(dy * x0 - dx * y0);
  const Numeric sqr    = std::sqrt(sqr_a - sqr_b);
  const Numeric term   = -a2 * dz * z0 - b2 * (dx * x0 + dy * y0);
  const Numeric invden = 1 / den;
  const Numeric t0     = (term + b * sqr) * invden;
  const Numeric t1     = (term - b * sqr) * invden;
  const Numeric& t     = min_geq0(t0, t1);

  if (sqr == 0) return {t, nan};
  return {t, (&t == &t0) ? t1 : t0};
}

namespace {
constexpr Vector3 ecef_vector_distance(const Vector3 ecef0,
                                       const Vector3 ecef1) {
  return {ecef1[0] - ecef0[0], ecef1[1] - ecef0[1], ecef1[2] - ecef0[2]};
}

constexpr std::pair<Numeric, Numeric> line_ellipsoid_latitude_intersect(
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
    const auto [ecefn, n] = geodetic_poslos2ecef({0, lat, 0}, {0, 0}, ell);
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

constexpr Numeric line_ellipsoid_longitude_intersect(const Numeric lon,
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
constexpr PropagationPathPoint path_at_distance(const Vector3 ecef,
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

Intersections pair_line_ellipsoid_intersect(
    const PropagationPathPoint& path,
    const AtmField& atm_field,
    const SurfaceField& surface_field,
    const Numeric safe_search_accuracy,
    const bool search_safe) {
  using enum PathPositionType;

  const auto [surf_h_min, surf_h_max] = minmax_surface_altitude(surface_field);
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
  const auto x =
      geodetic_poslos2ecef(path.pos, looking_los, surface_field.ellipsoid);

  const bool looking_down = path.los[0] < 90;

  const auto get_r_surface = [&, ecef = x.first, decef = x.second]() {
    return find_crossing_with_surface_z(path.pos,
                                        looking_los,
                                        ecef,
                                        decef,
                                        surface_field,
                                        safe_search_accuracy,
                                        search_safe);
  };

  const auto get_r_atm = [&, ecef = x.first, decef = x.second]() {
    return line_ellipsoid_altitude_intersect(
        atm_field.top_of_atmosphere, ecef, decef, surface_field.ellipsoid);
  };

  const auto get_point =
      [&, ecef = x.first, decef = x.second](
          Numeric r, PathPositionType start, PathPositionType end) {
        PropagationPathPoint p = path_at_distance(
            ecef, decef, surface_field.ellipsoid, r, start, end);
        if (end == surface)
          p.pos[0] = surface_altitude(surface_field, p.pos[1], p.pos[2]);
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
      if (not looking_down) {
        return error(space);
      }

      const auto [r_start_atm, r_end_atm] = get_r_atm();
      if (nonstd::isnan(r_start_atm) or r_start_atm < 0) {
        return error(space);
      }

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
      if (looking_down) {
        return error(surface);
      }
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

      if (r1 < 0) {
        return {.first           = get_point(r0, atm, space),
                .second          = path,
                .second_is_valid = false};
      }

      if (r0 < 0) {
        return {.first           = get_point(r1, atm, space),
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

constexpr bool not_looking_down(const Vector2 los) { return los[0] >= 90; }
}  // namespace

ArrayOfPropagationPathPoint& set_geometric_extremes(
    ArrayOfPropagationPathPoint& path,
    const AtmField& atm_field,
    const SurfaceField& surface_field,
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
      path.back(), atm_field, surface_field, safe_search_accuracy, search_safe);
  path.back().los_type = first.pos_type;
  path.push_back(first);
  if (second_is_valid) path.push_back(second);

  return path;
}

ArrayOfPropagationPathPoint& fill_geometric_stepwise(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surface_field,
    const Numeric max_step) {
  ARTS_USER_ERROR_IF(max_step <= 0, "Must move forward")
  if (path.size() == 0) return path;

  //! NOTE: grows path as it loops, so we cannot use iterators
  for (Size i = 0; i < path.size() - 1; ++i) {
    const auto& p1 = path[i];
    const auto& p2 = path[i + 1];
    if (p1.has(PathPositionType::atm) and p2.has(PathPositionType::atm)) {
      const auto [rad_start, rad_dir] =
          geodetic_poslos2ecef(p2.pos, p2.los, surface_field.ellipsoid);
      const Numeric distance = ecef_distance(
          rad_start, geodetic2ecef(p1.pos, surface_field.ellipsoid));
      if (distance <= max_step) continue;
      path.reserve(path.size() + static_cast<Size>(distance / max_step));
      for (Numeric d = distance - max_step; d > 0; d -= max_step) {
        path.insert(path.begin() + 1 + i,
                    path_at_distance<false>(rad_start,
                                            rad_dir,
                                            surface_field.ellipsoid,
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
    const SurfaceField& surface_field,
    const Numeric max_step) {
  ARTS_USER_ERROR_IF(max_step <= 0, "Must move forward")
  if (path.size() == 0) return path;

  //! NOTE: grows path as it loops, so we cannot use iterators
  for (Size i = 0; i < path.size() - 1; ++i) {
    const auto& p1 = path[i];
    const auto& p2 = path[i + 1];

    if (p1.has(PathPositionType::atm) and p2.has(PathPositionType::atm)) {
      const auto [rad_start, rad_dir] =
          geodetic_poslos2ecef(p2.pos, p2.los, surface_field.ellipsoid);
      const Numeric distance = ecef_distance(
          rad_start, geodetic2ecef(p1.pos, surface_field.ellipsoid));

      if (distance <= max_step) continue;

      path.insert(path.begin() + 1 + i,
                  path_at_distance<false>(rad_start,
                                          rad_dir,
                                          surface_field.ellipsoid,
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
    const SurfaceField& surface_field,
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
          geodetic_poslos2ecef(p2.pos, p2.los, surface_field.ellipsoid);
      const Numeric distance =
          ecef_distance(geodetic2ecef(p1.pos, surface_field.ellipsoid), ecef);

      potential_crossings.resize(0);
      for (auto alt : alt_grid) {
        const auto [l0, l1] = line_ellipsoid_altitude_intersect(
            alt, ecef, decef, surface_field.ellipsoid);
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
                                            surface_field.ellipsoid,
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
    const SurfaceField& surface_field,
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
          geodetic_poslos2ecef(p2.pos, p2.los, surface_field.ellipsoid);
      const Numeric distance =
          ecef_distance(geodetic2ecef(p1.pos, surface_field.ellipsoid), ecef);

      potential_crossings.resize(0);
      for (auto lat : lat_grid) {
        const auto [l0, l1] = line_ellipsoid_latitude_intersect(
            lat, ecef, decef, surface_field.ellipsoid);
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
                                            surface_field.ellipsoid,
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
    const SurfaceField& surface_field,
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
          geodetic_poslos2ecef(p2.pos, p2.los, surface_field.ellipsoid);
      const Numeric distance =
          ecef_distance(geodetic2ecef(p1.pos, surface_field.ellipsoid), ecef);

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
                                            surface_field.ellipsoid,
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
    const SurfaceField& surface_field,
    const Vector& alt_grid,
    const Vector& lat_grid,
    const Vector& lon_grid) {
  if (alt_grid.size() > 1)
    fill_geometric_altitude_crossings(path, surface_field, alt_grid);
  if (lat_grid.size() > 1)
    fill_geometric_latitude_crossings(path, surface_field, lat_grid);
  if (lon_grid.size() > 1)
    fill_geometric_longitude_crossings(path, surface_field, lon_grid);
  return path;
}

PropagationPathPoint find_geometric_limb(
    const ArrayOfPropagationPathPoint& path,
    const SurfaceField& surface_field) {
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
      geodetic2ecef(post_limb_point->pos, surface_field.ellipsoid);
  const auto pre_ecef =
      geodetic2ecef(pre_limb_point->pos, surface_field.ellipsoid);
  const Numeric distance = ecef_distance(pre_ecef, post_ecef);
  ARTS_USER_ERROR_IF(distance == 0,
                     "Distance between points is zero: {:B,} and {:B,}.",
                     post_ecef,
                     pre_ecef);

  const Vector3 decef{(post_ecef[0] - pre_ecef[0]) / distance,
                      (post_ecef[1] - pre_ecef[1]) / distance,
                      (post_ecef[2] - pre_ecef[2]) / distance};

  const auto get_limb_point = [ell   = surface_field.ellipsoid,
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
    ArrayOfPropagationPathPoint& x, const SurfaceField& surface_field) {
  const auto limb    = find_geometric_limb(x, surface_field);
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
                                           const SurfaceField& surface_field,
                                           const Numeric min_dist,
                                           const bool first) {
  const auto next = [&]() {
    return std::ranges::adjacent_find(
        path,
        [&](const Vector3& a, const Vector3& b) {
          return ecef_distance(geodetic2ecef(a, surface_field.ellipsoid),
                               geodetic2ecef(b, surface_field.ellipsoid)) <
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
                                    const SurfaceField& surface_field) {
  const auto in_atm = [](const PropagationPathPoint& p) {
    return p.has(PathPositionType::atm);
  };

  const auto first = std::ranges::find_if(path, in_atm);
  const auto last  = std::ranges::find_if_not(first, path.end(), in_atm) - 1;

  ARTS_USER_ERROR_IF(first == path.end(), "No path points in atmosphere")

  return ecef_distance(geodetic2ecef(first->pos, surface_field.ellipsoid),
                       geodetic2ecef(last->pos, surface_field.ellipsoid));
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
        geodetic_poslos2ecef(pos, mirror({za, azimuth}), ell);
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
}  // namespace path
