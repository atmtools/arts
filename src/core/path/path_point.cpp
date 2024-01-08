#include "path_point.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

#include "arts_constexpr_math.h"
#include "arts_conversions.h"
#include "artstime.h"
#include "atm.h"
#include "configtypes.h"
#include "debug.h"
#include "math_funcs.h"
#include "matpack_constexpr.h"
#include "nonstd.h"
#include "surf.h"

namespace path {
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

Vector3 los2enu(const Vector2 los) {
  const Numeric zarad = Conversion::deg2rad(los[0]);
  const Numeric aarad = Conversion::deg2rad(los[1]);
  const Numeric st = std::sin(zarad);
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

Vector3 geodetic2ecef(const Vector3 pos, const Vector2 refellipsoid) {
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
    const Numeric a2 = refellipsoid[0] * refellipsoid[0];
    const Numeric b2 = refellipsoid[1] * refellipsoid[1];
    const Numeric N =
        a2 / std::sqrt(a2 * coslat * coslat + b2 * sinlat * sinlat);
    const Numeric nhcos = (N + pos[0]) * coslat;
    ecef[0] = nhcos * std::cos(lonrad);
    ecef[1] = nhcos * std::sin(lonrad);
    ecef[2] = ((b2 / a2) * N + pos[0]) * sinlat;
  }

  return ecef;
}

std::pair<Vector3, Vector3> geodetic_poslos2ecef(const Vector3 pos,
                                                 const Vector2 los,
                                                 const Vector2 ell) {
  // lat = +-90
  // For lat = +- 90 the azimuth angle gives the longitude along which the
  // LOS goes
  // At the poles, no difference between geocentric and geodetic zenith
  Vector3 ecef, decef;
  if (nonstd::abs(pos[1]) > POLELAT) {
    const Numeric s = sign(pos[1]);
    const Numeric zarad = Conversion::deg2rad(los[0]);
    const Numeric aarad = Conversion::deg2rad(los[1]);
    ecef[0] = 0;
    ecef[1] = 0;
    ecef[2] = s * (pos[0] + ell[1]);
    decef[2] = s * std::cos(zarad);
    decef[0] = std::sin(zarad);
    decef[1] = decef[0] * std::sin(aarad);
    decef[0] = decef[0] * std::cos(aarad);
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
    const Numeric q = ecef[0] * ecef[0] + ecef[1] * ecef[1] +
                      ecef[2] * ecef[2] - ellipsoid[0] * ellipsoid[0];
    if (q > pp)
      l = l_min - 1.0;
    else {
      const Numeric sq = std::sqrt(pp - q);
      l = min_geq(-p - sq, -p + sq, l_min);
    }
  }

  // Ellipsoid case
  else {
    // Based on https://medium.com/@stephenhartzell/
    // satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6
    const Numeric a = ellipsoid[0];
    const Numeric b = ellipsoid[0];
    const Numeric c = ellipsoid[1];
    const Numeric a2 = a * a;
    const Numeric b2 = b * b;
    const Numeric c2 = c * c;
    const Numeric x2 = ecef[0] * ecef[0];
    const Numeric y2 = ecef[1] * ecef[1];
    const Numeric z2 = ecef[2] * ecef[2];
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
      l = min_geq((val - abc) / mag, (val + abc) / mag, l_min);
    }
  }
  return l;
}

Vector3 ecef2geocentric(const Vector3 ecef) {
  Vector3 pos;
  pos[0] = std::hypot(ecef[0], ecef[1], ecef[2]);
  pos[1] = Conversion::asind(ecef[2] / pos[0]);
  pos[2] = Conversion::atan2d(ecef[1], ecef[0]);
  return pos;
}

Vector3 ecef2geodetic(const Vector3 ecef, const Vector2 refellipsoid) {
  using Math::pow2, Math::pow3;

  Vector3 pos;
  // Use geocentric function if geoid is spherical
  if (is_ellipsoid_spherical(refellipsoid)) {
    pos = ecef2geocentric(ecef);
    pos[0] -= refellipsoid[0];

    // The general algorithm not stable for lat=+-90. Catch these cases
  } else if (ecef[0] == 0 && ecef[1] == 0) {
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

    const Numeric a2 = pow2(a);
    const Numeric b2 = pow2(b);
    const Numeric e2 = (a2 - b2) / a2;
    const Numeric DZ2 = (1 - e2) * Z * Z;
    const Numeric r = std::hypot(X, Y);
    const Numeric e2p = (a2 - b2) / b2;
    const Numeric F = 54 * pow2(b * Z);
    const Numeric G = pow2(r) + DZ2 - e2 * (a2 - b2);
    const Numeric c = pow2(e2) * F * pow2(r) / pow3(G);
    const Numeric s = std::cbrt(1 + c + std::sqrt(pow2(c) + 2 * c));
    const Numeric fP = F / (3 * pow2(G * (s + 1 / s + 1)));
    const Numeric Q = std::sqrt(1 + 2 * pow2(e2) * fP);
    const Numeric r0 = (-fP * e2 * r) / (1 + Q) +
                       std::sqrt(0.5 * a2 * (1 + 1 / Q) -
                                 fP * DZ2 / (Q * (1 + Q)) - 0.5 * fP * pow2(r));
    const Numeric U = std::hypot(r - e2 * r0, Z);
    const Numeric V = std::sqrt(pow2(r - e2 * r0) + DZ2);
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

Vector3 pos_at_distance(const Vector3 ecef,
                        const Vector3 decef,
                        const Vector2 refellipsoid,
                        const Numeric l) {
  return ecef2geodetic(ecef_at_distance(ecef, decef, l), refellipsoid);
}

Numeric ecef_distance(const Vector3 ecef1, const Vector3 ecef2) {
  return std::hypot(
      ecef2[0] - ecef1[0], ecef2[1] - ecef1[1], ecef2[2] - ecef1[2]);
}

Vector3 approx_geometrical_tangent_point(const Vector3 ecef,
                                         const Vector3 decef,
                                         const Vector2 refellipsoid) {
  // Spherical case (length simply obtained by dot product)
  if (is_ellipsoid_spherical(refellipsoid)) {
    return ecef_at_distance(ecef, decef, -(decef * ecef));
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

  zunit = cross3(decef, ecef);
  zunit /= std::hypot(zunit[0], zunit[1], zunit[2]);
  yunit = cross3(zunit, decef);
  yunit /= std::hypot(yunit[0], yunit[1], yunit[2]);

  const Numeric yr = ecef * yunit;
  const Numeric xr = ecef * decef;
  const Numeric B = 2.0 * ((decef[0] * yunit[0] + decef[1] * yunit[1]) / a2 +
                           (decef[2] * yunit[2]) / b2);
  Numeric xx;
  if (B == 0.0) {
    xx = 0.0;
  } else {
    const Numeric A = (decef[0] * decef[0] + decef[1] * decef[1]) / a2 +
                      decef[2] * decef[2] / b2;
    const Numeric C = (yunit[0] * yunit[0] + yunit[1] * yunit[1]) / a2 +
                      yunit[2] * yunit[2] / b2;
    const Numeric K = -2.0 * A / B;
    const Numeric factor = 1.0 / (A + (B + C * K) * K);
    xx = std::sqrt(factor);
    const Numeric yy = K * ecef[0];
    const Numeric dist1 = (xr - xx) * (xr - xx) + (yr - yy) * (yr - yy);
    const Numeric dist2 = (xr + xx) * (xr + xx) + (yr + yy) * (yr + yy);
    if (dist1 > dist2) xx = -xx;
  }

  return {decef[0] * xx + yunit[0] * yr,
          decef[1] * xx + yunit[1] * yr,
          decef[2] * xx + yunit[2] * yr};
}

Numeric surface_altitude(const SurfaceField& surface_field,
                         const Numeric lat,
                         const Numeric lon) {
  return surface_field.has(Surf::Key::h)
             ? surface_field.single_value(Surf::Key::h, lat, lon)
             : 0.0;
}

std::pair<Numeric, Numeric> minmax_surface_altitude(
    const SurfaceField& surface_field) {
  return surface_field.has(Surf::Key::h)
             ? surface_field.minmax_single_value(Surf::Key::h)
             : std::pair<Numeric, Numeric>{0.0, 0.0};
}

Numeric find_crossing_with_surface_z(const Vector3 pos,
                                     const Vector2 los,
                                     const Vector3 ecef,
                                     const Vector3 decef,
                                     const SurfaceField& surface_field,
                                     const Numeric& surface_search_accuracy,
                                     const bool surface_search_safe) {
  // Find min and max surface altitude
  const auto [z_min, z_max] = minmax_surface_altitude(surface_field);

  // Catch upward looking cases that can not have a surface intersection
  if (pos[0] >= z_max && los[0] <= 90) {
    return -1;
  }

  // Check that observation position is above ground
  if (pos[0] < z_max) {
    const Numeric z_surf = surface_altitude(surface_field, pos[1], pos[2]);
    if (pos[0] < z_surf - surface_search_accuracy)
      ARTS_USER_ERROR(
          "The sensor is below the surface. Not allowed!\n"
          "The sensor altitude is at ",
          pos[0],
          " m\n"
          "The surface altitude is ",
          z_surf,
          " m\n"
          "The position is (lat,lon): (",
          pos[1],
          ",",
          pos[2],
          ")");
  }

  // Constant surface altitude (in comparison to *surface_search_accuracy*)
  if (z_max - z_min < surface_search_accuracy / 100) {
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
  if (surface_search_safe) {
    Numeric l_test =
        l_min - surface_search_accuracy / 2;  // Remove l/2 to get exact result
    bool above_surface = true;                // if true l_test is 0
    while (above_surface && l_test < l_max) {
      l_test += surface_search_accuracy;
      Vector3 local_pos =
          pos_at_distance(ecef, decef, surface_field.ellipsoid, l_test);
      const Numeric z_surf =
          surface_altitude(surface_field, local_pos[1], local_pos[2]);
      if (local_pos[0] < z_surf) above_surface = false;
    }
    if (above_surface) {
      return -1;
    }
    return l_test - surface_search_accuracy / 2;
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
  while (l_max - l_min > 2 * surface_search_accuracy) {
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

constexpr bool up_or_down_los(const Vector2 los) {
  return los[0] == 0 or los[0] == 180;
}

Vector2 enu2los(const Vector3 enu) {
  // los[0] came out as Nan for a case as enu[2] was just below -1
  // So let's be safe and normalise enu[2], and get a cheap assert for free
  const Numeric twonorm = std::hypot(enu[0], enu[1], enu[2]);
  ARTS_ASSERT(nonstd::abs(twonorm - 1.0) < 1e-6);
  Vector2 out;
  out[0] = Conversion::acosd(enu[2] / twonorm);
  out[1] = (not up_or_down_los(out)) * Conversion::atan2d(enu[0], enu[1]);
  return out;
}

std::pair<Vector3, Vector2> ecef2geodetic_poslos(const Vector3 ecef,
                                                 const Vector3 decef,
                                                 const Vector2 refellipsoid) {
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

Vector2 mirror(const Vector2 los) {
  Vector2 los_mirrored;
  los_mirrored[0] = 180 - los[0];
  los_mirrored[1] = (not up_or_down_los(los)) * (los[1] + 180);
  if (los_mirrored[1] > 180) {
    los_mirrored[1] -= 360;
  }
  return los_mirrored;
}

std::pair<Vector3, Vector2> poslos_at_distance(const Vector3 ecef,
                                               const Vector3 decef,
                                               const Vector2 ell,
                                               const Numeric l) {
  const Vector3 new_ecef{
      ecef[0] + l * decef[0], ecef[1] + l * decef[1], ecef[2] + l * decef[2]};
  return ecef2geodetic_poslos(new_ecef, decef, ell);
}

ArrayOfPropagationPathPoint init(const Vector3& pos,
                                 const Vector2& los,
                                 const AtmField& atm_field,
                                 const SurfaceField& surface_field,
                                 bool as_sensor) {
  using enum PositionType;
  ARTS_USER_ERROR_IF(pos[1] > 90 or pos[1] < -90, "Non-polar coordinate")

  if (pos[0] >= atm_field.top_of_atmosphere) {
    return {PropagationPathPoint{.pos_type = space,
                                 .los_type = unknown,
                                 .pos = pos,
                                 .los = as_sensor ? mirror(los) : los}};
  }

  const Numeric surface_alt = surface_altitude(surface_field, pos[1], pos[2]);

  if (pos[0] < surface_alt) {
    return {PropagationPathPoint{.pos_type = subsurface,
                                 .los_type = unknown,
                                 .pos = pos,
                                 .los = as_sensor ? mirror(los) : los}};
  }

  return {PropagationPathPoint{.pos_type = pos[0] > surface_alt ? atm : surface,
                               .los_type = unknown,
                               .pos = pos,
                               .los = as_sensor ? mirror(los) : los}};
}

constexpr Numeric nan = std::numeric_limits<Numeric>::quiet_NaN();

constexpr const Numeric& min_geq0(const Numeric& a, const Numeric& b) noexcept {
  const bool at = a > 0;
  const bool bt = b > 0;
  return at and bt ? std::min(a, b) : at ? a : bt ? b : nan;
}

std::pair<Numeric, Numeric> line_ellipsoid_intersect(
    const Numeric alt,
    const Vector3 ecef,
    const Vector3 decef,
    const Vector2 ell) noexcept {
  using Math::pow2;

  const Numeric x0 = ecef[0];
  const Numeric y0 = ecef[1];
  const Numeric z0 = ecef[2];
  const Numeric dx = decef[0];
  const Numeric dy = decef[1];
  const Numeric dz = decef[2];
  const Numeric a = ell[0] + alt;
  const Numeric b = ell[1] + alt;
  const Numeric a2 = pow2(a);
  const Numeric b2 = pow2(b);
  const Numeric den = a2 * pow2(dz) + b2 * (pow2(dx) + pow2(dy));
  const Numeric sqr_a =
      a2 * (den - pow2(dy * z0 - dz * y0) - pow2(dx * z0 - dz * x0));
  const Numeric sqr_b = b2 * pow2(dy * x0 - dx * y0);
  const Numeric sqr = std::sqrt(sqr_a - sqr_b);
  const Numeric term = -a2 * dz * z0 - b2 * (dx * x0 + dy * y0);
  const Numeric invden = 1 / den;
  const Numeric t0 = (term + b * sqr) * invden;
  const Numeric t1 = (term - b * sqr) * invden;
  const Numeric& t = min_geq0(t0, t1);

  return {t, t == t0 ? t1 : t0};
}

template <bool mirror_los = true>
PropagationPathPoint path_at_distance(const Vector3 ecef,
                                      const Vector3 decef,
                                      const Vector2 ell,
                                      const Numeric l,
                                      const PositionType start,
                                      const PositionType end) {
  const auto [pos, los] = poslos_at_distance(ecef, decef, ell, l);
  if constexpr (mirror_los)
    return PropagationPathPoint{
        .pos_type = start, .los_type = end, .pos = pos, .los = mirror(los)};
  else
    return PropagationPathPoint{
        .pos_type = start, .los_type = end, .pos = pos, .los = los};
};

struct Intersections {
  PropagationPathPoint first;
  PropagationPathPoint second;
  bool second_is_valid;
};

Intersections pair_line_ellipsoid_intersect(
    const PropagationPathPoint& path,
    const AtmField& atm_field,
    const SurfaceField& surface_field,
    const Numeric surface_search_accuracy,
    const bool surface_search_safe) {
  using enum PositionType;

  const auto [surf_h_min, surf_h_max] = minmax_surface_altitude(surface_field);
  ARTS_USER_ERROR_IF(
      surf_h_max > atm_field.top_of_atmosphere,
      "The top of the atmosphere must be above the surface. "
      "This is not the case for the current surface and atmospheric fields. "
      "The maximum altitude of the surface is ",
      surf_h_max,
      " m, while the top of the atmosphere is at ",
      atm_field.top_of_atmosphere,
      " m. This is not allowed.");

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
                                        surface_search_accuracy,
                                        surface_search_safe);
  };

  const auto get_r_atm = [&, ecef = x.first, decef = x.second]() {
    return line_ellipsoid_intersect(
        atm_field.top_of_atmosphere, ecef, decef, surface_field.ellipsoid);
  };

  const auto get_point = [&, ecef = x.first, decef = x.second](
                             Numeric r, PositionType start, PositionType end) {
    PropagationPathPoint p =
        path_at_distance(ecef, decef, surface_field.ellipsoid, r, start, end);
    if (end == surface)
      p.pos[0] = surface_altitude(surface_field, p.pos[1], p.pos[2]);
    else
      p.pos[0] = atm_field.top_of_atmosphere;
    return p;
  };

  const auto error = [&path](PositionType end) {
    PropagationPathPoint p = path;
    p.los_type = end;
    return Intersections{p, p, false};
  };

  switch (path.pos_type) {
    case unknown:
      ARTS_USER_ERROR("Unknown start position");
    case space: {
      if (not looking_down) {
        return error(space);
      }

      const auto [r_start_atm, r_end_atm] = get_r_atm();
      if (nonstd::isnan(r_start_atm) or r_start_atm < 0) {
        return error(space);
      }

      const Numeric r_surface = get_r_surface();

      if (r_end_atm < r_surface or r_surface < 0) {
        return {get_point(r_start_atm, space, atm),
                get_point(r_end_atm, atm, space),
                true};
      }

      return {get_point(r_start_atm, space, atm),
              get_point(r_surface, atm, surface),
              true};
    }
    case surface:
      if (looking_down) {
        return error(surface);
      }
      [[fallthrough]];
    case atm:
      if (const Numeric r_surface = get_r_surface(); r_surface > 0) {
        return {get_point(r_surface, atm, surface), path, false};
      }
      return {get_point(get_r_atm().first, atm, space), path, false};
    case subsurface:
      ARTS_USER_ERROR("Unsupported subsurface start position")
    case FINAL:
      throw std::runtime_error("Invalid start position");
  }
}

ArrayOfPropagationPathPoint& set_geometric_extremes(
    ArrayOfPropagationPathPoint& path,
    const AtmField& atm_field,
    const SurfaceField& surface_field,
    const Numeric surface_search_accuracy,
    const bool surface_search_safe) {
  ARTS_USER_ERROR_IF(
      path.size() == 0,
      "Must have at least one path point, please call some init() first")

  const auto [first, second, second_is_valid] =
      pair_line_ellipsoid_intersect(path.back(),
                                    atm_field,
                                    surface_field,
                                    surface_search_accuracy,
                                    surface_search_safe);
  path.back().los_type = first.pos_type;
  path.push_back(first);
  if (second_is_valid) path.push_back(second);
  return path;
}

ArrayOfPropagationPathPoint& fill_geometric_atmosphere(
    ArrayOfPropagationPathPoint& path,
    const SurfaceField& surface_field,
    const Numeric max_step) {
  ARTS_USER_ERROR_IF(
      path.size() == 0,
      "Must have at least one path point, please call some init() first")

  //! NOTE: grows path as it loops, so we cannot use iterators
  for (Size i = 0; i < path.size() - 1; ++i) {
    const auto& p1 = path[i];
    const auto& p2 = path[i + 1];
    if (p1.has(PositionType::atm) and p2.has(PositionType::atm)) {
      const auto [ecef, decef] =
          geodetic_poslos2ecef(p2.pos, p2.los, surface_field.ellipsoid);
      const Numeric distance =
          ecef_distance(geodetic2ecef(p1.pos, surface_field.ellipsoid), ecef);
      path.reserve(path.size() + static_cast<Size>(distance / max_step));
      for (Numeric d = distance - max_step; d > 0; d -= max_step) {
        path.insert(path.begin() + 1 + i,
                    path_at_distance<false>(ecef,
                                            decef,
                                            surface_field.ellipsoid,
                                            d,
                                            PositionType::atm,
                                            PositionType::atm));
        ++i;
      }
    }
  }

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
  const Vector3 decef{post_ecef[0] - pre_ecef[0],
                      post_ecef[1] - pre_ecef[1],
                      post_ecef[2] - pre_ecef[2]};

  const auto get_limb_point = [ell = surface_field.ellipsoid,
                               ecef = pre_ecef,
                               decef = decef](Numeric dist) {
    return path_at_distance<false>(
        ecef, decef, ell, dist, PositionType::atm, PositionType::atm);
  };

  Numeric x0 = 0, x1 = 1.0, x = 0.5;
  PropagationPathPoint cur = *pre_limb_point, next = get_limb_point(x);
  while (std::nextafter(next.pos[0], cur.pos[0]) != cur.pos[0]) {
    cur = next;
    if (cur.los[0] >= 90) {
      x0 = x;
    } else {
      x1 = x;
    }
    x = std::midpoint<Numeric>(x0, x1);
    next = get_limb_point(x);
  }

  next.los[1] = next.los[1] - 180;
  return next;
}

Numeric total_geometric_path_length(const ArrayOfPropagationPathPoint& path,
                                    const SurfaceField& surface_field) {
  const auto in_atm = [](const PropagationPathPoint& p) {
    return p.has(PositionType::atm);
  };

  const auto first = std::ranges::find_if(path, in_atm);
  const auto last = std::ranges::find_if_not(first, path.end(), in_atm) - 1;

  ARTS_USER_ERROR_IF(first == path.end(), "No path points in atmosphere")

  return ecef_distance(geodetic2ecef(first->pos, surface_field.ellipsoid),
                       geodetic2ecef(last->pos, surface_field.ellipsoid));
}

std::ostream& operator<<(std::ostream& os, const PropagationPathPoint& p) {
  os << "pos: [" << p.pos << "], los: [" << p.los
     << "], pos-type: " << p.pos_type << ", los-type: " << p.los_type
     << ", nreal: " << p.n.nreal << ", ngroup: " << p.n.ngroup;
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const ArrayOfPropagationPathPoint& ps) {
  std::string_view sep = "";
  for (const auto& p : ps) {
    os << std::exchange(sep, "\n") << p;
  }
  return os;
}

bool valid_los_pos_pairs(const ArrayOfPropagationPathPoint& path) {
  return path.end() == std::adjacent_find(path.begin(),
                                          path.end(),
                                          [](const PropagationPathPoint& p1,
                                             const PropagationPathPoint& p2) {
                                            return p1.los_type != p2.pos_type;
                                          });
}

ArrayOfPropagationPathPoint& keep_only_atm(ArrayOfPropagationPathPoint& path) {
  using enum PositionType;
  path.erase(std::remove_if(
                 path.begin(),
                 path.end(),
                 [](const PropagationPathPoint& p) { return not p.has(atm); }),
             path.end());
  return path;
}

Numeric geometric_tangent_zenith(const Vector3 pos,
                                 const SurfaceField& surface_field,
                                 const Numeric alt,
                                 const Numeric azimuth) {
  ARTS_USER_ERROR_IF(alt > pos[0],
                     "Tangent altitude cannot be above the position")

  Numeric za0 = 0.0, za1 = 90.0, za = 45.0;
  Numeric t0 = 2 * (max(surface_field.ellipsoid) + alt),
          t1 = std::numeric_limits<Numeric>::quiet_NaN();

  while (std::nextafter(za0, za1) != za1) {
    const auto [ecef, decef] =
        geodetic_poslos2ecef(pos, mirror({za, azimuth}), surface_field.ellipsoid);
    const auto [l0, l1] =
        line_ellipsoid_intersect(alt, ecef, decef, surface_field.ellipsoid);
    if (l0 > 0 and l1 > 0) {
      const Numeric r = l1 - l0;
      if (t1 > t0) {
        t1 = r;
        za1 = za;
      } else {
        t0 = r;
        za0 = za;
      }
    } else {
      za1 = za;
    }
    za = std::midpoint(za0, za1);
  }

  return za;
}
}  // namespace path
