#include "sun_methods.h"

#include <path_point.h>
#include <surf.h>

#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

#include "arts_conversions.h"
#include "configtypes.h"
#include "debug.h"
#include "enums.h"
#include "sun.h"
#include "workspace_class.h"

Vector3 cart2sph_plain(const Vector3 cart) {
  using Conversion::asind;
  using Conversion::atan2d;

  const auto [x, y, z] = cart;
  const auto r = hypot(cart);
  return {r, asind(z / r), atan2d(y, x)};
}

/* 
 * \param[out]  h   Geodetic altitude
 * \param[out]  la  Geodetic latitude
 * \param[out]  lon Longitude
 * \param[in] x   x-coordinate (ECEF)
 * \param[in] y   y-coordinate (ECEF)
 * \param[in] z   z-coordinate (ECEF)
 * \param[in]  refellipsoid As the WSV with the same name.
 *
 * \author Patrick Eriksson
 * \date   2020-09-17
*/
Vector3 cart2geodetic(const Vector3 cart, const Vector2 ell) {
  using Conversion::atan2d;
  using Math::pow2;

  const auto [a, b] = ell;
  const Numeric e2 = 1 - pow2(b / a);

  // Use geocentric function if geoid is spherical
  if (e2 < 1e-14) {
    Vector3 sph = cart2sph_plain(cart);
    sph[0] -= a;  // Convert from geocentric to geodetic altitude
    return sph;
  }

  const auto [x, y, z] = cart;

  const Numeric sq = std::hypot(x, y);
  Numeric B0 = atan2(z, sq);
  Numeric B = B0 - 1, N{};

  Numeric h{};
  while (std::abs(B - B0) > 1e-10) {
    N = a / std::sqrt(1 - e2 * std::sin(B0) * std::sin(B0));
    h = sq / std::cos(B0) - N;
    B = B0;
    B0 = std::atan((z / sq) * 1 / (1 - e2 * N / (N + h)));
  }
  return {h, Conversion::rad2deg(B), atan2d(y, x)};
}

/** Converts ENU unit vector vector to zenith and azimuth

   This function and the sister function enu2zaaa handles transformation of
   line-of-sights, from and to ENU (east-north-up). The ENU vector is
   normalised to have length 1.

   @param[in]    enu_los [e, n, u]-part of LOS unit vector.

   @author Patrick Eriksson
   @date   2020-09-17
 */
Vector2 enu2zaaa(const Vector3 enu_los) {
  using Conversion::acosd;
  using Conversion::atan2d;

  const auto [de, dn, du] = enu_los;
  return {acosd(du), atan2d(de, dn)};
}

/*! 
   The inverse of *geodeticposlos2cart*.
   
   \param   cart     [x, y, z]
   \param   cart_los [dx, dy, dz]
   \param   ell      [a, b]
   \return  [za, aa]

   \author Patrick Eriksson
   \date   2020-09-17
*/
Vector2 cart2geodeticlos(const Vector3 cart,
                         const Vector3 cart_los,
                         const Vector2 ell) {
  using Conversion::cosd;
  using Conversion::sind;

  const auto [h, lat, lon] = cart2geodetic(cart, ell);

  const Numeric coslat = cosd(lat);
  const Numeric sinlat = sind(lat);
  const Numeric coslon = cosd(lon);
  const Numeric sinlon = sind(lon);

  // See
  // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
  const auto [dx, dy, dz] = cart_los;
  const Numeric de = -sinlon * dx + coslon * dy;
  const Numeric dn = -sinlat * coslon * dx - sinlat * sinlon * dy + coslat * dz;
  const Numeric du = coslat * coslon * dx + coslat * sinlon * dy + sinlat * dz;

  return enu2zaaa({de, dn, du});
}

/* Workspace method: Doxygen documentation will be auto-generated */
Vector3 geodetic2cart(const Vector3 pos, const Vector2 ell) {
  const auto [h, lat, lon] = pos;
  using Conversion::cosd;
  using Conversion::sind;
  using Math::pow2;

  const auto [a, b] = ell;

  const Numeric e2 = 1 - pow2(b / a);

  // Use geocentric function if geoid is spherical
  if (e2 < 1e-14) {
    return sph2cart({h + a, lat, lon});
  }

  const Numeric sinlat = sind(lat);
  const Numeric coslat = cosd(lat);
  const Numeric N = a / sqrt(1 - e2 * sinlat * sinlat);

  return {(N + h) * coslat * cosd(lon),
          (N + h) * coslat * sind(lon),
          (N * (1 - e2) + h) * sinlat};
}

Vector2 geometric_los(const Vector3 from, const Vector3 to, const Vector2 ell) {
  const Vector3 cart_from = geodetic2cart(from, ell);
  const Vector3 cart_to = geodetic2cart(to, ell);
  const Vector3 cart_los = normalized(cart_to - cart_from);

  return cart2geodeticlos(cart_from, cart_los, ell);
}

enum class angle_move { za_up, aa_up, aa_do, za_do };

angle_move next(angle_move x) {
  switch (x) {
    case angle_move::za_up:
      return angle_move::aa_up;
    case angle_move::aa_up:
      return angle_move::aa_do;
    case angle_move::aa_do:
      return angle_move::za_do;
    case angle_move::za_do:
      return angle_move::za_up;
  }
  std::unreachable();
}


void find_sun_path(const Workspace& ws,
                   ArrayOfPropagationPathPoint& sun_path,
                   const Sun& sun,
                   const Agenda& propagation_path_observer_agenda,
                   const SurfaceField& surface_field,
                   const Vector3& observer_pos,
                   const Numeric& angle_cut,
                   const bool just_hit) {
  ARTS_ASSERT(angle_cut >= 0.0)

  const Vector3 sun_pos{
      {sun.distance - surface_field.single_value(
                          SurfaceKey::h, observer_pos[1], observer_pos[2]),
       sun.latitude,
       sun.longitude}};
  auto los = geometric_los(observer_pos, sun_pos, surface_field.ellipsoid);
  auto best_los = los;

  constexpr Numeric za_min = 0.0;
   Numeric za_max = 180;

  Numeric best_beta = std::numeric_limits<Numeric>::infinity();
  angle_move move = angle_move::za_up;
  bool za_up{false}, aa_up{false}, aa_do{false}, za_do{false};
  do {
    propagation_path_observer_agendaExecute(
        ws, sun_path, observer_pos, los, propagation_path_observer_agenda);

    if (sun_path.back().los_type == PathPositionType::space) {
      auto [beta, hit] = hit_sun(sun,
                                 sun_path.back().pos,
                                 path::mirror(sun_path.back().los),
                                 surface_field.ellipsoid);
      if (hit and just_hit) return;

      beta = Conversion::rad2deg(std::abs(beta));
      if (beta < best_beta) {
        best_beta = beta;
        best_los = los;
        za_do = za_up = aa_up = aa_do = false;
      } else {
        los = best_los;

        switch (move) {
          case angle_move::za_up:
            za_up = true;
            break;
          case angle_move::aa_up:
            aa_up = true;
            break;
          case angle_move::aa_do:
            aa_do = true;
            break;
          case angle_move::za_do:
            za_do = true;
            break;
        }
      }
        
      move = next(move);
      switch (move) {
        case angle_move::za_up:
          los[0] = std::clamp(los[0] + 0.5 * best_beta, za_min, za_max);
          break;
        case angle_move::aa_up:
          los[1] += 0.5 * best_beta;
          break;
        case angle_move::aa_do:
          los[1] -= 0.5 * best_beta;
          break;
        case angle_move::za_do:
          los[0] = std::clamp(los[0] - 0.5 * best_beta, za_min, za_max);
          break;
      }
    } else {
      za_max = std::nextafter(los[0], za_min);
      los[0] = std::clamp(los[0] - 0.5 * best_beta, za_min, za_max);
      move = angle_move::za_up;
      za_do = za_up = aa_up = aa_do = false;
    }
  } while ((not za_do or not aa_do or not aa_up or not za_up) and
           best_beta >= angle_cut);

  // Add the best path, since it is possible we missed it in the loop (which only eliminates bad paths)
  propagation_path_observer_agendaExecute(
      ws, sun_path, observer_pos, best_los, propagation_path_observer_agenda);

  ARTS_USER_ERROR_IF(sun_path.back().los_type != PathPositionType::space,
                     "Path above the horizon is not found")
}