#include <arts_conversions.h>
#include <configtypes.h>
#include <geodetic.h>
#include <matpack.h>
#include <planet_data.h>

#include <cmath>
#include <format>
#include <print>

#include "time_test_util.h"

namespace {
struct Stats {
  Numeric max_dh   = 0;
  Numeric max_dlat = 0;
  Numeric max_dlon = 0;
  Numeric avg_dh   = 0;
  Numeric avg_dlat = 0;
  Numeric avg_dlon = 0;
  Index count      = 0;
};

void update_stats(Stats& s, const Vector3& old_pos, const Vector3& new_pos) {
  const Numeric dh   = std::abs(old_pos[0] - new_pos[0]);
  const Numeric dlat = std::abs(old_pos[1] - new_pos[1]);
  const Numeric dlon = std::abs(old_pos[2] - new_pos[2]);

  s.max_dh   = std::max(s.max_dh, dh);
  s.max_dlat = std::max(s.max_dlat, dlat);
  s.max_dlon = std::max(s.max_dlon, dlon);

  s.avg_dh   += dh;
  s.avg_dlat += dlat;
  s.avg_dlon += dlon;
  ++s.count;
}

void print_stats(const Stats& s, std::string_view label) {
  if (s.count == 0) {
    std::println("{}: no samples", label);
    return;
  }
  const auto n = static_cast<Numeric>(s.count);
  std::println(R"({}:
  max_dh   = {:.6e} m
  max_dlat = {:.6e} deg
  max_dlon = {:.6e} deg
  avg_dh   = {:.6e} m
  avg_dlat = {:.6e} deg
  avg_dlon = {:.6e} deg
  samples  = {})",
               label,
               s.max_dh,
               s.max_dlat,
               s.max_dlon,
               s.avg_dh / n,
               s.avg_dlat / n,
               s.avg_dlon / n,
               s.count);
}

inline constexpr Numeric RAD2DEG = Conversion::rad2deg(1);

/** Threshold for near-pole detection in ECEF to geodetic conversion

    When the horizontal distance from the z-axis (sqrt(x^2 + y^2)) is
    smaller than this threshold times the semi-major axis, the point
    is treated as being on the pole to avoid numerical instability.
*/
inline constexpr Numeric near_pole_threshold_ecef = 1e-15;

/** Maximum iterations for ECEF to geodetic conversion

     Prevents infinite loops in edge cases where the iterative
     algorithm fails to converge.
*/
inline constexpr Index max_ecef2geodetic_iter = 100;

/** Original iterative algorithm for ECEF to geodetic conversion */
Vector3 ecef2geodetic_old(Vector3 ecef, Vector2 refellipsoid) {
  Vector3 pos;

  // Use geocentric function if geoid is spherical
  if (is_ellipsoid_spherical(refellipsoid)) {
    pos     = ecef2geocentric(ecef);
    pos[0] -= refellipsoid[0];

  } else {
    // The general algorithm not stable for lat=+-90. Catch these cases
    // Also catch near-pole cases where numerical instability can occur
    const Numeric sq = std::hypot(ecef[0], ecef[1]);

    if (sq < near_pole_threshold_ecef * refellipsoid[0]) {
      // Near-pole case: use simplified formula (same as exact pole)
      pos[0] = fabs(ecef[2]) - refellipsoid[1];
      pos[1] = ecef[2] >= 0 ? 90 : -90;
      pos[2] = RAD2DEG * atan2(ecef[1], ecef[0]);

    } else {
      // General algorithm
      pos[2] = RAD2DEG * atan2(ecef[1], ecef[0]);

      Numeric B0       = atan2(ecef[2], sq);
      Numeric B        = B0 - 1, N;
      const Numeric e2 = 1 - (refellipsoid[1] * refellipsoid[1]) /
                                 (refellipsoid[0] * refellipsoid[0]);
      Index num_iter   = 0;
      // 1e-15 seems to give a accuracy of better than 2 cm
      while (fabs(B - B0) > 1e-15 && num_iter < max_ecef2geodetic_iter) {
        N      = refellipsoid[0] / sqrt(1 - e2 * sin(B0) * sin(B0));
        pos[0] = sq / cos(B0) - N;
        B      = B0;
        B0     = atan((ecef[2] / sq) * 1 / (1 - e2 * N / (N + pos[0])));
        ++num_iter;
      }
      ARTS_USER_ERROR_IF(num_iter == max_ecef2geodetic_iter,
                         "ECEF to geodetic conversion did not converge. "
                         "Input may be too close to singular values.");
      pos[1] = RAD2DEG * B;
    }
  }

  return pos;
}

}  // namespace

int main() {
  // WGS-84 ellipsoid (a, b)
  const Vector2 wgs84{Body::Earth::a, Body::Earth::b};
  const Vector2 sphere{6371000.0, 6371000.0};

  // Test grids
  const Vector alts{
      -1000.0,
      0.0,
      100.0,
      1000.0,
      10000.0,
      100000.0,
      250000.0,
      3.6e8,   // Distance to the Moon (m)
      3.3e11,  // Distance to Mars (m)
      8.7e11,  // Distance to Jupiter (m)
      1.5e12,  // Distance to Saturn (m)
      3.1e12,  // Distance to Uranus (m)
      5.2e12,  // Distance to Pluto (m)
  };
  const Vector lats{-90.0,
                    -89.99,
                    -80.0,
                    -60.0,
                    -45.0,
                    -30.0,
                    -15.0,
                    0.0,
                    15.0,
                    30.0,
                    45.0,
                    60.0,
                    80.0,
                    89.99,
                    90.0};
  const Vector lons{
      -180.0, -135.0, -90.0, -45.0, 0.0, 45.0, 90.0, 135.0, 180.0};

  Stats stats_wgs84_high;     // WGS-84 overall for high altitude (>100km)
  Stats stats_wgs84_low;      // WGS-84 overall for low altitude (<=100km)
  Stats stats_sphere_high;    // Sphere for high altitude (>100km)
  Stats stats_sphere_low;     // Sphere for low altitude (<=100km)
  Stats stats_poles_high;     // |lat| > 89.9 for high altitude (>100km)
  Stats stats_poles_low;      // |lat| > 89.9 for low altitude (<=100km)
  Stats stats_mid_high;       // |lat| < 80 for high altitude (>100km)
  Stats stats_mid_low;        // |lat| < 80 for low altitude (<=100km)
  Stats stats_wgs84_rt_low;   // WGS-84 round-trip for low altitude (<=100km)
  Stats stats_wgs84_rt_high;  // WGS-84 round-trip for high altitude (>100km)
  Stats stats_poles_rt_low;   // Near-pole round-trip for low altitude (<=100km)
  Stats stats_poles_rt_high;  // Near-pole round-trip for high altitude (>100km)
  Stats stats_mid_rt_low;     // Mid-lat round-trip for low altitude (<=100km)
  Stats stats_mid_rt_high;    // Mid-lat round-trip for high altitude (>100km)
  Stats stats_sphere_rt_low;  // Sphere round-trip for low altitude (<=100km)
  Stats stats_sphere_rt_high;  // Sphere round-trip for high altitude (>100km)

  // --- Accuracy test ---
  for (Numeric h : alts) {
    for (Numeric lat : lats) {
      for (Numeric lon : lons) {
        Vector3 geo{h, lat, lon};

        const bool is_high_altitude = std::abs(h) > 100e3;

        // Round-trip via old function
        Vector3 ecef    = geodetic2ecef(geo, wgs84);
        Vector3 old_pos = ecef2geodetic_old(ecef, wgs84);
        Vector3 new_pos = ecef2geodetic(ecef, wgs84);

        if (is_high_altitude)
          update_stats(stats_wgs84_high, old_pos, new_pos);
        else
          update_stats(stats_wgs84_low, old_pos, new_pos);

        if (is_high_altitude)
          update_stats(stats_wgs84_rt_high, geo, new_pos);
        else
          update_stats(stats_wgs84_rt_low, geo, new_pos);

        if (std::abs(lat) > 89.9) {
          if (is_high_altitude)
            update_stats(stats_poles_high, old_pos, new_pos);
          else
            update_stats(stats_poles_low, old_pos, new_pos);

          if (is_high_altitude)
            update_stats(stats_poles_rt_high, geo, new_pos);
          else
            update_stats(stats_poles_rt_low, geo, new_pos);
        } else if (std::abs(lat) < 80.0) {
          if (is_high_altitude)
            update_stats(stats_mid_high, old_pos, new_pos);
          else
            update_stats(stats_mid_low, old_pos, new_pos);

          if (is_high_altitude)
            update_stats(stats_mid_rt_high, geo, new_pos);
          else
            update_stats(stats_mid_rt_low, geo, new_pos);
        }

        // Spherical case
        Vector3 ecef_s = geodetic2ecef(geo, sphere);
        Vector3 old_s  = ecef2geodetic_old(ecef_s, sphere);
        Vector3 new_s  = ecef2geodetic(ecef_s, sphere);
        if (is_high_altitude)
          update_stats(stats_sphere_high, old_s, new_s);
        else
          update_stats(stats_sphere_low, old_s, new_s);

        if (is_high_altitude)
          update_stats(stats_sphere_rt_high, geo, new_s);
        else
          update_stats(stats_sphere_rt_low, geo, new_s);
      }
    }
  }

  std::println("=== Accuracy comparison (old vs new ecef2geodetic) ===");
  print_stats(stats_wgs84_low, "WGS-84 overall (low altitude <=100km)");
  print_stats(stats_wgs84_high, "WGS-84 overall (high altitude >100km)");
  print_stats(stats_poles_low, "Near-pole (low altitude <=100km)");
  print_stats(stats_poles_high, "Near-pole (high altitude >100km)");
  print_stats(stats_mid_low, "Mid-latitude (low altitude <=100km)");
  print_stats(stats_mid_high, "Mid-latitude (high altitude >100km)");
  print_stats(stats_sphere_low, "Sphere (low altitude <=100km)");
  print_stats(stats_sphere_high, "Sphere (high altitude >100km)");

  std::println("=== Round-trip error (original geodetic vs ecef2geodetic) ===");
  print_stats(stats_wgs84_rt_low, "WGS-84 round-trip (low altitude <=100km)");
  print_stats(stats_wgs84_rt_high, "WGS-84 round-trip (high altitude >100km)");
  print_stats(stats_poles_rt_low,
              "Near-pole round-trip (low altitude <=100km)");
  print_stats(stats_poles_rt_high,
              "Near-pole round-trip (high altitude >100km)");
  print_stats(stats_mid_rt_low,
              "Mid-latitude round-trip (low altitude <=100km)");
  print_stats(stats_mid_rt_high,
              "Mid-latitude round-trip (high altitude >100km)");
  print_stats(stats_sphere_rt_low, "Sphere round-trip (low altitude <=100km)");
  print_stats(stats_sphere_rt_high, "Sphere round-trip (high altitude >100km)");

  // --- Performance test ---
  constexpr Index n_perf = 1'000'000;

  // Build a list of ECEF points for benchmarking
  std::vector<Vector3> ecef_points;
  ecef_points.reserve(lats.size() * lons.size() * alts.size());
  for (Numeric h : alts) {
    for (Numeric lat : lats) {
      for (Numeric lon : lons) {
        ecef_points.push_back(geodetic2ecef(Vector3{h, lat, lon}, wgs84));
      }
    }
  }

  // Warm-up to reduce cache effects
  for (const auto& e : ecef_points) {
    volatile Vector3 p1 = ecef2geodetic_old(e, wgs84);
    volatile Vector3 p2 = ecef2geodetic(e, wgs84);
    (void)p1;
    (void)p2;
  }

  {
    test_timer_t timer("ecef2geodetic (old)");
    for (Index i = 0; i < n_perf; ++i) {
      volatile Vector3 p =
          ecef2geodetic_old(ecef_points[i % ecef_points.size()], wgs84);
      (void)p;
    }
  }

  {
    test_timer_t timer("ecef2geodetic_new");
    for (Index i = 0; i < n_perf; ++i) {
      volatile Vector3 p =
          ecef2geodetic(ecef_points[i % ecef_points.size()], wgs84);
      (void)p;
    }
  }

  print_time_points();

  return 0;
}
