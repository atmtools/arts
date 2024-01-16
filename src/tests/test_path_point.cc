#include <path_point.h>

#include <cstdlib>
#include <exception>
#include <stdexcept>

#include "debug.h"
#include "rng.h"

AtmField atm() {
  AtmField atm_field;
  atm_field.top_of_atmosphere = 1e5;
  return atm_field;
}

SurfaceField surf() {
  SurfaceField surface_field;
  surface_field.ellipsoid = {6378137, 6356752.314245};
  return surface_field;
}

void test_0_az_at_180_za(Size n) {
  const auto atm_field = atm();
  const auto surface_field = surf();

  RandomNumberGenerator<> rng;
  auto get_alt = rng.get<std::uniform_real_distribution>(200e3, 1200e3);
  auto get_lat = rng.get<std::uniform_real_distribution>(-90.0, 90.0);
  auto get_lon = rng.get<std::uniform_real_distribution>(-180.0, 180.0);
  auto get_az = rng.get<std::uniform_real_distribution>(-180.0, 180.0);

  while (n-- != 0) {
    const Vector3 pos{get_alt(), get_lat(), get_lon()};
    const Vector2 los{180, get_az()};
    ArrayOfPropagationPathPoint x{
        path::init(pos, los, atm_field, surface_field)};
    path::set_geometric_extremes(x, atm_field, surface_field);
    path::fill_geometric_stepwise(x, surface_field, 
    1e3);
    path::fix_updown_azimuth_to_first(x);
    path::keep_only_atm(x);
    for (auto p : x)
      ARTS_USER_ERROR_IF(p.los[1] != path::mirror(los)[1],
                         "los[1] == ",
                         p.los[1],
                         " for sensor pos = ",
                         pos,
                         " and sensor los = ",
                         path::mirror(los));
  }
}

void test_limb_finder(Size n) {
  const auto atm_field = atm();
  const auto surface_field = surf();

  RandomNumberGenerator<> rng;
  auto get_tangent_alt = rng.get<std::uniform_real_distribution>(50e3, 100e3);
  auto get_sensor_alt = rng.get<std::uniform_real_distribution>(500e3, 1200e3);
  auto get_lat = rng.get<std::uniform_real_distribution>(-80.0, 80.0);
  auto get_lon = rng.get<std::uniform_real_distribution>(-180.0, 180.0);
  auto get_az = rng.get<std::uniform_real_distribution>(0.0, 360.0);

  while (n-- != 0) {
    const Vector3 pos{get_sensor_alt(), get_lat(), get_lon()};
    const Numeric az = get_az();
    const Numeric tangent_altitude = get_tangent_alt();
    const Numeric za = path::geometric_tangent_zenith(
        pos, surface_field, tangent_altitude, az);
    const Vector2 los = path::mirror({za, az});
    ArrayOfPropagationPathPoint x{
        path::init(pos, los, atm_field, surface_field)};
    path::set_geometric_extremes(x, atm_field, surface_field);
    path::fill_geometric_stepwise(x, surface_field, 10e3);
    path::keep_only_atm(x);
    const Numeric computed_tangent_altitude =
        path::find_geometric_limb(x, surface_field).pos[0];
    ARTS_USER_ERROR_IF(
        nonstd::abs(tangent_altitude - computed_tangent_altitude) > 1.0,
        "More than 1 meter difference between set and estimated tangent altitude: ",
        tangent_altitude,
        " vs. ",
        computed_tangent_altitude,
        " for pos = ",
        pos,
        " and los = ",
        los);
  }
}

void test_geometric_fill(Size n) {
  const auto atm_field = atm();
  const auto surface_field = surf();

  const Vector alts = uniform_grid(1e3, 99, 1e3);
  const Vector lats = uniform_grid(-90, 181, 1.);
  const Vector lons = uniform_grid(-180, 361, 1.);

  const auto test = [&](const ArrayOfPropagationPathPoint& x,
                        const Numeric start_alt) {
    for (auto& p : x) {
      const bool toa = p.pos[0] == atm_field.top_of_atmosphere;
      const bool boa = p.pos[0] == 0;
      const bool at_alt = std::any_of(
          alts.begin(), alts.end(), [&](Numeric a) { return a == p.pos[0]; });
      const bool at_lat = std::any_of(
          lats.begin(), lats.end(), [&](Numeric a) { return a == p.pos[1]; });
      const bool at_lon = std::any_of(
          lons.begin(), lons.end(), [&](Numeric a) { return a == p.pos[2]; });
      const bool at_start = p.pos[0] == start_alt;
      const bool any_nan = nonstd::isnan(p.pos[0]) or nonstd::isnan(p.pos[1]) or
                           nonstd::isnan(p.pos[2]) or nonstd::isnan(p.los[0]) or
                           nonstd::isnan(p.los[1]);

      ARTS_USER_ERROR_IF(
          not(toa or boa or at_alt or at_lat or at_lon or at_start) or any_nan,
          "Bad path point: ",
          p,
          "\nIt is not at the top of the atmosphere or at the surface or at one of the grid points; or it contains NaNs.")
    }
  };

  {
    // Check crossing the pole
    const Vector3 pos{300e3, 75, 190};
    const Vector2 los{105, 0};
    ArrayOfPropagationPathPoint x{
        path::init(pos, los, atm_field, surface_field)};
    path::set_geometric_extremes(x, atm_field, surface_field);
    path::fill_geometric_altitude_crossings(x, surface_field, alts);
    path::fill_geometric_latitude_crossings(x, surface_field, lats);
    path::fill_geometric_longitude_crossings(x, surface_field, lons);
    test(x, pos[0]);
  }

  {
    // Check pole upwards
    const Vector3 pos{0.0, 90, 0};
    const Vector2 los{0.0, 0};
    ArrayOfPropagationPathPoint x{
        path::init(pos, los, atm_field, surface_field)};
    path::set_geometric_extremes(x, atm_field, surface_field);
    path::fill_geometric_altitude_crossings(x, surface_field, alts);
    path::fill_geometric_latitude_crossings(x, surface_field, lats);
    path::fill_geometric_longitude_crossings(x, surface_field, lons);
    test(x, pos[0]);
  }

  RandomNumberGenerator<> rng;
  auto get_alt = rng.get<std::uniform_real_distribution>(0.0, 200e3);
  auto get_lat = rng.get<std::uniform_real_distribution>(-90.0, 90.0);
  auto get_lon = rng.get<std::uniform_real_distribution>(-180.0, 180.0);
  auto get_za = rng.get<std::uniform_real_distribution>(0.0, 180.);
  auto get_az = rng.get<std::uniform_real_distribution>(0.0, 360.0);

  while (n-- != 0) {
    const Vector3 pos{get_alt(), get_lat(), get_lon()};
    const Vector2 los{get_za(), get_az()};
    ArrayOfPropagationPathPoint x{
        path::init(pos, los, atm_field, surface_field)};
    path::set_geometric_extremes(x, atm_field, surface_field);
    path::fill_geometric_altitude_crossings(x, surface_field, alts);
    path::fill_geometric_latitude_crossings(x, surface_field, lats);
    path::fill_geometric_longitude_crossings(x, surface_field, lons);
    test(x, pos[0]);
  }
}

int main() try {
  test_0_az_at_180_za(1'000);
  test_limb_finder(1'000);
  test_geometric_fill(1'000);

  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << "Caught error and rethrowing:\n" << e.what() << '\n';
  return EXIT_FAILURE;
}