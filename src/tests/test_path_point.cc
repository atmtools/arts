#include <path_point.h>

#include <iomanip>

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

void test_0_az_at_180_za() {
  const auto atm_field = atm();
  const auto surface_field = surf();

  RandomNumberGenerator<> rng;
  auto get_alt = rng.get<std::uniform_real_distribution>(200e3, 1200e3);
  auto get_lat = rng.get<std::uniform_real_distribution>(-90.0, 90.0);
  auto get_lon = rng.get<std::uniform_real_distribution>(-180.0, 180.0);
  auto get_az = rng.get<std::uniform_real_distribution>(0.0, 360.0);

  Index n = 1000;
  while (n-- > 0) {
    const Vector3 pos{get_alt(), get_lat(), get_lon()};
    const Vector2 los{180, get_az()};
    auto x = path::init(pos, los, atm_field, surface_field);
    path::set_geometric_extremes(x, atm_field, surface_field);
    path::fill_geometric_atmosphere(x, surface_field, 150e3);
    path::keep_only_atm(x);
    for (auto p : x)
      ARTS_USER_ERROR_IF(p.los[1] != 0,
                         "los[1] == ",
                         p.los[1],
                         " for pos = ",
                         pos,
                         " and los = ",
                         los);
  }
}

void test_limb_finder() {
  const auto atm_field = atm();
  const auto surface_field = surf();

  RandomNumberGenerator<> rng;
  auto get_tangent_alt = rng.get<std::uniform_real_distribution>(50e3, 100e3);
  auto get_sensor_alt = rng.get<std::uniform_real_distribution>(500e3, 1200e3);
  auto get_lat = rng.get<std::uniform_real_distribution>(-80.0, 80.0);
  auto get_lon = rng.get<std::uniform_real_distribution>(-180.0, 180.0);
  auto get_az = rng.get<std::uniform_real_distribution>(0.0, 360.0);

  Index n = 1000;
  while (n-- > 0) {
    const Vector3 pos{get_sensor_alt(), get_lat(), get_lon()};
    const Numeric az = get_az();
    const Numeric random_tangent_altitude = get_tangent_alt();
    const Numeric za = path::geometric_tangent_zenith(
        pos, surface_field, random_tangent_altitude, az);
    const Vector2 los = path::mirror({za, az});
    auto x = path::init(pos, los, atm_field, surface_field);
    path::set_geometric_extremes(x, atm_field, surface_field);
    path::fill_geometric_atmosphere(x, surface_field, 10e3);
    path::keep_only_atm(x);
    const Numeric computed_tangent_altitude =
        path::find_geometric_limb(x, surface_field).pos[0];
    ARTS_USER_ERROR_IF(
        nonstd::abs(random_tangent_altitude - computed_tangent_altitude) > 1.0,
        "More than 1 meter difference between random and computed tangent altitude: ",
        random_tangent_altitude,
        " vs. ",
        computed_tangent_altitude,
        " for pos = ",
        pos,
        " and los = ",
        los);
  }
}

int main() {
  test_0_az_at_180_za();
  test_limb_finder();
}