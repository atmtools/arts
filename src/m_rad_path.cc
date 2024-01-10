#include <path_point.h>

void rad_pathGeometric(ArrayOfPropagationPathPoint& rad_path,
                       const AtmField& atm_field,
                       const SurfaceField& surface_field,
                       const Vector3& pos,
                       const Vector2& los,
                       const Numeric& max_step,
                       const Index& as_sensor) {
  rad_path.resize(1);
  rad_path[0] = path::init(
      pos, los, atm_field, surface_field, static_cast<bool>(as_sensor));
  path::set_geometric_extremes(rad_path, atm_field, surface_field);
  path::fill_geometric_stepwise(rad_path, surface_field, max_step);
}

void rad_pathGeometricTangentAltitude(ArrayOfPropagationPathPoint& rad_path,
                                      const AtmField& atm_field,
                                      const SurfaceField& surface_field,
                                      const Vector3& pos,
                                      const Numeric& tangent_altitude,
                                      const Numeric& azimuth,
                                      const Numeric& max_step,
                                      const Index& as_sensor) {
  const Numeric aa =
      (as_sensor == 0) ? azimuth : path::mirror({90, azimuth})[1];

  const Numeric za =
      path::geometric_tangent_zenith(pos, surface_field, tangent_altitude, aa);
  const Vector2 los{za, aa};

  rad_path.resize(1);
  rad_path[0] = path::init(pos, los, atm_field, surface_field, false);
  path::set_geometric_extremes(rad_path, atm_field, surface_field);
  path::fill_geometric_stepwise(rad_path, surface_field, max_step);
  path::fill_geometric_limb(rad_path, surface_field);
}
