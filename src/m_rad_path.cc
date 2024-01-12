#include <path_point.h>

#include <algorithm>

void apply_options(ArrayOfPropagationPathPoint& rad_path,
                   const SurfaceField& surface_field,
                   const bool add_limb,
                   const bool remove_non_atm,
                   const bool fix_updown_azimuth) {
  if (fix_updown_azimuth) path::fix_updown_azimuth_to_first(rad_path);

  if (add_limb) path::fill_geometric_limb(rad_path, surface_field);

  if (remove_non_atm) path::keep_only_atm(rad_path);
}

void rad_pathGeometric(ArrayOfPropagationPathPoint& rad_path,
                       const AtmField& atm_field,
                       const SurfaceField& surface_field,
                       const Vector3& pos,
                       const Vector2& los,
                       const Numeric& max_step,
                       const Index& as_sensor,
                       const Index& add_limb,
                       const Index& remove_non_atm,
                       const Index& fix_updown_azimuth) {
  rad_path.resize(1);
  rad_path[0] = path::init(
      pos, los, atm_field, surface_field, static_cast<bool>(as_sensor));

  path::set_geometric_extremes(rad_path, atm_field, surface_field);

  path::fill_geometric_stepwise(rad_path, surface_field, max_step);

  apply_options(
      rad_path, surface_field, add_limb, remove_non_atm, fix_updown_azimuth);
}

void rad_pathGeometricTangentAltitude(ArrayOfPropagationPathPoint& rad_path,
                                      const AtmField& atm_field,
                                      const SurfaceField& surface_field,
                                      const Vector3& pos,
                                      const Numeric& tangent_altitude,
                                      const Numeric& azimuth,
                                      const Numeric& max_step,
                                      const Index& as_sensor,
                                      const Index& add_limb,
                                      const Index& remove_non_atm,
                                      const Index& fix_updown_azimuth) {
  const Numeric aa =
      static_cast<bool>(as_sensor) ? path::mirror({90, azimuth})[1] : azimuth;

  const Numeric za =
      path::geometric_tangent_zenith(pos, surface_field, tangent_altitude, aa);
  const Vector2 los{za, aa};

  rad_path.resize(1);
  rad_path[0] = path::init(pos, los, atm_field, surface_field, false);

  path::set_geometric_extremes(rad_path, atm_field, surface_field);

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(
          rad_path,
          [](const auto& p) { return p.has(path::PositionType::surface); }),
      "Path extremes contain surface point(s):\n",
      rad_path);

  path::fill_geometric_stepwise(rad_path, surface_field, max_step);

  apply_options(
      rad_path, surface_field, add_limb, remove_non_atm, fix_updown_azimuth);
}
