#include <path_point.h>

#include <algorithm>
#include <functional>

#include "atm.h"
#include "surf.h"
#include "workspace_class.h"

void apply_options(ArrayOfPropagationPathPoint& ray_path,
                   const SurfaceField& surface_field,
                   const bool add_limb,
                   const bool remove_non_atm,
                   const bool fix_updown_azimuth) {
  if (fix_updown_azimuth) path::fix_updown_azimuth_to_first(ray_path);

  if (add_limb) path::fill_geometric_limb(ray_path, surface_field);

  if (remove_non_atm) path::keep_only_atm(ray_path);
}

void ray_pathGeometric(ArrayOfPropagationPathPoint& ray_path,
                       const AtmField& atm_field,
                       const SurfaceField& surface_field,
                       const Vector3& pos,
                       const Vector2& los,
                       const Numeric& max_step,
                       const Numeric& surface_search_accuracy,
                       const Index& as_sensor,
                       const Index& add_limb,
                       const Index& remove_non_atm,
                       const Index& fix_updown_azimuth,
                       const Index& surface_safe_search) {
  ARTS_USER_ERROR_IF(any_nan(pos) or any_nan(los) or nonstd::isnan(max_step),
                     R"(There are NAN in the pos or los vector:
pos:      {:B,}
los:      {:B,}
max_step: {}
)",
                     pos,
                     los,
                     max_step);

  ray_path.resize(1);
  ray_path[0] = path::init(
      pos, los, atm_field, surface_field, static_cast<bool>(as_sensor));

  path::set_geometric_extremes(ray_path,
                               atm_field,
                               surface_field,
                               surface_search_accuracy,
                               static_cast<bool>(surface_safe_search));

  path::fill_geometric_stepwise(ray_path, surface_field, max_step);

  apply_options(
      ray_path, surface_field, add_limb, remove_non_atm, fix_updown_azimuth);
}

void ray_pathGeometricTangentAltitude(ArrayOfPropagationPathPoint& ray_path,
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

  ray_path.resize(1);
  ray_path[0] = path::init(pos, los, atm_field, surface_field, false);

  path::set_geometric_extremes(ray_path, atm_field, surface_field);

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(
          ray_path,
          [](const auto& p) { return p.has(PathPositionType::surface); }),
      "Path extremes contain surface point(s):\n{}",
      ray_path);

  path::fill_geometric_stepwise(ray_path, surface_field, max_step);

  apply_options(
      ray_path, surface_field, add_limb, remove_non_atm, fix_updown_azimuth);
}

void ray_path_pointBackground(PropagationPathPoint& ray_path_point,
                              const ArrayOfPropagationPathPoint& ray_path) {
  ARTS_USER_ERROR_IF(ray_path.size() == 0, "Empty propagation path.")
  ray_path_point = ray_path.back();
}

void ray_path_pointForeground(PropagationPathPoint& ray_path_point,
                              const ArrayOfPropagationPathPoint& ray_path) {
  ARTS_USER_ERROR_IF(ray_path.size() == 0, "Empty propagation path.")
  ray_path_point = ray_path.front();
}

void ray_path_pointLowestFromPath(PropagationPathPoint& ray_path_point,
                                  const ArrayOfPropagationPathPoint& ray_path) {
  ARTS_USER_ERROR_IF(ray_path.size() == 0, "Empty propagation path.")

  ray_path_point = *std::ranges::min_element(
      ray_path,
      [](const PropagationPathPoint& a, const PropagationPathPoint& b) {
        return a.altitude() < b.altitude();
      });
}

void ray_path_pointHighestFromPath(
    PropagationPathPoint& ray_path_point,
    const ArrayOfPropagationPathPoint& ray_path) {
  ARTS_USER_ERROR_IF(ray_path.size() == 0, "Empty propagation path.")

  ray_path_point = *std::ranges::max_element(
      ray_path,
      [](const PropagationPathPoint& a, const PropagationPathPoint& b) {
        return a.altitude() < b.altitude();
      });
}

void ray_pathGeometricUplooking(ArrayOfPropagationPathPoint& ray_path,
                                const AtmField& atmospheric_field,
                                const SurfaceField& surface_field,
                                const Numeric& latitude,
                                const Numeric& longitude,
                                const Numeric& max_step) {
  ray_pathGeometric(
      ray_path,
      atmospheric_field,
      surface_field,
      {surface_field.single_value(SurfaceKey::h, latitude, longitude),
       latitude,
       longitude},
      {0, 0},
      max_step,
      1.0,
      true,
      false,
      true,
      true,
      false);
}

void ray_pathGeometricDownlooking(ArrayOfPropagationPathPoint& ray_path,
                                  const AtmField& atmospheric_field,
                                  const SurfaceField& surface_field,
                                  const Numeric& latitude,
                                  const Numeric& longitude,
                                  const Numeric& max_step) {
  ray_pathGeometric(ray_path,
                    atmospheric_field,
                    surface_field,
                    {atmospheric_field.top_of_atmosphere, latitude, longitude},
                    {180, 0},
                    max_step,
                    1.0,
                    true,
                    false,
                    true,
                    true,
                    false);
}
