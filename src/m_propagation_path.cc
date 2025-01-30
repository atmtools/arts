#include <path_point.h>

#include <algorithm>
#include <functional>

#include "atm.h"
#include "auto_wsa_options.h"
#include "surf.h"
#include "workspace_agenda_creator.h"
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

void ray_pathInit(ArrayOfPropagationPathPoint& ray_path,
                  const AtmField& atmospheric_field,
                  const SurfaceField& surface_field,
                  const Vector3& pos,
                  const Vector2& los,
                  const Index& as_sensor) {
  ARTS_USER_ERROR_IF(any_nan(pos) or any_nan(los),
                     R"(There are NAN in the pos or los vector:
pos:      {:B,}
los:      {:B,}
)",
                     pos,
                     los);

  ray_path.resize(1);
  ray_path[0] = path::init(
      pos, los, atmospheric_field, surface_field, static_cast<bool>(as_sensor));
}

void ray_pathRemoveNearby(ArrayOfPropagationPathPoint& ray_path,
                          const SurfaceField& surface_field,
                          const Numeric& min_dist,
                          const Index& first) {
  path::erase_closeby(ray_path, surface_field, min_dist, first);
}

void ray_pathSetGeometricExtremes(ArrayOfPropagationPathPoint& ray_path,
                                  const AtmField& atmospheric_field,
                                  const SurfaceField& surface_field,
                                  const Numeric& surface_search_accuracy,
                                  const Index& surface_safe_search) {
  path::set_geometric_extremes(ray_path,
                               atmospheric_field,
                               surface_field,
                               surface_search_accuracy,
                               static_cast<bool>(surface_safe_search));
}

void ray_pathRemoveNonGeometricGridCrossings(
    ArrayOfPropagationPathPoint& ray_path,
    const AtmField& atmospheric_field,
    const AtmKey& atm_key) {
  const auto& data = atmospheric_field[atm_key];
  ARTS_USER_ERROR_IF(not std::holds_alternative<GriddedField3>(data.data),
                     "The data for key {} is not a GriddedField3",
                     atm_key);
  const auto& field = std::get<GriddedField3>(data.data);
  const auto& alt   = field.grid<0>();
  const auto& lat   = field.grid<1>();
  const auto& lon   = field.grid<2>();

  std::erase_if(ray_path, [&alt, &lat, &lon](auto& p) {
    return not(stdr::contains(alt, p.altitude()) or
               stdr::contains(lat, p.latitude()) or
               stdr::contains(lon, p.longitude()));
  });
}

void ray_pathAddGeometricGridCrossings(ArrayOfPropagationPathPoint& ray_path,
                                       const AtmField& atmospheric_field,
                                       const SurfaceField& surface_field,
                                       const AtmKey& atm_key) {
  const auto& data = atmospheric_field[atm_key];
  ARTS_USER_ERROR_IF(not std::holds_alternative<GriddedField3>(data.data),
                     "The data for key {} is not a GriddedField3",
                     atm_key);
  const auto& field = std::get<GriddedField3>(data.data);
  const auto& alt   = field.grid<0>();
  const auto& lat   = field.grid<1>();
  const auto& lon   = field.grid<2>();

  path::fill_geometric_crossings(ray_path, surface_field, alt, lat, lon);
}

void ray_pathFillGeometricHalfStep(ArrayOfPropagationPathPoint& ray_path,
                                   const SurfaceField& surface_field,
                                   const Numeric& max_step) {
  path::fill_geometric_by_half_steps(ray_path, surface_field, max_step);
}

void ray_pathFillGeometricStepwise(ArrayOfPropagationPathPoint& ray_path,
                                   const SurfaceField& surface_field,
                                   const Numeric& max_step) {
  path::fill_geometric_stepwise(ray_path, surface_field, max_step);
}

void ray_pathFixUpdownAzimuth(ArrayOfPropagationPathPoint& ray_path) {
  path::fix_updown_azimuth_to_first(ray_path);
}

void ray_pathAddLimbPoint(ArrayOfPropagationPathPoint& ray_path,
                          const SurfaceField& surface_field) {
  path::fill_geometric_limb(ray_path, surface_field);
}

void ray_pathRemoveNonAtm(ArrayOfPropagationPathPoint& ray_path) {
  path::keep_only_atm(ray_path);
}

void ray_path_observer_agendaSetGeometric(
    Agenda& ray_path_observer_agenda,
    const String& max_step_option,
    const Numeric& surface_search_accuracy,
    const Numeric& max_step,
    const Numeric& remove_nearby,
    const AtmKey& atm_key,
    const Index& surface_safe_search,
    const Index& remove_nearby_first,
    const Index& add_crossings,
    const Index& remove_non_crossings,
    const Index& fix_updown_azimuth,
    const Index& add_limb,
    const Index& remove_non_atm) {
  AgendaCreator creator("ray_path_observer_agenda");

  creator.add("ray_pathInit",
              SetWsv{"pos", "spectral_radiance_observer_position"},
              SetWsv{"los", "spectral_radiance_observer_line_of_sight"},
              SetWsv("as_sensor", Index{1}));

  creator.add("ray_pathSetGeometricExtremes",
              SetWsv("surface_search_accuracy", surface_search_accuracy),
              SetWsv("surface_safe_search", surface_safe_search));

  if (add_crossings) {
    creator.add("ray_pathAddGeometricGridCrossings",
                SetWsv("atm_key", atm_key));
  }

  switch (to<ray_path_observer_agendaSetGeometricMaxStep>(max_step_option)) {
    case ray_path_observer_agendaSetGeometricMaxStep::half:
      creator.add("ray_pathFillGeometricMaxStep", SetWsv("max_step", max_step));
      break;
    case ray_path_observer_agendaSetGeometricMaxStep::step:
      creator.add("ray_pathFillGeometricStepwise",
                  SetWsv("max_step", max_step));
      break;
    case ray_path_observer_agendaSetGeometricMaxStep::None: break;
  }

  if (remove_non_crossings) {
    creator.add("ray_pathRemoveNonGeometricGridCrossings",
                SetWsv("atm_key", atm_key));
  }

  if (fix_updown_azimuth) creator.add("ray_pathFixUpdownAzimuth");

  if (add_limb) creator.add("ray_pathAddLimbPoint");

  if (remove_non_atm) creator.add("ray_pathRemoveNonAtm");

  if (remove_nearby > 0.0) {
    creator.add("ray_pathRemoveNearby",
                SetWsv("min_dist", remove_nearby),
                SetWsv("first", remove_nearby_first));
  }

  ray_path_observer_agenda = std::move(creator).finalize(false);
}

void ray_pathGeometric(ArrayOfPropagationPathPoint& ray_path,
                       const AtmField& atmospheric_field,
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
  ray_pathInit(ray_path, atmospheric_field, surface_field, pos, los, as_sensor);
  ray_pathSetGeometricExtremes(ray_path,
                               atmospheric_field,
                               surface_field,
                               surface_search_accuracy,
                               surface_safe_search);
  ray_pathFillGeometricStepwise(ray_path, surface_field, max_step);
  if (add_limb) ray_pathAddLimbPoint(ray_path, surface_field);
  if (fix_updown_azimuth) ray_pathFixUpdownAzimuth(ray_path);
  if (remove_non_atm) ray_pathRemoveNonAtm(ray_path);
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
