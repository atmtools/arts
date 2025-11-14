#include <workspace.h>

#include <algorithm>

#include "atm_field.h"
#include "enumsSurfaceKey.h"
#include "geodetic.h"
#include "path_point.h"

void ray_pathInit(ArrayOfPropagationPathPoint& ray_path,
                  const AtmField& atm_field,
                  const SurfaceField& surf_field,
                  const Vector3& pos,
                  const Vector2& los,
                  const Index& as_sensor) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  ARTS_USER_ERROR_IF(any_nan(pos) or any_nan(los),
                     R"(There are NAN in the pos or los vector:
pos:      {:B,}
los:      {:B,}
)",
                     pos,
                     los);

  ray_path.resize(1);
  ray_path[0] =
      path::init(pos, los, atm_field, surf_field, static_cast<bool>(as_sensor));
}

void ray_pathRemoveNearby(ArrayOfPropagationPathPoint& ray_path,
                          const SurfaceField& surf_field,
                          const Numeric& min_dist,
                          const Index& first) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  path::erase_closeby(ray_path, surf_field, min_dist, first);
}

void ray_pathSetGeometricExtremes(ArrayOfPropagationPathPoint& ray_path,
                                  const AtmField& atm_field,
                                  const SurfaceField& surf_field,
                                  const Numeric& surf_search_accuracy,
                                  const Index& surf_safe_search) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  path::set_geometric_extremes(ray_path,
                               atm_field,
                               surf_field,
                               surf_search_accuracy,
                               static_cast<bool>(surf_safe_search));
}

void ray_pathRemoveNonGeometricGridCrossings(
    ArrayOfPropagationPathPoint& ray_path,
    const AtmField& atm_field,
    const AtmKey& atm_key) {
  ARTS_TIME_REPORT

  const auto& data = atm_field[atm_key];
  ARTS_USER_ERROR_IF(not std::holds_alternative<GeodeticField3>(data.data),
                     "The data for key {} is not a GeodeticField3",
                     atm_key);
  const auto& field = std::get<GeodeticField3>(data.data);
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
                                       const AtmField& atm_field,
                                       const SurfaceField& surf_field,
                                       const AtmKey& atm_key) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  const auto& data = atm_field[atm_key];
  ARTS_USER_ERROR_IF(not std::holds_alternative<GeodeticField3>(data.data),
                     "The data for key {} is not a GeodeticField3",
                     atm_key);
  const auto& field = std::get<GeodeticField3>(data.data);
  const auto& alt   = field.grid<0>();
  const auto& lat   = field.grid<1>();
  const auto& lon   = field.grid<2>();

  path::fill_geometric_crossings(ray_path, surf_field, alt, lat, lon);
}

void ray_pathFillGeometricHalfStep(ArrayOfPropagationPathPoint& ray_path,
                                   const SurfaceField& surf_field,
                                   const Numeric& max_step) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  path::fill_geometric_by_half_steps(ray_path, surf_field, max_step);
}

void ray_pathFillGeometricStepwise(ArrayOfPropagationPathPoint& ray_path,
                                   const SurfaceField& surf_field,
                                   const Numeric& max_step) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  path::fill_geometric_stepwise(ray_path, surf_field, max_step);
}

void ray_pathFixUpdownAzimuth(ArrayOfPropagationPathPoint& ray_path) {
  ARTS_TIME_REPORT

  path::fix_updown_azimuth_to_first(ray_path);
}

void ray_pathAddLimbPoint(ArrayOfPropagationPathPoint& ray_path,
                          const SurfaceField& surf_field) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  path::fill_geometric_limb(ray_path, surf_field);
}

void ray_pathRemoveNonAtm(ArrayOfPropagationPathPoint& ray_path) {
  ARTS_TIME_REPORT

  path::keep_only_atm(ray_path);
}

void ray_path_observer_agendaSetGeometric(Agenda& ray_path_observer_agenda,
                                          const String& max_step_option,
                                          const Numeric& surf_search_accuracy,
                                          const Numeric& remove_nearby,
                                          const AtmKey& atm_key,
                                          const Index& surf_safe_search,
                                          const Index& remove_nearby_first,
                                          const Index& add_crossings,
                                          const Index& remove_non_crossings,
                                          const Index& fix_updown_azimuth,
                                          const Index& add_limb,
                                          const Index& remove_non_atm) {
  ARTS_TIME_REPORT

  AgendaCreator creator("ray_path_observer_agenda");

  creator.add("ray_pathInit",
              SetWsv{"pos", "spectral_radiance_observer_position"},
              SetWsv{"los", "spectral_radiance_observer_line_of_sight"},
              SetWsv("as_sensor", Index{1}));

  creator.add("ray_pathSetGeometricExtremes",
              SetWsv("surf_search_accuracy", surf_search_accuracy),
              SetWsv("surf_safe_search", surf_safe_search));

  if (add_crossings) {
    creator.add("ray_pathAddGeometricGridCrossings",
                SetWsv("atm_key", atm_key));
  }

  switch (to<ray_path_observer_agendaSetGeometricMaxStep>(max_step_option)) {
    case ray_path_observer_agendaSetGeometricMaxStep::half:
      creator.add("ray_pathFillGeometricMaxStep");
      break;
    case ray_path_observer_agendaSetGeometricMaxStep::step:
      creator.add("ray_pathFillGeometricStepwise");
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
                       const AtmField& atm_field,
                       const SurfaceField& surf_field,
                       const Numeric& max_step,
                       const Vector3& pos,
                       const Vector2& los,
                       const Numeric& surf_search_accuracy,
                       const Index& as_sensor,
                       const Index& add_limb,
                       const Index& remove_non_atm,
                       const Index& fix_updown_azimuth,
                       const Index& surf_safe_search) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  ray_pathInit(ray_path, atm_field, surf_field, pos, los, as_sensor);
  ray_pathSetGeometricExtremes(
      ray_path, atm_field, surf_field, surf_search_accuracy, surf_safe_search);
  ray_pathFillGeometricStepwise(ray_path, surf_field, max_step);
  if (add_limb) ray_pathAddLimbPoint(ray_path, surf_field);
  if (fix_updown_azimuth) ray_pathFixUpdownAzimuth(ray_path);
  if (remove_non_atm) ray_pathRemoveNonAtm(ray_path);
}

void ray_path_pointBackground(PropagationPathPoint& ray_path_point,
                              const ArrayOfPropagationPathPoint& ray_path) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(ray_path.size() == 0, "Empty propagation path.")
  ray_path_point = ray_path.back();
}

void ray_path_pointForeground(PropagationPathPoint& ray_path_point,
                              const ArrayOfPropagationPathPoint& ray_path) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(ray_path.size() == 0, "Empty propagation path.")
  ray_path_point = ray_path.front();
}

void ray_path_pointLowestFromPath(PropagationPathPoint& ray_path_point,
                                  const ArrayOfPropagationPathPoint& ray_path) {
  ARTS_TIME_REPORT

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
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(ray_path.size() == 0, "Empty propagation path.")

  ray_path_point = *std::ranges::max_element(
      ray_path,
      [](const PropagationPathPoint& a, const PropagationPathPoint& b) {
        return a.altitude() < b.altitude();
      });
}

void ray_pathGeometricUplooking(ArrayOfPropagationPathPoint& ray_path,
                                const AtmField& atm_field,
                                const SurfaceField& surf_field,
                                const Numeric& latitude,
                                const Numeric& longitude,
                                const Numeric& max_step) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  ray_pathGeometric(
      ray_path,
      atm_field,
      surf_field,
      max_step,
      {surf_field.single_value(SurfaceKey::h, latitude, longitude),
       latitude,
       longitude},
      {0, 0},
      1.0,
      true,
      false,
      true,
      true,
      false);
}

void ray_pathGeometricDownlooking(ArrayOfPropagationPathPoint& ray_path,
                                  const AtmField& atm_field,
                                  const SurfaceField& surf_field,
                                  const Numeric& latitude,
                                  const Numeric& longitude,
                                  const Numeric& max_step) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  ray_pathGeometric(ray_path,
                    atm_field,
                    surf_field,
                    max_step,
                    {atm_field.top_of_atmosphere, latitude, longitude},
                    {180, 0},
                    1.0,
                    true,
                    false,
                    true,
                    true,
                    false);
}

void ray_path_pointPastGeometric(PropagationPathPoint& ray_path_point,
                                 const ArrayOfPropagationPathPoint& ray_path,
                                 const AtmField& atm_field,
                                 const SurfaceField& surf_field,
                                 const Numeric& max_stepsize,
                                 const Numeric& safe_search_accuracy,
                                 const Index& search_safe) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid);
  ARTS_USER_ERROR_IF(ray_path.size() == 0, "Empty propagation path.");

  ray_path_point = past_geometric(ray_path.back(),
                                  atm_field,
                                  surf_field,
                                  max_stepsize,
                                  safe_search_accuracy,
                                  search_safe);
}

void ray_path_pointPastRefractive(PropagationPathPoint& ray_path_point,
                                  const ArrayOfPropagationPathPoint& ray_path,
                                  const AtmField& atm_field,
                                  const SurfaceField& surf_field,
                                  const Numeric& max_stepsize,
                                  const Numeric& single_dispersion,
                                  const Numeric& safe_search_accuracy,
                                  const Index& search_safe) {
  ARTS_TIME_REPORT

  using Conversion::cosd, Conversion::acosd;

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid);
  ARTS_USER_ERROR_IF(ray_path.size() == 0, "Empty propagation path.");

  PropagationPathPoint future = ray_path.back();
  future.zenith() =
      acosd((future.nreal / (1.0 + single_dispersion)) * cosd(future.zenith()));

  future.nreal = 1 + single_dispersion;

  ray_path_point = past_geometric(future,
                                  atm_field,
                                  surf_field,
                                  max_stepsize,
                                  safe_search_accuracy,
                                  search_safe);
}
