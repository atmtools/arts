#include <arts_omp.h>
#include <workspace.h>

namespace {
Numeric surface_tangent_zenithFromAgenda(const Workspace& ws,
                                         const Vector3& pos,
                                         const Agenda& ray_path_observer_agenda,
                                         const Numeric azimuth) {
  ARTS_TIME_REPORT

  ArrayOfPropagationPathPoint ray_path;
  Numeric za0 = 180.0, za1 = 90.0;
  while (std::nextafter(za0, za1) != za1) {
    const Numeric za = std::midpoint(za0, za1);
    ray_path_observer_agendaExecute(
        ws, ray_path, pos, {za, azimuth}, ray_path_observer_agenda);

    // za0 points at surface, za1 points at space
    ((ray_path.back().los_type == PathPositionType::surface) ? za0 : za1) = za;
  }

  return za1;
}
}  // namespace

void ray_path_observersFieldProfilePseudo2D(
    const Workspace& ws,
    ArrayOfPropagationPathPoint& ray_path_observers,
    const AtmField& atmospheric_field,
    const SurfaceField& surface_field,
    const Agenda& ray_path_observer_agenda,
    const Numeric& latitude,
    const Numeric& longitude,
    const Numeric& azimuth,
    const Index& nup,
    const Index& nlimb,
    const Index& ndown) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(nup < 2 or nlimb < 2 or ndown < 2,
                     "Must have at least 2 observers per meta-direction.")

  const Vector3 top_pos = {
      atmospheric_field.top_of_atmosphere, latitude, longitude};
  const Vector3 bot_pos = {surface_field[SurfaceKey::h].at(latitude, longitude),
                           latitude,
                           longitude};

  const Numeric za_surface_limb = surface_tangent_zenithFromAgenda(
      ws, top_pos, ray_path_observer_agenda, azimuth);

  const Index N = nlimb + nup + ndown;
  arr::resize(0, ray_path_observers);
  arr::reserve(N, ray_path_observers);

  const Numeric miss_surf = std::nextafter(za_surface_limb, 0.0);
  const Numeric hit_surf  = std::nextafter(za_surface_limb, 180.0);
  const AscendingGrid upp = nlinspace(0.0, 90.0, nup);
  const AscendingGrid lmb = nlinspace(90.0, miss_surf, nlimb + 1);
  const AscendingGrid dwn = nlinspace(hit_surf, 180.0, ndown);

  for (auto za : upp) {
    ray_path_observers.push_back({.pos_type = PathPositionType::surface,
                                  .los_type = PathPositionType::space,
                                  .pos      = bot_pos,
                                  .los      = {za, azimuth}});
  }

  // Skip 90.0
  for (auto za : lmb | stdv::drop(1)) {
    ray_path_observers.push_back({.pos_type = PathPositionType::atm,
                                  .los_type = PathPositionType::space,
                                  .pos      = top_pos,
                                  .los      = {za, azimuth}});
  }

  for (auto za : dwn) {
    ray_path_observers.push_back({.pos_type = PathPositionType::atm,
                                  .los_type = PathPositionType::surface,
                                  .pos      = top_pos,
                                  .los      = {za, azimuth}});
  }
}

void ray_path_fieldFromObserverAgenda(
    const Workspace& ws,
    ArrayOfArrayOfPropagationPathPoint& ray_path_field,
    const ArrayOfPropagationPathPoint& ray_path_observers,
    const Agenda& ray_path_observer_agenda) {
  ARTS_TIME_REPORT

  const Size N = ray_path_observers.size();
  ARTS_USER_ERROR_IF(N == 0, "Must have at least one observer.")

  ray_path_field.resize(N);

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < N; i++) {
    try {
      ray_path_observer_agendaExecute(ws,
                                      ray_path_field[i],
                                      ray_path_observers[i].pos,
                                      ray_path_observers[i].los,
                                      ray_path_observer_agenda);
    } catch (const std::exception& e) {
#pragma omp critical
      error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(not error.empty(), "{}", error);
}
