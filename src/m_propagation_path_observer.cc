#include <array_algo.h>
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
    const AtmField& atm_field,
    const SurfaceField& surface_field,
    const Agenda& ray_path_observer_agenda,
    const Numeric& latitude,
    const Numeric& longitude,
    const Numeric& azimuth,
    const Index& nup,
    const Index& nlimb,
    const Index& ndown) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surface_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surface_field.ellipsoid)

  ARTS_USER_ERROR_IF(nup < 2 or nlimb < 2 or ndown < 2,
                     "Must have at least 2 observers per meta-direction.")

  const Vector3 top_pos = {atm_field.top_of_atmosphere, latitude, longitude};
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
                                  .los_type = PathPositionType::atm,
                                  .pos      = bot_pos,
                                  .los      = {za, azimuth}});
  }

  // Skip 90.0
  for (auto za : lmb | stdv::drop(1)) {
    ray_path_observers.push_back({.pos_type = PathPositionType::atm,
                                  .los_type = PathPositionType::atm,
                                  .pos      = top_pos,
                                  .los      = {za, azimuth}});
  }

  for (auto za : dwn) {
    ray_path_observers.push_back({.pos_type = PathPositionType::atm,
                                  .los_type = PathPositionType::atm,
                                  .pos      = top_pos,
                                  .los      = {za, azimuth}});
  }
}

void ray_path_observersFluxProfile(
    ArrayOfPropagationPathPoint& ray_path_observers,
    const AtmField& atm_field,
    const Numeric& azimuth,
    const Index& n,
    const AtmKey& atm_key) {
  ray_path_observers.clear();

  ARTS_USER_ERROR_IF(
      n < 3 or (n % 2) == 0,
      "Must have at least 3 observers, and an uneven number of them.")

  const auto& data = atm_field[atm_key].get<GeodeticField3>();

  const auto& alt_g  = data.grid<0>();
  const auto& lat    = data.grid<1>()[0];
  const auto& lon    = data.grid<2>()[0];
  const Vector zas_g = nlinspace(0.0, 180.0, n + 2);

  ARTS_USER_ERROR_IF(data.data.size() != alt_g.size(),
                     "Data size does not match altitude grid size")

  ray_path_observers.reserve(n * alt_g.size() + 2);
  for (Index i = 0; i < n; i++) {
    for (PropagationPathPoint p :
         alt_g | stdv::transform([za = zas_g[i + 1], lat, lon, azimuth](
                                     Numeric alt) {
           return PropagationPathPoint{.pos_type = PathPositionType::atm,
                                       .los_type = PathPositionType::atm,
                                       .pos      = {alt, lat, lon},
                                       .los      = {za, azimuth}};
         })) {
      ray_path_observers.push_back(p);
    }
  }
  ray_path_observers.push_back({.pos_type = PathPositionType::atm,
                                .los_type = PathPositionType::atm,
                                .pos      = {alt_g.front(), lat, lon},
                                .los      = {0.0, azimuth}});
  ray_path_observers.push_back({.pos_type = PathPositionType::atm,
                                .los_type = PathPositionType::atm,
                                .pos      = {alt_g.back(), lat, lon},
                                .los      = {180.0, azimuth}});
}

namespace {
Vector half_grid(const Numeric x0, const Numeric x1, const Numeric dx) {
  assert(dx >= 0.0);

  Numeric d = x1 - x0;

  if (d < dx) return {x0, x1};

  d      *= 0.5;
  Size s = 1, dn = 1;
  while (d > dx) {
    d  *= 0.5;
    dn *= 2;
    s  += dn;
  }

  return nlinspace(x0, x1, s + 2);
}
}  // namespace

void ray_path_fieldFluxProfile(
    const Workspace& ws,
    ArrayOfArrayOfPropagationPathPoint& ray_path_field,
    const AtmField& atm_field,
    const Agenda& ray_path_observer_agenda,
    const Numeric& azimuth,
    const Numeric& dza,
    const AtmKey& atm_key) try {
  ARTS_USER_ERROR_IF(dza <= 0.0, "Zenith angle step must be positive")

  ray_path_field.clear();

  const auto& data  = atm_field[atm_key].get<GeodeticField3>();
  const auto& alt_g = data.grid<0>();
  ARTS_USER_ERROR_IF(data.data.size() != alt_g.size(),
                     "Data size does not match altitude grid size")
  const auto& lat = data.grid<1>()[0];
  const auto& lon = data.grid<2>()[0];

  const Numeric za_limb = surface_tangent_zenithFromAgenda(
      ws, {alt_g.back(), lat, lon}, ray_path_observer_agenda, azimuth);
  const Numeric za_limb_miss = std::nextafter(za_limb, 0.0);
  const Numeric za_limb_hit  = std::nextafter(za_limb, 180.0);

  const Vector looking_down = half_grid(za_limb_hit, 180.0, dza);
  const Vector looking_up   = half_grid(0.0, 90, dza);

  // 4 extra points already added by
  ray_path_field.resize(looking_down.size() + looking_up.size());

  ray_path_observer_agendaExecute(ws,
                                  ray_path_field.front(),
                                  {alt_g.front(), lat, lon},
                                  {0, azimuth},
                                  ray_path_observer_agenda);

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 1; i < looking_up.size() - 1; i++) {
    try {
      ray_path_observer_agendaExecute(ws,
                                      ray_path_field[i],
                                      {alt_g.front(), lat, lon},
                                      {looking_up[i], azimuth},
                                      ray_path_observer_agenda);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }
  if (not error.empty()) throw std::runtime_error(error);

  ray_path_observer_agendaExecute(ws,
                                  ray_path_field[looking_up.size() - 1],
                                  {alt_g.back(), lat, lon},
                                  {za_limb_miss, azimuth},
                                  ray_path_observer_agenda);

  ray_path_observer_agendaExecute(ws,
                                  ray_path_field[looking_up.size()],
                                  {alt_g.back(), lat, lon},
                                  {za_limb_hit, azimuth},
                                  ray_path_observer_agenda);

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 1; i < looking_down.size() - 1; i++) {
    try {
      ray_path_observer_agendaExecute(ws,
                                      ray_path_field[looking_up.size() + i],
                                      {alt_g.back(), lat, lon},
                                      {looking_down[i], azimuth},
                                      ray_path_observer_agenda);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }
  if (not error.empty()) throw std::runtime_error(error);

  ray_path_observer_agendaExecute(ws,
                                  ray_path_field.back(),
                                  {alt_g.back(), lat, lon},
                                  {180, azimuth},
                                  ray_path_observer_agenda);

  std::erase_if(ray_path_field, [](const auto& x) { return x.size() < 2; });
}
ARTS_METHOD_ERROR_CATCH

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

  if (not error.empty()) throw std::runtime_error(error);
}
