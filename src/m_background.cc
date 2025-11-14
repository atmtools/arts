#include <debug.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <sun.h>
#include <surf.h>
#include <workspace.h>

#include <algorithm>

void spectral_radiance_backgroundAgendasAtEndOfPath(
    const Workspace& ws,
    StokvecVector& spectral_radiance_background,
    StokvecMatrix& spectral_radiance_background_jacobian,
    const AscendingGrid& freq_grid,
    const JacobianTargets& jacobian_targets,
    const PropagationPathPoint& ray_path_point,
    const SurfaceField& surface_field,
    const SubsurfaceField& subsurface_field,
    const Agenda& spectral_radiance_space_agenda,
    const Agenda& spectral_radiance_surface_agenda) try {
  ARTS_TIME_REPORT

  using enum PathPositionType;
  switch (ray_path_point.los_type) {
    case atm:
      ARTS_USER_ERROR("Undefined what to do with an atmospheric background")
      break;
    case unknown: ARTS_USER_ERROR("Undefined background type"); break;
    case space:
      spectral_radiance_space_agendaExecute(
          ws,
          spectral_radiance_background,
          spectral_radiance_background_jacobian,
          freq_grid,
          jacobian_targets,
          ray_path_point,
          spectral_radiance_space_agenda);
      break;
    case surface:
      spectral_radiance_surface_agendaExecute(
          ws,
          spectral_radiance_background,
          spectral_radiance_background_jacobian,
          freq_grid,
          jacobian_targets,
          ray_path_point,
          surface_field,
          subsurface_field,
          spectral_radiance_surface_agenda);
      break;
  }
}
ARTS_METHOD_ERROR_CATCH

namespace {
StokvecVector from_temp(const ConstVectorView& freq_grid, const Numeric t) {
  StokvecVector v(freq_grid.size(), 0.0);
  std::transform(
      freq_grid.begin(), freq_grid.end(), v.begin(), [t](auto f) -> Stokvec {
        return {planck(f, t), 0, 0, 0};
      });
  return v;
}
}  // namespace

void spectral_radianceUniformCosmicBackground(StokvecVector& spectral_radiance,
                                              const AscendingGrid& freq_grid) {
  ARTS_TIME_REPORT

  constexpr auto t  = Constant::cosmic_microwave_background_temperature;
  spectral_radiance = from_temp(freq_grid, t);
}

void spectral_radianceSunOrCosmicBackground(
    StokvecVector& spectral_radiance,
    const AscendingGrid& freq_grid,
    const ArrayOfPropagationPathPoint& sun_path,
    const Sun& sun,
    const SurfaceField& surface_field) try {
  ARTS_TIME_REPORT

  spectral_radiance.resize(freq_grid.size());

  if (set_spectral_radiance_if_sun_intersection(
          spectral_radiance, sun, sun_path.back(), surface_field))
    return;

  spectral_radianceUniformCosmicBackground(spectral_radiance, freq_grid);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceSunsOrCosmicBackground(
    StokvecVector& spectral_radiance,
    const AscendingGrid& freq_grid,
    const PropagationPathPoint& ray_path_point,
    const ArrayOfSun& suns,
    const SurfaceField& surface_field) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surface_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surface_field.ellipsoid)

  spectral_radiance.resize(freq_grid.size());

  for (auto& sun : suns) {
    if (set_spectral_radiance_if_sun_intersection(
            spectral_radiance, sun, ray_path_point, surface_field)) {
      return;
    }
  }

  spectral_radianceUniformCosmicBackground(spectral_radiance, freq_grid);
}

void spectral_radianceSurfaceBlackbody(
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_jacobian,
    const AscendingGrid& freq_grid,
    const SurfaceField& surface_field,
    const JacobianTargets& jacobian_targets,
    const PropagationPathPoint& ray_path_point) {
  ARTS_TIME_REPORT

  constexpr auto key = SurfaceKey::t;

  ARTS_USER_ERROR_IF(not surface_field.contains(key),
                     "Surface field does not contain temperature")

  const auto t = surface_field.single_value(
      key, ray_path_point.pos[1], ray_path_point.pos[2]);
  spectral_radiance = from_temp(freq_grid, t);

  spectral_radiance_jacobianEmpty(
      spectral_radiance_jacobian, freq_grid, jacobian_targets);

  for (auto& target : jacobian_targets.surf) {
    if (target.type == SurfaceKey::t) {
      const auto& data = surface_field[target.type];

      const auto ws =
          data.flat_weights(ray_path_point.pos[1], ray_path_point.pos[2]);

      for (Size i = 0; i < freq_grid.size(); i++) {
        const Numeric dBdt = dplanck_dt(freq_grid[i], t);
        for (const auto& w : ws) {
          spectral_radiance_jacobian[w.first + target.x_start, i] +=
              w.second * dBdt;
        }
      }
    }
  }
}

void transmission_matrix_backgroundFromPathPropagationBack(
    MuelmatVector& transmission_matrix_background,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix_cumulative) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(ray_path_transmission_matrix_cumulative.size() == 0,
                     "Cannot extract from empty list.")
  transmission_matrix_background =
      ray_path_transmission_matrix_cumulative.back();
}
ARTS_METHOD_ERROR_CATCH

void transmission_matrix_backgroundFromPathPropagationFront(
    MuelmatVector& transmission_matrix_background,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix_cumulative) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(ray_path_transmission_matrix_cumulative.size() == 0,
                     "Cannot extract from empty list.")
  transmission_matrix_background =
      ray_path_transmission_matrix_cumulative.front();
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceDefaultTransmission(
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_background,
    const AscendingGrid& freq_grid,
    const JacobianTargets& jacobian_targets) {
  ARTS_TIME_REPORT

  const Index nf = freq_grid.size();
  const Index nq = jacobian_targets.x_size();

  spectral_radiance_background.resize(nq, nf);
  spectral_radiance_background = 0.0;

  spectral_radiance.resize(nf);
  spectral_radiance = 1;
}

void single_rad_backgroundAgendasAtEndOfPath(
    const Workspace& ws,
    Stokvec& single_rad_background,
    StokvecVector& single_rad_background_jacobian,
    const Numeric& frequency,
    const JacobianTargets& jacobian_targets,
    const PropagationPathPoint& ray_path_point,
    const SurfaceField& surface_field,
    const SubsurfaceField& subsurface_field,
    const Agenda& single_rad_space_agenda,
    const Agenda& single_rad_surface_agenda) try {
  ARTS_TIME_REPORT

  using enum PathPositionType;
  switch (ray_path_point.los_type) {
    case atm:
      ARTS_USER_ERROR("Undefined what to do with an atmospheric background")
      break;
    case unknown: ARTS_USER_ERROR("Undefined background type"); break;
    case space:
      single_rad_space_agendaExecute(ws,
                                     single_rad_background,
                                     single_rad_background_jacobian,
                                     frequency,
                                     jacobian_targets,
                                     ray_path_point,
                                     single_rad_space_agenda);
      break;
    case surface:
      single_rad_surface_agendaExecute(ws,
                                       single_rad_background,
                                       single_rad_background_jacobian,
                                       frequency,
                                       jacobian_targets,
                                       ray_path_point,
                                       surface_field,
                                       subsurface_field,
                                       single_rad_surface_agenda);
      break;
  }
}
ARTS_METHOD_ERROR_CATCH
