#include <debug.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <sun.h>
#include <surf.h>
#include <workspace.h>

#include <algorithm>

void spectral_rad_bkgAgendasAtEndOfPath(
    const Workspace& ws,
    StokvecVector& spectral_rad_bkg,
    StokvecMatrix& spectral_rad_bkg_jac,
    const AscendingGrid& freq_grid,
    const JacobianTargets& jac_targets,
    const PropagationPathPoint& ray_path_point,
    const SurfaceField& surf_field,
    const SubsurfaceField& subsurf_field,
    const Agenda& spectral_rad_space_agenda,
    const Agenda& spectral_rad_surface_agenda) try {
  ARTS_TIME_REPORT

  using enum PathPositionType;
  switch (ray_path_point.los_type) {
    case atm:
      ARTS_USER_ERROR("Undefined what to do with an atmospheric background")
      break;
    case unknown: ARTS_USER_ERROR("Undefined background type"); break;
    case space:
      spectral_rad_space_agendaExecute(ws,
                                       spectral_rad_bkg,
                                       spectral_rad_bkg_jac,
                                       freq_grid,
                                       jac_targets,
                                       ray_path_point,
                                       spectral_rad_space_agenda);
      break;
    case surface:
      spectral_rad_surface_agendaExecute(ws,
                                         spectral_rad_bkg,
                                         spectral_rad_bkg_jac,
                                         freq_grid,
                                         jac_targets,
                                         ray_path_point,
                                         surf_field,
                                         subsurf_field,
                                         spectral_rad_surface_agenda);
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

void spectral_radUniformCosmicBackground(StokvecVector& spectral_rad,
                                         const AscendingGrid& freq_grid) {
  ARTS_TIME_REPORT

  constexpr auto t = Constant::cosmic_microwave_background_temperature;
  spectral_rad     = from_temp(freq_grid, t);
}

void spectral_radSunOrCosmicBackground(
    StokvecVector& spectral_rad,
    const AscendingGrid& freq_grid,
    const ArrayOfPropagationPathPoint& sun_path,
    const Sun& sun,
    const SurfaceField& surf_field) try {
  ARTS_TIME_REPORT

  spectral_rad.resize(freq_grid.size());

  if (set_spectral_rad_if_sun_intersection(
          spectral_rad, sun, sun_path.back(), surf_field))
    return;

  spectral_radUniformCosmicBackground(spectral_rad, freq_grid);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radSunsOrCosmicBackground(
    StokvecVector& spectral_rad,
    const AscendingGrid& freq_grid,
    const PropagationPathPoint& ray_path_point,
    const ArrayOfSun& suns,
    const SurfaceField& surf_field) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  spectral_rad.resize(freq_grid.size());

  for (auto& sun : suns) {
    if (set_spectral_rad_if_sun_intersection(
            spectral_rad, sun, ray_path_point, surf_field)) {
      return;
    }
  }

  spectral_radUniformCosmicBackground(spectral_rad, freq_grid);
}

void spectral_radSurfaceBlackbody(StokvecVector& spectral_rad,
                                  StokvecMatrix& spectral_rad_jac,
                                  const AscendingGrid& freq_grid,
                                  const SurfaceField& surf_field,
                                  const JacobianTargets& jac_targets,
                                  const PropagationPathPoint& ray_path_point) {
  ARTS_TIME_REPORT

  constexpr auto key = SurfaceKey::t;

  ARTS_USER_ERROR_IF(not surf_field.contains(key),
                     "Surface field does not contain temperature")

  const auto t = surf_field.single_value(
      key, ray_path_point.pos[1], ray_path_point.pos[2]);
  spectral_rad = from_temp(freq_grid, t);

  spectral_rad_jacEmpty(spectral_rad_jac, freq_grid, jac_targets);

  for (auto& target : jac_targets.surf) {
    if (target.type == SurfaceKey::t) {
      const auto& data = surf_field[target.type];

      const auto ws =
          data.flat_weights(ray_path_point.pos[1], ray_path_point.pos[2]);

      for (Size i = 0; i < freq_grid.size(); i++) {
        const Numeric dBdt = dplanck_dt(freq_grid[i], t);
        for (const auto& w : ws) {
          spectral_rad_jac[w.first + target.x_start, i] += w.second * dBdt;
        }
      }
    }
  }
}

void spectral_tramat_bkgFromPathPropagationBack(
    MuelmatVector& spectral_tramat_bkg,
    const ArrayOfMuelmatVector& spectral_tramat_cumulative_path) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(spectral_tramat_cumulative_path.size() == 0,
                     "Cannot extract from empty list.")
  spectral_tramat_bkg = spectral_tramat_cumulative_path.back();
}
ARTS_METHOD_ERROR_CATCH

void spectral_tramat_bkgFromPathPropagationFront(
    MuelmatVector& spectral_tramat_bkg,
    const ArrayOfMuelmatVector& spectral_tramat_cumulative_path) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(spectral_tramat_cumulative_path.size() == 0,
                     "Cannot extract from empty list.")
  spectral_tramat_bkg = spectral_tramat_cumulative_path.front();
}
ARTS_METHOD_ERROR_CATCH

void spectral_radDefaultTransmission(StokvecVector& spectral_rad,
                                     StokvecMatrix& spectral_rad_bkg,
                                     const AscendingGrid& freq_grid,
                                     const JacobianTargets& jac_targets) {
  ARTS_TIME_REPORT

  const Index nf = freq_grid.size();
  const Index nq = jac_targets.x_size();

  spectral_rad_bkg.resize(nq, nf);
  spectral_rad_bkg = 0.0;

  spectral_rad.resize(nf);
  spectral_rad = 1;
}

void single_rad_backgroundAgendasAtEndOfPath(
    const Workspace& ws,
    Stokvec& single_rad_background,
    StokvecVector& single_rad_background_jacobian,
    const Numeric& frequency,
    const JacobianTargets& jac_targets,
    const PropagationPathPoint& ray_path_point,
    const SurfaceField& surf_field,
    const SubsurfaceField& subsurf_field,
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
                                     jac_targets,
                                     ray_path_point,
                                     single_rad_space_agenda);
      break;
    case surface:
      single_rad_surface_agendaExecute(ws,
                                       single_rad_background,
                                       single_rad_background_jacobian,
                                       frequency,
                                       jac_targets,
                                       ray_path_point,
                                       surf_field,
                                       subsurf_field,
                                       single_rad_surface_agenda);
      break;
  }
}
ARTS_METHOD_ERROR_CATCH
