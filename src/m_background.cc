#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <surf.h>
#include <workspace.h>

#include <algorithm>

#include "auto_wsa.h"
#include "auto_wsm.h"
#include "debug.h"

void spectral_radiance_backgroundAgendasAtEndOfPath(
    const Workspace& ws,
    StokvecVector& spectral_radiance_background,
    StokvecMatrix& spectral_radiance_background_jacobian,
    const Vector& f_grid,
    const JacobianTargets& jacobian_targets,
    const ArrayOfPropagationPathPoint& rad_path,
    const Agenda& spectral_radiance_background_space_agenda,
    const Agenda& spectral_radiance_background_surface_agenda) {
  ARTS_USER_ERROR_IF(rad_path.size() == 0, "Empty propagation path.")

  const auto& path_point = rad_path.back();

  using enum path::PositionType;
  switch (path_point.los_type) {
    case atm:
      ARTS_USER_ERROR(
          "Undefined what to do with an atmospheric background in *rad_path*")
      break;
    case subsurface:
      ARTS_USER_ERROR(
          "Undefined what to do with a subsurface background in *rad_path*")
      break;
    case unknown:
      ARTS_USER_ERROR("Undefined background type in *rad_path*");
      break;
    case space:
      spectral_radiance_background_space_agendaExecute(
          ws,
          spectral_radiance_background,
          spectral_radiance_background_jacobian,
          f_grid,
          jacobian_targets,
          path_point,
          spectral_radiance_background_space_agenda);
      break;
    case surface:
      spectral_radiance_background_surface_agendaExecute(
          ws,
          spectral_radiance_background,
          spectral_radiance_background_jacobian,
          f_grid,
          jacobian_targets,
          path_point,
          spectral_radiance_background_surface_agenda);
      break;
    case FINAL:
      ARTS_USER_ERROR("Unkown background type in *rad_path*");
  }

  ARTS_USER_ERROR_IF(spectral_radiance_background.nelem() not_eq f_grid.nelem(),
                     "Bad size spectral_radiance_background (",
                     spectral_radiance_background.nelem(),
                     ").  It should have the same size as f_grid (",
                     f_grid.nelem(),
                     ")");

  ARTS_USER_ERROR_IF(
      static_cast<Size>(spectral_radiance_background_jacobian.nrows()) not_eq
              jacobian_targets.x_size() or
          spectral_radiance_background_jacobian.ncols() not_eq f_grid.nelem(),
      "Bad size of spectral_radiance_background_jacobian (",
      spectral_radiance_background_jacobian.nrows(),
      'x',
      spectral_radiance_background_jacobian.ncols(),
      ").  It should have the same size as jacobian_targets (",
      jacobian_targets.x_size(),
      ") and f_grid (",
      f_grid.nelem(),
      ")");
}

namespace detail {
StokvecVector from_temp(const ExhaustiveConstVectorView& f_grid,
                        const Numeric t) {
  StokvecVector v(f_grid.nelem(), 0.0);
  std::transform(
      f_grid.begin(), f_grid.end(), v.begin(), [t](auto f) -> Stokvec {
        return {planck(f, t), 0, 0, 0};
      });
  return v;
}
}  // namespace detail

void spectral_radiance_background_jacobianEmpty(
    StokvecMatrix& spectral_radiance_background_jacobian,
    const Vector& f_grid,
    const JacobianTargets& jacobian_targets) {
  spectral_radiance_jacobianEmpty(spectral_radiance_background_jacobian, f_grid, jacobian_targets);
}

void spectral_radiance_backgroundUniformCosmicBackground(
    StokvecVector& spectral_radiance_background, const Vector& f_grid) {
  constexpr auto t = Constant::cosmic_microwave_background_temperature;
  spectral_radiance_background = detail::from_temp(f_grid, t);
}

void spectral_radiance_backgroundSurfaceBlackbody(
    StokvecVector& spectral_radiance_background,
    StokvecMatrix& spectral_radiance_background_jacobian,
    const Vector& f_grid,
    const SurfaceField& surface_field,
    const JacobianTargets& jacobian_targets,
    const PropagationPathPoint& path_point) {
  constexpr auto key = Surf::Key::t;

  ARTS_USER_ERROR_IF(not surface_field.contains(key),
                     "Surface field does not contain temperature")

  const auto t = surface_field.single_value(key, path_point.pos[1], path_point.pos[2]);
  spectral_radiance_background = detail::from_temp(f_grid, t);

  spectral_radiance_background_jacobianEmpty(
      spectral_radiance_background_jacobian, f_grid, jacobian_targets);

  if (auto pair =
          jacobian_targets.find<Jacobian::SurfaceTarget>(SurfaceKeyVal{key});
      pair.first) {
    const auto x_start = pair.second->x_start;
    const auto& surf_data = surface_field[key];
    const auto weights = surf_data.flat_weights(path_point.pos[1], path_point.pos[2]);
    for (auto& w : weights) {
      if (w.second != 0.0) {
        const auto i = w.first + x_start;
        ARTS_ASSERT(i < static_cast<Size>(
                            spectral_radiance_background_jacobian.nrows()))

        std::transform(f_grid.begin(),
                       f_grid.end(),
                       spectral_radiance_background_jacobian[i].begin(),
                       [t, x = w.second](auto f) -> Stokvec {
                         return {x * dplanck_dt(f, t), 0, 0, 0};
                       });
      }
    }
  }
}

void background_transmittanceFromPathPropagationBack(
    MuelmatVector& background_transmittance,
    const ArrayOfMuelmatVector& ppvar_cumtramat) {
  ARTS_USER_ERROR_IF(ppvar_cumtramat.size() == 0,
                     "Cannot extract from empty list.")
  background_transmittance = ppvar_cumtramat.back();
}

void background_transmittanceFromPathPropagationFront(
    MuelmatVector& background_transmittance,
    const ArrayOfMuelmatVector& ppvar_cumtramat) {
  ARTS_USER_ERROR_IF(ppvar_cumtramat.size() == 0,
                     "Cannot extract from empty list.")
  background_transmittance = ppvar_cumtramat.front();
}
