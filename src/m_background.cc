#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <surf.h>
#include <workspace.h>

#include <algorithm>

#include "auto_wsm.h"
#include "debug.h"

void background_radFromPath2(const Workspace& ws,
                             StokvecVector& background_rad,
                             StokvecMatrix& background_drad,
                             const Vector& f_grid,
                             const JacobianTargets& jacobian_targets,
                             const ArrayOfPropagationPathPoint& rad_path,
                             const Agenda& space_radiation_agenda,
                             const Agenda& surface_radiation_agenda) {
  ARTS_USER_ERROR_IF(rad_path.size() == 0, "Empty propagation path.")

  using enum path::PositionType;
  switch (rad_path.back().los_type) {
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
      space_radiation_agendaExecute(ws,
                                    background_rad,
                                    background_drad,
                                    f_grid,
                                    jacobian_targets,
                                    Vector{rad_path.back().pos},
                                    Vector{rad_path.back().los},
                                    space_radiation_agenda);
      break;
    case surface:
      surface_radiation_agendaExecute(ws,
                                      background_rad,
                                      background_drad,
                                      f_grid,
                                      jacobian_targets,
                                      Vector{rad_path.back().pos},
                                      Vector{rad_path.back().los},
                                      surface_radiation_agenda);
      break;
    case FINAL:
      ARTS_USER_ERROR("Unkown background type in *rad_path*");
  }

  ARTS_USER_ERROR_IF(
      background_rad.nelem() not_eq f_grid.nelem(),
      "Bad size of background_rad, it's dimension should match the size of f_grid");

  ARTS_USER_ERROR_IF(
      static_cast<Size>(background_drad.nrows()) not_eq
          jacobian_targets.x_size(),
      "Bad size of background_drad, it's inner dimension should match the size of jacobian_targets");

  ARTS_USER_ERROR_IF(
      background_drad.ncols() not_eq f_grid.nelem(),
      "Bad size of background_drad, it's outer dimension should match the size of f_grid");
}

void background_radFromPath(const Workspace& ws,
                            StokvecVector& background_rad,
                            StokvecMatrix& background_drad,
                            const Vector& f_grid,
                            const JacobianTargets& jacobian_targets,
                            const Ppath& ppath,
                            const Agenda& space_radiation_agenda,
                            const Agenda& surface_radiation_agenda,
                            const Agenda& stop_distance_radiation_agenda) {
  using enum Options::PpathBackground;
  switch (ppath.background) {
    case Undefined:
      ARTS_USER_ERROR("Undefined background type in ppath.background");
      break;
    case Space:
      space_radiation_agendaExecute(ws,
                                    background_rad,
                                    background_drad,
                                    f_grid,
                                    jacobian_targets,
                                    ppath.end_pos,
                                    ppath.end_los,
                                    space_radiation_agenda);
      break;
    case Surface:
      surface_radiation_agendaExecute(ws,
                                      background_rad,
                                      background_drad,
                                      f_grid,
                                      jacobian_targets,
                                      ppath.end_pos,
                                      ppath.end_los,
                                      surface_radiation_agenda);
      break;
    case Cloudbox:
      ARTS_USER_ERROR("Cloudbox background type not implemented");
      break;
    case Transmitter:
      ARTS_USER_ERROR("Transmitter background type not implemented");
      break;
    case StopDistance:
      stop_distance_radiation_agendaExecute(ws,
                                            background_rad,
                                            background_drad,
                                            f_grid,
                                            jacobian_targets,
                                            ppath.end_pos,
                                            ppath.end_los,
                                            stop_distance_radiation_agenda);
      break;
    case FINAL:
      ARTS_USER_ERROR("Unkown background type in ppath.background");
  }

  ARTS_USER_ERROR_IF(
      background_rad.nelem() not_eq f_grid.nelem(),
      "Bad size of background_rad, it's dimension should match the size of f_grid");

  ARTS_USER_ERROR_IF(
      static_cast<Size>(background_drad.nrows()) not_eq
          jacobian_targets.x_size(),
      "Bad size of background_drad, it's inner dimension should match the size of jacobian_targets");

  ARTS_USER_ERROR_IF(
      background_drad.ncols() not_eq f_grid.nelem(),
      "Bad size of background_drad, it's outer dimension should match the size of f_grid");
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

void background_dradEmpty(StokvecMatrix& background_drad,
                          const Vector& f_grid,
                          const JacobianTargets& jacobian_targets) {
  dradEmpty(background_drad, f_grid, jacobian_targets);
}

void background_radCosmicBackground(StokvecVector& background_rad,
                                    const Vector& f_grid) {
  constexpr auto t = Constant::cosmic_microwave_background_temperature;
  background_rad = detail::from_temp(f_grid, t);
}

void background_radSurfaceFieldEmission(StokvecVector& background_rad,
                                        StokvecMatrix& background_drad,
                                        const Vector& f_grid,
                                        const SurfaceField& surface_field,
                                        const JacobianTargets& jacobian_targets,
                                        const Vector& rtp_pos) {
  constexpr auto key = Surf::Key::t;

  ARTS_USER_ERROR_IF(rtp_pos.size() != 3, "Bad size of rtp_pos, must be 3-long")
  ARTS_USER_ERROR_IF(not surface_field.contains(key),
                     "Surface field does not contain temperature")

  const auto t = surface_field.single_value(key, rtp_pos[1], rtp_pos[2]);
  background_rad = detail::from_temp(f_grid, t);

  background_dradEmpty(background_drad, f_grid, jacobian_targets);

  if (auto pair =
          jacobian_targets.find<Jacobian::SurfaceTarget>(SurfaceKeyVal{key});
      pair.first) {
    const auto x_start = pair.second->x_start;
    const auto& surf_data = surface_field[key];
    const auto weights = surf_data.flat_weights(rtp_pos[1], rtp_pos[2]);
    for (auto& w : weights) {
      if (w.second != 0.0) {
        const auto i = w.first + x_start;
        ARTS_ASSERT(i < static_cast<Size>(background_drad.nrows()))

        std::transform(f_grid.begin(),
                       f_grid.end(),
                       background_drad[i].begin(),
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
