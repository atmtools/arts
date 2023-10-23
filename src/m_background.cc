#include <physics_funcs.h>
#include <rtepack.h>
#include <surf.h>
#include <workspace.h>

void background_radFromPropagation(
    const Workspace& ws,
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
      static_cast<Size>(background_drad.nrows()) not_eq jacobian_targets.size(),
      "Bad size of background_drad, it's inner dimension should match the size of jacobian_targets");

  ARTS_USER_ERROR_IF(
      background_drad.ncols() not_eq f_grid.nelem(),
      "Bad size of background_drad, it's outer dimension should match the size of f_grid");
}

void dradAdaptBackground(StokvecMatrix& drad,
                         MuelmatVector& background_transmittance) {
  ARTS_USER_ERROR_IF(
      background_transmittance.nelem() not_eq drad.ncols(),
      "Bad size of cumtramat and drad, their frequency domain disagree");

  for (Index i = 0; i < drad.nrows(); ++i) {
    for (Index j = 0; j < drad.ncols(); ++j) {
      drad(i, j) = background_transmittance[j] * drad(i, j);
    }
  }
}

StokvecVector from_temp(const ExhaustiveConstVectorView& f_grid,
                        const Numeric t) {
  StokvecVector v(f_grid.nelem(), 0.0);
  std::transform(
      f_grid.begin(), f_grid.end(), v.begin(), [t](auto f) -> Stokvec {
        return {planck(f, t), 0, 0, 0};
      });
  return v;
}

void background_dradEmpty(StokvecMatrix& background_drad,
                          const Vector& f_grid,
                          const JacobianTargets& jacobian_targets) {
  dradEmpty(background_drad, f_grid, jacobian_targets);
}

void background_radCosmicBackground(StokvecVector &background_rad,
                                    const Vector &f_grid) {
  constexpr auto t = Constant::cosmic_microwave_background_temperature;
  background_rad = from_temp(f_grid, t);
}

void background_radSurfaceFieldEmission(StokvecVector &background_rad,
                                        const Vector &f_grid,
                                        const SurfaceField &surface_field,
                                        const Vector &rtp_pos) {
  ARTS_USER_ERROR_IF(rtp_pos.size() != 3, "Bad size of rtp_pos, must be 3-long")
  ARTS_USER_ERROR_IF(not surface_field.contains(Surf::Key::t),
                     "Surface field does not contain temperature")

  const auto t =
      surface_field.single_value(Surf::Key::t, rtp_pos[1], rtp_pos[2]);
  background_rad = from_temp(f_grid, t);

  //! FIXME: Surface frequency and surface temperature should be included in derivatives
}
