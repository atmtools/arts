#include <rtepack.h>
#include <workspace.h>

#include "debug.h"

void radBackground(const Workspace& ws,
                   StokvecVector& rad,
                   StokvecMatrix& drad,
                   const Vector& f_grid,
                   const JacobianTargets& jacobian_targets,
                   const Ppath& ppath,
                   const Agenda& space_radiation_agenda,
                   const Agenda& surface_radiation_agenda,
                   const Agenda& atm_radiation_agenda) {
  using enum Options::PpathBackground;
  switch (ppath.background) {
    case Undefined:
      ARTS_USER_ERROR("Undefined background type in ppath.background");
      break;
    case Space:
      space_radiation_agendaExecute(ws,
                                    rad,
                                    drad,
                                    f_grid,
                                    jacobian_targets,
                                    ppath.end_pos,
                                    ppath.end_los,
                                    space_radiation_agenda);
      break;
    case Surface:
      surface_radiation_agendaExecute(ws,
                                      rad,
                                      drad,
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
      atm_radiation_agendaExecute(ws,
                                  rad,
                                  drad,
                                  f_grid,
                                  jacobian_targets,
                                  ppath.end_pos,
                                  ppath.end_los,
                                  atm_radiation_agenda);
      break;
    case FINAL:
      ARTS_USER_ERROR("Unkown background type in ppath.background");
  }

  ARTS_USER_ERROR_IF(
      rad.nelem() not_eq f_grid.nelem(),
      "Bad size of rad, it's dimension should match the size of f_grid");

  ARTS_USER_ERROR_IF(
      static_cast<Size>(drad.nrows()) not_eq jacobian_targets.size(),
      "Bad size of drad, it's inner dimension should match the size of jacobian_targets");

  ARTS_USER_ERROR_IF(
      drad.ncols() not_eq f_grid.nelem(),
      "Bad size of drad, it's outer dimension should match the size of f_grid");
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
