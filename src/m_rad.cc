#include <atm.h>
#include <jacobian.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <surf.h>
#include <workspace.h>

#include <algorithm>
#include <exception>

#include "arts_omp.h"
#include "auto_wsa.h"
#include "auto_wsm.h"
#include "debug.h"
#include "fwd.h"
#include "rtepack_stokes_vector.h"
#include "sorted_grid.h"
#include "workspace_agenda_class.h"

void spectral_radiance_jacobianEmpty(
    StokvecMatrix &spectral_radiance_jacobian,
    const AscendingGrid &frequency_grid,
    const JacobianTargets &jacobian_targets) try {
  spectral_radiance_jacobian.resize(jacobian_targets.x_size(),
                                    frequency_grid.size());
  spectral_radiance_jacobian = Stokvec{0.0, 0.0, 0.0, 0.0};
}
ARTS_METHOD_ERROR_CATCH

void spectral_radiance_jacobianFromBackground(
    StokvecMatrix &spectral_radiance_jacobian,
    const StokvecMatrix &spectral_radiance_background_jacobian,
    const MuelmatVector &background_transmittance) try {
  ARTS_USER_ERROR_IF(
      spectral_radiance_background_jacobian.ncols() !=
          background_transmittance.nelem(),
      "spectral_radiance_background_jacobian must have same number of rows as the "
      "size of jacobian_targets")

  //! The radiance derivative shape is the background shape
  spectral_radiance_jacobian.resize(
      spectral_radiance_background_jacobian.shape());

  //! Set the background radiance derivative as that which is seen after "this" swath
  for (Index i = 0; i < spectral_radiance_jacobian.nrows(); i++) {
    std::transform(background_transmittance.begin(),
                   background_transmittance.end(),
                   spectral_radiance_background_jacobian[i].begin(),
                   spectral_radiance_jacobian[i].begin(),
                   std::multiplies<>());
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_radiance_jacobianAddPathPropagation(
    StokvecMatrix &spectral_radiance_jacobian,
    const ArrayOfStokvecMatrix &ray_path_spectral_radiance_jacobian,
    const JacobianTargets &jacobian_targets,
    const AtmField &atmospheric_field,
    const ArrayOfPropagationPathPoint &ray_path) try {
  const auto np = ray_path_spectral_radiance_jacobian.size();
  const auto nj = spectral_radiance_jacobian.nrows();
  const auto nf = spectral_radiance_jacobian.ncols();
  const auto nt = jacobian_targets.target_count();

  if (nt == 0) return;

  ARTS_USER_ERROR_IF(
      static_cast<Size>(spectral_radiance_jacobian.nrows()) !=
          jacobian_targets.x_size(),
      "Bad size of spectral_radiance_jacobian, it's inner dimension should match the size of jacobian_targets. Sizes: ",
      spectral_radiance_jacobian.nrows(),
      " != ",
      jacobian_targets.x_size())

  ARTS_USER_ERROR_IF(
      ray_path.size() != np,
      "ray_path must have same size as the size of ray_path_spectral_radiance_jacobian.  Sizes: ",
      ray_path.size(),
      " != ",
      np)

  for (auto &dr : ray_path_spectral_radiance_jacobian) {
    ARTS_USER_ERROR_IF(
        dr.ncols() != nf or dr.nrows() != static_cast<Index>(nt),
        "ray_path_spectral_radiance_jacobian elements must have same number of rows as the size of "
        "jacobian_targets.  Sizes: ",
        dr.shape(),
        " != ",
        nt,
        ", ",
        nf)
  }

  //! Checks that the jacobian_targets can be used and throws if not
  jacobian_targets.throwing_check(nj);

  //! The derivative part from the atmosphere
  for (auto &atm_block : jacobian_targets.atm()) {
    ARTS_USER_ERROR_IF(not atmospheric_field.contains(atm_block.type),
                       "No ",
                       atm_block.type,
                       " in atmospheric_field but in jacobian_targets")
    const auto &data = atmospheric_field[atm_block.type];
    for (Size ip = 0; ip < np; ip++) {
      const auto weights = data.flat_weight(ray_path[ip].pos);
      for (auto &w : weights) {
        if (w.second != 0.0) {
          const auto i = w.first + atm_block.x_start;
          ARTS_ASSERT(i < static_cast<Size>(nj))
          std::transform(
              ray_path_spectral_radiance_jacobian[ip]
                                                         [atm_block.target_pos]
                                                             .begin(),
              ray_path_spectral_radiance_jacobian[ip]
                                                         [atm_block.target_pos]
                                                             .end(),
              spectral_radiance_jacobian[i].begin(),
              spectral_radiance_jacobian[i].begin(),
              [x = w.second](auto &a, auto &b) { return x * a + b; });
        }
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceFromPathPropagation(
    StokvecVector &spectral_radiance,
    const ArrayOfStokvecVector &ray_path_spectral_radiance) try {
  ARTS_USER_ERROR_IF(ray_path_spectral_radiance.empty(),
                     "Empty ray_path_spectral_radiance")
  spectral_radiance = ray_path_spectral_radiance.front();
}
ARTS_METHOD_ERROR_CATCH

void spectral_radiance_jacobianApplyUnit(
    StokvecMatrix &spectral_radiance_jacobian,
    const StokvecVector &spectral_radiance,
    const AscendingGrid &frequency_grid,
    const PropagationPathPoint &ray_path_point,
    const String &spectral_radiance_unit) try {
  ARTS_USER_ERROR_IF(spectral_radiance.size() != frequency_grid.size(),
                     "spectral_radiance must have same size as frequency_grid")
  ARTS_USER_ERROR_IF(
      spectral_radiance_jacobian.size() != 0 and
          spectral_radiance_jacobian.ncols() != frequency_grid.size(),
      "spectral_radiance must have same size as frequency_grid")

  const auto dF = rtepack::dunit_converter(
      to<SpectralRadianceUnitType>(spectral_radiance_unit),
      ray_path_point.nreal);

  //! Must apply the unit to the spectral radiance jacobian first
  for (Index i = 0; i < spectral_radiance_jacobian.nrows(); i++) {
    for (Index j = 0; j < spectral_radiance_jacobian.ncols(); j++) {
      spectral_radiance_jacobian(i, j) = dF(spectral_radiance_jacobian(i, j),
                                            spectral_radiance[j],
                                            frequency_grid[j]);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceApplyUnit(
    StokvecVector &spectral_radiance,
    const AscendingGrid &frequency_grid,
    const PropagationPathPoint &ray_path_point,
    const String &spectral_radiance_unit) try {
  ARTS_USER_ERROR_IF(spectral_radiance.size() != frequency_grid.size(),
                     "spectral_radiance must have same size as frequency_grid")
  const auto F = rtepack::unit_converter(
      to<SpectralRadianceUnitType>(spectral_radiance_unit),
      ray_path_point.nreal);

  std::transform(spectral_radiance.begin(),
                 spectral_radiance.end(),
                 frequency_grid.begin(),
                 spectral_radiance.begin(),
                 [&F](const Stokvec &v, const Numeric f) { return F(v, f); });
}
ARTS_METHOD_ERROR_CATCH

void measurement_vectorFromSensor(
    const Workspace &ws,
    Vector &measurement_vector,
    Matrix &measurement_vector_jacobian,
    const ArrayOfSensorObsel &measurement_vector_sensor,
    const JacobianTargets &jacobian_targets,
    const AtmField &atmospheric_field,
    const SurfaceField &surface_field,
    const String &spectral_radiance_unit,
    const Agenda &spectral_radiance_observer_agenda,
    const Index &exhaustive_) try {
  //! Flag whether or not all frequency and pos-los grids are to be assumed the same
  const bool exhaustive = static_cast<bool>(exhaustive_);

  ARTS_USER_ERROR_IF(not all_ok(measurement_vector_sensor),
                     "Measurement vector sensor infromation is not OK")
  ARTS_USER_ERROR_IF(
      exhaustive and not sensor::is_exhaustive_like(measurement_vector_sensor),
      "Measurement vector sensor infromation is not exhaustive-like despite exhaustive flag")

  measurement_vector.resize(measurement_vector_sensor.size());
  measurement_vector = 0.0;

  measurement_vector_jacobian.resize(measurement_vector_sensor.size(), jacobian_targets.x_size()
                                     );
  measurement_vector_jacobian = 0.0;

  if (measurement_vector_sensor.empty()) return;

  auto frequency_grid =
      exhaustive ? measurement_vector_sensor.front().f_grid : AscendingGrid{};
  StokvecVector spectral_radiance;
  StokvecMatrix spectral_radiance_jacobian;

  for (const auto poslos :
       (exhaustive ? measurement_vector_sensor.front().poslos_grid
                   : sensor::collect_poslos(measurement_vector_sensor))) {
    if (not exhaustive) {
      sensor::collect_f_grid(frequency_grid, measurement_vector_sensor, poslos);
      spectral_radiance.resize(frequency_grid.size());
    }

    spectral_radiance_observer_agendaExecute(ws,
                                             spectral_radiance,
                                             spectral_radiance_jacobian,
                                             frequency_grid,
                                             jacobian_targets,
                                             poslos.pos,
                                             poslos.los,
                                             atmospheric_field,
                                             surface_field,
                                             spectral_radiance_unit,
                                             spectral_radiance_observer_agenda);

    if (exhaustive) {
      sensor::exhaustive_sumup(measurement_vector,
                               measurement_vector_jacobian,
                               spectral_radiance,
                               spectral_radiance_jacobian,
                               measurement_vector_sensor,
                               poslos);
    } else {
      sensor::sumup(measurement_vector,
                    measurement_vector_jacobian,
                    spectral_radiance,
                    spectral_radiance_jacobian,
                    frequency_grid,
                    measurement_vector_sensor,
                    poslos);
    }
  }
}
ARTS_METHOD_ERROR_CATCH
