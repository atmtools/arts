#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <sun.h>
#include <surf.h>
#include <workspace.h>

#include <algorithm>

#include "debug.h"

void spectral_radiance_backgroundAgendasAtEndOfPath(
    const Workspace& ws,
    StokvecVector& spectral_radiance_background,
    StokvecMatrix& spectral_radiance_background_jacobian,
    const AscendingGrid& frequency_grid,
    const JacobianTargets& jacobian_targets,
    const PropagationPathPoint& ray_path_point,
    const SurfaceField& surface_field,
    const Agenda& spectral_radiance_space_agenda,
    const Agenda& spectral_radiance_surface_agenda) try {
  using enum PathPositionType;
  switch (ray_path_point.los_type) {
    case atm:
      ARTS_USER_ERROR("Undefined what to do with an atmospheric background")
      break;
    case subsurface:
      ARTS_USER_ERROR("Undefined what to do with a subsurface background")
      break;
    case unknown:
      ARTS_USER_ERROR("Undefined background type");
      break;
    case space:
      spectral_radiance_space_agendaExecute(
          ws,
          spectral_radiance_background,
          spectral_radiance_background_jacobian,
          frequency_grid,
          jacobian_targets,
          ray_path_point,
          spectral_radiance_space_agenda);
      break;
    case surface:
      spectral_radiance_surface_agendaExecute(
          ws,
          spectral_radiance_background,
          spectral_radiance_background_jacobian,
          frequency_grid,
          jacobian_targets,
          ray_path_point,
          surface_field,
          spectral_radiance_surface_agenda);
      break;
  }

  ARTS_USER_ERROR_IF(
      spectral_radiance_background.nelem() not_eq frequency_grid.nelem(),
      "Bad size spectral_radiance_background ({}"
      ").  It should have the same size as frequency_grid ({})",
      spectral_radiance_background.nelem(),
      frequency_grid.nelem());

  ARTS_USER_ERROR_IF(
      static_cast<Size>(spectral_radiance_background_jacobian.nrows()) not_eq
              jacobian_targets.x_size() or
          spectral_radiance_background_jacobian.ncols() not_eq
              frequency_grid.nelem(),
      "Bad size of spectral_radiance_background_jacobian ({}x{}"
      ").  It should have the same size as jacobian_targets ({}"
      ") and frequency_grid ({})",
      spectral_radiance_background_jacobian.nrows(),
      spectral_radiance_background_jacobian.ncols(),
      jacobian_targets.x_size(),
      frequency_grid.nelem());
}
ARTS_METHOD_ERROR_CATCH

namespace detail {
StokvecVector from_temp(const ExhaustiveConstVectorView& frequency_grid,
                        const Numeric t) {
  StokvecVector v(frequency_grid.nelem(), 0.0);
  std::transform(frequency_grid.begin(),
                 frequency_grid.end(),
                 v.begin(),
                 [t](auto f) -> Stokvec { return {planck(f, t), 0, 0, 0}; });
  return v;
}
}  // namespace detail

void spectral_radianceUniformCosmicBackground(
    StokvecVector& spectral_radiance, const AscendingGrid& frequency_grid) {
  constexpr auto t = Constant::cosmic_microwave_background_temperature;
  spectral_radiance = detail::from_temp(frequency_grid, t);
}

void spectral_radianceSunOrCosmicBackground(
    StokvecVector& spectral_radiance,
    const AscendingGrid& frequency_grid,
    const ArrayOfPropagationPathPoint& sun_path,
    const Sun& sun,
    const SurfaceField& surface_field) try {
  spectral_radiance.resize(frequency_grid.size());

  if (set_spectral_radiance_if_sun_intersection(
          spectral_radiance, sun, sun_path.back(), surface_field))
    return;

  spectral_radianceUniformCosmicBackground(spectral_radiance, frequency_grid);
} ARTS_METHOD_ERROR_CATCH

void spectral_radianceSunsOrCosmicBackground(
    StokvecVector& spectral_radiance,
    const AscendingGrid& frequency_grid,
    const PropagationPathPoint& ray_path_point,
    const ArrayOfSun& suns,
    const SurfaceField& surface_field) {
  spectral_radiance.resize(frequency_grid.size());

  for (auto& sun : suns) {
    if (set_spectral_radiance_if_sun_intersection(
            spectral_radiance, sun, ray_path_point, surface_field)) {
      return;
    }
  }

  spectral_radianceUniformCosmicBackground(spectral_radiance, frequency_grid);
}

void spectral_radianceSurfaceBlackbody(
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_jacobian,
    const AscendingGrid& frequency_grid,
    const SurfaceField& surface_field,
    const JacobianTargets& jacobian_targets,
    const PropagationPathPoint& ray_path_point) {
  constexpr auto key = SurfaceKey::t;

  ARTS_USER_ERROR_IF(not surface_field.contains(key),
                     "Surface field does not contain temperature")

  const auto t = surface_field.single_value(
      key, ray_path_point.pos[1], ray_path_point.pos[2]);
  spectral_radiance = detail::from_temp(frequency_grid, t);

  spectral_radiance_jacobianEmpty(
      spectral_radiance_jacobian, frequency_grid, jacobian_targets);

  if (auto pair =
          jacobian_targets.find<Jacobian::SurfaceTarget>(SurfaceKeyVal{key});
      pair.first) {
    const auto x_start = pair.second->x_start;
    const auto& surf_data = surface_field[key];
    const auto weights =
        surf_data.flat_weights(ray_path_point.pos[1], ray_path_point.pos[2]);
    for (auto& w : weights) {
      if (w.second != 0.0) {
        const auto i = w.first + x_start;
        ARTS_ASSERT(i < static_cast<Size>(spectral_radiance_jacobian.nrows()))

        std::transform(frequency_grid.begin(),
                       frequency_grid.end(),
                       spectral_radiance_jacobian[i].begin(),
                       [t, x = w.second](auto f) -> Stokvec {
                         return {x * dplanck_dt(f, t), 0, 0, 0};
                       });
      }
    }
  }
}

void transmission_matrix_backgroundFromPathPropagationBack(
    MuelmatVector& transmission_matrix_background,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix_cumulative) try {
  ARTS_USER_ERROR_IF(ray_path_transmission_matrix_cumulative.size() == 0,
                     "Cannot extract from empty list.")
  transmission_matrix_background = ray_path_transmission_matrix_cumulative.back();
}
ARTS_METHOD_ERROR_CATCH

void transmission_matrix_backgroundFromPathPropagationFront(
    MuelmatVector& transmission_matrix_background,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix_cumulative) try {
  ARTS_USER_ERROR_IF(ray_path_transmission_matrix_cumulative.size() == 0,
                     "Cannot extract from empty list.")
  transmission_matrix_background = ray_path_transmission_matrix_cumulative.front();
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceDefaultTransmission(
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_background,
    const AscendingGrid& frequency_grid,
    const JacobianTargets& jacobian_targets) {
  const Index nf = frequency_grid.nelem();
  const Index nq = jacobian_targets.x_size();

  spectral_radiance_background.resize(nq, nf);
  spectral_radiance_background = 0.0;

  spectral_radiance.resize(nf);
  spectral_radiance = 1;
}
