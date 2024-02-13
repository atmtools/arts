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
#include "sorted_grid.h"
#include "workspace_agenda_class.h"

void spectral_radiance_jacobianEmpty(StokvecMatrix &spectral_radiance_jacobian,
                                     const AscendingGrid &frequency_grid,
                                     const JacobianTargets &jacobian_targets) {
  spectral_radiance_jacobian.resize(jacobian_targets.x_size(),
                                    frequency_grid.size());
  spectral_radiance_jacobian = Stokvec{0.0, 0.0, 0.0, 0.0};
}

void spectral_radiance_jacobianFromBackground(
    StokvecMatrix &spectral_radiance_jacobian,
    const StokvecMatrix &spectral_radiance_background_jacobian,
    const MuelmatVector &background_transmittance) {
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

void spectral_radiance_jacobianAddPathPropagation(
    StokvecMatrix &spectral_radiance_jacobian,
    const ArrayOfStokvecMatrix &propagation_path_spectral_radiance_jacobian,
    const JacobianTargets &jacobian_targets,
    const AtmField &atmospheric_field,
    const ArrayOfPropagationPathPoint &propagation_path) {
  const auto np = propagation_path_spectral_radiance_jacobian.size();
  const auto nj = spectral_radiance_jacobian.nrows();
  const auto nf = spectral_radiance_jacobian.ncols();

  ARTS_USER_ERROR_IF(
      static_cast<Size>(spectral_radiance_jacobian.nrows()) !=
          jacobian_targets.x_size(),
      "Bad size of spectral_radiance_jacobian, it's inner dimension should match the size of jacobian_targets")

  ARTS_USER_ERROR_IF(
      propagation_path.size() != np,
      "propagation_path must have same size as the size of propagation_path_spectral_radiance_jacobian")

  for (auto &dr : propagation_path_spectral_radiance_jacobian) {
    ARTS_USER_ERROR_IF(
        dr.nrows() != nj,
        "propagation_path_spectral_radiance_jacobian elements must have same number of rows as the size of "
        "jacobian_targets")
    ARTS_USER_ERROR_IF(
        dr.ncols() != nf,
        "propagation_path_spectral_radiance_jacobian elements must have same number of columns as the size of "
        "cumulative_transmission")
  }

  //! Checks that the jacobian_targets can be used and throws if not
  jacobian_targets.throwing_check(nj);

  //! The altitude, latitude and longitude vectors must be copied because of how atmospheric_field works
  const auto [alt, lat, lon] = [&]() {
    std::array<Vector, 3> g{Vector(np), Vector(np), Vector(np)};
    for (Size i = 0; i < np; i++) {
      g[0][i] = propagation_path[i].pos[0];
      g[1][i] = propagation_path[i].pos[1];
      g[2][i] = propagation_path[i].pos[2];
    }
    return g;
  }();

  //! The derivative part from the atmosphere
  for (auto &atm_block : jacobian_targets.atm()) {
    ARTS_USER_ERROR_IF(not atmospheric_field.contains(atm_block.type),
                       "No ",
                       atm_block.type,
                       " in atmospheric_field but in jacobian_targets")
    const auto &data = atmospheric_field[atm_block.type];
    const auto weights = data.flat_weights(alt, lat, lon);
    ARTS_ASSERT(weights.size() == np)

    for (Size j = 0; j < np; j++) {
      for (auto &w : weights[j]) {
        if (w.second != 0.0) {
          const auto i = w.first + atm_block.x_start;
          ARTS_ASSERT(i < static_cast<Size>(nj))
          std::transform(
              propagation_path_spectral_radiance_jacobian[j]
                                                         [atm_block.target_pos]
                                                             .begin(),
              propagation_path_spectral_radiance_jacobian[j]
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

void spectral_radianceFromPathPropagation(
    StokvecVector &spectral_radiance,
    const ArrayOfStokvecVector &propagation_path_spectral_radiance) {
  ARTS_USER_ERROR_IF(propagation_path_spectral_radiance.empty(),
                     "Empty propagation_path_spectral_radiance")
  spectral_radiance = propagation_path_spectral_radiance.front();
}

ENUMCLASS(
    SpectralRadianceUnitType, char, RJBT, PlanckBT, W_m2_m_sr, W_m2_m1_sr, unit)

void spectral_radianceApplyUnit(StokvecVector &spectral_radiance_with_unit,
                                const StokvecVector &spectral_radiance,
                                const AscendingGrid &frequency_grid,
                                const String &spectral_radiance_unit) try {
  const SpectralRadianceUnitType unit =
      toSpectralRadianceUnitTypeOrThrow(spectral_radiance_unit);

  ARTS_USER_ERROR_IF(spectral_radiance.size() != frequency_grid.size(),
                     "spectral_radiance must have same size as frequency_grid")

  spectral_radiance_with_unit.resize(spectral_radiance.size());

  ARTS_USER_ERROR_IF(
      unit == SpectralRadianceUnitType::unit,
      "No need to use this method with *spectral_radiance_unit* = \"unit\".");

  ARTS_USER_ERROR_IF(
      std::ranges::any_of(spectral_radiance,
                          [](const Stokvec &v) { return v.I() > 1e-3; }),
      "The spectrum matrix *spectral_radiance* is required to have original radiance\n"
      "unit, but this seems not to be the case. This as a value above\n"
      "1e-3 is found in *spectral_radiance*.")

  switch (unit) {
    case SpectralRadianceUnitType::unit:
      break;
    case SpectralRadianceUnitType::RJBT: {
      std::transform(spectral_radiance.begin(),
                     spectral_radiance.end(),
                     frequency_grid.begin(),
                     spectral_radiance_with_unit.begin(),
                     [&](const Stokvec &v, const Numeric f) {
                       const Numeric scfac = invrayjean(1.0, f);
                       return Stokvec{v.I() * scfac,
                                      2 * v.Q() * scfac,
                                      2 * v.U() * scfac,
                                      2 * v.V() * scfac};
                     });
    } break;
    case SpectralRadianceUnitType::PlanckBT: {
      std::transform(spectral_radiance.begin(),
                     spectral_radiance.end(),
                     frequency_grid.begin(),
                     spectral_radiance_with_unit.begin(),
                     [&](const Stokvec &v, const Numeric f) {
                       return Stokvec{invplanck(v.I(), f),
                                      invplanck(0.5 * (v.I() + v.Q()), f) -
                                          invplanck(0.5 * (v.I() - v.Q()), f),
                                      invplanck(0.5 * (v.I() + v.U()), f) -
                                          invplanck(0.5 * (v.I() - v.U()), f),
                                      invplanck(0.5 * (v.I() + v.V()), f) -
                                          invplanck(0.5 * (v.I() - v.V()), f)};
                     });
    } break;
    case SpectralRadianceUnitType::W_m2_m_sr: {
      std::transform(
          spectral_radiance.begin(),
          spectral_radiance.end(),
          frequency_grid.begin(),
          spectral_radiance_with_unit.begin(),
          [&](const Stokvec &v, const Numeric f) {
            const Numeric scfac = f * (f / Constant::c);
            return Stokvec{
                v.I() * scfac, v.Q() * scfac, v.U() * scfac, v.V() * scfac};
          });
    } break;
    case SpectralRadianceUnitType::W_m2_m1_sr: {
      std::transform(spectral_radiance.begin(),
                     spectral_radiance.end(),
                     spectral_radiance_with_unit.begin(),
                     [&](const Stokvec &v) { return v * Constant::c; });
    } break;
    case SpectralRadianceUnitType::FINAL: {
      ARTS_ASSERT(false)
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void sensor_radianceFromObservers(
    const Workspace &ws,
    StokvecVector &sensor_radiance,
    StokvecMatrix &sensor_radiance_jacobian,
    const AscendingGrid &frequency_grid,
    const JacobianTargets &jacobian_targets,
    const AtmField &atmospheric_field,
    const SurfaceField &surface_field,
    const Agenda &spectral_radiance_observer_agenda,
    const ArrayOfVector3 &pos,
    const ArrayOfVector2 &los) try {
  const Index nx = jacobian_targets.x_size();
  const Index nf = frequency_grid.size();
  const Size np = pos.size();

  ARTS_USER_ERROR_IF(np != los.size(),
                     "pos, "
                     "los, and "
                     "sensor_radiance_observer_weight must all have the same "
                     "size.\nThey do not, but have instead the sizes: ",
                     np,
                     ", ",
                     pos.size(),
                     ", and ",
                     los.size(),
                     ", respectively.")

  sensor_radiance.resize(nf);
  sensor_radiance_jacobian.resize(nx, nf);
  sensor_radiance = 0.;
  sensor_radiance_jacobian = 0.;

  if (np == 0) return;

  const auto summing_up = [w = 1.0 / static_cast<Numeric>(np)](
                              const Stokvec &v1, const Stokvec &v2) {
    return v1 * w + v2;
  };

  StokvecVector spectral_radiance(nf);
  StokvecMatrix spectral_radiance_jacobian(nx, nf);

  if (arts_omp_in_parallel() or arts_omp_get_max_threads() == 1) {
    for (Size ipos = 0; ipos < np; ipos++) {
      spectral_radiance_observer_agendaExecute(
          ws,
          spectral_radiance,
          spectral_radiance_jacobian,
          frequency_grid,
          jacobian_targets,
          pos[ipos],
          los[ipos],
          atmospheric_field,
          surface_field,
          spectral_radiance_observer_agenda);

      ARTS_USER_ERROR_IF(spectral_radiance.shape() != sensor_radiance.shape(),
                         "spectral_radiance has wrong shape")
      ARTS_USER_ERROR_IF(spectral_radiance_jacobian.shape() !=
                             sensor_radiance_jacobian.shape(),
                         "spectral_radiance has wrong shape")

      std::transform(spectral_radiance.begin(),
                     spectral_radiance.end(),
                     sensor_radiance.begin(),
                     sensor_radiance.begin(),
                     summing_up);
      std::transform(spectral_radiance_jacobian.elem_begin(),
                     spectral_radiance_jacobian.elem_end(),
                     sensor_radiance_jacobian.elem_begin(),
                     sensor_radiance_jacobian.elem_begin(),
                     summing_up);
    }
  } else {
    ArrayOfString fail_msg;
    bool do_abort = false;

#pragma omp parallel for firstprivate(spectral_radiance, \
                                          spectral_radiance_jacobian)
    for (Size ipos = 0; ipos < np; ipos++) {
      if (do_abort) continue;
      try {
        spectral_radiance_observer_agendaExecute(
            ws,
            spectral_radiance,
            spectral_radiance_jacobian,
            frequency_grid,
            jacobian_targets,
            pos[ipos],
            los[ipos],
            atmospheric_field,
            surface_field,
            spectral_radiance_observer_agenda);

        ARTS_USER_ERROR_IF(spectral_radiance.shape() != sensor_radiance.shape(),
                           "spectral_radiance has wrong shape")
        ARTS_USER_ERROR_IF(spectral_radiance_jacobian.shape() !=
                               sensor_radiance_jacobian.shape(),
                           "spectral_radiance has wrong shape")
#pragma omp critical
        {
          std::transform(spectral_radiance.begin(),
                         spectral_radiance.end(),
                         sensor_radiance.begin(),
                         sensor_radiance.begin(),
                         summing_up);
          std::transform(spectral_radiance_jacobian.elem_begin(),
                         spectral_radiance_jacobian.elem_end(),
                         sensor_radiance_jacobian.elem_begin(),
                         sensor_radiance_jacobian.elem_begin(),
                         summing_up);
        }
      } catch (const std::exception &e) {
#pragma omp critical
        {
          fail_msg.emplace_back(var_string(
              "Error executing agenda at position ", ipos, ":\n", e.what()));
          do_abort = true;
        }
      }
    }

    ARTS_USER_ERROR_IF(
        do_abort, "Error messages from failed cases:\n", fail_msg)
  }
}
ARTS_METHOD_ERROR_CATCH
