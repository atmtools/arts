/*===========================================================================
  ===  File description
  ===========================================================================*/

#include <workspace.h>

#include "arts_omp.h"
#include "atm.h"
#include "configtypes.h"
#include "debug.h"
#include "jacobian.h"
#include "path_point.h"
#include "physics_funcs.h"
#include "rtepack.h"
#include "sorted_grid.h"
#include "sun.h"
#include "sun_methods.h"
#include "workspace_class.h"

/*!
  \file   m_sun.cc
  \author Jon Petersen  <jon.petersen@studium.uni-hamburg.de>
          Manfred Brath  <manfred.brath@.uni-hamburg.de>
  \date   2021-02-08

  \brief  Workspace functions related to simulation of radiation fluxes.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

using Constant::pi;
using Math::pow2;

/*===========================================================================
  === The functions
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void sunBlackbody(Sun& sun,
                  // Inputs:
                  const AscendingGrid& frequency_grid,
                  const Numeric& radius,
                  const Numeric& distance,
                  const Numeric& temperature,
                  const Numeric& latitude,
                  const Numeric& longitude) {
  // some sanity checks
  ARTS_USER_ERROR_IF(distance < radius,
                     "The distance to the center of the sun (",
                     distance,
                     " m) \n"
                     " is smaller than the radius of the sun (",
                     radius,
                     " m )")

  // spectrum
  sun.spectrum = Matrix(frequency_grid.nelem(), 4, 0.);

  planck(sun.spectrum(joker, 0), frequency_grid, temperature);
  sun.spectrum *= pi;  // outgoing flux at the surface of the sun.

  sun.description = "Blackbody sun";
  sun.radius = radius;
  sun.distance = distance;
  sun.latitude = latitude;
  sun.longitude = longitude;
}

void sunsAddSun(ArrayOfSun& suns, const Sun& sun) { suns.push_back(sun); }

void sun_pathFromObserverAgenda(const Workspace& ws,
                                ArrayOfPropagationPathPoint& sun_path,
                                const SurfaceField& surface_field,
                                const Agenda& ray_path_observer_agenda,
                                const Sun& sun,
                                const Vector3& observer_pos,
                                const Numeric& angle_cut,
                                const Index& refinements,
                                const Index& just_hit) {
  find_sun_path(ws,
                sun_path,
                sun,
                ray_path_observer_agenda,
                surface_field,
                observer_pos,
                angle_cut,
                refinements,
                just_hit);
}

void ray_path_sun_pathFromPathObserver(
    const Workspace& ws,
    ArrayOfArrayOfPropagationPathPoint& ray_path_sun_path,
    const SurfaceField& surface_field,
    const Agenda& ray_path_observer_agenda,
    const ArrayOfPropagationPathPoint& ray_path,
    const Sun& sun,
    const Numeric& angle_cut,
    const Index& refinements,
    const Index& just_hit) {
  ARTS_USER_ERROR_IF(angle_cut < 0.0, "angle_cut must be positive")

  const Size np = ray_path.size();

  ray_path_sun_path.resize(np);
  if (arts_omp_in_parallel() and
      static_cast<Index>(np) >= arts_omp_get_max_threads()) {
    for (Size i = 0; i < np; ++i) {
      sun_pathFromObserverAgenda(ws,
                                 ray_path_sun_path[i],
                                 surface_field,
                                 ray_path_observer_agenda,
                                 sun,
                                 ray_path[i].pos,
                                 angle_cut,
                                 refinements,
                                 just_hit);
    }
  } else {
    String error{};

#pragma omp parallel for
    for (Size i = 0; i < np; ++i) {
      try {
        sun_pathFromObserverAgenda(ws,
                                   ray_path_sun_path[i],
                                   surface_field,
                                   ray_path_observer_agenda,
                                   sun,
                                   ray_path[i].pos,
                                   angle_cut,
                                   refinements,
                                   just_hit);
      } catch (const std::exception& e) {
#pragma omp critical
        error += e.what();
      }
    }

    ARTS_USER_ERROR_IF(error.size(), error)
  }
}

void ray_path_suns_pathFromPathObserver(
    const Workspace& ws,
    ArrayOfArrayOfArrayOfPropagationPathPoint& ray_path_suns_path,
    const SurfaceField& surface_field,
    const Agenda& ray_path_observer_agenda,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfSun& suns,
    const Numeric& angle_cut,
    const Index& refinements,
    const Index& just_hit) {
  ARTS_USER_ERROR_IF(angle_cut < 0.0, "angle_cut must be positive")

  const Size np = ray_path.size();
  const Size nsuns = suns.size();

  ray_path_suns_path.resize(np);
  for (auto& p : ray_path_suns_path) p.resize(nsuns);

  if (arts_omp_in_parallel() and
      static_cast<Index>(np * nsuns) >= arts_omp_get_max_threads()) {
    for (Size i = 0; i < np; ++i) {
      for (Size j = 0; j < nsuns; ++j) {
        sun_pathFromObserverAgenda(ws,
                                   ray_path_suns_path[i][j],
                                   surface_field,
                                   ray_path_observer_agenda,
                                   suns[j],
                                   ray_path[i].pos,
                                   angle_cut,
                                   refinements,
                                   just_hit);
      }
    }
  } else {
    String error{};

#pragma omp parallel for collapse(2)
    for (Size i = 0; i < np; ++i) {
      for (Size j = 0; j < nsuns; ++j) {
        try {
          sun_pathFromObserverAgenda(ws,
                                     ray_path_suns_path[i][j],
                                     surface_field,
                                     ray_path_observer_agenda,
                                     suns[j],
                                     ray_path[i].pos,
                                     angle_cut,
                                     refinements,
                                     just_hit);
        } catch (const std::exception& e) {
#pragma omp critical
          error += e.what();
        }
      }
    }

    ARTS_USER_ERROR_IF(error.size(), error)
  }
}

void gas_scattering_coefAirSimple(Vector& gas_scattering_coef,
                                  const AscendingGrid& frequency_grid,
                                  const AtmPoint& atmospheric_point) {
  static constexpr std::array coefficients{
      3.9729066, 4.6547659e-2, 4.5055995e-4, 2.3229848e-5};

  gas_scattering_coef.resize(frequency_grid.nelem());
  for (Index f = 0; f < frequency_grid.nelem(); f++) {
    const Numeric wavelen = Conversion::freq2wavelen(frequency_grid[f]) * 1e6;
    Numeric sum = 0;
    Numeric pows = 1;
    for (auto& coef : coefficients) {
      sum += coef * pows;
      pows /= Math::pow2(wavelen);
    }
    gas_scattering_coef[f] = 1e-32 * sum / Math::pow4(wavelen);
  }

  gas_scattering_coef *=
      number_density(atmospheric_point.pressure, atmospheric_point.temperature);
}

void ray_path_spectral_radiance_sourceAddSunsFirstOrderRayleighScattering(
    const Workspace& ws,
    ArrayOfStokvecVector& ray_path_spectral_radiance_source,
    const ArrayOfPropmatVector& ray_path_propagation_matrix,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const ArrayOfArrayOfArrayOfPropagationPathPoint& ray_path_suns_path,
    const ArrayOfSun& suns,
    const JacobianTargets& jacobian_targets,
    const AscendingGrid& frequency_grid,
    const AtmField& atmospheric_field,
    const SurfaceField& surface_field,
    const Agenda& propagation_matrix_agenda,
    const Numeric& depolarization_factor,
    const Numeric& rte_alonglos_v,
    const Index& hse_derivative) try {
  ARTS_USER_ERROR_IF(jacobian_targets.x_size(),
                     "Cannot have any Jacobian targets")

  const Size np = ray_path.size();
  ARTS_USER_ERROR_IF(
      np != ray_path_spectral_radiance_source.size(),
      "Bad ray_path_spectral_radiance_source: incorrect number of path points")
  ARTS_USER_ERROR_IF(
      np != ray_path_propagation_matrix.size(),
      "Bad ray_path_propagation_matrix: incorrect number of path points")
  ARTS_USER_ERROR_IF(np != ray_path_suns_path.size(),
                     "Bad ray_path_suns_path: incorrect number of path points")
  ARTS_USER_ERROR_IF(ray_path_atmospheric_point.size() != np,
                     "Bad ray_path_atmospheric_point: incorrect number of path points")

  const Size nsuns = suns.size();
  ARTS_USER_ERROR_IF(std::ranges::any_of(ray_path_suns_path,
                                         [nsuns](auto& ray_path_sun_path) {
                                           return ray_path_sun_path.size() !=
                                                  nsuns;
                                         }),
                     "Bad ray_path_suns_path: incorrect number of suns")

  const Index nv = frequency_grid.size();
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(ray_path_spectral_radiance_source,
                          [nv](auto& spectral_radiance_source) {
                            return nv != spectral_radiance_source.size();
                          }),
      "Bad ray_path_spectral_radiance_source: incorrect number of frequencies")
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(ray_path_propagation_matrix,
                          [nv](auto& propagation_matrix) {
                            return nv != propagation_matrix.size();
                          }),
      "Bad ray_path_propagation_matrix: incorrect number of frequencies")

  StokvecVector spectral_radiance{};
  StokvecMatrix spectral_radiance_jacobian{};
  StokvecVector spectral_radiance_background{};
  Vector gas_scattering_coef;
  const StokvecMatrix spectral_radiance_background_jacobian(0, nv);

  for (Size ip = 0; ip < np; ip++) {
    auto& spectral_radiance_source = ray_path_spectral_radiance_source[ip];
    const auto& propagation_matrix = ray_path_propagation_matrix[ip];
    const auto& ray_path_point = ray_path[ip];
    const auto& atmospheric_point = ray_path_atmospheric_point[ip];
    gas_scattering_coefAirSimple(
        gas_scattering_coef, frequency_grid, atmospheric_point);

    for (Size isun = 0; isun < nsuns; isun++) {
      const auto& sun_path = ray_path_suns_path[ip][isun];
      const auto& sun = suns[isun];

      spectral_radianceSunOrCosmicBackground(spectral_radiance_background,
                                             frequency_grid,
                                             sun_path,
                                             sun,
                                             surface_field);

      spectral_radianceClearskyBackgroundTransmission(
          ws,
          spectral_radiance,
          spectral_radiance_jacobian,
          atmospheric_field,
          frequency_grid,
          jacobian_targets,
          propagation_matrix_agenda,
          sun_path,
          spectral_radiance_background,
          spectral_radiance_background_jacobian,
          surface_field,
          rte_alonglos_v,
          hse_derivative);

      // irradiance ratio
      const Numeric radiance_2_irradiance =
          pi * suns[isun].sin_alpha_squared(sun_path.back().pos,
                                            surface_field.ellipsoid);

      ARTS_USER_ERROR_IF(spectral_radiance.size() != nv,
                         "Bad size spectral_radiance (",
                         spectral_radiance.size(),
                         ").  It should have the same size as frequency_grid (",
                         nv,
                         ")")

      const Muelmat scatmat = rtepack::rayleigh_scattering(
          sun_path.front().los, ray_path_point.los, depolarization_factor);

      // Scatter the source into the target direction
      // Direction the radiation is arriving: sun_path.front().los [MOVING TOWARDS POINT]
      // Direction the radiation is leaving: ray_path[ip].los [MOVING AWAY FROM POINT]

      // Add the source to the target
      for (Index iv = 0; iv < nv; iv++) {
        spectral_radiance_source += (radiance_2_irradiance / (4 * pi)) *
                                    inv(propagation_matrix[iv]) * scatmat *
                                    spectral_radiance[iv];
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH
