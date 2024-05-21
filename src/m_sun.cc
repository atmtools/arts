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

/*===========================================================================
  === The functions
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void sunBlackbody(Sun& sun,
                  // Inputs:
                  const AscendingGrid& f_grid,
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
  sun.spectrum = Matrix(f_grid.nelem(), 4, 0.);

  planck(sun.spectrum(joker, 0), f_grid, temperature);
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

void ray_path_spectral_radiance_sourceAddBackgroundSuns(
    const Workspace& ws,
    ArrayOfStokvecVector& ray_path_spectral_radiance_source,
    const ArrayOfPropmatVector& ray_path_propagation_matrix,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfArrayOfArrayOfPropagationPathPoint& ray_path_suns_path,
    const ArrayOfSun& suns,
    const JacobianTargets& jacobian_targets,
    const AscendingGrid& frequency_grid,
    const AtmField& atmospheric_field,
    const SurfaceField& surface_field,
    const Agenda& propagation_matrix_agenda,
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
  const StokvecMatrix spectral_radiance_background_jacobian(0, nv);

  for (Size ip = 0; ip < np; ip++) {
    for (Size isun = 0; isun < nsuns; isun++) {
      const auto& sun_path = ray_path_suns_path[ip][isun];
      const auto& sun = suns[isun];
      auto& spectral_radiance_source = ray_path_spectral_radiance_source[ip];
      const auto& propagation_matrix = ray_path_propagation_matrix[ip];

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

      ARTS_USER_ERROR_IF(spectral_radiance.size() != nv,
                         "Bad size spectral_radiance (",
                         spectral_radiance.size(),
                         ").  It should have the same size as frequency_grid (",
                         nv,
                         ")")

      // Scatter the source into the target direction
      // Direction the radiation is arriving: sun_path.front().los [MOVING TOWARDS POINT]
      // Direction the radiation is leaving: ray_path[ip].los [MOVING AWAY FROM POINT]
      spectral_radiance *= 0;  // placeholder for the actual scattering

      // Add the source to the target
      for (Index iv = 0; iv < nv; iv++) {
        spectral_radiance_source +=
            inv(propagation_matrix[iv]) * spectral_radiance[iv];
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH