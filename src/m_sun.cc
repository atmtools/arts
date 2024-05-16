/*===========================================================================
  ===  File description
  ===========================================================================*/

#include <workspace.h>

#include "arts_omp.h"
#include "atm.h"
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
void sunsAddSingleBlackbody(ArrayOfSun& suns,
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

  Sun& new_sun = suns.emplace_back();

  // spectrum
  new_sun.spectrum = Matrix(f_grid.nelem(), 4, 0.);

  planck(new_sun.spectrum(joker, 0), f_grid, temperature);
  new_sun.spectrum *= pi;  // outgoing flux at the surface of the sun.

  new_sun.description = "Blackbody sun";
  new_sun.radius = radius;
  new_sun.distance = distance;
  new_sun.latitude = latitude;
  new_sun.longitude = longitude;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void sunsAddSingleFromGrid(ArrayOfSun& suns,
                           // Inputs:
                           const AscendingGrid& f_grid,
                           const GriddedField2& sun_spectrum_raw,
                           const Numeric& radius,
                           const Numeric& distance,
                           const Numeric& temperature,
                           const Numeric& latitude,
                           const Numeric& longitude,
                           const String& description) {
  // some sanity checks
  ARTS_USER_ERROR_IF(distance < radius,
                     "The distance to the center of the sun (",
                     distance,
                     " m) \n"
                     " is smaller than the radius of the sun (",
                     radius,
                     " m )")

  // init sun
  Sun& new_sun = suns.emplace_back();
  new_sun.spectrum = regrid_sun_spectrum(sun_spectrum_raw,
                                         f_grid,
                                         temperature);  // set spectrum
  new_sun.description = description;
  new_sun.radius = radius;
  new_sun.distance = distance;
  new_sun.latitude = latitude;
  new_sun.longitude = longitude;
}

void sunsOff(ArrayOfSun& suns) {
  // create empty Array of Matrix for the sun_spectrum
  suns.resize(0);
}

void sun_pathFromObserverAgenda(const Workspace& ws,
                                ArrayOfPropagationPathPoint& sun_path,
                                const SurfaceField& surface_field,
                                const Agenda& propagation_path_observer_agenda,
                                const Sun& sun,
                                const Vector3& observer_pos,
                                const Numeric& za_cut,
                                const Index& just_hit) {
  find_sun_path(ws,
                sun_path,
                sun,
                propagation_path_observer_agenda,
                surface_field,
                observer_pos,
                za_cut,
                just_hit);
}

void sun_pathsFromPathObserver(
    const Workspace& ws,
    ArrayOfArrayOfPropagationPathPoint& sun_paths,
    const SurfaceField& surface_field,
    const Agenda& propagation_path_observer_agenda,
    const ArrayOfPropagationPathPoint& propagation_path,
    const Sun& sun,
    const Numeric& za_cut,
    const Index& just_hit) {
  const Size np = propagation_path.size();

  sun_paths.resize(np);
  if (arts_omp_in_parallel()) {
    for (Size i = 0; i < np; ++i) {
      find_sun_path(ws,
                    sun_paths[i],
                    sun,
                    propagation_path_observer_agenda,
                    surface_field,
                    propagation_path[i].pos,
                    za_cut,
                    just_hit);
    }
  } else {
    String error{};

#pragma omp parallel for
    for (Size i = 0; i < np; ++i) {
      try {
        find_sun_path(ws,
                      sun_paths[i],
                      sun,
                      propagation_path_observer_agenda,
                      surface_field,
                      propagation_path[i].pos,
                      za_cut,
                      just_hit);
      } catch (const std::exception& e) {
#pragma omp critical
        error += e.what();
      }
    }

    ARTS_USER_ERROR_IF(error.size(), error)
  }
}

void propagation_path_spectral_radiance_solarClearskyTransmissionFromBackground(
    const Workspace& ws,
    ArrayOfStokvecVector& propagation_path_spectral_radiance_solar,
    ArrayOfStokvecMatrix& propagation_path_spectral_radiance_solar_jacobian,
    const ArrayOfArrayOfPropagationPathPoint& sun_paths,
    const AtmField& atmospheric_field,
    const AscendingGrid& frequency_grid,
    // const JacobianTargets &jacobian_targets,
    const Agenda& propagation_matrix_agenda,
    const SurfaceField& surface_field,
    const Numeric& rte_alonglos_v,
    const Index& hse_derivative,
    const Sun& sun) {
  const Size np = sun_paths.size();
  propagation_path_spectral_radiance_solar.resize(np);
  propagation_path_spectral_radiance_solar_jacobian.resize(np);

  const ArrayOfSun suns{sun};
  StokvecVector spectral_radiance_background(frequency_grid.nelem(), 0.);
  const StokvecMatrix spectral_radiance_background_jacobian(
      0, frequency_grid.nelem());
  const String spectral_radiance_unit{"1"};
  const JacobianTargets jacobian_targets{};

  if (arts_omp_in_parallel()) {
    for (Size i = 0; i < np; ++i) {
      spectral_radianceSunOrCosmicBackground(spectral_radiance_background,
                                             frequency_grid,
                                             sun_paths[i].back(),
                                             surface_field,
                                             suns);

      spectral_radianceClearskyBackgroundTransmission(
          ws,
          propagation_path_spectral_radiance_solar[i],
          propagation_path_spectral_radiance_solar_jacobian[i],
          atmospheric_field,
          frequency_grid,
          jacobian_targets,
          propagation_matrix_agenda,
          sun_paths[i],
          spectral_radiance_background,
          spectral_radiance_background_jacobian,
          spectral_radiance_unit,
          surface_field,
          rte_alonglos_v,
          hse_derivative);
    }
  } else {
    String error{};

#pragma omp parallel for firstprivate(spectral_radiance_background)
    for (Size i = 0; i < np; ++i) {
      try {
        spectral_radianceSunOrCosmicBackground(spectral_radiance_background,
                                               frequency_grid,
                                               sun_paths[i].back(),
                                               surface_field,
                                               suns);

        spectral_radianceClearskyBackgroundTransmission(
            ws,
            propagation_path_spectral_radiance_solar[i],
            propagation_path_spectral_radiance_solar_jacobian[i],
            atmospheric_field,
            frequency_grid,
            jacobian_targets,
            propagation_matrix_agenda,
            sun_paths[i],
            spectral_radiance_background,
            spectral_radiance_background_jacobian,
            spectral_radiance_unit,
            surface_field,
            rte_alonglos_v,
            hse_derivative);
      } catch (const std::exception& e) {
#pragma omp critical
        error += e.what();
      }
    }

    ARTS_USER_ERROR_IF(error.size(), error)
  }
}
