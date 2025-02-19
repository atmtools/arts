/*===========================================================================
  ===  File description
  ===========================================================================*/

#include <workspace.h>

#include <algorithm>

#include "arts_omp.h"
#include "atm.h"
#include "compare.h"
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
void sunFromGrid(Sun& sun,
                 // Inputs:
                 const AscendingGrid& f_grid,
                 const Numeric& latitude,
                 const Numeric& longitude,
                 const GriddedField2& sun_spectrum_raw,
                 const Numeric& radius,
                 const Numeric& distance,
                 const Numeric& temperature,
                 const String& description) {
  ARTS_TIME_REPORT

  // some sanity checks
  ARTS_USER_ERROR_IF(distance < radius,
                     "The distance to the center of the sun (",
                     distance,
                     " m) \n"
                     " is smaller than the radius of the sun (",
                     radius,
                     " m )")

  // init sun
  sun.spectrum = regrid_sun_spectrum(
      sun_spectrum_raw, f_grid, temperature);  // set spectrum
  sun.description = description;
  sun.radius      = radius;
  sun.distance    = distance;
  sun.latitude    = latitude;
  sun.longitude   = longitude;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void sunBlackbody(Sun& sun,
                  // Inputs:
                  const AscendingGrid& frequency_grid,
                  const Numeric& latitude,
                  const Numeric& longitude,
                  const Numeric& radius,
                  const Numeric& distance,
                  const Numeric& temperature) {
  ARTS_TIME_REPORT

  // some sanity checks
  ARTS_USER_ERROR_IF(distance < radius,
                     "The distance to the center of the sun ({} m) \n"
                     " is smaller than the radius of the sun ({} m)",
                     distance,
                     radius)

  // spectrum
  sun.spectrum = Matrix(frequency_grid.size(), 4, 0.);

  planck(sun.spectrum[joker, 0], frequency_grid, temperature);
  sun.spectrum *= pi;  // outgoing flux at the surface of the sun.

  sun.description = "Blackbody sun";
  sun.radius      = radius;
  sun.distance    = distance;
  sun.latitude    = latitude;
  sun.longitude   = longitude;
}

void sunsAddSun(ArrayOfSun& suns, const Sun& sun) {
  ARTS_TIME_REPORT

  suns.push_back(sun);
}

void sun_pathFromObserverAgenda(const Workspace& ws,
                                ArrayOfPropagationPathPoint& sun_path,
                                const SurfaceField& surface_field,
                                const Agenda& ray_path_observer_agenda,
                                const Sun& sun,
                                const Vector3& observer_pos,
                                const Numeric& angle_cut,
                                const Index& refinements,
                                const Index& just_hit) {
  ARTS_TIME_REPORT

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
  ARTS_TIME_REPORT

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

    ARTS_USER_ERROR_IF(error.size(), "{}", error)
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
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(angle_cut < 0.0, "angle_cut must be positive")

  const Size np    = ray_path.size();
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

    ARTS_USER_ERROR_IF(error.size(), "{}", error)
  }
}

void propagation_matrix_scatteringInit(
    PropmatVector& propagation_matrix_scattering,
    const AscendingGrid& frequency_grid) {
  ARTS_TIME_REPORT

  propagation_matrix_scattering.resize(frequency_grid.size());
  propagation_matrix_scattering = 0.0;
}

void propagation_matrix_scatteringAirSimple(
    PropmatVector& propagation_matrix_scattering,
    const AscendingGrid& frequency_grid,
    const AtmPoint& atmospheric_point) {
  const Size nf = frequency_grid.size();
  ARTS_USER_ERROR_IF(
      propagation_matrix_scattering.size() != nf,
      "Mismatch in size of propagation_matrix_scattering and frequency_grid")

  static constexpr std::array coefficients{
      3.9729066, 4.6547659e-2, 4.5055995e-4, 2.3229848e-5};

  const Numeric nd =
      number_density(atmospheric_point.pressure, atmospheric_point.temperature);
  for (Size f = 0; f < nf; f++) {
    const Numeric wavelen = Conversion::freq2wavelen(frequency_grid[f]) * 1e6;
    Numeric sum           = 0;
    Numeric pows          = 1;
    for (auto& coef : coefficients) {
      sum  += coef * pows;
      pows /= Math::pow2(wavelen);
    }
    propagation_matrix_scattering[f].A() +=
        1e-32 * nd * sum / Math::pow4(wavelen);
  }
}

void ray_path_propagation_matrix_scatteringFromPath(
    const Workspace& ws,
    ArrayOfPropmatVector& ray_path_propagation_matrix_scattering,
    const Agenda& propagation_matrix_scattering_agenda,
    const ArrayOfAscendingGrid& ray_path_frequency_grid,
    const ArrayOfAtmPoint& ray_path_atmospheric_point) {
  ARTS_TIME_REPORT

  const Size np = ray_path_frequency_grid.size();
  ARTS_USER_ERROR_IF(
      np != ray_path_atmospheric_point.size(),
      "Bad ray_path_atmospheric_point: incorrect number of path points")

  ray_path_propagation_matrix_scattering.resize(np);
  if (arts_omp_in_parallel()) {
    for (Size ip = 0; ip < np; ip++) {
      propagation_matrix_scattering_agendaExecute(
          ws,
          ray_path_propagation_matrix_scattering[ip],
          ray_path_frequency_grid[ip],
          ray_path_atmospheric_point[ip],
          propagation_matrix_scattering_agenda);
    }
  } else {
    String error{};
#pragma omp parallel for
    for (Size ip = 0; ip < np; ip++) {
      try {
        propagation_matrix_scattering_agendaExecute(
            ws,
            ray_path_propagation_matrix_scattering[ip],
            ray_path_frequency_grid[ip],
            ray_path_atmospheric_point[ip],
            propagation_matrix_scattering_agenda);
      } catch (const std::exception& e) {
#pragma omp critical
        error += e.what();
      }
    }

    ARTS_USER_ERROR_IF(error.size(), "{}", error)
  }
}

void ray_path_spectral_radiance_sourceAddScattering(
    ArrayOfStokvecVector& ray_path_spectral_radiance_source,
    const ArrayOfStokvecVector& ray_path_spectral_radiance_scattering,
    const ArrayOfPropmatVector& ray_path_propagation_matrix) {
  ARTS_TIME_REPORT

  const Size np = ray_path_propagation_matrix.size();
  ARTS_USER_ERROR_IF(
      np != ray_path_spectral_radiance_source.size(),
      "Bad ray_path_propagation_matrix_scattering: incorrect number of path points")
  ARTS_USER_ERROR_IF(
      np != ray_path_spectral_radiance_scattering.size(),
      "Bad ray_path_propagation_matrix_scattering: incorrect number of path points")

  if (np == 0) return;
  const Size nf = ray_path_propagation_matrix.front().size();
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(ray_path_spectral_radiance_source,
                          Cmp::ne(nf),
                          [](auto& v) { return v.size(); }),
      "Mismatch frequency size of ray_path_propagation_matrix and ray_path_spectral_radiance_source")
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(ray_path_spectral_radiance_scattering,
                          Cmp::ne(nf),
                          [](auto& v) { return v.size(); }),
      "Mismatch frequency size of ray_path_propagation_matrix and ray_path_spectral_radiance_scattering")

  if (arts_omp_in_parallel()) {
    for (Size ip = 0; ip < np; ip++) {
      for (Size iv = 0; iv < nf; iv++) {
        ray_path_spectral_radiance_source[ip][iv] +=
            inv(ray_path_propagation_matrix[ip][iv]) *
            ray_path_spectral_radiance_scattering[ip][iv];
      }
    }
  } else {
#pragma omp parallel for collapse(2)
    for (Size ip = 0; ip < np; ip++) {
      for (Size iv = 0; iv < nf; iv++) {
        ray_path_spectral_radiance_source[ip][iv] +=
            inv(ray_path_propagation_matrix[ip][iv]) *
            ray_path_spectral_radiance_scattering[ip][iv];
      }
    }
  }
}

void ray_path_spectral_radiance_scatteringSunsFirstOrderRayleigh(
    const Workspace& ws,
    // [np, nf]:
    ArrayOfStokvecVector& ray_path_spectral_radiance_scattering,
    // [np, nf]:
    const ArrayOfPropmatVector& ray_path_propagation_matrix_scattering,
    // [np]:
    const ArrayOfPropagationPathPoint& ray_path,
    // [np, suns, np2]:
    const ArrayOfArrayOfArrayOfPropagationPathPoint& ray_path_suns_path,
    // [nsuns]:
    const ArrayOfSun& suns,
    // [njac]:
    const JacobianTargets& jacobian_targets,
    // [nf]:
    const AscendingGrid& frequency_grid,
    const AtmField& atmospheric_field,
    const SurfaceField& surface_field,
    const Agenda& propagation_matrix_agenda,
    const Numeric& depolarization_factor,
    const Index& hse_derivative) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(jacobian_targets.x_size(),
                     "Cannot have any Jacobian targets")

  const Size np = ray_path.size();
  ARTS_USER_ERROR_IF(
      np != ray_path_propagation_matrix_scattering.size(),
      "Bad ray_path_propagation_matrix_scattering: incorrect number of path points")
  ARTS_USER_ERROR_IF(np != ray_path_suns_path.size(),
                     "Bad ray_path_suns_path: incorrect number of path points")

  const Size nsuns = suns.size();
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(ray_path_suns_path,
                          Cmp::ne(nsuns),
                          &ArrayOfArrayOfPropagationPathPoint::size),
      "Bad ray_path_suns_path: incorrect number of suns")

  const Size nf = frequency_grid.size();
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(ray_path_spectral_radiance_scattering,
                          Cmp::ne(nf),
                          [](auto& v) { return v.size(); }),
      "Bad ray_path_spectral_radiance_source: incorrect number of frequencies")

  StokvecVector spectral_radiance{};
  StokvecMatrix spectral_radiance_jacobian{};
  StokvecVector spectral_radiance_background{};
  const StokvecMatrix spectral_radiance_background_jacobian(0, nf);

  ray_path_spectral_radiance_scattering.resize(np);
  for (auto& p : ray_path_spectral_radiance_scattering) {
    p.resize(nf);
    p = 0;
  }

  String error{};

#pragma omp parallel for firstprivate(    \
        spectral_radiance,                \
            spectral_radiance_jacobian,   \
            spectral_radiance_background, \
            spectral_radiance_background_jacobian) if (not arts_omp_in_parallel())
  for (Size ip = 0; ip < np; ip++) {
    try {
      auto& spectral_radiance_scattered =
          ray_path_spectral_radiance_scattering[ip];

      const auto& propagation_matrix_scattering =
          ray_path_propagation_matrix_scattering[ip];
      const auto& ray_path_point = ray_path[ip];
      const auto& suns_path      = ray_path_suns_path[ip];

      for (Size isun = 0; isun < nsuns; isun++) {
        const auto& sun_path = suns_path[isun];
        const auto& sun      = suns[isun];

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
            hse_derivative);

        ARTS_USER_ERROR_IF(
            spectral_radiance.size() != nf,
            "Bad size spectral_radiance (",
            spectral_radiance.size(),
            ").  It should have the same size as frequency_grid (",
            nf,
            ")")

        // irradiance ratio
        const Numeric radiance_2_irradiance =
            pi * suns[isun].sin_alpha_squared(sun_path.back().pos,
                                              surface_field.ellipsoid);

        const Muelmat scatmat =
            rtepack::rayleigh_scattering(sun_path.front().los,
                                         ray_path_point.los,
                                         depolarization_factor) /
            (4 * pi);

        // Add the source to the target
        for (Size iv = 0; iv < nf; iv++) {
          spectral_radiance_scattered[iv] += propagation_matrix_scattering[iv] *
                                             scatmat * radiance_2_irradiance *
                                             spectral_radiance[iv];
        }
      }
    } catch (const std::exception& e) {
#pragma omp critical
      error += e.what();
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "{}", error)
}
ARTS_METHOD_ERROR_CATCH
