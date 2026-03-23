/*===========================================================================
  ===  File description
  ===========================================================================*/

#include <arts_omp.h>
#include <atm.h>
#include <compare.h>
#include <configtypes.h>
#include <debug.h>
#include <jacobian.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <scattering_species.h>
#include <sun.h>
#include <sun_methods.h>
#include <workspace.h>

#include <algorithm>

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
  sun.spectrum    = regrid_sun_spectrum(sun_spectrum_raw, f_grid, temperature);
  sun.description = description;
  sun.radius      = radius;
  sun.distance    = distance;
  sun.latitude    = latitude;
  sun.longitude   = longitude;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void sunBlackbody(Sun& sun,
                  // Inputs:
                  const AscendingGrid& freq_grid,
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
  sun.spectrum = Matrix(freq_grid.size(), 4, 0.);

  planck(sun.spectrum[joker, 0], freq_grid, temperature);
  sun.spectrum *= pi;  // outgoing flux at the surface of the sun.

  sun.description = "Blackbody sun";
  sun.radius      = radius;
  sun.distance    = distance;
  sun.latitude    = latitude;
  sun.longitude   = longitude;
}

void sun_pathFromObserverAgenda(const Workspace& ws,
                                ArrayOfPropagationPathPoint& sun_path,
                                const SurfaceField& surf_field,
                                const Agenda& ray_path_observer_agenda,
                                const Sun& sun,
                                const Vector3& observer_pos,
                                const Numeric& angle_cut,
                                const Index& refinements,
                                const Index& just_hit) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  find_sun_path(ws,
                sun_path,
                sun,
                ray_path_observer_agenda,
                surf_field,
                observer_pos,
                angle_cut,
                refinements,
                just_hit);
}

void ray_path_suns_pathFromPathObserver(
    const Workspace& ws,
    ArrayOfArrayOfArrayOfPropagationPathPoint& ray_path_suns_path,
    const SurfaceField& surf_field,
    const Agenda& ray_path_observer_agenda,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfSun& suns,
    const Numeric& angle_cut,
    const Index& refinements,
    const Index& just_hit) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

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
                                   surf_field,
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
                                     surf_field,
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

void spectral_propmat_scatInit(PropmatVector& spectral_propmat_scat,
                               const AscendingGrid& freq_grid) {
  ARTS_TIME_REPORT

  spectral_propmat_scat.resize(freq_grid.size());
  spectral_propmat_scat = 0.0;
}

void spectral_propmat_scatAirSimple(PropmatVector& spectral_propmat_scat,
                                    const AscendingGrid& freq_grid,
                                    const AtmPoint& atm_point) {
  const Size nf = freq_grid.size();
  ARTS_USER_ERROR_IF(spectral_propmat_scat.size() != nf,
                     "Mismatch in size of spectral_propmat_scat and freq_grid")

  constexpr scattering::RayleighScatterer rayleigh{};
  for (Size f = 0; f < nf; f++) {
    const auto [sca, abs] = rayleigh.cross_sections(freq_grid[f], atm_point);
    spectral_propmat_scat[f].A() += sca + abs;
  }
}

void spectral_propmat_scat_pathFromPath(
    const Workspace& ws,
    ArrayOfPropmatVector& spectral_propmat_scat_path,
    const Agenda& spectral_propmat_scat_agenda,
    const ArrayOfAscendingGrid& freq_grid_path,
    const ArrayOfAtmPoint& atm_path) {
  ARTS_TIME_REPORT

  const Size np = freq_grid_path.size();
  ARTS_USER_ERROR_IF(np != atm_path.size(),
                     "Bad atm_path: incorrect number of path points")

  spectral_propmat_scat_path.resize(np);
  if (arts_omp_in_parallel()) {
    for (Size ip = 0; ip < np; ip++) {
      spectral_propmat_scat_agendaExecute(ws,
                                          spectral_propmat_scat_path[ip],
                                          freq_grid_path[ip],
                                          atm_path[ip],
                                          spectral_propmat_scat_agenda);
    }
  } else {
    String error{};
#pragma omp parallel for
    for (Size ip = 0; ip < np; ip++) {
      try {
        spectral_propmat_scat_agendaExecute(ws,
                                            spectral_propmat_scat_path[ip],
                                            freq_grid_path[ip],
                                            atm_path[ip],
                                            spectral_propmat_scat_agenda);
      } catch (const std::exception& e) {
#pragma omp critical
        error += e.what();
      }
    }

    ARTS_USER_ERROR_IF(error.size(), "{}", error)
  }
}

void spectral_rad_srcvec_pathAddScattering(
    SourceVector& spectral_rad_srcvec_path,
    const ArrayOfStokvecVector& spectral_rad_scat_path,
    const ArrayOfPropmatVector& spectral_propmat_path) {
  ARTS_TIME_REPORT

  const Size np = spectral_propmat_path.size();
  ARTS_USER_ERROR_IF(
      np != static_cast<Size>(spectral_rad_srcvec_path.J.ncols()),
      "Bad spectral_propmat_scat_path: incorrect number of path points")

  ARTS_USER_ERROR_IF(
      np != spectral_rad_scat_path.size(),
      "Bad spectral_propmat_scat_path: incorrect number of path points")

  if (np == 0) return;

  const Size nf = spectral_propmat_path.front().size();
  ARTS_USER_ERROR_IF(
      nf != static_cast<Size>(spectral_rad_srcvec_path.J.nrows()),
      "Mismatch frequency size of spectral_propmat_path and spectral_rad_srcvec_path")

  ARTS_USER_ERROR_IF(
      stdr::any_of(spectral_rad_scat_path,
                   Cmp::ne(nf),
                   [](auto& v) { return v.size(); }),
      "Mismatch frequency size of spectral_propmat_path and spectral_rad_scat_path")

#pragma omp parallel for collapse(2) if (!arts_omp_in_parallel())
  for (Size ip = 0; ip < np; ip++) {
    for (Size iv = 0; iv < nf; iv++) {
      spectral_rad_srcvec_path.J[iv, ip] +=
          inv(spectral_propmat_path[ip][iv]) * spectral_rad_scat_path[ip][iv];
    }
  }
}

void spectral_rad_scat_pathSunsFirstOrder(
    const Workspace& ws,
    // [np, nf]:
    ArrayOfStokvecVector& spectral_rad_scat_path,
    // [np]:
    const ArrayOfScatteringSpecies& scat_species,
    // [np]:
    const ArrayOfAtmPoint& atm_path,
    // [np]:
    const ArrayOfPropagationPathPoint& ray_path,
    // [np, suns, np2]:
    const ArrayOfArrayOfArrayOfPropagationPathPoint& ray_path_suns_path,
    // [nsuns]:
    const ArrayOfSun& suns,
    // [njac]:
    const JacobianTargets& jac_targets,
    // [nf]:
    const AscendingGrid& freq_grid,
    const AtmField& atm_field,
    const SurfaceField& surf_field,
    const Agenda& spectral_propmat_agenda,
    const TransmittanceOption& rte_option,
    const Index& hse_derivative) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  ARTS_USER_ERROR_IF(jac_targets.x_size(), "Cannot have any Jacobian targets")

  const Size np = ray_path.size();
  ARTS_USER_ERROR_IF(np != atm_path.size(),
                     "Bad atm_path: incorrect number of path points")
  ARTS_USER_ERROR_IF(np != ray_path_suns_path.size(),
                     "Bad ray_path_suns_path: incorrect number of path points")

  const Size nsuns = suns.size();
  ARTS_USER_ERROR_IF(stdr::any_of(ray_path_suns_path,
                                  Cmp::ne(nsuns),
                                  &ArrayOfArrayOfPropagationPathPoint::size),
                     "Bad ray_path_suns_path: incorrect number of suns")

  const Size nf = freq_grid.size();
  ARTS_USER_ERROR_IF(
      stdr::any_of(spectral_rad_scat_path,
                   Cmp::ne(nf),
                   [](auto& v) { return v.size(); }),
      "Bad spectral_rad_scat_path: incorrect number of frequencies")

  StokvecVector spectral_rad{};
  StokvecMatrix spectral_rad_jac{};
  StokvecVector spectral_rad_bkg{};
  const StokvecMatrix spectral_rad_bkg_jac(0, nf);

  spectral_rad_scat_path.resize(np);
  for (auto& p : spectral_rad_scat_path) {
    p.resize(nf);
    p = 0;
  }

  String error{};

#pragma omp parallel for firstprivate( \
        spectral_rad,                  \
            spectral_rad_jac,          \
            spectral_rad_bkg,          \
            spectral_rad_bkg_jac) if (not arts_omp_in_parallel())
  for (Size ip = 0; ip < np; ip++) {
    try {
      auto& spectral_rad_scattered = spectral_rad_scat_path[ip];

      const auto& ray_point = ray_path[ip];
      const auto& suns_path = ray_path_suns_path[ip];

      for (Size isun = 0; isun < nsuns; isun++) {
        const auto& sun_path = suns_path[isun];
        const auto& sun      = suns[isun];

        spectral_radSunOrCosmicBackground(
            spectral_rad_bkg, freq_grid, sun_path, sun, surf_field);

        spectral_radClearskyBackgroundTransmission(ws,
                                                   spectral_rad,
                                                   spectral_rad_jac,
                                                   atm_field,
                                                   freq_grid,
                                                   jac_targets,
                                                   sun_path,
                                                   rte_option,
                                                   spectral_propmat_agenda,
                                                   spectral_rad_bkg,
                                                   spectral_rad_bkg_jac,
                                                   surf_field,
                                                   hse_derivative);

        ARTS_USER_ERROR_IF(spectral_rad.size() != nf,
                           "Bad size spectral_rad (",
                           spectral_rad.size(),
                           ").  It should have the same size as freq_grid (",
                           nf,
                           ")")

        // irradiance ratio
        const Numeric radiance_2_irradiance =
            pi * suns[isun].sin_alpha_squared(sun_path.back().pos,
                                              surf_field.ellipsoid);

        // Scattering phase matrix at the angle between sun and ray.
        const auto Z = scat_species.get_phase_matrix_at_angle(
            atm_path[ip], freq_grid, sun_path.front().los, ray_point.los);

        for (Size iv = 0; iv < nf; iv++) {
          spectral_rad_scattered[iv] +=
              Z[iv] * radiance_2_irradiance * spectral_rad[iv];
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
