#include <arts_conversions.h>
#include <arts_omp.h>
#include <disort.h>
#include <workspace.h>
#pragma clang - format off
#if not ARTS_LGPL
#include <cdisort.h>
#endif
#pragma clang - format on

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

#include "arts_constants.h"
#include "configtypes.h"
#include "debug.h"
#include "path_point.h"

////////////////////////////////////////////////////////////////////////
// Core Disort
////////////////////////////////////////////////////////////////////////

#if not ARTS_LGPL
namespace {
void setup_cdisort(disort_state& ds,
                   const Vector& phis,
                   const DisortSettings& disort_settings,
                   const disort::main_data& dis,
                   const Numeric surface_temperature,
                   const bool intensity_correction = false) {
  // solar dependent properties if no sun is present
  // Number of azimuth angles
  Index nphi = 1;
  //local zenith angle of sun
  Numeric umu0 = 0.;
  //local azimuth angle of sun
  Numeric phi0 = 0.;

  // FIXME OLE Move into frequency loop? deciding based on first frequency here
  if (disort_settings.solar_source[0] >= 0) {
    nphi = phis.size();
    umu0 = Conversion::cosd(dis.solar_zenith());
    phi0 = dis.beam_azimuth();
    if (phi0 < 0) {
      phi0 = phi0 + 360.;
    }
  }

  ds.accur        = 0.005;
  ds.flag.prnt[0] = FALSE;
  ds.flag.prnt[1] = FALSE;
  ds.flag.prnt[2] = FALSE;
  ds.flag.prnt[3] = FALSE;
  ds.flag.prnt[4] = FALSE;

  ds.flag.usrtau         = FALSE;
  ds.flag.usrang         = TRUE;  // Toggle for custom zenith angles
  ds.flag.spher          = FALSE;
  ds.flag.general_source = FALSE;
  ds.flag.output_uum     = FALSE;

  ds.nlyr = static_cast<int>(disort_settings.nlay);  // pressure.nelem() - 1

  ds.flag.brdf_type = BRDF_NONE;

  ds.flag.ibcnd = GENERAL_BC;

  if (disort_settings.source_polynomial.ncols() > 0) {
    ds.flag.planck = TRUE;
  } else {
    ds.flag.planck = FALSE;
  }
  ds.flag.onlyfl = FALSE;
  ds.flag.lamber = TRUE;
  ds.flag.quiet  = FALSE;
  if (intensity_correction) {
    ds.flag.intensity_correction     = TRUE;
    ds.flag.old_intensity_correction = FALSE;
  } else {
    ds.flag.intensity_correction     = FALSE;
    ds.flag.old_intensity_correction = FALSE;
  }

  ds.nstr   = static_cast<int>(disort_settings.quadrature_dimension);
  ds.nphase = ds.nstr;
  ds.nmom   = ds.nstr;
  //ds.ntau = ds.nlyr + 1;   // With ds.flag.usrtau = FALSE; set by cdisort
  ds.numu = static_cast<int>(
      disort_settings.quadrature_dimension);  // za_grid.nelem();
  ds.nphi = static_cast<int>(nphi);

  // Looking direction of solar beam
  ds.bc.umu0 = umu0;
  ds.bc.phi0 = phi0;

  // Intensity of bottom-boundary isotropic illumination
  ds.bc.fluor = 0.;

  // Top of the atmosphere temperature and emissivity
  ds.bc.ttemp = Constant::cosmic_microwave_background_temperature;
  ds.bc.btemp = surface_temperature;
  ds.bc.temis = 1.;
}

void freq_setup_cdisort(disort_state& ds,
                        const disort::main_data& dis,
                        const Vector& phis,
                        const ArrayOfAtmPoint& ray_path_atmospheric_point,
                        const Numeric& frequency) {
  // fill up azimuth angle and temperature array
  for (Index i = 0; i < ds.nphi; i++) ds.phi[i] = phis[i];  // aa_grid[i];

  if (ds.flag.planck == TRUE) {
    for (Index i = 0; i <= ds.nlyr; i++)
      ds.temper[i] = ray_path_atmospheric_point[i].temperature;
  }

  // for (Index i = 0; i < ds.numu; i++) {
  //   ds.umu[i] = dis.mu()[i];  // -cos(za_grid[i] * PI / 180);
  // }

  for (Index i = 0; i < ds.numu / 2; i++) {
    ds.umu[i] = dis.mu()[ds.numu - i - 1];  // -cos(za_grid[i] * PI / 180);
    ds.umu[i + ds.numu / 2] = dis.mu()[i];  // -cos(za_grid[i] * PI / 180);
  }


  // Cosmic background
  // we use temis*ttemp as upper boundary specification, hence CBR set to 0.
  ds.bc.fisot = 0;

  std::vector<Numeric> dtauc(ds.nlyr);
  const auto& tau = dis.tau();
  std::adjacent_difference(tau.begin(), tau.end(), dtauc.begin());

  std::memcpy(ds.dtauc, dtauc.data(), sizeof(Numeric) * ds.nlyr);
  std::memcpy(ds.ssalb, dis.omega().data_handle(), sizeof(Numeric) * ds.nlyr);

  // Wavenumber in [1/cm]
  ds.wvnmhi = ds.wvnmlo  = frequency / (100. * Constant::c);
  ds.wvnmhi             += ds.wvnmhi * 1e-7;
  ds.wvnmlo             -= ds.wvnmlo * 1e-7;

  // set surface albedo
  if (dis.brdf_modes().size() > 0)
    ds.bc.albedo = dis.brdf_modes()[0]({0}, {0})[0, 0];
  else
    ds.bc.albedo = 0.0;

  ds.bc.fbeam = dis.beam_source();

  std::memcpy(ds.pmom,
              dis.all_legendre_coeffs().data_handle(),
              sizeof(Numeric) * dis.all_legendre_coeffs().size());
}

void run_cdisort(Tensor3View disort_spectral_radiance_field,
                 disort_state& ds,
                 disort_output& out) {
  Numeric umu0 = 0.;
  enum class Status { FIRST_TRY, RETRY, SUCCESS };
  Status tries      = Status::FIRST_TRY;
  const Numeric eps = 2e-4;  //two times the value defined in cdisort.c:3653
  do {
    try {
      c_disort(&ds, &out);
      tries = Status::SUCCESS;
    } catch (const std::runtime_error& e) {
      //catch cases if solar zenith angle=quadrature angle
      if (tries == Status::FIRST_TRY) {
        // change angle
        if (umu0 < 1 - eps) {
          umu0 += eps;
        } else if (umu0 > 1 - eps) {
          umu0 -= eps;
        }

        const Numeric shift =
            abs(Conversion::acosd(umu0) - Conversion::acosd(ds.bc.umu0));
        std::cerr
            << "Solar zenith angle coincided with one of the quadrature angles\n"
            << "We needed to shift the solar sun angle by " << shift
            << "deg.\n";

        ds.bc.umu0 = umu0;
        tries      = Status::RETRY;
      } else
        throw e;
    }
  } while (tries != Status::SUCCESS);

  for (Index k = 1; k <= ds.nlyr; k++) {
    for (Index i = 0; i < ds.nphi; i++) {
      for (Index j = 0; j < ds.numu / 2; j++) {
        disort_spectral_radiance_field[k - 1, i, ds.numu - j - 1] =
            out.uu[j + (k + i * (ds.nlyr + 1)) * ds.numu] /
            (ds.wvnmhi - ds.wvnmlo) / (100 * Constant::c);
        disort_spectral_radiance_field[k - 1, i, j] =
            out.uu[j + ds.numu / 2 + (k + i * (ds.nlyr + 1)) * ds.numu] /
            (ds.wvnmhi - ds.wvnmlo) / (100 * Constant::c);
      }
    }
  }
}

}  // namespace
#endif  // ARTS_LGPL

void cdisort_spectral_radiance_fieldCalc(
    Tensor4& disort_spectral_radiance_field,
    Vector& disort_quadrature_angles,
    Vector& disort_quadrature_weights,
    const DisortSettings& disort_settings,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const ArrayOfAscendingGrid& ray_path_frequency_grid,
    const ArrayOfPropagationPathPoint& ray_path,
    const SurfaceField& surface_field,
    const Vector& phis) {
  ARTS_TIME_REPORT
#if not ARTS_LGPL

  using Conversion::acosd;

  const Index nv    = disort_settings.nfreq;
  const Index np    = disort_settings.nlay;
  const Index nquad = disort_settings.quadrature_dimension;

  disort_spectral_radiance_field.resize(nv, np, phis.size(), nquad);

  disort::main_data dis = disort_settings.init();

  const Numeric surface_temperature =
      surface_field.single_value(SurfaceKey::t,
                                 ray_path[ray_path.size() - 1].latitude(),
                                 ray_path[ray_path.size() - 1].longitude());

  //! Supplementary outputs
  disort_quadrature_weights = dis.weights();
  disort_quadrature_angles.resize(nquad);

  std::transform(dis.mu().begin(),
                 dis.mu().end(),
                 disort_quadrature_angles.begin(),
                 [](const Numeric& mu) { return acosd(mu); });

  disort_state ds;
  setup_cdisort(ds, phis, disort_settings, dis, surface_temperature);

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(dis, ds)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      disort_output out;
      /* Allocate memory */
      c_disort_state_alloc(&ds);
      c_disort_out_alloc(&ds, &out);

      disort_settings.set_cdisort(dis, iv);

      freq_setup_cdisort(ds,
                         dis,
                         phis,
                         ray_path_atmospheric_point,
                         ray_path_frequency_grid[0][iv]);

      run_cdisort(disort_spectral_radiance_field[iv], ds, out);

      /* Free allocated memory */
      c_disort_out_free(&ds, &out);
      c_disort_state_free(&ds);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(
      error.size(), "Error occurred in disort-spectral:\n{}", error);
#else
  ARTS_USER_ERROR("Did not compile with -DENABLE_ARTS_LGPL=0")
#endif  // ARTS_LGPL
}
