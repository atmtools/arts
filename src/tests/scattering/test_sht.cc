#include <scattering/sht.h>

#include <complex>
#include <iostream>
#include <random>

#include "test_utils.h"

using namespace scattering;

///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Test functions
///////////////////////////////////////////////////////////////////////////////

bool test_initialize_sht() {
  // SHT with l_max = 5, m_max = 5, n_lon = 32, n_lat = 32.
  auto sht   = sht::SHT(5, 5, 32, 32);
  auto sht_p = *sht::provider.get_instance({5, 5, 32, 32});

  if (sht.get_n_spectral_coeffs() != sht_p.get_n_spectral_coeffs()) {
    return false;
  }
  if (sht.get_n_latitudes() != sht_p.get_n_latitudes()) {
    return false;
  }
  if (sht.get_n_longitudes() != sht_p.get_n_longitudes()) {
    return false;
  }

  auto sht_lonlat = *sht::provider.get_instance_lonlat(32, 32);
  auto sht_lm     = *sht::provider.get_instance_lm(15, 15);
  if (sht_lonlat.get_n_latitudes() != sht_lm.get_n_latitudes()) {
    return false;
  }
  if (sht_lonlat.get_n_longitudes() != sht_lm.get_n_longitudes()) {
    return false;
  }

  return true;
}

/** Calculates transforms of specific spherical harmonics and checks
 * that the resulting coefficient are zero almost everywhere except
 * at the coefficients that correspond to the given SH function.
 */
bool test_transform_harmonics() {
  auto sht              = *sht::provider.get_instance({5, 5, 32, 32});
  auto lons             = sht.get_longitude_grid(true);
  auto lats             = sht.get_latitude_grid(true);
  Matrix spatial_coeffs = static_cast<Matrix>(sht.get_spatial_coeffs());
  ComplexVector spectral_coeffs =
      static_cast<ComplexVector>(sht.get_spectral_coeffs());

  for (int l = 0; l < 3; ++l) {
    for (int m = -l; m <= l; ++m) {
      spatial_coeffs   = evaluate_spherical_harmonic(l, m, lons, lats);
      spectral_coeffs  = sht.transform(spatial_coeffs);
      auto coeff_index = sht.get_coeff_index(l, m);
      for (int i = 0; i < spectral_coeffs.size(); ++i) {
        double diff = 0.0;
        if (i == coeff_index) {
          if (m < 0) {
            diff = spectral_coeffs[i].real() - 1.0 * pow(-1.0, m);
          } else {
            diff = spectral_coeffs[i].real() - 1.0;
          }
        } else {
          diff = std::abs(spectral_coeffs[i]);
        }
        if (std::abs(diff) > 1e-6) {
          return false;
        }
      }
    }
  }
  return true;
}

/** Test symmetry.
 *
 * Ensure that applying the transform twice doesn't change
 * coefficients.
 */
bool test_symmetry(int n_trials) {
  for (int i = 0; i < n_trials; ++i) {
    auto sht_v          = *sht::provider.get_instance({16, 16, 34, 34});
    ComplexVector v     = random_spectral_coeffs(16, 16);
    auto v_spat         = sht_v.synthesize(v);
    ComplexVector v_ref = sht_v.transform(v_spat);
    auto diff           = max_error(v, v_ref);
    if (std::abs(diff) > 1e-6) {
      return false;
    }
  }
  return true;
}

/** Test addition of spectral coefficients.
 *
 * Calculates the sum of 2D fields in spatial and spectral
 * space and ensures that the spectral representation of
 * the sum is the same.
 */
bool test_add_coefficients(int n_trials) {
  for (int i = 0; i < n_trials; ++i) {
    auto sht_w = *sht::provider.get_instance({8, 8, 34, 34});
    auto sht_v = *sht::provider.get_instance({16, 16, 34, 34});

    ComplexVector v  = random_spectral_coeffs(16, 16);
    ComplexVector w  = random_spectral_coeffs(8, 8);
    ComplexVector vw = sht::add_coeffs(sht_v, v, sht_w, w);

    auto vw_spat  = sht_v.synthesize(v);
    vw_spat      += sht_w.synthesize(w);

    ComplexVector vw_ref = sht_v.transform(vw_spat);

    auto diff = max_error(vw_ref, vw);
    if (std::abs(diff) > 1e-6) {
      return false;
    }
  }
  return true;
}

bool test_grids() {
  auto sht = sht::provider.get_instance_lonlat(64, 64);

  auto lat_grid     = sht->get_latitude_grid();
  Numeric max_angle = max<Vector>(lat_grid);
  if (max_angle < 2.0 * scattering::sht::pi_v<Numeric>) {
    return false;
  }
  lat_grid  = sht->get_latitude_grid(true);
  max_angle = max<Vector>(lat_grid);
  if (max_angle > 2.0 * scattering::sht::pi_v<Numeric>) {
    return false;
  }

  auto lon_grid = sht->get_longitude_grid();
  max_angle     = max<Vector>(lon_grid);
  if (max_angle < 2.0 * scattering::sht::pi_v<Numeric>) {
    return false;
  }
  lon_grid  = sht->get_longitude_grid(true);
  max_angle = max<Vector>(lon_grid);
  if (max_angle > 2.0 * scattering::sht::pi_v<Numeric>) {
    return false;
  }

  max_angle = max<Vector>(*sht->get_za_grid_ptr());
  if (max_angle < 2.0 * scattering::sht::pi_v<Numeric>) {
    return false;
  }

  max_angle = max<Vector>(*sht->get_aa_grid_ptr());
  if (max_angle < 2.0 * scattering::sht::pi_v<Numeric>) {
    return false;
  }

  return true;
}

int main(int /*nargs*/, char** /*argv*/) {
#ifndef ARTS_NO_SHTNS
  bool passed = test_initialize_sht();
  std::cout << "test_initialization: ";
  if (passed) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
    return 1;
  }

  passed = test_transform_harmonics();
  std::cout << "test_transform_harmonics: ";
  if (passed) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
    return 1;
  }

  passed = test_symmetry(10);
  std::cout << "test_symmetry: ";
  if (passed) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
    return 1;
  }

  passed = test_add_coefficients(10);
  std::cout << "test_add_coefficients: ";
  if (passed) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
    return 1;
  }

  passed = test_grids();
  std::cout << "test_grids: ";
  if (passed) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
    return 1;
  }
#endif

  return 0;
}
