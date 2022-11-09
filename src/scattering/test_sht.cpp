#include <iostream>
#include <math.h>

#include <scattering/sht.h>


using namespace scattering::sht;
using namespace scattering::math;

///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

// Helper function to calculate factorial.
template <typename Scalar>
Scalar fac(int n) {
    if (n <= 0) {
        return 1;
    } else if (n == 1) {
        return 1;
    }
    return static_cast<Scalar>(n) * fac<Scalar>(n - 1);
}

/** Evaluate spherical harmonic.
 *
 * @param l The l parameter
 * @param m The m parameter
 * @param lons The longitude grid.
 * @param lats The latitude grid.
 *
 * @return A matrix containing the evaluated spherical harmonics
 * function with the latitudes oriented along columns and the
 * longitudes along rows.
 */
GridCoeffs evaluate_spherical_harmonic(
    int l,
    int m,
    scattering::math::Vector<double> lons,
    scattering::math::Vector<double> lats
) {

    auto n_lon = lons.cols();
    auto n_lat = lats.cols();

    double norm = sqrt((2 * l + 1.0) / (4.0 * M_PI));
    norm *= sqrt(fac<double>(l - m) / fac<double>(l + m));

    scattering::math::Matrix<double> results{n_lon, n_lat};
    for (int i = 0; i < n_lon; ++i) {
        for (int j = 0; j < n_lat; ++j) {
            auto clon = cos(m * lons[i]);

            auto clat = cos(lats[j]);
            double p = 0.0;
            if (l == 0) {
                p = 1.0;
            } else if (l == 1) {
                if (m == -1) {
                    p = 1.0 * sqrt(1.0 - clat * clat);
                } else if (m == 0) {
                    p = clat;
                } else if (m == 1) {
                    p = - 2.0 * sqrt(1.0 - clat * clat);
                } else {
                    throw std::runtime_error(
                        "m must be within [-l, l]."
                        );
                }
            } else if (l == 2) {
                if (m == -2) {
                    p = (1 - clat * clat) / 4;
                } else if (m == -1) {
                    p = 1.0 * clat * sqrt(1.0 - clat * clat);
                } else if (m == 0) {
                    p = 0.5 * (3 * clat * clat - 1);
                } else if (m == 1) {
                    p = -6 * clat * sqrt(1.0 - clat * clat);
                } else if (m == 2) {
                    p = 6 * (1 - clat * clat);
                }
            }
            if (l > 2) {
                throw std::runtime_error("l must be less than or equal to 1.");
            }
            results(i, j) = norm * p * clon;
        }
    }
    return results;
}

SpectralCoeffs random_spectral_coeffs(Index l_max, Index m_max) {
    SHT sht{l_max, m_max};
    auto n_coeffs = sht.get_n_spectral_coeffs();
    SpectralCoeffs coeffs = SpectralCoeffs::Random(1, n_coeffs);
    // Imaginary part of coeffs with m == 0 must be 0.
    for (int i = 0; i < l_max + 1; ++i) {
        coeffs[i] = coeffs[i].real();
    }
    return coeffs;
}

///////////////////////////////////////////////////////////////////////////////
// Test functions
///////////////////////////////////////////////////////////////////////////////

/** Calculates transforms of specific spherical harmonics and checks
 * that the resulting coefficient are zero almost everywhere except
 * at the coefficients that correspond to the given SH function.
 */
bool test_transform_harmonics() {

    auto sht = SHT{5, 5, 32, 32};
    auto lons = sht.get_longitude_grid();
    auto lats = sht.get_latitude_grid();
    auto spatial_coeffs = sht.get_spatial_coeffs();
    auto spectral_coeffs = sht.get_spectral_coeffs();

    for (int l = 0; l < 3; ++l) {
        for (int m = -l; m <= l; ++m) {
            spatial_coeffs = evaluate_spherical_harmonic(
                l, m, lons, lats
                );
            spectral_coeffs = sht.transform(spatial_coeffs);
            auto coeff_index = sht.get_coeff_index(l, m);
            for (int i = 0; i < spectral_coeffs.cols(); ++i) {
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
                if (!small(diff)) {
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

        auto sht_v = SHT{16, 16, 34, 34};
        SpectralCoeffs v = random_spectral_coeffs(16, 16);
        auto v_spat = sht_v.synthesize(v);
        SpectralCoeffs v_ref = sht_v.transform(v_spat);
        auto diff = (v_ref - v).cwiseAbs().maxCoeff();
        if (!small(diff)) {
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

      auto sht_w = SHT{8, 8, 34, 34};
      auto sht_v = SHT{16, 16, 34, 34};

      SpectralCoeffs v = random_spectral_coeffs(16, 16);
      SpectralCoeffs w = random_spectral_coeffs(8, 8);
      SpectralCoeffs vw = SHT::add_coeffs(sht_v, v, sht_w, w);

      auto v_spat = sht_v.synthesize(v);
      auto w_spat = sht_w.synthesize(w);
      SpectralCoeffs vw_ref = sht_v.transform(v_spat + w_spat);

      auto diff = (vw_ref - vw).cwiseAbs().maxCoeff();
      if (!small(diff)) {
          return false;
      }
    }
    return true;
}


int main(int /*nargs*/, char **/*argv*/) {

    auto passed = test_transform_harmonics();
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


    return 0;

}
