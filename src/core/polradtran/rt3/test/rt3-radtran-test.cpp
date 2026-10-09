// The C++ RADTRAN (radtran3.cc) against the Fortran RADTRAN (radtran3.f) it
// ports: both get the same inputs, and every output must agree to
// tolerance, relative to the largest magnitude in that output.  The test
// reports how many cases are bit-identical and the largest difference.
//
// Mathematically equivalent evaluations are accepted: the tolerances allow
// a few ulp (rounding, 16 epsilon) times what a computation amplifies
// rounding by, so that FMA contraction, vectorised libm functions and other
// BLAS kernels pass and an error of the port does not.  n doublings
// amplify rounding by 2^n, so the radiances and fluxes must agree to 1e-10
// (for the replaced quadratures and Planck function) plus rounding times
// 2^n for the layer doubled most; they differ by up to 1.2e-9 on AMD x86_64
// with MKL.
//
// Each porting step is first checked bit for bit against the previous one.
// The port was bit-identical to the Fortran until RT3's quadratures were
// replaced by ARTS's (polradtran::get_quadrature), which differ by rounding: the
// weights by up to 2.4e-12 relative, RT3's being the less accurate.  The
// radiances mostly change by 1e-15 to 1e-14 of the largest one, but by
// 1.7e-11 in an optically thick, strongly scattering case: a node 3 ulp
// lower flips a pivot of LINPACK's DGEFA in the doubling (1 or 2 ulp do
// not), and the inverse is ill-conditioned there.
//
// It also checks fft1dr's documented format against direct sums (what an
// FFTW build must also give), polradtran::get_quadrature against RT3's quadrature
// routines, thermal_radiance against RT3's (to 2e-12: RT3's Planck
// function loses digits at small h nu / k T), ground_surface against the
// Fortran grounds that RADTRAN used to make itself, doubling_integration,
// combine_layers and internal_radiance against RT3's (to rounding times
// 2^n kappa for n doublings that each invert a 1 - R R of condition number
// kappa, 1e-13 and 1e-13: LAPACK's inverse in place of LINPACK's, and the
// products of BLAS), check_norm's limit, and each ported routine against its Fortran:
// to rounding (1e-14 of the largest value; 1e-13 for scattering and
// direct_scattering, whose Legendre sums amplify the rounding of the
// scattering angle; 2e-12 for the ground radiances, with ARTS's Planck
// function), and exactly for the integer and copying ones (number_sums,
// matrix_symmetry, get_scattering, scatter_symmetry, get_direct).  The
// port is not bitwise: the C++ evaluates in its natural order and lets the
// compiler contract multiply-adds into FMAs.  The test reports how many
// cases are bit-identical for information.  Each RADTRAN case also runs
// again with one rt3_workdata reused over all cases and with one whose
// every array is NaN, which must not change a bit.
//
// The inputs are random, in RADTRAN's own layouts, and cover each branch of
// RADTRAN: 1 to 4 Stokes components, the quadratures (with extra angles,
// QUAD_TYPE 'E'), azimuth orders, delta-M, the source codes (none, solar,
// thermal, both), the Lambertian and Fresnel grounds, non-scattering and
// scattering layers (none, one and many doublings, and layers sharing a
// set), Rayleigh-, Mie- and general phase matrices (each of RT3's
// summation cases), negative gas extinction, output levels in any order,
// and the size limits.
#include <arts_constants.h>
#include <lin_alg.h>
#include <radintg.h>
#include <radintg3.h>
#include <radscat3.h>
#include <radtran3.h>
#include <radutil.h>
#include <radutil3.h>
#include <rt3.h>
#include <rt3_c_interface.h>
#include <rt3_fft.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace rt3 = polradtran::rt3;

namespace {
using polradtran::quadrature_type;

constexpr Numeric nan = std::numeric_limits<Numeric>::quiet_NaN();

//! The number of elements that differ, and the largest difference relative
//! to the largest magnitude in either array; NaN counts as infinitely far.
//! The elements are Numeric or Index.
std::pair<Index, Numeric> differ(const auto& a, const auto& b) {
  Index   count = 0;
  Numeric diff = 0.0, scale = 0.0;
  for (std::size_t k = 0; k < a.size(); k++) {
    const auto x = static_cast<Numeric>(a.data_handle()[k]);
    const auto y = static_cast<Numeric>(b.data_handle()[k]);
    scale        = std::max({scale, std::abs(x), std::abs(y)});
    if (x == y) continue;
    count++;
    const Numeric d = std::abs(x - y);
    diff            = std::max(diff, std::isnan(d) ? std::numeric_limits<Numeric>::infinity() : d);
  }
  return {count, count == 0 ? 0.0 : diff / scale};
}

/* polradtran::get_quadrature against RT3's quadrature routines, which it replaces
   in RADTRAN: the same rules, to rounding.  For nmu up to 64 the nodes
   differ by 4.4e-16 and the weights by 2.4e-12 relative (Apple arm64); the
   test allows 1e-15 and 1e-11 for other compilers and libms. */
void check_quadratures() {
  using fortran_rule = void (*)(std::int64_t, double*, double*);
  Numeric dmu = 0.0, dw = 0.0;
  for (Index n = 1; n <= 64; n++) {
    for (auto [type, fortran] :
         {std::pair<quadrature_type, fortran_rule>{quadrature_type::double_gauss, rt3_double_gauss_quadrature},
          std::pair<quadrature_type, fortran_rule>{quadrature_type::gauss, rt3_gauss_legendre_quadrature},
          std::pair<quadrature_type, fortran_rule>{quadrature_type::lobatto, rt3_lobatto_quadrature}}) {
      const auto q = polradtran::get_quadrature(n, type);
      Vector     mu(n), w(n);
      fortran(n, mu.data_handle(), w.data_handle());
      for (Index i = 0; i < n; i++) {
        dmu = std::max(dmu, std::abs(q.mu[i] - mu[i]));
        dw  = std::max(dw, std::abs(q.weights[i] - w[i]) / w[i]);
      }
    }
  }
  std::cout << std::format(
      "polradtran::get_quadrature against RT3's D, G and L routines, nmu 1 to 64: nodes within {:.2e}, weights within "
      "{:.2e} relative\n",
      dmu,
      dw);
  if (dmu > 1e-15 or dw > 1e-11) throw std::runtime_error("the quadratures differ by more than rounding");
}

/* The rounding of a few operations, per unit of amplification.  An
   evaluation that is mathematically the same but rounds differently (FMA
   contraction, a vectorised libm exp, another BLAS kernel or summation
   order, LAPACK's inverse for LINPACK's) changes a result by a few ulp
   times what the computation amplifies rounding by. */
constexpr Numeric rounding = 16 * std::numeric_limits<Numeric>::epsilon();

//! The differ() of one output, and how much its computation amplifies
//! rounding (1: not at all)
struct output_difference {
  Index   count{};
  Numeric rel{}, amplification{1.0};

  output_difference(std::pair<Index, Numeric> d, Numeric amp = 1.0)
      : count{d.first}, rel{d.second}, amplification{amp} {}
};

//! A ported routine against its Fortran, case by case
struct routine_tally {
  Index   ncase{}, identical{};
  Numeric worst{}, worst_per_amplification{};
  bool    amplified{false};

  //! One case, from the differ() of each output
  void add(std::initializer_list<output_difference> outputs) {
    Index ndiffer = 0;
    for (const auto& out : outputs) {
      ndiffer                 += out.count;
      worst                    = std::max(worst, out.rel);
      worst_per_amplification  = std::max(worst_per_amplification, out.rel / out.amplification);
      amplified                = amplified or out.amplification != 1.0;
    }
    identical += ndiffer == 0;
    ncase++;
  }

  //! Throws if an output differs by more than tolerance times its
  //! amplification
  void report(std::string_view cpp, std::string_view fortran, Numeric tolerance = 1e-14) const {
    std::cout << std::format("{} against {}: {} of {} cases bit-identical (largest relative difference {:.2e}{})\n",
                             cpp,
                             fortran,
                             identical,
                             ncase,
                             worst,
                             amplified ? std::format(", {:.2e} per unit of rounding amplification against {:.2e}",
                                                     worst_per_amplification,
                                                     tolerance)
                                       : std::string{});
    if (worst_per_amplification > tolerance)
      throw std::runtime_error(std::format("{} differs from {} by more than {:.2e}{}",
                                           cpp,
                                           fortran,
                                           tolerance,
                                           amplified ? " times the rounding amplification" : ""));
  }
};

/* The routines ported as is against their Fortran.  The C++ outputs start
   as NaN, the Fortran ones as 0, so an element a routine does not write
   differs. */
/* rt3::initialize against RT3_INITIALIZE, for 1 to 4 Stokes parameters,
   gauss nodes and the extra angle 1, albedos 0, 0.6 and 1, and random phase
   functions and thicknesses. */
void check_initialize() {
  std::mt19937_64                         gen(2013);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  const auto                              q = polradtran::get_quadrature(4, quadrature_type::gauss);
  Vector                                  mu(5);
  mu[Range{0, 4}] = q.mu;
  mu[4]           = 1.0;
  for (Index nstokes : {1, 2, 3, 4}) {
    const Index nummu = mu.size(), n = nstokes * nummu;
    for (Numeric albedo : {0.0, 0.6, 1.0}) {
      for (Index rep = 0; rep < 4; rep++) {
        const Numeric delta_z = 1e-3 * u(gen), extinction = 1e-3 * u(gen);
        Tensor5       phase(4, nummu, nstokes, nummu, nstokes);
        for (std::size_t e = 0; e < phase.size(); e++) phase.data_handle()[e] = 0.2 * u(gen) - 0.05;
        Tensor5 r(2, nummu, nstokes, nummu, nstokes, nan), t = r;
        Tensor3 rf(2, n, n, 0.0), tf(2, n, n, 0.0);
        rt3::initialize(delta_z, mu, extinction, albedo, phase, r, t);
        rt3_initialize(nstokes,
                       nummu,
                       n,
                       delta_z,
                       mu.data_handle(),
                       extinction,
                       albedo,
                       phase.data_handle(),
                       rf.data_handle(),
                       tf.data_handle());
        tally.add({differ(r, rf), differ(t, tf)});
      }
    }
  }
  tally.report("initialize", "RT3_INITIALIZE");
}

/* rt3::initial_source against RT3_INITIAL_SOURCE, for 1 to 4 Stokes
   parameters, gauss nodes and the extra angle 1, and random source vectors
   and thicknesses. */
void check_initial_source() {
  std::mt19937_64                         gen(2012);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  const auto                              q = polradtran::get_quadrature(4, quadrature_type::gauss);
  Vector                                  mu(5);
  mu[Range{0, 4}] = q.mu;
  mu[4]           = 1.0;
  for (Index nstokes : {1, 2, 3, 4}) {
    const Index nummu = mu.size(), n = nstokes * nummu;
    for (Index rep = 0; rep < 8; rep++) {
      const Numeric delta_z = 1e-3 * u(gen), extinction = 1e-3 * u(gen);
      Tensor3       vec(2, nummu, nstokes);
      for (std::size_t e = 0; e < vec.size(); e++) vec.data_handle()[e] = 2.0 * u(gen) - 1.0;
      Tensor3 src(2, nummu, nstokes, nan), srcf(2, nummu, nstokes, 0.0);
      rt3::initial_source(delta_z, mu, extinction, vec, src);
      rt3_initial_source(
          nstokes, nummu, n, delta_z, mu.data_handle(), extinction, vec.data_handle(), srcf.data_handle());
      tally.add({differ(src, srcf)});
    }
  }
  tally.report("initial_source", "RT3_INITIAL_SOURCE");
}

/* polradtran::nonscatter_layer against RT3_NONSCATTER_LAYER, for 1 to 4 Stokes
   parameters, gauss nodes and the extra angle 1, modes 0 and 1, optical
   depths from 0 to 30, and Planck functions rising, falling, equal and 0. */
void check_nonscatter_layer() {
  routine_tally tally;
  const auto    q = polradtran::get_quadrature(4, quadrature_type::gauss);
  Vector        mu(5);
  mu[Range{0, 4}] = q.mu;
  mu[4]           = 1.0;
  for (Index nstokes : {1, 2, 3, 4}) {
    const Index nummu = mu.size(), n = nstokes * nummu;
    for (Index mode : {0, 1}) {
      for (Numeric deltatau : {0.0, 1e-6, 0.7, 30.0}) {
        for (auto [planck0, planck1] : {std::pair{2.0e-15, 3.1e-15},
                                        std::pair{3.1e-15, 2.0e-15},
                                        std::pair{2.5e-15, 2.5e-15},
                                        std::pair{0.0, 0.0}}) {
          Tensor5 r(2, nummu, nstokes, nummu, nstokes, nan), t = r;
          Tensor3 src(2, nummu, nstokes, nan);
          Tensor3 rf(2, n, n, 0.0), tf(2, n, n, 0.0);
          Matrix  srcf(2, n, 0.0);
          polradtran::nonscatter_layer(mode, deltatau, mu, planck0, planck1, r, t, src);
          rt3_nonscatter_layer(nstokes,
                               nummu,
                               mode,
                               deltatau,
                               mu.data_handle(),
                               planck0,
                               planck1,
                               rf.data_handle(),
                               tf.data_handle(),
                               srcf.data_handle());
          tally.add({differ(r, rf), differ(t, tf), differ(src, srcf)});
        }
      }
    }
  }
  tally.report("nonscatter_layer", "RT3_NONSCATTER_LAYER");
}

/* polradtran::thermal_radiance (SI, ARTS's planck()) against RT3_THERMAL_RADIANCE
   (per micrometre, RT3's PLANCK_FUNCTION on the same exact constants),
   converted, from 1 GHz to 100 THz and 0 to 330 K, for both modes.  RT3
   evaluates exp(x) - 1 where planck() has expm1(x), so it loses digits at
   small x = h nu / k T. */
void check_thermal_radiance() {
  routine_tally tally;
  for (Numeric frequency : {1e9, 89e9, 1e12, 1e14}) {
    const Numeric wavelength       = 1e6 * Constant::c / frequency;
    const Numeric per_um_to_per_hz = wavelength / frequency;
    for (Numeric temperature : {0.0, 2.73, 150.0, 287.5, 330.0}) {
      for (Numeric albedo : {0.0, 0.27}) {
        for (Index mode : {0, 1}) {
          for (auto [nstokes, nummu] : {std::pair{Index{1}, Index{1}}, std::pair{Index{4}, Index{3}}}) {
            Tensor3 rad(2, nummu, nstokes, nan), radf(2, nummu, nstokes, nan);
            polradtran::thermal_radiance(mode, temperature, albedo, frequency, rad);
            rt3_thermal_radiance(nstokes, nummu, mode, temperature, albedo, wavelength, radf.data_handle());
            radf *= per_um_to_per_hz;
            tally.add({differ(rad, radf)});
          }
        }
      }
    }
  }
  tally.report("thermal_radiance", "RT3_THERMAL_RADIANCE", 2e-12);

  bool threw = false;
  try {
    Tensor3 rad(2, 1, 1);
    polradtran::thermal_radiance(0, -1.0, 0.0, 89e9, rad);
  } catch (const std::exception&) { threw = true; }
  if (not threw) throw std::runtime_error("thermal_radiance of a negative temperature did not throw");
}

/* polradtran::lambert_surface_layer against RT3_LAMBERT_SURFACE, for 1 to 4
   Stokes components, mode 0 (the reflecting one) and modes 1 and 2 (no
   reflection), random streams with a zero-weight extra angle, and the
   N = 64 limit. */
void check_lambert_surface() {
  std::mt19937_64                         gen(1991);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (Index nummu : {Index{1}, Index{5}, 64 / nstokes}) {
      for (Index mode : {0, 1, 2}) {
        Vector mu(nummu), w(nummu);
        for (Index j = 0; j < nummu; j++) {
          mu[j] = 0.02 + 0.98 * u(gen);
          w[j]  = nummu > 1 and j == nummu - 1 ? 0.0 : u(gen) / static_cast<Numeric>(nummu);
        }
        const Numeric albedo = u(gen);

        Tensor5 r(2, nummu, nstokes, nummu, nstokes, nan), t(2, nummu, nstokes, nummu, nstokes, nan);
        Tensor3 src(2, nummu, nstokes, nan);
        polradtran::lambert_surface_layer(mode, mu, w, albedo, r, t, src);

        Tensor5 rf(2, nummu, nstokes, nummu, nstokes), tf(2, nummu, nstokes, nummu, nstokes);
        Tensor3 srcf(2, nummu, nstokes);
        rt3_lambert_surface(nstokes,
                            nummu,
                            mode,
                            mu.data_handle(),
                            w.data_handle(),
                            albedo,
                            rf.data_handle(),
                            tf.data_handle(),
                            srcf.data_handle());
        tally.add({differ(r, rf), differ(t, tf), differ(src, srcf)});
      }
    }
  }
  tally.report("lambert_surface_layer", "RT3_LAMBERT_SURFACE");
}

/* polradtran::lambert_radiance (SI, ARTS's planck()) against
   RT3_LAMBERT_RADIANCE (per micrometre, RT3's PLANCK_FUNCTION), converted,
   from 1 GHz to 100 THz and 0 to 330 K, for modes 0 and 1 and every source
   code: they differ as the two Planck functions do (check_thermal_radiance),
   the reflected direct beam to rounding. */
void check_lambert_radiance() {
  std::mt19937_64                         gen(1991);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Numeric frequency : {1e9, 89e9, 1e12, 1e14}) {
    const Numeric wavelength       = 1e6 * Constant::c / frequency;
    const Numeric per_um_to_per_hz = wavelength / frequency;
    for (Numeric ground_temp : {0.0, 2.73, 150.0, 287.5, 330.0}) {
      for (Index mode : {0, 1}) {
        for (Index src_code : {0, 1, 2, 3}) {
          for (auto [nstokes, nummu] : {std::pair{Index{1}, Index{1}}, std::pair{Index{4}, Index{3}}}) {
            const Numeric albedo          = u(gen);
            const Numeric direct_sfc_flux = 1e-4 * u(gen);
            Matrix        rad(nummu, nstokes, nan), radf(nummu, nstokes);
            polradtran::lambert_radiance(mode, src_code, albedo, ground_temp, frequency, direct_sfc_flux, rad);
            rt3_lambert_radiance(nstokes,
                                 nummu,
                                 mode,
                                 src_code,
                                 albedo,
                                 ground_temp,
                                 wavelength,
                                 direct_sfc_flux / per_um_to_per_hz,
                                 radf.data_handle());
            radf *= per_um_to_per_hz;
            tally.add({differ(rad, radf)});
          }
        }
      }
    }
  }
  tally.report("lambert_radiance", "RT3_LAMBERT_RADIANCE", 2e-12);

  // planck() is negative below 0 K, where PLANCK_FUNCTION gave 0
  bool threw = false;
  try {
    Matrix rad(3, 2);
    polradtran::lambert_radiance(0, 2, 0.3, -1.0, 89e9, 0.0, rad);
  } catch (const std::exception&) { threw = true; }
  if (not threw) throw std::runtime_error("lambert_radiance at -1 K did not throw");
}

//! The refractive indices of the Fresnel checks: lossless to water-like
const std::vector<Complex> fresnel_indices{{1.33, 0.0}, {1.5, 0.0}, {3.0, 0.2}, {3.1, 0.4}, {5.5, 2.9}, {7.0, 2.5}};

/* polradtran::fresnel_surface_layer against RT3_FRESNEL_SURFACE, for 1 to 4 Stokes
   components (the [U, V] block with R3 and R4 too), random streams and
   mu = 1.  Not bit-identical: the amplitudes are ARTS's fresnel() at
   acos(mu) in degrees, with the transmitted cosine sqrt(1 - sin^2 / n^2)
   where RT3 has sqrt(n^2 - sin^2), and rtepack::fresnel_reflectance makes
   the Mueller matrix. */
void check_fresnel_surface() {
  std::mt19937_64                         gen(1991);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (const Complex& index : fresnel_indices) {
      const Index nummu = 9;
      Vector      mu(nummu);
      for (Index j = 0; j < nummu - 1; j++) mu[j] = 0.01 + 0.98 * u(gen);
      mu[nummu - 1] = 1.0;

      Tensor5 r(2, nummu, nstokes, nummu, nstokes, nan), t(2, nummu, nstokes, nummu, nstokes, nan);
      Tensor3 src(2, nummu, nstokes, nan);
      polradtran::fresnel_surface_layer(mu, index, r, t, src);

      Tensor5 rf(2, nummu, nstokes, nummu, nstokes), tf(2, nummu, nstokes, nummu, nstokes);
      Tensor3 srcf(2, nummu, nstokes);
      rt3_fresnel_surface(nstokes,
                          nummu,
                          mu.data_handle(),
                          index.real(),
                          index.imag(),
                          rf.data_handle(),
                          tf.data_handle(),
                          srcf.data_handle());
      tally.add({differ(r, rf), differ(t, tf), differ(src, srcf)});
    }
  }
  tally.report("fresnel_surface_layer", "RT3_FRESNEL_SURFACE", 1e-13);
}

/* polradtran::fresnel_radiance (SI, ARTS's planck()) against RT3_FRESNEL_RADIANCE
   (per micrometre, RT3's PLANCK_FUNCTION), converted, for modes 0 and 1:
   they differ as the two Planck functions do (check_thermal_radiance), and
   in the reflection to rounding. */
void check_fresnel_radiance() {
  std::mt19937_64                         gen(1991);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2, 4}) {
    for (const Complex& index : fresnel_indices) {
      for (Numeric frequency : {1e9, 89e9, 3e12}) {
        for (Numeric ground_temp : {0.0, 150.0, 300.0}) {
          for (Index mode : {0, 1}) {
            const Index nummu = 9;
            Vector      mu(nummu);
            for (Index j = 0; j < nummu - 1; j++) mu[j] = 0.01 + 0.98 * u(gen);
            mu[nummu - 1] = 1.0;

            Matrix rad(nummu, nstokes, nan), radf(nummu, nstokes);
            polradtran::fresnel_radiance(mode, mu, index, ground_temp, frequency, rad);

            const Numeric wavelength = 1e6 * Constant::c / frequency;
            rt3_fresnel_radiance(nstokes,
                                 nummu,
                                 mode,
                                 mu.data_handle(),
                                 index.real(),
                                 index.imag(),
                                 ground_temp,
                                 wavelength,
                                 radf.data_handle());
            radf *= wavelength / frequency;
            tally.add({differ(rad, radf)});
          }
        }
      }
    }
  }
  tally.report("fresnel_radiance", "RT3_FRESNEL_RADIANCE", 2e-12);

  // planck() is negative below 0 K, where PLANCK_FUNCTION gave 0
  bool threw = false;
  try {
    Matrix rad(2, 2);
    polradtran::fresnel_radiance(0, Vector{0.3, 0.8}, Complex{1.5, 0.0}, -1.0, 89e9, rad);
  } catch (const std::exception&) { threw = true; }
  if (not threw) throw std::runtime_error("fresnel_radiance at -1 K did not throw");
}

/* rt3::ground_surface against the Fortran grounds it gives RADTRAN, for
   the Lambertian and the Fresnel ground, 1 to 4 Stokes components, azimuth
   orders 0 and 2, and every source code the ground allows.  In each mode
   the surface layer that polradtran::external_surface_layer makes of surf_reflect
   must be the layer of LAMBERT_SURFACE or FRESNEL_SURFACE, and
   gnd_radiance + F direct_reflect (F only with the solar source) the
   radiance of LAMBERT_RADIANCE or FRESNEL_RADIANCE with the direct flux F
   on the ground, converted to SI: to rounding, but for ARTS's planck() in
   place of RT3's Planck function (2e-12, as in check_thermal_radiance). */
void check_ground_surface() {
  std::mt19937_64                         gen(1991);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  const Numeric frequency        = 89e9;
  const Numeric wavelength       = 1e6 * Constant::c / frequency;
  const Numeric per_um_to_per_hz = wavelength / frequency;
  const Numeric ground_temp = 287.5, direct_sfc_flux = 1.7e-4;
  const Index   nummu = 5;

  routine_tally tally;
  for (char type : {'L', 'F'}) {
    for (Index nstokes : {1, 2, 3, 4}) {
      for (Index aziorder : {0, 2}) {
        for (Index src_code : {0, 1, 2, 3}) {
          const bool solar = src_code == 1 or src_code == 3;
          if (type == 'F' and solar) continue;
          Vector mu(nummu), w(nummu);
          for (Index j = 0; j < nummu; j++) {
            mu[j] = 0.02 + 0.98 * u(gen);
            w[j]  = j == nummu - 1 ? 0.0 : u(gen) / static_cast<Numeric>(nummu);
          }
          const Numeric      albedo = u(gen);
          const Complex      index{1.5 + 4.0 * u(gen), 3.0 * u(gen)};
          const rt3::surface ground = type == 'F' ? rt3::surface{polradtran::fresnel_surface{.refractive_index = index}}
                                                  : rt3::surface{polradtran::lambertian_surface{.albedo = albedo}};

          Tensor5 surf_reflect(aziorder + 1, nummu, nstokes, nummu, nstokes, nan);
          Tensor3 gnd_radiance(aziorder + 1, nummu, nstokes, nan), direct_reflect(aziorder + 1, nummu, nstokes, nan);
          rt3::ground_surface(
              ground, src_code, mu, w, frequency, ground_temp, surf_reflect, gnd_radiance, direct_reflect);

          for (Index mode = 0; mode <= aziorder; mode++) {
            Tensor5 r(2, nummu, nstokes, nummu, nstokes, nan), t(2, nummu, nstokes, nummu, nstokes, nan);
            Tensor3 src(2, nummu, nstokes, nan);
            polradtran::external_surface_layer(surf_reflect[mode], r, t, src);
            Matrix rad(nummu, nstokes), direct(nummu, nstokes);
            rad = gnd_radiance[mode];
            if (solar) {
              direct  = direct_reflect[mode];
              direct *= direct_sfc_flux;
              rad    += direct;
            }

            Tensor5 rf(2, nummu, nstokes, nummu, nstokes), tf(2, nummu, nstokes, nummu, nstokes);
            Tensor3 srcf(2, nummu, nstokes);
            Matrix  radf(nummu, nstokes);
            if (type == 'L') {
              rt3_lambert_surface(nstokes,
                                  nummu,
                                  mode,
                                  mu.data_handle(),
                                  w.data_handle(),
                                  albedo,
                                  rf.data_handle(),
                                  tf.data_handle(),
                                  srcf.data_handle());
              rt3_lambert_radiance(nstokes,
                                   nummu,
                                   mode,
                                   src_code,
                                   albedo,
                                   ground_temp,
                                   wavelength,
                                   direct_sfc_flux / per_um_to_per_hz,
                                   radf.data_handle());
            } else {
              rt3_fresnel_surface(nstokes,
                                  nummu,
                                  mu.data_handle(),
                                  index.real(),
                                  index.imag(),
                                  rf.data_handle(),
                                  tf.data_handle(),
                                  srcf.data_handle());
              rt3_fresnel_radiance(nstokes,
                                   nummu,
                                   mode,
                                   mu.data_handle(),
                                   index.real(),
                                   index.imag(),
                                   ground_temp,
                                   wavelength,
                                   radf.data_handle());
            }
            radf *= per_um_to_per_hz;
            tally.add({differ(r, rf), differ(t, tf), differ(src, srcf), differ(rad, radf)});
          }
        }
      }
    }
  }
  tally.report("ground_surface", "LAMBERT_SURFACE, LAMBERT_RADIANCE, FRESNEL_SURFACE and FRESNEL_RADIANCE", 2e-12);

  // RT3 cannot reflect the direct beam from a Fresnel ground
  bool threw = false;
  try {
    Tensor5 surf_reflect(1, 2, 1, 2, 1);
    Tensor3 gnd_radiance(1, 2, 1), direct_reflect(1, 2, 1);
    rt3::ground_surface(polradtran::fresnel_surface{.refractive_index = Complex{1.5, 0.0}},
                        1,
                        Vector{0.3, 0.8},
                        Vector{0.5, 0.5},
                        frequency,
                        ground_temp,
                        surf_reflect,
                        gnd_radiance,
                        direct_reflect);
  } catch (const std::exception&) { threw = true; }
  if (not threw) throw std::runtime_error("ground_surface of a Fresnel ground with a solar source did not throw");
}

/* rt3::get_scat_set against GET_SCAT_SET, with and without delta-M, for
   series shorter and longer than the delta-M order 2 nummu, and random
   coefficients. */
void check_get_scat_set() {
  std::mt19937_64                         gen(1993);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  for (bool delta_m : {false, true}) {
    for (Index nummu : {1, 3, 6, 9}) {
      for (Index nlegin : {Index{0}, Index{2}, 2 * nummu - 1, 2 * nummu, 2 * nummu + 5}) {
        for (Index rep = 0; rep < 4; rep++) {
          Matrix coefin(nlegin + 1, 6);
          for (Index l = 0; l <= nlegin; l++)
            for (Index k = 0; k < 6; k++) coefin[l, k] = static_cast<Numeric>(2 * l + 1) * (u(gen) - 0.3);
          coefin[0, 0]         = 1.0;
          const Numeric extin  = 1e-4 * (0.5 + u(gen));
          const Numeric scatin = extin * u(gen);
          const Index   nrows  = std::max(nlegin + 1, 2 * nummu);
          Matrix        coef(nrows, 6, nan), coeff(nrows, 6, 0.0);
          Vector        out(3, nan), outf(3, 0.0);
          Index         nlegen  = -1;
          std::int64_t  nlegenf = -1;
          rt3::get_scat_set(delta_m, nummu, coefin, extin, scatin, nlegen, coef, out[1], out[2]);
          rt3_get_scat_set(delta_m ? 'Y' : 'N',
                           nummu,
                           nlegin,
                           coefin.data_handle(),
                           extin,
                           scatin,
                           &nlegenf,
                           coeff.data_handle(),
                           &outf[1],
                           &outf[2]);
          out[0]  = static_cast<Numeric>(nlegen);
          outf[0] = static_cast<Numeric>(nlegenf);
          tally.add({differ(coef, coeff), differ(out, outf)});
        }
      }
    }
  }
  tally.report("get_scat_set", "GET_SCAT_SET");

  // A delta-M scaling that divides by zero throws
  bool threw = false;
  try {
    Index   nlegen     = 0;
    Numeric extinction = 0.0, scatter = 0.0;
    Matrix  coef(2, 6);
    rt3::get_scat_set(true, 1, Matrix(1, 6, 1.0), 0.0, 0.0, nlegen, coef, extinction, scatter);
  } catch (const std::exception&) { threw = true; }
  if (not threw) throw std::runtime_error("get_scat_set did not throw for delta-M with zero extinction");
}

enum class layout { mixed, gas, thin, thick, shared };

struct case_spec {
  Index           nstokes, nquad, nuummu;
  quadrature_type quad;
  Index           aziorder, src_code;
  bool            delta_m;
  char            ground;
  Index           nlay;
  layout          lay;
  Numeric         max_delta_tau;
};

std::string describe(const case_spec& c) {
  constexpr const char* layouts[] = {"mixed", "gas", "thin", "thick", "shared"};
  constexpr const char* quads[]   = {"gauss", "double_gauss", "lobatto"};
  return std::format(
      "nstokes {}, nmu {} + {}, {}, aziorder {}, src_code {}, delta_m {}, ground '{}', {} layers {}, "
      "max_delta_tau {:.0e}",
      c.nstokes,
      c.nquad,
      c.nuummu,
      quads[static_cast<int>(c.quad)],
      c.aziorder,
      c.src_code,
      c.delta_m,
      c.ground,
      c.nlay,
      layouts[static_cast<int>(c.lay)],
      c.max_delta_tau);
}

//! RADTRAN's inputs
struct inputs {
  case_spec    spec;
  Index        nummu{};
  Vector       height, temperatures, gas_extinct;
  Vector       scat_extinct, scat_scatter;
  ArrayOfIndex scat_nlegen;
  Tensor3      scat_coef;
  ArrayOfIndex scatlayers, outlevels;
  Vector       extra_mu;
  Numeric      direct_flux{3e-4}, direct_mu{0.6}, ground_temp{287.5}, ground_albedo{0.27}, sky_temp{2.73};
  Numeric      wavelength{3370.0};  // 89 GHz
  Complex      ground_index{3.1, 0.4};
};

struct outputs {
  Vector  mu_values;
  Matrix  up_flux, down_flux;
  Tensor4 up_rad, down_rad;
};

//! [degree + 1, 6] Legendre coefficients (F11, F12, F33, F34, F22, F44) of
//! RT3's three summation cases: Rayleigh (F22 = F11, F44 = F33, no F34),
//! Mie (with F34) and general
Matrix legendre_set(int kind, Numeric g, Index degree) {
  Matrix c(degree + 1, 6, 0.0);
  if (kind == 0) {
    c[0, 0] = 1.0;
    c[0, 1] = -0.5;
    c[0, 4] = 1.0;
    if (degree >= 1) c[1, 2] = c[1, 5] = 1.5;
    if (degree >= 2) c[2, 0] = c[2, 1] = c[2, 4] = 0.5;
    return c;
  }
  for (Index l = 0; l <= degree; l++) {
    const Numeric hg = static_cast<Numeric>(2 * l + 1) * std::pow(g, l);
    c[l, 0]          = hg;
    c[l, 1]          = l > 0 ? -0.1 * hg : 0.0;
    c[l, 2]          = 0.9 * hg;
    c[l, 3]          = l > 0 ? 0.05 * hg : 0.0;
    c[l, 4]          = kind == 1 ? hg : 0.97 * hg;
    c[l, 5]          = kind == 1 ? 0.9 * hg : 0.85 * hg;
  }
  c[0, 0] = 1.0;
  return c;
}

/* rt3::get_scattering against GET_SCATTERING, for 1 to 4 Stokes
   parameters, 1 to 5 angles, azimuth orders 0 to 3, every mode and both of
   two sets.  It only copies and negates, so the two must agree exactly. */
void check_get_scattering() {
  std::mt19937_64                         gen(2008);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (Index nummu : {1, 3, 5}) {
      for (Index aziorder : {0, 1, 3}) {
        constexpr Index nsl = 2;
        Tensor7         scatbuf(nsl, aziorder + 1, 2, nummu, nummu, nstokes, nstokes);
        for (std::size_t e = 0; e < scatbuf.size(); e++) scatbuf.data_handle()[e] = 2.0 * u(gen) - 1.0;
        for (Index set = 0; set < nsl; set++) {
          for (Index mode = 0; mode <= aziorder; mode++) {
            Tensor5 sm(4, nummu, nstokes, nummu, nstokes, nan), smf(4, nummu, nstokes, nummu, nstokes, 0.0);
            rt3::get_scattering(mode, scatbuf[set], sm);
            rt3_get_scattering(nstokes, nummu, mode, aziorder, set + 1, scatbuf.data_handle(), smf.data_handle());
            tally.add({differ(sm, smf)});
          }
        }
      }
    }
  }
  tally.report("get_scattering", "GET_SCATTERING", 0.0);
}

/* rt3::check_norm against RT3_CHECK_NORM: on scattering matrices of
   normalised series both pass (the Fortran would stop the test otherwise),
   for 1 to 4 Stokes parameters, the three quadratures with extra angles
   (weight 0, skipped) and the three kinds of series.  The I-I term scaled
   by 1 + 5e-8 still passes, by 1 + 2e-7 throws (the Fortran's limit is
   1e-7; it is not run there, as it stops), and NaN throws. */
void check_check_norm() {
  std::mt19937_64                         gen(2010);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  Index                                   ncase = 0;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (auto [quad, nquad, nextra] : {std::tuple{quadrature_type::gauss, Index{6}, Index{2}},
                                       std::tuple{quadrature_type::double_gauss, Index{5}, Index{0}},
                                       std::tuple{quadrature_type::lobatto, Index{4}, Index{0}}}) {
      for (int kind : {0, 1, 2}) {
        const Index nummu = nquad + nextra;
        const auto  q     = polradtran::get_quadrature(nquad, quad);
        Vector      mu(nummu), w(nummu, 0.0);
        mu[Range{0, nquad}] = q.mu;
        w[Range{0, nquad}]  = q.weights;
        for (Index i = 0; i < nextra; i++) mu[nquad + i] = 0.97 - 0.5 * static_cast<Numeric>(i);
        const Index  degree = std::min<Index>(rt3::max_legendre_degree(nquad, quad), 7);
        const Matrix coef   = legendre_set(kind, 0.3 + 0.5 * u(gen), kind == 0 ? 2 : degree);

        Tensor6           buf(1, 2, nummu, nummu, nstokes, nstokes);
        rt3::rt3_workdata work;
        rt3::scattering(mu, w, coef, buf, work);
        Tensor5 sm(4, nummu, nstokes, nummu, nstokes);
        rt3::get_scattering(0, buf, sm);

        rt3::check_norm(w, sm);
        rt3_check_norm(nstokes, nummu, w.data_handle(), sm.data_handle());

        const auto scaled = [&](Numeric factor) {
          Tensor5 s2                     = sm;
          s2[joker, joker, 0, joker, 0] *= factor;
          return s2;
        };
        rt3::check_norm(w, scaled(1.0 + 5e-8));
        for (Numeric factor : {1.0 + 2e-7, nan}) {
          bool threw = false;
          try {
            rt3::check_norm(w, scaled(factor));
          } catch (const std::exception&) { threw = true; }
          if (not threw) throw std::runtime_error(std::format("check_norm passed an I-I term scaled by {}", factor));
        }
        ncase++;
      }
    }
  }
  std::cout << std::format(
      "check_norm against RT3_CHECK_NORM: {} normalised cases pass both, 1e-7 limit and NaN throw\n", ncase);
}

/* rt3::scatter_symmetry against SCATTER_SYMMETRY, for 1 to 4 Stokes
   parameters and 1 to 5 angles, with signed zeros among P++ and P+-.  It
   only copies and negates, so the two must agree exactly, signs of zeros
   included. */
void check_scatter_symmetry() {
  std::mt19937_64                         gen(2009);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (Index nummu : {1, 3, 5}) {
      Tensor5 scat(4, nummu, nstokes, nummu, nstokes);
      for (std::size_t e = 0; e < scat.size(); e++) scat.data_handle()[e] = 2.0 * u(gen) - 1.0;
      scat.data_handle()[0] = -0.0;
      scat.data_handle()[1] = 0.0;
      Tensor5 scatf         = scat;
      rt3::scatter_symmetry(scat);
      rt3_scatter_symmetry(nstokes, nummu, scatf.data_handle());
      Index nsign = 0;
      for (std::size_t e = 0; e < scat.size(); e++)
        nsign += std::signbit(scat.data_handle()[e]) != std::signbit(scatf.data_handle()[e]);
      tally.add({differ(scat, scatf), std::pair{nsign, static_cast<Numeric>(nsign)}});
    }
  }
  tally.report("scatter_symmetry", "SCATTER_SYMMETRY", 0.0);
}

/* rt3::get_direct against GET_DIRECT, for 1 to 4 Stokes parameters, 1 to 5
   angles, azimuth orders 0 to 3, every mode and both of two sets.  It only
   copies, so the two must agree exactly. */
void check_get_direct() {
  std::mt19937_64                         gen(2011);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (Index nummu : {1, 3, 5}) {
      for (Index aziorder : {0, 1, 3}) {
        constexpr Index nsl = 2;
        Tensor5         directbuf(nsl, aziorder + 1, 2, nummu, nstokes);
        for (std::size_t e = 0; e < directbuf.size(); e++) directbuf.data_handle()[e] = 2.0 * u(gen) - 1.0;
        for (Index set = 0; set < nsl; set++) {
          for (Index mode = 0; mode <= aziorder; mode++) {
            Tensor3 dv(2, nummu, nstokes, nan), dvf(2, nummu, nstokes, 0.0);
            rt3::get_direct(mode, directbuf[set], dv);
            rt3_get_direct(nstokes, nummu, mode, aziorder, set + 1, directbuf.data_handle(), dvf.data_handle());
            tally.add({differ(dv, dvf)});
          }
        }
      }
    }
  }
  tally.report("get_direct", "GET_DIRECT", 0.0);
}

/* rt3::number_sums against NUMBER_SUMS, for 1 to 4 Stokes parameters and
   Rayleigh, Mie and general series, and Rayleigh series but for one
   element of one row changed by 1 ulp (F34, F22 or F44, first or last
   row), -0.0 against 0.0 and NaN. */
void check_number_sums() {
  std::mt19937_64                         gen(1996);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  std::vector<Matrix>                     sets;
  for (int kind : {0, 1, 2})
    for (Index degree : {Index{0}, Index{2}, Index{9}})
      sets.push_back(legendre_set(kind, 0.3 + 0.5 * u(gen), kind == 0 ? std::min<Index>(degree, 2) : degree));
  for (Index col : {3, 4, 5}) {
    for (bool last : {false, true}) {
      Matrix      c = legendre_set(0, 0.0, 2);
      const Index l = last ? 2 : 0;
      c[l, col]     = std::nextafter(c[l, col], 2.0);
      sets.push_back(c);
    }
  }
  Matrix negzero = legendre_set(0, 0.0, 2);
  negzero[1, 3]  = -0.0;
  sets.push_back(negzero);
  for (Index col : {0, 3, 5}) {
    Matrix c  = legendre_set(1, 0.5, 4);
    c[2, col] = nan;
    sets.push_back(c);
  }

  routine_tally tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (const auto& c : sets) {
      const rt3::IndexVector6 dosum = rt3::number_sums(nstokes, c);
      rt3::IndexVector6       dosumf{-1, -1, -1, -1, -1, -1};
      rt3_number_sums(nstokes, c.nrows() - 1, c.data_handle(), dosumf.data_handle());
      tally.add({differ(dosum, dosumf)});
    }
  }
  tally.report("number_sums", "NUMBER_SUMS", 0.0);
}

/* rt3::sum_legendre against SUM_LEGENDRE, for series of degree 0 to 120
   at scattering angles from -1 to 1, with each DOSUM of NUMBER_SUMS and
   none.  Both phase matrices start with the same distinct values, so an
   element that only one of them writes differs.  Not bit-identical: the
   polynomials are ARTS's (Boost's recurrence), RT3's had two divisions
   per step.  A cosine rounded to just outside [-1, 1] is clamped (RT3
   summed its series there). */
void check_sum_legendre() {
  std::mt19937_64                         gen(1997);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  const std::vector<rt3::IndexVector6>    dosums{{1, 0, 0, 0, 0, 0},
                                                 {1, 1, 1, 0, 0, 0},
                                                 {1, 1, 1, 1, 0, 0},
                                                 {1, 1, 1, 0, 1, 0},
                                                 {1, 1, 1, 1, 1, 1},
                                                 {0, 0, 0, 0, 0, 0}};
  routine_tally                           tally;
  rt3::rt3_workdata                       work;
  for (Index nlegen : {0, 1, 2, 7, 30, 120}) {
    Matrix coef(nlegen + 1, 6);
    for (Index l = 0; l <= nlegen; l++)
      for (Index k = 0; k < 6; k++)
        coef[l, k] = static_cast<Numeric>(2 * l + 1) * std::pow(0.9, static_cast<Numeric>(l)) * (u(gen) - 0.4);
    for (Numeric x : {-1.0, -0.73, 0.0, 0.31, 0.999, 1.0, 2.0 * u(gen) - 1.0}) {
      for (const auto& dosum : dosums) {
        Matrix44 pm, pmf;
        for (Index k = 0; k < 16; k++) pm.data_handle()[k] = pmf.data_handle()[k] = 1000.0 + static_cast<Numeric>(k);
        rt3::sum_legendre(coef, x, dosum, pm, work);
        rt3_sum_legendre(nlegen, coef.data_handle(), x, dosum.data_handle(), pmf.data_handle());
        tally.add({differ(pm, pmf)});
      }
    }
  }
  tally.report("sum_legendre", "SUM_LEGENDRE");

  // One ulp outside [-1, 1] is the value at -1 and 1
  Matrix coef(31, 6);
  for (Index l = 0; l <= 30; l++)
    for (Index k = 0; k < 6; k++) coef[l, k] = static_cast<Numeric>(2 * l + 1) * std::pow(0.9, static_cast<Numeric>(l));
  for (Numeric x : {-1.0, 1.0}) {
    Matrix44 at{}, outside{};
    rt3::sum_legendre(coef, x, {1, 1, 1, 1, 1, 1}, at, work);
    rt3::sum_legendre(coef, std::nextafter(x, 2.0 * x), {1, 1, 1, 1, 1, 1}, outside, work);
    if (differ(at, outside).first != 0)
      throw std::runtime_error(std::format("sum_legendre does not clamp a cosine one ulp outside {}", x));
  }
}

/* rt3::rotate_phase_matrix against ROTATE_PHASE_MATRIX, for 1 to 4 Stokes
   parameters and the direction pairs of SCATTERING: gauss nodes and the
   extra angle 1, both signs of the incoming one, and the azimuths of a
   16-point FFT (forward and backward scattering among them), and random
   pairs.  Both outputs start with the same distinct values. */
void check_rotate_phase_matrix() {
  std::mt19937_64                         gen(1998);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  const auto                              q = polradtran::get_quadrature(4, quadrature_type::gauss);
  std::vector<Numeric>                    mus(q.mu.begin(), q.mu.end());
  mus.push_back(1.0);
  for (Index i = 0; i < 4; i++) mus.push_back(u(gen));

  routine_tally tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (Numeric mu1a : mus) {
      for (Numeric mu2 : mus) {
        for (Numeric sign : {1.0, -1.0}) {
          for (Index k = 1; k <= 9; k++) {
            const Numeric mu1 = sign * mu1a;
            const Numeric delphi =
                k == 9 ? Constant::two_pi * u(gen) : (Constant::two_pi * static_cast<Numeric>(k - 1)) / 16.0;
            const Numeric cos_scat = std::sqrt((1.0 - mu1 * mu1) * (1.0 - mu2 * mu2)) * std::cos(delphi) + mu1 * mu2;
            Matrix44      pm1;
            for (Index e = 0; e < 16; e++) pm1.data_handle()[e] = 2.0 * u(gen) - 1.0;
            Matrix44 pm2, pm2f;
            for (Index e = 0; e < 16; e++)
              pm2.data_handle()[e] = pm2f.data_handle()[e] = 1000.0 + static_cast<Numeric>(e);
            rt3::rotate_phase_matrix(pm1, mu1, mu2, delphi, cos_scat, pm2.view()[Range{0, nstokes}, Range{0, nstokes}]);
            rt3_rotate_phase_matrix(pm1.data_handle(), mu1, mu2, delphi, cos_scat, pm2f.data_handle(), nstokes);
            tally.add({differ(pm2, pm2f)});
          }
        }
      }
    }
  }
  tally.report("rotate_phase_matrix", "ROTATE_PHASE_MATRIX");
}

/* rt3::matrix_symmetry against MATRIX_SYMMETRY, for 1 to 4 Stokes
   parameters, into another matrix and in place (as SCATTERING calls it at
   delphi = pi), with signed zeros among the elements.  The outputs start
   with the same distinct values. */
void check_matrix_symmetry() {
  std::mt19937_64                         gen(1999);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (bool in_place : {false, true}) {
      for (Index rep = 0; rep < 8; rep++) {
        Matrix44 m1;
        for (Index e = 0; e < 16; e++) m1.data_handle()[e] = 2.0 * u(gen) - 1.0;
        m1.data_handle()[rep % 16]       = 0.0;
        m1.data_handle()[(rep + 8) % 16] = -0.0;
        Matrix44 m2 = m1, m2f = m1;
        if (not in_place)
          for (Index e = 0; e < 16; e++) m2.data_handle()[e] = m2f.data_handle()[e] = 1000.0 + static_cast<Numeric>(e);
        const Range stokes{0, nstokes};
        Matrix44    m1f = m1;
        if (in_place) {
          rt3::matrix_symmetry(m2.view()[stokes, stokes], m2.view()[stokes, stokes]);
          rt3_matrix_symmetry(nstokes, m2f.data_handle(), m2f.data_handle());
        } else {
          rt3::matrix_symmetry(m1.view()[stokes, stokes], m2.view()[stokes, stokes]);
          rt3_matrix_symmetry(nstokes, m1f.data_handle(), m2f.data_handle());
        }
        // A zero must keep its sign as the Fortran gives it
        Index nsign = 0;
        for (Index e = 0; e < 16; e++) nsign += std::signbit(m2.data_handle()[e]) != std::signbit(m2f.data_handle()[e]);
        tally.add({differ(m2, m2f), std::pair{nsign, static_cast<Numeric>(nsign)}});
      }
    }
  }
  tally.report("matrix_symmetry", "MATRIX_SYMMETRY", 0.0);
}

/* rt3::makephase against MAKEPHASE, for nmax a power of two from 1 to 512
   (FFT1DR's MN).  Both tables start with the same distinct values, so an
   entry that only one of them writes differs.  Another nmax throws: the
   halves of MAKEPHASE's table then overlap and run past 4 nmax (the
   Fortran's must not be called with one). */
void check_makephase() {
  routine_tally tally;
  for (Index nmax = 1; nmax <= 512; nmax *= 2) {
    Vector phase(4 * nmax), phasef(4 * nmax);
    for (Index e = 0; e < 4 * nmax; e++) phase[e] = phasef[e] = 1000.0 + static_cast<Numeric>(e);
    rt3::makephase(phase);
    rt3_makephase(phasef.data_handle(), nmax);
    tally.add({differ(phase, phasef)});
  }
  tally.report("makephase", "MAKEPHASE");

  for (Index nmax : {3, 5, 12}) {
    bool threw = false;
    try {
      Vector phase(4 * nmax);
      rt3::makephase(phase);
    } catch (const std::exception&) { threw = true; }
    if (not threw) throw std::runtime_error(std::format("makephase with nmax {} did not throw", nmax));
  }
}

/* rt3::fftc against FFTC, for 1 to 256 complex values, with either half of
   MAKEPHASE's table for 512.  Lengths that are no power of two and too
   short a table throw. */
void check_fftc() {
  std::mt19937_64                         gen(2005);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  Vector                                  table(4 * 512);
  rt3::makephase(table);
  routine_tally tally;
  for (Index n = 1; n <= 256; n *= 2) {
    for (Index half : {0, 1}) {
      const auto phase = table[Range{half * 1024, 1024}];
      Vector     data(2 * n);
      for (auto& v : data) v = 2.0 * u(gen) - 1.0;
      Vector dataf = data;
      rt3::fftc(data, phase);
      rt3_fftc(dataf.data_handle(), n, phase.data_handle());
      tally.add({differ(data, dataf)});
    }
  }
  tally.report("fftc", "FFTC");

  for (auto [nvalues, nphase] :
       {std::pair{Index{6}, Index{100}}, std::pair{Index{5}, Index{100}}, std::pair{Index{16}, Index{13}}}) {
    bool threw = false;
    try {
      Vector data(nvalues, 1.0), phase(nphase, 1.0);
      rt3::fftc(data, phase);
    } catch (const std::exception&) { threw = true; }
    if (not threw)
      throw std::runtime_error(std::format("fftc of {} values with {} phases did not throw", nvalues, nphase));
  }
}

/* rt3::fixreal against FIXREAL, both directions, for 1 to 256 complex
   values with the matching half of MAKEPHASE's table for 512. */
void check_fixreal() {
  std::mt19937_64                         gen(2006);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  Vector                                  table(4 * 512);
  rt3::makephase(table);
  routine_tally tally;
  for (Index n = 1; n <= 256; n *= 2) {
    for (auto isign : {rt3::fft_direction::forward, rt3::fft_direction::inverse}) {
      const bool forward = isign == rt3::fft_direction::forward;
      const auto phase   = table[Range{forward ? 0 : 1024, 1024}];
      Vector     data(2 * n);
      for (auto& v : data) v = 2.0 * u(gen) - 1.0;
      Vector  dataf = data;
      Vector2 nyquist{2.0 * u(gen) - 1.0, 0.25};
      Vector2 nyquistf = nyquist;
      rt3::fixreal(data, nyquist, isign, phase);
      rt3_fixreal(dataf.data_handle(), nyquistf.data_handle(), n, forward ? +1 : -1, phase.data_handle());
      tally.add({differ(data, dataf), differ(nyquist, nyquistf)});
    }
  }
  tally.report("fixreal", "FIXREAL");
}

/* rt3::fft1dr against FFT1DR, both directions, for lengths 2 to 512 in an
   order that grows and shrinks, with one fft_workdata throughout (as
   FFT1DR keeps its phase table) and with a fresh one per call.  Lengths
   that are no power of two or above 512 throw and leave the state as it
   was. */
void check_fft1dr() {
  std::mt19937_64                         gen(2002);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  rt3::fft_workdata                       kept;
  for (Index n : {8, 2, 512, 16, 4, 64, 128, 32, 256}) {
    for (auto isign : {rt3::fft_direction::forward, rt3::fft_direction::inverse}) {
      for (bool fresh : {false, true}) {
        Vector x(n);
        for (auto& v : x) v = 2.0 * u(gen) - 1.0;
        Vector            xf = x;
        rt3::fft_workdata own;
        rt3::fft1dr(x, isign, fresh ? own : kept);
        rt3_fft1dr(xf.data_handle(), n, isign == rt3::fft_direction::forward ? +1 : -1);
        tally.add({differ(x, xf)});
      }
    }
  }
  tally.report("fft1dr", "FFT1DR");

  for (Index n : {1, 12, 1024}) {
    bool threw = false;
    try {
      Vector x(n, 1.0);
      rt3::fft1dr(x, rt3::fft_direction::forward, kept);
    } catch (const std::exception&) { threw = true; }
    if (not threw) throw std::runtime_error(std::format("fft1dr of {} values did not throw", n));
  }
  if (kept.mn != 512 or kept.phase.size() != 4 * 512)
    throw std::runtime_error("fft1dr changed its state before throwing");
}

/* fft1dr's format, which an FFTW build must also give: the forward
   transform against the direct sums X_k = sum_j x_j exp(+2 pi i j k / n),
   packed [X_0, X_{n/2}, Re X_1, Im X_1, ...], and the inverse against
   x_j = X_0 + (-1)^j X_{n/2} + 2 sum_k Re(X_k exp(-2 pi i j k / n)), to
   1e-12 of the largest value. */
void check_fft1dr_format() {
  std::mt19937_64                         gen(2003);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  Numeric                                 worst = 0.0;
  rt3::fft_workdata                       fft;
  for (Index n : {2, 4, 8, 64, 512}) {
    const auto angle = [n](Index j, Index k) {
      return Constant::two_pi * static_cast<Numeric>(j * k % n) / static_cast<Numeric>(n);
    };
    Vector x(n), packed(n, 0.0);
    for (auto& v : x) v = 2.0 * u(gen) - 1.0;
    for (Index j = 0; j < n; j++) {
      packed[0] += x[j];
      packed[1] += x[j] * std::cos(angle(j, n / 2));
      for (Index k = 1; k < n / 2; k++) {
        packed[2 * k]     += x[j] * std::cos(angle(j, k));
        packed[2 * k + 1] += x[j] * std::sin(angle(j, k));
      }
    }
    Vector forward = x;
    rt3::fft1dr(forward, rt3::fft_direction::forward, fft);
    worst = std::max(worst, differ(forward, packed).second);

    Vector direct(n, 0.0);
    for (Index j = 0; j < n; j++) {
      direct[j] = packed[0] + (j % 2 == 0 ? 1.0 : -1.0) * packed[1];
      for (Index k = 1; k < n / 2; k++)
        direct[j] += 2.0 * (packed[2 * k] * std::cos(angle(j, k)) + packed[2 * k + 1] * std::sin(angle(j, k)));
    }
    Vector inverse = packed;
    rt3::fft1dr(inverse, rt3::fft_direction::inverse, fft);
    worst = std::max(worst, differ(inverse, direct).second);
  }
  std::cout << std::format("fft1dr against the direct sums of its format: within {:.2e} of the largest value\n", worst);
  if (worst > 1e-12) throw std::runtime_error("fft1dr does not give its documented format");
}

/* rt3::fourier_basis against FOURIER_BASIS, both directions, with and
   without the sines, for orders 0 (any numpts) and above (numpts a power of
   two from 2 to 512), up to and beyond numpts / 2 - 1.  Both vectors are
   compared, as FFT1DR overwrites the real one going to the basis. */
void check_fourier_basis() {
  std::mt19937_64                         gen(2001);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  for (auto direction : {rt3::fourier_direction::to_basis, rt3::fourier_direction::to_real}) {
    for (auto [numpts, order] : {std::pair{Index{1}, Index{0}},
                                 std::pair{Index{6}, Index{0}},
                                 std::pair{Index{30}, Index{0}},
                                 std::pair{Index{2}, Index{1}},
                                 std::pair{Index{8}, Index{1}},
                                 std::pair{Index{16}, Index{3}},
                                 std::pair{Index{16}, Index{12}},
                                 std::pair{Index{64}, Index{12}},
                                 std::pair{Index{512}, Index{40}}}) {
      for (bool even : {false, true}) {
        const Index numbasis = even ? order + 1 : 2 * order + 1;
        Vector      basis(numbasis), real(numpts);
        for (auto& x : basis) x = 2.0 * u(gen) - 1.0;
        for (auto& x : real) x = 2.0 * u(gen) - 1.0;
        Vector            basisf = basis, realf = real;
        rt3::fft_workdata fft;
        rt3::fourier_basis(order, direction, basis, real, fft);
        rt3_fourier_basis(numbasis,
                          order,
                          numpts,
                          direction == rt3::fourier_direction::to_real ? -1 : +1,
                          basisf.data_handle(),
                          realf.data_handle());
        tally.add({differ(basis, basisf), differ(real, realf)});
      }
    }
  }
  tally.report("fourier_basis", "FOURIER_BASIS");
}

/* rt3::fourier_matrix against FOURIER_MATRIX, for 1 to 4 Stokes
   parameters, the NUMPTS of SCATTERING (powers of two up to 512 with
   aziorder > 0, any even number with aziorder 0) and azimuth orders up to
   and beyond numpts / 2 - 1 (the modes FOURIER_BASIS sets to 0).  The
   outputs start with the same distinct values. */
void check_fourier_matrix() {
  std::mt19937_64                         gen(2000);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (auto [numpts, aziorder] : {std::pair{Index{6}, Index{0}},
                                    std::pair{Index{30}, Index{0}},
                                    std::pair{Index{8}, Index{1}},
                                    std::pair{Index{16}, Index{3}},
                                    std::pair{Index{16}, Index{12}},
                                    std::pair{Index{64}, Index{12}},
                                    std::pair{Index{512}, Index{40}}}) {
      const Index numazi = 2 * aziorder + 1;
      const Range stokes{0, nstokes};
      Tensor3     real_matrix(numpts, 4, 4);
      for (std::size_t e = 0; e < real_matrix.size(); e++) real_matrix.data_handle()[e] = 2.0 * u(gen) - 1.0;
      Tensor3 basis(numazi, 4, 4), basisf(numazi, 4, 4);
      for (std::size_t e = 0; e < basis.size(); e++)
        basis.data_handle()[e] = basisf.data_handle()[e] = 1000.0 + static_cast<Numeric>(e);
      rt3::rt3_workdata work;
      rt3::fourier_matrix(real_matrix[joker, stokes, stokes], basis[joker, stokes, stokes], work);
      rt3_fourier_matrix(aziorder, numpts, nstokes, real_matrix.data_handle(), basisf.data_handle());
      tally.add({differ(basis, basisf)});
    }
  }
  tally.report("fourier_matrix", "FOURIER_MATRIX");
}

/* rt3::combine_phase_modes against COMBINE_PHASE_MODES, for 1 to 4 Stokes
   parameters, azimuth orders 0 to 5 and every mode, with signed zeros
   among the Fourier modes.  OUT_MATRIX is a contiguous nstokes x nstokes
   matrix, as is the C++ output here; both start with the same distinct
   values. */
void check_combine_phase_modes() {
  std::mt19937_64                         gen(2004);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (Index aziorder : {0, 1, 5}) {
      const Range stokes{0, nstokes};
      Tensor3     basis(2 * aziorder + 1, 4, 4);
      for (std::size_t e = 0; e < basis.size(); e++) basis.data_handle()[e] = 2.0 * u(gen) - 1.0;
      basis.data_handle()[5] = -0.0;
      basis.data_handle()[6] = 0.0;
      for (Index m = 0; m <= aziorder; m++) {
        const Numeric tmp = u(gen);
        Matrix        out(nstokes, nstokes), outf(nstokes, nstokes);
        for (std::size_t e = 0; e < out.size(); e++)
          out.data_handle()[e] = outf.data_handle()[e] = 1000.0 + static_cast<Numeric>(e);
        rt3::combine_phase_modes(m, tmp, basis[joker, stokes, stokes], out);
        rt3_combine_phase_modes(nstokes, aziorder, m, tmp, basis.data_handle(), outf.data_handle());
        Index nsign = 0;
        for (std::size_t e = 0; e < out.size(); e++)
          nsign += std::signbit(out.data_handle()[e]) != std::signbit(outf.data_handle()[e]);
        tally.add({differ(out, outf), std::pair{nsign, static_cast<Numeric>(nsign)}});
      }
    }
  }
  tally.report("combine_phase_modes", "COMBINE_PHASE_MODES");
}

/* rt3::scattering against SCATTERING, for 1 to 4 Stokes parameters, the
   three quadratures with and without extra angles, azimuth orders 0 to 12,
   and Rayleigh, Mie and general series of several degrees (each summation
   case of NUMBER_SUMS). */
void check_scattering() {
  std::mt19937_64                         gen(1994);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (auto [quad, nquad, nextra] : {std::tuple{quadrature_type::gauss, Index{1}, Index{0}},
                                       std::tuple{quadrature_type::gauss, Index{6}, Index{2}},
                                       std::tuple{quadrature_type::double_gauss, Index{5}, Index{0}},
                                       std::tuple{quadrature_type::lobatto, Index{4}, Index{0}}}) {
      for (Index aziorder : {0, 1, 3, 12}) {
        for (int kind : {0, 1, 2}) {
          for (Index degree : {Index{2}, Index{7}, Index{24}}) {
            const Index nummu = nquad + nextra;
            const auto  q     = polradtran::get_quadrature(nquad, quad);
            Vector      mu(nummu), w(nummu, 0.0);
            mu[Range{0, nquad}] = q.mu;
            w[Range{0, nquad}]  = q.weights;
            for (Index i = 0; i < nextra; i++) mu[nquad + i] = 0.97 - 0.5 * static_cast<Numeric>(i);
            const Matrix coef = legendre_set(kind, 0.3 + 0.5 * u(gen), kind == 0 ? 2 : degree);

            Tensor6           buf(aziorder + 1, 2, nummu, nummu, nstokes, nstokes, nan);
            Tensor6           buff(aziorder + 1, 2, nummu, nummu, nstokes, nstokes, 0.0);
            rt3::rt3_workdata work;
            rt3::scattering(mu, w, coef, buf, work);
            rt3_scattering(nummu,
                           aziorder,
                           nstokes,
                           mu.data_handle(),
                           w.data_handle(),
                           coef.nrows() - 1,
                           coef.data_handle(),
                           1,
                           buff.data_handle());
            tally.add({differ(buf, buff)});
          }
        }
      }
    }
  }
  // A different rounding of COS_SCAT by 1 ulp moves the degree-24 phase
  // function by up to 1.3e-14 of its largest value
  tally.report("scattering", "SCATTERING", 1e-13);
}

/* rt3::direct_scattering against DIRECT_SCATTERING, for 1 to 4 Stokes
   parameters, the three quadratures with and without extra angles (1 among
   them), azimuth orders 0 to 12, Rayleigh, Mie and general series, and the
   sun at mu 0.6, 0.13 and 1 (the zenith). */
void check_direct_scattering() {
  std::mt19937_64                         gen(2007);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (auto [quad, nquad, nextra] : {std::tuple{quadrature_type::gauss, Index{1}, Index{0}},
                                       std::tuple{quadrature_type::gauss, Index{6}, Index{2}},
                                       std::tuple{quadrature_type::double_gauss, Index{5}, Index{0}},
                                       std::tuple{quadrature_type::lobatto, Index{4}, Index{0}}}) {
      for (Index aziorder : {0, 1, 3, 12}) {
        for (int kind : {0, 1, 2}) {
          for (Numeric direct_mu : {0.6, 0.13, 1.0}) {
            const Index nummu = nquad + nextra;
            const auto  q     = polradtran::get_quadrature(nquad, quad);
            Vector      mu(nummu);
            mu[Range{0, nquad}] = q.mu;
            for (Index i = 0; i < nextra; i++) mu[nquad + i] = 1.0 - 0.5 * static_cast<Numeric>(i);
            const Matrix coef = legendre_set(kind, 0.3 + 0.5 * u(gen), kind == 0 ? 2 : 24);

            Tensor4           buf(aziorder + 1, 2, nummu, nstokes, nan);
            Tensor4           buff(aziorder + 1, 2, nummu, nstokes, 0.0);
            rt3::rt3_workdata work;
            rt3::direct_scattering(mu, coef, direct_mu, buf, work);
            rt3_direct_scattering(nummu,
                                  aziorder,
                                  nstokes,
                                  mu.data_handle(),
                                  coef.nrows() - 1,
                                  coef.data_handle(),
                                  direct_mu,
                                  1,
                                  buff.data_handle());
            tally.add({differ(buf, buff)});
          }
        }
      }
    }
  }
  // As for scattering: the Legendre sums amplify the rounding of COS_SCAT
  tally.report("direct_scattering", "DIRECT_SCATTERING", 1e-13);
}

//! The condition number (infinity norm) of 1 - R R, the larger of plus and
//! minus, which the doubling inverts: it amplifies the inverse's rounding
Numeric inverse_condition(ConstTensor3View reflect) {
  const Index n     = reflect.extent(1);
  Numeric     kappa = 1.0;
  for (Index h = 0; h < 2; h++) {
    Matrix a(n, n), ainv(n, n);
    identity(a);
    mult(a, reflect[h], reflect[1 - h], -1.0, 1.0);
    inv(ainv, a);
    kappa = std::max(kappa, norm_inf(a) * norm_inf(ainv));
  }
  return kappa;
}

/* polradtran::doubling_integration against RT3_DOUBLING_INTEGRATION, on the thin
   starting layer of a scattering set (rt3::scattering, get_scattering,
   initialize) for 1 to 4 Stokes parameters (symmetric for 1 and 2), modes
   0 and 1, every source code and 0 to 20 doublings, with random source
   vectors.  All outputs are compared, and the overwritten inputs.  LAPACK's
   inverse replaces LINPACK's and DGEMM's alpha and beta absorb the matrix
   additions, so the two agree to rounding amplified by the doublings: each
   about squares T, doubling its relative error, and inverts 1 - R R,
   whose condition number kappa amplifies the inverse's rounding.  n
   doublings amplify rounding by 2^n kappa.  With 0 and 1 doublings the
   starting layer is optically thick (tau 4 and 2), R is near 1 and kappa
   up to 1e2; with 20 it is thin and kappa 1.  The C++ and the Fortran
   differ by up to 0.8 epsilon times 2^n kappa: 7e-15 for 1 doubling (kappa
   99) and 1.8e-10 for 20 (AMD x86_64, MKL). */
void check_doubling_integration() {
  std::mt19937_64                         gen(2014);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  routine_tally                           tally;
  const auto                              q = polradtran::get_quadrature(4, quadrature_type::gauss);
  Vector                                  mu(5), w(5, 0.0);
  mu[Range{0, 4}] = q.mu;
  w[Range{0, 4}]  = q.weights;
  mu[4]           = 1.0;
  for (Index nstokes : {1, 2, 3, 4}) {
    const Index nummu = mu.size(), n = nstokes * nummu;
    const bool  symmetric = nstokes <= 2;
    for (Index mode : {0, 1}) {
      for (Index num_doubles : {0, 1, 5, 20}) {
        for (Index src_code : {0, 1, 2, 3}) {
          const Matrix      coef = legendre_set(2, 0.3 + 0.5 * u(gen), 7);
          Tensor6           buf(2, 2, nummu, nummu, nstokes, nstokes);
          rt3::rt3_workdata work(nstokes, nummu, 0, 0, 0, 0);
          rt3::scattering(mu, w, coef, buf, work);
          Tensor5 sm(4, nummu, nstokes, nummu, nstokes);
          rt3::get_scattering(mode, buf, sm);

          const Numeric extinction = 2e-3, albedo = 0.3 + 0.65 * u(gen);
          const Numeric delta_z = 2000.0 / std::pow(2.0, static_cast<Numeric>(num_doubles));
          Tensor3       reflect(2, n, n), trans(2, n, n);
          rt3::initialize(delta_z,
                          mu,
                          extinction,
                          albedo,
                          sm,
                          reflect.view_as(2, nummu, nstokes, nummu, nstokes),
                          trans.view_as(2, nummu, nstokes, nummu, nstokes));
          Matrix exp_source(2, n), lin_source(2, n);
          for (std::size_t e = 0; e < exp_source.size(); e++) {
            exp_source.data_handle()[e] = 1e-15 * u(gen);
            lin_source.data_handle()[e] = 1e-15 * u(gen);
          }
          const Numeric expfactor = std::exp(-extinction * delta_z / 0.6), linfactor = 0.3 * delta_z / 2000.0;

          const Numeric amp = num_doubles == 0
                                  ? 1.0
                                  : std::pow(2.0, static_cast<Numeric>(num_doubles)) * inverse_condition(reflect);

          Tensor3 reflectf = reflect, transf = trans;
          Matrix  exp_sourcef = exp_source, lin_sourcef = lin_source;
          Tensor3 tr(2, n, n, nan), tt(2, n, n, nan), trf(2, n, n, 0.0), ttf(2, n, n, 0.0);
          Matrix  ts(2, n, nan), tsf(2, n, 0.0);
          polradtran::doubling_integration(num_doubles,
                                           src_code,
                                           symmetric,
                                           reflect,
                                           trans,
                                           exp_source,
                                           expfactor,
                                           lin_source,
                                           linfactor,
                                           tr,
                                           tt,
                                           ts,
                                           work);
          rt3_doubling_integration(n,
                                   num_doubles,
                                   src_code,
                                   symmetric,
                                   reflectf.data_handle(),
                                   transf.data_handle(),
                                   exp_sourcef.data_handle(),
                                   expfactor,
                                   lin_sourcef.data_handle(),
                                   linfactor,
                                   trf.data_handle(),
                                   ttf.data_handle(),
                                   tsf.data_handle());
          tally.add({{differ(tr, trf), amp},
                     {differ(tt, ttf), amp},
                     {differ(ts, tsf), amp},
                     {differ(reflect, reflectf), amp},
                     {differ(trans, transf), amp},
                     {differ(exp_source, exp_sourcef), amp},
                     {differ(lin_source, lin_sourcef), amp}});
        }
      }
    }
  }
  tally.report("doubling_integration", "RT3_DOUBLING_INTEGRATION", rounding);
}

//! A layer as RADTRAN combines them: [2, n, n] reflection and transmission
//! and [2, n] source
struct slab {
  Tensor3 reflect, trans;
  Matrix  source;
};

/* A scattering layer of optical depth tau for the adding checks, made as
   RADTRAN makes it with the ported routines: SCATTERING, GET_SCATTERING and
   INITIALIZE on a sublayer of optical depth at most 1e-6, doubled up to
   tau, with random thermal and solar sources. */
slab scattering_slab(
    std::mt19937_64& gen, ConstVectorView mu, ConstVectorView w, Index nstokes, Index mode, Numeric tau) {
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);
  const Index                             nummu = mu.size(), n = nstokes * nummu;
  const Matrix                            coef = legendre_set(2, 0.3 + 0.5 * u(gen), 7);
  Tensor6                                 buf(2, 2, nummu, nummu, nstokes, nstokes);
  rt3::rt3_workdata                       work(nstokes, nummu, 0, 0, 0, 0);
  rt3::scattering(mu, w, coef, buf, work);
  Tensor5 sm(4, nummu, nstokes, nummu, nstokes);
  rt3::get_scattering(mode, buf, sm);

  const Index   num_doubles = static_cast<Index>(std::ceil(std::log2(tau / 1e-6)));
  const Numeric extinction = 1e-3, zdiff = tau / extinction;
  const Numeric delta_z = zdiff / std::pow(2.0, static_cast<Numeric>(num_doubles));
  Tensor3       reflect(2, n, n), trans(2, n, n);
  rt3::initialize(delta_z,
                  mu,
                  extinction,
                  0.3 + 0.65 * u(gen),
                  sm,
                  reflect.view_as(2, nummu, nstokes, nummu, nstokes),
                  trans.view_as(2, nummu, nstokes, nummu, nstokes));
  Matrix exp_source(2, n), lin_source(2, n);
  for (std::size_t e = 0; e < exp_source.size(); e++) {
    exp_source.data_handle()[e] = 1e-15 * u(gen);
    lin_source.data_handle()[e] = 1e-15 * u(gen);
  }

  slab out{.reflect = Tensor3(2, n, n), .trans = Tensor3(2, n, n), .source = Matrix(2, n)};
  polradtran::doubling_integration(num_doubles,
                                   3,
                                   nstokes <= 2,
                                   reflect,
                                   trans,
                                   exp_source,
                                   std::exp(-extinction * delta_z / 0.6),
                                   lin_source,
                                   0.3 / std::pow(2.0, static_cast<Numeric>(num_doubles)),
                                   out.reflect,
                                   out.trans,
                                   out.source,
                                   work);
  return out;
}

/* polradtran::combine_layers against RT3_COMBINE_LAYERS: pairs of the layers
   RADTRAN combines, for 1 to 4 Stokes components, modes 0 and 1 and two
   numbers of streams (with an extra angle): thin (tau 1e-3) and thick
   (tau 4) scattering layers, a gas layer, and the Lambertian and Fresnel
   ground layers below.  Not as is, like doubling_integration: DGEMM's and
   DGEMV's alpha and beta absorb the MIDENTITY, MSUB and MADD, and the
   inverse is LAPACK's. */
void check_combine_layers() {
  std::mt19937_64                         gen(2015);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nquad : {2, 8}) {
    const auto q = polradtran::get_quadrature(nquad, quadrature_type::gauss);
    Vector     mu(nquad + 1), w(nquad + 1, 0.0);
    mu[Range{0, nquad}] = q.mu;
    w[Range{0, nquad}]  = q.weights;
    mu[nquad]           = 1.0;
    const Index nummu   = mu.size();
    for (Index nstokes : {1, 2, 3, 4}) {
      const Index n = nstokes * nummu;
      for (Index mode : {0, 1}) {
        // The thin and thick scattering layers, the gas layer, and the
        // Lambertian and Fresnel grounds
        const auto layer = [&](int kind) {
          if (kind < 2) return scattering_slab(gen, mu, w, nstokes, mode, kind == 0 ? 1e-3 : 4.0);
          Tensor5 r(2, nummu, nstokes, nummu, nstokes), t(2, nummu, nstokes, nummu, nstokes);
          Tensor3 src(2, nummu, nstokes);
          if (kind == 2) polradtran::nonscatter_layer(mode, 0.7, mu, 3e-15, 4e-15, r, t, src);
          if (kind == 3) polradtran::lambert_surface_layer(mode, mu, w, 0.3, r, t, src);
          if (kind == 4) polradtran::fresnel_surface_layer(mu, Complex{3.1, 0.4}, r, t, src);
          return slab{.reflect = Tensor3{r.view_as(2, n, n)},
                      .trans   = Tensor3{t.view_as(2, n, n)},
                      .source  = Matrix{src.view_as(2, n)}};
        };
        for (int top_kind : {0, 1, 2}) {
          for (int bottom_kind : {0, 1, 2, 3, 4}) {
            const slab top = layer(top_kind), bottom = layer(bottom_kind);

            Tensor3           r(2, n, n, nan), t(2, n, n, nan), rf(2, n, n), tf(2, n, n);
            Matrix            src(2, n, nan), srcf(2, n);
            rt3::rt3_workdata work(nstokes, nummu, 0, 0, 0, 0);
            polradtran::combine_layers(
                top.reflect, top.trans, top.source, bottom.reflect, bottom.trans, bottom.source, r, t, src, work);

            // The Fortran declares no intent: give it copies
            slab a = top, b = bottom;
            rt3_combine_layers(n,
                               a.reflect.data_handle(),
                               a.trans.data_handle(),
                               a.source.data_handle(),
                               b.reflect.data_handle(),
                               b.trans.data_handle(),
                               b.source.data_handle(),
                               rf.data_handle(),
                               tf.data_handle(),
                               srcf.data_handle());
            tally.add({differ(r, rf), differ(t, tf), differ(src, srcf)});
          }
        }
      }
    }
  }
  tally.report("combine_layers", "RT3_COMBINE_LAYERS", 1e-13);
}

/* polradtran::internal_radiance against RT3_INTERNAL_RADIANCE, for 1 to 4 Stokes
   components, modes 0 and 1 and two numbers of streams (with an extra
   angle): above the level a thin or thick scattering layer, or a gas
   layer, below it a scattering layer on a Lambertian or Fresnel ground
   (combined by polradtran::combine_layers), with random incident radiances at the
   top and bottom.  Not as is: BLAS's alpha and beta absorb the MIDENTITY,
   MSUB and MADD, and the inverse is LAPACK's. */
void check_internal_radiance() {
  std::mt19937_64                         gen(2016);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nquad : {2, 8}) {
    const auto q = polradtran::get_quadrature(nquad, quadrature_type::gauss);
    Vector     mu(nquad + 1), w(nquad + 1, 0.0);
    mu[Range{0, nquad}] = q.mu;
    w[Range{0, nquad}]  = q.weights;
    mu[nquad]           = 1.0;
    const Index nummu   = mu.size();
    for (Index nstokes : {1, 2, 3, 4}) {
      const Index n = nstokes * nummu;
      for (Index mode : {0, 1}) {
        for (int up_kind : {0, 1, 2}) {
          for (bool fresnel : {false, true}) {
            slab up;
            if (up_kind < 2) {
              up = scattering_slab(gen, mu, w, nstokes, mode, up_kind == 0 ? 1e-3 : 4.0);
            } else {
              Tensor5 r(2, nummu, nstokes, nummu, nstokes), t(2, nummu, nstokes, nummu, nstokes);
              Tensor3 src(2, nummu, nstokes);
              polradtran::nonscatter_layer(mode, 0.7, mu, 3e-15, 4e-15, r, t, src);
              up = {.reflect = Tensor3{r.view_as(2, n, n)},
                    .trans   = Tensor3{t.view_as(2, n, n)},
                    .source  = Matrix{src.view_as(2, n)}};
            }

            const slab layer = scattering_slab(gen, mu, w, nstokes, mode, 1.5);
            Tensor5    r(2, nummu, nstokes, nummu, nstokes), t(2, nummu, nstokes, nummu, nstokes);
            Tensor3    src(2, nummu, nstokes);
            if (fresnel)
              polradtran::fresnel_surface_layer(mu, Complex{3.1, 0.4}, r, t, src);
            else
              polradtran::lambert_surface_layer(mode, mu, w, 0.3, r, t, src);
            rt3::rt3_workdata work(nstokes, nummu, 0, 0, 0, 0);
            slab              down{.reflect = Tensor3(2, n, n), .trans = Tensor3(2, n, n), .source = Matrix(2, n)};
            polradtran::combine_layers(layer.reflect,
                                       layer.trans,
                                       layer.source,
                                       r.view_as(2, n, n),
                                       t.view_as(2, n, n),
                                       src.view_as(2, n),
                                       down.reflect,
                                       down.trans,
                                       down.source,
                                       work);

            Vector intoprad(n), inbottomrad(n);
            for (Index k = 0; k < n; k++) {
              intoprad[k]    = 1e-15 * u(gen);
              inbottomrad[k] = 1e-15 * u(gen);
            }

            Vector uprad(n, nan), downrad(n, nan), upradf(n), downradf(n);
            polradtran::internal_radiance(up.reflect,
                                          up.trans,
                                          up.source,
                                          down.reflect,
                                          down.trans,
                                          down.source,
                                          intoprad,
                                          inbottomrad,
                                          uprad,
                                          downrad,
                                          work);

            // The Fortran declares no intent: give it copies
            slab   a = up, b = down;
            Vector top = intoprad, bottom = inbottomrad;
            rt3_internal_radiance(n,
                                  a.reflect.data_handle(),
                                  a.trans.data_handle(),
                                  a.source.data_handle(),
                                  b.reflect.data_handle(),
                                  b.trans.data_handle(),
                                  b.source.data_handle(),
                                  top.data_handle(),
                                  bottom.data_handle(),
                                  upradf.data_handle(),
                                  downradf.data_handle());
            tally.add({differ(uprad, upradf), differ(downrad, downradf)});
          }
        }
      }
    }
  }
  tally.report("internal_radiance", "RT3_INTERNAL_RADIANCE", 1e-13);
}

inputs make_inputs(const case_spec& c, std::mt19937_64& gen) {
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  inputs in;
  in.spec         = c;
  in.nummu        = c.nquad + c.nuummu;
  in.height       = Vector(c.nlay + 1);
  in.temperatures = Vector(c.nlay + 1);
  in.gas_extinct  = Vector(c.nlay);
  for (Index l = 0; l <= c.nlay; l++) {
    in.height[l]       = 750.0 * static_cast<Numeric>(c.nlay - l);
    in.temperatures[l] = 205.0 + 86.0 * static_cast<Numeric>(l) / static_cast<Numeric>(c.nlay) + 6.0 * (u(gen) - 0.5);
  }
  for (Index l = 0; l < c.nlay; l++) in.gas_extinct[l] = c.lay == layout::thin ? 1e-12 + 1e-9 * u(gen) : 4e-4 * u(gen);
  if (c.nlay > 2) in.gas_extinct[2] = -1e-5;

  const Index nleglim = rt3::max_legendre_degree(c.nquad, c.quad);
  // Delta-M needs a series beyond 2 nummu to do anything
  const Index degree = c.delta_m ? 2 * in.nummu + 3 : std::min<Index>(6, nleglim);

  std::vector<Matrix> sets;
  in.scatlayers = ArrayOfIndex(c.nlay, 0);
  for (Index l = 0; l < c.nlay; l++) {
    if (c.lay == layout::gas or (c.lay == layout::mixed and l % 3 == 1)) continue;
    if (c.lay == layout::shared and sets.size() >= 2) {
      in.scatlayers[l] = 1 + l % 2;
      continue;
    }
    int kind = c.delta_m ? 2 : static_cast<int>(sets.size() % 3);
    if (kind == 0 and nleglim < 2) kind = 2;
    sets.push_back(legendre_set(kind, 0.3 + 0.5 * u(gen), kind == 0 ? 2 : degree));
    in.scatlayers[l] = static_cast<Index>(sets.size());
  }

  const Index nsl    = static_cast<Index>(sets.size());
  Index       ldcoef = 1;
  for (const auto& s : sets) ldcoef = std::max(ldcoef, s.nrows());
  in.scat_extinct     = Vector(nsl);
  in.scat_scatter     = Vector(nsl);
  in.scat_nlegen      = ArrayOfIndex(nsl);
  in.scat_coef        = Tensor3(nsl, ldcoef, 6, 0.0);
  const Numeric scale = c.lay == layout::thin ? 1e-11 : c.lay == layout::thick ? 5e-3 : 2e-4;
  for (Index s = 0; s < nsl; s++) {
    in.scat_extinct[s]                         = scale * (0.5 + u(gen));
    in.scat_scatter[s]                         = in.scat_extinct[s] * (0.3 + 0.65 * u(gen));
    in.scat_nlegen[s]                          = sets[s].nrows() - 1;
    in.scat_coef[s, Range{0, sets[s].nrows()}] = sets[s];
  }

  // Every level, the bottom one first
  in.outlevels = ArrayOfIndex(c.nlay + 1);
  for (Index l = 0; l <= c.nlay; l++) in.outlevels[l] = c.nlay + 1 - l;

  in.extra_mu = Vector(c.nuummu);
  for (Index i = 0; i < c.nuummu; i++) in.extra_mu[i] = 1.0 - 0.6 * static_cast<Numeric>(i);
  return in;
}

//! The C++ RADTRAN, with its ground made by rt3::ground_surface; all of
//! mu_values is output.  It works in SI at the frequency, the Fortran per
//! micrometre at the wavelength: its inputs and outputs are converted
//! (B_nu = B_lambda[um^-1] lambda[um] / f)
//! 2^n for the most doublings n of a scattering layer, as RADTRAN counts
//! them (delta-M, which only lowers the extinction, aside): the
//! amplification of rounding (see check_doubling_integration)
Numeric doubling_amplification(const inputs& in) {
  Numeric amp = 1.0;
  for (std::size_t l = 0; l < in.scatlayers.size(); l++) {
    const Index set = in.scatlayers[l];
    if (set < 1 or in.scat_scatter[set - 1] == 0.0) continue;
    const Numeric extinct = in.scat_extinct[set - 1] + std::max(in.gas_extinct[l], 0.0);
    const Numeric f =
        std::log2(std::max(extinct * std::abs(in.height[l] - in.height[l + 1]), 1e-7) / in.spec.max_delta_tau);
    if (f > 0.0) amp = std::max(amp, std::pow(2.0, std::floor(f) + 1.0));
  }
  return amp;
}

outputs run_cpp(const inputs& in, rt3::rt3_workdata& work) {
  const Numeric frequency        = 1e6 * Constant::c / in.wavelength;
  const Numeric per_um_to_per_hz = in.wavelength / frequency;
  const Index   no = static_cast<Index>(in.outlevels.size()), ns = in.spec.nstokes;
  outputs       out{.mu_values = Vector(in.nummu, nan),
                    .up_flux   = Matrix(no, ns, nan),
                    .down_flux = Matrix(no, ns, nan),
                    .up_rad    = Tensor4(no, in.spec.aziorder + 1, in.nummu, ns, nan),
                    .down_rad  = Tensor4(no, in.spec.aziorder + 1, in.nummu, ns, nan)};

  // The ground, made by rt3::ground_surface on the streams RADTRAN makes
  const Index nquad = in.spec.nquad, nmode = in.spec.aziorder + 1;
  const auto  q = polradtran::get_quadrature(nquad, in.spec.quad);
  Vector      mu(in.nummu), w(in.nummu, 0.0);
  mu[Range{0, nquad}]              = q.mu;
  mu[Range{nquad, in.spec.nuummu}] = in.extra_mu;
  w[Range{0, nquad}]               = q.weights;
  const rt3::surface ground = in.spec.ground == 'F'
                                  ? rt3::surface{polradtran::fresnel_surface{.refractive_index = in.ground_index}}
                                  : rt3::surface{polradtran::lambertian_surface{.albedo = in.ground_albedo}};
  Tensor5            surf_reflect(nmode, in.nummu, ns, in.nummu, ns);
  Tensor3            gnd_radiance(nmode, in.nummu, ns), direct_reflect(nmode, in.nummu, ns);
  rt3::ground_surface(
      ground, in.spec.src_code, mu, w, frequency, in.ground_temp, surf_reflect, gnd_radiance, direct_reflect);

  rt3::radtran(in.spec.max_delta_tau,
               in.spec.src_code,
               in.spec.quad,
               in.spec.delta_m,
               in.direct_flux * per_um_to_per_hz,
               in.direct_mu,
               surf_reflect,
               gnd_radiance,
               direct_reflect,
               in.sky_temp,
               frequency,
               in.height,
               in.temperatures,
               in.gas_extinct,
               in.scat_extinct,
               in.scat_scatter,
               in.scat_nlegen,
               in.scat_coef,
               in.scatlayers,
               in.outlevels,
               in.extra_mu,
               out.mu_values,
               out.up_flux,
               out.down_flux,
               out.up_rad,
               out.down_rad,
               work);
  out.up_flux   /= per_um_to_per_hz;
  out.down_flux /= per_um_to_per_hz;
  out.up_rad    /= per_um_to_per_hz;
  out.down_rad  /= per_um_to_per_hz;
  return out;
}

/* A work data sized for RADTRAN's inputs, as rt3::radtran sizes it, with
   every array NaN, so that a read before a write shows in the results.  The
   scratch of the scattering routines, which they size themselves, is NaN
   at its largest size.  FFT1DR's table is kept between calls, so it stays
   empty. */
rt3::rt3_workdata poisoned_workdata(const inputs& in) {
  Index legendre_rows = 2 * in.nummu;
  for (Index l : in.scat_nlegen) legendre_rows = std::max(legendre_rows, l + 1);
  rt3::rt3_workdata w(
      in.spec.nstokes, in.nummu, in.spec.aziorder, in.spec.nlay, in.scat_extinct.extent(0), legendre_rows);
  w.scat_matrix.resize(1025, 4, 4);
  w.basis_matrix.resize(1025, 4, 4);
  w.legendre_p.resize(1024);
  w.real_vector.resize(1024);
  w.basis_vector.resize(1025);
  for (Vector* v : {&w.quad_weights,
                    &w.set_extinct,
                    &w.set_scatter,
                    &w.extinctions,
                    &w.albedos,
                    &w.direct_level_flux,
                    &w.legendre_p,
                    &w.real_vector,
                    &w.basis_vector,
                    &w.xv,
                    &w.yv})
    *v = nan;
  for (Matrix* m : {&w.legendre_coef,
                    &w.direct_vector,
                    &w.thermal_vector,
                    &w.exp_source,
                    &w.lin_source,
                    &w.source1,
                    &w.upsource,
                    &w.downsource,
                    &w.sky_radiance,
                    &w.ground_radiance,
                    &w.direct_radiance,
                    &w.t_lin,
                    &w.cnst,
                    &w.t_const,
                    &w.x,
                    &w.y,
                    &w.gamma})
    *m = nan;
  for (Tensor3* t : {&w.scat_matrix,
                     &w.basis_matrix,
                     &w.source,
                     &w.reflect1,
                     &w.upreflect,
                     &w.downreflect,
                     &w.trans1,
                     &w.uptrans,
                     &w.downtrans})
    *t = nan;
  w.reflect        = nan;
  w.trans          = nan;
  w.scatter_matrix = nan;
  w.directbuf      = nan;
  w.scatbuf        = nan;
  for (Index& s : w.scat_nums) s = -1;
  return w;
}

//! The Fortran RADTRAN, which marks the extra angles by non-zero MU_VALUES
//! (QUAD_TYPE 'E') and declares no intent: it gets copies
outputs run_fortran(inputs in) {
  const Index no = static_cast<Index>(in.outlevels.size()), ns = in.spec.nstokes;
  outputs     out{.mu_values = Vector(in.nummu, 0.0),
                  .up_flux   = Matrix(no, ns, 0.0),
                  .down_flux = Matrix(no, ns, 0.0),
                  .up_rad    = Tensor4(no, in.spec.aziorder + 1, in.nummu, ns, 0.0),
                  .down_rad  = Tensor4(no, in.spec.aziorder + 1, in.nummu, ns, 0.0)};
  out.mu_values[Range{in.spec.nquad, in.spec.nuummu}] = in.extra_mu;

  char quad = 'G';
  if (in.spec.quad == quadrature_type::double_gauss) quad = 'D';
  if (in.spec.quad == quadrature_type::lobatto) quad = 'L';
  if (in.spec.nuummu > 0) quad = 'E';

  const Index nsl = in.scat_extinct.extent(0);
  // At least one entry, so that no array is empty
  Vector extinct(std::max<Index>(nsl, 1), 0.0), scatter(std::max<Index>(nsl, 1), 0.0);
  extinct[Range{0, nsl}] = in.scat_extinct;
  scatter[Range{0, nsl}] = in.scat_scatter;
  Tensor3 coef(std::max<Index>(nsl, 1), in.scat_coef.extent(1), 6, 0.0);
  coef[Range{0, nsl}] = in.scat_coef;
  ArrayOfIndex nlegen = in.scat_nlegen;
  nlegen.resize(std::max<Index>(nsl, 1));

  rt3_radtran(ns,
              in.nummu,
              in.spec.aziorder,
              in.spec.max_delta_tau,
              in.spec.src_code,
              quad,
              in.spec.delta_m ? 'Y' : 'N',
              in.direct_flux,
              in.direct_mu,
              in.ground_temp,
              in.spec.ground,
              in.ground_albedo,
              in.ground_index.real(),
              in.ground_index.imag(),
              in.sky_temp,
              in.wavelength,
              in.spec.nlay,
              in.height.data_handle(),
              in.temperatures.data_handle(),
              in.gas_extinct.data_handle(),
              nsl,
              extinct.data_handle(),
              scatter.data_handle(),
              nlegen.data(),
              in.scat_coef.extent(1),
              coef.data_handle(),
              in.scatlayers.data(),
              no,
              in.outlevels.data(),
              out.mu_values.data_handle(),
              out.up_flux.data_handle(),
              out.down_flux.data_handle(),
              out.up_rad.data_handle(),
              out.down_rad.data_handle());
  return out;
}

std::vector<case_spec> cases() {
  constexpr auto G = quadrature_type::gauss, D = quadrature_type::double_gauss, L = quadrature_type::lobatto;

  std::vector<case_spec> out;
  // Stokes, quadrature, extra angles, azimuth order, sources and ground
  for (Index ns : {1, 2, 3, 4})
    for (auto [quad, extra] :
         {std::pair{G, Index{0}}, std::pair{G, Index{2}}, std::pair{D, Index{0}}, std::pair{L, Index{0}}})
      for (Index aziorder : {0, 3})
        for (auto [src, ground] : {std::pair{Index{2}, 'L'},
                                   std::pair{Index{3}, 'L'},
                                   std::pair{Index{1}, 'L'},
                                   std::pair{Index{0}, 'L'},
                                   std::pair{Index{2}, 'F'}})
          out.push_back({ns, 6, extra, quad, aziorder, src, false, ground, 5, layout::mixed, 1e-6});
  // Delta-M, with the quadratures whose NLEGLIM keeps the scaled series
  for (auto quad : {G, L})
    for (Index src : {2, 3}) out.push_back({4, 6, 0, quad, 2, src, true, 'L', 5, layout::mixed, 1e-6});
  // Layouts and doubling depths
  for (layout lay : {layout::gas, layout::thin, layout::thick, layout::shared})
    for (Index src : {2, 3}) out.push_back({3, 8, 1, G, 5, src, false, 'L', 7, lay, 1e-6});
  for (Numeric mdt : {1e-3, 1e-8}) out.push_back({2, 8, 0, D, 2, 3, false, 'L', 4, layout::mixed, mdt});
  // A high azimuth order, the N = 64 limit, many layers, one layer
  out.push_back({4, 4, 0, G, 12, 3, false, 'L', 6, layout::mixed, 1e-6});
  out.push_back({4, 16, 0, D, 1, 3, false, 'L', 3, layout::mixed, 1e-6});
  out.push_back({2, 8, 1, G, 4, 3, false, 'L', 40, layout::mixed, 1e-6});
  out.push_back({1, 1, 0, G, 0, 2, false, 'L', 1, layout::thick, 1e-6});
  return out;
}
}  // namespace

int main() try {
  check_quadratures();
  check_thermal_radiance();
  check_lambert_surface();
  check_lambert_radiance();
  check_fresnel_surface();
  check_fresnel_radiance();
  check_ground_surface();
  check_nonscatter_layer();
  check_initial_source();
  check_initialize();
  check_doubling_integration();
  check_combine_layers();
  check_internal_radiance();
  check_get_scat_set();
  check_scatter_symmetry();
  check_check_norm();
  check_get_scattering();
  check_get_direct();
  check_number_sums();
  check_sum_legendre();
  check_rotate_phase_matrix();
  check_matrix_symmetry();
  check_makephase();
  check_fftc();
  check_fixreal();
  check_fft1dr();
  check_fft1dr_format();
  check_fourier_basis();
  check_fourier_matrix();
  check_combine_phase_modes();
  check_scattering();
  check_direct_scattering();

  // ARTS's quadratures and Planck function for RT3's; to this is added the
  // rounding amplified by the doublings
  constexpr Numeric tolerance = 1e-10;

  std::mt19937_64 gen(20261008);
  Index           identical = 0, failed = 0, reuse_differ = 0;
  Numeric         worst = 0.0, worst_of_tolerance = 0.0;
  const auto      all   = cases();
  // One work data over all cases, whose sizes differ, and one with every
  // array NaN, against a fresh one per case: neither may change a bit
  rt3::rt3_workdata shared;
  for (const auto& c : all) {
    const auto        in = make_inputs(c, gen);
    rt3::rt3_workdata fresh, poisoned = poisoned_workdata(in);
    const auto        cpp = run_cpp(in, fresh);
    const auto        f77 = run_fortran(in);
    for (rt3::rt3_workdata* work : {&shared, &poisoned}) {
      const auto again  = run_cpp(in, *work);
      reuse_differ     += differ(again.up_rad, cpp.up_rad).first + differ(again.down_rad, cpp.down_rad).first +
                          differ(again.up_flux, cpp.up_flux).first + differ(again.down_flux, cpp.down_flux).first +
                          differ(again.mu_values, cpp.mu_values).first;
    }

    std::string   bad;
    Index         ndiffer        = 0;
    const Numeric case_tolerance = tolerance + rounding * doubling_amplification(in);
    const auto    check          = [&](const char* name, const auto& a, const auto& b) {
      const auto [count, rel]  = differ(a, b);
      ndiffer                 += count;
      worst                    = std::max(worst, rel);
      worst_of_tolerance       = std::max(worst_of_tolerance, rel / case_tolerance);
      if (rel > case_tolerance) bad += std::format(" {}: {} differ, max rel {:.3e};", name, count, rel);
    };
    check("up_rad", cpp.up_rad, f77.up_rad);
    check("down_rad", cpp.down_rad, f77.down_rad);
    check("up_flux", cpp.up_flux, f77.up_flux);
    check("down_flux", cpp.down_flux, f77.down_flux);
    check("mu_values", cpp.mu_values, f77.mu_values);
    identical += ndiffer == 0;
    if (not bad.empty()) {
      failed++;
      std::cout << std::format("Differs by more than {:.2e}: {}:{}\n", case_tolerance, describe(c), bad);
    }
  }

  std::cout << std::format(
      "One rt3_workdata reused over all {} cases (of different sizes), and one with every array NaN, against a "
      "fresh one per case: {} values differ\n",
      all.size(),
      reuse_differ);
  if (reuse_differ > 0) throw std::runtime_error("reusing an rt3_workdata changes the results");

  // The routines refuse a work data that is not sized for their streams
  {
    rt3::rt3_workdata work(2, 4, 0, 0, 0, 0);  // 8 streams
    Tensor3           r(2, 2, 2, 0.0), t(2, 2, 2, 0.0);
    Matrix            src(2, 2, 0.0);
    Vector            top(2, 0.0), bottom(2, 0.0), up(2), down(2);
    bool              refused = false;
    try {
      polradtran::internal_radiance(r, t, src, r, t, src, top, bottom, up, down, work);
    } catch (const std::exception&) { refused = true; }
    if (not refused) throw std::runtime_error("internal_radiance accepted an rt3_workdata for 8 streams with 2");
  }

  // A STOP of the Fortran is an error of the port
  bool threw = false;
  try {
    auto in = make_inputs({2, 6, 0, quadrature_type::gauss, 0, 3, false, 'F', 3, layout::mixed, 1e-6}, gen);
    rt3::rt3_workdata work;
    run_cpp(in, work);
  } catch (const std::exception&) { threw = true; }
  if (not threw) throw std::runtime_error("a solar source over a Fresnel ground did not throw");

  std::cout << std::format(
      "C++ against Fortran RADTRAN: {} of {} cases within {:.0e} plus the rounding amplified by their doublings "
      "(largest relative difference {:.2e}, at most {:.2f} of a case's tolerance), {} bit-identical\n",
      all.size() - failed,
      all.size(),
      tolerance,
      worst,
      worst_of_tolerance,
      identical);
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
} catch (const std::exception& e) {
  std::cerr << "rt3-radtran-test failed: " << e.what() << '\n';
  return EXIT_FAILURE;
}
