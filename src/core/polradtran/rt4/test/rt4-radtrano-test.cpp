// The C++ RADTRANO (radtran4.cc) against the Fortran RADTRANO (radtran4.f)
// it ports: both get the same inputs, and every output must agree to
// tolerance, relative to the largest magnitude in that output.
//
// Each porting step is first checked bit for bit against the previous one.
// The port was bit-identical to the Fortran until RT4's quadratures were
// replaced by ARTS's (rt4::get_quadrature), which differ by rounding: the
// weights by up to 2.4e-12 relative, RT4's being the less accurate.  Then
// ARTS's planck() replaced PLANCK_FUNCTION for the layers, the more accurate
// at small h nu / k T.  Then the port of DOUBLING_INTEGRATION folded the
// MADD, MSUB and MIDENTITY after a product into DGEMM's alpha and beta, and
// inverts with LAPACK instead of LINPACK.
//
// Mathematically equivalent evaluations are accepted: the tolerances allow
// a few ulp (rounding, 16 epsilon) times what a computation amplifies
// rounding by, so that FMA contraction, vectorised libm functions and other
// BLAS kernels pass and an error of the port does not.  n doublings
// amplify rounding by 2^n, so the radiances must agree to 1e-11 (for the
// replaced quadratures and planck()) plus rounding times 2^n for the
// layer doubled most; they differ by up to 4e-13 of the largest radiance
// on Apple arm64 with OpenBLAS and 1.5e-9 on AMD x86_64 with MKL.  The test
// reports how many cases are still bit-identical and the largest
// difference.
//
// It also checks each replaced or ported routine against its Fortran: the
// quadratures, planck() and the ground and sky radiances (lambert_,
// fresnel_, specular_ and thermal_radiance, which use planck()) to rounding,
// fresnel_surface_layer (ARTS's fresnel()), combine_layers and
// internal_radiance to 1e-13, doubling_integration to rounding times 2^n,
// and the routines ported as is (initialize, initial_source,
// nonscatter_layer, lambert_surface_layer, specular_surface_layer,
// external_surface_layer) to 1e-14, nonscatter_layer's source to rounding
// times the cancellation of its terms, reporting how many cases are
// bit-identical (all of them with Apple clang and gfortran; on glibc,
// gfortran vectorises NONSCATTER_LAYER's exp to libmvec's).
//
// The inputs are random, in RADTRANO's own layouts, and cover each branch of
// RADTRANO: the three quadratures, extra angles, the four ground types,
// non-scattering and scattering layers (none, one and many doublings, and
// layers sharing an optics set), a 0 K layer top (LINFACTOR = 0), negative
// gas extinction (clipped in place), and the size limits.  The C++ RADTRANO
// works in SI at the frequency, the Fortran one per micrometre at the
// wavelength; the test converts.
#include <arts_constants.h>
#include <physics_funcs.h>
#include <radintg.h>
#include <radintg4.h>
#include <radtran4.h>
#include <radutil.h>
#include <radutil4.h>
#include <rt4.h>
#include <rt4_c_interface.h>

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
#include <utility>
#include <vector>

namespace rt4 = polradtran::rt4;

namespace {
struct inputs {
  Index                nstokes{}, nummu{}, nuummu{};
  Numeric              max_delta_tau{};
  rt4::quadrature_type quad_type{};
  char                 ground_type{};
  Numeric              ground_temp{}, ground_albedo{}, sky_temp{}, frequency{};
  Complex              ground_index{};
  Matrix               ground_reflec;
  Tensor4              surf_reflect;
  Matrix               gnd_radiance;
  Vector               height, temperatures, gas_extinct, scatlayers;
  Tensor5              extinct_matrix;
  Tensor4              emis_vector;
  Tensor6              scatter_matrix;
  Vector               mu_values;
};

struct outputs {
  Matrix  gnd_radiance;
  Vector  gas_extinct, mu_values;
  Tensor3 up_rad, down_rad;
};

enum class layout { mixed, thin, thick, shared };

struct case_spec {
  Index                nstokes, nquad, nuummu;
  rt4::quadrature_type quad;
  char                 ground;
  Index                nlay;
  layout               lay;
  Numeric              max_delta_tau;
  bool                 zero_kelvin_top{false};
};

//! RADTRANO's QUAD_TYPE, for the Fortran
char fortran_quad_type(rt4::quadrature_type type) {
  switch (type) {
    case rt4::quadrature_type::double_gauss: return 'D';
    case rt4::quadrature_type::gauss:        return 'G';
    case rt4::quadrature_type::lobatto:      return 'L';
  }
  throw std::runtime_error("unknown quadrature type");
}

std::string describe(const case_spec& c) {
  constexpr const char* names[] = {"mixed", "thin", "thick", "shared"};
  return std::format("nstokes {}, nmu {} + {}, quad '{}', ground '{}', {} layers {}, max_delta_tau {:.0e}{}",
                     c.nstokes,
                     c.nquad,
                     c.nuummu,
                     fortran_quad_type(c.quad),
                     c.ground,
                     c.nlay,
                     names[static_cast<int>(c.lay)],
                     c.max_delta_tau,
                     c.zero_kelvin_top ? ", 0 K top" : "");
}

inputs make_inputs(const case_spec& c, std::mt19937_64& gen) {
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  const Index ns = c.nstokes, nmu = c.nquad + c.nuummu, nlay = c.nlay;
  // Layers with optics: every one except every third (non-scattering) in
  // "mixed"; "shared" layers alternate between the first two sets
  ArrayOfIndex set(nlay, 0);
  Index        nsl = 0;
  for (Index l = 0; l < nlay; l++) {
    if (c.lay == layout::mixed and l % 3 == 1) continue;
    if (c.lay == layout::shared and nsl >= 2) {
      set[l] = 1 + l % 2;
      continue;
    }
    set[l] = ++nsl;
  }

  inputs in{.nstokes        = ns,
            .nummu          = nmu,
            .nuummu         = c.nuummu,
            .max_delta_tau  = c.max_delta_tau,
            .quad_type      = c.quad,
            .ground_type    = c.ground,
            .ground_temp    = 287.5,
            .ground_albedo  = 0.27,
            .sky_temp       = 2.73,
            .frequency      = 89e9,
            .ground_index   = Complex{3.1, 0.4},
            .ground_reflec  = Matrix(ns, ns),
            .surf_reflect   = Tensor4(nmu, ns, nmu, ns),
            .gnd_radiance   = Matrix(nmu, ns, 0.0),
            .height         = Vector(nlay + 1),
            .temperatures   = Vector(nlay + 1),
            .gas_extinct    = Vector(nlay),
            .scatlayers     = Vector(nlay),
            .extinct_matrix = Tensor5(std::max<Index>(nsl, 1), 2, nmu, ns, ns, 0.0),
            .emis_vector    = Tensor4(std::max<Index>(nsl, 1), 2, nmu, ns, 0.0),
            .scatter_matrix = Tensor6(std::max<Index>(nsl, 1), 4, nmu, ns, nmu, ns, 0.0),
            .mu_values      = Vector(nmu, 0.0)};

  const Numeric scale = c.lay == layout::thin ? 1e-11 : c.lay == layout::thick ? 5e-3 : 2e-4;
  for (Index l = 0; l <= nlay; l++) {
    in.height[l]       = 750.0 * static_cast<Numeric>(nlay - l);
    in.temperatures[l] = 205.0 + 86.0 * static_cast<Numeric>(l) / static_cast<Numeric>(nlay) + 6.0 * (u(gen) - 0.5);
  }
  if (c.zero_kelvin_top) in.temperatures[0] = 0.0;
  for (Index l = 0; l < nlay; l++) {
    in.gas_extinct[l] = c.lay == layout::thin ? 1e-12 + 1e-9 * u(gen) : 4e-4 * u(gen);
    in.scatlayers[l]  = static_cast<Numeric>(set[l]);
  }
  if (nlay > 2) in.gas_extinct[2] = -1e-5;

  // Mirror-symmetric optics (RADTRANO's SYMMETRIC), K11 above what is scattered
  for (Index s = 0; s < nsl; s++) {
    for (Index i = 0; i < nmu; i++) {
      for (Index si = 0; si < ns; si++) {
        for (Index j = 0; j < nmu; j++) {
          for (Index so = 0; so < ns; so++) {
            const Numeric fwd = scale * (so == si ? u(gen) : 0.1 * (u(gen) - 0.5));
            const Numeric bwd = scale * (so == si ? u(gen) : 0.1 * (u(gen) - 0.5));
            // q = 1 +<-+, 2 +<--, 3 -<-+, 4 -<--: [set, q, in mu, in s, out mu, out s]
            in.scatter_matrix[s, 0, i, si, j, so] = in.scatter_matrix[s, 3, i, si, j, so] = fwd;
            in.scatter_matrix[s, 1, i, si, j, so] = in.scatter_matrix[s, 2, i, si, j, so] = bwd;
          }
        }
      }
      const Numeric k11 = scale * (2.0 * Constant::pi * 2.0 + 0.5 + 2.5 * u(gen));
      const Numeric k12 = 0.05 * k11 * (u(gen) - 0.5);
      const Numeric a1  = scale * (0.5 + 2.5 * u(gen));
      const Numeric a2  = 0.05 * a1 * (u(gen) - 0.5);
      for (Index h = 0; h < 2; h++) {
        for (Index si = 0; si < ns; si++) {
          in.extinct_matrix[s, h, i, si, si] = k11;
          in.emis_vector[s, h, i, si]        = si == 0 ? a1 : a2;
        }
        if (ns > 1) in.extinct_matrix[s, h, i, 0, 1] = in.extinct_matrix[s, h, i, 1, 0] = k12;
      }
    }
  }

  const Numeric R[2][2] = {{0.3, 0.05}, {-0.1, 0.25}};
  for (Index a = 0; a < ns; a++)
    for (Index b = 0; b < ns; b++) in.ground_reflec[a, b] = R[a][b];
  for (auto& x : in.surf_reflect | by_elem) x = 0.02 * u(gen);
  if (c.ground != 'F' and c.ground != 'L' and c.ground != 'S')
    for (Index i = 0; i < nmu; i++)
      for (Index s = 0; s < ns; s++) in.gnd_radiance[i, s] = s == 0 ? 3e-15 * u(gen) : 1e-15 * (u(gen) - 0.5);

  // The extra angles follow the quadrature nodes
  for (Index i = c.nquad; i < nmu; i++) in.mu_values[i] = 1.0 - 0.6 * static_cast<Numeric>(i - c.nquad);
  return in;
}

//! RADTRANO's ground inputs as an rt4::surface
rt4::surface ground_of(const inputs& in) {
  switch (in.ground_type) {
    case 'L': return rt4::lambertian_surface{.albedo = in.ground_albedo};
    case 'F': return rt4::fresnel_surface{.refractive_index = in.ground_index};
    case 'S': return rt4::specular_surface{.reflectivity = in.ground_reflec};
    default:  break;
  }
  // SURF_REFLECT is [in mu, in s, out mu, out s]
  rt4::discrete_surface d{.reflection = Tensor4(in.nummu, in.nummu, in.nstokes, in.nstokes),
                          .emission   = in.gnd_radiance};
  for (Index io = 0; io < in.nummu; io++)
    for (Index ii = 0; ii < in.nummu; ii++)
      for (Index so = 0; so < in.nstokes; so++)
        for (Index si = 0; si < in.nstokes; si++) d.reflection[io, ii, so, si] = in.surf_reflect[ii, si, io, so];
  return d;
}

//! The C++ RADTRANO, with its ground made by rt4::ground_surface on the
//! streams RADTRANO makes
outputs run_cpp(inputs in, polradtran::workdata& work) {
  Tensor3 up_rad(in.height.extent(0), in.nummu, in.nstokes, 0.0);
  Tensor3 down_rad(in.height.extent(0), in.nummu, in.nstokes, 0.0);
  // The extra angles go in on their own; all of mu_values is output
  const Index  nquad = in.nummu - in.nuummu;
  const Vector extra_mu{in.mu_values[Range{nquad, in.nuummu}]};
  const auto   q = rt4::get_quadrature(nquad, in.quad_type);
  Vector       mu(in.nummu), w(in.nummu, 0.0);
  mu[Range{0, nquad}]         = q.mu;
  mu[Range{nquad, in.nuummu}] = extra_mu;
  w[Range{0, nquad}]          = q.weights;

  Tensor4 surf_reflect(in.nummu, in.nstokes, in.nummu, in.nstokes);
  Matrix  gnd_radiance(in.nummu, in.nstokes);
  rt4::ground_surface(ground_of(in), mu, w, in.frequency, in.ground_temp, surf_reflect, gnd_radiance);

  in.mu_values = std::numeric_limits<Numeric>::quiet_NaN();
  rt4::radtrano(in.max_delta_tau,
                in.quad_type,
                surf_reflect,
                gnd_radiance,
                in.sky_temp,
                in.frequency,
                in.height,
                in.temperatures,
                in.gas_extinct,
                in.scatlayers,
                in.extinct_matrix,
                in.emis_vector,
                in.scatter_matrix,
                extra_mu,
                in.mu_values,
                up_rad,
                down_rad,
                work);
  return {.gnd_radiance = std::move(gnd_radiance),
          .gas_extinct  = std::move(in.gas_extinct),
          .mu_values    = std::move(in.mu_values),
          .up_rad       = std::move(up_rad),
          .down_rad     = std::move(down_rad)};
}

//! The Fortran RADTRANO works per micrometre at the wavelength in um; its
//! radiances are converted to SI with B_nu = B_lambda * lambda / nu
outputs run_fortran(inputs in, Index nsl) {
  Tensor3       up_rad(in.height.extent(0), in.nummu, in.nstokes, 0.0);
  Tensor3       down_rad(in.height.extent(0), in.nummu, in.nstokes, 0.0);
  const Numeric wavelength        = 1e6 * Constant::c / in.frequency;
  const Numeric per_um_to_per_hz  = wavelength / in.frequency;
  in.gnd_radiance                /= per_um_to_per_hz;
  rt4_radtrano(in.nstokes,
               in.nummu,
               in.nuummu,
               in.max_delta_tau,
               fortran_quad_type(in.quad_type),
               in.ground_temp,
               in.ground_type,
               in.ground_albedo,
               in.ground_index.real(),
               in.ground_index.imag(),
               in.ground_reflec.data_handle(),
               in.surf_reflect.data_handle(),
               in.gnd_radiance.data_handle(),
               in.sky_temp,
               wavelength,
               in.gas_extinct.extent(0),
               in.height.data_handle(),
               in.temperatures.data_handle(),
               in.gas_extinct.data_handle(),
               nsl,
               in.scatlayers.data_handle(),
               in.extinct_matrix.data_handle(),
               in.emis_vector.data_handle(),
               in.scatter_matrix.data_handle(),
               in.mu_values.data_handle(),
               up_rad.data_handle(),
               down_rad.data_handle());
  in.gnd_radiance *= per_um_to_per_hz;
  up_rad          *= per_um_to_per_hz;
  down_rad        *= per_um_to_per_hz;
  return {.gnd_radiance = std::move(in.gnd_radiance),
          .gas_extinct  = std::move(in.gas_extinct),
          .mu_values    = std::move(in.mu_values),
          .up_rad       = std::move(up_rad),
          .down_rad     = std::move(down_rad)};
}

//! 2^n for the most doublings n of a scattering layer, as RADTRANO counts
//! them: the amplification of rounding (see check_doubling_integration)
Numeric doubling_amplification(const inputs& in) {
  Numeric amp = 1.0;
  for (Index l = 0; l < in.scatlayers.extent(0); l++) {
    const Index set = std::lround(in.scatlayers[l]);
    if (set < 1) continue;
    const Numeric extinct = in.extinct_matrix[set - 1, 0, 0, 0, 0] + std::max(in.gas_extinct[l], 0.0);
    const Numeric f = std::log2(std::max(extinct * std::abs(in.height[l] - in.height[l + 1]), 1e-7) / in.max_delta_tau);
    if (f > 0.0) amp = std::max(amp, std::pow(2.0, std::floor(f) + 1.0));
  }
  return amp;
}

//! The number of elements that differ, and the largest difference relative
//! to the largest magnitude in either array (Q can be near 0 beside I)
std::pair<Index, Numeric> differ(const auto& a, const auto& b) {
  Index          count = 0;
  Numeric        diff = 0.0, scale = 0.0;
  const Numeric* x = a.data_handle();
  const Numeric* y = b.data_handle();
  for (std::size_t k = 0; k < a.size(); k++) {
    scale = std::max({scale, std::abs(x[k]), std::abs(y[k])});
    if (x[k] == y[k]) continue;
    count++;
    const Numeric d = std::abs(x[k] - y[k]);
    diff            = std::max(diff, std::isnan(d) ? std::numeric_limits<Numeric>::infinity() : d);
  }
  return {count, count == 0 ? 0.0 : diff / scale};
}

/* rt4::get_quadrature against RT4's quadrature routines, which it replaces
   in RADTRANO: the same rules, to rounding.  For nmu up to 64 the nodes
   differ by 4.4e-16 and the weights by 2.4e-12 relative (Apple arm64); the
   test allows 1e-15 and 1e-11 for other compilers and libms. */
void check_quadratures() {
  using fortran_rule = void (*)(std::int64_t, double*, double*);
  Numeric dmu = 0.0, dw = 0.0;
  for (Index n = 1; n <= 64; n++) {
    for (auto [type, fortran] :
         {std::pair<rt4::quadrature_type, fortran_rule>{rt4::quadrature_type::double_gauss,
                                                        rt4_double_gauss_quadrature},
          std::pair<rt4::quadrature_type, fortran_rule>{rt4::quadrature_type::gauss, rt4_gauss_legendre_quadrature},
          std::pair<rt4::quadrature_type, fortran_rule>{rt4::quadrature_type::lobatto, rt4_lobatto_quadrature}}) {
      const auto q = rt4::get_quadrature(n, type);
      Vector     mu(n), w(n);
      fortran(n, mu.data_handle(), w.data_handle());
      for (Index i = 0; i < n; i++) {
        dmu = std::max(dmu, std::abs(q.mu[i] - mu[i]));
        dw  = std::max(dw, std::abs(q.weights[i] - w[i]) / w[i]);
      }
    }
  }
  std::cout << std::format(
      "rt4::get_quadrature against RT4's D, G and L routines, nmu 1 to 64: nodes within {:.2e}, weights within "
      "{:.2e} relative\n",
      dmu,
      dw);
  if (dmu > 1e-15 or dw > 1e-11) throw std::runtime_error("the quadratures differ by more than rounding");
}

/* ARTS's planck() against RT4's PLANCK_FUNCTION, which it replaces for
   the layers in RADTRANO: the same function on the same exact SI
   constants, but RT4 evaluates exp(x) - 1 where planck() has expm1(x), so
   RT4 loses digits at small x = h nu / k T. */
void check_planck() {
  Numeric worst = 0.0, worst_f = 0.0, worst_t = 0.0;
  for (Numeric f : {1e9, 10e9, 89e9, 183e9, 664e9, 3e12}) {
    for (Numeric t : {Constant::cosmic_microwave_background_temperature, 50.0, 150.0, 250.0, 330.0}) {
      const Numeric wavelength = 1e6 * Constant::c / f;
      const Numeric rt4        = rt4_planck_function(t, wavelength) * wavelength / f;
      const Numeric diff       = std::abs(planck(f, t) - rt4) / planck(f, t);
      if (diff > worst) {
        worst   = diff;
        worst_f = f;
        worst_t = t;
      }
    }
  }
  std::cout << std::format(
      "planck() against RT4's PLANCK_FUNCTION, 1 GHz to 3 THz and 2.7 to 330 K: within {:.2e} relative (at {:.0f} "
      "GHz, {} K)\n",
      worst,
      worst_f / 1e9,
      worst_t);
  if (worst > 1e-12) throw std::runtime_error("planck() and PLANCK_FUNCTION differ by more than rounding");
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

constexpr Numeric nan = std::numeric_limits<Numeric>::quiet_NaN();

/* The routines ported as is against their Fortran.  The C++ outputs start
   as NaN, so an element a routine does not write differs. */
/* polradtran::lambert_surface_layer against LAMBERT_SURFACE, ported as is, for mode 0
   (RT4's) and 1 (no reflection), with zero-weight extra angles. */
void check_lambert_surface() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2}) {
    for (Index nummu : {Index{1}, Index{5}, 64 / nstokes}) {
      for (Index mode : {0, 1}) {
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
        rt4_lambert_surface(nstokes,
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
  tally.report("lambert_surface_layer", "LAMBERT_SURFACE");
}

/* polradtran::lambert_radiance against LAMBERT_RADIANCE, whose PLANCK_FUNCTION per
   micrometre is converted to SI: they differ as planck() and
   PLANCK_FUNCTION do (check_planck). */
void check_lambert_radiance() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2}) {
    for (Numeric f : {1e9, 89e9, 664e9, 3e12}) {
      for (Numeric t : {0.0, Constant::cosmic_microwave_background_temperature, 150.0, 330.0}) {
        const Index   nummu  = 7;
        const Numeric albedo = u(gen);

        Matrix rad(nummu, nstokes, nan), radf(nummu, nstokes);
        polradtran::lambert_radiance(0, 2, albedo, t, f, 0.0, rad);

        const Numeric wavelength = 1e6 * Constant::c / f;
        rt4_lambert_radiance(nstokes, nummu, albedo, t, wavelength, radf.data_handle());
        radf *= wavelength / f;
        tally.add({differ(rad, radf)});
      }
    }
  }
  tally.report("lambert_radiance", "LAMBERT_RADIANCE", 2e-12);

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

/* polradtran::fresnel_surface_layer against FRESNEL_SURFACE, for up to 4 Stokes
   components (the [U, V] block too).  Not bit-identical: the amplitudes
   are ARTS's fresnel() at acos(mu) in degrees, with the transmitted cosine
   sqrt(1 - sin^2 / n^2) where RT4 has sqrt(n^2 - sin^2), and C++ complex
   arithmetic. */
void check_fresnel_surface() {
  std::mt19937_64                         gen(1995);
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
      rt4_fresnel_surface(nstokes,
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
  tally.report("fresnel_surface_layer", "FRESNEL_SURFACE", 1e-13);
}

/* polradtran::fresnel_radiance against FRESNEL_RADIANCE, whose PLANCK_FUNCTION per
   micrometre is converted to SI: they differ as planck() and
   PLANCK_FUNCTION do (6e-13 at 1 GHz and 150 K, where RT4's exp(x) - 1
   loses digits), and in the reflection to rounding. */
void check_fresnel_radiance() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2}) {
    for (const Complex& index : fresnel_indices) {
      for (Numeric f : {1e9, 89e9, 3e12}) {
        for (Numeric t : {0.0, 150.0, 300.0}) {
          const Index nummu = 9;
          Vector      mu(nummu);
          for (Index j = 0; j < nummu - 1; j++) mu[j] = 0.01 + 0.98 * u(gen);
          mu[nummu - 1] = 1.0;

          Matrix rad(nummu, nstokes, nan), radf(nummu, nstokes);
          polradtran::fresnel_radiance(0, mu, index, t, f, rad);

          const Numeric wavelength = 1e6 * Constant::c / f;
          rt4_fresnel_radiance(
              nstokes, nummu, mu.data_handle(), index.real(), index.imag(), t, wavelength, radf.data_handle());
          radf *= wavelength / f;
          tally.add({differ(rad, radf)});
        }
      }
    }
  }
  tally.report("fresnel_radiance", "FRESNEL_RADIANCE", 2e-12);
}

/* rt4::specular_surface_layer against SPECULAR_SURFACE, ported as is, for
   up to 4 Stokes components and a non-symmetric reflectivity R(out, in)
   (passed row-major, as SPECULAR_SURFACE reads it transposed). */
void check_specular_surface() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (Index nummu : {Index{1}, Index{5}, 64 / nstokes}) {
      Matrix ref(nstokes, nstokes);
      for (auto& x : ref | by_elem) x = 0.3 * (u(gen) - 0.2);

      Tensor5 r(2, nummu, nstokes, nummu, nstokes, nan), t(2, nummu, nstokes, nummu, nstokes, nan);
      Tensor3 src(2, nummu, nstokes, nan);
      rt4::specular_surface_layer(ref, r, t, src);

      Tensor5 rf(2, nummu, nstokes, nummu, nstokes), tf(2, nummu, nstokes, nummu, nstokes);
      Tensor3 srcf(2, nummu, nstokes);
      rt4_specular_surface(nstokes, nummu, ref.data_handle(), rf.data_handle(), tf.data_handle(), srcf.data_handle());
      tally.add({differ(r, rf), differ(t, tf), differ(src, srcf)});
    }
  }
  tally.report("specular_surface_layer", "SPECULAR_SURFACE");
}

/* rt4::specular_radiance against SPECULAR_RADIANCE (I and Q, RT4's Stokes
   components), whose PLANCK_FUNCTION per micrometre is converted to SI:
   they differ as planck() and PLANCK_FUNCTION do. */
void check_specular_radiance() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2}) {
    for (Numeric f : {1e9, 89e9, 3e12}) {
      for (Numeric t : {0.0, 150.0, 300.0}) {
        const Index nummu = 7;
        Matrix      ref(nstokes, nstokes);
        for (auto& x : ref | by_elem) x = 0.3 * (u(gen) - 0.2);

        Matrix rad(nummu, nstokes, nan), radf(nummu, nstokes);
        rt4::specular_radiance(ref, t, f, rad);

        const Numeric wavelength = 1e6 * Constant::c / f;
        rt4_specular_radiance(nstokes, nummu, ref.data_handle(), t, wavelength, radf.data_handle());
        radf *= wavelength / f;
        tally.add({differ(rad, radf)});
      }
    }
  }
  tally.report("specular_radiance", "SPECULAR_RADIANCE", 2e-12);
}

/* polradtran::external_surface_layer against EXTERNAL_SURFACE, ported as is (its
   unused RADIANCE argument dropped). */
void check_external_surface() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2, 3, 4}) {
    for (Index nummu : {Index{1}, Index{5}, 64 / nstokes}) {
      Tensor4 refl(nummu, nstokes, nummu, nstokes);
      for (auto& x : refl | by_elem) x = 0.05 * (u(gen) - 0.2);
      Matrix radiance(nummu, nstokes, 1.0);

      Tensor5 r(2, nummu, nstokes, nummu, nstokes, nan), t(2, nummu, nstokes, nummu, nstokes, nan);
      Tensor3 src(2, nummu, nstokes, nan);
      polradtran::external_surface_layer(refl, r, t, src);

      Tensor5 rf(2, nummu, nstokes, nummu, nstokes), tf(2, nummu, nstokes, nummu, nstokes);
      Tensor3 srcf(2, nummu, nstokes);
      rt4_external_surface(nstokes,
                           nummu,
                           refl.data_handle(),
                           radiance.data_handle(),
                           rf.data_handle(),
                           tf.data_handle(),
                           srcf.data_handle());
      tally.add({differ(r, rf), differ(t, tf), differ(src, srcf)});
    }
  }
  tally.report("external_surface_layer", "EXTERNAL_SURFACE");
}

/* polradtran::thermal_radiance against THERMAL_RADIANCE, whose PLANCK_FUNCTION per
   micrometre is converted to SI: they differ as planck() and
   PLANCK_FUNCTION do. */
void check_thermal_radiance() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2}) {
    for (Numeric f : {1e9, 89e9, 3e12}) {
      for (Numeric t : {0.0, Constant::cosmic_microwave_background_temperature, 150.0, 300.0}) {
        for (bool zero_albedo : {true, false}) {
          const Index   nummu  = 7;
          const Numeric albedo = zero_albedo ? 0.0 : u(gen);

          Tensor3 rad(2, nummu, nstokes, nan), radf(2, nummu, nstokes);
          polradtran::thermal_radiance(0, t, albedo, f, rad);

          const Numeric wavelength = 1e6 * Constant::c / f;
          rt4_thermal_radiance(nstokes, nummu, t, albedo, wavelength, radf.data_handle());
          radf *= wavelength / f;
          tally.add({differ(rad, radf)});
        }
      }
    }
  }
  tally.report("thermal_radiance", "THERMAL_RADIANCE", 2e-12);

  bool threw = false;
  try {
    Tensor3 rad(2, 3, 2);
    polradtran::thermal_radiance(0, -1.0, 0.0, 89e9, rad);
  } catch (const std::exception&) { threw = true; }
  if (not threw) throw std::runtime_error("thermal_radiance at -1 K did not throw");
}

/* polradtran::doubling_integration against DOUBLING_INTEGRATION, from the thin
   initial sublayer of random, physical optics (single-scattering albedo 0.5
   to 0.95) made by rt4::initialize and rt4::initial_source.  It is not as
   is: DGEMM's alpha and beta absorb the MIDENTITY, MSUB and MADD after a
   product (one rounding fewer each), and the inverse is LAPACK's instead of
   LINPACK's, so the two agree to rounding amplified by the doublings.  Each
   doubling about squares T, which doubles its relative error, so n
   doublings amplify rounding by 2^n: perturbing the Fortran's inputs by
   1 ulp changes its output by 3.7e-9 for 24 doublings, and in quad
   precision the C++ and the Fortran are both 1.9e-9 off.  The C++ and the
   Fortran differ by 2e-14 for 6 doublings and 1.2e-9 for 24 (AMD x86_64,
   MKL; 1.3e-11 on Apple arm64 with OpenBLAS). */
void check_doubling_integration() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2}) {
    for (Index nummu : {Index{1}, Index{4}, 32 / nstokes}) {
      for (Index num_doubles : {0, 1, 6, 24}) {
        for (bool symmetric : {true, false}) {
          const Index n = nstokes * nummu;
          const auto  q = rt4::get_quadrature(nummu, rt4::quadrature_type::double_gauss);

          // Extinction k, scattering omega k spread over the streams; with
          // symmetric, the minus hemisphere mirrors the plus one
          Tensor4 ext(2, nummu, nstokes, nstokes, 0.0);
          Tensor5 sca(4, nummu, nstokes, nummu, nstokes);
          Tensor3 emis(2, nummu, nstokes, 0.0);
          for (Index h = 0; h < 2; h++) {
            for (Index j = 0; j < nummu; j++) {
              const Numeric k = 0.5 + u(gen), omega = 0.5 + 0.45 * u(gen);
              for (Index s = 0; s < nstokes; s++) ext[h, j, s, s] = k;
              if (nstokes > 1) ext[h, j, 0, 1] = ext[h, j, 1, 0] = 0.05 * k * (u(gen) - 0.5);
              emis[h, j, 0] = (1.0 - omega) * k;
            }
          }
          for (auto& z : sca | by_elem) z = (0.5 + 0.45 * u(gen)) * (0.5 + u(gen)) / (4 * Constant::pi) * u(gen);
          if (symmetric) {
            ext[1]  = ext[0];
            emis[1] = emis[0];
            sca[3]  = sca[0];
            sca[2]  = sca[1];
          }

          Tensor3 reflect(2, n, n), trans(2, n, n);
          Matrix  lin(2, n);
          rt4::initialize(1e-6,
                          q.mu,
                          q.weights,
                          1e-2,
                          ext,
                          sca,
                          reflect.view_as(2, nummu, nstokes, nummu, nstokes),
                          trans.view_as(2, nummu, nstokes, nummu, nstokes));
          rt4::initial_source(1e-6, q.mu, 1e-15, emis, 1e-2, lin.view_as(2, nummu, nstokes));
          const Numeric linfactor = 0.3 / std::pow(2.0, num_doubles);

          Tensor3              rf{reflect}, tf{trans};
          Matrix               lf{lin};
          Tensor3              tr(2, n, n, nan), tt(2, n, n, nan), trf(2, n, n), ttf(2, n, n);
          Matrix               ts(2, n, nan), tsf(2, n);
          polradtran::workdata work(nstokes, nummu, 0);
          polradtran::doubling_integration(num_doubles, symmetric, reflect, trans, lin, linfactor, tr, tt, ts, work);
          rt4_doubling_integration(n,
                                   num_doubles,
                                   symmetric,
                                   rf.data_handle(),
                                   tf.data_handle(),
                                   lf.data_handle(),
                                   linfactor,
                                   trf.data_handle(),
                                   ttf.data_handle(),
                                   tsf.data_handle());
          const Numeric amp = std::pow(2.0, num_doubles);
          tally.add({{differ(tr, trf), amp},
                     {differ(tt, ttf), amp},
                     {differ(ts, tsf), amp},
                     {differ(reflect, rf), amp},
                     {differ(trans, tf), amp},
                     {differ(lin, lf), amp}});
        }
      }
    }
  }
  tally.report("doubling_integration", "DOUBLING_INTEGRATION", rounding);
}

//! A physical slab: R, T and S of random optics (single-scattering albedo
//! 0.5 to 0.95), num_doubles doublings of a 1e-6 thick initial sublayer
struct slab {
  Tensor3 reflect, trans;
  Matrix  source;
};

slab random_slab(std::mt19937_64& gen, Index nstokes, Index nummu, Index num_doubles, bool symmetric) {
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  const Index n = nstokes * nummu;
  const auto  q = rt4::get_quadrature(nummu, rt4::quadrature_type::double_gauss);
  Tensor4     ext(2, nummu, nstokes, nstokes, 0.0);
  Tensor5     sca(4, nummu, nstokes, nummu, nstokes);
  Tensor3     emis(2, nummu, nstokes, 0.0);
  for (Index h = 0; h < 2; h++) {
    for (Index j = 0; j < nummu; j++) {
      const Numeric k = 0.5 + u(gen), omega = 0.5 + 0.45 * u(gen);
      for (Index s = 0; s < nstokes; s++) ext[h, j, s, s] = k;
      if (nstokes > 1) ext[h, j, 0, 1] = ext[h, j, 1, 0] = 0.05 * k * (u(gen) - 0.5);
      emis[h, j, 0] = (1.0 - omega) * k;
    }
  }
  for (auto& z : sca | by_elem) z = (0.5 + 0.45 * u(gen)) * (0.5 + u(gen)) / (4 * Constant::pi) * u(gen);
  if (symmetric) {
    ext[1]  = ext[0];
    emis[1] = emis[0];
    sca[3]  = sca[0];
    sca[2]  = sca[1];
  }

  Tensor3 r1(2, n, n), t1(2, n, n);
  Matrix  lin(2, n);
  rt4::initialize(1e-6,
                  q.mu,
                  q.weights,
                  1e-2,
                  ext,
                  sca,
                  r1.view_as(2, nummu, nstokes, nummu, nstokes),
                  t1.view_as(2, nummu, nstokes, nummu, nstokes));
  rt4::initial_source(1e-6, q.mu, 1e-15 * (0.5 + u(gen)), emis, 1e-2, lin.view_as(2, nummu, nstokes));

  slab                 out{.reflect = Tensor3(2, n, n), .trans = Tensor3(2, n, n), .source = Matrix(2, n)};
  polradtran::workdata work(nstokes, nummu, 0);
  polradtran::doubling_integration(
      num_doubles, symmetric, r1, t1, lin, 0.3 / std::pow(2.0, num_doubles), out.reflect, out.trans, out.source, work);
  return out;
}

/* polradtran::combine_layers against COMBINE_LAYERS: pairs of physical slabs,
   thin and thick, symmetric or not, and a slab on a Lambertian ground's
   surface layer.  Not as is, like doubling_integration: DGEMM's beta
   absorbs the MADD and MSUB after a product, and the inverse is LAPACK's. */
void check_combine_layers() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2}) {
    for (Index nummu : {Index{1}, Index{4}, 32 / nstokes}) {
      const Index n = nstokes * nummu;
      for (Index d1 : {0, 6, 20}) {
        for (Index d2 : {0, 6, 20, -1}) {  // -1: the ground
          for (bool symmetric : {true, false}) {
            const slab top = random_slab(gen, nstokes, nummu, d1, symmetric);
            slab       bottom;
            if (d2 < 0) {
              const auto q = rt4::get_quadrature(nummu, rt4::quadrature_type::double_gauss);
              Tensor5    r(2, nummu, nstokes, nummu, nstokes), t(2, nummu, nstokes, nummu, nstokes);
              Tensor3    src(2, nummu, nstokes);
              polradtran::lambert_surface_layer(0, q.mu, q.weights, 0.3, r, t, src);
              bottom = {.reflect = Tensor3{r.view_as(2, n, n)},
                        .trans   = Tensor3{t.view_as(2, n, n)},
                        .source  = Matrix{src.view_as(2, n)}};
            } else {
              bottom = random_slab(gen, nstokes, nummu, d2, symmetric);
            }

            Tensor3              r(2, n, n, nan), t(2, n, n, nan), rf(2, n, n), tf(2, n, n);
            Matrix               src(2, n, nan), srcf(2, n);
            polradtran::workdata work(nstokes, nummu, 0);
            polradtran::combine_layers(
                top.reflect, top.trans, top.source, bottom.reflect, bottom.trans, bottom.source, r, t, src, work);

            // The Fortran declares no intent: give it copies
            slab a = top, b = bottom;
            rt4_combine_layers(n,
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
  tally.report("combine_layers", "COMBINE_LAYERS", 1e-13);
}

/* polradtran::internal_radiance against INTERNAL_RADIANCE: the radiances at a level
   between physical slabs, symmetric or not, and at the top, under vacuum
   (R = 0, T = 1, S = 0, as RADTRANO has it).  Not as is, like
   combine_layers. */
void check_internal_radiance() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2}) {
    for (Index nummu : {Index{1}, Index{4}, 32 / nstokes}) {
      const Index n = nstokes * nummu;
      for (Index d1 : {-1, 0, 6, 20}) {  // -1: vacuum above (the top level)
        for (Index d2 : {0, 6, 20}) {
          for (bool symmetric : {true, false}) {
            slab above;
            if (d1 < 0) {
              above = {.reflect = Tensor3(2, n, n, 0.0), .trans = Tensor3(2, n, n), .source = Matrix(2, n, 0.0)};
              for (Index h = 0; h < 2; h++) identity(above.trans[h]);
            } else {
              above = random_slab(gen, nstokes, nummu, d1, symmetric);
            }
            const slab below = random_slab(gen, nstokes, nummu, d2, symmetric);
            Vector     top(n), bottom(n);
            for (auto& x : top) x = 1e-17 * u(gen);
            for (auto& x : bottom) x = 1e-15 * u(gen);

            Vector               up(n, nan), down(n, nan), upf(n), downf(n);
            polradtran::workdata work(nstokes, nummu, 0);
            polradtran::internal_radiance(above.reflect,
                                          above.trans,
                                          above.source,
                                          below.reflect,
                                          below.trans,
                                          below.source,
                                          top,
                                          bottom,
                                          up,
                                          down,
                                          work);

            // The Fortran declares no intent: give it copies
            slab   a = above, b = below;
            Vector tf{top}, bf{bottom};
            rt4_internal_radiance(n,
                                  a.reflect.data_handle(),
                                  a.trans.data_handle(),
                                  a.source.data_handle(),
                                  b.reflect.data_handle(),
                                  b.trans.data_handle(),
                                  b.source.data_handle(),
                                  tf.data_handle(),
                                  bf.data_handle(),
                                  upf.data_handle(),
                                  downf.data_handle());
            tally.add({differ(up, upf), differ(down, downf)});
          }
        }
      }
    }
  }
  tally.report("internal_radiance", "INTERNAL_RADIANCE", 1e-13);
}

void check_initialize() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2}) {
    for (Index nummu : {Index{1}, Index{5}, 64 / nstokes}) {
      for (Numeric delta_z : {0.0, 1e-9, 1e-3, 1.0}) {
        for (bool gas : {false, true}) {
          // The last stream is an extra angle of weight 0 when there are several
          Vector mu(nummu), w(nummu);
          for (Index j = 0; j < nummu; j++) {
            mu[j] = 0.02 + 0.98 * u(gen);
            w[j]  = nummu > 1 and j == nummu - 1 ? 0.0 : u(gen) / static_cast<Numeric>(nummu);
          }
          Tensor4 ext(2, nummu, nstokes, nstokes);
          Tensor5 sca(4, nummu, nstokes, nummu, nstokes);
          for (auto& e : ext | by_elem) e = 1e-3 * (u(gen) - 0.1);
          for (auto& z : sca | by_elem) z = 1e-4 * (u(gen) - 0.1);
          const Numeric gas_extinct = gas ? 1e-3 * u(gen) : 0.0;

          Tensor5 r(2, nummu, nstokes, nummu, nstokes, nan), t(2, nummu, nstokes, nummu, nstokes, nan);
          rt4::initialize(delta_z, mu, w, gas_extinct, ext, sca, r, t);

          Tensor5 rf(2, nummu, nstokes, nummu, nstokes), tf(2, nummu, nstokes, nummu, nstokes);
          rt4_initialize(nstokes,
                         nummu,
                         delta_z,
                         mu.data_handle(),
                         w.data_handle(),
                         gas_extinct,
                         ext.data_handle(),
                         sca.data_handle(),
                         rf.data_handle(),
                         tf.data_handle());
          tally.add({differ(r, rf), differ(t, tf)});
        }
      }
    }
  }
  tally.report("initialize", "INITIALIZE");
}

void check_initial_source() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2}) {
    for (Index nummu : {Index{1}, Index{5}, 64 / nstokes}) {
      for (Numeric delta_z : {0.0, 1e-9, 1.0, 750.0}) {
        for (bool zero_planck : {false, true}) {
          for (bool gas : {false, true}) {
            Vector mu(nummu);
            for (auto& m : mu) m = 0.02 + 0.98 * u(gen);
            Tensor3 emis(2, nummu, nstokes);
            for (auto& e : emis | by_elem) e = 1e-4 * (u(gen) - 0.2);
            const Numeric planck      = zero_planck ? 0.0 : 1e-15 * u(gen);
            const Numeric gas_extinct = gas ? 1e-3 * u(gen) : 0.0;

            Tensor3 src(2, nummu, nstokes, nan), srcf(2, nummu, nstokes);
            rt4::initial_source(delta_z, mu, planck, emis, gas_extinct, src);
            rt4_initial_source(
                nstokes, nummu, delta_z, mu.data_handle(), planck, emis.data_handle(), gas_extinct, srcf.data_handle());
            tally.add({differ(src, srcf)});
          }
        }
      }
    }
  }
  tally.report("initial_source", "INITIAL_SOURCE");
}

void check_nonscatter_layer() {
  std::mt19937_64                         gen(1995);
  std::uniform_real_distribution<Numeric> u(0.0, 1.0);

  routine_tally tally;
  for (Index nstokes : {1, 2}) {
    for (Index nummu : {Index{1}, Index{5}, 64 / nstokes}) {
      for (Numeric deltatau : {0.0, 1e-12, 1e-4, 0.3, 7.0, 60.0}) {
        for (int b = 0; b < 3; b++) {
          Vector mu(nummu);
          for (auto& m : mu) m = 0.02 + 0.98 * u(gen);
          const Numeric planck0 = b == 2 ? 0.0 : 1e-15 * u(gen);
          const Numeric planck1 = b == 1 ? planck0 : 1e-15 * u(gen);

          Tensor5 r(2, nummu, nstokes, nummu, nstokes, nan), t(2, nummu, nstokes, nummu, nstokes, nan);
          Tensor3 src(2, nummu, nstokes, nan);
          polradtran::nonscatter_layer(0, deltatau, mu, planck0, planck1, r, t, src);

          Tensor5 rf(2, nummu, nstokes, nummu, nstokes), tf(2, nummu, nstokes, nummu, nstokes);
          Tensor3 srcf(2, nummu, nstokes);
          rt4_nonscatter_layer(nstokes,
                               nummu,
                               deltatau,
                               mu.data_handle(),
                               planck0,
                               planck1,
                               rf.data_handle(),
                               tf.data_handle(),
                               srcf.data_handle());
          // The source's terms, of size planck + slope, cancel to about
          // planck * path, so at a small path an ulp of exp(-path) (gfortran
          // vectorises it to libmvec's, up to 3.5 ulp off) is amplified by
          // the sum of their magnitudes over the largest source
          Numeric terms = 0.0, largest = 0.0;
          if (deltatau > 0.0) {
            for (Index j = 0; j < nummu; j++) {
              const Numeric path = deltatau / mu[j], slope = std::abs(planck1 - planck0) / path,
                            p = std::max(planck0, planck1);
              terms           = std::max(terms, p + slope + (p + slope * (1.0 + path)) * std::exp(-path));
            }
          }
          for (auto s : src | by_elem) largest = std::max(largest, std::abs(s));
          const Numeric amp = largest > 0.0 ? std::max(1.0, terms / largest) : 1.0;
          tally.add({differ(r, rf), differ(t, tf), {differ(src, srcf), amp}});
        }
      }
    }
  }
  tally.report("nonscatter_layer", "NONSCATTER_LAYER", rounding);
}

std::vector<case_spec> cases() {
  constexpr auto D = rt4::quadrature_type::double_gauss, G = rt4::quadrature_type::gauss,
                 L = rt4::quadrature_type::lobatto;

  std::vector<case_spec> out;
  for (Index ns : {1, 2})
    for (auto quad : {D, G, L})
      for (Index extra : {0, 1, 2})
        for (char ground : {'L', 'F', 'S', 'A'}) out.push_back({ns, 6, extra, quad, ground, 5, layout::mixed, 1e-6});
  for (layout lay : {layout::thin, layout::thick, layout::shared})
    for (char ground : {'L', 'A'}) out.push_back({2, 8, 1, D, ground, 7, lay, 1e-6});
  for (Numeric mdt : {1e-3, 1e-8}) out.push_back({2, 8, 0, D, 'F', 4, layout::mixed, mdt});
  out.push_back({2, 8, 1, D, 'L', 4, layout::mixed, 1e-6, true});
  out.push_back({2, 30, 2, D, 'F', 4, layout::mixed, 1e-6});   // N = 64 = MAXV
  out.push_back({1, 1, 0, D, 'L', 400, layout::mixed, 1e-6});  // NUM_LAYERS = MAXLAY
  out.push_back({1, 1, 0, G, 'S', 1, layout::thick, 1e-6});    // one layer
  return out;
}
}  // namespace

int main() try {
  // ARTS's quadratures and planck() for RT4's, which differ by 2.3e-12 and
  // 4.5e-13; to this is added the rounding amplified by the doublings
  constexpr Numeric tolerance = 1e-11;

  check_quadratures();
  check_planck();
  check_lambert_surface();
  check_lambert_radiance();
  check_fresnel_surface();
  check_fresnel_radiance();
  check_specular_surface();
  check_specular_radiance();
  check_external_surface();
  check_thermal_radiance();
  check_doubling_integration();
  check_combine_layers();
  check_internal_radiance();
  check_initialize();
  check_initial_source();
  check_nonscatter_layer();

  std::mt19937_64 gen(20261007);
  Index           identical = 0, failed = 0, reuse_differ = 0;
  Numeric         worst = 0.0, worst_of_tolerance = 0.0;
  const auto      all = cases();

  // One work data over all cases, whose sizes differ, against a fresh one
  // per case: reuse must not change a bit
  polradtran::workdata shared;
  for (const auto& c : all) {
    const auto in  = make_inputs(c, gen);
    Index      nsl = 0;
    for (auto t : in.scatlayers) nsl = std::max(nsl, static_cast<Index>(t));
    polradtran::workdata fresh;
    const auto           cpp    = run_cpp(in, fresh);
    const auto           reused = run_cpp(in, shared);
    const auto           f77    = run_fortran(in, nsl);
    for (auto [count, rel] : {differ(cpp.up_rad, reused.up_rad),
                              differ(cpp.down_rad, reused.down_rad),
                              differ(cpp.gnd_radiance, reused.gnd_radiance),
                              differ(cpp.mu_values, reused.mu_values)})
      reuse_differ += count;

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
    check("gnd_radiance", cpp.gnd_radiance, f77.gnd_radiance);
    check("gas_extinct", cpp.gas_extinct, f77.gas_extinct);
    check("mu_values", cpp.mu_values, f77.mu_values);
    identical += ndiffer == 0;
    if (not bad.empty()) {
      failed++;
      std::cout << std::format("Differs by more than {:.2e}: {}:{}\n", case_tolerance, describe(c), bad);
    }
  }

  // A STOP of the Fortran is an error of the port
  bool threw = false;
  try {
    auto in = make_inputs({1, 1, 0, rt4::quadrature_type::double_gauss, 'L', 401, layout::mixed, 1e-6}, gen);
    polradtran::workdata work;
    run_cpp(in, work);
  } catch (const std::exception&) { threw = true; }
  if (not threw) throw std::runtime_error("NUM_LAYERS = 401 > MAXLAY did not throw");

  std::cout << std::format(
      "One workdata reused over all {} cases (of different sizes) against a fresh one per case: {} values "
      "differ\n",
      all.size(),
      reuse_differ);
  if (reuse_differ > 0) throw std::runtime_error("reusing a workdata changes the results");

  // The routines refuse a work data that is not sized for their streams
  {
    bool                 refused = false;
    polradtran::workdata work(2, 4, 0);  // 8 streams
    Tensor3              r(2, 2, 2, 0.0), t(2, 2, 2, 0.0);
    Matrix               src(2, 2, 0.0);
    Vector               top(2, 0.0), bottom(2, 0.0), up(2), down(2);
    try {
      polradtran::internal_radiance(r, t, src, r, t, src, top, bottom, up, down, work);
    } catch (const std::exception&) { refused = true; }
    if (not refused) throw std::runtime_error("internal_radiance accepted a workdata for 8 streams with 2");
  }

  std::cout << std::format(
      "C++ against Fortran RADTRANO: {} of {} cases within {:.0e} plus the rounding amplified by their doublings "
      "(largest relative difference {:.2e}, at most {:.2f} of a case's tolerance), {} bit-identical\n",
      all.size() - failed,
      all.size(),
      tolerance,
      worst,
      worst_of_tolerance,
      identical);
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
} catch (const std::exception& e) {
  std::cerr << "rt4-radtrano-test failed: " << e.what() << '\n';
  return EXIT_FAILURE;
}
