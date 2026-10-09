// The C++ RADTRAN (rt3::radtran, radtran3.cc) against the outputs of the
// Fortran RADTRAN it ports (radtran3.f), kept as constants
// (rt3-radtran-reference.h) since the Fortran was removed from ARTS.  Every
// output must agree to tolerance, relative to the largest magnitude in that
// output.
//
// The port is not bitwise: the C++ evaluates in its natural order and lets
// the compiler contract multiply-adds into FMAs, and it replaced RT3's
// quadratures and Planck function with ARTS's, LINPACK's inverse with
// LAPACK's and Evans' matrix products with BLAS.  Mathematically equivalent
// evaluations are accepted: the tolerances allow a few ulp (rounding, 16
// epsilon) times what a computation amplifies rounding by, so that FMA
// contraction, vectorised libm functions and other BLAS kernels pass and an
// error of the port does not.  n doublings amplify rounding by 2^n, so the
// radiances and fluxes must agree to 1e-10 (for the replaced quadratures
// and Planck function) plus rounding times 2^n for the layer doubled most.
// Each case also runs again with one rt3_workdata reused over all cases and
// with one whose every array is NaN, which must not change a bit.
//
// The cases (reference_cases) cover each branch of RADTRAN: 1 to 4 Stokes
// components, the quadratures (with an extra angle, QUAD_TYPE 'E'),
// azimuth orders, delta-M, the source codes (none, solar, thermal, both),
// the Lambertian and Fresnel grounds, non-scattering and scattering layers
// (none, one and many doublings, and layers sharing a set), Rayleigh-, Mie-
// and general phase matrices (each of RT3's summation cases), negative gas
// extinction and output levels in any order.  Their inputs are made from
// mt19937_64, whose output is the same on every platform.
//
// It also checks fft1dr's documented format against direct sums (what an
// FFTW build must also give), up to 4096 values, beyond Evans' 512, that
// SCATTERING and DIRECT_SCATTERING give the same m = 0 mode through it
// for 1024 and 4096 azimuths as without it, and that makephase refuses a
// table that is not a power of two.
#include <arts_constants.h>
#include <radintg.h>
#include <radscat3.h>
#include <radtran3.h>
#include <radutil3.h>
#include <rt3.h>
#include <rt3_fft.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <format>
#include <iostream>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "rt3-radtran-reference.h"

namespace rt3 = polradtran::rt3;

namespace {
using polradtran::quadrature_type;

constexpr Numeric nan = std::numeric_limits<Numeric>::quiet_NaN();

//! The values of a matpack array, row-major
std::span<const Numeric> values(const auto& a) { return {a.data_handle(), a.size()}; }

//! The number of elements that differ, and the largest difference relative
//! to the largest magnitude in either; NaN counts as infinitely far
std::pair<Index, Numeric> differ(std::span<const Numeric> a, std::span<const Numeric> b) {
  if (a.size() != b.size())
    return {static_cast<Index>(std::max(a.size(), b.size())), std::numeric_limits<Numeric>::infinity()};
  Index   count = 0;
  Numeric diff = 0.0, scale = 0.0;
  for (std::size_t k = 0; k < a.size(); k++) {
    scale = std::max({scale, std::abs(a[k]), std::abs(b[k])});
    if (a[k] == b[k]) continue;
    count++;
    const Numeric d = std::abs(a[k] - b[k]);
    diff            = std::max(diff, std::isnan(d) ? std::numeric_limits<Numeric>::infinity() : d);
  }
  return {count, count == 0 ? 0.0 : diff / scale};
}

/* The rounding of a few operations, per unit of amplification.  An
   evaluation that is mathematically the same but rounds differently (FMA
   contraction, a vectorised libm exp, another BLAS kernel or summation
   order, LAPACK's inverse for LINPACK's) changes a result by a few ulp
   times what the computation amplifies rounding by. */
constexpr Numeric rounding = 16 * std::numeric_limits<Numeric>::epsilon();

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

//! Uniform in [0, 1) from the 53 high bits of the generator's output, the
//! same on every platform (std::uniform_real_distribution is not)
Numeric uniform(std::mt19937_64& gen) { return static_cast<Numeric>(gen() >> 11) * 0x1.0p-53; }

inputs make_inputs(const case_spec& c, std::mt19937_64& gen) {
  const auto u = uniform;

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

//! 2^n for the most doublings n of a scattering layer, as RADTRAN counts
//! them (delta-M, which only lowers the extinction, aside): the
//! amplification of rounding
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

//! The C++ RADTRAN, with its ground made by rt3::ground_surface; all of
//! mu_values is output.  It works in SI at the frequency, the Fortran per
//! micrometre at the wavelength: its inputs and outputs are converted
//! (B_nu = B_lambda[um^-1] lambda[um] / f)
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
  w.scat_matrix.resize(1025);
  w.basis_matrix.resize(1025);
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
  for (Tensor3* t : {&w.source, &w.reflect1, &w.upreflect, &w.downreflect, &w.trans1, &w.uptrans, &w.downtrans})
    *t = nan;
  w.reflect        = nan;
  w.trans          = nan;
  w.scatter_matrix = nan;
  for (MuelmatVector* v : {&w.scat_matrix, &w.basis_matrix}) *v = Muelmat::constant(nan);
  w.scatbuf   = Muelmat::constant(nan);
  w.directbuf = Stokvec{nan, nan, nan, nan};
  for (Index& s : w.scat_nums) s = -1;
  return w;
}

/* fft1dr's format, which an FFTW build must also give: the forward
   transform against the direct sums X_k = sum_j x_j exp(+2 pi i j k / n),
   packed [X_0, X_{n/2}, Re X_1, Im X_1, ...], and the inverse against
   x_j = X_0 + (-1)^j X_{n/2} + 2 sum_k Re(X_k exp(-2 pi i j k / n)), to
   1e-12 of the largest value. */
void check_fft1dr_format() {
  std::mt19937_64   gen(2003);
  const auto        u     = uniform;
  Numeric           worst = 0.0;
  rt3::fft_workdata fft;
  for (Index n : {2, 4, 8, 64, 512, 1024, 4096}) {
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
    worst = std::max(worst, differ(values(forward), values(packed)).second);

    Vector direct(n, 0.0);
    for (Index j = 0; j < n; j++) {
      direct[j] = packed[0] + (j % 2 == 0 ? 1.0 : -1.0) * packed[1];
      for (Index k = 1; k < n / 2; k++)
        direct[j] += 2.0 * (packed[2 * k] * std::cos(angle(j, k)) + packed[2 * k + 1] * std::sin(angle(j, k)));
    }
    Vector inverse = packed;
    rt3::fft1dr(inverse, rt3::fft_direction::inverse, fft);
    worst = std::max(worst, differ(values(inverse), values(direct)).second);
  }
  std::cout << std::format("fft1dr against the direct sums of its format: within {:.2e} of the largest value\n", worst);
  if (worst > 1e-12) throw std::runtime_error("fft1dr does not give its documented format");
}

/* SCATTERING and DIRECT_SCATTERING above the 512 azimuth samples of
   Evans' FFT1DR, which the port does not have: with aziorder > 0 they
   sample 2 * 2^int(log2(degree + 4) + 1) azimuths (1024 for degree 300,
   4096 for degree 1100, also above Evans' 1023) and transform them with
   fft1dr; with aziorder 0 they average 2 int((degree + 1) / 2) + 4 samples
   without it.  P11 is a polynomial of the degree in cos(phi), so both
   sample its m = 0 mode without aliasing and must agree to 1e-12 of the
   largest value.  (The polarized elements also depend on the rotation
   angles, which are not band-limited: RT3 samples them with aliasing.) */
void check_large_fft() {
  const Vector mu{0.15, 0.55, 0.95}, w{0.3, 0.4, 0.3};
  Numeric      worst = 0.0;
  for (Index degree : {300, 1100}) {
    const Matrix      coef = legendre_set(2, 0.6, degree);
    rt3::rt3_workdata work;
    MuelmatTensor4    s0(1, 2, 3, 3), s2(3, 2, 3, 3);
    rt3::scattering(mu, w, coef, 4, s0, work);
    rt3::scattering(mu, w, coef, 4, s2, work);
    Tensor3 p0(2, 3, 3), p2(2, 3, 3);
    for (Size i = 0; i < p0.size(); i++) {
      p0.elem_at(i) = s0[0].elem_at(i)[0, 0];
      p2.elem_at(i) = s2[0].elem_at(i)[0, 0];
    }
    worst = std::max(worst, differ(values(p0), values(p2)).second);
    StokvecTensor3 d0(1, 2, 3), d2(3, 2, 3);
    rt3::direct_scattering(mu, coef, 0.6, 4, d0, work);
    rt3::direct_scattering(mu, coef, 0.6, 4, d2, work);
    Matrix q0(2, 3), q2(2, 3);
    for (Size i = 0; i < q0.size(); i++) {
      q0.elem_at(i) = d0[0].elem_at(i).I();
      q2.elem_at(i) = d2[0].elem_at(i).I();
    }
    worst = std::max(worst, differ(values(q0), values(q2)).second);
  }
  std::cout << std::format(
      "SCATTERING and DIRECT_SCATTERING with 1024 and 4096 azimuths: the m = 0 P11 through fft1dr within {:.2e} "
      "of the mean\n",
      worst);
  if (not(worst <= 1e-12)) throw std::runtime_error("the m = 0 mode through fft1dr differs from the mean");
}

//! makephase's table fits in 4 nmax values only for nmax a power of two,
//! as FFT1DR uses it; for another nmax the Fortran wrote past it
void check_makephase_limits() {
  for (Index nmax : {3, 5, 12}) {
    bool threw = false;
    try {
      Vector phase(4 * nmax);
      rt3::makephase(phase);
    } catch (const std::exception&) { threw = true; }
    if (not threw) throw std::runtime_error(std::format("makephase with nmax {} did not throw", nmax));
  }
}

//! The cases whose Fortran RADTRAN outputs are the reference, each made
//! by make_inputs from its own generator, seeded 20261009 + its index
std::vector<case_spec> reference_cases() {
  constexpr auto G = quadrature_type::gauss, D = quadrature_type::double_gauss, L = quadrature_type::lobatto;
  return {
      {1, 2, 0, G, 0, 2, false, 'L', 3, layout::mixed, 1e-6},   // [I], the thermal source
      {2, 3, 0, D, 1, 3, false, 'L', 3, layout::mixed, 1e-6},   // both sources
      {3, 3, 0, L, 2, 1, false, 'L', 3, layout::mixed, 1e-6},   // the solar source, Lobatto
      {4, 2, 1, G, 2, 3, false, 'L', 3, layout::mixed, 1e-6},   // an extra angle (QUAD_TYPE 'E')
      {4, 3, 0, G, 1, 2, false, 'F', 3, layout::mixed, 1e-6},   // a Fresnel ground
      {4, 3, 0, G, 0, 0, false, 'L', 3, layout::mixed, 1e-6},   // no source but the sky
      {4, 3, 0, G, 1, 3, true, 'L', 3, layout::mixed, 1e-6},    // delta-M
      {2, 3, 1, G, 2, 3, false, 'L', 3, layout::thick, 1e-6},   // many doublings
      {2, 3, 0, G, 1, 3, false, 'L', 3, layout::thin, 1e-6},    // no doubling
      {3, 2, 0, G, 1, 3, false, 'L', 4, layout::shared, 1e-6},  // layers sharing sets
      {2, 3, 0, D, 1, 3, false, 'L', 3, layout::gas, 1e-6},     // gas only
      {2, 3, 0, D, 2, 3, false, 'L', 3, layout::mixed, 1e-3},   // a coarse max_delta_tau
  };
}
}  // namespace

int main() try {
  check_fft1dr_format();
  check_makephase_limits();
  check_large_fft();

  // ARTS's quadratures and Planck function for RT3's; to this is added the
  // rounding amplified by the doublings
  constexpr Numeric tolerance = 1e-10;

  const auto all = reference_cases();
  if (all.size() != rt3_reference::fortran.size())
    throw std::runtime_error("rt3-radtran-reference.h has the outputs of another number of cases");

  Index   failed = 0, reuse_differ = 0;
  Numeric worst = 0.0, worst_of_tolerance = 0.0;
  // One work data over all cases, whose sizes differ, and one with every
  // array NaN, against a fresh one per case: neither may change a bit
  rt3::rt3_workdata shared;
  for (std::size_t i = 0; i < all.size(); i++) {
    std::mt19937_64   gen(20261009 + i);
    const auto        in = make_inputs(all[i], gen);
    rt3::rt3_workdata fresh, poisoned = poisoned_workdata(in);
    const auto        cpp = run_cpp(in, fresh);
    for (rt3::rt3_workdata* work : {&shared, &poisoned}) {
      const auto again  = run_cpp(in, *work);
      reuse_differ     += differ(values(again.up_rad), values(cpp.up_rad)).first +
                          differ(values(again.down_rad), values(cpp.down_rad)).first +
                          differ(values(again.up_flux), values(cpp.up_flux)).first +
                          differ(values(again.down_flux), values(cpp.down_flux)).first +
                          differ(values(again.mu_values), values(cpp.mu_values)).first;
    }

    const auto&   f77 = rt3_reference::fortran[i];
    std::string   bad;
    const Numeric case_tolerance = tolerance + rounding * doubling_amplification(in);
    const auto    check          = [&](const char* name, const auto& a, const std::vector<Numeric>& b) {
      const auto [count, rel] = differ(values(a), b);
      worst                   = std::max(worst, rel);
      worst_of_tolerance      = std::max(worst_of_tolerance, rel / case_tolerance);
      if (not(rel <= case_tolerance)) bad += std::format(" {}: {} differ, max rel {:.3e};", name, count, rel);
    };
    check("up_rad", cpp.up_rad, f77.up_rad);
    check("down_rad", cpp.down_rad, f77.down_rad);
    check("up_flux", cpp.up_flux, f77.up_flux);
    check("down_flux", cpp.down_flux, f77.down_flux);
    check("mu_values", cpp.mu_values, f77.mu_values);
    if (not bad.empty()) {
      failed++;
      std::cout << std::format("Differs by more than {:.2e}: {}:{}\n", case_tolerance, describe(all[i]), bad);
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
    std::mt19937_64 gen(1);
    auto            in = make_inputs({2, 6, 0, quadrature_type::gauss, 0, 3, false, 'F', 3, layout::mixed, 1e-6}, gen);
    rt3::rt3_workdata work;
    run_cpp(in, work);
  } catch (const std::exception&) { threw = true; }
  if (not threw) throw std::runtime_error("a solar source over a Fresnel ground did not throw");

  std::cout << std::format(
      "C++ RADTRAN against the Fortran's outputs: {} of {} cases within {:.0e} plus the rounding amplified by "
      "their doublings (largest relative difference {:.2e}, at most {:.2f} of a case's tolerance)\n",
      all.size() - failed,
      all.size(),
      tolerance,
      worst,
      worst_of_tolerance);
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
} catch (const std::exception& e) {
  std::cerr << "rt3-radtran-test failed: " << e.what() << '\n';
  return EXIT_FAILURE;
}
