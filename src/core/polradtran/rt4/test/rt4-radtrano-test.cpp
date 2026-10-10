// The C++ RADTRANO (rt4::radtrano, radtran4.cc) against the outputs of the
// Fortran RADTRANO it ports (radtran4.f), kept as constants
// (rt4-radtrano-reference.h) since the Fortran was removed from ARTS.  Every
// output must agree to tolerance, relative to the largest magnitude in that
// output.
//
// The port is not bitwise: it replaced RT4's quadratures and Planck function
// with ARTS's, folded the MADD, MSUB and MIDENTITY after a product into
// DGEMM's alpha and beta, and inverts with LAPACK instead of LINPACK.
// Mathematically equivalent evaluations are accepted: the tolerances allow
// a few ulp (rounding, 16 epsilon) times what a computation amplifies
// rounding by, so that FMA contraction, vectorised libm functions and other
// BLAS kernels pass and an error of the port does not.  n doublings amplify
// rounding by 2^n, so the radiances must agree to 1e-11 (for the replaced
// quadratures and planck()) plus rounding times 2^n for the layer doubled
// most.  Each case also runs again with one workdata reused over all cases,
// which must not change a bit.
//
// The cases (reference_cases) cover each branch of RADTRANO: the three
// quadratures, extra angles, the four ground types, non-scattering and
// scattering layers (none, one and many doublings, and layers sharing an
// optics set), a 0 K layer top (LINFACTOR = 0) and negative gas extinction
// (clipped in place).  Their inputs are made from mt19937_64, whose output
// is the same on every platform.
#include <arts_constants.h>
#include <radintg.h>
#include <radtran4.h>
#include <radutil4.h>
#include <rt4.h>

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

#include "rt4-radtrano-reference.h"

namespace rt4 = polradtran::rt4;

namespace {
struct inputs {
  Index                       nstokes{}, nummu{}, nuummu{};
  Numeric                     max_delta_tau{};
  polradtran::quadrature_type quad_type{};
  char                        ground_type{};
  Numeric                     ground_temp{}, ground_albedo{}, sky_temp{}, frequency{};
  Complex                     ground_index{};
  Matrix                      ground_reflec;
  Tensor4                     surf_reflect;
  Matrix                      gnd_radiance;
  Vector                      height, temperatures, gas_extinct, scatlayers;
  Tensor5                     extinct_matrix;
  Tensor4                     emis_vector;
  Tensor6                     scatter_matrix;
  Vector                      mu_values;
};

struct outputs {
  Matrix  gnd_radiance;
  Vector  gas_extinct, mu_values;
  Tensor3 up_rad, down_rad;
};

enum class layout { mixed, thin, thick, shared };

struct case_spec {
  Index                       nstokes, nquad, nuummu;
  polradtran::quadrature_type quad;
  char                        ground;
  Index                       nlay;
  layout                      lay;
  Numeric                     max_delta_tau;
  bool                        zero_kelvin_top{false};
};

//! RADTRANO's QUAD_TYPE letter of a quadrature
char quad_letter(polradtran::quadrature_type type) {
  switch (type) {
    case polradtran::quadrature_type::double_gauss: return 'D';
    case polradtran::quadrature_type::gauss:        return 'G';
    case polradtran::quadrature_type::lobatto:      return 'L';
  }
  throw std::runtime_error("unknown quadrature type");
}

std::string describe(const case_spec& c) {
  constexpr const char* names[] = {"mixed", "thin", "thick", "shared"};
  return std::format("nstokes {}, nmu {} + {}, quad '{}', ground '{}', {} layers {}, max_delta_tau {:.0e}{}",
                     c.nstokes,
                     c.nquad,
                     c.nuummu,
                     quad_letter(c.quad),
                     c.ground,
                     c.nlay,
                     names[static_cast<int>(c.lay)],
                     c.max_delta_tau,
                     c.zero_kelvin_top ? ", 0 K top" : "");
}

//! Uniform in [0, 1) from the 53 high bits of the generator's output, the
//! same on every platform (std::uniform_real_distribution is not)
Numeric uniform(std::mt19937_64& gen) { return static_cast<Numeric>(gen() >> 11) * 0x1.0p-53; }

inputs make_inputs(const case_spec& c, std::mt19937_64& gen) {
  const auto u = uniform;

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
    case 'L': return polradtran::lambertian_surface{.albedo = in.ground_albedo};
    case 'F': return polradtran::fresnel_surface{.refractive_index = in.ground_index};
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
  const auto   q = polradtran::get_quadrature(nquad, in.quad_type);
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

//! 2^n for the most doublings n of a scattering layer, as RADTRANO counts
//! them: the amplification of rounding
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

/* The rounding of a few operations, per unit of amplification.  An
   evaluation that is mathematically the same but rounds differently (FMA
   contraction, a vectorised libm exp, another BLAS kernel or summation
   order, LAPACK's inverse for LINPACK's) changes a result by a few ulp
   times what the computation amplifies rounding by. */
constexpr Numeric rounding = 16 * std::numeric_limits<Numeric>::epsilon();

//! The values of a matpack array, row-major
std::span<const Numeric> values(const auto& a) { return {a.data_handle(), a.size()}; }

//! The number of elements that differ, and the largest difference relative
//! to the largest magnitude in either (Q can be near 0 beside I); NaN
//! counts as infinitely far
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

//! The cases whose Fortran RADTRANO outputs are the reference, each made
//! by make_inputs from its own generator, seeded 20261009 + its index
std::vector<case_spec> reference_cases() {
  constexpr auto D = polradtran::quadrature_type::double_gauss, G = polradtran::quadrature_type::gauss,
                 L = polradtran::quadrature_type::lobatto;
  return {
      {1, 3, 0, D, 'L', 3, layout::mixed, 1e-6},        // [I], a Lambertian ground
      {2, 3, 0, G, 'F', 3, layout::mixed, 1e-6},        // a Fresnel ground, Gauss
      {2, 3, 1, L, 'S', 3, layout::mixed, 1e-6},        // a specular ground, Lobatto, an extra angle
      {2, 3, 2, D, 'A', 3, layout::mixed, 1e-6},        // an external ground, two extra angles
      {2, 3, 1, D, 'L', 3, layout::thick, 1e-6},        // many doublings
      {2, 3, 0, D, 'F', 3, layout::thin, 1e-6},         // no doubling
      {2, 3, 0, D, 'A', 4, layout::shared, 1e-6},       // layers sharing sets
      {2, 3, 0, D, 'F', 3, layout::mixed, 1e-3},        // a coarse max_delta_tau
      {2, 3, 1, D, 'L', 3, layout::mixed, 1e-6, true},  // a 0 K top (LINFACTOR = 0)
  };
}
}  // namespace

int main() try {
  // ARTS's quadratures and planck() for RT4's, which differ by 2.3e-12 and
  // 4.5e-13; to this is added the rounding amplified by the doublings
  constexpr Numeric tolerance = 1e-11;

  const auto all = reference_cases();
  if (all.size() != rt4_reference::fortran.size())
    throw std::runtime_error("rt4-radtrano-reference.h has the outputs of another number of cases");

  Index   failed = 0, reuse_differ = 0;
  Numeric worst = 0.0, worst_of_tolerance = 0.0;
  // One work data over all cases, whose sizes differ, against a fresh one
  // per case: reuse must not change a bit
  polradtran::workdata shared;
  for (std::size_t i = 0; i < all.size(); i++) {
    std::mt19937_64      gen(20261009 + i);
    const auto           in = make_inputs(all[i], gen);
    polradtran::workdata fresh;
    const auto           cpp    = run_cpp(in, fresh);
    const auto           reused = run_cpp(in, shared);
    for (auto [count, rel] : {differ(values(cpp.up_rad), values(reused.up_rad)),
                              differ(values(cpp.down_rad), values(reused.down_rad)),
                              differ(values(cpp.gnd_radiance), values(reused.gnd_radiance)),
                              differ(values(cpp.mu_values), values(reused.mu_values))})
      reuse_differ += count;

    const auto&   f77 = rt4_reference::fortran[i];
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
    check("gnd_radiance", cpp.gnd_radiance, f77.gnd_radiance);
    check("gas_extinct", cpp.gas_extinct, f77.gas_extinct);
    check("mu_values", cpp.mu_values, f77.mu_values);
    if (not bad.empty()) {
      failed++;
      std::cout << std::format("Differs by more than {:.2e}: {}:{}\n", case_tolerance, describe(all[i]), bad);
    }
  }

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
      "C++ RADTRANO against the Fortran's outputs: {} of {} cases within {:.0e} plus the rounding amplified by "
      "their doublings (largest relative difference {:.2e}, at most {:.2f} of a case's tolerance)\n",
      all.size() - failed,
      all.size(),
      tolerance,
      worst,
      worst_of_tolerance);
  return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
} catch (const std::exception& e) {
  std::cerr << "rt4-radtrano-test failed: " << e.what() << '\n';
  return EXIT_FAILURE;
}
