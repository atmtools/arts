/* RT4, RT3 and VDISORT on one ARTS atmosphere, through their ARTS-native
   path builders (rt4::problem_from_path, rt3::problem_from_path and
   vdisort::main_data_from_path).

   Where the expected answer comes from.  RT4 and RT3 are Evans'
   doubling-adding solvers; VDISORT is a discrete-ordinate eigen-solver
   (see vdisort-rt4-comparison.cpp and vdisort-rt3-comparison.cpp for why
   that makes them external references for each other).  Here the three
   solvers get their inputs from the same ARTS data by three different
   routes:
   - RT4: the azimuthal mean of ARTS's laboratory-frame phase matrix
     (get_bulk_scattering_properties_aro_fourier, Mishchenko-style
     rotations) on RT4's streams;
   - RT3: the species' Legendre series, which RT3 rotates and
     Fourier-transforms itself;
   - VDISORT: the Fourier modes of the same laboratory-frame phase matrix on
     VDISORT's streams, in VDISORT's vector-geometry basis.
   With the same double-Gauss streams, the solvers then solve the same
   problem, and RT3 and RT4 differ from VDISORT by their first-order
   doubling error, below 10 max_delta_tau / mu0 of max I (mu0 = 1 without a
   beam).  A convention error in a builder (hemisphere, stream or Stokes
   order, sign of Q, U or V, azimuth sense, normalisation, layer averaging,
   source or surface) would not.

   The routes evaluate the same scattering matrices: the GasScatterer's
   closed form, and the cloud's Legendre series (exact from its tabulated
   matrix on the nodes of a Gauss-Legendre rule).  test_thermal() checks that RT4's layer phase matrices equal
   VDISORT's to round-off.

   The atmosphere: an AtmField with exponential pressure, a temperature
   profile from 290 K at the surface to 230 K at 6 km, Rayleigh scattering
   by a GasScatterer with a constant cross-section, and a cloud of 1.5 mm
   liquid spheres (ARTS's Mie code with Ellison's water permittivity,
   size parameter 1.4 at 89 GHz, F12 and F34 non-zero) with a Gaussian
   number-density profile, plus an unpolarized gas absorption profile. */
#include <arts_constants.h>
#include <arts_conversions.h>
#include <atm_path.h>
#include <legendre.h>
#include <physics_funcs.h>
#include <rt3_arts.h>
#include <rt4_arts.h>
#include <vdisort_arts.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string>

namespace rt3 = polradtran::rt3;
namespace rt4 = polradtran::rt4;

namespace {
constexpr Numeric frequency     = 89e9;
constexpr Index   nmu           = 8;
constexpr Index   native_angles = 4000;  // the cloud's scattering-angle grid, the nodes of its projection rule
constexpr Index   cloud_degree  = 64;    // the degree of the cloud's Legendre series, converged to rounding

void require(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error(what);
}

const ScatteringSpeciesProperty drops{"drops", ParticulateProperty::NumberDensity};

struct atmosphere {
  ArrayOfPropagationPathPoint ray_path;
  ArrayOfAtmPoint             atm_path;
  ArrayOfPropmatVector        propmat;
  AscendingGrid               freq_grid{Vector{frequency}};
  ArrayOfScatteringSpecies    species;
};

atmosphere make_atmosphere() {
  AtmField field;
  field.top_of_atmosphere = 6000.0;
  field[AtmKey::p] = Atm::FunctionalData{[](Numeric alt, Numeric, Numeric) { return 1e5 * std::exp(-alt / 8000.0); }};
  field[AtmKey::t] = Atm::FunctionalData{[](Numeric alt, Numeric, Numeric) { return 290.0 - 0.01 * alt; }};
  field[drops]     = Atm::FunctionalData{
      [](Numeric alt, Numeric, Numeric) { return 30.0 * std::exp(-std::pow((alt - 2500.0) / 1500.0, 2)); }};

  atmosphere a;
  for (Numeric alt : {6000.0, 4500.0, 3200.0, 2000.0, 1000.0, 0.0}) {
    PropagationPathPoint p;
    p.pos = {alt, 0.0, 0.0};
    p.los = {180.0, 0.0};
    a.ray_path.push_back(p);
    PropmatVector pm(1);
    pm[0].A() = 5e-5 * std::exp(-alt / 2000.0);
    a.propmat.push_back(pm);
  }
  a.atm_path = forward_atm_path(a.ray_path, field);

  a.species.add(GasScatterer{scattering::ConstantGasScattering{2e-31}, scattering::RayleighGasScattering{0.0}});

  // The cloud's native grid is the Gauss-Legendre rule of the projection, in ascending scattering angle
  Vector x(native_angles), w(native_angles), angles(native_angles);
  Legendre::GaussLegendre(x, w);
  stdr::reverse(x);
  for (Index i = 0; i < native_angles; i++) angles[i] = Conversion::rad2deg(std::acos(x[i]));
  Vector     t_grid{220.0, 250.0, 280.0, 310.0}, f_grid{frequency}, diameter{1.5e-3};
  const auto habit = ParticleHabit::liquid_sphere(
      t_grid, f_grid, diameter, scattering::ZenithAngleGrid{scattering::IrregularZenithAngleGrid(angles)});

  // The Legendre series, a_l = 2 pi sqrt((2 l + 1) / 4 pi) sum_i w_i F(x_i) P_l(x_i), exact for data on the nodes
  using Gridded =
      scattering::SingleScatteringData<Numeric, scattering::Format::TRO, scattering::Representation::Gridded>;
  using Spectral =
      scattering::SingleScatteringData<Numeric, scattering::Format::TRO, scattering::Representation::Spectral>;
  const auto& gridded = std::get<Gridded>(habit[0]);
  scattering::PhaseMatrixData<Numeric, scattering::Format::TRO, scattering::Representation::Spectral> series(
      gridded.phase_matrix->get_t_grid(), gridded.phase_matrix->get_f_grid(), cloud_degree);
  Vector p(cloud_degree + 1);
  for (Index i = 0; i < native_angles; i++) {
    Legendre::legendre_polynomials(p, x[i]);
    for (Index l = 0; l <= cloud_degree; l++) {
      // The Gauss-Legendre weights are symmetric, so w is that of the reversed x too
      const Numeric f =
          2.0 * Constant::pi * std::sqrt((2.0 * static_cast<Numeric>(l) + 1.0) / (4.0 * Constant::pi)) * w[i] * p[l];
      for (Index it = 0; it < series.extent(0); it++)
        for (Index k = 0; k < 6; k++) series[it, 0, l, k] += f * (*gridded.phase_matrix)[it, 0, i, k];
    }
  }
  const Spectral spectral(gridded.properties,
                          series,
                          gridded.extinction_matrix.to_spectral(),
                          gridded.absorption_vector.to_spectral(),
                          gridded.backscatter_matrix,
                          gridded.forwardscatter_matrix);
  a.species.add(ScatteringHabit{
      ParticleHabit{std::vector<Spectral>{spectral}}, scattering::PSD{scattering::MonodispersePSD{drops}}, 1.0, 3.0});
  return a;
}

struct radiance {
  Tensor4 up;    // [level, azimuth, stream, Stokes]
  Tensor4 down;  // [level, azimuth, stream, Stokes]
};

//! max |a - b| over levels, azimuths, streams and both directions, per Stokes, relative to max |I_b|
Vector4 deviation(const radiance& a, const radiance& b) {
  Numeric scale = 0.0;
  for (const auto* t : {&b.up, &b.down})
    for (Index x = 0; x < t->extent(0); x++)
      for (Index k = 0; k < t->extent(1); k++)
        for (Index i = 0; i < t->extent(2); i++) scale = std::max(scale, std::abs((*t)[x, k, i, 0]));
  Vector4 d{};
  for (Index l = 0; l < a.up.extent(0); l++)
    for (Index k = 0; k < a.up.extent(1); k++)
      for (Index i = 0; i < a.up.extent(2); i++)
        for (Index s = 0; s < a.up.extent(3); s++)
          d[s] = std::max({d[s],
                           std::abs(a.up[l, k, i, s] - b.up[l, k, i, s]) / scale,
                           std::abs(a.down[l, k, i, s] - b.down[l, k, i, s]) / scale});
  return d;
}

//! VDISORT at every level and azimuth phi0 + psi, the first ns Stokes components
radiance vdisort_radiance(const vdisort::main_data& v, const Vector& psi, Numeric phi0, Index ns) {
  const Index N = static_cast<Index>(v.weights().size()), nlev = static_cast<Index>(v.tau().size()) + 1;
  const Index na = static_cast<Index>(psi.size());
  radiance    r{.up = Tensor4(nlev, na, N, ns), .down = Tensor4(nlev, na, N, ns)};
  for (Index l = 0; l < nlev; l++) {
    for (Index k = 0; k < na; k++) {
      vdisort::u_data u;
      Numeric         phi = std::fmod(phi0 + psi[k], 2 * Constant::pi);
      if (phi < 0) phi += 2 * Constant::pi;
      v.u(u, l == 0 ? 0.0 : v.tau()[l - 1], phi);
      for (Index i = 0; i < N; i++) {
        for (Index s = 0; s < ns; s++) {
          r.up[l, k, i, s]   = u.intensities[i][s];
          r.down[l, k, i, s] = u.intensities[N + i][s];
        }
      }
    }
  }
  return r;
}

radiance rt3_radiance(const rt3::result& res, const Vector& psi) {
  return {.up = rt3::azimuth_radiance(res.up, psi), .down = rt3::azimuth_radiance(res.down, psi)};
}

radiance rt4_radiance(const rt4::result& res) {
  const Index nlev = res.up.extent(0), n = res.up.extent(1), ns = res.up.extent(2);
  radiance    r{.up = Tensor4(nlev, 1, n, ns), .down = Tensor4(nlev, 1, n, ns)};
  for (Index l = 0; l < nlev; l++) {
    for (Index i = 0; i < n; i++) {
      for (Index s = 0; s < ns; s++) {
        r.up[l, 0, i, s]   = res.up[l, i, s];
        r.down[l, 0, i, s] = res.down[l, i, s];
      }
    }
  }
  return r;
}

Numeric report(std::string_view what, const Vector4& d, Index ns, Numeric tol, Numeric scale) {
  Numeric worst = 0.0;
  for (Index s = 0; s < ns; s++) worst = std::max(worst, d[s]);
  std::cout << std::format("{:<48} I {:9.3e}  Q {:9.3e}", what, d[0], d[1]);
  if (ns > 2) std::cout << std::format("  U {:9.3e}  V {:9.3e}", d[2], d[3]);
  std::cout << std::format("  = {:5.2f} max_delta_tau / mu0; tolerance {:.2e}\n", worst / scale, tol);
  require(worst <= tol, std::format("{}: the solvers must agree to {:.2e} of max I, got {:.3e}", what, tol, worst));
  return worst;
}

/* The largest relative difference between RT4's layer phase matrices and
   the [I, Q] block of sigma / (4 pi) C^0 from vdisort::scattering_optics,
   both layer means, on the shared streams: the input difference of the RT4
   route. */
Numeric rt4_input_difference(const atmosphere& a, const rt4::problem& p) {
  Vector mu(2 * nmu), inv(2 * nmu), w(nmu);
  disort_common::initialize_streams(mu, inv, w);
  Numeric worst = 0.0;
  for (Index l = 0; l < static_cast<Index>(p.layer_optics_index.size()); l++) {
    const auto&                            o = p.optics[p.layer_optics_index[l]];
    std::array<vdisort::fourier_optics, 2> v;
    for (Index j = 0; j < 2; j++)
      v[j] = vdisort::scattering_optics(a.species, a.atm_path[l + j], frequency, mu, mu, 1, 1e-6);
    Numeric scale = 0.0, diff = 0.0;
    for (Index ho = 0; ho < 2; ho++) {
      for (Index hi = 0; hi < 2; hi++) {
        for (Index io = 0; io < nmu; io++) {
          for (Index ii = 0; ii < nmu; ii++) {
            const Index so = ho == rt4::up ? io : nmu + io, si = hi == rt4::up ? ii : nmu + ii;
            for (Index s = 0; s < 2; s++) {
              for (Index t = 0; t < 2; t++) {
                Numeric ref = 0.0;
                for (Index j = 0; j < 2; j++)
                  ref += 0.5 * v[j].scattering / (4 * Constant::pi) * v[j].cosine[0, so, si][s, t];
                scale = std::max(scale, std::abs(ref));
                diff  = std::max(diff, std::abs(o.phase[ho, hi, io, ii, s, t] - ref));
              }
            }
          }
        }
      }
    }
    worst = std::max(worst, diff / scale);
  }
  return worst;
}

struct thermal_result {
  Vector4 rt3_vdisort, rt4_vdisort, rt4_rt3;
  Numeric input;  // the relative difference of the RT4 layer phase matrices from VDISORT's
  Numeric q;      // max |Q_up| / max |I_up|
};

thermal_result thermal(const atmosphere& a, Numeric max_delta_tau, bool fresnel) {
  const rt4::surface     g4 = fresnel ? rt4::surface{polradtran::fresnel_surface{.refractive_index = Complex{3.0, 0.2}}}
                                      : rt4::surface{polradtran::lambertian_surface{.albedo = 0.3}};
  const rt3::surface     g3 = fresnel ? rt3::surface{polradtran::fresnel_surface{.refractive_index = Complex{3.0, 0.2}}}
                                      : rt3::surface{polradtran::lambertian_surface{.albedo = 0.3}};
  const vdisort::surface gv = fresnel
                                  ? vdisort::surface{vdisort::fresnel_surface{.refractive_index = Complex{3.0, 0.2}}}
                                  : vdisort::surface{vdisort::lambertian_surface{.albedo = 0.3}};

  const auto p4 = rt4::problem_from_path(
      a.ray_path,
      a.atm_path,
      a.propmat,
      a.freq_grid,
      0,
      a.species,
      {.nstokes = 2, .nmu = nmu, .quad = polradtran::quadrature_type::double_gauss, .max_delta_tau = max_delta_tau},
      g4,
      288.0,
      2.7);
  const auto p3 = rt3::problem_from_path(a.ray_path,
                                         a.atm_path,
                                         a.propmat,
                                         a.freq_grid,
                                         0,
                                         a.species,
                                         {.nstokes                 = 2,
                                          .nmu                     = nmu,
                                          .quad                    = polradtran::quadrature_type::double_gauss,
                                          .max_delta_tau           = max_delta_tau,
                                          .normalisation_tolerance = 1e-6},
                                         g3,
                                         288.0,
                                         2.7);
  const auto v  = vdisort::main_data_from_path(a.ray_path,
                                               a.atm_path,
                                               a.propmat,
                                               a.freq_grid,
                                               0,
                                               a.species,
                                               {.nquad = 2 * nmu, .nfourier = 1, .normalisation_tolerance = 1e-6},
                                               gv,
                                               288.0,
                                               2.7);

  const Vector psi{0.0};
  const auto   rv = vdisort_radiance(v, psi, 0.0, 2);
  const auto   r3 = rt3_radiance(rt3::solve(p3), psi);
  const auto   r4 = rt4_radiance(rt4::solve(p4));

  Numeric qmax = 0.0, imax = 0.0;
  for (Index x = 0; x < static_cast<Index>(rv.up.size()) / 2; x++) {
    imax = std::max(imax, std::abs(rv.up.data_handle()[2 * x]));
    qmax = std::max(qmax, std::abs(rv.up.data_handle()[2 * x + 1]));
  }
  return {.rt3_vdisort = deviation(r3, rv),
          .rt4_vdisort = deviation(r4, rv),
          .rt4_rt3     = deviation(r4, r3),
          .input       = rt4_input_difference(a, p4),
          .q           = qmax / imax};
}

/* Thermal emission, Lambertian and Fresnel surfaces, nstokes 2.
   - RT4's layer phase matrices must equal VDISORT's to round-off (1e-12
     relative).
   - RT3 and RT4 run Evans' identical doubling scheme on the same problem,
     so their doubling errors cancel and RT4 - RT3 must be round-off (1e-8
     of max I).
   - RT3 and RT4 vs VDISORT: 10 max_delta_tau.
   - The RT3 - VDISORT difference must fall with max_delta_tau like the
     first-order doubling error (a factor 5 to 20 per decade). */
void test_thermal(const atmosphere& a) {
  constexpr Numeric mdt = 1e-7;
  const Numeric     tol = 10 * mdt;

  for (bool fresnel : {false, true}) {
    const std::string surf = fresnel ? "Fresnel 3+0.2i" : "Lambertian 0.3";
    const auto        r    = thermal(a, mdt, fresnel);
    std::cout << std::format(
        "Thermal, {}, max_delta_tau {:.0e}: RT4 layer phase matrices differ from VDISORT's by {:.2e} (relative); "
        "max |Q_up| / max |I_up| = {:.2e}\n",
        surf,
        mdt,
        r.input,
        r.q);
    require(r.input <= 1e-12,
            std::format("RT4's layer phase matrices must equal VDISORT's to 1e-12 (relative), got {:.2e}", r.input));
    report("    RT3 vs VDISORT", r.rt3_vdisort, 2, tol, mdt);
    report("    RT4 vs RT3 (identical doubling)", r.rt4_rt3, 2, 1e-8, mdt);
    report("    RT4 vs VDISORT", r.rt4_vdisort, 2, tol, mdt);
    require(r.q > 1e-3, "The thermal comparison must exercise Q");
  }

  const auto    coarse = thermal(a, 1e-6, false);
  const auto    fine   = thermal(a, 1e-7, false);
  const Numeric ratio  = coarse.rt3_vdisort[0] / fine.rt3_vdisort[0];
  std::cout << std::format("Thermal, Lambertian: RT3 - VDISORT ratio for max_delta_tau 1e-6 / 1e-7: {:.2f}\n", ratio);
  require(ratio > 5.0 and ratio < 20.0,
          std::format("The RT3 - VDISORT difference must fall like RT3's first-order doubling error, by a factor 5 to "
                      "20 per decade of max_delta_tau; got {:.2f}",
                      ratio));
}

void test_solar(const atmosphere& a, Numeric max_delta_tau) {
  constexpr Numeric mu0 = 0.6, flux = 1e-15;  // flux on the horizontal, of the order of pi B at 89 GHz
  constexpr Index   aziorder = 7;
  const auto        p3       = [&] {
    auto p        = rt3::problem_from_path(a.ray_path,
                                           a.atm_path,
                                           a.propmat,
                                           a.freq_grid,
                                           0,
                                           a.species,
                                           {.nstokes                 = 4,
                                            .nmu                     = nmu,
                                            .quad                    = polradtran::quadrature_type::double_gauss,
                                            .aziorder                = aziorder,
                                            .max_delta_tau           = max_delta_tau,
                                            .normalisation_tolerance = 1e-6},
                                           polradtran::lambertian_surface{.albedo = 0.3},
                                           288.0,
                                           2.7);
    p.direct_flux = flux;
    p.direct_mu   = mu0;
    return p;
  }();
  const auto v = vdisort::main_data_from_path(
      a.ray_path,
      a.atm_path,
      a.propmat,
      a.freq_grid,
      0,
      a.species,
      {.nquad = 2 * nmu, .nfourier = aziorder + 1, .normalisation_tolerance = 1e-6, .beam_flux = flux, .beam_mu = mu0},
      vdisort::lambertian_surface{.albedo = 0.3},
      288.0,
      2.7);

  Vector psi(6);
  for (Index k = 0; k < 6; k++) psi[k] = Conversion::deg2rad(std::array{0.0, 30.0, 75.0, 135.0, 180.0, 250.0}[k]);
  const auto rv = vdisort_radiance(v, psi, 0.0, 4);
  const auto r3 = rt3_radiance(rt3::solve(p3), psi);

  const auto d = deviation(r3, rv);
  std::cout << std::format("Solar (mu0 {}) and thermal, Lambertian 0.3, {} Fourier modes, max_delta_tau {:.0e}\n",
                           mu0,
                           aziorder + 1,
                           max_delta_tau);
  report("    RT3 vs VDISORT, 6 azimuths", d, 4, 10 * max_delta_tau / mu0, max_delta_tau / mu0);

  Vector4 top{};
  for (Index x = 0; x < static_cast<Index>(rv.up.size()) / 4; x++)
    for (Index s = 0; s < 4; s++) top[s] = std::max(top[s], std::abs(rv.up.data_handle()[4 * x + s]));
  std::cout << std::format("    max |Q|, |U|, |V| / max |I| (upward) = {:.2e}, {:.2e}, {:.2e}\n",
                           top[1] / top[0],
                           top[2] / top[0],
                           top[3] / top[0]);
  require(top[1] > 1e-3 * top[0] and top[2] > 1e-3 * top[0] and top[3] > 1e-6 * top[0],
          "The solar comparison must exercise Q, U and V");
}
}  // namespace

int main() try {
  require(rt3::available() and rt4::available(), "This test requires ENABLE_RT3=ON and ENABLE_RT4=ON");
  const auto a = make_atmosphere();
  test_thermal(a);
  test_solar(a, 1e-7);
  std::cout << "vdisort-arts comparison passed\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
