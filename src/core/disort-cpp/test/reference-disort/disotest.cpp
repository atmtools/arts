#include <disort-brdf.h>
#include <legendre.h>

#include <limits>

#include "../reference-data.h"
#include "../test-helpers.h"

namespace legacy_disotest {
namespace {
void expect_reference(const std::string_view name,
                      const Numeric          actual,
                      const Numeric          expected,
                      const Numeric          relative_tolerance,
                      const Numeric          absolute_tolerance = 0.0) {
  if (std::abs(actual - expected) > absolute_tolerance + relative_tolerance * std::max(1.0, std::abs(expected)))
    arts::user_error("{}: expected {}, got {} (difference {})", name, expected, actual, actual - expected);
}

void expect_small_reference(const std::string_view name,
                            const Numeric          actual,
                            const Numeric          expected,
                            const Numeric          relative_tolerance = 2e-3,
                            const Numeric          absolute_tolerance = 1e-10) {
  if (std::abs(actual - expected) > absolute_tolerance + relative_tolerance * std::abs(expected))
    arts::user_error("{}: expected {}, got {} (difference {})", name, expected, actual, actual - expected);
}

Numeric blackbody_radiance(const Numeric temperature) {
  return Constant::sigma * std::pow(temperature, 4) * Constant::inv_pi;
}

Numeric band_blackbody_radiance(const Numeric temperature,
                                const Numeric wavenumber_low,
                                const Numeric wavenumber_high) {
  if (temperature == 0.0 || wavenumber_high == wavenumber_low) return 0.0;
  constexpr Index intervals = 4096;
  const Numeric   scale     = Constant::h * Constant::c * 100.0 / (Constant::k * temperature);
  const Numeric   x0        = scale * wavenumber_low;
  const Numeric   x1        = scale * wavenumber_high;
  const Numeric   dx        = (x1 - x0) / intervals;
  const auto      integrand = [](const Numeric x) {
    if (x == 0.0) return 0.0;
    if (x > 700.0) return 0.0;
    return x * x * x / std::expm1(x);
  };
  Numeric integral = integrand(x0) + integrand(x1);
  for (Index i = 1; i < intervals; ++i) integral += (i % 2 ? 4.0 : 2.0) * integrand(x0 + static_cast<Numeric>(i) * dx);
  integral *= dx / 3.0;
  return blackbody_radiance(temperature) * integral / (Math::pow2(Constant::pi) * Math::pow2(Constant::pi) / 15.0);
}

Matrix henyey_greenstein(const Vector& g, const Index nmom) {
  Matrix out(g.size(), nmom + 1, 1.0);
  for (Index l = 1; l <= nmom; ++l) {
    for (Size i = 0; i < g.size(); ++i) out[i, l] = std::pow(g[i], l);
  }
  return out;
}

Matrix linear_polynomial(const AscendingGrid& tau, const Vector& values) {
  if (values.size() != tau.size() + 1) arts::user_error("A value is required at every layer interface");

  Matrix  out(tau.size(), 2);
  Numeric tau0 = 0;
  for (Size i = 0; i < tau.size(); ++i) {
    const Numeric b0    = values[i];
    const Numeric b1    = values[i + 1];
    const Numeric slope = (b1 - b0) / (tau[i] - tau0);
    // CPP-DISORT applies the LTE absorption factor (1 - omega).  Store the
    // Planck radiance polynomial itself, as in the physical Test 6 input.
    out[i, 0] = b0 - slope * tau0;
    out[i, 1] = slope;
    tau0      = tau[i];
  }
  return out;
}

//! The cosines of the upward streams of nquad streams
Vector upward_streams(const Index nquad) {
  Vector mu(nquad / 2), weights(nquad / 2);
  Legendre::PositiveDoubleGaussLegendre(mu, weights);
  return mu;
}

//! The beam cosine a surface reflects: mu0 with a beam, NaN without
Numeric beam_cosine(const Numeric beam, const Numeric mu0) {
  return beam == 0.0 ? std::numeric_limits<Numeric>::quiet_NaN() : mu0;
}

disort::main_data make_disort(const Index                NQuad,
                              AscendingGrid              tau,
                              Vector                     omega,
                              Matrix                     legendre,
                              const Numeric              mu0,
                              const Numeric              I0,
                              Matrix                     boundary_up   = {},
                              Matrix                     boundary_down = {},
                              Matrix                     source        = {},
                              disort::surface_reflection surface       = {},
                              Vector                     f             = {}) {
  const Index nfourier = NQuad;
  const Index nlayers  = static_cast<Index>(tau.size());
  if (boundary_up.size() == 0) {
    boundary_up.resize(nfourier, NQuad / 2);
    boundary_up = 0;
  }
  if (boundary_down.size() == 0) {
    boundary_down.resize(nfourier, NQuad / 2);
    boundary_down = 0;
  }
  if (source.size() == 0) source.resize(tau.size(), 0);
  if (f.empty()) {
    f.resize(nlayers);
    f = 0.0;
  }

  return disort::main_data(NQuad,
                           NQuad,
                           nfourier,
                           std::move(tau),
                           std::move(omega),
                           std::move(legendre),
                           std::move(boundary_up),
                           std::move(boundary_down),
                           std::move(f),
                           std::move(source),
                           std::move(surface),
                           mu0,
                           I0,
                           0);
}

void run_case(const std::string_view   name,
              const disort::main_data& dis,
              const Vector&            taus,
              const Vector&            phis = Vector{0.0}) {
  disort::u_data    u;
  disort::flux_data flux;
  for (const Numeric tau : taus) {
    for (const Numeric phi : phis) {
      dis.u(u, tau, phi);
      for (const Numeric value : u.intensities) {
        if (not std::isfinite(value)) arts::user_error("{} produced a non-finite intensity", name);
      }
    }

    const auto values = dis.flux(flux, tau);
    if (not(std::isfinite(values.up) and std::isfinite(values.down_diffuse) and std::isfinite(values.down_direct) and
            std::isfinite(values.dfdt)))
      arts::user_error("{} produced a non-finite flux", name);
  }
  std::cout << std::format("{} completed\n", name);
}

void disort_test06() {
  using namespace disort_test::reference;
  constexpr Index nquad = 16;

  for (Index case_index = 0; case_index < static_cast<Index>(problem_6.size()); ++case_index) {
    const auto&         test      = problem_6[case_index];
    const auto&         reference = problem_6_flux[case_index];
    const Numeric       depth     = test.optical_depth == 0.0 ? 1e-12 : test.optical_depth;
    const AscendingGrid tau{depth};
    Matrix              up(nquad, nquad / 2, 0.0), down(nquad, nquad / 2, 0.0);

    disort::brdf::RawFunction raw;
    if (test.surface == surface_type::lambertian)
      raw = [albedo = test.lambertian_albedo](Numeric, Numeric, Numeric) { return albedo * Constant::inv_pi; };
    else if (test.surface == surface_type::hapke)
      raw = disort::brdf::Hapke{.opposition_amplitude     = test.hapke_parameters[0],
                                .opposition_width         = test.hapke_parameters[1],
                                .single_scattering_albedo = test.hapke_parameters[2]};

    const Numeric top_planck =
        band_blackbody_radiance(test.top_temperature, problem_6_wavenumber_low, problem_6_wavenumber_high);
    down[0] = test.top_isotropic + test.top_emissivity * top_planck;

    // The surface's Planck radiance; DISORT applies its emissivity, as DISORT's BEMST
    up[0] = band_blackbody_radiance(test.bottom_temperature, problem_6_wavenumber_low, problem_6_wavenumber_high);

    Matrix source;
    if (test.interface_temperature.size() == 2) {
      source = linear_polynomial(
          tau,
          Vector{band_blackbody_radiance(
                     test.interface_temperature[0], problem_6_wavenumber_low, problem_6_wavenumber_high),
                 band_blackbody_radiance(
                     test.interface_temperature[1], problem_6_wavenumber_low, problem_6_wavenumber_high)});
    }

    auto dis = make_disort(
        nquad,
        tau,
        Vector{std::max(test.single_scattering_albedo, 1e-8)},
        henyey_greenstein(Vector{0.0}, nquad),
        test.beam_mu,
        test.beam,
        std::move(up),
        std::move(down),
        std::move(source),
        raw ? disort::brdf::fourier_modes(
                  raw, nquad, upward_streams(nquad), upward_streams(nquad), beam_cosine(test.beam, test.beam_mu))
            : disort::surface_reflection{});

    disort::flux_data flux_data;
    for (Index level = 0; level < static_cast<Index>(test.output_tau.size()); ++level) {
      const Numeric output_tau = test.optical_depth == 0.0 ? 0.0 : test.output_tau[level];
      const auto    flux       = dis.flux(flux_data, output_tau);
      expect_reference(
          std::format("{} direct flux [{}]", test.name, level), flux.down_direct, reference.direct[level], 7e-5);
      expect_reference(std::format("{} diffuse-down flux [{}]", test.name, level),
                       flux.down_diffuse,
                       reference.diffuse_down[level],
                       1e-3);
      // The solver's emissivity comes from the streams, cDISORT's from a fine quadrature; with the Hapke
      // hotspot, 6h's DFDT at the surface (a small difference of large terms) differs by 0.4%
      const Numeric thermal_brdf_tolerance = case_index >= 5 ? 1e-2 : 1e-3;
      expect_reference(
          std::format("{} up flux [{}]", test.name, level), flux.up, reference.up[level], thermal_brdf_tolerance);
      expect_reference(
          std::format("{} DFDT [{}]", test.name, level), flux.dfdt, reference.dfdt[level], thermal_brdf_tolerance);
    }
    std::cout << std::format("test_{} completed\n", test.name);
  }
}

void disort_test07() {
  for (Index case_index = 0; case_index < static_cast<Index>(disort_test::reference::problem_7.size()); ++case_index) {
    const auto&         c         = disort_test::reference::problem_7[case_index];
    const auto&         reference = disort_test::reference::problem_7_flux[case_index];
    const AscendingGrid tau{c.optical_depth};
    const Vector        omega{c.single_scattering_albedo};
    Matrix              up(c.streams, c.streams / 2, 0), down(c.streams, c.streams / 2, 0);
    const Numeric       bottom_planck =
        band_blackbody_radiance(c.bottom_boundary_temperature, c.wavenumber_low, c.wavenumber_high);
    down[0] =
        c.top_isotropic + band_blackbody_radiance(c.top_boundary_temperature, c.wavenumber_low, c.wavenumber_high);

    // The surface's Planck radiance; DISORT applies its emissivity, as DISORT's BEMST
    disort::brdf::RawFunction raw;
    up[0] = bottom_planck;
    if (c.surface == disort_test::reference::surface_type::lambertian) {
      raw = [albedo = c.lambertian_albedo](Numeric, Numeric, Numeric) { return albedo * Constant::inv_pi; };
    } else if (c.surface == disort_test::reference::surface_type::hapke) {
      raw = disort::brdf::Hapke{};
    }

    const Vector streams = upward_streams(c.streams);
    auto         surface = raw ? disort::brdf::fourier_modes(raw, c.streams, streams, streams, beam_cosine(c.beam, 0.5))
                               : disort::surface_reflection{};

    auto dis = make_disort(
        c.streams,
        tau,
        omega,
        henyey_greenstein(Vector{c.asymmetry_parameter}, c.streams),
        0.5,
        c.beam,
        up,
        down,
        linear_polynomial(
            tau,
            Vector{band_blackbody_radiance(c.atmosphere_top_temperature, c.wavenumber_low, c.wavenumber_high),
                   band_blackbody_radiance(c.atmosphere_bottom_temperature, c.wavenumber_low, c.wavenumber_high)}),
        std::move(surface));
    const Vector      output_tau = case_index < 2 ? Vector{0.0, c.optical_depth} : Vector{0.0, 0.5, 1.0};
    disort::flux_data flux_data;
    for (Index level = 0; level < static_cast<Index>(output_tau.size()); ++level) {
      const auto flux = dis.flux(flux_data, output_tau[level]);
      if (case_index == 1) {
        expect_small_reference(
            std::format("{} direct flux [{}]", c.name, level), flux.down_direct, reference.direct[level]);
        expect_small_reference(
            std::format("{} diffuse-down flux [{}]", c.name, level), flux.down_diffuse, reference.diffuse_down[level]);
        expect_small_reference(std::format("{} up flux [{}]", c.name, level), flux.up, reference.up[level]);
        expect_small_reference(std::format("{} DFDT [{}]", c.name, level), flux.dfdt, reference.dfdt[level]);
      } else {
        constexpr Numeric tolerance = 2e-3;
        expect_reference(
            std::format("{} direct flux [{}]", c.name, level), flux.down_direct, reference.direct[level], tolerance);
        expect_reference(std::format("{} diffuse-down flux [{}]", c.name, level),
                         flux.down_diffuse,
                         reference.diffuse_down[level],
                         tolerance);
        expect_reference(std::format("{} up flux [{}]", c.name, level), flux.up, reference.up[level], tolerance);
        expect_reference(std::format("{} DFDT [{}]", c.name, level), flux.dfdt, reference.dfdt[level], tolerance);
      }
    }
    std::cout << std::format("test_{} completed\n", c.name);
  }
}

disort::main_data problem10(const Index nquad) {
  const auto&         test = disort_test::reference::problem_9[2];
  const AscendingGrid tau{Vector(test.cumulative_tau)};
  Matrix              up(nquad, nquad / 2, 0), down(nquad, nquad / 2, 0);
  up[0] = band_blackbody_radiance(test.bottom_temperature, test.wavenumber_low, test.wavenumber_high);
  down[0] =
      test.top_isotropic + band_blackbody_radiance(test.top_temperature, test.wavenumber_low, test.wavenumber_high);
  Vector source_radiance(test.interface_temperature.size());
  for (Size i = 0; i < test.interface_temperature.size(); ++i)
    source_radiance[i] =
        band_blackbody_radiance(test.interface_temperature[i], test.wavenumber_low, test.wavenumber_high);
  const Numeric mu0     = nquad == 2 ? 0.500001 : 0.5;
  const Vector  streams = upward_streams(nquad);
  return make_disort(nquad,
                     tau,
                     test.single_scattering_albedo,
                     henyey_greenstein(disort_test::reference::problem_9c_asymmetry, nquad),
                     mu0,
                     Constant::pi,
                     up,
                     down,
                     linear_polynomial(tau, source_radiance),
                     disort::brdf::lambertian_fourier_modes(test.surface_albedo, 1, streams, streams, mu0));
}

void disort_test09c() {
  // Problem 10 is the four-stream form of the original fully general 9c.
  // Retain 9c itself with eight streams and its five requested optical depths.
  const auto dis = problem10(8);
  run_case(
      "test_9c", dis, Vector{0.0, 1.05, 2.1, 6.0, 21.0}, Vector{Constant::pi / 3, 2 * Constant::pi / 3, Constant::pi});
}

void disort_test10() {
  const auto   dis = problem10(4);
  const Vector taus{0.0, 2.1, 21.0};
  const Vector phis{Constant::pi / 3, 2 * Constant::pi / 3};

  // The C++ solver always returns the quadrature streams.  Exercise both its
  // scalar and bulk evaluation paths, the counterparts of USRANG true/false.
  const Tensor3 pointwise = compute_u(dis, taus, phis, false);
  Tensor3       bulk(taus.size(), phis.size(), dis.mu().size());
  dis.ungridded_u(bulk, AscendingGrid{taus}, phis);
  Tensor3 bulk_in_pointwise_order(phis.size(), taus.size(), dis.mu().size());
  for (Size iphi = 0; iphi < phis.size(); ++iphi) {
    for (Size itau = 0; itau < taus.size(); ++itau) { bulk_in_pointwise_order[iphi, itau] = bulk[itau, iphi]; }
  }
  if (not is_good(pointwise, bulk_in_pointwise_order)) arts::user_error("Test 10 scalar and bulk angle paths differ");
  run_case("test_10a (scalar angle path)", dis, taus, phis);
  std::cout << "test_10b (quadrature/bulk angle path) completed\n";
}

void disort_test12() {
  constexpr Index nquad = disort_test::reference::problem_12_streams;
  const Vector    g1{disort_test::reference::problem_12_asymmetry};
  auto            one_layer = make_disort(nquad,
                                          AscendingGrid{disort_test::reference::problem_12_output_tau.back()},
                                          Vector{disort_test::reference::problem_12_omega},
                                          henyey_greenstein(g1, nquad),
                                          disort_test::reference::problem_12_beam_mu,
                                          disort_test::reference::problem_12_beam);
  auto            split = make_disort(nquad,
                                      AscendingGrid{Vector(disort_test::reference::problem_12_subdivided_tau)},
                                      Vector(3, disort_test::reference::problem_12_omega),
                                      henyey_greenstein(Vector(3, disort_test::reference::problem_12_asymmetry), nquad),
                                      disort_test::reference::problem_12_beam_mu,
                                      disort_test::reference::problem_12_beam);
  run_case("test_12a", one_layer, disort_test::reference::problem_12_output_tau);
  run_case("test_12b", split, disort_test::reference::problem_12_output_tau);
}

void disort_test13() {
  constexpr Index nquad = 16;
  for (const bool multilayer : {false, true}) {
    const AscendingGrid tau   = multilayer ? AscendingGrid{0.5, 1.0} : AscendingGrid{1.0};
    const Vector        omega = multilayer ? Vector{0.99, 0.50} : Vector{0.99};
    const Vector        g(tau.size(), 0.8);
    const Vector        streams = upward_streams(nquad);
    auto dis = make_disort(nquad,
                           tau,
                           omega,
                           henyey_greenstein(g, nquad),
                           0.5,
                           2.0,
                           {},
                           {},
                           {},
                           disort::brdf::lambertian_fourier_modes(
                               disort_test::reference::problem_13_surface_albedo, 1, streams, streams, 0.5),
                           Vector(tau.size(), std::pow(0.8, nquad)));
    if (multilayer) {
      run_case("test_13c (special-boundary equivalent)", dis, Vector{0.0, 1.0});
      run_case("test_13d (regular method)", dis, Vector{0.0, 1.0});
    } else {
      run_case("test_13a (special-boundary equivalent)", dis, Vector{0.0, 1.0});
      run_case("test_13b (regular method)", dis, Vector{0.0, 1.0});
    }
  }
}

void disort_test14() {
  // twostr() has no separate C++ implementation.  Keep both original sides
  // executable as the two-stream (NQuad=2) and four-stream configurations.
  const Vector taus{0.0, 2.1, 21.0};
  const auto   four_stream = problem10(4);
  const auto   two_stream  = problem10(2);
  run_case("test_14a (four stream)", four_stream, taus);
  run_case("test_14b (two stream)", two_stream, taus);
}
}  // namespace
}  // namespace legacy_disotest

int main() try {
  std::cout << std::setprecision(16);

  legacy_disotest::disort_test06();
  legacy_disotest::disort_test07();
  legacy_disotest::disort_test09c();

  legacy_disotest::disort_test10();

  legacy_disotest::disort_test12();
  legacy_disotest::disort_test13();
  legacy_disotest::disort_test14();

  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::cerr << "Error in disotest:\n" << e.what() << '\n';
  return EXIT_FAILURE;
}
