#include <legendre.h>

#include "../reference-data.h"
#include "../test-helpers.h"

namespace legacy_disotest {
namespace {
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

Numeric hapke(const Numeric mu,
              const Numeric mup,
              const Numeric dphi,
              const Numeric b0 = 1.0,
              const Numeric hh = 0.06,
              const Numeric w  = 0.6) {
  const Numeric ctheta =
      std::clamp(mu * mup + std::sqrt((1.0 - mu * mu) * (1.0 - mup * mup)) * std::cos(dphi), -1.0, 1.0);
  const Numeric theta      = std::acos(ctheta);
  const Numeric opposition = b0 * hh / (hh + std::tan(0.5 * theta));
  const Numeric gamma      = std::sqrt(1.0 - w);
  const Numeric h0         = (1.0 + 2.0 * mup) / (1.0 + 2.0 * gamma * mup);
  const Numeric h          = (1.0 + 2.0 * mu) / (1.0 + 2.0 * gamma * mu);
  return 0.25 * w * ((1.0 + opposition) * (1.0 + 0.5 * ctheta) + h0 * h - 1.0) / (mu + mup);
}

Matrix henyey_greenstein(const Vector& g, const Index nmom) {
  Matrix out(g.size(), nmom + 1, 1.0);
  for (Index l = 1; l <= nmom; ++l) {
    for (Size i = 0; i < g.size(); ++i) out[i, l] = std::pow(g[i], l);
  }
  return out;
}

Matrix linear_polynomial(const AscendingGrid& tau, const Vector& values) {
  ARTS_USER_ERROR_IF(values.size() != tau.size() + 1, "A value is required at every layer interface");

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

Matrix linear_source(const AscendingGrid& tau, const Vector& temperature) {
  Vector radiance(temperature.size());
  std::ranges::transform(temperature, radiance.begin(), blackbody_radiance);
  return linear_polynomial(tau, radiance);
}

disort::main_data make_disort(const Index               NQuad,
                              AscendingGrid             tau,
                              Vector                    omega,
                              Matrix                    legendre,
                              const Numeric             mu0,
                              const Numeric             I0,
                              Matrix                    boundary_up   = {},
                              Matrix                    boundary_down = {},
                              Matrix                    source        = {},
                              std::vector<disort::BDRF> brdf          = {},
                              Vector                    f             = {}) {
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
                           std::move(brdf),
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
        ARTS_USER_ERROR_IF(not std::isfinite(value), "{} produced a non-finite intensity", name);
      }
    }

    const auto [down_diffuse, down_direct] = dis.flux_down(flux, tau);
    const Numeric up                       = dis.flux_up(flux, tau);
    ARTS_USER_ERROR_IF(not(std::isfinite(up) and std::isfinite(down_diffuse) and std::isfinite(down_direct)),
                       "{} produced a non-finite flux",
                       name);
  }
  std::cout << std::format("{} completed\n", name);
}

void disort_test06() {
  constexpr Index NQuad = 16;
  const Vector    phis{Constant::pi / 2};

  // 6a: transparent medium and beam source.  A strictly zero optical
  // thickness is not a valid main_data grid, so use its floating-point limit.
  {
    const AscendingGrid tau{1e-12};
    auto                dis = make_disort(NQuad, tau, Vector{1e-8}, henyey_greenstein(Vector{0.0}, NQuad), 0.5, 200.0);
    run_case("test_6a", dis, Vector{0.0, 1e-12}, phis);
  }

  // 6b: add optical depth.
  {
    const AscendingGrid tau{1.0};
    auto                dis = make_disort(NQuad, tau, Vector{1e-8}, henyey_greenstein(Vector{0.0}, NQuad), 0.5, 200.0);
    run_case("test_6b", dis, Vector{0.0, 0.5, 1.0}, phis);
  }

  // 6c and 6d: Lambertian and non-Lambertian reflection.  The C++ API
  // represents both through Fourier BRDF callbacks.
  for (const auto& [name, albedo] : {std::pair{"test_6c", 0.5}, std::pair{"test_6d", 0.25}}) {
    const AscendingGrid       tau{1.0};
    std::vector<disort::BDRF> brdf{disort::BDRF{[=](auto c, auto&, auto&) { c = albedo; }}};
    auto                      dis = make_disort(
        NQuad, tau, Vector{1e-8}, henyey_greenstein(Vector{0.0}, NQuad), 0.5, 200.0, {}, {}, {}, std::move(brdf));
    run_case(name, dis, Vector{0.0, 0.5, 1.0}, phis);
  }

  // 6e: bottom emission.
  {
    const AscendingGrid tau{1.0};
    Matrix              up(NQuad, NQuad / 2, 0);
    up[0]    = blackbody_radiance(300.0);
    auto dis = make_disort(NQuad, tau, Vector{1e-8}, henyey_greenstein(Vector{0.0}, NQuad), 0.5, 200.0, up);
    run_case("test_6e", dis, Vector{0.0, 0.5, 1.0}, phis);
  }

  // 6f: prescribed plus thermally emitted top incidence and bottom emission.
  {
    const AscendingGrid tau{1.0};
    Matrix              up(NQuad, NQuad / 2, 0), down(NQuad, NQuad / 2, 0);
    up[0]    = blackbody_radiance(300.0);
    down[0]  = 100.0 * Constant::inv_pi + blackbody_radiance(250.0);
    auto dis = make_disort(NQuad, tau, Vector{1e-8}, henyey_greenstein(Vector{0.0}, NQuad), 0.5, 200.0, up, down);
    run_case("test_6f", dis, Vector{0.0, 0.5, 1.0}, phis);
  }

  // 6g and 6h: add a linear internal thermal source, then increase depth.
  for (const auto& [name, depth] : {std::pair{"test_6g", 1.0}, std::pair{"test_6h", 10.0}}) {
    const AscendingGrid tau{depth};
    const Vector        omega{1e-8};
    Matrix              up(NQuad, NQuad / 2, 0), down(NQuad, NQuad / 2, 0);
    up[0]    = blackbody_radiance(300.0);
    down[0]  = blackbody_radiance(250.0);
    auto dis = make_disort(NQuad,
                           tau,
                           omega,
                           henyey_greenstein(Vector{0.0}, NQuad),
                           0.5,
                           200.0,
                           up,
                           down,
                           linear_source(tau, Vector{250.0, 300.0}));
    run_case(name, dis, depth == 1.0 ? Vector{0.0, 0.5, 1.0} : Vector{0.0, 1.0, 10.0}, phis);
  }
}

void disort_test07() {
  for (const auto& c : disort_test::reference::problem_7) {
    const AscendingGrid tau{c.optical_depth};
    const Vector        omega{c.single_scattering_albedo};
    Matrix              up(c.streams, c.streams / 2, 0), down(c.streams, c.streams / 2, 0);
    const Numeric       bottom_planck =
        band_blackbody_radiance(c.bottom_boundary_temperature, c.wavenumber_low, c.wavenumber_high);
    down[0] =
        c.top_isotropic + band_blackbody_radiance(c.top_boundary_temperature, c.wavenumber_low, c.wavenumber_high);

    std::vector<disort::BDRF> brdf;
    if (c.surface != disort_test::reference::surface_type::hapke) up[0] = (1.0 - c.lambertian_albedo) * bottom_planck;
    if (c.surface == disort_test::reference::surface_type::lambertian) {
      brdf.push_back(disort::BDRF{[albedo = c.lambertian_albedo](auto value, auto&, auto&) { value = albedo; }});
    } else if (c.surface == disort_test::reference::surface_type::hapke) {
      constexpr Index naz = 512;
      for (Index m = 0; m < c.streams; ++m) {
        brdf.push_back(disort::BDRF{[m](auto value, const auto& outgoing, const auto& incoming) {
          const Index nout = outgoing.size();
          const Index nin  = incoming.size();
          for (Index i = 0; i < nout; ++i)
            for (Index j = 0; j < nin; ++j) {
              Numeric coefficient = 0.0;
              for (Index k = 0; k < naz; ++k) {
                const Numeric phi = Constant::two_pi * (static_cast<Numeric>(k) + 0.5) / static_cast<Numeric>(naz);
                coefficient += hapke(outgoing[i], std::abs(incoming[j]), phi) * std::cos(static_cast<Numeric>(m) * phi);
              }
              value[i, j] = coefficient / naz;
            }
        }});
      }
      // Directional Kirchhoff emissivity at the computational angles.
      Vector surface_mu(c.streams / 2), surface_weight(c.streams / 2);
      Legendre::PositiveDoubleGaussLegendre(surface_mu, surface_weight);
      for (Index i = 0; i < c.streams / 2; ++i) {
        constexpr Index nmu         = 256;
        Numeric         reflectance = 0.0;
        for (Index j = 0; j < nmu; ++j) {
          const Numeric mup             = (static_cast<Numeric>(j) + 0.5) / static_cast<Numeric>(nmu);
          Numeric       azimuth_average = 0.0;
          for (Index k = 0; k < naz; ++k)
            azimuth_average += hapke(
                surface_mu[i], mup, Constant::two_pi * (static_cast<Numeric>(k) + 0.5) / static_cast<Numeric>(naz));
          reflectance += 2.0 * mup * azimuth_average / (nmu * naz);
        }
        up[0, i] = (1.0 - reflectance) * bottom_planck;
      }
    }

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
        std::move(brdf));
    run_case(c.name,
             dis,
             c.optical_depth == 100.0 ? Vector{0.0, 100.0} : Vector{0.0, 0.5, 1.0},
             Vector{0.0, Constant::pi / 2});
  }
}

disort::main_data problem10(const Index nquad) {
  const auto&         test = disort_test::reference::problem_9[2];
  const AscendingGrid tau{Vector(test.cumulative_tau)};
  Matrix              up(nquad, nquad / 2, 0), down(nquad, nquad / 2, 0);
  up[0] = (1.0 - test.surface_albedo) *
          band_blackbody_radiance(test.bottom_temperature, test.wavenumber_low, test.wavenumber_high);
  down[0] =
      test.top_isotropic + band_blackbody_radiance(test.top_temperature, test.wavenumber_low, test.wavenumber_high);
  Vector source_radiance(test.interface_temperature.size());
  for (Size i = 0; i < test.interface_temperature.size(); ++i)
    source_radiance[i] =
        band_blackbody_radiance(test.interface_temperature[i], test.wavenumber_low, test.wavenumber_high);
  std::vector<disort::BDRF> brdf{
      disort::BDRF{[albedo = test.surface_albedo](auto value, auto&, auto&) { value = albedo; }}};
  return make_disort(nquad,
                     tau,
                     test.single_scattering_albedo,
                     henyey_greenstein(disort_test::reference::problem_9c_asymmetry, nquad),
                     nquad == 2 ? 0.500001 : 0.5,
                     Constant::pi,
                     up,
                     down,
                     linear_polynomial(tau, source_radiance),
                     std::move(brdf));
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
  ARTS_USER_ERROR_IF(not is_good(pointwise, bulk_in_pointwise_order), "Test 10 scalar and bulk angle paths differ");
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
    const AscendingGrid       tau   = multilayer ? AscendingGrid{0.5, 1.0} : AscendingGrid{1.0};
    const Vector              omega = multilayer ? Vector{0.99, 0.50} : Vector{0.99};
    const Vector              g(tau.size(), 0.8);
    std::vector<disort::BDRF> brdf{
        disort::BDRF{[](auto value, auto&, auto&) { value = disort_test::reference::problem_13_surface_albedo; }}};
    auto dis = make_disort(nquad,
                           tau,
                           omega,
                           henyey_greenstein(g, nquad),
                           0.5,
                           2.0,
                           {},
                           {},
                           {},
                           std::move(brdf),
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
