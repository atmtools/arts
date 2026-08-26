#define DISORT_TEST_NO_MAIN

// Keep the scalar DISORT reference suite in the same problem order as
// 3rdparty/cdisort/disotest.cc.  The case bodies remain separately includable
// so this translation unit can later be mirrored by the VDISORT test driver.
#include "disort-test-1.cpp"
#include "disort-test-2.cpp"
#include "disort-test-3.cpp"
#include "disort-test-4.cpp"
#define Leg_coeffs_ALL Leg_coeffs_ALL_5
#include "disort-test-5.cpp"
#undef Leg_coeffs_ALL
#include "disort-test-8.cpp"
#include "disort-test-9.cpp"
#include "disort-test-11.cpp"

namespace legacy_disotest {

Numeric blackbody_radiance(const Numeric temperature) {
  return Constant::sigma * std::pow(temperature, 4) * Constant::inv_pi;
}

Matrix henyey_greenstein(const Vector& g, const Index nmom) {
  Matrix out(g.size(), nmom + 1, 1.0);
  for (Index l = 1; l <= nmom; ++l) {
    for (Size i = 0; i < g.size(); ++i) out[i, l] = std::pow(g[i], l);
  }
  return out;
}

Matrix linear_source(const AscendingGrid& tau,
                     const Vector& temperature,
                     const Vector& omega) {
  ARTS_USER_ERROR_IF(temperature.size() != tau.size() + 1,
                     "A temperature is required at every layer interface");

  Matrix out(tau.size(), 2);
  Numeric tau0 = 0;
  for (Size i = 0; i < tau.size(); ++i) {
    const Numeric b0 = blackbody_radiance(temperature[i]);
    const Numeric b1 = blackbody_radiance(temperature[i + 1]);
    const Numeric slope = (b1 - b0) / (tau[i] - tau0);
    out[i, 0] = (b0 - slope * tau0) / (1 - omega[i]);
    out[i, 1] = slope / (1 - omega[i]);
    tau0 = tau[i];
  }
  return out;
}

disort::main_data make_disort(const Index NQuad,
                              AscendingGrid tau,
                              Vector omega,
                              Matrix legendre,
                              const Numeric mu0,
                              const Numeric I0,
                              Matrix boundary_up = {},
                              Matrix boundary_down = {},
                              Matrix source = {},
                              std::vector<disort::BDRF> brdf = {}) {
  const Index nfourier = NQuad;
  const Index nlayers = static_cast<Index>(tau.size());
  if (boundary_up.size() == 0) {
    boundary_up.resize(nfourier, NQuad / 2);
    boundary_up = 0;
  }
  if (boundary_down.size() == 0) {
    boundary_down.resize(nfourier, NQuad / 2);
    boundary_down = 0;
  }
  if (source.size() == 0) source.resize(tau.size(), 0);

  return disort::main_data(NQuad,
                           NQuad,
                           nfourier,
                           std::move(tau),
                           std::move(omega),
                           std::move(legendre),
                           std::move(boundary_up),
                           std::move(boundary_down),
                           Vector(nlayers, 0),
                           std::move(source),
                           std::move(brdf),
                           mu0,
                           I0,
                           0);
}

void run_case(const std::string_view name,
              const disort::main_data& dis,
              const Vector& taus,
              const Vector& phis = Vector{0.0}) {
  disort::u_data u;
  disort::flux_data flux;
  for (const Numeric tau : taus) {
    for (const Numeric phi : phis) {
      dis.u(u, tau, phi);
      for (const Numeric value : u.intensities) {
        ARTS_USER_ERROR_IF(not std::isfinite(value), "{} produced a non-finite intensity", name);
      }
    }

    const auto [down_diffuse, down_direct] = dis.flux_down(flux, tau);
    const Numeric up = dis.flux_up(flux, tau);
    ARTS_USER_ERROR_IF(not(std::isfinite(up) and std::isfinite(down_diffuse) and std::isfinite(down_direct)),
                       "{} produced a non-finite flux",
                       name);
  }
  std::cout << std::format("{} completed\n", name);
}

void disort_test06() {
  constexpr Index NQuad = 16;
  const Vector phis{Constant::pi / 2};

  // 6a: transparent medium and beam source.  A strictly zero optical
  // thickness is not a valid main_data grid, so use its floating-point limit.
  {
    const AscendingGrid tau{1e-12};
    auto dis = make_disort(NQuad, tau, Vector{1e-8}, henyey_greenstein(Vector{0.0}, NQuad), 0.5, 200.0);
    run_case("test_6a", dis, Vector{0.0, 1e-12}, phis);
  }

  // 6b: add optical depth.
  {
    const AscendingGrid tau{1.0};
    auto dis = make_disort(NQuad, tau, Vector{1e-8}, henyey_greenstein(Vector{0.0}, NQuad), 0.5, 200.0);
    run_case("test_6b", dis, Vector{0.0, 0.5, 1.0}, phis);
  }

  // 6c and 6d: Lambertian and non-Lambertian reflection.  The C++ API
  // represents both through Fourier BRDF callbacks.
  for (const auto& [name, albedo] : {std::pair{"test_6c", 0.5}, std::pair{"test_6d", 0.25}}) {
    const AscendingGrid tau{1.0};
    std::vector<disort::BDRF> brdf{disort::BDRF{[=](auto c, auto&, auto&) { c = albedo; }}};
    auto dis = make_disort(NQuad,
                           tau,
                           Vector{1e-8},
                           henyey_greenstein(Vector{0.0}, NQuad),
                           0.5,
                           200.0,
                           {},
                           {},
                           {},
                           std::move(brdf));
    run_case(name, dis, Vector{0.0, 0.5, 1.0}, phis);
  }

  // 6e: bottom emission.
  {
    const AscendingGrid tau{1.0};
    Matrix up(NQuad, NQuad / 2, 0);
    up[0] = blackbody_radiance(300.0);
    auto dis = make_disort(NQuad, tau, Vector{1e-8}, henyey_greenstein(Vector{0.0}, NQuad), 0.5, 200.0, up);
    run_case("test_6e", dis, Vector{0.0, 0.5, 1.0}, phis);
  }

  // 6f: prescribed plus thermally emitted top incidence and bottom emission.
  {
    const AscendingGrid tau{1.0};
    Matrix up(NQuad, NQuad / 2, 0), down(NQuad, NQuad / 2, 0);
    up[0] = blackbody_radiance(300.0);
    down[0] = 100.0 * Constant::inv_pi + blackbody_radiance(250.0);
    auto dis = make_disort(
        NQuad, tau, Vector{1e-8}, henyey_greenstein(Vector{0.0}, NQuad), 0.5, 200.0, up, down);
    run_case("test_6f", dis, Vector{0.0, 0.5, 1.0}, phis);
  }

  // 6g and 6h: add a linear internal thermal source, then increase depth.
  for (const auto& [name, depth] : {std::pair{"test_6g", 1.0}, std::pair{"test_6h", 10.0}}) {
    const AscendingGrid tau{depth};
    const Vector omega{1e-8};
    Matrix up(NQuad, NQuad / 2, 0), down(NQuad, NQuad / 2, 0);
    up[0] = blackbody_radiance(300.0);
    down[0] = blackbody_radiance(250.0);
    auto dis = make_disort(NQuad,
                           tau,
                           omega,
                           henyey_greenstein(Vector{0.0}, NQuad),
                           0.5,
                           200.0,
                           up,
                           down,
                           linear_source(tau, Vector{250.0, 300.0}, omega));
    run_case(name, dis, depth == 1.0 ? Vector{0.0, 0.5, 1.0} : Vector{0.0, 1.0, 10.0}, phis);
  }
}

void disort_test07() {
  struct Case {
    std::string_view name;
    Index nquad;
    Numeric depth;
    Numeric omega;
    Numeric g;
    Numeric top_temperature;
    Numeric bottom_temperature;
    Numeric beam;
    Numeric top_isotropic;
    Numeric bottom_albedo;
  };

  const std::array cases{
      Case{"test_7a", 16, 1.0, 0.10, 0.05, 200.0, 300.0, 0.0, 0.0, 0.0},
      Case{"test_7b", 16, 100.0, 0.95, 0.75, 200.0, 300.0, 0.0, 0.0, 0.0},
      Case{"test_7c", 12, 1.0, 0.50, 0.80, 300.0, 200.0, 200.0, 100.0, 0.0},
      Case{"test_7d", 12, 1.0, 0.50, 0.80, 300.0, 200.0, 200.0, 100.0, 1.0},
      Case{"test_7e", 12, 1.0, 0.50, 0.80, 300.0, 200.0, 200.0, 100.0, 0.5},
  };

  for (const auto& c : cases) {
    const AscendingGrid tau{c.depth};
    const Vector omega{c.omega};
    Matrix up(c.nquad, c.nquad / 2, 0), down(c.nquad, c.nquad / 2, 0);
    up[0] = c.bottom_temperature == 0 ? 0 : blackbody_radiance(c.bottom_temperature);
    down[0] = c.top_isotropic + blackbody_radiance(c.top_temperature);

    std::vector<disort::BDRF> brdf;
    if (c.bottom_albedo != 0) {
      brdf.push_back(
          disort::BDRF{[albedo = c.bottom_albedo](auto value, auto&, auto&) { value = albedo; }});
    }

    auto dis = make_disort(c.nquad,
                           tau,
                           omega,
                           henyey_greenstein(Vector{c.g}, c.nquad),
                           0.5,
                           c.beam,
                           up,
                           down,
                           linear_source(tau, Vector{c.top_temperature, c.bottom_temperature}, omega),
                           std::move(brdf));
    run_case(c.name,
             dis,
             c.depth == 100.0 ? Vector{0.0, 100.0} : Vector{0.0, 0.5, 1.0},
             Vector{0.0, Constant::pi / 2});
  }
}

disort::main_data problem10(const Index nquad) {
  const AscendingGrid tau{1.0, 3.0, 6.0, 10.0, 15.0, 21.0};
  const Vector omega{0.65, 0.70, 0.75, 0.80, 0.85, 0.90};
  const Vector g{1.0 / 7, 2.0 / 7, 3.0 / 7, 4.0 / 7, 5.0 / 7, 6.0 / 7};
  const Vector temperature{600.0, 610.0, 620.0, 630.0, 640.0, 650.0, 660.0};
  Matrix up(nquad, nquad / 2, 0), down(nquad, nquad / 2, 0);
  up[0] = blackbody_radiance(700.0);
  down[0] = 1.0 + blackbody_radiance(550.0);
  std::vector<disort::BDRF> brdf{disort::BDRF{[](auto value, auto&, auto&) { value = 0.5; }}};
  return make_disort(nquad,
                     tau,
                     omega,
                     henyey_greenstein(g, nquad),
                     nquad == 2 ? 0.500001 : 0.5,
                     Constant::pi,
                     up,
                     down,
                     linear_source(tau, temperature, omega),
                     std::move(brdf));
}

void disort_test09c() {
  // Problem 10 is the four-stream form of the original fully general 9c.
  // Retain 9c itself with eight streams and its five requested optical depths.
  const auto dis = problem10(8);
  run_case("test_9c",
           dis,
           Vector{0.0, 1.05, 2.1, 6.0, 21.0},
           Vector{Constant::pi / 3, 2 * Constant::pi / 3, Constant::pi});
}

void disort_test10() {
  const auto dis = problem10(4);
  const Vector taus{0.0, 2.1, 21.0};
  const Vector phis{Constant::pi / 3, 2 * Constant::pi / 3};

  // The C++ solver always returns the quadrature streams.  Exercise both its
  // scalar and bulk evaluation paths, the counterparts of USRANG true/false.
  const Tensor3 pointwise = compute_u(dis, taus, phis, false);
  Tensor3 bulk(taus.size(), phis.size(), dis.mu().size());
  dis.ungridded_u(bulk, AscendingGrid{taus}, phis);
  Tensor3 bulk_in_pointwise_order(phis.size(), taus.size(), dis.mu().size());
  for (Size iphi = 0; iphi < phis.size(); ++iphi) {
    for (Size itau = 0; itau < taus.size(); ++itau) {
      bulk_in_pointwise_order[iphi, itau] = bulk[itau, iphi];
    }
  }
  ARTS_USER_ERROR_IF(not is_good(pointwise, bulk_in_pointwise_order),
                     "Test 10 scalar and bulk angle paths differ");
  run_case("test_10a (scalar angle path)", dis, taus, phis);
  std::cout << "test_10b (quadrature/bulk angle path) completed\n";
}

void disort_test12() {
  const Index nquad = 20;
  const Vector g1{0.9};
  auto one_layer = make_disort(nquad,
                               AscendingGrid{20.1},
                               Vector{0.5},
                               henyey_greenstein(g1, nquad),
                               1.0,
                               1.0);
  auto split = make_disort(nquad,
                           AscendingGrid{10.0, 19.9, 20.1},
                           Vector(3, 0.5),
                           henyey_greenstein(Vector(3, 0.9), nquad),
                           1.0,
                           1.0);
  const Vector taus{0.0, 10.0, 19.9, 20.1};
  run_case("test_12a", one_layer, taus);
  run_case("test_12b", split, taus);
}

void disort_test13() {
  constexpr Index nquad = 16;
  for (const bool multilayer : {false, true}) {
    const AscendingGrid tau = multilayer ? AscendingGrid{0.5, 1.0} : AscendingGrid{1.0};
    const Vector omega = multilayer ? Vector{0.99, 0.50} : Vector{0.99};
    const Vector g(tau.size(), 0.8);
    auto dis = make_disort(nquad, tau, omega, henyey_greenstein(g, nquad), 0.5, 2.0);
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
  const auto four_stream = problem10(4);
  const auto two_stream = problem10(2);
  run_case("test_14a (four stream)", four_stream, taus);
  run_case("test_14b (two stream)", two_stream, taus);
}

}  // namespace legacy_disotest

int main() try {
  std::cout << std::setprecision(16);

  test_1a();
  test_1b();
  test_1c();
  test_1d();
  test_1e();
  test_1f();

  test_2a();
  test_2b();
  test_2c();
  test_2d();

  test_3a();
  test_3b();

  test_4a();
  test_4b();
  test_4c();

  test_5a();
  test_5b();
  test_5BDRF();

  legacy_disotest::disort_test06();
  legacy_disotest::disort_test07();

  test_8a();
  test_8b();
  test_8c();

  test_9a();
  test_9b();
  legacy_disotest::disort_test09c();

  legacy_disotest::disort_test10();

  test_11a_1layer();
  test_11a_multilayer();

  legacy_disotest::disort_test12();
  legacy_disotest::disort_test13();
  legacy_disotest::disort_test14();

  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::cerr << "Error in disotest:\n" << e.what() << '\n';
  return EXIT_FAILURE;
}
