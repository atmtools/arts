#include <arts_constants.h>
#include <disort-brdf.h>
#include <disort.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <utility>

#include "../reference-data.h"

namespace {
void expect_reference(const std::string_view name,
                      const Numeric          actual,
                      const Numeric          expected,
                      const Numeric          tolerance = 7e-5) {
  ARTS_USER_ERROR_IF(std::abs(actual - expected) > tolerance * std::max(1.0, std::abs(expected)),
                     "{}: expected {}, got {} (difference {})",
                     name,
                     expected,
                     actual,
                     actual - expected);
}

void run_isotropic_case(const disort_test::reference::single_layer_case& test) {
  constexpr Index nquad = disort_test::reference::problem_1_streams;
  Matrix          legendre(1, 17, 0.0);
  legendre[0, 0] = 1.0;
  Matrix down(nquad, nquad / 2, 0.0);
  if (not test.beam) down[0] = 1.0;

  const Numeric           mu0 = disort_test::reference::problem_1_beam_mu;
  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{test.depth},
                              Vector{test.omega},
                              legendre,
                              Matrix(nquad, nquad / 2, 0.0),
                              down,
                              Vector{0.0},
                              Matrix(1, 0),
                              {},
                              mu0,
                              test.beam ? Constant::pi / mu0 : 0.0,
                              0.0);

  const Vector        tau{0.0, test.depth};
  disort::user_u_data user;
  disort::flux_data   flux;
  for (Index level = 0; level < 2; ++level) {
    dis.u_user(user, tau[level], 0.0, disort_test::reference::problem_1_user_mu);
    for (Index angle = 0; angle < 6; ++angle)
      expect_reference(std::format("{} radiance [{}, {}]", test.name, level, angle),
                       user.intensities[angle],
                       test.radiance[level, angle]);

    const auto [diffuse_down, direct] = dis.flux_down(flux, tau[level]);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), direct, test.direct[level]);
    expect_reference(
        std::format("{} diffuse-down flux [{}]", test.name, level), diffuse_down, test.diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), dis.flux_up(flux, tau[level]), test.up[level]);
  }
}

void test_problem_1() {
  for (const auto& test : disort_test::reference::problem_1) run_isotropic_case(test);
}

void run_rayleigh_case(const disort_test::reference::single_layer_case& test) {
  constexpr Index nquad = disort_test::reference::problem_2_streams;
  Matrix          legendre(1, 17, 0.0);
  legendre[0, 0]              = 1.0;
  legendre[0, 2]              = 0.1;
  constexpr Numeric       mu0 = disort_test::reference::problem_2_beam_mu;
  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{test.depth},
                              Vector{test.omega},
                              legendre,
                              Matrix(nquad, nquad / 2, 0.0),
                              Matrix(nquad, nquad / 2, 0.0),
                              Vector{0.0},
                              Matrix(1, 0),
                              {},
                              mu0,
                              Constant::pi,
                              0.0);
  disort::user_u_data     user;
  disort::flux_data       flux;
  for (Index level = 0; level < 2; ++level) {
    const Numeric tau = level == 0 ? 0.0 : test.depth;
    dis.u_user(user, tau, 0.0, disort_test::reference::problem_2_user_mu);
    for (Index angle = 0; angle < 6; ++angle)
      expect_reference(std::format("{} radiance [{}, {}]", test.name, level, angle),
                       user.intensities[angle],
                       test.radiance[level, angle]);
    const auto [down, beam] = dis.flux_down(flux, tau);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), beam, test.direct[level]);
    expect_reference(std::format("{} diffuse-down flux [{}]", test.name, level), down, test.diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), dis.flux_up(flux, tau), test.up[level]);
  }
}

void test_problem_2() {
  for (const auto& test : disort_test::reference::problem_2) run_rayleigh_case(test);
}

void run_henyey_greenstein_case(const disort_test::reference::single_layer_case& test) {
  constexpr Index nquad = disort_test::reference::problem_3_streams;
  Matrix          legendre(1, disort_test::reference::problem_3_moments, 0.0);
  for (Index i = 0; i < disort_test::reference::problem_3_moments; ++i)
    legendre[0, i] = std::pow(disort_test::reference::problem_3_asymmetry, i);

  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{test.depth},
                              Vector{1.0},
                              legendre,
                              Matrix(nquad, nquad / 2, 0.0),
                              Matrix(nquad, nquad / 2, 0.0),
                              Vector{legendre[0, nquad]},
                              Matrix(1, 0),
                              {},
                              1.0,
                              Constant::pi,
                              0.0);
  disort::user_u_data     user;
  disort::tms_data        tms;
  disort::flux_data       flux;
  Vector                  ims;
  for (Index level = 0; level < 2; ++level) {
    const Numeric tau = level == 0 ? 0.0 : test.depth;
    dis.u_user_corr(user, ims, tms, tau, 0.0, disort_test::reference::problem_3_user_mu);
    for (Index angle = 0; angle < 6; ++angle)
      expect_reference(std::format("{} radiance [{}, {}]", test.name, level, angle),
                       user.intensities[angle],
                       test.radiance[level, angle]);
    const auto [down, beam] = dis.flux_down(flux, tau);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), beam, test.direct[level]);
    expect_reference(std::format("{} diffuse-down flux [{}]", test.name, level), down, test.diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), dis.flux_up(flux, tau), test.up[level]);
  }
}

void test_problem_3() {
  for (const auto& test : disort_test::reference::problem_3) run_henyey_greenstein_case(test);
}

void run_haze_l_case(const disort_test::reference::haze_l_case& test) {
  constexpr Index         nquad   = disort_test::reference::problem_4_streams;
  const Matrix            moments = disort_test::reference::haze_l_moments();
  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{disort_test::reference::problem_4_total_tau},
                              Vector{test.omega},
                              moments,
                              Matrix(nquad, nquad / 2, 0.0),
                              Matrix(nquad, nquad / 2, 0.0),
                              Vector{moments[0, nquad]},
                              Matrix(1, 0),
                              {},
                              test.beam_mu,
                              Constant::pi,
                              0.0);
  disort::user_u_data     user;
  disort::tms_data        tms;
  disort::flux_data       flux;
  Vector                  ims;
  const Index             nphi = test.azimuth.size();
  for (Index p = 0; p < nphi; ++p)
    for (Index level = 0; level < 3; ++level) {
      dis.u_user_corr(user,
                      ims,
                      tms,
                      disort_test::reference::problem_4_output_tau[level],
                      test.azimuth[p],
                      disort_test::reference::problem_4_user_mu);
      for (Index angle = 0; angle < 6; ++angle)
        expect_reference(std::format("{} radiance [{}, {}, {}]", test.name, p, level, angle),
                         user.intensities[angle],
                         test.radiance[p, level, angle]);
    }
  for (Index level = 0; level < 3; ++level) {
    const auto values = dis.flux(flux, disort_test::reference::problem_4_output_tau[level]);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), values.down_direct, test.direct[level]);
    expect_reference(
        std::format("{} diffuse-down flux [{}]", test.name, level), values.down_diffuse, test.diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), values.up, test.up[level]);
  }
}

void test_problem_4() {
  for (const auto& test : disort_test::reference::problem_4) run_haze_l_case(test);
}

void run_cloud_c1_case(const disort_test::reference::scalar_case& test) {
  constexpr Index         nquad   = disort_test::reference::problem_5_streams;
  const Matrix            moments = disort_test::reference::cloud_c1_moments();
  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{disort_test::reference::problem_5_total_tau},
                              Vector{test.omega},
                              moments,
                              Matrix(nquad, nquad / 2, 0.0),
                              Matrix(nquad, nquad / 2, 0.0),
                              Vector{moments[0, nquad]},
                              Matrix(1, 0),
                              {},
                              1.0,
                              Constant::pi,
                              0.0);
  disort::user_u_data     user;
  disort::tms_data        tms;
  disort::flux_data       flux;
  Vector                  ims;
  for (Index level = 0; level < 3; ++level) {
    dis.u_user_corr(user, ims, tms, test.tau[level], 0.0, disort_test::reference::problem_5_user_mu);
    for (Index angle = 0; angle < 6; ++angle)
      expect_reference(std::format("{} radiance [{}, {}]", test.name, level, angle),
                       user.intensities[angle],
                       test.radiance[level, angle],
                       2e-4);
    const auto [down, beam] = dis.flux_down(flux, test.tau[level]);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), beam, test.direct[level]);
    expect_reference(std::format("{} diffuse-down flux [{}]", test.name, level), down, test.diffuse_down[level]);
    expect_reference(
        std::format("{} up flux [{}]", test.name, level), dis.flux_up(flux, test.tau[level]), test.up[level]);
  }
}

void test_problem_5() {
  run_cloud_c1_case(disort_test::reference::problem_5a);
  run_cloud_c1_case(disort_test::reference::problem_5b);
}

void run_problem_8_case(const disort_test::reference::layered_isotropic_case& test) {
  constexpr Index nquad = disort_test::reference::problem_8_streams;
  Matrix          moments(2, nquad + 1, 0.0);
  moments[joker, 0] = 1.0;

  Matrix boundary_down(nquad, nquad / 2, 0.0);
  boundary_down[0] = Constant::inv_pi;

  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{test.cumulative_tau[0], test.cumulative_tau[1]},
                              test.single_scattering_albedo,
                              moments,
                              Matrix(nquad, nquad / 2, 0.0),
                              boundary_down,
                              Vector(2, 0.0),
                              Matrix(2, 0),
                              {},
                              0.5,
                              0.0,
                              0.0);

  disort::user_u_data user;
  disort::flux_data   flux;
  for (Index level = 0; level < 3; ++level) {
    const Numeric tau = test.output_tau[level];
    dis.u_user(user, tau, disort_test::reference::problem_8_azimuth, disort_test::reference::problem_8_user_mu);
    for (Index angle = 0; angle < 4; ++angle)
      expect_reference(std::format("{} radiance [{}, {}]", test.name, level, angle),
                       user.intensities[angle],
                       test.radiance[level, angle]);

    const auto [down, direct] = dis.flux_down(flux, tau);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), direct, test.direct[level]);
    expect_reference(std::format("{} diffuse-down flux [{}]", test.name, level), down, test.diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), dis.flux_up(flux, tau), test.up[level]);
  }
}

void test_problem_8() {
  for (const auto& test : disort_test::reference::problem_8) run_problem_8_case(test);
}

Numeric band_blackbody_radiance(const Numeric temperature,
                                const Numeric wavenumber_low,
                                const Numeric wavenumber_high) {
  if (temperature == 0.0 || wavenumber_low == wavenumber_high) return 0.0;
  constexpr Index intervals = 4096;
  const Numeric   scale     = Constant::h * Constant::c * 100.0 / (Constant::k * temperature);
  const Numeric   x0        = scale * wavenumber_low;
  const Numeric   x1        = scale * wavenumber_high;
  const Numeric   dx        = (x1 - x0) / intervals;
  const auto      integrand = [](const Numeric x) {
    if (x == 0.0 || x > 700.0) return 0.0;
    return x * x * x / std::expm1(x);
  };
  Numeric integral = integrand(x0) + integrand(x1);
  for (Index i = 1; i < intervals; ++i) integral += (i % 2 ? 4.0 : 2.0) * integrand(x0 + static_cast<Numeric>(i) * dx);
  integral *= dx / 3.0;
  return Constant::sigma * std::pow(temperature, 4) * Constant::inv_pi * integral /
         (Math::pow2(Constant::pi) * Math::pow2(Constant::pi) / 15.0);
}

Matrix problem_9_moments(const disort_test::reference::phase_type phase,
                         const Index                              nquad = disort_test::reference::problem_9_streams) {
  Matrix moments(6, nquad + 1, 0.0);
  if (phase == disort_test::reference::phase_type::isotropic) {
    moments[joker, 0] = 1.0;
  } else if (phase == disort_test::reference::phase_type::problem_9b) {
    for (Index layer = 0; layer < 6; ++layer) moments[layer] = disort_test::reference::problem_9b_moments;
  } else {
    for (Index layer = 0; layer < 6; ++layer)
      for (Index moment = 0; moment <= nquad; ++moment)
        moments[layer, moment] = std::pow(disort_test::reference::problem_9c_asymmetry[layer], moment);
  }
  return moments;
}

Matrix problem_9_source(const disort_test::reference::general_multilayer_case& test) {
  if (not test.thermal) return Matrix(6, 0);
  Vector planck(7);
  for (Index level = 0; level < 7; ++level)
    planck[level] =
        band_blackbody_radiance(test.interface_temperature[level], test.wavenumber_low, test.wavenumber_high);
  Matrix  source(6, 2);
  Numeric tau0 = 0.0;
  for (Index layer = 0; layer < 6; ++layer) {
    const Numeric tau1  = test.cumulative_tau[layer];
    const Numeric slope = (planck[layer + 1] - planck[layer]) / (tau1 - tau0);
    source[layer, 0]    = planck[layer] - slope * tau0;
    source[layer, 1]    = slope;
    tau0                = tau1;
  }
  return source;
}

void run_problem_9_case(const disort_test::reference::general_multilayer_case& test) {
  constexpr Index nquad = disort_test::reference::problem_9_streams;
  Matrix          up(nquad, nquad / 2, 0.0), down(nquad, nquad / 2, 0.0);
  down[0] = test.top_isotropic;
  std::vector<disort::BDRF> brdf;
  if (test.thermal) {
    down[0] += band_blackbody_radiance(test.top_temperature, test.wavenumber_low, test.wavenumber_high);
    up[0]    = (1.0 - test.surface_albedo) *
               band_blackbody_radiance(test.bottom_temperature, test.wavenumber_low, test.wavenumber_high);
  }
  if (test.surface_albedo != 0.0)
    brdf.push_back(disort::BDRF{[albedo = test.surface_albedo](auto value, auto&, auto&) { value = albedo; }});

  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{Vector(test.cumulative_tau)},
                              test.single_scattering_albedo,
                              problem_9_moments(test.phase),
                              up,
                              down,
                              Vector(6, 0.0),
                              problem_9_source(test),
                              std::move(brdf),
                              test.beam_mu,
                              test.beam,
                              0.0);

  disort::user_u_data user;
  disort::tms_data    tms;
  disort::flux_data   flux;
  Vector              ims;
  for (Index azimuth = 0; azimuth < static_cast<Index>(test.azimuth.size()); ++azimuth)
    for (Index level = 0; level < 5; ++level) {
      if (test.phase == disort_test::reference::phase_type::henyey_greenstein)
        dis.u_user_corr(
            user, ims, tms, test.output_tau[level], test.azimuth[azimuth], disort_test::reference::problem_9_user_mu);
      else
        dis.u_user(user, test.output_tau[level], test.azimuth[azimuth], disort_test::reference::problem_9_user_mu);
      for (Index angle = 0; angle < 4; ++angle)
        expect_reference(std::format("{} radiance [{}, {}, {}]", test.name, azimuth, level, angle),
                         user.intensities[angle],
                         test.radiance[azimuth, level, angle],
                         test.thermal ? 1e-3 : 7e-5);
    }

  for (Index level = 0; level < 5; ++level) {
    const Numeric tau              = test.output_tau[level];
    const auto [down_flux, direct] = dis.flux_down(flux, tau);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), direct, test.direct[level]);
    expect_reference(std::format("{} diffuse-down flux [{}]", test.name, level),
                     down_flux,
                     test.diffuse_down[level],
                     test.thermal ? 1e-3 : 7e-5);
    expect_reference(std::format("{} up flux [{}]", test.name, level),
                     dis.flux_up(flux, tau),
                     test.up[level],
                     test.thermal ? 1e-3 : 7e-5);
  }
}

void test_problem_9() {
  for (const auto& test : disort_test::reference::problem_9) run_problem_9_case(test);
}

void test_problem_10() {
  constexpr Index nquad = disort_test::reference::problem_10_streams;
  const auto&     test  = disort_test::reference::problem_9[2];
  Matrix          up(nquad, nquad / 2, 0.0), down(nquad, nquad / 2, 0.0);
  down[0] =
      test.top_isotropic + band_blackbody_radiance(test.top_temperature, test.wavenumber_low, test.wavenumber_high);
  up[0] = (1.0 - test.surface_albedo) *
          band_blackbody_radiance(test.bottom_temperature, test.wavenumber_low, test.wavenumber_high);
  std::vector<disort::BDRF> brdf{
      disort::BDRF{[albedo = test.surface_albedo](auto value, auto&, auto&) { value = albedo; }}};

  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{Vector(test.cumulative_tau)},
                              test.single_scattering_albedo,
                              problem_9_moments(test.phase, nquad),
                              up,
                              down,
                              Vector(6, 0.0),
                              problem_9_source(test),
                              std::move(brdf),
                              test.beam_mu,
                              test.beam,
                              0.0);

  disort::u_data      quadrature;
  disort::user_u_data formal;
  // CPP stores +mu from grazing to vertical, followed by the corresponding
  // -mu values; the original DISORT user list is monotonically increasing.
  constexpr std::array<Index, nquad> user_to_quadrature{3, 2, 0, 1};
  for (const Numeric phi : disort_test::reference::problem_10_azimuth)
    for (const Numeric tau : disort_test::reference::problem_10_output_tau) {
      dis.u(quadrature, tau, phi);
      dis.u_user(formal, tau, phi, disort_test::reference::problem_10_user_mu);
      for (Index angle = 0; angle < nquad; ++angle)
        expect_reference(std::format("Problem 10 formal/quadrature [{}, {}, {}]", phi, tau, angle),
                         formal.intensities[angle],
                         quadrature.intensities[user_to_quadrature[angle]],
                         2e-6);
    }

  Tensor3 bulk(
      disort_test::reference::problem_10_output_tau.size(), disort_test::reference::problem_10_azimuth.size(), nquad);
  dis.ungridded_u(bulk,
                  AscendingGrid{Vector(disort_test::reference::problem_10_output_tau)},
                  disort_test::reference::problem_10_azimuth);
  for (Index level = 0; level < 3; ++level)
    for (Index azimuth = 0; azimuth < 2; ++azimuth) {
      dis.u(quadrature,
            disort_test::reference::problem_10_output_tau[level],
            disort_test::reference::problem_10_azimuth[azimuth]);
      for (Index angle = 0; angle < nquad; ++angle)
        expect_reference(std::format("Problem 10 pointwise/bulk [{}, {}, {}]", level, azimuth, angle),
                         bulk[level, azimuth, angle],
                         quadrature.intensities[angle],
                         2e-12);
    }
}

disort::main_data problem_11_atmosphere(const bool subdivided) {
  constexpr Index nquad  = disort_test::reference::problem_11_streams;
  const Index     layers = subdivided ? 3 : 1;
  Matrix          moments(layers, nquad + 1, 0.0);
  moments[joker, 0] = 1.0;
  Matrix down(nquad, nquad / 2, 0.0);
  down[0] = disort_test::reference::problem_11_top_isotropic;
  std::vector<disort::BDRF> brdf{
      disort::BDRF{[](auto value, auto&, auto&) { value = disort_test::reference::problem_11_surface_albedo; }}};
  return disort::main_data(nquad,
                           nquad,
                           nquad,
                           subdivided ? AscendingGrid{Vector(disort_test::reference::problem_11_subdivided_tau)}
                                      : AscendingGrid{disort_test::reference::problem_11_output_tau.back()},
                           Vector(layers, disort_test::reference::problem_11_omega),
                           std::move(moments),
                           Matrix(nquad, nquad / 2, 0.0),
                           std::move(down),
                           Vector(layers, 0.0),
                           Matrix(layers, 0),
                           std::move(brdf),
                           disort_test::reference::problem_11_beam_mu,
                           disort_test::reference::problem_11_beam / disort_test::reference::problem_11_beam_mu,
                           0.0);
}

void test_problem_11() {
  const auto          one_layer  = problem_11_atmosphere(false);
  const auto          subdivided = problem_11_atmosphere(true);
  disort::user_u_data one_user, subdivided_user;
  disort::flux_data   one_flux, subdivided_flux;

  for (const Numeric phi : disort_test::reference::problem_11_azimuth)
    for (const Numeric tau : disort_test::reference::problem_11_output_tau) {
      one_layer.u_user(one_user, tau, phi, disort_test::reference::problem_11_user_mu);
      subdivided.u_user(subdivided_user, tau, phi, disort_test::reference::problem_11_user_mu);
      for (Index angle = 0; angle < 4; ++angle)
        expect_reference(std::format("Problem 11 radiance [{}, {}, {}]", phi, tau, angle),
                         subdivided_user.intensities[angle],
                         one_user.intensities[angle],
                         2e-10);
    }

  for (const Numeric tau : disort_test::reference::problem_11_output_tau) {
    const auto one  = one_layer.flux(one_flux, tau);
    const auto many = subdivided.flux(subdivided_flux, tau);
    expect_reference(std::format("Problem 11 upward flux [{}]", tau), many.up, one.up, 2e-10);
    expect_reference(std::format("Problem 11 diffuse-down flux [{}]", tau), many.down_diffuse, one.down_diffuse, 2e-10);
    expect_reference(std::format("Problem 11 direct-down flux [{}]", tau), many.down_direct, one.down_direct, 2e-10);
    expect_reference(std::format("Problem 11 DFDT [{}]", tau), many.dfdt, one.dfdt, 2e-10);
  }
}

disort::main_data problem_12_atmosphere(const bool subdivided) {
  constexpr Index nquad  = disort_test::reference::problem_12_streams;
  const Index     layers = subdivided ? 3 : 1;
  Matrix          moments(layers, nquad + 1);
  for (Index layer = 0; layer < layers; ++layer)
    for (Index moment = 0; moment <= nquad; ++moment)
      moments[layer, moment] = std::pow(disort_test::reference::problem_12_asymmetry, moment);
  std::vector<disort::BDRF> brdf{
      disort::BDRF{[](auto value, auto&, auto&) { value = disort_test::reference::problem_12_surface_albedo; }}};
  return disort::main_data(nquad,
                           nquad,
                           nquad,
                           subdivided ? AscendingGrid{Vector(disort_test::reference::problem_12_subdivided_tau)}
                                      : AscendingGrid{disort_test::reference::problem_12_output_tau.back()},
                           Vector(layers, disort_test::reference::problem_12_omega),
                           std::move(moments),
                           Matrix(nquad, nquad / 2, 0.0),
                           Matrix(nquad, nquad / 2, 0.0),
                           Vector(layers, std::pow(disort_test::reference::problem_12_asymmetry, nquad)),
                           Matrix(layers, 0),
                           std::move(brdf),
                           disort_test::reference::problem_12_beam_mu,
                           disort_test::reference::problem_12_beam / disort_test::reference::problem_12_beam_mu,
                           0.0);
}

void test_problem_12() {
  const auto          one_layer  = problem_12_atmosphere(false);
  const auto          subdivided = problem_12_atmosphere(true);
  disort::user_u_data one_user, subdivided_user;
  disort::tms_data    one_tms, subdivided_tms;
  Vector              one_ims, subdivided_ims;
  disort::flux_data   one_flux, subdivided_flux;

  for (const Numeric tau : disort_test::reference::problem_12_output_tau) {
    one_layer.u_user_corr(one_user,
                          one_ims,
                          one_tms,
                          tau,
                          disort_test::reference::problem_12_azimuth,
                          disort_test::reference::problem_12_user_mu);
    subdivided.u_user_corr(subdivided_user,
                           subdivided_ims,
                           subdivided_tms,
                           tau,
                           disort_test::reference::problem_12_azimuth,
                           disort_test::reference::problem_12_user_mu);
    for (Index angle = 0; angle < 4; ++angle)
      expect_reference(std::format("Problem 12 radiance [{}, {}]", tau, angle),
                       subdivided_user.intensities[angle],
                       one_user.intensities[angle],
                       2e-9);

    const auto one  = one_layer.flux(one_flux, tau);
    const auto many = subdivided.flux(subdivided_flux, tau);
    expect_reference(std::format("Problem 12 upward flux [{}]", tau), many.up, one.up, 2e-9);
    expect_reference(std::format("Problem 12 diffuse-down flux [{}]", tau), many.down_diffuse, one.down_diffuse, 2e-9);
    expect_reference(std::format("Problem 12 direct-down flux [{}]", tau), many.down_direct, one.down_direct, 2e-12);
    expect_reference(std::format("Problem 12 DFDT [{}]", tau), many.dfdt, one.dfdt, 2e-9);
  }
}

void test_problem_13_boundary_limits() {
  using namespace disort_test::reference;
  constexpr Numeric depth = 0.25;
  Matrix            moments(1, problem_13_streams + 1, 0.0);
  moments[0, 0] = 1.0;
  std::vector<disort::BDRF> brdf{disort::BDRF{[](auto value, auto&, auto&) { value = problem_13_surface_albedo; }}};
  const disort::main_data   absorbing(problem_13_streams,
                                      problem_13_streams,
                                      1,
                                      AscendingGrid{depth},
                                      Vector{0.0},
                                      std::move(moments),
                                      Matrix(1, problem_13_streams / 2, 0.0),
                                      Matrix(1, problem_13_streams / 2, 0.0),
                                      Vector{0.0},
                                      Matrix(1, 0),
                                      std::move(brdf),
                                      problem_13_beam_mu,
                                      problem_13_beam,
                                      0.0);

  const Numeric incident_flux       = problem_13_beam * problem_13_beam_mu;
  const Numeric bottom_beam         = incident_flux * std::exp(-depth / problem_13_beam_mu);
  Numeric       upward_transmission = 0.0;
  for (Index i = 0; i < problem_13_streams / 2; ++i)
    upward_transmission += absorbing.weights()[i] * absorbing.mu()[i] * std::exp(-depth / absorbing.mu()[i]);
  const Numeric expected_up = 2.0 * problem_13_surface_albedo * bottom_beam * upward_transmission;

  disort::flux_data flux;
  const auto        top    = absorbing.flux(flux, 0.0);
  const auto        bottom = absorbing.flux(flux, depth);
  expect_reference("Problem 13 absorbing Lambertian top reflection", top.up, expected_up, 2e-12);
  expect_reference("Problem 13 absorbing direct transmission", bottom.down_direct, bottom_beam, 2e-12);
  expect_reference("Problem 13 absorbing diffuse transmission", bottom.down_diffuse, 0.0, 2e-12);
}

void test_problem_13() {
  using namespace disort_test::reference;
  for (const auto& test : problem_13) {
    const Index layers = test.cumulative_tau.size();
    Matrix      moments(layers, problem_13_streams + 1);
    for (Index layer = 0; layer < layers; ++layer)
      for (Index moment = 0; moment <= problem_13_streams; ++moment)
        moments[layer, moment] = std::pow(problem_13_asymmetry, moment);
    std::vector<disort::BDRF> brdf{disort::BDRF{[](auto value, auto&, auto&) { value = problem_13_surface_albedo; }}};
    const disort::main_data   dis(problem_13_streams,
                                  problem_13_streams,
                                  1,
                                  AscendingGrid{test.cumulative_tau},
                                  test.single_scattering_albedo,
                                  std::move(moments),
                                  Matrix(1, problem_13_streams / 2, 0.0),
                                  Matrix(1, problem_13_streams / 2, 0.0),
                                  Vector(layers, std::pow(problem_13_asymmetry, problem_13_streams)),
                                  Matrix(layers, 0),
                                  std::move(brdf),
                                  problem_13_beam_mu,
                                  problem_13_beam,
                                  0.0);

    disort::flux_data flux;
    const Numeric     incident_flux = problem_13_beam * problem_13_beam_mu;
    const auto        top           = dis.flux(flux, 0.0);
    const auto        bottom        = dis.flux(flux, test.cumulative_tau.back());
    expect_reference(std::format("{} albedo", test.name), top.up / incident_flux, test.albedo);
    expect_reference(std::format("{} transmission", test.name),
                     (bottom.down_direct + bottom.down_diffuse) / incident_flux,
                     test.transmission);
  }
}

void test_raw_brdfs_and_fourier_helper() {
  constexpr Numeric mu0 = 0.8660254037844386;
  using namespace disort::brdf;
  expect_reference("Hapke raw BRDF", mu0 * Hapke{}(0.1, mu0, 0.0), 0.0519252, 2e-5);
  expect_reference("Cox-Munk raw BRDF", mu0 * CoxMunk{}(0.5, mu0, Constant::pi / 4.0), 0.00374155, 2e-5);
  expect_reference("RPV raw BRDF", mu0 * RPV{}(0.2, mu0, Constant::pi), 0.0668223, 2e-5);
  expect_reference("Ross-Li raw BRDF", mu0 * RossLi{}(0.5, mu0, Constant::pi / 2.0), 0.0660729, 2e-5);

  constexpr Numeric reflectance = 0.37;
  const auto   modes = fourier_modes([](Numeric, Numeric, Numeric) { return reflectance * Constant::inv_pi; }, 4, 32);
  const Vector outgoing{0.2, 0.8};
  const Vector incoming{-0.3, -0.9};
  for (Index mode = 0; mode < 4; ++mode) {
    const Matrix coefficient = modes[mode](outgoing, incoming);
    for (Index i = 0; i < coefficient.nrows(); ++i)
      for (Index j = 0; j < coefficient.ncols(); ++j)
        expect_reference(std::format("constant BRDF Fourier mode {}", mode),
                         coefficient[i, j],
                         mode == 0 ? reflectance : 0.0,
                         2e-13);
  }
}

disort::brdf::RawFunction problem_14_raw(const disort_test::reference::brdf_type type) {
  using enum disort_test::reference::brdf_type;
  switch (type) {
    case hapke:    return disort::brdf::Hapke{};
    case cox_munk: return disort::brdf::CoxMunk{};
    case rpv:      return disort::brdf::RPV{};
    case ross_li:  return disort::brdf::RossLi{};
  }
  std::unreachable();
}

disort::brdf::RawFunction problem_15_raw(const disort_test::reference::brdf_type type) {
  if (type == disort_test::reference::brdf_type::cox_munk) return disort::brdf::CoxMunk{.shadowing = true};
  return problem_14_raw(type);
}

void test_problem_14() {
  using namespace disort_test::reference;
  constexpr Numeric transparent_depth      = 1e-12;
  constexpr Numeric transparent_scattering = 1e-8;

  for (const auto& test : problem_14) {
    const auto raw = problem_14_raw(test.type);
    for (Index azimuth = 0; azimuth < 3; ++azimuth)
      for (Index angle = 0; angle < 4; ++angle)
        expect_reference(std::format("{} raw radiance [{}, {}]", test.name, azimuth, angle),
                         problem_14_beam_mu * raw(problem_14_user_mu[angle], problem_14_beam_mu, test.azimuth[azimuth]),
                         test.radiance[azimuth, angle],
                         2e-5);

    Matrix moments(1, problem_14_streams + 1, 0.0);
    moments[0, 0] = 1.0;
    const disort::main_data dis(problem_14_streams,
                                problem_14_streams,
                                problem_14_streams,
                                AscendingGrid{transparent_depth},
                                Vector{transparent_scattering},
                                std::move(moments),
                                Matrix(problem_14_streams, problem_14_streams / 2, 0.0),
                                Matrix(problem_14_streams, problem_14_streams / 2, 0.0),
                                Vector{0.0},
                                Matrix(1, 0),
                                disort::brdf::fourier_modes(raw, problem_14_streams),
                                problem_14_beam_mu,
                                problem_14_beam,
                                0.0);

    disort::user_u_data user;
    for (Index azimuth = 0; azimuth < 3; ++azimuth) {
      dis.u_user(user, 0.0, test.azimuth[azimuth], problem_14_user_mu);
      for (Index angle = 0; angle < 4; ++angle)
        expect_reference(std::format("{} solver radiance [{}, {}]", test.name, azimuth, angle),
                         user.intensities[angle],
                         test.radiance[azimuth, angle],
                         2e-5);
    }

    disort::flux_data flux;
    const auto        values = dis.flux(flux, 0.0);
    expect_reference(std::format("{} direct flux", test.name), values.down_direct, test.direct);
    expect_reference(std::format("{} diffuse-down flux", test.name), values.down_diffuse, test.diffuse_down);
    expect_reference(std::format("{} upward flux", test.name), values.up, test.up, 2e-5);
    expect_reference(std::format("{} DFDT", test.name), values.dfdt, test.dfdt, 2e-5);
  }
}

void test_problem_15() {
  using namespace disort_test::reference;

  constexpr Index problem_15_moment_count = 600;
  Matrix          moments(2, problem_15_moment_count, 0.0);
  moments[0, 0] = 1.0;
  moments[0, 2] = 0.1;
  moments[1]    = kokhanovsky_aerosol_moments[Range{0, problem_15_moment_count}];

  for (const auto& test : problem_15) {
    const auto              raw = problem_15_raw(test.type);
    const disort::main_data dis(problem_15_streams,
                                problem_15_streams,
                                problem_15_streams,
                                AscendingGrid{problem_15_tau},
                                Vector{1.0, 1.0},
                                moments,
                                Matrix(problem_15_streams, problem_15_streams / 2, 0.0),
                                Matrix(problem_15_streams, problem_15_streams / 2, 0.0),
                                Vector{0.0, kokhanovsky_aerosol_moments[problem_15_streams]},
                                Matrix(2, 0),
                                disort::brdf::fourier_modes(raw, problem_15_streams),
                                problem_15_beam_mu,
                                problem_15_beam,
                                0.0);

    disort::user_u_data user;
    disort::tms_data    tms;
    Vector              ims;
    for (Index azimuth = 0; azimuth < std::ssize(problem_15_azimuth); ++azimuth) {
      for (Index level = 0; level < std::ssize(problem_15_output_tau); ++level) {
        dis.u_user_corr(user, ims, tms, problem_15_output_tau[level], problem_15_azimuth[azimuth], problem_15_user_mu);
        for (Index angle = 0; angle < std::ssize(problem_15_user_mu); ++angle)
          expect_reference(std::format("{} radiance [{}, {}, {}]", test.name, azimuth, level, angle),
                           user.intensities[angle],
                           test.radiance[azimuth, level, angle]);
      }
    }

    disort::flux_data flux;
    for (Index level = 0; level < std::ssize(problem_15_output_tau); ++level) {
      const auto values = dis.flux(flux, problem_15_output_tau[level]);
      expect_reference(std::format("{} direct flux [{}]", test.name, level), values.down_direct, test.direct[level]);
      expect_reference(
          std::format("{} diffuse-down flux [{}]", test.name, level), values.down_diffuse, test.diffuse_down[level]);
      expect_reference(std::format("{} upward flux [{}]", test.name, level), values.up, test.up[level]);
      expect_reference(std::format("{} DFDT [{}]", test.name, level), values.dfdt, test.dfdt[level], 2e-6);
    }
  }
}

void test_problem_17() {
  using namespace disort_test::reference;

  const std::array<const Vector*, 2> phase_moments{
      &kokhanovsky_aerosol_moments,
      &kokhanovsky_cloud_moments,
  };
  const std::array expected_fraction{0.24306, 0.44398};

  for (Index case_index = 0; case_index < std::ssize(problem_17); ++case_index) {
    const auto& test = problem_17[case_index];
    Matrix      moments(1, phase_moments[case_index]->size());
    moments[0]         = *phase_moments[case_index];
    const auto scaling = disort::delta_m_plus(moments, problem_17_streams);

    expect_reference(
        std::format("{} delta-M-plus fraction", test.name), scaling.fraction[0], expected_fraction[case_index], 2e-5);

    const disort::main_data dis(problem_17_streams,
                                problem_17_streams,
                                problem_17_streams,
                                AscendingGrid{test.depth},
                                Vector{1.0},
                                std::move(moments),
                                Matrix(problem_17_streams, problem_17_streams / 2, 0.0),
                                Matrix(problem_17_streams, problem_17_streams / 2, 0.0),
                                scaling.fraction,
                                Matrix(1, 0),
                                {},
                                problem_17_beam_mu,
                                problem_17_beam,
                                0.0,
                                scaling.moments);

    expect_reference(std::format("{} scaled depth", test.name),
                     dis.scaled_tau().back(),
                     test.depth * (1.0 - (1.0 - 1e-8) * scaling.fraction[0]),
                     2e-13);

    disort::user_u_data user;
    const Vector        output_tau{0.0, test.depth};
    for (Index level = 0; level < 2; ++level) {
      for (Index azimuth = 0; azimuth < std::ssize(problem_17_azimuth); ++azimuth) {
        dis.u_user(user, output_tau[level], problem_17_azimuth[azimuth], problem_17_user_mu);
        for (Index angle = 0; angle < std::ssize(problem_17_user_mu); ++angle)
          expect_reference(std::format("{} radiance [{}, {}, {}]", test.name, level, angle, azimuth),
                           user.intensities[angle],
                           test.radiance[level, angle, azimuth]);
      }
    }
  }
}

}  // namespace

int main() try {
  test_problem_1();
  test_problem_2();
  test_problem_3();
  test_problem_4();
  test_problem_5();
  test_problem_8();
  test_problem_9();
  test_problem_10();
  test_problem_11();
  test_problem_12();
  test_problem_13_boundary_limits();
  test_problem_13();
  test_raw_brdfs_and_fourier_helper();
  test_problem_14();
  test_problem_15();
  test_problem_17();
  return EXIT_SUCCESS;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return EXIT_FAILURE;
}
