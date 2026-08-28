#include <arts_constants.h>
#include <disort.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

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

struct isotropic_case {
  std::string_view name;
  Numeric          depth;
  Numeric          omega;
  bool             beam;
  Vector           direct;
  Vector           diffuse_down;
  Vector           up;
  Matrix           radiance;
};

void run_isotropic_case(const isotropic_case& test) {
  constexpr Index nquad = 16;
  Matrix          legendre(1, 17, 0.0);
  legendre[0, 0] = 1.0;
  Matrix down(nquad, nquad / 2, 0.0);
  if (not test.beam) down[0] = 1.0;

  const Numeric mu0 = 0.1;
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

  const Vector user_mu{-1.0, -0.5, -0.1, 0.1, 0.5, 1.0};
  const Vector tau{0.0, test.depth};
  disort::user_u_data user;
  disort::flux_data   flux;
  for (Index level = 0; level < 2; ++level) {
    dis.u_user(user, tau[level], 0.0, user_mu);
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

void test_1a_user_angles() {
  constexpr Index nquad = 16;
  Matrix          legendre(1, 17, 0.0);
  legendre[0, 0] = 1.0;

  const Numeric mu0 = 0.1;
  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{0.03125},
                              Vector{0.2},
                              legendre,
                              Matrix(nquad, nquad / 2, 0.0),
                              Matrix(nquad, nquad / 2, 0.0),
                              Vector{0.0},
                              Matrix(1, 0),
                              {},
                              mu0,
                              Constant::pi / mu0,
                              0.0);

  const Vector user_mu{-1.0, -0.5, -0.1, 0.1, 0.5, 1.0};
  const Matrix expected{Vector{0.0,
                               0.0,
                               0.0,
                               1.17771e-1,
                               2.64170e-2,
                               1.34041e-2,
                               1.33826e-2,
                               2.63324e-2,
                               1.15898e-1,
                               0.0,
                               0.0,
                               0.0}
                            .reshape(2, 6)};

  disort::user_u_data data;
  Matrix              actual(2, 6);
  dis.u_user(data, 0.0, 0.0, user_mu);
  actual[0] = data.intensities;
  dis.u_user(data, 0.03125, 0.0, user_mu);
  actual[1] = data.intensities;

  auto value     = actual.elem_begin();
  auto reference = expected.elem_begin();
  for (Size i = 0; i < actual.size(); ++i, ++value, ++reference) {
    ARTS_USER_ERROR_IF(std::abs(*value - *reference) > 7e-6 * std::max(1.0, std::abs(*reference)),
                       "DISORT 4.0.99 test 1a user-angle mismatch at {}: {} vs {}",
                       i,
                       *value,
                       *reference);
  }
}

void test_problem_1() {
  const auto radiance = [](std::initializer_list<Numeric> values) {
    Vector data(values.size());
    std::ranges::copy(values, data.begin());
    return Matrix{std::move(data).reshape(2, 6)};
  };
  run_isotropic_case({"1a",
                      0.03125,
                      0.2,
                      true,
                      {3.14159, 2.29844},
                      {0.0, 7.94108e-2},
                      {7.99451e-2, 0.0},
                      radiance({0, 0, 0, 1.17771e-1, 2.64170e-2, 1.34041e-2,
                                1.33826e-2, 2.63324e-2, 1.15898e-1, 0, 0, 0})});
  run_isotropic_case({"1b",
                      0.03125,
                      1.0,
                      true,
                      {3.14159, 2.29844},
                      {0.0, 4.20233e-1},
                      {4.22922e-1, 0.0},
                      radiance({0, 0, 0, 6.22884e-1, 1.39763e-1, 7.09192e-2,
                                7.08109e-2, 1.39337e-1, 6.13458e-1, 0, 0, 0})});
  run_isotropic_case({"1c",
                      0.03125,
                      0.99,
                      false,
                      {0, 0},
                      {3.14159, 3.04897},
                      {9.06556e-2, 0},
                      radiance({1, 1, 1, 1.33177e-1, 2.99879e-2, 1.52233e-2,
                                9.84447e-1, 9.69363e-1, 8.63946e-1, 0, 0, 0})});
  run_isotropic_case({"1d",
                      32.0,
                      0.2,
                      true,
                      {3.14159, 0},
                      {0, 0},
                      {2.59686e-1, 0},
                      radiance({0, 0, 0, 2.62972e-1, 9.06967e-2, 5.02853e-2,
                                1.22980e-15, 1.30698e-17, 6.88840e-18, 0, 0, 0})});
  run_isotropic_case({"1e",
                      32.0,
                      1.0,
                      true,
                      {3.14159, 0},
                      {0, 6.76954e-2},
                      {3.07390, 0},
                      radiance({0, 0, 0, 1.93321, 1.02732, 7.97199e-1,
                                2.71316e-2, 1.87805e-2, 1.16385e-2, 0, 0, 0})});
  run_isotropic_case({"1f",
                      32.0,
                      0.99,
                      false,
                      {0, 0},
                      {3.14159, 4.60048e-3},
                      {2.49618, 0},
                      radiance({1, 1, 1, 8.77510e-1, 8.15136e-1, 7.52715e-1,
                                1.86840e-3, 1.26492e-3, 7.79280e-4, 0, 0, 0})});
}

void run_rayleigh_case(const std::string_view name,
                       const Numeric          depth,
                       const Numeric          omega,
                       const Vector&          direct,
                       const Vector&          diffuse_down,
                       const Vector&          up,
                       const Matrix&          radiance) {
  constexpr Index nquad = 16;
  Matrix          legendre(1, 17, 0.0);
  legendre[0, 0] = 1.0;
  legendre[0, 2] = 0.1;
  constexpr Numeric mu0 = 0.080442;
  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{depth},
                              Vector{omega},
                              legendre,
                              Matrix(nquad, nquad / 2, 0.0),
                              Matrix(nquad, nquad / 2, 0.0),
                              Vector{0.0},
                              Matrix(1, 0),
                              {},
                              mu0,
                              Constant::pi,
                              0.0);
  const Vector user_mu{-0.981986, -0.538263, -0.018014, 0.018014, 0.538263, 0.981986};
  disort::user_u_data user;
  disort::flux_data   flux;
  for (Index level = 0; level < 2; ++level) {
    const Numeric tau = level == 0 ? 0.0 : depth;
    dis.u_user(user, tau, 0.0, user_mu);
    for (Index angle = 0; angle < 6; ++angle)
      expect_reference(std::format("{} radiance [{}, {}]", name, level, angle),
                       user.intensities[angle],
                       radiance[level, angle]);
    const auto [down, beam] = dis.flux_down(flux, tau);
    expect_reference(std::format("{} direct flux [{}]", name, level), beam, direct[level]);
    expect_reference(std::format("{} diffuse-down flux [{}]", name, level), down, diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", name, level), dis.flux_up(flux, tau), up[level]);
  }
}

void test_problem_2() {
  const auto matrix = [](std::initializer_list<Numeric> values) {
    Vector data(values.size());
    std::ranges::copy(values, data.begin());
    return Matrix{std::move(data).reshape(2, 6)};
  };
  run_rayleigh_case("2a", 0.2, 0.5, {2.52716e-1, 2.10311e-2}, {0, 4.41791e-2}, {5.35063e-2, 0},
                    matrix({0, 0, 0, 1.61796e-1, 2.11501e-2, 7.86713e-3,
                            7.71897e-3, 2.00778e-2, 2.57685e-2, 0, 0, 0}));
  run_rayleigh_case("2b", 0.2, 1.0, {2.52716e-1, 2.10311e-2}, {0, 1.06123e-1}, {1.25561e-1, 0},
                    matrix({0, 0, 0, 3.47678e-1, 4.87120e-2, 1.89387e-2,
                            1.86027e-2, 4.64061e-2, 6.77603e-2, 0, 0, 0}));
  run_rayleigh_case("2c", 5.0, 0.5, {2.52716e-1, 2.56077e-28}, {0, 2.51683e-4}, {6.24730e-2, 0},
                    matrix({0, 0, 0, 1.62566e-1, 2.45786e-2, 1.01498e-2,
                            1.70004e-4, 3.97168e-5, 1.32472e-5, 0, 0, 0}));
  run_rayleigh_case("2d", 5.0, 1.0, {2.52716e-1, 0}, {0, 2.68008e-2}, {2.25915e-1, 0},
                    matrix({0, 0, 0, 3.64010e-1, 8.26993e-2, 4.92370e-2,
                            1.05950e-2, 7.69337e-3, 3.79276e-3, 0, 0, 0}));
}

void run_henyey_greenstein_case(const std::string_view name,
                                const Numeric depth,
                                const Vector& direct,
                                const Vector& diffuse_down,
                                const Vector& up,
                                const Matrix& radiance) {
  constexpr Index nquad = 16;
  Matrix legendre(1, 33, 0.0);
  for (Index i = 0; i < 33; ++i) legendre[0, i] = std::pow(0.75, i);

  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{depth},
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
  const Vector user_mu{-1.0, -0.5, -0.1, 0.1, 0.5, 1.0};
  disort::user_u_data user;
  disort::tms_data    tms;
  disort::flux_data   flux;
  Vector              ims;
  for (Index level = 0; level < 2; ++level) {
    const Numeric tau = level == 0 ? 0.0 : depth;
    dis.u_user_corr(user, ims, tms, tau, 0.0, user_mu);
    for (Index angle = 0; angle < 6; ++angle)
      expect_reference(std::format("{} radiance [{}, {}]", name, level, angle),
                       user.intensities[angle],
                       radiance[level, angle]);
    const auto [down, beam] = dis.flux_down(flux, tau);
    expect_reference(std::format("{} direct flux [{}]", name, level), beam, direct[level]);
    expect_reference(std::format("{} diffuse-down flux [{}]", name, level), down, diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", name, level), dis.flux_up(flux, tau), up[level]);
  }
}

void test_problem_3() {
  const auto matrix = [](std::initializer_list<Numeric> values) {
    Vector data(values.size());
    std::ranges::copy(values, data.begin());
    return Matrix{std::move(data).reshape(2, 6)};
  };
  run_henyey_greenstein_case("3a", 1.0,
                             {3.14159, 1.15573}, {0.0, 1.73849}, {0.247374, 0.0},
                             matrix({0, 0, 0, 0.151159, 0.101103, 0.0395460,
                                     3.05855, 0.266648, 0.213750, 0, 0, 0}));
  run_henyey_greenstein_case("3b", 8.0,
                             {3.14159, 0.00105389}, {0.0, 1.54958}, {1.59096, 0.0},
                             matrix({0, 0, 0, 0.379740, 0.519598, 0.493302,
                                     0.669581, 0.422350, 0.236362, 0, 0, 0}));
}

Matrix haze_l_moments() {
  const Vector haze_l{2.41260, 3.23047, 3.37296, 3.23150, 2.89350, 2.49594, 2.11361, 1.74812,
                      1.44692, 1.17714, 0.96643, 0.78237, 0.64114, 0.51966, 0.42563, 0.34688,
                      0.28351, 0.23317, 0.18963, 0.15788, 0.12739, 0.10762, 0.08597, 0.07381,
                      0.05828, 0.05089, 0.03971, 0.03524, 0.02720, 0.02451, 0.01874, 0.01711};
  Matrix moments(1, 33, 0.0);
  moments[0, 0] = 1.0;
  for (Index k = 1; k <= 32; ++k)
    moments[0, k] = haze_l[k - 1] / static_cast<Numeric>(2 * k + 1);
  return moments;
}

void run_haze_l_case(const std::string_view name,
                     const Numeric omega,
                     const Numeric mu0,
                     const Vector& phi,
                     const Vector& direct,
                     const Vector& diffuse_down,
                     const Vector& up,
                     const Tensor3& radiance) {
  constexpr Index nquad = 32;
  const Matrix moments = haze_l_moments();
  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{1.0},
                              Vector{omega},
                              moments,
                              Matrix(nquad, nquad / 2, 0.0),
                              Matrix(nquad, nquad / 2, 0.0),
                              Vector{moments[0, nquad]},
                              Matrix(1, 0),
                              {},
                              mu0,
                              Constant::pi,
                              0.0);
  const Vector user_mu{-1.0, -0.5, -0.1, 0.1, 0.5, 1.0};
  const Vector tau{0.0, 0.5, 1.0};
  disort::user_u_data user;
  disort::tms_data tms;
  disort::flux_data flux;
  Vector ims;
  const Index nphi = phi.size();
  for (Index p = 0; p < nphi; ++p)
    for (Index level = 0; level < 3; ++level) {
      dis.u_user_corr(user, ims, tms, tau[level], phi[p], user_mu);
      for (Index angle = 0; angle < 6; ++angle)
        expect_reference(std::format("{} radiance [{}, {}, {}]", name, p, level, angle),
                         user.intensities[angle],
                         radiance[p, level, angle]);
    }
  for (Index level = 0; level < 3; ++level) {
    const auto [down, beam] = dis.flux_down(flux, tau[level]);
    expect_reference(std::format("{} direct flux [{}]", name, level), beam, direct[level]);
    expect_reference(std::format("{} diffuse-down flux [{}]", name, level), down, diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", name, level), dis.flux_up(flux, tau[level]), up[level]);
  }
}

void test_problem_4() {
  const auto tensor = [](const Index nphi, std::initializer_list<Numeric> values) {
    Vector data(values.size());
    std::ranges::copy(values, data.begin());
    return Tensor3{std::move(data).reshape(nphi, 3, 6)};
  };
  run_haze_l_case("4a", 1.0, 1.0, {0.0},
                  {3.14159, 1.90547, 1.15573}, {0.0, 1.17401, 1.81264}, {0.173223, 0.111113, 0.0},
                  tensor(1, {0, 0, 0, 0.0926837, 0.0659569, 0.0364755,
                             2.51608, 0.119287, 0.134962, 0.123887, 0.0402058, 0.0177746,
                             3.37302, 0.219835, 0.156893, 0, 0, 0}));
  run_haze_l_case("4b", 0.9, 1.0, {0.0},
                  {3.14159, 1.90547, 1.15573}, {0.0, 1.01517, 1.51554}, {0.123665, 0.0788690, 0.0},
                  tensor(1, {0, 0, 0, 0.0653056, 0.0455144, 0.0282693,
                             2.24258, 0.0966049, 0.0961335, 0.0843278, 0.0279473, 0.0138835,
                             2.97057, 0.167698, 0.108115, 0, 0, 0}));
  run_haze_l_case("4c", 0.9, 0.5, {0.0, Constant::pi / 2.0, Constant::pi},
                  {1.57080, 0.577864, 0.212584}, {0.0, 0.702764, 0.803294}, {0.225487, 0.123848, 0.0},
                  tensor(3, {
                      0, 0, 0, 0.870812, 0.224960, 0.0227572,
                      0.0477016, 3.02631, 1.41195, 0.697692, 0.109130, 0.00932861,
                      0.0838488, 2.70538, 0.876523, 0, 0, 0,
                      0, 0, 0, 0.0888117, 0.0577411, 0.0227572,
                      0.0477016, 0.0580971, 0.104502, 0.0916071, 0.0295842, 0.00932861,
                      0.0838488, 0.0942187, 0.0895457, 0, 0, 0,
                      0, 0, 0, 0.0698247, 0.0502877, 0.0227572,
                      0.0477016, 0.0258544, 0.0625954, 0.0591273, 0.0247702, 0.00932861,
                      0.0838488, 0.0399383, 0.0467155, 0, 0, 0}));
}

void run_cloud_c1_case(const disort_test::reference::scalar_case& test) {
  constexpr Index nquad = 48;
  const Matrix moments = disort_test::reference::cloud_c1_moments();
  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{64.0},
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
  const Vector user_mu{-1.0, -0.5, -0.1, 0.1, 0.5, 1.0};
  disort::user_u_data user;
  disort::tms_data tms;
  disort::flux_data flux;
  Vector ims;
  for (Index level = 0; level < 3; ++level) {
    dis.u_user_corr(user, ims, tms, test.tau[level], 0.0, user_mu);
    for (Index angle = 0; angle < 6; ++angle)
      expect_reference(std::format("{} radiance [{}, {}]", test.name, level, angle),
                       user.intensities[angle],
                       test.radiance[level, angle],
                       2e-4);
    const auto [down, beam] = dis.flux_down(flux, test.tau[level]);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), beam, test.direct[level]);
    expect_reference(std::format("{} diffuse-down flux [{}]", test.name, level), down, test.diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), dis.flux_up(flux, test.tau[level]), test.up[level]);
  }
}

void test_problem_5() {
  run_cloud_c1_case(disort_test::reference::problem_5a);
  run_cloud_c1_case(disort_test::reference::problem_5b);
}

void run_problem_8_case(const disort_test::reference::layered_isotropic_case& test) {
  constexpr Index nquad = 8;
  Matrix moments(2, nquad + 1, 0.0);
  moments[joker, 0] = 1.0;

  Matrix boundary_down(nquad, nquad / 2, 0.0);
  boundary_down[0] = Constant::inv_pi;

  const disort::main_data dis(
      nquad,
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
  disort::flux_data flux;
  for (Index level = 0; level < 3; ++level) {
    const Numeric tau = test.output_tau[level];
    dis.u_user(user,
               tau,
               disort_test::reference::problem_8_azimuth,
               disort_test::reference::problem_8_user_mu);
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
  const Numeric scale = Constant::h * Constant::c * 100.0 / (Constant::k * temperature);
  const Numeric x0 = scale * wavenumber_low;
  const Numeric x1 = scale * wavenumber_high;
  const Numeric dx = (x1 - x0) / intervals;
  const auto integrand = [](const Numeric x) {
    if (x == 0.0 || x > 700.0) return 0.0;
    return x * x * x / std::expm1(x);
  };
  Numeric integral = integrand(x0) + integrand(x1);
  for (Index i = 1; i < intervals; ++i)
    integral += (i % 2 ? 4.0 : 2.0) * integrand(x0 + static_cast<Numeric>(i) * dx);
  integral *= dx / 3.0;
  return Constant::sigma * std::pow(temperature, 4) * Constant::inv_pi * integral /
         (Math::pow2(Constant::pi) * Math::pow2(Constant::pi) / 15.0);
}

Matrix problem_9_moments(const disort_test::reference::phase_type phase, const Index nquad = 8) {
  Matrix moments(6, nquad + 1, 0.0);
  if (phase == disort_test::reference::phase_type::isotropic) {
    moments[joker, 0] = 1.0;
  } else if (phase == disort_test::reference::phase_type::problem_9b) {
    const Vector coefficients{1.0, 2.00916 / 3.0, 1.56339 / 5.0, 0.67407 / 7.0,
                              0.22215 / 9.0, 0.04725 / 11.0, 0.00671 / 13.0,
                              0.00068 / 15.0, 0.00005 / 17.0};
    for (Index layer = 0; layer < 6; ++layer) moments[layer] = coefficients;
  } else {
    for (Index layer = 0; layer < 6; ++layer)
      for (Index moment = 0; moment <= nquad; ++moment)
        moments[layer, moment] = std::pow(static_cast<Numeric>(layer + 1) / 7.0, moment);
  }
  return moments;
}

Matrix problem_9_source(const disort_test::reference::general_multilayer_case& test) {
  if (not test.thermal) return Matrix(6, 0);
  Vector planck(7);
  for (Index level = 0; level < 7; ++level)
    planck[level] = band_blackbody_radiance(test.interface_temperature[level],
                                           test.wavenumber_low,
                                           test.wavenumber_high);
  Matrix source(6, 2);
  Numeric tau0 = 0.0;
  for (Index layer = 0; layer < 6; ++layer) {
    const Numeric tau1 = test.cumulative_tau[layer];
    const Numeric slope = (planck[layer + 1] - planck[layer]) / (tau1 - tau0);
    source[layer, 0] = planck[layer] - slope * tau0;
    source[layer, 1] = slope;
    tau0 = tau1;
  }
  return source;
}

void run_problem_9_case(const disort_test::reference::general_multilayer_case& test) {
  constexpr Index nquad = 8;
  Matrix up(nquad, nquad / 2, 0.0), down(nquad, nquad / 2, 0.0);
  down[0] = test.top_isotropic;
  std::vector<disort::BDRF> brdf;
  if (test.thermal) {
    down[0] += band_blackbody_radiance(test.top_temperature, test.wavenumber_low, test.wavenumber_high);
    up[0] = (1.0 - test.surface_albedo) *
            band_blackbody_radiance(test.bottom_temperature, test.wavenumber_low, test.wavenumber_high);
  }
  if (test.surface_albedo != 0.0)
    brdf.push_back(disort::BDRF{[albedo = test.surface_albedo](auto value, auto&, auto&) { value = albedo; }});

  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{1.0, 3.0, 6.0, 10.0, 15.0, 21.0},
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
  disort::tms_data tms;
  disort::flux_data flux;
  Vector ims;
  for (Index azimuth = 0; azimuth < static_cast<Index>(test.azimuth.size()); ++azimuth)
    for (Index level = 0; level < 5; ++level) {
      if (test.phase == disort_test::reference::phase_type::henyey_greenstein)
        dis.u_user_corr(user,
                        ims,
                        tms,
                        test.output_tau[level],
                        test.azimuth[azimuth],
                        disort_test::reference::problem_9_user_mu);
      else
        dis.u_user(user,
                   test.output_tau[level],
                   test.azimuth[azimuth],
                   disort_test::reference::problem_9_user_mu);
      for (Index angle = 0; angle < 4; ++angle)
        expect_reference(std::format("{} radiance [{}, {}, {}]", test.name, azimuth, level, angle),
                         user.intensities[angle],
                         test.radiance[azimuth, level, angle],
                         test.thermal ? 1e-3 : 7e-5);
    }

  for (Index level = 0; level < 5; ++level) {
    const Numeric tau = test.output_tau[level];
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
  constexpr Index nquad = 4;
  const auto& test = disort_test::reference::problem_9[2];
  Matrix up(nquad, nquad / 2, 0.0), down(nquad, nquad / 2, 0.0);
  down[0] = test.top_isotropic +
            band_blackbody_radiance(test.top_temperature, test.wavenumber_low, test.wavenumber_high);
  up[0] = (1.0 - test.surface_albedo) *
          band_blackbody_radiance(test.bottom_temperature, test.wavenumber_low, test.wavenumber_high);
  std::vector<disort::BDRF> brdf{
      disort::BDRF{[albedo = test.surface_albedo](auto value, auto&, auto&) { value = albedo; }}};

  const disort::main_data dis(nquad,
                              nquad,
                              nquad,
                              AscendingGrid{1.0, 3.0, 6.0, 10.0, 15.0, 21.0},
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

  disort::u_data quadrature;
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

  Tensor3 bulk(disort_test::reference::problem_10_output_tau.size(),
               disort_test::reference::problem_10_azimuth.size(),
               nquad);
  dis.ungridded_u(bulk,
                  AscendingGrid{0.0, 2.1, 21.0},
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
}  // namespace

int main() try {
  test_1a_user_angles();
  test_problem_1();
  test_problem_2();
  test_problem_3();
  test_problem_4();
  test_problem_5();
  test_problem_8();
  test_problem_9();
  test_problem_10();
  return EXIT_SUCCESS;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return EXIT_FAILURE;
}
