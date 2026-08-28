#include <arts_constants.h>
#include <disort.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

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
}  // namespace

int main() try {
  test_1a_user_angles();
  test_problem_1();
  test_problem_2();
  return EXIT_SUCCESS;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return EXIT_FAILURE;
}
