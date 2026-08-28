#include <arts_constants.h>
#include <disort.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace {
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
}  // namespace

int main() try {
  test_1a_user_angles();
  return EXIT_SUCCESS;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return EXIT_FAILURE;
}
