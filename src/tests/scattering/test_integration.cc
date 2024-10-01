/** \file scattering/test_integration.cc
 *
 * Unit tests for the quadrature classes defined in scattering/integration.h
 *
 * @author Simon Pfreundschuh, 2023
 */
#include <iostream>
#include <numbers>

#include "integration.h"
#include "matpack/matpack_math.h"

using std::numbers::pi_v;

template <matpack::any_matpack_type IN>
Numeric max_error(const IN& a, const IN& b) {
  auto delta = a;
  delta -= b;
  std::transform(delta.begin(), delta.end(), delta.begin(), [](Numeric x) {
    return fabs(x);
  });
  return max(delta);
}

bool test_gauss_legendre() {
  scattering::GaussLegendreQuadrature quad(3);

  Vector nodes = quad.get_nodes();
  Vector weights = quad.get_weights();
  Vector nodes_ref = {-sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0)};
  Vector weights_ref = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

  auto error = max_error(nodes, nodes_ref);
  if (error > 1e-6) return false;
  error = max_error(weights, weights_ref);
  if (error > 1e-6) return false;

  return true;
}

/** Test Gauss Legendre quadrature.
 *
 * Ensure that Gauss-Legendre nodes and weights are correct for the
 * case n = 3.
 *
 */
bool test_double_gauss() {
  scattering::DoubleGaussQuadrature quad(6);
  Vector nodes = quad.get_nodes();
  Vector weights = quad.get_weights();
  Vector nodes_ref = {-1.0 + 0.5 * (-sqrt(3.0 / 5.0) + 1.0),
                      -0.5,
                      -1.0 + 0.5 * (sqrt(3.0 / 5.0) + 1.0),
                      1.0 - 0.5 * (sqrt(3.0 / 5.0) + 1.0),
                      0.5,
                      1.0 - 0.5 * (-sqrt(3.0 / 5.0) + 1.0)};
  Vector weights_ref = {0.5 * 5.0 / 9.0,
                        0.5 * 8.0 / 9.0,
                        0.5 * 5.0 / 9.0,
                        0.5 * 5.0 / 9.0,
                        0.5 * 8.0 / 9.0,
                        0.5 * 5.0 / 9.0};

  auto error = max_error(nodes, nodes_ref);
  if (error > 1e-6) return false;
  error = max_error(weights, weights_ref);
  if (error > 1e-6) return false;

  return true;
}

bool test_lobatto_quadrature() {
  scattering::LobattoQuadrature quad(4);
  Vector nodes = quad.get_nodes();
  Vector weights = quad.get_weights();
  Vector nodes_ref = {-1.0, -sqrt(5.0 / 25.0), sqrt(5.0 / 25.0), 1.0};
  Vector weights_ref = {1.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0};

  auto error = max_error(nodes, nodes_ref);
  if (error > 1e-6) return false;
  error = max_error(weights, weights_ref);
  if (error > 1e-6) return false;

  return true;
}

bool test_clenshaw_curtis_quadrature() {
  scattering::ClenshawCurtisQuadrature quad(4);
  Vector nodes = quad.get_nodes();
  Vector weights = quad.get_weights();

  Vector nodes_ref = {-1.0, -0.5, 0.5, 1.0};
  auto error = max_error(nodes, nodes_ref);
  if (error > 1e-6) return false;

  return true;
}

bool test_fejer_quadrature() {
  scattering::FejerQuadrature quad(4);
  Vector nodes = quad.get_nodes();
  Vector weights = quad.get_weights();

  Vector nodes_ref = {-cos(0.125 * pi_v<Numeric>),
                      -cos(0.375 * pi_v<Numeric>),
                      -cos(0.625 * pi_v<Numeric>),
                      -cos(0.875 * pi_v<Numeric>)};

  auto error = max_error(nodes, nodes_ref);
  if (error > 1e-6) return false;

  return true;
}

bool test_zenith_angle_integration() {
  scattering::GaussLegendreGrid grid{2};
  auto nodes = grid.get_angle_cosines();

  // GL quadrature of degree 2 must be exact for integrationg cos(x)^2.
  Vector y = {pow(nodes[0], 2.0), pow(nodes[1], 2.0)};

  Numeric y_int = scattering::integrate_zenith_angle(y, grid);
  if (fabs(y_int - 2.0 / 3.0) > 1e-6) return false;
  return true;
}

bool test_calculate_downsampling_weights() {
  Vector old_grid{{0.0, 0.3, 0.7, 1.3, 1.9, 2.1, 2.22, 2.8, pi_v<Numeric>}};
  Vector new_grid{{0.0, 0.9, 1.7, pi_v<Numeric>}};

  auto weights =
      scattering::calculate_downsampling_matrix<Numeric>(old_grid, new_grid);

  Vector y(old_grid);
  std::transform(
      y.begin(), y.end(), y.begin(), [](Numeric x) { return sin(x); });
  Vector y_sub(4);
  mult(y_sub, weights, y);

  Numeric y_int_ref = 0.0;
  for (Index i = 0; i < y.size() - 1; ++i) {
    y_int_ref += 0.5 * (y[i] + y[i + 1]) * (old_grid[i + 1] - old_grid[i]);
  }

  Numeric y_int = 0.0;
  for (Index i = 0; i < y_sub.size() - 1; ++i) {
    y_int += 0.5 * (y_sub[i] + y_sub[i + 1]) * (new_grid[i + 1] - new_grid[i]);
  }

  if (std::abs(y_int - y_int_ref) > 1e-6) {
    return false;
  }
  return true;
}

int main(int /*argc*/, const char** /*argv*/) {
#ifndef ARTS_NO_SHTNS

  bool passed = true;

  std::cout << "Testing Gauss-Legendre quadrature: ";
  passed &= test_gauss_legendre();
  if (passed) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
    return 1;
  }

  std::cout << "Testing double Gauss-Legendre quadrature: ";
  passed &= test_double_gauss();
  if (passed) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
    return 1;
  }

  std::cout << "Testing Lobatto quadrature: ";
  passed &= test_lobatto_quadrature();
  if (passed) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
    return 1;
  }

  std::cout << "Testing Clenshaw-Curtis quadrature: ";
  passed &= test_clenshaw_curtis_quadrature();
  if (passed) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
    return 1;
  }

  std::cout << "Testing Fejer quadrature: ";
  passed &= test_fejer_quadrature();
  if (passed) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
    return 1;
  }

  std::cout << "Testing zenith-angle integration: ";
  passed &= test_zenith_angle_integration();
  if (passed) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
    return 1;
  }

  std::cout << "Testing zenith-angle integration: ";
  passed &= test_calculate_downsampling_weights();
  if (passed) {
    std::cout << "PASSED" << std::endl;
  } else {
    std::cout << "FAILED" << std::endl;
    return 1;
  }
#endif

  return 0;
};
