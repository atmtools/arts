#include <arts_constants.h>
#include <arts_conversions.h>

#include <cmath>
#include <iostream>
#include <stdexcept>

#include "tro_legendre.h"

namespace {
using namespace scattering;

void check(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error("FAILED: " + what);
  std::cout << "PASSED: " << what << '\n';
}

Numeric ylm0_norm(Index l) { return std::sqrt(static_cast<Numeric>(2 * l + 1) / (4.0 * Constant::pi)); }

/* Data that are linear in the scattering angle are their own interpolant, so
   the projection must give the closed-form coefficients to rounding:
   F = Theta on the grid {0, 180} has b_0 = int_0^pi Theta sin(Theta) dTheta = pi
   and b_1 = int_0^pi Theta cos(Theta) sin(Theta) dTheta = -pi / 4. */
void test_linear_is_exact() {
  const Vector angles{0.0, 180.0};
  Matrix       values(2, 1);
  values[0, 0] = 0.0;
  values[1, 0] = Constant::pi;
  const Matrix a = tro_legendre::project(angles, values, 1);
  const Numeric b0 = a[0, 0] / (2.0 * Constant::pi * ylm0_norm(0)), b1 = a[1, 0] / (2.0 * Constant::pi * ylm0_norm(1));
  check(std::abs(b0 - Constant::pi) < 1e-14 and std::abs(b1 + 0.25 * Constant::pi) < 1e-14,
        "projection of a function linear in the scattering angle is exact");
}

/* Constant ends: data only on [30, 150] deg are held constant beyond, so a
   constant F = 1 there is F = 1 everywhere: a_0 = sqrt(4 pi), all else 0. */
void test_constant_ends() {
  const Vector angles{30.0, 90.0, 150.0};
  Matrix       values(3, 2, 1.0);
  const Matrix a  = tro_legendre::project(angles, values, 12);
  Numeric      e  = std::abs(a[0, 0] - std::sqrt(4.0 * Constant::pi));
  for (Index l = 1; l <= 12; l++) e = std::max(e, std::abs(a[l, 1]));
  check(e < 1e-14, "the constant extension beyond the first and last angle");
}

/* evaluate() is the series sum_l a_l Y_l0: Henyey-Greenstein has a_l =
   sqrt((2 l + 1) / 4 pi) g^l, so its truncated series is a known sum. */
void test_evaluate() {
  const Index   L = 40;
  const Numeric g = 0.6;
  Matrix        a(L + 1, 1);
  for (Index l = 0; l <= L; l++) a[l, 0] = ylm0_norm(l) * std::pow(g, l);
  const Vector angles{0.0, 33.0, 90.0, 147.0, 180.0};
  const Matrix f = tro_legendre::evaluate(a, angles);
  Numeric      e = 0.0;
  for (Size i = 0; i < angles.size(); i++) {
    const Numeric x = std::cos(Conversion::deg2rad(angles[i]));
    Numeric       p0 = 1.0, p1 = x, ref = 1.0 / (4.0 * Constant::pi) * (1.0 + 3.0 * g * x);
    for (Index l = 1; l < L; l++) {
      const auto    dl = static_cast<Numeric>(l);
      const Numeric p2 = ((2.0 * dl + 1.0) * x * p1 - dl * p0) / (dl + 1.0);
      ref += (2.0 * dl + 3.0) / (4.0 * Constant::pi) * std::pow(g, l + 1) * p2;
      p0 = p1;
      p1 = p2;
    }
    e = std::max(e, std::abs(f[i, 0] - ref) / ref);
  }
  check(e < 1e-13, "evaluation of a Legendre series");
}

/* Henyey-Greenstein on a fine grid: the projection converges to the analytic
   coefficients as the grid refines (second order in the spacing), and the
   report sees the truncation. */
void test_hg_convergence() {
  const Numeric g = 0.7;
  const Index   L = 16;
  Numeric       previous = 0.0;
  for (const Index n : {181, 721}) {
    Vector angles(n);
    Matrix values(n, 6, 0.0);
    for (Index i = 0; i < n; i++) {
      angles[i]     = 180.0 * static_cast<Numeric>(i) / static_cast<Numeric>(n - 1);
      const Numeric x = std::cos(Conversion::deg2rad(angles[i]));
      values[i, 0]  = (1.0 - g * g) / (4.0 * Constant::pi * std::pow(1.0 + g * g - 2.0 * g * x, 1.5));
    }
    const Matrix a = tro_legendre::project(angles, values, L);
    Numeric      e = 0.0;
    for (Index l = 0; l <= L; l++) e = std::max(e, std::abs(a[l, 0] - ylm0_norm(l) * std::pow(g, l)));
    if (n == 721) check(e < previous / 10.0, "projection of gridded HG converges at second order");
    previous = e;

    if (n == 721) {
      const auto r = tro_legendre::assess(a, angles, values);
      check(std::abs(r.asymmetry - g) < 1e-5, "the report's asymmetry parameter");
      check(std::abs(r.tail[0] - std::pow(g, L) * ylm0_norm(L) / ylm0_norm(0)) < 1e-5, "the report's tail");
      check(r.reconstruction_error[0] > 1e-4, "the report sees the truncation of a forward-peaked function");
      check(r.reconstruction_error[1] == 0.0, "the report of an element that is zero");
    }
  }
}

void test_errors() {
  const auto throws = [](auto&& f) {
    try {
      f();
    } catch (const std::exception&) {
      return true;
    }
    return false;
  };
  Matrix v(2, 1, 1.0);
  check(throws([&] { (void)tro_legendre::project(Vector{0.0, 180.0}, v, -1); }), "a negative degree throws");
  check(throws([&] { (void)tro_legendre::project(Vector{90.0, 90.0}, v, 2); }), "repeated angles throw");
  check(throws([&] { (void)tro_legendre::project(Vector{0.0, 190.0}, v, 2); }), "angles beyond 180 deg throw");
  check(throws([&] { (void)tro_legendre::project(Vector{0.0, 90.0, 180.0}, v, 2); }), "a shape mismatch throws");
}
}  // namespace

int main() try {
  test_linear_is_exact();
  test_constant_ends();
  test_evaluate();
  test_hg_convergence();
  test_errors();
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
