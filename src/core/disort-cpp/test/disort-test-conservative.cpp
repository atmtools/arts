#include <arts_constants.h>
#include <disort.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {
constexpr Numeric top_downward  = 0.2;
constexpr Numeric bottom_upward = 0.7;
constexpr Numeric total_depth   = 3.0;

void expect_close(const std::string_view label,
                  const Numeric          actual,
                  const Numeric          expected,
                  const Numeric          tolerance = 5.0e-10) {
  if (not std::isfinite(actual) or std::abs(actual - expected) > tolerance * std::max(1.0, std::abs(expected)))
    throw std::runtime_error(std::format(
        "{}: expected {:.17g}, got {:.17g} (difference {:.3g})", label, expected, actual, actual - expected));
}

Numeric sinhc(const Numeric x) {
  if (std::abs(x) > 1.0e-5) return std::sinh(x) / x;
  const Numeric x2 = x * x;
  return 1.0 + x2 * (1.0 / 6.0 + x2 * (1.0 / 120.0 + x2 / 5040.0));
}

std::array<Numeric, 2> analytic_two_stream(const Numeric tau, const Numeric epsilon) {
  // CPP-DISORT uses double Gauss-Legendre quadrature on each hemisphere;
  // its one positive node is mu=1/2 rather than the full-range node 1/sqrt(3).
  const Numeric p      = 2.0;
  const Numeric lambda = p * std::sqrt(epsilon);
  const Numeric cT     = std::cosh(lambda * total_depth);
  const Numeric hT     = total_depth * sinhc(lambda * total_depth);
  const Numeric aT     = cT + p * epsilon * hT;
  const Numeric bT     = p * hT + cT;
  const Numeric q0     = (2.0 * bottom_upward - 2.0 * top_downward * aT) / (aT + bT);
  const Numeric s0     = q0 + 2.0 * top_downward;

  const Numeric c = std::cosh(lambda * tau);
  const Numeric h = tau * sinhc(lambda * tau);
  const Numeric s = c * s0 + p * h * q0;
  const Numeric q = p * epsilon * h * s0 + c * q0;
  return {0.5 * (s + q), 0.5 * (s - q)};
}

disort::main_data make_model(const AscendingGrid& depths, const Numeric omega, Matrix source = {}) {
  const Index layers = depths.size();
  Matrix      moments(layers, 2, 0.0);
  moments[joker, 0] = 1.0;
  Matrix up(1, 1, bottom_upward), down(1, 1, top_downward);
  if (source.empty()) source = Matrix(layers, 0);
  return disort::main_data(2,
                           2,
                           1,
                           AscendingGrid{depths},
                           Vector(layers, omega),
                           std::move(moments),
                           std::move(up),
                           std::move(down),
                           Vector(layers, 0.0),
                           std::move(source),
                           {},
                           1.0,
                           0.0,
                           0.0);
}

void check_analytic_sweep() {
  constexpr std::array epsilon{0.0, 1.0e-14, 1.0e-12, 1.0e-10, 1.0e-8, 2.0e-8};
  constexpr std::array depths{0.0, 0.3, 1.2, total_depth};
  for (const Numeric eps : epsilon) {
    const auto        dis = make_model(AscendingGrid{total_depth}, 1.0 - eps);
    disort::u0_data   u0;
    disort::u_data    u;
    disort::flux_data flux;
    for (const Numeric tau : depths) {
      const auto expected = analytic_two_stream(tau, eps);
      dis.u0(u0, tau);
      expect_close("u0 upward", u0.u0[0], expected[0]);
      expect_close("u0 downward", u0.u0[1], expected[1]);

      dis.u(u, tau, 0.37);
      expect_close("u upward", u.intensities[0], expected[0]);
      expect_close("u downward", u.intensities[1], expected[1]);

      const auto    values      = dis.flux(flux, tau);
      const Numeric flux_factor = Constant::pi;
      expect_close("upward flux", values.up, flux_factor * expected[0]);
      expect_close("downward flux", values.down_diffuse, flux_factor * expected[1]);
      expect_close("direct flux", values.down_direct, 0.0);
    }

    Vector output_depths(depths.size());
    stdr::copy(depths, output_depths.begin());
    Tensor3 field(depths.size(), 1, 2);
    dis.ungridded_u(field, AscendingGrid{std::move(output_depths)}, Vector{0.37});
    for (Index i = 0; i < static_cast<Index>(depths.size()); ++i) {
      const auto expected = analytic_two_stream(depths[static_cast<std::size_t>(i)], eps);
      expect_close("ungridded upward", field[i, 0, 0], expected[0]);
      expect_close("ungridded downward", field[i, 0, 1], expected[1]);
    }
  }
}

void check_layer_splitting() {
  const auto      one   = make_model(AscendingGrid{total_depth}, 1.0);
  const auto      split = make_model(AscendingGrid{0.4, 1.1, total_depth}, 1.0);
  disort::u0_data one_u0, split_u0;
  for (const Numeric tau : {0.0, 0.3, 0.8, 1.7, total_depth}) {
    one.u0(one_u0, tau);
    split.u0(split_u0, tau);
    for (Index state = 0; state < 2; ++state)
      expect_close("exact-conservative layer splitting", split_u0.u0[state], one_u0.u0[state], 2.0e-11);
  }
}

void check_zero_emission_and_repeat_update() {
  Matrix source(1, 2);
  source[0, 0]                = 123.0;
  source[0, 1]                = -45.0;
  auto            with_source = make_model(AscendingGrid{total_depth}, 1.0, std::move(source));
  auto            no_source   = make_model(AscendingGrid{total_depth}, 1.0);
  disort::u0_data source_u0, reference_u0;
  for (const Numeric tau : {0.0, 0.9, total_depth}) {
    with_source.u0(source_u0, tau);
    no_source.u0(reference_u0, tau);
    for (Index state = 0; state < 2; ++state)
      expect_close("zero conservative emission", source_u0.u0[state], reference_u0.u0[state], 2.0e-11);
  }

  no_source.u0(reference_u0, 1.3);
  const Vector before = reference_u0.u0;
  no_source.update_all();
  expect_close("omega retained by update", no_source.omega()[0], 1.0, 0.0);
  no_source.u0(reference_u0, 1.3);
  for (Index state = 0; state < 2; ++state)
    expect_close("repeat update", reference_u0.u0[state], before[state], 2.0e-11);
}
}  // namespace

int main() try {
  check_analytic_sweep();
  check_layer_splitting();
  check_zero_emission_and_repeat_update();
  return 0;
} catch (const std::exception& error) {
  std::cerr << "Conservative DISORT test failed: " << error.what() << '\n';
  return 1;
}
