#include <arts_constants.h>
#include <lbl_lineshape_voigt_ecs_hadded.h>
#include <wigner_functions.h>

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {
namespace nh3 = lbl::voigt::ecs::hadded;
using enum nh3::inversion;

void require(bool condition, std::string_view message) {
  if (not condition) throw std::runtime_error(std::string(message));
}

void near(
    Numeric actual, Numeric expected, std::string_view message, Numeric relative = 2e-13, Numeric absolute = 2e-14) {
  if (not std::isfinite(actual) or std::abs(actual - expected) > absolute + relative * std::abs(expected))
    throw std::runtime_error(std::format("{}: got {}, expected {}", message, actual, expected));
}

nh3::rotational_line line(Index Ju, Index Jl, Index K, nh3::inversion lower = symmetric) {
  return {{Ju, K, lower == symmetric ? antisymmetric : symmetric}, {Jl, K, lower}};
}

template <class F> void rejects(F&& function) {
  try {
    function();
  } catch (const std::exception&) { return; }
  throw std::runtime_error("Expected invalid NH3 core input to be rejected");
}

void energy_and_dipoles() {
  const Numeric B = 9.9 * Constant::h * Constant::c * 100;
  const Numeric C = 6.3 * Constant::h * Constant::c * 100;
  near(nh3::rotational_energy(3, 2, B, C) / B, 8 + 4 * C / B, "Rigid symmetric-top energy");
  near(nh3::rotational_energy(0, 0, B, C), 0, "Rotational energy zero");

  // Paper I, Eq. 9, evaluated independently as exact Wigner coefficients.
  near(nh3::reduced_dipole(line(1, 1, 1)), -std::sqrt(0.5), "Signed Q(1,1) dipole");
  near(nh3::reduced_dipole(line(2, 1, 1)), std::sqrt(0.5), "Signed R(1,1) dipole");
  near(nh3::reduced_dipole(line(1, 2, 1)), -std::sqrt(0.3), "Signed P(2,1) dipole");
  near(nh3::reduced_dipole(line(4, 4, 3)), -std::sqrt(9.0 / 20), "Signed Q(4,3) dipole");
  near(nh3::reduced_dipole(line(3, 3, 3, antisymmetric)), -std::sqrt(0.75), "Inversion partner dipole");
  rejects([] { (void)nh3::reduced_dipole(line(1, 1, 0)); });
  rejects([] {
    auto forbidden           = line(2, 1, 1);
    forbidden.upper.symmetry = forbidden.lower.symmetry;
    (void)nh3::reduced_dipole(forbidden);
  });
}

void adiabatic_factors() {
  constexpr Numeric duration = 2e-12;
  const Numeric     unit     = Constant::h_bar / duration;
  const Vector      gap{0, 2 * unit, 3 * unit};
  const Vector      dduration{0.3 * duration, -0.2 * duration, 0};
  Matrix            dgap(3, 3);
  dgap[0] = Vector{0.1 * unit, 0.2 * unit, -0.4 * unit};
  dgap[1] = Vector{0, 0, 0};
  dgap[2] = Vector{0.5 * unit, -0.1 * unit, 0.3 * unit};
  Vector omega(3);
  Matrix domega(3, 3);
  nh3::adiabatic_factors(omega, gap, duration, domega, dgap, dduration);
  for (Index i = 0; i < 3; ++i) {
    const Numeric x = gap[i] * duration / Constant::h_bar;
    near(omega[i], std::pow(1 + x * x / 24, 2), "Paper-I Eq. 19");
    for (Index q = 0; q < 3; ++q) {
      const Numeric dx = (dgap[q, i] * duration + gap[i] * dduration[q]) / Constant::h_bar;
      near(domega[q, i], (1 + x * x / 24) * x * dx / 6, "Adiabatic factor derivative");
    }
  }
  Vector sudden(3);
  nh3::adiabatic_factors(sudden, gap, 0);
  for (Numeric value : sudden) near(value, 1, "Zero-duration IOS limit");
}

Matrix matrix(std::span<const nh3::rotational_line> lines,
              const nh3::basis_data&                basis,
              const Vector&                         e0,
              const Vector&                         omega,
              Numeric                               T = 240) {
  Matrix W(lines.size(), lines.size(), -99);
  for (Size i = 0; i < lines.size(); ++i) W[i, i] = Numeric(5 + i);
  nh3::relaxation_matrix_offdiagonal(W, lines, basis, e0, omega, T);
  for (Size i = 0; i < lines.size(); ++i) near(W[i, i], Numeric(5 + i), "Kernel preserves supplied widths", 0, 0);
  return W;
}

void exact_angular_fixtures() {
  const Vector e0{20 * Constant::k, 0}, omega{1, 1};
  {
    const std::array lines{line(1, 1, 1), line(1, 1, 1, antisymmetric)};
    nh3::basis_data  basis{{{1, 0, 0}}, Vector{2.6}, Vector{1}};
    const auto       W = matrix(lines, basis, e0, omega);
    near(W[1, 0], -2.6 / 4, "Same-J inversion coupling, no extra (2L+1)");
    const std::array same_symmetry{line(1, 1, 1), line(2, 2, 1)};
    const auto       forbidden = matrix(same_symmetry, basis, e0, omega);
    near(forbidden[0, 1], 0, "Odd-L same-inversion Q-branch selection rule", 0, 0);
    const std::array spin_species{line(1, 1, 1), line(3, 3, 3)};
    const auto       separate = matrix(spin_species, basis, e0, omega);
    near(separate[0, 1], 0, "Ortho/para blocks do not mix", 0, 0);
    near(separate[1, 0], 0, "Ortho/para reverse block remains zero", 0, 0);
  }

  // Independent factorial Racah sums for Eq. 10. Each test enables exactly
  // one signed channel, including both off-diagonal dynamical factors.
  const auto check_channels = [&](const auto& lines, const auto& channels, const auto& expected) {
    for (Size c = 0; c < channels.size(); ++c) {
      const nh3::basis_data basis{{channels[c]}, Vector{1}, Vector{1}};
      const auto            W = matrix(lines, basis, e0, omega);
      near(W[1, 0], expected[c], "Independent four-term angular coefficient");
    }
  };
  check_channels(
      std::array{line(4, 4, 3), line(3, 3, 3)},
      std::array<nh3::collision_channel, 4>{{{6, 0, 0}, {6, 6, 0}, {6, 0, 6}, {6, 6, 6}}},
      std::array{
          std::sqrt(245.0 / 15704832), -std::sqrt(35.0 / 118976), std::sqrt(35.0 / 118976), -std::sqrt(15.0 / 2704)});
  check_channels(std::array{line(7, 7, 3), line(7, 7, 6)},
                 std::array<nh3::collision_channel, 4>{{{9, -3, -3}, {9, 9, -3}, {9, -3, 9}, {9, 9, 9}}},
                 std::array{-2805.0 / 97888,
                            -std::sqrt(50139375.0 / 67074423808),
                            std::sqrt(50139375.0 / 67074423808),
                            17875.0 / 685216});

  // K=0 has N=1 and epsilon=0, rather than the K>0 normalization.
  check_channels(std::array{line(5, 4, 0, antisymmetric), line(4, 3, 0)},
                 std::array<nh3::collision_channel, 1>{{{1, 0, 0}}},
                 std::array{-2 * std::sqrt(35.0) / 27});
  const Numeric zero_to_three = -std::sqrt(5.0) / 88;
  check_channels(std::array{line(5, 4, 0, antisymmetric), line(4, 3, 3)},
                 std::array<nh3::collision_channel, 4>{{{3, -3, -3}, {3, 3, -3}, {3, -3, 3}, {3, 3, 3}}},
                 std::array{zero_to_three, zero_to_three, zero_to_three, zero_to_three});

  // Equal lower energies: the two optical branches must give the same result
  // regardless of which row appears first in the caller's matrix.
  const std::array      tied{line(4, 3, 3), line(3, 3, 3)};
  const std::array      reversed{tied[1], tied[0]};
  const nh3::basis_data basis{{{6, 0, 0}, {6, 6, 0}, {6, 0, 6}, {6, 6, 6}}, Vector{1, 2, 3, 4}, Vector{1, 1, 1, 1}};
  const Numeric expected = std::sqrt(49.0 / 15704832) + 2 * std::sqrt(343.0 / 118976) + 3 * std::sqrt(7.0 / 118976) +
                           4 * std::sqrt(147.0 / 2704);
  const auto    W        = matrix(tied, basis, Vector{0, 0}, omega);
  const auto    swapped  = matrix(reversed, basis, Vector{0, 0}, omega);
  near(W[0, 1], expected, "Exact-energy tie angular coefficient");
  near(W[0, 1], W[1, 0], "Equal-population detailed balance");
  near(W[0, 1], swapped[1, 0], "Exact-energy tie permutation invariance");
}

void balance_and_batched_derivatives() {
  constexpr Index       n = 4, nq = 6;
  constexpr Numeric     T = 240, step = 1e-4;
  const std::array      lines{line(1, 1, 1), line(1, 1, 1, antisymmetric), line(2, 2, 1), line(3, 2, 2, antisymmetric)};
  const Vector          e0{100 * Constant::k, 102 * Constant::k, 135 * Constant::k, 120 * Constant::k};
  const Vector          omega{1.2, 1.1, 1.7, 1.35};
  const nh3::basis_data basis{
      {{1, 0, 0}, {2, 0, 0}, {3, 3, 3}, {3, 0, 0}}, Vector{1.1, 0.7, 0.4, 0}, Vector{1.03, 1.07, 1.11, 1.2}};
  Vector dT(nq, 0.0);
  Matrix de0(nq, n, 0.0), dQ(nq, 4, 0.0), dOmega_basis(nq, 4, 0.0), dOmega_line(nq, n, 0.0);
  dT[0]           = 1;
  de0[1]          = Vector{0.1 * Constant::k, -0.2 * Constant::k, 0.3 * Constant::k, 0.4 * Constant::k};
  dQ[2, 3]        = 1;  // Derivative must survive when the primal basis rate is zero.
  dOmega_basis[3] = Vector{0.1, -0.2, 0.2, 0.3};
  dOmega_line[4]  = Vector{0.2, -0.1, 0.3, -0.2};
  dT[5]           = -0.5;
  de0[5]          = de0[1];
  dQ[5]           = Vector{0.2, -0.3, 0.1, 0.4};
  dOmega_basis[5] = dOmega_basis[3];
  dOmega_line[5]  = dOmega_line[4];
  const nh3::derivative_data derivatives{dT, de0, dQ, dOmega_basis, dOmega_line};
  Matrix                     W(n, n, 0.0);
  Tensor3                    dW(nq, n, n, 0.0);
  for (Index i = 0; i < n; ++i) {
    W[i, i] = Numeric(5 + i);
    for (Index q = 0; q < nq; ++q) dW[q, i, i] = Numeric((q + 1) * (i + 1)) * 0.01;
  }
  const Tensor3 seeded = dW;
  nh3::relaxation_matrix_offdiagonal(W, lines, basis, e0, omega, T, dW, derivatives);
  const auto plain = matrix(lines, basis, e0, omega, T);
  for (Index i = 0; i < n; ++i) {
    for (Index j = 0; j < n; ++j) near(W[i, j], plain[i, j], "Derivatives preserve primal matrix", 0, 0);
    for (Index q = 0; q < nq; ++q) near(dW[q, i, i], seeded[q, i, i], "Diagonal derivative preserved", 0, 0);
    const Numeric rho_i = Numeric(2 * (2 * lines[i].lower.J + 1)) * std::exp(-e0[i] / (Constant::k * T));
    for (Index j = i + 1; j < n; ++j) {
      const Numeric rho_j = Numeric(2 * (2 * lines[j].lower.J + 1)) * std::exp(-e0[j] / (Constant::k * T));
      // W[i, j] couples into i from j, so the source population symmetrises.
      near(rho_j * W[i, j], rho_i * W[j, i], "Detailed balance with lower-state degeneracy");
    }
  }
  require(std::abs(dW[2, 1, 2]) > 1e-5, "A zero-Q channel must have a nonzero coupling derivative");

  for (Index q = 0; q < nq; ++q) {
    auto evaluate = [&](Numeric amount) {
      auto shifted_basis  = basis;
      auto shifted_energy = e0;
      auto shifted_omega  = omega;
      for (Index i = 0; i < n; ++i) {
        shifted_energy[i] += amount * de0[q, i];
        shifted_omega[i]  += amount * dOmega_line[q, i];
      }
      for (Index c = 0; c < 4; ++c) {
        shifted_basis.Q[c]     += amount * dQ[q, c];
        shifted_basis.Omega[c] += amount * dOmega_basis[q, c];
      }
      return matrix(lines, shifted_basis, shifted_energy, shifted_omega, T + amount * dT[q]);
    };
    const auto high = evaluate(step), low = evaluate(-step);
    for (Index i = 0; i < n; ++i)
      for (Index j = 0; j < n; ++j)
        if (i != j)
          near(
              dW[q, i, j], (high[i, j] - low[i, j]) / (2 * step), "Batched directional finite difference", 3e-7, 2e-10);
  }

  const std::array<Index, n>          order{2, 0, 3, 1};
  std::array<nh3::rotational_line, n> permuted;
  Vector                              pe(n), po(n);
  for (Index i = 0; i < n; ++i) {
    permuted[i] = lines[order[i]];
    pe[i]       = e0[order[i]];
    po[i]       = omega[order[i]];
  }
  const auto permuted_W = matrix(permuted, basis, pe, po, T);
  for (Index i = 0; i < n; ++i)
    for (Index j = 0; j < n; ++j)
      if (i != j) near(permuted_W[i, j], W[order[i], order[j]], "Energy-selected matrix permutation invariance");

  Vector dipr(n);
  for (Index i = 0; i < n; ++i) dipr[i] = nh3::reduced_dipole(lines[i]);
  const auto original             = W;
  const auto original_derivatives = dW;
  nh3::sum_rule_diagonal(W, dipr, dW);
  for (Index i = 0; i < n; ++i) {
    Numeric residual = 0;
    for (Index j = 0; j < n; ++j) {
      residual += dipr[j] * W[i, j];
      if (i != j) near(W[i, j], original[i, j], "Width estimate preserves off-diagonals", 0, 0);
    }
    near(residual, 0, "Optional optical sum-rule row closure");
    for (Index q = 0; q < nq; ++q) {
      Numeric derivative_residual = 0;
      for (Index j = 0; j < n; ++j) {
        derivative_residual += dipr[j] * dW[q, i, j];
        if (i != j)
          near(dW[q, i, j], original_derivatives[q, i, j], "Width estimate preserves coupling derivatives", 0, 0);
      }
      near(derivative_residual, 0, "Differentiated optical sum-rule closure");
    }
  }
}
}  // namespace

int main() try {
  WignerInformation wigner(100, 0, true, true);
  energy_and_dipoles();
  adiabatic_factors();
  exact_angular_fixtures();
  balance_and_batched_derivatives();
  wigner.unload();
  std::cout << "NH3 Hadded core equation tests passed (synthetic collision basis)\n";
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return 1;
}
