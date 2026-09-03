#include <lin_alg.h>
#include <rtepack.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

namespace {
constexpr Numeric tol           = 2e-10;
constexpr Numeric stability_tol = 5e-13;

struct MatrixFunctionReference {
  Matrix exp;
  Matrix expm1;
  Matrix phi1;
};

MatrixFunctionReference matrix_function_reference(const Propmat& generator) {
  const Matrix g = to_matrix(generator);

  Matrix exp(4, 4, 0.0);
  Matrix expm1(4, 4, 0.0);
  Matrix phi1(4, 4, 0.0);
  Matrix exp_term(4, 4, 0.0);
  Matrix phi1_term(4, 4, 0.0);
  Matrix product(4, 4);
  for (Index i = 0; i < 4; ++i) {
    exp[i, i] = exp_term[i, i] = 1.0;
    phi1[i, i] = phi1_term[i, i] = 1.0;
  }

  // All callers keep ||generator|| small enough that this direct series is
  // an accurate reference independent of rtepack's spectral coefficients.
  for (Index n = 1; n <= 80; ++n) {
    mult(product, exp_term, g);
    exp_term  = product;
    exp_term /= static_cast<Numeric>(n);
    exp      += exp_term;
    expm1    += exp_term;

    mult(product, phi1_term, g);
    phi1_term  = product;
    phi1_term /= static_cast<Numeric>(n + 1);
    phi1      += phi1_term;
  }

  return {std::move(exp), std::move(expm1), std::move(phi1)};
}

Specmat to_specmat(const Propmat& k) {
  const Muelmat m = rtepack::to_muelmat(k);
  return {m[0, 0],
          m[0, 1],
          m[0, 2],
          m[0, 3],
          m[1, 0],
          m[1, 1],
          m[1, 2],
          m[1, 3],
          m[2, 0],
          m[2, 1],
          m[2, 2],
          m[2, 3],
          m[3, 0],
          m[3, 1],
          m[3, 2],
          m[3, 3]};
}

Numeric max_abs_diff(const Specmat& a, const Specmat& b) {
  Numeric out{};
  for (Size i = 0; i < 4; ++i) {
    for (Size j = 0; j < 4; ++j) {
      const Numeric error = std::abs(a[i, j] - b[i, j]);
      if (not std::isfinite(error)) return std::numeric_limits<Numeric>::infinity();
      out = std::max(out, error);
    }
  }
  return out;
}

Numeric max_abs(const Specmat& a) {
  Numeric out{};
  for (Size i = 0; i < 4; ++i) {
    for (Size j = 0; j < 4; ++j) out = std::max(out, std::abs(a[i, j]));
  }
  return out;
}

Numeric max_abs_diff(const Matrix& a, const Matrix& b) {
  Numeric out{};
  for (Index i = 0; i < a.nrows(); ++i) {
    for (Index j = 0; j < a.ncols(); ++j) out = std::max(out, std::abs(a[i, j] - b[i, j]));
  }
  return out;
}

Numeric max_abs(const Matrix& a) {
  Numeric out{};
  for (Index i = 0; i < a.nrows(); ++i) {
    for (Index j = 0; j < a.ncols(); ++j) out = std::max(out, std::abs(a[i, j]));
  }
  return out;
}

bool close(const Muelmat& a, const Muelmat& b, const Numeric tolerance = tol) {
  return max_abs_diff(Matrix{a}, Matrix{b}) < tolerance;
}

Numeric max_abs_diff(const Stokvec& a, const Stokvec& b) {
  return std::max({std::abs(a.I() - b.I()), std::abs(a.Q() - b.Q()), std::abs(a.U() - b.U()), std::abs(a.V() - b.V())});
}

Muelmat numerical_magnus_deriv(
    const Propmat& k0, const Propmat& k1, const Propmat& dk, const Numeric r, const Numeric dr, const bool k0_deriv) {
  constexpr Numeric eps = 1e-6;
  const Propmat     dk0 = k0_deriv ? dk : Propmat{};
  const Propmat     dk1 = k0_deriv ? Propmat{} : dk;
  const Muelmat     tp = rtepack::tran{k0 + eps * dk0, k1 + eps * dk1, r + eps * dr, rtepack::MagnusOperator::magnus}();
  const Muelmat     tm = rtepack::tran{k0 - eps * dk0, k1 - eps * dk1, r - eps * dr, rtepack::MagnusOperator::magnus}();
  return (tp - tm) / (2.0 * eps);
}

Muelmat numerical_magnus_linsrc_deriv(
    const Propmat& k0, const Propmat& k1, const Propmat& dk, const Numeric r, const Numeric dr, const bool k0_deriv) {
  constexpr Numeric   eps = 1e-6;
  const Propmat       dk0 = k0_deriv ? dk : Propmat{};
  const Propmat       dk1 = k0_deriv ? Propmat{} : dk;
  const Propmat       k0p = k0 + eps * dk0;
  const Propmat       k1p = k1 + eps * dk1;
  const Propmat       k0m = k0 - eps * dk0;
  const Propmat       k1m = k1 - eps * dk1;
  const Numeric       rp  = r + eps * dr;
  const Numeric       rm  = r - eps * dr;
  const rtepack::tran trp{k0p, k1p, rp, rtepack::MagnusOperator::magnus};
  const rtepack::tran trm{k0m, k1m, rm, rtepack::MagnusOperator::magnus};
  return (trp.magnus_linsrc(k0p, k1p, rp) - trm.magnus_linsrc(k0m, k1m, rm)) / (2.0 * eps);
}

Muelmat numerical_linprop_linsrc_deriv(
    const Propmat& k0, const Propmat& k1, const Propmat& dk, const Numeric r, const Numeric dr, const bool k0_deriv) {
  constexpr Numeric   eps = 1e-6;
  const Propmat       dk0 = k0_deriv ? dk : Propmat{};
  const Propmat       dk1 = k0_deriv ? Propmat{} : dk;
  const Propmat       k0p = k0 + eps * dk0;
  const Propmat       k1p = k1 + eps * dk1;
  const Propmat       k0m = k0 - eps * dk0;
  const Propmat       k1m = k1 - eps * dk1;
  const Numeric       rp  = r + eps * dr;
  const Numeric       rm  = r - eps * dr;
  const rtepack::tran trp{k0p, k1p, rp};
  const rtepack::tran trm{k0m, k1m, rm};
  const Muelmat       tp = trp();
  const Muelmat       tm = trm();
  return (trp.linsrc_linprop(tp, k0p, k1p, rp) - trm.linsrc_linprop(tm, k0m, k1m, rm)) / (2.0 * eps);
}

bool test_exponential() {
  constexpr Numeric r = 0.7;
  const Propmat     k0{0.8, 0.03, -0.02, 0.01, 0.04, -0.015, 0.025};
  const Propmat     k1{1.1, -0.01, 0.05, 0.02, -0.03, 0.035, -0.02};

  const Matrix m0 = to_matrix(k0);
  const Matrix m1 = to_matrix(k1);
  Matrix       m1m0(4, 4);
  Matrix       m0m1(4, 4);
  Matrix       omega(4, 4);
  mult(m1m0, m1, m0);
  mult(m0m1, m0, m1);
  for (Index i = 0; i < 4; ++i) {
    for (Index j = 0; j < 4; ++j) {
      omega[i, j] = -0.5 * r * (m0[i, j] + m1[i, j]) + (r * r / 12.0) * (m1m0[i, j] - m0m1[i, j]);
    }
  }
  Matrix reference(4, 4, 0.0);
  Matrix term(4, 4, 0.0);
  Matrix product(4, 4);
  for (Index i = 0; i < 4; ++i) reference[i, i] = term[i, i] = 1.0;
  for (Index n = 1; n <= 30; ++n) {
    mult(product, term, omega);
    term       = product;
    term      /= static_cast<Numeric>(n);
    reference += term;
  }

  const rtepack::tran tr{k0, k1, r, rtepack::MagnusOperator::magnus};
  const Muelmat       actual          = tr();
  const Matrix        generator       = to_matrix(Propmat{tr.a, tr.b, tr.c, tr.d, tr.u, tr.v, tr.w});
  const Numeric       generator_error = max_abs_diff(generator, omega);
  if (generator_error >= tol) {
    std::cerr << "Magnus generator does not match the matrix commutator: " << generator_error << '\n';
    return false;
  }
  const Numeric exponential_error = max_abs_diff(Matrix{actual}, reference);
  if (exponential_error >= tol) {
    std::cerr << "Magnus exponential does not match the general matrix exponential: " << exponential_error << '\n';
    return false;
  }

  // A matrix commutes with its own perturbation direction.  These identities
  // exercise the analytical Frechet derivatives without a finite-difference
  // step or another matrix-function implementation:
  //   D exp(G)[G]   = G exp(G)
  //   D phi_1(G)[G] = exp(G) - phi_1(G).
  const Propmat exponent{tr.a, tr.b, tr.c, tr.d, tr.u, tr.v, tr.w};
  const Muelmat phi1 = tr.linsrc();
  if (not close(tr.deriv(actual, exponent), rtepack::to_muelmat(exponent) * actual, stability_tol) or
      not close(tr.linsrc_deriv(exponent), actual - phi1, stability_tol)) {
    std::cerr << "Matrix-function derivatives violate their commuting-direction identities\n";
    return false;
  }

  constexpr Size ordered_steps = 4096;
  Muelmat        ordered       = Muelmat::id();
  for (Size i = 0; i < ordered_steps; ++i) {
    const Numeric t    = (static_cast<Numeric>(i) + 0.5) / static_cast<Numeric>(ordered_steps);
    const Propmat kmid = k0 + t * (k1 - k0);
    ordered            = rtepack::tran{kmid, kmid, r / static_cast<Numeric>(ordered_steps)}() * ordered;
  }
  const Muelmat reverse_magnus = rtepack::tran{k1, k0, r, rtepack::MagnusOperator::magnus}();
  const Numeric ordered_error  = max_abs_diff(Matrix{actual}, Matrix{ordered});
  const Numeric reverse_error  = max_abs_diff(Matrix{reverse_magnus}, Matrix{ordered});
  if (ordered_error >= reverse_error) {
    std::cerr << "Magnus commutator has the wrong propagation order\n";
    return false;
  }

  const Stokvec j0{0.8, 0.02, 0.01, -0.01};
  const Stokvec j1{1.1, -0.01, 0.03, 0.015};
  const Stokvec incoming{2.0, 0.1, -0.05, 0.02};
  Stokvec       source_reference = incoming;
  for (Size i = 0; i < ordered_steps; ++i) {
    const Numeric t    = (static_cast<Numeric>(i) + 0.5) / static_cast<Numeric>(ordered_steps);
    const Propmat kmid = k0 + t * (k1 - k0);
    const Stokvec jmid = j0 + t * (j1 - j0);
    const Muelmat step = rtepack::tran{kmid, kmid, r / static_cast<Numeric>(ordered_steps)}();
    source_reference   = step * (source_reference - jmid) + jmid;
  }
  const Muelmat       magnus_l      = tr.magnus_linsrc(k0, k1, r);
  const Stokvec       source_magnus = j1 + actual * (incoming - j0) + magnus_l * (j0 - j1);
  const rtepack::tran reverse_tr{k1, k0, r, rtepack::MagnusOperator::magnus};
  const Stokvec       source_reverse =
      j1 + reverse_magnus * (incoming - j0) + reverse_tr.magnus_linsrc(k1, k0, r) * (j0 - j1);
  if (max_abs_diff(source_magnus, source_reference) >= max_abs_diff(source_reverse, source_reference)) {
    std::cerr << "Magnus linear-source correction has the wrong propagation order\n";
    return false;
  }

  const Muelmat constant        = rtepack::tran{k0, k0, r}();
  const Muelmat magnus_constant = rtepack::tran{k0, k0, r, rtepack::MagnusOperator::magnus}();
  if (not close(constant, magnus_constant)) {
    std::cerr << "Magnus does not reduce to the constant case\n";
    return false;
  }

  const Propmat proportional     = 1.7 * k0;
  const Muelmat average          = rtepack::tran{k0, proportional, r}();
  const Muelmat magnus_commuting = rtepack::tran{k0, proportional, r, rtepack::MagnusOperator::magnus}();
  if (not close(average, magnus_commuting)) {
    std::cerr << "Magnus does not reduce to the averaged exponential for commuting matrices\n";
    return false;
  }

  // Clustered polarized eigenvalues make the closed divided differences in
  // phi_1 prone to cancellation.  Compare the stable fallback to an
  // independent power series.
  const Propmat                 clustered{0.3, 1.1e-4, -0.9e-4, 0.7e-4, -1.0e-4, 0.8e-4, -0.6e-4};
  const MatrixFunctionReference clustered_reference = matrix_function_reference(-clustered);
  const rtepack::tran           clustered_tr{clustered, clustered, 1.0};
  if (max_abs_diff(Matrix{clustered_tr()}, clustered_reference.exp) >= stability_tol or
      max_abs_diff(Matrix{clustered_tr.expm1()}, clustered_reference.expm1) >= stability_tol or
      max_abs_diff(Matrix{clustered_tr.linsrc()}, clustered_reference.phi1) >= stability_tol) {
    std::cerr << "Matrix functions are unstable near clustered eigenvalues\n";
    return false;
  }

  // Check both sides of the small-eigenvalue and coefficient-series
  // boundaries.  The b/w cases have x=|b| and y=|w| exactly, while the
  // b/u cases are strongly non-normal and have nearly coalescing eigenvalues.
  constexpr Numeric          boundary_delta    = 1e-6;
  constexpr Numeric          eigen_boundary    = 1e-4;
  const Numeric              coefficient_below = std::sqrt(0.5e-4 * (1.0 - boundary_delta));
  const Numeric              coefficient_above = std::sqrt(0.5e-4 * (1.0 + boundary_delta));
  const std::vector<Propmat> boundary_cases{
      {0.3, eigen_boundary * (1.0 - boundary_delta), 0.0, 0.0, 0.0, 0.0, 0.1},
      {0.3, eigen_boundary * (1.0 + boundary_delta), 0.0, 0.0, 0.0, 0.0, 0.1},
      {0.3, 0.1, 0.0, 0.0, 0.0, 0.0, eigen_boundary * (1.0 - boundary_delta)},
      {0.3, 0.1, 0.0, 0.0, 0.0, 0.0, eigen_boundary * (1.0 + boundary_delta)},
      {0.3, coefficient_below, 0.0, 0.0, 0.0, 0.0, coefficient_below},
      {0.3, coefficient_above, 0.0, 0.0, 0.0, 0.0, coefficient_above},
      {1.0, 1.0, 0.0, 0.0, std::sqrt(1.0 + 1e-6), 0.0, 0.0},
      {1.0, 1.0, 0.0, 0.0, std::sqrt(1.0 + 2e-3), 0.0, 0.0},
  };

  for (Size i = 0; i < boundary_cases.size(); ++i) {
    const Propmat&                k   = boundary_cases[i];
    const MatrixFunctionReference ref = matrix_function_reference(-k);
    const rtepack::tran           boundary_tr{k, k, 1.0};
    const Numeric                 exp_error   = max_abs_diff(Matrix{boundary_tr()}, ref.exp);
    const Numeric                 expm1_error = max_abs_diff(Matrix{boundary_tr.expm1()}, ref.expm1);
    const Numeric                 phi1_error  = max_abs_diff(Matrix{boundary_tr.linsrc()}, ref.phi1);
    if (exp_error >= stability_tol or expm1_error >= stability_tol or phi1_error >= stability_tol) {
      std::cerr << "Matrix-function boundary case " << i << " is unstable: " << exp_error << ", " << expm1_error << ", "
                << phi1_error << '\n';
      return false;
    }
  }

  // An absolute tolerance would allow an implementation to return zero for
  // exp(G)-I when G is tiny.  Require accuracy relative to the remainder.
  const Propmat                 tiny{3e-12, 1.1e-12, -0.9e-12, 0.7e-12, -1.0e-12, 0.8e-12, -0.6e-12};
  const MatrixFunctionReference tiny_reference = matrix_function_reference(-tiny);
  const rtepack::tran           tiny_tr{tiny, tiny, 1.0};
  const Numeric                 tiny_expm1_error = max_abs_diff(Matrix{tiny_tr.expm1()}, tiny_reference.expm1);
  if (tiny_expm1_error >= stability_tol * max_abs(tiny_reference.expm1)) {
    std::cerr << "Matrix expm1 loses its tiny polarized remainder: " << tiny_expm1_error << '\n';
    return false;
  }

  // At high optical depth exp(a) and the polarized hyperbolic coefficients
  // can be individually extreme even though their product is tiny.  In this
  // regime exp(G)-I must approach -I without cancellation between two large
  // terms in its diagonal construction.
  const Propmat opaque{115.91402277, 46.18724517, -43.08403181, 58.86233400, -39.00367487, -25.67156451, 25.61306909};
  const rtepack::tran opaque_tr{opaque, opaque, 1.0};
  Matrix              opaque_reference{opaque_tr()};
  for (Index i = 0; i < 4; ++i) opaque_reference[i, i] -= 1.0;
  if (max_abs_diff(Matrix{opaque_tr.expm1()}, opaque_reference) >= stability_tol) {
    std::cerr << "Matrix expm1 is unstable at high optical depth\n";
    return false;
  }

  // This generator previously selected the cancellation-safe diagonal
  // expression merely because exp(a) C0 - 1 was small.  Its individual a and
  // C0-1 terms are large, however, so that expression rounds to zero or -1.
  // The result must agree with exp(G)-I irrespective of which diagonal form is
  // selected.
  const Propmat       expm1_selector_generator{-67.16886217313498,
                                               34.558419206478604,
                                               82.16181435011583,
                                               33.043707618338715,
                                               -130.31572316043608,
                                               90.53558666731178,
                                               44.637457236401126};
  const rtepack::tran expm1_selector_tr{-expm1_selector_generator, -expm1_selector_generator, 1.0};
  Matrix              expm1_selector_reference{expm1_selector_tr()};
  for (Index i = 0; i < 4; ++i) expm1_selector_reference[i, i] -= 1.0;
  if (max_abs_diff(Matrix{expm1_selector_tr.expm1()}, expm1_selector_reference) >= stability_tol) {
    std::cerr << "Matrix expm1 selects an unstable diagonal expression\n";
    return false;
  }

  // exp(-720 I + 719 N), N^2=I in the leading 2x2 block, has a finite
  // transmitted mode exp(-1) even though exp(-720) underflows on some
  // platforms.  Check exp, expm1, and phi_1 against their analytic modes.
  const Propmat       extreme_generator{-720.0, 719.0};
  const rtepack::tran extreme_tr{-extreme_generator, -extreme_generator, 1.0};
  const Numeric       exp_plus            = std::exp(-1.0);
  const Numeric       exp_minus           = std::exp(-1439.0);
  const Numeric       exp_scalar          = std::exp(-720.0);
  const Numeric       exp_half_sum        = 0.5 * (exp_plus + exp_minus);
  const Numeric       exp_half_difference = 0.5 * (exp_plus - exp_minus);
  const Muelmat       extreme_exp_reference{exp_half_sum,
                                            exp_half_difference,
                                            0.0,
                                            0.0,
                                            exp_half_difference,
                                            exp_half_sum,
                                            0.0,
                                            0.0,
                                            0.0,
                                            0.0,
                                            exp_scalar,
                                            0.0,
                                            0.0,
                                            0.0,
                                            0.0,
                                            exp_scalar};
  Muelmat             extreme_expm1_reference = extreme_exp_reference;
  for (Size i = 0; i < 4; ++i) extreme_expm1_reference[i, i] -= 1.0;

  const auto    phi1_scalar_reference = [](const Numeric z) { return std::expm1(z) / z; };
  const Numeric phi_plus              = phi1_scalar_reference(-1.0);
  const Numeric phi_minus             = phi1_scalar_reference(-1439.0);
  const Numeric phi_scalar            = phi1_scalar_reference(-720.0);
  const Numeric phi_half_sum          = 0.5 * (phi_plus + phi_minus);
  const Numeric phi_half_difference   = 0.5 * (phi_plus - phi_minus);
  const Muelmat extreme_phi1_reference{phi_half_sum,
                                       phi_half_difference,
                                       0.0,
                                       0.0,
                                       phi_half_difference,
                                       phi_half_sum,
                                       0.0,
                                       0.0,
                                       0.0,
                                       0.0,
                                       phi_scalar,
                                       0.0,
                                       0.0,
                                       0.0,
                                       0.0,
                                       phi_scalar};

  if (not close(extreme_tr(), extreme_exp_reference, stability_tol) or
      not close(extreme_tr.expm1(), extreme_expm1_reference, stability_tol) or
      not close(extreme_tr.linsrc(), extreme_phi1_reference, stability_tol)) {
    std::cerr << "Matrix functions lose the finite mode of an extreme dichroic generator\n";
    return false;
  }

  // Repeated squaring must not erase a finite eigenmode merely because the
  // other mode is many orders of magnitude more opaque.  This scale is
  // intentionally far beyond the point where exp(a) and cosh(x) can be
  // formed separately.
  constexpr Numeric   huge_depth = 1e12;
  const Propmat       huge_generator{-huge_depth, huge_depth - 1.0};
  const rtepack::tran huge_tr{-huge_generator, -huge_generator, 1.0};
  const Numeric       huge_fast_mode = std::exp(-1.0);
  const Numeric       huge_slow_phi  = -std::expm1(-2.0 * huge_depth + 1.0) / (2.0 * huge_depth - 1.0);
  const Numeric       huge_fast_phi  = -std::expm1(-1.0);
  const Muelmat       huge_dexp_reference{-0.5 * huge_fast_mode,
                                          -0.5 * huge_fast_mode,
                                          0.0,
                                          0.0,
                                          -0.5 * huge_fast_mode,
                                          -0.5 * huge_fast_mode,
                                          0.0,
                                          0.0,
                                          0.0,
                                          0.0,
                                          0.0,
                                          0.0,
                                          0.0,
                                          0.0,
                                          0.0,
                                          0.0};
  const Muelmat       huge_exp   = huge_tr();
  const Muelmat       huge_phi   = huge_tr.linsrc();
  const Muelmat       huge_dexp  = huge_tr.deriv(huge_exp, huge_generator);
  const Muelmat       huge_dphi1 = huge_tr.linsrc_deriv(huge_generator);
  if (std::abs(huge_exp[0, 0] - 0.5 * huge_fast_mode) > 2e-14 or
      std::abs(huge_exp[0, 1] - 0.5 * huge_fast_mode) > 2e-14 or
      std::abs(huge_phi[0, 0] - 0.5 * (huge_fast_phi + huge_slow_phi)) > 2e-14 or
      std::abs(huge_phi[0, 1] - 0.5 * (huge_fast_phi - huge_slow_phi)) > 2e-14 or
      not close(huge_dexp, huge_dexp_reference, 2e-12) or not close(huge_dphi1, huge_exp - huge_phi, 2e-12)) {
    std::cerr << "Matrix functions erase a finite mode in a very opaque dichroic layer: " << huge_exp[0, 0] << ", "
              << huge_exp[0, 1] << ", " << huge_phi[0, 0] << ", " << huge_phi[0, 1] << "; derivative errors "
              << max_abs_diff(Matrix{huge_dexp}, Matrix{huge_dexp_reference}) << ", "
              << max_abs_diff(Matrix{huge_dphi1}, Matrix{huge_exp - huge_phi}) << '\n';
    return false;
  }

  // Reconstructing the +x eigenspace perturbation as da+dx loses digits
  // when both terms are O(depth) but their surviving eigenvalue is -10.
  // Applying the perturbation between projectors retains this gap directly.
  constexpr Numeric   finite_gap = 10.0;
  const Propmat       gap_generator{-huge_depth, huge_depth - finite_gap};
  const rtepack::tran gap_tr{-gap_generator, -gap_generator, 1.0};
  const Muelmat       gap_exp  = gap_tr();
  const Muelmat       gap_phi  = gap_tr.linsrc();
  const Numeric       gap_mode = -0.5 * finite_gap * std::exp(-finite_gap);
  const Muelmat       gap_dexp_reference{
      gap_mode, gap_mode, 0.0, 0.0, gap_mode, gap_mode, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  if (not close(gap_tr.deriv(gap_exp, gap_generator), gap_dexp_reference, 5e-13) or
      not close(gap_tr.linsrc_deriv(gap_generator), gap_exp - gap_phi, 5e-13)) {
    std::cerr << "Opaque Frechet derivatives lose a finite eigenvalue gap\n";
    return false;
  }

  // A unit rotation is not absolutely small, but its eigenvalue is tiny
  // relative to this layer's 1e12 real separation.  The near-real expansion
  // must use that relative scale and preserve the surviving real mode.
  const Propmat       near_real_generator{-huge_depth, huge_depth - 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
  const rtepack::tran near_real_tr{-near_real_generator, -near_real_generator, 1.0};
  const Muelmat       near_real_exp   = near_real_tr();
  const Muelmat       near_real_phi   = near_real_tr.linsrc();
  const Muelmat       near_real_dexp  = near_real_tr.deriv(near_real_exp, near_real_generator);
  const Muelmat       near_real_dphi1 = near_real_tr.linsrc_deriv(near_real_generator);
  if (not close(near_real_dexp, huge_dexp_reference, 5e-13) or
      not close(near_real_dphi1, near_real_exp - near_real_phi, 5e-13)) {
    std::cerr << "Opaque nearly-real Frechet derivatives lose the finite eigenmode: "
              << max_abs_diff(Matrix{near_real_dexp}, Matrix{huge_dexp_reference}) << ", "
              << max_abs_diff(Matrix{near_real_dphi1}, Matrix{near_real_exp - near_real_phi}) << '\n';
    return false;
  }

  // Also use a noncommuting direction.  In the zero-splitting limit its C
  // component couples the +x and -x real modes to the repeated a mode.  The
  // divided-difference reference is analytic and independent of rtepack's
  // Cayley-Hamilton derivative coefficients.
  const Propmat near_real_noncommuting_direction{0.0, 0.0, 1.0};
  const Numeric real_separation  = huge_depth - 1.0;
  const Numeric phi_at_fast      = phi1_scalar_reference(-1.0);
  const Numeric phi_at_center    = phi1_scalar_reference(-huge_depth);
  const Numeric phi_at_slow      = phi1_scalar_reference(-2.0 * huge_depth + 1.0);
  const Numeric exp_plus_center  = huge_fast_mode / real_separation;
  const Numeric phi_plus_center  = (phi_at_fast - phi_at_center) / real_separation;
  const Numeric phi_minus_center = (phi_at_slow - phi_at_center) / -real_separation;
  const Muelmat near_real_noncommuting_dexp_reference{0.0,
                                                      0.0,
                                                      0.5 * exp_plus_center,
                                                      0.0,
                                                      0.0,
                                                      0.0,
                                                      0.5 * exp_plus_center,
                                                      0.0,
                                                      0.5 * exp_plus_center,
                                                      0.5 * exp_plus_center,
                                                      0.0,
                                                      0.0,
                                                      0.0,
                                                      0.0,
                                                      0.0,
                                                      0.0};
  const Muelmat near_real_noncommuting_dphi_reference{0.0,
                                                      0.0,
                                                      0.5 * (phi_plus_center + phi_minus_center),
                                                      0.0,
                                                      0.0,
                                                      0.0,
                                                      0.5 * (phi_plus_center - phi_minus_center),
                                                      0.0,
                                                      0.5 * (phi_plus_center + phi_minus_center),
                                                      0.5 * (phi_plus_center - phi_minus_center),
                                                      0.0,
                                                      0.0,
                                                      0.0,
                                                      0.0,
                                                      0.0,
                                                      0.0};
  const Muelmat near_real_noncommuting_dexp = near_real_tr.deriv(near_real_exp, near_real_noncommuting_direction);
  const Muelmat near_real_noncommuting_dphi = near_real_tr.linsrc_deriv(near_real_noncommuting_direction);
  if (not close(near_real_noncommuting_dexp, near_real_noncommuting_dexp_reference, 5e-24) or
      not close(near_real_noncommuting_dphi, near_real_noncommuting_dphi_reference, 5e-24)) {
    std::cerr << "Opaque nearly-real noncommuting Frechet derivatives are unstable: "
              << max_abs_diff(Matrix{near_real_noncommuting_dexp}, Matrix{near_real_noncommuting_dexp_reference})
              << ", "
              << max_abs_diff(Matrix{near_real_noncommuting_dphi}, Matrix{near_real_noncommuting_dphi_reference})
              << '\n';
    return false;
  }

  constexpr Numeric eigenvalue = 0.25;
  const Propmat     dichroism_generator{0.0, eigenvalue};
  const Numeric     ch = std::cosh(eigenvalue);
  const Numeric     sh = std::sinh(eigenvalue);
  const Muelmat     dichroism_reference{ch, sh, 0, 0, sh, ch, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
  if (not close(rtepack::tran{-dichroism_generator, -dichroism_generator, 1.0}(), dichroism_reference)) {
    std::cerr << "Pure-dichroism exponential uses the wrong eigenvalue\n";
    return false;
  }

  const Propmat rotation_generator{0.0, 0.0, 0.0, 0.0, eigenvalue};
  const Numeric co = std::cos(eigenvalue);
  const Numeric si = std::sin(eigenvalue);
  const Muelmat rotation_reference{1, 0, 0, 0, 0, co, si, 0, 0, -si, co, 0, 0, 0, 0, 1};
  if (not close(rtepack::tran{-rotation_generator, -rotation_generator, 1.0}(), rotation_reference)) {
    std::cerr << "Pure-rotation exponential uses the wrong eigenvalue\n";
    return false;
  }

  return true;
}

bool test_sqrt_stability() {
  // With N represented by {0, 1, 0, 0, 1, 0, 0}, N^3=0.  Consequently
  // sqrt(2 I + N) is a finite quadratic polynomial, despite all four
  // eigenvalues being exactly coalesced.  This catches implementations that
  // mistake a zero spectral separation for a scalar matrix.
  const Propmat k{2.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
  const Specmat q        = rtepack::sqrt(k);
  const Numeric residual = max_abs_diff(q * q, to_specmat(k));
  if (residual >= stability_tol) {
    std::cerr << "Propagation-matrix square root is unstable for a shifted nilpotent matrix: " << residual << '\n';
    return false;
  }

  // A tiny nonzero rotation must not be classified as the zero matrix.  Use
  // a relative residual because the entries of the squared root are itself
  // tiny; an absolute tolerance would not detect a root that was all zero.
  const Propmat tiny_rotation{0.0, 0.0, 0.0, 0.0, 1e-300, -2e-300, 3e-300};
  const Specmat tiny_rotation_matrix = to_specmat(tiny_rotation);
  const Specmat tiny_rotation_root   = rtepack::sqrt(tiny_rotation);
  const Numeric tiny_rotation_residual =
      max_abs_diff(tiny_rotation_root * tiny_rotation_root, tiny_rotation_matrix) / max_abs(tiny_rotation_matrix);
  if (not std::isfinite(tiny_rotation_residual) or tiny_rotation_residual >= 5e-14) {
    std::cerr << "Propagation-matrix square root loses a tiny rotational matrix: " << tiny_rotation_residual << '\n';
    return false;
  }

  // The spectral invariants contain products and squares of the propagation
  // coefficients.  They must remain finite far above the range where a
  // naive q*q calculation would overflow.
  const Propmat large{2e78, 3e77, -2e77, 1e77, 2.5e77, -1.5e77, 1.2e77};
  const Specmat large_matrix   = to_specmat(large);
  const Specmat large_root     = rtepack::sqrt(large);
  const Numeric large_residual = max_abs_diff(large_root * large_root, large_matrix) / max_abs(large_matrix);
  if (not std::isfinite(large_residual) or large_residual >= 5e-13) {
    std::cerr << "Propagation-matrix square root is unstable at large finite scale: " << large_residual << '\n';
    return false;
  }

  return true;
}

bool test_scalar_linprop_stability() {
  const auto source_operator = [](const Numeric start, const Numeric end, const Numeric r) {
    const Propmat       ks{start};
    const Propmat       ke{end};
    const rtepack::tran tr{ks, ke, r};
    return tr.linsrc_linprop(tr(), ks, ke, r)[0, 0];
  };

  const std::array cases{
      std::array{10.0, 0.0, 1.0, 0.39571230961051357},
      std::array{0.0, 10.0, 1.0, 0.11570508894007744},
      std::array{1.0, 1.0, 1.0, 0.6321205588285577},
      std::array{100.0, 99.9, 1.0, 0.010009909712423902},
      std::array{1000.0, 999.0, 1.0, 0.001000999997998006},
  };
  for (const auto& [start, end, r, reference] : cases) {
    const Numeric actual = source_operator(start, end, r);
    if (std::abs(actual - reference) > 2e-14 * std::max(1.0, std::abs(reference))) {
      std::cerr << "Scalar linprop source operator is unstable for endpoints " << start << ", " << end << ": " << actual
                << " instead of " << reference << '\n';
      return false;
    }
  }

  // At a zero gradient the value reduces to phi_1(-tau), but derivatives
  // with respect to the two endpoints are different because the within-layer
  // gradient changes in opposite directions.
  const Propmat       k{1.0};
  const Propmat       dk{1.0};
  const rtepack::tran tr{k, k, 1.0};
  const Muelmat       t      = tr();
  const Muelmat       l      = tr.linsrc_linprop(t, k, k, 1.0);
  const Numeric       dstart = tr.linsrc_linprop_deriv(l, t, k, k, dk, Muelmat{}, 1.0, 0.0, true)[0, 0];
  const Numeric       dend   = tr.linsrc_linprop_deriv(l, t, k, k, dk, Muelmat{}, 1.0, 0.0, false)[0, 0];
  if (std::abs(dstart - (-0.08030139707139419)) > 2e-14 or std::abs(dend - (-0.18393972058572117)) > 2e-14) {
    std::cerr << "Scalar linprop endpoint derivatives are discontinuous at zero gradient: " << dstart << ", " << dend
              << '\n';
    return false;
  }

  // r=0 is a regular limit and must not divide by the layer length.
  const Propmat       start{2.0};
  const Propmat       end{5.0};
  const rtepack::tran zero_tr{start, end, 0.0};
  const Muelmat       zero_t = zero_tr();
  const Muelmat       zero_l = zero_tr.linsrc_linprop(zero_t, start, end, 0.0);
  const Numeric       zero_dr =
      zero_tr.linsrc_linprop_deriv(zero_l, zero_t, start, end, Propmat{}, Muelmat{}, 0.0, 1.0, true)[0, 0];
  if (zero_l[0, 0] != 1.0 or std::abs(zero_dr - (-2.0)) > 2e-14) {
    std::cerr << "Scalar linprop has an incorrect zero-length limit\n";
    return false;
  }

  return true;
}

bool test_linprop_helper_consistency() {
  constexpr Numeric r = 0.8;

  const auto check_path =
      [](const Propmat& start, const Propmat& end, const Propmat& dk_start, const Propmat& dk_end, const char* label) {
        // TransmittanceMatrix paths are ordered observer-to-background, whereas
        // tran receives the background (start) endpoint first.
        const std::vector<Propmat>       k{end, start};
        const std::vector<PropmatVector> dk{PropmatVector(1, dk_end), PropmatVector(1, dk_start)};
        const Vector                     distance{r};
        Tensor3                          dr(2, 1, 1, 0.0);

        rtepack::TransmittanceMatrix tm;
        tm.init(k, dk, distance, dr, TransmittanceOption::linprop);

        const rtepack::tran direct{start, end, r};
        const Muelmat       t        = direct();
        const Muelmat       l        = direct.linsrc_linprop(t, start, end, r);
        const Muelmat       dt_end   = direct.deriv(t, start, end, dk_end, r, 0.0);
        const Muelmat       dt_start = direct.deriv(t, start, end, dk_start, r, 0.0);
        const Muelmat       dl_end   = direct.linsrc_linprop_deriv(l, t, start, end, dk_end, dt_end, r, 0.0, false);
        const Muelmat       dl_start = direct.linsrc_linprop_deriv(l, t, start, end, dk_start, dt_start, r, 0.0, true);

        if (not close(tm.T[0, 1], t, stability_tol) or not close(tm.L[0, 1], l, stability_tol) or
            not close(tm.dL[0, 0, 0, 0], dl_end, stability_tol) or
            not close(tm.dL[1, 0, 1, 0], dl_start, stability_tol)) {
          std::cerr << label << " linprop direct and TransmittanceMatrix paths disagree\n";
          return false;
        }
        return true;
      };

  const Propmat scalar_start{1.4};
  const Propmat scalar_end{0.3};
  const Propmat polarized_dstart{0.02, 0.03, -0.01, 0.015, -0.025, 0.02, 0.01};
  const Propmat polarized_dend{-0.01, -0.02, 0.025, 0.01, 0.015, -0.01, 0.02};
  if (not check_path(scalar_start, scalar_end, polarized_dstart, polarized_dend, "Scalar-base")) return false;

  const Propmat polarized_start{1.4, 0.08, -0.05, 0.03, 0.06, -0.04, 0.02};
  const Propmat polarized_end{0.3, -0.03, 0.07, 0.04, -0.02, 0.05, -0.06};
  if (not check_path(polarized_start, polarized_end, polarized_dstart, polarized_dend, "Polarized")) return false;

  // A polarized direction at an exactly scalar base exercises the
  // degenerate Frechet derivative.  Comparing against a centered difference
  // also makes sure the fast scalar value path does not discard the Q/U/V
  // derivative components.
  const rtepack::tran scalar_tr{scalar_start, scalar_end, r};
  const Muelmat       scalar_t   = scalar_tr();
  const Muelmat       scalar_l   = scalar_tr.linsrc_linprop(scalar_t, scalar_start, scalar_end, r);
  const Muelmat       scalar_dt  = scalar_tr.deriv(scalar_t, scalar_start, scalar_end, polarized_dstart, r, 0.0);
  const Muelmat       analytical = scalar_tr.linsrc_linprop_deriv(
      scalar_l, scalar_t, scalar_start, scalar_end, polarized_dstart, scalar_dt, r, 0.0, true);
  const Muelmat numerical = numerical_linprop_linsrc_deriv(scalar_start, scalar_end, polarized_dstart, r, 0.0, true);
  if (not close(analytical, numerical, 2e-8)) {
    std::cerr << "Polarized linprop derivative at a scalar base disagrees with centered differences\n";
    return false;
  }

  // The polarized construction must approach the exact scalar linprop
  // integral continuously, rather than falling back to phi_1 of the average
  // propagation matrix.
  const Propmat       tiny_start{10.0, 1e-14, -2e-14, 1.5e-14, -1e-14, 0.5e-14, 2.5e-14};
  const Propmat       tiny_end{0.0, -1.5e-14, 0.5e-14, -1e-14, 2e-14, -2.5e-14, 1e-14};
  const rtepack::tran tiny_tr{tiny_start, tiny_end, 1.0};
  const Numeric       tiny_l           = tiny_tr.linsrc_linprop(tiny_tr(), tiny_start, tiny_end, 1.0)[0, 0];
  constexpr Numeric   scalar_reference = 0.39571230961051357;
  if (std::abs(tiny_l - scalar_reference) >= 5e-14) {
    std::cerr << "Polarized linprop is discontinuous at the scalar limit: " << tiny_l << '\n';
    return false;
  }

  return true;
}

bool test_extreme_construction() {
  constexpr Numeric r = 1e-308;

  // Forming k1+k2 overflows here even though the scaled average is O(1).
  // Both constructors must average before scaling.
  const Propmat       scalar_k{1e308};
  const Numeric       tau          = r * scalar_k.A();
  const Numeric       transmission = std::exp(-tau);
  const rtepack::tran ordinary_scalar{scalar_k, scalar_k, r};
  const rtepack::tran magnus_scalar{scalar_k, scalar_k, r, rtepack::MagnusOperator::magnus};
  if (not close(ordinary_scalar(), Muelmat{transmission}, 5e-15) or
      not close(magnus_scalar(), Muelmat{transmission}, 5e-15)) {
    std::cerr << "Ordinary or Magnus construction overflows a finite scaled scalar average\n";
    return false;
  }

  const Muelmat scalar_t   = ordinary_scalar();
  const Muelmat scalar_l   = ordinary_scalar.linsrc_linprop(scalar_t, scalar_k, scalar_k, r);
  const Numeric expected_l = -std::expm1(-tau) / tau;
  if (std::abs(scalar_l[0, 0] - expected_l) >= 5e-15) {
    std::cerr << "Scalar linprop loses a finite value for huge unscaled endpoints\n";
    return false;
  }

  // Perturb the path length by its own magnitude.  The optical-depth
  // derivative is finite, although an intermediate ks+ke is not.
  const Numeric dr        = r;
  const Muelmat scalar_dt = ordinary_scalar.deriv(scalar_t, scalar_k, scalar_k, Propmat{}, r, dr);
  const Numeric actual_dl = ordinary_scalar.linsrc_linprop_deriv(
      scalar_l, scalar_t, scalar_k, scalar_k, Propmat{}, scalar_dt, r, dr, true)[0, 0];
  const Numeric expected_dl = ((tau + 1.0) * std::exp(-tau) - 1.0) / tau;
  if (not std::isfinite(actual_dl) or std::abs(actual_dl - expected_dl) >= 2e-14) {
    std::cerr << "Scalar linprop derivative overflows a finite scaled endpoint average: " << actual_dl << '\n';
    return false;
  }

  // A noncommuting Magnus pair also has to form its commutator from the
  // scaled endpoint matrices.  Computing either unscaled matrix product
  // would overflow by hundreds of orders of magnitude.
  const Propmat k0{1e308, 2e307, -1e307, 1.5e307, 0.7e307, -0.5e307, 0.3e307};
  const Propmat k1{9e307, -1.5e307, 2e307, -0.5e307, -0.2e307, 0.8e307, -0.4e307};
  const Matrix  unscaled_m0 = to_matrix(k0);
  const Matrix  unscaled_m1 = to_matrix(k1);
  const Matrix  m0          = to_matrix(r * k0);
  const Matrix  m1          = to_matrix(r * k1);
  Matrix        m1m0(4, 4);
  Matrix        m0m1(4, 4);
  Matrix        omega(4, 4);
  mult(m1m0, m1, m0);
  mult(m0m1, m0, m1);
  for (Index i = 0; i < 4; ++i) {
    for (Index j = 0; j < 4; ++j) {
      omega[i, j] = -r * std::midpoint(unscaled_m0[i, j], unscaled_m1[i, j]) + (m1m0[i, j] - m0m1[i, j]) / 12.0;
    }
  }
  const Propmat expected_generator{
      omega[0, 0], omega[0, 1], omega[0, 2], omega[0, 3], omega[1, 2], omega[1, 3], omega[2, 3]};
  const MatrixFunctionReference expected = matrix_function_reference(expected_generator);
  const rtepack::tran           magnus{k0, k1, r, rtepack::MagnusOperator::magnus};
  const Matrix                  actual_generator =
      to_matrix(Propmat{magnus.a, magnus.b, magnus.c, magnus.d, magnus.u, magnus.v, magnus.w});
  if (max_abs_diff(actual_generator, omega) >= 5e-14 or max_abs_diff(Matrix{magnus()}, expected.exp) >= 5e-13) {
    std::cerr << "Magnus construction overflows a finite scaled commutator\n";
    return false;
  }

  return true;
}

bool test_transmittance_matrix() {
  constexpr Numeric r = 0.7;
  const Propmat     k0{0.8, 0.03, -0.02, 0.01, 0.04, -0.015, 0.025};
  const Propmat     k1{1.1, -0.01, 0.05, 0.02, -0.03, 0.035, -0.02};
  const Propmat     dk0{0.1, -0.02, 0.01, 0.03, 0.02, -0.01, 0.015};
  const Propmat     dk1{-0.04, 0.03, 0.02, -0.01, 0.015, 0.025, -0.02};

  const std::vector<PropmatVector> K{PropmatVector(1, k0), PropmatVector(1, k1)};
  const std::vector<PropmatMatrix> dK{PropmatMatrix(1, 1, dk0), PropmatMatrix(1, 1, dk1)};
  const Vector                     distance{r};
  Tensor3                          dr(2, 1, 1, 0.0);
  dr[0, 0, 0] = 0.04;
  dr[1, 0, 0] = -0.03;

  rtepack::TransmittanceMatrix tm;
  tm.init(K, dK, distance, dr, TransmittanceOption::magop);
  tm.check(2, 1, 1, "test_rtepack_magnus");

  // The path is stored observer-to-background, while radiation propagates
  // from k1 to k0.
  const rtepack::tran tr{k1, k0, r, rtepack::MagnusOperator::magnus};
  if (not close(tm.T[0, 1], tr())) {
    std::cerr << "TransmittanceMatrix did not store the Magnus transmission\n";
    return false;
  }

  if (not close(tm.dT[0, 0, 0, 0], tr.magnus_deriv(k1, k0, dk0, r, dr[0, 0, 0], false)) or
      not close(tm.dT[1, 0, 1, 0], tr.magnus_deriv(k1, k0, dk1, r, dr[1, 0, 0], true))) {
    std::cerr << "TransmittanceMatrix did not store the Magnus endpoint derivatives\n";
    return false;
  }

  const std::vector<Propmat>       scalar_K{k0, k1};
  const std::vector<PropmatVector> scalar_dK{PropmatVector(1, dk0), PropmatVector(1, dk1)};
  rtepack::TransmittanceMatrix     scalar_tm;
  scalar_tm.init(scalar_K, scalar_dK, distance, dr, TransmittanceOption::magop);
  scalar_tm.check(2, 1, 1, "test_rtepack_magnus_single_frequency");
  if (not close(scalar_tm.T[0, 1], tr()) or
      not close(scalar_tm.dT[0, 0, 0, 0], tr.magnus_deriv(k1, k0, dk0, r, dr[0, 0, 0], false)) or
      not close(scalar_tm.dT[1, 0, 1, 0], tr.magnus_deriv(k1, k0, dk1, r, dr[1, 0, 0], true))) {
    std::cerr << "The single-frequency TransmittanceMatrix Magnus path is inconsistent\n";
    return false;
  }

  rtepack::TransmittanceMatrix linsrc_tm;
  linsrc_tm.init(K, dK, distance, dr, TransmittanceOption::magop_linsrc);
  linsrc_tm.check(2, 1, 1, "test_rtepack_magnus_linsrc");

  const Muelmat expected_l = tr.linsrc() * (Muelmat::id() - (r / 12.0) * rtepack::to_muelmat(k0 - k1));
  if (not close(linsrc_tm.T[0, 1], tr()) or not close(linsrc_tm.L[0, 1], expected_l) or
      not close(linsrc_tm.dT[0, 0, 0, 0], tr.magnus_deriv(k1, k0, dk0, r, dr[0, 0, 0], false)) or
      not close(linsrc_tm.dT[1, 0, 1, 0], tr.magnus_deriv(k1, k0, dk1, r, dr[1, 0, 0], true)) or
      not close(linsrc_tm.dL[0, 0, 0, 0], tr.magnus_linsrc_deriv(k1, k0, dk0, r, dr[0, 0, 0], false)) or
      not close(linsrc_tm.dL[1, 0, 1, 0], tr.magnus_linsrc_deriv(k1, k0, dk1, r, dr[1, 0, 0], true))) {
    std::cerr << "The Magnus linear-source TransmittanceMatrix path is inconsistent\n";
    return false;
  }

  rtepack::TransmittanceMatrix scalar_linsrc_tm;
  scalar_linsrc_tm.init(scalar_K, scalar_dK, distance, dr, TransmittanceOption::magop_linsrc);
  scalar_linsrc_tm.check(2, 1, 1, "test_rtepack_magnus_linsrc_single_frequency");
  if (not close(scalar_linsrc_tm.T[0, 1], linsrc_tm.T[0, 1]) or
      not close(scalar_linsrc_tm.L[0, 1], linsrc_tm.L[0, 1]) or
      not close(scalar_linsrc_tm.dT[0, 0, 0, 0], linsrc_tm.dT[0, 0, 0, 0]) or
      not close(scalar_linsrc_tm.dL[1, 0, 1, 0], linsrc_tm.dL[1, 0, 1, 0])) {
    std::cerr << "The single-frequency Magnus linear-source path is inconsistent\n";
    return false;
  }

  rtepack::TransmittanceMatrix linprop_tm;
  linprop_tm.init(K, dK, distance, dr, TransmittanceOption::linprop);
  const rtepack::tran ordinary_tr{k1, k0, r};
  const Muelmat       ordinary_t = ordinary_tr();
  const Muelmat       ordinary_l = ordinary_tr.linsrc_linprop(ordinary_t, k1, k0, r);
  if (not close(linprop_tm.T[0, 1], ordinary_t) or not close(linprop_tm.L[0, 1], ordinary_l) or
      not close(linprop_tm.dL[0, 0, 0, 0],
                ordinary_tr.linsrc_linprop_deriv(
                    ordinary_l, ordinary_t, k1, k0, dk0, linprop_tm.dT[0, 0, 0, 0], r, dr[0, 0, 0], false)) or
      not close(linprop_tm.dL[1, 0, 1, 0],
                ordinary_tr.linsrc_linprop_deriv(
                    ordinary_l, ordinary_t, k1, k0, dk1, linprop_tm.dT[1, 0, 1, 0], r, dr[1, 0, 0], true))) {
    std::cerr << "The linprop path reverses endpoints or does not use its analytic polarized source operator\n";
    return false;
  }

  StokvecMatrix radiance(1, 2, Stokvec{});
  SourceVector  source;
  source.J.resize(1, 2);
  radiance[0, 1]                  = Stokvec{2.0, 0.1, -0.05, 0.02};
  source.J[0, 0]                  = Stokvec{0.8, 0.02, 0.01, -0.01};
  source.J[0, 1]                  = Stokvec{1.1, -0.01, 0.03, 0.015};
  const Stokvec expected_radiance = linsrc_tm.T[0, 1] * (radiance[0, 1] - source.J[0, 1]) +
                                    linsrc_tm.L[0, 1] * (source.J[0, 1] - source.J[0, 0]) + source.J[0, 0];
  rtepack::rte_emission_path(radiance, linsrc_tm, source);
  const Numeric stokes_error = std::max({std::abs(radiance[0, 0].I() - expected_radiance.I()),
                                         std::abs(radiance[0, 0].Q() - expected_radiance.Q()),
                                         std::abs(radiance[0, 0].U() - expected_radiance.U()),
                                         std::abs(radiance[0, 0].V() - expected_radiance.V())});
  if (stokes_error >= tol) {
    std::cerr << "Magnus linear-source emission routing is inconsistent\n";
    return false;
  }

  return true;
}

bool test_analytical_derivatives() {
  constexpr Numeric   r   = 0.65;
  constexpr Numeric   dr0 = 0.04;
  constexpr Numeric   dr1 = -0.03;
  const Propmat       k0{0.9, 0.18, -0.13, 0.09, 0.14, -0.11, 0.08};
  const Propmat       k1{1.2, -0.12, 0.16, 0.11, -0.09, 0.13, -0.15};
  const Propmat       dk0{0.07, -0.04, 0.03, 0.02, -0.05, 0.01, 0.025};
  const Propmat       dk1{-0.03, 0.05, 0.015, -0.025, 0.02, 0.04, -0.035};
  const rtepack::tran tr{k0, k1, r, rtepack::MagnusOperator::magnus};

  const Muelmat dt0 = tr.magnus_deriv(k0, k1, dk0, r, dr0, true);
  const Muelmat dt1 = tr.magnus_deriv(k0, k1, dk1, r, dr1, false);
  const Muelmat nt0 = numerical_magnus_deriv(k0, k1, dk0, r, dr0, true);
  const Muelmat nt1 = numerical_magnus_deriv(k0, k1, dk1, r, dr1, false);
  if (not close(dt0, nt0, 2e-8) or not close(dt1, nt1, 2e-8)) {
    std::cerr << "Analytical Magnus transmission derivatives disagree with centered differences\n";
    return false;
  }

  const Muelmat dl0 = tr.magnus_linsrc_deriv(k0, k1, dk0, r, dr0, true);
  const Muelmat dl1 = tr.magnus_linsrc_deriv(k0, k1, dk1, r, dr1, false);
  const Muelmat nl0 = numerical_magnus_linsrc_deriv(k0, k1, dk0, r, dr0, true);
  const Muelmat nl1 = numerical_magnus_linsrc_deriv(k0, k1, dk1, r, dr1, false);
  if (not close(dl0, nl0, 2e-8) or not close(dl1, nl1, 2e-8)) {
    std::cerr << "Analytical Magnus linear-source derivatives disagree with centered differences\n";
    return false;
  }

  const Propmat       linprop_k0{0.4};
  const Propmat       linprop_k1{1.1};
  const Propmat       linprop_dk0{0.07};
  const Propmat       linprop_dk1{-0.03};
  const rtepack::tran linprop_tr{linprop_k0, linprop_k1, r};
  const Muelmat       linprop_t   = linprop_tr();
  const Muelmat       linprop_l   = linprop_tr.linsrc_linprop(linprop_t, linprop_k0, linprop_k1, r);
  const Muelmat       linprop_dt0 = linprop_tr.deriv(linprop_t, linprop_k0, linprop_k1, linprop_dk0, r, dr0);
  const Muelmat       linprop_dt1 = linprop_tr.deriv(linprop_t, linprop_k0, linprop_k1, linprop_dk1, r, dr1);
  if (not close(linprop_tr.linsrc_linprop_deriv(
                    linprop_l, linprop_t, linprop_k0, linprop_k1, linprop_dk0, linprop_dt0, r, dr0, true),
                numerical_linprop_linsrc_deriv(linprop_k0, linprop_k1, linprop_dk0, r, dr0, true),
                2e-8) or
      not close(linprop_tr.linsrc_linprop_deriv(
                    linprop_l, linprop_t, linprop_k0, linprop_k1, linprop_dk1, linprop_dt1, r, dr1, false),
                numerical_linprop_linsrc_deriv(linprop_k0, linprop_k1, linprop_dk1, r, dr1, false),
                2e-8)) {
    std::cerr << "Analytical linear-propagation source derivatives disagree with centered differences\n";
    return false;
  }

  // Exercise the degenerate fallback: scalar base propagation with a
  // polarized derivative direction.
  const Propmat       scalar0{0.8};
  const Propmat       scalar1{1.1};
  const Propmat       polarized_dk{0.02, 0.03, -0.01, 0.015, -0.025, 0.02, 0.01};
  const rtepack::tran scalar_tr{scalar0, scalar1, r, rtepack::MagnusOperator::magnus};
  if (not close(scalar_tr.magnus_deriv(scalar0, scalar1, polarized_dk, r, dr0, true),
                numerical_magnus_deriv(scalar0, scalar1, polarized_dk, r, dr0, true),
                2e-8) or
      not close(scalar_tr.magnus_linsrc_deriv(scalar0, scalar1, polarized_dk, r, dr0, true),
                numerical_magnus_linsrc_deriv(scalar0, scalar1, polarized_dk, r, dr0, true),
                2e-8)) {
    std::cerr << "Degenerate Magnus analytical derivatives disagree with centered differences\n";
    return false;
  }

  const auto check_degenerate_case = [&](const Propmat& base0, const Propmat& base1, const char* label) {
    const rtepack::tran degenerate_tr{base0, base1, r, rtepack::MagnusOperator::magnus};
    if (not close(degenerate_tr.magnus_deriv(base0, base1, polarized_dk, r, dr0, true),
                  numerical_magnus_deriv(base0, base1, polarized_dk, r, dr0, true),
                  2e-8) or
        not close(degenerate_tr.magnus_linsrc_deriv(base0, base1, polarized_dk, r, dr0, true),
                  numerical_magnus_linsrc_deriv(base0, base1, polarized_dk, r, dr0, true),
                  2e-8)) {
      std::cerr << label << " Magnus analytical derivatives disagree with centered differences\n";
      return false;
    }
    return true;
  };

  if (not check_degenerate_case(Propmat{0.8, 0.2}, Propmat{1.1, -0.15}, "Pure-dichroism") or
      not check_degenerate_case(Propmat{0.8, 0.0, 0.0, 0.0, 0.2, -0.1, 0.15},
                                Propmat{1.1, 0.0, 0.0, 0.0, -0.12, 0.18, -0.08},
                                "Pure-rotation")) {
    return false;
  }

  return true;
}
}  // namespace

int main() {
  return test_exponential() and test_sqrt_stability() and test_scalar_linprop_stability() and
                 test_linprop_helper_consistency() and test_extreme_construction() and test_analytical_derivatives() and
                 test_transmittance_matrix()
             ? 0
             : 1;
}
