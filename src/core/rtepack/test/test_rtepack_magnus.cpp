#include <lin_alg.h>
#include <rtepack.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
constexpr Numeric tol = 2e-10;

Numeric max_abs_diff(const Matrix& a, const Matrix& b) {
  Numeric out{};
  for (Index i = 0; i < a.nrows(); ++i) {
    for (Index j = 0; j < a.ncols(); ++j) out = std::max(out, std::abs(a[i, j] - b[i, j]));
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
  const Propmat clustered{0.3, 1.1e-4, -0.9e-4, 0.7e-4, -1.0e-4, 0.8e-4, -0.6e-4};
  const Matrix  clustered_g = to_matrix(-clustered);
  Matrix        phi1_reference(4, 4, 0.0);
  Matrix        phi1_term(4, 4, 0.0);
  Matrix        phi1_product(4, 4);
  for (Index i = 0; i < 4; ++i) phi1_reference[i, i] = phi1_term[i, i] = 1.0;
  for (Index n = 1; n <= 30; ++n) {
    mult(phi1_product, phi1_term, clustered_g);
    phi1_term       = phi1_product;
    phi1_term      /= static_cast<Numeric>(n + 1);
    phi1_reference += phi1_term;
  }
  const Matrix phi1_actual = Matrix{rtepack::tran{clustered, clustered, 1.0}.linsrc()};
  if (max_abs_diff(phi1_actual, phi1_reference) >= tol) {
    std::cerr << "Linear-source operator is unstable near clustered eigenvalues\n";
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

int main() { return test_exponential() and test_analytical_derivatives() and test_transmittance_matrix() ? 0 : 1; }
