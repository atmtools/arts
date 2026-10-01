#include <lin_alg.h>

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>

namespace {
void require(const bool condition, const std::string_view message) {
  if (not condition) throw std::runtime_error(std::string(message));
}

void close(const Complex actual, const Complex expected, const std::string_view message) {
  require(std::abs(actual - expected) <= 2e-12 * std::max(1.0, std::abs(expected)), message);
}

void rejects(const std::function<void()>& operation, const std::string_view message) {
  try {
    operation();
  } catch (const std::exception& error) {
    require(not std::string_view(error.what()).empty(), "Validation error has no explanation");
    return;
  }
  throw std::runtime_error(std::string(message));
}

/** A nonfinite input that is not searched for must still be visible.
 *
 * LAPACK reports what it detects itself; anything it does not report has to
 * reach the output, where the caller detects it once for a whole result.
 */
void rejects_or_propagates(const std::function<bool()>& operation, const std::string_view message) {
  try {
    require(operation(), message);
  } catch (const std::exception& error) {
    require(not std::string_view(error.what()).empty(), "Validation error has no explanation");
  }
}

bool nonfinite(const Complex value) { return not(std::isfinite(value.real()) and std::isfinite(value.imag())); }

template <class T> bool any_nonfinite(const T& values) {
  for (const Complex value : matpack::elemwise_range(values))
    if (nonfinite(value)) return true;
  return false;
}

ComplexMatrix nonsymmetric() {
  ComplexMatrix A(3, 3);
  A[0, 0] = {3, 1};
  A[0, 1] = {2, -1};
  A[0, 2] = {-1, 2};
  A[1, 0] = {-1, 0};
  A[1, 1] = {4, -2};
  A[1, 2] = {0, 1};
  A[2, 0] = {2, 1};
  A[2, 1] = {-2, -1};
  A[2, 2] = {5, 3};
  return A;
}

ComplexVector product(StridedConstComplexMatrixView A, StridedConstComplexVectorView x) {
  ComplexVector out(A.nrows(), 0);
  for (Index i = 0; i < A.nrows(); ++i)
    for (Index j = 0; j < A.ncols(); ++j) out[i] += A[i, j] * x[j];
  return out;
}

void same(StridedConstComplexMatrixView A, StridedConstComplexMatrixView B) {
  require(A.shape() == B.shape(), "Changed matrix dimensions");
  for (Index i = 0; i < A.nrows(); ++i)
    for (Index j = 0; j < A.ncols(); ++j) require(A[i, j] == B[i, j], "Changed input matrix");
}

void eigen_residual(StridedConstComplexMatrixView A, StridedConstComplexMatrixView P, StridedConstComplexVectorView W) {
  for (Index j = 0; j < A.ncols(); ++j) {
    Numeric norm = 0;
    for (Index i = 0; i < A.nrows(); ++i) {
      Complex ap{};
      for (Index k = 0; k < A.ncols(); ++k) ap += A[i, k] * P[k, j];
      close(ap, P[i, j] * W[j], "Right-eigenvector residual is too large");
      norm += std::norm(P[i, j]);
    }
    close(norm, 1, "Eigenvector is not normalized");
  }
}

void ordinary() {
  const auto          A = nonsymmetric();
  const ComplexVector expected{{1, 2}, {-2, 1}, {3, -1}};
  const auto          b        = product(A, expected);
  const auto          original = A;
  ComplexVector       x(3);
  const Numeric       rcond = solve(x, A, b, 1e-12);
  require(rcond > 0 and rcond <= 1, "Invalid reciprocal condition estimate");
  for (Size i = 0; i < x.size(); ++i) close(x[i], expected[i], "Incorrect complex solution");
  same(A, original);

  auto rhs = b;
  solve(rhs, A, rhs);
  for (Size i = 0; i < rhs.size(); ++i) close(rhs[i], expected[i], "Aliased solve is incorrect");

  ComplexMatrix                P(3, 3);
  ComplexVector                W(3);
  complex_diagonalize_workdata work(3);
  diagonalize(P, W, A, work);
  eigen_residual(A, P, W);
  same(A, original);
  require(work.work.size() >= 6, "Eigenvalue workspace query was not retained");
  diagonalize(P, W, A, work);
  eigen_residual(A, P, W);

  auto alias = A;
  diagonalize(alias, W, alias);
  eigen_residual(A, alias, W);
}

void strided() {
  constexpr Complex   guard{17, -31};
  ComplexMatrix       backing(6, 7, guard), vectors(6, 7, guard);
  ComplexVector       right(7, guard), result(7, guard), values(7, guard);
  auto                A = backing[StridedRange(0, 3, 2), StridedRange(1, 3, 2)];
  auto                P = vectors[StridedRange(0, 3, 2), StridedRange(1, 3, 2)];
  auto                b = right[StridedRange(1, 3, 2)];
  auto                x = result[StridedRange(1, 3, 2)];
  auto                W = values[StridedRange(1, 3, 2)];
  const ComplexVector expected{{1, -2}, {2, 1}, {-1, 3}};
  A                   = nonsymmetric();
  b                   = product(A, expected);
  const auto original = backing;
  solve(x, A, b);
  for (Size i = 0; i < x.size(); ++i) close(x[i], expected[i], "Incorrect strided solution");
  diagonalize(P, W, A);
  eigen_residual(A, P, W);
  same(backing, original);
  for (Index i = 0; i < 7; i += 2) {
    require(result[i] == guard, "Strided solve overwrote an unselected value");
    require(values[i] == guard, "Strided eigenvalues overwrote an unselected value");
  }
  for (Index i = 0; i < 6; ++i)
    for (Index j = 0; j < 7; ++j)
      if (i % 2 != 0 or j % 2 == 0) require(vectors[i, j] == guard, "Strided eigenvectors overwrote padding");
}

void multiple_rhs_solve() {
  const auto    A = nonsymmetric();
  ComplexMatrix expected(3, 4), rhs(3, 4, 0), solution(3, 4);
  for (Index i = 0; i < 3; ++i)
    for (Index j = 0; j < 4; ++j) expected[i, j] = Complex{Numeric(1 + i + 2 * j), Numeric(2 * i - j)};
  for (Index i = 0; i < 3; ++i)
    for (Index j = 0; j < 4; ++j)
      for (Index k = 0; k < 3; ++k) rhs[i, j] += A[i, k] * expected[k, j];
  const auto    original_rhs = rhs;
  const Numeric rcond        = solve(solution, A, rhs, 1e-12);
  require(rcond > 0 and rcond <= 1, "Invalid multiple-RHS condition estimate");
  same(rhs, original_rhs);
  for (Index i = 0; i < 3; ++i)
    for (Index j = 0; j < 4; ++j) close(solution[i, j], expected[i, j], "Incorrect multiple-RHS solution");
  auto alias = rhs;
  close(solve(alias, A, alias, 1e-12), rcond, "Condition estimate changed with aliased matrix RHS");
  for (Index i = 0; i < 3; ++i)
    for (Index j = 0; j < 4; ++j) close(alias[i, j], expected[i, j], "Incorrect aliased matrix-RHS solve");

  constexpr Complex guard{17, -31};
  ComplexMatrix     rhs_storage(6, 9, guard), solution_storage(6, 9, guard);
  auto              strided_rhs      = rhs_storage[StridedRange(0, 3, 2), StridedRange(1, 4, 2)];
  auto              strided_solution = solution_storage[StridedRange(0, 3, 2), StridedRange(1, 4, 2)];
  strided_rhs                        = rhs;
  close(solve(strided_solution, A, strided_rhs), rcond, "Strided matrix-RHS condition estimate changed");
  for (Index i = 0; i < 3; ++i)
    for (Index j = 0; j < 4; ++j) close(strided_solution[i, j], expected[i, j], "Incorrect strided matrix-RHS solve");
  for (Index i = 0; i < 6; ++i)
    for (Index j = 0; j < 9; ++j)
      if (i % 2 != 0 or j % 2 == 0)
        require(solution_storage[i, j] == guard, "Matrix-RHS solve overwrote strided padding");

  solution  = guard;
  rhs[2, 3] = {0, std::numeric_limits<Numeric>::infinity()};
  rejects_or_propagates(
      [&] {
        solve(solution, A, rhs);
        return any_nonfinite(solution);
      },
      "Nonfinite later RHS column neither rejected nor propagated");
  solution = guard;
  ComplexMatrix wrong_output(3, 2);
  rejects([&] { solve(wrong_output, A, original_rhs); }, "Accepted mismatched multiple-RHS output dimensions");
  ComplexMatrix no_rhs(3, 0), no_solution(3, 0);
  close(solve(no_solution, A, no_rhs), rcond, "Zero-RHS condition estimate is incorrect");
  ComplexMatrix empty_A(0, 0), empty_rhs(0, 4), empty_solution(0, 4);
  require(solve(empty_solution, empty_A, empty_rhs) == 1, "Empty matrix-RHS system has incorrect condition estimate");
}

void invalid_and_empty() {
  ComplexMatrix singular(2, 2, 1);
  ComplexVector b(2, 1), x(2, Complex{17, -31});
  rejects([&] { solve(x, singular, b); }, "Accepted a singular complex matrix");
  for (Complex value : x) require(value == Complex{17, -31}, "Failed solve changed its output");

  ComplexMatrix ill(2, 2, 0);
  ill[0, 0] = 1;
  ill[1, 1] = 1e-14;
  rejects([&] { solve(x, ill, b, 1e-12); }, "Accepted an ill-conditioned complex matrix");
  const Numeric rcond = solve(x, ill, b);
  require(std::abs(rcond / 1e-14 - 1) < 1e-12, "Incorrect condition estimate for diagonal matrix");
  close(x[1], 1e14, "Incorrect unthresholded solution");

  for (const Numeric threshold : {-1.0, 1.1, std::numeric_limits<Numeric>::quiet_NaN()})
    rejects([&] { solve(x, ill, b, threshold); }, "Accepted an invalid condition threshold");

  ComplexMatrix P(2, 2, 7), rectangular(2, 3);
  ComplexVector W(2, 7), short_vector(1);
  rejects([&] { solve(x, rectangular, b); }, "Accepted a nonsquare solve matrix");
  rejects([&] { solve(short_vector, ill, b); }, "Accepted incorrect solution dimensions");
  rejects([&] { diagonalize(P, W, rectangular); }, "Accepted a nonsquare eigenproblem");
  rejects([&] { diagonalize(P, short_vector, ill); }, "Accepted incorrect eigenvalue dimensions");
  complex_diagonalize_workdata wrong_work(1);
  rejects([&] { diagonalize(P, W, ill, wrong_work); }, "Accepted incorrect workspace dimensions");
  rejects([] { complex_diagonalize_workdata invalid(-1); }, "Accepted a negative workspace size");

  // The 1-norm needed for the condition estimate rejects a nonfinite matrix
  // for free.  A nonfinite right-hand side has no such aggregate and is left
  // to propagate, as is a nonfinite eigenproblem that LAPACK accepts.
  ill[0, 0] = {1, std::numeric_limits<Numeric>::infinity()};
  rejects([&] { solve(x, ill, b); }, "Accepted nonfinite solve matrix input");
  rejects_or_propagates(
      [&] {
        diagonalize(P, W, ill);
        return any_nonfinite(W) or any_nonfinite(P);
      },
      "Nonfinite eigenproblem neither rejected nor propagated");
  P         = Complex{7, 0};
  W         = Complex{7, 0};
  ill[0, 0] = 1;
  b[0]      = {std::numeric_limits<Numeric>::quiet_NaN(), 0};
  rejects_or_propagates(
      [&] {
        solve(x, ill, b);
        return any_nonfinite(x);
      },
      "Nonfinite right-hand side neither rejected nor propagated");
  b[0] = 1;

  // A valid eigenproblem can have a singular eigenvector basis. The caller's
  // condition threshold must reject it before computing equivalent strengths.
  ComplexMatrix jordan(2, 2, 0);
  jordan[0, 0] = jordan[1, 1] = jordan[0, 1] = 1;
  diagonalize(P, W, jordan);
  b = 1;
  rejects([&] { solve(x, P, b, 1e-12); }, "Accepted a defective eigenvector basis");

  ComplexMatrix empty(0, 0);
  ComplexVector empty_vector(0);
  require(solve(empty_vector, empty, empty_vector) == 1, "Empty system has incorrect condition estimate");
  diagonalize(empty, empty_vector, empty);
  complex_diagonalize_workdata empty_work;
  diagonalize(empty, empty_vector, empty, empty_work);
}

void analytic_derivatives() {
  // A(t)=S(t) D(t) S(t)^-1 with S'(0)=K S(0). The exact mode
  // derivatives are D' and K*p (after unit-norm/phase projection).
  // A is deliberately nonnormal and complex, so P^H cannot replace P^-1.
  ComplexMatrix S(3, 3, 0), inverse_S(3, 3, 0), K(3, 3);
  for (Index i = 0; i < 3; ++i) S[i, i] = inverse_S[i, i] = 1;
  S[0, 1]         = {0.8, -0.4};
  S[0, 2]         = {-0.3, 0.6};
  S[1, 2]         = {0.5, 0.2};
  inverse_S[0, 1] = -S[0, 1];
  inverse_S[1, 2] = -S[1, 2];
  inverse_S[0, 2] = S[0, 1] * S[1, 2] - S[0, 2];
  K[0, 0]         = {0.1, 0.2};
  K[0, 1]         = {0.3, -0.4};
  K[0, 2]         = {-0.2, 0.1};
  K[1, 0]         = {-0.5, 0.2};
  K[1, 1]         = {0.2, -0.1};
  K[1, 2]         = {0.4, 0.3};
  K[2, 0]         = {0.1, 0.4};
  K[2, 1]         = {-0.3, 0.2};
  K[2, 2]         = {-0.2, 0.3};
  const ComplexVector expected_values{{1, 0.4}, {3, 1.2}, {6, -0.5}};
  const ComplexVector expected_derivatives{{0.2, -0.1}, {-0.3, 0.4}, {0.5, 0.2}};
  ComplexMatrix       A(3, 3, 0), dA(3, 3, 0);
  for (Index i = 0; i < 3; ++i)
    for (Index j = 0; j < 3; ++j)
      for (Index k = 0; k < 3; ++k) {
        A[i, j]  += S[i, k] * expected_values[k] * inverse_S[k, j];
        dA[i, j] += S[i, k] * expected_derivatives[k] * inverse_S[k, j];
      }
  for (Index i = 0; i < 3; ++i)
    for (Index j = 0; j < 3; ++j)
      for (Index k = 0; k < 3; ++k) dA[i, j] += K[i, k] * A[k, j] - A[i, k] * K[k, j];
  const auto                   original_A = A, original_dA = dA;
  ComplexMatrix                P(3, 3), dP(3, 3);
  ComplexVector                W(3), dW(3);
  complex_diagonalize_workdata work(3);
  diagonalize(P, W, dP, dW, A, dA, work);
  same(A, original_A);
  same(dA, original_dA);
  eigen_residual(A, P, W);

  const auto matching = [](const Complex value, StridedConstComplexVectorView candidates) {
    Index best = 0;
    for (Index i = 1; i < static_cast<Index>(candidates.size()); ++i)
      if (std::abs(candidates[i] - value) < std::abs(candidates[best] - value)) best = i;
    return best;
  };
  for (Index mode = 0; mode < 3; ++mode) {
    const Index exact = matching(W[mode], expected_values);
    close(W[mode], expected_values[exact], "Wrong eigenvalue in analytic derivative fixture");
    close(dW[mode], expected_derivatives[exact], "Wrong analytic eigenvalue derivative");
    ComplexVector expected(3, 0);
    Complex       projection = 0, gauge = 0;
    for (Index i = 0; i < 3; ++i) {
      for (Index j = 0; j < 3; ++j) expected[i] += K[i, j] * P[j, mode];
      projection += std::conj(P[i, mode]) * expected[i];
      gauge      += std::conj(P[i, mode]) * dP[i, mode];
    }
    close(gauge, 0, "Eigenvector derivative violates its normalization/phase gauge");
    for (Index i = 0; i < 3; ++i) {
      expected[i] -= P[i, mode] * projection;
      close(dP[i, mode], expected[i], "Wrong analytic eigenvector derivative");
      Complex residual = -dW[mode] * P[i, mode] - W[mode] * dP[i, mode];
      for (Index j = 0; j < 3; ++j) residual += dA[i, j] * P[j, mode] + A[i, j] * dP[j, mode];
      close(residual, 0, "Differentiated eigenproblem has a nonzero residual");
    }
  }

  // Several Jacobian targets share the eigendecomposition: a general complex
  // perturbation, a carrier shift, and a pure similarity transformation.
  ComplexTensor3 directions(3, 3, 3, 0), batch_dP(3, 3, 3);
  ComplexMatrix  batch_dW(3, 3), batch_P(3, 3);
  ComplexVector  batch_W(3);
  directions[0] = dA;
  constexpr Complex carrier_derivative{2, -0.7};
  for (Index i = 0; i < 3; ++i) directions[1, i, i] = carrier_derivative;
  // K2 has its only nonzero entry K2[0,1]=1.
  for (Index i = 0; i < 3; ++i)
    for (Index j = 0; j < 3; ++j) directions[2, i, j] = (i == 0 ? A[1, j] : Complex{}) - (j == 1 ? A[i, 0] : Complex{});
  diagonalize(batch_P, batch_W, batch_dP, batch_dW, A, directions, work);
  for (Index mode = 0; mode < 3; ++mode) {
    close(batch_W[mode], W[mode], "Batched primal eigenvalues differ");
    close(batch_dW[0, mode], dW[mode], "Batched general eigenvalue derivative differs");
    close(batch_dW[1, mode], carrier_derivative, "Batched carrier derivative is incorrect");
    close(batch_dW[2, mode], 0, "A similarity transformation changed an eigenvalue");
    const Complex projection = std::conj(P[0, mode]) * P[1, mode];
    for (Index row = 0; row < 3; ++row) {
      close(batch_P[row, mode], P[row, mode], "Batched primal eigenvectors differ");
      close(batch_dP[0, row, mode], dP[row, mode], "Batched general eigenvector derivative differs");
      close(batch_dP[1, row, mode], 0, "A carrier shift changed an eigenvector");
      const Complex expected = (row == 0 ? P[1, mode] : Complex{}) - P[row, mode] * projection;
      close(batch_dP[2, row, mode], expected, "Batched similarity eigenvector derivative is incorrect");
    }
  }

  // Verify target-page strides as well as matrix/vector strides.
  constexpr Complex batch_guard{19, -23};
  ComplexTensor3    direction_storage(6, 6, 7, batch_guard), derivative_storage(6, 6, 7, batch_guard);
  ComplexMatrix     eigenvalue_storage(6, 7, batch_guard);
  auto strided_directions  = direction_storage[StridedRange(0, 3, 2), StridedRange(0, 3, 2), StridedRange(1, 3, 2)];
  auto strided_derivatives = derivative_storage[StridedRange(0, 3, 2), StridedRange(0, 3, 2), StridedRange(1, 3, 2)];
  auto strided_value_derivatives = eigenvalue_storage[StridedRange(0, 3, 2), StridedRange(1, 3, 2)];
  strided_directions             = directions;
  diagonalize(batch_P, batch_W, strided_derivatives, strided_value_derivatives, A, strided_directions);
  for (Index q = 0; q < 3; ++q)
    for (Index i = 0; i < 3; ++i) {
      close(strided_value_derivatives[q, i], batch_dW[q, i], "Incorrect strided batched eigenvalue derivative");
      for (Index j = 0; j < 3; ++j)
        close(strided_derivatives[q, i, j], batch_dP[q, i, j], "Incorrect strided batched eigenvector derivative");
    }
  for (Index q = 0; q < 6; ++q)
    for (Index i = 0; i < 6; ++i)
      for (Index j = 0; j < 7; ++j)
        if (q % 2 != 0 or i % 2 != 0 or j % 2 == 0)
          require(derivative_storage[q, i, j] == batch_guard, "Batched derivative overwrote tensor padding");
  for (Index q = 0; q < 6; ++q)
    for (Index i = 0; i < 7; ++i)
      if (q % 2 != 0 or i % 2 == 0)
        require(eigenvalue_storage[q, i] == batch_guard, "Batched derivative overwrote eigenvalue padding");

  // Phase-aligned central differences independently check the derivative of
  // the actual LAPACK result. No finite difference enters the analytic path.
  Numeric coarse_error = 0, fine_error = 0;
  for (const Numeric h : {1e-3, 2e-4, 4e-5}) {
    ComplexMatrix plus_A(A), minus_A(A), plus_P(3, 3), minus_P(3, 3);
    ComplexVector plus_W(3), minus_W(3);
    for (Index i = 0; i < 3; ++i)
      for (Index j = 0; j < 3; ++j) {
        plus_A[i, j]  += h * dA[i, j];
        minus_A[i, j] -= h * dA[i, j];
      }
    diagonalize(plus_P, plus_W, plus_A);
    diagonalize(minus_P, minus_W, minus_A);
    Numeric error = 0;
    for (Index mode = 0; mode < 3; ++mode) {
      const Index plus = matching(W[mode], plus_W), minus = matching(W[mode], minus_W);
      error                = std::max(error, std::abs((plus_W[plus] - minus_W[minus]) / (2 * h) - dW[mode]));
      Complex plus_overlap = 0, minus_overlap = 0;
      for (Index i = 0; i < 3; ++i) {
        plus_overlap  += std::conj(P[i, mode]) * plus_P[i, plus];
        minus_overlap += std::conj(P[i, mode]) * minus_P[i, minus];
      }
      const Complex plus_phase  = std::conj(plus_overlap) / std::abs(plus_overlap);
      const Complex minus_phase = std::conj(minus_overlap) / std::abs(minus_overlap);
      for (Index i = 0; i < 3; ++i) {
        const Complex fd = (plus_phase * plus_P[i, plus] - minus_phase * minus_P[i, minus]) / (2 * h);
        error            = std::max(error, std::abs(fd - dP[i, mode]));
      }
    }
    if (h == 1e-3) coarse_error = error;
    fine_error = error;
  }
  require(fine_error < 2e-8 and fine_error < 0.01 * coarse_error,
          "Analytic eigendecomposition derivatives do not converge against central differences");

  // Input/output aliases remain safe because all four outputs are committed last.
  auto          alias_A = A, alias_dA = dA;
  ComplexVector alias_W(3), alias_dW(3);
  diagonalize(alias_A, alias_W, alias_dA, alias_dW, alias_A, alias_dA);
  for (Index i = 0; i < 3; ++i) {
    close(alias_W[i], W[i], "Aliased derivative eigenvalues changed");
    close(alias_dW[i], dW[i], "Aliased eigenvalue derivatives changed");
    for (Index j = 0; j < 3; ++j) close(alias_dA[i, j], dP[i, j], "Aliased eigenvector derivatives changed");
  }

  constexpr Complex guard{17, -31};
  ComplexMatrix     backing_A(6, 7, guard), backing_dA(6, 7, guard), backing_P(6, 7, guard), backing_dP(6, 7, guard);
  ComplexVector     backing_W(7, guard), backing_dW(7, guard);
  auto              strided_A  = backing_A[StridedRange(0, 3, 2), StridedRange(1, 3, 2)];
  auto              strided_dA = backing_dA[StridedRange(0, 3, 2), StridedRange(1, 3, 2)];
  auto              strided_P  = backing_P[StridedRange(0, 3, 2), StridedRange(1, 3, 2)];
  auto              strided_dP = backing_dP[StridedRange(0, 3, 2), StridedRange(1, 3, 2)];
  auto              strided_W  = backing_W[StridedRange(1, 3, 2)];
  auto              strided_dW = backing_dW[StridedRange(1, 3, 2)];
  strided_A                    = A;
  strided_dA                   = dA;
  diagonalize(strided_P, strided_W, strided_dP, strided_dW, strided_A, strided_dA, work);
  for (Index i = 0; i < 3; ++i) {
    close(strided_dW[i], dW[i], "Incorrect strided eigenvalue derivative");
    for (Index j = 0; j < 3; ++j) close(strided_dP[i, j], dP[i, j], "Incorrect strided eigenvector derivative");
  }
  for (Index i = 0; i < 6; ++i)
    for (Index j = 0; j < 7; ++j)
      if (i % 2 != 0 or j % 2 == 0)
        require(backing_P[i, j] == guard and backing_dP[i, j] == guard,
                "Strided eigendecomposition derivative overwrote matrix padding");
  for (Index i = 0; i < 7; i += 2)
    require(backing_W[i] == guard and backing_dW[i] == guard,
            "Strided eigendecomposition derivative overwrote vector padding");
}

void derivative_edge_cases() {
  constexpr Complex guard{17, -31};
  ComplexMatrix     P(3, 3, guard), dP(3, 3, guard), A(3, 3, 0), dA(3, 3, 0);
  ComplexVector     W(3, guard), dW(3, guard);
  A[0, 0] = A[1, 1] = 1;
  A[2, 2]           = 3;
  for (Index i = 0; i < 3; ++i) dA[i, i] = 1;
  rejects([&] { diagonalize(P, W, dP, dW, A, dA); }, "Accepted repeated eigenvalue derivatives");
  A[1, 1] = 1 + 1e-15;
  rejects([&] { diagonalize(P, W, dP, dW, A, dA); }, "Accepted unresolved eigenvalue derivatives");
  A[1, 1] = 2;
  A[0, 1] = 1e14;
  rejects([&] { diagonalize(P, W, dP, dW, A, dA); }, "Accepted an unstable derivative eigenvector basis");
  A[0, 1] = 0;
  ComplexMatrix wrong_direction(2, 2);
  rejects([&] { diagonalize(P, W, dP, dW, A, wrong_direction); }, "Accepted incorrect direction dimensions");
  complex_diagonalize_workdata wrong_work(2);
  rejects([&] { diagonalize(P, W, dP, dW, A, dA, wrong_work); }, "Accepted incorrect derivative workspace");
  for (Index i = 0; i < 3; ++i) {
    require(W[i] == guard and dW[i] == guard, "Failed derivative calculation changed output eigenvalues");
    for (Index j = 0; j < 3; ++j)
      require(P[i, j] == guard and dP[i, j] == guard, "Failed derivative calculation changed output eigenvectors");
  }

  // A nonfinite direction is not searched for; it has to reach dP and dW.
  dA[1, 0] = {0, std::numeric_limits<Numeric>::infinity()};
  rejects_or_propagates(
      [&] {
        diagonalize(P, W, dP, dW, A, dA);
        return any_nonfinite(dW) or any_nonfinite(dP);
      },
      "Nonfinite eigendecomposition direction neither rejected nor propagated");
  dA[1, 0] = 0;

  ComplexTensor3 directions(2, 3, 3, 0), batch_dP(2, 3, 3, guard);
  ComplexMatrix  batch_dW(2, 3, guard), wrong_batch_dW(1, 3);
  rejects([&] { diagonalize(P, W, batch_dP, wrong_batch_dW, A, directions); },
          "Accepted mismatched derivative target counts");
  for (Index q = 0; q < 2; ++q)
    for (Index i = 0; i < 3; ++i) {
      require(batch_dW[q, i] == guard, "Failed batch changed eigenvalue derivatives");
      for (Index j = 0; j < 3; ++j) require(batch_dP[q, i, j] == guard, "Failed batch changed eigenvector derivatives");
    }
  directions[1, 1, 1] = {std::numeric_limits<Numeric>::quiet_NaN(), 0};
  rejects_or_propagates(
      [&] {
        diagonalize(P, W, batch_dP, batch_dW, A, directions);
        return any_nonfinite(batch_dW) or any_nonfinite(batch_dP);
      },
      "Nonfinite later derivative target neither rejected nor propagated");
  // An empty target set has primal semantics, including repeated eigenvalues.
  ComplexTensor3 no_directions(0, 3, 3), no_dP(0, 3, 3);
  ComplexMatrix  no_dW(0, 3);
  A[1, 1] = A[0, 0];
  diagonalize(P, W, no_dP, no_dW, A, no_directions);
  eigen_residual(A, P, W);

  ComplexMatrix single_A(1, 1, Complex{6e10, 2e7}), single_dA(1, 1, Complex{3, -4});
  ComplexMatrix single_P(1, 1), single_dP(1, 1);
  ComplexVector single_W(1), single_dW(1);
  diagonalize(single_P, single_W, single_dP, single_dW, single_A, single_dA);
  close(single_dW[0], single_dA[0, 0], "Incorrect single-mode eigenvalue derivative");
  close(single_dP[0, 0], 0, "Single-mode eigenvector derivative is not zero");
  ComplexMatrix empty_A(0, 0), empty_P(0, 0), empty_dP(0, 0);
  ComplexVector empty_W(0), empty_dW(0);
  diagonalize(empty_P, empty_W, empty_dP, empty_dW, empty_A, empty_A);
}

}  // namespace

int main() try {
  ordinary();
  strided();
  multiple_rhs_solve();
  invalid_and_empty();
  analytic_derivatives();
  derivative_edge_cases();
  std::cout << "Complex LAPACK solve and diagonalization checks passed.\n";
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return 1;
}
