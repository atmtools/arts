#include <lin_alg.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <string_view>

namespace {

/** Construct an exactly singular two-by-two matrix. */
Matrix singular_matrix() {
  Matrix matrix(2, 2);
  matrix[0, 0] = 1.0;
  matrix[0, 1] = 2.0;
  matrix[1, 0] = 2.0;
  matrix[1, 1] = 4.0;
  return matrix;
}

/** Return whether an operation reports the expected singular LU factorization. */
template <typename Operation> bool reports_singular_matrix(Operation&& operation) {
  try {
    operation();
  } catch (const std::exception& error) {
    const std::string_view message{error.what()};
    return message.contains("DGETRF") and message.contains("singular matrix");
  }
  return false;
}

/** Check that a regular in-place solve is unaffected by LAPACK status checks. */
bool test_regular_solve() {
  Matrix matrix(2, 2);
  matrix[0, 0] = 3.0;
  matrix[0, 1] = 1.0;
  matrix[1, 0] = 1.0;
  matrix[1, 1] = 2.0;

  Vector solution{9.0, 8.0};
  solve_inplace(solution, matrix);

  constexpr Numeric tolerance = 1e-13;
  return std::abs(solution[0] - 2.0) < tolerance and std::abs(solution[1] - 3.0) < tolerance;
}

/** Check that the allocating and in-place solve paths reject singular input. */
bool test_singular_solves() {
  const Matrix singular = singular_matrix();
  const Vector rhs{1.0, 2.0};

  const bool allocating = reports_singular_matrix([&] {
    Vector solution(2);
    solve(solution, singular, rhs);
  });

  const bool inplace = reports_singular_matrix([&] {
    Matrix matrix   = singular;
    Vector solution = rhs;
    solve_inplace(solution, matrix);
  });

  return allocating and inplace;
}

/** Check that matrix inversion rejects a singular LU factorization. */
bool test_singular_inverse() {
  return reports_singular_matrix([] {
    Matrix singular = singular_matrix();
    inv_inplace(singular);
  });
}

/** Check B A^-1 in place against a known solution, for a non-square B */
bool test_right_solve() {
  ComplexMatrix A(3, 3), X(2, 3);
  for (Index i = 0; i < 3; ++i)
    for (Index j = 0; j < 3; ++j) A[i, j] = Complex{i == j ? 5.0 : 1.0, 0.5 * static_cast<Numeric>(i - j)};
  for (Index i = 0; i < 2; ++i)
    for (Index j = 0; j < 3; ++j) X[i, j] = Complex{static_cast<Numeric>(i + j), static_cast<Numeric>(i - 2 * j)};

  ComplexMatrix B(2, 3, Complex{});
  for (Index i = 0; i < 2; ++i)
    for (Index j = 0; j < 3; ++j)
      for (Index k = 0; k < 3; ++k) B[i, j] += X[i, k] * A[k, j];

  solve_workdata wo;
  right_solve_inplace(B, A, wo);
  Numeric error = 0.0;
  for (Index i = 0; i < 2; ++i)
    for (Index j = 0; j < 3; ++j) error = std::max(error, std::abs(B[i, j] - X[i, j]));
  return error < 1e-13;
}

/** matrix_exp against analytic exponentials: exp([[0, -t], [t, 0]]) = [[cos t, -sin t], [sin t, cos t]] and that of
 * a diagonal.  matrix_exp scales A by 2^-j, the binary exponent of its norm, to a norm below 1, sums its Pade
 * approximant to a few epsilon, and squares j times, each squaring about doubling the relative error, so the bound
 * is 16 epsilon 2^j. */
bool test_matrix_exp() {
  const auto check = [](const Matrix& A, const Matrix& expected) {
    Matrix F(A.shape());
    matrix_exp(F, A);
    int j = 0;
    std::frexp(norm_inf(A), &j);
    const Numeric bound = std::ldexp(16 * std::numeric_limits<Numeric>::epsilon(), std::max(0, j));
    Numeric       worst = 0.0;
    for (Index r = 0; r < A.nrows(); ++r)
      for (Index c = 0; c < A.ncols(); ++c) worst = std::max(worst, std::abs(F[r, c] - expected[r, c]));
    return worst <= bound;
  };
  bool ok = true;
  for (const Numeric t : {0.3, 2.5}) {
    Matrix A(2, 2, 0.0), E(2, 2);
    A[0, 1] = -t;
    A[1, 0] = t;
    E[0, 0] = E[1, 1] = std::cos(t);
    E[0, 1]           = -std::sin(t);
    E[1, 0]           = std::sin(t);
    ok                = ok and check(A, E);
  }
  Matrix D(2, 2, 0.0), E(2, 2, 0.0);
  D[0, 0] = -1.2;
  D[1, 1] = 0.7;
  E[0, 0] = std::exp(-1.2);
  E[1, 1] = std::exp(0.7);
  return ok and check(D, E);
}
}  // namespace

int main() {
  if (not test_matrix_exp()) {
    std::cerr << "The matrix exponential differs from the analytic one beyond its rounding.\n";
    return 1;
  }
  if (not test_right_solve()) {
    std::cerr << "The right division B A^-1 returned the wrong solution.\n";
    return 1;
  }
  if (not test_regular_solve()) {
    std::cerr << "Regular LU solve returned the wrong solution.\n";
    return 1;
  }
  if (not test_singular_solves()) {
    std::cerr << "A singular LU solve was not rejected by DGETRF.\n";
    return 1;
  }
  if (not test_singular_inverse()) {
    std::cerr << "A singular matrix inverse was not rejected by DGETRF.\n";
    return 1;
  }
}
