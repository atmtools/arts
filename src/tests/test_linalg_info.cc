#include <lin_alg.h>

#include <cmath>
#include <exception>
#include <iostream>
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

}  // namespace

int main() {
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
