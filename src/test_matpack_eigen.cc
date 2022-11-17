#include <cstdlib>
#include <iomanip>
#include <ios>
#define EIGEN_RUNTIME_NO_MALLOC

#include <string>
#include <vector>

#include "matpack_eigen.h"

class TestDiag {
  struct Unit {
    std::string name;
    bool error;
  };

  std::vector<std::string> test;
  std::vector<std::vector<Unit>> unit_result;

public:
  void add_test(const char *test_name) {
    test.emplace_back(test_name);
    unit_result.emplace_back();
  }

  void run_test(const char *test_desc, bool no_error) {
    unit_result.back().emplace_back(Unit{test_desc, not no_error});
  }

  bool print_and_return_error() {
    bool any_error{false};

    for (std::size_t i = 0; i < test.size(); i++) {
      std::cout << std::setfill('#') << std::setw(80) << '\n';
      std::cout << std::setfill(' ') << std::setw(0) << test[i] << ":\n";
      for (auto &[name, error] : unit_result[i]) {
        std::cout << "    " << name
                  << std::setw(56 - static_cast<int>(name.size())) << " : "
                  << (error ? "Fail" : "Pass") << '\n';
        any_error = any_error or error;
      }
      if (i == test.size() - 1)
        std::cout << std::setfill('#') << std::setw(80) << '\n';
    }

    return any_error;
  }
};

using namespace matpack::eigen;

void test_col_and_row_vec(TestDiag &tests) {
  tests.add_test("test_col_and_row_vec");

  {
    Vector arts({1, 2, 3, 4, 5});
    Eigen::VectorXd eigen(5);
    eigen << 1, 2, 3, 4, 5;
    tests.run_test("Real    eigen == row_vec(arts)", eigen == row_vec(arts));
    tests.run_test("Real    eigen == col_vec(arts)", eigen == col_vec(arts));
  }

  {
    ComplexVector arts({Complex{1, 2}, Complex{3, 4}, Complex{5, 6},
                        Complex{7, 8}, Complex{9, 10}});
    Eigen::VectorXcd eigen(5);
    eigen << Complex{1, 2}, Complex{3, 4}, Complex{5, 6}, Complex{7, 8},
        Complex{9, 10};
    tests.run_test("Complex eigen == row_vec(arts)", eigen == row_vec(arts));
    tests.run_test("Complex eigen == col_vec(arts)", eigen == col_vec(arts));
  }
}

void test_multiplication_and_setting(TestDiag &tests) {
  tests.add_test("test_multiplication_and_setting");

  Matrix A(3, 4, 1);
  Vector x(4, 2);

  Matrix XTX;
  Vector y;

  y = A * x;
  XTX = A * transpose(A);

  Eigen::internal::set_is_malloc_allowed(false);
  as_eigen(y).noalias() = A * x;
  mat(XTX).noalias() = A * transpose(A);
  Eigen::internal::set_is_malloc_allowed(true);

  tests.run_test("y = A * x; y.size() == 3", y.size() == 3);
  tests.run_test("y[0] == 8", (A * x)[0] == 8);
  tests.run_test("A * x == as_eigen(Vector({8, 8, 8}))",
                 A * x == as_eigen(Vector({8, 8, 8})));
  tests.run_test("XTX = A * transpose(A); XTX.shape() == Shape<2>{3, 3}",
                 XTX.shape() == Shape<2>{3, 3});
  tests.run_test("A * transpose(A) == as_eigen(Matrix(3, 3, 4))",
                 A * transpose(A) == as_eigen(Matrix(3, 3, 4)));
  tests.run_test("No alloc (as_eigen(y).noalias() = A * x)", true);
  tests.run_test("No alloc (mat(XTX).noalias() = A * transpose(A))", true);
}

void test_chaining_expressions(TestDiag &tests) {
  tests.add_test("test_chaining_expressions");

  Matrix A(3, 4, 1);
  Vector x({1, 2, 3, 4}), b({5, 6, 7}), y{A * x + b};

  tests.run_test("A * x + b == as_eigen(Vector{15, 16, 17})",
                 A * x + b == as_eigen(Vector{15, 16, 17}));
  tests.run_test("b + A * x == as_eigen(Vector{15, 16, 17})",
                 b + A * x == as_eigen(Vector{15, 16, 17}));
  tests.run_test("A * x - b == as_eigen(Vector{5, 4, 3})",
                 A * x - b == as_eigen(Vector{5, 4, 3}));
  tests.run_test("b - A * x == as_eigen(Vector{-5, -4, -3})",
                 b - A * x == as_eigen(Vector{-5, -4, -3}));
  tests.run_test("5*y == 5*as_eigen(Vector{15, 16, 17})",
                 5 * y == 5 * as_eigen(Vector{15, 16, 17}));
  tests.run_test("-1.2*y == -1.2*as_eigen(Vector{15, 16, 17})",
                 -1.2 * y == -1.2 * as_eigen(Vector{15, 16, 17}));
  tests.run_test("5*A == 5*as_eigen(Matrix{3, 4, 1})",
                 5 * A == 5 * as_eigen(Matrix{3, 4, 1}));
  tests.run_test("-1.2*A == -1.2*as_eigen(Matrix{3, 4, 1})",
                 -1.2 * A == -1.2 * as_eigen(Matrix{3, 4, 1}));
}

int main() {
  TestDiag tests;
  test_col_and_row_vec(tests);
  test_multiplication_and_setting(tests);
  test_chaining_expressions(tests);
  if (tests.print_and_return_error())
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
