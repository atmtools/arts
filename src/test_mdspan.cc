#include "debug.h"
#include "matpack_mdspan.h"
#include <iostream>

#include <cstdlib>
#include <cxxabi.h>
#include <memory>
#include <typeinfo>
std::string demangle(const char *name) {
  int status = -4; // some arbitrary value to eliminate the compiler warning

  // enable c++11 by passing the flag -std=c++11 to g++
  std::unique_ptr<char, void (*)(void *)> res{
      abi::__cxa_demangle(name, nullptr, nullptr, &status), std::free};

  return (status == 0) ? res.get() : name;
}

template <class T> std::string type(const T &t) {
  return demangle(typeid(t).name());
}

Matrix build_test_matrix(Index rows, Index cols) {
  Matrix a(rows, cols);
  for (Index i = 0; i < rows; i++) {
    for (Index j = 0; j < cols; j++) {
      a(i, j) = static_cast<Numeric>(10 * i + j + 1);
    }
  }
  return a;
}

void test_mult() {
  std::cout << "# MULT TEST ######################################\n";
  {
    std::cout << "TEST 1 (y = A * x; NOT STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2, 3});
    Matrix A (3, 3, 1);
    mult_fast(y, A, x);
    std::cout << "A (in):\n" << A << '\n'
              << "x (in):\n" << x << '\n'
              << "y (out):\n" << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "6 6 6", "Bad mult")
  }
  {
    std::cout << "TEST 2 (y = A * x; NOT STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2, 3});
    Matrix A = build_test_matrix(3, 3);
    mult_fast(y, A, x);
    std::cout << "A (in):\n" << A << '\n'
              << "x (in):\n" << x << '\n'
              << "y (out):\n" << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "14 74 134", "Bad mult")
  }
  {
    std::cout << "TEST 3 (y = A * x; NOT STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2});
    Matrix A = build_test_matrix(3, 2);
    mult_fast(y, A, x);
    std::cout << "A (in):\n" << A << '\n'
              << "x (in):\n" << x << '\n'
              << "y (out):\n" << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "5 35 65", "Bad mult")
  }
    {
    std::cout << "TEST 4 (y = A * x; STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2, 3});
    Matrix A (3, 3, 1);
    mult(y, A, x);
    std::cout << "A (in):\n" << A << '\n'
              << "x (in):\n" << x << '\n'
              << "y (out):\n" << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "6 6 6", "Bad mult")
  }
  {
    std::cout << "TEST 5 (y = A * x; STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2, 3});
    Matrix A = build_test_matrix(3, 3);
    mult(y, A, x);
    std::cout << "A (in):\n" << A << '\n'
              << "x (in):\n" << x << '\n'
              << "y (out):\n" << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "14 74 134", "Bad mult")
  }
  {
    std::cout << "TEST 6 (y = A * x; STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2});
    Matrix A = build_test_matrix(3, 2);
    mult(y, A, x);
    std::cout << "A (in):\n" << A << '\n'
              << "x (in):\n" << x << '\n'
              << "y (out):\n" << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "5 35 65", "Bad mult")
  }
  {
    std::cout << "TEST 7 (C = A * B; NON STRIDED):\n";
    Matrix C(3, 3);
    Matrix B(3, 3, 1);
    Matrix A (3, 3, 1);
    mult_fast(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "3 3 3\n3 3 3\n3 3 3", "Bad mult")
  }
  {
    std::cout << "TEST 8 (C = A * B; NON STRIDED):\n";
    Matrix C(3, 4);
    Matrix A (3, 2, 1);
    Matrix B(2, 4, 1);
    mult_fast(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "2 2 2 2\n2 2 2 2\n2 2 2 2", "Bad mult")
  }
  {
    std::cout << "TEST 9 (C = A * B; NON STRIDED):\n";
    Matrix C(3, 4);
    Matrix A = build_test_matrix(3, 2);
    Matrix B = build_test_matrix(2, 4);
    mult_fast(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "23 26 29 32\n143 166 189 212\n263 306 349 392", "Bad mult")
  }
  {
    std::cout << "TEST 10 (C = A * B; NON STRIDED):\n";
    Matrix C(4, 3);
    Matrix A = build_test_matrix(4, 2);
    Matrix B = build_test_matrix(2, 3);
    mult_fast(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "23 26 29\n143 166 189\n263 306 349\n383 446 509", "Bad mult")
  }
  {
    std::cout << "TEST 11 (C = A * B; STRIDED):\n";
    Matrix C(3, 3);
    Matrix B(3, 3, 1);
    Matrix A (3, 3, 1);
    mult(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "3 3 3\n3 3 3\n3 3 3", "Bad mult")
  }
  {
    std::cout << "TEST 12 (C = A * B; STRIDED):\n";
    Matrix C(3, 4);
    Matrix A (3, 2, 1);
    Matrix B(2, 4, 1);
    mult(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "2 2 2 2\n2 2 2 2\n2 2 2 2", "Bad mult")
  }
  {
    std::cout << "TEST 13 (C = A * B; STRIDED):\n";
    Matrix C(3, 4);
    Matrix A = build_test_matrix(3, 2);
    Matrix B = build_test_matrix(2, 4);
    mult(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "23 26 29 32\n143 166 189 212\n263 306 349 392", "Bad mult")
  }
  {
    std::cout << "TEST 14 (C = A * B; STRIDED):\n";
    Matrix C(4, 3);
    Matrix A = build_test_matrix(4, 2);
    Matrix B = build_test_matrix(2, 3);
    mult(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "23 26 29\n143 166 189\n263 306 349\n383 446 509", "Bad mult")
  }
  std::cout << "#/MULT TEST ######################################\n";
}

ComplexMatrix build_test_complex_matrix(Index rows, Index cols) {
  ComplexMatrix a(rows, cols);
  for (Index i = 0; i < rows; i++) {
    for (Index j = 0; j < cols; j++) {
      auto n = static_cast<Numeric>(10 * i + j + 1);
      a(i, j) = Complex{n, 2*n};
    }
  }
  return a;
}

void test_complex_mult() {
  std::cout << "# COMPLEX MULT TEST ######################################\n";
  {
    std::cout << "TEST 1 (y = A * x; NOT STRIDED):\n";
    ComplexVector y(3);
    ComplexVector x(std::vector<Complex>{Complex{1, 1}, Complex{2, 2}, Complex{3, 3}});
    ComplexMatrix A (3, 3, Complex{1, 1});
    mult_fast(y, A, x);
    std::cout << "A (in):\n" << A << '\n'
              << "x (in):\n" << x << '\n'
              << "y (out):\n" << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "(0,12) (0,12) (0,12)", "Bad mult")
  }
  {
    std::cout << "TEST 2 (y = A * x; NOT STRIDED):\n";
    ComplexVector y(3);
    ComplexVector x(std::vector<Complex>{Complex{1, 1}, Complex{2, 2}, Complex{3, 3}});
    ComplexMatrix A = build_test_complex_matrix(3, 3);
    mult_fast(y, A, x);
    std::cout << "A (in):\n" << A << '\n'
              << "x (in):\n" << x << '\n'
              << "y (out):\n" << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "(-14,42) (-74,222) (-134,402)", "Bad mult")
  }
  {
    std::cout << "TEST 3 (y = A * x; NOT STRIDED):\n";
    ComplexVector y(3);
    ComplexVector x(std::vector<Complex>{Complex{1, 1}, Complex{2, 2}});
    ComplexMatrix A = build_test_complex_matrix(3, 2);
    mult_fast(y, A, x);
    std::cout << "A (in):\n" << A << '\n'
              << "x (in):\n" << x << '\n'
              << "y (out):\n" << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "(-5,15) (-35,105) (-65,195)", "Bad mult")
  }
    {
    std::cout << "TEST 4 (y = A * x; STRIDED):\n";
    ComplexVector y(3);
    ComplexVector x(std::vector<Complex>{Complex{1, 1}, Complex{2, 2}, Complex{3, 3}});
    ComplexMatrix A (3, 3, Complex{1, 1});
    mult(y, A, x);
    std::cout << "A (in):\n" << A << '\n'
              << "x (in):\n" << x << '\n'
              << "y (out):\n" << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "(0,12) (0,12) (0,12)", "Bad mult")
  }
  {
    std::cout << "TEST 5 (y = A * x; STRIDED):\n";
    ComplexVector y(3);
    ComplexVector x(std::vector<Complex>{Complex{1, 1}, Complex{2, 2}, Complex{3, 3}});
    ComplexMatrix A = build_test_complex_matrix(3, 3);
    mult(y, A, x);
    std::cout << "A (in):\n" << A << '\n'
              << "x (in):\n" << x << '\n'
              << "y (out):\n" << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "(-14,42) (-74,222) (-134,402)", "Bad mult")
  }
  {
    std::cout << "TEST 6 (y = A * x; STRIDED):\n";
    ComplexVector y(3);
    ComplexVector x(std::vector<Complex>{Complex{1, 1}, Complex{2, 2}});
    ComplexMatrix A = build_test_complex_matrix(3, 2);
    mult(y, A, x);
    std::cout << "A (in):\n" << A << '\n'
              << "x (in):\n" << x << '\n'
              << "y (out):\n" << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "(-5,15) (-35,105) (-65,195)", "Bad mult")
  }
  {
    std::cout << "TEST 7 (C = A * B; NON STRIDED):\n";
    ComplexMatrix C(3, 3);
    ComplexMatrix B(3, 3, Complex{1, 1});
    ComplexMatrix A (3, 3, Complex{1, 1});
    mult_fast(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "(0,6) (0,6) (0,6)\n(0,6) (0,6) (0,6)\n(0,6) (0,6) (0,6)", "Bad mult")
  }
  {
    std::cout << "TEST 8 (C = A * B; NON STRIDED):\n";
    ComplexMatrix C(3, 4);
    ComplexMatrix A (3, 2, Complex{1, 1});
    ComplexMatrix B(2, 4, Complex{1, 1});
    mult_fast(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "(0,4) (0,4) (0,4) (0,4)\n(0,4) (0,4) (0,4) (0,4)\n(0,4) (0,4) (0,4) (0,4)", "Bad mult")
  }
  {
    std::cout << "TEST 9 (C = A * B; NON STRIDED):\n";
    ComplexMatrix C(3, 4);
    ComplexMatrix A = build_test_complex_matrix(3, 2);
    ComplexMatrix B = build_test_complex_matrix(2, 4);
    mult_fast(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "(-69,92) (-78,104) (-87,116) (-96,128)\n(-429,572) (-498,664) (-567,756) (-636,848)\n(-789,1052) (-918,1224) (-1047,1396) (-1176,1568)", "Bad mult")
  }
  {
    std::cout << "TEST 10 (C = A * B; NON STRIDED):\n";
    ComplexMatrix C(4, 3);
    ComplexMatrix A = build_test_complex_matrix(4, 2);
    ComplexMatrix B = build_test_complex_matrix(2, 3);
    mult_fast(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "(-69,92) (-78,104) (-87,116)\n(-429,572) (-498,664) (-567,756)\n(-789,1052) (-918,1224) (-1047,1396)\n(-1149,1532) (-1338,1784) (-1527,2036)", "Bad mult")
  }
  {
    std::cout << "TEST 11 (C = A * B; STRIDED):\n";
    ComplexMatrix C(3, 3);
    ComplexMatrix B(3, 3, Complex{1, 1});
    ComplexMatrix A (3, 3, Complex{1, 1});
    mult(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "(0,6) (0,6) (0,6)\n(0,6) (0,6) (0,6)\n(0,6) (0,6) (0,6)", "Bad mult")
  }
  {
    std::cout << "TEST 12 (C = A * B; STRIDED):\n";
    ComplexMatrix C(3, 4);
    ComplexMatrix A (3, 2, Complex{1, 1});
    ComplexMatrix B(2, 4, Complex{1, 1});
    mult(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "(0,4) (0,4) (0,4) (0,4)\n(0,4) (0,4) (0,4) (0,4)\n(0,4) (0,4) (0,4) (0,4)", "Bad mult")
  }
  {
    std::cout << "TEST 13 (C = A * B; STRIDED):\n";
    ComplexMatrix C(3, 4);
    ComplexMatrix A = build_test_complex_matrix(3, 2);
    ComplexMatrix B = build_test_complex_matrix(2, 4);
    mult(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "(-69,92) (-78,104) (-87,116) (-96,128)\n(-429,572) (-498,664) (-567,756) (-636,848)\n(-789,1052) (-918,1224) (-1047,1396) (-1176,1568)", "Bad mult")
  }
  {
    std::cout << "TEST 14 (C = A * B; STRIDED):\n";
    ComplexMatrix C(4, 3);
    ComplexMatrix A = build_test_complex_matrix(4, 2);
    ComplexMatrix B = build_test_complex_matrix(2, 3);
    mult(C, A, B);
    std::cout << "A (in):\n" << A << '\n'
              << "B (in):\n" << B << '\n'
              << "C (out):\n" << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "(-69,92) (-78,104) (-87,116)\n(-429,572) (-498,664) (-567,756)\n(-789,1052) (-918,1224) (-1047,1396)\n(-1149,1532) (-1338,1784) (-1527,2036)", "Bad mult")
  }
  std::cout << "#/COMPLEX MULT TEST ######################################\n";
}

void test_transpose() {
  std::cout << "# TRANSPOSE TEST ######################################\n";
  {
    std::cout << "TEST 1:\n";
    auto A = build_test_matrix(3, 3);
    std::cout << A << '\n';
    ConstMatrixView X=A;
    std::cout << X.transpose() << '\n';
    std::cout << transpose(X).transpose() << '\n';
  }
  {
    std::cout << "TEST 2:\n";
    auto A = build_test_matrix(2, 3);
    std::cout << A << '\n';
    ConstMatrixView X=A;
    std::cout << X.transpose() << '\n';
    std::cout << transpose(X).transpose() << '\n';
  }
  {
    std::cout << "TEST 3:\n";
    auto A = build_test_matrix(4, 3);
    std::cout << A << '\n';
    ConstMatrixView X=A;
    std::cout << X.transpose() << '\n';
    std::cout << transpose(X).transpose() << '\n';
  }
  std::cout << "#/TRANSPOSE TEST ######################################\n";
}

int main() {
  using namespace matpack::md;

  auto x = Vector(4);
  const auto y = Matrix(4, 3, 4);
  auto yc = Matrix{};
  auto z = Tensor3(4, 3, 2, 4);

  z(1, 2, 1) += 1;
  z(1, 1, 1) += 1;

  for (auto v : z) {
    v(0, 0) += 3;
    for (auto t : v) {
      t[1] += 5;
      for (auto &s : t)
        s *= 2;
    }
  }
  for (auto &v : x)
    v += 2;

  x = std::move(z).flatten();
  std::cout << "SHOULD BE SAME (TOP IS STORED CONST CHAR *):\n14 18 8 18 8 18 "
               "14 18 8 20 8 20 14 18 8 18 8 "
               "18 14 18 8 18 8 18\n";
  std::cout << x << '\n';
  std::cout << z << '\n';
  z = std::move(std::move(x)).reshape(2, 3, 4);
  std::cout << x << '\n';
  std::cout << z << '\n';
  z = std::move(std::move(z).flatten()).reshape(4, 3, 2);
  std::cout << x << '\n';
  std::cout << z << '\n';
  yc = std::move(std::move(z).flatten()).reshape(4, 6);
  std::cout << x << '\n';
  std::cout << yc << '\n';

  z = simple_view<double, 3, true>{std::move(yc).reshape(4,3,2)};
  std::cout << z << '\n';

  x = std::move(z).flatten();
  std::cout << x.sum() << ' ' << x * x << '\n';

  std::cout << ConstVectorView{x} * x << ' ' << x * ConstVectorView{x} << '\n';

  Matrix G{4, 5};
  for (Index i = 0; i < 4; i++)
    for (Index j = 0; j < 5; j++)
      G(i, j) = static_cast<Numeric>(10 * i + j);
  VectorView a = G(matpack::md::joker, 0);
  FastVectorView b = G(0, matpack::md::joker);
  std::cout << "THE NEXT TWO ROWS SHOULD BE IDENTICAL (STRIDED):\n"
            << a << '\n'
            << "0 10 20 30\n";
  std::cout << "THE NEXT TWO ROWS SHOULD BE IDENTICAL (NON-STRIDED):\n"
            << b << '\n'
            << "0 1 2 3 4\n";
  std::cout << "TEST GOOD (IF 6-LINES ABOVE LOOKS GOOD): Could get strided and non-strided data\n";
  try {
    auto f = FastVectorView{a};
    std::cout << "THE NEXT TWO ROWS SHOULD BE IDENTICAL:\n"
              << f << '\n'
              << "0 10 20 30\n";
    std::cout << "ERROR: You should not be able to do this" << '\n';
    std::exit(1);
  } catch (...) {
    std::cout << "TEST GOOD: Could not convert between non-strided and strided "
                 "view\n";
  }

  test_mult();
  test_complex_mult();
  test_transpose();
}
