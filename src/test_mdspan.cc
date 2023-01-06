#include "debug.h"
#include "matpack_mdspan.h"
#include "matpack_mdspan_eigen.h"

#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::runtime_error;

Numeric by_reference(const Numeric &x) { return x + 1; }

Numeric by_value(Numeric x) { return x + 1; }

void fill_with_junk(VectorView x) { x = 999; }

void fill_with_junk(MatrixView x) { x = 888; }

int test1() {
  Vector v(20);

  cout << "v.nelem() = " << v.nelem() << "\n";

  for (Index i = 0; i < v.nelem(); ++i)
    v[i] = (Numeric)i;

  cout << "v.begin() = " << *v.begin() << "\n";

  cout << "v = \n" << v << "\n";

  fill_with_junk(v[Range(1, 8, 2)][Range(2, joker)]);
  //  fill_with_junk(v);

  Vector v2 = v[Range(2, 4)];

  cout << "v2 = \n" << v2 << "\n";

  for (Index i = 0; i < 1000; ++i) {
    Vector v3(1000);
    v3 = (Numeric)i;
  }

  v2[Range(joker)] = 88;

  v2[Range(0, 2)] = 77;

  cout << "v = \n" << v << "\n";
  cout << "v2 = \n" << v2 << "\n";
  cout << "v2.nelem() = \n" << v2.nelem() << "\n";

  Vector v3;
  v3.resize(v2.nelem());
  v3 = v2;

  cout << "\nv3 = \n" << v3 << "\n";
  fill_with_junk((VectorView)v2);
  cout << "\nv3 after junking v2 = \n" << v3 << "\n";
  v3 *= 2;
  cout << "\nv3 after *2 = \n" << v3 << "\n";

  Matrix M(10, 15);
  {
    Numeric n = 0;
    for (Index i = 0; i < M.nrows(); ++i)
      for (Index j = 0; j < M.ncols(); ++j)
        M(i, j) = ++n;
  }

  cout << "\nM =\n" << M << "\n";

  cout << "\nM(Range(2,4),Range(2,4)) =\n"
       << M(Range(2, 4), Range(2, 4)) << "\n";

  cout << "\nM(Range(2,4),Range(2,4))(Range(1,2),Range(1,2)) =\n"
       << M(Range(2, 4), Range(2, 4))(Range(1, 2), Range(1, 2)) << "\n";

  cout << "\nM(1,Range(joker)) =\n" << M(1, Range(joker)) << "\n";

  cout << "\nFilling M(1,Range(1,2)) with junk.\n";
  fill_with_junk(M(1, Range(1, 2)));

  cout << "\nM(Range(0,4),Range(0,4)) =\n"
       << M(Range(0, 4), Range(0, 4)) << "\n";

  cout << "\nFilling M(Range(4,2,2),Range(6,3)) with junk.\n";

  MatrixView s = M(Range(4, 2, 2), Range(6, 3));
  fill_with_junk(s);

  cout << "\nM =\n" << M << "\n";

  const Matrix C = M;

  cout << "\nC(Range(3,4,2),Range(2,3,3)) =\n"
       << C(Range(3, 4, 2), Range(2, 3, 3)) << "\n";

  cout << "\nC(Range(3,4,2),Range(2,3,3)).transpose() =\n"
       << transpose(C(Range(3, 4, 2), Range(2, 3, 3))) << "\n";

  return 0;
}

void test2() {
  Vector v(50000000);

  cout << "v.nelem() = " << v.nelem() << "\n";

  cout << "Filling\n";
  //   for (Index i=0; i<v.nelem(); ++i )
  //     v[i] = sqrt(i);
  v = 1.;
  cout << "Done\n";
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
    Matrix A(3, 3, 1);
    mult_fast(y, A, x);
    std::cout << "A (in):\n"
              << A << '\n'
              << "x (in):\n"
              << x << '\n'
              << "y (out):\n"
              << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "6 6 6", "Bad mult")
  }
  {
    std::cout << "TEST 2 (y = A * x; NOT STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2, 3});
    Matrix A = build_test_matrix(3, 3);
    mult_fast(y, A, x);
    std::cout << "A (in):\n"
              << A << '\n'
              << "x (in):\n"
              << x << '\n'
              << "y (out):\n"
              << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "14 74 134", "Bad mult")
  }
  {
    std::cout << "TEST 3 (y = A * x; NOT STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2});
    Matrix A = build_test_matrix(3, 2);
    mult_fast(y, A, x);
    std::cout << "A (in):\n"
              << A << '\n'
              << "x (in):\n"
              << x << '\n'
              << "y (out):\n"
              << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "5 35 65", "Bad mult")
  }
  {
    std::cout << "TEST 4 (y = A * x; STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2, 3});
    Matrix A(3, 3, 1);
    mult(y, A, x);
    std::cout << "A (in):\n"
              << A << '\n'
              << "x (in):\n"
              << x << '\n'
              << "y (out):\n"
              << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "6 6 6", "Bad mult")
  }
  {
    std::cout << "TEST 5 (y = A * x; STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2, 3});
    Matrix A = build_test_matrix(3, 3);
    mult(y, A, x);
    std::cout << "A (in):\n"
              << A << '\n'
              << "x (in):\n"
              << x << '\n'
              << "y (out):\n"
              << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "14 74 134", "Bad mult")
  }
  {
    std::cout << "TEST 6 (y = A * x; STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2});
    Matrix A = build_test_matrix(3, 2);
    mult(y, A, x);
    std::cout << "A (in):\n"
              << A << '\n'
              << "x (in):\n"
              << x << '\n'
              << "y (out):\n"
              << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "5 35 65", "Bad mult")
  }
  {
    std::cout << "TEST 7 (C = A * B; NON STRIDED):\n";
    Matrix C(3, 3);
    Matrix B(3, 3, 1);
    Matrix A(3, 3, 1);
    mult_fast(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "3 3 3\n3 3 3\n3 3 3", "Bad mult")
  }
  {
    std::cout << "TEST 8 (C = A * B; NON STRIDED):\n";
    Matrix C(3, 4);
    Matrix A(3, 2, 1);
    Matrix B(2, 4, 1);
    mult_fast(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "2 2 2 2\n2 2 2 2\n2 2 2 2",
                       "Bad mult")
  }
  {
    std::cout << "TEST 9 (C = A * B; NON STRIDED):\n";
    Matrix C(3, 4);
    Matrix A = build_test_matrix(3, 2);
    Matrix B = build_test_matrix(2, 4);
    mult_fast(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq
                           "23 26 29 32\n143 166 189 212\n263 306 349 392",
                       "Bad mult")
  }
  {
    std::cout << "TEST 10 (C = A * B; NON STRIDED):\n";
    Matrix C(4, 3);
    Matrix A = build_test_matrix(4, 2);
    Matrix B = build_test_matrix(2, 3);
    mult_fast(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq
                           "23 26 29\n143 166 189\n263 306 349\n383 446 509",
                       "Bad mult")
  }
  {
    std::cout << "TEST 11 (C = A * B; STRIDED):\n";
    Matrix C(3, 3);
    Matrix B(3, 3, 1);
    Matrix A(3, 3, 1);
    mult(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "3 3 3\n3 3 3\n3 3 3", "Bad mult")
  }
  {
    std::cout << "TEST 12 (C = A * B; STRIDED):\n";
    Matrix C(3, 4);
    Matrix A(3, 2, 1);
    Matrix B(2, 4, 1);
    mult(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq "2 2 2 2\n2 2 2 2\n2 2 2 2",
                       "Bad mult")
  }
  {
    std::cout << "TEST 13 (C = A * B; STRIDED):\n";
    Matrix C(3, 4);
    Matrix A = build_test_matrix(3, 2);
    Matrix B = build_test_matrix(2, 4);
    mult(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq
                           "23 26 29 32\n143 166 189 212\n263 306 349 392",
                       "Bad mult")
  }
  {
    std::cout << "TEST 14 (C = A * B; STRIDED):\n";
    Matrix C(4, 3);
    Matrix A = build_test_matrix(4, 2);
    Matrix B = build_test_matrix(2, 3);
    mult(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq
                           "23 26 29\n143 166 189\n263 306 349\n383 446 509",
                       "Bad mult")
  }
  std::cout << "#/MULT TEST ######################################\n";
}

ComplexMatrix build_test_complex_matrix(Index rows, Index cols) {
  ComplexMatrix a(rows, cols);
  for (Index i = 0; i < rows; i++) {
    for (Index j = 0; j < cols; j++) {
      auto n = static_cast<Numeric>(10 * i + j + 1);
      a(i, j) = Complex{n, 2 * n};
    }
  }
  return a;
}

void test_complex_mult() {
  std::cout << "# COMPLEX MULT TEST ######################################\n";
  {
    std::cout << "TEST 1 (y = A * x; NOT STRIDED):\n";
    ComplexVector y(3);
    ComplexVector x(
        std::vector<Complex>{Complex{1, 1}, Complex{2, 2}, Complex{3, 3}});
    ComplexMatrix A(3, 3, Complex{1, 1});
    mult_fast(y, A, x);
    std::cout << "A (in):\n"
              << A << '\n'
              << "x (in):\n"
              << x << '\n'
              << "y (out):\n"
              << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "(0,12) (0,12) (0,12)", "Bad mult")
  }
  {
    std::cout << "TEST 2 (y = A * x; NOT STRIDED):\n";
    ComplexVector y(3);
    ComplexVector x(
        std::vector<Complex>{Complex{1, 1}, Complex{2, 2}, Complex{3, 3}});
    ComplexMatrix A = build_test_complex_matrix(3, 3);
    mult_fast(y, A, x);
    std::cout << "A (in):\n"
              << A << '\n'
              << "x (in):\n"
              << x << '\n'
              << "y (out):\n"
              << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "(-14,42) (-74,222) (-134,402)",
                       "Bad mult")
  }
  {
    std::cout << "TEST 3 (y = A * x; NOT STRIDED):\n";
    ComplexVector y(3);
    ComplexVector x(std::vector<Complex>{Complex{1, 1}, Complex{2, 2}});
    ComplexMatrix A = build_test_complex_matrix(3, 2);
    mult_fast(y, A, x);
    std::cout << "A (in):\n"
              << A << '\n'
              << "x (in):\n"
              << x << '\n'
              << "y (out):\n"
              << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "(-5,15) (-35,105) (-65,195)",
                       "Bad mult")
  }
  {
    std::cout << "TEST 4 (y = A * x; STRIDED):\n";
    ComplexVector y(3);
    ComplexVector x(
        std::vector<Complex>{Complex{1, 1}, Complex{2, 2}, Complex{3, 3}});
    ComplexMatrix A(3, 3, Complex{1, 1});
    mult(y, A, x);
    std::cout << "A (in):\n"
              << A << '\n'
              << "x (in):\n"
              << x << '\n'
              << "y (out):\n"
              << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "(0,12) (0,12) (0,12)", "Bad mult")
  }
  {
    std::cout << "TEST 5 (y = A * x; STRIDED):\n";
    ComplexVector y(3);
    ComplexVector x(
        std::vector<Complex>{Complex{1, 1}, Complex{2, 2}, Complex{3, 3}});
    ComplexMatrix A = build_test_complex_matrix(3, 3);
    mult(y, A, x);
    std::cout << "A (in):\n"
              << A << '\n'
              << "x (in):\n"
              << x << '\n'
              << "y (out):\n"
              << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "(-14,42) (-74,222) (-134,402)",
                       "Bad mult")
  }
  {
    std::cout << "TEST 6 (y = A * x; STRIDED):\n";
    ComplexVector y(3);
    ComplexVector x(std::vector<Complex>{Complex{1, 1}, Complex{2, 2}});
    ComplexMatrix A = build_test_complex_matrix(3, 2);
    mult(y, A, x);
    std::cout << "A (in):\n"
              << A << '\n'
              << "x (in):\n"
              << x << '\n'
              << "y (out):\n"
              << y << '\n';
    ARTS_USER_ERROR_IF(var_string(y) not_eq "(-5,15) (-35,105) (-65,195)",
                       "Bad mult")
  }
  {
    std::cout << "TEST 7 (C = A * B; NON STRIDED):\n";
    ComplexMatrix C(3, 3);
    ComplexMatrix B(3, 3, Complex{1, 1});
    ComplexMatrix A(3, 3, Complex{1, 1});
    mult_fast(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(
        var_string(C) not_eq
            "(0,6) (0,6) (0,6)\n(0,6) (0,6) (0,6)\n(0,6) (0,6) (0,6)",
        "Bad mult")
  }
  {
    std::cout << "TEST 8 (C = A * B; NON STRIDED):\n";
    ComplexMatrix C(3, 4);
    ComplexMatrix A(3, 2, Complex{1, 1});
    ComplexMatrix B(2, 4, Complex{1, 1});
    mult_fast(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq
                           "(0,4) (0,4) (0,4) (0,4)\n(0,4) (0,4) (0,4) "
                           "(0,4)\n(0,4) (0,4) (0,4) (0,4)",
                       "Bad mult")
  }
  {
    std::cout << "TEST 9 (C = A * B; NON STRIDED):\n";
    ComplexMatrix C(3, 4);
    ComplexMatrix A = build_test_complex_matrix(3, 2);
    ComplexMatrix B = build_test_complex_matrix(2, 4);
    mult_fast(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq
                           "(-69,92) (-78,104) (-87,116) (-96,128)\n(-429,572) "
                           "(-498,664) (-567,756) (-636,848)\n(-789,1052) "
                           "(-918,1224) (-1047,1396) (-1176,1568)",
                       "Bad mult")
  }
  {
    std::cout << "TEST 10 (C = A * B; NON STRIDED):\n";
    ComplexMatrix C(4, 3);
    ComplexMatrix A = build_test_complex_matrix(4, 2);
    ComplexMatrix B = build_test_complex_matrix(2, 3);
    mult_fast(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(
        var_string(C) not_eq
            "(-69,92) (-78,104) (-87,116)\n(-429,572) (-498,664) "
            "(-567,756)\n(-789,1052) (-918,1224) (-1047,1396)\n(-1149,1532) "
            "(-1338,1784) (-1527,2036)",
        "Bad mult")
  }
  {
    std::cout << "TEST 11 (C = A * B; STRIDED):\n";
    ComplexMatrix C(3, 3);
    ComplexMatrix B(3, 3, Complex{1, 1});
    ComplexMatrix A(3, 3, Complex{1, 1});
    mult(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(
        var_string(C) not_eq
            "(0,6) (0,6) (0,6)\n(0,6) (0,6) (0,6)\n(0,6) (0,6) (0,6)",
        "Bad mult")
  }
  {
    std::cout << "TEST 12 (C = A * B; STRIDED):\n";
    ComplexMatrix C(3, 4);
    ComplexMatrix A(3, 2, Complex{1, 1});
    ComplexMatrix B(2, 4, Complex{1, 1});
    mult(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq
                           "(0,4) (0,4) (0,4) (0,4)\n(0,4) (0,4) (0,4) "
                           "(0,4)\n(0,4) (0,4) (0,4) (0,4)",
                       "Bad mult")
  }
  {
    std::cout << "TEST 13 (C = A * B; STRIDED):\n";
    ComplexMatrix C(3, 4);
    ComplexMatrix A = build_test_complex_matrix(3, 2);
    ComplexMatrix B = build_test_complex_matrix(2, 4);
    mult(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(var_string(C) not_eq
                           "(-69,92) (-78,104) (-87,116) (-96,128)\n(-429,572) "
                           "(-498,664) (-567,756) (-636,848)\n(-789,1052) "
                           "(-918,1224) (-1047,1396) (-1176,1568)",
                       "Bad mult")
  }
  {
    std::cout << "TEST 14 (C = A * B; STRIDED):\n";
    ComplexMatrix C(4, 3);
    ComplexMatrix A = build_test_complex_matrix(4, 2);
    ComplexMatrix B = build_test_complex_matrix(2, 3);
    mult(C, A, B);
    std::cout << "A (in):\n"
              << A << '\n'
              << "B (in):\n"
              << B << '\n'
              << "C (out):\n"
              << C << '\n';
    ARTS_USER_ERROR_IF(
        var_string(C) not_eq
            "(-69,92) (-78,104) (-87,116)\n(-429,572) (-498,664) "
            "(-567,756)\n(-789,1052) (-918,1224) (-1047,1396)\n(-1149,1532) "
            "(-1338,1784) (-1527,2036)",
        "Bad mult")
  }
  std::cout << "#/COMPLEX MULT TEST ######################################\n";
}

void test_transpose() {
  std::cout << "# TRANSPOSE TEST ######################################\n";
  {
    std::cout << "TEST 1:\n";
    auto A = build_test_matrix(3, 3);
    std::cout << A << '\n';
    ConstMatrixView X = A;
    std::cout << X.transpose() << '\n';
    std::cout << A.transpose() << '\n';
  }
  {
    std::cout << "TEST 2:\n";
    Matrix A = build_test_matrix(2, 3);
    std::cout << A << '\n';
    ConstMatrixView X = A;
    std::cout << X.transpose() << '\n';
    std::cout << A.transpose() << '\n';
  }
  {
    std::cout << "TEST 3:\n";
    auto A = build_test_matrix(4, 3);
    std::cout << A << '\n';
    ConstMatrixView X = A;
    std::cout << X.transpose() << '\n';
    std::cout << X.transpose().transpose() << '\n';
  }
  {
    std::cout << "TEST 3:\n";
    auto A = build_test_matrix(4, 3);
    ConstMatrixView X = A;
    std::cout << X.transpose() << '\n';
    std::cout << transpose(X) << '\n';
    std::cout << FastConstMatrixView{X.transpose().transpose()} << '\n';
  }
  std::cout << "#/TRANSPOSE TEST ######################################\n";
}

void test_complex_view() {
  std::cout << "# COMPLEX VIEWS TEST ######################################\n";
  {
    ComplexMatrix A = build_test_complex_matrix(3, 2);
    std::cout << A << '\n';
    std::cout << FastComplexMatrixView{A}.real() << '\n';
    std::cout << FastConstComplexMatrixView{A}.real() << '\n';
    std::cout << FastComplexMatrixView{A}.imag() << '\n';
    std::cout << FastConstComplexMatrixView{A}.imag() << '\n';
    std::cout << '\n';
    std::cout << A.transpose() << '\n';
    std::cout << FastComplexMatrixView{A}.real().transpose() << '\n';
    std::cout << FastConstComplexMatrixView{A}.real().transpose() << '\n';
    std::cout << FastComplexMatrixView{A}.imag().transpose() << '\n';
    std::cout << FastConstComplexMatrixView{A}.imag().transpose() << '\n';
  }
  {
    ComplexMatrix A = build_test_complex_matrix(3, 2);
    std::cout << A << '\n';
    std::cout << ComplexMatrixView{A}.real() << '\n';
    std::cout << ConstComplexMatrixView{A}.real() << '\n';
    std::cout << ComplexMatrixView{A}.imag() << '\n';
    std::cout << ConstComplexMatrixView{A}.imag() << '\n';
    std::cout << '\n';
    std::cout << A.transpose() << '\n';
    std::cout << ComplexMatrixView{A}.real().transpose() << '\n';
    std::cout << ConstComplexMatrixView{A}.real().transpose() << '\n';
    std::cout << ComplexMatrixView{A}.imag().transpose() << '\n';
    std::cout << ConstComplexMatrixView{A}.imag().transpose() << '\n';
    std::cout << '\n';
    std::cout << A << '\n';
    std::cout << ComplexMatrixView{A.transpose()}.real().transpose() << '\n';
    std::cout << ConstComplexMatrixView{A.transpose()}.real().transpose()
              << '\n';
    std::cout << ComplexMatrixView{A.transpose()}.imag().transpose() << '\n';
    std::cout << ConstComplexMatrixView{A.transpose()}.imag().transpose()
              << '\n';
    std::cout << '\n';
    std::cout << A.transpose() << '\n';
    std::cout << ComplexMatrixView{A.transpose()}.real() << '\n';
    std::cout << ConstComplexMatrixView{A.transpose()}.real() << '\n';
    std::cout << ComplexMatrixView{A.transpose()}.imag() << '\n';
    std::cout << ConstComplexMatrixView{A.transpose()}.imag() << '\n';
  }
  std::cout << "#/COMPLEX VIEWS TEST ######################################\n";
}

void test_range_view() {
  Vector VEC{std::vector<double>{0, 1, 2, 3, 4, 5}};
  std::cout << "VEC\n";
  std::cout << VEC << '\n';
  std::cout << VEC[Range(0, 2, 2)] << '\n';
  std::cout << VEC[Range(0, 3, 2)] << '\n';
  std::cout << VEC[Range(1, 2, 2)] << '\n';
  std::cout << VEC[Range(3, 2, 2)] << '\n';
  std::cout << VEC[Range(1, joker)] << '\n';
  std::cout << VEC[Range(1, joker, 2)] << '\n';

  Matrix MAT = build_test_matrix(6, 5);
  std::cout << "MAT\n";
  std::cout << MAT << '\n';
  std::cout << MAT[Range(0, 2, 2)] << '\n';
  std::cout << MAT[Range(0, 3, 2)] << '\n';
  std::cout << MAT[Range(1, 2, 2)] << '\n';
  std::cout << MAT[Range(3, 2, 2)] << '\n';

  std::cout << MAT(0, Range(0, 2, 2)) << '\n';
  std::cout << MAT(1, Range(0, 3, 2)) << '\n';
  std::cout << MAT(2, Range(1, 2, 2)) << '\n';
  std::cout << MAT(3, Range(3, 2, 1)) << '\n';

  std::cout << MAT(Range(0, 2, 2), Range(0, 2, 2)) << '\n';
  std::cout << MAT(Range(0, 2, 2), Range(0, 3, 2)) << '\n';
  std::cout << MAT(Range(0, 3, 2), Range(0, 2, 2)) << '\n';

  std::cout << MAT(Range(0, 2, 2), 1) << '\n';

  std::cout << MAT << '\n';
  std::cout << MAT[Range(joker)] << '\n';
  std::cout << MAT[joker] << '\n';
  std::cout << MAT[Range(1, 2)] << '\n';
  std::cout << MAT(Range(1, 2), joker) << '\n';
}

void test_impl() {
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

  z = std::move(yc).reshape(4, 3, 2);
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
  std::cout << "TEST GOOD (IF 6-LINES ABOVE LOOKS GOOD): Could get strided and "
               "non-strided data\n";
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
}

void test4() {
  Vector a(10);
  Vector b(a.nelem());

  for (Index i = 0; i < a.nelem(); ++i) {
    a[i] = (Numeric)(i + 1);
    b[i] = (Numeric)(a.nelem() - i);
  }

  cout << "a = \n" << a << "\n";
  cout << "b = \n" << b << "\n";
  cout << "a*b \n= " << a * b << "\n";

  Matrix A(11, 6);
  Matrix B(10, 20);
  Matrix C(20, 5);

  B = 2;
  C = 3;
  mult(A(Range(1, joker), Range(1, joker)), B, C);

  //  cout << "\nB =\n" << B << "\n";
  //  cout << "\nC =\n" << C << "\n";
  cout << "\nB*C =\n" << A << "\n";
}

void test5() {
  Vector a(10);
  Vector b(20);
  Matrix M(10, 20);

  // Fill b and M with a constant number:
  b = 1;
  M = 2;

  cout << "b = \n" << b << "\n";
  cout << "M =\n" << M << "\n";

  mult(a, M, b);  // a = M*b
  cout << "\na = M*b = \n" << a << "\n";

  mult(transpose((MatrixView)b), transpose((MatrixView)a), M);  // b^t = a^t * M
  cout << "\nb^t = a^t * M = \n" << transpose((MatrixView)b) << "\n";
}

void test6() {
  Index n = 5000;
  Vector x(1, n, 1), y(n);
  Matrix M(n, n);
  M = 1;
  //  cout << "x = \n" << x << "\n";

  cout << "Transforming.\n";
  //  transform(x,sin,x);
  // transform(transpose(y),sin,transpose(x));
  //  cout << "sin(x) =\n" << y << "\n";
  for (Index i = 0; i < 1000; ++i) {
    //      mult(y,M,x);
    auto tmp = static_cast<MatrixView>(y);  // FIXME: I have to have a temporary!
    transform(tmp, [](auto a){return sin(a);}, static_cast<MatrixView>(x));
    x += 1;
  }
  //  cout << "y =\n" << y << "\n";

  cout << "Done.\n";
}

void test7() {
  Vector x(1, 20000000, 1);
  Vector y(x.nelem());
  transform(y, [](auto a){return sin(a);}, x);
  cout << "min(sin(x)), max(sin(x)) = " << min(y) << ", " << max(y) << "\n";
}

void test8() {
  Vector x(80000000);
  for (Index i = 0; i < x.nelem(); ++i) x[i] = (Numeric)i;
  cout << "Done."
       << "\n";
}

void test9() {
  // Initialization of Matrix with view of other Matrix:
  Matrix A(4, 8);
  Matrix B(A(Range(joker), Range(0, 3)));
  cout << "B = " << B << "\n";
}

void test10() {
  // Initialization of Matrix with a vector (giving a 1 column Matrix).

  Vector v(1, 8, 1);
  Matrix M((MatrixView)v);  //FIXME: You have to view the vector as (MatrixView)
  cout << "M = " << M << "\n";
}

void test11() {
  // Assignment between Vector and Matrix:

  // At the moment doing this with a non-const Vector will result in a
  // warning message.
  Vector v(1, 8, 1);
  Matrix M(v.nelem(), 1);
  M = MatrixView{v};  //FIXME: You have to view the vector as MatrixView
  cout << "M = " << M << "\n";
}

void test13() {
  // Mix vector and one-column matrix in += operator.
  const Vector v(1, 8, 1);  // The const is necessary here to
                            // avoid compiler warnings about
                            // different conversion paths.
  Matrix M(ConstMatrixView{v});
  M += ConstMatrixView{v};  // FIXME: Should we allow implicit conversions to larger types???
  cout << "M = \n" << M << "\n";
}

void test17() {
  // Test Sum.
  Vector a(1, 10, 1);
  cout << "a.sum() = " << a.sum() << "\n";
}

void test18() {
  // Test elementvise square of a vector:
  Vector a(1, 10, 1);
  a *= a;
  cout << "a *= a =\n" << a << "\n";
}

void test19() {
  // There exists no explicit filling constructor of the form
  // Vector a(3,1.7).
  // But you can use the more general filling constructor with 3 arguments.

  Vector a(1, 10, 1);
  Vector b(5.3, 10, 0);
  cout << "a =\n" << a << "\n";
  cout << "b =\n" << b << "\n";
}

void test20() {
  // Test initialization list constructor:
  Vector a{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  cout << "a =\n" << a << "\n";
}

void test23() {
  // Test constructors that fills with constant:
  Vector a(10, 3.5);
  cout << "a =\n" << a << "\n";
  Matrix b(10, 10, 4.5);
  cout << "b =\n" << b << "\n";
}

void test24() {
  // Try element-vise multiplication of Matrix and Vector:
  Matrix a(5, 1, 2.5);
  Vector b(1, 5, 1);
  a *= MatrixView{b};  // FIXME: Must view as MatrixView
  cout << "a*=b =\n" << a << "\n";
  a /= MatrixView{b};  // FIXME: Must view as MatrixView
  cout << "a/=b =\n" << a << "\n";
  a += MatrixView{b};  // FIXME: Must view as MatrixView
  cout << "a+=b =\n" << a << "\n";
  a -= MatrixView{b};  // FIXME: Must view as MatrixView
  cout << "a-=b =\n" << a << "\n";
}

void test28() {
  cout << "Test default constructor for Matrix:\n";
  Matrix a;
  Matrix b(a);
  cout << "b =\n" << b << "\n";
}

void test30() {
  cout << "Test Matrices of size 0:\n";
  Matrix a(0, 0);
  //  cout << "a(0,0) =\n" << a(0,0) << "\n";
  a.resize(2, 2);
  a = 1;
  cout << "a =\n" << a << "\n";

  Matrix b(3, 0);
  //  cout << "b(0,0) =\n" << b(0,0) << "\n";
  b.resize(b.nrows(), b.ncols() + 3);
  b = 2;
  cout << "b =\n" << b << "\n";

  Matrix c(0, 3);
  //  cout << "c(0,0) =\n" << c(0,0) << "\n";
  c.resize(c.nrows() + 3, c.ncols());
  c = 3;
  cout << "c =\n" << c << "\n";
}

void test31() {
  cout << "Test Tensor3:\n";

  Tensor3 a(2, 3, 4, 1.0);

  Index fill = 0;

  // Fill with some numbers
  for (Index i = 0; i < a.npages(); ++i)
    for (Index j = 0; j < a.nrows(); ++j)
      for (Index k = 0; k < a.ncols(); ++k) a(i, j, k) = (Numeric)(++fill);

  cout << "a =\n" << a << "\n";

  cout << "Taking out first row of first page:\n"
       << a(0, 0, Range(joker)) << "\n";

  cout << "Taking out last column of second page:\n"
       << a(1, Range(joker), a.ncols() - 1) << "\n";

  cout << "Taking out the first letter on every page:\n"
       << a(Range(joker), 0, 0) << "\n";

  cout << "Taking out first page:\n"
       << a(0, Range(joker), Range(joker)) << "\n";

  cout << "Taking out last row of all pages:\n"
       << a(Range(joker), a.nrows() - 1, Range(joker)) << "\n";

  cout << "Taking out second column of all pages:\n"
       << a(Range(joker), Range(joker), 1) << "\n";

  a *= 2;

  cout << "After element-vise multiplication with 2:\n" << a << "\n";

  transform(a, [](auto x){return sqrt(x);}, a);  // FIXME: Cannot call just sqrt, have to lambda-wrap-it

  cout << "After taking the square-root:\n" << a << "\n";

  Index s = 200;
  cout << "Let's allocate a large tensor, "
       << (Numeric)(s * s * s * 8) / 1024. / 1024. << " MB...\n";

  a.resize(s, s, s);

  cout << "Set it to 1...\n";

  a = 1;

  cout << "a(900,900,900) = " << a(90, 90, 90) << "\n";

  fill = 0;

  cout << "Fill with running numbers, using for loops...\n";
  for (Index i = 0; i < a.npages(); ++i)
    for (Index j = 0; j < a.nrows(); ++j)
      for (Index k = 0; k < a.ncols(); ++k) a(i, j, k) = (Numeric)(++fill);

  cout << "Max(a) = ...\n";

  cout << max(a) << "\n";
}

void test32() {
  cout << "Test von X = A*X:\n";
  Matrix X(3, 3), A(3, 3), B(3, 3);

  for (Index j = 0; j < A.nrows(); ++j)
    for (Index k = 0; k < A.ncols(); ++k) {
      X(j, k) = 1;
      A(j, k) = (Numeric)(j + k);
    }
  cout << "A:\n" << A << "\n";
  cout << "X:\n" << X << "\n";

  mult(B, A, X);
  cout << "B = A*X:\n" << B << "\n";
  mult(X, A, X);
  cout << "X = A*X:\n" << X << "\n";

  cout << "This is not the same, and should not be, because you\n"
       << "are not allowed to use the same matrix as input and output\n"
       << "for mult!\n";
}

std::string describe(ConstTensor7View x) {
  std::string out;
  out += "Tensor 7 [";
  for (auto s: x.shape()) out += var_string(s, ", ");
  return out + "]";
}
std::string describe(ConstTensor6View x) {
  std::string out;
  out += "Tensor 6 [";
  for (auto s: x.shape()) out += var_string(s, ", ");
  return out + "]";
}
std::string describe(ConstTensor5View x) {
  std::string out;
  out += "Tensor 5 [";
  for (auto s: x.shape()) out += var_string(s, ", ");
  return out + "]";
}
std::string describe(ConstTensor4View x) {
  std::string out;
  out += "Tensor 4 [";
  for (auto s: x.shape()) out += var_string(s, ", ");
  return out + "]";
}
std::string describe(ConstTensor3View x) {
  std::string out;
  out += "Tensor 3 [";
  for (auto s: x.shape()) out += var_string(s, ", ");
  return out + "]";
}
std::string describe(ConstMatrixView x) {
  std::string out;
  out += "Matrix [";
  for (auto s: x.shape()) out += var_string(s, ", ");
  return out + "]";
}
std::string describe(ConstVectorView x) {
  std::string out;
  out += "Vector [";
  for (auto s: x.shape()) out += var_string(s, ", ");
  return out + "]";
}
std::string describe(const Numeric& x) {
  return var_string("Scalar (", x, ")");
}

void test33() {
  cout << "Making things look bigger than they are...\n";

  {
    cout << "1. Make a scalar look like a vector:\n";
    Numeric a = 3.1415;  // Just any number here.
    VectorView av(a);
    cout << "a, viewed as a vector: " << av << "\n";
    cout << "Describe a: " << describe(a) << "\n";
    av[0] += 1;
    cout << "a, after the first element\n"
         << "of the vector has been increased by 1: " << a << "\n";
  }

  {
    cout
        << "\n2. Make a vector look like a matrix:\n"
        << "This is an exception, because the new dimension is added at the end.\n";
    Vector a{1, 2, 3, 4, 5};
    MatrixView am = MatrixView{a};  // FIXME: Has to explicitly cast
    cout << "a, viewed as a matrix:\n" << am << "\n";
    cout << "Trasnpose view:\n" << transpose(am) << "\n";
    cout << "Describe transpose(am): " << describe(transpose(am)) << "\n";

    cout << "\n3. Make a matrix look like a tensor3:\n";
    Tensor3View at3 = Tensor3View{am};  // FIXME: Has to explicitly cast
    cout << "at3 = \n" << at3 << "\n";
    cout << "Describe at3: " << describe(at3) << "\n";
    at3(2, 0, 0) += 1;  // FIXME: The order of dimensions have changed
    cout << "a after Increasing element at3(0,2,0) by 1: \n" << a << "\n\n";

    Tensor4View at4 = Tensor4View{at3};  // FIXME: Has to explicitly cast
    cout << "at4 = \n" << at4 << "\n";
    cout << "Describe at4: " << describe(at4) << "\n";

    Tensor5View at5 = Tensor5View{at4};  // FIXME: Has to explicitly cast
    cout << "at5 = \n" << at5 << "\n";
    cout << "Describe at5: " << describe(at5) << "\n";

    Tensor6View at6 = Tensor6View{at5};  // FIXME: Has to explicitly cast
    cout << "at6 = \n" << at6 << "\n";
    cout << "Describe at6: " << describe(at6) << "\n";

    Tensor7View at7 = Tensor7View{at6};  // FIXME: Has to explicitly cast
    cout << "at7 = \n" << at7 << "\n";
    cout << "Describe at7: " << describe(at7) << "\n";

    at7(2, 0, 0, 0, 0, 0, 0) -= 1;  // FIXME: 

    cout << "After subtracting 1 from at7(0,0,0,0,0,2,0)\n"
         << "a = " << a << "\n";

    cout << "\nAll in one go:\n";
    Numeric b = 3.1415;  // Just any number here.
    Tensor7View bt7 = Tensor7View(b);  // FIXME: Can perform the shorter cast
    cout << "bt7:\n" << bt7 << "\n";
    cout << "Describe bt7: " << describe(bt7) << "\n";
  }
}

void junk4(Tensor4View a) { cout << "Describe a: " << describe(a) << "\n"; }

void junk2(ConstVectorView a) { cout << "Describe a: " << describe(a) << "\n"; }

void test34() {
  cout << "Test, if dimension expansion works implicitly.\n";

  Tensor3 t3(2, 3, 4);
  junk4(Tensor4View{t3});  // FIXME: Has to explicitly cast

  Numeric x;
  junk2(ConstVectorView(x));
}

void test35() {
  cout << "Test the new copy semantics.\n";
  Vector a(1, 4, 1);
  Vector b;

  b = a;
  cout << "b = " << b << "\n";

  Vector aa(1, 5, 1);
  ConstVectorView c = aa;
  b = c;
  cout << "b = " << b << "\n";

  Vector aaa(1, 6, 1);
  VectorView d = aaa;
  b = d;
  cout << "b = " << b << "\n";
}

void test36() {
  cout << "Test using naked joker on Vector.\n";
  Vector a(1, 4, 1);
  VectorView b = a[joker];
  cout << "a = " << a << "\n";
  cout << "b = " << b << "\n";
}

void test37(const Index& i) {
  Vector v1(5e-15, 10, 0.42e-15 / 11);
  Vector v2 = v1;
  //  Index i = 10000000;

  v1 /= (Numeric)i;
  v2 /= (Numeric)i;
  cout.precision(12);
  //  cout.setf(ios_base::scientific,ios_base::floatfield);
  v1 *= v1;
  v2 *= v2;
  cout << v1 << endl;
  cout << v2 << endl;
}

void test38() {
  Vector v(5, 0.);
  Numeric* const a = v.data_handle();

  a[4] = 5.;

  cout << v << endl;
  cout << endl << "========================" << endl << endl;

  Matrix m(5, 5, 0.);
  Numeric* const b = m.data_handle();

  b[4] = 5.;

  cout << m << endl;
  cout << endl << "========================" << endl << endl;

  Tensor3 t3(5, 6, 7, 0.);
  Numeric* const c = t3.data_handle();

  c[6] = 5.;

  cout << t3 << endl;
}

void test39() {
  Vector v1(1, 5, 1), v2(5);

  // v2 = v1 * 2;  // FIXME: This no longer compiles, should it?
  // Unfortunately, this thing compiles, but at least it gives an
  // assertion failure at runtime. I think what happens is that it
  // implicitly creates a one element vector out of the "2", then
  // tries to do a scalar product with v1.

  cout << v2 << endl;
}

void test40() {
  Vector v(100);

  v = 5;

  cout << v << endl;
}

void test41() {
  const Vector v1(10, 1);

  ConstVectorView vv = v1[Range(0, 5)];

  cout << "Vector:     " << v1 << endl;
  cout << "VectorView: " << vv << endl;

  try {
    // vv = 2;  //FIXME: This no longer compiles, should it?
  } catch (const std::runtime_error& e) {
    std::cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  cout << "VectorView: " << vv << endl;
  cout << "Vector:     " << v1 << endl;
}

void test42() {
  cout << "test42\n";
  Vector x{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  cout << "x: " << x << endl;

  VectorView y = x[Range(2, 4, 2)];
  cout << "y: " << y << endl;

  ConstMatrixView z((MatrixView)y);
  cout << "z:\n" << z << endl;

  cout << "Every other z:\n" << z(Range(1, 2, 2), joker) << endl;

  // Try write access also:
  MatrixView zz(y);
  zz(Range(1, 2, 2), joker) = 0;
  cout << "zz:\n" << zz << endl;
  cout << "New x: " << x << endl;
}

void test44() {
#define docheck(fn, val, expect)                                           \
  cout << #fn << "(" << val << ") = " << fn(x) << " (expected " << #expect \
       << ")" << std::endl;

  Vector x{1, 2, 3};
  docheck(is_increasing, x, 1) docheck(is_decreasing, x, 0)
      docheck(is_sorted, x, 1)

          x = {3, 2, 1};
  docheck(is_increasing, x, 0) docheck(is_decreasing, x, 1)
      docheck(is_sorted, x, 0)

          x = {1, 2, 2};
  docheck(is_increasing, x, 0) docheck(is_decreasing, x, 0)
      docheck(is_sorted, x, 1)

          x = {2, 2, 1};
  docheck(is_increasing, x, 0) docheck(is_decreasing, x, 0)
      docheck(is_sorted, x, 0)

          x = {1, NAN, 2};
  docheck(is_increasing, x, 0) docheck(is_decreasing, x, 0)
      docheck(is_sorted, x, 0)

          x = {2, NAN, 1};
  docheck(is_increasing, x, 0) docheck(is_decreasing, x, 0)
      docheck(is_sorted, x, 0)

          x = {NAN, NAN, NAN};
  docheck(is_increasing, x, 0) docheck(is_decreasing, x, 0)
      docheck(is_sorted, x, 0)

#undef docheck
}

void test46() {
  Vector v(5, 0.);
  nlinspace(v, 1, 10, 10);
  VectorView v1 = v;
  auto compare_func = [](Numeric n) { return n != 0; };
  cout << std::any_of(v1.begin(), v1.end(), compare_func) << endl;
  v1 = 0.;
  cout << std::any_of(v1.begin(), v1.end(), compare_func) << endl;
}

void test47() {
  // Selecting empty matpack dimensions with Joker shouldn't fail
  Vector v;
  Matrix m;

  VectorView vv = v[joker];
  MatrixView mv = m(joker, joker);

  std::cout << "vv.nelem: " << vv.nelem() << std::endl;
  std::cout << "mv.nrows: " << mv.nrows() << " ";
  std::cout << "mv.ncols: " << mv.ncols() << std::endl;

  Matrix m2(0, 5);
  MatrixView mv2 = MatrixView{m2(joker, 3)};  // FIXME: Must explicitly cast
  std::cout << "mv2.nrows: " << mv2.nrows() << " ";
  std::cout << "mv2.ncols: " << mv2.ncols() << std::endl;
}

void test48() {
  // Test each simple reduction
  std::cout << Tensor7(2, 2, 2, 2, 2, 2, 1, 8).reduce_rank<0, 1, 2, 3, 4, 5>() << '\n';
  std::cout << Tensor7(2, 2, 2, 2, 2, 1, 1, 9).reduce_rank<0, 1, 2, 3, 4>() << '\n';
  std::cout << Tensor7(2, 2, 2, 2, 1, 1, 1, 10).reduce_rank<0, 1, 2, 3>() << '\n';
  std::cout << Tensor7(2, 2, 2, 1, 1, 1, 1, 11).reduce_rank<0, 1, 2>() << '\n';
  std::cout << Tensor7(2, 2, 1, 1, 1, 1, 1, 12).reduce_rank<0, 1>() << '\n';
  std::cout << Tensor7(2, 1, 1, 1, 1, 1, 1, 13).reduce_rank<0>() << '\n';
  std::cout << Tensor7(1, 1, 1, 1, 1, 1, 1, 14) << '\n';
  std::cout << Tensor6(2, 2, 2, 2, 2, 1, 15).reduce_rank<0, 1, 2, 3, 4>() << '\n';
  std::cout << Tensor6(2, 2, 2, 2, 1, 1, 16).reduce_rank<0, 1, 2, 3>() << '\n';
  std::cout << Tensor6(2, 2, 2, 1, 1, 1, 17).reduce_rank<0, 1, 2>() << '\n';
  std::cout << Tensor6(2, 2, 1, 1, 1, 1, 18).reduce_rank<0, 1>() << '\n';
  std::cout << Tensor6(2, 1, 1, 1, 1, 1, 19).reduce_rank<0>() << '\n';
  std::cout << Tensor6(1, 1, 1, 1, 1, 1, 20) << '\n';
  std::cout << Tensor5(2, 2, 2, 2, 1, 21).reduce_rank<0, 1, 2, 3>() << '\n';
  std::cout << Tensor5(2, 2, 2, 1, 1, 22).reduce_rank<0, 1, 2>() << '\n';
  std::cout << Tensor5(2, 2, 1, 1, 1, 23).reduce_rank<0, 1>() << '\n';
  std::cout << Tensor5(2, 1, 1, 1, 1, 24).reduce_rank<0>() << '\n';
  std::cout << Tensor5(1, 1, 1, 1, 1, 25) << '\n';
  std::cout << Tensor4(2, 2, 2, 1, 26).reduce_rank<0, 1, 2>() << '\n';
  std::cout << Tensor4(2, 2, 1, 1, 27).reduce_rank<0, 1>() << '\n';
  std::cout << Tensor4(2, 1, 1, 1, 28).reduce_rank<0>() << '\n';
  std::cout << Tensor4(1, 1, 1, 1, 29) << '\n';
  std::cout << Tensor3(2, 2, 1, 30).reduce_rank<0, 1>() << '\n';
  std::cout << Tensor3(2, 1, 1, 31).reduce_rank<0>() << '\n';
  std::cout << Tensor3(1, 1, 1, 32) << '\n';
  std::cout << Matrix(2, 1, 33).reduce_rank<0>() << '\n';
  std::cout << Matrix(1, 1, 34) << '\n';
  
  // Test that the reductions work along different axis for Vector
  std::cout << matpack::eigen::mat(Tensor7(2, 1, 1, 1, 1, 1, 1, 35).reduce_rank<0, 1>()).transpose() << '\n'
            << matpack::eigen::mat(Tensor7(2, 1, 1, 1, 1, 1, 1, 35).reduce_rank<0, 2>()).transpose() << '\n'
            << matpack::eigen::mat(Tensor7(2, 1, 1, 1, 1, 1, 1, 35).reduce_rank<0, 3>()).transpose() << '\n'
            << matpack::eigen::mat(Tensor7(2, 1, 1, 1, 1, 1, 1, 35).reduce_rank<0, 4>()).transpose() << '\n'
            << matpack::eigen::mat(Tensor7(2, 1, 1, 1, 1, 1, 1, 35).reduce_rank<0, 5>()).transpose() << '\n'
            << matpack::eigen::mat(Tensor7(2, 1, 1, 1, 1, 1, 1, 35).reduce_rank<0, 6>()).transpose() << '\n';
  
  // Test that the reductions work along different axis for Matrix
  std::cout << matpack::eigen::mat(Tensor7(2, 2, 1, 1, 1, 1, 1, 36).reduce_rank<0, 1>()) << '\n'
            << matpack::eigen::mat(Tensor7(2, 1, 2, 1, 1, 1, 1, 36).reduce_rank<0, 2>()) << '\n'
            << matpack::eigen::mat(Tensor7(2, 1, 1, 2, 1, 1, 1, 36).reduce_rank<0, 3>()) << '\n'
            << matpack::eigen::mat(Tensor7(2, 1, 1, 1, 2, 1, 1, 36).reduce_rank<0, 4>()) << '\n'
            << matpack::eigen::mat(Tensor7(2, 1, 1, 1, 1, 2, 1, 36).reduce_rank<0, 5>()) << '\n'
            << matpack::eigen::mat(Tensor7(2, 1, 1, 1, 1, 1, 2, 36).reduce_rank<0, 6>()) << '\n';

  // Test that the right elements are accessed
  {
    Index val=1;
    Tensor7 testvar(9, 2, 4, 3, 5, 7, 11);for (Index i=0; i<9; i++)
    for (Index j=0; j<2; j++) {
      for (Index k=0; k<4; k++) {
        for (Index l=0; l<3; l++) {
          for (Index m=0; m<5; m++) {
            for (Index n=0; n<7; n++) {
              for (Index o=0; o<11; o++) {
                testvar(i, j, k, l, m, n, o) = Numeric(val++);
              }
            }
          }
        }
      }
    }
    
    //  Vector test
    for (Index j=0; j<2; j++) {
      for (Index k=0; k<4; k++) {
        for (Index l=0; l<3; l++) {
          for (Index m=0; m<5; m++) {
            for (Index n=0; n<7; n++) {
              for (Index o=0; o<11; o++) {
                std::cout << Tensor7(testvar(joker, Range(j, 1), Range(k, 1), Range(l, 1), Range(m, 1), Range(n, 1), Range(o, 1))).reduce_rank<0>() << '\n';
                std::cout << testvar(joker, j, k, l, m, n, o) << '\n';
                std::cout << '\n';
              }
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index k=0; k<4; k++) {
        for (Index l=0; l<3; l++) {
          for (Index m=0; m<5; m++) {
            for (Index n=0; n<7; n++) {
              for (Index o=0; o<11; o++) {
                std::cout << Tensor7(testvar(Range(i, 1), joker, Range(k, 1), Range(l, 1), Range(m, 1), Range(n, 1), Range(o, 1))).reduce_rank<1>() << '\n';
                std::cout << testvar(i, joker, k, l, m, n, o) << '\n';
                std::cout << '\n';
              }
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index l=0; l<3; l++) {
          for (Index m=0; m<5; m++) {
            for (Index n=0; n<7; n++) {
              for (Index o=0; o<11; o++) {
                std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), joker, Range(l, 1), Range(m, 1), Range(n, 1), Range(o, 1))).reduce_rank<2>() << '\n';
                std::cout << testvar(i, j, joker, l, m, n, o) << '\n';
                std::cout << '\n';
              }
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index k=0; k<4; k++) {
          for (Index m=0; m<5; m++) {
            for (Index n=0; n<7; n++) {
              for (Index o=0; o<11; o++) {
                std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), Range(k, 1), joker, Range(m, 1), Range(n, 1), Range(o, 1))).reduce_rank<3>() << '\n';
                std::cout << testvar(i, j, k, joker, m, n, o) << '\n';
                std::cout << '\n';
              }
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index k=0; k<4; k++) {
          for (Index l=0; l<3; l++) {
            for (Index n=0; n<7; n++) {
              for (Index o=0; o<11; o++) {
                std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), Range(k, 1), Range(l, 1), joker, Range(n, 1), Range(o, 1))).reduce_rank<4>() << '\n';
                std::cout << testvar(i, j, k, l, joker, n, o) << '\n';
                std::cout << '\n';
              }
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index k=0; k<4; k++) {
          for (Index l=0; l<3; l++) {
            for (Index m=0; m<5; m++) {
              for (Index o=0; o<11; o++) {
                std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), Range(k, 1), Range(l, 1), Range(m, 1), joker, Range(o, 1))).reduce_rank<5>() << '\n';
                std::cout << testvar(i, j, k, l, m, joker, o) << '\n';
                std::cout << '\n';
              }
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index k=0; k<4; k++) {
          for (Index l=0; l<3; l++) {
            for (Index m=0; m<5; m++) {
              for (Index n=0; n<7; n++) {
                std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), Range(k, 1), Range(l, 1), Range(m, 1), Range(n, 1), joker)).reduce_rank<6>() << '\n';
                std::cout << testvar(i, j, k, l, m, n, joker) << '\n';
                std::cout << '\n';
              }
            }
          }
        }
      }
    }
    
    // Matrix test
    for (Index k=0; k<4; k++) {
      for (Index l=0; l<3; l++) {
        for (Index m=0; m<5; m++) {
          for (Index n=0; n<7; n++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(joker, joker, Range(k, 1), Range(l, 1), Range(m, 1), Range(n, 1), Range(o, 1))).reduce_rank<0, 1>() << '\n';
              std::cout << testvar(joker, joker, k, l, m, n, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index j=0; j<2; j++) {
      for (Index l=0; l<3; l++) {
        for (Index m=0; m<5; m++) {
          for (Index n=0; n<7; n++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(joker, Range(j, 1), joker, Range(l, 1), Range(m, 1), Range(n, 1), Range(o, 1))).reduce_rank<0, 2>() << '\n';
              std::cout << testvar(joker, j, joker, l, m, n, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index j=0; j<2; j++) {
      for (Index k=0; k<4; k++) {
        for (Index m=0; m<5; m++) {
          for (Index n=0; n<7; n++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(joker, Range(j, 1), Range(k, 1), joker, Range(m, 1), Range(n, 1), Range(o, 1))).reduce_rank<0, 3>() << '\n';
              std::cout << testvar(joker, j, k, joker, m, n, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index j=0; j<2; j++) {
      for (Index k=0; k<4; k++) {
        for (Index l=0; l<3; l++) {
          for (Index n=0; n<7; n++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(joker, Range(j, 1), Range(k, 1), Range(l, 1), joker, Range(n, 1), Range(o, 1))).reduce_rank<0, 4>() << '\n';
              std::cout << testvar(joker, j, k, l, joker, n, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index j=0; j<2; j++) {
      for (Index k=0; k<4; k++) {
        for (Index l=0; l<3; l++) {
          for (Index m=0; m<5; m++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(joker, Range(j, 1), Range(k, 1), Range(l, 1), Range(m, 1), joker, Range(o, 1))).reduce_rank<0, 5>() << '\n';
              std::cout << testvar(joker, j, k, l, m, joker, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index j=0; j<2; j++) {
      for (Index k=0; k<4; k++) {
        for (Index l=0; l<3; l++) {
          for (Index m=0; m<5; m++) {
            for (Index n=0; n<7; n++) {
              std::cout << Tensor7(testvar(joker, Range(j, 1), Range(k, 1), Range(l, 1), Range(m, 1), Range(n, 1), joker)).reduce_rank<0, 6>() << '\n';
              std::cout << testvar(joker, j, k, l, m, n, joker) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index l=0; l<3; l++) {
        for (Index m=0; m<5; m++) {
          for (Index n=0; n<7; n++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(Range(i, 1), joker, joker, Range(l, 1), Range(m, 1), Range(n, 1), Range(o, 1))).reduce_rank<1, 2>() << '\n';
              std::cout << testvar(i, joker, joker, l, m, n, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index k=0; k<4; k++) {
        for (Index m=0; m<5; m++) {
          for (Index n=0; n<7; n++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(Range(i, 1), joker, Range(k, 1), joker, Range(m, 1), Range(n, 1), Range(o, 1))).reduce_rank<1, 3>() << '\n';
              std::cout << testvar(i, joker, k, joker, m, n, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index k=0; k<4; k++) {
        for (Index l=0; l<3; l++) {
          for (Index n=0; n<7; n++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(Range(i, 1), joker, Range(k, 1), Range(l, 1), joker, Range(n, 1), Range(o, 1))).reduce_rank<1, 4>() << '\n';
              std::cout << testvar(i, joker, k, l, joker, n, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index k=0; k<4; k++) {
        for (Index l=0; l<3; l++) {
          for (Index m=0; m<5; m++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(Range(i, 1), joker, Range(k, 1), Range(l, 1), Range(m, 1), joker, Range(o, 1))).reduce_rank<1, 5>() << '\n';
              std::cout << testvar(i, joker, k, l, m, joker, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index k=0; k<4; k++) {
        for (Index l=0; l<3; l++) {
          for (Index m=0; m<5; m++) {
            for (Index n=0; n<7; n++) {
              std::cout << Tensor7(testvar(Range(i, 1), joker, Range(k, 1), Range(l, 1), Range(m, 1), Range(n, 1), joker)).reduce_rank<1, 6>() << '\n';
              std::cout << testvar(i, joker, k, l, m, n, joker) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index m=0; m<5; m++) {
          for (Index n=0; n<7; n++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), joker, joker, Range(m, 1), Range(n, 1), Range(o, 1))).reduce_rank<2, 3>() << '\n';
              std::cout << testvar(i, j, joker, joker, m, n, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index l=0; l<3; l++) {
          for (Index n=0; n<7; n++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), joker, Range(l, 1), joker, Range(n, 1), Range(o, 1))).reduce_rank<2, 4>() << '\n';
              std::cout << testvar(i, j, joker, l, joker, n, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index l=0; l<3; l++) {
          for (Index m=0; m<5; m++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), joker, Range(l, 1), Range(m, 1), joker, Range(o, 1))).reduce_rank<2, 5>() << '\n';
              std::cout << testvar(i, j, joker, l, m, joker, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index l=0; l<3; l++) {
          for (Index m=0; m<5; m++) {
            for (Index n=0; n<7; n++) {
              std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), joker, Range(l, 1), Range(m, 1), Range(n, 1), joker)).reduce_rank<2, 6>() << '\n';
              std::cout << testvar(i, j, joker, l, m, n, joker) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index k=0; k<4; k++) {
          for (Index n=0; n<7; n++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), Range(k, 1), joker, joker, Range(n, 1), Range(o, 1))).reduce_rank<3, 4>() << '\n';
              std::cout << testvar(i, j, k, joker, joker, n, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index k=0; k<4; k++) {
          for (Index m=0; m<5; m++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), Range(k, 1), joker, Range(m, 1), joker, Range(o, 1))).reduce_rank<3, 5>() << '\n';
              std::cout << testvar(i, j, k, joker, m, joker, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index k=0; k<4; k++) {
          for (Index m=0; m<5; m++) {
            for (Index n=0; n<7; n++) {
              std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), Range(k, 1), joker, Range(m, 1), Range(n, 1), joker)).reduce_rank<3, 6>() << '\n';
              std::cout << testvar(i, j, k, joker, m, n, joker) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index k=0; k<4; k++) {
          for (Index l=0; l<3; l++) {
            for (Index o=0; o<11; o++) {
              std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), Range(k, 1), Range(l, 1), joker, joker, Range(o, 1))).reduce_rank<4, 5>() << '\n';
              std::cout << testvar(i, j, k, l, joker, joker, o) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index k=0; k<4; k++) {
          for (Index l=0; l<3; l++) {
            for (Index n=0; n<7; n++) {
              std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), Range(k, 1), Range(l, 1), joker, Range(n, 1), joker)).reduce_rank<4, 6>() << '\n';
              std::cout << testvar(i, j, k, l, joker, n, joker) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
    for (Index i=0; i<9; i++) {
      for (Index j=0; j<2; j++) {
        for (Index k=0; k<4; k++) {
          for (Index l=0; l<3; l++) {
            for (Index m=0; m<5; m++) {
              std::cout << Tensor7(testvar(Range(i, 1), Range(j, 1), Range(k, 1), Range(l, 1), Range(m, 1), joker, joker)).reduce_rank<5, 6>() << '\n';
              std::cout << testvar(i, j, k, l, m, joker, joker) << '\n';
              std::cout << '\n';
            }
          }
        }
      }
    }
  }
}

void test_eigen_conv() {
  Matrix a = build_test_matrix(3, 5);
  std::cout << a << '\n';
  std::cout << matpack::eigen::mat(a) << '\n';

  Vector b({1, 2, 3, 4, 5});
  std::cout << b << '\n';
  std::cout << matpack::eigen::row_vec(b) << '\n';
  std::cout << matpack::eigen::col_vec(b) << '\n';

  {
    FastMatrixView c = a;
    std::cout << c << '\n';
    std::cout << matpack::eigen::mat(c) << '\n';

    FastConstMatrixView d = c;
    std::cout << d << '\n';
    std::cout << matpack::eigen::mat(d) << '\n';

    FastVectorView e = b;
    std::cout << e << '\n';
    std::cout << matpack::eigen::row_vec(e) << '\n';
    std::cout << matpack::eigen::col_vec(e) << '\n';

    FastConstVectorView f = e;
    std::cout << f << '\n';
    std::cout << matpack::eigen::row_vec(f) << '\n';
    std::cout << matpack::eigen::col_vec(f) << '\n';
  }
  {
    Vector g(3);
    matpack::eigen::as_eigen(g).noalias() =
        matpack::eigen::as_eigen(a) * matpack::eigen::as_eigen(b);
    std::cout << g << '\n';
    std::cout << matpack::eigen::as_eigen(g) << '\n';
  }

  {
    MatrixView c = a;
    std::cout << c << '\n';
    std::cout << matpack::eigen::mat(c) << '\n';

    ConstMatrixView d = c;
    std::cout << d << '\n';
    std::cout << matpack::eigen::mat(d) << '\n';

    VectorView e = b;
    std::cout << e << '\n';
    std::cout << matpack::eigen::row_vec(e) << '\n';
    std::cout << matpack::eigen::col_vec(e) << '\n';

    ConstVectorView f = e;
    std::cout << f << '\n';
    std::cout << matpack::eigen::row_vec(f) << '\n';
    std::cout << matpack::eigen::col_vec(f) << '\n';
  }
}

void test_eigen_complex_conv() {
  ComplexMatrix a = build_test_complex_matrix(3, 5);
  std::cout << a << '\n';
  std::cout << matpack::eigen::mat(a) << '\n';

  ComplexVector b({{1,2}, {3,4}, {5,6}, {7,8}, {9,10}});
  std::cout << b << '\n';
  std::cout << matpack::eigen::row_vec(b) << '\n';
  std::cout << matpack::eigen::col_vec(b) << '\n';

  {
    FastComplexMatrixView c = a;
    std::cout << c << '\n';
    std::cout << matpack::eigen::mat(c) << '\n';

    MatrixView d = a.real();
    std::cout << "REAL\n";
    std::cout << d << '\n';
    std::cout << matpack::eigen::mat(d) << '\n';

    MatrixView e = a.imag();
    std::cout << "IMAG\n";
    std::cout << e << '\n';
    std::cout << matpack::eigen::mat(e) << '\n';
  }
}

void test_eigen_base_set() {
  {
    Vector a = {1,2,3,4,5};
    std::cout << a << '\n';

    auto b = matpack::eigen::row_vec(a);
    std::cout << b << '\n';
    std::cout << Vector(b) << '\n';

    auto c = matpack::eigen::col_vec(a);
    std::cout << c << '\n';
    std::cout << Vector(c) << '\n';
  }

  {
    ComplexVector a = {{1,2},{2,3},{3,4},{4,5},{5,6}};
    std::cout << a << '\n';

    auto b = matpack::eigen::row_vec(a);
    std::cout << b << '\n';
    std::cout << ComplexVector(b) << '\n';

    auto c = matpack::eigen::col_vec(a);
    std::cout << c << '\n';
    std::cout << ComplexVector(c) << '\n';
  }

  {
    Matrix a = build_test_matrix(3, 4);
    std::cout << a << '\n';

    auto b = matpack::eigen::mat(a);
    std::cout << b << '\n';
    std::cout << Matrix(b) << '\n';

    auto c = matpack::eigen::mat(a.transpose());
    std::cout << c << '\n';
    std::cout << Matrix(c) << '\n';
  }

  {
    ComplexMatrix a = build_test_complex_matrix(3, 4);
    std::cout << a << '\n';

    auto b = matpack::eigen::mat(a);
    std::cout << b << '\n';
    std::cout << ComplexMatrix(b) << '\n';

    auto c = matpack::eigen::mat(a.transpose());
    std::cout << c << '\n';
    std::cout << ComplexMatrix(c) << '\n';

    auto d = matpack::eigen::mat(a.real());
    std::cout << d << '\n';
    std::cout << Matrix(d) << '\n';

    auto e = matpack::eigen::mat(a.imag());
    std::cout << e << '\n';
    std::cout << Matrix(e) << '\n';

    auto f = matpack::eigen::mat(a.transpose().real());
    std::cout << f << '\n';
    std::cout << Matrix(f) << '\n';

    auto g = matpack::eigen::mat(a.transpose().imag());
    std::cout << g << '\n';
    std::cout << Matrix(g) << '\n';

    auto h = matpack::eigen::mat(a.real().transpose());
    std::cout << h << '\n';
    std::cout << Matrix(h) << '\n';

    auto i = matpack::eigen::mat(a.imag().transpose());
    std::cout << i << '\n';
    std::cout << Matrix(i) << '\n';
  }
}

void test_eigen_base_equal() {
  {
    Vector a = {1,2,3,4,5}, b, c;
    b = a;
    auto d = matpack::eigen::row_vec(a);
    c = d;
    std::cout << a << '\n';
    std::cout << b << '\n';
    std::cout << c << '\n';
    std::cout << d << '\n';
  }

  {
    Vector a = {1,2,3,4,5}, b, c;
    b = a;
    auto d = matpack::eigen::col_vec(a);
    c = d;
    std::cout << a << '\n';
    std::cout << b << '\n';
    std::cout << c << '\n';
    std::cout << d << '\n';
  }

  {
    Matrix a = build_test_matrix(3, 4), b, c;
    b = a;
    auto d = matpack::eigen::mat(a);
    c = d;
    std::cout << a << '\n';
    std::cout << b << '\n';
    std::cout << c << '\n';
    std::cout << d << '\n';
  }

  {
    ComplexMatrix a = build_test_complex_matrix(3, 4), b, c;
    b = a;
    auto d = matpack::eigen::mat(a);
    c = d;
    std::cout << a << '\n';
    std::cout << b << '\n';
    std::cout << c << '\n';
    std::cout << d << '\n';
  }
}

int main() {
  // test_impl();
  // test_mult();
  // test_complex_mult();
  // test_transpose();
  // test_complex_view();
  // test_range_view();
  // test1();
  // test2();
  // test4();
  // test5();
  // test6();
  // test7();
  // test8();
  // test9();
  // test10();
  // test11();
  // test13();
  // test17();
  // test18();
  // test19();
  // test20();
  // test23();
  // test24();
  // test28();
  // test30();
  // test31();
  // test32();
  // test33();
  // test34();
  // test35();
  // test36();
  // test37(5);
  // test38();
  // test39();
  // test40();
  // test41();
  // test42();
  // test44();
  // test46();
  // test47();
  // test48();
  // test_eigen_conv();
  // test_eigen_complex_conv();
  // test_eigen_base_set();
  test_eigen_base_equal();
}
