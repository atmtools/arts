#include "debug.h"
#include "matpack_mdspan.h"

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
  test10();
}
