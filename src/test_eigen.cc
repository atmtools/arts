// Define this so that the Eigen::internal::set_is_malloc_allowed testing command is available
#define EIGEN_RUNTIME_NO_MALLOC

#include "matpack_eigen.h"

using namespace matpack::eigen;

void test_col_and_row_vec() {
  Vector a({1,2,3,4,5});
  ComplexVector b({Complex{1, 2},Complex{3, 4},Complex{5,6},Complex{7,8},Complex{9,10}});

  std::cout << col_vec(a) << '\n';
  std::cout << col_vec(b) << '\n';
  std::cout << col_vec(b.real()) << '\n';
  std::cout << col_vec(b).real() << '\n';
  std::cout << col_vec(b.imag()) << '\n';
  std::cout << col_vec(b).imag() << '\n';

  std::cout << row_vec(a) << '\n';
  std::cout << row_vec(b) << '\n';
  std::cout << row_vec(b.real()) << '\n';
  std::cout << row_vec(b).real() << '\n';
  std::cout << row_vec(b.imag()) << '\n';
  std::cout << row_vec(b).imag() << '\n';
}

void test_multiplication_and_setting() {
  Matrix A(3, 4, 1);
  Vector x(4, 2);

  std::cout << A << "\n*\n" << x << "\n=\n" << A * x << '\n';
  std::cout << Vector(A * x) << '\n';
  Vector y;
  y = A * x;
  std::cout << y << '\n' << '\n';

  std::cout << A << "\n*\n" << transpose(A) << "\n=\n" <<  A * transpose(A) << '\n';
  std::cout << Matrix(A * transpose(A)) << '\n';
  Matrix XTX;
  XTX = A * transpose(A);
  std::cout << XTX << '\n' << '\n';

  // Neither of these allocate!
  Eigen::internal::set_is_malloc_allowed(false);
  row_vec(y).noalias() = A * x;
  mat(XTX).noalias() = A * transpose(A);
  Eigen::internal::set_is_malloc_allowed(true);
}

void test_chaining_expressions() {
  Matrix A(3, 4, 1);
  Vector x({1, 2, 3, 4}), b({5, 6, 7}), y(3);

  std::cout << A << "\n*\n" << x << "\n+\n" << b << "\n=\n" << A * x + b << '\n';
  row_vec(y).noalias() = A * x + b;
  std::cout << y << '\n' << '\n';

  std::cout << b << "\n+\n" << A << "\n*\n" << x << "\n=\n" << b + A * x << '\n';
  row_vec(y).noalias() = b + A * x;
  std::cout << y << '\n' << '\n';

  std::cout << A << "\n*\n" << x << "\n-\n" << b << "\n=\n" << A * x - b << '\n';
  row_vec(y).noalias() = A * x - b;
  std::cout << y << '\n' << '\n';

  std::cout << b << "\n-\n" << A << "\n*\n" << x << "\n=\n" << b - A * x << '\n';
  row_vec(y).noalias() = b - A * x;
  std::cout << y << '\n' << '\n';

  std::cout << 5*y << '\n';
  std::cout << -1.2*y << '\n';
  std::cout << 5*A << '\n';
  std::cout << -1.2*A << '\n';
}

int main() {
  test_col_and_row_vec();
  test_multiplication_and_setting();
  test_chaining_expressions();
}
