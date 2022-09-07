#include "matpack_eigen.h"

int main() {
  using namespace matpack::eigen;
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