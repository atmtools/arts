/*!
  \file   test_cia.cc
  \author  <sbuehler@ltu.se>
  \date   2012-11-30
  
  \brief  Test Collission Induced Absorption (CIA) functions.
*/

#include "cia.h"
#include "matpack_math.h"

void test01() {
  std::cout << "Testing CIA Interpolation.\n";
  GriddedField2 cia_data;

  Matrix A(5, 3, 0.);
  A(2, 1) = 1;
  //    cout << "A:" << A << std::endl;

  cia_data.data = A;
  cia_data.set_grid(0, {1, 2, 3, 4, 5});
  cia_data.set_grid(1, {100, 200, 300});

  std::cout << "cia_data:" << cia_data << std::endl;

  // Output frequencies and temperature:
  Vector f_out=uniform_grid(1, 9, 0.5);
  std::cout << "f_out:" << f_out << std::endl;
  Numeric T_out = 150;
  std::cout << "T_out:" << T_out << std::endl;

  Vector result(9);
  cia_interpolation(result, f_out, T_out, cia_data, 0.5, 0);
  std::cout << "result:" << result << std::endl;
}

int main() {
  test01();
  return 0;
}
