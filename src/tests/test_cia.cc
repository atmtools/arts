/*!
  \file   test_cia.cc
  \author  <sbuehler@ltu.se>
  \date   2012-11-30
  
  \brief  Test Collission Induced Absorption (CIA) functions.
*/

#include "cia.h"
#include "matpack.h"
#include <iostream>

void test01() {
  std::cout << "Testing CIA Interpolation.\n";
  GriddedField2 cia_data;

  Matrix A(5, 3, 0.);
  A[2, 1] = 1;
  //    cout << "A:" << A << '\n';

  cia_data.data = A;
  cia_data.grid<0>() = {1, 2, 3, 4, 5};
  cia_data.grid<1>() = {100, 200, 300};

  std::cout << "cia_data:" << std::format("{}", cia_data) << '\n';

  // Output frequencies and temperature:
  Vector f_out=matpack::uniform_grid(1, 9, 0.5);
  std::cout << "f_out:" << std::format("{}", f_out) << '\n';
  Numeric T_out = 150;
  std::cout << "T_out:" << std::format("{}", T_out) << '\n';

  Vector result(9);
  cia_interpolation(result, f_out, T_out, cia_data, 0.5, 0);
  std::cout << "result:" << std::format("{}", result) << '\n';
}

int main() {
  test01();
  return 0;
}
