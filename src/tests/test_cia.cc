/*!
  \file   test_cia.cc
  \author  <sbuehler@ltu.se>
  \date   2012-11-30
  
  \brief  Test Collission Induced Absorption (CIA) functions.
*/

#include <iostream>

#include "cia.h"
#include "matpack.h"

void test01() {
  std::cout << "Testing CIA Interpolation.\n";
  GriddedField2 cia_data;

  Matrix A(5, 3, 0.);
  A[2, 1] = 1;
  //    cout << "A:" << A << '\n';

  cia_data.data      = A;
  cia_data.grid<0>() = {1, 2, 3, 4, 5};
  cia_data.grid<1>() = {100, 200, 300};

  arts::println("cia_data:{}", cia_data);

  // Output frequencies and temperature:
  Vector f_out = matpack::uniform_grid(1, 9, 0.5);
  arts::println("f_out:{}", f_out);
  Numeric T_out = 150;
  arts::println("T_out:{}", T_out);

  Vector result(9);
  cia_interpolation(result, f_out, T_out, cia_data, 0.5, 0);
  arts::println("result:{}", result);
}

int main() {
  test01();
  return 0;
}
