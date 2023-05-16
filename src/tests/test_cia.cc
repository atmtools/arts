/*!
  \file   test_cia.cc
  \author  <sbuehler@ltu.se>
  \date   2012-11-30
  
  \brief  Test Collission Induced Absorption (CIA) functions.
*/

#include "cia.h"
#include "matpack_math.h"

void test01() {
  cout << "Testing CIA Interpolation.\n";
  GriddedField2 cia_data;

  Matrix A(5, 3, 0.);
  A(2, 1) = 1;
  //    cout << "A:" << A << endl;

  cia_data.data = A;
  cia_data.set_grid(0, {1, 2, 3, 4, 5});
  cia_data.set_grid(1, {100, 200, 300});

  cout << "cia_data:" << cia_data << endl;

  // Output frequencies and temperature:
  Vector f_out=uniform_grid(1, 9, 0.5);
  cout << "f_out:" << f_out << endl;
  Numeric T_out = 150;
  cout << "T_out:" << T_out << endl;

  Vector result(9);
  cia_interpolation(result, f_out, T_out, cia_data, 0.5, 0, Verbosity(3, 3, 0));
  cout << "result:" << result << endl;
}

int main() {
  test01();
  return 0;
}
