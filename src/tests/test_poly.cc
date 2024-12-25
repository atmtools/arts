#include <poly_roots.h>

#include <iostream>

int main(void) {
  Vector v(9, 0);
  Matrix s(8, 2);

  v[0] = 2;
  v[4] = 1;
  v[8] = 8;

  int status = poly_root_solve(s, v);

  std::cout << status << '\n';
  std::cout << std::format("{}", s) << '\n';

  return (0);
}
