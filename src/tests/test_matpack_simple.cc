#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "matpack_data.h"

// Simple element access operator benchmark
void test1() {
  Matrix m(500, 500);
  for (Index k = 0; k < 1000; k++)
    for (Index i = 0; i < m.nrows(); i++)
      for (Index j = 0; j < m.ncols(); j++) m(i, j) = 2.;
}

int main() {
  test1();

  return 0;
}
