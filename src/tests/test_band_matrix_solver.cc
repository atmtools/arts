#include <matpack.h>

#include <iostream>
#include <limits>

#include "debug.h"
#include "matpack_band_matrix.h"

int main() {
  const Matrix ex = []() {
    Matrix out(4, 4);
    out[0] = {-0.23, 2.54, -3.66, 0.0};
    out[1] = {-6.98, 2.46, -2.73, -2.13};
    out[2] = {0.0, 2.56, 2.46, 4.07};
    out[3] = {0.0, 0.0, -4.78, -3.82};
    return out;
  }();

  matpack::band_matrix bd(ex);
  const Vector b{4.42, 27.13, -6.14, 10.50};
  Vector sparse_b{b};
  Vector dense_y{b};

  //! Solves inline
  bd.solve(sparse_b);

  //! Solves out-of-place
  solve(dense_y, ex, b);

  //! Diff should be just 0s
  dense_y -= sparse_b;

  //! Ensure that the difference is within the machine epsilon
  for (auto& x : dense_y) {
    ARTS_USER_ERROR_IF(std::abs(x) > 10*std::numeric_limits<Numeric>::epsilon(),
                       "Error in band matrix solver!\nOutput supposed to be: ",
                       sparse_b,
                       "\nBut diff between dense and banded matrix solutions are: ",
                       dense_y)
  }

  return 0;
}
