#include <matpack.h>

#include <iostream>
#include <limits>

#include "debug.h"

int main() {
  const Matrix ex = [] {
    Matrix out(4, 4);
    out[0] = std::array{-0.23, 2.54, -3.66, 0.0};
    out[1] = std::array{-6.98, 2.46, -2.73, -2.13};
    out[2] = std::array{0.0, 2.56, 2.46, 4.07};
    out[3] = std::array{0.0, 0.0, -4.78, -3.82};
    return out;
  }();

  matpack::band_matrix bd(ex);
  const Vector         b{4.42, 27.13, -6.14, 10.50};
  Vector               sparse_b{b};
  Vector               dense_y{b};

  //! Solves inline
  bd.solve(sparse_b);

  //! Solves out-of-place
  solve(dense_y, ex, b);

  //! Diff should be just 0s
  dense_y -= sparse_b;

  //! Ensure that the difference is within the machine epsilon
  for (auto& x : dense_y) {
    ARTS_USER_ERROR_IF(std::abs(x) > 1000 * std::numeric_limits<Numeric>::epsilon(),
                       "Error in band matrix solver!\nOutput supposed to be: {}"
                       "\nBut diff between dense and banded matrix solutions are: {}",
                       sparse_b,
                       dense_y)
  }

  const ComplexMatrix complex_a = [] {
    ComplexMatrix out(3, 3, 0.0);
    out[0, 0] = Complex{3.0, 1.0};
    out[0, 1] = Complex{-1.0, 0.5};
    out[1, 0] = Complex{2.0, -1.0};
    out[1, 1] = Complex{4.0, 0.0};
    out[1, 2] = Complex{0.5, 2.0};
    out[2, 1] = Complex{-2.0, 0.25};
    out[2, 2] = Complex{1.0, -3.0};
    return out;
  }();
  const ComplexVector          complex_b{Complex{1.0, 2.0}, Complex{-3.0, 0.5}, Complex{2.0, -1.0}};
  ComplexVector                complex_x{complex_b};
  matpack::complex_band_matrix complex_bd(complex_a);
  ARTS_USER_ERROR_IF(complex_bd.solve(complex_x) != 0, "Complex band matrix solver failed");

  for (Index i = 0; i < complex_a.nrows(); ++i) {
    Complex residual = -complex_b[i];
    for (Index j = 0; j < complex_a.ncols(); ++j) residual += complex_a[i, j] * complex_x[j];
    ARTS_USER_ERROR_IF(std::abs(residual) > 1000 * std::numeric_limits<Numeric>::epsilon(),
                       "Complex band matrix solver residual is {}",
                       residual);
  }

  return 0;
}
