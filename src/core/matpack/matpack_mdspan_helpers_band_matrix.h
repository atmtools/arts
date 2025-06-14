#pragma once

#include "matpack_mdspan_data_t.h"

namespace matpack {
class band_matrix {
  Index KU;
  Index KL;
  Index M;
  Index N;
  Matrix AB;
  std::vector<int> ipiv;

  static Matrix mat(Index KL, Index KU, Index N);

 public:
  // Empty matrix of known size
  band_matrix(Index ku, Index kl, Index m, Index n);

  band_matrix();
  band_matrix(const band_matrix&);
  band_matrix(band_matrix&&) noexcept;
  band_matrix& operator=(const band_matrix&);
  band_matrix& operator=(band_matrix&&) noexcept;

  [[nodiscard]] Index end_row(Index j) const;

  [[nodiscard]] Index start_row(Index j) const;

  band_matrix(const Matrix& ab);

  Numeric& operator[](Index i, Index j);

  //! Solves the system of equations A * x = b destructively
  int solve(Vector& bx);
};
}  // namespace matpack