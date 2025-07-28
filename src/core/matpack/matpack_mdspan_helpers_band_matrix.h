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

  [[nodiscard]] constexpr Index end_row(Index j) const {
    return std::min<Index>(M, j + KL + 1);
  }

  [[nodiscard]] constexpr Index start_row(Index j) const {
    return std::max<Index>(0, j - KU);
  }

  band_matrix(const Matrix& ab);

  template <std::integral A, std::integral B>
  [[nodiscard]] constexpr Numeric& operator[](const A& i, const B& j) {
    assert(static_cast<Index>(i) >= start_row(static_cast<Index>(j)));
    assert(static_cast<Index>(i) < end_row(static_cast<Index>(j)));
    return AB[j, KU + KL + i - j];
  }

  //! Solves the system of equations A * x = b destructively
  int solve(Vector& bx);
};
}  // namespace matpack
