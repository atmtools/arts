#pragma once

#include "matpack_mdspan_data_t.h"

namespace matpack {
template <typename T> class basic_band_matrix {
  Index            KU;
  Index            KL;
  Index            M;
  Index            N;
  data_t<T, 2>     AB;
  std::vector<int> ipiv;

  static data_t<T, 2> mat(Index KL, Index KU, Index N) { return data_t<T, 2>(N, 2 * KL + KU + 1); }

 public:
  // Empty matrix of known size
  basic_band_matrix(Index ku, Index kl, Index m, Index n) : KU(ku), KL(kl), M(m), N(n), AB(mat(KL, KU, N)), ipiv(N) {}

  basic_band_matrix()                                        = default;
  basic_band_matrix(const basic_band_matrix&)                = default;
  basic_band_matrix(basic_band_matrix&&) noexcept            = default;
  basic_band_matrix& operator=(const basic_band_matrix&)     = default;
  basic_band_matrix& operator=(basic_band_matrix&&) noexcept = default;

  [[nodiscard]] constexpr Index end_row(Index j) const { return std::min<Index>(M, j + KL + 1); }

  [[nodiscard]] constexpr Index start_row(Index j) const { return std::max<Index>(0, j - KU); }

  explicit basic_band_matrix(const data_t<T, 2>& ab) : KU(0), KL(0), M(ab.nrows()), N(ab.ncols()), AB(0, 0), ipiv(N) {
    for (Index i = 0; i < M; ++i) {
      for (Index j = 0; j < N; ++j) {
        if (i == j) continue;
        if (ab[i, j] != T{}) {
          if (j > i)
            KU = std::max(j - i, KU);
          else
            KL = std::max(i - j, KL);
        }
      }
    }

    AB = mat(KL, KU, N);
    for (Index j = 0; j < N; ++j)
      for (Index i = start_row(j); i < end_row(j); ++i) operator[](i, j) = ab[i, j];
  }

  template <std::integral A, std::integral B> [[nodiscard]] constexpr T& operator[](const A& i, const B& j) {
    assert(static_cast<Index>(i) >= start_row(static_cast<Index>(j)));
    assert(static_cast<Index>(i) < end_row(static_cast<Index>(j)));
    return AB[j, KU + KL + i - j];
  }

  template <typename U> void zero(this U&& self) { std::forward<U>(self).AB = T{}; }

  //! Solves the system of equations A * x = b destructively
  int solve(data_t<T, 1>& bx);
};

using band_matrix         = basic_band_matrix<Numeric>;
using complex_band_matrix = basic_band_matrix<Complex>;
}  // namespace matpack
