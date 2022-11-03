#ifndef template_partfun_h
#define template_partfun_h

#include <algorithm>
#include <array>
#include <cmath>

#include "enums.h"
#include "matpack.h"
#include "matpackI.h"

namespace PartitionFunctions {
ENUMCLASS(Type, Index, Interp, Coeff, StaticInterp)

struct Data {
  Type type;
  Matrix data;
  
  friend std::ostream& operator<<(std::ostream& os, const Data& d) {
    return os << d.data << '\n';
  }
};

enum class Derivatives : bool {No, Yes};

template <Derivatives deriv, std::size_t N> 
constexpr Numeric linterp(const std::array<Numeric, N>& Tg,
                          const std::array<Numeric, N>& Qg,
                          const Numeric T) noexcept {
  static_assert(N > 1);

  // First position (note that we are only at left-most grid point when T<=Tg[0] because of the logical operand)
  const auto i_low = std::distance(Tg.cbegin(), std::lower_bound(Tg.cbegin(), Tg.cend(), T));
  const auto i = std::min<std::size_t>(i_low - (i_low > 0), N - 2);

  if constexpr (Derivatives::No == deriv) {
    return Qg[i] + (T - Tg[i]) * (Qg[i + 1] - Qg[i]) / (Tg[i + 1] - Tg[i]);
  } else {
    return (Qg[i + 1] - Qg[i]) / (Tg[i + 1] - Tg[i]);
  }
}

#if __cpp_nontype_template_args >= 201911L
#define STATIC_LINTERP(deriv, Q, T, dT, T0) static_linterp<deriv, dT, T0>(Q, T)
template <Derivatives deriv, Numeric dT, Numeric T0, std::size_t N>
constexpr Numeric static_linterp(const std::array<Numeric, N> &Q,
                                 const Numeric T) noexcept {
  constexpr auto r_dT = 1.0 / dT;
#else
#define STATIC_LINTERP(deriv, Q, T, dT, T0) static_linterp<deriv>(Q, T, dT, T0)
template <Derivatives deriv, std::size_t N>
constexpr Numeric static_linterp(const std::array<Numeric, N> &Q,
                                 const Numeric T, Numeric dT,
                                 Numeric T0) noexcept {
  const auto r_dT = 1.0 / dT;
#endif
  static_assert(N > 1);

  const auto Tx = (T - T0) * r_dT;
  const auto iTx = static_cast<std::size_t>(Tx);
  const auto i = std::min<std::size_t>(iTx, N - 2);

  if constexpr (Derivatives::No == deriv) {
    const auto To = Tx - static_cast<Numeric>(i);
    return Q[i] + To * (Q[i + 1] - Q[i]);
  } else {
    return (Q[i + 1] - Q[i]) * r_dT;
  }
}

template <Derivatives deriv, std::size_t N> constexpr
Numeric polynom(const std::array<Numeric, N>& coeffs, const Numeric T) {
  static_assert(N > 0);
  
  Numeric result = 0.;
  Numeric TN = 1;
  
  if constexpr (Derivatives::No == deriv) {
    for (auto& c: coeffs) {
      result += TN * c;
      TN *= T;
    }
  } else {
    for (std::size_t i = 1; i < N; i++) {
      result += Numeric(i) * TN * coeffs[i];
      TN *= T;
    }
  }
  
  return result;
}
} // namespace PartitionFunctions

using PartitionFunctionsType = PartitionFunctions::Type;
using PartitionFunctionsData = PartitionFunctions::Data;

#endif  // template_partfun_h
