#ifndef template_partfun_h
#define template_partfun_h

#include <algorithm>
#include <array>

#include "enums.h"
#include "matpackI.h"

namespace PartitionFunctions {
ENUMCLASS(Type, Index, Interp, Coeff)

struct Data {
  Type type;
  Matrix data;
  
  void print_data() const;
  
  void print_method() const;
  
  friend std::ostream& operator<<(std::ostream& os, const Data& d) {
    return os << d.data << '\n';
  }
};

enum class Derivatives : bool {No, Yes};

template <Derivatives deriv, std::size_t N> 
Numeric linterp(const std::array<Numeric, N>& Tg,
                const std::array<Numeric, N>& Qg,
                const Numeric T) noexcept {
  static_assert(N > 1);
  
  // First position
  const std::size_t i = std::distance(Tg.cbegin(),
    std::min(std::lower_bound(Tg.cbegin(), Tg.cend(), T), Tg.cend() - 2));
  
  if constexpr (Derivatives::No == deriv) {
    return Qg[i] + (T - Tg[i]) * (Qg[i + 1] - Qg[i]) / (Tg[i + 1] - Tg[i]);
  } else {
    return (Qg[i + 1] - Qg[i]) / (Tg[i + 1] - Tg[i]);
  }
}

template <Derivatives deriv, std::size_t N> constexpr
Numeric polynom(const std::array<Numeric, N>& coeffs, const Numeric T) noexcept {
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
}

using PartitionFunctionsType = PartitionFunctions::Type;
using PartitionFunctionsData = PartitionFunctions::Data;

#endif  // template_partfun_h
