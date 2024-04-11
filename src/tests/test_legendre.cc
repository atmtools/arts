#include "legendre.h"

#include <iostream>
#include <ranges>

Numeric factorial(Index i) {
  Index result = 1;
  for (Index j = 1; j <= i; j++) {
    result *= j;
  }
  return static_cast<Numeric>(result);

}

int main() {
  for (int l = 0; l < 3; l++) {
    for (int m =-l; m <= l; m++) {
      std::cout << "l = " << l << ", m = " << m << ": " << Legendre::assoc_legendre(l, m, 0.5) << '\n';
    }
  }

  std::cout << std::ranges::iota_view(1, 5)[2] << '\n';

  for (Index x=1; x<5; x++){
    for (Index y = 1; y < 5; y++) {
      std::cout << Legendre::tgamma_ratio(static_cast<Numeric>(x) + 1,
                                          static_cast<Numeric>(y) + 1)
                << " vs " << factorial(x) / factorial(y) << '\n';
    }
  }
}
