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
   std::vector ms {0, 1, 2};
   std::vector ls{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
   Numeric x = .080442;
   for (auto m: ms) {
    for (auto l: ls | std::ranges::views::drop(m)) {
      std::cout << "P(" << l << "," << m << "): " << Legendre::assoc_legendre(l, m, -x) << std::endl;
    }
   }
}
