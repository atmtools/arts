#include <matpack.h>

#include "interp.h"

#include <iostream>

int main() {
  const GriddedField1 f{.data = Vector{1, 2, 3, 4, 5, 6},
                        .grids = Vector{2, 4, 6, 8, 10, 12}};
  std::cout << f << '\n';

  std::cout << f.reinterp<LagrangeInterpolation>({2, 3}, 1) << '\n';
  std::cout << f.reinterp<FixedLagrangeInterpolation<1>>({2, 3, 4, 1}) << '\n';
  std::cout << f.interp<LagrangeInterpolation>(20, 1) << '\n';
  std::cout << f.interp<FixedLagrangeInterpolation<1>>(-1) << '\n';

  const GriddedField2 g{.data = Vector{1, 2, 3, 4, 5, 6}.reshape(2, 3),
                        .grids = {Vector{1, 2}, Vector{1, 2, 3}}};
  std::cout << g << '\n';

  std::cout << g.reinterp<LagrangeInterpolation, LagrangeInterpolation>(
                   {1.5}, {1.5}, 1, 20)
            << '\n';
  std::cout << g.reinterp<FixedLagrangeInterpolation<1>,
                          FixedLagrangeInterpolation<1>>(
                   {2, 3, 4, 1}, {-5, 3}, 20)
            << '\n';
  std::cout << g.interp<LagrangeInterpolation, LagrangeInterpolation>(0, 0, 1)
            << '\n';
  std::cout
      << g.interp<FixedLagrangeInterpolation<1>, FixedLagrangeInterpolation<1>>(
             -1, -2)
      << '\n';
}
