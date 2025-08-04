#include <matpack.h>

#include <iostream>

#include "interp.h"

int main() {
  const GriddedField1 f{.data_name  = "TestData-1D",
                        .data       = Vector{1, 2, 3, 4, 5, 6},
                        .grid_names = {"x"},
                        .grids      = {Vector{2, 4, 6, 8, 10, 12}}};
  std::print(std::cout, "Matpack data:\n{:B,Ns}\n", f);

  std::print(std::cout,
             "Matpack data:\n{:B,Ns}\n",
             f.reinterp<LagrangeInterpolation>({2, 3}, 1));
  std::print(std::cout,
             "Matpack data:\n{:B,Ns}\n",
             f.reinterp<FixedLagrangeInterpolation<1>>({2, 3, 4, 1}));
  std::print(
      std::cout, "Numeric: {}\n", f.interp<LagrangeInterpolation>(20, 1));
  std::print(
      std::cout, "Numeric: {}\n", f.interp<FixedLagrangeInterpolation<1>>(-1));

  const GriddedField2 g{.data_name  = "TestData-2D",
                        .data       = Vector{1, 2, 3, 4, 5, 6}.reshape(2, 3),
                        .grid_names = {"x", "y"},
                        .grids      = {Vector{1, 2}, Vector{1, 2, 3}}};
  std::print(std::cout, "Matpack data:\n{:B,Ns}\n", g);

  std::print(std::cout,
             "Matpack data:\n{:B,Ns}\n",
             g.reinterp<LagrangeInterpolation, LagrangeInterpolation>(
                 {1.5}, {1.5}, 1, 20));
  std::print(
      std::cout,
      "Matpack data:\n{:B,Ns}\n",
      g.reinterp<FixedLagrangeInterpolation<1>, FixedLagrangeInterpolation<1>>(
          {2, 3, 4, 1}, {-5, 3}, 20));
  std::print(std::cout,
             "Numeric: {}\n",
             g.interp<LagrangeInterpolation, LagrangeInterpolation>(0, 0, 1));
  std::print(
      std::cout,
      "Numeric: {}\n",
      g.interp<FixedLagrangeInterpolation<1>, FixedLagrangeInterpolation<1>>(
          -1, -2));
}
