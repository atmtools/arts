#include <rng.h>
#include <iostream>
#include <random>

int main() {
  auto x = random_numbers<std::uniform_int_distribution>(10, 0, 20);
  std::cout << std::format("{}", x) << '\n';
  auto y = random_numbers<2>({2, 3}, 0.0, 1.0);
  std::cout << std::format("{}", y) << '\n';
}
