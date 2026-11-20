#include <rng.h>

#include <iostream>
#include <random>

int main() {
  auto x = random_numbers<1, std::uniform_int_distribution>(10, 0, 20);
  arts::println("{}", x);
  Matrix y = random_numbers<2>({2, 3}, 0.0, 1.0);
  arts::println("{}", y);
}
