#include "rng.h"
#include <random>

int main() {
  auto x = random_numbers<std::uniform_int_distribution>(10, 0, 20);
  std::cout << x << '\n';
  auto y = random_numbers<2>({2, 3}, 0.0, 1.0);
  std::cout << y << '\n';
}
