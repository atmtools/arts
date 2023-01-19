#include "artstime.h"
#include "matpack_arrays.h"
#include "matpack_math.h"

#include <cstdlib>

void test_sum(Index N) {
  Numeric X;
  {
    DebugTime timer{"Vector"};
    Vector a(N, 1);
    X = sum(a);
  }
  std::cout << "Got " << X << " expected something close to " << static_cast<Numeric>(N) << '\n';
}

int main(int argc, char** c) {
  Index N = 100'000;
  if (argc == 2) N = static_cast<Index>(std::atoll(c[1]));

  test_sum(N);
}