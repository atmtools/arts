#include "artstime.h"
#include "matpack_arrays.h"
#include "matpack_math.h"
#include "matpack_view.h"

#include <cstdlib>

void test_sum(Index N) {
  Numeric X;
  Vector a(N, 1);
  {
    DebugTime timer{"Summing a Vector"};
    X = sum(a);
  }
  const ConstVectorView b = a;
  {
    DebugTime timer{"Summing a const ConstVectorView"};
    X = sum(b);
  }
}

int main(int argc, char** c) {
  Index N = 100'000;
  if (argc == 2) N = static_cast<Index>(std::atoll(c[1]));

  test_sum(N);
}