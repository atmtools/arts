#include "artstime.h"
#include "matpack_arrays.h"
#include "matpack_math.h"
#include "matpack_view.h"

#include <cstdlib>

void test_sum(Index N) {
  {
    Numeric X;
    Vector a(N, 1);
    const ConstVectorView b = a;
    VectorView c = a;

    {
      DebugTime timer{"STARTUP"};
      X = sum(a) + sum(a) + sum(a) + sum(a) + sum(a);
    }
    std::cout << X << '\n';

    {
      DebugTime timer{"Summing a Vector                         "};
      X = sum(a);
    }
    std::cout << X << '\n';

    {
      DebugTime timer{"Summing a const ConstVectorView          "};
      X = sum(b);
    }
    std::cout << X << '\n';

    {
      DebugTime timer{"Summing a VectorView                     "};
      X = sum(c);
    }
    std::cout << X << '\n';

    {
      DebugTime timer{"STARTUP"};
      X = mean(a) + mean(a) + mean(a) + mean(a) + mean(a);
    }
    std::cout << X << '\n';

    {
      DebugTime timer{"Mean a Vector                         "};
      X = mean(a);
    }
    std::cout << X << '\n';

    {
      DebugTime timer{"Mean a const ConstVectorView          "};
      X = mean(b);
    }
    std::cout << X << '\n';

    {
      DebugTime timer{"Mean a VectorView                     "};
      X = mean(c);
    }
    std::cout << X << '\n';
  }
}

int main(int argc, char** c) {
  Index N = 100'000;
  if (argc == 2) N = static_cast<Index>(std::atoll(c[1]));

  test_sum(N);
}