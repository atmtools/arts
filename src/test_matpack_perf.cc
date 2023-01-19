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
    const ExhaustiveConstVectorView d = a;
    ExhaustiveVectorView e = a;

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
      DebugTime timer{"Summing a const ExhaustiveConstVectorView"};
      X = sum(d);
    }
    std::cout << X << '\n';

    {
      DebugTime timer{"Summing a ExhaustiveVectorView           "};
      X = sum(e);
    }
    std::cout << X << '\n';
  }
  {
    Numeric X;
    Matrix a(N/20, 20, 1);
    const ConstMatrixView b = a;
    MatrixView c = a;
    const ExhaustiveConstMatrixView d = a;
    ExhaustiveMatrixView e = a;

    {
      DebugTime timer{"Summing a Matrix                         "};
      X = sum(a);
    }
    std::cout << X << '\n';

    {
      DebugTime timer{"Summing a const ConstMatrixView          "};
      X = sum(b);
    }
    std::cout << X << '\n';

    {
      DebugTime timer{"Summing a MatrixView                     "};
      X = sum(c);
    }
    std::cout << X << '\n';

    {
      DebugTime timer{"Summing a const ExhaustiveConstMatrixView"};
      X = sum(d);
    }
    std::cout << X << '\n';

    {
      DebugTime timer{"Summing a ExhaustiveMatrixView           "};
      X = sum(e);
    }
    std::cout << X << '\n';
  }
}

int main(int argc, char** c) {
  Index N = 100'000;
  if (argc == 2) N = static_cast<Index>(std::atoll(c[1]));

  test_sum(N);
}