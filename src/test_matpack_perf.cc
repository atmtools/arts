#include "artstime.h"
#include "matpack_arrays.h"
#include "matpack_eigen.h"
#include "matpack_math.h"
#include "matpack_view.h"

#include <cstdlib>

void test_sum(Index N) {
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

  {
    DebugTime timer{"STARTUP"};
    X = nanmean(a) + nanmean(a) + nanmean(a) + nanmean(a) + nanmean(a);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"NAN-Mean a Vector                         "};
    X = nanmean(a);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"NAN-Mean a const ConstVectorView          "};
    X = nanmean(b);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"NAN-Mean a VectorView                     "};
    X = nanmean(c);
  }
  std::cout << X << '\n';
}

void test_dot(Index N) {
  Vector a(N, 1);
  Vector b(N, 1);

  VectorView av=a;
  VectorView bv=b;

  Numeric X;
  {
    DebugTime timer{"STARTUP"};
    X = a * b + a * b + a * b + a * b + a * b + a * b;
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"Vector dot product    "};
    X = a * b;
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"VectorView dot product"};
    X = av * bv;
  }
  std::cout << X << '\n';
}

void test_vec_mult(Index N) {
  Vector x(N, 1);
  Matrix A(N, N, 1);
  Vector y(N, 1);

  VectorView xv=x;
  MatrixView Av=A;
  VectorView yv=y;

  Numeric X;
  {
    DebugTime timer{"STARTUP"};
    mult(x, A, y);
    X = x[0];
    mult(x, A, y);
    X += x[0];
    mult(x, A, y);
    X += x[0];
    mult(x, A, y);
    X += x[0];
    mult(x, A, y);
    X += x[0];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"Vector vector mult product    "};
    mult(x, A, y);
    X = sum(x);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"VectorView vector mult product"};
    mult(xv, Av, yv);
    X = sum(xv);
  }
  std::cout << X << '\n';
}

int main(int argc, char** c) {
  Index N = 100'000;
  if (argc == 2) N = static_cast<Index>(std::atoll(c[1]));

  //test_sum(N);
  //test_dot(N);
  test_vec_mult(N);
}