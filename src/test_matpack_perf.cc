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

void test_mat_multiply(Index N) {
  Matrix A(20, 20), C(N, 20, 1), B(20, N, 1);
  MatrixView Av=A, Bv=B, Cv=C;

  Numeric X;
  {
    DebugTime timer{"STARTUP"};
    mult(A, B, C);
    mult(A, B, C);
    mult(A, B, C);
    mult(A, B, C);
    X = A.elem_at(0);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"Matrix mult    "};
    mult(A, B, C);
    X = A(0, 0);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"MatrixView mult"};
    mult(Av, Bv, Cv);
    X = Av(0, 0);
  }
  std::cout << X << '\n';
}

void test_elementary_ops(Index N) {
  Vector A(N, 1);
  VectorView Av = A;

  Numeric X;
  {
    DebugTime timer{"STARTUP"};
    A = 1;
    A += 1;
    A *= 2;
    A /= 2;
    A -= 1;
    X = A[0] + A[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"Vector Set    "};
    A = 1;
    X = A[0] + A[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"VectorView Set"};
    Av = 1;
    X = Av[0] + Av[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"Vector PlusEq    "};
    A += 1.33;
    X = A[0] + A[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"VectorView PlusEq"};
    Av += 1.33;
    X = Av[0] + Av[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"Vector MinusEq    "};
    A -= 1.5;
    X = A[0] + A[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"VectorView MinusEq"};
    Av -= 1.5;
    X = Av[0] + Av[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"Vector TimesEq    "};
    A *= 3.5;
    X = A[0] + A[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"VectorView TimesEq"};
    Av *= 3.5;
    X = Av[0] + Av[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"Vector DivEq    "};
    A /= 7.77;
    X = A[0] + A[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"VectorView DivEq"};
    A /= 7.77;
    X = Av[0] + Av[N-1];
  }
  std::cout << X << '\n';
}

int main(int argc, char** c) {
  Index N = 100'000;
  if (argc == 2) N = static_cast<Index>(std::atoll(c[1]));

  //test_sum(N);
  //test_dot(N);
  //test_vec_mult(N);
  //test_mat_multiply(N);
  test_elementary_ops(N);
}