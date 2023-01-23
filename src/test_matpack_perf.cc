#include "artstime.h"
#include "matpack_arrays.h"
#include "matpack_eigen.h"
#include "matpack_math.h"
#include "matpack_view.h"

#include <cstdlib>
#include <stdexcept>

void test_sum(Index N) {
  Numeric X;
  Vector a(N, 1);
  auto b = static_cast<ExhaustiveVectorView>(a);
  VectorView c = a;

  {
    DebugTime timer{"STARTUP"};
    X = sum(a) + sum(a) + sum(a) + sum(a) + sum(a);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"sum(vec)"};
    X = sum(a);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"sum(exview)"};
    X = sum(b);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"sum(view)"};
    X = sum(c);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"STARTUP"};
    X = mean(a) + mean(a) + mean(a) + mean(a) + mean(a);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mean(vec)"};
    X = mean(a);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mean(exview)"};
    X = mean(b);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mean(view)"};
    X = mean(c);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"STARTUP"};
    X = nanmean(a) + nanmean(a) + nanmean(a) + nanmean(a) + nanmean(a);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"nanmean(vec)"};
    X = nanmean(a);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"nanmean(exview)"};
    X = nanmean(b);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"nanmean(view)"};
    X = nanmean(c);
  }
  std::cout << X << '\n';
}

void test_dot(Index N) {
  Vector a(N, 1);
  Vector b(N, 1);

  auto ae=static_cast<ExhaustiveVectorView>(a);
  auto be=static_cast<ExhaustiveVectorView>(b);

  VectorView av=a;
  VectorView bv=b;

  Numeric X;
  {
    DebugTime timer{"STARTUP"};
    X = a * b + a * b + a * b + a * b + a * b + a * b;
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"vec * vec"};
    X = a * b;
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exview * exview"};
    X = ae * be;
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"view * view"};
    X = av * bv;
  }
  std::cout << X << '\n';
}

void test_vec_mult(Index N) {
  Vector x(N, 1);
  Matrix A(N, N, 1);
  Vector y(N, 1);

  auto xe=static_cast<ExhaustiveVectorView>(x);
  MatrixView Ae=static_cast<ExhaustiveMatrixView>(A);
  auto ye=static_cast<ExhaustiveVectorView>(y);

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
    DebugTime timer{"vec = mat * vec"};
    mult(x, A, y);
    X = sum(x);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exview = exmview * exview"};
    mult(xe, Ae, ye);
    X = sum(x);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"view = mview * view"};
    mult(xv, Av, yv);
    X = sum(xv);
  }
  std::cout << X << '\n';
}

void test_mat_multiply(Index N) {
  Matrix A(20, 20), C(N, 20, 1), B(20, N, 1);
  auto Ae=static_cast<ExhaustiveMatrixView>(A), Be=static_cast<ExhaustiveMatrixView>(B), Ce=static_cast<ExhaustiveMatrixView>(C);
  MatrixView Av=A, Bv=B, Cv=C;

  Numeric X;
  {
    DebugTime timer{"STARTUP"};
    mult(A, B, C);
    mult(A, B, C);
    mult(A, B, C);
    mult(A, B, C);
    X = A.elem_at(0, 0);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mat = mat * mat"};
    mult(A, B, C);
    X = A(0, 0);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exmview = exmview * exmview"};
    mult(Ae, Be, Ce);
    X = Ae(0, 0);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mview = mview * mview"};
    mult(Av, Bv, Cv);
    X = Av(0, 0);
  }
  std::cout << X << '\n';
}

void test_elementary_ops_vec(Index N) {
  Vector A(N, 1);
  auto Ae = static_cast<ExhaustiveVectorView>(A);
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
    DebugTime timer{"vec = 1"};
    A = 1;
    X = A[0] + A[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exview = 1"};
    Ae = 1;
    X = Ae[0] + Ae[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"view = 1"};
    Av = 1;
    X = Av[0] + Av[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"vec += 1.33"};
    A += 1.33;
    X = A[0] + A[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exview += 1.33"};
    Ae += 1.33;
    X = Ae[0] + Ae[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"view += 1.33"};
    Av += 1.33;
    X = Av[0] + Av[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"vec -= 1.5"};
    A -= 1.5;
    X = A[0] + A[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exview -= 1.5"};
    Ae -= 1.5;
    X = Ae[0] + Ae[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"view -= 1.5"};
    Av -= 1.5;
    X = Av[0] + Av[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"vec *= 3.5"};
    A *= 3.5;
    X = A[0] + A[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exview *= 3.5"};
    Ae *= 3.5;
    X = Ae[0] + Ae[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"view *= 3.5"};
    Av *= 3.5;
    X = Av[0] + Av[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"vec /= 7.77"};
    A /= 7.77;
    X = A[0] + A[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exview /= 7.77"};
    Ae /= 7.77;
    X = Ae[0] + Ae[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"view /= 7.77"};
    Av /= 7.77;
    X = Av[0] + Av[N-1];
  }
  std::cout << X << '\n';
}

void test_elementary_ops_mat(Index N) {
  Matrix A(N, N, 1);
  auto Ae = static_cast<ExhaustiveMatrixView>(A);
  MatrixView Av = A;

  Numeric X;
  {
    DebugTime timer{"STARTUP"};
    A = 1;
    A += 1;
    A *= 2;
    A /= 2;
    A -= 1;
    X = A(0, 0) + A(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mat = 1"};
    A = 1;
    X = A(0, 0) + A(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exmview = 1"};
    Ae = 1;
    X = Ae(0, 0) + Ae(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mview = 1"};
    Av = 1;
    X = Av(0, 0) + Av(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mat += 1.33"};
    A += 1.33;
    X = A(0, 0) + A(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exmview += 1.33"};
    Ae += 1.33;
    X = Ae(0, 0) + Ae(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mview += 1.33"};
    Av += 1.33;
    X = Av(0, 0) + Av(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mat -= 1.5"};
    A -= 1.5;
    X = A(0, 0) + A(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exmview -= 1.5"};
    Ae -= 1.5;
    X = Ae(0, 0) + Ae(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mview -= 1.5"};
    Av -= 1.5;
    X = Av(0, 0) + Av(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mat *= 3.5"};
    A *= 3.5;
    X = A(0, 0) + A(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exmview *= 3.5"};
    Ae *= 3.5;
    X = Ae(0, 0) + Ae(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mview *= 3.5"};
    Av *= 3.5;
    X = Av(0, 0) + Av(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mat /= 7.77"};
    A /= 7.77;
    X = A(0, 0) + A(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exmview /= 7.77"};
    Ae /= 7.77;
    X = Ae(0, 0) + Ae(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mview /= 7.77"};
    Av /= 7.77;
    X = Av(0, 0) + Av(N-1, N-1);
  }
  std::cout << X << '\n';
}

void test_ops_vec(Index N) {
  Vector a(N, 1.1);
  Vector b(N, 2.2);
  auto ae=static_cast<ExhaustiveVectorView>(a);
  auto be=static_cast<ExhaustiveVectorView>(b);
  VectorView av=a;
  VectorView bv=b;

  Numeric X;
  {
    DebugTime timer{"STARTUP"};
    a += b;
    a -= b;
    a *= b;
    a /= b;
    X = a[0] + a[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"vec += vec"};
    a += b;
    X = a[0] + a[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exview += exview"};
    ae += be;
    X = ae[0] + ae[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"view += view"};
    av += bv;
    X = av[0] + av[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"vec -= vec"};
    a -= b;
    X = a[0] + a[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exview -= exview"};
    ae -= be;
    X = ae[0] + ae[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"view -= view"};
    av -= bv;
    X = av[0] + av[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"vec *= vec"};
    a *= b;
    X = a[0] + a[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exview *= exview"};
    ae *= be;
    X = ae[0] + ae[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"view *= exview"};
    av *= bv;
    X = av[0] + av[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"vec /= vec"};
    a /= b;
    X = a[0] + a[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exview /= exview"};
    ae /= be;
    X = ae[0] + ae[N-1];
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"view /= view"};
    av /= bv;
    X = av[0] + av[N-1];
  }
  std::cout << X << '\n';
}

void test_ops_mat(Index N) {
  Matrix a(N, N, 1.1);
  Matrix b(N, N, 2.2);
  auto ae=static_cast<ExhaustiveMatrixView>(a);
  auto be=static_cast<ExhaustiveMatrixView>(b);
  MatrixView av=a;
  MatrixView bv=b;

  Numeric X;
  {
    DebugTime timer{"STARTUP"};
    a += b;
    a -= b;
    a *= b;
    a /= b;
    X = a(0, 0) + a(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mat += mat"};
    a += b;
    X = a(0, 0) + a(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exmview += exmview"};
    ae += be;
    X = ae(0, 0) + ae(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mview += mview"};
    av += bv;
    X = av(0, 0) + av(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mat -= mat"};
    a -= b;
    X = a(0, 0) + a(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exmview -= exmview"};
    ae -= be;
    X = ae(0, 0) + ae(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mview -= mview"};
    av -= bv;
    X = av(0, 0) + av(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mat *= mat"};
    a *= b;
    X = a(0, 0) + a(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exmview *= exmview"};
    ae *= be;
    X = ae(0, 0) + ae(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mview *= mview"};
    av *= bv;
    X = av(0, 0) + av(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mat /= mat"};
    a /= b;
    X = a(0, 0) + a(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"exmview /= exmview"};
    ae /= be;
    X = ae(0, 0) + ae(N-1, N-1);
  }
  std::cout << X << '\n';

  {
    DebugTime timer{"mview /= mview"};
    av /= bv;
    X = av(0, 0) + av(N-1, N-1);
  }
  std::cout << X << '\n';
}

int main(int argc, char** c) {
  std::array <Index, 8> N;
  if (static_cast<std::size_t>(argc) < 1 + 1 + N.size()) {
    std::cerr << "Expects PROGNAME NREPEAT NSIZE..., wehere NSIZE is " << N.size() << " indices\n";
    return EXIT_FAILURE;
  }

  const auto n = static_cast<Index>(std::atoll(c[1]));
  for (std::size_t i=0; i<N.size(); i++)  N[i] = static_cast<Index>(std::atoll(c[2 + i]));

  for (Index i=0; i<n; i++) {
    test_sum(N[0]);
    test_dot(N[1]);
    test_vec_mult(N[2]);
    test_mat_multiply(N[3]);
    test_elementary_ops_vec(N[4]);
    test_elementary_ops_mat(N[5]);
    test_ops_vec(N[6]);
    test_ops_mat(N[7]);
  }
}