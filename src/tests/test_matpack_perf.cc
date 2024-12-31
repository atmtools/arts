#include <artstime.h>
#include <matpack.h>

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "test_perf.h"

Array<Timing> test_sum(Index N) {
  Numeric X;
  Vector a(N, 1);
  auto b = static_cast<VectorView>(a);
  VectorView c = a;

  Array<Numeric> results_;
  Array<Timing> out;

  { X = sum(a) + sum(a) + sum(a) + sum(a) + sum(a); };
  results_.push_back(X);

  out.emplace_back("sum(Vector)")([&]() { X = sum(a); });
  results_.push_back(X);

  out.emplace_back("sum(VectorView)")([&]() { X = sum(b); });
  results_.push_back(X);

  out.emplace_back("sum(VectorView)")([&]() { X = sum(c); });
  results_.push_back(X);

  out.emplace_back("mean(Vector)")([&]() { X = mean(a); });
  results_.push_back(X);

  out.emplace_back("mean(VectorView)")([&]() { X = mean(b); });
  results_.push_back(X);

  out.emplace_back("mean(VectorView)")([&]() { X = mean(c); });
  results_.push_back(X);

  out.emplace_back("nanmean(Vector)")([&]() { X = nanmean(a); });
  results_.push_back(X);

  out.emplace_back("nanmean(VectorView)")([&]() { X = nanmean(b); });
  results_.push_back(X);

  out.emplace_back("nanmean(VectorView)")([&]() { X = nanmean(c); });
  results_.push_back(X);

  out.emplace_back("dummy")([results_]() { return results_[results_.size() - 1]; });

  return out;
}

Array<Timing> test_dot(Index N) {
  Vector a(N, 1);
  Vector b(N, 1);

  auto ae = static_cast<VectorView>(a);
  auto be = static_cast<VectorView>(b);

  VectorView av = a;
  VectorView bv = b;

  Array<Numeric> results_;
  Array<Timing> out;

  Numeric X;
  { X = dot(a , b) + dot(a, b) + dot(a, b) + dot(a, b) + dot(a, b) + dot(a, b); };
  results_.push_back(X);

  out.emplace_back("Vector*Vector")([&]() { X = dot(a, b); });
  results_.push_back(X);

  out.emplace_back("VectorView*VectorView")([&]() { X = dot(ae, be); });
  results_.push_back(X);

  out.emplace_back("VectorView*VectorView")([&]() { X = dot(av, bv); });
  results_.push_back(X);

  out.emplace_back("dummy")([results_]() { return results_[results_.size() - 1]; });

  return out;
}

Array<Timing> test_Vector_mult(Index N) {
  Vector x(N, 1);
  Matrix A(N, N, 1);
  Vector y(N, 1);

  auto xe = static_cast<VectorView>(x);
  MatrixView Ae = static_cast<MatrixView>(A);
  auto ye = static_cast<VectorView>(y);

  VectorView xv = x;
  MatrixView Av = A;
  VectorView yv = y;

  Array<Numeric> results_;
  Array<Timing> out;

  Numeric X;
  {
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
  results_.push_back(X);

  out.emplace_back("mult(Vector,Matrix,Vector)")([&]() {
    mult(x, A, y);
    X = x[0] + x[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("mult(VectorView,MatrixView,VectorView)")([&]() {
    mult(xe, Ae, ye);
    X = xe[0] + xe[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("mult(VectorView,MatrixView,VectorView)")([&]() {
    mult(xv, Av, yv);
    X = xv[0] + xv[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("dummy")([results_]() { return results_[results_.size() - 1]; });

  return out;
}

Array<Timing> test_Matrix_multiply(Index N) {
  Matrix A(N, N), C(N, N, 1), B(N, N, 1);
  auto Ae = static_cast<MatrixView>(A),
       Be = static_cast<MatrixView>(B),
       Ce = static_cast<MatrixView>(C);
  MatrixView Av = A, Bv = B, Cv = C;

  Array<Numeric> results_;
  Array<Timing> out;

  Numeric X;
  {
    mult(A, B, C);
    mult(A, B, C);
    mult(A, B, C);
    mult(A, B, C);
    X = A[0, 0];
  }
  results_.push_back(X);

  out.emplace_back("mult(Matrix,Matrix,Matrix)")([&]() {
    mult(A, B, C);
    X = A[0, 0] + A[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("mult(MatrixView,MatrixView,MatrixView)")([&]() {
    mult(Ae, Be, Ce);
    X = Ae[0, 0] + Ae[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("mult(MatrixView,MatrixView,MatrixView)")([&]() {
    mult(Av, Bv, Cv);
    X = Av[0, 0] + Ae[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("dummy")([results_]() { return results_[results_.size() - 1]; });

  return out;
}

Array<Timing> test_elementary_ops_Vector(Index N) {
  Vector A(N, 1);
  auto Ae = static_cast<VectorView>(A);
  VectorView Av = A;

  Array<Numeric> results_;
  Array<Timing> out;

  Numeric X;
  {
    A = 1;
    A += 1;
    A *= 2;
    A /= 2;
    A -= 1;
    X = A[0] + A[N - 1];
  }
  results_.push_back(X);

  out.emplace_back("Vector::operator=(Numeric)")([&]() {
    A = 1;
    X = A[0] + A[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator=(Numeric)")([&]() {
    Ae = 1;
    X = Ae[0] + Ae[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator=(Numeric)")([&]() {
    Av = 1;
    X = Av[0] + Av[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Vector::operator+=(Numeric)")([&]() {
    A += 1.33;
    X = A[0] + A[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator+=(Numeric)")([&]() {
    Ae += 1.33;
    X = Ae[0] + Ae[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator+=(Numeric)")([&]() {
    Av += 1.33;
    X = Av[0] + Av[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Vector::operator-=(Numeric)")([&]() {
    A -= 1.5;
    X = A[0] + A[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator-=(Numeric)")([&]() {
    Ae -= 1.5;
    X = Ae[0] + Ae[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator-=(Numeric)")([&]() {
    Av -= 1.5;
    X = Av[0] + Av[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Vector::operator*=(Numeric)")([&]() {
    A *= 3.5;
    X = A[0] + A[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator*=(Numeric)")([&]() {
    Ae *= 3.5;
    X = Ae[0] + Ae[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator*=(Numeric)")([&]() {
    Av *= 3.5;
    X = Av[0] + Av[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Vector::operator/=(Numeric)")([&]() {
    A /= 7.77;
    X = A[0] + A[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator/=(Numeric)")([&]() {
    Ae /= 7.77;
    X = Ae[0] + Ae[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator/=(Numeric)")([&]() {
    Av /= 7.77;
    X = Av[0] + Av[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("dummy")([results_]() { return results_[results_.size() - 1]; });

  return out;
}

Array<Timing> test_elementary_ops_Matrix(Index N) {
  Matrix A(N, N, 1);
  auto Ae = static_cast<MatrixView>(A);
  MatrixView Av = A;

  Array<Numeric> results_;
  Array<Timing> out;

  Numeric X;
  {
    A = 1;
    A += 1;
    A *= 2;
    A /= 2;
    A -= 1;
    X = A[0, 0] + A[N - 1, N - 1];
  }
  results_.push_back(X);

  out.emplace_back("Matrix::operator=(Numeric)")([&]() {
    A = 1;
    X = A[0, 0] + A[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator=(Numeric)")([&]() {
    Ae = 1;
    X = Ae[0, 0] + Ae[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator=(Numeric)")([&]() {
    Av = 1;
    X = Av[0, 0] + Av[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Matrix::operator+=(Numeric)")([&]() {
    A += 1.33;
    X = A[0, 0] + A[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator+=(Numeric)")([&]() {
    Ae += 1.33;
    X = Ae[0, 0] + Ae[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator+=(Numeric)")([&]() {
    Av += 1.33;
    X = Av[0, 0] + Av[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Matrix::operator-=(Numeric)")([&]() {
    A -= 1.5;
    X = A[0, 0] + A[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator-=(Numeric)")([&]() {
    Ae -= 1.5;
    X = Ae[0, 0] + Ae[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator-=(Numeric)")([&]() {
    Av -= 1.5;
    X = Av[0, 0] + Av[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Matrix::operator*=(Numeric)")([&]() {
    A *= 3.5;
    X = A[0, 0] + A[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator*=(Numeric)")([&]() {
    Ae *= 3.5;
    X = Ae[0, 0] + Ae[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator*=(Numeric)")([&]() {
    Av *= 3.5;
    X = Av[0, 0] + Av[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Matrix::operator/=(Numeric)")([&]() {
    A /= 7.77;
    X = A[0, 0] + A[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator/=(Numeric)")([&]() {
    Ae /= 7.77;
    X = Ae[0, 0] + Ae[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator/=(Numeric)")([&]() {
    Av /= 7.77;
    X = Av[0, 0] + Av[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("dummy")([results_]() { return results_[results_.size() - 1]; });

  return out;
}

Array<Timing> test_ops_Vector(Index N) {
  Vector a(N, 1.1);
  Vector b(N, 2.2);
  auto ae = static_cast<VectorView>(a);
  auto be = static_cast<VectorView>(b);
  VectorView av = a;
  VectorView bv = b;

  Array<Numeric> results_;
  Array<Timing> out;

  Numeric X;
  {
    a += b;
    a -= b;
    a *= b;
    a /= b;
    X = a[0] + a[N - 1];
  }
  results_.push_back(X);

  out.emplace_back("Vector::operator+=(Vector)")([&]() {
    a += b;
    X = a[0] + a[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator+=(VectorView)")([&]() {
    ae += be;
    X = ae[0] + ae[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator+=(VectorView)")([&]() {
    av += bv;
    X = av[0] + av[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Vector::operator-=(Vector)")([&]() {
    a -= b;
    X = a[0] + a[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator-=(VectorView)")([&]() {
    ae -= be;
    X = ae[0] + ae[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator-=(VectorView)")([&]() {
    av -= bv;
    X = av[0] + av[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Vector::operator*=(Vector)")([&]() {
    a *= b;
    X = a[0] + a[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator*=(VectorView)")([&]() {
    ae *= be;
    X = ae[0] + ae[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator*=(VectorView)")([&]() {
    av *= bv;
    X = av[0] + av[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Vector::operator/=(Vector)")([&]() {
    a /= b;
    X = a[0] + a[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator/=(VectorView)")([&]() {
    ae /= be;
    X = ae[0] + ae[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("VectorView::operator/=(VectorView)")([&]() {
    av /= bv;
    X = av[0] + av[N - 1];
  });
  results_.push_back(X);

  out.emplace_back("dummy")([results_]() { return results_[results_.size() - 1]; });

  return out;
}

Array<Timing> test_ops_Matrix(Index N) {
  Matrix a(N, N, 1.1);
  Matrix b(N, N, 2.2);
  auto ae = static_cast<MatrixView>(a);
  auto be = static_cast<MatrixView>(b);
  MatrixView av = a;
  MatrixView bv = b;

  Array<Numeric> results_;
  Array<Timing> out;

  Numeric X;
  {
    a += b;
    a -= b;
    a *= b;
    a /= b;
    X = a[0, 0] + a[N - 1, N - 1];
  }
  results_.push_back(X);

  out.emplace_back("Matrix::operator+=(Matrix)")([&]() {
    a += b;
    X = a[0, 0] + a[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator+=(MatrixView)")([&]() {
    ae += be;
    X = ae[0, 0] + ae[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator+=(MatrixView)")([&]() {
    av += bv;
    X = av[0, 0] + av[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Matrix::operator-=(Matrix)")([&]() {
    a -= b;
    X = a[0, 0] + a[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator-=(MatrixView)")([&]() {
    ae -= be;
    X = ae[0, 0] + ae[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator-=(MatrixView)")([&]() {
    av -= bv;
    X = av[0, 0] + av[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Matrix::operator*=(Matrix)")([&]() {
    a *= b;
    X = a[0, 0] + a[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator*=(MatrixView)")([&]() {
    ae *= be;
    X = ae[0, 0] + ae[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator*=(MatrixView)")([&]() {
    av *= bv;
    X = av[0, 0] + av[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("Matrix::operator/=(Matrix)")([&]() {
    a /= b;
    X = a[0, 0] + a[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator/=(MatrixView)")([&]() {
    ae /= be;
    X = ae[0, 0] + ae[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("MatrixView::operator/=(MatrixView)")([&]() {
    av /= bv;
    X = av[0, 0] + av[N - 1, N - 1];
  });
  results_.push_back(X);

  out.emplace_back("dummy")([results_]() { return results_[results_.size() - 1]; });

  return out;
}

int main(int argc, char** c) {
  std::array<Index, 8> N;
  if (static_cast<std::size_t>(argc) < 1 + 1 + N.size()) {
    std::cerr << "Expects PROGNAME NREPEAT NSIZE..., wehere NSIZE is "
              << N.size() << " indices\n";
    return EXIT_FAILURE;
  }

  const auto n = static_cast<Index>(std::atoll(c[1]));
  for (std::size_t i = 0; i < N.size(); i++)
    N[i] = static_cast<Index>(std::atoll(c[2 + i]));

  std::cout << n << " matpack-performance-tests\n\n";
  for (Index i = 0; i < n; i++) {
    std::cout << N[0] << " sums\n" << test_sum(N[0]) << '\n';
    std::cout << N[1] << " dots\n" << test_dot(N[1]) << '\n';
    std::cout << N[2] << " vector_mult\n"
              << test_Vector_mult(N[2]) << '\n';
    std::cout << N[3] << " matrix_mult\n"
              << test_Matrix_multiply(N[3]) << '\n';
    std::cout << N[4] << " numeric_ops_vector\n"
              << test_elementary_ops_Vector(N[4]) << '\n';
    std::cout << N[5] << " numeric_ops_matrix\n"
              << test_elementary_ops_Matrix(N[5]) << '\n';
    std::cout << N[6] << " vector_ops_vector\n" << test_ops_Vector(N[6]) << '\n';
    std::cout << N[7] << " matrix_ops_matrix\n" << test_ops_Matrix(N[7]) << '\n';
  }
}
