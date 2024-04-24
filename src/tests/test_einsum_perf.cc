#include <matpack.h>

#include <cstdlib>
#include <iostream>

#include "test_perf.h"

void matvec(int N) {
  constexpr int NTEST = 10;

  constexpr Index n = 1000;
  constexpr Index m = 100000;

  Array<Timing> ts;
  ts.reserve(N * NTEST);

  Array<Numeric> some_results;
  some_results.reserve(ts.capacity());

  for (Index i = 0; i < N; i++) {
    ts.emplace_back("matvec-naive");
    ts.back()([&] {
      const Matrix A(n, m, 1);
      const Vector x(m, 1);
      Vector y(n);
      for (Index i = 0; i < n; i++) {
        y[i] = 0;
        for (Index j = 0; j < m; j++) {
          y[i] += A(i, j) * x[j];
        }
      }
      some_results.push_back(y[0]);
    });

    ts.emplace_back("matvec-partial-einsum");
    ts.back()([&] {
      const Matrix A(n, m, 1);
      const Vector x(m, 1);
      Vector y(n);
      for (Index i = 0; i < n; i++) {
        y[i] = einsum<Numeric, "", "m", "m">({}, A[i], x);
      }
      some_results.push_back(y[0]);
    });

    ts.emplace_back("matvec-partial-dotprod");
    ts.back()([&] {
      const Matrix A(n, m, 1);
      const Vector x(m, 1);
      Vector y(n);
      for (Index i = 0; i < n; i++) {
        y[i] = A[i] * x;
      }
      some_results.push_back(y[0]);
    });

    ts.emplace_back("matvec-lapack");
    ts.back()([&] {
      const Matrix A(n, m, 1);
      const Vector x(m, 1);
      Vector y(n);
      mult(y, A, x);
      some_results.push_back(y[0]);
    });

    ts.emplace_back("matvec-einsum-view");
    ts.back()([&] {
      const Matrix A(n, m, 1);
      const Vector x(m, 1);
      Vector y(n, 1);
      einsum<"n", "nm", "m">(y, A, x);
      some_results.push_back(y[0]);
    });

    ts.emplace_back("matvec-einsum-ret");
    ts.back()([&] {
      const Matrix A(n, m, 1);
      const Vector x(m, 1);
      const auto y = einsum<Vector, "n", "nm", "m">({n}, A, x);
      some_results.push_back(y[0]);
    });

    ts.emplace_back("matvec-einsum-view-reverse");
    ts.back()([&] {
      const Matrix A(n, m, 1);
      Vector x(m, 1);
      const Vector y(n, 1);
      einsum<"m", "nm", "n">(x, A, y);
      some_results.push_back(x[0]);
    });

    ts.emplace_back("matvec-einsum-ret-reverse");
    ts.back()([&] {
      const Matrix A(n, m, 1);
      const Vector y(n, 1);
      const auto x = einsum<Vector, "m", "nm", "n">({m}, A, y);
      some_results.push_back(x[0]);
    });
  }

  std::cout << ts;
}

void matmat(int N) {
  constexpr int NTEST = 10;

  constexpr Index n = 100;
  constexpr Index m = 1000;
  constexpr Index p = 10000;

  Array<Timing> ts;
  ts.reserve(N * NTEST);

  Array<Numeric> some_results;
  some_results.reserve(ts.capacity());

  for (Index i = 0; i < N; i++) {
    ts.emplace_back("matmat-naive");
    ts.back()([&] {
      const Matrix A(m, n, 1);
      const Matrix B(n, p, 1);
      Matrix C(m, p);
      for (Index i = 0; i < m; i++) {
        for (Index j = 0; j < p; j++) {
          C(i, j) = 0;
          for (Index k = 0; k < n; k++) {
            C(i, j) += A(i, k) * B(k, j);
          }
        }
      }
      some_results.push_back(C(0, 0));
    });

    ts.emplace_back("matmat-partial-einsum");
    ts.back()([&] {
      const Matrix A(m, n, 1);
      const Matrix B(n, p, 1);
      Matrix C(m, p);
      for (Index i = 0; i < m; i++) {
        einsum<"p", "n", "np">(C[i], A[i], B);
      }
      some_results.push_back(C(0, 0));
    });

    ts.emplace_back("matmat-naive-dotprod");
    ts.back()([&] {
      const Matrix A(m, n, 1);
      const Matrix B(n, p, 1);
      Matrix C(m, p);
      for (Index i = 0; i < m; i++) {
        for (Index j = 0; j < p; j++) {
          C(i, j) = A[i] * B(joker, j);
        }
      }
      some_results.push_back(C(0, 0));
    });

    ts.emplace_back("matmat-naive-optimized");
    ts.back()([&] {
      const Matrix A(m, n, 1);
      const Matrix B(n, p, 1);
      Matrix C(m, p);
      C = 0;
      for (Index i = 0; i < m; i++) {
        for (Index k = 0; k < n; k++) {
          for (Index j = 0; j < p; j++) {
            C(i, j) += A(i, k) * B(k, j);
          }
        }
      }
      some_results.push_back(C(0, 0));
    });

    ts.emplace_back("matmat-partial-einsum-matvec");
    ts.back()([&] {
      const Matrix A(m, n, 1);
      const Matrix B(n, p, 1);
      Matrix C(m, p);
      C = 0;
      for (Index i = 0; i < p; i++) {
        einsum<"m", "mn", "n">(C(joker, i), A, B(joker, i));
      }
      some_results.push_back(C(0, 0));
    });

    ts.emplace_back("matmat-partial-lapack-matvec");
    ts.back()([&] {
      const Matrix A(m, n, 1);
      const Matrix B(n, p, 1);
      Matrix C(m, p);
      C = 0;
      for (Index i = 0; i < p; i++) {
        mult(C(joker, i), A, B(joker, i));
      }
      some_results.push_back(C(0, 0));
    });

    ts.emplace_back("matmat-lapack");
    ts.back()([&] {
      const Matrix A(m, n, 1);
      const Matrix B(n, p, 1);
      Matrix C(m, p);
      mult(C, A, B);
      some_results.push_back(C(0, 0));
    });

    ts.emplace_back("matmat-einsum-view");
    ts.back()([&] {
      const Matrix A(m, n, 1);
      const Matrix B(n, p, 1);
      Matrix C(m, p);
      einsum<"mp", "mn", "np">(C, A, B);
      some_results.push_back(C(0, 0));
    });

    ts.emplace_back("matmat-einsum-retval");
    ts.back()([&] {
      const Matrix A(m, n, 1);
      const Matrix B(n, p, 1);
      const auto C = einsum<Matrix, "mp", "mn", "np">({m, p}, A, B);
      some_results.push_back(C(0, 0));
    });
  }

  std::cout << ts;
}

void pathing(int N) {
  constexpr int NTEST = 10;

  constexpr Index i = 10;
  constexpr Index j = 100;
  constexpr Index k = 1000;
  constexpr Index l = 1000;

  Array<Timing> ts;
  ts.reserve(N * NTEST);

  Array<Numeric> some_results;
  some_results.reserve(ts.capacity());

  for (Index n = 0; n < N; n++) {
    ts.emplace_back("ij,jk,kl->il-direct-einsum");
    ts.back()([&] {
      const Matrix A(i, j, .1);
      const Matrix B(j, k, .2);
      Matrix C(k, l, .3);
      C(0, 0) = 2;
      Matrix D(i, l);
      einsum<"il", "ij", "jk", "kl">(D, A, B, C);
      some_results.push_back(D(0, 0));
    });

    ts.emplace_back("ij,jk,kl->il-np.einsum_path-einsum");
    ts.back()([&] {
      const Matrix A(i, j, .1);
      const Matrix B(j, k, .2);
      Matrix C(k, l, .3);
      C(0, 0) = 2;
      Matrix CB(j, l);
      Matrix D(i, l);
      einsum<"jl", "kl", "jk">(CB, C, B);
      einsum<"il", "jl", "ij">(D, CB, A);
      some_results.push_back(D(0, 0));
    });

    ts.emplace_back("ea,fb,abcd,gc,hd->efgh-direct-einsum");
    ts.back()([&]() {
      const Tensor4 I(i, i, i, i, 1);
      Matrix C(i, i, 1);
      C(0, 0) = 2;
      Tensor4 D(i, i, i, i);
      einsum<"efgh", "ea", "fb", "abcd", "gc", "hd">(D, C, C, I, C, C);
      some_results.push_back(D(0, 0, 0, 0));
    });

    ts.emplace_back("ea,fb,abcd,gc,hd->efgh-np.einsum_path-einsum");
    ts.back()([&]() {
      const Tensor4 I(i, i, i, i, 1);
      Matrix C(i, i, 1);
      C(0, 0) = 2;
      Tensor4 B(i, i, i, i);
      Tensor4 D(i, i, i, i);
      einsum<"bcde", "abcd", "ea">(B, I, C);
      einsum<"cdef", "bcde", "fb">(D, B, C);
      einsum<"defg", "cdef", "gc">(B, D, C);
      einsum<"efgh", "defg", "hd">(D, B, C);
      some_results.push_back(D(0, 0, 0, 0));
    });
  }

  std::cout << ts;
}

int main() try {
  std::cout << "einsum-perf-test\n";

  std::cout << "matrix-vector-mult\n";
  matvec(20);

  std::cout << "matrix-matrix-mult\n";
  matmat(10);

  std::cout << "trying-np.einsum_path\n";
  pathing(10);

  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::cerr << "Error: " << e.what() << '\n';
  return EXIT_FAILURE;
}