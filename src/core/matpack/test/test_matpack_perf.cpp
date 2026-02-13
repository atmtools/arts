#include <matpack.h>
#include <rng.h>
#include <time_report.h>

#include <cstdlib>
#include <iostream>

namespace {
Numeric dot_vector(const ConstVectorView& a, const ConstVectorView& b) {
  ARTS_NAMED_TIME_REPORT("dot_vector");

  return dot(a, b);
}

Numeric dot_tensor(const ConstTensor7View& a, const ConstTensor7View& b) {
  ARTS_NAMED_TIME_REPORT("dot_tensor");

  return dot(a, b);
}

Numeric vector_matrix_mult(VectorView& x,
                           const ConstMatrixView& A,
                           const ConstVectorView& y) {
  ARTS_NAMED_TIME_REPORT(
      std::format("vector_matrix_mult; {:B,} from {:B,} x {:B,}",
                  x.shape(),
                  A.shape(),
                  y.shape()));

  mult(x, A, y);

  return x.front();
}

Numeric matrix_matrix_mult(MatrixView& C,
                           const ConstMatrixView& A,
                           const ConstMatrixView& B) {
  ARTS_NAMED_TIME_REPORT(
      std::format("matrix_matrix_mult; {:B,} from {:B,} x {:B,}",
                  C.shape(),
                  A.shape(),
                  B.shape()));

  mult(C, A, B);

  return C.front();
}

Numeric plus_vector(VectorView& a, const ConstVectorView& b) {
  ARTS_NAMED_TIME_REPORT("vector +=");

  a += b;

  return a.front();
}

Numeric plus_tensor(Tensor7View& a, const ConstTensor7View& b) {
  ARTS_NAMED_TIME_REPORT("tensor7 +=");

  a += b;

  return a.front();
}

Numeric mult_vector(VectorView& a, const ConstVectorView& b) {
  ARTS_NAMED_TIME_REPORT("vector *=");

  a *= b;

  return a.front();
}

Numeric mult_tensor(Tensor7View& a, const ConstTensor7View& b) {
  ARTS_NAMED_TIME_REPORT("tensor7 *=");

  a *= b;

  return a.front();
}

Numeric sub_vector(VectorView& a, const ConstVectorView& b) {
  ARTS_NAMED_TIME_REPORT("vector -=");

  a -= b;

  return a.front();
}

Numeric sub_tensor(Tensor7View& a, const ConstTensor7View& b) {
  ARTS_NAMED_TIME_REPORT("tensor7 -=");

  a -= b;

  return a.front();
}

Numeric div_vector(VectorView& a, const ConstVectorView& b) {
  ARTS_NAMED_TIME_REPORT("vector /=");

  a /= b;

  return a.front();
}

Numeric div_tensor(Tensor7View& a, const ConstTensor7View& b) {
  ARTS_NAMED_TIME_REPORT("tensor7 /=");

  a /= b;

  return a.front();
}

Numeric sort_vector(VectorView& a) {
  ARTS_NAMED_TIME_REPORT("sort_vector");

  matpack::sort(a);

  return a.front();
}

Numeric sort_tensor(Tensor7View& a) {
  ARTS_NAMED_TIME_REPORT("sort_tensor");

  matpack::sort(a, {}, [](auto&& A) { return A.front(); });

  return a.front();
}

Numeric dot_sum(const ConstMatrixView& a) {
  ARTS_NAMED_TIME_REPORT(std::format("dot_sum; {:B,}", a.shape()));

  return std::transform_reduce(
      a.begin(), a.end(), Numeric{0}, std::plus<>{}, [](auto&& v) {
        return dot(v, v);
      });
}
}  // namespace

int main() {
  Numeric buf{};

  {
    constexpr Size N = 100'000'000;
    const Vector b   = random_numbers<1>({N});
    Vector a         = random_numbers<1>({N});

    buf += plus_vector(a, b);
    buf += sub_vector(a, b);
    buf += mult_vector(a, b);
    buf += div_vector(a, b);
    buf += dot_vector(a, b);

    buf += sort_vector(a);
  }

  {
    constexpr Size N = 10;
    const Tensor7 b  = random_numbers<7>({N, N, N, N, N, N, N * N});
    Tensor7 a        = random_numbers<7>({N, N, N, N, N, N, N * N});

    buf += plus_tensor(a, b);
    buf += sub_tensor(a, b);
    buf += mult_tensor(a, b);
    buf += div_tensor(a, b);
    buf += dot_tensor(a, b);

    buf += sort_tensor(a);
  }

  {
    constexpr Size N = 10'000;
    constexpr Size M = 100;
    constexpr Size P = 10'000;
    const Matrix B   = random_numbers<2>({M, P});
    const Matrix A   = random_numbers<2>({N, M});
    Matrix C         = random_numbers<2>({N, P});

    buf += matrix_matrix_mult(C, A, B);
  }

  {
    constexpr Size N = 10'000;
    constexpr Size M = 10'000;
    constexpr Size P = 100;
    const Matrix B   = random_numbers<2>({M, P});
    const Matrix A   = random_numbers<2>({N, M});
    Matrix C         = random_numbers<2>({N, P});

    buf += matrix_matrix_mult(C, A, B);
  }

  {
    constexpr Size N = 100;
    constexpr Size M = 10'000;
    constexpr Size P = 10'000;
    const Matrix B   = random_numbers<2>({M, P});
    const Matrix A   = random_numbers<2>({N, M});
    Matrix C         = random_numbers<2>({N, P});

    buf += matrix_matrix_mult(C, A, B);
  }

  {
    constexpr Size N = 100'000;
    constexpr Size M = 1'000;
    const Vector y   = random_numbers<1>({M});
    const Matrix A   = random_numbers<2>({N, M});
    Vector x         = random_numbers<1>({N});

    buf += vector_matrix_mult(x, A, y);
  }

  {
    constexpr Size N = 1'000;
    constexpr Size M = 100'000;
    const Vector y   = random_numbers<1>({M});
    const Matrix A   = random_numbers<2>({N, M});
    Vector x         = random_numbers<1>({N});

    buf += vector_matrix_mult(x, A, y);
  }

  {
    constexpr Size N = 100'000;
    constexpr Size M = 1'000;
    const Matrix A   = random_numbers<2>({N, M});

    buf += dot_sum(A);
  }

  {
    constexpr Size M = 100'000;
    constexpr Size N = 1'000;
    const Matrix A   = random_numbers<2>({N, M});

    buf += dot_sum(A);
  }

  std::println(std::cerr, "Prevent optimizing away: {}", buf);
  arts::print_report();
}
