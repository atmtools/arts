#include <debug.h>

#include <algorithm>
#include <complex>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "matpack.h"
#include "matpack_mdspan_helpers.h"
#include "matpack_mdspan_helpers_eigen.h"

namespace {
void test_strided_view() {
  Vector x{1, 2, 3, 4, 5, 6, 7, 8};

  ARTS_ASSERT(x[0] != 0 and x[2] != 0 and x[4] != 0, "{:B}", x);

  x[StridedRange(0, 3, 2)] = 0;

  ARTS_ASSERT(x[0] == 0 and x[2] == 0 and x[4] == 0, "{:B}", x);
}

void test_view() {
  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1, 2, 3, 4, 5, 6, 7, 8};
    ConstTensor3View y{{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i = 0; i < 2; i++) {
      [[maybe_unused]] auto g = y[i, joker, joker][i, i];
      [[maybe_unused]] auto h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric>);
    }

    Numeric d = 1.0;
    for (auto a : y)
      for (auto b : a)
        for (auto c : b) {
          ARTS_USER_ERROR_IF(c != d, "{} {}", c, d);
          d += 1.0;
        }
  }

  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1, 2, 3, 4, 5, 6, 7, 8};
    Tensor3View y{{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i = 0; i < 2; i++) {
      [[maybe_unused]] auto& g = y[i, joker, joker][i, i];
      [[maybe_unused]] auto& h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric&>);
    }

    Numeric d = 1.0;
    for (const auto& a : y)
      for (auto b : a)
        for (auto c : b) {
          ARTS_USER_ERROR_IF(c != d, "{} {}", c, d);
          d += 1.0;
        }
  }

  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1, 2, 3, 4, 5, 6, 7, 8};
    ConstTensor3View y{{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i = 0; i < 2; i++) {
      [[maybe_unused]] auto g = y[i, joker, joker][i, i];
      [[maybe_unused]] auto h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric>);
    }

    Numeric d = 1.0;
    for (auto a : y)
      for (auto b : a)
        for (auto c : b) {
          ARTS_USER_ERROR_IF(c != d, "{} {}", c, d);
          d += 1.0;
        }
  }

  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1, 2, 3, 4, 5, 6, 7, 8};
    Tensor3View y{{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i = 0; i < 2; i++) {
      [[maybe_unused]] auto& g = y[i, joker, joker][i, i];
      [[maybe_unused]] auto& h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric&>);
    }

    Numeric d = 1.0;
    for (const auto& a : y)
      for (auto b : a)
        for (auto c : b) {
          ARTS_USER_ERROR_IF(c != d, "{} {}", c, d);
          d += 1.0;
        }
  }

  //! Allow assignment between views
  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> xr{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> xi{0, 0, 0, 0, 0, 0, 0, 0};
    ComplexTensor3View y{{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const Tensor3View yr{{xr.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const Tensor3View yi{{xi.data(), std::array{Index{2}, Index{2}, Index{2}}}};

    ARTS_ASSERT(yr != yi, "{:B,} != {:B,}", yr, yi)
    ARTS_ASSERT(y.real() != y.imag(), "{:B,} != {:B,}", y.real(), y.imag())
    ARTS_ASSERT(y.real() == yr, "{:B,} == {:B,}", y.real(), yr)
    ARTS_ASSERT(y.imag() == yi, "{:B,} == {:B,}", y.imag(), yi)
    ARTS_ASSERT(y.imag() != yr, "{:B,} != {:B,}", y.imag(), yr)

    y.imag() = y.real();

    ARTS_ASSERT(yr != yi)
    ARTS_ASSERT(y.real() == y.imag())
    ARTS_ASSERT(y.real() == yr)
    ARTS_ASSERT(y.imag() != yi)
    ARTS_ASSERT(y.imag() == yr)
  }
  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> xr{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> xi{0, 0, 0, 0, 0, 0, 0, 0};
    ComplexTensor3View y{{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const Tensor3View yr{{xr.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const Tensor3View yi{{xi.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    auto z = y[joker, StridedRange(0, 2, 2), joker];

    std::vector<Numeric> reszr{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> reszi{1, 2, 0, 0, 5, 6, 0, 0};
    const Tensor3View zr{
        {reszr.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const Tensor3View zi{
        {reszi.data(), std::array{Index{2}, Index{2}, Index{2}}}};

    ARTS_ASSERT(yr != yi)
    ARTS_ASSERT(y.real() != y.imag())
    ARTS_ASSERT(z.real() != z.imag())
    ARTS_ASSERT(y.real() == yr)
    ARTS_ASSERT(y.imag() == yi)
    ARTS_ASSERT(y.real() == zr)
    ARTS_ASSERT(y.imag() != zi)

    z.imag() = z.real();

    ARTS_ASSERT(yr != yi)
    for(auto x: elemwise_range(y.real())) std::print("{} ", x);
    std::print("\n");
    for(auto x: elemwise_range(y.imag())) std::print("{} ", x);
    std::print("\n");
    std::print("\n");
    ARTS_ASSERT(y.real() != y.imag(), "{:B,} != {:B,}", y.real(), y.imag())
    ARTS_ASSERT(z.real() == z.imag())
    ARTS_ASSERT(y.real() == yr)
    ARTS_ASSERT(y.imag() != yi)
    ARTS_ASSERT(y.real() == zr)
    ARTS_ASSERT(y.imag() == zi)
  }

  //! Test basic +=, -=, *=, /=, = for arithmetic change or assignment
  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexTensor3View y{{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};

    std::vector<std::complex<Numeric>> copy = x;
    const ComplexTensor3View test{
        {copy.data(), std::array{Index{2}, Index{2}, Index{2}}}};

    ARTS_ASSERT(y == test)

    y += 3.;
    for (auto& a : copy) a += 3;
    ARTS_ASSERT(y == test)

    y /= 3.;
    for (auto& a : copy) a /= 3;
    ARTS_ASSERT(y == test)

    y.imag() -= 2;
    for (auto& a : copy) a -= std::complex<Numeric>(0, 2);
    ARTS_ASSERT(y == test)

    y -= 5.;
    for (auto& a : copy) a -= 5;
    ARTS_ASSERT(y == test)

    y *= 5.;
    for (auto& a : copy) a *= 5;
    ARTS_ASSERT(y == test)

    y = 42;
    for (auto& a : copy) a = 42;
    ARTS_ASSERT(y == test)
  }

  //! Test external assignment
  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexVectorView y{{x.data(), std::array<Index, 1>{8}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexVectorView test{{copy.data(), std::array<Index, 1>{8}}};

    for (auto& a : copy) a += std::complex<Numeric>(3, 2);

    ARTS_ASSERT(y != test)

    y = copy;

    ARTS_ASSERT(y == test)
  }
}

void test_eigen() {
  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexMatrixView y{{x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{{copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_ASSERT(y == test)

    y = 2 * y;

    ARTS_ASSERT(y != test)
    for (auto& a : copy) a *= 2;
    ARTS_ASSERT(y == test)

    ComplexMatrix SSS;
    SSS = 2 * y;
  }

  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexVectorView y{{x.data(), std::array<Index, 1>{8}}};

    const std::complex<Numeric> z = dot(y, y);
    ARTS_USER_ERROR_IF(
        z != std::transform_reduce(
                 x.begin(), x.end(), x.begin(), std::complex<Numeric>{0}),
        "{}",
        z)
  }

  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexMatrixView y{{x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{{copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_ASSERT(y == test)

    y = 2 * y + y;

    ARTS_ASSERT(y != test)
    for (auto& a : copy) a = 2.0 * a + a;
    ARTS_ASSERT(y == test)
  }

  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexMatrixView y{{x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{{copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_ASSERT(y == test)

    y = 2 * y - y;

    ARTS_ASSERT(y == test)
  }
}

void test_data() {
  matpack::data_t<Numeric, 3> x(4, 2, 3, 2.1);
  std::print("{}\n\n", x);
  std::print("{}\n\n", x[0]);
  std::print("{}\n\n", x[0][0]);
  std::print("{}\n\n", x[0][0][0]);
  x[0][0][0] += 1;
  std::print("{}\n\n", x);

  auto y = std::move(x).flatten();
  std::print("{}\n\n", y);
  std::print("{}\n\n", x);

  matpack::data_t<Numeric, 1> z;
  z.swap(y);
}

void test_complex() {
  {
    Complex x{0, 0};
    const Complex y{0, 0};
    ARTS_USER_ERROR_IF(x != y, "{} {}", x, y)
    ARTS_ASSERT(real_val(x) == real_val(y))
    ARTS_ASSERT(imag_val(x) == imag_val(y))
  }

  {
    Complex x{1, 1};
    x = 2 + x;
    ARTS_ASSERT(x == (Complex{3, 1}))
    x = 2 - x;
    ARTS_ASSERT(x == (Complex{-1, -1}))
    x = 2 * x;
    ARTS_ASSERT(x == (Complex{-2, -2}))
    x = 2 / x;
    ARTS_ASSERT(x == (Complex{-0.5, 0.5}))
  }

  {
    Complex x{1, 1};
    x = x + 2;
    ARTS_ASSERT(x == (Complex{3, 1}))
    x = x - 2;
    ARTS_ASSERT(x == (Complex{1, 1}))
    x = x * 2;
    ARTS_ASSERT(x == (Complex{2, 2}))
    x = x / 2;
    ARTS_ASSERT(x == (Complex{1, 1}))
  }

  {
    Complex x{1, 1};
    x = 2 + x;
    ARTS_ASSERT(x == (Complex{3, 1}))
    x = 2 - x;
    ARTS_ASSERT(x == (Complex{-1, -1}))
    x = 2 * x;
    ARTS_ASSERT(x == (Complex{-2, -2}))
    x = 2 / x;
    ARTS_ASSERT(x == (Complex{-0.5, 0.5}))
  }

  {
    Complex x{1, 1};
    x = x + 2;
    ARTS_ASSERT(x == (Complex{3, 1}))
    x = x - 2;
    ARTS_ASSERT(x == (Complex{1, 1}))
    x = x * 2;
    ARTS_ASSERT(x == (Complex{2, 2}))
    x = x / 2;
    ARTS_ASSERT(x == (Complex{1, 1}))
  }
}

void test_math() {
  {
    Vector x{0, 1, 2, 3, 4};

    ARTS_ASSERT(max(x) == 4)
    ARTS_ASSERT(max(VectorView{x}) == 4)
    ARTS_ASSERT(max(ConstVectorView{x}) == 4)
    ARTS_ASSERT(max(StridedVectorView{x}) == 4)
    ARTS_ASSERT(max(StridedConstVectorView{x}) == 4)

    ARTS_ASSERT(min(x) == 0, "min(x) == {}\n {:B,}", min(x), x)
    ARTS_ASSERT(min(VectorView{x}) == 0)
    ARTS_ASSERT(min(ConstVectorView{x}) == 0)
    ARTS_ASSERT(min(StridedVectorView{x}) == 0)
    ARTS_ASSERT(min(StridedConstVectorView{x}) == 0)

    static_assert(std::random_access_iterator<decltype(x.elem_begin())>);
    static_assert(
        std::random_access_iterator<decltype(VectorView{x}.elem_begin())>);
    static_assert(
        std::random_access_iterator<decltype(ConstVectorView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<
                  decltype(StridedVectorView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<
                  decltype(StridedConstVectorView{x}.elem_begin())>);

    static_assert(std::random_access_iterator<decltype(x.begin())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.begin())>);
    static_assert(
        std::random_access_iterator<decltype(ConstVectorView{x}.begin())>);
    static_assert(
        std::random_access_iterator<decltype(StridedVectorView{x}.begin())>);
    static_assert(std::random_access_iterator<
                  decltype(StridedConstVectorView{x}.begin())>);

    static_assert(std::random_access_iterator<decltype(x.elem_end())>);
    static_assert(
        std::random_access_iterator<decltype(VectorView{x}.elem_end())>);
    static_assert(
        std::random_access_iterator<decltype(ConstVectorView{x}.elem_end())>);
    static_assert(
        std::random_access_iterator<decltype(StridedVectorView{x}.elem_end())>);
    static_assert(std::random_access_iterator<
                  decltype(StridedConstVectorView{x}.elem_end())>);

    static_assert(std::random_access_iterator<decltype(x.end())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.end())>);
    static_assert(
        std::random_access_iterator<decltype(ConstVectorView{x}.end())>);
    static_assert(
        std::random_access_iterator<decltype(StridedVectorView{x}.end())>);
    static_assert(
        std::random_access_iterator<decltype(StridedConstVectorView{x}.end())>);
  }

  {
    Matrix x(3, 3, 0);
    ARTS_ASSERT(min(x) == 0)
    ARTS_ASSERT(min(MatrixView{x}) == 0)
    ARTS_ASSERT(min(ConstMatrixView{x}) == 0)

    static_assert(std::random_access_iterator<decltype(x.elem_begin())>);
    static_assert(
        std::random_access_iterator<decltype(MatrixView{x}.elem_begin())>);
    static_assert(
        std::random_access_iterator<decltype(ConstMatrixView{x}.elem_begin())>);
    static_assert(
        std::random_access_iterator<decltype(MatrixView{x}.elem_begin())>);
    static_assert(
        std::random_access_iterator<decltype(ConstMatrixView{x}.elem_begin())>);

    static_assert(std::random_access_iterator<decltype(x.begin())>);
    static_assert(std::random_access_iterator<decltype(MatrixView{x}.begin())>);
    static_assert(
        std::random_access_iterator<decltype(ConstMatrixView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(MatrixView{x}.begin())>);
    static_assert(
        std::random_access_iterator<decltype(ConstMatrixView{x}.begin())>);

    static_assert(std::random_access_iterator<decltype(x.elem_end())>);
    static_assert(
        std::random_access_iterator<decltype(MatrixView{x}.elem_end())>);
    static_assert(
        std::random_access_iterator<decltype(ConstMatrixView{x}.elem_end())>);
    static_assert(
        std::random_access_iterator<decltype(MatrixView{x}.elem_end())>);
    static_assert(
        std::random_access_iterator<decltype(ConstMatrixView{x}.elem_end())>);

    static_assert(std::random_access_iterator<decltype(x.end())>);
    static_assert(std::random_access_iterator<decltype(MatrixView{x}.end())>);
    static_assert(
        std::random_access_iterator<decltype(ConstMatrixView{x}.end())>);
    static_assert(std::random_access_iterator<decltype(MatrixView{x}.end())>);
    static_assert(
        std::random_access_iterator<decltype(ConstMatrixView{x}.end())>);
  }
}

void test_mult() {
  {
    Matrix A(4, 4, 1);
    for (Index i = 0; i < 4; i++)
      for (Index j = 0; j < 4; j++) A[i, j] = static_cast<Numeric>(i);
    const Vector y({1, 2, 3, 4});
    Vector x(4);
    Vector eig_x{A * y};
    mult(x, A, y);
    ARTS_USER_ERROR_IF(eig_x not_eq x, "{} vs {}", eig_x, x)
  }
  {
    ComplexMatrix A(4, 4, 1);
    for (Index i = 0; i < 4; i++)
      for (Index j = 0; j < 4; j++)
        A[i, j] = {static_cast<Numeric>(i), 2 * static_cast<Numeric>(i)};
    ComplexVector y({1, 2, 3, 4}), x(4);

    ComplexVector eig_x{A * y};
    mult(x, A, y);
    ARTS_USER_ERROR_IF(eig_x not_eq x, "{} vs {}", eig_x, x)
  }
}

void test_const_data() {
  {
    constexpr matpack::cdata_t<Numeric, 3> x{1, 2, 3};
    Vector y(x);

    ARTS_USER_ERROR_IF(x != y, "Bad")
    ARTS_USER_ERROR_IF((x[2] != y.back()), "Bad")
    ARTS_USER_ERROR_IF((x[0] != y.front()), "Bad")
  }
  {
    constexpr matpack::cdata_t<Numeric, 3, 3> x{2, 3, 4};
    Matrix y(x);

    ARTS_USER_ERROR_IF(x != y, "Bad")
    ARTS_USER_ERROR_IF((x[2, 2] != y[2, 2]), "Bad")
    ARTS_USER_ERROR_IF((x[0, 0] != y[0, 0]), "Bad")
  }
}

void test_my_interp() {
  my_interp::Lagrange<1, true> x(
      0, 3.5, matpack::cdata_t<Numeric, 7>{1, 2, 3, 4, 5, 6, 7});
  my_interp::Lagrange<1, false> x2(
      0, 3.5, matpack::cdata_t<Numeric, 7>{1, 2, 3, 4, 5, 6, 7});
  auto iw = interpweights(std::array{x, x}, std::array{x2, x2});

  interp(Matrix(7, 7, 1), interpweights(x, x2), x, x2);

  [[maybe_unused]] auto f =
      reinterp(Matrix(7, 7, 1), iw, std::array{x, x}, std::array{x2, x2});

  Matrix Z  = matpack::uniform_grid(-5, 49, 0.2).reshape(7, 7);
  Vector X  = matpack::uniform_grid(1, 7, 1);
  Vector Y  = matpack::uniform_grid(1, 7, 1);
  Vector XN = matpack::uniform_grid(1, 14, 0.5);
  Vector YN = matpack::uniform_grid(1, 28, 0.25);

  for (auto mx : XN) {
    for (auto my : YN) {
      auto lagx = my_interp::Lagrange<1, false>(0, mx, X);
      auto lagy = my_interp::Lagrange<1, false>(0, my, Y);
      auto iw2  = my_interp::interpweights(lagx, lagy);
      ARTS_USER_ERROR_IF(
          std::abs(interp(Z, iw2, lagx, lagy) - interp(Z, lagx, lagy)) > 1e-16,
          "Bad");
    }
  }
}

void test_sorted_grid() {
  {
    const AscendingGrid x{1, 2, 3, 4, 5, 6, 7, 8};
    ARTS_USER_ERROR_IF(x.size() != 8,
                       "Should be working, just testing existence of size()");
  }
  {
    bool correctly_fails;
    try {
      const AscendingGrid x{1, 2, 3, 4, 5, 6, 7, 8, 8};
      correctly_fails = false;
    } catch (...) {
      correctly_fails = true;
    }
    ARTS_USER_ERROR_IF(
        not correctly_fails,
        "Should not be able to initialize with non-ascending grid")
  }
}

void test_lapack_vector_mult() {
  Vector x{1, 2, 3};
  Matrix A = Vector{1, 2, 3, 4, 5, 6, 7, 8, 9}.reshape(3, 3);
  Vector y{0, 0, 0};

  mult(y, A, x);
  ARTS_USER_ERROR_IF((y != Vector{14., 32., 50.}), "Bad values:\n", y);

  mult(y, A, x, 0.5);
  ARTS_USER_ERROR_IF((y != Vector{7., 16., 25.}), "Bad values:\n", y);

  mult(y, A, x, 1.0, 1.0);
  ARTS_USER_ERROR_IF((y != Vector{21., 48., 75.}), "Bad values:\n", y);
}

void test_grid() {
  AscendingGrid x{1, 2, 3, 4, 5, 6, 7, 8};

  AscendingGridView y{x};

  std::print("{}\n", x);
  std::print("{}\n", y);

  y = {1, 2, 3, 4, 5, 6, 7, 9};
  std::print("{}\n", x);
  std::print("{}\n", y);

  x.emplace_back(10);
  std::print("{}\n", x);
  std::print("{}\n", y);

  bool set_low;
  try {
    x.emplace_back(8);
    set_low = true;
  } catch (...) {
    set_low = false;
  }
  ARTS_USER_ERROR_IF(set_low, "Should not be able to set low value");

  bool set_same;
  try {
    x.emplace_back(x.back());
    set_same = true;
  } catch (...) {
    {};
    set_same = false;
  }
  ARTS_USER_ERROR_IF(set_same, "Should not be able to set same value");

  Vector z{1, 2, 3, 4, 5, 6, 7, 9, 20};
  y = z;
  std::print("{}\n", x);
  std::print("{}\n", y);

  y = VectorView{z};
  std::print("{}\n", x);
  std::print("{}\n", y);
}
}  // namespace

#define EXECUTE_TEST(X)                                                       \
  std::cout << "#########################################################\n"; \
  std::cout << "Executing test: " #X << '\n';                                 \
  std::cout << "#########################################################\n"; \
  X();                                                                        \
  std::cout << "#########################################################\n";

int main() {
  EXECUTE_TEST(test_strided_view)
  EXECUTE_TEST(test_view)
  EXECUTE_TEST(test_eigen)
  EXECUTE_TEST(test_data)
  EXECUTE_TEST(test_complex)
  EXECUTE_TEST(test_math)
  EXECUTE_TEST(test_mult)
  EXECUTE_TEST(test_const_data)
  EXECUTE_TEST(test_my_interp)
  EXECUTE_TEST(test_sorted_grid)
  EXECUTE_TEST(test_lapack_vector_mult)
  EXECUTE_TEST(test_grid)
}
