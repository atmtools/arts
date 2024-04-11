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

#include "artstime/artstime.h"
#include "debug.h"
#include "interp.h"
#include "lin_alg.h"
#include "logic.h"
#include "matpack_algo.h"
#include "matpack_arrays.h"
#include "matpack_complex.h"
#include "matpack_concepts.h"
#include "matpack_constexpr.h"
#include "matpack_data.h"
#include "matpack_eigen.h"
#include "matpack_einsum.h"
#include "matpack_iter.h"
#include "matpack_math.h"
#include "matpack_view.h"
#include "sorted_grid.h"

void test_view() {
  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1, 2, 3, 4, 5, 6, 7, 8};
    ConstTensor3View y{matpack::exhaustive_mdspan<Numeric, 3>{
        x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i = 0; i < 2; i++) {
      [[maybe_unused]] auto g =
          y(i, joker, matpack::matpack_strided_access{joker})(i, i);
      [[maybe_unused]] auto h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric>);
    }

    Numeric d = 1.0;
    for (auto a : y)
      for (auto b : a)
        for (auto c : b) {
          ARTS_USER_ERROR_IF(c != d, c, ' ', d);
          d += 1.0;
        }
  }

  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1, 2, 3, 4, 5, 6, 7, 8};
    Tensor3View y{matpack::exhaustive_mdspan<Numeric, 3>{
        x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i = 0; i < 2; i++) {
      [[maybe_unused]] auto& g =
          y(i, joker, matpack::matpack_strided_access{joker})(i, i);
      [[maybe_unused]] auto& h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric&>);
    }

    Numeric d = 1.0;
    for (const auto& a : y)
      for (auto b : a)
        for (auto c : b) {
          ARTS_USER_ERROR_IF(c != d, c, ' ', d);
          d += 1.0;
        }
  }

  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1, 2, 3, 4, 5, 6, 7, 8};
    ExhaustiveConstTensor3View y{matpack::exhaustive_mdspan<Numeric, 3>{
        x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i = 0; i < 2; i++) {
      [[maybe_unused]] auto g =
          y(i, joker, matpack::matpack_strided_access{joker})(i, i);
      [[maybe_unused]] auto h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric>);
    }

    Numeric d = 1.0;
    for (auto a : y)
      for (auto b : a)
        for (auto c : b) {
          ARTS_USER_ERROR_IF(c != d, c, ' ', d);
          d += 1.0;
        }
  }

  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1, 2, 3, 4, 5, 6, 7, 8};
    ExhaustiveTensor3View y{matpack::exhaustive_mdspan<Numeric, 3>{
        x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i = 0; i < 2; i++) {
      [[maybe_unused]] auto& g =
          y(i, joker, matpack::matpack_strided_access{joker})(i, i);
      [[maybe_unused]] auto& h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric&>);
    }

    Numeric d = 1.0;
    for (const auto& a : y)
      for (auto b : a)
        for (auto c : b) {
          ARTS_USER_ERROR_IF(c != d, c, ' ', d);
          d += 1.0;
        }
  }

  //! Allow assignment between views
  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> xr{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> xi{0, 0, 0, 0, 0, 0, 0, 0};
    ExhaustiveComplexTensor3View y{
        matpack::exhaustive_mdspan<std::complex<Numeric>, 3>{
            x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View yr{matpack::exhaustive_mdspan<Numeric, 3>{
        xr.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View yi{matpack::exhaustive_mdspan<Numeric, 3>{
        xi.data(), std::array{Index{2}, Index{2}, Index{2}}}};

    ARTS_ASSERT(yr != yi)
    ARTS_ASSERT(y.real() != y.imag())
    ARTS_ASSERT(y.real() == yr)
    ARTS_ASSERT(y.imag() == yi)
    ARTS_ASSERT(y.imag() != yr)

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
    ExhaustiveComplexTensor3View y{
        matpack::exhaustive_mdspan<std::complex<Numeric>, 3>{
            x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View yr{matpack::exhaustive_mdspan<Numeric, 3>{
        xr.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View yi{matpack::exhaustive_mdspan<Numeric, 3>{
        xi.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    auto z = y(joker, matpack::matpack_strided_access(0, -1, 2), joker);

    std::vector<Numeric> reszr{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> reszi{1, 2, 0, 0, 5, 6, 0, 0};
    const ExhaustiveTensor3View zr{matpack::exhaustive_mdspan<Numeric, 3>{
        reszr.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View zi{matpack::exhaustive_mdspan<Numeric, 3>{
        reszi.data(), std::array{Index{2}, Index{2}, Index{2}}}};

    ARTS_ASSERT(yr != yi)
    ARTS_ASSERT(y.real() != y.imag())
    ARTS_ASSERT(z.real() != z.imag())
    ARTS_ASSERT(y.real() == yr)
    ARTS_ASSERT(y.imag() == yi)
    ARTS_ASSERT(y.real() == zr)
    ARTS_ASSERT(y.imag() != zi)

    z.imag() = z.real();

    ARTS_ASSERT(yr != yi)
    ARTS_ASSERT(y.real() != y.imag())
    ARTS_ASSERT(z.real() == z.imag())
    ARTS_ASSERT(y.real() == yr)
    ARTS_ASSERT(y.imag() != yi)
    ARTS_ASSERT(y.real() == zr)
    ARTS_ASSERT(y.imag() == zi)
  }

  //! Test basic +=, -=, *=, /=, = for arithmetic change or assignment
  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ExhaustiveComplexTensor3View y{
        matpack::exhaustive_mdspan<std::complex<Numeric>, 3>{
            x.data(), std::array{Index{2}, Index{2}, Index{2}}}};

    std::vector<std::complex<Numeric>> copy = x;
    const ExhaustiveComplexTensor3View test{
        matpack::exhaustive_mdspan<std::complex<Numeric>, 3>{
            copy.data(), std::array{Index{2}, Index{2}, Index{2}}}};

    ARTS_ASSERT(y == test)

    y += 3;
    for (auto& a : copy) a += 3;
    ARTS_ASSERT(y == test)

    y /= 3;
    for (auto& a : copy) a /= 3;
    ARTS_ASSERT(y == test)

    y.imag() -= 2;
    for (auto& a : copy) a -= std::complex<Numeric>(0, 2);
    ARTS_ASSERT(y == test)

    y -= 5;
    for (auto& a : copy) a -= 5;
    ARTS_ASSERT(y == test)

    y *= 5;
    for (auto& a : copy) a *= 5;
    ARTS_ASSERT(y == test)

    y = 42;
    for (auto& a : copy) a = 42;
    ARTS_ASSERT(y == test)
  }

  //! Test external assignment
  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexVectorView y{
        matpack::exhaustive_mdspan<std::complex<Numeric>, 1>{x.data(), 8}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexVectorView test{
        matpack::exhaustive_mdspan<std::complex<Numeric>, 1>{copy.data(), 8}};

    for (auto& a : copy) a += std::complex<Numeric>(3, 2);

    ARTS_ASSERT(y != test)

    y = copy;

    ARTS_ASSERT(y == test)
  }
}

void test_eigen() {
  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexMatrixView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{
        x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{
        matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{
            copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_ASSERT(y == test)

    y = 2 * matpack::eigen::as_eigen(y);

    ARTS_ASSERT(y != test)
    for (auto& a : copy) a *= 2;
    ARTS_ASSERT(y == test)
  }

  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexVectorView y{
        matpack::exhaustive_mdspan<std::complex<Numeric>, 1>{x.data(), 8}};

    const std::complex<Numeric> z = y * y;
    ARTS_USER_ERROR_IF(
        z != std::transform_reduce(
                 x.begin(), x.end(), x.begin(), std::complex<Numeric>{0}),
        z)
  }

  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexMatrixView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{
        x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{
        matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{
            copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_ASSERT(y == test)

    y = 2 * matpack::eigen::as_eigen(y) + y;

    ARTS_ASSERT(y != test)
    for (auto& a : copy) a = 2.0 * a + a;
    ARTS_ASSERT(y == test)
  }

  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexMatrixView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{
        x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{
        matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{
            copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_ASSERT(y == test)

    y = 2 * matpack::eigen::as_eigen(y) - y;

    ARTS_ASSERT(y == test)
  }
}

void test_data() {
  matpack::matpack_data<Numeric, 3> x(4, 2, 3, 2.1);
  std::cout << x << '\n' << '\n';
  std::cout << x[0] << '\n' << '\n';
  std::cout << x[0][0] << '\n' << '\n';
  std::cout << x[0][0][0] << '\n' << '\n';
  x[0][0][0] += 1;
  std::cout << x << '\n' << '\n';

  auto y = std::move(x).flatten();
  std::cout << y << '\n';
  std::cout << x << '\n';

  matpack::matpack_data<Numeric, 1> z;
  z.swap(y);
}

void test_complex() {
  {
    Complex x{0, 0};
    const Complex y{0, 0};
    ARTS_USER_ERROR_IF(x != y, x, ' ', y)
    ARTS_ASSERT(real_val(x) == real_val(y))
    ARTS_ASSERT(imag_val(x) == imag_val(y))
  }

  {
    Complex x{1, 1};
    x = 2. + x;
    ARTS_ASSERT(x == (Complex{3, 1}))
    x = 2. - x;
    ARTS_ASSERT(x == (Complex{-1, -1}))
    x = 2. * x;
    ARTS_ASSERT(x == (Complex{-2, -2}))
    x = 2. / x;
    ARTS_ASSERT(x == (Complex{-0.5, 0.5}))
  }

  {
    Complex x{1, 1};
    x = x + 2.;
    ARTS_ASSERT(x == (Complex{3, 1}))
    x = x - 2.;
    ARTS_ASSERT(x == (Complex{1, 1}))
    x = x * 2.;
    ARTS_ASSERT(x == (Complex{2, 2}))
    x = x / 2.;
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

  {
    Complex x{3, 3};
    const std::complex<int> y{3, 3};
    x = x + y;
    ARTS_ASSERT(x == (Complex{6, 6}))
    x = x - y;
    ARTS_ASSERT(x == (Complex{3, 3}))
    x = x * y;
    ARTS_ASSERT(x == (Complex{0, 18}))
    x = x / y;
    ARTS_ASSERT(x == (Complex{3, 3}))
  }
}

void test_math() {
  {
    Vector x{0, 1, 2, 3, 4};
    ARTS_ASSERT(min(x) == 0)
    ARTS_ASSERT(min(ExhaustiveVectorView{x}) == 0)
    ARTS_ASSERT(min(ExhaustiveConstVectorView{x}) == 0)
    ARTS_ASSERT(min(VectorView{x}) == 0)
    ARTS_ASSERT(min(ConstVectorView{x}) == 0)

    static_assert(std::random_access_iterator<decltype(x.elem_begin())>);
    static_assert(std::random_access_iterator<
                  decltype(ExhaustiveVectorView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<
                  decltype(ExhaustiveConstVectorView{x}.elem_begin())>);
    static_assert(
        std::random_access_iterator<decltype(VectorView{x}.elem_begin())>);
    static_assert(
        std::random_access_iterator<decltype(ConstVectorView{x}.elem_begin())>);

    static_assert(std::random_access_iterator<decltype(x.begin())>);
    static_assert(
        std::random_access_iterator<decltype(ExhaustiveVectorView{x}.begin())>);
    static_assert(std::random_access_iterator<
                  decltype(ExhaustiveConstVectorView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.begin())>);
    static_assert(
        std::random_access_iterator<decltype(ConstVectorView{x}.begin())>);

    static_assert(std::random_access_iterator<decltype(x.elem_end())>);
    static_assert(std::random_access_iterator<
                  decltype(ExhaustiveVectorView{x}.elem_end())>);
    static_assert(std::random_access_iterator<
                  decltype(ExhaustiveConstVectorView{x}.elem_end())>);
    static_assert(
        std::random_access_iterator<decltype(VectorView{x}.elem_end())>);
    static_assert(
        std::random_access_iterator<decltype(ConstVectorView{x}.elem_end())>);

    static_assert(std::random_access_iterator<decltype(x.end())>);
    static_assert(
        std::random_access_iterator<decltype(ExhaustiveVectorView{x}.end())>);
    static_assert(std::random_access_iterator<
                  decltype(ExhaustiveConstVectorView{x}.end())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.end())>);
    static_assert(
        std::random_access_iterator<decltype(ConstVectorView{x}.end())>);
  }

  {
    Matrix x(3, 3, 0);
    ARTS_ASSERT(min(x) == 0)
    ARTS_ASSERT(min(ExhaustiveMatrixView{x}) == 0)
    ARTS_ASSERT(min(ExhaustiveConstMatrixView{x}) == 0)

    static_assert(std::random_access_iterator<decltype(x.elem_begin())>);
    static_assert(std::random_access_iterator<
                  decltype(ExhaustiveMatrixView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<
                  decltype(ExhaustiveConstMatrixView{x}.elem_begin())>);
    static_assert(
        std::random_access_iterator<decltype(MatrixView{x}.elem_begin())>);
    static_assert(
        std::random_access_iterator<decltype(ConstMatrixView{x}.elem_begin())>);

    static_assert(std::random_access_iterator<decltype(x.begin())>);
    static_assert(
        std::random_access_iterator<decltype(ExhaustiveMatrixView{x}.begin())>);
    static_assert(std::random_access_iterator<
                  decltype(ExhaustiveConstMatrixView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(MatrixView{x}.begin())>);
    static_assert(
        std::random_access_iterator<decltype(ConstMatrixView{x}.begin())>);

    static_assert(std::random_access_iterator<decltype(x.elem_end())>);
    static_assert(std::random_access_iterator<
                  decltype(ExhaustiveMatrixView{x}.elem_end())>);
    static_assert(std::random_access_iterator<
                  decltype(ExhaustiveConstMatrixView{x}.elem_end())>);
    static_assert(
        std::random_access_iterator<decltype(MatrixView{x}.elem_end())>);
    static_assert(
        std::random_access_iterator<decltype(ConstMatrixView{x}.elem_end())>);

    static_assert(std::random_access_iterator<decltype(x.end())>);
    static_assert(
        std::random_access_iterator<decltype(ExhaustiveMatrixView{x}.end())>);
    static_assert(std::random_access_iterator<
                  decltype(ExhaustiveConstMatrixView{x}.end())>);
    static_assert(std::random_access_iterator<decltype(MatrixView{x}.end())>);
    static_assert(
        std::random_access_iterator<decltype(ConstMatrixView{x}.end())>);
  }
}

void test_mult() {
  {
    Matrix A(4, 4, 1);
    for (Index i = 0; i < 4; i++)
      for (Index j = 0; j < 4; j++) A(i, j) = static_cast<Numeric>(i);
    Vector y({1, 2, 3, 4}), x(4);

    Vector eig_x{A * y};
    mult(x, A, y);
    ARTS_USER_ERROR_IF(eig_x not_eq x, eig_x, " vs ", x)
  }
  {
    ComplexMatrix A(4, 4, 1);
    for (Index i = 0; i < 4; i++)
      for (Index j = 0; j < 4; j++)
        A(i, j) = {static_cast<Numeric>(i), 2 * static_cast<Numeric>(i)};
    ComplexVector y({1, 2, 3, 4}), x(4);

    ComplexVector eig_x{A * y};
    mult(x, A, y);
    ARTS_USER_ERROR_IF(eig_x not_eq x, eig_x, " vs ", x)
  }
}

void test_const_view() {
  std::array<double, 3 * 3> x{1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::array<double, 3> x3{1, 2, 3};

  {
    matpack::matpack_constant_view<double, false, 3, 3> xv(x);
    ARTS_USER_ERROR_IF(xv != xv, "Bad")
    ARTS_USER_ERROR_IF(xv(2, 2) != x.back(), "Bad")
    ARTS_USER_ERROR_IF(xv(0, 0) != x.front(), "Bad")
    ARTS_USER_ERROR_IF(
        xv[0] != (matpack::matpack_constant_view<double, false, 3>{x3}), "Bad")
    ARTS_USER_ERROR_IF(xv[0][0] != x.front(), "Bad")
  }

  {
    matpack::matpack_constant_view<double, true, 3, 3> xv(x);
    ARTS_USER_ERROR_IF(xv != xv, "Bad")
    ARTS_USER_ERROR_IF(xv(2, 2) != x.back(), "Bad")
    ARTS_USER_ERROR_IF(xv(0, 0) != x.front(), "Bad")
    ARTS_USER_ERROR_IF(
        xv[0] != (matpack::matpack_constant_view<double, true, 3>{x3}), "Bad")
    ARTS_USER_ERROR_IF(xv[0][0] != x.front(), "Bad")
  }
}

void test_const_data() {
  constexpr matpack::matpack_constant_data<Numeric, 3> x{1, 2, 3};
  Vector y(x);
  constexpr matpack::matpack_constant_data<Numeric, 3, 3> z{2, 3, 4};
  Matrix t(z);
}

void test_my_interp() {
  my_interp::Lagrange<1, true> x(
      0, 3.5, matpack::matpack_constant_data<Numeric, 7>{1, 2, 3, 4, 5, 6, 7});
  my_interp::Lagrange<1, false> x2(
      0, 3.5, matpack::matpack_constant_data<Numeric, 7>{1, 2, 3, 4, 5, 6, 7});
  auto iw = interpweights(std::array{x, x}, std::array{x2, x2});

  matpack::shape_help{iw.shape()};

  interp(Matrix(7, 7, 1), interpweights(x, x2), x, x2);

  reinterp(Matrix(7, 7, 1), iw, std::array{x, x}, std::array{x2, x2});

  Matrix Z = uniform_grid(-5, 49, 0.2).reshape(7, 7);
  Vector X = uniform_grid(1, 7, 1);
  Vector Y = uniform_grid(1, 7, 1);
  Vector XN = uniform_grid(1, 14, 0.5);
  Vector YN = uniform_grid(1, 28, 0.25);

  for (auto mx : XN) {
    for (auto my : YN) {
      auto lagx = my_interp::Lagrange<1, false>(0, mx, X);
      auto lagy = my_interp::Lagrange<1, false>(0, my, Y);
      auto iw2 = my_interp::interpweights(lagx, lagy);
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

void test_einsum_arr() {
  const Matrix Z = uniform_grid(2, 490000, 0.2).reshape(700, 700);

  {
    Matrix xm(700, 700, 1.0), ym(700, 700, 0.0);
    matpack::detail::einsum_arr<std::array{'i', 'j'},
                                std::array{'i', 'j'},
                                std::array{'j', 'i'},
                                std::array{'i', 'j'}>(xm, Z, Z, Z);

    for (Size i = 0; i < 700; i++) {
      for (Size j = 0; j < 700; j++) {
        ym[i][j] += Z(i, j) * Z(j, i) * Z(i, j);
      }
    }

    xm /= ym;
    xm -= 1.0;
    ARTS_USER_ERROR_IF(std::max(max(xm), std::abs(min(xm))) > 1e-10,
                       "Large discrepancy (main): ",
                       max(xm),
                       " or ",
                       min(xm));
  }

  {
    Matrix xm(700, 700, 1.0), ym(700, 700, 0.0);
    matpack::detail::einsum_arr<std::array{'j', 'i'},
                                std::array{'i', 'j'},
                                std::array{'j', 'i'},
                                std::array{'i', 'j'}>(xm, Z, Z, Z);

    for (Size i = 0; i < 700; i++) {
      for (Size j = 0; j < 700; j++) {
        ym[j][i] += Z(i, j) * Z(j, i) * Z(i, j);
      }
    }

    xm /= ym;
    xm -= 1.0;
    ARTS_USER_ERROR_IF(std::max(max(xm), std::abs(min(xm))) > 1e-10,
                       "Large discrepancy (transpose): ",
                       max(xm),
                       " or ",
                       min(xm));
  }

  {
    Vector xv(700, 1.0), yv(700, 0.0);
    matpack::detail::einsum_arr<std::array{'i'},
                                std::array{'i', 'j'},
                                std::array{'j', 'i'},
                                std::array{'i', 'j'}>(xv, Z, Z, Z);

    for (Size i = 0; i < 700; i++) {
      for (Size j = 0; j < 700; j++) {
        yv[i] += Z(i, j) * Z(j, i) * Z(i, j);
      }
    }

    xv /= yv;
    xv -= 1.0;
    ARTS_USER_ERROR_IF(std::max(max(xv), std::abs(min(xv))) > 1e-10,
                       "Large discrepancy (main): ",
                       max(xv),
                       " or ",
                       min(xv));
  }

  {
    Vector xv(700, 1.0), yv(700, 0.0);
    matpack::detail::einsum_arr<std::array{'j'},
                                std::array{'i', 'j'},
                                std::array{'j', 'i'},
                                std::array{'i', 'j'}>(xv, Z, Z, Z);

    for (Size i = 0; i < 700; i++) {
      for (Size j = 0; j < 700; j++) {
        yv[j] += Z(i, j) * Z(j, i) * Z(i, j);
      }
    }

    xv /= yv;
    xv -= 1.0;
    ARTS_USER_ERROR_IF(std::max(max(xv), std::abs(min(xv))) > 1e-10,
                       "Large discrepancy (transpose): ",
                       max(xv),
                       " or ",
                       min(xv));
  }

  {
    Numeric xn{0.0}, yn{0.0};
    matpack::detail::einsum_arr<std::array<char, 0>{},
                                std::array{'i', 'j'},
                                std::array{'j', 'i'},
                                std::array{'i', 'j'}>(xn, Z, Z, Z);

    for (Size i = 0; i < 700; i++) {
      for (Size j = 0; j < 700; j++) {
        yn += Z(i, j) * Z(j, i) * Z(i, j);
      }
    }

    xn /= yn;
    xn -= 1.0;
    ARTS_USER_ERROR_IF(std::abs(xn) > 1e-10, "Large discrepancy: ", xn);
  }

  {
    const Index L = 30, M = 20, N = 10, P = 25, C = 32;
    const Tensor3 LMN = uniform_grid(2, L * M * N, 0.2).reshape(L, M, N);
    const Tensor3 LMP = uniform_grid(2, L * M * P, 0.2).reshape(L, M, P);
    const Tensor4 LMNC =
        uniform_grid(2, L * M * N * C, 0.2).reshape(L, M, N, C);
    const Tensor4 LMPC =
        uniform_grid(2, L * M * P * C, 0.2).reshape(L, M, P, C);
    const Tensor5 LMNCP =
        uniform_grid(2, L * M * N * C * P, 0.2).reshape(L, M, N, C, P);
    Matrix x(L, M, 0.0), y(L, M, 0.0);

    matpack::detail::einsum_arr<std::array{'l', 'm'},
                                std::array{'l', 'm', 'n'},
                                std::array{'l', 'm', 'p'},
                                std::array{'l', 'm', 'n', 'c'},
                                std::array{'l', 'm', 'p', 'c'},
                                std::array{'l', 'm', 'n', 'c', 'p'}>(
        x, LMN, LMP, LMNC, LMPC, LMNCP);

    for (Index l = 0; l < L; l++) {
      for (Index m = 0; m < M; m++) {
        Numeric sum = 0.0;
        for (Index n = 0; n < N; n++) {
          for (Index p = 0; p < P; p++) {
            for (Index c = 0; c < C; c++) {
              sum += LMN(l, m, n) * LMP(l, m, p) * LMNC(l, m, n, c) *
                     LMPC(l, m, p, c) * LMNCP(l, m, n, c, p);
            }
          }
        }
        y(l, m) = sum;
      }
    }

    x /= y;
    x -= 1.0;
    ARTS_USER_ERROR_IF(std::max(max(x), std::abs(min(x))) > 1e-10,
                       "Large discrepancy (main): ",
                       max(x),
                       " or ",
                       min(x));
  }

  {
    const Index L = 30, M = 20, N = 10, P = 25, C = 32;
    const Tensor3 LMN = uniform_grid(2, L * M * N, 0.2).reshape(L, M, N);
    const Tensor3 LMP = uniform_grid(2, L * M * P, 0.2).reshape(L, M, P);
    const Tensor4 LMNC =
        uniform_grid(2, L * M * N * C, 0.2).reshape(L, M, N, C);
    const Tensor4 LMPC =
        uniform_grid(2, L * M * P * C, 0.2).reshape(L, M, P, C);
    const Tensor5 LMNCP =
        uniform_grid(2, L * M * N * C * P, 0.2).reshape(L, M, N, C, P);
    Matrix x(M, L, 0.0), y(M, L, 0.0), z(M, L, 0.0);

    matpack::detail::einsum_arr<std::array{'m', 'l'},
                                std::array{'l', 'm', 'n'},
                                std::array{'l', 'm', 'p'},
                                std::array{'l', 'm', 'n', 'c'},
                                std::array{'l', 'm', 'p', 'c'},
                                std::array{'l', 'm', 'n', 'c', 'p'}>(
        x, LMN, LMP, LMNC, LMPC, LMNCP);
    einsum<"ml", "lmn", "lmp", "lmnc", "lmpc", "lmncp">(
        z, LMN, LMP, LMNC, LMPC, LMNCP);

    for (Index l = 0; l < L; l++) {
      for (Index m = 0; m < M; m++) {
        Numeric sum = 0.0;
        for (Index n = 0; n < N; n++) {
          for (Index p = 0; p < P; p++) {
            for (Index c = 0; c < C; c++) {
              sum += LMN(l, m, n) * LMP(l, m, p) * LMNC(l, m, n, c) *
                     LMPC(l, m, p, c) * LMNCP(l, m, n, c, p);
            }
          }
        }
        y(m, l) = sum;
      }
    }

    z /= y;
    z -= 1.0;
    ARTS_USER_ERROR_IF(std::max(max(z), std::abs(min(z))) > 1e-10,
                       "Large discrepancy (transpose): ",
                       max(z),
                       " or ",
                       min(z));
    x /= y;
    x -= 1.0;
    ARTS_USER_ERROR_IF(std::max(max(x), std::abs(min(x))) > 1e-10,
                       "Large discrepancy (transpose): ",
                       max(x),
                       " or ",
                       min(x));
  }

  {
    Vector x{1, 2, 3};
    Matrix y(3, 3, 0.0);
    einsum<"ii", "i">(y, x);
    ARTS_USER_ERROR_IF(
        (Matrix{y}.reshape(9) != Vector{1., 0., 0., 0., 2., 0., 0., 0., 3.}),
        "Bad values:\n",
        y);
  }

  {
    Vector x{1, 2, 3};
    Matrix y(3, 3, 0.0);
    einsum<"ij", "i", "">(y, x, 2.0);
    ARTS_USER_ERROR_IF(
        (Matrix{y}.reshape(9) != Vector{2., 2., 2., 4., 4., 4., 6., 6., 6.}),
        "Bad values:\n",
        y);
  }

  {
    Vector x{1, 2, 3};
    const Vector y = einsum<Vector, "i", "i", "">({3}, x, -1);
    x += y;
    ARTS_USER_ERROR_IF(std::max(std::abs(max(x)), std::abs(min(x))) > 1e-10,
                       "Bad values:\n",
                       x);
  }

  {
    const Index N = 10'000;
    Matrix A(N, N, 1);
    A(0, 0) = 2;
    A(N - 1, N - 1) = -2;
    Vector x(N, 1), y(N, 0);

    y = 0;
    {
      DebugTime t{"vec loop 1"};
      for (Index j = 0; j < N; j++) {
        for (Index i = 0; i < N; i++) {
          y[i] += A(i, j) * x[j];
        }
      }
    }

    {
      DebugTime t{"vec mult"};
      mult(y, A, x);
    }

    {
      DebugTime t{"vec einsum"};
      einsum<"i", "ij", "j">(y, A, x);
    }

    y = 0;
    {
      DebugTime t{"vec loop 2"};
      for (Index i = 0; i < N; i++) {
        for (Index j = 0; j < N; j++) {
          y[i] += A(i, j) * x[j];
        }
      }
    }
  }

  {
    const Index N = 700;
    Matrix A(N, N), B(N, N, 1), C(N, N, 1);

    A = 0;
    {
      DebugTime t{"mat loop 1"};
      for (Index i = 0; i < N; i++) {
        for (Index k = 0; k < N; k++) {
          for (Index j = 0; j < N; j++) {
            A(i, k) += B(i, j) * C(j, k);
          }
        }
      }
    }

    {
      DebugTime t{"mat mult"};
      mult(A, B, C);
    }

    {
      DebugTime t{"mat einsum"};
      einsum<"ik", "ij", "jk">(A, B, C);
    }

    A = 0;
    {
      DebugTime t{"mat loop 2"};
      for (Index i = 0; i < N; i++) {
        for (Index j = 0; j < N; j++) {
          for (Index k = 0; k < N; k++) {
            A(i, k) += B(i, j) * C(j, k);
          }
        }
      }
    }
    {
      DebugTime t{"mat loop 3"};
      for (Index i = 0; i < N; i++) {
        for (Index j = 0; j < N; j++) {
          einsum<"k", "", "k">(A[i], B(i, j), C[j]);
        }
      }
    }
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

void test_solve2() {
  ComplexMatrix A(4, 4);
  A[0] = ComplexVector{1.8 + 1i, 2.88, 2.05, -0.89};
  A[1] = ComplexVector{5.25, -2.95, -0.95, -3.8};
  A[2] = ComplexVector{1.58, -2.69, -2.9, -1.04};
  A[3] = ComplexVector{-1.11, -0.66, -0.59, 0.8};
  ComplexVector b(4);
  ComplexMatrix C(4, 4);
  diagonalize(C, b, A);
  std::cout << "A:\n" << A << '\n';
  std::cout  << "b:\n" << b << '\n';
  std::cout  << "C:\n" << C << '\n';
}

#define EXECUTE_TEST(X)                                                       \
  std::cout << "#########################################################\n"; \
  std::cout << "Executing test: " #X << '\n';                                 \
  std::cout << "#########################################################\n"; \
  X();                                                                        \
  std::cout << "#########################################################\n";

int main() {
  // EXECUTE_TEST(test_view)
  // EXECUTE_TEST(test_eigen)
  // EXECUTE_TEST(test_data)
  // EXECUTE_TEST(test_complex)
  // EXECUTE_TEST(test_math)
  // EXECUTE_TEST(test_mult)
  // EXECUTE_TEST(test_const_view)
  // EXECUTE_TEST(test_const_data)
  // EXECUTE_TEST(test_my_interp)
  // EXECUTE_TEST(test_sorted_grid)
  // EXECUTE_TEST(test_einsum_arr)
  // EXECUTE_TEST(test_lapack_vector_mult)
  test_solve2();
}
