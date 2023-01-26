#include <algorithm>
#include <complex>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

#include "debug.h"
#include "experimental/__p0009_bits/full_extent_t.hpp"
#include "experimental/__p0009_bits/layout_stride.hpp"
#include "logic.h"

#include "matpack_complex.h"
#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_eigen.h"
#include "matpack_iter.h"
#include "matpack_math.h"
#include "matpack_view.h"

void test_view() {
  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1,2,3,4,5, 6, 7, 8};
    ConstTensor3View y{matpack::exhaustive_mdspan<Numeric, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i=0; i<2; i++) {
      [[maybe_unused]] auto g = y(i, joker, matpack::matpack_strided_access{joker})(i, i);
      [[maybe_unused]] auto h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric>);
    }
    
    Numeric d=1.0;
    for (auto a: y) for (auto b: a) for (auto c: b) {ARTS_USER_ERROR_IF(c != d, c, ' ', d); d+=1.0;}
  }

  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1,2,3,4,5, 6, 7, 8};
    Tensor3View y{matpack::exhaustive_mdspan<Numeric, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i=0; i<2; i++) {
      [[maybe_unused]] auto& g = y(i, joker, matpack::matpack_strided_access{joker})(i, i);
      [[maybe_unused]] auto& h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric&>);
    }
    
    Numeric d=1.0;
    for (const auto& a: y) for (auto b: a) for (auto c: b) {ARTS_USER_ERROR_IF(c != d, c, ' ', d); d+=1.0;}
  }

  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1,2,3,4,5, 6, 7, 8};
    ExhaustiveConstTensor3View y{matpack::exhaustive_mdspan<Numeric, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i=0; i<2; i++) {
      [[maybe_unused]] auto g = y(i, joker, matpack::matpack_strided_access{joker})(i, i);
      [[maybe_unused]] auto h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric>);
    }
    
    Numeric d=1.0;
    for (auto a: y) for (auto b: a) for (auto c: b) {ARTS_USER_ERROR_IF(c != d, c, ' ', d); d+=1.0;}
  }

  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1,2,3,4,5, 6, 7, 8};
    ExhaustiveTensor3View y{matpack::exhaustive_mdspan<Numeric, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i=0; i<2; i++) {
      [[maybe_unused]] auto& g = y(i, joker, matpack::matpack_strided_access{joker})(i, i);
      [[maybe_unused]] auto& h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric&>);
    }
    
    Numeric d=1.0;
    for (const auto& a: y) for (auto b: a) for (auto c: b) {ARTS_USER_ERROR_IF(c != d, c, ' ', d); d+=1.0;}
  }

  //! Allow assignment between views
  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> xr{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> xi{0, 0, 0, 0, 0, 0, 0, 0};
    ExhaustiveComplexTensor3View y{matpack::exhaustive_mdspan<std::complex<Numeric>, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View yr{matpack::exhaustive_mdspan<Numeric, 3>{xr.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View yi{matpack::exhaustive_mdspan<Numeric, 3>{xi.data(), std::array{Index{2}, Index{2}, Index{2}}}};

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
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    std::vector<Numeric> xr{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> xi{0, 0, 0, 0, 0, 0, 0, 0};
    ExhaustiveComplexTensor3View y{matpack::exhaustive_mdspan<std::complex<Numeric>, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View yr{matpack::exhaustive_mdspan<Numeric, 3>{xr.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View yi{matpack::exhaustive_mdspan<Numeric, 3>{xi.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    auto z = y(joker, matpack::matpack_strided_access(0,-1,2), joker);

    std::vector<Numeric> reszr{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> reszi{1, 2, 0, 0, 5, 6, 0, 0};
    const ExhaustiveTensor3View zr{matpack::exhaustive_mdspan<Numeric, 3>{reszr.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View zi{matpack::exhaustive_mdspan<Numeric, 3>{reszi.data(), std::array{Index{2}, Index{2}, Index{2}}}};

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
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    ExhaustiveComplexTensor3View y{matpack::exhaustive_mdspan<std::complex<Numeric>, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    
    std::vector<std::complex<Numeric>> copy = x;
    const ExhaustiveComplexTensor3View test{matpack::exhaustive_mdspan<std::complex<Numeric>, 3>{copy.data(), std::array{Index{2}, Index{2}, Index{2}}}};

    ARTS_ASSERT(y == test)

    y += 3;
    for (auto& a: copy) a += 3;
    ARTS_ASSERT(y == test)

    y /= 3;
    for (auto& a: copy) a /= 3;
    ARTS_ASSERT(y == test)

    y.imag() -= 2;
    for (auto& a: copy) a -= std::complex<Numeric>(0, 2);
    ARTS_ASSERT(y == test)

    y -= 5;
    for (auto& a: copy) a -= 5;
    ARTS_ASSERT(y == test)

    y *= 5;
    for (auto& a: copy) a *= 5;
    ARTS_ASSERT(y == test)

    y = 42;
    for (auto& a: copy) a = 42;
    ARTS_ASSERT(y == test)
  }

  //! Test external assignment
  {
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    ComplexVectorView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 1>{x.data(), 8}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexVectorView test{matpack::exhaustive_mdspan<std::complex<Numeric>, 1>{copy.data(), 8}};

    for(auto& a: copy) a+=std::complex<Numeric>(3, 2);

    ARTS_ASSERT(y != test)
    
    y = copy;

    ARTS_ASSERT(y == test)
  }
}

void test_eigen() {
  {
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    ComplexMatrixView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_ASSERT(y == test)

    y = 2 * matpack::eigen::as_eigen(y);

    ARTS_ASSERT(y != test)
    for (auto& a: copy) a *= 2;
    ARTS_ASSERT(y == test)
  }

  {
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    ComplexVectorView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 1>{x.data(), 8}};

    const std::complex<Numeric> z = y * y;
    ARTS_USER_ERROR_IF(z != std::transform_reduce(x.begin(), x.end(), x.begin(), std::complex<Numeric>{0}), z)
  }

  {
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    ComplexMatrixView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_ASSERT(y == test)

    y = 2 * matpack::eigen::as_eigen(y) + y;

    ARTS_ASSERT(y != test)
    for (auto& a: copy) a = 2.0 * a + a;
    ARTS_ASSERT(y == test)
  }

  {
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    ComplexMatrixView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_ASSERT(y == test)

    y = 2 * matpack::eigen::as_eigen(y) - y;
    
    ARTS_ASSERT(y == test)
  }
}

void test_data() {
  matpack::matpack_data<Numeric, 3> x(4, 2, 3, 2.1);
  std::cout << x << '\n' << '\n';
  std::cout << x[0]  << '\n' << '\n';
  std::cout << x[0][0]  << '\n' << '\n';
  std::cout << x[0][0][0]  << '\n' << '\n';
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
    static_assert(std::random_access_iterator<decltype(ExhaustiveVectorView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstVectorView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(ConstVectorView{x}.elem_begin())>);

    static_assert(std::random_access_iterator<decltype(x.begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveVectorView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstVectorView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(ConstVectorView{x}.begin())>);

    static_assert(std::random_access_iterator<decltype(x.elem_end())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveVectorView{x}.elem_end())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstVectorView{x}.elem_end())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.elem_end())>);
    static_assert(std::random_access_iterator<decltype(ConstVectorView{x}.elem_end())>);

    static_assert(std::random_access_iterator<decltype(x.end())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveVectorView{x}.end())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstVectorView{x}.end())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.end())>);
    static_assert(std::random_access_iterator<decltype(ConstVectorView{x}.end())>);
  }

  {
    Matrix x(3, 3, 0);
    ARTS_ASSERT(min(x) == 0)
    ARTS_ASSERT(min(ExhaustiveMatrixView{x}) == 0)
    ARTS_ASSERT(min(ExhaustiveConstMatrixView{x}) == 0)

    static_assert(std::random_access_iterator<decltype(x.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveMatrixView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstMatrixView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(MatrixView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(ConstMatrixView{x}.elem_begin())>);

    static_assert(std::random_access_iterator<decltype(x.begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveMatrixView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstMatrixView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(MatrixView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(ConstMatrixView{x}.begin())>);

    static_assert(std::random_access_iterator<decltype(x.elem_end())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveMatrixView{x}.elem_end())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstMatrixView{x}.elem_end())>);
    static_assert(std::random_access_iterator<decltype(MatrixView{x}.elem_end())>);
    static_assert(std::random_access_iterator<decltype(ConstMatrixView{x}.elem_end())>);

    static_assert(std::random_access_iterator<decltype(x.end())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveMatrixView{x}.end())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstMatrixView{x}.end())>);
    static_assert(std::random_access_iterator<decltype(MatrixView{x}.end())>);
    static_assert(std::random_access_iterator<decltype(ConstMatrixView{x}.end())>);
  }
}

void test_mult() {
  {
    Matrix A(4,4,1);
    for (Index i=0; i<4; i++) for (Index j=0; j<4; j++) A(i, j) = static_cast<Numeric>(i);
    Vector y({1,2,3,4}), x(4);

    Vector eig_x=A*y;
    mult(x, A, y);
    ARTS_USER_ERROR_IF(eig_x not_eq x, eig_x, " vs ", x)
  }
  {
    ComplexMatrix A(4,4,1);
    for (Index i=0; i<4; i++) for (Index j=0; j<4; j++) A(i, j) = {static_cast<Numeric>(i), 2*static_cast<Numeric>(i)};
    ComplexVector y({1,2,3,4}), x(4);

    ComplexVector eig_x=A*y;
    mult(x, A, y);
    ARTS_USER_ERROR_IF(eig_x not_eq x, eig_x, " vs ", x)
  }
}

#define EXECUTE_TEST(X) \
std::cout << "#########################################################\n";\
std::cout << "Executing test: " #X << '\n'; \
std::cout << "#########################################################\n";\
X();\
std::cout << "#########################################################\n";

int main() {
  EXECUTE_TEST(test_view)
  EXECUTE_TEST(test_eigen)
  EXECUTE_TEST(test_data)
  EXECUTE_TEST(test_complex)
  EXECUTE_TEST(test_math)
  EXECUTE_TEST(test_mult)
}

