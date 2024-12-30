#include <debug.h>

#include <algorithm>
#include <complex>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "matpack.h"
#include "matpack_mdspan_helpers.h"
#include "matpack_mdspan_helpers_eigen.h"
#include "matpack_mdspan_helpers_vector.h"

namespace {
void test_view() {
  {
    Vector x{1, 2, 3, 4, 5, 6, 7, 8};

    ARTS_USER_ERROR_IF(not(x[0] != 0 and x[2] != 0 and x[4] != 0), "{:B}", x);
    ARTS_USER_ERROR_IF(not(x[1] != 0 and x[3] != 0 and x[5] != 0), "{:B}", x);

    VectorView v = x[Range(0, 3)];
    v            = 0;

    ARTS_USER_ERROR_IF(not(x[0] == 0 and x[2] == 0 and x[4] != 0), "{:B}", x);
    ARTS_USER_ERROR_IF(not(x[1] == 0 and x[3] != 0 and x[5] != 0), "{:B}", x);

    ARTS_USER_ERROR_IF(&x[2] != &v[2], "Pointing error (Vector and views)");
  }

  {
    Vector x{1, 2, 3, 4, 5, 6, 7, 8};

    ARTS_USER_ERROR_IF(not(x[0] != 0 and x[2] != 0 and x[4] != 0), "{:B}", x);
    ARTS_USER_ERROR_IF(not(x[1] != 0 and x[3] != 0 and x[5] != 0), "{:B}", x);

    StridedVectorView v = x[StridedRange(0, 3, 2)];
    v                   = 0;

    ARTS_USER_ERROR_IF(not(x[0] == 0 and x[2] == 0 and x[4] == 0), "{:B}", x);
    ARTS_USER_ERROR_IF(not(x[1] != 0 and x[3] != 0 and x[5] != 0), "{:B}", x);

    ARTS_USER_ERROR_IF(&x[2] != &v[1],
                       "Pointing error (Vector and strided views)");
  }

  {
    std::array<Index, 3> exts{5, 5, 5};
    const Size N = matpack::mdsize(exts);
    Tensor3 x    = matpack::uniform_grid(1, N, 1.0).reshape(exts);

    ARTS_USER_ERROR_IF(not(x[0] != 0 and x[2] != 0), "{:B}", x);
    ARTS_USER_ERROR_IF(not(x[1] != 0 and x[3] != 0), "{:B}", x);

    Tensor3View v = x[Range(0, 3)];
    v             = 0;

    ARTS_USER_ERROR_IF(not(x[0] == 0 and x[2] == 0), "{:B}", x);
    ARTS_USER_ERROR_IF(not(x[1] == 0 and x[3] != 0), "{:B}", x);

    ARTS_USER_ERROR_IF(x[2].base::data_handle() != v[2].base::data_handle(),
                       "Pointing error (Tensor3 and views)");
    ARTS_USER_ERROR_IF(x[2].base::data_handle() == x[1].base::data_handle(),
                       "Pointing error (Tensor3 and views)");
  }

  {
    std::array<Index, 3> exts{5, 5, 5};
    const Size N = matpack::mdsize(exts);
    Tensor3 x    = matpack::uniform_grid(1, N, 1.0).reshape(exts);

    ARTS_USER_ERROR_IF(not(x[0] != 0 and x[2] != 0), "{:B}", x);
    ARTS_USER_ERROR_IF(not(x[1] != 0 and x[3] != 0), "{:B}", x);

    StridedTensor3View v = x[StridedRange(0, 3, 2)];
    v                    = 0;

    ARTS_USER_ERROR_IF(not(x[0] == 0 and x[2] == 0), "{:B}", x);
    ARTS_USER_ERROR_IF(not(x[1] != 0 and x[3] != 0), "{:B}", x);

    ARTS_USER_ERROR_IF(x[2].base::data_handle() != v[1].base::data_handle(),
                       "Pointing error (Tensor3 and strided views)");
    ARTS_USER_ERROR_IF(x[2].base::data_handle() == x[1].base::data_handle(),
                       "Pointing error (Tensor3 and strided views)");
  }

  {
    std::array<Index, 3> exts{5, 5, 5};
    const Size N = matpack::mdsize(exts);
    Tensor3 xt   = matpack::uniform_grid(1, N, 1.0).reshape(exts);
    ComplexTensor3 x(xt.shape());
    stdr::transform(elemwise_range(xt), x.elem_begin(), [](auto a) { return Complex{a, a + 1.0}; });

    ARTS_USER_ERROR_IF(not(x[0] != 0 and x[2] != 0), "{:B}", x);
    ARTS_USER_ERROR_IF(not(x[1] != 0 and x[3] != 0), "{:B}", x);

    ComplexTensor3View v = x[Range(0, 2)];
    v                    = 0;

    ARTS_USER_ERROR_IF(not(x[0] == 0 and x[2] != 0), "{:B}", x);
    ARTS_USER_ERROR_IF(not(x[1] == 0 and x[3] != 0), "{:B}", x);

    ARTS_USER_ERROR_IF(x[2].base::data_handle() != v[2].base::data_handle(),
                       "Pointing error (ComplexTensor3 and views)");
    ARTS_USER_ERROR_IF(x[2].base::data_handle() == x[1].base::data_handle(),
                       "Pointing error (ComplexTensor3 and views)");
  }

  {
    std::array<Index, 3> exts{5, 5, 5};
    const Size N = matpack::mdsize(exts);
    Tensor3 xt   = matpack::uniform_grid(1, N, 1.0).reshape(exts);
    ComplexTensor3 x(xt.shape());
    stdr::transform(elemwise_range(xt), x.elem_begin(), [](auto a) { return Complex{a, a + 1.0}; });

    ARTS_USER_ERROR_IF(not(x[0] != 0 and x[2] != 0), "{:B}", x);
    ARTS_USER_ERROR_IF(not(x[1] != 0 and x[3] != 0), "{:B}", x);

    StridedComplexTensor3View v = x[StridedRange(0, 2, 2)];
    v                           = 0;

    ARTS_USER_ERROR_IF(not(x[0] == 0 and x[2] == 0), "{:B}", x);
    ARTS_USER_ERROR_IF(not(x[1] != 0 and x[3] != 0), "{:B}", x);

    ARTS_USER_ERROR_IF(x[2].base::data_handle() != v[1].base::data_handle(),
                       "Pointing error (ComplexTensor3 and views)");
    ARTS_USER_ERROR_IF(x[2].base::data_handle() == x[1].base::data_handle(),
                       "Pointing error (ComplexTensor3 and views)");
  }

  {
    std::array<Index, 3> exts{5, 5, 5};
    const Size N = matpack::mdsize(exts);
    ComplexTensor3 x(matpack::uniform_grid(1, N, 1.0).reshape(exts),
                     matpack::uniform_grid(2, N, 1.0).reshape(exts));
    StridedTensor3View xr{x.real()};
    StridedTensor3View xi{x.imag()};

    ARTS_USER_ERROR_IF(not(xr[0] != 0 and xr[2] != 0), "real {:B}", x);
    ARTS_USER_ERROR_IF(not(xr[1] != 0 and xr[3] != 0), "real {:B}", x);
    ARTS_USER_ERROR_IF(not(xi[0] != 0 and xi[2] != 0), "imag {:B}", x);
    ARTS_USER_ERROR_IF(not(xi[1] != 0 and xi[3] != 0), "imag {:B}", x);

    xr[StridedRange(0, 2, 2)] = 0;
    xi[StridedRange(1, 2, 2)] = 0;

    ARTS_USER_ERROR_IF(not(xr[0] == 0 and xr[2] == 0), "real {:B}", x);
    ARTS_USER_ERROR_IF(not(xr[1] != 0 and xr[3] != 0), "real {:B}", x);
    ARTS_USER_ERROR_IF(not(xi[0] != 0 and xi[2] != 0), "imag {:B}", x);
    ARTS_USER_ERROR_IF(not(xi[1] == 0 and xi[3] == 0), "imag {:B}", x);

    ARTS_USER_ERROR_IF(xr[2].base::data_handle() !=
                           xr[StridedRange(0, 3, 2)][1].base::data_handle(),
                       "Pointing error (ComplexTensor3 real views)");

    ARTS_USER_ERROR_IF(xi[2].base::data_handle() !=
                           xi[StridedRange(0, 3, 2)][1].base::data_handle(),
                       "Pointing error (ComplexTensor3 imag views)");
  }

  for (Index RIJ = 0; RIJ < 5; RIJ++) {
    for (Index RIK = 0; RIK < 5; RIK++) {
      for (Index IIJ = 0; IIJ < 5; IIJ++) {
        for (Index IIK = 0; IIK < 5; IIK++) {
          std::array<Index, 3> exts{5, 5, 5};
          const Size N = matpack::mdsize(exts);
          Tensor3 xt   = matpack::uniform_grid(1, N, 1.0).reshape(exts);
          ComplexTensor3 x(xt.shape());
          stdr::transform(elemwise_range(xt), x.elem_begin(), [](auto a) { return Complex{a, a + 1.0}; });

          StridedTensor3View xr{x.real()};
          StridedTensor3View xi{x.imag()};
          ComplexTensor3 y(xr, xi);

          xr[joker, RIJ, RIK] = 0;
          xi[joker, IIJ, IIK] = 0;

          for (Index i = 0; i < 5; i++) {
            for (Index j = 0; j < 5; j++) {
              for (Index k = 0; k < 5; k++) {
                const bool r_zero = (j == RIJ and k == RIK);
                const bool i_zero = (j == IIJ and k == IIK);

                Complex c{r_zero ? 0 : y[i, j, k].real(),
                          i_zero ? 0 : y[i, j, k].imag()};
                ARTS_USER_ERROR_IF(
                    (x[i, j, k] != c),
                    "Failed at ({0}, {1}, {2}). x[i,j,k] != c. x[{0}, {1}, {2}] := {3}, c := {4}",
                    i,
                    j,
                    k,
                    x[i, j, k],
                    c);
              }
            }
          }
        }
      }
    }
  }

  for (Index RII = 0; RII < 5; RII++) {
    for (Index RIK = 0; RIK < 5; RIK++) {
      for (Index III = 0; III < 5; III++) {
        for (Index IIK = 0; IIK < 5; IIK++) {
          std::array<Index, 3> exts{5, 5, 5};
          const Size N = matpack::mdsize(exts);
          Tensor3 xt   = matpack::uniform_grid(1, N, 1.0).reshape(exts);
          ComplexTensor3 x(xt.shape());
          stdr::transform(elemwise_range(xt), x.elem_begin(), [](auto a) { return Complex{a, a + 1.0}; });

          auto xr{x.real()};
          auto xi{x.imag()};
          ComplexTensor3 y(xr, xi);

          xr[RII, joker, RIK] = 0;
          xi[III, joker, IIK] = 0;

          for (Index i = 0; i < 5; i++) {
            for (Index j = 0; j < 5; j++) {
              for (Index k = 0; k < 5; k++) {
                const bool r_zero = (i == RII and k == RIK);
                const bool i_zero = (i == III and k == IIK);

                Complex c{r_zero ? 0 : y[i, j, k].real(),
                          i_zero ? 0 : y[i, j, k].imag()};
                ARTS_USER_ERROR_IF(
                    (x[i, j, k] != c),
                    "Failed at ({0}, {1}, {2}). x[i,j,k] != c. x[{0}, {1}, {2}] := {3}, c := {4}",
                    i,
                    j,
                    k,
                    x[i, j, k],
                    c);
              }
            }
          }
        }
      }
    }
  }

  for (Index RII = 0; RII < 5; RII++) {
    for (Index RIJ = 0; RIJ < 5; RIJ++) {
      for (Index III = 0; III < 5; III++) {
        for (Index IIJ = 0; IIJ < 5; IIJ++) {
          std::array<Index, 3> exts{5, 5, 5};
          const Size N = matpack::mdsize(exts);
          Tensor3 xt   = matpack::uniform_grid(1, N, 1.0).reshape(exts);
          ComplexTensor3 x(xt.shape());
          stdr::transform(elemwise_range(xt), x.elem_begin(), [](auto a) { return Complex{a, a + 1.0}; });

          auto xr{x.real()};
          auto xi{x.imag()};
          ComplexTensor3 y(xr, xi);

          xr[RII, RIJ, joker] = 0;
          xi[III, IIJ, joker] = 0;

          for (Index i = 0; i < 5; i++) {
            for (Index j = 0; j < 5; j++) {
              for (Index k = 0; k < 5; k++) {
                const bool r_zero = (i == RII and j == RIJ);
                const bool i_zero = (i == III and j == IIJ);

                Complex c{r_zero ? 0 : y[i, j, k].real(),
                          i_zero ? 0 : y[i, j, k].imag()};
                ARTS_USER_ERROR_IF(
                    (x[i, j, k] != c),
                    "Failed at ({0}, {1}, {2}). x[i,j,k] != c. x[{0}, {1}, {2}] := {3}, c := {4}",
                    i,
                    j,
                    k,
                    x[i, j, k],
                    c);
              }
            }
          }
        }
      }
    }
  }

  for (Index RII = 0; RII < 5; RII++) {
    for (Index III = 0; III < 5; III++) {
      std::array<Index, 3> exts{5, 5, 5};
      const Size N = matpack::mdsize(exts);
      Tensor3 xt   = matpack::uniform_grid(1, N, 1.0).reshape(exts);
      ComplexTensor3 x(xt.shape());
      stdr::transform(elemwise_range(xt), x.elem_begin(), [](auto a) { return Complex{a, a + 1.0}; });

      auto xr{x.real()};
      auto xi{x.imag()};
      ComplexTensor3 y(xr, xi);

      xr[RII, joker, joker] = 0;
      xi[III, joker, joker] = 0;

      for (Index i = 0; i < 5; i++) {
        for (Index j = 0; j < 5; j++) {
          for (Index k = 0; k < 5; k++) {
            const bool r_zero = (i == RII);
            const bool i_zero = (i == III);

            Complex c{r_zero ? 0 : y[i, j, k].real(),
                      i_zero ? 0 : y[i, j, k].imag()};
            ARTS_USER_ERROR_IF(
                (x[i, j, k] != c),
                "Failed at ({0}, {1}, {2}). x[i,j,k] != c. x[{0}, {1}, {2}] := {3}, c := {4}",
                i,
                j,
                k,
                x[i, j, k],
                c);
          }
        }
      }
    }
  }

  for (Index RIJ = 0; RIJ < 5; RIJ++) {
    for (Index IIJ = 0; IIJ < 5; IIJ++) {
      std::array<Index, 3> exts{5, 5, 5};
      const Size N = matpack::mdsize(exts);
      Tensor3 xt   = matpack::uniform_grid(1, N, 1.0).reshape(exts);
      ComplexTensor3 x(xt.shape());
      stdr::transform(elemwise_range(xt), x.elem_begin(), [](auto a) { return Complex{a, a + 1.0}; });

      auto xr{x.real()};
      auto xi{x.imag()};
      ComplexTensor3 y(xr, xi);

      xr[joker, RIJ, joker] = 0;
      xi[joker, IIJ, joker] = 0;

      for (Index i = 0; i < 5; i++) {
        for (Index j = 0; j < 5; j++) {
          for (Index k = 0; k < 5; k++) {
            const bool r_zero = (j == RIJ);
            const bool i_zero = (j == IIJ);

            Complex c{r_zero ? 0 : y[i, j, k].real(),
                      i_zero ? 0 : y[i, j, k].imag()};
            ARTS_USER_ERROR_IF(
                (x[i, j, k] != c),
                "Failed at ({0}, {1}, {2}). x[i,j,k] != c. x[{0}, {1}, {2}] := {3}, c := {4}",
                i,
                j,
                k,
                x[i, j, k],
                c);
          }
        }
      }
    }
  }

  for (Index RIK = 0; RIK < 5; RIK++) {
    for (Index IIK = 0; IIK < 5; IIK++) {
      std::array<Index, 3> exts{5, 5, 5};
      const Size N = matpack::mdsize(exts);
      Tensor3 xt   = matpack::uniform_grid(1, N, 1.0).reshape(exts);
      ComplexTensor3 x(xt.shape());
      stdr::transform(elemwise_range(xt), x.elem_begin(), [](auto a) { return Complex{a, a + 1.0}; });

      auto xr{x.real()};
      auto xi{x.imag()};
      ComplexTensor3 y(xr, xi);

      xr[joker, joker, RIK] = 0;
      xi[joker, joker, IIK] = 0;

      for (Index i = 0; i < 5; i++) {
        for (Index j = 0; j < 5; j++) {
          for (Index k = 0; k < 5; k++) {
            const bool r_zero = (k == RIK);
            const bool i_zero = (k == IIK);

            Complex c{r_zero ? 0 : y[i, j, k].real(),
                      i_zero ? 0 : y[i, j, k].imag()};
            ARTS_USER_ERROR_IF(
                (x[i, j, k] != c),
                "Failed at ({0}, {1}, {2}). x[i,j,k] != c. x[{0}, {1}, {2}] := {3}, c := {4}",
                i,
                j,
                k,
                x[i, j, k],
                c);
          }
        }
      }
    }
  }
}

void test_eigen() {
  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexMatrixView y{{x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{{copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_USER_ERROR_IF(not(y == test), "Error")

    y = 2 * y;

    ARTS_USER_ERROR_IF(not(y != test), "Error")
    for (auto& a : copy) a *= 2;
    ARTS_USER_ERROR_IF(not(y == test), "Error")

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

    ARTS_USER_ERROR_IF(not(y == test), "Error")

    y = 2 * y + y;

    ARTS_USER_ERROR_IF(not(y != test), "Error")
    for (auto& a : copy) a = 2.0 * a + a;
    ARTS_USER_ERROR_IF(not(y == test), "Error")
  }

  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    ComplexMatrixView y{{x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{{copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_USER_ERROR_IF(not(y == test), "Error")

    y = 2 * y - y;

    ARTS_USER_ERROR_IF(not(y == test), "Error")
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
    ARTS_USER_ERROR_IF(not(real_val(x) == real_val(y)), "Error")
    ARTS_USER_ERROR_IF(not(imag_val(x) == imag_val(y)), "Error")
  }

  {
    Complex x{1, 1};
    x = 2 + x;
    ARTS_USER_ERROR_IF(not(x == (Complex{3, 1})), "Error")
    x = 2 - x;
    ARTS_USER_ERROR_IF(not(x == (Complex{-1, -1})), "Error")
    x = 2 * x;
    ARTS_USER_ERROR_IF(not(x == (Complex{-2, -2})), "Error")
    x = 2 / x;
    ARTS_USER_ERROR_IF(not(x == (Complex{-0.5, 0.5})), "Error")
  }

  {
    Complex x{1, 1};
    x = x + 2;
    ARTS_USER_ERROR_IF(not(x == (Complex{3, 1})), "Error")
    x = x - 2;
    ARTS_USER_ERROR_IF(not(x == (Complex{1, 1})), "Error")
    x = x * 2;
    ARTS_USER_ERROR_IF(not(x == (Complex{2, 2})), "Error")
    x = x / 2;
    ARTS_USER_ERROR_IF(not(x == (Complex{1, 1})), "Error")
  }

  {
    Complex x{1, 1};
    x = 2 + x;
    ARTS_USER_ERROR_IF(not(x == (Complex{3, 1})), "Error")
    x = 2 - x;
    ARTS_USER_ERROR_IF(not(x == (Complex{-1, -1})), "Error")
    x = 2 * x;
    ARTS_USER_ERROR_IF(not(x == (Complex{-2, -2})), "Error")
    x = 2 / x;
    ARTS_USER_ERROR_IF(not(x == (Complex{-0.5, 0.5})), "Error")
  }

  {
    Complex x{1, 1};
    x = x + 2;
    ARTS_USER_ERROR_IF(not(x == (Complex{3, 1})), "Error")
    x = x - 2;
    ARTS_USER_ERROR_IF(not(x == (Complex{1, 1})), "Error")
    x = x * 2;
    ARTS_USER_ERROR_IF(not(x == (Complex{2, 2})), "Error")
    x = x / 2;
    ARTS_USER_ERROR_IF(not(x == (Complex{1, 1})), "Error")
  }
}

void test_math() {
  {
    Vector x{0, 1, 2, 3, 4};

    ARTS_USER_ERROR_IF(not(max(x) == 4), "Error")
    ARTS_USER_ERROR_IF(not(max(VectorView{x}) == 4), "Error")
    ARTS_USER_ERROR_IF(not(max(ConstVectorView{x}) == 4), "Error")
    ARTS_USER_ERROR_IF(not(max(StridedVectorView{x}) == 4), "Error")
    ARTS_USER_ERROR_IF(not(max(StridedConstVectorView{x}) == 4), "Error")

    ARTS_USER_ERROR_IF(not(min(x) == 0), "min(x) == {}\n {:B,}", min(x), x)
    ARTS_USER_ERROR_IF(not(min(VectorView{x}) == 0), "Error")
    ARTS_USER_ERROR_IF(not(min(ConstVectorView{x}) == 0), "Error")
    ARTS_USER_ERROR_IF(not(min(StridedVectorView{x}) == 0), "Error")
    ARTS_USER_ERROR_IF(not(min(StridedConstVectorView{x}) == 0), "Error")

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
    ARTS_USER_ERROR_IF(not(min(x) == 0), "Error")
    ARTS_USER_ERROR_IF(not(min(MatrixView{x}) == 0), "Error")
    ARTS_USER_ERROR_IF(not(min(ConstMatrixView{x}) == 0), "Error")

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

void test_einsum() {
  const Size n = 5, m = 3;
  const Matrix A = matpack::uniform_grid(1.0, n * m, 1.0).reshape(n, m);
  const Vector x = matpack::uniform_grid(1.0, m, 1.0);

  // mat-vec mult
  {
    Vector y0(n), y1(n);

    einsum<"i", "ij", "j">(y0, A, x);
    mult(y1, A, x);
    ARTS_USER_ERROR_IF(y0 != y1, "{:B,} != {:B,}", y0, y1);
  }

  // mat-vec collaps
  {
    Vector y0(m, 0), y1(m, 0);
    einsum<"j", "ij", "j">(y0, A, x);

    for (Size i = 0; i < n; i++) {
      for (Size j = 0; j < m; j++) {
        y1[j] += A[i, j] * x[j];
      }
    }

    ARTS_USER_ERROR_IF(y0 != y1, "{:B,} != {:B,}", y0, y1);
  }
}
}  // namespace

#define EXECUTE_TEST(X)                                                       \
  std::cout << "#########################################################\n"; \
  std::cout << "Executing test: " #X << '\n';                                 \
  std::cout << "#########################################################\n"; \
  X();                                                                        \
  std::cout << "#########################################################\n";

int main() try {
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
  EXECUTE_TEST(test_einsum)

  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::print(
      std::cerr, "ERROR: See below\n\n{}\n\nERROR: See above\n", e.what());
  return EXIT_FAILURE;
}
