#pragma once

#include <python_interface.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "configtypes.h"
#include "py_macros.h"

namespace Python {
using Scalar = std::variant<Numeric, Index>;

using numpy_array = py::array;

template <typename T, Index... sz>
auto& fix_matpack_constant_data(
    py::class_<matpack::matpack_constant_data<T, sz...>>& r)
  requires(sizeof...(sz) > 0)
{
  using matpackT = matpack::matpack_constant_data<T, sz...>;
  using arrT = py::array_t<T, py::array::forcecast>;

  static constexpr std::size_t dim = sizeof...(sz);
  static constexpr std::size_t size = (sz * ...);
  static constexpr std::array<Index, dim> shape = matpackT::shape();

  r.def(py::init([](const arrT& x) -> std::shared_ptr<matpackT> {
          if (static_cast<std::size_t>(x.ndim()) != dim) {
            throw std::runtime_error(
                var_string("Bad rank: ", dim, " vs ", x.ndim()));
          }

          std::array<Index, dim> arr_shape;
          for (std::size_t i = 0; i < dim; i++) {
            arr_shape[i] = static_cast<Index>(x.shape(i));
          }

          if (arr_shape != shape) {
            throw std::runtime_error(
                var_string("Bad shape: ",
                           matpack::shape_help<dim>(shape),
                           " vs ",
                           matpack::shape_help<dim>(arr_shape)));
          }

          return std::make_shared<matpackT>(
              x.template cast<std::array<T, size>>());
        }),
        "Cast from numeric :class:`~numpy.array`")
      .def(py::init([](const py::list& x) {
             return py::cast<matpackT>(x.cast<arrT>());
           }),
           "Cast from :class:`list`")
      .def(py::init([](const py::array& x) {
             return py::cast<matpackT>(x.cast<arrT>());
           }),
           "Cast from any :class:`~numpy.array`")
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](matpackT& x) -> py::buffer_info {
        std::vector<ssize_t> shpe{static_cast<ssize_t>(sz)...};
        std::vector<ssize_t> strides(shpe.size(),
                                     static_cast<ssize_t>(sizeof(T)));
        for (std::size_t j = 0; j < sizeof...(sz) - 1; j++) {
          strides[j] *= std::reduce(
              shpe.begin() + j + 1, shpe.end(), 1, std::multiplies<>());
        }
        return py::buffer_info(x.data.data(),
                               sizeof(T),
                               py::format_descriptor<T>::format(),
                               sizeof...(sz),
                               std::move(shpe),
                               std::move(strides));
      })
      .def_prop_rw(
          "value",
          py::cpp_function(
              [](matpackT& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](matpackT& x, matpackT& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"));
  py::implicitly_convertible<arrT, matpackT>();
  py::implicitly_convertible<py::list, matpackT>();
  py::implicitly_convertible<py::array, matpackT>();

  return r;
}

template <typename T, Index... sz>
auto register_matpack_constant_data(py::module_& m, const char* const name)
  requires(sizeof...(sz) > 0)
{
  using matpackT = matpack::matpack_constant_data<T, sz...>;

  static constexpr std::size_t dim = sizeof...(sz);
  static constexpr std::array<Index, dim> shape = matpackT::shape();

  py::class_<matpackT> r(m, name, py::buffer_protocol());
  r.def(py::init([]() { return std::make_shared<matpackT>(); }),
        "Default constant data")
      .PythonInterfaceCopyValue(matpackT)
      .PythonInterfaceBasicRepresentation(matpackT);
  fix_matpack_constant_data(r);

  return r;
}

template <typename T>
void test_correct_size(const std::vector<T>& x) {
  if constexpr (not std::is_same_v<Scalar, T>) {
    // T is a vector!
    ARTS_USER_ERROR_IF(
        x.size() and std::any_of(x.begin() + 1,
                                 x.end(),
                                 [n = x.front().size()](auto& v) {
                                   return v.size() != n;
                                 }),
        "Bad size")
    for (auto& y : x) test_correct_size(y);
  }
}

template <Size N, typename T>
std::array<Index, N> shape(const numpy_array& x)
  requires(N > 0)
{
  const Size D = x.ndim();
  const Size S = N - D;

  ARTS_USER_ERROR_IF(
      N < D, "Cannot cast rank ", D, " input into rank ", N, " output")

  std::array<Index, N> out;
  std::ranges::fill(out | std::views::take(S), 1);
  for (Size i = S; i < N; i++) {
    out[i] = static_cast<Index>(x.shape(i - S));
  }

  return out;
}

template <Size N, typename T>
std::shared_ptr<matpack::matpack_data<T, N>> copy(const numpy_array& x) {
  auto out = std::make_shared<matpack::matpack_data<T, N>>(shape<N, T>(x));
  py::object np = py::module_::import("numpy");
  np.attr("copyto")(py::cast(out).attr("value"), x);
  return out;
}
}  // namespace Python