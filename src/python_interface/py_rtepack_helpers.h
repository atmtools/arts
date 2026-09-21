#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>
#include <rtepack.h>

#include <algorithm>

#include "hpy_numpy.h"

namespace Python {
template <typename T, Index M, size_t... N, typename Array>
void rtepack_array_from_ndarray(matpack::data_t<T, M> *y, const Array &x) {
  constexpr std::array<Size, sizeof...(N)> component_shape{N...};
  constexpr Size                           rank = M + sizeof...(N);
  const Size                               sz   = T{}.size();

  std::array<Index, rank> input_shape{};
  for (Size i = 0; i < rank; ++i) input_shape[i] = static_cast<Index>(x.shape(i));

  bool valid_component_shape = true;
  for (Size i = 0; i < component_shape.size(); ++i) {
    valid_component_shape &= static_cast<Size>(x.shape(M + i)) == component_shape[i];
  }
  if (x.size() != 0 and not valid_component_shape) {
    throw std::invalid_argument(
        std::format("Trailing dimensions of input must match the underlying rtepack type. "
                    "The shape is {:B,} and the required trailing shape is {:B,}",
                    input_shape,
                    component_shape));
  }

  std::array<Index, M> shape{};
  std::copy_n(input_shape.begin(), M, shape.begin());
  new (y) matpack::data_t<T, M>(shape);

  for (Size i = 0; i < y->size(); ++i) {
    std::array<typename T::value_type, (N * ...)> value{};
    std::copy_n(x.data() + i * sz, sz, value.begin());
    y->elem_begin()[i] = T{value};
  }
}

template <typename T, Index M, size_t... N> void rtepack_array(py::class_<matpack::data_t<T, M>> &c) {
  using U = T::value_type;
  static_assert(sizeof...(N) > 0);
  static_assert((N * ...) == T{}.size());

  constexpr Size rank = M + sizeof...(N);
  using nd            = py::ndarray<py::numpy, U, py::ndim<rank>, py::c_contig>;
  using const_nd      = py::ndarray<py::numpy, const U, py::ndim<rank>, py::c_contig>;

  c.def("__init__", [](matpack::data_t<T, M> *y, const nd &x) { rtepack_array_from_ndarray<T, M, N...>(y, x); }, "x"_a);
  c.def(
      "__init__",
      [](matpack::data_t<T, M> *y, const const_nd &x) { rtepack_array_from_ndarray<T, M, N...>(y, x); },
      "x"_a);
  py::implicitly_convertible<nd, matpack::data_t<T, M>>();
  py::implicitly_convertible<const_nd, matpack::data_t<T, M>>();

  c.def(
      "__init__",
      [](matpack::data_t<T, M> *y, const matpack::data_t<U, M> &x) {
        new (y) matpack::data_t<T, M>(x.shape());
        std::transform(x.elem_begin(), x.elem_end(), y->elem_begin(), [](const U &z) { return T(z); });
      },
      "x"_a);
  py::implicitly_convertible<matpack::data_t<U, M>, matpack::data_t<T, M>>();
  c.def(
      "__init__",
      [](matpack::data_t<T, M> *y, const matpack::data_t<U, M + sizeof...(N)> &x) {
        constexpr std::array<Size, sizeof...(N)> component_shape{N...};
        const Size                               sz = T{}.size();

        bool valid_component_shape = true;
        for (Size i = 0; i < component_shape.size(); ++i) {
          valid_component_shape &= static_cast<Size>(x.extent(M + i)) == component_shape[i];
        }
        if (x.size() != 0 and not valid_component_shape) {
          throw std::invalid_argument(
              std::format("Trailing dimensions of input must match the underlying rtepack type. "
                          "The shape is {:B,} and the required trailing shape is {:B,}",
                          x.shape(),
                          component_shape));
        }

        std::array<Index, M> shape{};
        for (Size i = 0; i < M; i++) shape[i] = x.extent(i);

        new (y) matpack::data_t<T, M>(shape);

        auto outview = y->view_as(y->size());
        auto inview  = x.view_as(y->size(), sz);
        for (Size i = 0; i < y->size(); i++) {
          for (Size j = 0; j < sz; j++) { outview[i][j] = inview[i, j]; }
        }
      },
      "x"_a);
  py::implicitly_convertible<matpack::data_t<U, M + sizeof...(N)>, matpack::data_t<T, M>>();

  c.def(
      "__array__",
      [](matpack::data_t<T, M> &v, py::object dtype, py::object copy) -> std::variant<nd, py::object> {
        constexpr auto        n = M + sizeof...(N);
        std::array<size_t, n> shape{};
        stdr::copy(v.shape(), shape.begin());
        stdr::copy(std::array{N...}, shape.begin() + M);
        auto np = py::module_::import_("numpy");
        auto x  = nd(v.data_handle(), n, shape.data(), py::cast(&v));

        if (not dtype.is_none()) { return np.attr("asarray")(x, "dtype"_a = dtype, "copy"_a = copy); }

        if (copy.is_none() or not py::bool_(copy)) { return x.cast(py::rv_policy::automatic_reference); }
        return x.cast(py::rv_policy::copy);
      },
      "dtype"_a.none() = py::none(),
      "copy"_a.none()  = py::none(),
      "Returns a :class:`~numpy.ndarray` of the object.  The last dimensions are fixed to the underlying rtepack type's size.");

  c.def_prop_rw(
      "value",
      [](py::object &x) { return x.attr("__array__")(); },
      [](matpack::data_t<T, M> &x, matpack::data_t<T, M> &y) { x = y; },
      "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`~numpy.ndarray`");

  common_ndarray(c);
}
}  // namespace Python
