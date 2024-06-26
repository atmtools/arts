#include <nanobind/operators.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "mystring.h"
#include "nanobind/nanobind.h"
#include "python_interface_value_type.h"
#include "supergeneric.h"

namespace Python {
void py_basic(py::module_& m) try {
  py::class_<ValueHolder<String>> str(m, "String");
  value_holder_interface(str);
  workspace_group_interface(str);
  str.def("__init__", [](ValueHolder<String>* s, const py::bytes& b) {
    new (s) ValueHolder<String>{b.c_str()};
  });
  py::implicitly_convertible<py::bytes, ValueHolder<String>>();

  py::class_<ValueHolder<Numeric>> num(m, "Numeric");
  value_holder_interface(num);
  workspace_group_interface(num);
  num.def(
         "__array__",
         [](ValueHolder<Numeric>& n, py::object dtype, py::object copy) {
           using nd =
               py::ndarray<py::numpy, Numeric, py::ndim<0>, py::c_contig>;
           std::array<size_t, 0> shape{};

           auto np = py::module_::import_("numpy");
           auto x  = nd(n.val.get(), 0, shape.data(), py::handle());
           return np.attr("asarray")(
               x, py::arg("dtype") = dtype, py::arg("copy") = copy);
         },
         "dtype"_a.none() = py::none(),
         "copy"_a.none()  = py::none())
      .def_prop_rw(
          "value",
          [](py::object& v) {
            auto np = py::module_::import_("numpy");
            return np.attr("asarray")(v, py::arg("copy") = false);
          },
          [](ValueHolder<Numeric>& a, const ValueHolder<Numeric>& b) {
            a = b;
          });
  common_ndarray(num);

  py::class_<ValueHolder<Index>> ind(m, "Index");
  value_holder_interface(ind);
  workspace_group_interface(ind);
  ind.def(
         "__array__",
         [](ValueHolder<Index>& n, py::object dtype, py::object copy) {
           using nd = py::ndarray<py::numpy, Index, py::ndim<0>, py::c_contig>;
           std::array<size_t, 0> shape{};

           auto np = py::module_::import_("numpy");
           auto x  = nd(n.val.get(), 0, shape.data(), py::handle());
           return np.attr("asarray")(
               x, py::arg("dtype") = dtype, py::arg("copy") = copy);
         },
         "dtype"_a.none() = py::none(),
         "copy"_a.none()  = py::none())
      .def_prop_rw(
          "value",
          [](py::object& v) {
            auto np = py::module_::import_("numpy");
            return np.attr("asarray")(v, py::arg("copy") = false);
          },
          [](ValueHolder<Index>& a, const ValueHolder<Index>& b) { a = b; });
  common_ndarray(ind);

  auto aos = py::class_<ArrayOfString>(m, "ArrayOfString");
  value_holder_vector_interface(aos);
  workspace_group_interface(aos);

  auto aoi = py::class_<ArrayOfIndex>(m, "ArrayOfIndex");
  aoi.def(
         "__array__",
         [](ArrayOfIndex& v, py::object dtype, py::object copy) {
           using nd = py::ndarray<py::numpy, Index, py::ndim<1>, py::c_contig>;
           std::array<size_t, 1> shape{static_cast<size_t>(v.size())};

           auto np = py::module_::import_("numpy");
           auto x  = nd(v.data(), 1, shape.data(), py::handle());
           return np.attr("asarray")(
               x, py::arg("dtype") = dtype, py::arg("copy") = copy);
         },
         "dtype"_a.none() = py::none(),
         "copy"_a.none()  = py::none())
      .def_prop_rw(
          "value",
          [](ArrayOfIndex& v) {
            auto np = py::module_::import_("numpy");
            return np.attr("asarray")(v, py::arg("copy") = false);
          },
          [](ArrayOfIndex& a, const ArrayOfIndex& b) { a = b; });
  value_holder_vector_interface(aoi);
  common_ndarray(aoi);
  workspace_group_interface(aoi);

  auto aon = py::class_<ArrayOfNumeric>(m, "ArrayOfNumeric");
  aon.def(
         "__array__",
         [](ArrayOfNumeric& v, py::object dtype, py::object copy) {
           using nd =
               py::ndarray<py::numpy, Numeric, py::ndim<1>, py::c_contig>;
           std::array<size_t, 1> shape{static_cast<size_t>(v.size())};

           auto np = py::module_::import_("numpy");
           auto x  = nd(v.data(), 1, shape.data(), py::handle());
           return np.attr("asarray")(
               x, py::arg("dtype") = dtype, py::arg("copy") = copy);
         },
         "dtype"_a.none() = py::none(),
         "copy"_a.none()  = py::none())
      .def_prop_rw(
          "value",
          [](ArrayOfNumeric& v) {
            auto np = py::module_::import_("numpy");
            return np.attr("asarray")(v, py::arg("copy") = false);
          },
          [](ArrayOfNumeric& a, const ArrayOfNumeric& b) { a = b; });
  value_holder_vector_interface(aon);
  common_ndarray(aon);

  auto a1 =
      py::bind_vector<ArrayOfArrayOfIndex, py::rv_policy::reference_internal>(
          m, "ArrayOfArrayOfIndex");
  vector_interface(a1);
  workspace_group_interface(a1);

  auto a2 =
      py::bind_vector<ArrayOfArrayOfString, py::rv_policy::reference_internal>(
          m, "ArrayOfArrayOfString");
  vector_interface(a2);
  workspace_group_interface(a2);

  py::class_<Any> any(m, "Any");
  workspace_group_interface(any);
  any.def("__init__",
          [](Any* a, const py::args&, const py::kwargs&) { new (a) Any{}; })
      .def("__getstate__", [](const Any&) { return std::tuple<>(); })
      .def("__setstate__", [](Any* a, const std::tuple<>&) { new (a) Any{}; });
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize basic\n", e.what()));
}
}  // namespace Python
