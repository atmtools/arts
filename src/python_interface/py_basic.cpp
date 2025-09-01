#include <nanobind/operators.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/variant.h>
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
  generic_interface(str);
  str.def("__init__", [](ValueHolder<String>* s, const py::bytes& b) {
    new (s) ValueHolder<String>{b.c_str()};
  });
  py::implicitly_convertible<py::bytes, ValueHolder<String>>();

  py::class_<ValueHolder<Numeric>> num(m, "Numeric");
  value_holder_interface(num);
  generic_interface(num);
  using nd_numeric = py::ndarray<py::numpy, Numeric, py::ndim<0>, py::c_contig>;
  num.def(
         "__array__",
         [](ValueHolder<Numeric>& n,
            py::object dtype,
            py::object copy) -> std::variant<nd_numeric, py::object> {
           std::array<size_t, 0> shape{};

           auto np = py::module_::import_("numpy");
           auto x  = nd_numeric(n.val.get(), 0, shape.data(), py::cast(n));

           if (not dtype.is_none()) {
             return np.attr("asarray")(x, "dtype"_a = dtype, "copy"_a = copy);
           }

           return x.cast((copy.is_none() or not py::bool_(copy))
                             ? py::rv_policy::automatic_reference
                             : py::rv_policy::copy);
         },
         "dtype"_a.none() = py::none(),
         "copy"_a.none()  = py::none())
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")(); },
          [](ValueHolder<Numeric>& a, const ValueHolder<Numeric>& b) { a = b; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`numpy.ndarray`");
  common_ndarray(num);

  py::class_<ValueHolder<Index>> ind(m, "Index");
  value_holder_interface(ind);
  generic_interface(ind);
  using nd_index = py::ndarray<py::numpy, Index, py::ndim<0>, py::c_contig>;
  ind.def(
         "__array__",
         [](ValueHolder<Index>& n,
            py::object dtype,
            py::object copy) -> std::variant<nd_index, py::object> {
           std::array<size_t, 0> shape{};

           auto np = py::module_::import_("numpy");
           auto x  = nd_index(n.val.get(), 0, shape.data(), py::cast(n));

           if (not dtype.is_none()) {
             return np.attr("asarray")(x, "dtype"_a = dtype, "copy"_a = copy);
           }

           return x.cast((copy.is_none() or not py::bool_(copy))
                             ? py::rv_policy::automatic_reference
                             : py::rv_policy::copy);
         },
         "dtype"_a.none() = py::none(),
         "copy"_a.none()  = py::none(),
         "Returns a :class:`~numpy.ndarray` of the object.")
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")(); },
          [](ValueHolder<Index>& a, const ValueHolder<Index>& b) { a = b; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`numpy.ndarray`");
  common_ndarray(ind);

  auto aos = py::class_<ArrayOfString>(m, "ArrayOfString");
  value_holder_vector_interface(aos);
  generic_interface(aos);

  auto aoi = py::class_<ArrayOfIndex>(m, "ArrayOfIndex");
  using nd_index_array =
      py::ndarray<py::numpy, Index, py::ndim<1>, py::c_contig>;
  aoi.def(
         "__array__",
         [](ArrayOfIndex& v,
            py::object dtype,
            py::object copy) -> std::variant<nd_index_array, py::object> {
           std::array<size_t, 1> shape{static_cast<size_t>(v.size())};

           auto np = py::module_::import_("numpy");
           auto x  = nd_index_array(v.data(), 1, shape.data(), py::cast(v));

           if (not dtype.is_none()) {
             return np.attr("asarray")(x, "dtype"_a = dtype, "copy"_a = copy);
           }

           return x.cast((copy.is_none() or not py::bool_(copy))
                             ? py::rv_policy::automatic_reference
                             : py::rv_policy::copy);
         },
         "dtype"_a.none() = py::none(),
         "copy"_a.none()  = py::none())
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")(); },
          [](ArrayOfIndex& a, const ArrayOfIndex& b) { a = b; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`numpy.ndarray`");
  common_ndarray(aoi);
  value_holder_vector_interface(aoi);
  generic_interface(aoi);

  auto aon = py::class_<ArrayOfNumeric>(m, "ArrayOfNumeric");
  using nd_numeric_array =
      py::ndarray<py::numpy, Numeric, py::ndim<1>, py::c_contig>;
  aon.def(
         "__array__",
         [](ArrayOfNumeric& v,
            py::object dtype,
            py::object copy) -> std::variant<nd_numeric_array, py::object> {
           std::array<size_t, 1> shape{static_cast<size_t>(v.size())};

           auto np = py::module_::import_("numpy");
           auto x  = nd_numeric_array(v.data(), 1, shape.data(), py::cast(v));

           if (not dtype.is_none()) {
             return np.attr("asarray")(x, "dtype"_a = dtype, "copy"_a = copy);
           }

           return x.cast((copy.is_none() or not py::bool_(copy))
                             ? py::rv_policy::automatic_reference
                             : py::rv_policy::copy);
         },
         "dtype"_a.none() = py::none(),
         "copy"_a.none()  = py::none(),
         "Returns a :class:`~numpy.ndarray` of the object.")
      .def_prop_rw(
          "value",
          [](py::object& x) { return x.attr("__array__")(); },
          [](ArrayOfNumeric& a, const ArrayOfNumeric& b) { a = b; },
          "A :class:`~numpy.ndarray` of the object.\n\n.. :class:`numpy.ndarray`")
      .doc() = "A list of :class:`~pyarts3.arts.Numeric`";
  common_ndarray(aon);
  value_holder_vector_interface(aon);
  generic_interface(aon);

  auto a1 =
      py::bind_vector<ArrayOfArrayOfIndex, py::rv_policy::reference_internal>(
          m, "ArrayOfArrayOfIndex");
  vector_interface(a1);
  generic_interface(a1);

  auto a2 =
      py::bind_vector<ArrayOfArrayOfString, py::rv_policy::reference_internal>(
          m, "ArrayOfArrayOfString");
  vector_interface(a2);
  generic_interface(a2);

  py::class_<Any> any(m, "Any");
  generic_interface(any);
  any.def(
         "__init__",
         [](Any* a, const py::arg&, const py::kwargs&) { new (a) Any{}; },
         "Anything is allowed.")
      .def("__getstate__", [](const Any&) { return std::tuple<>(); })
      .def("__setstate__", [](Any* a, const std::tuple<>&) { new (a) Any{}; });
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize basic\n{}", e.what()));
}
}  // namespace Python
