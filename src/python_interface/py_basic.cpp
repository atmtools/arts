#include <nanobind/operators.h>
#include <nanobind/stl/string.h>
#include <python_interface.h>

#include <algorithm>
#include <iomanip>
#include <memory>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "mystring.h"
#include "python_interface_value_type.h"
#include "supergeneric.h"

namespace Python {
void py_basic(py::module_& m) try {
  py::class_<ValueHolder<String>> str(m, "String");
  value_holder_interface(str);

  py::class_<ValueHolder<Numeric>> num(m, "Numeric");
  value_holder_interface(num);

  py::class_<ValueHolder<Index>> ind(m, "Index");
  value_holder_interface(ind);

  auto aos = py::bind_vector<ArrayOfString>(m, "ArrayOfString");
  vector_interface(aos);

  auto aoi = py::bind_vector<ArrayOfIndex>(m, "ArrayOfIndex");
  vector_interface(aoi);
  common_ndarray(aoi);

  auto aon = py::bind_vector<ArrayOfNumeric>(m, "ArrayOfNumeric");
  vector_interface(aon);
  common_ndarray(aon);

  py::class_<Any>(m, "Any")
      .def("__init__",
           [](Any* a, const py::args&, const py::kwargs&) { new (a) Any{}; })
      .def("__getstate__", [](const Any&) { return std::tuple<>(); })
      .def("__setstate__", [](Any* a, const std::tuple<>&) { new (a) Any{}; });
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize basic\n", e.what()));
}
}  // namespace Python
