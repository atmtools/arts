#include <nanobind/operators.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "hpy_numpy.h"
#include "hpy_vector.h"
#include "mystring.h"
#include "python_interface_value_type.h"
#include "supergeneric.h"

namespace Python {
void py_basic(py::module_& m) try {
  py::class_<ValueHolder<String>> str(m, "String");
  value_holder_interface(str);
  workspace_group_interface(str);

  py::class_<ValueHolder<Numeric>> num(m, "Numeric");
  value_holder_interface(num);
  workspace_group_interface(num);

  py::class_<ValueHolder<Index>> ind(m, "Index");
  value_holder_interface(ind);
  workspace_group_interface(ind);

  auto aos = py::bind_vector<ArrayOfString>(m, "ArrayOfString");
  vector_interface(aos);
  workspace_group_interface(aos);

  auto aoi = py::bind_vector<ArrayOfIndex>(m, "ArrayOfIndex");
  vector_interface(aoi);
  common_ndarray(aoi);
  workspace_group_interface(aoi);

  auto aon = py::bind_vector<ArrayOfNumeric>(m, "ArrayOfNumeric");
  vector_interface(aon);
  common_ndarray(aon);

  auto a1 = py::bind_vector<ArrayOfArrayOfIndex>(m, "ArrayOfArrayOfIndex");
  workspace_group_interface(a1);

  auto a2 = py::bind_vector<ArrayOfArrayOfString>(m, "ArrayOfArrayOfString");
  workspace_group_interface(a2);

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
