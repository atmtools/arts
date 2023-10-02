#include <operators.h>

#include "python_interface.h"

namespace Python {
void py_operators(py::module_& m) {
  artsclass<NumericUnaryOperator>(m, "NumericUnaryOperator")
      .def(py::init<>())
      .def(py::init<NumericUnaryOperator::func_t>())
      .def("__call__",
           py::vectorize(&NumericUnaryOperator::operator()),
           py::arg("x"),
           "Apply the operator to a value");

  artsclass<NumericTernaryOperator>(m, "NumericTernaryOperator")
      .def(py::init<>())
      .def(py::init<NumericTernaryOperator::func_t>())
      .def("__call__",
           py::vectorize(&NumericTernaryOperator::operator()),
           py::arg("x"),
           py::arg("y"),
           py::arg("z"),
           "Apply the operator to a value");
}
}  // namespace Python