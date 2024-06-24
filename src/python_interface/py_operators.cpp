#include <operators.h>

#include "hpy_arts.h"
#include "nanobind/nanobind.h"
#include "python_interface.h"

namespace Python {
void py_operators(py::module_& m) {
  py::class_<NumericUnaryOperator> nuop(m,"NumericUnaryOperator");
      nuop.def(py::init_implicit<NumericUnaryOperator::func_t>())
      .def("__call__",
           &NumericUnaryOperator::operator(),
           py::arg("x"),
           "Apply the operator to a value");
           workspace_group_interface(nuop);

  py::class_<NumericTernaryOperator> ntop(m,"NumericTernaryOperator");
     ntop .def(py::init_implicit<NumericTernaryOperator::func_t>())
      .def("__call__",
           &NumericTernaryOperator::operator(),
           py::arg("x"),
           py::arg("y"),
           py::arg("z"),
           "Apply the operator to a value");
           workspace_group_interface(ntop);
}
}  // namespace Python