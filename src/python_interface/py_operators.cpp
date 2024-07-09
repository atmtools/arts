#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/vector.h>
#include <operators.h>

#include "hpy_arts.h"
#include "hpy_numpy.h"

namespace Python {
void py_operators(py::module_& m) {
  py::class_<NumericUnaryOperator> nuop(m, "NumericUnaryOperator");
  nuop.def(py::init_implicit<NumericUnaryOperator::func_t>())
      .def(
          "__call__",
          [](NumericUnaryOperator& f, py::object x) {
            return vectorize(f.f, x);
          },
          "x"_a);
  workspace_group_interface(nuop);

  py::class_<NumericTernaryOperator> ntop(m, "NumericTernaryOperator");
  ntop.def(py::init_implicit<NumericTernaryOperator::func_t>())
      .def(
          "__call__",
          [](NumericTernaryOperator& f,
             py::object x,
             py::object y,
             py::object z) { return vectorize(f.f, x, y, z); },
          "x"_a,
          "y"_a,
          "z"_a);
  workspace_group_interface(ntop);
}
}  // namespace Python