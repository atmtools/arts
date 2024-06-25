#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/vector.h>
#include <operators.h>

#include "hpy_arts.h"

namespace Python {
void py_operators(py::module_& m) {
  py::class_<NumericUnaryOperator> nuop(m, "NumericUnaryOperator");
  nuop.def(py::init_implicit<NumericUnaryOperator::func_t>())
      .def(
          "__call__",
          [](NumericUnaryOperator& f, Numeric x) { return f(x); },
          py::arg("x"))
      .def(
          "__call__",
          [](NumericUnaryOperator& f, py::ndarray<Numeric> x) {
            auto u = py::ndarray<Numeric>(x);
            std::transform(x.data(), x.data() + x.size(), u.data(), f);
            return u;
          },
          py::arg("x"));
  workspace_group_interface(nuop);

  py::class_<NumericTernaryOperator> ntop(m, "NumericTernaryOperator");
  ntop.def(py::init_implicit<NumericTernaryOperator::func_t>())
      .def(
          "__call__",
          [](NumericTernaryOperator& f, Numeric x, Numeric y, Numeric z) {
            return f(x, y, z);
          },
          py::arg("x"),
          py::arg("y"),
          py::arg("z"));
  workspace_group_interface(ntop);
}
}  // namespace Python