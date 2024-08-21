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
  nuop.def("__init__",
           [](NumericUnaryOperator* op, NumericUnaryOperator::func_t f) {
             new (op) NumericUnaryOperator([f](Numeric x) {
               py::gil_scoped_acquire gil{};
               return f(x);
             });
           })
      .def(
          "__call__",
          [](NumericUnaryOperator& f, py::object x) {
            return vectorize(f.f, x);
          },
          "x"_a);
  workspace_group_interface(nuop);
  py::implicitly_convertible<NumericUnaryOperator::func_t, NumericUnaryOperator>();

  py::class_<NumericBinaryOperator> nbop(m, "NumericBinaryOperator");
  nbop.def("__init__",
           [](NumericBinaryOperator* op, NumericBinaryOperator::func_t f) {
             new (op) NumericBinaryOperator([f](Numeric x, Numeric y) {
               py::gil_scoped_acquire gil{};
               return f(x, y);
             });
           })
      .def(
          "__call__",
          [](NumericBinaryOperator& f, py::object x, py::object y) {
            return vectorize(f.f, x, y);
          },
          "x"_a,
          "y"_a);
  workspace_group_interface(nbop);
  py::implicitly_convertible<NumericBinaryOperator::func_t, NumericBinaryOperator>();

  py::class_<NumericTernaryOperator> ntop(m, "NumericTernaryOperator");
  ntop.def("__init__",
           [](NumericTernaryOperator* op, NumericTernaryOperator::func_t f) {
             new (op)
                 NumericTernaryOperator([f](Numeric x, Numeric y, Numeric z) {
                   py::gil_scoped_acquire gil{};
                   return f(x, y, z);
                 });
           })
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
  py::implicitly_convertible<NumericTernaryOperator::func_t, NumericTernaryOperator>();
}
}  // namespace Python
