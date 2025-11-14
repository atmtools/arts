#include <enumsSpectralRadianceUnitType.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
#include <operators.h>

#include "henyey_greenstein.h"
#include "hpy_arts.h"
#include "hpy_numpy.h"
namespace Python {
void py_operators(py::module_& m) {
  py::class_<NumericUnaryOperator> nuop(m, "NumericUnaryOperator");
  nuop.def("__init__",
           [](NumericUnaryOperator* op, NumericUnaryOperator::func_t f) {
             new (op) NumericUnaryOperator([f = std::move(f)](Numeric x) {
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
  generic_interface(nuop);
  py::implicitly_convertible<NumericUnaryOperator::func_t,
                             NumericUnaryOperator>();

  py::class_<ExtSSACallback> esop(m, "ExtSSACallback");
  esop.def("__init__",
           [](ExtSSACallback* op, ExtSSACallback::func_t f) {
             new (op) ExtSSACallback(
                 [f = std::move(f)](Numeric x, const AtmPoint& atm_point) {
                   py::gil_scoped_acquire gil{};
                   return f(x, atm_point);
                 });
           })
      .def(
          "__call__",
          [](ExtSSACallback& f, Numeric x, const AtmPoint& y) {
            return f.f(x, y);
          },
          "frequency"_a,
          "atm_point"_a)
      .doc() =
      "A callback to get [extinction, ssa] from a single scattering albedo and extinction field.";
  generic_interface(esop);
  py::implicitly_convertible<ExtSSACallback::func_t, ExtSSACallback>();

  py::class_<NumericBinaryOperator> nbop(m, "NumericBinaryOperator");
  nbop.def("__init__",
           [](NumericBinaryOperator* op, NumericBinaryOperator::func_t f) {
             new (op) NumericBinaryOperator(
                 [f = std::move(f)](Numeric x, Numeric y) {
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
  generic_interface(nbop);
  py::implicitly_convertible<NumericBinaryOperator::func_t,
                             NumericBinaryOperator>();

  py::class_<NumericTernaryOperator> ntop(m, "NumericTernaryOperator");
  ntop.def("__init__",
           [](NumericTernaryOperator* op, NumericTernaryOperator::func_t f) {
             new (op) NumericTernaryOperator(
                 [f = std::move(f)](Numeric x, Numeric y, Numeric z) {
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
  generic_interface(ntop);
  py::implicitly_convertible<NumericTernaryOperator::func_t,
                             NumericTernaryOperator>();

  py::class_<ScatteringGeneralSpectralTROFunc> sgstro(
      m, "ScatteringGeneralSpectralTROFunc");
  sgstro
      .def("__init__",
           [](ScatteringGeneralSpectralTROFunc* op,
              ScatteringGeneralSpectralTROFunc::func_t f) {
             new (op) ScatteringGeneralSpectralTROFunc(
                 [f = std::move(f)](
                     const AtmPoint& x, const Vector& y, Index z) {
                   py::gil_scoped_acquire gil{};
                   return f(x, y, z);
                 });
           })
      .def(
          "__call__",
          [](ScatteringGeneralSpectralTROFunc& f,
             py::object x,
             py::object y,
             py::object z) { return vectorize(f.f, x, y, z); },
          "x"_a,
          "y"_a,
          "z"_a);
  sgstro.doc() =
      "A callback to get the bulk scattering properties of a species.";
  generic_interface(sgstro);
  py::implicitly_convertible<ScatteringGeneralSpectralTROFunc::func_t,
                             ScatteringGeneralSpectralTROFunc>();

  py::class_<SpectralRadianceTransformOperator> srtop(
      m, "SpectralRadianceTransformOperator");
  srtop.def(py::init_implicit<const std::string_view&>())
      .def("__init__",
           [](SpectralRadianceTransformOperator* op,
              SpectralRadianceTransformOperator::Op::func_t f) {
             new (op) SpectralRadianceTransformOperator(
                 SpectralRadianceTransformOperator::Op(
                     [f = std::move(f)](StokvecVector& x,
                                        StokvecMatrix& y,
                                        const AscendingGrid& z,
                                        const PropagationPathPoint& c) {
                       py::gil_scoped_acquire gil{};
                       f(x, y, z, c);
                     }));
           })
      .def(
          "__call__",
          [](const SpectralRadianceTransformOperator& f,
             StokvecVector& x,
             StokvecMatrix& y,
             const AscendingGrid& z,
             const PropagationPathPoint& c) { return f(x, y, z, c); },
          "spectral_radiance"_a,
          "spectral_radiance_jacobian"_a,
          "freq_grid"_a,
          "ray_path_point"_a);
  generic_interface(srtop);
  py::implicitly_convertible<SpectralRadianceTransformOperator::Op::func_t,
                             SpectralRadianceTransformOperator>();
}
}  // namespace Python
