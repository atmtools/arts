#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/vector.h>
#include <workspace.h>

NB_MAKE_OPAQUE(
    std::vector<std::pair<LineShapeModelVariable, lbl::temperature::data>>);
NB_MAKE_OPAQUE(std::vector<lbl::line_shape::species_model>);
NB_MAKE_OPAQUE(std::vector<lbl::line>)
NB_MAKE_OPAQUE(Array<LagrangeInterpolation>);
NB_MAKE_OPAQUE(Array<AbsorptionSingleLine>);

namespace Python {
namespace py = nanobind;
using namespace py::literals;

template <typename T>
void vector_interface(py::class_<Array<T>>& c) {
  using Vec = Array<T>;

  c.def("__copy__", [](const Vec& v) { return Vec(v); });

  c.def("__deepcopy__", [](const Vec& v, const py::dict&) { return Vec(v); });

  c.def("__getstate__", [](const Vec& v) { return std::tuple<Vec>{v}; });

  c.def("__setstate__",
        [](Vec& v, const std::tuple<Vec>& x) { new (&v) Vec{std::get<0>(x)}; });

  py::implicitly_convertible<py::list, Vec>();
}
}  // namespace Python
