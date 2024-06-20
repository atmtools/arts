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
  using nd  = py::ndarray<py::numpy, T, py::ndim<1>, py::c_contig>;

  c.def("__copy__", [](const Vec& v) { return Vec(v); });

  c.def("__deepcopy__", [](const Vec& v, const py::dict&) { return Vec(v); });

  c.def("__getstate__", [](const Vec& v) { return std::tuple<Vec>{v}; });

  c.def("__setstate__",
        [](Vec& v, const std::tuple<Vec>& x) { new (&v) Vec{std::get<0>(x)}; });
  c.def(
      "__array__",
      [](Vec& v) {
        std::array<size_t, 1> shape{static_cast<size_t>(v.size())};

        return nd(v.data(), 1, shape.data(), py::handle());
      },
      py::rv_policy::reference_internal);

  c.def_prop_rw(
      "value",
      py::cpp_function(
          [](Vec& x) {
            py::object np = py::module_::import_("numpy");
            return np.attr("array")(x, py::arg("copy") = false);
          },
          py::keep_alive<0, 1>()),
      [](Vec& x, Vec& y) { x = y; },
      "Operate on type as if :class:`numpy.ndarray` type");

  py::implicitly_convertible<py::list, Vec>();
  py::implicitly_convertible<nd, Vec>();
}
}  // namespace Python
