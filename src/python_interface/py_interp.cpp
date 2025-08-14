#include <interp.h>
#include <lagrange_interp.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>

#include <variant>

namespace Python {
using PythonLags = std::variant<LagrangeInterpolation,
                                LagrangeCyclic0to360Interpolation,
                                LagrangeCyclicPM180Interpolation>;

void py_interp(py::module_& m) {
  auto interp = m.def_submodule("interp", "Interpolation methods");

  py::class_<my_interp::Empty>(interp, "Empty")
      .def(py::init<>())
      .def("__repr__", [](const my_interp::Empty&) { return "Empty()"; })
      .doc() =
      "An empty interpolation object, it is just a filler and has no use but naming itself";

  py::class_<LagrangeInterpolation>(interp, "PolyGrid")
      .def(
          "__init__",
          [](LagrangeInterpolation* l,
             const Numeric x,
             const Vector& xi,
             const Index polyorder) {
            new (l) LagrangeInterpolation(0, x, xi, polyorder);
          },
          "val"_a,
          "vec"_a,
          "polyorder"_a = Index{1},
          "Construct a Lagrange interpolation object")
      .def_prop_ro(
          "order",
          [](const LagrangeInterpolation& l) { return l.size() - 1; },
          "The order of interpolation")
      .def_rw("pos", &LagrangeInterpolation::pos, "The interpolation positions")
      .def_rw(
          "weights", &LagrangeInterpolation::lx, "The interpolation weights")
      .def("__getstate__",
           [](const LagrangeInterpolation& self) {
             return std::make_tuple(self.pos, self.lx, self.dlx);
           })
      .def("__setstate__",
           [](LagrangeInterpolation* l,
              const std::tuple<Index, Array<Numeric>, Array<Numeric>>& state) {
             new (l) LagrangeInterpolation{};

             l->pos = std::get<0>(state);
             l->lx  = std::get<1>(state);
             l->dlx = std::get<2>(state);
           })
      .doc() = "Polynomial interpolation object";

  py::bind_vector<ArrayOfLagrangeInterpolation,
                  py::rv_policy::reference_internal>(
      interp, "ArrayOfLagrangeInterpolation")
      .doc() = "List of :class:`~pyarts.arts.LagrangeInterpolation`";

  py::class_<LagrangeCyclic0to360Interpolation>(interp, "CyclicGrid0to360")
      .def(
          "__init__",
          [](LagrangeCyclic0to360Interpolation* l,
             const Numeric x,
             const Vector& xi,
             const Index polyorder) {
            new (l) LagrangeCyclic0to360Interpolation(0, x, xi, polyorder);
          },
          "val"_a,
          "vec"_a,
          "polyorder"_a = Index{1},
          "Construct a Lagrange interpolation object")
      .def_prop_ro(
          "order",
          [](const LagrangeCyclic0to360Interpolation& l) {
            return l.size() - 1;
          },
          "The order of the polynomial interpolation")
      .def_rw("pos",
              &LagrangeCyclic0to360Interpolation::pos,
              "The interpolation positions")
      .def_rw("weights",
              &LagrangeCyclic0to360Interpolation::lx,
              "The interpolation weights")
      .def("__getstate__",
           [](const LagrangeCyclic0to360Interpolation& self) {
             return std::make_tuple(self.pos, self.lx, self.dlx);
           })
      .def("__setstate__",
           [](LagrangeCyclic0to360Interpolation* l,
              const std::tuple<Index, Array<Numeric>, Array<Numeric>>& state) {
             new (l) LagrangeCyclic0to360Interpolation{};

             l->pos = std::get<0>(state);
             l->lx  = std::get<1>(state);
             l->dlx = std::get<2>(state);
           })
      .doc() = "Cyclic [0, 360] polynomial interpolation object";

  py::class_<LagrangeCyclicPM180Interpolation>(interp, "CyclicGridPM180")
      .def(
          "__init__",
          [](LagrangeCyclicPM180Interpolation* l,
             const Numeric x,
             const Vector& xi,
             const Index polyorder) {
            new (l) LagrangeCyclicPM180Interpolation(0, x, xi, polyorder);
          },
          "val"_a,
          "vec"_a,
          "polyorder"_a = Index{1},
          "Construct a Lagrange interpolation object")
      .def_prop_ro(
          "order",
          [](const LagrangeCyclicPM180Interpolation& l) {
            return l.size() - 1;
          },
          "The order of the polynomial interpolation")
      .def_rw("pos",
              &LagrangeCyclicPM180Interpolation::pos,
              "The interpolation positions")
      .def_rw("weights",
              &LagrangeCyclicPM180Interpolation::lx,
              "The interpolation weights")
      .def("__getstate__",
           [](const LagrangeCyclicPM180Interpolation& self) {
             return std::make_tuple(self.pos, self.lx, self.dlx);
           })
      .def("__setstate__",
           [](LagrangeCyclicPM180Interpolation* l,
              const std::tuple<Index, Array<Numeric>, Array<Numeric>>& state) {
             new (l) LagrangeCyclicPM180Interpolation{};

             l->pos = std::get<0>(state);
             l->lx  = std::get<1>(state);
             l->dlx = std::get<2>(state);
           })
      .doc() = "Cyclic [-180, 180] polynomial interpolation object";

  interp.def(
      "interp",
      [](const Vector& yi, const PythonLags& lag1) {
        return std::visit([&yi](auto& l1) { return my_interp::interp(yi, l1); },
                          lag1);
      },
      "yi"_a,
      "lag1"_a,
      "Interpolate a 1D vector");

  interp.def(
      "interp",
      [](const Matrix& yi, const PythonLags& lag1, const PythonLags& lag2) {
        return std::visit(
            [&yi](auto& l1, auto& l2) { return my_interp::interp(yi, l1, l2); },
            lag1,
            lag2);
      },
      "yi"_a,
      "lag1"_a,
      "lag2"_a,
      "Interpolate a 2D matrix");

  interp.def(
      "interp",
      [](const Tensor3& yi,
         const PythonLags& lag1,
         const PythonLags& lag2,
         const PythonLags& lag3) {
        return std::visit(
            [&yi](auto& l1, auto& l2, auto& l3) {
              return my_interp::interp(yi, l1, l2, l3);
            },
            lag1,
            lag2,
            lag3);
      },
      "yi"_a,
      "lag1"_a,
      "lag2"_a,
      "lag3"_a,
      "Interpolate a 3D tensor");

  interp.def(
      "interp1d",
      [](const Vector& yi,
         const Vector& xi,
         const Vector& x,
         const Index polyorder,
         bool cyclic) {
        if (cyclic)
          return reinterp(yi,
                          lagrange_interp::make_lags<lagrange_interp::loncross>(
                              xi, x, polyorder, 0.0));
        return reinterp(yi,
                        lagrange_interp::make_lags<>(xi, x, polyorder, 0.0));
      },
      "yi"_a,
      "xi"_a,
      "x"_a,
      "polyorder"_a = 1,
      "cyclic"_a    = false,
      "Interpolate a 1D vector");
}
}  // namespace Python
