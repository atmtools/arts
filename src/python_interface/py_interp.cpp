#include <interp.h>
#include <python_interface.h>

#include <variant>

#include "py_macros.h"

namespace Python {
using PythonLags =
    std::variant<LagrangeInterpolation, LagrangeCyclic0to360Interpolation, LagrangeCyclicPM180Interpolation>;

void py_interp(py::module_& m) {
  auto interp = m.def_submodule("interp", "Interpolation methods");

  artsclass<LagrangeInterpolation>(interp, "PolyGrid")
      .def(py::init(
               [](const Numeric x, const Vector& xi, const Index polyorder) {
                 return LagrangeInterpolation(0, x, xi, polyorder);
               }),
           py::arg("val"),
           py::arg("vec"),
           py::arg("polyorder") = Index{1},
           py::doc("Construct a Lagrange interpolation object"))
      .def_property_readonly(
          "order", [](const LagrangeInterpolation& l) { return l.size() - 1; })
      .def_readwrite("pos", &LagrangeInterpolation::pos)
      .def_readwrite("weights", &LagrangeInterpolation::lx)
      .PythonInterfaceBasicRepresentation(LagrangeInterpolation)
      .def(py::pickle(
          [](const LagrangeInterpolation& self) {
            return py::make_tuple(self.pos, self.lx, self.dlx);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 3, "Invalid state!")

            auto out = std::make_shared<LagrangeInterpolation>();

            out->pos = t[0].cast<Index>();
            out->lx = t[1].cast<Array<Numeric>>();
            out->dlx = t[2].cast<Array<Numeric>>();
            return out;
          })).doc() = "Polynomial interpolation object";

  artsarray<ArrayOfLagrangeInterpolation>(interp, "ArrayOfLagrangeInterpolation")
      .doc() = "List of :class:`~pyarts.arts.LagrangeInterpolation`";

  artsclass<LagrangeCyclic0to360Interpolation>(interp, "CyclicGrid0to360")
      .def(py::init(
               [](const Numeric x, const Vector& xi, const Index polyorder) {
                 return LagrangeCyclic0to360Interpolation(0, x, xi, polyorder);
               }),
           py::arg("val"),
           py::arg("vec"),
           py::arg("polyorder") = Index{1},
           py::doc("Construct a Lagrange interpolation object"))
      .def_property_readonly("order",
                             [](const LagrangeCyclic0to360Interpolation& l) {
                               return l.size() - 1;
                             })
      .def_readwrite("pos", &LagrangeCyclic0to360Interpolation::pos)
      .def_readwrite("weights", &LagrangeCyclic0to360Interpolation::lx)
      .PythonInterfaceBasicRepresentation(LagrangeCyclic0to360Interpolation)
      .PythonInterfaceBasicRepresentation(LagrangeCyclic0to360Interpolation)
      .def(py::pickle(
          [](const LagrangeCyclic0to360Interpolation& self) {
            return py::make_tuple(self.pos, self.lx, self.dlx);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 3, "Invalid state!")

            auto out = std::make_shared<LagrangeCyclic0to360Interpolation>();

            out->pos = t[0].cast<Index>();
            out->lx = t[1].cast<Array<Numeric>>();
            out->dlx = t[2].cast<Array<Numeric>>();
            return out;
          })).doc() = "Cyclic [0, 360] polynomial interpolation object";

  artsclass<LagrangeCyclicPM180Interpolation>(interp, "CyclicGridPM180")
      .def(py::init(
               [](const Numeric x, const Vector& xi, const Index polyorder) {
                 return LagrangeCyclicPM180Interpolation(0, x, xi, polyorder);
               }),
           py::arg("val"),
           py::arg("vec"),
           py::arg("polyorder") = Index{1},
           py::doc("Construct a Lagrange interpolation object"))
      .def_property_readonly("order",
                             [](const LagrangeCyclicPM180Interpolation& l) {
                               return l.size() - 1;
                             })
      .def_readwrite("pos", &LagrangeCyclicPM180Interpolation::pos)
      .def_readwrite("weights", &LagrangeCyclicPM180Interpolation::lx)
      .PythonInterfaceBasicRepresentation(LagrangeCyclicPM180Interpolation)
      .PythonInterfaceBasicRepresentation(LagrangeCyclicPM180Interpolation)
      .def(py::pickle(
          [](const LagrangeCyclicPM180Interpolation& self) {
            return py::make_tuple(self.pos, self.lx, self.dlx);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 3, "Invalid state!")

            auto out = std::make_shared<LagrangeCyclicPM180Interpolation>();

            out->pos = t[0].cast<Index>();
            out->lx = t[1].cast<Array<Numeric>>();
            out->dlx = t[2].cast<Array<Numeric>>();
            return out;
          })).doc() = "Cyclic [-180, 180] polynomial interpolation object";

  interp.def("interp", [](const Vector& yi, const PythonLags& lag1) {
    return std::visit([&yi](auto& l1) { return my_interp::interp(yi, l1); },
                      lag1);
  });

  interp.def(
      "interp",
      [](const Matrix& yi, const PythonLags& lag1, const PythonLags& lag2) {
        return std::visit(
            [&yi](auto& l1, auto& l2) { return my_interp::interp(yi, l1, l2); },
            lag1,
            lag2);
      });

  interp.def("interp",
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
             });
}
}  // namespace Python
