#include <py_auto_interface.h>
#include <pybind11/detail/common.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#include <type_traits>
#include <variant>

#include "debug.h"
#include "py_macros.h"

namespace Python {
void py_rte(py::module_& m) {
  py::class_<LagrangeInterpolation>(m, "LagrangeInterpolation")
      .def(py::init([]() { return std::make_unique<LagrangeInterpolation>(); }))
      .PythonInterfaceCopyValue(LagrangeInterpolation)
      .def(py::pickle(
          [](const LagrangeInterpolation& self) {
            return py::make_tuple(self.pos, self.lx, self.dlx);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 3, "Invalid state!")

            auto out = std::make_unique<LagrangeInterpolation>();

            out->pos = t[0].cast<Index>();
            out->lx = t[1].cast<Array<Numeric>>();
            out->dlx = t[2].cast<Array<Numeric>>();
            return out;
          }));

  py::class_<ArrayOfLagrangeInterpolation>(m, "ArrayOfLagrangeInterpolation")
      .PythonInterfaceArrayDefault(LagrangeInterpolation)
      .PythonInterfaceBasicRepresentation(ArrayOfLagrangeInterpolation);
  py::implicitly_convertible<std::vector<LagrangeInterpolation>,
                             ArrayOfLagrangeInterpolation>();

  py::class_<GasAbsLookup>(m, "GasAbsLookup")
      .def(py::init([]() { return std::make_unique<GasAbsLookup>(); }))
      .PythonInterfaceCopyValue(GasAbsLookup)
      .PythonInterfaceWorkspaceVariableConversion(GasAbsLookup)
      .PythonInterfaceFileIO(GasAbsLookup)
      .PythonInterfaceBasicRepresentation(GasAbsLookup)
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, species, Species, Species)
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, non_linear_species, NonLinearSpecies, NonLinearSpecies)
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, f_grid, Fgrid, Fgrid)
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, flag_default, FLAGDefault, FLAGDefault)
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, p_grid, Pgrid, Pgrid)
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, log_p_grid, LogPgrid, LogPgrid)
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, vmrs, VMRs, VMRs)
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, t_ref, Tref, Tref)
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, t_pert, Tpert, Tpert)
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, nls_pert, NLSPert, NLSPert)
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, xsec, Xsec, Xsec)
      .def(py::pickle(
          [](GasAbsLookup& self) {
            return py::make_tuple(self.Species(),
                                  self.NonLinearSpecies(),
                                  self.Fgrid(),
                                  self.FLAGDefault(),
                                  self.Pgrid(),
                                  self.LogPgrid(),
                                  self.VMRs(),
                                  self.Tref(),
                                  self.Tpert(),
                                  self.NLSPert(),
                                  self.Xsec());
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 11, "Invalid state!")

            auto out = std::make_unique<GasAbsLookup>();

            out->Species() = t[0].cast<ArrayOfArrayOfSpeciesTag>();
            out->NonLinearSpecies() = t[1].cast<ArrayOfIndex>();
            out->Fgrid() = t[2].cast<Vector>();
            out->FLAGDefault() = t[3].cast<ArrayOfLagrangeInterpolation>();
            out->Pgrid() = t[4].cast<Vector>();
            out->LogPgrid() = t[5].cast<Vector>();
            out->VMRs() = t[6].cast<Matrix>();
            out->Tref() = t[7].cast<Vector>();
            out->Tpert() = t[8].cast<Vector>();
            out->NLSPert() = t[9].cast<Vector>();
            out->Xsec() = t[10].cast<Tensor4>();

            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(GasAbsLookup);
}
}  // namespace Python
