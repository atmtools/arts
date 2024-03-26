#include <python_interface.h>

#include <type_traits>
#include <variant>

#include "debug.h"
#include "py_macros.h"

namespace Python {
void py_rte(py::module_& m) try {
  py_staticGasAbsLookup(m)
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, species, Species, Species, ":class:`~pyarts.arts.ArrayOfArrayOfSpeciesTag` Active species")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, non_linear_species, NonLinearSpecies, NonLinearSpecies, ":class:`~pyarts.arts.ArrayOfIndex` Non-linear species")
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, f_grid, Fgrid, Fgrid, ":class:`~pyarts.arts.ArrayOfIndex` Trained frequencies")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, flag_default, FLAGDefault, FLAGDefault, ":class:`~pyarts.arts.ArrayOfLagrangeInterpolation` Active interpolation")
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, p_grid, Pgrid, Pgrid, ":class:`~pyarts.arts.Vector` Trained pressures")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, log_p_grid, LogPgrid, LogPgrid, ":class:`~pyarts.arts.Vector` Trained pressure in log")
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, vmrs, VMRs, VMRs, ":class:`~pyarts.arts.Matrix` Training VMRs")
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, t_ref, Tref, Tref, ":class:`~pyarts.arts.Vector` Trained temperatures")
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, t_pert, Tpert, Tpert, ":class:`~pyarts.arts.Vector` Temperature perturbations")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, nls_pert, NLSPert, NLSPert, ":class:`~pyarts.arts.Vector` Non-linear perturbations")
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, xsec, Xsec, Xsec, ":class:`~pyarts.arts.Tensor4` Cross-section data")
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

            auto out = std::make_shared<GasAbsLookup>();

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
          }));
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize rte\n", e.what()));
}
}  // namespace Python
