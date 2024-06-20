#include <python_interface.h>

#include <type_traits>
#include <variant>

#include "debug.h"
#include "py_macros.h"

namespace Python {
void py_rte(py::module_& m) try {
  py::class_<GasAbsLookup>(m, "GasAbsLookup")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup,
          species,
          Species,
          Species,
          ":class:`~pyarts.arts.ArrayOfArrayOfSpeciesTag` Active species")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup,
          non_linear_species,
          NonLinearSpecies,
          NonLinearSpecies,
          ":class:`~pyarts.arts.ArrayOfIndex` Non-linear species")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup,
          f_grid,
          Fgrid,
          Fgrid,
          ":class:`~pyarts.arts.ArrayOfIndex` Trained frequencies")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup,
          flag_default,
          FLAGDefault,
          FLAGDefault,
          ":class:`~pyarts.arts.ArrayOfLagrangeInterpolation` Active interpolation")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup,
          p_grid,
          Pgrid,
          Pgrid,
          ":class:`~pyarts.arts.Vector` Trained pressures")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup,
          log_p_grid,
          LogPgrid,
          LogPgrid,
          ":class:`~pyarts.arts.Vector` Trained pressure in log")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup,
          vmrs,
          VMRs,
          VMRs,
          ":class:`~pyarts.arts.Matrix` Training VMRs")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup,
          t_ref,
          Tref,
          Tref,
          ":class:`~pyarts.arts.Vector` Trained temperatures")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup,
          t_pert,
          Tpert,
          Tpert,
          ":class:`~pyarts.arts.Vector` Temperature perturbations")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup,
          nls_pert,
          NLSPert,
          NLSPert,
          ":class:`~pyarts.arts.Vector` Non-linear perturbations")
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup,
          xsec,
          Xsec,
          Xsec,
          ":class:`~pyarts.arts.Tensor4` Cross-section data")
      .def("__getstate__",
           [](GasAbsLookup& self) {
             return std::tuple<ArrayOfArrayOfSpeciesTag,
                               ArrayOfIndex,
                               Vector,
                               ArrayOfLagrangeInterpolation,
                               Vector,
                               Matrix,
                               Vector,
                               Vector,
                               Vector,
                               Tensor4>{self.Species(),
                                        self.NonLinearSpecies(),
                                        self.Fgrid(),
                                        self.FLAGDefault(),
                                        self.Pgrid(),
                                        self.LogPgrid(),
                                        self.VMRs(),
                                        self.Tref(),
                                        self.Tpert(),
                                        self.NLSPert(),
                                        self.Xsec()};
           })
      .def("__setstate__",
           [](GasAbsLookup* self,
              const std::tuple<ArrayOfArrayOfSpeciesTag,
                               ArrayOfIndex,
                               Vector,
                               ArrayOfLagrangeInterpolation,
                               Vector,
                               Matrix,
                               Vector,
                               Vector,
                               Vector,
                               Tensor4>& state) {
             new (self) GasAbsLookup();

             self->Species()          = std::get<0>(state);
             self->NonLinearSpecies() = std::get<1>(state);
             self->Fgrid()            = std::get<2>(state);
             self->FLAGDefault()      = std::get<3>(state);
             self->Pgrid()            = std::get<4>(state);
             self->LogPgrid()         = std::get<5>(state);
             self->VMRs()             = std::get<6>(state);
             self->Tref()             = std::get<7>(state);
             self->Tpert()            = std::get<8>(state);
             self->NLSPert()          = std::get<9>(state);
             self->Xsec()             = std::get<10>(state);
           });
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize rte\n", e.what()));
}
}  // namespace Python
