#include <debug.h>
#include <python_interface.h>

#include "hpy_arts.h"

namespace Python {
void py_rte(py::module_& m) try {
  py::class_<GasAbsLookup> galu(m, "GasAbsLookup");
  workspace_group_interface(galu);
  galu.def_rw("species",
              &GasAbsLookup::species,
              ":class:`~pyarts.arts.ArrayOfArrayOfSpeciesTag` Active species")
      .def_rw("non_linear_species",
              &GasAbsLookup::nonlinear_species,
              ":class:`~pyarts.arts.ArrayOfIndex` Non-linear species")
      .def_rw("f_grid",
              &GasAbsLookup::f_grid,
              ":class:`~pyarts.arts.ArrayOfIndex` Trained frequencies")
      .def_rw(
          "flag_default",
          &GasAbsLookup::flag_default,
          ":class:`~pyarts.arts.ArrayOfLagrangeInterpolation` Active interpolation")
      .def_rw("p_grid",
              &GasAbsLookup::p_grid,
              ":class:`~pyarts.arts.Vector` Trained pressures")
      .def_rw("log_p_grid",
              &GasAbsLookup::log_p_grid,
              ":class:`~pyarts.arts.Vector` Trained pressure in log")
      .def_rw("vmrs",
              &GasAbsLookup::vmrs_ref,
              ":class:`~pyarts.arts.Matrix` Training VMRs")
      .def_rw("t_ref",
              &GasAbsLookup::t_ref,
              ":class:`~pyarts.arts.Vector` Trained temperatures")
      .def_rw("t_pert",
              &GasAbsLookup::t_pert,
              ":class:`~pyarts.arts.Vector` Temperature perturbations")
      .def_rw("nls_pert",
              &GasAbsLookup::nls_pert,
              ":class:`~pyarts.arts.Vector` Non-linear perturbations")
      .def_rw("xsec",
              &GasAbsLookup::xsec,
              ":class:`~pyarts.arts.Tensor4` Cross-section data")
      .def("__getstate__",
           [](GasAbsLookup& self) {
             return std::tuple<ArrayOfArrayOfSpeciesTag,
                               ArrayOfIndex,
                               Vector,
                               ArrayOfLagrangeInterpolation,
                               Vector,
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
