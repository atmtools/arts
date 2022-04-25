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
  py::class_<TransmissionMatrix>(m, "TransmissionMatrix", py::buffer_protocol())
      .def(py::init([](Index nf, Index ns) {
             ARTS_USER_ERROR_IF(nf < 0, "Bad frequency size")
             ARTS_USER_ERROR_IF(ns < 1 or ns > 4, "Bad stokes_dim")
             return new TransmissionMatrix(nf, ns);
           }),
           py::arg("nf") = 0,
           py::arg("stokes_dim") = 1)
      .PythonInterfaceCopyValue(TransmissionMatrix)
      .PythonInterfaceWorkspaceVariableConversion(TransmissionMatrix)
      .PythonInterfaceFileIO(TransmissionMatrix)
      .PythonInterfaceBasicRepresentation(TransmissionMatrix)
      .def_buffer([](TransmissionMatrix& t) -> py::buffer_info {
        Numeric* ptr = t.stokes_dim == 1
                           ? t.T1.data()->data()
                           : (t.stokes_dim == 2
                                  ? t.T2.data()->data()
                                  : (t.stokes_dim == 3 ? t.T3.data()->data()
                                                       : t.T4.data()->data()));
        return py::buffer_info(ptr,
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               3,
                               {t.Frequencies(), t.stokes_dim, t.stokes_dim},
                               {sizeof(Numeric) * t.stokes_dim * t.stokes_dim,
                                sizeof(Numeric) * t.stokes_dim,
                                sizeof(Numeric)});
      })
      .def_property("value",
                    py::cpp_function(
                        [](TransmissionMatrix& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](TransmissionMatrix& x, TransmissionMatrix& y) { x = y; })
      .PythonInterfaceValueOperators.def(py::pickle(
          [](const TransmissionMatrix& self) {
            return py::make_tuple(
                self.stokes_dim, self.T1, self.T2, self.T3, self.T4);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 5, "Invalid state!")

            auto* out = new TransmissionMatrix{};
            out->stokes_dim = t[0].cast<Index>();
            out->T1 = t[1].cast<decltype(out->T1)>();
            out->T2 = t[2].cast<decltype(out->T2)>();
            out->T3 = t[3].cast<decltype(out->T3)>();
            out->T4 = t[4].cast<decltype(out->T4)>();
            return out;
          }));

  py::class_<RadiationVector>(m, "RadiationVector", py::buffer_protocol())
      .def(py::init([](Index nf, Index ns) {
             ARTS_USER_ERROR_IF(nf < 0, "Bad requency size, valid: [0, ...)")
             ARTS_USER_ERROR_IF(ns < 1 or ns > 4,
                                "Bad stokes dimensionality, valid: [1, 4]")
             return new RadiationVector(nf, ns);
           }),
           py::arg("nf") = 0,
           py::arg("ns") = 1,
           py::doc("Init by frequency size and stokes dimensionality"))
      .PythonInterfaceCopyValue(RadiationVector)
      .PythonInterfaceWorkspaceVariableConversion(RadiationVector)
      .PythonInterfaceFileIO(RadiationVector)
      .PythonInterfaceBasicRepresentation(RadiationVector)
      .def_buffer([](RadiationVector& t) -> py::buffer_info {
        Numeric* ptr = t.stokes_dim == 1
                           ? t.R1.data()->data()
                           : (t.stokes_dim == 2
                                  ? t.R2.data()->data()
                                  : (t.stokes_dim == 3 ? t.R3.data()->data()
                                                       : t.R4.data()->data()));
        return py::buffer_info(
            ptr,
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            2,
            {t.Frequencies(), t.stokes_dim},
            {sizeof(Numeric) * t.stokes_dim, sizeof(Numeric)});
      })
      .def_property("value",
                    py::cpp_function(
                        [](RadiationVector& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](RadiationVector& x, RadiationVector& y) { x = y; })
      .PythonInterfaceValueOperators.def(py::pickle(
          [](const RadiationVector& self) {
            return py::make_tuple(
                self.stokes_dim, self.R1, self.R2, self.R3, self.R4);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 5, "Invalid state!")

            auto* out = new RadiationVector{};
            out->stokes_dim = t[0].cast<Index>();
            out->R1 = t[1].cast<decltype(out->R1)>();
            out->R2 = t[2].cast<decltype(out->R2)>();
            out->R3 = t[3].cast<decltype(out->R3)>();
            out->R4 = t[4].cast<decltype(out->R4)>();
            return out;
          }));

  py::class_<PropagationMatrix>(m, "PropagationMatrix")
      .def(py::init([](Index nf, Index ns, Index nza, Index naa, Numeric v) {
             ARTS_USER_ERROR_IF(ns < 1 or ns > 4,
                                "Bad Stokes Dimension ",
                                ns,
                                " should be [1, 4]")
             ARTS_USER_ERROR_IF(nf < 0 or nza < 0 or naa < 0,
                                "Negative size index")
             return new PropagationMatrix(nf, ns, nza, naa, v);
           }),
           py::arg("nf") = 0,
           py::arg("ns") = 1,
           py::arg("nza") = 1,
           py::arg("naa") = 1,
           py::arg("v") = 0)
      .PythonInterfaceCopyValue(PropagationMatrix)
      .PythonInterfaceWorkspaceVariableConversion(PropagationMatrix)
      .PythonInterfaceFileIO(PropagationMatrix)
      .PythonInterfaceBasicRepresentation(PropagationMatrix)
      .def_property(
          "data",
          py::cpp_function(
              [](PropagationMatrix& s) -> Tensor4& { return s.Data(); },
              py::return_value_policy::reference_internal),
          [](PropagationMatrix& s, const Tensor4& d) {
            ARTS_USER_ERROR_IF(s.Data().nbooks() not_eq d.nbooks() or
                                   s.Data().npages() not_eq d.npages() or
                                   s.Data().nrows() not_eq d.nrows() or
                                   s.Data().ncols() not_eq d.ncols(),
                               "Expect shape: (",
                               s.Data().shape(),
                               ')',
                               "\nGot shape: (",
                               d.shape(),
                               ')')
            s.Data() = d;
          })
      .def(py::pickle(
          [](const PropagationMatrix& self) {
            return py::make_tuple(self.NumberOfFrequencies(),
                                  self.StokesDimensions(),
                                  self.NumberOfZenithAngles(),
                                  self.NumberOfAzimuthAngles(),
                                  self.Data());
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 5, "Invalid state!")

            auto* out = new PropagationMatrix{t[0].cast<Index>(),
                                              t[1].cast<Index>(),
                                              t[2].cast<Index>(),
                                              t[3].cast<Index>()};
            out->Data() = t[4].cast<Tensor4>();
            return out;
          }));

  py::class_<StokesVector>(m, "StokesVector")
      .def(py::init([](Index nf, Index ns, Index nza, Index naa, Numeric v) {
             ARTS_USER_ERROR_IF(ns < 1 or ns > 4,
                                "Bad Stokes Dimension ",
                                ns,
                                " should be [1, 4]")
             ARTS_USER_ERROR_IF(nf < 0 or nza < 0 or naa < 0,
                                "Negative size index")
             return new StokesVector(nf, ns, nza, naa, v);
           }),
           py::arg_v("nf", Index(0), "Index(0)"),
           py::arg_v("ns", Index(1), "Index(1)"),
           py::arg_v("nza", Index(1), "Index(1)"),
           py::arg_v("naa", Index(1), "Index(1)"),
           py::arg_v("v", Numeric(0), "Numeric(0)"))
      .PythonInterfaceCopyValue(StokesVector)
      .PythonInterfaceWorkspaceVariableConversion(StokesVector)
      .PythonInterfaceFileIO(StokesVector)
      .PythonInterfaceBasicRepresentation(StokesVector)
      .def_property(
          "data",
          py::cpp_function([](StokesVector& s) -> Tensor4& { return s.Data(); },
                           py::return_value_policy::reference_internal),
          [](StokesVector& s, const Tensor4& d) {
            ARTS_USER_ERROR_IF(s.Data().nbooks() not_eq d.nbooks() or
                                   s.Data().npages() not_eq d.npages() or
                                   s.Data().nrows() not_eq d.nrows() or
                                   s.Data().ncols() not_eq d.ncols(),
                               "\nExpect shape: (",
                               s.Data().shape(),
                               ')',
                               "\nGot shape: (",
                               d.shape(),
                               ')')
            s.Data() = d;
          })
      .def(py::pickle(
          [](const StokesVector& self) {
            return py::make_tuple(self.NumberOfFrequencies(),
                                  self.StokesDimensions(),
                                  self.NumberOfZenithAngles(),
                                  self.NumberOfAzimuthAngles(),
                                  self.Data());
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 5, "Invalid state!")

            auto* out = new StokesVector{t[0].cast<Index>(),
                                         t[1].cast<Index>(),
                                         t[2].cast<Index>(),
                                         t[3].cast<Index>()};
            out->Data() = t[4].cast<Tensor4>();
            return out;
          }));

  PythonInterfaceWorkspaceArray(TransmissionMatrix);
  PythonInterfaceWorkspaceArray(PropagationMatrix);
  PythonInterfaceWorkspaceArray(RadiationVector);
  PythonInterfaceWorkspaceArray(StokesVector);

  PythonInterfaceWorkspaceArray(ArrayOfTransmissionMatrix);
  PythonInterfaceWorkspaceArray(ArrayOfPropagationMatrix);
  PythonInterfaceWorkspaceArray(ArrayOfRadiationVector);
  PythonInterfaceWorkspaceArray(ArrayOfStokesVector);

  py::class_<LagrangeInterpolation>(m, "LagrangeInterpolation")
      .def(py::init([]() { return new LagrangeInterpolation{}; }))
      .PythonInterfaceCopyValue(LagrangeInterpolation)
      .def(py::pickle(
          [](const LagrangeInterpolation& self) {
            return py::make_tuple(self.pos, self.lx, self.dlx);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 3, "Invalid state!")

            auto* out = new LagrangeInterpolation{};

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
      .def(py::init([]() { return new GasAbsLookup{}; }))
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

            auto* out = new GasAbsLookup{};

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
}
}  // namespace Python
