#include <py_auto_interface.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "debug.h"
#include "py_macros.h"

namespace Python {
void py_rte(py::module_& m) {
  py::class_<TransmissionMatrix>(m, "TransmissionMatrix")
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
      .def(
          "__getitem__",
          [](TransmissionMatrix& t,
             Index i) -> std::variant<Eigen::Ref<Eigen::Matrix4d>,
                                      Eigen::Ref<Eigen::Matrix3d>,
                                      Eigen::Ref<Eigen::Matrix2d>,
                                      Eigen::Ref<Eigen::Matrix<double, 1, 1>>> {
            if (i < 0 or i > t.Frequencies())
              throw std::out_of_range(
                  var_string(i, " in range [0", t.Frequencies(), ')'));
            switch (t.StokesDim()) {
              case 1:
                return Eigen::Ref<Eigen::Matrix<double, 1, 1>>(t.Mat1(i));
              case 2:
                return Eigen::Ref<Eigen::Matrix2d>(t.Mat2(i));
              case 3:
                return Eigen::Ref<Eigen::Matrix3d>(t.Mat3(i));
              case 4:
                return Eigen::Ref<Eigen::Matrix4d>(t.Mat4(i));
            }
            throw std::out_of_range("bad stokes dim");
          })
      .def("__setitem__",
           [](TransmissionMatrix& t, Index i, const py::array_t<Numeric>& arr) {
             if (i < 0 or i > t.Frequencies())
               throw std::out_of_range(
                   var_string(i, " in range [0", t.Frequencies(), ')'));
             ARTS_USER_ERROR_IF(arr.request().ndim not_eq 2, "Bad size array")

             auto arr_val = arr.unchecked<2>();
             ARTS_USER_ERROR_IF(arr_val.shape(0) not_eq arr_val.shape(1) or
                                    arr_val.shape(0) not_eq t.StokesDim(),
                                "Bad input size!")
             const Index n = t.StokesDim();

             switch (n) {
               case 1:
                 for (Index r = 0; r < n; r++)
                   for (Index c = 0; c < n; c++)
                     t.Mat1(i)(r, c) = arr_val(r, c);
                 break;
               case 2:
                 for (Index r = 0; r < n; r++)
                   for (Index c = 0; c < n; c++)
                     t.Mat2(i)(r, c) = arr_val(r, c);
                 break;
               case 3:
                 for (Index r = 0; r < n; r++)
                   for (Index c = 0; c < n; c++)
                     t.Mat3(i)(r, c) = arr_val(r, c);
                 break;
               case 4:
                 for (Index r = 0; r < n; r++)
                   for (Index c = 0; c < n; c++)
                     t.Mat4(i)(r, c) = arr_val(r, c);
                 break;
             }
           });

  py::class_<RadiationVector>(m, "RadiationVector")
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
      .def(
          "__getitem__",
          [](RadiationVector& t,
             Index i) -> std::variant<Eigen::Ref<Eigen::Vector4d>,
                                      Eigen::Ref<Eigen::Vector3d>,
                                      Eigen::Ref<Eigen::Vector2d>,
                                      Eigen::Ref<Eigen::Matrix<double, 1, 1>>> {
            if (i < 0 or i > t.Frequencies())
              throw std::out_of_range(
                  var_string(i, " in range [0", t.Frequencies(), ')'));
            switch (t.StokesDim()) {
              case 1:
                return Eigen::Ref<Eigen::Matrix<double, 1, 1>>(t.Vec1(i));
              case 2:
                return Eigen::Ref<Eigen::Vector2d>(t.Vec2(i));
              case 3:
                return Eigen::Ref<Eigen::Vector3d>(t.Vec3(i));
              case 4:
                return Eigen::Ref<Eigen::Vector4d>(t.Vec4(i));
            }
            throw std::out_of_range("bad stokes dim");
          })
      .def(
          "__setitem__",
          [](RadiationVector& t, Index i, const py::array_t<Numeric>& arr) {
            if (i < 0 or i > t.Frequencies())
              throw std::out_of_range(
                  var_string(i, " in range [0", t.Frequencies(), ')'));
            ARTS_USER_ERROR_IF(arr.request().ndim not_eq 1, "Bad size array")

            auto arr_val = arr.unchecked<1>();
            ARTS_USER_ERROR_IF(arr_val.shape(0) not_eq t.StokesDim(),
                               "Bad input size!")
            const Index n = t.StokesDim();

            switch (n) {
              case 1:
                for (Index r = 0; r < n; r++) t.Vec1(i)[r] = *arr_val.data(r);
                break;
              case 2:
                for (Index r = 0; r < n; r++) t.Vec2(i)[r] = *arr_val.data(r);
                break;
              case 3:
                for (Index r = 0; r < n; r++) t.Vec3(i)[r] = *arr_val.data(r);
                break;
              case 4:
                for (Index r = 0; r < n; r++) t.Vec4(i)[r] = *arr_val.data(r);
                break;
            }
          },
          py::arg("i"),
          py::arg("vec"),
          py::doc("Sets a single vector"));

  py::class_<StokesVector>(m, "StokesVector")
      .def(py::init([](Index nf, Index ns, Index nza, Index naa, Numeric v) {
             ARTS_USER_ERROR_IF(ns < 1 or ns > 4,
                                "Bad Stokes Dimension ",
                                ns,
                                " should be [1, 4]")
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
                               s.Data().nbooks(),
                               ", ",
                               s.Data().npages(),
                               ", ",
                               s.Data().nrows(),
                               ", ",
                               s.Data().ncols(),
                               ')',
                               "\nGot shape: (",
                               d.nbooks(),
                               ", ",
                               d.npages(),
                               ", ",
                               d.nrows(),
                               ", ",
                               d.ncols(),
                               ')')
            s.Data() = d;
          });

  py::class_<PropagationMatrix>(m, "PropagationMatrix")
      .def(py::init([](Index nf, Index ns, Index nza, Index naa, Numeric v) {
             ARTS_USER_ERROR_IF(ns < 1 or ns > 4,
                                "Bad Stokes Dimension ",
                                ns,
                                " should be [1, 4]")
             ARTS_USER_ERROR_IF(nf < 0 or ns < 0 or nza < 0 or naa < 0,
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
                               s.Data().nbooks(),
                               ", ",
                               s.Data().npages(),
                               ", ",
                               s.Data().nrows(),
                               ", ",
                               s.Data().ncols(),
                               ')')
            s.Data() = d;
          });

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
      .PythonInterfaceCopyValue(LagrangeInterpolation);

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
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, xsec, Xsec, Xsec);
}
}  // namespace Python