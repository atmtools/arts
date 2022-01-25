#include <auto_md.h>
#include <xml_io.h>

#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_rte(py::module_& m) {
  py::class_<TransmissionMatrix>(m, "TransmissionMatrix")
      .def(py::init<>())
      .def(py::init<Index>())
      .def(py::init<Index, Index>())
      .PythonInterfaceFileIO(TransmissionMatrix)
      .PythonInterfaceBasicRepresentation(TransmissionMatrix)
      .def("__getitem__",
           [](const TransmissionMatrix& t, Index i) {
             if (i < 0 or i > t.Frequencies())
               throw std::out_of_range(var_string("Out of bounds"));
             return t.Mat(i);
           })
      .def("__setitem__",
           [](TransmissionMatrix& t, Index i, const py::array_t<Numeric>& arr) {
             if (i < 0 or i > t.Frequencies())
               throw std::out_of_range(var_string("Out of bounds"));
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
      .def(py::init<>())
      .def(py::init<Index>())
      .def(py::init<Index, Index>())
      .PythonInterfaceFileIO(RadiationVector)
      .PythonInterfaceBasicRepresentation(RadiationVector)
      .def("__getitem__",
           [](const RadiationVector& t, Index i) {
             ARTS_USER_ERROR_IF(i < 0 or i > t.Frequencies(), "Out of bounds")
             return t.Vec(i);
           })
      .def("__setitem__",
           [](RadiationVector& t, Index i, const py::array_t<Numeric>& arr) {
             ARTS_USER_ERROR_IF(i < 0 or i > t.Frequencies(), "Out of bounds")
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
           });

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
      .PythonInterfaceFileIO(StokesVector)
      .PythonInterfaceBasicRepresentation(StokesVector)
      .def_property(
          "data",
          [](const StokesVector& s) { return s.Data(); },
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
             return new PropagationMatrix(nf, ns, nza, naa, v);
           }),
           py::arg_v("nf", Index(0), "Index(0)"),
           py::arg_v("ns", Index(1), "Index(1)"),
           py::arg_v("nza", Index(1), "Index(1)"),
           py::arg_v("naa", Index(1), "Index(1)"),
           py::arg_v("v", Numeric(0), "Numeric(0)"))
      .PythonInterfaceFileIO(PropagationMatrix)
      .PythonInterfaceBasicRepresentation(PropagationMatrix)
      .def_property(
          "data",
          [](const PropagationMatrix& s) { return s.Data(); },
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

#define PythonInterfaceBasicReferenceProperty(       \
    Type, PropertyName, ReadFunction, WriteFunction) \
  def_property(                                      \
      #PropertyName,                                 \
      [](Type& x) { return x.ReadFunction(); },      \
      [](Type& x, decltype(x.ReadFunction()) y) {    \
        return x.WriteFunction() = std::move(y);     \
      })

  py::class_<LagrangeInterpolation>(m, "LagrangeInterpolation")
      .def(py::init<>());

  py::class_<ArrayOfLagrangeInterpolation>(m, "ArrayOfLagrangeInterpolation")
      .PythonInterfaceArrayDefault(LagrangeInterpolation)
      .PythonInterfaceBasicRepresentation(ArrayOfLagrangeInterpolation);
  py::implicitly_convertible<std::vector<LagrangeInterpolation>,
                             ArrayOfLagrangeInterpolation>();

  py::class_<GasAbsLookup>(m, "GasAbsLookup")
      .def(py::init<>())
      .PythonInterfaceFileIO(GasAbsLookup)
      .PythonInterfaceBasicRepresentation(GasAbsLookup)
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, species, Species, Species)
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, nonlin_species, NonLinearSpecies, NonLinearSpecies)
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, f_grid, Fgrid, Fgrid)
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, flag_default, FLAGDefault, FLAGDefault)
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, p_grid, Pgrid, Pgrid)
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, log_p_grid, LogPgrid, LogPgrid)
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, vmrs_ref, VMRs, VMRs)
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, t_ref, Tref, Tref)
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, t_pert, Tpert, Tpert)
      .PythonInterfaceBasicReferenceProperty(
          GasAbsLookup, nls_pert, NLSPert, NLSPert)
      .PythonInterfaceBasicReferenceProperty(GasAbsLookup, xsec, Xsec, Xsec);
}
}  // namespace Python