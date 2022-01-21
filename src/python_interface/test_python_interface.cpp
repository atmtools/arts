#include "array.h"
#include "auto_md.h"
#include "debug.h"
#include "matpackIII.h"
#include "ppath.h"
#include "python_interface.h"
#include "xml_io.h"
#include <cstdint>
#include <functional>
#include <pybind11/attr.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <tuple>
#include <vector>

#include "py_macros.h"

namespace Python {
PYBIND11_MODULE(pyarts_classes, m) {

  py::class_<GridPos>(m, "GridPos")
      .def(py::init<>())
      .def(py::init<Index, std::array<Numeric, 2>>())
      .PythonInterfaceFileIO(GridPos)
      .PythonInterfaceBasicRepresentation(GridPos)
      .def_readwrite("idx", &GridPos::idx)
      .def_readwrite("fd", &GridPos::fd);

  py::class_<ArrayOfGridPos>(m, "ArrayOfGridPos")
      .PythonInterfaceArrayDefault(ArrayOfGridPos)
      .PythonInterfaceBasicRepresentation(ArrayOfGridPos);

  py::class_<Ppath>(m, "Ppath")
      .def(py::init<>())
      .def(py::init<Index,
                    Index,
                    Numeric,
                    String,
                    Vector,
                    Vector,
                    Numeric,
                    Matrix,
                    Matrix,
                    Vector,
                    Vector,
                    Vector,
                    Vector,
                    Numeric,
                    Vector,
                    Vector,
                    ArrayOfGridPos,
                    ArrayOfGridPos,
                    ArrayOfGridPos>())
      .PythonInterfaceFileIO(Ppath)
      .def_readwrite("dim", &Ppath::dim)
      .def_readwrite("np", &Ppath::np)
      .def_readwrite("constant", &Ppath::constant)
      .def_readwrite("background", &Ppath::background)
      .def_readwrite("start_pos", &Ppath::start_pos)
      .def_readwrite("start_los", &Ppath::start_los)
      .def_readwrite("start_lstep", &Ppath::start_lstep)
      .def_readwrite("pos", &Ppath::pos)
      .def_readwrite("los", &Ppath::los)
      .def_readwrite("r", &Ppath::r)
      .def_readwrite("lstep", &Ppath::lstep)
      .def_readwrite("end_pos", &Ppath::end_pos)
      .def_readwrite("end_los", &Ppath::end_los)
      .def_readwrite("end_lstep", &Ppath::end_lstep)
      .def_readwrite("nreal", &Ppath::nreal)
      .def_readwrite("ngroup", &Ppath::ngroup)
      .def_readwrite("gp_p", &Ppath::gp_p)
      .def_readwrite("gp_lat", &Ppath::gp_lat)
      .def_readwrite("gp_lon", &Ppath::gp_lon);

  py::class_<ArrayOfPpath>(m, "ArrayOfPpath")
      .PythonInterfaceFileIO(ArrayOfPpath)
      .PythonInterfaceArrayDefault(ArrayOfPpath);

  py::class_<Rational>(m, "Rational")
      .def(py::init<>())
      .def(py::init<Index>())
      .def(py::init<Index, Index>())
      .PythonInterfaceFileIO(Rational)
      .PythonInterfaceBasicRepresentation(Rational)
      .PythonInterfaceInPlaceMathOperators(Rational, Rational)
      .PythonInterfaceInPlaceMathOperators(Rational, Index)
      .PythonInterfaceMathOperators(Rational, Rational)
      .PythonInterfaceMathOperators(Rational, Index);
  py::implicitly_convertible<Index, Rational>();

  py::class_<Time>(m, "Time")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceFileIO(Time)
      .PythonInterfaceBasicRepresentation(Time)
      .def("sec", [](const Time& t){return t.Seconds();});

  py::class_<GriddedField1>(m, "GriddedField1")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceFileIO(GriddedField1)
      .PythonInterfaceBasicRepresentation(GriddedField1)
      .PythonInterfaceGriddedField(GriddedField1);

  py::class_<GriddedField2>(m, "GriddedField2")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceFileIO(GriddedField2)
      .PythonInterfaceBasicRepresentation(GriddedField2)
      .PythonInterfaceGriddedField(GriddedField2);

  py::class_<GriddedField3>(m, "GriddedField3")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceFileIO(GriddedField3)
      .PythonInterfaceBasicRepresentation(GriddedField3)
      .PythonInterfaceGriddedField(GriddedField3);

  py::class_<GriddedField4>(m, "GriddedField4")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceFileIO(GriddedField4)
      .PythonInterfaceBasicRepresentation(GriddedField4)
      .PythonInterfaceGriddedField(GriddedField4);

  py::class_<GriddedField5>(m, "GriddedField5")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceFileIO(GriddedField5)
      .PythonInterfaceBasicRepresentation(GriddedField5)
      .PythonInterfaceGriddedField(GriddedField5);

  py::class_<GriddedField6>(m, "GriddedField6")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceFileIO(GriddedField6)
      .PythonInterfaceBasicRepresentation(GriddedField6)
      .PythonInterfaceGriddedField(GriddedField6);

  py::class_<TransmissionMatrix>(m, "TransmissionMatrix")
      .def(py::init<>())
      .def(py::init<Index>())
      .def(py::init<Index, Index>())
      .PythonInterfaceFileIO(TransmissionMatrix)
      .PythonInterfaceBasicRepresentation(TransmissionMatrix)
      .def("__getitem__",
           [](const TransmissionMatrix& t, Index i) {
             ARTS_USER_ERROR_IF(i < 0 or i > t.Frequencies(), "Out of bounds")
             return t.Mat(i);
           })
      .def("__setitem__",
           [](TransmissionMatrix& t, Index i, const py::array_t<Numeric>& arr) {
             ARTS_USER_ERROR_IF(i < 0 or i > t.Frequencies(), "Out of bounds")
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

  py::class_<TessemNN>(m, "TessemNN")
      .def(py::init<>())
      .PythonInterfaceFileIO(TessemNN)
      .def_readwrite("nb_inputs", &TessemNN::nb_inputs)
      .def_readwrite("nb_outputs", &TessemNN::nb_outputs)
      .def_readwrite("nb_cache", &TessemNN::nb_cache)
      .def_readwrite("b1", &TessemNN::b1)
      .def_readwrite("b2", &TessemNN::b2)
      .def_readwrite("w1", &TessemNN::w1)
      .def_readwrite("w2", &TessemNN::w2)
      .def_readwrite("x_min", &TessemNN::x_min)
      .def_readwrite("x_max", &TessemNN::x_max)
      .def_readwrite("y_min", &TessemNN::y_min)
      .def_readwrite("y_max", &TessemNN::y_max);

  py::class_<Timer>(m, "Timer")
      .def(py::init<>())
      .PythonInterfaceFileIO(Timer)
      .def_readwrite("running", &Timer::running)
      .def_readwrite("finished", &Timer::finished)
#ifdef TIME_SUPPORT
      .def_readwrite("tms cputime_start", &Timer::cputime_start)
      .def_readwrite("realtime_start", &Timer::realtime_start)
      .def_readwrite("tms cputime_end", &Timer::cputime_end)
      .def_readwrite("realtime_end", &Timer::realtime_end)
#endif
      ;
  
  py::class_<TelsemAtlas>(m, "TelsemAtlas")
      .def(py::init<>())
      .PythonInterfaceFileIO(TelsemAtlas)
      .PythonInterfaceBasicRepresentation(TelsemAtlas);

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
             ARTS_USER_ERROR_IF(arr_val.shape(0) not_eq t.StokesDim(), "Bad input size!")
             const Index n = t.StokesDim();

             switch (n) {
               case 1:
                 for (Index r = 0; r < n; r++)
                     t.Vec1(i)[r] = *arr_val.data(r);
                 break;
               case 2:
                 for (Index r = 0; r < n; r++)
                     t.Vec2(i)[r] = *arr_val.data(r);
                 break;
               case 3:
                 for (Index r = 0; r < n; r++)
                     t.Vec3(i)[r] = *arr_val.data(r);
                 break;
               case 4:
                 for (Index r = 0; r < n; r++)
                     t.Vec4(i)[r] = *arr_val.data(r);
                 break;
             }
           });

  py::class_<StokesVector>(m, "StokesVector")
      .def(py::init([](Index nf, Index ns, Index nza, Index naa, Numeric v) {
             ARTS_USER_ERROR_IF(ns < 1 or ns > 4, "Bad Stokes Dimension ", ns, " should be [1, 4]")
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
                               s.Data().nbooks(), ", ",
                               s.Data().npages(), ", ",
                               s.Data().nrows(), ", ",
                               s.Data().ncols(), ')',
                               "\nGot shape: (",
                               d.nbooks(), ", ",
                               d.npages(), ", ",
                               d.nrows(), ", ",
                               d.ncols(), ')')
            s.Data() = d;
          });

  py::class_<PropagationMatrix>(m, "PropagationMatrix")
      .def(py::init([](Index nf, Index ns, Index nza, Index naa, Numeric v) {
             ARTS_USER_ERROR_IF(ns < 1 or ns > 4, "Bad Stokes Dimension ", ns, " should be [1, 4]")
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
                               s.Data().nbooks(), ", ",
                               s.Data().npages(), ", ",
                               s.Data().nrows(), ", ",
                               s.Data().ncols(), ')')
            s.Data() = d;
          });

  PythonInterfaceWorkspaceArray(ArrayOfIndex);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfIndex);
  PythonInterfaceWorkspaceArray(ArrayOfTime);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfTime);
  PythonInterfaceWorkspaceArray(ArrayOfString);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfString);
  PythonInterfaceWorkspaceArray(ArrayOfVector);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfVector);
  PythonInterfaceWorkspaceArray(ArrayOfMatrix);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfMatrix);
  PythonInterfaceWorkspaceArray(ArrayOfTensor3);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfTensor3);
  PythonInterfaceWorkspaceArray(ArrayOfTensor4);
  PythonInterfaceWorkspaceArray(ArrayOfTensor5);
  PythonInterfaceWorkspaceArray(ArrayOfTensor6);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfTensor6);
  PythonInterfaceWorkspaceArray(ArrayOfTensor7);
  PythonInterfaceWorkspaceArray(ArrayOfGriddedField1);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfGriddedField1);
  PythonInterfaceWorkspaceArray(ArrayOfGriddedField2);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfGriddedField2);
  PythonInterfaceWorkspaceArray(ArrayOfGriddedField3);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfGriddedField3);
  PythonInterfaceWorkspaceArray(ArrayOfGriddedField4);
  PythonInterfaceWorkspaceArray(ArrayOfTransmissionMatrix);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfTransmissionMatrix);
  PythonInterfaceWorkspaceArray(ArrayOfRadiationVector);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfRadiationVector);
  PythonInterfaceWorkspaceArray(ArrayOfStokesVector);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfStokesVector);
  PythonInterfaceWorkspaceArray(ArrayOfPropagationMatrix);
  PythonInterfaceWorkspaceArray(ArrayOfArrayOfPropagationMatrix);
}
}  // namespace Python

#undef PythonInterfaceFileIO
#undef PythonInterfaceIndexItemAccess
#undef PythonInterfaceBasicRepresentation
#undef PythonInterfaceArrayDefault
