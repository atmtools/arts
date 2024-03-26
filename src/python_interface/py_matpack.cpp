#include "py_matpack.h"

#include <pybind11/cast.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include <python_interface.h>

#include <memory>
#include <stdexcept>

#include "configtypes.h"

namespace Python {
void py_matpack(py::module_& m) try {
  fix_matpack_constant_data(py_staticVector2(m));
  fix_matpack_constant_data(py_staticVector3(m));

  artsclass<Range>(m, "Range")
      .def(py::init([](Index a, Index b, Index c) {
             ARTS_USER_ERROR_IF(0 > a, "Bad offset")
             ARTS_USER_ERROR_IF(0 > b, "Bad extent")
             return std::make_shared<Range>(a, b, c);
           }),
           py::arg("offset"),
           py::arg("extent"),
           py::arg("stride") = 1,
           "Valued initialization")
      .PythonInterfaceBasicRepresentation(Range)
      .def(py::pickle(
          [](const Range& self) {
            return py::make_tuple(self.offset, self.extent, self.stride);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 3, "Invalid state!")
            return std::make_shared<Range>(
                t[0].cast<Index>(), t[1].cast<Index>(), t[2].cast<Index>());
          }))
      .doc() = "A range, used to select parts of a matpack type";

  artsclass<IndexVector>(m, "IndexVector", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<const IndexVector&>())
      .def(py::init([](const numpy_array& x) { return copy<1, Index>(x); }),
           py::arg("vec").none(false),
           py::doc("From :class:`numpy.ndarray` equivalent"))
      .def(py::init([](const py::list& x) {
             return py::cast<IndexVector>(x.cast<numpy_array>());
           }),
           py::arg("lst"),
           py::doc("From :class:`list` equivalent via numpy"))
      .PythonInterfaceCopyValue(IndexVector)
      .PythonInterfaceBasicRepresentation(IndexVector)
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](IndexVector& x) -> py::buffer_info {
        return py::buffer_info(x.data_handle(),
                               sizeof(Index),
                               py::format_descriptor<Index>::format(),
                               1,
                               {x.size()},
                               {sizeof(Index)});
      })
      .def_property(
          "value",
          py::cpp_function(
              [](IndexVector& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](IndexVector& x, IndexVector& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"));
  py::implicitly_convertible<numpy_array, IndexVector>();
  py::implicitly_convertible<py::list, IndexVector>();

  py_staticVector(m)
      .def(py::init([](const numpy_array& x) { return copy<1, Numeric>(x); }),
           py::arg("vec").none(false),
           py::doc("From :class:`numpy.ndarray` equivalent"))
      .def(py::init([](const py::list& x) {
             return py::cast<Vector>(x.cast<numpy_array>());
           }),
           py::arg("lst"),
           py::doc("From :class:`list` equivalent via numpy"))
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](Vector& x) -> py::buffer_info {
        return py::buffer_info(x.data_handle(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               1,
                               {x.size()},
                               {sizeof(Numeric)});
      })
      .def_property(
          "value",
          py::cpp_function(
              [](Vector& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](Vector& x, Vector& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"));

  py_staticMatrix(m)
      .def(py::init([](const numpy_array& x) { return copy<2, Numeric>(x); }),
           py::arg("mat").none(false),
           py::doc("From :class:`numpy.ndarray` equivalent"))
      .def(py::init([](const py::list& x) {
             return py::cast<Matrix>(x.cast<numpy_array>());
           }),
           py::arg("lst"),
           py::doc("From :class:`list` equivalent via numpy"))
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](Matrix& x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            2,
            {static_cast<ssize_t>(x.nrows()), static_cast<ssize_t>(x.ncols())},
            {sizeof(Numeric) * static_cast<ssize_t>(x.ncols()),
             sizeof(Numeric)});
      })
      .def_property_readonly(
          "T",
          [](const Matrix& x) { return Matrix(transpose(x)); },
          "Non-trivial transpose")
      .def_property(
          "value",
          py::cpp_function(
              [](Matrix& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](Matrix& x, Matrix& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"));

  py_staticTensor3(m)
      .def(py::init([](const numpy_array& x) { return copy<3, Numeric>(x); }),
           py::arg("ten").none(false),
           py::doc("From :class:`numpy.ndarray` equivalent"))
      .def(py::init([](const py::list& x) {
             return py::cast<Tensor3>(x.cast<numpy_array>());
           }),
           py::arg("lst"),
           py::doc("From :class:`list` equivalent via numpy"))
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](Tensor3& x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            3,
            {static_cast<ssize_t>(x.npages()),
             static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols())},
            {sizeof(Numeric) * static_cast<ssize_t>(x.nrows() * x.ncols()),
             sizeof(Numeric) * static_cast<ssize_t>(x.ncols()),
             sizeof(Numeric)});
      })
      .def_property(
          "value",
          py::cpp_function(
              [](Tensor3& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](Tensor3& x, Tensor3& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"))
      .def(py::pickle(
          [](const py::object& self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<Tensor3>()(t[0]).cast<Tensor3>();
          }));

  py_staticTensor4(m)
      .def(py::init([](const numpy_array& x) { return copy<4, Numeric>(x); }),
           py::arg("ten").none(false),
           py::doc("From :class:`numpy.ndarray` equivalent"))
      .def(py::init([](const py::list& x) {
             return py::cast<Tensor4>(x.cast<numpy_array>());
           }),
           py::arg("lst"),
           py::doc("From :class:`list` equivalent via numpy"))
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](Tensor4& x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            4,
            {static_cast<ssize_t>(x.nbooks()),
             static_cast<ssize_t>(x.npages()),
             static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols())},
            {sizeof(Numeric) *
                 static_cast<ssize_t>(x.npages() * x.ncols() * x.nrows()),
             sizeof(Numeric) * static_cast<ssize_t>(x.ncols() * x.nrows()),
             sizeof(Numeric) * static_cast<ssize_t>(x.ncols()),
             sizeof(Numeric)});
      })
      .def_property(
          "value",
          py::cpp_function(
              [](Tensor4& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](Tensor4& x, Tensor4& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"));

  py_staticTensor5(m)
      .def(py::init([](const numpy_array& x) { return copy<5, Numeric>(x); }),
           py::arg("ten").none(false),
           py::doc("From :class:`numpy.ndarray` equivalent"))
      .def(py::init([](const py::list& x) {
             return py::cast<Tensor5>(x.cast<numpy_array>());
           }),
           py::arg("lst"),
           py::doc("From :class:`list` equivalent via numpy"))
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](Tensor5& x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            5,
            {static_cast<ssize_t>(x.nshelves()),
             static_cast<ssize_t>(x.nbooks()),
             static_cast<ssize_t>(x.npages()),
             static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols())},
            {sizeof(Numeric) * static_cast<ssize_t>(x.nbooks() * x.npages() *
                                                    x.ncols() * x.nrows()),
             sizeof(Numeric) *
                 static_cast<ssize_t>(x.npages() * x.ncols() * x.nrows()),
             sizeof(Numeric) * static_cast<ssize_t>(x.ncols() * x.nrows()),
             sizeof(Numeric) * static_cast<ssize_t>(x.ncols()),
             sizeof(Numeric)});
      })
      .def_property(
          "value",
          py::cpp_function(
              [](Tensor5& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](Tensor5& x, Tensor5& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"));

  py_staticTensor6(m)
      .def(py::init([](const numpy_array& x) { return copy<6, Numeric>(x); }),
           py::arg("ten").none(false),
           py::doc("From :class:`numpy.ndarray` equivalent"))
      .def(py::init([](const py::list& x) {
             return py::cast<Tensor6>(x.cast<numpy_array>());
           }),
           py::arg("lst"),
           py::doc("From :class:`list` equivalent via numpy"))
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](Tensor6& x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            6,
            {static_cast<ssize_t>(x.nvitrines()),
             static_cast<ssize_t>(x.nshelves()),
             static_cast<ssize_t>(x.nbooks()),
             static_cast<ssize_t>(x.npages()),
             static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols())},
            {sizeof(Numeric) *
                 static_cast<ssize_t>(x.nshelves() * x.nbooks() * x.npages() *
                                      x.ncols() * x.nrows()),
             sizeof(Numeric) * static_cast<ssize_t>(x.nbooks() * x.npages() *
                                                    x.ncols() * x.nrows()),
             sizeof(Numeric) *
                 static_cast<ssize_t>(x.npages() * x.ncols() * x.nrows()),
             sizeof(Numeric) * static_cast<ssize_t>(x.ncols() * x.nrows()),
             sizeof(Numeric) * static_cast<ssize_t>(x.ncols()),
             sizeof(Numeric)});
      })
      .def_property(
          "value",
          py::cpp_function(
              [](Tensor6& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](Tensor6& x, Tensor6& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"));

  py_staticTensor7(m)
      .def(py::init([](const numpy_array& x) { return copy<7, Numeric>(x); }),
           py::arg("ten").none(false),
           py::doc("From :class:`numpy.ndarray` equivalent"))
      .def(py::init([](const py::list& x) {
             return py::cast<Tensor7>(x.cast<numpy_array>());
           }),
           py::arg("lst"),
           py::doc("From :class:`list` equivalent via numpy"))
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](Tensor7& x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            7,
            {static_cast<ssize_t>(x.nlibraries()),
             static_cast<ssize_t>(x.nvitrines()),
             static_cast<ssize_t>(x.nshelves()),
             static_cast<ssize_t>(x.nbooks()),
             static_cast<ssize_t>(x.npages()),
             static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols())},
            {sizeof(Numeric) * static_cast<ssize_t>(
                                   x.nvitrines() * x.nshelves() * x.nbooks() *
                                   x.npages() * x.ncols() * x.nrows()),
             sizeof(Numeric) *
                 static_cast<ssize_t>(x.nshelves() * x.nbooks() * x.npages() *
                                      x.ncols() * x.nrows()),
             sizeof(Numeric) * static_cast<ssize_t>(x.nbooks() * x.npages() *
                                                    x.ncols() * x.nrows()),
             sizeof(Numeric) *
                 static_cast<ssize_t>(x.npages() * x.ncols() * x.nrows()),
             sizeof(Numeric) * static_cast<ssize_t>(x.ncols() * x.nrows()),
             sizeof(Numeric) * static_cast<ssize_t>(x.ncols()),
             sizeof(Numeric)});
      })
      .def_property(
          "value",
          py::cpp_function(
              [](Tensor7& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](Tensor7& x, Tensor7& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"));

  py::implicitly_convertible<numpy_array, Vector>();
  py::implicitly_convertible<numpy_array, Matrix>();
  py::implicitly_convertible<numpy_array, Tensor3>();
  py::implicitly_convertible<numpy_array, Tensor4>();
  py::implicitly_convertible<numpy_array, Tensor5>();
  py::implicitly_convertible<numpy_array, Tensor6>();
  py::implicitly_convertible<numpy_array, Tensor7>();
  py::implicitly_convertible<py::list, Vector>();
  py::implicitly_convertible<py::list, Matrix>();
  py::implicitly_convertible<py::list, Tensor3>();
  py::implicitly_convertible<py::list, Tensor4>();
  py::implicitly_convertible<py::list, Tensor5>();
  py::implicitly_convertible<py::list, Tensor6>();
  py::implicitly_convertible<py::list, Tensor7>();

  py_staticRational(m)
      .def(py::init([](Index n, Index d) {
             ARTS_USER_ERROR_IF(d < 1, "Must be positive")
             return Rational(n, d);
           }),
           py::arg("n") = 0,
           py::arg("d") = 1,
           py::doc("Default rational"))
      .def(py::init(
               [](const String& s) { return std::make_shared<Rational>(s); }),
           "From :class:`str`")
      .def(py::init([](Numeric n) { return Rational(std::to_string(n)); }),
           "From :class:`float`")
      .def("__float__", [](const Rational& x) { return Numeric(x); })
      .def("__int__", [](const Rational& x) { return Index(x); })
      .PythonInterfaceCommonMath(Index)
      .PythonInterfaceCommonMathSelf
      .def_readwrite("n", &Rational::numer, ":class:`int` Numerator")
      .def_readwrite("d", &Rational::denom, ":class:`int` Denominator")
      .def(py::pickle(
          [](const Rational& self) {
            return py::make_tuple(self.numer, self.denom);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 2, "Invalid state!")

            return std::make_shared<Rational>(t[0].cast<Index>(),
                                              t[1].cast<Index>());
          }));
  py::implicitly_convertible<Index, Rational>();

  artsclass<ComplexVector>(m, "ComplexVector", py::buffer_protocol())
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_property(
          "value",
          py::cpp_function(
              [](ComplexVector& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](ComplexVector& x, ComplexVector& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"))
      .def_buffer([](ComplexVector& x) -> py::buffer_info {
        return py::buffer_info(x.data_handle(),
                               sizeof(Complex),
                               py::format_descriptor<Complex>::format(),
                               1,
                               {static_cast<ssize_t>(x.nelem())},
                               {sizeof(Complex)});
      })
      .doc() = R"--(Holds complex vector data.
)--";

  artsclass<ComplexMatrix>(m, "ComplexMatrix", py::buffer_protocol())
      .def(py::init([]() { return std::make_shared<ComplexMatrix>(); }),
           "Default matrix")
      .def(py::init([](const numpy_array& x) { return copy<2, Complex>(x); }))
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_property(
          "value",
          py::cpp_function(
              [](ComplexMatrix& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](ComplexMatrix& x, ComplexMatrix& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"))
      .def_buffer([](ComplexMatrix& x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Complex),
            py::format_descriptor<Complex>::format(),
            2,
            {static_cast<ssize_t>(x.nrows()), static_cast<ssize_t>(x.ncols())},
            {sizeof(Complex) * static_cast<ssize_t>(x.ncols()),
             sizeof(Complex)});
      })
      .def(py::pickle(
          [](const py::object& self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<ComplexMatrix>()(t[0]).cast<ComplexMatrix>();
          }))
      .doc() = R"--(Holds complex matrix data.
)--";

  artsclass<ComplexTensor3>(m, "ComplexTensor3", py::buffer_protocol())
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_property(
          "value",
          py::cpp_function(
              [](ComplexTensor3& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](ComplexTensor3& x, ComplexTensor3& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"))
      .def_buffer([](ComplexTensor3& x) -> py::buffer_info {
        return py::buffer_info(
            x.data_handle(),
            sizeof(Complex),
            py::format_descriptor<Complex>::format(),
            3,
            {static_cast<ssize_t>(x.npages()),
             static_cast<ssize_t>(x.nrows()),
             static_cast<ssize_t>(x.ncols())},
            {sizeof(Complex) * static_cast<ssize_t>(x.nrows() * x.ncols()),
             sizeof(Complex) * static_cast<ssize_t>(x.ncols()),
             sizeof(Complex)});
      })
      .doc() = R"--(Holds complex tensor3 data.
)--";

  py_staticAscendingGrid(m)
      .def(py::init<Vector>(), py::arg("x"), "From :class:`Vector`")
      .def(py::init([](const numpy_array& x) -> AscendingGrid {
             return AscendingGrid{*copy<1, Numeric>(x)};
           }),
           py::arg("vec").none(false),
           py::doc("From :class:`numpy.ndarray` equivalent"))
      .def(py::init([](const py::list& x) {
             return py::cast<AscendingGrid>(x.cast<numpy_array>());
           }),
           py::arg("lst"),
           py::doc("From :class:`list` equivalent via numpy"))
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](AscendingGrid& x) -> py::buffer_info {
        return py::buffer_info(x.data_handle(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               1,
                               {x.size()},
                               {sizeof(Numeric)},
                               true);
      })
      .def_property(
          "value",
          py::cpp_function(
              [](AscendingGrid& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](AscendingGrid& x, Vector& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"));

  py::implicitly_convertible<numpy_array, AscendingGrid>();
  py::implicitly_convertible<py::list, AscendingGrid>();
  py::implicitly_convertible<Vector, AscendingGrid>();

  py_staticVector(m).def(py::init([](const AscendingGrid& x) {
                           return std::make_shared<Vector>(x);
                         }),
                         py::arg("grid"),
                         py::doc("From :class:`~pyarts.arts.AscendingGrid`"));
  py::implicitly_convertible<AscendingGrid, Vector>();
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize matpack\n", e.what()));
}
}  // namespace Python
