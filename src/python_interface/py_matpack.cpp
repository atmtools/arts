#include <matpackVII.h>
#include <pybind11/attr.h>
#include <xml_io.h>

#include <functional>
#include <stdexcept>

#include "debug.h"
#include "matpackI.h"
#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_matpack(py::module_& m) {
  py::class_<Vector>(m, "Vector", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index>())
      .def(py::init<Index, Numeric>())
      .def(py::init<const std::vector<Numeric>&>())
      .def(py::init<Numeric, Index, Numeric>())
      .def_property_readonly(
          "shape",
          [](const Vector& x) { return std::array{x.nelem()}; },
          "The shape of the data")
      .PythonInterfaceBasicRepresentation(Vector)
      .PythonInterfaceFileIO(Vector)
      .PythonInterfaceIndexItemAccess(Vector)
      .def_buffer([](Vector& x) -> py::buffer_info {
        return py::buffer_info(x.get_c_array(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               1,
                               {x.nelem()},
                               {sizeof(Numeric)});
      })
      .doc() =
      "The Arts Vector class\n"
      "\n"
      "This class is compatible with numpy arrays.  The data can "
      "be accessed without copy using np.array(x, copy=False), "
      "with x as an instance of this class.  Note that access to "
      "any and all of the mathematical operations are only available "
      "via the numpy interface.  Also note that the equality operation "
      "only checks pointer equality and not element-wise equality\n"
      "\n"
      "Initialization:\n"
      "    Vector(): for basic initialization\n\n"
      "    Vector(Index): for constant size, unknown value\n\n"
      "    Vector(Index, Numeric): for constant size, constant value\n\n"
      "    Vector(List or Array): to copy elements\n\n"
      "    Vector(Numeric x, Index n, Numeric dx): as Vector([x+i*dx for i in range(n)])\n\n";

  py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index, Index>())
      .def(py::init<Index, Index, Numeric>())
      .def(py::init([](const py::array_t<Numeric>& arr) {
        ARTS_USER_ERROR_IF(arr.request().ndim not_eq 2, "Bad size array")

        auto arr_val = arr.unchecked<2>();

        Matrix out(arr_val.shape(0), arr_val.shape(1));

        for (Index r = 0; r < out.nrows(); r++) {
          for (Index c = 0; c < out.ncols(); c++) {
            out(r, c) = arr_val(r, c);
          }
        }

        return out;
      }))
      .PythonInterfaceBasicRepresentation(Matrix)
      .PythonInterfaceFileIO(Matrix)
      .def_property_readonly(
          "shape",
          [](const Matrix& x) {
            return std::array{x.nrows(), x.ncols()};
          },
          "The shape of the data")
      .def("__getitem__",
           [](const Matrix& x, std::tuple<Index, Index> inds) {
             auto [r, c] = inds;
             if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0)
               throw std::out_of_range("Out of bounds");
             return x(r, c);
           })
      .def("__setitem__",
           [](Matrix& x, std::tuple<Index, Index> inds, Numeric y) {
             auto [r, c] = inds;
             if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0)
               throw std::out_of_range("Out of bounds");
             x(r, c) = y;
           })
      .def_buffer([](Matrix& x) -> py::buffer_info {
        return py::buffer_info(x.get_c_array(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               2,
                               {x.nrows(), x.ncols()},
                               {sizeof(Numeric) * x.ncols(), sizeof(Numeric)});
      })
      .doc() =
      "The Arts Matrix class\n"
      "\n"
      "This class is compatible with numpy arrays.  The data can "
      "be accessed without copy using np.array(x, copy=False), "
      "with x as an instance of this class.  Note that access to "
      "any and all of the mathematical operations are only available "
      "via the numpy interface.  Also note that the equality operation "
      "only checks pointer equality and not element-wise equality\n"
      "\n"
      "Initialization:\n"
      "    Matrix(): for basic initialization\n\n"
      "    Matrix(Index, Index): for constant size, unknown value\n\n"
      "    Matrix(Index, Index, Numeric): for constant size, constant value\n\n"
      "    Matrix(List or Array): to copy elements\n\n";

  py::class_<Tensor3>(m, "Tensor3", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index, Index, Index>())
      .def(py::init<Index, Index, Index, Numeric>())
      .def(py::init([](const py::array_t<Numeric>& arr) {
        ARTS_USER_ERROR_IF(arr.request().ndim not_eq 3, "Bad size array")

        auto arr_val = arr.unchecked<3>();

        Tensor3 out(arr_val.shape(0), arr_val.shape(1), arr_val.shape(2));

        for (Index p = 0; p < out.npages(); p++) {
          for (Index r = 0; r < out.nrows(); r++) {
            for (Index c = 0; c < out.ncols(); c++) {
              out(p, r, c) = arr_val(p, r, c);
            }
          }
        }

        return out;
      }))
      .PythonInterfaceBasicRepresentation(Tensor3)
      .PythonInterfaceFileIO(Tensor3)
      .def_property_readonly(
          "shape",
          [](const Tensor3& x) {
            return std::array{x.npages(), x.nrows(), x.ncols()};
          },
          "The shape of the data")
      .def(
          "__getitem__",
          [](const Tensor3& x, std::tuple<Index, Index, Index> inds) {
            auto [p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0)
              throw std::out_of_range("Out of bounds");
            return x(p, r, c);
          },
          py::return_value_policy::reference_internal)
      .def("__setitem__",
           [](Tensor3& x, std::tuple<Index, Index, Index> inds, Numeric y) {
             auto [p, r, c] = inds;
             if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0)
               throw std::out_of_range("Out of bounds");
             x(p, r, c) = y;
           })
      .def_buffer([](Tensor3& x) -> py::buffer_info {
        return py::buffer_info(x.get_c_array(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               3,
                               {x.npages(), x.nrows(), x.ncols()},
                               {sizeof(Numeric) * x.nrows() * x.ncols(),
                                sizeof(Numeric) * x.ncols(),
                                sizeof(Numeric)});
      })
      .doc() =
      "The Arts Tensor3 class\n"
      "\n"
      "This class is compatible with numpy arrays.  The data can "
      "be accessed without copy using np.array(x, copy=False), "
      "with x as an instance of this class.  Note that access to "
      "any and all of the mathematical operations are only available "
      "via the numpy interface.  Also note that the equality operation "
      "only checks pointer equality and not element-wise equality\n"
      "\n"
      "Initialization:\n"
      "    Tensor3(): for basic initialization\n\n"
      "    Tensor3(Index, Index, Index): for constant size, unknown value\n\n"
      "    Tensor3(Index, Index, Index, Numeric): for constant size, constant value\n\n"
      "    Tensor3(List or Array): to copy elements\n\n";

  py::class_<Tensor4>(m, "Tensor4", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index, Index, Index, Index>())
      .def(py::init<Index, Index, Index, Index, Numeric>())
      .def(py::init([](const py::array_t<Numeric>& arr) {
        ARTS_USER_ERROR_IF(arr.request().ndim not_eq 4, "Bad size array")

        auto arr_val = arr.unchecked<4>();

        Tensor4 out(arr_val.shape(0),
                    arr_val.shape(1),
                    arr_val.shape(2),
                    arr_val.shape(3));

        for (Index b = 0; b < out.nbooks(); b++) {
          for (Index p = 0; p < out.npages(); p++) {
            for (Index r = 0; r < out.nrows(); r++) {
              for (Index c = 0; c < out.ncols(); c++) {
                out(b, p, r, c) = arr_val(b, p, r, c);
              }
            }
          }
        }

        return out;
      }))
      .PythonInterfaceBasicRepresentation(Tensor4)
      .PythonInterfaceFileIO(Tensor4)
      .def_property_readonly(
          "shape",
          [](const Tensor4& x) {
            return std::array{x.nbooks(), x.npages(), x.nrows(), x.ncols()};
          },
          "The shape of the data")
      .def(
          "__getitem__",
          [](const Tensor4& x, std::tuple<Index, Index, Index, Index> inds) {
            auto [b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0)
              throw std::out_of_range("Out of bounds");
            return x(b, p, r, c);
          },
          py::return_value_policy::reference_internal)
      .def("__setitem__",
           [](Tensor4& x,
              std::tuple<Index, Index, Index, Index> inds,
              Numeric y) {
             auto [b, p, r, c] = inds;
             if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0)
               throw std::out_of_range("Out of bounds");
             x(b, p, r, c) = y;
           })
      .def_buffer([](Tensor4& x) -> py::buffer_info {
        return py::buffer_info(
            x.get_c_array(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            4,
            {x.nbooks(), x.npages(), x.nrows(), x.ncols()},
            {sizeof(Numeric) * x.npages() * x.ncols() * x.nrows(),
             sizeof(Numeric) * x.ncols() * x.nrows(),
             sizeof(Numeric) * x.ncols(),
             sizeof(Numeric)});
      })
      .doc() =
      "The Arts Tensor4 class\n"
      "\n"
      "This class is compatible with numpy arrays.  The data can "
      "be accessed without copy using np.array(x, copy=False), "
      "with x as an instance of this class.  Note that access to "
      "any and all of the mathematical operations are only available "
      "via the numpy interface.  Also note that the equality operation "
      "only checks pointer equality and not element-wise equality\n"
      "\n"
      "Initialization:\n"
      "    Tensor4(): for basic initialization\n\n"
      "    Tensor4(Index, Index, Index, Index): for constant size, unknown value\n\n"
      "    Tensor4(Index, Index, Index, Index, Numeric): for constant size, constant value\n\n"
      "    Tensor4(List or Array): to copy elements\n\n";

  py::class_<Tensor5>(m, "Tensor5", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index, Index, Index, Index, Index>())
      .def(py::init<Index, Index, Index, Index, Index, Numeric>())
      .def(py::init([](const py::array_t<Numeric>& arr) {
        ARTS_USER_ERROR_IF(arr.request().ndim not_eq 5, "Bad size array")

        auto arr_val = arr.unchecked<5>();

        Tensor5 out(arr_val.shape(0),
                    arr_val.shape(1),
                    arr_val.shape(2),
                    arr_val.shape(3),
                    arr_val.shape(4));

        for (Index s = 0; s < out.nshelves(); s++) {
          for (Index b = 0; b < out.nbooks(); b++) {
            for (Index p = 0; p < out.npages(); p++) {
              for (Index r = 0; r < out.nrows(); r++) {
                for (Index c = 0; c < out.ncols(); c++) {
                  out(s, b, p, r, c) = arr_val(s, b, p, r, c);
                }
              }
            }
          }
        }

        return out;
      }))
      .PythonInterfaceBasicRepresentation(Tensor5)
      .PythonInterfaceFileIO(Tensor5)
      .def_property_readonly(
          "shape",
          [](const Tensor5& x) {
            return std::array{
                x.nshelves(), x.nbooks(), x.npages(), x.nrows(), x.ncols()};
          },
          "The shape of the data")
      .def(
          "__getitem__",
          [](const Tensor5& x,
             std::tuple<Index, Index, Index, Index, Index> inds) {
            auto [s, b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                x.nshelves() <= s or s < 0)
              throw std::out_of_range("Out of bounds");
            return x(s, b, p, r, c);
          },
          py::return_value_policy::reference_internal)
      .def("__setitem__",
           [](Tensor5& x,
              std::tuple<Index, Index, Index, Index, Index> inds,
              Numeric y) {
             auto [s, b, p, r, c] = inds;
             if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                 x.nshelves() <= s or s < 0)
               throw std::out_of_range("Out of bounds");
             x(s, b, p, r, c) = y;
           })
      .def_buffer([](Tensor5& x) -> py::buffer_info {
        return py::buffer_info(
            x.get_c_array(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            5,
            {x.nshelves(), x.nbooks(), x.npages(), x.nrows(), x.ncols()},
            {sizeof(Numeric) * x.nbooks() * x.npages() * x.ncols() * x.nrows(),
             sizeof(Numeric) * x.npages() * x.ncols() * x.nrows(),
             sizeof(Numeric) * x.ncols() * x.nrows(),
             sizeof(Numeric) * x.ncols(),
             sizeof(Numeric)});
      })
      .doc() =
      "The Arts Tensor5 class\n"
      "\n"
      "This class is compatible with numpy arrays.  The data can "
      "be accessed without copy using np.array(x, copy=False), "
      "with x as an instance of this class.  Note that access to "
      "any and all of the mathematical operations are only available "
      "via the numpy interface.  Also note that the equality operation "
      "only checks pointer equality and not element-wise equality\n"
      "\n"
      "Initialization:\n"
      "    Tensor5(): for basic initialization\n\n"
      "    Tensor5(Index, Index, Index, Index, Index): for constant size, unknown value\n\n"
      "    Tensor5(Index, Index, Index, Index, Index, Numeric): for constant size, constant value\n\n"
      "    Tensor5(List or Array): to copy elements\n\n";

  py::class_<Tensor6>(m, "Tensor6", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index, Index, Index, Index, Index, Index>())
      .def(py::init<Index, Index, Index, Index, Index, Index, Numeric>())
      .def(py::init([](const py::array_t<Numeric>& arr) {
        ARTS_USER_ERROR_IF(arr.request().ndim not_eq 6, "Bad size array")

        auto arr_val = arr.unchecked<6>();

        Tensor6 out(arr_val.shape(0),
                    arr_val.shape(1),
                    arr_val.shape(2),
                    arr_val.shape(3),
                    arr_val.shape(4),
                    arr_val.shape(5));

        for (Index v = 0; v < out.nvitrines(); v++) {
          for (Index s = 0; s < out.nshelves(); s++) {
            for (Index b = 0; b < out.nbooks(); b++) {
              for (Index p = 0; p < out.npages(); p++) {
                for (Index r = 0; r < out.nrows(); r++) {
                  for (Index c = 0; c < out.ncols(); c++) {
                    out(v, s, b, p, r, c) = arr_val(v, s, b, p, r, c);
                  }
                }
              }
            }
          }
        }

        return out;
      }))
      .PythonInterfaceBasicRepresentation(Tensor6)
      .PythonInterfaceFileIO(Tensor6)
      .def_property_readonly(
          "shape",
          [](const Tensor6& x) {
            return std::array{x.nvitrines(),
                              x.nshelves(),
                              x.nbooks(),
                              x.npages(),
                              x.nrows(),
                              x.ncols()};
          },
          "The shape of the data")
      .def(
          "__getitem__",
          [](const Tensor6& x,
             std::tuple<Index, Index, Index, Index, Index, Index> inds) {
            auto [v, s, b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                x.nshelves() <= s or s < 0 or x.nvitrines() <= v or v < 0)
              throw std::out_of_range("Out of bounds");
            return x(v, s, b, p, r, c);
          },
          py::return_value_policy::reference_internal)
      .def("__setitem__",
           [](Tensor6& x,
              std::tuple<Index, Index, Index, Index, Index, Index> inds,
              Numeric y) {
             auto [v, s, b, p, r, c] = inds;
             if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                 x.nshelves() <= s or s < 0 or x.nvitrines() <= v or v < 0)
               throw std::out_of_range("Out of bounds");
             x(v, s, b, p, r, c) = y;
           })
      .def_buffer([](Tensor6& x) -> py::buffer_info {
        return py::buffer_info(
            x.get_c_array(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            6,
            {x.nvitrines(),
             x.nshelves(),
             x.nbooks(),
             x.npages(),
             x.nrows(),
             x.ncols()},
            {sizeof(Numeric) * x.nshelves() * x.nbooks() * x.npages() *
                 x.ncols() * x.nrows(),
             sizeof(Numeric) * x.nbooks() * x.npages() * x.ncols() * x.nrows(),
             sizeof(Numeric) * x.npages() * x.ncols() * x.nrows(),
             sizeof(Numeric) * x.ncols() * x.nrows(),
             sizeof(Numeric) * x.ncols(),
             sizeof(Numeric)});
      })
      .doc() =
      "The Arts Tensor6 class\n"
      "\n"
      "This class is compatible with numpy arrays.  The data can "
      "be accessed without copy using np.array(x, copy=False), "
      "with x as an instance of this class.  Note that access to "
      "any and all of the mathematical operations are only available "
      "via the numpy interface.  Also note that the equality operation "
      "only checks pointer equality and not element-wise equality\n"
      "\n"
      "Initialization:\n"
      "    Tensor6(): for basic initialization\n\n"
      "    Tensor6(Index, Index, Index, Index, Index, Index): for constant size, unknown value\n\n"
      "    Tensor6(Index, Index, Index, Index, Index, Index, Numeric): for constant size, constant value\n\n"
      "    Tensor6(List or Array): to copy elements\n\n";

  py::class_<Tensor7>(m, "Tensor7", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index, Index, Index, Index, Index, Index, Index>())
      .def(py::init<Index, Index, Index, Index, Index, Index, Index, Numeric>())
      .def(py::init([](const py::array_t<Numeric>& arr) {
        ARTS_USER_ERROR_IF(arr.request().ndim not_eq 7, "Bad size array")

        auto arr_val = arr.unchecked<7>();

        Tensor7 out(arr_val.shape(0),
                    arr_val.shape(1),
                    arr_val.shape(2),
                    arr_val.shape(3),
                    arr_val.shape(4),
                    arr_val.shape(5),
                    arr_val.shape(6));

        for (Index l = 0; l < out.nlibraries(); l++) {
          for (Index v = 0; v < out.nvitrines(); v++) {
            for (Index s = 0; s < out.nshelves(); s++) {
              for (Index b = 0; b < out.nbooks(); b++) {
                for (Index p = 0; p < out.npages(); p++) {
                  for (Index r = 0; r < out.nrows(); r++) {
                    for (Index c = 0; c < out.ncols(); c++) {
                      out(l, v, s, b, p, r, c) = arr_val(l, v, s, b, p, r, c);
                    }
                  }
                }
              }
            }
          }
        }

        return out;
      }))
      .PythonInterfaceBasicRepresentation(Tensor7)
      .PythonInterfaceFileIO(Tensor7)
      .def_property_readonly(
          "shape",
          [](const Tensor7& x) {
            return std::array{x.nlibraries(),
                              x.nvitrines(),
                              x.nshelves(),
                              x.nbooks(),
                              x.npages(),
                              x.nrows(),
                              x.ncols()};
          },
          "The shape of the data")
      .def(
          "__getitem__",
          [](const Tensor7& x,
             std::tuple<Index, Index, Index, Index, Index, Index, Index> inds) {
            auto [l, v, s, b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                x.nshelves() <= s or s < 0 or x.nvitrines() <= v or v < 0 or
                x.nlibraries() <= l or l < 0)
              throw std::out_of_range("Out of bounds");
            return x(l, v, s, b, p, r, c);
          },
          py::return_value_policy::reference_internal)
      .def("__setitem__",
           [](Tensor7& x,
              std::tuple<Index, Index, Index, Index, Index, Index, Index> inds,
              Numeric y) {
             auto [l, v, s, b, p, r, c] = inds;
             if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                 x.nshelves() <= s or s < 0 or x.nvitrines() <= v or v < 0 or
                 x.nlibraries() <= l or l < 0)
               throw std::out_of_range("Out of bounds");
             x(l, v, s, b, p, r, c) = y;
           })
      .def_buffer([](Tensor7& x) -> py::buffer_info {
        return py::buffer_info(
            x.get_c_array(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            7,
            {x.nlibraries(),
             x.nvitrines(),
             x.nshelves(),
             x.nbooks(),
             x.npages(),
             x.nrows(),
             x.ncols()},
            {sizeof(Numeric) * x.nvitrines() * x.nshelves() * x.nbooks() *
                 x.npages() * x.ncols() * x.nrows(),
             sizeof(Numeric) * x.nshelves() * x.nbooks() * x.npages() *
                 x.ncols() * x.nrows(),
             sizeof(Numeric) * x.nbooks() * x.npages() * x.ncols() * x.nrows(),
             sizeof(Numeric) * x.npages() * x.ncols() * x.nrows(),
             sizeof(Numeric) * x.ncols() * x.nrows(),
             sizeof(Numeric) * x.ncols(),
             sizeof(Numeric)});
      })
      .doc() =
      "The Arts Tensor7 class\n"
      "\n"
      "This class is compatible with numpy arrays.  The data can "
      "be accessed without copy using np.array(x, copy=False), "
      "with x as an instance of this class.  Note that access to "
      "any and all of the mathematical operations are only available "
      "via the numpy interface.  Also note that the equality operation "
      "only checks pointer equality and not element-wise equality\n"
      "\n"
      "Initialization:\n"
      "    Tensor7(): for basic initialization\n\n"
      "    Tensor7(Index, Index, Index, Index, Index, Index, Index): for constant size, unknown value\n\n"
      "    Tensor7(Index, Index, Index, Index, Index, Index, Index, Numeric): for constant size, constant value\n\n"
      "    Tensor7(List or Array): to copy elements\n\n";

  PythonInterfaceWorkspaceArray(Vector);
  PythonInterfaceWorkspaceArray(Matrix);
  PythonInterfaceWorkspaceArray(Tensor3);
  PythonInterfaceWorkspaceArray(Tensor4);
  PythonInterfaceWorkspaceArray(Tensor5);
  PythonInterfaceWorkspaceArray(Tensor6);
  PythonInterfaceWorkspaceArray(Tensor7);

  PythonInterfaceWorkspaceArray(ArrayOfVector);
  PythonInterfaceWorkspaceArray(ArrayOfMatrix);
  PythonInterfaceWorkspaceArray(ArrayOfTensor3);
  PythonInterfaceWorkspaceArray(ArrayOfTensor6);

  py::implicitly_convertible<py::array, Vector>();
  py::implicitly_convertible<py::array, Matrix>();
  py::implicitly_convertible<py::array, Tensor3>();
  py::implicitly_convertible<py::array, Tensor4>();
  py::implicitly_convertible<py::array, Tensor5>();
  py::implicitly_convertible<py::array, Tensor6>();
  py::implicitly_convertible<py::array, Tensor7>();

  py::implicitly_convertible<py::list, Vector>();
  py::implicitly_convertible<py::list, Matrix>();
  py::implicitly_convertible<py::list, Tensor3>();
  py::implicitly_convertible<py::list, Tensor4>();
  py::implicitly_convertible<py::list, Tensor5>();
  py::implicitly_convertible<py::list, Tensor6>();
  py::implicitly_convertible<py::list, Tensor7>();

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
}
}  // namespace Python