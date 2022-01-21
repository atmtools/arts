#include <matpackVII.h>
#include <xml_io.h>

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
      });
  py::implicitly_convertible<py::array, Vector>();
  py::implicitly_convertible<py::list, Vector>();

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
      .def("ncols", [](const Matrix& x) { return x.ncols(); })
      .def("nrows", [](const Matrix& x) { return x.nrows(); })
      .def("__getitem__",
           [](const Matrix& x, std::tuple<Index, Index> inds) {
             auto [r, c] = inds;
             if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0)
               throw std::out_of_range(var_string("Bad index access ",
                                                  '[',
                                                  r,
                                                  ", ",
                                                  c,
                                                  ']',
                                                  " in object of shape ",
                                                  '(',
                                                  x.nrows(),
                                                  ", ",
                                                  x.ncols(),
                                                  ')'));
             return x(r, c);
           })
      .def("__setitem__",
           [](Matrix& x, std::tuple<Index, Index> inds, Numeric y) {
             auto [r, c] = inds;
             if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0)
               throw std::out_of_range(var_string("Bad index access ",
                                                  '[',
                                                  r,
                                                  ", ",
                                                  c,
                                                  ']',
                                                  " in object of shape ",
                                                  '(',
                                                  x.nrows(),
                                                  ", ",
                                                  x.ncols(),
                                                  ')'));
             x(r, c) = y;
           })
      .def_buffer([](Matrix& x) -> py::buffer_info {
        return py::buffer_info(x.get_c_array(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               2,
                               {x.nrows(), x.ncols()},
                               {sizeof(Numeric) * x.ncols(), sizeof(Numeric)});
      });

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
      .def("ncols", [](const Tensor3& x) { return x.ncols(); })
      .def("nrows", [](const Tensor3& x) { return x.nrows(); })
      .def("npages", [](const Tensor3& x) { return x.npages(); })
      .def(
          "__getitem__",
          [](const Tensor3& x, std::tuple<Index, Index, Index> inds) {
            auto [p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0)
              throw std::out_of_range(var_string("Bad index access ",
                                                 '[',
                                                 p,
                                                 ", ",
                                                 r,
                                                 ", ",
                                                 c,
                                                 ']',
                                                 " in object of shape ",
                                                 '(',
                                                 x.npages(),
                                                 ", ",
                                                 x.nrows(),
                                                 ", ",
                                                 x.ncols(),
                                                 ')'));
            return x(p, r, c);
          },
          py::return_value_policy::reference_internal)
      .def("__setitem__",
           [](Tensor3& x, std::tuple<Index, Index, Index> inds, Numeric y) {
             auto [p, r, c] = inds;
             if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0)
               throw std::out_of_range(var_string("Bad index access ",
                                                  '[',
                                                  p,
                                                  ", ",
                                                  r,
                                                  ", ",
                                                  c,
                                                  ']',
                                                  " in object of shape ",
                                                  '(',
                                                  x.npages(),
                                                  ", ",
                                                  x.nrows(),
                                                  ", ",
                                                  x.ncols(),
                                                  ')'));
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
      });

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
      .def("ncols", [](const Tensor4& x) { return x.ncols(); })
      .def("nrows", [](const Tensor4& x) { return x.nrows(); })
      .def("npages", [](const Tensor4& x) { return x.npages(); })
      .def("nbooks", [](const Tensor4& x) { return x.nbooks(); })
      .def(
          "__getitem__",
          [](const Tensor4& x, std::tuple<Index, Index, Index, Index> inds) {
            auto [b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0)
              throw std::out_of_range(var_string("Bad index access ",
                                                 '[',
                                                 b,
                                                 ", ",
                                                 p,
                                                 ", ",
                                                 r,
                                                 ", ",
                                                 c,
                                                 ']',
                                                 " in object of shape ",
                                                 '(',
                                                 x.nbooks(),
                                                 ", ",
                                                 x.npages(),
                                                 ", ",
                                                 x.nrows(),
                                                 ", ",
                                                 x.ncols(),
                                                 ')'));
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
               throw std::out_of_range(var_string("Bad index access ",
                                                  '[',
                                                  b,
                                                  ", ",
                                                  p,
                                                  ", ",
                                                  r,
                                                  ", ",
                                                  c,
                                                  ']',
                                                  " in object of shape ",
                                                  '(',
                                                  x.nbooks(),
                                                  ", ",
                                                  x.npages(),
                                                  ", ",
                                                  x.nrows(),
                                                  ", ",
                                                  x.ncols(),
                                                  ')'));
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
      });

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
      .def("ncols", [](const Tensor5& x) { return x.ncols(); })
      .def("nrows", [](const Tensor5& x) { return x.nrows(); })
      .def("npages", [](const Tensor5& x) { return x.npages(); })
      .def("nbooks", [](const Tensor5& x) { return x.nbooks(); })
      .def("nshelves", [](const Tensor5& x) { return x.nshelves(); })
      .def(
          "__getitem__",
          [](const Tensor5& x,
             std::tuple<Index, Index, Index, Index, Index> inds) {
            auto [s, b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                x.nshelves() <= s or s < 0)
              throw std::out_of_range(var_string("Bad index access ",
                                                 '[',
                                                 s,
                                                 ", ",
                                                 b,
                                                 ", ",
                                                 p,
                                                 ", ",
                                                 r,
                                                 ", ",
                                                 c,
                                                 ']',
                                                 " in object of shape ",
                                                 '(',
                                                 x.nshelves(),
                                                 ", ",
                                                 x.nbooks(),
                                                 ", ",
                                                 x.npages(),
                                                 ", ",
                                                 x.nrows(),
                                                 ", ",
                                                 x.ncols(),
                                                 ')'));
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
               throw std::out_of_range(var_string("Bad index access ",
                                                  '[',
                                                  s,
                                                  ", ",
                                                  b,
                                                  ", ",
                                                  p,
                                                  ", ",
                                                  r,
                                                  ", ",
                                                  c,
                                                  ']',
                                                  " in object of shape ",
                                                  '(',
                                                  x.nshelves(),
                                                  ", ",
                                                  x.nbooks(),
                                                  ", ",
                                                  x.npages(),
                                                  ", ",
                                                  x.nrows(),
                                                  ", ",
                                                  x.ncols(),
                                                  ')'));
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
      });

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
      .def("ncols", [](const Tensor6& x) { return x.ncols(); })
      .def("nrows", [](const Tensor6& x) { return x.nrows(); })
      .def("npages", [](const Tensor6& x) { return x.npages(); })
      .def("nbooks", [](const Tensor6& x) { return x.nbooks(); })
      .def("nshelves", [](const Tensor6& x) { return x.nshelves(); })
      .def("nvitrines", [](const Tensor6& x) { return x.nvitrines(); })
      .def(
          "__getitem__",
          [](const Tensor6& x,
             std::tuple<Index, Index, Index, Index, Index, Index> inds) {
            auto [v, s, b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                x.nshelves() <= s or s < 0 or x.nvitrines() <= v or v < 0)
              throw std::out_of_range(var_string("Bad index access ",
                                                 '[',
                                                 v,
                                                 ", ",
                                                 s,
                                                 ", ",
                                                 b,
                                                 ", ",
                                                 p,
                                                 ", ",
                                                 r,
                                                 ", ",
                                                 c,
                                                 ']',
                                                 " in object of shape ",
                                                 '(',
                                                 x.nvitrines(),
                                                 ", ",
                                                 x.nshelves(),
                                                 ", ",
                                                 x.nbooks(),
                                                 ", ",
                                                 x.npages(),
                                                 ", ",
                                                 x.nrows(),
                                                 ", ",
                                                 x.ncols(),
                                                 ')'));
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
               throw std::out_of_range(var_string("Bad index access ",
                                                  '[',
                                                  v,
                                                  ", ",
                                                  s,
                                                  ", ",
                                                  b,
                                                  ", ",
                                                  p,
                                                  ", ",
                                                  r,
                                                  ", ",
                                                  c,
                                                  ']',
                                                  " in object of shape ",
                                                  '(',
                                                  x.nvitrines(),
                                                  ", ",
                                                  x.nshelves(),
                                                  ", ",
                                                  x.nbooks(),
                                                  ", ",
                                                  x.npages(),
                                                  ", ",
                                                  x.nrows(),
                                                  ", ",
                                                  x.ncols(),
                                                  ')'));
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
      });

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
      .def("ncols", [](const Tensor7& x) { return x.ncols(); })
      .def("nrows", [](const Tensor7& x) { return x.nrows(); })
      .def("npages", [](const Tensor7& x) { return x.npages(); })
      .def("nbooks", [](const Tensor7& x) { return x.nbooks(); })
      .def("nshelves", [](const Tensor7& x) { return x.nshelves(); })
      .def("nvitrines", [](const Tensor7& x) { return x.nvitrines(); })
      .def("nlibraries", [](const Tensor7& x) { return x.nlibraries(); })
      .def(
          "__getitem__",
          [](const Tensor7& x,
             std::tuple<Index, Index, Index, Index, Index, Index, Index> inds) {
            auto [l, v, s, b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                x.nshelves() <= s or s < 0 or x.nvitrines() <= v or v < 0 or
                x.nlibraries() <= l or l < 0)
              throw std::out_of_range(var_string("Bad index access ",
                                                 '[',
                                                 l,
                                                 ", ",
                                                 v,
                                                 ", ",
                                                 s,
                                                 ", ",
                                                 b,
                                                 ", ",
                                                 p,
                                                 ", ",
                                                 r,
                                                 ", ",
                                                 c,
                                                 ']',
                                                 " in object of shape ",
                                                 '(',
                                                 x.nlibraries(),
                                                 ", ",
                                                 x.nvitrines(),
                                                 ", ",
                                                 x.nshelves(),
                                                 ", ",
                                                 x.nbooks(),
                                                 ", ",
                                                 x.npages(),
                                                 ", ",
                                                 x.nrows(),
                                                 ", ",
                                                 x.ncols(),
                                                 ')'));
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
               throw std::out_of_range(var_string("Bad index access ",
                                                  '[',
                                                  l,
                                                  ", ",
                                                  v,
                                                  ", ",
                                                  s,
                                                  ", ",
                                                  b,
                                                  ", ",
                                                  p,
                                                  ", ",
                                                  r,
                                                  ", ",
                                                  c,
                                                  ']',
                                                  " in object of shape ",
                                                  '(',
                                                  x.nlibraries(),
                                                  ", ",
                                                  x.nvitrines(),
                                                  ", ",
                                                  x.nshelves(),
                                                  ", ",
                                                  x.nbooks(),
                                                  ", ",
                                                  x.npages(),
                                                  ", ",
                                                  x.nrows(),
                                                  ", ",
                                                  x.ncols(),
                                                  ')'));
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
      });
}
}  // namespace Python