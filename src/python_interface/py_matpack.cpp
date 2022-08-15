#include <py_auto_interface.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#include <type_traits>
#include <vector>

#include "py_macros.h"

namespace Python {
using Scalar = std::variant<Numeric, Index>;

template <typename T>
void test_correct_size(const std::vector<T>& x) {
  if constexpr (not std::is_same_v<Scalar, T>) {
    // T is a vector!
    ARTS_USER_ERROR_IF(
        x.size() and std::any_of(x.begin() + 1,
                                 x.end(),
                                 [n = x.front().size()](auto& v) {
                                   return v.size() != n;
                                 }),
        "Bad size")
    for (auto& y : x) test_correct_size(y);
  }
}

void py_matpack(py::module_& m) {
  py::class_<Range>(m, "Range")
      .def(py::init([](Index a, Index b, Index c) {
             ARTS_USER_ERROR_IF(0 > a, "Bad start")
             ARTS_USER_ERROR_IF(0 > b, "Bad extent")
             return new Range{a, b, c};
           }),
           py::arg("start"),
           py::arg("extent"),
           py::arg("stride") = 1)
      .PythonInterfaceBasicRepresentation(Range)
      .def(py::pickle(
          [](const Range& self) {
            return py::make_tuple(
                self.get_start(), self.get_extent(), self.get_stride());
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 3, "Invalid state!")
            return new Range{
                t[0].cast<Index>(), t[1].cast<Index>(), t[2].cast<Index>()};
          }));

  py::class_<ConstVectorView>(m, "ConstVectorView");
  py::class_<ConstMatrixView>(m, "ConstMatrixView");
  py::class_<ConstTensor3View>(m, "ConstTensor3View");
  py::class_<ConstTensor4View>(m, "ConstTensor4View");
  py::class_<ConstTensor5View>(m, "ConstTensor5View");
  py::class_<ConstTensor6View>(m, "ConstTensor6View");
  py::class_<ConstTensor7View>(m, "ConstTensor7View");
  py::class_<VectorView, ConstVectorView>(m, "VectorView");
  py::class_<MatrixView, ConstMatrixView>(m, "MatrixView");
  py::class_<Tensor3View, ConstTensor3View>(m, "Tensor3View");
  py::class_<Tensor4View, ConstTensor4View>(m, "Tensor4View");
  py::class_<Tensor5View, ConstTensor5View>(m, "Tensor5View");
  py::class_<Tensor6View, ConstTensor6View>(m, "Tensor6View");
  py::class_<Tensor7View, ConstTensor7View>(m, "Tensor7View");

  py::class_<Vector, VectorView>(m, "Vector", py::buffer_protocol())
      .def(py::init([]() { return new Vector{}; }))
      .def(py::init([](const std::vector<Scalar>& v) {
             auto* out = new Vector(v.size());
             for (size_t i = 0; i < v.size(); i++)
               out->operator[](i) = std::visit(
                   [](auto&& x) { return static_cast<Numeric>(x); }, v[i]);
             return out;
           }),
           py::arg("vec").none(false))
      .PythonInterfaceCopyValue(Vector)
      .PythonInterfaceWorkspaceVariableConversion(Vector)
      .PythonInterfaceValueOperators.PythonInterfaceBasicRepresentation(Vector)
      .PythonInterfaceFileIO(Vector)
      .def_buffer([](Vector& x) -> py::buffer_info {
        return py::buffer_info(x.get_c_array(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               1,
                               {x.nelem()},
                               {sizeof(Numeric)});
      })
      .def_property("value",
                    py::cpp_function(
                        [](Vector& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Vector& x, Vector& y) { x = y; })
      .def(py::pickle(
          [](const py::object& self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<Vector>()(t[0]).cast<Vector>();
          }))
      .PythonInterfaceWorkspaceDocumentationExtra(Vector, R"--(

This class is mostly compatible with numpy arrays including numpy math.
The data can be accessed without copy using np.array(x, copy=False) or
via x.value)--");

  py::class_<Matrix, MatrixView>(m, "Matrix", py::buffer_protocol())
      .def(py::init([]() { return new Matrix{}; }))
      .def(py::init([](const std::vector<std::vector<Scalar>>& v) {
             test_correct_size(v);
             auto n1 = v.size();
             auto n2 = n1 > 0 ? v[0].size() : 0;
             auto* out = new Matrix(n1, n2);
             for (size_t i = 0; i < n1; i++) {
               for (size_t j = 0; j < n2; j++) {
                 out->operator()(i, j) = std::visit(
                     [](auto&& x) { return static_cast<Numeric>(x); }, v[i][j]);
               }
             }
             return out;
           }),
           py::arg("mat").none(false))
      .PythonInterfaceCopyValue(Matrix)
      .PythonInterfaceWorkspaceVariableConversion(Matrix)
      .PythonInterfaceBasicRepresentation(Matrix)
      .PythonInterfaceFileIO(Matrix)
      .PythonInterfaceValueOperators
      .def_property_readonly(
          "T",
          [](const Matrix& x) { return Matrix(transpose(x)); },
          "Non-trivial transpose")
      .def_buffer([](Matrix& x) -> py::buffer_info {
        return py::buffer_info(
            x.get_c_array(),
            sizeof(Numeric),
            py::format_descriptor<Numeric>::format(),
            2,
            {static_cast<ssize_t>(x.nrows()), static_cast<ssize_t>(x.ncols())},
            {sizeof(Numeric) * static_cast<ssize_t>(x.ncols()),
             sizeof(Numeric)});
      })
      .def_property("value",
                    py::cpp_function(
                        [](Matrix& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Matrix& x, Matrix& y) { x = y; })
      .def(py::pickle(
          [](const py::object& self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<Matrix>()(t[0]).cast<Matrix>();
          }))
      .PythonInterfaceWorkspaceDocumentationExtra(Matrix, R"--(

This class is mostly compatible with numpy arrays including numpy math.
The data can be accessed without copy using np.array(x, copy=False) or
via x.value)--");

  py::class_<Tensor3, Tensor3View>(m, "Tensor3", py::buffer_protocol())
      .def(py::init([]() { return new Tensor3{}; }))
      .def(py::init([](const std::vector<std::vector<std::vector<Scalar>>>& v) {
             test_correct_size(v);
             auto n1 = v.size();
             auto n2 = n1 > 0 ? v[0].size() : 0;
             auto n3 = n2 > 0 ? v[0][0].size() : 0;
             auto* out = new Tensor3(n1, n2, n3);
             for (size_t i1 = 0; i1 < n1; i1++) {
               for (size_t i2 = 0; i2 < n2; i2++) {
                 for (size_t i3 = 0; i3 < n3; i3++) {
                   out->operator()(i1, i2, i3) = std::visit(
                       [](auto&& x) { return static_cast<Numeric>(x); },
                       v[i1][i2][i3]);
                 }
               }
             }
             return out;
           }),
           py::arg("ten3").none(false))
      .PythonInterfaceCopyValue(Tensor3)
      .PythonInterfaceWorkspaceVariableConversion(Tensor3)
      .PythonInterfaceBasicRepresentation(Tensor3)
      .PythonInterfaceFileIO(Tensor3)
      .PythonInterfaceValueOperators
      .def_buffer([](Tensor3& x) -> py::buffer_info {
        return py::buffer_info(
            x.get_c_array(),
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor3& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor3& x, Tensor3& y) { x = y; })
      .def(py::pickle(
          [](const py::object& self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<Tensor3>()(t[0]).cast<Tensor3>();
          }))
      .PythonInterfaceWorkspaceDocumentationExtra(Tensor3, R"--(

This class is mostly compatible with numpy arrays including numpy math.
The data can be accessed without copy using np.array(x, copy=False) or
via x.value)--");

  py::class_<Tensor4, Tensor4View>(m, "Tensor4", py::buffer_protocol())
      .def(py::init([]() { return new Tensor4{}; }))
      .def(py::init([](const std::vector<
                        std::vector<std::vector<std::vector<Scalar>>>>& v) {
             test_correct_size(v);
             auto n1 = v.size();
             auto n2 = n1 > 0 ? v[0].size() : 0;
             auto n3 = n2 > 0 ? v[0][0].size() : 0;
             auto n4 = n3 > 0 ? v[0][0][0].size() : 0;
             auto* out = new Tensor4(n1, n2, n3, n4);
             for (size_t i1 = 0; i1 < n1; i1++) {
               for (size_t i2 = 0; i2 < n2; i2++) {
                 for (size_t i3 = 0; i3 < n3; i3++) {
                   for (size_t i4 = 0; i4 < n4; i4++) {
                     out->operator()(i1, i2, i3, i4) = std::visit(
                         [](auto&& x) { return static_cast<Numeric>(x); },
                         v[i1][i2][i3][i4]);
                   }
                 }
               }
             }
             return out;
           }),
           py::arg("ten4").none(false))
      .PythonInterfaceCopyValue(Tensor4)
      .PythonInterfaceWorkspaceVariableConversion(Tensor4)
      .PythonInterfaceBasicRepresentation(Tensor4)
      .PythonInterfaceFileIO(Tensor4)
      .PythonInterfaceValueOperators
      .def_buffer([](Tensor4& x) -> py::buffer_info {
        return py::buffer_info(
            x.get_c_array(),
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor4& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor4& x, Tensor4& y) { x = y; })
      .def(py::pickle(
          [](const py::object& self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<Tensor4>()(t[0]).cast<Tensor4>();
          }))
      .PythonInterfaceWorkspaceDocumentationExtra(Tensor4, R"--(

This class is mostly compatible with numpy arrays including numpy math.
The data can be accessed without copy using np.array(x, copy=False) or
via x.value)--");

  py::class_<Tensor5, Tensor5View>(m, "Tensor5", py::buffer_protocol())
      .def(py::init([]() { return new Tensor5{}; }))
      .def(py::init([](const std::vector<std::vector<
                           std::vector<std::vector<std::vector<Scalar>>>>>& v) {
             test_correct_size(v);
             auto n1 = v.size();
             auto n2 = n1 > 0 ? v[0].size() : 0;
             auto n3 = n2 > 0 ? v[0][0].size() : 0;
             auto n4 = n3 > 0 ? v[0][0][0].size() : 0;
             auto n5 = n4 > 0 ? v[0][0][0][0].size() : 0;
             auto* out = new Tensor5(n1, n2, n3, n4, n5);
             for (size_t i1 = 0; i1 < n1; i1++) {
               for (size_t i2 = 0; i2 < n2; i2++) {
                 for (size_t i3 = 0; i3 < n3; i3++) {
                   for (size_t i4 = 0; i4 < n4; i4++) {
                     for (size_t i5 = 0; i5 < n5; i5++) {
                       out->operator()(i1, i2, i3, i4, i5) = std::visit(
                           [](auto&& x) { return static_cast<Numeric>(x); },
                           v[i1][i2][i3][i4][i5]);
                     }
                   }
                 }
               }
             }
             return out;
           }),
           py::arg("ten5").none(false))
      .PythonInterfaceCopyValue(Tensor5)
      .PythonInterfaceWorkspaceVariableConversion(Tensor5)
      .PythonInterfaceBasicRepresentation(Tensor5)
      .PythonInterfaceFileIO(Tensor5)
      .PythonInterfaceValueOperators
      .def_buffer([](Tensor5& x) -> py::buffer_info {
        return py::buffer_info(
            x.get_c_array(),
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor5& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor5& x, Tensor5& y) { x = y; })
      .def(py::pickle(
          [](const py::object& self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<Tensor5>()(t[0]).cast<Tensor5>();
          }))
      .PythonInterfaceWorkspaceDocumentationExtra(Tensor5, R"--(

This class is mostly compatible with numpy arrays including numpy math.
The data can be accessed without copy using np.array(x, copy=False) or
via x.value)--");

  py::class_<Tensor6, Tensor6View>(m, "Tensor6", py::buffer_protocol())
      .def(py::init([]() { return new Tensor6{}; }))
      .def(
          py::init([](const std::vector<std::vector<std::vector<
                          std::vector<std::vector<std::vector<Scalar>>>>>>& v) {
            test_correct_size(v);
            auto n1 = v.size();
            auto n2 = n1 > 0 ? v[0].size() : 0;
            auto n3 = n2 > 0 ? v[0][0].size() : 0;
            auto n4 = n3 > 0 ? v[0][0][0].size() : 0;
            auto n5 = n4 > 0 ? v[0][0][0][0].size() : 0;
            auto n6 = n5 > 0 ? v[0][0][0][0][0].size() : 0;
            auto* out = new Tensor6(n1, n2, n3, n4, n5, n6);
            for (size_t i1 = 0; i1 < n1; i1++) {
              for (size_t i2 = 0; i2 < n2; i2++) {
                for (size_t i3 = 0; i3 < n3; i3++) {
                  for (size_t i4 = 0; i4 < n4; i4++) {
                    for (size_t i5 = 0; i5 < n5; i5++) {
                      for (size_t i6 = 0; i6 < n6; i6++) {
                        out->operator()(i1, i2, i3, i4, i5, i6) = std::visit(
                            [](auto&& x) { return static_cast<Numeric>(x); },
                            v[i1][i2][i3][i4][i5][i6]);
                      }
                    }
                  }
                }
              }
            }
            return out;
          }),
          py::arg("ten6").none(false))
      .PythonInterfaceCopyValue(Tensor6)
      .PythonInterfaceWorkspaceVariableConversion(Tensor6)
      .PythonInterfaceBasicRepresentation(Tensor6)
      .PythonInterfaceFileIO(Tensor6)
      .PythonInterfaceValueOperators
      .def_buffer([](Tensor6& x) -> py::buffer_info {
        return py::buffer_info(
            x.get_c_array(),
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor6& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor6& x, Tensor6& y) { x = y; })
      .def(py::pickle(
          [](const py::object& self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<Tensor6>()(t[0]).cast<Tensor6>();
          }))
      .PythonInterfaceWorkspaceDocumentationExtra(Tensor6, R"--(

This class is mostly compatible with numpy arrays including numpy math.
The data can be accessed without copy using np.array(x, copy=False) or
via x.value)--");

  py::class_<Tensor7, Tensor7View>(m, "Tensor7", py::buffer_protocol())
      .def(py::init([]() { return new Tensor7{}; }))
      .def(py::init(
               [](const std::vector<std::vector<std::vector<std::vector<
                      std::vector<std::vector<std::vector<Scalar>>>>>>>& v) {
                 test_correct_size(v);
                 auto n1 = v.size();
                 auto n2 = n1 > 0 ? v[0].size() : 0;
                 auto n3 = n2 > 0 ? v[0][0].size() : 0;
                 auto n4 = n3 > 0 ? v[0][0][0].size() : 0;
                 auto n5 = n4 > 0 ? v[0][0][0][0].size() : 0;
                 auto n6 = n5 > 0 ? v[0][0][0][0][0].size() : 0;
                 auto n7 = n6 > 0 ? v[0][0][0][0][0][0].size() : 0;
                 auto* out = new Tensor7(n1, n2, n3, n4, n5, n6, n7);
                 for (size_t i1 = 0; i1 < n1; i1++) {
                   for (size_t i2 = 0; i2 < n2; i2++) {
                     for (size_t i3 = 0; i3 < n3; i3++) {
                       for (size_t i4 = 0; i4 < n4; i4++) {
                         for (size_t i5 = 0; i5 < n5; i5++) {
                           for (size_t i6 = 0; i6 < n6; i6++) {
                             for (size_t i7 = 0; i7 < n7; i7++) {
                               out->operator()(i1, i2, i3, i4, i5, i6, i7) =
                                   std::visit(
                                       [](auto&& x) {
                                         return static_cast<Numeric>(x);
                                       },
                                       v[i1][i2][i3][i4][i5][i6][i7]);
                             }
                           }
                         }
                       }
                     }
                   }
                 }
                 return out;
               }),
           py::arg("ten7").none(false))
      .PythonInterfaceCopyValue(Tensor7)
      .PythonInterfaceWorkspaceVariableConversion(Tensor7)
      .PythonInterfaceFileIO(Tensor7)
      .PythonInterfaceValueOperators.PythonInterfaceBasicRepresentation(Tensor7)
      .def_buffer([](Tensor7& x) -> py::buffer_info {
        return py::buffer_info(
            x.get_c_array(),
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor7& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor7& x, Tensor7& y) { x = y; })
      .def(py::pickle(
          [](const py::object& self) {
            return py::make_tuple(self.attr("value"));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return py::type::of<Tensor7>()(t[0]).cast<Tensor7>();
          }))
      .PythonInterfaceWorkspaceDocumentationExtra(Tensor7, R"--(

This class is mostly compatible with numpy arrays including numpy math.
The data can be accessed without copy using np.array(x, copy=False) or
via x.value)--");

  py::implicitly_convertible<std::vector<Scalar>, Vector>();
  py::implicitly_convertible<std::vector<std::vector<Scalar>>, Matrix>();
  py::implicitly_convertible<std::vector<std::vector<std::vector<Scalar>>>, Tensor3>();
  py::implicitly_convertible<std::vector<std::vector<std::vector<std::vector<Scalar>>>>, Tensor4>();
  py::implicitly_convertible<std::vector<std::vector<std::vector<std::vector<std::vector<Scalar>>>>>, Tensor5>();
  py::implicitly_convertible<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Scalar>>>>>>, Tensor6>();
  py::implicitly_convertible<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Scalar>>>>>>>, Tensor7>();

  PythonInterfaceWorkspaceArray(Vector);
  PythonInterfaceWorkspaceArray(Matrix);
  PythonInterfaceWorkspaceArray(Tensor3);
  PythonInterfaceWorkspaceArray(Tensor4);
  PythonInterfaceWorkspaceArray(Tensor5);
  PythonInterfaceWorkspaceArray(Tensor6);
  PythonInterfaceWorkspaceArray(Tensor7);

  PythonInterfaceWorkspaceArray(ArrayOfVector)
      .def(py::init([](const std::vector<std::vector<Vector>>& x) {
             ArrayOfArrayOfVector y(x.size());
             std::copy(x.begin(), x.end(), y.begin());
             return y;
           }),
           py::arg("arr").none(false));
  py::implicitly_convertible<std::vector<std::vector<Vector>>,
                             ArrayOfArrayOfVector>();

  PythonInterfaceWorkspaceArray(ArrayOfMatrix)
      .def(py::init([](const std::vector<std::vector<Matrix>>& x) {
             ArrayOfArrayOfMatrix y(x.size());
             std::copy(x.begin(), x.end(), y.begin());
             return y;
           }),
           py::arg("arr").none(false));
  py::implicitly_convertible<std::vector<std::vector<Matrix>>,
                             ArrayOfArrayOfMatrix>();

  PythonInterfaceWorkspaceArray(ArrayOfTensor3)
      .def(py::init([](const std::vector<std::vector<Tensor3>>& x) {
             ArrayOfArrayOfTensor3 y(x.size());
             std::copy(x.begin(), x.end(), y.begin());
             return y;
           }),
           py::arg("arr").none(false));
  py::implicitly_convertible<std::vector<std::vector<Tensor3>>,
                             ArrayOfArrayOfTensor3>();

  PythonInterfaceWorkspaceArray(ArrayOfTensor6)
      .def(py::init([](const std::vector<std::vector<Tensor6>>& x) {
             ArrayOfArrayOfTensor6 y(x.size());
             std::copy(x.begin(), x.end(), y.begin());
             return y;
           }),
           py::arg("arr").none(false));
  py::implicitly_convertible<std::vector<std::vector<Tensor6>>,
                             ArrayOfArrayOfTensor6>();

  py::class_<Rational>(m, "Rational")
      .def(py::init([](Index n, Index d) {
             ARTS_USER_ERROR_IF(d < 1, "Must be positive")
             return Rational(n, d);
           }),
           py::arg("n") = 0,
           py::arg("d") = 1)
      .PythonInterfaceCopyValue(Rational)
      .PythonInterfaceWorkspaceVariableConversion(Rational)
      .def(py::init([](const String& s) { return new Rational{s}; }))
      .def(py::init([](Numeric n) { return Rational(std::to_string(n)); }))
      .PythonInterfaceFileIO(Rational)
      .PythonInterfaceBasicRepresentation(Rational)
      .PythonInterfaceInPlaceMathOperators(Rational, Rational)
      .PythonInterfaceInPlaceMathOperators(Rational, Index)
      .PythonInterfaceMathOperators(Rational, Rational)
      .PythonInterfaceMathOperators(Rational, Index)
      .PythonInterfaceComparisonOperators(Rational, Rational)
      .PythonInterfaceComparisonOperators(Rational, Index)
      .def(py::pickle(
          [](const Rational& self) {
            return py::make_tuple(self.numer, self.denom);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 2, "Invalid state!")

            return new Rational{t[0].cast<Index>(), t[1].cast<Index>()};
          }))
      .PythonInterfaceWorkspaceDocumentation(Rational);
  py::implicitly_convertible<Index, Rational>();
}
}  // namespace Python
