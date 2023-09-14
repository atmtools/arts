#include <pybind11/cast.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include <python_interface.h>

#include <array>
#include <cstddef>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "array.h"
#include "matpack_constexpr.h"
#include "matpack_iter.h"
#include "py_macros.h"

namespace Python {
using Scalar = std::variant<Numeric, Index>;

template <typename T, Index... sz>
auto register_matpack_constant_data(py::module_& m, const char* const name)
  requires(sizeof...(sz) > 0)
{
  using matpackT = matpack::matpack_constant_data<T, sz...>;
  using arrT = py::array_t<T, py::array::forcecast>;

  static constexpr std::size_t dim = sizeof...(sz);
  static constexpr std::size_t size = (sz * ...);
  static constexpr std::array<Index, dim> shape = matpackT::shape();

  artsclass<matpackT> r(m, name, py::buffer_protocol());
  r.def(py::init([]() { return std::make_shared<matpackT>(); }),
        "Default constant data")
      .def(py::init([](const arrT& x) -> std::shared_ptr<matpackT> {
             if (static_cast<std::size_t>(x.ndim()) != dim) {
               throw std::runtime_error(
                   var_string("Bad rank: ", dim, " vs ", x.ndim()));
             }

             std::array<Index, dim> arr_shape;
             for (std::size_t i = 0; i < dim; i++) {
               arr_shape[i] = static_cast<Index>(x.shape(i));
             }

             if (arr_shape != shape) {
               throw std::runtime_error(
                   var_string("Bad shape: ",
                              matpack::shape_help<dim>(shape),
                              " vs ",
                              matpack::shape_help<dim>(arr_shape)));
             }

             return std::make_shared<matpackT>(
                 x.template cast<std::array<T, size>>());
           }),
           "Cast from numeric :class:`~numpy.array`")
      .def(py::init([](const py::list& x) {
             return py::cast<matpackT>(x.cast<arrT>());
           }),
           "Cast from :class:`list`")
      .def(py::init([](const py::array& x) {
             return py::cast<matpackT>(x.cast<arrT>());
           }),
           "Cast from any :class:`~numpy.array`")
      .PythonInterfaceCopyValue(matpackT)
      .PythonInterfaceBasicRepresentation(matpackT)
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](matpackT& x) -> py::buffer_info {
        std::vector<ssize_t> shpe{static_cast<ssize_t>(sz)...};
        std::vector<ssize_t> strides(shpe.size(),
                                     static_cast<ssize_t>(sizeof(T)));
        for (std::size_t j = 0; j < sizeof...(sz) - 1; j++) {
          strides[j] *= std::reduce(
              shpe.begin() + j + 1, shpe.end(), 1, std::multiplies<>());
        }
        return py::buffer_info(x.data.data(),
                               sizeof(T),
                               py::format_descriptor<T>::format(),
                               sizeof...(sz),
                               std::move(shpe),
                               std::move(strides));
      })
      .def_property(
          "value",
          py::cpp_function(
              [](matpackT& x) {
                py::object np = py::module_::import("numpy");
                return np.attr("array")(x, py::arg("copy") = false);
              },
              py::keep_alive<0, 1>()),
          [](matpackT& x, matpackT& y) { x = y; },
          py::doc(":class:`~numpy.ndarray` Data array"));
  py::implicitly_convertible<arrT, matpackT>();
  py::implicitly_convertible<py::list, matpackT>();
  py::implicitly_convertible<py::array, matpackT>();

  return r;
}

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

void py_matpack(py::module_& m) try {
  register_matpack_constant_data<Numeric, 3>(m, "Vector3").doc() = "A 3D vector";

  artsclass<Range>(m, "Range")
      .def(py::init([](Index a, Index b, Index c) {
             ARTS_USER_ERROR_IF(0 > a, "Bad offset")
             ARTS_USER_ERROR_IF(0 > b, "Bad extent")
             return std::make_shared<Range>(a, b, c);
           }),
           py::arg("offset"),
           py::arg("extent"),
           py::arg("stride") = 1, "Valued initialization")
      .PythonInterfaceBasicRepresentation(Range)
      .def(py::pickle(
          [](const Range& self) {
            return py::make_tuple(
                self.offset, self.extent, self.stride);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 3, "Invalid state!")
            return std::make_shared<Range>(
                t[0].cast<Index>(), t[1].cast<Index>(), t[2].cast<Index>());
          })).doc() = "A range, used to select parts of a matpack type";


  artsclass<Vector>(m, "Vector", py::buffer_protocol())
      .def(py::init([]() { return std::make_shared<Vector>(); }), "Default vector")
      .def(py::init([](const std::vector<Scalar>& v) {
             auto out = std::make_shared<Vector>(v.size());;
             for (size_t i = 0; i < v.size(); i++)
               out->operator[](i) = std::visit(
                   [](auto&& x) { return static_cast<Numeric>(x); }, v[i]);
             return out;
           }),
           py::arg("vec").none(false), py::doc("From :class:`list` equivalent"))
      .PythonInterfaceCopyValue(Vector)
      .PythonInterfaceWorkspaceVariableConversion(Vector)
      .PythonInterfaceBasicRepresentation(Vector)
      .PythonInterfaceFileIO(Vector)
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](Vector& x) -> py::buffer_info {
        return py::buffer_info(x.data_handle(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               1,
                               {x.size()},
                               {sizeof(Numeric)});
      })
      .def_property("value",
                    py::cpp_function(
                        [](Vector& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Vector& x, Vector& y) { x = y; }, py::doc(":class:`~numpy.ndarray` Data array"))
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
The data can be accessed without copy using ``np.array(x, copy=False)`` or
via x.value
)--");

  artsclass<Matrix>(m, "Matrix", py::buffer_protocol())
      .def(py::init([]() { return std::make_shared<Matrix>(); }), "Default matrix")
      .def(py::init([](const std::vector<std::vector<Scalar>>& v) {
             test_correct_size(v);
             auto n1 = v.size();
             auto n2 = n1 > 0 ? v[0].size() : 0;
             auto out = std::make_shared<Matrix>(n1, n2);
             for (size_t i = 0; i < n1; i++) {
               for (size_t j = 0; j < n2; j++) {
                 out->operator()(i, j) = std::visit(
                     [](auto&& x) { return static_cast<Numeric>(x); }, v[i][j]);
               }
             }
             return out;
           }),
           py::arg("mat").none(false), py::doc("From :class:`list` equivalent"))
      .PythonInterfaceCopyValue(Matrix)
      .PythonInterfaceWorkspaceVariableConversion(Matrix)
      .PythonInterfaceBasicRepresentation(Matrix)
      .PythonInterfaceFileIO(Matrix)
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
      .def_property("value",
                    py::cpp_function(
                        [](Matrix& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Matrix& x, Matrix& y) { x = y; }, py::doc(":class:`~numpy.ndarray` Data array"))
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
The data can be accessed without copy using ``np.array(x, copy=False)`` or
via x.value
)--");

  artsclass<Tensor3>(m, "Tensor3", py::buffer_protocol())
      .def(py::init([]() { return std::make_shared<Tensor3>(); }), "Default tensor")
      .def(py::init([](const std::vector<std::vector<std::vector<Scalar>>>& v) {
             test_correct_size(v);
             auto n1 = v.size();
             auto n2 = n1 > 0 ? v[0].size() : 0;
             auto n3 = n2 > 0 ? v[0][0].size() : 0;
             auto out = std::make_shared<Tensor3>(n1, n2, n3);
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
           py::arg("ten3").none(false), py::doc("From :class:`list` equivalent"))
      .PythonInterfaceCopyValue(Tensor3)
      .PythonInterfaceWorkspaceVariableConversion(Tensor3)
      .PythonInterfaceBasicRepresentation(Tensor3)
      .PythonInterfaceFileIO(Tensor3)
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor3& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor3& x, Tensor3& y) { x = y; }, py::doc(":class:`~numpy.ndarray` Data array"))
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
The data can be accessed without copy using ``np.array(x, copy=False)`` or
via x.value
)--");

  artsclass<Tensor4>(m, "Tensor4", py::buffer_protocol())
      .def(py::init([]() { return std::make_shared<Tensor4>(); }), "Default tensor")
      .def(py::init([](const std::vector<
                        std::vector<std::vector<std::vector<Scalar>>>>& v) {
             test_correct_size(v);
             auto n1 = v.size();
             auto n2 = n1 > 0 ? v[0].size() : 0;
             auto n3 = n2 > 0 ? v[0][0].size() : 0;
             auto n4 = n3 > 0 ? v[0][0][0].size() : 0;
             auto out = std::make_shared<Tensor4>(n1, n2, n3, n4);
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
           py::arg("ten4").none(false), py::doc("From :class:`list` equivalent"))
      .PythonInterfaceCopyValue(Tensor4)
      .PythonInterfaceWorkspaceVariableConversion(Tensor4)
      .PythonInterfaceBasicRepresentation(Tensor4)
      .PythonInterfaceFileIO(Tensor4)
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor4& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor4& x, Tensor4& y) { x = y; }, py::doc(":class:`~numpy.ndarray` Data array"))
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
The data can be accessed without copy using ``np.array(x, copy=False)`` or
via x.value
)--");

  artsclass<Tensor5>(m, "Tensor5", py::buffer_protocol())
      .def(py::init([]() { return std::make_shared<Tensor5>(); }), "Default tensor")
      .def(py::init([](const std::vector<std::vector<
                           std::vector<std::vector<std::vector<Scalar>>>>>& v) {
             test_correct_size(v);
             auto n1 = v.size();
             auto n2 = n1 > 0 ? v[0].size() : 0;
             auto n3 = n2 > 0 ? v[0][0].size() : 0;
             auto n4 = n3 > 0 ? v[0][0][0].size() : 0;
             auto n5 = n4 > 0 ? v[0][0][0][0].size() : 0;
             auto out = std::make_shared<Tensor5>(n1, n2, n3, n4, n5);
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
           py::arg("ten5").none(false), py::doc("From :class:`list` equivalent"))
      .PythonInterfaceCopyValue(Tensor5)
      .PythonInterfaceWorkspaceVariableConversion(Tensor5)
      .PythonInterfaceBasicRepresentation(Tensor5)
      .PythonInterfaceFileIO(Tensor5)
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor5& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor5& x, Tensor5& y) { x = y; }, py::doc(":class:`~numpy.ndarray` Data array"))
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
The data can be accessed without copy using ``np.array(x, copy=False)`` or
via x.value
)--");

  artsclass<Tensor6>(m, "Tensor6", py::buffer_protocol())
      .def(py::init([]() { return std::make_shared<Tensor6>(); }), "Default tensor")
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
            auto out = std::make_shared<Tensor6>(n1, n2, n3, n4, n5, n6);
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
          py::arg("ten6").none(false), py::doc("From :class:`list` equivalent"))
      .PythonInterfaceCopyValue(Tensor6)
      .PythonInterfaceWorkspaceVariableConversion(Tensor6)
      .PythonInterfaceBasicRepresentation(Tensor6)
      .PythonInterfaceFileIO(Tensor6)
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor6& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor6& x, Tensor6& y) { x = y; }, py::doc(":class:`~numpy.ndarray` Data array"))
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
The data can be accessed without copy using ``np.array(x, copy=False)`` or
via x.value
)--");

  artsclass<Tensor7>(m, "Tensor7", py::buffer_protocol())
      .def(py::init([]() { return std::make_shared<Tensor7>(); }), "Default tensor")
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
                 auto out = std::make_shared<Tensor7>(n1, n2, n3, n4, n5, n6, n7);
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
           py::arg("ten7").none(false), py::doc("From :class:`list` equivalent"))
      .PythonInterfaceCopyValue(Tensor7)
      .PythonInterfaceWorkspaceVariableConversion(Tensor7)
      .PythonInterfaceBasicRepresentation(Tensor7)
      .PythonInterfaceFileIO(Tensor7)
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor7& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor7& x, Tensor7& y) { x = y; }, py::doc(":class:`~numpy.ndarray` Data array"))
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
The data can be accessed without copy using ``np.array(x, copy=False)`` or
via x.value
)--");

  py::implicitly_convertible<std::vector<Scalar>, Vector>();
  py::implicitly_convertible<std::vector<std::vector<Scalar>>, Matrix>();
  py::implicitly_convertible<std::vector<std::vector<std::vector<Scalar>>>, Tensor3>();
  py::implicitly_convertible<std::vector<std::vector<std::vector<std::vector<Scalar>>>>, Tensor4>();
  py::implicitly_convertible<std::vector<std::vector<std::vector<std::vector<std::vector<Scalar>>>>>, Tensor5>();
  py::implicitly_convertible<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Scalar>>>>>>, Tensor6>();
  py::implicitly_convertible<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Scalar>>>>>>>, Tensor7>();

  
  artsarray<ArrayOfVector>(m, "ArrayOfVector")
      .PythonInterfaceFileIO(ArrayOfVector)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfVector);
  
  artsarray<ArrayOfMatrix>(m, "ArrayOfMatrix")
      .PythonInterfaceFileIO(ArrayOfMatrix)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfMatrix);
  
  artsarray<ArrayOfTensor3>(m, "ArrayOfTensor3")
      .PythonInterfaceFileIO(ArrayOfTensor3)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfTensor3);
  
  artsarray<ArrayOfTensor4>(m, "ArrayOfTensor4")
      .PythonInterfaceFileIO(ArrayOfTensor4)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfTensor4);
  
  artsarray<ArrayOfTensor5>(m, "ArrayOfTensor5")
      .PythonInterfaceFileIO(ArrayOfTensor5)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfTensor5);
  
  artsarray<ArrayOfTensor6>(m, "ArrayOfTensor6")
      .PythonInterfaceFileIO(ArrayOfTensor6)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfTensor6);
  
  artsarray<ArrayOfTensor7>(m, "ArrayOfTensor7")
      .PythonInterfaceFileIO(ArrayOfTensor7)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfTensor7);
  
  artsarray<ArrayOfArrayOfVector>(m, "ArrayOfArrayOfVector")
      .PythonInterfaceFileIO(ArrayOfArrayOfVector)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfArrayOfVector);
  
  artsarray<ArrayOfArrayOfMatrix>(m, "ArrayOfArrayOfMatrix")
      .PythonInterfaceFileIO(ArrayOfArrayOfMatrix)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfArrayOfMatrix);
  
  artsarray<ArrayOfArrayOfTensor3>(m, "ArrayOfArrayOfTensor3")
      .PythonInterfaceFileIO(ArrayOfArrayOfTensor3)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfArrayOfTensor3);
  
  artsarray<ArrayOfArrayOfTensor6>(m, "ArrayOfArrayOfTensor6")
      .PythonInterfaceFileIO(ArrayOfArrayOfTensor6)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfArrayOfTensor6);

  artsclass<Rational>(m, "Rational")
      .def(py::init([](Index n, Index d) {
             ARTS_USER_ERROR_IF(d < 1, "Must be positive")
             return Rational(n, d);
           }),
           py::arg("n") = 0,
           py::arg("d") = 1, py::doc("Default rational"))
      .PythonInterfaceCopyValue(Rational)
      .PythonInterfaceWorkspaceVariableConversion(Rational)
      .def(py::init([](const String& s) { return std::make_shared<Rational>(s); }), "From :class:`str`")
      .def(py::init([](Numeric n) { return Rational(std::to_string(n)); }), "From :class:`float`")
      .def("__float__", [](const Rational& x) { return Numeric(x); })
      .def("__int__", [](const Rational& x) { return Index(x); })
      .PythonInterfaceFileIO(Rational)
      .PythonInterfaceBasicRepresentation(Rational)
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

            return std::make_shared<Rational>(t[0].cast<Index>(), t[1].cast<Index>());
          }))
      .PythonInterfaceWorkspaceDocumentation(Rational);
  py::implicitly_convertible<Index, Rational>();

  artsclass<ComplexVector>(m, "ComplexVector", py::buffer_protocol())
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_property("value",
                    py::cpp_function(
                        [](ComplexVector& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](ComplexVector& x, ComplexVector& y) { x = y; }, py::doc(":class:`~numpy.ndarray` Data array"))
      .doc() = R"--(Holds complex vector data.
)--";
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize matpack\n", e.what()));
}
}  // namespace Python
