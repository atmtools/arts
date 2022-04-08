#include <py_auto_interface.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include "py_macros.h"

#define PythonInterfaceMatpackMath(Type)                          \
  def(                                                            \
      "__pos__",                                                  \
      [](const Type& a) {                                         \
        Type c = a;                                               \
        return c;                                                 \
      },                                                          \
      py::is_operator())                                          \
      .def(                                                       \
          "__neg__",                                              \
          [](const Type& a) {                                     \
            Type c = a;                                           \
            c *= -1;                                              \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
                                                                  \
      .def(                                                       \
          "__rpow__",                                             \
          [](const Type& a, Numeric_ b) {                         \
            Type c = a;                                           \
            c.transform_elementwise(                              \
                [b](Numeric x) { return std::pow(b, x); });       \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__rmul__",                                             \
          [](const Type& a, Numeric_ b) {                         \
            Type c = a;                                           \
            c *= b;                                               \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__radd__",                                             \
          [](const Type& a, Numeric_ b) {                         \
            Type c = a;                                           \
            c += b;                                               \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__rtruediv__",                                         \
          [](const Type& a, Numeric_ b) {                         \
            Type c = a;                                           \
            c = b;                                                \
            c /= a;                                               \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__rsub__",                                             \
          [](Type& a, Numeric_ b) {                               \
            Type c = a;                                           \
            c = b;                                                \
            c -= a;                                               \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
                                                                  \
      .def(                                                       \
          "__pow__",                                              \
          [](const Type& a, Numeric_ b) {                         \
            Type c = a;                                           \
            c.transform_elementwise(                              \
                [b](Numeric x) { return std::pow(x, b); });       \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__mul__",                                              \
          [](const Type& a, Numeric_ b) {                         \
            Type c = a;                                           \
            c *= b;                                               \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__add__",                                              \
          [](const Type& a, Numeric_ b) {                         \
            Type c = a;                                           \
            c += b;                                               \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__truediv__",                                          \
          [](const Type& a, Numeric_ b) {                         \
            Type c = a;                                           \
            c /= b;                                               \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__sub__",                                              \
          [](const Type& a, Numeric_ b) {                         \
            Type c = a;                                           \
            c -= b;                                               \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
                                                                  \
      .def(                                                       \
          "__mul__",                                              \
          [](const Type& a, const Type& b) {                      \
            ARTS_USER_ERROR_IF(a.shape() not_eq b.shape(),        \
                               "Invalid operation with shapes (", \
                               a.shape(),                         \
                               ") and (",                         \
                               b.shape(),                         \
                               ')')                               \
            Type c = a;                                           \
            c *= b;                                               \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__add__",                                              \
          [](const Type& a, const Type& b) {                      \
            ARTS_USER_ERROR_IF(a.shape() not_eq b.shape(),        \
                               "Invalid operation with shapes (", \
                               a.shape(),                         \
                               ") and (",                         \
                               b.shape(),                         \
                               ')')                               \
            Type c = a;                                           \
            c += b;                                               \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__truediv__",                                          \
          [](const Type& a, const Type& b) {                      \
            ARTS_USER_ERROR_IF(a.shape() not_eq b.shape(),        \
                               "Invalid operation with shapes (", \
                               a.shape(),                         \
                               ") and (",                         \
                               b.shape(),                         \
                               ')')                               \
            Type c = a;                                           \
            c /= b;                                               \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__sub__",                                              \
          [](const Type& a, const Type& b) {                      \
            ARTS_USER_ERROR_IF(a.shape() not_eq b.shape(),        \
                               "Invalid operation with shapes (", \
                               a.shape(),                         \
                               ") and (",                         \
                               b.shape(),                         \
                               ')')                               \
            Type c = a;                                           \
            c -= b;                                               \
            return c;                                             \
          },                                                      \
          py::is_operator())                                      \
                                                                  \
      .PythonInterfaceInPlaceMathOperators(Type, Numeric_)        \
                                                                  \
      .def(                                                       \
          "__imul__",                                             \
          [](Type& a, const Type& b) {                            \
            ARTS_USER_ERROR_IF(a.shape() not_eq b.shape(),        \
                               "Invalid operation with shapes (", \
                               a.shape(),                         \
                               ") and (",                         \
                               b.shape(),                         \
                               ')')                               \
            a *= b;                                               \
            return a;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__iadd__",                                             \
          [](Type& a, const Type& b) {                            \
            ARTS_USER_ERROR_IF(a.shape() not_eq b.shape(),        \
                               "Invalid operation with shapes (", \
                               a.shape(),                         \
                               ") and (",                         \
                               b.shape(),                         \
                               ')')                               \
            a += b;                                               \
            return a;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__itruediv__",                                         \
          [](Type& a, const Type& b) {                            \
            ARTS_USER_ERROR_IF(a.shape() not_eq b.shape(),        \
                               "Invalid operation with shapes (", \
                               a.shape(),                         \
                               ") and (",                         \
                               b.shape(),                         \
                               ')')                               \
            a /= b;                                               \
            return a;                                             \
          },                                                      \
          py::is_operator())                                      \
      .def(                                                       \
          "__isub__",                                             \
          [](Type& a, const Type& b) {                            \
            ARTS_USER_ERROR_IF(a.shape() not_eq b.shape(),        \
                               "Invalid operation with shapes (", \
                               a.shape(),                         \
                               ") and (",                         \
                               b.shape(),                         \
                               ')')                               \
            a -= b;                                               \
            return a;                                             \
          },                                                      \
          py::is_operator())

namespace Python {
void py_matpack(py::module_& m) {
  py::class_<ConstVectorView>(m,
                              "Const"
                              "Vector"
                              "View");
  py::class_<ConstMatrixView>(m,
                              "Const"
                              "Matrix"
                              "View");
  py::class_<ConstTensor3View>(m,
                               "Const"
                               "Tensor3"
                               "View");
  py::class_<ConstTensor4View>(m,
                               "Const"
                               "Tensor4"
                               "View");
  py::class_<ConstTensor5View>(m,
                               "Const"
                               "Tensor5"
                               "View");
  py::class_<ConstTensor6View>(m,
                               "Const"
                               "Tensor6"
                               "View");
  py::class_<ConstTensor7View>(m,
                               "Const"
                               "Tensor7"
                               "View");
  py::class_<VectorView, ConstVectorView>(m,
                                          "Vector"
                                          "View");
  py::class_<MatrixView, ConstMatrixView>(m,
                                          "Matrix"
                                          "View");
  py::class_<Tensor3View, ConstTensor3View>(m,
                                            "Tensor3"
                                            "View");
  py::class_<Tensor4View, ConstTensor4View>(m,
                                            "Tensor4"
                                            "View");
  py::class_<Tensor5View, ConstTensor5View>(m,
                                            "Tensor5"
                                            "View");
  py::class_<Tensor6View, ConstTensor6View>(m,
                                            "Tensor6"
                                            "View");
  py::class_<Tensor7View, ConstTensor7View>(m,
                                            "Tensor7"
                                            "View");

  py::class_<Vector, VectorView>(m, "Vector", py::buffer_protocol())
      .def(py::init([]() { return new Vector{}; }))
      .def(
          py::init([](const std::vector<Numeric>& v) { return new Vector{v}; }))
      .PythonInterfaceCopyValue(Vector)
      .PythonInterfaceWorkspaceVariableConversion(Vector)
      .PythonInterfaceMatpackMath(Vector)
      .def(
          "__matmul__",
          [](const Vector& a, const Vector& b) {
            ARTS_USER_ERROR_IF(a.shape() not_eq b.shape(),
                               "Invalid operation with shapes (",
                               a.shape(),
                               ") and (",
                               b.shape(),
                               ')')
            return a * b;
          },
          py::is_operator())
      .def_property_readonly("size", [](Vector& x) { return x.size(); })
      .def_property_readonly(
          "shape",
          [](Vector& x) { return x.shape().data; },
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
      .def_property("value",
                    py::cpp_function(
                        [](Vector& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Vector& x, Vector& y) { x = y; })
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

  py::class_<Matrix, MatrixView>(m, "Matrix", py::buffer_protocol())
      .def(py::init([]() { return new Matrix{}; }))
      .def(py::init([](const py::array_t<Numeric>& arr) {
        auto info = arr.request();
        if (info.ndim == 1 and info.shape[0] == 0) return Matrix{};

        Index i = 2;
        ARTS_USER_ERROR_IF(info.ndim not_eq i,
                           "Can only initialize from empty or ",
                           i,
                           "D array")
        std::size_t nc = info.shape[--i];
        std::size_t nr = info.shape[--i];

        auto* val = reinterpret_cast<Numeric*>(info.ptr);
        Matrix out(nr, nc);
        std::copy(val, val + out.size(), out.get_raw_data());

        return out;
      }))
      .PythonInterfaceCopyValue(Matrix)
      .PythonInterfaceWorkspaceVariableConversion(Matrix)
      .PythonInterfaceBasicRepresentation(Matrix)
      .PythonInterfaceFileIO(Matrix)
      .PythonInterfaceMatpackMath(Matrix)
      .def(
          "__matmul__",
          [](const Matrix& B, const Matrix& C) {
            ARTS_USER_ERROR_IF(B.ncols() not_eq C.nrows(),
                               "Invalid operation with shapes (",
                               B.shape(),
                               ") and (",
                               C.shape(),
                               ')')
            Matrix A(B.nrows(), C.ncols());
            mult(A, B, C);
            return A;
          },
          py::is_operator())
      .def(
          "__matmul__",
          [](const Matrix& B, const Vector& C) {
            ARTS_USER_ERROR_IF(B.ncols() not_eq C.nelem(),
                               "Invalid operation with shapes (",
                               B.shape(),
                               ") and (",
                               C.shape(),
                               ')')
            Vector A(B.nrows());
            mult(A, B, C);
            return A;
          },
          py::is_operator())
      .def_property_readonly(
          "T",
          [](const Matrix& x) { return Matrix(transpose(x)); },
          "Non-trivial transpose")
      .def_property_readonly("size", [](Matrix& x) { return x.size(); })
      .def_property_readonly(
          "shape",
          [](Matrix& x) { return x.shape().data; },
          "The shape of the data")
      .def(
          "__getitem__",
          [](Matrix& x, std::tuple<Index, Index> inds) -> Numeric& {
            auto [r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0)
              throw std::out_of_range("Out of bounds");
            return as_ref(x(r, c));
          },
          py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](Matrix& x, std::tuple<Index, Index> inds, Numeric_ y) {
            auto [r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0)
              throw std::out_of_range("Out of bounds");
            x(r, c) = y;
          },
          py::return_value_policy::reference_internal)
      .def_buffer([](Matrix& x) -> py::buffer_info {
        return py::buffer_info(x.get_c_array(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               2,
                               {x.nrows(), x.ncols()},
                               {sizeof(Numeric) * x.ncols(), sizeof(Numeric)});
      })
      .def_property("value",
                    py::cpp_function(
                        [](Matrix& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Matrix& x, Matrix& y) { x = y; })
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

  py::class_<Tensor3, Tensor3View>(m, "Tensor3", py::buffer_protocol())
      .def(py::init([]() { return new Tensor3{}; }))
      .def(py::init([](const py::array_t<Numeric>& arr) {
        auto info = arr.request();
        if (info.ndim == 1 and info.shape[0] == 0) return Tensor3{};

        Index i = 3;
        ARTS_USER_ERROR_IF(info.ndim not_eq i,
                           "Can only initialize from empty or ",
                           i,
                           "D array")
        std::size_t nc = info.shape[--i];
        std::size_t nr = info.shape[--i];
        std::size_t np = info.shape[--i];

        auto* val = reinterpret_cast<Numeric*>(info.ptr);
        Tensor3 out(np, nr, nc);
        std::copy(val, val + out.size(), out.get_c_array());

        return out;
      }))
      .PythonInterfaceCopyValue(Tensor3)
      .PythonInterfaceWorkspaceVariableConversion(Tensor3)
      .PythonInterfaceBasicRepresentation(Tensor3)
      .PythonInterfaceFileIO(Tensor3)
      .PythonInterfaceMatpackMath(Tensor3)
      .def_property_readonly("size", [](Tensor3& x) { return x.size(); })
      .def_property_readonly(
          "shape",
          [](Tensor3& x) { return x.shape().data; },
          "The shape of the data")
      .def(
          "__getitem__",
          [](Tensor3& x, std::tuple<Index, Index, Index> inds) -> Numeric& {
            auto [p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0)
              throw std::out_of_range("Out of bounds");
            return as_ref(x(p, r, c));
          },
          py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](Tensor3& x, std::tuple<Index, Index, Index> inds, Numeric_ y) {
            auto [p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0)
              throw std::out_of_range("Out of bounds");
            x(p, r, c) = y;
          },
          py::return_value_policy::reference_internal)
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor3& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor3& x, Tensor3& y) { x = y; })
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

  py::class_<Tensor4, Tensor4View>(m, "Tensor4", py::buffer_protocol())
      .def(py::init([]() { return new Tensor4{}; }))
      .def(py::init([](const py::array_t<Numeric>& arr) {
        auto info = arr.request();
        if (info.ndim == 1 and info.shape[0] == 0) return Tensor4{};

        Index i = 4;
        ARTS_USER_ERROR_IF(info.ndim not_eq i,
                           "Can only initialize from empty or ",
                           i,
                           "D array")
        std::size_t nc = info.shape[--i];
        std::size_t nr = info.shape[--i];
        std::size_t np = info.shape[--i];
        std::size_t nb = info.shape[--i];

        auto* val = reinterpret_cast<Numeric*>(info.ptr);
        Tensor4 out(nb, np, nr, nc);
        std::copy(val, val + out.size(), out.get_c_array());

        return out;
      }))
      .PythonInterfaceCopyValue(Tensor4)
      .PythonInterfaceWorkspaceVariableConversion(Tensor4)
      .PythonInterfaceBasicRepresentation(Tensor4)
      .PythonInterfaceFileIO(Tensor4)
      .PythonInterfaceMatpackMath(Tensor4)
      .def_property_readonly("size", [](Tensor4& x) { return x.size(); })
      .def_property_readonly(
          "shape",
          [](Tensor4& x) { return x.shape().data; },
          "The shape of the data")
      .def(
          "__getitem__",
          [](Tensor4& x,
             std::tuple<Index, Index, Index, Index> inds) -> Numeric& {
            auto [b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0)
              throw std::out_of_range("Out of bounds");
            return as_ref(x(b, p, r, c));
          },
          py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](Tensor4& x,
             std::tuple<Index, Index, Index, Index> inds,
             Numeric_ y) {
            auto [b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0)
              throw std::out_of_range("Out of bounds");
            x(b, p, r, c) = y;
          },
          py::return_value_policy::reference_internal)
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor4& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor4& x, Tensor4& y) { x = y; })
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

  py::class_<Tensor5, Tensor5View>(m, "Tensor5", py::buffer_protocol())
      .def(py::init([]() { return new Tensor5{}; }))
      .def(py::init([](const py::array_t<Numeric>& arr) {
        auto info = arr.request();
        if (info.ndim == 1 and info.shape[0] == 0) return Tensor5{};

        Index i = 5;
        ARTS_USER_ERROR_IF(info.ndim not_eq i,
                           "Can only initialize from empty or ",
                           i,
                           "D array")
        std::size_t nc = info.shape[--i];
        std::size_t nr = info.shape[--i];
        std::size_t np = info.shape[--i];
        std::size_t nb = info.shape[--i];
        std::size_t ns = info.shape[--i];

        auto* val = reinterpret_cast<Numeric*>(info.ptr);
        Tensor5 out(ns, nb, np, nr, nc);
        std::copy(val, val + out.size(), out.get_c_array());

        return out;
      }))
      .PythonInterfaceCopyValue(Tensor5)
      .PythonInterfaceWorkspaceVariableConversion(Tensor5)
      .PythonInterfaceBasicRepresentation(Tensor5)
      .PythonInterfaceFileIO(Tensor5)
      .PythonInterfaceMatpackMath(Tensor5)
      .def_property_readonly("size", [](Tensor5& x) { return x.size(); })
      .def_property_readonly(
          "shape",
          [](Tensor5& x) { return x.shape().data; },
          "The shape of the data")
      .def(
          "__getitem__",
          [](Tensor5& x,
             std::tuple<Index, Index, Index, Index, Index> inds) -> Numeric& {
            auto [s, b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                x.nshelves() <= s or s < 0)
              throw std::out_of_range("Out of bounds");
            return as_ref(x(s, b, p, r, c));
          },
          py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](Tensor5& x,
             std::tuple<Index, Index, Index, Index, Index> inds,
             Numeric_ y) {
            auto [s, b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                x.nshelves() <= s or s < 0)
              throw std::out_of_range("Out of bounds");
            x(s, b, p, r, c) = y;
          },
          py::return_value_policy::reference_internal)
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor5& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor5& x, Tensor5& y) { x = y; })
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

  py::class_<Tensor6, Tensor6View>(m, "Tensor6", py::buffer_protocol())
      .def(py::init([]() { return new Tensor6{}; }))
      .def(py::init([](const py::array_t<Numeric>& arr) {
        auto info = arr.request();
        if (info.ndim == 1 and info.shape[0] == 0) return Tensor6{};

        Index i = 6;
        ARTS_USER_ERROR_IF(info.ndim not_eq i,
                           "Can only initialize from empty or ",
                           i,
                           "D array")
        std::size_t nc = info.shape[--i];
        std::size_t nr = info.shape[--i];
        std::size_t np = info.shape[--i];
        std::size_t nb = info.shape[--i];
        std::size_t ns = info.shape[--i];
        std::size_t nv = info.shape[--i];

        auto* val = reinterpret_cast<Numeric*>(info.ptr);
        Tensor6 out(nv, ns, nb, np, nr, nc);
        std::copy(val, val + out.size(), out.get_c_array());

        return out;
      }))
      .PythonInterfaceCopyValue(Tensor6)
      .PythonInterfaceWorkspaceVariableConversion(Tensor6)
      .PythonInterfaceBasicRepresentation(Tensor6)
      .PythonInterfaceFileIO(Tensor6)
      .PythonInterfaceMatpackMath(Tensor6)
      .def_property_readonly("size", [](Tensor6& x) { return x.size(); })
      .def_property_readonly(
          "shape",
          [](Tensor6& x) { return x.shape().data; },
          "The shape of the data")
      .def(
          "__getitem__",
          [](Tensor6& x,
             std::tuple<Index, Index, Index, Index, Index, Index> inds)
              -> Numeric& {
            auto [v, s, b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                x.nshelves() <= s or s < 0 or x.nvitrines() <= v or v < 0)
              throw std::out_of_range("Out of bounds");
            return as_ref(x(v, s, b, p, r, c));
          },
          py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](Tensor6& x,
             std::tuple<Index, Index, Index, Index, Index, Index> inds,
             Numeric_ y) {
            auto [v, s, b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                x.nshelves() <= s or s < 0 or x.nvitrines() <= v or v < 0)
              throw std::out_of_range("Out of bounds");
            x(v, s, b, p, r, c) = y;
          },
          py::return_value_policy::reference_internal)
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor6& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor6& x, Tensor6& y) { x = y; })
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

  py::class_<Tensor7, Tensor7View>(m, "Tensor7", py::buffer_protocol())
      .def(py::init([]() { return new Tensor7{}; }))
      .def(py::init([](const py::array_t<Numeric>& arr) {
        auto info = arr.request();
        if (info.ndim == 1 and info.shape[0] == 0) return Tensor7{};

        Index i = 7;
        ARTS_USER_ERROR_IF(info.ndim not_eq i,
                           "Can only initialize from empty or ",
                           i,
                           "D array")
        std::size_t nc = info.shape[--i];
        std::size_t nr = info.shape[--i];
        std::size_t np = info.shape[--i];
        std::size_t nb = info.shape[--i];
        std::size_t ns = info.shape[--i];
        std::size_t nv = info.shape[--i];
        std::size_t nl = info.shape[--i];

        auto* val = reinterpret_cast<Numeric*>(info.ptr);
        Tensor7 out(nl, nv, ns, nb, np, nr, nc);
        std::copy(val, val + out.size(), out.get_c_array());

        return out;
      }))
      .PythonInterfaceCopyValue(Tensor7)
      .PythonInterfaceWorkspaceVariableConversion(Tensor7)
      .PythonInterfaceBasicRepresentation(Tensor7)
      .PythonInterfaceFileIO(Tensor7)
      .PythonInterfaceMatpackMath(Tensor7)
      .def_property_readonly("size", [](Tensor7& x) { return x.size(); })
      .def_property_readonly(
          "shape",
          [](Tensor7& x) { return x.shape().data; },
          "The shape of the data")
      .def(
          "__getitem__",
          [](Tensor7& x,
             std::tuple<Index, Index, Index, Index, Index, Index, Index> inds)
              -> Numeric& {
            auto [l, v, s, b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                x.nshelves() <= s or s < 0 or x.nvitrines() <= v or v < 0 or
                x.nlibraries() <= l or l < 0)
              throw std::out_of_range("Out of bounds");
            return as_ref(x(l, v, s, b, p, r, c));
          },
          py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](Tensor7& x,
             std::tuple<Index, Index, Index, Index, Index, Index, Index> inds,
             Numeric_ y) {
            auto [l, v, s, b, p, r, c] = inds;
            if (x.ncols() <= c or c < 0 or x.nrows() <= r or r < 0 or
                x.npages() <= p or p < 0 or x.nbooks() <= b or b < 0 or
                x.nshelves() <= s or s < 0 or x.nvitrines() <= v or v < 0 or
                x.nlibraries() <= l or l < 0)
              throw std::out_of_range("Out of bounds");
            x(l, v, s, b, p, r, c) = y;
          },
          py::return_value_policy::reference_internal)
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
      .def_property("value",
                    py::cpp_function(
                        [](Tensor7& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](Tensor7& x, Tensor7& y) { x = y; })
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
      }));
  py::implicitly_convertible<std::vector<std::vector<Vector>>,
                             ArrayOfArrayOfVector>();

  PythonInterfaceWorkspaceArray(ArrayOfMatrix)
      .def(py::init([](const std::vector<std::vector<Matrix>>& x) {
        ArrayOfArrayOfMatrix y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }));
  py::implicitly_convertible<std::vector<std::vector<Matrix>>,
                             ArrayOfArrayOfMatrix>();

  PythonInterfaceWorkspaceArray(ArrayOfTensor3)
      .def(py::init([](const std::vector<std::vector<Tensor3>>& x) {
        ArrayOfArrayOfTensor3 y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }));
  py::implicitly_convertible<std::vector<std::vector<Tensor3>>,
                             ArrayOfArrayOfTensor3>();

  PythonInterfaceWorkspaceArray(ArrayOfTensor6)
      .def(py::init([](const std::vector<std::vector<Tensor6>>& x) {
        ArrayOfArrayOfTensor6 y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }));
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
      .def(
          "__eq__",
          [](Rational& a, Rational b) { return a == b; },
          py::is_operator())
      .def(
          "__neq__",
          [](Rational& a, Rational b) { return a != b; },
          py::is_operator())
      .def(
          "__lt__",
          [](Rational& a, Rational b) { return a < b; },
          py::is_operator())
      .def(
          "__le__",
          [](Rational& a, Rational b) { return a <= b; },
          py::is_operator())
      .def(
          "__gt__",
          [](Rational& a, Rational b) { return a > b; },
          py::is_operator())
      .def(
          "__ge__",
          [](Rational& a, Rational b) { return a >= b; },
          py::is_operator());
  py::implicitly_convertible<Index, Rational>();
}
}  // namespace Python
