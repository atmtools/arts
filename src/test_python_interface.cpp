#include <tuple>
#include "array.h"
#include "auto_md.h"
#include "debug.h"
#include "python_interface.h"
#include "xml_io.h"

#define PythonInterfaceFileIO(Type)                                        \
  def(                                                                     \
      "savexml",                                                           \
      [](const Type& x,                                                    \
         const char* const file,                                           \
         const char* const type,                                           \
         bool clobber) {                                                   \
        xml_write_to_file(                                                 \
            file, x, string2filetype(type), clobber ? 0 : 1, Verbosity()); \
      },                                                                   \
      py::arg("file"),                                                     \
      py::arg("type") = "ascii",                                           \
      py::arg("clobber") = true)                                           \
      .def("loadxml", [](Type& x, const char* const file) {                \
        xml_read_from_file(file, x, Verbosity());                          \
      })

#define PythonInterfaceIndexItemAccess(Type)                       \
  def("__len__", [](const Type& x) { return x.nelem(); })          \
      .def("__getitem__",                                          \
           [](const Type& x, Index i) {                            \
             ARTS_USER_ERROR_IF(x.nelem() <= i or i < 0,           \
                                "Bad index access: ",              \
                                i,                                 \
                                " in object of size [0, ",         \
                                x.size(),                          \
                                ")")                               \
             return x[i];                                          \
           })                                                      \
      .def("__setitem__", [](Type& x, Index i, decltype(x[i]) y) { \
        ARTS_USER_ERROR_IF(x.nelem() <= i or i < 0,                \
                           "Bad index access: ",                   \
                           i,                                      \
                           " in object of size [0, ",              \
                           x.size(),                               \
                           ")")                                    \
        x[i] = std::move(y);                                       \
      })

#define PythonInterfaceBasicRepresentation(Type) \
  def("__str__", [](const Type& x) {             \
    return var_string(x);                        \
  }).def("__repr__", [](const Type& x) { return var_string(x); })

/** Gives the basic Array<X> interface of Arts
 * 
 * The type should be the full name of the class
 *
 * Python operators provided includes:
 *
 * x = Array<X>()
 * x = Array<X>(INDEX)
 * x = Array<X>(INDEX, X())
 * x.append(X())
 * list(x)
 * x[INDEX]
 * len(x)
 * for a in x: DO_SOMETHING(a)
 */
#define PythonInterfaceArrayDefault(Type)                                \
  PythonInterfaceIndexItemAccess(Type)                                   \
      .def(py::init<>())                                                 \
      .def(py::init<Index>())                                            \
      .def(py::init<Index, decltype(Type{}[0])>())                       \
      .def(                                                              \
          "__iter__",                                                    \
          [](Type& x) { return py::make_iterator(x.begin(), x.end()); }, \
          py::keep_alive<0, 1>())                                        \
      .def("append",                                                     \
           [](Type& x, decltype(x[0]) y) { x.emplace_back(std::move(y)); })

/** Provides -=, +=, /=. and *= for all LHS Type
 * 
 * For interoperability with python, the result of the operator is moved into
 * another Type() instance, which is returned.  The reason for this is that
 * we do not want the matpack to return views or other temporary objects being
 * returned that are not core types supported by Arts.  For most common types,
 * this should be a cost-less operator but not for all
 */
#define PythonInterfaceInPlaceMathOperators(Type, Other)                   \
  def(                                                                     \
      "__imul__",                                                          \
      [](Type& x, const Other& y) { return Type(std::move(x *= y)); },     \
      py::is_operator())                                                   \
      .def(                                                                \
          "__itruediv__",                                                  \
          [](Type& x, const Other& y) { return Type(std::move(x /= y)); }, \
          py::is_operator())                                               \
      .def(                                                                \
          "__iadd__",                                                      \
          [](Type& x, const Other& y) { return Type(std::move(x += y)); }, \
          py::is_operator())                                               \
      .def(                                                                \
          "__isub__",                                                      \
          [](Type& x, const Other& y) { return Type(std::move(x -= y)); }, \
          py::is_operator())

/** Provides -, +, /. and * for Type
 */
#define PythonInterfaceOperators(Type, Other)                       \
  def(                                                              \
      "__mul__",                                                    \
      [](Type& x, const Other& y) { return std::move(x * y); },     \
      py::is_operator())                                            \
      .def(                                                         \
          "__truediv__",                                            \
          [](Type& x, const Other& y) { return std::move(x / y); }, \
          py::is_operator())                                        \
      .def(                                                         \
          "__add__",                                                \
          [](Type& x, const Other& y) { return std::move(x + y); }, \
          py::is_operator())                                        \
      .def(                                                         \
          "__sub__",                                                \
          [](Type& x, const Other& y) { return std::move(x - y); }, \
          py::is_operator())

namespace Python {
PYBIND11_MODULE(pyarts_classes, m) {
  py::class_<Vector>(m, "Vector", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index>())
      .def(py::init<Index, Numeric>())
      .def(py::init<const std::vector<Numeric>&>())
      .def(py::init<Numeric, Index, Numeric>())
      .PythonInterfaceBasicRepresentation(Vector)
      .PythonInterfaceFileIO(Vector)
      .PythonInterfaceIndexItemAccess(Vector)
      .PythonInterfaceInPlaceMathOperators(Vector, Vector)
      .PythonInterfaceInPlaceMathOperators(Vector, Numeric)
      .def("__matmul__", [](const Vector& x, const Vector& y){return x*y;}, py::is_operator())
      .def(
          "__iter__",
          [](Vector& x) { return py::make_iterator(x.begin(), x.end()); },
          py::keep_alive<0, 1>())
      .def_buffer([](Vector& x) -> py::buffer_info {
        return py::buffer_info(x.get_c_array(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               1,
                               {x.nelem()},
                               {sizeof(Numeric)});
      });

  py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index, Index>())
      .def(py::init<Index, Index, Numeric>())
      .PythonInterfaceBasicRepresentation(Matrix)
      .PythonInterfaceFileIO(Matrix)
      .PythonInterfaceInPlaceMathOperators(Matrix, Matrix)
      .PythonInterfaceInPlaceMathOperators(Matrix, Numeric)
      .def("ncols", [](const Matrix& x) { return x.ncols(); })
      .def("nrows", [](const Matrix& x) { return x.nrows(); })
      .def("__getitem__",
           [](const Matrix& x, std::tuple<Index, Index> inds) {
             auto [r, c] = inds;
             ARTS_USER_ERROR_IF(
                 x.ncols() <= c or c < 0 or
                 x.nrows() <= r or r < 0,
                 "Bad index access ",
                 '[', r, ", ", c, ']',
                 " in object of shape ",
                 '(', x.nrows(), ", ", x.ncols(), ')')
             return x(r, c);
           })
      .def("__setitem__",
           [](Matrix& x, std::tuple<Index, Index> inds, Numeric y) {
             auto [r, c] = inds;
             ARTS_USER_ERROR_IF(
                 x.ncols() <= c or c < 0 or
                 x.nrows() <= r or r < 0,
                 "Bad index access ",
                 '[', r, ", ", c, ']',
                 " in object of size ",
                 '[', x.nrows(), ", ", x.ncols(), ']')
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
      .PythonInterfaceBasicRepresentation(Tensor3)
      .PythonInterfaceFileIO(Tensor3)
      .PythonInterfaceInPlaceMathOperators(Tensor3, Tensor3)
      .PythonInterfaceInPlaceMathOperators(Tensor3, Numeric)
      .def("ncols", [](const Tensor3& x) { return x.ncols(); })
      .def("nrows", [](const Tensor3& x) { return x.nrows(); })
      .def("npages", [](const Tensor3& x) { return x.npages(); })
      .def("__getitem__",
           [](const Tensor3& x, std::tuple<Index, Index, Index> inds) {
             auto [p, r, c] = inds;
             ARTS_USER_ERROR_IF(
                 x.ncols() <= c or c < 0 or
                 x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0,
                 "Bad index access ",
                 '[', p, ", ", r, ", ", c, ']',
                 " in object of shape ",
                 '(', x.npages(), ", ", x.nrows(), ", ", x.ncols(), ')')
             return x(p, r, c);
           })
      .def("__setitem__",
           [](Tensor3& x, std::tuple<Index, Index, Index> inds, Numeric y) {
             auto [p, r, c] = inds;
             ARTS_USER_ERROR_IF(
                 x.ncols() <= c or c < 0 or
                 x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0,
                 "Bad index access ",
                 '[', p, ", ", r, ", ", c, ']',
                 " in object of shape ",
                 '(', x.npages(), ", ", x.nrows(), ", ", x.ncols(), ')')
             x(p, r, c) = y;
           })
      .def_buffer([](Tensor3& x) -> py::buffer_info {
        return py::buffer_info(x.get_c_array(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               3,
                               {x.npages(), x.nrows(), x.ncols()},
                               {sizeof(Numeric) * x.ncols() * x.nrows(), sizeof(Numeric) * x.ncols(), sizeof(Numeric)});
      });

  py::class_<Tensor4>(m, "Tensor4", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index, Index, Index, Index>())
      .def(py::init<Index, Index, Index, Index, Numeric>())
      .PythonInterfaceBasicRepresentation(Tensor4)
      .PythonInterfaceFileIO(Tensor4)
      .PythonInterfaceInPlaceMathOperators(Tensor4, Tensor4)
      .PythonInterfaceInPlaceMathOperators(Tensor4, Numeric)
      .def("ncols", [](const Tensor4& x) { return x.ncols(); })
      .def("nrows", [](const Tensor4& x) { return x.nrows(); })
      .def("npages", [](const Tensor4& x) { return x.npages(); })
      .def("nbooks", [](const Tensor4& x) { return x.nbooks(); })
      .def("__getitem__",
           [](const Tensor4& x, std::tuple<Index, Index, Index, Index> inds) {
             auto [b, p, r, c] = inds;
             ARTS_USER_ERROR_IF(
                 x.ncols() <= c or c < 0 or
                 x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0 or
                 x.nbooks() <= b or b < 0,
                 "Bad index access ",
                 '[', b, ", ", p, ", ", r, ", ", c, ']',
                 " in object of shape ",
                 '(', x.nbooks(), ", ", x.npages(), ", ", x.nrows(), ", ", x.ncols(), ')')
             return x(b, p, r, c);
           })
      .def("__setitem__",
           [](Tensor4& x, std::tuple<Index, Index, Index, Index> inds, Numeric y) {
             auto [b, p, r, c] = inds;
             ARTS_USER_ERROR_IF(
                 x.ncols() <= c or c < 0 or
                 x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0 or
                 x.nbooks() <= b or b < 0,
                 "Bad index access ",
                 '[', b, ", ", p, ", ", r, ", ", c, ']',
                 " in object of shape ",
                 '(', x.nbooks(), ", ", x.npages(), ", ", x.nrows(), ", ", x.ncols(), ')')
             x(b, p, r, c) = y;
           })
      .def_buffer([](Tensor4& x) -> py::buffer_info {
        return py::buffer_info(x.get_c_array(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               4,
                               {x.nbooks(), x.npages(), x.nrows(), x.ncols()},
                               {sizeof(Numeric) * x.npages() * x.ncols() * x.nrows(), sizeof(Numeric) * x.ncols() * x.nrows(), sizeof(Numeric) * x.ncols(), sizeof(Numeric)});
      });

  py::class_<Tensor5>(m, "Tensor5", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index, Index, Index, Index, Index>())
      .def(py::init<Index, Index, Index, Index, Index, Numeric>())
      .PythonInterfaceBasicRepresentation(Tensor5)
      .PythonInterfaceFileIO(Tensor5)
      .PythonInterfaceInPlaceMathOperators(Tensor5, Tensor5)
      .PythonInterfaceInPlaceMathOperators(Tensor5, Numeric)
      .def("ncols", [](const Tensor5& x) { return x.ncols(); })
      .def("nrows", [](const Tensor5& x) { return x.nrows(); })
      .def("npages", [](const Tensor5& x) { return x.npages(); })
      .def("nbooks", [](const Tensor5& x) { return x.nbooks(); })
      .def("nshelves", [](const Tensor5& x) { return x.nshelves(); })
      .def("__getitem__",
           [](const Tensor5& x, std::tuple<Index, Index, Index, Index, Index> inds) {
             auto [s, b, p, r, c] = inds;
             ARTS_USER_ERROR_IF(
                 x.ncols() <= c or c < 0 or
                 x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0 or
                 x.nbooks() <= b or b < 0 or
                 x.nshelves() <= s or s < 0,
                 "Bad index access ",
                 '[', s, ", ", b, ", ", p, ", ", r, ", ", c, ']',
                 " in object of shape ",
                 '(', x.nshelves(), ", ", x.nbooks(), ", ", x.npages(), ", ", x.nrows(), ", ", x.ncols(), ')')
             return x(s, b, p, r, c);
           })
      .def("__setitem__",
           [](Tensor5& x, std::tuple<Index, Index, Index, Index, Index> inds, Numeric y) {
             auto [s, b, p, r, c] = inds;
             ARTS_USER_ERROR_IF(
                 x.ncols() <= c or c < 0 or
                 x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0 or
                 x.nbooks() <= b or b < 0 or
                 x.nshelves() <= s or s < 0,
                 "Bad index access ",
                 '[', s, ", ", b, ", ", p, ", ", r, ", ", c, ']',
                 " in object of shape ",
                 '(', x.nshelves(), ", ", x.nbooks(), ", ", x.npages(), ", ", x.nrows(), ", ", x.ncols(), ')')
             x(s, b, p, r, c) = y;
           })
      .def_buffer([](Tensor5& x) -> py::buffer_info {
        return py::buffer_info(x.get_c_array(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               5,
                               {x.nshelves(), x.nbooks(), x.npages(), x.nrows(), x.ncols()},
                               {sizeof(Numeric) * x.nbooks() * x.npages() * x.ncols() * x.nrows(), sizeof(Numeric) * x.npages() * x.ncols() * x.nrows(), sizeof(Numeric) * x.ncols() * x.nrows(), sizeof(Numeric) * x.ncols(), sizeof(Numeric)});
      });

  py::class_<Tensor6>(m, "Tensor6", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index, Index, Index, Index, Index, Index>())
      .def(py::init<Index, Index, Index, Index, Index, Index, Numeric>())
      .PythonInterfaceBasicRepresentation(Tensor6)
      .PythonInterfaceFileIO(Tensor6)
      .PythonInterfaceInPlaceMathOperators(Tensor6, Tensor6)
      .PythonInterfaceInPlaceMathOperators(Tensor6, Numeric)
      .def("ncols", [](const Tensor6& x) { return x.ncols(); })
      .def("nrows", [](const Tensor6& x) { return x.nrows(); })
      .def("npages", [](const Tensor6& x) { return x.npages(); })
      .def("nbooks", [](const Tensor6& x) { return x.nbooks(); })
      .def("nshelves", [](const Tensor6& x) { return x.nshelves(); })
      .def("nvitrines", [](const Tensor6& x) { return x.nvitrines(); })
      .def("__getitem__",
           [](const Tensor6& x, std::tuple<Index, Index, Index, Index, Index, Index> inds) {
             auto [v, s, b, p, r, c] = inds;
             ARTS_USER_ERROR_IF(
                 x.ncols() <= c or c < 0 or
                 x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0 or
                 x.nbooks() <= b or b < 0 or
                 x.nshelves() <= s or s < 0 or
                 x.nvitrines() <= v or v < 0,
                 "Bad index access ",
                 '[', v, ", ", s, ", ", b, ", ", p, ", ", r, ", ", c, ']',
                 " in object of shape ",
                 '(', x.nvitrines(), ", ", x.nshelves(), ", ", x.nbooks(), ", ", x.npages(), ", ", x.nrows(), ", ", x.ncols(), ')')
             return x(v, s, b, p, r, c);
           })
      .def("__setitem__",
           [](Tensor6& x, std::tuple<Index, Index, Index, Index, Index, Index> inds, Numeric y) {
             auto [v, s, b, p, r, c] = inds;
             ARTS_USER_ERROR_IF(
                 x.ncols() <= c or c < 0 or
                 x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0 or
                 x.nbooks() <= b or b < 0 or
                 x.nshelves() <= s or s < 0 or
                 x.nvitrines() <= v or v < 0,
                 "Bad index access ",
                 '[', v, ", ", s, ", ", b, ", ", p, ", ", r, ", ", c, ']',
                 " in object of shape ",
                 '(', x.nvitrines(), ", ", x.nshelves(), ", ", x.nbooks(), ", ", x.npages(), ", ", x.nrows(), ", ", x.ncols(), ')')
             x(v, s, b, p, r, c) = y;
           })
      .def_buffer([](Tensor6& x) -> py::buffer_info {
        return py::buffer_info(x.get_c_array(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               6,
                               {x.nvitrines(), x.nshelves(), x.nbooks(), x.npages(), x.nrows(), x.ncols()},
                               {sizeof(Numeric) * x.nshelves() * x.nbooks() * x.npages() * x.ncols() * x.nrows(), sizeof(Numeric) * x.nbooks() * x.npages() * x.ncols() * x.nrows(), sizeof(Numeric) * x.npages() * x.ncols() * x.nrows(), sizeof(Numeric) * x.ncols() * x.nrows(), sizeof(Numeric) * x.ncols(), sizeof(Numeric)});
      });

  py::class_<Tensor7>(m, "Tensor7", py::buffer_protocol())
      .def(py::init<>())
      .def(py::init<Index, Index, Index, Index, Index, Index, Index>())
      .def(py::init<Index, Index, Index, Index, Index, Index, Index, Numeric>())
      .PythonInterfaceBasicRepresentation(Tensor7)
      .PythonInterfaceFileIO(Tensor7)
      .PythonInterfaceInPlaceMathOperators(Tensor7, Tensor7)
      .PythonInterfaceInPlaceMathOperators(Tensor7, Numeric)
      .def("ncols", [](const Tensor7& x) { return x.ncols(); })
      .def("nrows", [](const Tensor7& x) { return x.nrows(); })
      .def("npages", [](const Tensor7& x) { return x.npages(); })
      .def("nbooks", [](const Tensor7& x) { return x.nbooks(); })
      .def("nshelves", [](const Tensor7& x) { return x.nshelves(); })
      .def("nvitrines", [](const Tensor7& x) { return x.nvitrines(); })
      .def("nlibraries", [](const Tensor7& x) { return x.nlibraries(); })
      .def("__getitem__",
           [](const Tensor7& x, std::tuple<Index, Index, Index, Index, Index, Index, Index> inds) {
             auto [l, v, s, b, p, r, c] = inds;
             ARTS_USER_ERROR_IF(
                 x.ncols() <= c or c < 0 or
                 x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0 or
                 x.nbooks() <= b or b < 0 or
                 x.nshelves() <= s or s < 0 or
                 x.nvitrines() <= v or v < 0 or
                 x.nlibraries() <= l or l < 0,
                 "Bad index access ",
                 '[', l, ", ", v, ", ", s, ", ", b, ", ", p, ", ", r, ", ", c, ']',
                 " in object of shape ",
                 '(', x.nlibraries(), ", ", x.nvitrines(), ", ", x.nshelves(), ", ", x.nbooks(), ", ", x.npages(), ", ", x.nrows(), ", ", x.ncols(), ')')
             return x(l, v, s, b, p, r, c);
           })
      .def("__setitem__",
           [](Tensor7& x, std::tuple<Index, Index, Index, Index, Index, Index, Index> inds, Numeric y) {
             auto [l, v, s, b, p, r, c] = inds;
             ARTS_USER_ERROR_IF(
                 x.ncols() <= c or c < 0 or
                 x.nrows() <= r or r < 0 or
                 x.npages() <= p or p < 0 or
                 x.nbooks() <= b or b < 0 or
                 x.nshelves() <= s or s < 0 or
                 x.nvitrines() <= v or v < 0 or
                 x.nlibraries() <= l or l < 0,
                 "Bad index access ",
                 '[', l, ", ", v, ", ", s, ", ", b, ", ", p, ", ", r, ", ", c, ']',
                 " in object of shape ",
                 '(', x.nlibraries(), ", ", x.nvitrines(), ", ", x.nshelves(), ", ", x.nbooks(), ", ", x.npages(), ", ", x.nrows(), ", ", x.ncols(), ')')
             x(l, v, s, b, p, r, c) = y;
           })
      .def_buffer([](Tensor7& x) -> py::buffer_info {
        return py::buffer_info(x.get_c_array(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               7,
                               {x.nlibraries(), x.nvitrines(), x.nshelves(), x.nbooks(), x.npages(), x.nrows(), x.ncols()},
                               {sizeof(Numeric) * x.nvitrines() * x.nbooks() * x.npages() * x.ncols() * x.nrows(), sizeof(Numeric) * x.nshelves() * x.nbooks() * x.npages() * x.ncols() * x.nrows(), sizeof(Numeric) * x.nbooks() * x.npages() * x.ncols() * x.nrows(), sizeof(Numeric) * x.npages() * x.ncols() * x.nrows(), sizeof(Numeric) * x.ncols() * x.nrows(), sizeof(Numeric) * x.ncols(), sizeof(Numeric)});
      });

  py::class_<ArrayOfVector>(m, "ArrayOfVector")
      .PythonInterfaceFileIO(ArrayOfVector)
      .PythonInterfaceArrayDefault(ArrayOfVector)
      .PythonInterfaceBasicRepresentation(ArrayOfVector);

  py::class_<ArrayOfArrayOfVector>(m, "ArrayOfArrayOfVector")
      .PythonInterfaceFileIO(ArrayOfArrayOfVector)
      .PythonInterfaceArrayDefault(ArrayOfArrayOfVector)
      .PythonInterfaceBasicRepresentation(ArrayOfArrayOfVector);

  py::class_<Verbosity>(m, "Verbosity")
      .def(py::init<>())
      .def(py::init<Index, Index, Index>())
      .PythonInterfaceFileIO(Verbosity)
      .PythonInterfaceBasicRepresentation(Verbosity)
      .def("get_agenda_verbosity", &Verbosity::get_agenda_verbosity)
      .def("get_screen_verbosity", &Verbosity::get_screen_verbosity)
      .def("get_file_verbosity", &Verbosity::get_file_verbosity)
      .def("is_main_agenda", &Verbosity::is_main_agenda)
      .def("set_agenda_verbosity", &Verbosity::set_agenda_verbosity)
      .def("set_screen_verbosity", &Verbosity::set_screen_verbosity)
      .def("set_file_verbosity", &Verbosity::set_file_verbosity)
      .def("set_main_agenda", &Verbosity::set_main_agenda);

  py::class_<String>(m, "String")
      .def(py::init<>())
      .def(py::init<const char* const>())
      .PythonInterfaceFileIO(String)
      .def(
          "__iter__",
          [](String& x) { return py::make_iterator(x.begin(), x.end()); },
          py::keep_alive<0, 1>())
      .PythonInterfaceBasicRepresentation(String);
  py::implicitly_convertible<const char * const, String>();

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
}
}  // namespace Python

#undef PythonInterfaceFileIO
#undef PythonInterfaceIndexItemAccess
#undef PythonInterfaceBasicRepresentation
#undef PythonInterfaceArrayDefault
