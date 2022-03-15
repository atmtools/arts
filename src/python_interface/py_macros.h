#ifndef py_macros_h
#define py_macros_h

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <debug.h>
#include <xml_io.h>

#include "python_interface_value_type.h"

namespace Python { namespace py = pybind11; }

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
      py::arg("file").none(false),                                         \
      py::arg("type").none(false) = "ascii",                               \
      py::arg("clobber") = true,                                           \
      py::doc("Saves " #Type " to file\n"                                  \
              "\n"                                                         \
              "Parameters:\n"                                              \
              "    file (str): The path to which the file is written."     \
              " Note that several of the options might modify the"         \
              " name or write more files\n"                                \
              "    type (str): Type of file to save (ascii. zascii,"       \
              " or binary)\n"                                              \
              "    clobber (bool): Overwrite existing files or add new"    \
              " file with modified name?\n"                                \
              "\n"                                                         \
              "On Error:\n"                                                \
              "    Throws RuntimeError for any failure to save"))          \
      .def(                                                                \
          "readxml",                                                       \
          [](Type& x, const char* const file) {                            \
            xml_read_from_file(file, x, Verbosity());                      \
          },                                                               \
          py::arg("file").none(false),                                     \
          py::doc("Read " #Type " from file\n"                             \
                  "\n"                                                     \
                  "Parameters:\n"                                          \
                  "    file (str): A file that can be read\n"              \
                  "\n"                                                     \
                  "On Error:\n"                                            \
                  "    Throws RuntimeError for any failure to read"))

#define PythonInterfaceIndexItemAccess(Type)                                \
  def("__len__", [](const Type& x) { return x.nelem(); })                   \
      .def(                                                                 \
          "__getitem__",                                                    \
          [](Type& x, Index i) -> decltype(x[i])& {                         \
            if (x.nelem() <= i or i < 0)                                    \
              throw std::out_of_range(var_string("Bad index access: ",      \
                                                 i,                         \
                                                 " in object of size [0, ", \
                                                 x.size(),                  \
                                                 ")"));                     \
            return x[i];                                                    \
          },                                                                \
          py::return_value_policy::reference_internal)                      \
      .def(                                                                 \
          "__setitem__",                                                    \
          [](Type& x, Index i, decltype(x[i]) y) {                          \
            if (x.nelem() <= i or i < 0)                                    \
              throw std::out_of_range(var_string("Bad index access: ",      \
                                                 i,                         \
                                                 " in object of size [0, ", \
                                                 x.size(),                  \
                                                 ")"));                     \
            x[i] = std::move(y);                                            \
          },                                                                \
          py::return_value_policy::reference_internal)

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
 * y = x.pop()
 * list(x)
 * x[INDEX]
 * x[INDEX] = X()
 * len(x)
 * for a in x: DO_SOMETHING(a)
 */
#define PythonInterfaceArrayDefault(BaseType)                                 \
  PythonInterfaceIndexItemAccess(Array<BaseType>)                             \
      .def(py::init<>())                                                      \
      .def(py::init<Index>())                                                 \
      .def(py::init<Index, BaseType>())                                       \
      .def(py::init<std::vector<BaseType>>())                                 \
      .def(                                                                   \
          "append",                                                           \
          [](Array<BaseType>& x, BaseType y) {                                \
            x.emplace_back(std::move(y));                                     \
          },                                                                  \
          py::doc("Appends a " #BaseType " at the end of the Array"))         \
      .def(                                                                   \
          "pop",                                                              \
          [](Array<BaseType>& x) -> BaseType {                                \
            if (x.size() < 1) throw std::out_of_range("pop from empty list"); \
            BaseType copy = std::move(x.back());                              \
            x.pop_back();                                                     \
            return copy;                                                      \
          },                                                                  \
          py::doc("Pops a " #BaseType " from the Array"))

/** Provides -=, +=, /=. and *= for all LHS Type */
#define PythonInterfaceInPlaceMathOperators(Type, Other) \
  def(                                                   \
      "__imul__",                                        \
      [](Type& x, const Other& y) {                      \
        x *= y;                                          \
        return x;                                        \
      },                                                 \
      py::is_operator())                                 \
      .def(                                              \
          "__itruediv__",                                \
          [](Type& x, const Other& y) {                  \
            x /= y;                                      \
            return x;                                    \
          },                                             \
          py::is_operator())                             \
      .def(                                              \
          "__iadd__",                                    \
          [](Type& x, const Other& y) {                  \
            x += y;                                      \
            return x;                                    \
          },                                             \
          py::is_operator())                             \
      .def(                                              \
          "__isub__",                                    \
          [](Type& x, const Other& y) {                  \
            x -= y;                                      \
            return x;                                    \
          },                                             \
          py::is_operator())

/** Provides -, +, /. and * for Type (right and left)
 */
#define PythonInterfaceMathOperators(Type, Other)        \
  def(                                                   \
      "__mul__",                                         \
      [](Type& x, const Other& y) { return x * y; },     \
      py::is_operator())                                 \
      .def(                                              \
          "__truediv__",                                 \
          [](Type& x, const Other& y) { return x / y; }, \
          py::is_operator())                             \
      .def(                                              \
          "__add__",                                     \
          [](Type& x, const Other& y) { return x + y; }, \
          py::is_operator())                             \
      .def(                                              \
          "__sub__",                                     \
          [](Type& x, const Other& y) { return x - y; }, \
          py::is_operator())                             \
      .def(                                              \
          "__rmul__",                                    \
          [](Type& x, const Other& y) { return y * x; }, \
          py::is_operator())                             \
      .def(                                              \
          "__rtruediv__",                                \
          [](Type& x, const Other& y) { return y / x; }, \
          py::is_operator())                             \
      .def(                                              \
          "__radd__",                                    \
          [](Type& x, const Other& y) { return y + x; }, \
          py::is_operator())                             \
      .def(                                              \
          "__rsub__",                                    \
          [](Type& x, const Other& y) { return y - x; }, \
          py::is_operator())

#define PythonInterfaceWorkspaceVariableConversion(Type) \
  def(py::init([](const WorkspaceVariable* w) {          \
        Type& v = *w;                                    \
        return new Type(v);                              \
      }),                                                \
      py::arg("wsv").none(false).noconvert(),            \
      py::doc("Automatic conversion from a workspace variable"))

/*! The workspace array interface

It takes the base type as input, e.g., ArrayOfIndex (Array<Index>) takes
"Index" as input.  This allows setting up the default conversion from
std::vector<BaseType> to be created.

The main limitation for this is that the type must be names ArrayOfBaseType
in Arts.  Without this limitation, the representation and automatic documentation
that is generated by this macro would show the C++ type name rather than the
desired python name.  "ArrayOfBaseType" is the class exposed to python
*/
#define PythonInterfaceWorkspaceArray(BaseType)                           \
  auto auto_impl_name##BaseType =                                         \
      py::class_<ArrayOf##BaseType>(m, "ArrayOf" #BaseType);              \
  auto_impl_name##BaseType.doc() = "The Arts ArrayOf" #BaseType " class"; \
                                                                          \
  py::implicitly_convertible<std::vector<BaseType>, ArrayOf##BaseType>(); \
  py::implicitly_convertible<py::list, ArrayOf##BaseType>();              \
                                                                          \
  auto_impl_name##BaseType.PythonInterfaceFileIO(ArrayOf##BaseType)       \
      .PythonInterfaceBasicRepresentation(ArrayOf##BaseType)              \
      .PythonInterfaceArrayDefault(BaseType)                              \
      .PythonInterfaceWorkspaceVariableConversion(ArrayOf##BaseType)

#define PythonInterfaceGriddedField(Type)                                \
  def_readwrite("data", &Type::data)                                     \
      .def("get_dim", [](Type& gf) { return gf.get_dim(); })             \
      .def("get_name", [](Type& gf) { return gf.get_name(); })           \
      .def("set_name",                                                   \
           [](Type& gf, const String& str) { return gf.set_name(str); }) \
      .def("get_grid_name",                                              \
           [](Type& gf, Index i) {                                       \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.get_grid_name(i);                                 \
           })                                                            \
      .def("get_grid_size",                                              \
           [](Type& gf, Index i) {                                       \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.get_grid_size(i);                                 \
           })                                                            \
      .def("get_grid_type",                                              \
           [](Type& gf, Index i) {                                       \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.get_grid_type(i);                                 \
           })                                                            \
      .def("get_numeric_grid",                                           \
           [](Type& gf, Index i) {                                       \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.get_numeric_grid(i);                              \
           })                                                            \
      .def("get_string_grid",                                            \
           [](Type& gf, Index i) {                                       \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.get_string_grid(i);                               \
           })                                                            \
      .def("set_grid_name",                                              \
           [](Type& gf, Index i, const String& str) {                    \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.set_grid_name(i, str);                            \
           })                                                            \
      .def("set_grid",                                                   \
           [](Type& gf, Index i, const Vector& vec) {                    \
             ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")      \
             return gf.set_grid(i, vec);                                 \
           })                                                            \
      .def("set_grid", [](Type& gf, Index i, const ArrayOfString& str) { \
        ARTS_USER_ERROR_IF(i >= gf.get_dim(), "Out of bounds")           \
        return gf.set_grid(i, str);                                      \
      });

#define PythonInterfaceReadWriteData(Type, data) \
  def_readwrite(#data, &Type::data, py::return_value_policy::reference_internal)

#endif
