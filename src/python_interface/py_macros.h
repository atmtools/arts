#ifndef py_macros_h
#define py_macros_h

#include <python_interface.h>
#include <xml_io.h>

#include <functional>

#include "pydocs.h"
#include "python_interface_value_type.h"

namespace Python {
namespace py = pybind11;
}  // namespace Python

constexpr Index negative_clamp(const Index i, const Index n) noexcept {
  return (i < 0) ? i + n : i;
}

#define PythonInterfaceFileIO(Type)                                         \
  def(                                                                      \
      "savexml",                                                            \
      [](const Type& x,                                                     \
         const char* const file,                                            \
         const char* const type,                                            \
         bool clobber) {                                                    \
        xml_write_to_file(file, x, string2filetype(type), clobber ? 0 : 1); \
      },                                                                    \
      py::arg("file").none(false),                                          \
      py::arg("type").none(false) = "ascii",                                \
      py::arg("clobber") = true,                                            \
      py::doc("Saves :class:`" #Type "` to file\n"                          \
              "\n"                                                          \
              "Parameters:\n"                                               \
              "    file (str): The path to which the file is written."      \
              " Note that several of the options might modify the"          \
              " name or write more files\n"                                 \
              "    type (str): Type of file to save (ascii. zascii,"        \
              " or binary)\n"                                               \
              "    clobber (bool): Overwrite existing files or add new"     \
              " file with modified name?\n"                                 \
              "\n"                                                          \
              "On Error:\n"                                                 \
              "    Throws RuntimeError for any failure to save"))           \
      .def(                                                                 \
          "readxml",                                                        \
          [](Type& x, const char* const file) {                             \
            xml_read_from_file(file, x);                                    \
          },                                                                \
          py::arg("file").none(false),                                      \
          py::doc("Read :class:`" #Type "` from file\n"                     \
                  "\n"                                                      \
                  "Parameters:\n"                                           \
                  "    file (str): A file that can be read\n"               \
                  "\n"                                                      \
                  "On Error:\n"                                             \
                  "    Throws RuntimeError for any failure to read"))       \
      .def_static(                                                          \
          "fromxml",                                                        \
          [](const char* const file) {                                      \
            Type x;                                                         \
            xml_read_from_file(file, x);                                    \
            return x;                                                       \
          },                                                                \
          py::arg("file").none(false),                                      \
          py::doc("Create :class:`" #Type "` from file\n"                   \
                  "\n"                                                      \
                  "Parameters:\n"                                           \
                  "    file (str): A file that can be read\n"               \
                  "\n"                                                      \
                  "On Error:\n"                                             \
                  "    Throws RuntimeError for any failure to read"))

#define PythonInterfaceIndexItemAccess(Type)                                 \
  def("__len__", [](const Type& x) { return x.size(); })                     \
      .def(                                                                  \
          "__getitem__",                                                     \
          [](Type& x, Index i)                                               \
              -> std::shared_ptr<std::remove_cvref_t<decltype(x[i])>> {      \
            i = negative_clamp(i, x.size());                                 \
            if (x.size() <= static_cast<Size>(i) or i < 0)                   \
              throw std::out_of_range(var_string("Bad index access: ",       \
                                                 i,                          \
                                                 " in object of range [0, ", \
                                                 x.size(),                   \
                                                 ")"));                      \
            return std::shared_ptr<std::remove_cvref_t<decltype(x[i])>>(     \
                &x[i], [](void*) {});                                        \
          },                                                                 \
          py::return_value_policy::reference_internal,                       \
          py::keep_alive<0, 1>())                                            \
      .def("__setitem__", [](Type& x, Index i, decltype(x[i]) y) {           \
        i = negative_clamp(i, x.size());                                     \
        if (x.size() <= static_cast<Size>(i) or i < 0)                       \
          throw std::out_of_range(var_string("Bad index access: ",           \
                                             i,                              \
                                             " in object of range [0, ",     \
                                             x.size(),                       \
                                             ")"));                          \
        x[i] = std::move(y);                                                 \
      })

#define PythonInterfaceBasicRepresentation(Type) \
  def("__str__", [](const Type& x) {             \
    return var_string(x);                        \
  }).def("__repr__", [](const Type& x) { return var_string(x); })

#define PythonInterfaceCopyValue(Type)  \
  def("__copy__", [](Type& t) -> Type { \
    return t;                           \
  }).def("__deepcopy__", [](Type& t, py::dict&) -> Type { return t; })

#define PythonInterfaceCommonMath(Type) \
  def(py::self + Type())                \
      .def(py::self - Type())           \
      .def(py::self* Type())            \
      .def(py::self / Type())           \
      .def(py::self += Type())          \
      .def(py::self -= Type())          \
      .def(py::self *= Type())          \
      .def(py::self /= Type())          \
      .def(py::self == Type())          \
      .def(py::self != Type())          \
      .def(py::self <= Type())          \
      .def(py::self < Type())           \
      .def(py::self >= Type())          \
      .def(py::self > Type())           \
      .def(Type() + py::self)           \
      .def(Type() - py::self)           \
      .def(Type() * py::self)           \
      .def(Type() / py::self)           \
      .def(Type() == py::self)          \
      .def(Type() != py::self)          \
      .def(Type() <= py::self)          \
      .def(Type() < py::self)           \
      .def(Type() >= py::self)          \
      .def(Type() > py::self)

#define PythonInterfaceCommonMathSelf \
  def(+py::self)                      \
      .def(-py::self)                 \
      .def(py::self + py::self)       \
      .def(py::self - py::self)       \
      .def(py::self* py::self)        \
      .def(py::self / py::self)       \
      .def(py::self += py::self)      \
      .def(py::self -= py::self)      \
      .def(py::self *= py::self)      \
      .def(py::self /= py::self)      \
      .def(py::self == py::self)      \
      .def(py::self != py::self)      \
      .def(py::self <= py::self)      \
      .def(py::self < py::self)       \
      .def(py::self >= py::self)      \
      .def(py::self > py::self)

#define PythonInterfaceWorkspaceVariableConversion(Type)                 \
  def(py::init([](const Type& x) { return std::make_shared<Type>(x); }), \
      py::arg("val"),                                                    \
      py::doc("Copy instance"))

//! Place at the end!
#define PythonInterfaceWorkspaceDocumentation(Type)                           \
  doc() = unwrap_stars(var_string(internal_workspace_groups().at(#Type).desc, \
                                  '\n',                                       \
                                  group_generics_inout(#Type),                \
                                  group_workspace_types(#Type)))              \
              .c_str()

//! Place at the end!  Add a few line-breaks to fit in extra documentation!
#define PythonInterfaceWorkspaceDocumentationExtra(Type, extra)               \
  doc() = unwrap_stars(var_string(internal_workspace_groups().at(#Type).desc, \
                                  extra,                                      \
                                  '\n',                                       \
                                  group_generics_inout(#Type),                \
                                  group_workspace_types(#Type)))              \
              .c_str()

#define PythonInterfaceReadWriteData(Type, data, docstr)     \
  def_readwrite(#data,                                       \
                &Type::data,                                 \
                py::return_value_policy::reference_internal, \
                py::doc(docstr))

#define PythonInterfaceBasicReferenceProperty(                               \
    Type, PropertyName, ReadFunction, WriteFunction, docstr)                 \
  def_property(                                                              \
      #PropertyName,                                                         \
      py::cpp_function(                                                      \
          [](Type& x)                                                        \
              -> std::remove_cv_t<                                           \
                  std::add_lvalue_reference_t<decltype(x.ReadFunction())>> { \
            return x.ReadFunction();                                         \
          },                                                                 \
          py::return_value_policy::reference_internal),                      \
      [](Type& x, std::remove_reference_t<decltype(x.ReadFunction())> y) {   \
        x.WriteFunction() = std::move(y);                                    \
      },                                                                     \
      py::doc(docstr))

#define PythonInterfaceSelfAttribute(ATTR)                       \
  def_property_readonly(                                         \
      #ATTR,                                                     \
      [](py::object& x) { return x.attr("value").attr(#ATTR); }, \
      "As for :class:`numpy.ndarray`")

#define PythonInterfaceSelfOperator(ATTR)                          \
  def(                                                             \
      #ATTR,                                                       \
      [](py::object& x) { return x.attr("value").attr(#ATTR)(); }, \
      "As for :class:`numpy.ndarray`")

#define PythonInterfaceValueOperator(ATTR)     \
  def(                                         \
      #ATTR,                                   \
      [](py::object& x, py::object& y) {       \
        return x.attr("value").attr(#ATTR)(y); \
      },                                       \
      "As for :class:`numpy.ndarray`")

#define PythonInterfaceTwoValueOperator(ATTR)           \
  def(                                                  \
      #ATTR,                                            \
      [](py::object& x, py::object& y, py::object& z) { \
        return x.attr("value").attr(#ATTR)(y, z);       \
      },                                                \
      "As for :class:`numpy.ndarray`")

#define PythonInterfaceInternalValueOperator(ATTR) \
  def(                                             \
      #ATTR,                                       \
      [](py::object& x, py::object& y) {           \
        x.attr("value").attr(#ATTR)(y);            \
        return x;                                  \
      },                                           \
      "As for :class:`numpy.ndarray`")

#define PythonInterfaceInternalTwoValueOperator(ATTR)   \
  def(                                                  \
      #ATTR,                                            \
      [](py::object& x, py::object& y, py::object& z) { \
        x.attr("value").attr(#ATTR)(y, z);              \
        return x;                                       \
      },                                                \
      "As for :class:`numpy.ndarray`")

#define PythonInterfaceNumpyValueProperties \
  PythonInterfaceSelfAttribute(ndim)        \
      .PythonInterfaceSelfAttribute(shape)  \
      .PythonInterfaceSelfAttribute(size)   \
      .PythonInterfaceSelfAttribute(dtype)

#define PythonInterfaceValueOperators                      \
  PythonInterfaceSelfOperator(__len__)                     \
      .PythonInterfaceSelfOperator(__iter__)               \
      .PythonInterfaceSelfOperator(__pos__)                \
      .PythonInterfaceSelfOperator(__neg__)                \
      .PythonInterfaceSelfOperator(__abs__)                \
      .PythonInterfaceSelfOperator(__invert__)             \
      .PythonInterfaceSelfOperator(__complex__)            \
      .PythonInterfaceSelfOperator(__int__)                \
      .PythonInterfaceSelfOperator(__float__)              \
                                                           \
      .PythonInterfaceValueOperator(__add__)               \
      .PythonInterfaceValueOperator(__sub__)               \
      .PythonInterfaceValueOperator(__mul__)               \
      .PythonInterfaceValueOperator(__matmul__)            \
      .PythonInterfaceValueOperator(__truediv__)           \
      .PythonInterfaceValueOperator(__floordiv__)          \
      .PythonInterfaceValueOperator(__mod__)               \
      .PythonInterfaceValueOperator(__divmod__)            \
      .PythonInterfaceValueOperator(__pow__)               \
      .PythonInterfaceTwoValueOperator(__pow__)            \
                                                           \
      .PythonInterfaceValueOperator(__radd__)              \
      .PythonInterfaceValueOperator(__rsub__)              \
      .PythonInterfaceValueOperator(__rmul__)              \
      .PythonInterfaceValueOperator(__rmatmul__)           \
      .PythonInterfaceValueOperator(__rtruediv__)          \
      .PythonInterfaceValueOperator(__rfloordiv__)         \
      .PythonInterfaceValueOperator(__rmod__)              \
      .PythonInterfaceValueOperator(__rdivmod__)           \
      .PythonInterfaceValueOperator(__rpow__)              \
      .PythonInterfaceTwoValueOperator(__rpow__)           \
                                                           \
      .PythonInterfaceInternalValueOperator(__iadd__)      \
      .PythonInterfaceInternalValueOperator(__isub__)      \
      .PythonInterfaceInternalValueOperator(__imul__)      \
      .PythonInterfaceInternalValueOperator(__imatmul__)   \
      .PythonInterfaceInternalValueOperator(__itruediv__)  \
      .PythonInterfaceInternalValueOperator(__ifloordiv__) \
      .PythonInterfaceInternalValueOperator(__imod__)      \
      .PythonInterfaceInternalValueOperator(__idivmod__)   \
      .PythonInterfaceInternalValueOperator(__ipow__)      \
      .PythonInterfaceInternalTwoValueOperator(__ipow__)   \
                                                           \
      .PythonInterfaceValueOperator(__lt__)                \
      .PythonInterfaceValueOperator(__le__)                \
      .PythonInterfaceValueOperator(__eq__)                \
      .PythonInterfaceValueOperator(__ne__)                \
      .PythonInterfaceValueOperator(__ge__)                \
      .PythonInterfaceValueOperator(__gt__)                \
                                                           \
      .PythonInterfaceValueOperator(__contains__)          \
                                                           \
      .PythonInterfaceValueOperator(__getitem__)           \
      .PythonInterfaceTwoValueOperator(__setitem__)        \
                                                           \
      .PythonInterfaceSelfOperator(flatten)

#endif
