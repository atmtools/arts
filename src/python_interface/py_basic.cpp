#include <pybind11/pybind11.h>
#include <python_interface.h>

#include <iomanip>

#include "mystring.h"
#include "py_macros.h"
#include "python_interface_value_type.h"

#define PythonInterfaceValueHolder(Type)                                     \
  artsclass<ValueHolder<Type>>(m, #Type)                                     \
      .def(py::init([]() -> ValueHolder<Type> {                              \
             return std::make_shared<Type>();                                \
           }),                                                               \
           "Create default")                                                 \
      .def(py::init([](const Type& a) -> ValueHolder<Type> {                 \
             return std::make_shared<Type>(a);                               \
           }),                                                               \
           "Create from value")                                              \
      .def_property(                                                         \
          "value",                                                           \
          [](const ValueHolder<Type>& x) { return *x.val; },                 \
          [](ValueHolder<Type>& x, const ValueHolder<Type>& y) {             \
            *x.val = *y.val;                                                 \
          })                                                                 \
      .def("__str__",                                                        \
           [](py::object& x) { return x.attr("value").attr("__str__")(); })  \
      .def("__repr__",                                                       \
           [](py::object& x) { return x.attr("value").attr("__repr__")(); }) \
      .def("__add__",                                                        \
           [](const py::object& x, const py::object& y) {                    \
             return x.attr("value").attr("__add__")(y);                      \
           })                                                                \
      .def("__sub__",                                                        \
           [](const py::object& x, const py::object& y) {                    \
             return x.attr("value").attr("__sub__")(y);                      \
           })                                                                \
      .def("__mul__",                                                        \
           [](const py::object& x, const py::object& y) {                    \
             return x.attr("value").attr("__mul__")(y);                      \
           })                                                                \
      .def("__div__",                                                        \
           [](const py::object& x, const py::object& y) {                    \
             return x.attr("value").attr("__div__")(y);                      \
           })                                                                \
      .def("__truediv__",                                                    \
           [](const py::object& x, const py::object& y) {                    \
             return x.attr("value").attr("__truediv__")(y);                  \
           })                                                                \
      .PythonInterfaceWorkspaceDocumentation(Type);                          \
  py::implicitly_convertible<Type, ValueHolder<Type>>();

namespace Python {
void py_basic(py::module_& m) try {
  PythonInterfaceValueHolder(String);
  PythonInterfaceValueHolder(Numeric);
  PythonInterfaceValueHolder(Index);

  //PythonInterfaceWorkspaceArray(String);
  artsclass<ArrayOfString>(m, "ArrayOfString");
  artsclass<ArrayOfIndex>(m, "ArrayOfIndex");
  py::implicitly_convertible<std::vector<String>, ArrayOfString>();
  py::implicitly_convertible<std::vector<Index>, ArrayOfIndex>();

  PythonInterfaceWorkspaceArray(ArrayOfString);
  PythonInterfaceWorkspaceArray(ArrayOfIndex);

  artsclass<Any>(m, "Any")
      .def(py::init([]() { return std::make_shared<Any>(); }), "Create empty")
      .def(py::init([](const py::args&, const py::kwargs&) {
             return std::make_shared<Any>();
           }),
           "Create empty")
      .PythonInterfaceBasicRepresentation(Any)
      .PythonInterfaceFileIO(Any)
      .def(py::pickle([](const py::object&) { return py::make_tuple(); },
                      [](const py::tuple& t) {
                        ARTS_USER_ERROR_IF(t.size() != 0, "Invalid state!")
                        return std::make_shared<Any>();
                      }))
      .PythonInterfaceWorkspaceDocumentation(Any);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize basic\n", e.what()));
}
}  // namespace Python
