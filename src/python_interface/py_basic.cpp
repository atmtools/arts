#include <iomanip>

#include <python_interface.h>

#include "py_macros.h"

namespace Python {
void py_basic(py::module_& m) try {
  py::class_<String>(m, "String")
      .def(py::init([]() { return std::make_unique<String>(); }), "Create empty")
      .def(py::init([](const std::string& s) { return std::make_unique<String>(s); }), "Create from :class:`str`")
      .PythonInterfaceCopyValue(String)
      .PythonInterfaceWorkspaceVariableConversion(String)
      .PythonInterfaceFileIO(String)
      .PythonInterfaceIndexItemAccess(String)
      .def("__repr__", [](const String& s){return var_string(std::quoted(s));}, py::is_operator())
      .def("__str__", [](const String& s){return var_string(s);}, py::is_operator())
      .def(py::self + py::self)
      .def(py::self += py::self)
      .def(py::self == py::self)
      .def("__hash__", [](String& x) { return py::hash(py::str(x)); })
      .def(py::pickle(
          [](const String& self) { return py::make_tuple(std::string(self)); },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_unique<String>(t[0].cast<std::string>());
          }))
      .PythonInterfaceWorkspaceDocumentationExtra(String,
                                                  R"--(

This class is compatible with python strings.  The data can 
be accessed without copy using element-wise access operators)--");

  PythonInterfaceWorkspaceArray(String);
  PythonInterfaceWorkspaceArray(ArrayOfString);

  py::class_<ArrayOfIndex>(m, "ArrayOfIndex", py::buffer_protocol())
      .PythonInterfaceFileIO(ArrayOfIndex)
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfIndex)
      .PythonInterfaceBasicRepresentation(ArrayOfIndex)
      .PythonInterfaceArrayDefault(Index)
      .PythonInterfaceValueOperators
      .PythonInterfaceNumpyValueProperties
      .def_buffer([](ArrayOfIndex& x) -> py::buffer_info {
        return py::buffer_info(x.data(),
                               sizeof(Index),
                               py::format_descriptor<Index>::format(),
                               1,
                               {x.nelem()},
                               {sizeof(Index)});
      })
      .def_property("value",
                    py::cpp_function(
                        [](ArrayOfIndex& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>()),
                    [](ArrayOfIndex& x, ArrayOfIndex& y) { x = y; }, "Operate on type as if :class:`numpy.ndarray` type")
      .PythonInterfaceWorkspaceDocumentation(ArrayOfIndex);
  py::implicitly_convertible<std::vector<Index>, ArrayOfIndex>();

  PythonInterfaceWorkspaceArray(ArrayOfIndex);
  py::implicitly_convertible<std::vector<std::vector<Index>>,
                             ArrayOfArrayOfIndex>();

  py::class_<Numeric_>(m, "Numeric")
      .def(py::init([]() { return std::make_unique<Numeric_>(); }), "Create default")
      .def(py::init([](Index i) -> Numeric_ {
        return Numeric_{static_cast<Numeric>(i)};
      }), "Create from :class:`int`")
      .def(py::init([](Numeric n) { return std::make_unique<Numeric_>(n); }), "Create from :class:`float`")
      .PythonInterfaceCopyValue(Numeric_)
      .PythonInterfaceWorkspaceVariableConversion(Numeric_)
      .def("__hash__", [](Numeric_& x) { return py::hash(py::float_(x.val)); })
      .PythonInterfaceCommonMath(Numeric)
      .PythonInterfaceCommonMath(Index_)
      .PythonInterfaceCommonMath(Index)
      .PythonInterfaceCommonMathSelf
      .def_property(
          "value",
          [](Numeric_& x) { return x.val; },
          [](Numeric_& x, Numeric_ y) { x.val = y.val; }, "Value of instance as :class:`float`")
      .def("__repr__", [](const Numeric_& x){return var_string(x);}, py::is_operator())
      .def("__str__", [](const Numeric_& x){return var_string(x);}, py::is_operator())
      .def(
          "savexml",
          [](const Numeric_& x,
             const char* const file,
             const char* const type,
             bool clobber) {
            xml_write_to_file(file,
                              x.val,
                              string2filetype(type),
                              clobber ? 0 : 1);
          },
          py::arg("file").none(false),
          py::arg("type").none(false) = "ascii",
          py::arg("clobber") = true,
          py::doc("Saves Numeric to file\n"
                  "\n"
                  "Parameters:\n"
                  "    file (str): The path to which the file is written."
                  " Note that several of the options might modify the"
                  " name or write more files\n"
                  "    type (str): Type of file to save (ascii. zascii,"
                  " or binary)\n"
                  "    clobber (bool): Overwrite existing files or add new"
                  " file with modified name?\n"
                  "\n"
                  "On Error:\n"
                  "    Throws RuntimeError for any failure to save"))
      .def(
          "readxml",
          [](Numeric_& x, const char* const file) {
            xml_read_from_file(file, x.val);
          },
          py::arg("file").none(false),
          py::doc("Read Numeric from file\n"
                  "\n"
                  "Parameters:\n"
                  "    file (str): A file that can be read\n"
                  "\n"
                  "On Error:\n"
                  "    Throws RuntimeError for any failure to read"))
      .def(py::pickle(
          [](const Numeric_& self) { return py::make_tuple(self.val); },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_unique<Numeric_>(t[0].cast<Numeric>());
          }))
      .PythonInterfaceWorkspaceDocumentation(Numeric);

  py::class_<Index_>(m, "Index")
      .def(py::init([]() { return std::make_unique<Index_>(); }), "Create default")
      .def(py::init([](Index i) { return std::make_unique<Index_>(i); }), "Create from :class:`int`")
      .PythonInterfaceCopyValue(Index_)
      .PythonInterfaceWorkspaceVariableConversion(Index_)
      .def("__hash__", [](Index_& x) { return py::hash(py::int_(x.val)); })
      .PythonInterfaceCommonMath(Numeric)
      .PythonInterfaceCommonMath(Index)
      .PythonInterfaceCommonMathSelf
      .def_property(
          "value",
          [](Index_& x) { return x.val; },
          [](Index_& x, Index_ y) { x.val = y.val; }, "Value of instance as :class:`int`")
      .def("__repr__", [](const Index_& x){return var_string(x);}, py::is_operator())
      .def("__str__", [](const Index_& x){return var_string(x);}, py::is_operator())
      .def(
          "savexml",
          [](const Index_& x,
             const char* const file,
             const char* const type,
             bool clobber) {
            xml_write_to_file(file,
                              x.val,
                              string2filetype(type),
                              clobber ? 0 : 1);
          },
          py::arg("file").none(false),
          py::arg("type").none(false) = "ascii",
          py::arg("clobber") = true,
          py::doc("Saves Index to file\n"
                  "\n"
                  "Parameters:\n"
                  "    file (str): The path to which the file is written."
                  " Note that several of the options might modify the"
                  " name or write more files\n"
                  "    type (str): Type of file to save (ascii. zascii,"
                  " or binary)\n"
                  "    clobber (bool): Overwrite existing files or add new"
                  " file with modified name?\n"
                  "\n"
                  "On Error:\n"
                  "    Throws RuntimeError for any failure to save"))
      .def(
          "readxml",
          [](Index_& x, const char* const file) {
            xml_read_from_file(file, x.val);
          },
          py::arg("file").none(false),
          py::doc("Read Index from file\n"
                  "\n"
                  "Parameters:\n"
                  "    file (str): A file that can be read\n"
                  "\n"
                  "On Error:\n"
                  "    Throws RuntimeError for any failure to read"))
      .def(py::pickle(
          [](const Index_& self) { return py::make_tuple(self.val); },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_unique<Index_>(t[0].cast<Index>());
          }))
      .PythonInterfaceWorkspaceDocumentation(Index);

  py::implicitly_convertible<py::str, String>();
  py::implicitly_convertible<std::vector<py::str>, ArrayOfString>();
  py::implicitly_convertible<Index, Numeric_>();
  py::implicitly_convertible<Numeric, Numeric_>();
  py::implicitly_convertible<Index, Index_>();

  py::class_<Any>(m, "Any")
      .def(py::init([]() { return std::make_unique<Any>(); }), "Create empty")
      .def(py::init([](const py::args&, const py::kwargs&) { return std::make_unique<Any>(); }), "Create empty")
      .def("__repr__", [](Any&) { return "Any"; }, py::is_operator())
      .def("__str__", [](Any&) { return "Any"; }, py::is_operator())
      .def(py::pickle([](const py::object&) { return py::make_tuple(); },
                      [](const py::tuple& t) {
                        ARTS_USER_ERROR_IF(t.size() != 0, "Invalid state!")
                        return std::make_unique<Any>();
                      }))
      .PythonInterfaceWorkspaceDocumentation(Any);

  py::class_<ArrayOfNumeric>(m, "ArrayOfNumeric", py::buffer_protocol())
      // .PythonInterfaceFileIO(ArrayOfNumeric)
      .PythonInterfaceBasicRepresentation(ArrayOfNumeric)
      .PythonInterfaceArrayDefault(Numeric)
      .PythonInterfaceValueOperators.PythonInterfaceNumpyValueProperties
      .def_buffer([](ArrayOfNumeric& x) -> py::buffer_info {
        return py::buffer_info(x.data(),
                               sizeof(Numeric),
                               py::format_descriptor<Numeric>::format(),
                               1,
                               {x.nelem()},
                               {sizeof(Numeric)});
      })
      .def_property("value",
                    py::cpp_function(
                        [](ArrayOfNumeric& x) {
                          py::object np = py::module_::import("numpy");
                          return np.attr("array")(x, py::arg("copy") = false);
                        },
                        py::keep_alive<0, 1>(), "Value of instance as :class:`numpy.ndarray`"),
                    [](ArrayOfNumeric& x, ArrayOfNumeric& y) { x = y; })
      .doc() = R"--(A list of :class:`~pyarts.arts.Numeric`)--";
  py::implicitly_convertible<std::vector<Numeric>, ArrayOfNumeric>();
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize basic\n", e.what()));
}
}  // namespace Python
