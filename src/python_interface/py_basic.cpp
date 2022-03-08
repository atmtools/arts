#include <py_auto_interface.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <vector>

#include "py_macros.h"

namespace Python {
void py_basic(py::module_& m) {
  py::class_<Verbosity>(m, "Verbosity")
      .def(py::init<>())
      .PythonInterfaceWorkspaceVariableConversion(Verbosity)
      .def(py::init<Index, Index, Index>(),
           py::arg("agenda") = 0,
           py::arg("screen") = 0,
           py::arg("file") = 0)
      .PythonInterfaceFileIO(Verbosity)
      .PythonInterfaceBasicRepresentation(Verbosity)
      .def_property(
          "agenda",
          [](Verbosity& x) { return x.get_agenda_verbosity(); },
          [](Verbosity& x, Index i) { return x.set_agenda_verbosity(i); },
          py::doc("Agenda output level"))
      .def_property(
          "screen",
          [](Verbosity& x) { return x.get_screen_verbosity(); },
          [](Verbosity& x, Index i) { return x.set_screen_verbosity(i); },
          py::doc("Screen output level"))
      .def_property(
          "file",
          [](Verbosity& x) { return x.get_file_verbosity(); },
          [](Verbosity& x, Index i) { return x.set_file_verbosity(i); },
          py::doc("File output level"))
      .def_property(
          "main",
          [](Verbosity& x) { return x.is_main_agenda(); },
          [](Verbosity& x, bool i) { return x.set_main_agenda(i); },
          py::doc("Main class flag"));

  py::class_<String>(m, "String")
      .def(py::init<>())
      .def(py::init<const char*>(), py::arg("str").none(false))
      .PythonInterfaceWorkspaceVariableConversion(String)
      .PythonInterfaceFileIO(String)
      .PythonInterfaceIndexItemAccess(String)
      .PythonInterfaceBasicRepresentation(String)
      .def(py::self + py::self)
      .def(py::self += py::self)
      .def(py::self == py::self)
      .def("__hash__", [](String& x){return py::hash(py::str(x));})
      .doc() =
      R"--(
The Arts String class

This class is compatible with python strings.  The data can 
be accessed without copy using element-wise access operators.
)--";

  PythonInterfaceWorkspaceArray(String).def(
      "index", [](ArrayOfString& a, String& b) {
        auto ptr = std::find(a.begin(), a.end(), b);
        if (ptr == a.end())
          throw std::invalid_argument(var_string(b, " not in list"));
        return std::distance(a.begin(), ptr);
      });
  PythonInterfaceWorkspaceArray(ArrayOfString);

  py::class_<ArrayOfIndex>(m, "ArrayOfIndex", py::buffer_protocol())
      .PythonInterfaceFileIO(ArrayOfIndex)
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfIndex)
      .PythonInterfaceBasicRepresentation(ArrayOfIndex)
      .PythonInterfaceArrayDefault(Index)
      .def_buffer([](ArrayOfIndex& x) -> py::buffer_info {
        return py::buffer_info(x.data(),
                               sizeof(Index),
                               py::format_descriptor<Index>::format(),
                               1,
                               {x.nelem()},
                               {sizeof(Index)});
      })
      .doc() =
      "The Arts ArrayOfIndex class\n"
      "\n"
      "This class is compatible with numpy arrays.  The data can "
      "be accessed without copy using np.array(x, copy=False), "
      "with x as an instance of this class.  Note that access to "
      "any and all of the mathematical operations are only available "
      "via the numpy interface.  The main constern equality operations, "
      "which only checks pointer equality if both LHS and RHS of this "
      "object's instances (i.e., no element-wise comparisions)\n"
      "\n"
      "Initialization:\n"
      "    ArrayOfIndex(): for basic initialization\n\n"
      "    ArrayOfIndex(Index): for constant size, unknown value\n\n"
      "    ArrayOfIndex(Index, Index): for constant size, constant value\n\n"
      "    ArrayOfIndex(List or Array): to copy elements\n\n";
  py::implicitly_convertible<std::vector<Index>, ArrayOfIndex>();

  PythonInterfaceWorkspaceArray(ArrayOfIndex);
  py::implicitly_convertible<std::vector<std::vector<Index>>, ArrayOfArrayOfIndex>();

  py::class_<Numeric_>(m, "Numeric")
      .def(py::init<>())
      .def(py::init<Numeric>())
      .PythonInterfaceWorkspaceVariableConversion(Numeric_)
      .def(
          "__ne__",
          [](Numeric_& a, Numeric_ b) { return a.val != b.val; },
          py::is_operator())
      .def(
          "__eq__",
          [](Numeric_& a, Numeric_ b) { return a.val == b.val; },
          py::is_operator())
      .def(
          "__le__",
          [](Numeric_& a, Numeric_ b) { return a.val <= b.val; },
          py::is_operator())
      .def(
          "__ge__",
          [](Numeric_& a, Numeric_ b) { return a.val >= b.val; },
          py::is_operator())
      .def(
          "__lt__",
          [](Numeric_& a, Numeric_ b) { return a.val < b.val; },
          py::is_operator())
      .def(
          "__gt__",
          [](Numeric_& a, Numeric_ b) { return a.val > b.val; },
          py::is_operator())
      .def(
          "__mul__",
          [](Numeric_& a, Numeric_ b) { return a.val * b.val; },
          py::is_operator())
      .def(
          "__add__",
          [](Numeric_& a, Numeric_ b) { return a.val + b.val; },
          py::is_operator())
      .def(
          "__truediv__",
          [](Numeric_& a, Numeric_ b) { return a.val / b.val; },
          py::is_operator())
      .def(
          "__sub__",
          [](Numeric_& a, Numeric_ b) { return a.val - b.val; },
          py::is_operator())
      .def(
          "__pow__",
          [](Numeric_& a, Numeric_ b) { return std::pow(a.val, b.val); },
          py::is_operator())
      .def_property(
          "val",
          [](Numeric_& x) { return x.val; },
          [](Numeric_& x, Numeric y) { x.val = y; })
      .PythonInterfaceBasicRepresentation(Numeric_)
      .def(
          "savexml",
          [](const Numeric_& x,
             const char* const file,
             const char* const type,
             bool clobber) {
            xml_write_to_file(file,
                              x.val,
                              string2filetype(type),
                              clobber ? 0 : 1,
                              Verbosity());
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
            xml_read_from_file(file, x.val, Verbosity());
          },
          py::arg("file").none(false),
          py::doc("Read Numeric from file\n"
                  "\n"
                  "Parameters:\n"
                  "    file (str): A file that can be read\n"
                  "\n"
                  "On Error:\n"
                  "    Throws RuntimeError for any failure to read"))
      .doc() =
      R"--(
This is a wrapper class for Arts Numeric.

You can get copies and set the value by the \"val\" property
)--";

  py::class_<Index_>(m, "Index")
      .def(py::init<>())
      .def(py::init<Index>())
      .def(py::init<Numeric>())
      .PythonInterfaceWorkspaceVariableConversion(Index_)
      .def(
          "__ne__",
          [](Index_& a, Index_ b) { return a.val != b.val; },
          py::is_operator())
      .def(
          "__eq__",
          [](Index_& a, Index_ b) { return a.val == b.val; },
          py::is_operator())
      .def(
          "__le__",
          [](Index_& a, Index_ b) { return a.val <= b.val; },
          py::is_operator())
      .def(
          "__ge__",
          [](Index_& a, Index_ b) { return a.val >= b.val; },
          py::is_operator())
      .def(
          "__lt__",
          [](Index_& a, Index_ b) { return a.val < b.val; },
          py::is_operator())
      .def(
          "__gt__",
          [](Index_& a, Index_ b) { return a.val > b.val; },
          py::is_operator())
      .def(
          "__mul__",
          [](Index_& a, Index_ b) { return a.val * b.val; },
          py::is_operator())
      .def(
          "__add__",
          [](Index_& a, Index_ b) { return a.val + b.val; },
          py::is_operator())
      .def(
          "__div__",
          [](Index_& a, Index_ b) { return a.val / b.val; },
          py::is_operator())
      .def(
          "__truediv__",
          [](Index_& a, Index_ b) { return Numeric(Rational(a.val, b.val)); },
          py::is_operator())
      .def(
          "__sub__",
          [](Index_& a, Index_ b) { return a.val - b.val; },
          py::is_operator())
      .def(
          "__pow__",
          [](Index_& a, Index_ b) { return std::pow(a.val, b.val); },
          py::is_operator())
      .def_property(
          "val",
          [](Index_& x) { return x.val; },
          [](Index_& x, Index y) { x.val = y; })
      .PythonInterfaceBasicRepresentation(Index_)
      .def(
          "savexml",
          [](const Index_& x,
             const char* const file,
             const char* const type,
             bool clobber) {
            xml_write_to_file(file,
                              x.val,
                              string2filetype(type),
                              clobber ? 0 : 1,
                              Verbosity());
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
            xml_read_from_file(file, x.val, Verbosity());
          },
          py::arg("file").none(false),
          py::doc("Read Index from file\n"
                  "\n"
                  "Parameters:\n"
                  "    file (str): A file that can be read\n"
                  "\n"
                  "On Error:\n"
                  "    Throws RuntimeError for any failure to read"))
      .doc() =
      R"--(
This is a wrapper class for Arts Index.

You can get copies and set the value by the \"val\" property
)--";

  py::implicitly_convertible<py::str, String>();
  py::implicitly_convertible<std::vector<py::str>, ArrayOfString>();
  py::implicitly_convertible<Index, Numeric_>();
  py::implicitly_convertible<Numeric, Numeric_>();
  py::implicitly_convertible<Index, Index_>();
  py::implicitly_convertible<Numeric, Index_>();

  py::class_<Any>(m, "Any")
      .def("__repr__", [](Any&) { return "Any"; })
      .def("__str__", [](Any&) { return "Any"; });
}
}  // namespace Python
