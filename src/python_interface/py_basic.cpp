#include <auto_md.h>
#include <pybind11/cast.h>
#include <pybind11/pybind11.h>
#include <xml_io.h>

#include <initializer_list>

#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_basic(py::module_& m) {
  py::class_<Verbosity>(m, "Verbosity")
      .def(py::init<>())
      .def(py::init<Index, Index, Index>(),
      py::arg("agenda")=0,
      py::arg("screen")=0,
      py::arg("file")=0)
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
      .def(py::init<char*>())
      .PythonInterfaceFileIO(String)
      .PythonInterfaceIndexItemAccess(String)
      .PythonInterfaceBasicIteration(String)
      .PythonInterfaceBasicRepresentation(String)
      .def(py::self + py::self)
      .def(py::self += py::self)
      .doc() = "The Arts String class\n"
      "\n"
      "This class is compatible with python strings.  The data can "
      "be accessed without copy using element-wise access operators.";
  py::implicitly_convertible<py::str, String>();

  PythonInterfaceWorkspaceArray(String);
  PythonInterfaceWorkspaceArray(ArrayOfString);

  py::class_<ArrayOfIndex>(m, "ArrayOfIndex", py::buffer_protocol())
      .PythonInterfaceFileIO(ArrayOfIndex)
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
  py::implicitly_convertible<py::array, ArrayOfIndex>();
  py::implicitly_convertible<py::list, ArrayOfIndex>();

  py::implicitly_convertible<std::vector<Index>, ArrayOfIndex>();  
  PythonInterfaceWorkspaceArray(ArrayOfIndex);      
}
}  // namespace Python