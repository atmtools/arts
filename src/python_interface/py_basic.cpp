#include "python_interface.h"

#include "py_macros.h"
#include <auto_md.h>
#include <xml_io.h>

namespace Python {
void py_basic(py::module_& m) {
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
  py::implicitly_convertible<const char* const, String>();
}
}  // namespace Python