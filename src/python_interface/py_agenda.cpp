#include <auto_md.h>
#include <xml_io.h>
#include <global_data.h>
#include <algorithm>

#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_agenda(py::module_& m) {
  py::class_<MdRecord>(m, "MdRecord")
      .def(py::init<>())
      .def_property_readonly("name", &MdRecord::Name);

  py::class_<Array<MdRecord>>(m, "ArrayOfMdRecord")
      .PythonInterfaceArrayDefault(MdRecord);

  py::class_<Agenda>(m, "Agenda")
      .def(py::init<>())
      .PythonInterfaceBasicRepresentation(Agenda);
}

void py_agenda_methods(py::module_& m) {
  m.def("workspace_methods_list", [](){return global_data::md_data;});

  m.def("workspace_methods_map", []() {
    return global_data::MdMap;
  });

  m.def("workspace_methods_raw_list", [](){return global_data::md_data_raw;});

  m.def("workspace_methods_raw_map", []() {
    return global_data::MdRawMap;
  });
}
}  // namespace Python
