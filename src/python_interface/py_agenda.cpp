#include <auto_md.h>
#include <xml_io.h>
#include <global_data.h>
#include <algorithm>

#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_agenda(py::module_& m) {
  py::class_<CallbackFunction>(m, "CallbackFunction")
      .def(py::init<>())
      .def(py::init<std::function<void(void)>>())
      .def("__call__", [](CallbackFunction& f){f();});
  py::implicitly_convertible<std::function<void(void)>, CallbackFunction>();

  py::class_<MdRecord>(m, "MdRecord")
      .def(py::init<>())
      .def_property_readonly("name", &MdRecord::Name);

  py::class_<Array<MdRecord>>(m, "ArrayOfMdRecord")
      .PythonInterfaceArrayDefault(MdRecord);

  py::class_<Agenda>(m, "Agenda")
      .def(py::init<>())
      .def("add_method",
           [](Agenda& a, String method_name, py::args& args, py::kwargs& kwargs) {
             const auto nargs = args.size();
             const auto nkwargs = kwargs.size();
             //py::print(method_name);
             //auto method_index = global_data::md_data_raw.at(method_name);
             //auto method_data = global_data::md_data_raw.at(method_index);
             //py::print(method_data.Name());
           })
      .def("__repr__", [](Agenda& a) {
        std::string out = var_string("Arts Agenda ", a.name(), ":\n");
        std::ostringstream os;
        a.print(os, "   ");
        out += os.str();
        return out;
      })
      .def("__str__", [](Agenda& a) {
        std::string out = var_string("Arts Agenda ", a.name(), ":\n");
        std::ostringstream os;
        a.print(os, "   ");
        out += os.str();
        return out;
      });
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
