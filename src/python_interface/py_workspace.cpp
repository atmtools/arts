#include <auto_md.h>
#include <global_data.h>
#include <pybind11/attr.h>
#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/gil.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <xml_io.h>
#include <supergeneric.h>
#include <iterator>

#include "debug.h"
#include "jacobian.h"
#include "py_macros.h"
#include "python_interface.h"
#include "workspace_ng.h"

namespace Python {

void py_auto_workspace(py::class_<Workspace>&);
void py_workspace(py::module_& m) {
  define_wsv_group_names();
  Workspace::define_wsv_data();
  Workspace::define_wsv_map();
  define_md_data_raw();
  expand_md_data_raw_to_md_data();
  define_md_map();
  define_agenda_data();
  define_agenda_map();
  global_data::workspace_memory_handler.initialize();

  auto ws = py::class_<Workspace>(m, "Workspace")
                .def(py::init([]() {
                  Workspace w{};
                  w.initialize();
                  return w;
                }));

  ws.def_property_readonly("size", &Workspace::nelem);

  ws.def(
      "create_variable",
      [](Workspace& w, char* group, char* name, char* desc) {
        using global_data::workspace_memory_handler;

        ARTS_USER_ERROR_IF(std::find_if(Workspace::WsvMap.begin(),
                                        Workspace::WsvMap.end(),
                                        [name](auto& b) {
                                          return b.first == name;
                                        }) not_eq Workspace::WsvMap.end(),
                           "A variable of this name already exist: ",
                           name)

        const Index group_index = get_wsv_group_id(group);
        ARTS_USER_ERROR_IF(group_index < 0, "Cannot recognize group: ", group)

        WsvRecord x(name, desc ? desc : "Created by pybind11 API", group_index);

        Index out = w.add_wsv_inplace(x);
        w.push(out, workspace_memory_handler.allocate(group_index));
      },
      py::arg("group"),
      py::arg("name"),
      py::arg("desc") = nullptr,
      py::doc(
          R"--(
Creates a named variable of the given group, initializing it

The variable can be accessed with its property name upon completion
)--"));

  ws.def("__delattr__", [](Workspace& w, const char* name) {
    auto varpos = std::find_if(Workspace::WsvMap.begin(),
                               Workspace::WsvMap.end(),
                               [name](auto& b) { return b.first == name; });
    ARTS_USER_ERROR_IF(varpos == Workspace::WsvMap.end(),
                       "No workspace variable: ",
                       name,
                       "\n\nCustom workspace variables have to "
                       "be created using create_variable or by explicit set")

    if (w.is_initialized(varpos->second)) {
      w.pop_free(varpos->second);
    }
  });

  py_auto_workspace(ws);
}
}  // namespace Python