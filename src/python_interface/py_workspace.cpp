#include <auto_md.h>
#include <global_data.h>
#include <pybind11/cast.h>
#include <pybind11/pybind11.h>
#include <xml_io.h>

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

  auto ws = py::class_<Workspace>(m, "Workspace").def(py::init([]() {
    Workspace w{};
    w.initialize();
    return w;
  }));

  ws.def_property_readonly("size", &Workspace::nelem);

  py_auto_workspace(ws);
}
}  // namespace Python