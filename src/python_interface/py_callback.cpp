#include <callback.h>
#include <pybind11/attr.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "py_macros.h"

namespace Python {
void py_callback(py::module_& m) {
  py::class_<CallbackOperator>(m, "CallbackOperator")
      .def(py::init<>())
      .def(py::init<std::function<void(std::shared_ptr<Workspace>)>,
                    std::vector<std::string>,
                    std::vector<std::string>>(),
           py::arg("fn"),
           py::arg("inputs"),
           py::arg("outputs"))
      .PythonInterfaceCopyValue(CallbackOperator)
      .PythonInterfaceBasicRepresentation(CallbackOperator)
      .PythonInterfaceFileIO(CallbackOperator)
      .def(
          "__call__",
          [](const CallbackOperator& cb, Workspace& w) { cb(w); },
          py::arg("ws"),
          py::is_operator())
      .PythonInterfaceWorkspaceDocumentation(CallbackOperator);
}
}  // namespace Python
