#include <callback.h>
#include <pybind11/pybind11.h>

#include "py_macros.h"

namespace Python {
void py_callback(py::module_& m) {
  py::class_<CallbackOperator>(m, "CallbackOperator")
      .def(py::init<>())
      .def(py::init<std::function<void(Workspace&)>,
                    std::vector<std::string>,
                    std::vector<std::string>>(),
           py::arg("fn"),
           py::arg("inputs"),
           py::arg("outputs"))
      .PythonInterfaceCopyValue(CallbackOperator)
      .PythonInterfaceBasicRepresentation(CallbackOperator)
      .PythonInterfaceFileIO(CallbackOperator)
      .def("__call__", &CallbackOperator::operator(), py::arg("ws"))
      ;
}
}  // namespace Python
