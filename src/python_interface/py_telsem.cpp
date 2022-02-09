#include "py_macros.h"
#include <py_auto_interface.h>

namespace Python {
void py_telsem(py::module_& m) {
  py::class_<TelsemAtlas>(m, "TelsemAtlas")
      .def(py::init<>())
      .PythonInterfaceWorkspaceVariableConversion(TelsemAtlas)
      .PythonInterfaceFileIO(TelsemAtlas)
      .PythonInterfaceBasicRepresentation(TelsemAtlas);

  PythonInterfaceWorkspaceArray(TelsemAtlas);
}
}  // namespace Python