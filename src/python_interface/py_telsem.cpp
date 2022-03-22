#include <py_auto_interface.h>

#include "py_macros.h"

namespace Python {
void py_telsem(py::module_& m) {
  py::class_<TelsemAtlas>(m, "TelsemAtlas")
      .def(py::init<>())
      .PythonInterfaceCopyValue(TelsemAtlas)
      .PythonInterfaceWorkspaceVariableConversion(TelsemAtlas)
      .PythonInterfaceFileIO(TelsemAtlas)
      .PythonInterfaceBasicRepresentation(TelsemAtlas);

  PythonInterfaceWorkspaceArray(TelsemAtlas);
}
}  // namespace Python