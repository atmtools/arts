#include <py_auto_interface.h>

#include "py_macros.h"

namespace Python {
void py_mcantenna(py::module_& m) {
  py::class_<MCAntenna>(m, "MCAntenna")
      .def(py::init([]() { return new MCAntenna{}; }))
      .PythonInterfaceCopyValue(MCAntenna)
      .PythonInterfaceWorkspaceVariableConversion(MCAntenna)
      .PythonInterfaceFileIO(MCAntenna)
      .PythonInterfaceBasicRepresentation(MCAntenna);
}
}  // namespace Python
