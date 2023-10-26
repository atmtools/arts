#include <python_interface.h>

#include "py_macros.h"

#include "new_jacobian.h"

namespace Python {
void py_jac(py::module_& m) try {
  artsclass<JacobianTargets>(m, "JacobianTargets")
      .def(py::init([]() { return std::make_shared<JacobianTargets>(); }), "Default target")
      .PythonInterfaceCopyValue(JacobianTargets)
      .PythonInterfaceWorkspaceVariableConversion(JacobianTargets)
      .PythonInterfaceBasicRepresentation(JacobianTargets)
      .PythonInterfaceFileIO(JacobianTargets)
      .def_property_readonly("atm", &JacobianTargets::atm, "List of atmospheric targets")
      .def_property_readonly("surf", &JacobianTargets::surf, "List of atmospheric targets")
      .PythonInterfaceWorkspaceDocumentation(JacobianTargets);
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize jac\n", e.what()));
}
}  // namespace Python