#include <py_auto_interface.h>

#include "py_macros.h"

namespace Python {
void py_xsec(py::module_& m) {
  py::class_<XsecRecord>(m, "XsecRecord")
      .def(py::init([]() { return new XsecRecord{}; }))
      .PythonInterfaceBasicRepresentation(XsecRecord);

  PythonInterfaceWorkspaceArray(XsecRecord);
}
}  // namespace Python