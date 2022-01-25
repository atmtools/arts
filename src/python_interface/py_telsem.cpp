#include <auto_md.h>
#include <xml_io.h>

#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_telsem(py::module_& m) {
  py::class_<TelsemAtlas>(m, "TelsemAtlas")
      .def(py::init<>())
      .PythonInterfaceFileIO(TelsemAtlas)
      .PythonInterfaceBasicRepresentation(TelsemAtlas);
  
  PythonInterfaceWorkspaceArray(TelsemAtlas);
}
}  // namespace Python