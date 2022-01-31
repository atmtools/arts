#include <auto_md.h>
#include <xml_io.h>

#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_mcantenna(py::module_& m) {
  py::class_<MCAntenna>(m, "MCAntenna")
      .def(py::init<>())
      .PythonInterfaceFileIO(MCAntenna)
      .PythonInterfaceBasicRepresentation(MCAntenna);
}
}  // namespace Python
