#include <python_interface.h>

#include "py_macros.h"

namespace Python {
void py_scattering_species(py::module_& m) try {

  artsclass<HenyeyGreenstein>(m, "HenyeyGreenstein")
      .def(py::init([](Numeric g) { return std::make_shared<HenyeyGreenstein>(g); }), "Henyey-Greenstein Phase Function")
      .def(
           "evaluate_phase_function",
           static_cast<Vector (HenyeyGreenstein::*)(const Vector &)>(&HenyeyGreenstein::evaluate_phase_function),
           "Evaluate the Henyey-Greenstein phase function."
           )
      .def(
           "evaluate_phase_function",
           static_cast<Numeric (HenyeyGreenstein::*)(const Numeric &)>(&HenyeyGreenstein::evaluate_phase_function),
           "Evaluate the Henyey-Greenstein phase function."
           )
      .PythonInterfaceCopyValue(HenyeyGreenstein)
      .PythonInterfaceWorkspaceVariableConversion(HenyeyGreenstein)
      .PythonInterfaceFileIO(HenyeyGreenstein)
      .PythonInterfaceBasicRepresentation(HenyeyGreenstein)
      .PythonInterfaceWorkspaceDocumentation(HenyeyGreenstein);

} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize scattering species:\n", e.what()));
};
}
