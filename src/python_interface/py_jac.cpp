#include <python_interface.h>
#include <py_auto_wsg_init.h>

#include "py_macros.h"

#include "jacobian.h"

namespace Python {
void py_jac(py::module_& m) try {
  py_staticJacobianTargets(m)
      .def_property_readonly("atm", &JacobianTargets::atm, "List of atmospheric targets")
      .def_property_readonly("surf", &JacobianTargets::surf, "List of surface targets")
      .def_property_readonly("line", &JacobianTargets::line, "List of line targets");
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize jac\n", e.what()));
}
}  // namespace Python