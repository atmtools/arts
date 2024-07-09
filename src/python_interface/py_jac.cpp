#include <python_interface.h>

#include "hpy_arts.h"
#include "jacobian.h"

namespace Python {
void py_jac(py::module_& m) try {
  py::class_<JacobianTargets> jacs(m, "JacobianTargets");
  workspace_group_interface(jacs);
  jacs.def_prop_ro("atm", &JacobianTargets::atm, "List of atmospheric targets")
      .def_prop_ro("surf", &JacobianTargets::surf, "List of surface targets")
      .def_prop_ro("line", &JacobianTargets::line, "List of line targets");
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize jac\n", e.what()));
}
}  // namespace Python