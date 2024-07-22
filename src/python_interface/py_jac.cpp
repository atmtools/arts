#include <python_interface.h>

#include "hpy_arts.h"
#include "jacobian.h"

namespace Python {
void py_jac(py::module_& m) try {
  py::class_<JacobianTargets> jacs(m, "JacobianTargets");
  workspace_group_interface(jacs);
  jacs.def_prop_rw(
      "atm",
      [](JacobianTargets& j) { return j.atm(); },
      [](JacobianTargets& j, std::vector<Jacobian::AtmTarget> t) {
        j.atm() = std::move(t);
      },
      "List of atmospheric targets");
  jacs.def_prop_rw(
      "surf",
      [](JacobianTargets& j) { return j.surf(); },
      [](JacobianTargets& j, std::vector<Jacobian::SurfaceTarget> t) {
        j.surf() = std::move(t);
      },
      "List of surface targets");
  jacs.def_prop_rw(
      "line",
      [](JacobianTargets& j) { return j.line(); },
      [](JacobianTargets& j, std::vector<Jacobian::LineTarget> t) {
        j.line() = std::move(t);
      },
      "List of line targets");
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize jac\n", e.what()));
}
}  // namespace Python