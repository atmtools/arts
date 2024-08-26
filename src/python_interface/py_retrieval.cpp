#include "hpy_arts.h"
#include "python_interface.h"
#include <nanobind/stl/variant.h>

namespace Python {
void py_retrieval(py::module_& m) try {
  py::class_<JacobianTargetsDiagonalCovarianceMatrixMap> jtdcmm(
      m, "JacobianTargetsDiagonalCovarianceMatrixMap");
  workspace_group_interface(jtdcmm);

  py::class_<JacobianTargetType> jtt(m, "JacobianTargetType");
  jtt.def_rw("value", &JacobianTargetType::target);
  workspace_group_interface(jtt);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize retrieval\n", e.what()));
}
}  // namespace Python
