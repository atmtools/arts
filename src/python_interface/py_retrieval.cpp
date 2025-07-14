#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/variant.h>

#include "hpy_arts.h"
#include "python_interface.h"

namespace Python {
void py_retrieval(py::module_& m) try {
  auto jtdcmm = py::bind_map<JacobianTargetsDiagonalCovarianceMatrixMap,
                             py::rv_policy::reference_internal>(
      m, "JacobianTargetsDiagonalCovarianceMatrixMap");
  generic_interface(jtdcmm);
  
  py::class_<PairOfBlockMatrix> pobm(m, "PairOfBlockMatrix");
  generic_interface(pobm);
  pobm.def_rw("first", &PairOfBlockMatrix::first, "Matrix");
  pobm.def_rw("second", &PairOfBlockMatrix::second, "Inverse of Matrix");

  py::class_<JacobianTargetType> jtt(m, "JacobianTargetType");
  jtt.def_rw("value", &JacobianTargetType::target, "Target");
  generic_interface(jtt);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize retrieval\n{}", e.what()));
}
}  // namespace Python
