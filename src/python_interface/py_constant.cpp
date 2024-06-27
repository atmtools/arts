#include <arts_constants.h>
#include <debug.h>
#include <nanobind/nanobind.h>

namespace Python {
namespace py = nanobind;

struct ConstantDummy {};

#define PythonInterfaceConstant(name) constants.attr(#name) = Constant::name;

void py_constants(py::module_& m) try {
  auto constants = m.def_submodule(
      "constants", R"--(Contain copies of constants of Arts internals
)--");

  PythonInterfaceConstant(pi);
  PythonInterfaceConstant(c);
  PythonInterfaceConstant(h);
  PythonInterfaceConstant(h_bar);
  PythonInterfaceConstant(k);
  PythonInterfaceConstant(doppler_broadening_const_squared);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize constant\n", e.what()));
}
}  // namespace Python
