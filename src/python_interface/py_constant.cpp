#include <arts_constants.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace Python {
namespace py = pybind11;

struct ConstantDummy {};

#define PythonInterfaceConstant(name) \
  constant.def_readonly_static(       \
      #name, &Constant::name, py::doc(var_string("Value: ", Constant::name).c_str()))

void py_constants(py::module_& m) {
  auto constant = py::class_<ConstantDummy>(m, "constant");
  constant.doc() = R"--(Contains constant internal values of Arts

Note that this "class" simply contains static data and that it should NOT
be initialized
)--";
  PythonInterfaceConstant(pi);
  PythonInterfaceConstant(c);
  PythonInterfaceConstant(h);
  PythonInterfaceConstant(h_bar);
  PythonInterfaceConstant(k);
}
}  // namespace Python
