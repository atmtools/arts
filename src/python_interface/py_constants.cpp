#include <constants.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace Python {
namespace py = pybind11;

struct ConstantDummy {};

#define PythonInterfaceConstant(name) \
  constant.def_readonly_static(       \
      #name, &Constant::name, py::doc(var_string("Value: ", Constant::name).c_str()))

//! Wraps a conversion function with basic information
#define PythonInterfaceConvert(name, t, from, to)  \
  convert.def(#name,                               \
              py::vectorize(&Conversion::name<t>), \
              py::doc("Converts from " from " to " to))

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

  auto convert = m.def_submodule("convert");
  convert.doc() = R"--(Contains several unit conversion functions used in Arts

These all should work with any array-like type
)--";
  PythonInterfaceConvert(freq2kaycm, Numeric, "Freq. [Hz]", "Kayser [cm-1]");
  PythonInterfaceConvert(kaycm2freq, Numeric, "Kayser [cm-1]", "Freq. [Hz]");
  PythonInterfaceConvert(pa2torr, Numeric, "Pressure [Pa]", "Torr [Torr]");
  PythonInterfaceConvert(torr2pa, Numeric, "Torr [Torr]", "Pressure [Pa]");
}
}  // namespace Python
