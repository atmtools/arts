#include <constants.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace Python {
namespace py = pybind11;

struct ConstantDummy {};
struct ConversionDummy {};

#define PythonInterfaceConstant(name) \
  def_readonly_static(                \
      #name, &Constant::name, py::doc("Value of " #name " in Arts"))

#define PythonInterfaceConvert(name, t, from, to) \
  def_static(#name,                               \
             py::vectorize(&Conversion::name<t>), \
             py::doc("Converts from " from " to " to))

void py_constants(py::module_& m) {
  py::class_<ConstantDummy>(m, "constant")
      .PythonInterfaceConstant(pi)
      .PythonInterfaceConstant(c)
      .PythonInterfaceConstant(h)
      .PythonInterfaceConstant(h_bar)
      .PythonInterfaceConstant(k);

  py::class_<ConversionDummy>(m, "convert")
      .PythonInterfaceConvert(
          freq2kaycm, Numeric, "Frequency [Hz]", "Kayser [cm-1]")
      .PythonInterfaceConvert(
          kaycm2freq, Numeric, "Kayser [cm-1]", "Frequency [Hz]")
      .PythonInterfaceConvert(pa2torr, Numeric, "Pressure [Pa]", "Torr [Torr]")
      .PythonInterfaceConvert(torr2pa, Numeric, "Torr [Torr]", "Pressure [Pa]");
}
}  // namespace Python