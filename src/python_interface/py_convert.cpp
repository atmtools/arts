#include <arts_conversions.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace Python {
namespace py = pybind11;

//! Wraps a conversion function with basic information
#define PythonInterfaceConvert(name, t, from, to)  \
  convert.def(#name,                               \
              py::vectorize(&Conversion::name<t>), \
              py::doc("Converts from " from " to " to))

void py_conversions(py::module_& m) {
  auto convert = m.def_submodule("convert");
  convert.doc() = R"--(Contains several unit conversion functions used in Arts

These all should work with any array-like type
)--";
  PythonInterfaceConvert(freq2kaycm, Numeric, "Frequency [Hz]", "Kayser [cm-1]");
  PythonInterfaceConvert(kaycm2freq, Numeric, "Kayser [cm-1]", "Frequency [Hz]");
  PythonInterfaceConvert(freq2wavelen, Numeric, "Frequency [Hz]", "Wavelenth [m]");
  PythonInterfaceConvert(wavelen2freq, Numeric, "Wavelenth [m]", "Frequency [Hz]");
  PythonInterfaceConvert(pa2torr, Numeric, "Pressure [Pa]", "Torr [Torr]");
  PythonInterfaceConvert(torr2pa, Numeric, "Torr [Torr]", "Pressure [Pa]");
}
}  // namespace Python
