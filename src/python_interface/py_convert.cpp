#include <arts_conversions.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace Python {
namespace py = pybind11;

//! Wraps a conversion function with basic information
#define PythonInterfaceConvert(name, t, from, to, var)                                 \
  convert.def(                                                                         \
      #name,                                                                           \
      py::vectorize(&Conversion::name<t>),                                             \
      py::arg(var),                                                                    \
      py::doc(                                                                         \
          "Converts from " from " to " to "\n\nParameters:\n    " var                  \
          " (:class:`~float` or :class:`~list` or :class:`~numpy.ndarray`): " from \
          "\n\nReturn\n------\n:class:`~float` or :class:`~list` or :class:`~numpy.ndarray`\n    " to))

void py_conversions(py::module_& m) {
  auto convert = m.def_submodule("convert");
  convert.doc() = R"--(Contains several unit conversion functions used in Arts

These all should work with any array-like type
)--";
  PythonInterfaceConvert(freq2kaycm, Numeric, "Frequency [Hz]", "Kayser [cm-1]", "f");
  PythonInterfaceConvert(kaycm2freq, Numeric, "Kayser [cm-1]", "Frequency [Hz]", "v");
  PythonInterfaceConvert(freq2wavelen, Numeric, "Frequency [Hz]", "Wavelenth [m]", "f");
  PythonInterfaceConvert(wavelen2freq, Numeric, "Wavelenth [m]", "Frequency [Hz]", "wl");
  PythonInterfaceConvert(pa2torr, Numeric, "Pressure [Pa]", "Torr [Torr]", "pa");
  PythonInterfaceConvert(torr2pa, Numeric, "Torr [Torr]", "Pressure [Pa]", "torr");
  PythonInterfaceConvert(hz_per_msquared2kaycm_per_cmsquared, Numeric, "Spectral line intensity [Hz/m2]", "CGS spectral line intensity [cm-1/cm-2]", "s");
  PythonInterfaceConvert(kaycm_per_cmsquared2hz_per_msquared, Numeric, "CGS spectral line intensity [cm-1/cm-2]", "Spectral line intensity [Hz/m2]", "s_cgs");
  PythonInterfaceConvert(hz_per_pa2kaycm_per_atm, Numeric, "Line shape parameter [Hz/Pa]", "CGS line shape parameter [cm-1/atm]", "lsm");
  PythonInterfaceConvert(kaycm_per_atm2hz_per_pa, Numeric, "CGS line shape parameter [cm-1/atm]", "Line shape parameter [Hz/Pa]", "lsm_cgs");
  PythonInterfaceConvert(joule2kaycm, Numeric, "Energy [J]", "CGS Energy [cm-1]", "j");
  PythonInterfaceConvert(kaycm2joule, Numeric, "CGS Energy [cm-1]", "Energy [J]", "v");
}
}  // namespace Python
