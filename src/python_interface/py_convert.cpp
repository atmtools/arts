#include <arts_conversions.h>
#include <debug.h>
#include <nanobind/nanobind.h>

#include "hpy_numpy.h"

namespace Python {
namespace py = nanobind;

//! Wraps a conversion function with basic information
#define PythonInterfaceConvert(name, t, from, to, var)                         \
  convert.def(                                                                 \
      #name,                                                                   \
      [](py::object v) {                                                       \
        return vectorize([](t x) { return Conversion::name(x); }, v);          \
      },                                                                       \
      py::arg(var),                                                            \
                                                                               \
      "Converts from " from " to " to "\n\nParameters:\n    " var              \
      " (:class:`~float` or :class:`~list` or :class:`~numpy.ndarray`): " from \
      "\n\nReturn\n------\n:class:`~float` or :class:`~list` or :class:`~numpy.ndarray`\n    " to)

void py_conversions(py::module_& m) try {
  auto convert  = m.def_submodule("convert");
  convert.doc() = R"--(Contains several unit conversion functions used in Arts

These all should work with any array-like type
)--";
  PythonInterfaceConvert(
      freq2kaycm, Numeric, "Frequency [Hz]", "Kayser [cm-1]", "f");
  PythonInterfaceConvert(
      kaycm2freq, Numeric, "Kayser [cm-1]", "Frequency [Hz]", "v");
  PythonInterfaceConvert(
      freq2wavelen, Numeric, "Frequency [Hz]", "Wavelenth [m]", "f");
  PythonInterfaceConvert(
      wavelen2freq, Numeric, "Wavelenth [m]", "Frequency [Hz]", "wl");
  PythonInterfaceConvert(
      pa2torr, Numeric, "Pressure [Pa]", "Torr [Torr]", "pa");
  PythonInterfaceConvert(
      torr2pa, Numeric, "Torr [Torr]", "Pressure [Pa]", "torr");
  PythonInterfaceConvert(hz_per_msquared2kaycm_per_cmsquared,
                         Numeric,
                         "Spectral line intensity [Hz/m2]",
                         "CGS spectral line intensity [cm-1/cm-2]",
                         "s");
  PythonInterfaceConvert(kaycm_per_cmsquared2hz_per_msquared,
                         Numeric,
                         "CGS spectral line intensity [cm-1/cm-2]",
                         "Spectral line intensity [Hz/m2]",
                         "s_cgs");
  PythonInterfaceConvert(hz_per_pa2kaycm_per_atm,
                         Numeric,
                         "Line shape parameter [Hz/Pa]",
                         "CGS line shape parameter [cm-1/atm]",
                         "lsm");
  PythonInterfaceConvert(kaycm_per_atm2hz_per_pa,
                         Numeric,
                         "CGS line shape parameter [cm-1/atm]",
                         "Line shape parameter [Hz/Pa]",
                         "lsm_cgs");
  PythonInterfaceConvert(
      joule2kaycm, Numeric, "Energy [J]", "CGS Energy [cm-1]", "j");
  PythonInterfaceConvert(
      kaycm2joule, Numeric, "CGS Energy [cm-1]", "Energy [J]", "v");
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize convert\n", e.what()));
}
}  // namespace Python
