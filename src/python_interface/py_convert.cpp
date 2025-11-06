#include <arts_conversions.h>
#include <debug.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>

#include "hpy_numpy.h"

namespace Python {
namespace py = nanobind;

//! Wraps a conversion function with basic information
#define PythonInterfaceConvert(name, t, from, to, varin, varout)           \
  convert.def(                                                             \
      #name,                                                               \
      [](py::object v) {                                                   \
        return vectorize([](t x) { return Conversion::name(x); }, v);      \
      },                                                                   \
      py::arg(varin),                                                      \
                                                                           \
      "Converts from " from " to " to "\n\nParameters\n----------\n" varin \
      " : float | list | ~numpy.ndarray\n  " from                          \
      "\n\nReturn\n------\n" varout " : ~numpy.ndarray\n  " to)

void py_conversions(py::module_& m) try {
  auto convert  = m.def_submodule("convert");
  convert.doc() = R"--(Contains several unit conversion functions used in Arts

These all should work with any array-like type
)--";
  PythonInterfaceConvert(freq2kaycm,
                         Numeric,
                         "Frequency [Hz]",
                         "Kayser :math:`\\left[\\mathrm{cm}^{-1}\\right]`",
                         "f",
                         "v");
  PythonInterfaceConvert(kaycm2freq,
                         Numeric,
                         "Kayser :math:`\\left[\\mathrm{cm}^{-1}\\right]`",
                         "Frequency [Hz]",
                         "v",
                         "f");
  PythonInterfaceConvert(
      freq2wavelen, Numeric, "Frequency [Hz]", "Wavelenth [m]", "f", "wl");
  PythonInterfaceConvert(
      wavelen2freq, Numeric, "Wavelenth [m]", "Frequency [Hz]", "wl", "f");
  PythonInterfaceConvert(
      pa2torr, Numeric, "Pressure [Pa]", "Torr [Torr]", "pa", "torr");
  PythonInterfaceConvert(
      torr2pa, Numeric, "Torr [Torr]", "Pressure [Pa]", "torr", "pa");
  PythonInterfaceConvert(
      hz_per_msquared2kaycm_per_cmsquared,
      Numeric,
      "Spectral line intensity :math:`\\left[\\mathrm{Hz} \\middle/ \\mathrm{m}^{2}\\right]`",
      "CGS spectral line intensity :math:`\\left[\\mathrm{cm}^{-1} \\middle/ \\mathrm{cm}^{-2}\\right]`",
      "s",
      "s_cgs");
  PythonInterfaceConvert(
      kaycm_per_cmsquared2hz_per_msquared,
      Numeric,
      "CGS spectral line intensity :math:`\\left[\\mathrm{cm}^{-1} \\middle/ \\mathrm{cm}^{-2}\\right]`",
      "Spectral line intensity :math:`\\left[\\mathrm{Hz} \\middle/ \\mathrm{m}^{2}\\right]`",
      "s_cgs",
      "s");
  PythonInterfaceConvert(
      hz_per_pa2kaycm_per_atm,
      Numeric,
      "Line shape parameter [Hz/Pa]",
      "CGS line shape parameter :math:`\\left[\\mathrm{cm}^{-1} \\middle/ \\mathrm{atm}\\right]`",
      "lsm",
      "lsm_cgs");
  PythonInterfaceConvert(
      kaycm_per_atm2hz_per_pa,
      Numeric,
      "CGS line shape parameter :math:`\\left[\\mathrm{cm}^{-1} \\middle/ \\mathrm{atm}\\right]`",
      "Line shape parameter [Hz/Pa]",
      "lsm_cgs",
      "lsm");
  PythonInterfaceConvert(joule2kaycm,
                         Numeric,
                         "Energy [J]",
                         "CGS Energy :math:`\\left[\\mathrm{cm}^{-1}\\right]`",
                         "j",
                         "j_cgs");
  PythonInterfaceConvert(kaycm2joule,
                         Numeric,
                         "CGS Energy :math:`\\left[\\mathrm{cm}^{-1}\\right]`",
                         "Energy [J]",
                         "j_cgs",
                         "j");
  PythonInterfaceConvert(
      rad2deg, Numeric, "radians [rad]", "degrees [deg]", "rad", "deg");
  PythonInterfaceConvert(
      deg2rad, Numeric, "degrees [deg]", "radians [rad]", "deg", "rad");

  convert.def(
      "metric_prefix",
      [](Numeric x) { return Conversion::metric_prefix(x); },
      "x"_a,
      "Returns a suitable metric prefix for the given value x, also returns the scaled value");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize convert\n{}", e.what()));
}
}  // namespace Python
