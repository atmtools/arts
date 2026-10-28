#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <physics_funcs.h>

#include "hpy_numpy.h"

namespace Python {
namespace py = nanobind;

void py_physics(py::module_& m) try {
  auto physics  = m.def_submodule("physics");
  physics.doc() = R"--(Contains simple physics functions in arts
)--";

  physics.def(
      "refractive_index_water_visible_nir_harvey98",
      [](py::object frequency, py::object temperature, py::object density, bool check_validity) {
        return vectorize(
            [check_validity](Numeric f, Numeric t, Numeric rho) {
              return refractive_index_water_visible_nir_harvey98(f, t, rho, check_validity);
            },
            frequency,
            temperature,
            density);
      },
      "frequency"_a,
      "temperature"_a,
      "density"_a,
      "check_validity"_a = true,
      R"--(Real refractive index of water and steam in the visible and near infrared.

Implements Harvey et al. (1998), "Revised Formulation for the Refractive
Index of Water and Steam as a Function of Wavelength, Temperature and
Density".  The function broadcasts scalar and array inputs using NumPy rules.

Parameters
----------
frequency : Numeric or numpy.ndarray
    Frequency [Hz].
temperature : Numeric or numpy.ndarray
    Temperature [K].
density : Numeric or numpy.ndarray
    Water mass density [kg/m3].
check_validity : bool, optional
    Enforce the published ranges: 261.15 K < temperature < 773.15 K,
    0 kg/m3 < density < 1060 kg/m3, and 0.2 micrometres < wavelength <
    1.9 micrometres.

Returns
-------
n : Numeric or numpy.ndarray
    Real refractive index.
)--");

  physics.def(
      "permittivity_water_liebe93",
      [](py::object frequency, py::object temperature, bool check_validity) {
        return vectorize(
            [check_validity](Numeric f, Numeric t) { return permittivity_water_liebe93(f, t, check_validity); },
            frequency,
            temperature);
      },
      "frequency"_a,
      "temperature"_a,
      "check_validity"_a = true,
      R"--(Complex permittivity of liquid water in the microwave.

Implements the double-Debye model of Liebe, Hufford and Manabe (1991), "A
model for the complex permittivity of water at frequencies below 1 THz", as
ARTS2 and Atmlab's eps_water_liebe93.m.  The imaginary part is positive.  The
function broadcasts scalar and array inputs using NumPy rules.

Parameters
----------
frequency : Numeric or numpy.ndarray
    Frequency [Hz].
temperature : Numeric or numpy.ndarray
    Temperature [K].
check_validity : bool, optional
    Enforce the ranges 10 GHz <= frequency <= 1000 GHz and
    233.15 K <= temperature <= 373.15 K.

Returns
-------
eps : complex or numpy.ndarray
    Complex relative permittivity.
)--");

  physics.def(
      "number_density",
      [](py::object p, py::object t) {
        return vectorize([](Numeric a, Numeric b) { return number_density(a, b); }, p, t);
      },
      "P"_a,
      "T"_a,
      R"--(Calculates the atmospheric number density.

Parameters
----------
  P : Numeric or numpy.ndarray
    Pressure [Pa]

  T : Numeric or numpy.ndarray
    Temperature [K]

Returns
-------
  n : Numeric or numpy.ndarray
    Number density [1/m³]
)--");

  physics.def(
      "planck",
      [](py::object frequency, py::object temperature) {
        return vectorize([](Numeric f, Numeric t) { return planck(f, t); }, frequency, temperature);
      },
      "frequency"_a,
      "temperature"_a,
      R"--(Calculates the Planck function.

.. math::
    I =\frac{2h\nu^3}{c^2} \frac{1}{e^{\frac{h\nu}{kT}} - 1},

where :math:`\nu` is the frequency and :math:`T` is the temperature.

Parameters
----------
frequency : Numeric or numpy.ndarray
    Frequency [Hz]

temperature : Numeric or numpy.ndarray
    Temperature [K]

Returns
-------
B : Numeric or numpy.ndarray
    Planck function [W/(m² Hz sr)]
)--");
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize physics\n{}", e.what()));
}
}  // namespace Python
