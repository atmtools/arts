#include <nanobind/nanobind.h>
#include <physics_funcs.h>

#include "hpy_numpy.h"

namespace Python {
namespace py = nanobind;

void py_physics(py::module_& m) try {
  auto physics  = m.def_submodule("physics");
  physics.doc() = R"--(Contains simple physics functions in arts
)--";

  physics.def(
      "number_density",
      [](py::object p, py::object t) {
        return vectorize(
            [](Numeric a, Numeric b) { return number_density(a, b); }, p, t);
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
        return vectorize([](Numeric f, Numeric t) { return planck(f, t); },
                         frequency,
                         temperature);
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
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize physics\n{}", e.what()));
}
}  // namespace Python
