#include <physics_funcs.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace Python {
namespace py = pybind11;

void py_physics(py::module_& m) try {
  auto physics = m.def_submodule("physics");
  physics.doc() = R"--(Contains simple physics functions in arts
)--";

  physics.def("number_density",
              py::vectorize(&number_density),
              py::arg("P"),
              py::arg("T"),
              py::doc(R"--(Calculates the atmospheric number density.

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
)--"));

  physics.def("planck",
              py::vectorize([](Numeric frequency, Numeric temperature) {
                return planck(frequency, temperature);
              }),
              py::arg("frequency"),
              py::arg("temperature"),
              py::doc(R"--(Calculates the Planck function.

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
)--"));
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize physics\n", e.what()));
}
}  // namespace Python
