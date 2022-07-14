#include <physics_funcs.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace Python {
namespace py = pybind11;

void py_physics(py::module_& m) {
  auto physics = m.def_submodule("physics");
  physics.doc() = R"--(Contains simple physics functions in arts
)--";

  physics.def("number_density",
              py::vectorize(&number_density),
              py::doc(R"--(Calculates the atmospheric number density.

Parameters
----------
  P : Numeric
    Pressure [Pa]

  T : Numeric
    Temperature [K]

Returns
-------
  n : Numeric
    Number density [1/m³]
)--"),
              py::arg("P"),
              py::arg("T"));
}
}  // namespace Python
