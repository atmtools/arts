#include <py_auto_interface.h>
#include <pybind11/attr.h>
#include <pybind11/pybind11.h>

#include <hitran_species.h>

namespace Python {
void py_hitran(py::module_ &m) {
  auto hit = m.def_submodule("hitran");

  hit.def("quantumidentity", &Hitran::id_from_lookup, py::arg("mol"),
          py::arg("iso"),
          py::doc("Get the ARTS QuantumIdentifier of the HITRAN molecule"));

  hit.def("ratio", &Hitran::ratio_from_lookup, py::arg("mol"),
          py::arg("iso"),
          py::doc("Get the isotopologue ratio in HITRAN of HITRAN molecule"));
}
} // namespace Python
