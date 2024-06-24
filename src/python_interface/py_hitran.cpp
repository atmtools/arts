#include <python_interface.h>

#include <hitran_species.h>

namespace Python {
void py_hitran(py::module_ &m) try {
  auto hit = m.def_submodule("hitran");
  hit.doc() = "Helpers to interface with HITRAN";

  hit.def("quantumidentity",
          &Hitran::id_from_lookup,
          py::arg("mol"),
          py::arg("iso"),
          R"--(Get the Arts quantum identifier of the HITRAN molecule

Parameters
----------
mol : ~pyarts.arts.Index
    Moluecule index
iso : bytes
    Isotopologue character (single byte)

Return
------
:class:`~pyarts.arts.QuantumIdentifier`
    Identifier
)--");

  hit.def("ratio",
          &Hitran::ratio_from_lookup,
          py::arg("mol"),
          py::arg("iso"),
          R"--(Get the isotopologue ratio in HITRAN of HITRAN molecule

Parameters
----------
mol : ~pyarts.arts.Index
    Moluecule index
iso : bytes
    Isotopologue character (single byte)

Return
------
:class:`float`
    Atmospheric ratio of isotopologue
)--");
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize hitran\n", e.what()));
}
}  // namespace Python
