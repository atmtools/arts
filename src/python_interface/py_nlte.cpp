#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/unordered_map.h>
#include <python_interface.h>

#include "hpy_arts.h"

namespace Python {
void py_nlte(py::module_& m) try {
  auto vec =
      py::bind_map<QuantumIdentifierVectorMap>(m, "QuantumIdentifierVectorMap");
  generic_interface(vec);

  auto num =
      py::bind_map<QuantumIdentifierNumericMap>(m, "QuantumIdentifierNumericMap");
  generic_interface(num);

  auto gf1 =
      py::bind_map<QuantumIdentifierGriddedField1Map>(m, "QuantumIdentifierGriddedField1Map");
  generic_interface(gf1);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize nlte\n{}", e.what()));
}
}  // namespace Python
