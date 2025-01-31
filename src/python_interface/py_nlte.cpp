#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/unordered_map.h>
#include <python_interface.h>

#include "hpy_arts.h"

namespace Python {
void py_nlte(py::module_& m) try {
  auto nltelfp =
      py::bind_map<NonlteLineFluxProfile>(m, "NonlteLineFluxProfile");
  workspace_group_interface(nltelfp);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize nlte\n{}", e.what()));
}
}  // namespace Python
