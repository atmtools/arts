#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <python_interface.h>
#include <radar.h>

#include <algorithm>

#include "hpy_arts.h"
#include "hpy_vector.h"

namespace Python {
void py_radar(py::module_& m) try {
  auto aora = py::bind_vector<ArrayOfRadarAuxiliary, py::rv_policy::reference_internal>(m, "ArrayOfRadarAuxiliary");
  generic_interface(aora);
  vector_interface(aora);
  aora.def("__init__", [](ArrayOfRadarAuxiliary* self, const std::vector<std::string>& x) {
    new (self) ArrayOfRadarAuxiliary();
    self->reserve(x.size());
    std::transform(
        x.begin(), x.end(), std::back_inserter(*self), [](const std::string& s) { return to<RadarAuxiliary>(s); });
  });
  py::implicitly_convertible<std::vector<std::string>, ArrayOfRadarAuxiliary>();
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize radar\n{}", e.what()));
}
}  // namespace Python
