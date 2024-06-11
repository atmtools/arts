#include <python_interface.h>

#include <igrf13.h>

namespace Python {
void py_igrf(py::module_& m) try {
  m.def("igrf", &IGRF::igrf);
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize IGRF\n", e.what()));
}
}  // namespace Python
