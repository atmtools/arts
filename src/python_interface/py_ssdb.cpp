#include <python_interface.h>

#include "core/scattering/arts_ssdb.h"
#include "py_macros.h"

namespace Python {
void py_ssdb(py::module_& m) try {
} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "DEV ERROR:\nCannot initialize scattering species:\n", e.what()));
};
}  // namespace Python
