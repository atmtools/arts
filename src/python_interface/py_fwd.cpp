#include <fwd.h>

#include <python_interface.h>

#include "debug.h"
#include "py_macros.h"

namespace Python {
void py_fwd(py::module_& m) try {
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize fwd\n", e.what()));
}
}  // namespace Python