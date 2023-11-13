#include <python_interface.h>

#include "lbl_data.h"
#include "py_macros.h"

#include <lbl.h>

namespace Python {
void py_lbl(py::module_& m) try {
  artsclass<lbl::band>(m, "band")
      .def(py::init([]() { return std::make_shared<lbl::band>(); }), "Default target")
      .PythonInterfaceBasicRepresentation(lbl::band);

  artsarray<AbsorptionBands>(m, "AbsorptionBands");
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize lbl\n", e.what()));
}
}  // namespace Python