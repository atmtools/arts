#include <stdexcept>

#include <rtepack.h>

#include "debug.h"
#include "hpy_matpack.h"

namespace Python {
void py_griddedfield(py::module_& m) try {
  py::class_<GriddedField1> gf1(m, "GriddedField1");
  py::class_<GriddedField2> gf2(m, "GriddedField2");
  py::class_<GriddedField3> gf3(m, "GriddedField3");
  py::class_<GriddedField4> gf4(m, "GriddedField4");
  py::class_<GriddedField5> gf5(m, "GriddedField5");
  py::class_<GriddedField6> gf6(m, "GriddedField6");
  gridded_data_interface(gf1);
  gridded_data_interface(gf2);
  gridded_data_interface(gf3);
  gridded_data_interface(gf4);
  gridded_data_interface(gf5);
  gridded_data_interface(gf6);

  py::class_<NamedGriddedField2> ngf2(m, "NamedGriddedField2");
  py::class_<NamedGriddedField3> ngf3(m, "NamedGriddedField3");
  gridded_data_interface(ngf2);
  gridded_data_interface(ngf3);

  py::class_<GriddedField1Named> gf1n(m, "GriddedField1Named");
  gridded_data_interface(gf1n);

  py::class_<ComplexGriddedField2> gf2c(m, "ComplexGriddedField2");
  gridded_data_interface(gf2c);

  py::class_<StokvecGriddedField6> gf6sv(m, "StokvecGriddedField6");
  gridded_data_interface(gf6sv);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize gridded field\n", e.what()));
}
}  // namespace Python