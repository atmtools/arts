#include "hpy_arts.h"
#include "hpy_matpack.h"
#include "hpy_vector.h"
#include "python_interface.h"

namespace Python {
void py_griddedfield_stokvec(py::module_& m) {
  py::class_<StokvecSortedGriddedField1> gf1sv(m, "StokvecSortedGriddedField1");
  gridded_data_interface(gf1sv);
  generic_interface(gf1sv);

  py::class_<StokvecSortedGriddedField2> gf2sv(m, "StokvecSortedGriddedField2");
  gridded_data_interface(gf2sv);
  generic_interface(gf2sv);

  py::class_<StokvecSortedGriddedField3> gf3sv(m, "StokvecSortedGriddedField3");
  gridded_data_interface(gf3sv);
  generic_interface(gf3sv);

  py::class_<StokvecSortedGriddedField4> gf4sv(m, "StokvecSortedGriddedField4");
  gridded_data_interface(gf4sv);
  generic_interface(gf4sv);

  py::class_<StokvecSortedGriddedField5> gf5sv(m, "StokvecSortedGriddedField5");
  gridded_data_interface(gf5sv);
  generic_interface(gf5sv);

  py::class_<StokvecSortedGriddedField6> gf6sv(m, "StokvecSortedGriddedField6");
  gridded_data_interface(gf6sv);
  generic_interface(gf6sv);
}
}  // namespace Python
