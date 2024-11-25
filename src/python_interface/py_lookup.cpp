#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/shared_ptr.h>
#include <python_interface.h>

#include "hpy_arts.h"

namespace Python {
void py_lookup(py::module_& m) try {
  py::class_<AbsorptionLookupTable> alt(m, "AbsorptionLookupTable");
  workspace_group_interface(alt);
  alt.def_rw(
      "f_grid", &AbsorptionLookupTable::f_grid, "The frequency grid in Hz");
  alt.def_rw("log_p_grid",
             &AbsorptionLookupTable::log_p_grid,
             "The pressure grid in log Pa [same dimension as atm]");
  alt.def_rw(
      "t_pert",
      &AbsorptionLookupTable::t_pert,
      "The temperautre perturbation grid in K [any number of elements or empty for nothing]");
  alt.def_rw(
      "w_pert",
      &AbsorptionLookupTable::w_pert,

      "The humidity perturbation grid in fractional units [any number of elements or empty for nothing]");
  alt.def_rw("water_atmref",
             &AbsorptionLookupTable::water_atmref,
             "Local grids so that pressure interpolation may work");
  alt.def_rw("t_atmref",
             &AbsorptionLookupTable::t_atmref,
             "Local grids so that pressure interpolation may work");
  alt.def_rw("xsec",
             &AbsorptionLookupTable::xsec,
             "The absorption cross section table");

  auto alts = py::bind_map<AbsorptionLookupTables>(m, "AbsorptionLookupTables");
  workspace_group_interface(alts);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize lookup\n{}", e.what()));
}
}  // namespace Python
