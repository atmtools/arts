#include <path_point.h>
#include <python_interface.h>

#include "py_macros.h"

namespace Python {
void py_path(py::module_& m) try {
  py::class_<PropagationPathPoint>(m, "PropagationPathPoint")
      .def_rw(
          "pos_type",
          &PropagationPathPoint::pos_type,
          ":class:`~pyarts.arts.options.PathPositionType` Path position type")
      .def_rw(
          "los_type",
          &PropagationPathPoint::los_type,
          ":class:`~pyarts.arts.options.PathPositionType` Path line-of-sight type")
      .def_rw("pos",
              &PropagationPathPoint::pos,
              ":class:`~pyarts.arts.Vector3` Path position")
      .def_rw("los",
              &PropagationPathPoint::los,
              ":class:`~pyarts.arts.Vector2` Path line-of-sight")
      .def_rw("nreal",
              &PropagationPathPoint::nreal,
              ":class:`float` Path real refractive index")
      .def_rw("ngroup",
              &PropagationPathPoint::ngroup,
              ":class:`float` Path group refractive index");
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize ppath\n", e.what()));
}
}  // namespace Python