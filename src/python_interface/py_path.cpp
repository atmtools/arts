#include <path_point.h>
#include <pybind11/pybind11.h>
#include <python_interface.h>
#include <py_auto_wsg_init.h>

#include "py_macros.h"

namespace Python {
void py_path(py::module_& m) try {
  py_staticPropagationPathPoint(m)
      .def_readwrite(
          "pos_type",
          &PropagationPathPoint::pos_type,
          ":class:`~pyarts.arts.options.PathPositionType` Path position type")
      .def_readwrite(
          "los_type",
          &PropagationPathPoint::los_type,
          ":class:`~pyarts.arts.options.PathPositionType` Path line-of-sight type")
      .def_readwrite("pos",
                     &PropagationPathPoint::pos,
                     ":class:`~pyarts.arts.Vector3` Path position")
      .def_readwrite("los",
                     &PropagationPathPoint::los,
                     ":class:`~pyarts.arts.Vector2` Path line-of-sight")
      .def_readwrite("nreal",
                     &PropagationPathPoint::nreal,
                     ":class:`float` Path real refractive index")
      .def_readwrite("ngroup",
                     &PropagationPathPoint::ngroup,
                     ":class:`float` Path group refractive index");

  py_staticArrayOfPropagationPathPoint(m);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize ppath\n", e.what()));
}
}  // namespace Python