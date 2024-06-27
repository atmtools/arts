#include <nanobind/stl/bind_vector.h>
#include <path_point.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "hpy_vector.h"
#include "py_macros.h"

namespace Python {
void py_path(py::module_& m) try {
  py::class_<PropagationPathPoint> pppp(m, "PropagationPathPoint");
  workspace_group_interface(pppp);
  pppp.def_rw(
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

  auto a1 = py::bind_vector<ArrayOfPropagationPathPoint,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfPropagationPathPoint");
  workspace_group_interface(a1);
  vector_interface(a1);
  auto a2 = py::bind_vector<ArrayOfArrayOfPropagationPathPoint,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfPropagationPathPoint");
  workspace_group_interface(a2);
  vector_interface(a2);
  auto a3 = py::bind_vector<ArrayOfArrayOfArrayOfPropagationPathPoint,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfArrayOfPropagationPathPoint");
  workspace_group_interface(a3);
  vector_interface(a3);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize ppath\n", e.what()));
}
}  // namespace Python
