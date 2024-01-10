#include <path_point.h>
#include <pybind11/pybind11.h>
#include <python_interface.h>

#include "py_macros.h"

namespace Python {
void py_path(py::module_& m) try {
  artsclass<PropagationPathPoint>(m, "PropagationPathPoint")
      .def(py::init([]() { return std::make_shared<PropagationPathPoint>(); }),
           "Empty path point")
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
      .def_readwrite("pos",
                     &PropagationPathPoint::los,
                     ":class:`~pyarts.arts.Vector2` Path line-of-sight")
      .def_readwrite("pos",
                     &PropagationPathPoint::nreal,
                     ":class:`float` Path real refractive index")
      .def_readwrite("pos",
                     &PropagationPathPoint::ngroup,
                     ":class:`float` Path group refractive index")
      .PythonInterfaceCopyValue(PropagationPathPoint)
      .PythonInterfaceWorkspaceVariableConversion(PropagationPathPoint)
      .PythonInterfaceFileIO(PropagationPathPoint)
      .PythonInterfaceWorkspaceDocumentation(PropagationPathPoint);

  artsarray<ArrayOfPropagationPathPoint>(m, "ArrayOfPropagationPathPoint")
      .PythonInterfaceFileIO(ArrayOfPropagationPathPoint)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfPropagationPathPoint);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize ppath\n", e.what()));
}
}  // namespace Python