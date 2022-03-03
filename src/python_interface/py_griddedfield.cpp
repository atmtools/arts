#include <py_auto_interface.h>

#include "py_macros.h"

namespace Python {
void py_griddedfield(py::module_& m) {
  py::class_<GriddedField1>(m, "GriddedField1")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceWorkspaceVariableConversion(GriddedField1)
      .PythonInterfaceFileIO(GriddedField1)
      .PythonInterfaceBasicRepresentation(GriddedField1)
      .PythonInterfaceGriddedField(GriddedField1);

  py::class_<GriddedField2>(m, "GriddedField2")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceWorkspaceVariableConversion(GriddedField2)
      .PythonInterfaceFileIO(GriddedField2)
      .PythonInterfaceBasicRepresentation(GriddedField2)
      .PythonInterfaceGriddedField(GriddedField2);

  py::class_<GriddedField3>(m, "GriddedField3")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceWorkspaceVariableConversion(GriddedField3)
      .PythonInterfaceFileIO(GriddedField3)
      .PythonInterfaceBasicRepresentation(GriddedField3)
      .PythonInterfaceGriddedField(GriddedField3);

  py::class_<GriddedField4>(m, "GriddedField4")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceWorkspaceVariableConversion(GriddedField4)
      .PythonInterfaceFileIO(GriddedField4)
      .PythonInterfaceBasicRepresentation(GriddedField4)
      .PythonInterfaceGriddedField(GriddedField4);

  py::class_<GriddedField5>(m, "GriddedField5")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceWorkspaceVariableConversion(GriddedField5)
      .PythonInterfaceFileIO(GriddedField5)
      .PythonInterfaceBasicRepresentation(GriddedField5)
      .PythonInterfaceGriddedField(GriddedField5);

  py::class_<GriddedField6>(m, "GriddedField6")
      .def(py::init<>())
      .def(py::init<const String&>())
      .PythonInterfaceWorkspaceVariableConversion(GriddedField6)
      .PythonInterfaceFileIO(GriddedField6)
      .PythonInterfaceBasicRepresentation(GriddedField6)
      .PythonInterfaceGriddedField(GriddedField6);

  PythonInterfaceWorkspaceArray(GriddedField1);
  PythonInterfaceWorkspaceArray(GriddedField2);
  PythonInterfaceWorkspaceArray(GriddedField3);
  PythonInterfaceWorkspaceArray(GriddedField4);

  PythonInterfaceWorkspaceArray(ArrayOfGriddedField1)
      .def(py::init([](const std::vector<std::vector<GriddedField1>>& x) {
        ArrayOfArrayOfGriddedField1 y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }));
  py::implicitly_convertible<std::vector<std::vector<GriddedField1>>,
                             ArrayOfArrayOfGriddedField1>();

  PythonInterfaceWorkspaceArray(ArrayOfGriddedField2)
      .def(py::init([](const std::vector<std::vector<GriddedField2>>& x) {
        ArrayOfArrayOfGriddedField2 y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }));
  py::implicitly_convertible<std::vector<std::vector<GriddedField2>>,
                             ArrayOfArrayOfGriddedField2>();

  PythonInterfaceWorkspaceArray(ArrayOfGriddedField3)
      .def(py::init([](const std::vector<std::vector<GriddedField3>>& x) {
        ArrayOfArrayOfGriddedField3 y(x.size());
        std::copy(x.begin(), x.end(), y.begin());
        return y;
      }));
  py::implicitly_convertible<std::vector<std::vector<GriddedField3>>,
                             ArrayOfArrayOfGriddedField3>();
}
}  // namespace Python