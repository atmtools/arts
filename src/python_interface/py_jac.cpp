#include <py_auto_interface.h>

#include "jacobian.h"
#include "py_macros.h"

namespace Python {
void py_jac(py::module_& m) {
  py::class_<Jacobian::Type>(m, "JacobianType")
      .def(py::init<>())
      .def(py::init([](const char* c) { return Jacobian::toTypeOrThrow(c); }), py::arg("str").none(false))
      .PythonInterfaceCopyValue(Jacobian::Type)
      .PythonInterfaceBasicRepresentation(Jacobian::Type);
  py::implicitly_convertible<py::str, Jacobian::Type>();

  py::class_<Jacobian::Atm>(m, "JacobianAtm")
      .def(py::init<>())
      .def(py::init([](const char* c) { return Jacobian::toAtmOrThrow(c); }), py::arg("str").none(false))
      .PythonInterfaceCopyValue(Jacobian::Atm)
      .PythonInterfaceBasicRepresentation(Jacobian::Atm);
  py::implicitly_convertible<py::str, Jacobian::Atm>();

  py::class_<Jacobian::Line>(m, "JacobianLine")
      .def(py::init<>())
      .def(py::init([](const char* c) { return Jacobian::toLineOrThrow(c); }), py::arg("str").none(false))
      .PythonInterfaceCopyValue(Jacobian::Line)
      .PythonInterfaceBasicRepresentation(Jacobian::Line);
  py::implicitly_convertible<py::str, Jacobian::Line>();

  py::class_<Jacobian::Sensor>(m, "JacobianSensor")
      .def(py::init<>())
      .def(py::init([](const char* c) { return Jacobian::toSensorOrThrow(c); }), py::arg("str").none(false))
      .PythonInterfaceCopyValue(Jacobian::Sensor)
      .PythonInterfaceBasicRepresentation(Jacobian::Sensor);
  py::implicitly_convertible<py::str, Jacobian::Sensor>();

  py::class_<Jacobian::Special>(m, "JacobianSpecial")
      .def(py::init<>())
      .def(py::init([](const char* c) { return Jacobian::toSpecialOrThrow(c); }), py::arg("str").none(false))
      .PythonInterfaceCopyValue(Jacobian::Special)
      .PythonInterfaceBasicRepresentation(Jacobian::Special);
  py::implicitly_convertible<py::str, Jacobian::Special>();

  py::class_<JacobianTarget>(m, "JacobianTarget")
      .def(py::init<>())
      .PythonInterfaceCopyValue(JacobianTarget)
      .PythonInterfaceWorkspaceVariableConversion(JacobianTarget)
      .PythonInterfaceFileIO(JacobianTarget)
      .PythonInterfaceBasicRepresentation(JacobianTarget)
      .PythonInterfaceReadWriteData(JacobianTarget, type)
      .PythonInterfaceReadWriteData(JacobianTarget, atm)
      .PythonInterfaceReadWriteData(JacobianTarget, line)
      .PythonInterfaceReadWriteData(JacobianTarget, sensor)
      .PythonInterfaceReadWriteData(JacobianTarget, special)
      .PythonInterfaceReadWriteData(JacobianTarget, perturbation)
      .PythonInterfaceReadWriteData(JacobianTarget, qid)
      .PythonInterfaceReadWriteData(JacobianTarget, species_array_id)
      .PythonInterfaceReadWriteData(JacobianTarget, string_id)
      .PythonInterfaceReadWriteData(JacobianTarget, species_id);

  PythonInterfaceWorkspaceArray(JacobianTarget);

  py::class_<RetrievalQuantity>(m, "RetrievalQuantity");

  PythonInterfaceWorkspaceArray(RetrievalQuantity);
}
}  // namespace Python