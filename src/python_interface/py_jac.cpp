#include <py_auto_interface.h>

#include "jacobian.h"
#include "py_macros.h"
#include "quantum_numbers.h"
#include "species.h"
#include "species_tags.h"

namespace Python {
void py_jac(py::module_& m) {
  py::class_<Jacobian::Type>(m, "JacobianType")
      .def(py::init([]() { return new Jacobian::Type{}; }))
      .def(py::init(
               [](const std::string& c) { return Jacobian::toTypeOrThrow(c); }),
           py::arg("str").none(false))
      .PythonInterfaceCopyValue(Jacobian::Type)
      .PythonInterfaceBasicRepresentation(Jacobian::Type)
      .def(py::pickle(
          [](const Jacobian::Type& t) {
            return py::make_tuple(std::string(Jacobian::toString(t)));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return new Jacobian::Type{
                Jacobian::toType(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, Jacobian::Type>();

  py::class_<Jacobian::Atm>(m, "JacobianAtm")
      .def(py::init([]() { return new Jacobian::Atm{}; }))
      .def(py::init(
               [](const std::string& c) { return Jacobian::toAtmOrThrow(c); }),
           py::arg("str").none(false))
      .PythonInterfaceCopyValue(Jacobian::Atm)
      .PythonInterfaceBasicRepresentation(Jacobian::Atm)
      .def(py::pickle(
          [](const Jacobian::Atm& t) {
            return py::make_tuple(std::string(Jacobian::toString(t)));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return new Jacobian::Atm{
                Jacobian::toAtm(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, Jacobian::Atm>();

  py::class_<Jacobian::Line>(m, "JacobianLine")
      .def(py::init([]() { return new Jacobian::Line{}; }))
      .def(py::init(
               [](const std::string& c) { return Jacobian::toLineOrThrow(c); }),
           py::arg("str").none(false))
      .PythonInterfaceCopyValue(Jacobian::Line)
      .PythonInterfaceBasicRepresentation(Jacobian::Line)
      .def(py::pickle(
          [](const Jacobian::Line& t) {
            return py::make_tuple(std::string(Jacobian::toString(t)));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return new Jacobian::Line{
                Jacobian::toLine(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, Jacobian::Line>();

  py::class_<Jacobian::Sensor>(m, "JacobianSensor")
      .def(py::init([]() { return new Jacobian::Sensor{}; }))
      .def(py::init([](const std::string& c) {
             return Jacobian::toSensorOrThrow(c);
           }),
           py::arg("str").none(false))
      .PythonInterfaceCopyValue(Jacobian::Sensor)
      .PythonInterfaceBasicRepresentation(Jacobian::Sensor)
      .def(py::pickle(
          [](const Jacobian::Sensor& t) {
            return py::make_tuple(std::string(Jacobian::toString(t)));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return new Jacobian::Sensor{
                Jacobian::toSensor(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, Jacobian::Sensor>();

  py::class_<Jacobian::Special>(m, "JacobianSpecial")
      .def(py::init([]() { return new Jacobian::Special{}; }))
      .def(py::init([](const std::string& c) {
             return Jacobian::toSpecialOrThrow(c);
           }),
           py::arg("str").none(false))
      .PythonInterfaceCopyValue(Jacobian::Special)
      .PythonInterfaceBasicRepresentation(Jacobian::Special)
      .def(py::pickle(
          [](const Jacobian::Special& t) {
            return py::make_tuple(std::string(Jacobian::toString(t)));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return new Jacobian::Special{
                Jacobian::toSpecial(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, Jacobian::Special>();

  py::class_<JacobianTarget>(m, "JacobianTarget")
      .def(py::init([]() { return new JacobianTarget{}; }))
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
      .PythonInterfaceReadWriteData(JacobianTarget, species_id)
      .def(py::pickle(
          [](const JacobianTarget& t) {
            return py::make_tuple(t.type, t.atm, t.line, t.sensor, t.special, t.perturbation, t.qid, t.species_array_id, t.string_id, t.species_id);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 10, "Invalid state!")
            
            auto* out = new JacobianTarget{};
            
            out -> type = t[0].cast<Jacobian::Type>();
            out -> atm = t[1].cast<Jacobian::Atm>();
            out -> line = t[2].cast<Jacobian::Line>();
            out -> sensor = t[3].cast<Jacobian::Sensor>();
            out -> special = t[4].cast<Jacobian::Special>();
            out -> perturbation = t[5].cast<Numeric>();
            out -> qid = t[6].cast<QuantumIdentifier>();
            out -> species_array_id = t[7].cast<ArrayOfSpeciesTag>();
            out -> string_id = t[8].cast<String>();
            out -> species_id = t[9].cast<Species::Species>();

            return out;
          }));

  PythonInterfaceWorkspaceArray(JacobianTarget);

  py::class_<RetrievalQuantity>(m, "RetrievalQuantity").def(py::init([]() {
    return new RetrievalQuantity{};
  }));

  PythonInterfaceWorkspaceArray(RetrievalQuantity);
}
}  // namespace Python