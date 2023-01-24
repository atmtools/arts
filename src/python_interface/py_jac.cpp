#include <py_auto_interface.h>

#include "jacobian.h"
#include "py_macros.h"
#include "quantum_numbers.h"
#include "species.h"
#include "species_tags.h"

namespace Python {
void py_jac(py::module_& m) {
  py::class_<Jacobian::Type>(m, "JacobianType")
      .def(py::init([]() { return Jacobian::Type{}; }))
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
            return Jacobian::Type{
                Jacobian::toType(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, Jacobian::Type>();

  py::class_<Jacobian::Atm>(m, "JacobianAtm")
      .def(py::init([]() { return Jacobian::Atm{}; }))
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
            return Jacobian::Atm{
                Jacobian::toAtm(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, Jacobian::Atm>();

  py::class_<Jacobian::Line>(m, "JacobianLine")
      .def(py::init([]() { return Jacobian::Line{}; }))
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
            return Jacobian::Line{
                Jacobian::toLine(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, Jacobian::Line>();

  py::class_<Jacobian::Sensor>(m, "JacobianSensor")
      .def(py::init([]() { return Jacobian::Sensor{}; }))
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
            return Jacobian::Sensor{
                Jacobian::toSensor(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, Jacobian::Sensor>();

  py::class_<Jacobian::Special>(m, "JacobianSpecial")
      .def(py::init([]() { return Jacobian::Special{}; }))
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
            return Jacobian::Special{
                Jacobian::toSpecial(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, Jacobian::Special>();

  py::class_<JacobianTarget>(m, "JacobianTarget")
      .def(py::init([]() { return std::make_unique<JacobianTarget>(); }))
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
            
            auto out = std::make_unique<JacobianTarget>();
            
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
          }))
      .PythonInterfaceWorkspaceDocumentation(JacobianTarget);

  PythonInterfaceWorkspaceArray(JacobianTarget);

  py::class_<RetrievalQuantity>(m, "RetrievalQuantity")
      .def(py::init([]() { return std::make_unique<RetrievalQuantity>(); }))
      .def(py::pickle(
          [](const RetrievalQuantity& self) {
            return py::make_tuple(self.SubTag(),
                                  self.SubSubTag(),
                                  self.Mode(),
                                  self.Grids(),
                                  self.TransformationFunc(),
                                  self.TFuncParameters(),
                                  self.Transformation(),
                                  self.Offset(),
                                  self.Target());
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 9, "Invalid state!")

            auto out = std::make_unique<RetrievalQuantity>();

            out->SubTag() = t[0].cast<String>();
            out->SubSubTag() = t[1].cast<String>();
            out->Mode() = t[2].cast<String>();
            out->Grids() = t[3].cast<ArrayOfVector>();
            out->TransformationFunc() = t[4].cast<String>();
            out->TFuncParameters() = t[5].cast<Vector>();
            out->Transformation() = t[6].cast<Matrix>();
            out->Offset() = t[7].cast<Vector>();
            out->Target() = t[8].cast<JacobianTarget>();
            
            return out;
          }));

  PythonInterfaceWorkspaceArray(RetrievalQuantity);
}
}  // namespace Python