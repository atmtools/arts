#include <python_interface.h>

#include "jacobian.h"
#include "new_jacobian.h"
#include "py_macros.h"
#include "quantum_numbers.h"
#include "species.h"
#include "species_tags.h"

namespace Python {
void py_jac(py::module_& m) try {
  artsclass<JacobianTarget>(m, "JacobianTarget")
      .def(py::init([]() { return std::make_shared<JacobianTarget>(); }), "Default target")
      .PythonInterfaceCopyValue(JacobianTarget)
      .PythonInterfaceWorkspaceVariableConversion(JacobianTarget)
      .PythonInterfaceBasicRepresentation(JacobianTarget)
      .PythonInterfaceReadWriteData(JacobianTarget, type, ":class:`~pyarts.arts.options.JacobianType` Type of target")
      .PythonInterfaceReadWriteData(JacobianTarget, atm, ":class:`~pyarts.arts.options.JacobianAtm` Type of atmospheric target")
      .PythonInterfaceReadWriteData(JacobianTarget, line, ":class:`~pyarts.arts.options.JacobianLine` Type of line target")
      .PythonInterfaceReadWriteData(JacobianTarget, sensor, ":class:`~pyarts.arts.options.JacobianSensor` Type of sensor target")
      .PythonInterfaceReadWriteData(JacobianTarget, special, ":class:`~pyarts.arts.options.JacobianSpecial` Type of special target")
      .PythonInterfaceReadWriteData(JacobianTarget, perturbation, ":class:`float` Perturbation magnitude")
      .PythonInterfaceReadWriteData(JacobianTarget, qid, ":class:`~pyarts.arts.QuantumIdentifier` Identifier")
      .PythonInterfaceReadWriteData(JacobianTarget, species_array_id, ":class:`~pyarts.arts.ArrayOfSpeciesTag` Identifier")
      .PythonInterfaceReadWriteData(JacobianTarget, string_id, ":class:`~pyarts.arts.String` Identifier")
      .PythonInterfaceReadWriteData(JacobianTarget, species_id, ":class:`~pyarts.arts.Species` Identifier")
      .def(py::pickle(
          [](const JacobianTarget& t) {
            return py::make_tuple(t.type, t.atm, t.line, t.sensor, t.special, t.perturbation, t.qid, t.species_array_id, t.string_id, t.species_id);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 10, "Invalid state!")
            
            auto out = std::make_shared<JacobianTarget>();
            
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
          })).doc() = "A jacobian target";
  
  artsclass<RetrievalQuantity>(m, "RetrievalQuantity")
      .def(py::init([]() { return std::make_shared<RetrievalQuantity>(); }), "Empty quantity")
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

            auto out = std::make_shared<RetrievalQuantity>();

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
          }))
      .def(py::init<const RetrievalQuantity&>())
      .PythonInterfaceCopyValue(RetrievalQuantity)
      .PythonInterfaceBasicRepresentation(RetrievalQuantity)
      .PythonInterfaceFileIO(RetrievalQuantity)
      .PythonInterfaceWorkspaceDocumentation(RetrievalQuantity);

  artsarray<ArrayOfRetrievalQuantity>(m, "ArrayOfRetrievalQuantity")
      .PythonInterfaceFileIO(ArrayOfRetrievalQuantity)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfRetrievalQuantity);

  artsclass<JacobianTargets>(m, "JacobianTargets")
      .PythonInterfaceWorkspaceVariableConversion(JacobianTargets)
      .PythonInterfaceBasicRepresentation(JacobianTargets)
      .PythonInterfaceFileIO(JacobianTargets)
      .def_readwrite("atm", &JacobianTargets::atm, "List of atmospheric targets")
      .PythonInterfaceWorkspaceDocumentation(JacobianTargets);
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize jac\n", e.what()));
}
}  // namespace Python