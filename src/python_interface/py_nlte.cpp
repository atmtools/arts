
#include <py_auto_interface.h>

#include "py_macros.h"

namespace Python {
void py_nlte(py::module_& m) {
py::class_<EnergyLevelMapType>(m, "EnergyLevelMapType")
      .def(py::init([]() { return new EnergyLevelMapType{}; }))
      .def(py::init([](const std::string& c) {
        return toEnergyLevelMapTypeOrThrow(c);
      }))
      .PythonInterfaceCopyValue(EnergyLevelMapType)
      .PythonInterfaceBasicRepresentation(EnergyLevelMapType)
      .def(py::pickle(
          [](const EnergyLevelMapType& t) {
            return py::make_tuple(std::string(toString(t)));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return new EnergyLevelMapType{
                toEnergyLevelMapType(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, EnergyLevelMapType>();

  py::class_<EnergyLevelMap>(m, "EnergyLevelMap")
      .def(py::init([]() { return new EnergyLevelMap{}; }))
      .PythonInterfaceCopyValue(EnergyLevelMap)
      .PythonInterfaceWorkspaceVariableConversion(EnergyLevelMap)
      .PythonInterfaceBasicRepresentation(EnergyLevelMap)
      .PythonInterfaceFileIO(EnergyLevelMap)
      .PythonInterfaceReadWriteData(EnergyLevelMap, type)
      .PythonInterfaceReadWriteData(EnergyLevelMap, levels)
      .PythonInterfaceReadWriteData(EnergyLevelMap, vib_energy)
      .PythonInterfaceReadWriteData(EnergyLevelMap, value)
      .def(py::pickle(
          [](const EnergyLevelMap& t) {
            return py::make_tuple(t.type, t.levels, t.vib_energy, t.value);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 4, "Invalid state!")
            auto*  out = new EnergyLevelMap{};
            out->type = t[0].cast<EnergyLevelMapType>();
            out->levels = t[1].cast<ArrayOfQuantumIdentifier>();
            out->vib_energy = t[2].cast<Vector>();
            out->value = t[3].cast<Tensor4>();
            return out;
          }));
}
}  // namespace Python