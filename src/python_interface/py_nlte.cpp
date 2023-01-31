
#include <memory>
#include <py_auto_interface.h>
#include <pybind11/attr.h>
#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>

#include "debug.h"
#include "py_macros.h"
#include "python_interface_value_type.h"
#include "quantum_numbers.h"

namespace Python {
void py_nlte(py::module_& m) {
py::class_<EnergyLevelMapType>(m, "EnergyLevelMapType")
      .def(py::init([]() { return EnergyLevelMapType{}; }))
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
            return EnergyLevelMapType{
                toEnergyLevelMapType(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, EnergyLevelMapType>();

  py::class_<EnergyLevelMap>(m, "EnergyLevelMap")
      .def(py::init([]() { return std::make_unique<EnergyLevelMap>(); }))
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
            auto  out = std::make_unique<EnergyLevelMap>();
            out->type = t[0].cast<EnergyLevelMapType>();
            out->levels = t[1].cast<ArrayOfQuantumIdentifier>();
            out->vib_energy = t[2].cast<Vector>();
            out->value = t[3].cast<Tensor4>();
            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(EnergyLevelMap);

  py::class_<VibrationalEnergyLevels>(m, "VibrationalEnergyLevels")
      .def(py::init([](){return std::make_unique<VibrationalEnergyLevels>();}))
      .def(
          "__getitem__",
          [](VibrationalEnergyLevels &x, const QuantumIdentifier &q) {
            if (x.data.find(q) == x.end()) throw py::key_error(var_string(q));
            return x[q];
          }, py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](VibrationalEnergyLevels &x, const QuantumIdentifier &q,
             Numeric y) { x[q] = y; })
      .PythonInterfaceCopyValue(VibrationalEnergyLevels)
      .PythonInterfaceWorkspaceVariableConversion(VibrationalEnergyLevels)
      .PythonInterfaceBasicRepresentation(VibrationalEnergyLevels)
      .PythonInterfaceFileIO(VibrationalEnergyLevels)
      .PythonInterfaceWorkspaceDocumentation(VibrationalEnergyLevels);
}
}  // namespace Python