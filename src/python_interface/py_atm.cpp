
#include <memory>

#include <py_auto_interface.h>
#include <pybind11/attr.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "py_macros.h"

#include <atm.h>
#include <debug.h>
#include <quantum_numbers.h>
#include <species_tags.h>

namespace Python {

void py_atm(py::module_ &m) {
  py::class_<Atm::Data>(m, "AtmData")
      .def(py::init([]() { return std::make_unique<Atm::Data>(); }))
      .def_readwrite("data", &Atm::Data::data)
      .def_readwrite("alt_upp", &Atm::Data::alt_upp)
      .def_readwrite("alt_low", &Atm::Data::alt_low)
      .def_readwrite("lat_upp", &Atm::Data::lat_upp)
      .def_readwrite("lat_low", &Atm::Data::lat_low)
      .def_readwrite("lon_upp", &Atm::Data::lon_upp)
      .def_readwrite("lon_low", &Atm::Data::lon_low);

  auto pnt = py::class_<AtmPoint>(m, "AtmPoint").def(py::init([] {
    return std::make_unique<AtmPoint>();
  }));

  auto fld = py::class_<AtmField>(m, "AtmField").def(py::init([] {
    return std::make_unique<AtmField>();
  }));

  pnt.def_readwrite("temperature", &AtmPoint::temperature)
      .def_readwrite("pressure", &AtmPoint::pressure)
      .def_readwrite("mag", &AtmPoint::mag)
      .def_readwrite("wind", &AtmPoint::wind)
      .def(
          "__getitem__",
          [](AtmPoint &atm, Atm::Key x) {
            ARTS_USER_ERROR_IF(not atm.has(x), "No data for ", x) return atm[x];
          },
          py::is_operator(), py::return_value_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmPoint &atm, const QuantumIdentifier &x) {
            ARTS_USER_ERROR_IF(not atm.has(x), "No data for ", x)
            return atm[x];
          },
          py::is_operator(), py::return_value_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmPoint &atm, const ArrayOfSpeciesTag &x) {
            ARTS_USER_ERROR_IF(not atm.has(x), "No data for ", x)
            return atm[x];
          },
          py::is_operator(), py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](AtmPoint &atm, Atm::Key x, Numeric data) { atm[x] = data; },
          py::is_operator())
      .def(
          "__setitem__",
          [](AtmPoint &atm, const QuantumIdentifier &x, Numeric data) {
            atm[x] = data;
          },
          py::is_operator())
      .def(
          "__setitem__",
          [](AtmPoint &atm, const ArrayOfSpeciesTag &x, Numeric data) {
            atm[x] = data;
          },
          py::is_operator())
      .PythonInterfaceCopyValue(AtmPoint)
      .PythonInterfaceWorkspaceVariableConversion(AtmPoint)
      .PythonInterfaceFileIO(AtmPoint)
      .PythonInterfaceBasicRepresentation(AtmPoint)
      .PythonInterfaceWorkspaceDocumentation(AtmPoint);

  fld.def(
         "__getitem__",
         [](AtmField &atm, Atm::Key x) -> Atm::Data & {
           ARTS_USER_ERROR_IF(not atm.has(x), "No data for ", x) return atm[x];
         },
         py::is_operator(), py::return_value_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmField &atm, const QuantumIdentifier &x) -> Atm::Data & {
            ARTS_USER_ERROR_IF(not atm.has(x), "No data for ", x)
            return atm[x];
          },
          py::is_operator(), py::return_value_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmField &atm, const ArrayOfSpeciesTag &x) -> Atm::Data & {
            ARTS_USER_ERROR_IF(not atm.has(x), "No data for ", x)
            return atm[x];
          },
          py::is_operator(), py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](AtmField &atm, Atm::Key x, const Atm::Data &data) {
            atm[x] = data;
          },
          py::is_operator())
      .def(
          "__setitem__",
          [](AtmField &atm, const QuantumIdentifier &x, const Atm::Data &data) {
            atm[x] = data;
          },
          py::is_operator())
      .def(
          "__setitem__",
          [](AtmField &atm, const ArrayOfSpeciesTag &x, const Atm::Data &data) {
            atm[x] = data;
          },
          py::is_operator())
      .PythonInterfaceCopyValue(AtmField)
      .PythonInterfaceWorkspaceVariableConversion(AtmField)
      .PythonInterfaceFileIO(AtmField)
      .PythonInterfaceBasicRepresentation(AtmField)
      .PythonInterfaceWorkspaceDocumentation(AtmField);
}
} // namespace Python
