#include <py_auto_interface.h>

#include "debug.h"
#include "py_macros.h"
#include "quantum_numbers.h"

namespace Python {
void py_quantum(py::module_& m) {
  py::class_<QuantumNumberType>(m, "QuantumNumberType")
      .def(py::init<>())
      .def(py::init([](char* c) { return Quantum::Number::toTypeOrThrow(c); }))
      .PythonInterfaceBasicRepresentation(QuantumNumberType);
  py::implicitly_convertible<py::str, QuantumNumberType>();

  py::class_<QuantumNumberLocalState>(m, "QuantumNumberLocalState")
      .def(py::init<>())
      .def(py::init([](QuantumNumberValueList& ql) {
        auto out = QuantumNumberLocalState{};
        out.val = ql;
        return out;
      }))
      .PythonInterfaceBasicRepresentation(QuantumNumberLocalState)
      .def_readwrite("state", &QuantumNumberLocalState::val);

  py::class_<QuantumNumberValue>(m, "QuantumNumberValue")
      .def(py::init<>())
      .def(py::init<char*>())
      .PythonInterfaceBasicRepresentation(QuantumNumberValue)
      .def_readwrite("type", &QuantumNumberValue::type)
      .def_property("str_upp",
                    &QuantumNumberValue::str_upp,
                    [](QuantumNumberValue& x, String& y) { x.set(y, true); })
      .def_property("str_low",
                    &QuantumNumberValue::str_low,
                    [](QuantumNumberValue& x, String& y) { x.set(y, false); })
      .def_property("upp",
                    &QuantumNumberValue::upp,
                    [](QuantumNumberValue& x, Rational& y) {
                      x.set(var_string(y), true);
                    })
      .def_property("low",
                    &QuantumNumberValue::low,
                    [](QuantumNumberValue& x, Rational& y) {
                      x.set(var_string(y), false);
                    });
  py::implicitly_convertible<char*, QuantumNumberValue>();

  py::class_<QuantumNumberValueList>(m, "QuantumNumberValueList")
      .def(py::init<>())
      .def(py::init<char*>())
      .def("get",
           [](QuantumNumberValueList& x, QuantumNumberType y) {
             ARTS_USER_ERROR_IF(not x.has(y), "Out of range: ", y) return x[y];
           })
      .def("set",
           [](QuantumNumberValueList& x, QuantumNumberValue y) { x.set(y); })
      .PythonInterfaceBasicRepresentation(QuantumNumberValueList);
  py::implicitly_convertible<char*, QuantumNumberValueList>();

  py::class_<QuantumIdentifier>(m, "QuantumIdentifier")
      .def(py::init<>())
      .def(py::init<char*>())
      .PythonInterfaceWorkspaceVariableConversion(QuantumIdentifier)
      .PythonInterfaceFileIO(QuantumIdentifier)
      .PythonInterfaceBasicRepresentation(QuantumIdentifier)
      .def_readonly("isotopologue_index",
                    &QuantumIdentifier::isotopologue_index)
      .def_readwrite("state", &QuantumIdentifier::val)
      .def_property(
          "isotopologue",
          [](QuantumIdentifier& qid) {
            return Species::Isotopologues.at(qid.isotopologue_index);
          },
          [](QuantumIdentifier& qid, SpeciesIsotopeRecord& iso) {
            Index res = Species::find_species_index(iso);
            ARTS_USER_ERROR_IF(res < 0, "Bad species: ", iso)
            qid.isotopologue_index = res;
          });
  py::implicitly_convertible<char*, QuantumIdentifier>();

  PythonInterfaceWorkspaceArray(QuantumIdentifier);
}
}  // namespace Python