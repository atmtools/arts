#include <py_auto_interface.h>
#include <quantum_numbers.h>
#include <quantum_term_symbol.h>

#include "py_macros.h"

#include <pybind11/operators.h>

namespace Python {
void py_quantum(py::module_& m) {
  py::class_<QuantumNumberType>(m, "QuantumNumberType")
      .def(py::init([]() { return QuantumNumberType{}; }))
      .def(py::init([](const std::string& c) {
        return Quantum::Number::toTypeOrThrow(c);
      }))
      .PythonInterfaceCopyValue(QuantumNumberType)
      .PythonInterfaceBasicRepresentation(QuantumNumberType)
      .def(py::pickle(
          [](const QuantumNumberType& t) {
            return py::make_tuple(std::string(Quantum::Number::toString(t)));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return QuantumNumberType{
                Quantum::Number::toType(t[0].cast<std::string>())};
          }));
  py::implicitly_convertible<std::string, QuantumNumberType>();

  py::class_<QuantumNumberValue>(m, "QuantumNumberValue")
      .def(py::init([]() { return std::make_unique<QuantumNumberValue>(); }))
      .def(py::init(
          [](const std::string& s) { return std::make_unique<QuantumNumberValue>(s); }))
      .PythonInterfaceCopyValue(QuantumNumberValue)
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
                    })
      .def(py::pickle(
          [](const QuantumNumberValue& t) {
            return py::make_tuple(var_string(t));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_unique<QuantumNumberValue>(t[0].cast<std::string>());
          }));
  py::implicitly_convertible<std::string, QuantumNumberValue>();

  py::class_<QuantumNumberValueList>(m, "QuantumNumberValueList")
      .def(py::init([]() { return std::make_unique<QuantumNumberValueList>(); }))
      .def(py::init(
          [](const std::string& s) { return std::make_unique<QuantumNumberValueList>(s); }))
      .PythonInterfaceCopyValue(QuantumNumberValueList)
      .def("get",
           [](QuantumNumberValueList& x, QuantumNumberType y) {
             ARTS_USER_ERROR_IF(not x.has(y), "Out of range: ", y) return x[y];
           })
      .def("set",
           [](QuantumNumberValueList& x, QuantumNumberValue y) { x.set(y); })
      .PythonInterfaceBasicRepresentation(QuantumNumberValueList)
      .def(py::pickle(
          [](const QuantumNumberValueList& t) {
            return py::make_tuple(var_string(t));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_unique<QuantumNumberValueList>(t[0].cast<std::string>());
          }));
  py::implicitly_convertible<std::string, QuantumNumberValueList>();

  py::class_<QuantumNumberLocalState>(m, "QuantumNumberLocalState")
      .def(py::init([]() { return std::make_unique<QuantumNumberLocalState>(); }))
      .def(py::init([](const std::string& str) {
        auto out = QuantumNumberLocalState{};
        out.val = QuantumNumberValueList{str};
        return out;
      }))
      .def(py::init([](QuantumNumberValueList& ql) {
        auto out = QuantumNumberLocalState{};
        out.val = ql;
        return out;
      }))
      .PythonInterfaceCopyValue(QuantumNumberLocalState)
      .PythonInterfaceBasicRepresentation(QuantumNumberLocalState)
      .def_readwrite("state", &QuantumNumberLocalState::val)
      .def(py::pickle(
          [](const QuantumNumberLocalState& t) {
            return py::make_tuple(t.val);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            auto out = std::make_unique<QuantumNumberLocalState>();
            out->val = t[0].cast<QuantumNumberValueList>();
            return out;
          }));
  py::implicitly_convertible<std::string, QuantumNumberLocalState>();

  py::class_<QuantumIdentifier>(m, "QuantumIdentifier")
      .def(py::init([]() { return std::make_unique<QuantumIdentifier>(); }))
      .def(py::init(
          [](const std::string& s) { return std::make_unique<QuantumIdentifier>(s); }))
      .PythonInterfaceCopyValue(QuantumIdentifier)
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
          })
      .def(py::self == py::self)
      .def("__hash__", [](QuantumIdentifier& x) { return py::hash(py::str(var_string(x))); })
      .def("as_symbol", &Quantum::Helpers::molecular_term_symbol)
      .def(py::pickle(
          [](const QuantumIdentifier& t) {
            return py::make_tuple(t.isotopologue_index, t.val);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 2, "Invalid state!")
            auto out = std::make_unique<QuantumIdentifier>();

            out->isotopologue_index = t[0].cast<Index>();
            out->val = t[1].cast<Quantum::Number::ValueList>();

            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(QuantumIdentifier);
  py::implicitly_convertible<std::string, QuantumIdentifier>();

  PythonInterfaceWorkspaceArray(QuantumIdentifier);
}
}  // namespace Python