#include <python_interface.h>
#include <quantum_numbers.h>
#include <quantum_term_symbol.h>

#include "py_macros.h"

#include <pybind11/operators.h>

namespace Python {
void py_quantum(py::module_& m) {
  py::class_<QuantumNumberValue>(m, "QuantumNumberValue")
      .def(py::init([]() { return std::make_unique<QuantumNumberValue>(); }), "Default value")
      .def(py::init(
          [](const std::string& s) { return std::make_unique<QuantumNumberValue>(s); }), "From :class:`str`")
      .PythonInterfaceCopyValue(QuantumNumberValue)
      .PythonInterfaceBasicRepresentation(QuantumNumberValue)
      .def_readwrite("type", &QuantumNumberValue::type, ":class:`~pyarts.arts.options.QuantumNumberType` Type of number")
      .def_property("str_upp",
                    &QuantumNumberValue::str_upp,
                    [](QuantumNumberValue& x, String& y) { x.set(y, true); }, py::doc(":class:`~pyarts.arts.String` Upper value"))
      .def_property("str_low",
                    &QuantumNumberValue::str_low,
                    [](QuantumNumberValue& x, String& y) { x.set(y, false); }, py::doc(":class:`~pyarts.arts.String` Lower value"))
      .def_property("upp",
                    &QuantumNumberValue::upp,
                    [](QuantumNumberValue& x, Rational& y) {
                      x.set(var_string(y), true);
                    }, py::doc(":class:`~pyarts.arts.Rational` Upper value"))
      .def_property("low",
                    &QuantumNumberValue::low,
                    [](QuantumNumberValue& x, Rational& y) {
                      x.set(var_string(y), false);
                    }, py::doc(":class:`~pyarts.arts.Rational` Lower value"))
      .def(py::pickle(
          [](const QuantumNumberValue& t) {
            return py::make_tuple(var_string(t));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_unique<QuantumNumberValue>(t[0].cast<std::string>());
          })).doc() = "A single quantum number with a value";
  py::implicitly_convertible<std::string, QuantumNumberValue>();

  py::class_<QuantumNumberValueList>(m, "QuantumNumberValueList")
      .def(py::init([]() { return std::make_unique<QuantumNumberValueList>(); }), "Default list")
      .def(py::init(
          [](const std::string& s) { return std::make_unique<QuantumNumberValueList>(s); }), "From :class:`str`")
      .PythonInterfaceCopyValue(QuantumNumberValueList)
      .def("get",
           [](QuantumNumberValueList& x, QuantumNumberType y) {
             ARTS_USER_ERROR_IF(not x.has(y), "Out of range: ", y) return x[y];
           }, py::arg("qt"), R"(Set a quantum number value

Parameters
----------
qt : ~pyarts.arts.QuantumNumberType
    The type to get

Returns
-------
qn : ~pyarts.arts.QuantumNumberValue
    The value
)")
      .def("set",
           [](QuantumNumberValueList& x, QuantumNumberValue y) { x.set(y); }, py::arg("qn"), R"(Set a quantum number value

Parameters
----------
qn : ~pyarts.arts.QuantumNumberValue
    The value to set
)")
      .PythonInterfaceBasicRepresentation(QuantumNumberValueList)
      .def(py::pickle(
          [](const QuantumNumberValueList& t) {
            return py::make_tuple(var_string(t));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_unique<QuantumNumberValueList>(t[0].cast<std::string>());
          })).doc() = "A list of unique :class:`~pyarts.arts.QuantumNumberValue`";
  py::implicitly_convertible<std::string, QuantumNumberValueList>();

  py::class_<QuantumNumberLocalState>(m, "QuantumNumberLocalState")
      .def(py::init([]() { return std::make_unique<QuantumNumberLocalState>(); }), "Default state")
      .def(py::init([](const std::string& str) {
        auto out = QuantumNumberLocalState{};
        out.val = QuantumNumberValueList{str};
        return out;
      }), "From :class:`str`")
      .def(py::init([](QuantumNumberValueList& ql) {
        auto out = QuantumNumberLocalState{};
        out.val = ql;
        return out;
      }), "From :class:`~pyarts.arts.QuantumNumberValueList`")
      .PythonInterfaceCopyValue(QuantumNumberLocalState)
      .PythonInterfaceBasicRepresentation(QuantumNumberLocalState)
      .def_readwrite("state", &QuantumNumberLocalState::val, ":class:`~pyarts.arts.QuantumNumberValueList` The values that make up the state")
      .def(py::pickle(
          [](const QuantumNumberLocalState& t) {
            return py::make_tuple(t.val);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            auto out = std::make_unique<QuantumNumberLocalState>();
            out->val = t[0].cast<QuantumNumberValueList>();
            return out;
          })).doc() = "A local state of quantum numbers";
  py::implicitly_convertible<std::string, QuantumNumberLocalState>();

  py::class_<QuantumIdentifier>(m, "QuantumIdentifier")
      .def(py::init([]() { return std::make_unique<QuantumIdentifier>(); }), "Default ID")
      .def(py::init(
          [](const std::string& s) { return std::make_unique<QuantumIdentifier>(s); }), "From :class:`str`")
      .PythonInterfaceCopyValue(QuantumIdentifier)
      .PythonInterfaceWorkspaceVariableConversion(QuantumIdentifier)
      .PythonInterfaceFileIO(QuantumIdentifier)
      .PythonInterfaceBasicRepresentation(QuantumIdentifier)
      .def_readonly("isotopologue_index",
                    &QuantumIdentifier::isotopologue_index, ":class:`int` The isotopologue index")
      .def_readwrite("state", &QuantumIdentifier::val, ":class:`~pyarts.arts.QuantumNumberValueList` The values that make up the state")
      .def_property(
          "isotopologue",
          [](QuantumIdentifier& qid) {
            return Species::Isotopologues.at(qid.isotopologue_index);
          },
          [](QuantumIdentifier& qid, SpeciesIsotopeRecord& iso) {
            Index res = Species::find_species_index(iso);
            ARTS_USER_ERROR_IF(res < 0, "Bad species: ", iso)
            qid.isotopologue_index = res;
          }, ":class:`SpeciesIsotopeRecord` The isotopologue")
      .def(py::self == py::self)
      .def("__hash__", [](QuantumIdentifier& x) { return py::hash(py::str(var_string(x))); })
      .def("as_symbol", &Quantum::Helpers::molecular_term_symbol, R"(Get the molecular symbol as often seen in literature
Returns
-------
symbol : str
    The symbol representation
)")
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