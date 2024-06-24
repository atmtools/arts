#include <nanobind/operators.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <python_interface.h>
#include <quantum_numbers.h>
#include <quantum_term_symbol.h>

#include "py_macros.h"

namespace Python {
void py_quantum(py::module_& m) try {
  py::class_<QuantumNumberValue>(m, "QuantumNumberValue")
      .def(py::init<>())
      .def(py::init<std::string>())
      .PythonInterfaceCopyValue(QuantumNumberValue)
      .PythonInterfaceBasicRepresentation(QuantumNumberValue)
      .def_rw("type",
              &QuantumNumberValue::type,
              ":class:`~pyarts.arts.options.QuantumNumberType` Type of number")
      .def_prop_rw(
          "str_upp",
          &QuantumNumberValue::str_upp,
          [](QuantumNumberValue& x, String& y) { x.set(y, true); },
          ":class:`~pyarts.arts.String` Upper value")
      .def_prop_rw(
          "str_low",
          &QuantumNumberValue::str_low,
          [](QuantumNumberValue& x, String& y) { x.set(y, false); },
          ":class:`~pyarts.arts.String` Lower value")
      .def_prop_rw(
          "upp",
          &QuantumNumberValue::upp,
          [](QuantumNumberValue& x, Rational& y) {
            x.set(var_string(y), true);
          },
          ":class:`~pyarts.arts.Rational` Upper value")
      .def_prop_rw(
          "low",
          &QuantumNumberValue::low,
          [](QuantumNumberValue& x, Rational& y) {
            x.set(var_string(y), false);
          },
          ":class:`~pyarts.arts.Rational` Lower value")
      .def("__getstate__",
           [](const QuantumNumberValue& qnv) {
             return std::tuple<std::string>{var_string(qnv)};
           })
      .def("__setstate__",
           [](QuantumNumberValue* qnv, const std::tuple<std::string>& state) {
             new (qnv) QuantumNumberValue(std::get<0>(state));
           })
      .doc() = "A single quantum number with a value";
  py::implicitly_convertible<std::string, QuantumNumberValue>();

  py::class_<QuantumNumberValueList>(m, "QuantumNumberValueList")
      .def(py::init<>())
      .def(py::init<std::string>())
      .PythonInterfaceCopyValue(QuantumNumberValueList)
      .def(
          "get",
          [](QuantumNumberValueList& x, QuantumNumberType y) {
            ARTS_USER_ERROR_IF(not x.has(y), "Out of range: ", y) return x[y];
          },
          py::arg("qt"),
          R"(Set a quantum number value

Parameters
----------
qt : ~pyarts.arts.QuantumNumberType
    The type to get

Returns
-------
qn : ~pyarts.arts.QuantumNumberValue
    The value
)")
      .def(
          "set",
          [](QuantumNumberValueList& x, QuantumNumberValue y) { x.set(y); },
          py::arg("qn"),
          R"(Set a quantum number value

Parameters
----------
qn : ~pyarts.arts.QuantumNumberValue
    The value to set
)")
      .PythonInterfaceBasicRepresentation(QuantumNumberValueList)
      .def("__getstate__",
           [](const QuantumNumberValueList& qnv) {
             return std::tuple<std::string>{var_string(qnv)};
           })
      .def("__setstate__",
           [](QuantumNumberValueList* qnv,
              const std::tuple<std::string>& state) {
             new (qnv) QuantumNumberValueList(std::get<0>(state));
           })
      .doc() = "A list of unique :class:`~pyarts.arts.QuantumNumberValue`";
  py::implicitly_convertible<std::string, QuantumNumberValueList>();

  py::class_<QuantumNumberLocalState>(m, "QuantumNumberLocalState")
      .def(py::init<>())
      .def(py::init<std::string>())
      .PythonInterfaceCopyValue(QuantumNumberLocalState)
      .PythonInterfaceBasicRepresentation(QuantumNumberLocalState)
      .def_rw(
          "state",
          &QuantumNumberLocalState::val,
          ":class:`~pyarts.arts.QuantumNumberValueList` The values that make up the state")
      .def("__getstate__",
           [](const QuantumNumberLocalState& qnv) {
             return std::tuple<std::string>{var_string(qnv)};
           })
      .def("__setstate__",
           [](QuantumNumberLocalState* qnv,
              const std::tuple<std::string>& state) {
             new (qnv) QuantumNumberLocalState(std::get<0>(state));
           })
      .doc() = "A local state of quantum numbers";
  py::implicitly_convertible<std::string, QuantumNumberLocalState>();

  py::class_<QuantumIdentifier>(m, "QuantumIdentifier")
      .def(py::init<std::string>())
      .def_ro("isotopologue_index",
              &QuantumIdentifier::isotopologue_index,
              ":class:`int` The isotopologue index")
      .def_rw(
          "state",
          &QuantumIdentifier::val,
          ":class:`~pyarts.arts.QuantumNumberValueList` The values that make up the state")
      .def_prop_rw(
          "isotopologue",
          [](QuantumIdentifier& qid) {
            return Species::Isotopologues.at(qid.isotopologue_index);
          },
          [](QuantumIdentifier& qid, SpeciesIsotope& iso) {
            Index res = Species::find_species_index(iso);
            ARTS_USER_ERROR_IF(res < 0, "Bad species: ", iso)
            qid.isotopologue_index = res;
          },
          ":class:`SpeciesIsotope` The isotopologue")
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self <= py::self)
      .def(py::self >= py::self)
      .def(py::self < py::self)
      .def(py::self > py::self)
      .def(py::hash(py::self))
      .def("as_symbol",
           &Quantum::Helpers::molecular_term_symbol,
           R"(Get the molecular symbol as often seen in literature
Returns
-------
symbol : str
    The symbol representation
)")
      .def("__getstate__",
           [](const QuantumIdentifier& qnv) {
             return std::tuple<std::string>{var_string(qnv)};
           })
      .def("__setstate__",
           [](QuantumIdentifier* qnv, const std::tuple<std::string>& state) {
             new (qnv) QuantumIdentifier(std::get<0>(state));
           });
  py::implicitly_convertible<std::string, QuantumIdentifier>();
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize quantum\n", e.what()));
}
}  // namespace Python