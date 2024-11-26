#include <nanobind/operators.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <python_interface.h>
#include <quantum_numbers.h>
#include <quantum_term_symbol.h>

#include <stdexcept>

#include "hpy_arts.h"
#include "hpy_vector.h"
#include "nanobind/nanobind.h"
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
              ":class:`~pyarts.arts.QuantumNumberType` Type of number")
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
            x.set(std::format("{}", y), true);
          },
          ":class:`~pyarts.arts.Rational` Upper value")
      .def_prop_rw(
          "low",
          &QuantumNumberValue::low,
          [](QuantumNumberValue& x, Rational& y) {
            x.set(std::format("{}", y), false);
          },
          ":class:`~pyarts.arts.Rational` Lower value")
      .def("__getstate__",
           [](const QuantumNumberValue& qnv) {
             return std::tuple<std::string>{std::format("{}", qnv)};
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
      .def(
          "__len__",
          [](QuantumNumberValueList& x) { return x.values.size(); },
          "Number of values")
      .def(
          "__iter__",
          [](QuantumNumberValueList& x) {
            return py::make_iterator<py::rv_policy::reference_internal>(
                py::type<QuantumNumberValueList>(),
                "iterator",
                x.values.begin(),
                x.values.end());
          },
          py::keep_alive<0, 1>(),
          "Iterate over the values")
      .def("__getitem__",
           [](QuantumNumberValueList& x, QuantumNumberType y) {
             if (not x.has(y)) throw std::out_of_range(std::format("{}", y));
             return x[y];
           })
      .def("__getitem__",
           [](QuantumNumberValueList& x, Size y) {
             if (x.values.size() <= y)
               throw std::out_of_range(std::format("{}", y));
             return x.values[y];
           })
      .def("__setitem__",
           [](QuantumNumberValueList& x, Size y, QuantumNumberValue z) {
             if (x.values.size() <= y)
               throw std::out_of_range(std::format("{}", y));
             return x.values[y] = z;
           })
      .PythonInterfaceCopyValue(QuantumNumberValueList)
      .def(
          "get",
          [](QuantumNumberValueList& x, QuantumNumberType y) {
            ARTS_USER_ERROR_IF(not x.has(y), "Out of range: {}", y) return x[y];
          },
          "qt"_a,
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
          "qn"_a,
          R"(Set a quantum number value

Parameters
----------
qn : ~pyarts.arts.QuantumNumberValue
    The value to set
)")
      .PythonInterfaceBasicRepresentation(QuantumNumberValueList)
      .def("__getstate__",
           [](const QuantumNumberValueList& qnv) {
             return std::tuple<std::string>{std::format("{}", qnv)};
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
      .def(
          "__len__",
          [](QuantumNumberLocalState& x) { return x.val.values.size(); },
          "Number of values")
      .def(
          "__iter__",
          [](QuantumNumberLocalState& x) {
            return py::make_iterator<py::rv_policy::reference_internal>(
                py::type<QuantumNumberLocalState>(),
                "iterator",
                x.val.values.begin(),
                x.val.values.end());
          },
          py::keep_alive<0, 1>(),
          "Iterate over the values")
      .def("__getitem__",
           [](QuantumNumberLocalState& x, QuantumNumberType y) {
             if (not x.val.has(y))
               throw std::out_of_range(std::format("{}", y));
             return x.val[y];
           })
      .def("__getitem__",
           [](QuantumNumberLocalState& x, Size y) {
             if (x.val.values.size() <= y)
               throw std::out_of_range(std::format("{}", y));
             return x.val.values[y];
           })
      .def("__setitem__",
           [](QuantumNumberLocalState& x, Size y, QuantumNumberValue z) {
             if (x.val.values.size() <= y)
               throw std::out_of_range(std::format("{}", y));
             return x.val.values[y] = z;
           })
      .def_rw(
          "state",
          &QuantumNumberLocalState::val,
          ":class:`~pyarts.arts.QuantumNumberValueList` The values that make up the state")
      .def("__getstate__",
           [](const QuantumNumberLocalState& qnv) {
             return std::tuple<std::string>{std::format("{}", qnv)};
           })
      .def("__setstate__",
           [](QuantumNumberLocalState* qnv,
              const std::tuple<std::string>& state) {
             new (qnv) QuantumNumberLocalState(std::get<0>(state));
           })
      .doc() = "A local state of quantum numbers";
  py::implicitly_convertible<std::string, QuantumNumberLocalState>();

  py::class_<QuantumIdentifier> qids(m, "QuantumIdentifier");
  workspace_group_interface(qids);
  qids.def(py::init<std::string>())
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
            ARTS_USER_ERROR_IF(res < 0, "Bad species: {}", iso)
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
             return std::tuple<std::string>{std::format("{}", qnv)};
           })
      .def("__setstate__",
           [](QuantumIdentifier* qnv, const std::tuple<std::string>& state) {
             new (qnv) QuantumIdentifier(std::get<0>(state));
           });
  py::implicitly_convertible<std::string, QuantumIdentifier>();

  auto a1 = py::bind_vector<ArrayOfQuantumIdentifier,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfQuantumIdentifier");
  workspace_group_interface(a1);
  vector_interface(a1);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize quantum\n{}", e.what()));
}
}  // namespace Python
