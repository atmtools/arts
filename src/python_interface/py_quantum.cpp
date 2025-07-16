#include <nanobind/operators.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>
#include <quantum.h>
#include <quantum_term_symbol.h>

#include <stdexcept>
#include <variant>

#include "enumsQuantumNumberType.h"
#include "hpy_arts.h"
#include "hpy_vector.h"
#include "nanobind/nanobind.h"
#include "nanobind/stl/bind_map.h"
#include "rational.h"

namespace Python {
void py_quantum(py::module_& m) try {
  py::class_<Quantum::Value> qval(m, "QuantumValue");
  qval.def(py::init<String>())
      .def(py::init<Rational>())
      .def(py::init<QuantumNumberType>())
      .def_prop_rw(
          "value",
          [](Quantum::Value& x) { return x.value; },
          [](Quantum::Value& x, std::variant<Rational, String> y) {
            x = std::visit(
                [](auto& v) -> Quantum::Value { return Quantum::Value{v}; }, y);
          },
          ":class:`~pyarts.arts.String` Upper value")
      .def("__getstate__",
           [](const Quantum::Value& qnv) {
             return std::tuple<std::string>{std::format("{}", qnv)};
           })
      .def("__setstate__",
           [](Quantum::Value* qnv, const std::tuple<std::string>& state) {
             new (qnv) Quantum::Value(std::get<0>(state));
           })
      .doc() = "A single quantum number value";
  generic_interface(qval);

  py::class_<Quantum::UpperLower> qul(m, "QuantumUpperLower");
  qul.def_rw("upper", &Quantum::UpperLower::upper, "Upper state");
  qul.def_rw("lower", &Quantum::UpperLower::lower, "Lower state");
  qul.doc() = "Uppler and lower quantum number values";
  generic_interface(qul);

  auto qlvl = py::bind_map<QuantumLevel, py::rv_policy::reference_internal>(
      m, "QuantumLevel");
  generic_interface(qlvl);

  auto qstate = py::bind_map<QuantumState, py::rv_policy::reference_internal>(
      m, "QuantumState");
  generic_interface(qstate);

  py::class_<QuantumIdentifier> qid(m, "QuantumIdentifier");
  qid.def(py::init_implicit<const std::string_view>());
  qid.def_rw("isot", &QuantumIdentifier::isot, "Isotopologue");
  qid.def_rw("state", &QuantumIdentifier::state, "State");
  generic_interface(qid);

  py::class_<QuantumLevelIdentifier> qlid(m, "QuantumLevelIdentifier");
  qlid.def(py::init_implicit<const std::string_view>());
  qlid.def_rw("isot", &QuantumLevelIdentifier::isot, "Isotopologue");
  qlid.def_rw("state", &QuantumLevelIdentifier::state, "State");
  generic_interface(qlid);

  qid.def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self <= py::self)
      .def(py::self >= py::self)
      .def(py::self < py::self)
      .def(py::self > py::self)
      .def("__hash__",
           [](const QuantumIdentifier& x) {
             return std::hash<QuantumIdentifier>{}(x);
           })
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

  qlid.def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self <= py::self)
      .def(py::self >= py::self)
      .def(py::self < py::self)
      .def(py::self > py::self)
      .def("__hash__",
           [](const QuantumLevelIdentifier& x) {
             return std::hash<QuantumLevelIdentifier>{}(x);
           })
      .def("__getstate__",
           [](const QuantumLevelIdentifier& qnv) {
             return std::tuple<std::string>{std::format("{}", qnv)};
           })
      .def("__setstate__",
           [](QuantumLevelIdentifier* qnv,
              const std::tuple<std::string>& state) {
             new (qnv) QuantumLevelIdentifier(std::get<0>(state));
           });

  auto a1 = py::bind_vector<ArrayOfQuantumIdentifier,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfQuantumIdentifier");
  generic_interface(a1);
  vector_interface(a1);

  auto a2 = py::bind_vector<ArrayOfQuantumLevelIdentifier,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfQuantumLevelIdentifier");
  generic_interface(a1);
  vector_interface(a1);
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize quantum\n{}", e.what()));
}
}  // namespace Python
