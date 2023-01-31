
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
      .def(py::init([](const GriddedField3 &x) {
        return std::make_unique<Atm::Data>(x);
      }))
      .def(py::init(
          [](const Tensor3 &x) { return std::make_unique<Atm::Data>(x); }))
      .def(py::init(
          [](const Numeric &x) { return std::make_unique<Atm::Data>(x); }))
      .def(py::init([](const Index &x) {
        return std::make_unique<Atm::Data>(static_cast<Numeric>(x));
      }))
      .def(py::init([](const Atm::FunctionalData &x) {
        return std::make_unique<Atm::Data>(x);
      }))
      .def_readwrite("data", &Atm::Data::data)
      .def_readwrite("alt_upp", &Atm::Data::alt_upp)
      .def_readwrite("alt_low", &Atm::Data::alt_low)
      .def_readwrite("lat_upp", &Atm::Data::lat_upp)
      .def_readwrite("lat_low", &Atm::Data::lat_low)
      .def_readwrite("lon_upp", &Atm::Data::lon_upp)
      .def_readwrite("lon_low", &Atm::Data::lon_low);
  py::implicitly_convertible<GriddedField3, Atm::Data>();
  py::implicitly_convertible<Tensor3, Atm::Data>();
  py::implicitly_convertible<Numeric, Atm::Data>();
  py::implicitly_convertible<Index, Atm::Data>();
  py::implicitly_convertible<Atm::FunctionalData, Atm::Data>();

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
            if (not atm.has(x))
              throw py::key_error(var_string(x));
            return atm[x];
          },
          py::return_value_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmPoint &atm, const QuantumIdentifier &x) {
            if (not atm.has(x))
              throw py::key_error(var_string(x));
            return atm[x];
          },
          py::return_value_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmPoint &atm, const ArrayOfSpeciesTag &x) {
            if (not atm.has(x))
              throw py::key_error(var_string(x));
            return atm[x];
          },
          py::return_value_policy::reference_internal)
      .def("__setitem__",
           [](AtmPoint &atm, Atm::Key x, Numeric data) { atm[x] = data; })
      .def("__setitem__", [](AtmPoint &atm, const QuantumIdentifier &x,
                             Numeric data) { atm[x] = data; })
      .def("__setitem__", [](AtmPoint &atm, const ArrayOfSpeciesTag &x,
                             Numeric data) { atm[x] = data; })
      .PythonInterfaceCopyValue(AtmPoint)
      .PythonInterfaceWorkspaceVariableConversion(AtmPoint)
      .PythonInterfaceFileIO(AtmPoint)
      .PythonInterfaceBasicRepresentation(AtmPoint)
      .PythonInterfaceWorkspaceDocumentation(AtmPoint);

  fld.def(
         "__getitem__",
         [](AtmField &atm, Atm::Key x) -> Atm::Data & {
           if (not atm.has(x))
             throw py::key_error(var_string(x));
           return atm[x];
         }, py::return_value_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmField &atm, const QuantumIdentifier &x) -> Atm::Data & {
            if (not atm.has(x))
              throw py::key_error(var_string(x));
            return atm[x];
          },
          py::return_value_policy::reference_internal)
      .def(
          "__getitem__",
          [](AtmField &atm, const ArrayOfSpeciesTag &x) -> Atm::Data & {
            if (not atm.has(x))
              throw py::key_error(var_string(x));
            return atm[x];
          },
          py::return_value_policy::reference_internal)
      .def("__setitem__", [](AtmField &atm, Atm::Key x,
                             const Atm::Data &data) { atm[x] = data; })
      .def("__setitem__", [](AtmField &atm, const QuantumIdentifier &x,
                             const Atm::Data &data) { atm[x] = data; })
      .def("__setitem__", [](AtmField &atm, const ArrayOfSpeciesTag &x,
                             const Atm::Data &data) { atm[x] = data; })
      .def("regularize", &AtmField::regularize)
      .def("regularized_shape", &AtmField::regularized_shape)
      .def_readwrite("top_of_atmosphere", &AtmField::top_of_atmosphere)
      .PythonInterfaceCopyValue(AtmField)
      .PythonInterfaceWorkspaceVariableConversion(AtmField)
      .PythonInterfaceFileIO(AtmField)
      .PythonInterfaceBasicRepresentation(AtmField)
      .PythonInterfaceWorkspaceDocumentation(AtmField);
}
} // namespace Python
