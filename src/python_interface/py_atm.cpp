
#include <memory>

#include <python_interface.h>

#include "py_macros.h"

#include <atm.h>
#include <debug.h>
#include <quantum_numbers.h>
#include <species_tags.h>

namespace Python {

void py_atm(py::module_ &m) try {
  artsclass<Atm::Data>(m, "AtmData")
      .def(py::init([]() { return std::make_shared<Atm::Data>(); }))
      .def(py::init([](const GriddedField3 &x) {
        return std::make_shared<Atm::Data>(x);
      }))
      .def(py::init(
          [](const Numeric &x) { return std::make_shared<Atm::Data>(x); }))
      .def(py::init([](const Index &x) {
        return std::make_shared<Atm::Data>(static_cast<Numeric>(x));
      }))
      .def(py::init([](const Atm::FunctionalData &x) {
        return std::make_shared<Atm::Data>(x);
      }))
      .def_readwrite("data", &Atm::Data::data)
      .def_readwrite("alt_upp", &Atm::Data::alt_upp)
      .def_readwrite("alt_low", &Atm::Data::alt_low)
      .def_readwrite("lat_upp", &Atm::Data::lat_upp)
      .def_readwrite("lat_low", &Atm::Data::lat_low)
      .def_readwrite("lon_upp", &Atm::Data::lon_upp)
      .def_readwrite("lon_low", &Atm::Data::lon_low)
      .def(py::pickle(
          [](const Atm::Data &t) {
            return py::make_tuple(t.data, t.alt_low, t.alt_upp, t.lat_low,
                                  t.lat_upp, t.lon_low, t.lon_upp);
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 7, "Invalid state!")

            auto out = std::make_shared<Atm::Data>();
            out->data = t[0].cast<Atm::FieldData>();
            out->alt_low = t[1].cast<Atm::Extrapolation>();
            out->alt_upp = t[2].cast<Atm::Extrapolation>();
            out->lat_low = t[3].cast<Atm::Extrapolation>();
            out->lat_upp = t[4].cast<Atm::Extrapolation>();
            out->lon_low = t[5].cast<Atm::Extrapolation>();
            out->lon_upp = t[6].cast<Atm::Extrapolation>();
            return out;
          }));
  py::implicitly_convertible<GriddedField3, Atm::Data>();
  py::implicitly_convertible<Numeric, Atm::Data>();
  py::implicitly_convertible<Index, Atm::Data>();
  py::implicitly_convertible<Atm::FunctionalData, Atm::Data>();

  auto pnt = artsclass<AtmPoint>(m, "AtmPoint").def(py::init([] {
    return std::make_shared<AtmPoint>();
  }));

  auto fld = artsclass<AtmField>(m, "AtmField").def(py::init([] {
    return std::make_shared<AtmField>();
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
      .def("keys", &AtmPoint::keys)
      .PythonInterfaceCopyValue(AtmPoint)
      .PythonInterfaceWorkspaceVariableConversion(AtmPoint)
      .PythonInterfaceFileIO(AtmPoint)
      .PythonInterfaceBasicRepresentation(AtmPoint)
      .def(py::pickle(
          [](const AtmPoint &t) {
            auto k = t.keys();
            std::vector<Numeric> v;
            v.reserve(k.size());
            for (auto &kn : k)
              v.emplace_back(
                  std::visit([&](auto &&key) { return t[key]; }, kn));
            return py::make_tuple(k, v);
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 2, "Invalid state!")

            auto k = t[0].cast<std::vector<Atm::KeyVal>>();
            auto v = t[1].cast<std::vector<Numeric>>();
            ARTS_USER_ERROR_IF(k.size() != v.size(), "Invalid state!")

            auto out = std::make_shared<AtmPoint>();
            for (std::size_t i = 0; i < k.size(); i++)
              std::visit(
                  [&](auto &&key) -> Numeric & { return out->operator[](key); },
                  k[i]) = v[i];

            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(AtmPoint);

  fld.def(
         "__getitem__",
         [](AtmField &atm, Atm::Key x) -> Atm::Data & {
           if (not atm.has(x))
             throw py::key_error(var_string(x));
           return atm[x];
         },
         py::return_value_policy::reference_internal)
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
      .def("at", [](const AtmField &atm, const Vector& h, const Vector& lat,
                    const Vector& lon) { return atm.at(h, lat, lon); })
      .def_readwrite("top_of_atmosphere", &AtmField::top_of_atmosphere)
      .PythonInterfaceCopyValue(AtmField)
      .PythonInterfaceWorkspaceVariableConversion(AtmField)
      .PythonInterfaceFileIO(AtmField)
      .PythonInterfaceBasicRepresentation(AtmField)
      .def(py::pickle(
          [](const AtmField &t) {
            auto k = t.keys();
            std::vector<Atm::Data> v;
            v.reserve(k.size());
            for (auto &kn : k)
              v.emplace_back(
                  std::visit([&](auto &&key) { return t[key]; }, kn));
            return py::make_tuple(k, v, t.top_of_atmosphere);
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 3, "Invalid state!")

            auto k = t[0].cast<std::vector<Atm::KeyVal>>();
            auto v = t[1].cast<std::vector<Atm::Data>>();
            ARTS_USER_ERROR_IF(k.size() != v.size(), "Invalid state!")

            auto out = std::make_shared<AtmField>();
            out->top_of_atmosphere = t[2].cast<Numeric>();

            for (std::size_t i = 0; i < k.size(); i++)
              std::visit(
                  [&](auto &&key) -> Atm::Data & {
                    return out->operator[](key);
                  },
                  k[i]) = v[i];

            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(AtmField);

  PythonInterfaceWorkspaceArray(AtmPoint);
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize atm\n", e.what()));
}
} // namespace Python
