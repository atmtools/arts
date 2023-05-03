
#include <memory>

#include <py_auto_interface.h>
#include <pybind11/attr.h>
#include <pybind11/detail/common.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "py_macros.h"

#include <surf.h>
#include <debug.h>
#include <quantum_numbers.h>
#include <species_tags.h>

namespace Python {

void py_surf(py::module_ &m) {
  py::class_<Surf::Data>(m, "SurfData")
      .def(py::init([]() { return std::make_unique<Surf::Data>(); }))
      .def(py::init([](const GriddedField2 &x) {
        return std::make_unique<Surf::Data>(x);
      }))
      .def(py::init(
          [](const Numeric &x) { return std::make_unique<Surf::Data>(x); }))
      .def(py::init([](const Index &x) {
        return std::make_unique<Surf::Data>(static_cast<Numeric>(x));
      }))
      .def(py::init([](const Surf::FunctionalData &x) {
        return std::make_unique<Surf::Data>(x);
      }))
      .def_readwrite("data", &Surf::Data::data)
      .def_readwrite("lat_upp", &Surf::Data::lat_upp)
      .def_readwrite("lat_low", &Surf::Data::lat_low)
      .def_readwrite("lon_upp", &Surf::Data::lon_upp)
      .def_readwrite("lon_low", &Surf::Data::lon_low)
      .def(py::pickle(
          [](const Surf::Data &t) {
            return py::make_tuple(t.data, t.lat_low,
                                  t.lat_upp, t.lon_low, t.lon_upp);
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 7, "Invalid state!")

            auto out = std::make_unique<Surf::Data>();
            out->data = t[0].cast<Surf::FieldData>();
            out->lat_low = t[3].cast<Surf::Extrapolation>();
            out->lat_upp = t[4].cast<Surf::Extrapolation>();
            out->lon_low = t[5].cast<Surf::Extrapolation>();
            out->lon_upp = t[6].cast<Surf::Extrapolation>();
            return out;
          }));
  py::implicitly_convertible<GriddedField3, Surf::Data>();
  py::implicitly_convertible<Tensor3, Surf::Data>();
  py::implicitly_convertible<Numeric, Surf::Data>();
  py::implicitly_convertible<Index, Surf::Data>();
  py::implicitly_convertible<Surf::FunctionalData, Surf::Data>();

  auto pnt = py::class_<SurfPoint>(m, "SurfPoint").def(py::init([] {
    return std::make_unique<SurfPoint>();
  }));

  auto fld = py::class_<SurfField>(m, "SurfField").def(py::init([] {
    return std::make_unique<SurfField>();
  }));

  pnt.def_readwrite("temperature", &SurfPoint::temperature)
      .def_readwrite("altitude", &SurfPoint::altitude)
      .def_readwrite("normal", &SurfPoint::normal)
      .def_readwrite("wind", &SurfPoint::wind)
      .def(
          "__getitem__",
          [](SurfPoint &surf, Surf::Key x) {
            if (not surf.has(x))
              throw py::key_error(var_string(x));
            return surf[x];
          },
          py::return_value_policy::reference_internal)
      .def("__setitem__",
           [](SurfPoint &surf, Surf::Key x, Numeric data) { surf[x] = data; })
      .PythonInterfaceCopyValue(SurfPoint)
      .PythonInterfaceWorkspaceVariableConversion(SurfPoint)
      .PythonInterfaceFileIO(SurfPoint)
      .PythonInterfaceBasicRepresentation(SurfPoint)
      .def(py::pickle(
          [](const SurfPoint &t) {
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

            auto k = t[0].cast<std::vector<Surf::KeyVal>>();
            auto v = t[1].cast<std::vector<Numeric>>();
            ARTS_USER_ERROR_IF(k.size() != v.size(), "Invalid state!")

            auto out = std::make_unique<SurfPoint>();
            for (std::size_t i = 0; i < k.size(); i++)
              std::visit(
                  [&](auto &&key) -> Numeric & { return out->operator[](key); },
                  k[i]) = v[i];

            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(SurfPoint);

  fld.def(
         "__getitem__",
         [](SurfField &surf, Surf::Key x) -> Surf::Data & {
           if (not surf.has(x))
             throw py::key_error(var_string(x));
           return surf[x];
         },
         py::return_value_policy::reference_internal)
      .def("__setitem__", [](SurfField &surf, Surf::Key x,
                             const Surf::Data &data) { surf[x] = data; })
      .def("at", [](const SurfField &surf, Numeric lat, Numeric lon, Vector2 ell) { return surf.at(lat, lon, ell); })
      .PythonInterfaceCopyValue(SurfField)
      .PythonInterfaceWorkspaceVariableConversion(SurfField)
      .PythonInterfaceFileIO(SurfField)
      .PythonInterfaceBasicRepresentation(SurfField)
      .def(py::pickle(
          [](const SurfField &t) {
            auto k = t.keys();
            std::vector<Surf::Data> v;
            v.reserve(k.size());
            for (auto &kn : k)
              v.emplace_back(
                  std::visit([&](auto &&key) { return t[key]; }, kn));
            return py::make_tuple(k, v);
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 2, "Invalid state!")

            auto k = t[0].cast<std::vector<Surf::KeyVal>>();
            auto v = t[1].cast<std::vector<Surf::Data>>();
            ARTS_USER_ERROR_IF(k.size() != v.size(), "Invalid state!")

            auto out = std::make_unique<SurfField>();

            for (std::size_t i = 0; i < k.size(); i++)
              std::visit(
                  [&](auto &&key) -> Surf::Data & {
                    return out->operator[](key);
                  },
                  k[i]) = v[i];

            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(SurfField);
}
} // namespace Python
