
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

  auto pnt = py::class_<SurfacePoint>(m, "SurfacePoint").def(py::init([] {
    return std::make_unique<SurfacePoint>();
  }));

  auto fld = py::class_<SurfaceField>(m, "SurfaceField").def(py::init([] {
    return std::make_unique<SurfaceField>();
  }));

  pnt.def_readwrite("temperature", &SurfacePoint::temperature)
      .def_readwrite("elevation", &SurfacePoint::elevation)
      .def_readwrite("normal", &SurfacePoint::normal)
      .def_readwrite("wind", &SurfacePoint::wind)
      .def(
          "__getitem__",
          [](SurfacePoint &surf, Surf::Key x) {
            if (not surf.has(x))
              throw py::key_error(var_string(x));
            return surf[x];
          },
          py::return_value_policy::reference_internal)
      .def("__setitem__",
           [](SurfacePoint &surf, Surf::Key x, Numeric data) { surf[x] = data; })
      .PythonInterfaceCopyValue(SurfacePoint)
      .PythonInterfaceWorkspaceVariableConversion(SurfacePoint)
      .PythonInterfaceFileIO(SurfacePoint)
      .PythonInterfaceBasicRepresentation(SurfacePoint)
      .def(py::pickle(
          [](const SurfacePoint &t) {
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

            auto out = std::make_unique<SurfacePoint>();
            for (std::size_t i = 0; i < k.size(); i++)
              std::visit(
                  [&](auto &&key) -> Numeric & { return out->operator[](key); },
                  k[i]) = v[i];

            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(SurfacePoint);

  fld.def(
         "__getitem__",
         [](SurfaceField &surf, Surf::Key x) -> Surf::Data & {
           if (not surf.has(x))
             throw py::key_error(var_string(x));
           return surf[x];
         },
         py::return_value_policy::reference_internal)
      .def("__setitem__", [](SurfaceField &surf, Surf::Key x,
                             const Surf::Data &data) { surf[x] = data; })
      .def("at", [](const SurfaceField &surf, Numeric lat, Numeric lon, Vector2 ell) { return surf.at(lat, lon, ell); })
      .PythonInterfaceCopyValue(SurfaceField)
      .PythonInterfaceWorkspaceVariableConversion(SurfaceField)
      .PythonInterfaceFileIO(SurfaceField)
      .PythonInterfaceBasicRepresentation(SurfaceField)
      .def(py::pickle(
          [](const SurfaceField &t) {
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

            auto out = std::make_unique<SurfaceField>();

            for (std::size_t i = 0; i < k.size(); i++)
              std::visit(
                  [&](auto &&key) -> Surf::Data & {
                    return out->operator[](key);
                  },
                  k[i]) = v[i];

            return out;
          }))
      .PythonInterfaceWorkspaceDocumentation(SurfaceField);
}
} // namespace Python
