
#include <debug.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>
#include <quantum_numbers.h>
#include <species_tags.h>
#include <surf.h>

#include "hpy_arts.h"
namespace Python {
void py_surf(py::module_ &m) try {
  py::class_<SurfacePropertyTag> spt =
      py::class_<SurfacePropertyTag>(m, "SurfacePropertyTag");
  workspace_group_interface(spt);
  spt.def_rw("name", &SurfacePropertyTag::name, "Name of property");
  spt.def(py::init_implicit<String>());

  py::class_<SurfaceTypeTag> stt =
      py::class_<SurfaceTypeTag>(m, "SurfaceTypeTag");
  workspace_group_interface(stt);
  stt.def_rw("name", &SurfaceTypeTag::name, "Name of type");
  stt.def(py::init_implicit<String>());

  py::class_<Surf::Data>(m, "SurfData")
      .def(py::init<>())
      .def(py::init_implicit<GriddedField2>())
      .def(py::init_implicit<Numeric>())
      .def(py::init_implicit<Surf::FunctionalData>())
      .def_rw("data", &Surf::Data::data, "Data")
      .def_rw("lat_upp", &Surf::Data::lat_upp, "Latitude upper bound")
      .def_rw("lat_low", &Surf::Data::lat_low, "Latitude lower bound")
      .def_rw("lon_upp", &Surf::Data::lon_upp, "Longitude upper bound")
      .def_rw("lon_low", &Surf::Data::lon_low, "Longitude lower bound")
      .def("__getstate__",
           [](const Surf::Data &t) {
             return std::tuple<Surf::FieldData,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation>(
                 t.data, t.lat_low, t.lat_upp, t.lon_low, t.lon_upp);
           })
      .def("__setstate__",
           [](Surf::Data *s,
              const std::tuple<Surf::FieldData,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation> &state) {
             new (s) Surf::Data{};
             s->data    = std::get<0>(state);
             s->lat_low = std::get<1>(state);
             s->lat_upp = std::get<2>(state);
             s->lon_low = std::get<3>(state);
             s->lon_upp = std::get<4>(state);
           })
      .doc() = "Surface data";
  py::implicitly_convertible<GriddedField2, Surf::Data>();
  py::implicitly_convertible<Numeric, Surf::Data>();
  py::implicitly_convertible<Index, Surf::Data>();
  py::implicitly_convertible<Surf::FunctionalData::func_t, Surf::Data>();

  auto pnt = py::class_<SurfacePoint>(m, "SurfacePoint");
  workspace_group_interface(pnt);

  auto fld = py::class_<SurfaceField>(m, "SurfaceField");
  workspace_group_interface(fld);

  pnt.def_rw("temperature", &SurfacePoint::temperature, "Temperature [K]")
      .def_rw("elevation", &SurfacePoint::elevation, "Surface elevation [m]")
      .def_rw("normal", &SurfacePoint::normal, "Surface normal vector")
      .def_rw("wind", &SurfacePoint::wind, "Wind vector [m/s]")
      .def(
          "__getitem__",
          [](SurfacePoint &surf, SurfaceKey x) {
            if (not surf.has(x)) throw py::key_error(std::format("{}", x).c_str());
            return surf[x];
          },
          py::rv_policy::reference_internal)
      .def("__setitem__",
           [](SurfacePoint &surf, SurfaceKey x, Numeric data) {
             surf[x] = data;
           })
      .def(
          "__getstate__",
          [](const SurfacePoint &t) {
            auto k = t.keys();
            std::vector<Numeric> v;
            v.reserve(k.size());
            for (auto &kn : k) {
              v.emplace_back(
                  std::visit([&](auto &&key) { return t[key]; }, kn));
            }
            return std::tuple<std::vector<SurfaceKeyVal>, std::vector<Numeric>>(
                k, v);
          })
      .def("__setstate__",
           [](SurfacePoint *s,
              const std::tuple<std::vector<SurfaceKeyVal>, std::vector<Numeric>>
                  &state) {
             auto k = std::get<0>(state);
             auto v = std::get<1>(state);
             ARTS_USER_ERROR_IF(k.size() != v.size(), "Invalid state!")

             new (s) SurfacePoint{};
             for (std::size_t i = 0; i < k.size(); i++) {
               std::visit(
                   [&](auto &&key) -> Numeric & { return s->operator[](key); },
                   k[i]) = v[i];
             }
           });

  fld.def(
         "__getitem__",
         [](SurfaceField &surf, SurfaceKey x) -> Surf::Data & {
           if (not surf.has(x)) throw py::key_error(std::format("{}", x).c_str());
           return surf[x];
         },
         py::rv_policy::reference_internal)
      .def("__setitem__",
           [](SurfaceField &surf, SurfaceKey x, const Surf::Data &data) {
             surf[x] = data;
           })
      .def("at",
           [](const SurfaceField &surf, Numeric lat, Numeric lon) {
             return surf.at(lat, lon);
           }, "lat"_a, "lon"_a, "Get the data at a point")
      .def("__getstate__",
           [](const SurfaceField &t) {
             auto k = t.keys();
             std::vector<Surf::Data> v;
             v.reserve(k.size());
             for (auto &kn : k) {
               v.emplace_back(
                   std::visit([&](auto &&key) { return t[key]; }, kn));
             }
             return std::tuple<std::vector<SurfaceKeyVal>,
                               std::vector<Surf::Data>>(k, v);
           })
      .def("__setstate__",
           [](SurfaceField *s,
              const std::tuple<std::vector<SurfaceKeyVal>,
                               std::vector<Surf::Data>> &state) {
             auto k = std::get<0>(state);
             auto v = std::get<1>(state);
             ARTS_USER_ERROR_IF(k.size() != v.size(), "Invalid state!")

             new (s) SurfaceField{};

             for (std::size_t i = 0; i < k.size(); i++) {
               std::visit(
                   [&](auto &&key) -> Surf::Data & {
                     return s->operator[](key);
                   },
                   k[i]) = v[i];
             }
           });
} catch (std::exception &e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize surf\n{}", e.what()));
}
}  // namespace Python
