
#include <debug.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>
#include <species_tags.h>
#include <surf.h>

#include "hpy_arts.h"

namespace Python {
void py_surf(py::module_ &m) try {
  py::class_<Surf::Data> surfdata(m, "SurfaceData");
  surfdata.def(py::init_implicit<SortedGriddedField2>())
      .def(py::init_implicit<Numeric>())
      .def(py::init_implicit<Surf::FunctionalData>())
      .def(
          "__init__",
          [](Surf::Data *a, const GriddedField2 &v) {
            new (a) Surf::Data(SortedGriddedField2(v));
          },
          "v"_a,
          "Initialize with a sorted field")
      .def_rw("data", &Surf::Data::data, "The data")
      .def_rw("lat_upp", &Surf::Data::lat_upp, "Upper latitude limit")
      .def_rw("lat_low", &Surf::Data::lat_low, "Lower latitude limit")
      .def_rw("lon_upp", &Surf::Data::lon_upp, "Upper longitude limit")
      .def_rw("lon_low", &Surf::Data::lon_low, "Lower longitude limit")
      .def(
          "set_extrapolation",
          [](Surf::Data &self, InterpolationExtrapolation x) {
            self.lat_upp = x;
            self.lat_low = x;
            self.lon_upp = x;
            self.lon_low = x;
          },
          "extrapolation"_a,
          "Set the extrapolation for all dimensions")
      .def(
          "__call__",
          [](const Surf::Data &d, Numeric lat, Numeric lon) {
            return d.at(lat, lon);
          },
          "lat"_a,
          "lon"_a,
          "Get a point of data at the position")
      .def(
          "__call__",
          [](const SurfaceData &surf, const Vector &latv, const Vector &lonv) {
            const Size N = latv.size();
            if (N != lonv.size())
              throw std::logic_error(std::format(R"(Not same size:
  lat: {:B,} (size: {})
  lon: {:B,} (size: {})
)",
                                                 latv,
                                                 latv.size(),
                                                 lonv,
                                                 lonv.size()));

            Vector out(N);
            for (Size i = 0; i < N; i++) {
              out[i] = surf.at(latv[i], lonv[i]);
            }
            return out;
          },
          "lat"_a,
          "lon"_a,
          "Get the data at a list of points")
      .def(
          "ws",
          [](const Surf::Data &d, Numeric lat, Numeric lon) {
            return d.flat_weights(lat, lon);
          },
          "lat"_a,
          "lon"_a,
          "Get the weights of neighbors at a position")
      .def_prop_ro("data_type", &Surf::Data::data_type, "The data type")
      .def("__getstate__",
           [](const Surf::Data &t) {
             return std::make_tuple(
                 t.data, t.lat_low, t.lat_upp, t.lon_low, t.lon_upp);
           })
      .def("__setstate__",
           [](Surf::Data *a,
              const std::tuple<Surf::FieldData,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation> &state) {
             new (a) Surf::Data();
             a->data    = std::get<0>(state);
             a->lat_low = std::get<1>(state);
             a->lat_upp = std::get<2>(state);
             a->lon_low = std::get<3>(state);
             a->lon_upp = std::get<4>(state);
           });
  surfdata.doc() = "Surface data";
  py::implicitly_convertible<Surf::FunctionalData::func_t, Surf::Data>();
  py::implicitly_convertible<GriddedField2, Surf::Data>();

  py::class_<SurfacePropertyTag> spt =
      py::class_<SurfacePropertyTag>(m, "SurfacePropertyTag");
  generic_interface(spt);
  spt.def_rw("name", &SurfacePropertyTag::name, "Name of property");
  spt.def(py::init_implicit<String>());

  auto pnt = py::class_<SurfacePoint>(m, "SurfacePoint");
  generic_interface(pnt);

  auto asp =
      py::bind_vector<ArrayOfSurfacePoint, py::rv_policy::reference_internal>(
          m, "ArrayOfSurfacePoint");
  asp.doc() = "Array of SurfacePoint";

  auto fld = py::class_<SurfaceField>(m, "SurfaceField");
  fld.def(
         "__init__",
         [](SurfaceField *sf, const String &planet) {
           new (sf) SurfaceField();
           surface_fieldPlanet(*sf, planet, 0.0);
         },
         "planet"_a)
      .def_rw("ellipsoid",
              &SurfaceField::ellipsoid,
              "Ellipsoid parameters (semi-major axis, semi-minor axis)");
  generic_interface(fld);
  py::implicitly_convertible<String, SurfaceField>();

  pnt.def_rw("temperature", &SurfacePoint::temperature, "Temperature [K]")
      .def_rw("elevation", &SurfacePoint::elevation, "Surface elevation [m]")
      .def_rw("normal", &SurfacePoint::normal, "Surface normal vector")
      .def(
          "__getitem__",
          [](SurfacePoint &surf, SurfaceKey x) {
            if (not surf.has(x)) {
              const auto error_message = std::format("{}", x);
              throw py::key_error(error_message.c_str());
            }
            return surf[x];
          },
          py::rv_policy::reference_internal)
      .def("__setitem__",
           [](SurfacePoint &surf, SurfaceKey x, Numeric data) {
             surf[x] = data;
           })
      .def(
          "__getitem__",
          [](SurfacePoint &surf, const SurfacePropertyTag &x) {
            if (not surf.has(x)) {
              const auto error_message = std::format("{}", x);
              throw py::key_error(error_message.c_str());
            }
            return surf[x];
          },
          py::rv_policy::reference_internal)
      .def("__setitem__",
           [](SurfacePoint &surf, const SurfacePropertyTag &x, Numeric data) {
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
           if (not surf.contains(x)) {
             const auto error_message = std::format("{}", x);
             throw py::key_error(error_message.c_str());
           }
           return surf[x];
         },
         py::rv_policy::reference_internal)
      .def("__setitem__",
           [](SurfaceField &surf, SurfaceKey x, const Surf::Data &data) {
             surf[x] = data;
           })
      .def(
          "__getitem__",
          [](SurfaceField &surf, const SurfacePropertyTag &x) -> Surf::Data & {
            if (not surf.contains(x)) {
              const auto error_message = std::format("{}", x);
              throw py::key_error(error_message.c_str());
            }
            return surf[x];
          },
          py::rv_policy::reference_internal)
      .def("__setitem__",
           [](SurfaceField &surf,
              const SurfacePropertyTag &x,
              const Surf::Data &data) { surf[x] = data; })
      .def(
          "__call__",
          [](const SurfaceField &surf, Numeric lat, Numeric lon) {
            return surf.at(lat, lon);
          },
          "lat"_a,
          "lon"_a,
          "Get the data at a point")
      .def(
          "__call__",
          [](const SurfaceField &surf, const Vector &latv, const Vector &lonv) {
            const Size N = latv.size();
            if (N != lonv.size())
              throw std::logic_error(std::format(R"(Not same size:
  lat: {:B,} (size: {})
  lon: {:B,} (size: {})
)",
                                                 latv,
                                                 latv.size(),
                                                 lonv,
                                                 lonv.size()));

            Array<SurfacePoint> out(N);
            for (Size i = 0; i < N; i++) {
              out[i] = surf.at(latv[i], lonv[i]);
            }
            return out;
          },
          "lat"_a,
          "lon"_a,
          "Get the data at a list of points")
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
  fld.def_rw("other", &SurfaceField::other, "Other data in the surface field")
      .def_rw("props", &SurfaceField::props, "Properties of the surface field");

  py::class_<SubsurfaceField> ssf(m, "SubsurfaceField");
  ssf.def_rw(
      "other", &SubsurfaceField::other, "Other data in the subsurface field");
  ssf.def_rw("bottom_depth",
             &SubsurfaceField::bottom_depth,
             "The depth of the bottom of the subsurface [m]");
  generic_interface(ssf);

  py::class_<SubsurfacePoint> ssp(m, "SubsurfacePoint");
  generic_interface(ssp);
  auto assp = py::bind_vector<Array<SubsurfacePoint>,
                              py::rv_policy::reference_internal>(
      m, "ArrayOfSubsurfacePoint");
  generic_interface(assp);
  vector_interface(assp);
} catch (std::exception &e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize surf\n{}", e.what()));
}
}  // namespace Python
