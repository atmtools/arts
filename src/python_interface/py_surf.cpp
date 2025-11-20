
#include <debug.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <python_interface.h>
#include <species_tags.h>
#include <subsurf_field.h>
#include <surf.h>

#include "hpy_arts.h"

namespace Python {
void py_surf(py::module_ &m) try {
  py::class_<Surf::Data> surfdata(m, "SurfaceData");
  surfdata.def(py::init_implicit<GeodeticField2>())
      .def(py::init_implicit<Numeric>())
      .def(py::init_implicit<Surf::FunctionalData>())
      .def(
          "__init__",
          [](Surf::Data *a, const GriddedField2 &v) {
            new (a) Surf::Data(GeodeticField2(v));
          },
          "v"_a,
          "Initialize with a sorted field")
      .def_rw(
          "data",
          &Surf::Data::data,
          "The data\n\n.. :class:`GeodeticField2`\n\n.. :class:`Numeric`\n\n.. :class:`NumericBinaryOperator`")
      .def_rw("lat_upp",
              &Surf::Data::lat_upp,
              "Upper latitude limit\n\n.. :class:`InterpolationExtrapolation`")
      .def_rw("lat_low",
              &Surf::Data::lat_low,
              "Lower latitude limit\n\n.. :class:`InterpolationExtrapolation`")
      .def_rw("lon_upp",
              &Surf::Data::lon_upp,
              "Upper longitude limit\n\n.. :class:`InterpolationExtrapolation`")
      .def_rw("lon_low",
              &Surf::Data::lon_low,
              "Lower longitude limit\n\n.. :class:`InterpolationExtrapolation`")
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
      .def_prop_ro("data_type",
                   &Surf::Data::data_type,
                   "The data type\n\n.. :class:`String`")
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
  generic_interface(surfdata);

  auto aosd =
      py::bind_vector<Array<SurfaceData>, py::rv_policy::reference_internal>(
          m, "ArrayOfSurfaceData");
  aosd.doc() = "A list of surface data";
  vector_interface(aosd);
  generic_interface(aosd);

  py::class_<SurfacePropertyTag> spt(m, "SurfacePropertyTag");
  generic_interface(spt);
  spt.def_rw("name",
             &SurfacePropertyTag::name,
             "Name of property\n\n.. :class:`String`");
  spt.def(py::init_implicit<String>());

  auto pnt = py::class_<SurfacePoint>(m, "SurfacePoint");
  generic_interface(pnt);

  auto asp =
      py::bind_vector<ArrayOfSurfacePoint, py::rv_policy::reference_internal>(
          m, "ArrayOfSurfacePoint");
  asp.doc() = "Array of SurfacePoint";
  generic_interface(asp);
  vector_interface(asp);

  auto fld = py::class_<SurfaceField>(m, "SurfaceField");
  fld.def(
         "__init__",
         [](SurfaceField *sf, const String &planet) {
           new (sf) SurfaceField();
           surf_fieldPlanet(*sf, planet, 0.0);
         },
         "planet"_a)
      .def_rw(
          "ellipsoid",
          &SurfaceField::ellipsoid,
          "Ellipsoid parameters (semi-major axis, semi-minor axis)\n\n.. :class:`Vector2`");
  fld.def("keys", &SurfaceField::keys, "Available keys");
  generic_interface(fld);
  py::implicitly_convertible<String, SurfaceField>();

  pnt.def_rw("temperature",
             &SurfacePoint::temperature,
             "Temperature [K]\n\n.. :class:`Numeric`")
      .def_rw("elevation",
              &SurfacePoint::elevation,
              "Surface elevation [m]\n\n.. :class:`Numeric`")
      .def_rw("normal",
              &SurfacePoint::normal,
              "Surface normal vector\n\n.. :class:`Vector2`")
      .def(
          "__getitem__",
          [](SurfacePoint &surf, const SurfaceKeyVal &x) {
            if (not surf.contains(x)) {
              const auto error_message = std::format("{}", x);
              throw py::key_error(error_message.c_str());
            }
            return surf[x];
          },
          py::rv_policy::reference_internal)
      .def("__setitem__",
           [](SurfacePoint &surf, const SurfaceKeyVal &x, Numeric data) {
             surf[x] = data;
           })
      .def("__contains__",
           [](const SurfacePoint &self, const SurfaceKeyVal &x) {
             return self.contains(x);
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
           })
      .def("keys", &SurfacePoint::keys, "Available keys");

  fld.def(
         "__getitem__",
         [](SurfaceField &surf, const SurfaceKeyVal &x) -> Surf::Data & {
           if (not surf.contains(x)) {
             const auto error_message = std::format("{}", x);
             throw py::key_error(error_message.c_str());
           }
           return surf[x];
         },
         py::rv_policy::reference_internal)
      .def("__setitem__",
           [](SurfaceField &surf,
              const SurfaceKeyVal &x,
              const Surf::Data &data) { surf[x] = data; })
      .def("__contains__",
           [](const SurfaceField &self, const SurfaceKeyVal &x) {
             return self.contains(x);
           })
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
  fld.def_rw(
         "other",
         &SurfaceField::other,
         "Other data in the surface field\n\n.. :class:`dict[SurfaceKey, SurfaceData]`")
      .def_rw(
          "props",
          &SurfaceField::props,
          "Properties of the surface field\n\n.. :class:`dict[SurfacePropertyTag, SurfaceData]`");

  py::class_<SubsurfacePropertyTag> sptag(m, "SubsurfacePropertyTag");
  sptag.def_rw("name",
               &SubsurfacePropertyTag::name,
               "Name of the subsurface property\n\n.. :class:`String`");
  sptag.def(py::init_implicit<String>());
  generic_interface(sptag);

  py::class_<SubsurfaceField> ssf(m, "SubsurfaceField");
  ssf.def(
      "__init__",
      [](SubsurfaceField *sf, Numeric bottom_depth) {
        new (sf) SubsurfaceField();
        sf->bottom_depth = bottom_depth;
      },
      "bottom_depth"_a);
  ssf.def_rw(
      "other",
      &SubsurfaceField::other,
      "Other data in the subsurface field\n\n.. :class:`dict[SubsurfaceKey, SubsurfaceData]`");
  ssf.def_rw(
      "prop",
      &SubsurfaceField::prop,
      "Properties of the subsurface field\n\n.. :class:`dict[SubsurfacePropertyTag, SubsurfaceData]`");
  ssf.def_rw(
      "bottom_depth",
      &SubsurfaceField::bottom_depth,
      "The depth of the bottom of the subsurface [m]\n\n.. :class:`Numeric`");
  ssf.def(
      "__call__",
      [](const SubsurfaceField &d, Numeric alt, Numeric lat, Numeric lon) {
        return d.at(alt, lat, lon);
      },
      "alt"_a,
      "lat"_a,
      "lon"_a,
      "Get a point of data at the position");
  ssf.def(
      "__call__",
      [](const SubsurfaceField &atm,
         const Vector &hv,
         const Vector &latv,
         const Vector &lonv) {
        const Size N = hv.size();
        if (latv.size() != lonv.size() or N != latv.size())
          throw std::logic_error(std::format(R"(Not same size:
  h:   {:B,} (size: {})
  lat: {:B,} (size: {})
  lon: {:B,} (size: {})
)",
                                             hv,
                                             hv.size(),
                                             latv,
                                             latv.size(),
                                             lonv,
                                             lonv.size()));
        ArrayOfSubsurfacePoint out;
        out.reserve(N);
        for (Size i = 0; i < N; i++)
          out.emplace_back(atm.at(hv[i], latv[i], lonv[i]));
        return out;
      },
      "h"_a,
      "lat"_a,
      "lon"_a,
      "Get the data as a list");
  ssf.def("__setitem__",
          [](SubsurfaceField &surf,
             const SubsurfaceKeyVal &x,
             const SubsurfaceData &data) { surf[x] = data; })
      .def(
          "__getitem__",
          [](SubsurfaceField &surf,
             const SubsurfaceKeyVal &x) -> SubsurfaceData & {
            if (not surf.contains(x)) {
              const auto error_message = std::format("{}", x);
              throw py::key_error(error_message.c_str());
            }
            return surf[x];
          },
          py::rv_policy::reference_internal)
      .def("__contains__",
           [](const SubsurfaceField &surf, const SubsurfaceKeyVal &x) {
             return surf.contains(x);
           });
  ssf.def("keys", &SubsurfaceField::keys, "Available keys");
  generic_interface(ssf);

  py::class_<SubsurfacePoint> ssp(m, "SubsurfacePoint");
  ssp.def_rw("temperature",
             &SubsurfacePoint::temperature,
             "Temperature [K]\n\n.. :class:`Numeric`");
  ssp.def_rw("density",
             &SubsurfacePoint::density,
             "Density [kg/m^3]\n\n.. :class:`Numeric`");
  ssp.def_rw(
      "prop",
      &SubsurfacePoint::prop,
      "Properties of the subsurface point\n\n.. :class:`dict[SubsurfacePropertyTag, Numeric]`");
  ssp.def("__getitem__",
          [](SubsurfacePoint &self, const SubsurfaceKeyVal &key) {
            if (self.contains(key)) return self[key];
            throw py::key_error(std::format("{}", key).c_str());
          });
  ssp.def("__setitem__",
          [](SubsurfacePoint &self, const SubsurfaceKeyVal &key, Numeric x) {
            self[key] = x;
          });
  ssp.def("__contains__",
          [](const SubsurfacePoint &self, const SubsurfaceKeyVal &x) {
            return self.contains(x);
          });
  ssp.def("keys", &SubsurfacePoint::keys, "Available keys");
  generic_interface(ssp);

  py::class_<SubsurfaceData> ssd(m, "SubsurfaceData");
  ssd.def(py::init_implicit<GeodeticField3>())
      .def(py::init_implicit<Numeric>())
      .def(py::init_implicit<Subsurface::FunctionalData>())
      .def(
          "__init__",
          [](Subsurface::Data *a, const GriddedField3 &v) {
            new (a) Subsurface::Data(GeodeticField3(v));
          },
          "v"_a,
          "Initialize with a gridded field")
      .def_rw(
          "data",
          &Subsurface::Data::data,
          "The data.\n\n.. :class:`~pyarts3.arts.GeodeticField3`\n\n.. :class:`~pyarts3.arts.Numeric`\n\n.. :class:`~pyarts3.arts.NumericTernaryOperator`")
      .def_rw(
          "alt_upp",
          &Subsurface::Data::alt_upp,
          "Upper altitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw(
          "alt_low",
          &Subsurface::Data::alt_low,
          "Lower altitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw(
          "lat_upp",
          &Subsurface::Data::lat_upp,
          "Upper latitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw(
          "lat_low",
          &Subsurface::Data::lat_low,
          "Lower latitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw(
          "lon_upp",
          &Subsurface::Data::lon_upp,
          "Upper longitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw(
          "lon_low",
          &Subsurface::Data::lon_low,
          "Lower longitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def(
          "set_extrapolation",
          [](Subsurface::Data &self, InterpolationExtrapolation x) {
            self.alt_upp = x;
            self.alt_low = x;
            self.lat_upp = x;
            self.lat_low = x;
            self.lon_upp = x;
            self.lon_low = x;
          },
          "extrapolation"_a,
          "Set the extrapolation for all dimensions")
      .def(
          "__call__",
          [](const Subsurface::Data &d, Numeric alt, Numeric lat, Numeric lon) {
            return d.at(alt, lat, lon);
          },
          "alt"_a,
          "lat"_a,
          "lon"_a,
          "Get a point of data at the position")
      .def(
          "ws",
          [](const Subsurface::Data &d, Numeric alt, Numeric lat, Numeric lon) {
            return d.flat_weight(alt, lat, lon);
          },
          "alt"_a,
          "lat"_a,
          "lon"_a,
          "Get the weights of neighbors at a position")
      .def_prop_ro("data_type",
                   &Subsurface::Data::data_type,
                   "The data type\n\n.. :class:`~pyarts3.arts.String`")
      .def("__getstate__",
           [](const Subsurface::Data &t) {
             return std::make_tuple(t.data,
                                    t.alt_low,
                                    t.alt_upp,
                                    t.lat_low,
                                    t.lat_upp,
                                    t.lon_low,
                                    t.lon_upp);
           })
      .def("__setstate__",
           [](Subsurface::Data *a,
              const std::tuple<Subsurface::FieldData,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation,
                               InterpolationExtrapolation> &state) {
             new (a) Subsurface::Data();
             a->data    = std::get<0>(state);
             a->alt_low = std::get<1>(state);
             a->alt_upp = std::get<2>(state);
             a->lat_low = std::get<3>(state);
             a->lat_upp = std::get<4>(state);
             a->lon_low = std::get<5>(state);
             a->lon_upp = std::get<6>(state);
           });
  py::implicitly_convertible<Subsurface::FunctionalData::func_t,
                             Subsurface::Data>();
  py::implicitly_convertible<GriddedField3, Subsurface::Data>();
  generic_interface(ssd);

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
