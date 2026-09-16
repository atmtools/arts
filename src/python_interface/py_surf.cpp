
#include <debug.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <python_interface.h>
#include <species_tags.h>
#include <subsurf_field.h>
#include <surf.h>

#include "hpy_arts.h"

namespace Python {
void py_surf(py::module_ &m) try {
  auto tessem = py::class_<TessemNN>(m, "TessemNN");
  tessem.def(py::init<>())
      .def_rw("nb_inputs", &TessemNN::nb_inputs, "Number of neural-network inputs\n\n.. :class:`~pyarts3.arts.Index`")
      .def_rw("nb_outputs", &TessemNN::nb_outputs, "Number of neural-network outputs\n\n.. :class:`~pyarts3.arts.Index`")
      .def_rw("nb_cache", &TessemNN::nb_cache, "Number of hidden neural-network nodes\n\n.. :class:`~pyarts3.arts.Index`")
      .def_rw("b1", &TessemNN::b1, "Hidden-layer biases\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("b2", &TessemNN::b2, "Output-layer biases\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("w1", &TessemNN::w1, "Hidden-layer weights\n\n.. :class:`~pyarts3.arts.Matrix`")
      .def_rw("w2", &TessemNN::w2, "Output-layer weights\n\n.. :class:`~pyarts3.arts.Matrix`")
      .def_rw("x_min", &TessemNN::x_min, "Minimum values used to scale the inputs\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("x_max", &TessemNN::x_max, "Maximum values used to scale the inputs\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("y_min", &TessemNN::y_min, "Minimum values used to scale the outputs\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("y_max", &TessemNN::y_max, "Maximum values used to scale the outputs\n\n.. :class:`~pyarts3.arts.Vector`")
      .def(
          "__call__",
          [](const TessemNN &self, const Vector &input) { return tessem_emissivity(self, input); },
          "input"_a,
          "Evaluate the neural network")
      .def_static(
          "from_ascii",
          [](const String &filename) {
            TessemNN out;
            tessem_read_ascii(filename, out);
            return out;
          },
          "filename"_a,
          "Read an original TESSEM2 neural-network parameter file");
  generic_interface(tessem);

  auto telsem = py::class_<TelsemAtlas>(m, "TelsemAtlas");
  telsem.def(py::init<>())
      .def_ro("ndat", &TelsemAtlas::ndat, "Number of populated atlas cells\n\n.. :class:`~pyarts3.arts.Index`")
      .def_ro_static("nchan", &TelsemAtlas::nchan, "Number of atlas channels\n\n.. :class:`~pyarts3.arts.Index`")
      .def_ro("name", &TelsemAtlas::name, "Atlas name\n\n.. :class:`~pyarts3.arts.String`")
      .def_ro("month", &TelsemAtlas::month, "Atlas month\n\n.. :class:`~pyarts3.arts.Index`")
      .def_ro("dlat", &TelsemAtlas::dlat, "Atlas latitude resolution [degrees]\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def("contains", &TelsemAtlas::contains, "cell_number"_a, "Check whether an atlas cell contains data")
      .def("cell_number", &TelsemAtlas::calc_cellnum, "lat"_a, "lon"_a, "Return the atlas cell number at a position")
      .def("coordinates", &TelsemAtlas::get_coordinates, "cell_number"_a, "Return the coordinates of an atlas cell")
      .def(
          "emissivity",
          [](const TelsemAtlas &self,
             Numeric            lat,
             Numeric            lon,
             Numeric            incidence_angle,
             Numeric            frequency,
             Numeric            max_distance) {
            const auto out = self.emissivity(lat, lon, incidence_angle, frequency, max_distance);
            return std::pair{out[0], out[1]};
          },
          "lat"_a,
          "lon"_a,
          "incidence_angle"_a,
          "frequency"_a,
          "max_distance"_a = -1,
          "Evaluate vertical and horizontal emissivity")
      .def_static(
          "from_ascii",
          [](const String &filename, Index month) {
            TelsemAtlas out;
            telsem_read_ascii(filename, out, month);
            return out;
          },
          "filename"_a,
          "month"_a = 0,
          "Read an original TELSEM2 monthly atlas file");
  generic_interface(telsem);

  py::class_<Surf::Data> surfdata(m, "SurfaceData");
  surfdata.def(py::init_implicit<GeodeticField2>())
      .def(py::init_implicit<Numeric>())
      .def(py::init_implicit<Surf::FunctionalData>())
      .def(
          "__init__",
          [](Surf::Data *a, const GriddedField2 &v) { new (a) Surf::Data(GeodeticField2(v)); },
          "v"_a,
          "Initialize with a sorted field")
      .def_rw("data",
              &Surf::Data::data,
              "The data\n\n.. :class:`~pyarts3.arts.GeodeticField2`\n\n.. :class:`~pyarts3.arts.Numeric`\n\n.. :class:`~pyarts3.arts.NumericBinaryOperator`")
      .def_rw("lat_upp", &Surf::Data::lat_upp, "Upper latitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("lat_low", &Surf::Data::lat_low, "Lower latitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("lon_upp", &Surf::Data::lon_upp, "Upper longitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("lon_low", &Surf::Data::lon_low, "Lower longitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
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
          [](const Surf::Data &d, Numeric lat, Numeric lon) { return d.at(lat, lon); },
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
            for (Size i = 0; i < N; i++) { out[i] = surf.at(latv[i], lonv[i]); }
            return out;
          },
          "lat"_a,
          "lon"_a,
          "Get the data at a list of points")
      .def(
          "ws",
          [](const Surf::Data &d, Numeric lat, Numeric lon) { return d.flat_weights(lat, lon); },
          "lat"_a,
          "lon"_a,
          "Get the weights of neighbors at a position")
      .def_prop_ro("data_type", &Surf::Data::data_type, "The data type\n\n.. :class:`~pyarts3.arts.String`");
  surfdata.doc() = "Surface data";
  py::implicitly_convertible<Surf::FunctionalData::func_t, Surf::Data>();
  py::implicitly_convertible<GriddedField2, Surf::Data>();
  generic_interface(surfdata);

  auto aosd  = py::bind_vector<Array<SurfaceData>, py::rv_policy::reference_internal>(m, "ArrayOfSurfaceData");
  aosd.doc() = "A list of surface data";
  vector_interface(aosd);
  generic_interface(aosd);

  py::class_<SurfacePropertyTag> spt(m, "SurfacePropertyTag");
  generic_interface(spt);
  spt.def_rw("name", &SurfacePropertyTag::name, "Name of property\n\n.. :class:`~pyarts3.arts.String`");
  spt.def(py::init_implicit<String>());

  auto pnt = py::class_<SurfacePoint>(m, "SurfacePoint");
  generic_interface(pnt);

  auto asp  = py::bind_vector<ArrayOfSurfacePoint, py::rv_policy::reference_internal>(m, "ArrayOfSurfacePoint");
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
      .def_rw("ellipsoid",
              &SurfaceField::ellipsoid,
              "Ellipsoid parameters (semi-major axis, semi-minor axis)\n\n.. :class:`~pyarts3.arts.Vector2`");
  fld.def("keys", &SurfaceField::keys, "Available keys");
  fld.def("single_value", &SurfaceField::single_value, "key"_a, "lat"_a, "lon"_a, "Get a single value at a position");
  generic_interface(fld);
  py::implicitly_convertible<String, SurfaceField>();

  pnt.def_rw("temperature", &SurfacePoint::temperature, "Temperature [K]\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("elevation", &SurfacePoint::elevation, "Surface elevation [m]\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("normal", &SurfacePoint::normal, "Surface normal vector\n\n.. :class:`~pyarts3.arts.Vector2`")
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
      .def("__setitem__", [](SurfacePoint &surf, const SurfaceKeyVal &x, Numeric data) { surf[x] = data; })
      .def("__contains__", [](const SurfacePoint &self, const SurfaceKeyVal &x) { return self.contains(x); })

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
      .def("__setitem__", [](SurfaceField &surf, const SurfaceKeyVal &x, const Surf::Data &data) { surf[x] = data; })
      .def("__contains__", [](const SurfaceField &self, const SurfaceKeyVal &x) { return self.contains(x); })
      .def(
          "__call__",
          [](const SurfaceField &surf, Numeric lat, Numeric lon) { return surf.at(lat, lon); },
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
            for (Size i = 0; i < N; i++) { out[i] = surf.at(latv[i], lonv[i]); }
            return out;
          },
          "lat"_a,
          "lon"_a,
          "Get the data at a list of points");
  fld.def_rw(
         "other", &SurfaceField::other, "Other data in the surface field\n\n.. :class:`dict[SurfaceKey, SurfaceData]`")
      .def_rw("props",
              &SurfaceField::props,
              "Properties of the surface field\n\n.. :class:`dict[SurfacePropertyTag, SurfaceData]`");

  py::class_<SubsurfacePropertyTag> sptag(m, "SubsurfacePropertyTag");
  sptag.def_rw("name", &SubsurfacePropertyTag::name, "Name of the subsurface property\n\n.. :class:`~pyarts3.arts.String`");
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
  ssf.def_rw("other",
             &SubsurfaceField::other,
             "Other data in the subsurface field\n\n.. :class:`dict[SubsurfaceKey, SubsurfaceData]`");
  ssf.def_rw("props",
             &SubsurfaceField::props,
             "Properties of the subsurface field\n\n.. :class:`dict[SubsurfacePropertyTag, SubsurfaceData]`");
  ssf.def_rw("bottom_depth",
             &SubsurfaceField::bottom_depth,
             "The depth of the bottom of the subsurface [m]\n\n.. :class:`~pyarts3.arts.Numeric`");
  ssf.def(
      "__call__",
      [](const SubsurfaceField &d, Numeric alt, Numeric lat, Numeric lon) { return d.at(alt, lat, lon); },
      "alt"_a,
      "lat"_a,
      "lon"_a,
      "Get a point of data at the position");
  ssf.def(
      "__call__",
      [](const SubsurfaceField &atm, const Vector &hv, const Vector &latv, const Vector &lonv) {
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
        for (Size i = 0; i < N; i++) out.emplace_back(atm.at(hv[i], latv[i], lonv[i]));
        return out;
      },
      "h"_a,
      "lat"_a,
      "lon"_a,
      "Get the data as a list");
  ssf.def("__setitem__",
          [](SubsurfaceField &surf, const SubsurfaceKeyVal &x, const SubsurfaceData &data) { surf[x] = data; })
      .def(
          "__getitem__",
          [](SubsurfaceField &surf, const SubsurfaceKeyVal &x) -> SubsurfaceData & {
            if (not surf.contains(x)) {
              const auto error_message = std::format("{}", x);
              throw py::key_error(error_message.c_str());
            }
            return surf[x];
          },
          py::rv_policy::reference_internal)
      .def("__contains__", [](const SubsurfaceField &surf, const SubsurfaceKeyVal &x) { return surf.contains(x); });
  ssf.def("keys", &SubsurfaceField::keys, "Available keys");
  generic_interface(ssf);

  py::class_<SubsurfacePoint> ssp(m, "SubsurfacePoint");
  ssp.def_rw("temperature", &SubsurfacePoint::temperature, "Temperature [K]\n\n.. :class:`~pyarts3.arts.Numeric`");
  ssp.def_rw("density", &SubsurfacePoint::density, "Density [kg/m^3]\n\n.. :class:`~pyarts3.arts.Numeric`");
  ssp.def_rw("props",
             &SubsurfacePoint::props,
             "Properties of the subsurface point\n\n.. :class:`dict[SubsurfacePropertyTag, Numeric]`");
  ssp.def("__getitem__", [](SubsurfacePoint &self, const SubsurfaceKeyVal &key) {
    if (self.contains(key)) return self[key];
    throw py::key_error(std::format("{}", key).c_str());
  });
  ssp.def("__setitem__", [](SubsurfacePoint &self, const SubsurfaceKeyVal &key, Numeric x) { self[key] = x; });
  ssp.def("__contains__", [](const SubsurfacePoint &self, const SubsurfaceKeyVal &x) { return self.contains(x); });
  ssp.def("keys", &SubsurfacePoint::keys, "Available keys");
  generic_interface(ssp);

  py::class_<SubsurfaceData> ssd(m, "SubsurfaceData");
  ssd.def(py::init_implicit<GeodeticField3>())
      .def(py::init_implicit<Numeric>())
      .def(py::init_implicit<Subsurface::FunctionalData>())
      .def(
          "__init__",
          [](Subsurface::Data *a, const GriddedField3 &v) { new (a) Subsurface::Data(GeodeticField3(v)); },
          "v"_a,
          "Initialize with a gridded field")
      .def_rw(
          "data",
          &Subsurface::Data::data,
          "The data.\n\n.. :class:`~pyarts3.arts.GeodeticField3`\n\n.. :class:`~pyarts3.arts.Numeric`\n\n.. :class:`~pyarts3.arts.NumericTernaryOperator`")
      .def_rw("alt_upp",
              &Subsurface::Data::alt_upp,
              "Upper altitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("alt_low",
              &Subsurface::Data::alt_low,
              "Lower altitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("lat_upp",
              &Subsurface::Data::lat_upp,
              "Upper latitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("lat_low",
              &Subsurface::Data::lat_low,
              "Lower latitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("lon_upp",
              &Subsurface::Data::lon_upp,
              "Upper longitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("lon_low",
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
          [](const Subsurface::Data &d, Numeric alt, Numeric lat, Numeric lon) { return d.at(alt, lat, lon); },
          "alt"_a,
          "lat"_a,
          "lon"_a,
          "Get a point of data at the position")
      .def(
          "ws",
          [](const Subsurface::Data &d, Numeric alt, Numeric lat, Numeric lon) { return d.flat_weight(alt, lat, lon); },
          "alt"_a,
          "lat"_a,
          "lon"_a,
          "Get the weights of neighbors at a position")
      .def_prop_ro("data_type", &Subsurface::Data::data_type, "The data type\n\n.. :class:`~pyarts3.arts.String`");
  py::implicitly_convertible<Subsurface::FunctionalData::func_t, Subsurface::Data>();
  py::implicitly_convertible<GriddedField3, Subsurface::Data>();
  generic_interface(ssd);

  auto assp = py::bind_vector<Array<SubsurfacePoint>, py::rv_policy::reference_internal>(m, "ArrayOfSubsurfacePoint");
  generic_interface(assp);
  vector_interface(assp);
} catch (std::exception &e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize surf\n{}", e.what()));
}
}  // namespace Python
