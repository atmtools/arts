#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>
#include <surface_scattering/lambertian.h>
#include <surface_scattering/surface_scattering_model.h>
#include <surface_scattering/surface_scattering_properties.h>

#include "hpy_arts.h"

namespace Python {

namespace {
template <typename T> void model_interface(py::class_<T>& c, const char* doc) {
  c.def(
      "get_surface_scattering_model_properties",
      [](const T&            self,
         const SurfacePoint& surf_point,
         Numeric             lat,
         Numeric             lon,
         const Vector&       f_grid,
         const Vector&       za_inc_grid,
         const Vector&       aa_inc_grid,
         const Vector&       za_scat_grid,
         const Vector&       aa_scat_grid) {
        return self.get_surface_scattering_model_properties(
            surf_point, lat, lon, f_grid, za_inc_grid, aa_inc_grid, za_scat_grid, aa_scat_grid);
      },
      "surf_point"_a,
      "lat"_a,
      "lon"_a,
      "f_grid"_a,
      "za_inc_grid"_a,
      "aa_inc_grid"_a,
      "za_scat_grid"_a,
      "aa_scat_grid"_a,
      doc);
}
}  // namespace

void py_surface_scattering(py::module_& m) try {
  //
  // SurfaceScatteringModelProperties
  //
  py::class_<surface_scattering::SurfaceScatteringModelProperties>(
      m, "SurfaceScatteringModelProperties")
      .def(py::init<>())
      .def_rw("brdf_matrix_diffuse",
              &surface_scattering::SurfaceScatteringModelProperties::brdf_matrix_diffuse,
              "Diffuse BRDF Mueller matrix: dims [nf, nza_inc, naa_inc, nza_scat, naa_scat, 4, 4]\n\nNone if no model provides a diffuse component.\n\n.. :class:`~pyarts3.arts.MuelmatTensor5 | None`")
      .def_rw("emissivity_vector_diffuse",
              &surface_scattering::SurfaceScatteringModelProperties::emissivity_vector_diffuse,
              "Diffuse emissivity vector: dims [nf, nza_scat, naa_scat, 4, 4]\n\nNone if no model provides a diffuse component.\n\n.. :class:`~pyarts3.arts.MuelmatTensor3 | None`")
      .def_rw("brdf_matrix_specular",
              &surface_scattering::SurfaceScatteringModelProperties::brdf_matrix_specular,
              "Specular BRDF Mueller matrix: dims [nf, nza_inc, naa_inc, nza_scat, naa_scat, 4, 4]\n\nNone if no model provides a specular component.\n\n.. :class:`~pyarts3.arts.MuelmatTensor5 | None`")
      .def_rw("emissivity_vector_specular",
              &surface_scattering::SurfaceScatteringModelProperties::emissivity_vector_specular,
              "Specular emissivity vector: dims [nf, nza_scat, naa_scat, 4, 4]\n\nNone if no model provides a specular component.\n\n.. :class:`~pyarts3.arts.MuelmatTensor3 | None`")
      .doc() = "Bulk surface scattering properties (BRDF matrix + emissivity vector).\n\nComponents that no contributing model physically provides are None.";

  //
  // SurfaceSpectralData
  //
  py::class_<SurfaceSpectralData> ssd(m, "SurfaceSpectralData");
  generic_interface(ssd);
  ssd.def(py::init_implicit<SpectralGeodeticField3>())
      .def(py::init_implicit<Numeric>())
      .def(py::init_implicit<surface_scattering::FunctionalData>())
      .def_rw(
          "data",
          &SurfaceSpectralData::data,
          "The data.\n\n.. :class:`~pyarts3.arts.SpectralGeodeticField3`\n\n.. :class:`~pyarts3.arts.Numeric`\n\n.. :class:`~pyarts3.arts.NumericTernaryOperator`")
      .def_rw("frq_upp",
              &SurfaceSpectralData::frq_upp,
              "Upper frequency limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("frq_low",
              &SurfaceSpectralData::frq_low,
              "Lower frequency limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("lat_upp",
              &SurfaceSpectralData::lat_upp,
              "Upper latitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("lat_low",
              &SurfaceSpectralData::lat_low,
              "Lower latitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("lon_upp",
              &SurfaceSpectralData::lon_upp,
              "Upper longitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def_rw("lon_low",
              &SurfaceSpectralData::lon_low,
              "Lower longitude limit\n\n.. :class:`~pyarts3.arts.InterpolationExtrapolation`")
      .def(
          "set_extrapolation",
          [](SurfaceSpectralData &self, InterpolationExtrapolation x) {
            self.frq_upp = x;
            self.frq_low = x;
            self.lat_upp = x;
            self.lat_low = x;
            self.lon_upp = x;
            self.lon_low = x;
          },
          "extrapolation"_a,
          "Set the extrapolation for all dimensions")
      .def(
          "__call__",
          [](const SurfaceSpectralData &d, Numeric frq, Numeric lat, Numeric lon) { return d.at(frq, lat, lon); },
          "frq"_a,
          "lat"_a,
          "lon"_a,
          "Get a point of data at the frequency and position")
      .def(
          "ws",
          [](const SurfaceSpectralData &d, Numeric frq, Numeric lat, Numeric lon) {
            return d.flat_weights(frq, lat, lon);
          },
          "frq"_a,
          "lat"_a,
          "lon"_a,
          "Get the weights of neighbors at a frequency and position")
      .def_prop_ro("data_type", &SurfaceSpectralData::data_type, "The data type\n\n.. :class:`~pyarts3.arts.String`");
  py::implicitly_convertible<surface_scattering::FunctionalData::func_t, SurfaceSpectralData>();
  ssd.doc() = R"(Surface data over frequency [Hz], latitude [deg], and longitude [deg].

The data is a :class:`~pyarts3.arts.SpectralGeodeticField3` with frequency,
latitude, and longitude grids, a constant :class:`~pyarts3.arts.Numeric`,
or a function of (frequency, latitude, longitude).  Single-point grids and
constant or functional data extrapolate as the nearest value; otherwise the
extrapolation outside each grid is set per dimension.
)";

  //
  // LambertianSurfaceScatterer
  //
  py::class_<LambertianSurfaceScatterer> lss(m, "LambertianSurfaceScatterer");
  lss.def(py::init<SurfaceSpectralData>(),
           "reflectivity"_a,
           R"(Create a Lambertian surface scatterer.

Parameters
----------
reflectivity : SurfaceSpectralData
    Reflectivity over frequency [Hz], latitude [deg], and longitude [deg].
    Values are expected in [0, 1]; values outside this range are clamped
    at runtime.
)")
      .def_rw("reflectivity",
              &LambertianSurfaceScatterer::reflectivity,
              "Reflectivity over frequency, latitude, and longitude\n\n.. :class:`~pyarts3.arts.SurfaceSpectralData`");
  model_interface(lss, "Compute bulk surface scattering properties for this Lambertian model");
  generic_interface(lss);
  lss.doc() = R"(Lambertian (isotropic) surface scattering model.

The reflectivity is stored as :class:`~pyarts3.arts.SurfaceSpectralData`.
At runtime it is interpolated onto the simulation's ``f_grid`` at the surface
position, so the stored spectral and spatial resolutions are fully
independent of the simulation.
)";

  //
  // FresnelSurfaceScatterer
  //
  py::class_<FresnelSurfaceScatterer> fss(m, "FresnelSurfaceScatterer");
  fss.def(py::init<SurfaceSpectralData, Numeric>(),
           "refractive_index"_a,
           "n1"_a = 1.0,
           R"(Create a Fresnel surface scatterer.

Parameters
----------
refractive_index : SurfaceSpectralData
    Refractive index of the surface (medium 2) over frequency [Hz],
    latitude [deg], and longitude [deg].  Values must be positive.
n1 : float
    Refractive index of the medium the radiation propagates in (must be > 0).
)")
      .def_rw("refractive_index",
              &FresnelSurfaceScatterer::refractive_index,
              "Surface (medium 2) refractive index over frequency, latitude, and longitude\n\n.. :class:`~pyarts3.arts.SurfaceSpectralData`")
      .def_rw("n1",
              &FresnelSurfaceScatterer::n1,
              "Refractive index of the medium the radiation propagates in\n\n.. :class:`~pyarts3.arts.Numeric`");
  model_interface(fss, "Compute bulk surface scattering properties for this Fresnel model");
  generic_interface(fss);
  fss.doc() = R"(Fresnel (specular) surface scattering model.

Stores the surface refractive index (medium 2) as
:class:`~pyarts3.arts.SurfaceSpectralData`, plus the refractive index ``n1``
of the incident medium (default 1.0).  At runtime the refractive index is
interpolated onto the simulation ``f_grid`` at the surface position.  The
model fills the specular BRDF and the complementary emissivity; the diffuse
tensors are None.
)";

  //
  // FlatScalarSurfaceScatterer
  //
  py::class_<FlatScalarSurfaceScatterer> fsc(m, "FlatScalarSurfaceScatterer");
  fsc.def(py::init<SurfaceSpectralData>(),
           "reflectivity"_a,
           R"(Create a flat-scalar surface scatterer.

Parameters
----------
reflectivity : SurfaceSpectralData
    Scalar reflectivity over frequency [Hz], latitude [deg], and longitude
    [deg].  Values are expected in [0, 1]; values outside this range are
    clamped at runtime.
)")
      .def_rw("reflectivity",
              &FlatScalarSurfaceScatterer::reflectivity,
              "Scalar reflectivity over frequency, latitude, and longitude\n\n.. :class:`~pyarts3.arts.SurfaceSpectralData`");
  model_interface(fsc, "Compute bulk surface scattering properties for this flat-scalar model");
  generic_interface(fsc);
  fsc.doc() = R"(Flat-scalar (specular) surface scattering model.

The scalar reflectivity is stored as :class:`~pyarts3.arts.SurfaceSpectralData`.
At runtime it is interpolated onto the simulation's ``f_grid`` at the surface
position.  The model fills the specular BRDF (``Muelmat{R}``) and emissivity
(``Muelmat{1 - R}``); the diffuse tensors are None.
)";

  //
  // MapOfSurfaceScatteringModel
  //
  py::class_<MapOfSurfaceScatteringModel> mossm(m, "MapOfSurfaceScatteringModel");
  mossm.def(py::init<>())
      .def(
          "add",
          [](MapOfSurfaceScatteringModel& self,
             const std::string& name,
             const surface_scattering::SurfaceScatteringModel& model) {
            self.add(name, model);
          },
          "name"_a,
          "model"_a,
          "Insert or replace a named surface scattering model")
      .def(
          "__getitem__",
          [](const MapOfSurfaceScatteringModel& self, const std::string& name)
              -> const surface_scattering::SurfaceScatteringModel& {
            auto it = self.models.find(name);
            if (it == self.models.end())
              throw py::key_error(name.c_str());
            return it->second;
          },
          py::rv_policy::reference_internal,
          "name"_a)
      .def(
          "__setitem__",
          [](MapOfSurfaceScatteringModel& self,
             const std::string& name,
             const surface_scattering::SurfaceScatteringModel& model) {
            self.models[name] = model;
          },
          "name"_a,
          "model"_a)
      .def(
          "__contains__",
          [](const MapOfSurfaceScatteringModel& self, const std::string& name) {
            return self.models.contains(name);
          },
          "name"_a)
      .def(
          "__len__",
          [](const MapOfSurfaceScatteringModel& self) {
            return self.models.size();
          })
      .def(
          "__iter__",
          [](const MapOfSurfaceScatteringModel& self) {
            return py::make_iterator(py::type<MapOfSurfaceScatteringModel>(),
                                     "mossm-iterator",
                                     self.models.begin(),
                                     self.models.end());
          },
          py::rv_policy::reference_internal)
      .def(
          "get_surface_scattering_model_properties",
          [](const MapOfSurfaceScatteringModel& self,
             const SurfacePoint& surf_point,
             Numeric lat,
             Numeric lon,
             const Vector& f_grid,
             const Vector& za_inc_grid,
             const Vector& aa_inc_grid,
             const Vector& za_scat_grid,
             const Vector& aa_scat_grid) {
            return self.get_surface_scattering_model_properties(
                surf_point, lat, lon, f_grid, za_inc_grid, aa_inc_grid, za_scat_grid, aa_scat_grid);
          },
          "surf_point"_a,
          "lat"_a,
          "lon"_a,
          "f_grid"_a,
          "za_inc_grid"_a,
          "aa_inc_grid"_a,
          "za_scat_grid"_a,
          "aa_scat_grid"_a,
          "Accumulate bulk surface scattering properties from all stored models")
      .def_rw(
          "weighting_option",
          &MapOfSurfaceScatteringModel::weighting_option,
          "Weighting option for combining multiple models: Maximum or Average\n\n.. :class:`~pyarts3.arts.MapOfSurfaceScatteringModel.Weighting`");
  generic_interface(mossm);
  mossm.doc() = "Named map of surface scattering models";

  // Weighting enum for combining multiple models
  py::enum_<MapOfSurfaceScatteringModel::Weighting>(mossm, "Weighting")
      .value("Maximum", MapOfSurfaceScatteringModel::Weighting::Maximum,
             "Winner takes all: the model(s) with the largest surface mask value at the\n"
             "point get weight 1 (split equally on ties), all others 0.  All-zero masks\n"
             "give every model weight 1/N.")
      .value("Average", MapOfSurfaceScatteringModel::Weighting::Average,
             "Mask-weighted average: model weights are the surface mask values normalized\n"
             "to sum 1.  All-zero masks give all-zero weights (no scattering).");

} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize surface scattering:\n{}", e.what()));
}

}  // namespace Python
