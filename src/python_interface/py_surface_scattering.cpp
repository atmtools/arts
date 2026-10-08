#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>
#include <surface_scattering/lambertian.h>
#include <surface_scattering/surface_scattering_model.h>
#include <surface_scattering/surface_scattering_properties.h>

#include "hpy_arts.h"

namespace Python {

void py_surface_scattering(py::module_& m) try {
  //
  // SurfaceScatteringModelProperties
  //
  py::class_<surface_scattering::SurfaceScatteringModelProperties>(
      m, "SurfaceScatteringModelProperties")
      .def(py::init<>())
      .def_rw("brdf_matrix_diffuse",
              &surface_scattering::SurfaceScatteringModelProperties::brdf_matrix_diffuse,
              "Diffuse BRDF Mueller matrix: dims [nf, nza_inc, naa_inc, nza_scat, naa_scat, 4, 4]\n\n.. :class:`class-information`")
      .def_rw("emissivity_vector_diffuse",
              &surface_scattering::SurfaceScatteringModelProperties::emissivity_vector_diffuse,
              "Diffuse emissivity vector: dims [nf, nza_scat, naa_scat, 4, 4]\n\n.. :class:`class-information`")
      .def_rw("brdf_matrix_specular",
              &surface_scattering::SurfaceScatteringModelProperties::brdf_matrix_specular,
              "Specular BRDF Mueller matrix: dims [nf, nza_inc, naa_inc, nza_scat, naa_scat, 4, 4]\n\n.. :class:`class-information`")
      .def_rw("emissivity_vector_specular",
              &surface_scattering::SurfaceScatteringModelProperties::emissivity_vector_specular,
              "Specular emissivity vector: dims [nf, nza_scat, naa_scat, 4, 4]\n\n.. :class:`class-information`")
      .doc() = "Bulk surface scattering properties (BRDF matrix + emissivity vector).";

  //
  // LambertianSurfaceScatterer
  //
  py::class_<LambertianSurfaceScatterer> lss(m, "LambertianSurfaceScatterer");
  lss.def(py::init<>())
      .def(py::init<SortedGriddedField1>(),
           "reflectivity_spectrum"_a,
           R"(Create a Lambertian surface scatterer from a spectral reflectivity field.

Parameters
----------
reflectivity_spectrum : SortedGriddedField1
    Reflectivity as a function of frequency [Hz].  The frequency grid must
    be sorted in ascending order.  Values are expected in [0, 1]; values
    outside this range are clamped at runtime.
)")
      .def_prop_rw(
          "reflectivity_spectrum",
          [](const LambertianSurfaceScatterer& self) {
            return self.get_reflectivity_spectrum();
          },
          [](LambertianSurfaceScatterer& self, const SortedGriddedField1& s) {
            self.set_reflectivity_spectrum(s);
          },
          R"(Spectral reflectivity field on an arbitrary sorted frequency grid.

The frequency axis must be in Hz (ascending).  Values should lie in [0, 1];
they are clamped when the BRDF is computed.

.. :class:`SortedGriddedField1`
)")
      .def(
          "get_surface_scattering_model_properties",
          [](const LambertianSurfaceScatterer& self,
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
          "Compute bulk surface scattering properties for this Lambertian model");
  generic_interface(lss);
  lss.doc() = R"(Lambertian (isotropic) surface scattering model.

The reflectivity is stored as a :class:`~pyarts3.arts.SortedGriddedField1`
on an arbitrary sorted frequency grid.  At runtime the spectrum is linearly
interpolated onto the simulation's ``f_grid``, so the stored spectral
resolution is fully independent of the simulation grid.
)";

  //
  // LambertianSurfaceScattererField
  //
  py::class_<LambertianSurfaceScattererField> lssf(m, "LambertianSurfaceScattererField");
  lssf.def(py::init<>())
      .def(py::init<SortedGriddedField3>(),
           "reflectivity_field"_a,
           R"(Create a spatially-varying Lambertian surface scatterer.

Parameters
----------
reflectivity_field : SortedGriddedField3
    Reflectivity as a function of latitude [deg], longitude [deg], and
    frequency [Hz].  All three grids must be sorted in ascending order.
    Values are expected in [0, 1]; values outside this range are clamped
    at runtime.
)")
      
      .def_prop_rw(
          "reflectivity_field",
          [](const LambertianSurfaceScattererField& self) {
            return self.get_reflectivity_field();
          },
          [](LambertianSurfaceScattererField& self, const SortedGriddedField3& f) {
            self.set_reflectivity_field(f);
          },
          R"(Spatially-varying reflectivity field on sorted (lat, lon, freq) grids.

Grid dimensions: latitude [deg], longitude [deg], frequency [Hz] — all ascending.
Values should lie in [0, 1]; they are clamped when the BRDF is computed.

.. :class:`SortedGriddedField3`
)")
      .def_rw(
          "interp_extrapolation",
          &LambertianSurfaceScattererField::interp_extrapolation,
          "Interpolation and extrapolation method for the frequency dimension\n\n.. :class:`class-information`")
      .def(
          "get_surface_scattering_model_properties",
          [](const LambertianSurfaceScattererField& self,
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
          "Compute bulk surface scattering properties for this spatially-varying Lambertian model");
  generic_interface(lssf);
  lssf.doc() = R"(Spatially-varying Lambertian (isotropic) surface scattering model.

The reflectivity is stored as a :class:`~pyarts3.arts.SortedGriddedField3`
on sorted (latitude [deg], longitude [deg], frequency [Hz]) grids.
At runtime the field is bilinearly interpolated in the geographic dimensions
and linearly interpolated onto the simulation's ``f_grid``, so the stored
spatial and spectral resolutions are fully independent of the simulation.
)";

  //
  // FresnelSurfaceScatterer
  //
  {
    py::class_<FresnelSurfaceScatterer> fss(m, "FresnelSurfaceScatterer");
    fss.def(py::init<>())
        .def(py::init<SortedGriddedField1>(),
             "refractive_index_spectrum"_a,
             R"(Create a Fresnel surface scatterer from a spectral refractive index field.

Parameters
----------
refractive_index_spectrum : SortedGriddedField1
    Refractive index of the surface (medium 2) as a function of frequency
    [Hz].  The frequency grid must be sorted in ascending order.  Values
    must be positive.
)")
        .def(py::init<SortedGriddedField1, Numeric>(),
             "refractive_index_spectrum"_a,
             "n1"_a,
             R"(Create a Fresnel surface scatterer from a spectral refractive index field.

Parameters
----------
refractive_index_spectrum : SortedGriddedField1
    Refractive index of the surface (medium 2) as a function of frequency [Hz].
n1 : float
    Refractive index of the medium the radiation propagates in (must be > 0).
)")
        .def_prop_rw(
            "refractive_index_spectrum",
            [](const FresnelSurfaceScatterer& self) {
              return self.get_refractive_index_spectrum();
            },
            [](FresnelSurfaceScatterer& self, const SortedGriddedField1& s) {
              self.set_refractive_index_spectrum(s);
            },
            R"(Refractive index of the surface (medium 2) on a sorted frequency grid [Hz].

Must be positive.

.. :class:`SortedGriddedField1`
)")
        .def_rw(
            "n1",
            &FresnelSurfaceScatterer::n1,
            "Refractive index of medium 1 (the medium the radiation propagates in)\n\n.. :class:`class-information`")
        .def_rw(
            "interp_extrapolation",
            &FresnelSurfaceScatterer::interp_extrapolation,
            "Interpolation and extrapolation method for frequency\n\n.. :class:`class-information`")
        .def(
            "get_surface_scattering_model_properties",
            [](const FresnelSurfaceScatterer& self,
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
            "Compute bulk surface scattering properties for this Fresnel model");
    generic_interface(fss);
    fss.doc() = R"(Fresnel (specular) surface scattering model.

Stores the surface refractive index (medium 2) as a
:class:`~pyarts3.arts.SortedGriddedField1` on an arbitrary sorted frequency
grid, plus the refractive index ``n1`` of the incident medium (default 1.0).
At runtime the spectrum is linearly interpolated onto the simulation
``f_grid``.  The model fills the specular BRDF and the complementary
emissivity; the diffuse tensors remain zero.
)";
  }

  //
  // FresnelSurfaceScattererField
  //
  {
    py::class_<FresnelSurfaceScattererField> fssf(m, "FresnelSurfaceScattererField");
    fssf.def(py::init<>())
        .def(py::init<SortedGriddedField3>(),
             "refractive_index_field"_a,
             R"(Create a spatially-varying Fresnel surface scatterer.

Parameters
----------
refractive_index_field : SortedGriddedField3
    Refractive index of the surface as a function of latitude [deg],
    longitude [deg], and frequency [Hz].  All three grids must be sorted in
    ascending order.  Values must be positive.
)")
        .def(py::init<SortedGriddedField3, Numeric>(),
             "refractive_index_field"_a,
             "n1"_a,
             R"(Create a spatially-varying Fresnel surface scatterer.

Parameters
----------
refractive_index_field : SortedGriddedField3
    Refractive index of the surface as a function of latitude [deg],
    longitude [deg], and frequency [Hz].
n1 : float
    Refractive index of the medium the radiation propagates in (must be > 0).
)")
        .def_prop_rw(
            "refractive_index_field",
            [](const FresnelSurfaceScattererField& self) {
              return self.get_refractive_index_field();
            },
            [](FresnelSurfaceScattererField& self, const SortedGriddedField3& f) {
              self.set_refractive_index_field(f);
            },
            R"(Spatially-varying surface refractive index on sorted (lat, lon, freq) grids.

Grid dimensions: latitude [deg], longitude [deg], frequency [Hz] — all
ascending.  Values must be positive.

.. :class:`SortedGriddedField3`
)")
        .def_rw(
            "n1",
            &FresnelSurfaceScattererField::n1,
            "Refractive index of medium 1 (the medium the radiation propagates in)\n\n.. :class:`class-information`")
        .def_rw(
            "interp_extrapolation",
            &FresnelSurfaceScattererField::interp_extrapolation,
            "Interpolation and extrapolation method for latitude, longitude, and frequency dimensions\n\n.. :class:`class-information`")
        .def(
            "get_surface_scattering_model_properties",
            [](const FresnelSurfaceScattererField& self,
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
            "Compute bulk surface scattering properties for this spatially-varying Fresnel model");
    generic_interface(fssf);
    fssf.doc() = R"(Spatially-varying Fresnel (specular) surface scattering model.

Stores the surface refractive index as a :class:`~pyarts3.arts.SortedGriddedField3`
on sorted (latitude [deg], longitude [deg], frequency [Hz]) grids, plus the
refractive index ``n1`` of the incident medium (default 1.0).  At runtime the
field is bilinearly interpolated in the geographic dimensions and linearly
interpolated onto the simulation ``f_grid``.  The model fills the specular
BRDF and the complementary emissivity; the diffuse tensors remain zero.
)";
  }

  //
  // FlatScalarSurfaceScatterer
  //
  {
    py::class_<FlatScalarSurfaceScatterer> fsc(m, "FlatScalarSurfaceScatterer");
    fsc.def(py::init<>())
        .def(py::init<SortedGriddedField1>(),
             "reflectivity_spectrum"_a,
             R"(Create a flat-scalar surface scatterer from a spectral reflectivity field.

Parameters
----------
reflectivity_spectrum : SortedGriddedField1
    Scalar reflectivity as a function of frequency [Hz].  The frequency grid
    must be sorted in ascending order.  Values are expected in [0, 1]; values
    outside this range are clamped at runtime.
)")
        .def_prop_rw(
            "reflectivity_spectrum",
            [](const FlatScalarSurfaceScatterer& self) {
              return self.get_reflectivity_spectrum();
            },
            [](FlatScalarSurfaceScatterer& self, const SortedGriddedField1& s) {
              self.set_reflectivity_spectrum(s);
            },
            R"(Spectral reflectivity field on an arbitrary sorted frequency grid.

The frequency axis must be in Hz (ascending).  Values should lie in [0, 1];
they are clamped when the BRDF is computed.

.. :class:`SortedGriddedField1`
)")
        .def_rw(
            "interp_extrapolation",
            &FlatScalarSurfaceScatterer::interp_extrapolation,
            "Interpolation and extrapolation method for frequency\n\n.. :class:`class-information`")
        .def(
            "get_surface_scattering_model_properties",
            [](const FlatScalarSurfaceScatterer& self,
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
            "Compute bulk surface scattering properties for this flat-scalar model");
    generic_interface(fsc);
    fsc.doc() = R"(Flat-scalar (specular) surface scattering model.

The scalar reflectivity is stored as a :class:`~pyarts3.arts.SortedGriddedField1`
on an arbitrary sorted frequency grid.  At runtime the spectrum is linearly
interpolated onto the simulation's ``f_grid``, so the stored spectral
resolution is fully independent of the simulation grid.  The model fills the
specular BRDF (``Muelmat{R}``) and emissivity (``Muelmat{1 - R}``); the
diffuse tensors remain zero.
)";
  }

  //
  // FlatScalarSurfaceScattererField
  //
  {
    py::class_<FlatScalarSurfaceScattererField> fscf(m, "FlatScalarSurfaceScattererField");
    fscf.def(py::init<>())
        .def(py::init<SortedGriddedField3>(),
             "reflectivity_field"_a,
             R"(Create a spatially-varying flat-scalar surface scatterer.

Parameters
----------
reflectivity_field : SortedGriddedField3
    Scalar reflectivity as a function of latitude [deg], longitude [deg], and
    frequency [Hz].  All three grids must be sorted in ascending order.
    Values are expected in [0, 1]; values outside this range are clamped
    at runtime.
)")
        .def_prop_rw(
            "reflectivity_field",
            [](const FlatScalarSurfaceScattererField& self) {
              return self.get_reflectivity_field();
            },
            [](FlatScalarSurfaceScattererField& self, const SortedGriddedField3& f) {
              self.set_reflectivity_field(f);
            },
            R"(Spatially-varying scalar reflectivity field on sorted (lat, lon, freq) grids.

Grid dimensions: latitude [deg], longitude [deg], frequency [Hz] — all ascending.
Values should lie in [0, 1]; they are clamped when the BRDF is computed.

.. :class:`SortedGriddedField3`
)")
        .def_rw(
            "interp_extrapolation",
            &FlatScalarSurfaceScattererField::interp_extrapolation,
            "Interpolation and extrapolation method for latitude, longitude, and frequency dimensions\n\n.. :class:`class-information`")
        .def(
            "get_surface_scattering_model_properties",
            [](const FlatScalarSurfaceScattererField& self,
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
            "Compute bulk surface scattering properties for this spatially-varying flat-scalar model");
    generic_interface(fscf);
    fscf.doc() = R"(Spatially-varying flat-scalar (specular) surface scattering model.

Stores the scalar reflectivity as a :class:`~pyarts3.arts.SortedGriddedField3`
on sorted (latitude [deg], longitude [deg], frequency [Hz]) grids.  At runtime
the field is bilinearly interpolated in the geographic dimensions and linearly
interpolated onto the simulation's ``f_grid``.  The model fills the specular
BRDF (``Muelmat{R}``) and emissivity (``Muelmat{1 - R}``); the diffuse tensors
remain zero.
)";
  }

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
          "Weighting option for combining multiple models: Maximum or Average\n\n.. :class:`class-information`");
  generic_interface(mossm);
  mossm.doc() = "Named map of surface scattering models";

  // Weighting enum for combining multiple models
  py::enum_<MapOfSurfaceScatteringModel::Weighting>(mossm, "Weighting")
      .value("Maximum", MapOfSurfaceScatteringModel::Weighting::Maximum,
             "Take the maximum value across models for each BRDF/emissivity element")
      .value("Average", MapOfSurfaceScatteringModel::Weighting::Average,
             "Take the average value across models for each BRDF/emissivity element");

} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize surface scattering:\n{}", e.what()));
}

}  // namespace Python
