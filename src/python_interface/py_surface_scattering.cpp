#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <python_interface.h>
#include <surface_scattering/lambertian.h>
#include <surface_scattering/surface_scattering_model.h>

#include "hpy_arts.h"

namespace Python {

void py_surface_scattering(py::module_& m) try {
  //
  // LambertianSurfaceScatterer
  //
  py::class_<LambertianSurfaceScatterer> lss(m, "LambertianSurfaceScatterer");
  lss.def(py::init<>())
      .def(py::init<SurfacePropertyTag, Vector>(),
           "reflectivity_tag"_a,
           "reflectivity"_a,
           "Create a Lambertian surface scatterer from a tag and reflectivity vector")
      .def_rw("reflectivity_tag",
              &LambertianSurfaceScatterer::reflectivity_tag,
              "Surface property tag identifying this model\n\n.. :class:`SurfacePropertyTag`")
      .def_prop_rw(
          "reflectivity",
          [](const LambertianSurfaceScatterer& self) { return self.get_reflectivity(); },
          [](LambertianSurfaceScatterer& self, const Vector& r) { self.set_reflectivity(r); },
          "Reflectivity vector over f_grid\n\n.. :class:`Vector`")
      .def(
          "get_bulk_surface_scattering_properties",
          [](const LambertianSurfaceScatterer& self,
             const SurfacePoint& surf_point,
             const Vector& f_grid,
             const Vector& za_inc_grid,
             const Vector& aa_inc_grid,
             const Vector& za_scat_grid,
             const Vector& aa_scat_grid) {
            return self.get_surface_scattering_model_properties(
                surf_point, f_grid, za_inc_grid, aa_inc_grid, za_scat_grid, aa_scat_grid);
          },
          "surf_point"_a,
          "f_grid"_a,
          "za_inc_grid"_a,
          "aa_inc_grid"_a,
          "za_scat_grid"_a,
          "aa_scat_grid"_a,
          "Compute bulk surface scattering properties for this Lambertian model");
  generic_interface(lss);
  lss.doc() = "Lambertian surface scattering model";

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
          "get_bulk_surface_scattering_properties",
          [](const MapOfSurfaceScatteringModel& self,
             const SurfacePoint& surf_point,
             const Vector& f_grid,
             const Vector& za_inc_grid,
             const Vector& aa_inc_grid,
             const Vector& za_scat_grid,
             const Vector& aa_scat_grid) {
            return self.get_surface_scattering_model_properties(
                surf_point, f_grid, za_inc_grid, aa_inc_grid, za_scat_grid, aa_scat_grid);
          },
          "surf_point"_a,
          "f_grid"_a,
          "za_inc_grid"_a,
          "aa_inc_grid"_a,
          "za_scat_grid"_a,
          "aa_scat_grid"_a,
          "Accumulate bulk surface scattering properties from all stored models");
  generic_interface(mossm);
  mossm.doc() = "Named map of surface scattering models";

} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize surface scattering:\n{}", e.what()));
}

}  // namespace Python

