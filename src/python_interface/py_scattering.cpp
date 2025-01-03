#include <nanobind/stl/bind_vector.h>
#include <optproperties.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "hpy_vector.h"
#include "py_macros.h"

namespace Python {
void py_scattering(py::module_& m) try {
  py::enum_<PType>(m, "PType")
      .value("PTYPE_GENERAL", PType::PTYPE_GENERAL, "As general")
      .value(
          "PTYPE_AZIMUTH_RND", PType::PTYPE_AZIMUTH_RND, "Azimuthally random")
      .value("PTYPE_TOTAL_RND", PType::PTYPE_TOTAL_RND, "Totally random")
      .PythonInterfaceCopyValue(PType)
      .def("__getstate__",
           [](const PType& self) {
             return py::make_tuple(static_cast<Index>(self));
           })
      .def("__setstate__", [](PType* self, const py::tuple& t) {
        ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")

        new (self) PType{py::cast<Index>(t[0])};
      });

  py::class_<SingleScatteringData> ssdb(m, "SingleScatteringData");
  workspace_group_interface(ssdb);
  ssdb.def_rw("ptype",
              &SingleScatteringData::ptype,
              ":class:`~pyarts.arts.PType` The type")
      .def_rw("description",
              &SingleScatteringData::description,
              ":class:`~pyarts.arts.String` The description")
      .def_rw("f_grid",
              &SingleScatteringData::f_grid,
              ":class:`~pyarts.arts.Vector` The frequency grid")
      .def_rw("T_grid",
              &SingleScatteringData::T_grid,
              ":class:`~pyarts.arts.Vector` The temperature grid")
      .def_rw("za_grid",
              &SingleScatteringData::za_grid,
              ":class:`~pyarts.arts.Vector` The zenith grid")
      .def_rw("aa_grid",
              &SingleScatteringData::aa_grid,
              ":class:`~pyarts.arts.Vector` The azimuth grid")
      .def_rw("pha_mat_data",
              &SingleScatteringData::pha_mat_data,
              ":class:`~pyarts.arts.Tensor7` The phase matrix")
      .def_rw("ext_mat_data",
              &SingleScatteringData::ext_mat_data,
              ":class:`~pyarts.arts.Tensor5` The extinction matrix")
      .def_rw("abs_vec_data",
              &SingleScatteringData::abs_vec_data,
              ":class:`~pyarts.arts.Tensor5` The absorption vector")
      .def("__getstate__",
           [](const SingleScatteringData& self) {
             return py::make_tuple(self.ptype,
                                   self.description,
                                   self.f_grid,
                                   self.T_grid,
                                   self.za_grid,
                                   self.aa_grid,
                                   self.pha_mat_data,
                                   self.ext_mat_data,
                                   self.abs_vec_data);
           })
      .def("__setstate__", [](SingleScatteringData* self, const py::tuple& t) {
        ARTS_USER_ERROR_IF(t.size() != 9, "Invalid state!")

        new (self)
            SingleScatteringData{.ptype        = py::cast<PType>(t[0]),
                                 .description  = py::cast<String>(t[1]),
                                 .f_grid       = py::cast<Vector>(t[2]),
                                 .T_grid       = py::cast<Vector>(t[3]),
                                 .za_grid      = py::cast<Vector>(t[4]),
                                 .aa_grid      = py::cast<Vector>(t[5]),
                                 .pha_mat_data = py::cast<Tensor7>(t[6]),
                                 .ext_mat_data = py::cast<Tensor5>(t[7]),
                                 .abs_vec_data = py::cast<Tensor5>(t[8])};
      });

  py::class_<ScatteringMetaData> smdb(m, "ScatteringMetaData");
  workspace_group_interface(smdb);
  smdb.def_rw("description",
              &ScatteringMetaData::description,
              ":class:`~pyarts.arts.String` The description")
      .def_rw("source",
              &ScatteringMetaData::source,
              ":class:`~pyarts.arts.String` The source")
      .def_rw("refr_index",
              &ScatteringMetaData::refr_index,
              ":class:`~pyarts.arts.String` The refractive index")
      .def_rw("mass", &ScatteringMetaData::mass, ":class:`float` The mass")
      .def_rw("diameter_max",
              &ScatteringMetaData::diameter_max,
              ":class:`float` The max diameter")
      .def_rw("diameter_volume_equ",
              &ScatteringMetaData::diameter_volume_equ,
              ":class:`float` The volume equivalent diameter")
      .def_rw("diameter_area_equ_aerodynamical",
              &ScatteringMetaData::diameter_area_equ_aerodynamical,
              ":class:`float` The diameter area equivalent")
      .def("__getstate__",
           [](const ScatteringMetaData& self) {
             return py::make_tuple(self.description,
                                   self.source,
                                   self.refr_index,
                                   self.mass,
                                   self.diameter_max,
                                   self.diameter_volume_equ,
                                   self.diameter_area_equ_aerodynamical);
           })
      .def("__setstate__", [](ScatteringMetaData* self, const py::tuple& t) {
        ARTS_USER_ERROR_IF(t.size() != 7, "Invalid state!")

        new (self) ScatteringMetaData{
            .description                     = py::cast<String>(t[0]),
            .source                          = py::cast<String>(t[1]),
            .refr_index                      = py::cast<String>(t[2]),
            .mass                            = py::cast<Numeric>(t[3]),
            .diameter_max                    = py::cast<Numeric>(t[4]),
            .diameter_volume_equ             = py::cast<Numeric>(t[5]),
            .diameter_area_equ_aerodynamical = py::cast<Numeric>(t[6])};
      });

  auto a1 = py::bind_vector<ArrayOfScatteringMetaData,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfScatteringMetaData");
  workspace_group_interface(a1);
  vector_interface(a1);
  auto a2 = py::bind_vector<ArrayOfArrayOfScatteringMetaData,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfScatteringMetaData");
  workspace_group_interface(a2);
  vector_interface(a2);
  auto a3 = py::bind_vector<ArrayOfSingleScatteringData,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfSingleScatteringData");
  workspace_group_interface(a3);
  vector_interface(a3);
  auto a4 = py::bind_vector<ArrayOfArrayOfSingleScatteringData,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfSingleScatteringData");
  workspace_group_interface(a4);
  vector_interface(a4);

} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize scattering\n{}", e.what()));
}
}  // namespace Python
