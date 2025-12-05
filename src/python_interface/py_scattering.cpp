#include <nanobind/stl/bind_vector.h>
#include <optproperties.h>
#include <python_interface.h>

#include "hpy_arts.h"
#include "hpy_vector.h"

namespace Python {
void py_scattering(py::module_& m) try {
  py::enum_<PType>(m, "PType")
      .value("PTYPE_GENERAL", PType::PTYPE_GENERAL, "As general")
      .value(
          "PTYPE_AZIMUTH_RND", PType::PTYPE_AZIMUTH_RND, "Azimuthally random")
      .value("PTYPE_TOTAL_RND", PType::PTYPE_TOTAL_RND, "Totally random");

  py::class_<SingleScatteringData> ssdb(m, "SingleScatteringData");
  generic_interface(ssdb);
  ssdb.def_rw("ptype",
              &SingleScatteringData::ptype,
              "The type\n\n.. :class:`~pyarts3.arts.PType`")
      .def_rw("description",
              &SingleScatteringData::description,
              "The description\n\n.. :class:`~pyarts3.arts.String`")
      .def_rw("f_grid",
              &SingleScatteringData::f_grid,
              "The frequency grid\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("T_grid",
              &SingleScatteringData::T_grid,
              "The temperature grid\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("za_grid",
              &SingleScatteringData::za_grid,
              "The zenith grid\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("aa_grid",
              &SingleScatteringData::aa_grid,
              "The azimuth grid\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_rw("pha_mat_data",
              &SingleScatteringData::pha_mat_data,
              "The phase matrix\n\n.. :class:`~pyarts3.arts.Tensor7`")
      .def_rw("ext_mat_data",
              &SingleScatteringData::ext_mat_data,
              "The extinction matrix\n\n.. :class:`~pyarts3.arts.Tensor5`")
      .def_rw("abs_vec_data",
              &SingleScatteringData::abs_vec_data,
              "The absorption vector\n\n.. :class:`~pyarts3.arts.Tensor5`");

  py::class_<ScatteringMetaData> smdb(m, "ScatteringMetaData");
  generic_interface(smdb);
  smdb.def_rw("description",
              &ScatteringMetaData::description,
              "The description\n\n.. :class:`~pyarts3.arts.String`")
      .def_rw("source",
              &ScatteringMetaData::source,
              "The source\n\n.. :class:`~pyarts3.arts.String`")
      .def_rw("refr_index",
              &ScatteringMetaData::refr_index,
              "The refractive index\n\n.. :class:`~pyarts3.arts.String`")
      .def_rw(
          "mass", &ScatteringMetaData::mass, "The mass\n\n.. :class:`float`")
      .def_rw("diameter_max",
              &ScatteringMetaData::diameter_max,
              "The max diameter\n\n.. :class:`float`")
      .def_rw("diameter_volume_equ",
              &ScatteringMetaData::diameter_volume_equ,
              "The volume equivalent diameter\n\n.. :class:`float`")
      .def_rw("diameter_area_equ_aerodynamical",
              &ScatteringMetaData::diameter_area_equ_aerodynamical,
              "The diameter area equivalent\n\n.. :class:`float`");

  auto a1 = py::bind_vector<ArrayOfScatteringMetaData,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfScatteringMetaData");
  generic_interface(a1);
  vector_interface(a1);
  auto a2 = py::bind_vector<ArrayOfArrayOfScatteringMetaData,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfScatteringMetaData");
  generic_interface(a2);
  vector_interface(a2);
  auto a3 = py::bind_vector<ArrayOfSingleScatteringData,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfSingleScatteringData");
  generic_interface(a3);
  vector_interface(a3);
  auto a4 = py::bind_vector<ArrayOfArrayOfSingleScatteringData,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfSingleScatteringData");
  generic_interface(a4);
  vector_interface(a4);

} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize scattering\n{}", e.what()));
}
}  // namespace Python
