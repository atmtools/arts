#include <python_interface.h>

#include "py_macros.h"

namespace Python {
void py_scattering(py::module_& m) try {
  py::enum_<PType>(m, "PType")
      .value("PTYPE_GENERAL", PType::PTYPE_GENERAL, "As general")
      .value("PTYPE_AZIMUTH_RND", PType::PTYPE_AZIMUTH_RND, "Azimuthally random")
      .value("PTYPE_TOTAL_RND", PType::PTYPE_TOTAL_RND, "Totallt random")
      .PythonInterfaceCopyValue(PType)
      .def(py::pickle(
          [](const PType& self) {
            return py::make_tuple(static_cast<Index>(self));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            
            return static_cast<PType>(t[0].cast<Index>());
          }));

  artsclass<SingleScatteringData>(m, "SingleScatteringData")
      .def(py::init([]() { return std::make_shared<SingleScatteringData>(); }), "Empty scattering data")
      .PythonInterfaceCopyValue(SingleScatteringData)
      .PythonInterfaceWorkspaceVariableConversion(SingleScatteringData)
      .PythonInterfaceFileIO(SingleScatteringData)
      .PythonInterfaceBasicRepresentation(SingleScatteringData)
      .def_readwrite("ptype", &SingleScatteringData::ptype, ":class:`~pyarts.arts.PType` The type")
      .def_readwrite("description", &SingleScatteringData::description, ":class:`~pyarts.arts.String` The description")
      .def_readwrite("f_grid", &SingleScatteringData::f_grid, ":class:`~pyarts.arts.Vector` The frequency grid")
      .def_readwrite("T_grid", &SingleScatteringData::T_grid, ":class:`~pyarts.arts.Vector` The temperature grid")
      .def_readwrite("za_grid", &SingleScatteringData::za_grid, ":class:`~pyarts.arts.Vector` The zenith grid")
      .def_readwrite("aa_grid", &SingleScatteringData::aa_grid, ":class:`~pyarts.arts.Vector` The azimuth grid")
      .def_readwrite("pha_mat_data", &SingleScatteringData::pha_mat_data, ":class:`~pyarts.arts.Tensor7` The phase matrix")
      .def_readwrite("ext_mat_data", &SingleScatteringData::ext_mat_data, ":class:`~pyarts.arts.Tensor5` The extinction matrix")
      .def_readwrite("abs_vec_data", &SingleScatteringData::abs_vec_data, ":class:`~pyarts.arts.Tensor5` The absorption vector")
      .def(py::pickle(
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
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 9, "Invalid state!")

            return std::make_shared<SingleScatteringData>(SingleScatteringData{t[0].cast<PType>(),
                                            t[1].cast<String>(),
                                            t[2].cast<Vector>(),
                                            t[3].cast<Vector>(),
                                            t[4].cast<Vector>(),
                                            t[5].cast<Vector>(),
                                            t[6].cast<Tensor7>(),
                                            t[7].cast<Tensor5>(),
                                            t[8].cast<Tensor5>()});
          }))
      .PythonInterfaceWorkspaceDocumentation(SingleScatteringData);

  artsclass<ScatteringMetaData>(m, "ScatteringMetaData")
      .def(py::init([]() { return std::make_shared<ScatteringMetaData>(); }), "Empty meta data")
      .PythonInterfaceCopyValue(ScatteringMetaData)
      .PythonInterfaceWorkspaceVariableConversion(ScatteringMetaData)
      .PythonInterfaceFileIO(ScatteringMetaData)
      .PythonInterfaceBasicRepresentation(ScatteringMetaData)
      .def_readwrite("description", &ScatteringMetaData::description, ":class:`~pyarts.arts.String` The description")
      .def_readwrite("source", &ScatteringMetaData::source, ":class:`~pyarts.arts.String` The source")
      .def_readwrite("refr_index", &ScatteringMetaData::refr_index, ":class:`~pyarts.arts.String` The refractive index")
      .def_readwrite("mass", &ScatteringMetaData::mass, ":class:`float` The mass")
      .def_readwrite("diameter_max", &ScatteringMetaData::diameter_max, ":class:`float` The max diameter")
      .def_readwrite("diameter_volume_equ",
                     &ScatteringMetaData::diameter_volume_equ, ":class:`float` The volume equivalent diameter")
      .def_readwrite("diameter_area_equ_aerodynamical",
                     &ScatteringMetaData::diameter_area_equ_aerodynamical, ":class:`float` The diameter area equivalent")
      .def(py::pickle(
          [](const ScatteringMetaData& self) {
            return py::make_tuple(self.description,
                                  self.source,
                                  self.refr_index,
                                  self.mass,
                                  self.diameter_max,
                                  self.diameter_volume_equ,
                                  self.diameter_area_equ_aerodynamical);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 7, "Invalid state!")

            return std::make_shared<ScatteringMetaData>(ScatteringMetaData{t[0].cast<String>(),
                                          t[1].cast<String>(),
                                          t[2].cast<String>(),
                                          t[3].cast<Numeric>(),
                                          t[4].cast<Numeric>(),
                                          t[5].cast<Numeric>(),
                                          t[6].cast<Numeric>()});
          }))
      .PythonInterfaceWorkspaceDocumentation(ScatteringMetaData);

  PythonInterfaceWorkspaceArray(ScatteringMetaData);
  PythonInterfaceWorkspaceArray(SingleScatteringData);
  PythonInterfaceWorkspaceArray(ArrayOfScatteringMetaData);
  PythonInterfaceWorkspaceArray(ArrayOfSingleScatteringData);
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize scattering\n", e.what()));
}
}  // namespace Python
