#include <py_auto_interface.h>

#include "py_macros.h"

namespace Python {
void py_scattering(py::module_& m) {
  py::enum_<PType>(m, "PType")
      .value("PTYPE_GENERAL", PType::PTYPE_GENERAL)
      .value("PTYPE_AZIMUTH_RND", PType::PTYPE_AZIMUTH_RND)
      .value("PTYPE_TOTAL_RND", PType::PTYPE_TOTAL_RND)
      .PythonInterfaceCopyValue(PType)
      .def(py::pickle(
          [](const PType& self) {
            return py::make_tuple(static_cast<Index>(self));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            
            return static_cast<PType>(t[0].cast<Index>());
          }));

  py::class_<SingleScatteringData>(m, "SingleScatteringData")
      .def(py::init([]() { return new SingleScatteringData{}; }))
      .PythonInterfaceCopyValue(SingleScatteringData)
      .PythonInterfaceWorkspaceVariableConversion(SingleScatteringData)
      .PythonInterfaceFileIO(SingleScatteringData)
      .PythonInterfaceBasicRepresentation(SingleScatteringData)
      .def_readwrite("ptype", &SingleScatteringData::ptype)
      .def_readwrite("description", &SingleScatteringData::description)
      .def_readwrite("f_grid", &SingleScatteringData::f_grid)
      .def_readwrite("T_grid", &SingleScatteringData::T_grid)
      .def_readwrite("za_grid", &SingleScatteringData::za_grid)
      .def_readwrite("aa_grid", &SingleScatteringData::aa_grid)
      .def_readwrite("pha_mat_data", &SingleScatteringData::pha_mat_data)
      .def_readwrite("ext_mat_data", &SingleScatteringData::ext_mat_data)
      .def_readwrite("abs_vec_data", &SingleScatteringData::abs_vec_data)
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

            return new SingleScatteringData{t[0].cast<PType>(),
                                            t[1].cast<String>(),
                                            t[2].cast<Vector>(),
                                            t[3].cast<Vector>(),
                                            t[4].cast<Vector>(),
                                            t[5].cast<Vector>(),
                                            t[6].cast<Tensor7>(),
                                            t[7].cast<Tensor5>(),
                                            t[8].cast<Tensor5>()};
          }))
      .PythonInterfaceWorkspaceDocumentation(SingleScatteringData);

  py::class_<ScatteringMetaData>(m, "ScatteringMetaData")
      .def(py::init([]() { return new ScatteringMetaData{}; }))
      .PythonInterfaceCopyValue(ScatteringMetaData)
      .PythonInterfaceWorkspaceVariableConversion(ScatteringMetaData)
      .PythonInterfaceFileIO(ScatteringMetaData)
      .PythonInterfaceBasicRepresentation(ScatteringMetaData)
      .def_readwrite("description", &ScatteringMetaData::description)
      .def_readwrite("source", &ScatteringMetaData::source)
      .def_readwrite("refr_index", &ScatteringMetaData::refr_index)
      .def_readwrite("mass", &ScatteringMetaData::mass)
      .def_readwrite("diameter_max", &ScatteringMetaData::diameter_max)
      .def_readwrite("diameter_volume_equ",
                     &ScatteringMetaData::diameter_volume_equ)
      .def_readwrite("diameter_area_equ_aerodynamical",
                     &ScatteringMetaData::diameter_area_equ_aerodynamical)
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

            return new ScatteringMetaData{t[0].cast<String>(),
                                          t[1].cast<String>(),
                                          t[2].cast<String>(),
                                          t[3].cast<Numeric>(),
                                          t[4].cast<Numeric>(),
                                          t[5].cast<Numeric>(),
                                          t[6].cast<Numeric>()};
          }))
      .PythonInterfaceWorkspaceDocumentation(ScatteringMetaData);

  PythonInterfaceWorkspaceArray(ScatteringMetaData);
  PythonInterfaceWorkspaceArray(SingleScatteringData);
  PythonInterfaceWorkspaceArray(ArrayOfScatteringMetaData);
  PythonInterfaceWorkspaceArray(ArrayOfSingleScatteringData);
}
}  // namespace Python
