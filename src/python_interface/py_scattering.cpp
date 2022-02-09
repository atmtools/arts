#include "py_macros.h"
#include <py_auto_interface.h>

namespace Python {
void py_scattering(py::module_& m) {
  py::enum_<PType>(m, "PType")
      .value("PTYPE_GENERAL", PType::PTYPE_GENERAL)
      .value("PTYPE_AZIMUTH_RND", PType::PTYPE_AZIMUTH_RND)
      .value("PTYPE_TOTAL_RND", PType::PTYPE_TOTAL_RND);

  py::class_<SingleScatteringData>(m, "SingleScatteringData")
      .def(py::init<>())
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
      .def_readwrite("abs_vec_data", &SingleScatteringData::abs_vec_data);

  py::class_<ScatteringMetaData>(m, "ScatteringMetaData")
      .def(py::init<>())
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
                     &ScatteringMetaData::diameter_area_equ_aerodynamical);

  PythonInterfaceWorkspaceArray(ScatteringMetaData);
  PythonInterfaceWorkspaceArray(SingleScatteringData);
  PythonInterfaceWorkspaceArray(ArrayOfScatteringMetaData);
  PythonInterfaceWorkspaceArray(ArrayOfSingleScatteringData);
}
}  // namespace Python
