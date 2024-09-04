#include <python_interface.h>

#include "hpy_arts.h"
#include "py_macros.h"

namespace Python {
void py_scattering_species(py::module_& m) try {
  //
  // ScatSpeciesProperty
  //

  py::class_<ScatteringSpeciesProperty> ssp(m, "ScatteringSpeciesProperty");


  //
  // Modified gamma PSD
  //

  py::class_<MGDSingleMoment>(m, "MGDSingleMoment");
  //py::arg("n0") = std::nullopt,
  //py::arg("mu") = std::nullopt,
  //py::arg("la") = std::nullopt,
  //py::arg("ga") = std::nullopt,
  //py::arg("t_min") = std::nullopt,
  //py::arg("t_max") = std::nullopt);

  py::class_<ParticleHabit>(m, "ParticleHabit");

  py::class_<ScatteringHabit>(m, "ScatteringHabit");

  py::class_<HenyeyGreenstein>(m, "HenyeyGreenstein")
      .def("evaluate_phase_function",
           static_cast<Vector (HenyeyGreenstein::*)(const Vector&)>(
               &HenyeyGreenstein::evaluate_phase_function),
           "Evaluate the Henyey-Greenstein phase function.")
      .def("evaluate_phase_function",
           static_cast<Numeric (HenyeyGreenstein::*)(const Numeric&)>(
               &HenyeyGreenstein::evaluate_phase_function),
           "Evaluate the Henyey-Greenstein phase function.");
  //.PythonInterfaceCopyValue(HenyeyGreenstein)
  //.PythonInterfaceWorkspaceVariableConversion(HenyeyGreenstein)
  //.PythonInterfaceFileIO(HenyeyGreenstein)
  //.PythonInterfaceBasicRepresentation(HenyeyGreenstein)
  //.PythonInterfaceWorkspaceDocumentation(HenyeyGreenstein);
  py::class_<ArrayOfScatteringSpecies>(m, "ArrayOfScatteringSpecies");

} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "DEV ERROR:\nCannot initialize scattering species:\n", e.what()));
};
}  // namespace Python
