#include <python_interface.h>

#include "py_macros.h"

namespace Python {
void py_scattering_species(py::module_& m) try {
  //
  // Bulk property
  //

  artsclass<BulkProperty>(m, "BulkProperty")
      .def(py::init([]() { return std::make_shared<BulkProperty>(); }),
           "Default value")
      .def(py::init(
               [](const std::string& c) { return toBulkPropertyOrThrow(c); }),
           "From :class:`str`");

  //
  // Modified gamma PSD
  //

  py::implicitly_convertible<std::string, SpeciesEnum>();
  artsclass<ModifiedGamma>(m, "ModifiedGamma")
      .def(py::init([](std::string second_moment,
                       std::optional<Numeric> n0,
                       std::optional<Numeric> mu,
                       std::optional<Numeric> la,
                       std::optional<Numeric> ga,
                       std::optional<Numeric> t_min,
                       std::optional<Numeric> t_max) {
             return std::make_shared<ModifiedGamma>(
                 second_moment, n0, mu, la, ga, t_min, t_max);
           }),
           "A modified gamma PSD.",
           py::arg("second_moment"),
           py::arg("n0") = std::nullopt,
           py::arg("mu") = std::nullopt,
           py::arg("la") = std::nullopt,
           py::arg("ga") = std::nullopt,
           py::arg("t_min") = std::nullopt,
           py::arg("t_max") = std::nullopt);

  //
  // ScatSpeciesProperty
  //

  artsclass<ParticleHabit>(m, "ParticleHabit")
      .def(
          py::init([](std::string path) {
            return std::make_shared<ParticleHabit>(path);
          }),
          "A collection scattering properties for a range of particle of difference sizes.");

  artsclass<ScatteringHabit>(m, "ScatteringHabit")
      .def(
          py::init([](ParticleHabit particle_habit, PSD psd) {
            return std::make_shared<ScatteringHabit>(particle_habit, psd);
          }),
          "A particle habit with associated PSD representing a distribution of scattering particles in the atmosphere.");

  artsclass<HenyeyGreenstein>(m, "HenyeyGreenstein")
      .def(py::init(
               [](Numeric g) { return std::make_shared<HenyeyGreenstein>(g); }),
           "Henyey-Greenstein Phase Function")
      .def("evaluate_phase_function",
           static_cast<Vector (HenyeyGreenstein::*)(const Vector&)>(
               &HenyeyGreenstein::evaluate_phase_function),
           "Evaluate the Henyey-Greenstein phase function.")
      .def("evaluate_phase_function",
           static_cast<Numeric (HenyeyGreenstein::*)(const Numeric&)>(
               &HenyeyGreenstein::evaluate_phase_function),
           "Evaluate the Henyey-Greenstein phase function.")
      .PythonInterfaceCopyValue(HenyeyGreenstein)
      .PythonInterfaceWorkspaceVariableConversion(HenyeyGreenstein)
      .PythonInterfaceFileIO(HenyeyGreenstein)
      .PythonInterfaceBasicRepresentation(HenyeyGreenstein)
      .PythonInterfaceWorkspaceDocumentation(HenyeyGreenstein);
  artsclass<ArrayOfScatteringSpecies>(m, "ArrayOfScatteringSpecies")
      .def(py::init(
          []() { return std::make_shared<ArrayOfScatteringSpecies>(); }));

} catch (std::exception& e) {
  throw std::runtime_error(var_string(
      "DEV ERROR:\nCannot initialize scattering species:\n", e.what()));
};
}  // namespace Python
