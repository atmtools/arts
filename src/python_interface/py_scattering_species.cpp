#include <python_interface.h>

#include "py_macros.h"

namespace Python {
void py_scattering_species(py::module_& m) try {
  //
  // ScatSpeciesProperty
  //

  artsclass<ParticulateProperty>(m, "ParticulateProperty")
      .def(py::init([]() { return std::make_shared<ParticulateProperty>(); }),
           "Default value")
      .def(py::init([](const std::string& c) {
             return toParticulatePropertyOrThrow(c);
           }),
           "From :class:`str`");

  artsclass<ScatteringSpeciesProperty>(m, "ScatteringSpeciesProperty")
      .def(py::init(
               []() { return std::make_shared<ScatteringSpeciesProperty>(); }),
           "Default value")
      .def(py::init([](const std::string& species,
                       const ParticulateProperty& pprop) {
             return ScatteringSpeciesProperty(species, pprop);
           }),
           "From :class:`str`");

  //
  // Modified gamma PSD
  //

  artsclass<MGDSingleMoment>(m, "MGDSingleMoment")
      .def(py::init([](ScatteringSpeciesProperty moment,
                       Numeric n_alpha,
                       Numeric n_b,
                       Numeric mu,
                       Numeric gamma,
                       Numeric t_min,
                       Numeric t_max,
                       bool picky) {
             return std::make_shared<MGDSingleMoment>(
                 moment, n_alpha, n_b, mu, gamma, t_min, t_max, picky);
           }),
           "A modified gamma PSD.")
      .def(py::init([](ScatteringSpeciesProperty moment,
                       std::string name,
                       Numeric t_min,
                       Numeric t_max,
                       bool picky) {
             return std::make_shared<MGDSingleMoment>(
                 moment, name, t_min, t_max, picky);
           }),
           "A modified gamma PSD.")
      .def("evaluate",
           [](const MGDSingleMoment& psd,
              const AtmPoint& point,
              const Vector& sizes,
              Numeric scat_species_a,
              Numeric scat_species_b) {
             return psd.evaluate(point, sizes, scat_species_a, scat_species_b);
           });
  //py::arg("n0") = std::nullopt,
  //py::arg("mu") = std::nullopt,
  //py::arg("la") = std::nullopt,
  //py::arg("ga") = std::nullopt,
  //py::arg("t_min") = std::nullopt,
  //py::arg("t_max") = std::nullopt);

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
