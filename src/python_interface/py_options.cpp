#include <arts_options.h>
#include <lineshapemodel.h>
#include "atm.h"
#include "ppath_struct.h"

#include "py_macros.h"
#include <python_interface.h>

//! See DeclareOption macro, but this may rename the python class
#define DeclareOptionRenamed(opt_rename, opt_namespace, opt_localname)             \
  [&]() {                                                                          \
    auto cls =                                                                     \
        artsclass<opt_namespace::opt_localname>(opt, #opt_rename)                  \
            .def(py::init([]() { return opt_namespace::opt_localname{}; }),        \
                 "Default value")                                                  \
            .def(py::init([](const std::string& s) {                               \
                   return opt_namespace::to##opt_localname##OrThrow(s);            \
                 }),                                                               \
                 py::arg("str"),                                                   \
                 "From :class:`str`")                                              \
            .def("__hash__",                                                       \
                 [](opt_namespace::opt_localname& x) {                             \
                   return std::hash<opt_namespace::opt_localname>{}(x);            \
                 })                                                                \
            .PythonInterfaceCopyValue(opt_namespace::opt_localname)                \
            .PythonInterfaceBasicRepresentation(opt_namespace::opt_localname)      \
            .def(py::self == py::self)                                             \
            .def(py::self != py::self)                                             \
            .def(py::pickle(                                                       \
                [](const opt_namespace::opt_localname& t) {                        \
                  return py::make_tuple(                                           \
                      std::string(opt_namespace::toString(t)));                    \
                },                                                                 \
                [](const py::tuple& t) {                                           \
                  ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")              \
                  return opt_namespace::opt_localname{                             \
                      opt_namespace::to##opt_localname##OrThrow(                   \
                          t[0].cast<std::string>())};                              \
                }))                                                                \
            .def_static(                                                           \
                "get_options",                                                     \
                []() {                                                             \
                  return opt_namespace::enumtyps::opt_localname##Types;            \
                },                                                                 \
                py::doc(":class:`list` of full set of options available"))         \
            .def_static(                                                           \
                "get_options_as_strings",                                          \
                []() {                                                             \
                  return opt_namespace::enumstrs::opt_localname##Names;            \
                },                                                                 \
                py::doc(                                                           \
                    ":class:`list` of full set of options available as strings")); \
    cls.doc() = "Options for " #opt_rename;                                        \
    for (auto& x : opt_namespace::enumtyps::opt_localname##Types) {                \
      String str = String{toString(x)};                                            \
      if (str == "None") str += "_";                                               \
      cls.def_property_readonly_static(                                            \
          str.c_str(),                                                             \
          [x](const py::object&) { return x; },                                    \
          py::doc(":class:`~pyarts.options." #opt_rename                           \
                  "` static value as named"));                                     \
    }                                                                              \
    py::implicitly_convertible<std::string, opt_namespace::opt_localname>();       \
    return cls;                                                                    \
  }();

//! Exposes and option defined by the ARTS internal ENUMCLASS macro to pyarts
#define DeclareOption(opt_namespace, opt_localname) \
  DeclareOptionRenamed(opt_localname, opt_namespace, opt_localname)
namespace Python {
void py_options(py::module_& m) try {
  auto opt = m.def_submodule("options");
  opt.doc() = "Various named options of Arts";

  // Default agenda options:
  DeclareOption(Options, iy_main_agendaDefaultOptions)
  DeclareOption(Options, iy_loop_freqs_agendaDefaultOptions)
  DeclareOption(Options, iy_space_agendaDefaultOptions)
  DeclareOption(Options, iy_surface_agendaDefaultOptions)
  DeclareOption(Options, iy_cloudbox_agendaDefaultOptions)
  DeclareOption(Options, ppath_agendaDefaultOptions)
  DeclareOption(Options, ppath_step_agendaDefaultOptions)
  DeclareOption(Options, refr_index_air_agendaDefaultOptions)
  DeclareOption(Options, water_p_eq_agendaDefaultOptions)
  DeclareOption(Options, gas_scattering_agendaDefaultOptions)
  DeclareOption(Options, surface_rtprop_agendaDefaultOptions)
  DeclareOption(Options, g0_agendaDefaultOptions)
  DeclareOption(Options, test_agendaDefaultOptions)
  DeclareOption(Options, dobatch_calc_agendaDefaultOptions)
  DeclareOption(Options, ybatch_calc_agendaDefaultOptions)
  DeclareOption(Options, spt_calc_agendaDefaultOptions)
  DeclareOption(Options, sensor_response_agendaDefaultOptions)
  DeclareOption(Options, propmat_clearsky_agendaDefaultOptions)
  DeclareOption(Options, pha_mat_spt_agendaDefaultOptions)
  DeclareOption(Options, met_profile_calc_agendaDefaultOptions)
  DeclareOption(Options, main_agendaDefaultOptions)
  DeclareOption(Options, jacobian_agendaDefaultOptions)
  DeclareOption(Options, iy_radar_agendaDefaultOptions)
  DeclareOption(Options, iy_independent_beam_approx_agendaDefaultOptions)
  DeclareOption(Options, inversion_iterate_agendaDefaultOptions)
  DeclareOption(Options, forloop_agendaDefaultOptions)
  DeclareOption(Options, doit_scat_field_agendaDefaultOptions)
  DeclareOption(Options, doit_rte_agendaDefaultOptions)
  DeclareOption(Options, doit_mono_agendaDefaultOptions)
  DeclareOption(Options, doit_conv_test_agendaDefaultOptions)
  DeclareOption(Options, ppvar_rtprop_agendaDefaultOptions)
  DeclareOption(Options, rte_background_agendaDefaultOptions)

  // Default multiple-choice options:
  DeclareOption(Options, planetDefaultOptions)

  // Enum options relating to spectroscopy
  DeclareOptionRenamed(LineShapeTemperatureModel, LineShape, TemperatureModel)
  DeclareOptionRenamed(AbsorptionCutoffType, Absorption, CutoffType)
  DeclareOptionRenamed(AbsorptionMirroringType, Absorption, MirroringType)
  DeclareOptionRenamed(AbsorptionPopulationType, Absorption,PopulationType)
  DeclareOptionRenamed(AbsorptionNormalizationType, Absorption, NormalizationType)
  DeclareOptionRenamed(LineShapeType, LineShape, Type)
  DeclareOptionRenamed(LineShapeVariable, LineShape, Variable)

  // Ppath
  DeclareOption(Options, PpathBackground)

  // Atm
  DeclareOptionRenamed(AtmExtrapolation, Atm, Extrapolation)
  DeclareOptionRenamed(AtmKey, Atm, Key)


  // Jacobian enums
  DeclareOptionRenamed(JacobianType, Jacobian, Type)
  DeclareOptionRenamed(JacobianAtm, Jacobian, Atm)
  DeclareOptionRenamed(JacobianLine, Jacobian, Line)
  DeclareOptionRenamed(JacobianSensor, Jacobian, Sensor)
  DeclareOptionRenamed(JacobianSpecial, Jacobian, Special)

  // Predef enums
  DeclareOptionRenamed(PredefinedModelDataKey, Absorption::PredefinedModel, DataKey)

  // Predef enums
  DeclareOptionRenamed(QuantumNumberType, Quantum::Number, Type)

  // Species enums
  DeclareOptionRenamed(SpeciesTagType, Species, TagType)
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize options\n", e.what()));
}
}  // namespace Python
