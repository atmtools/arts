#include <arts_options.h>
#include <lineshapemodel.h>

#include "py_macros.h"
#include "py_auto_interface.h"
#include "python_interface.h"

//! See DeclareOption macro, but this may rename the python class
#define DeclareOptionRenamed(opt_rename, opt_namespace, opt_localname)       \
  py::class_<opt_namespace::opt_localname>(opt, #opt_rename)                 \
      .def(py::init([]() { return opt_namespace::opt_localname{}; }))        \
      .def(py::init([](const std::string& s) {                               \
        return opt_namespace::to##opt_localname##OrThrow(s);                 \
      }))                                                                    \
      .def_property(                                                         \
          "value",                                                           \
          [](const opt_namespace::opt_localname& x) { return toString(x); }, \
          [](opt_namespace::opt_localname& x, const std::string& s) {        \
            x = opt_namespace::to##opt_localname##OrThrow(s);                \
          })                                                                 \
      .PythonInterfaceCopyValue(opt_namespace::opt_localname)                \
      .PythonInterfaceBasicRepresentation(opt_namespace::opt_localname)      \
      .def(py::self == py::self)                                             \
      .def(py::self != py::self)                                             \
      .def(py::pickle(                                                       \
          [](const opt_namespace::opt_localname& t) {                        \
            return py::make_tuple(std::string(opt_namespace::toString(t)));  \
          },                                                                 \
          [](const py::tuple& t) {                                           \
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")              \
            return opt_namespace::opt_localname{                             \
                opt_namespace::to##opt_localname(t[0].cast<std::string>())}; \
          }))                                                                \
      .def_static(                                                           \
          "get_options",                                                     \
          []() { return opt_namespace::enumtyps::opt_localname##Types; })    \
      .def_static(                                                           \
          "get_options_as_strings",                                          \
          []() { return opt_namespace::enumstrs::opt_localname##Names; })    \
      .doc() = "Options for " #opt_rename;                                   \
  py::implicitly_convertible<std::string, opt_namespace::opt_localname>();

//! Exposes and option defined by the ARTS internal ENUMCLASS macro to pyarts
#define DeclareOption(opt_namespace, opt_localname) \
  DeclareOptionRenamed(opt_localname, opt_namespace, opt_localname)
namespace Python {
void py_options(py::module_& m) {
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
}
}  // namespace Python
