#include <arts_options.h>
#include <lineshapemodel.h>
#include <python_interface.h>

#include "atm.h"
#include <path_point.h>
#include "py_macros.h"

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
  }()

//! Exposes and option defined by the ARTS internal ENUMCLASS macro to pyarts
#define DeclareOption(opt_namespace, opt_localname) \
  DeclareOptionRenamed(opt_localname, opt_namespace, opt_localname)
namespace Python {
void py_options(py::module_& m) try {
  auto opt = m.def_submodule("options");
  opt.doc() = "Various named options of Arts";

  // Default multiple-choice options:
  DeclareOption(Options, planetOption);

  // Enum options relating to spectroscopy
  DeclareOptionRenamed(LineShapeTemperatureModel, LineShape, TemperatureModel);
  DeclareOptionRenamed(AbsorptionCutoffType, Absorption, CutoffType);
  DeclareOptionRenamed(AbsorptionMirroringType, Absorption, MirroringType);
  DeclareOptionRenamed(AbsorptionPopulationType, Absorption, PopulationType);
  DeclareOptionRenamed(
      AbsorptionNormalizationType, Absorption, NormalizationType);
  DeclareOptionRenamed(LineShapeType, LineShape, Type);
  DeclareOptionRenamed(LineShapeVariableOLDOLD, LineShape, Variable);

  // Atm
  DeclareOptionRenamed(AtmExtrapolation, Atm, Extrapolation);
  DeclareOptionRenamed(AtmKey, Atm, Key);

  // Surface
  DeclareOptionRenamed(SurfaceKey, Surf, Key);

  // Predef enums
  DeclareOptionRenamed(
      PredefinedModelDataKey, Absorption::PredefinedModel, DataKey);

  // Position types
  DeclareOptionRenamed(PathPositionType, path, PositionType);

  // Quantum numbers
  DeclareOptionRenamed(QuantumNumberType, Quantum::Number, Type);

  // Line shape types
  DeclareOptionRenamed(TemperatureModelType, lbl::temperature, model_type);
  DeclareOptionRenamed(LineShapeVariable, lbl::line_shape, variable);
  DeclareOptionRenamed(LineshapeNEWNEW, lbl, Lineshape);
  DeclareOptionRenamed(CutoffTypeNEWNEW, lbl, CutoffType);

  // Species enums
  DeclareOptionRenamed(SpeciesTagType, Species, TagType);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize options\n", e.what()));
}
}  // namespace Python
