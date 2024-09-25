#include <lineshapemodel.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/bind_vector.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <python_interface.h>
#include <zeemandata.h>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "absorptionlines.h"
#include "cia.h"
#include "debug.h"
#include "hpy_arts.h"
#include "hpy_vector.h"
#include "isotopologues.h"
#include "physics_funcs.h"
#include "py_macros.h"
#include "quantum_numbers.h"
#include "species_tags.h"

NB_MAKE_OPAQUE(Zeeman::Polarization)

namespace Python {
void py_spectroscopy(py::module_& m) try {
  static_assert(LineShapeModelParameters::N == 4);
  py::class_<LineShapeModelParameters>(m, "LineShapeModelParameters")
      .def(py::init<LineShapeTemperatureModelOld,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric>(),
           "type"_a = LineShapeTemperatureModelOld::None,
           "X0"_a   = 0,
           "X1"_a   = 0,
           "X2"_a   = 0,
           "X3"_a   = 0,
           "From values")
      .PythonInterfaceCopyValue(LineShapeModelParameters)
      .def_rw(
          "type",
          &LineShapeModelParameters::type,
          ":class:`~pyarts.arts.options.LineShapeTemperatureModel` The temperature model")
      .def_rw(
          "X0", &LineShapeModelParameters::X0, ":class:`float` 1st coefficient")
      .def_rw(
          "X1", &LineShapeModelParameters::X1, ":class:`float` 2nd coefficient")
      .def_rw(
          "X2", &LineShapeModelParameters::X2, ":class:`float` 3rd coefficient")
      .def_rw(
          "X3", &LineShapeModelParameters::X3, ":class:`float` 4th coefficient")
      .PythonInterfaceBasicRepresentation(LineShapeModelParameters)
      .def("__getstate__",
           [](const LineShapeModelParameters& t) {
             return std::tuple<LineShapeTemperatureModelOld,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric>(t.type, t.X0, t.X1, t.X2, t.X3);
           })
      .def("__setstate__",
           [](LineShapeModelParameters* lsmp,
              const std::tuple<LineShapeTemperatureModelOld,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric>& state) {
             new (lsmp) LineShapeModelParameters{std::get<0>(state),
                                                 std::get<1>(state),
                                                 std::get<2>(state),
                                                 std::get<3>(state),
                                                 std::get<4>(state)};
           })
      .doc() = "A temperature model calculator";

  auto ZeemanPolarizationStringGetter =
      [](Zeeman::Polarization c) -> std::string {
    if (c == Zeeman::Polarization::SigmaMinus) return "SigmaMinus";
    if (c == Zeeman::Polarization::Pi) return "Pi";
    if (c == Zeeman::Polarization::SigmaPlus) return "SigmaPlus";
    if (c == Zeeman::Polarization::None) return "None";
    ARTS_USER_ERROR("Bad enum state")
  };
  auto ZeemanPolarizationEnumGetter =
      [](const std::string& c) -> Zeeman::Polarization {
    if (c == "SigmaMinus") return Zeeman::Polarization::SigmaMinus;
    if (c == "Pi") return Zeeman::Polarization::Pi;
    if (c == "SigmaPlus") return Zeeman::Polarization::SigmaPlus;
    if (c == "None") return Zeeman::Polarization::None;
    ARTS_USER_ERROR("Bad enum value {}", c);
  };
  py::class_<Zeeman::Polarization>(m, "ZeemanPolarization")
      .def(
          "__init__",
          [ZeemanPolarizationEnumGetter](Zeeman::Polarization* z,
                                         const std::string& c) {
            new (z) Zeeman::Polarization{ZeemanPolarizationEnumGetter(c)};
          },
          "str"_a = std::string{"None"},
          "From :class:`str`")
      // .PythonInterfaceCopyValue(Zeeman::Polarization)
      .def("__repr__",
           [ZeemanPolarizationStringGetter](const Zeeman::Polarization& c) {
             return ZeemanPolarizationStringGetter(c);
           })
      .def("__str__",
           [ZeemanPolarizationStringGetter](const Zeeman::Polarization& c) {
             return ZeemanPolarizationStringGetter(c);
           })
      .def("__getstate__",
           [ZeemanPolarizationStringGetter](const Zeeman::Polarization& c) {
             return std::tuple<std::string>(ZeemanPolarizationStringGetter(c));
           })
      .def("__setstate__",
           [ZeemanPolarizationEnumGetter](
               Zeeman::Polarization* z, const std::tuple<std::string>& state) {
             new (z) Zeeman::Polarization{
                 ZeemanPolarizationEnumGetter(std::get<0>(state))};
           })
      .doc() = "Options for ZeemanPolarization";
  py::implicitly_convertible<std::string, Zeeman::Polarization>();

  py::class_<Zeeman::Model>(m, "ZeemanModel")
      .def(py::init<>())
      .def(py::init<Numeric, Numeric>(),
           "gu"_a,
           "gl"_a,
           "From two numeric values")
      .def(
          "__init__",
          [](Zeeman::Model* z, std::array<Numeric, 2> a) {
            new (z) Zeeman::Model(a[0], a[1]);
          },
          "gs"_a,
          "From list of two values")
      .def_prop_rw(
          "gu",
          [](Zeeman::Model z) { return z.gu(); },
          [](Zeeman::Model& z, Numeric g) { z.gu() = g; },
          ":class:`~float` The upper state g")
      .def_prop_rw(
          "gl",
          [](Zeeman::Model z) { return z.gl(); },
          [](Zeeman::Model& z, Numeric g) { z.gl() = g; },
          ":class:`~float` The lower state g")
      .def("__getstate__",
           [](const Zeeman::Model& t) {
             return std::tuple<Numeric, Numeric>(t.gu(), t.gl());
           })
      .def("__setstate__",
           [](Zeeman::Model* z, const std::tuple<Numeric, Numeric>& state) {
             new (z) Zeeman::Model(std::get<0>(state), std::get<1>(state));
           })
      .doc() = "A Zeeman model";
  py::implicitly_convertible<std::array<Numeric, 2>, Zeeman::Model>();

  py::class_<LineShape::Output>(m, "LineShapeOutput")
      .def(py::init<Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric>(),
           "g0"_a  = 0,
           "d0"_a  = 0,
           "g2"_a  = 0,
           "d2"_a  = 0,
           "fvc"_a = 0,
           "eta"_a = 0,
           "y"_a   = 0,
           "g"_a   = 0,
           "dv"_a  = 0,
           "From values")
      .PythonInterfaceCopyValue(LineShape::Output)
      .def_rw("G0",
              &LineShape::Output::G0,
              ":class:`float`: Pressure broadening speed-independent")
      .def_rw("D0",
              &LineShape::Output::D0,
              ":class:`float`: Pressure f-shifting speed-independent")
      .def_rw("G2",
              &LineShape::Output::G2,
              ":class:`float`: Pressure broadening speed-dependent")
      .def_rw("D2",
              &LineShape::Output::D2,
              ":class:`float`: Pressure f-shifting speed-dependent")
      .def_rw("ETA", &LineShape::Output::ETA, ":class:`float`: Correlation")
      .def_rw("FVC",
              &LineShape::Output::FVC,
              ":class:`float`: Frequency of velocity-changing collisions")
      .def_rw("Y",
              &LineShape::Output::Y,
              ":class:`float`: First order line mixing coefficient")
      .def_rw("G",
              &LineShape::Output::G,
              ":class:`float`: Second order line mixing coefficient")
      .def_rw("DV",
              &LineShape::Output::DV,
              ":class:`float`: Second order line mixing f-shifting")
      // .PythonInterfaceBasicRepresentation(LineShape::Output)
      .def("__getstate__",
           [](const LineShape::Output& t) {
             return std::tuple<Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric>(
                 t.G0, t.D0, t.G2, t.D2, t.FVC, t.ETA, t.Y, t.G, t.DV);
           })
      .def("__setstate__",
           [](LineShape::Output* lo,
              const std::tuple<Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric>& state) {
             new (lo) LineShape::Output(std::get<0>(state),
                                        std::get<1>(state),
                                        std::get<2>(state),
                                        std::get<3>(state),
                                        std::get<4>(state),
                                        std::get<5>(state),
                                        std::get<6>(state),
                                        std::get<7>(state),
                                        std::get<8>(state));
           })
      .doc() = "Derived line shape parameters";

  py::class_<LineShapeSingleSpeciesModel>(m, "LineShapeSingleSpeciesModel")
      .def(py::init<LineShape::ModelParameters,
                    LineShape::ModelParameters,
                    LineShape::ModelParameters,
                    LineShape::ModelParameters,
                    LineShape::ModelParameters,
                    LineShape::ModelParameters,
                    LineShape::ModelParameters,
                    LineShape::ModelParameters,
                    LineShape::ModelParameters>(),
           "G0"_a  = LineShape::ModelParameters{},
           "D0"_a  = LineShape::ModelParameters{},
           "G2"_a  = LineShape::ModelParameters{},
           "D2"_a  = LineShape::ModelParameters{},
           "FVC"_a = LineShape::ModelParameters{},
           "ETA"_a = LineShape::ModelParameters{},
           "Y"_a   = LineShape::ModelParameters{},
           "G"_a   = LineShape::ModelParameters{},
           "DV"_a  = LineShape::ModelParameters{},
           "From values")
      .PythonInterfaceCopyValue(LineShapeSingleSpeciesModel)
      .PythonInterfaceBasicRepresentation(LineShapeSingleSpeciesModel)
      .def_prop_rw(
          "G0",
          [](const LineShapeSingleSpeciesModel& x) { return x.G0(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.G0() = y;
          },
          ":class:`~pyarts.arts.LineShapeModelParameters`: Pressure broadening speed-independent")
      .def_prop_rw(
          "D0",
          [](const LineShapeSingleSpeciesModel& x) { return x.D0(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.D0() = y;
          },
          ":class:`~pyarts.arts.LineShapeModelParameters`: Pressure f-shifting speed-independent")
      .def_prop_rw(
          "G2",
          [](const LineShapeSingleSpeciesModel& x) { return x.G2(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.G2() = y;
          },
          ":class:`~pyarts.arts.LineShapeModelParameters`: Pressure broadening speed-dependent")
      .def_prop_rw(
          "D2",
          [](const LineShapeSingleSpeciesModel& x) { return x.D2(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.D2() = y;
          },
          ":class:`~pyarts.arts.LineShapeModelParameters`: Pressure f-shifting speed-dependent")
      .def_prop_rw(
          "FVC",
          [](const LineShapeSingleSpeciesModel& x) { return x.FVC(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.FVC() = y;
          },
          ":class:`~pyarts.arts.LineShapeModelParameters`: Frequency of velocity-changing collisions")
      .def_prop_rw(
          "ETA",
          [](const LineShapeSingleSpeciesModel& x) { return x.ETA(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.ETA() = y;
          },
          ":class:`~pyarts.arts.LineShapeModelParameters`: Correlation")
      .def_prop_rw(
          "Y",
          [](const LineShapeSingleSpeciesModel& x) { return x.Y(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.Y() = y;
          },
          ":class:`~pyarts.arts.LineShapeModelParameters`: First order line mixing coefficient")
      .def_prop_rw(
          "G",
          [](const LineShapeSingleSpeciesModel& x) { return x.G(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.G() = y;
          },
          ":class:`~pyarts.arts.LineShapeModelParameters`: Second order line mixing coefficient")
      .def_prop_rw(
          "DV",
          [](const LineShapeSingleSpeciesModel& x) { return x.DV(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.DV() = y;
          },
          ":class:`~pyarts.arts.LineShapeModelParameters`: Second order line mixing f-shifting")
      .def("__getstate__",
           [](const LineShapeSingleSpeciesModel& t) {
             return std::tuple<LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters>(t.G0(),
                                                           t.D0(),
                                                           t.G2(),
                                                           t.D2(),
                                                           t.FVC(),
                                                           t.ETA(),
                                                           t.Y(),
                                                           t.G(),
                                                           t.DV());
           })
      .def("__setstate__",
           [](LineShapeSingleSpeciesModel* l,
              const std::tuple<LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters,
                               LineShape::ModelParameters>& state) {
             new (l) LineShapeSingleSpeciesModel(std::get<0>(state),
                                                 std::get<1>(state),
                                                 std::get<2>(state),
                                                 std::get<3>(state),
                                                 std::get<4>(state),
                                                 std::get<5>(state),
                                                 std::get<6>(state),
                                                 std::get<7>(state),
                                                 std::get<8>(state));
           })
      .doc() = "Single species line shape model";

  py::class_<LineShapeModel>(m, "LineShapeModel")
      .def(py::init<>())
      .def(py::init<std::vector<LineShapeSingleSpeciesModel>>())
      .PythonInterfaceCopyValue(LineShapeModel)
      .PythonInterfaceBasicRepresentation(LineShapeModel)
      .PythonInterfaceIndexItemAccess(LineShapeModel)
      .def_prop_rw(
          "data",
          [](const LineShapeModel& x) { return x.Data(); },
          [](LineShapeModel& x, std::vector<LineShapeSingleSpeciesModel> y) {
            x.Data() = std::move(y);
          },
          ":class:`list` of single species models")
      .def("__getstate__",
           [](const LineShapeModel& t) {
             return std::tuple<std::vector<LineShapeSingleSpeciesModel>>(
                 t.Data());
           })
      .def("__setstate__",
           [](LineShapeModel* l,
              const std::tuple<std::vector<LineShapeSingleSpeciesModel>>&
                  state) { new (l) LineShapeModel(std::get<0>(state)); })
      .doc() = "Multi-species line shape model";
  py::implicitly_convertible<std::vector<LineShapeSingleSpeciesModel>,
                             LineShapeModel>();

  py::class_<AbsorptionSingleLine>(m, "AbsorptionSingleLine")
      .def(py::init<Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    Numeric,
                    ZeemanModel,
                    LineShapeModel,
                    Quantum::Number::LocalState>(),
           "F0"_a          = 0,
           "I0"_a          = 0,
           "E0"_a          = 0,
           "glow"_a        = 0,
           "gupp"_a        = 0,
           "A"_a           = 0,
           "zeeman"_a      = Zeeman::Model(),
           "lineshape"_a   = LineShape::Model(),
           "localquanta"_a = Quantum::Number::LocalState{},
           "From values")
      .PythonInterfaceCopyValue(AbsorptionSingleLine)
      .PythonInterfaceBasicRepresentation(AbsorptionSingleLine)
      .def_rw(
          "F0", &AbsorptionSingleLine::F0, ":class:`float` Central frequency")
      .def_rw(
          "I0", &AbsorptionSingleLine::I0, ":class:`float` Reference intensity")
      .def_rw(
          "E0", &AbsorptionSingleLine::E0, ":class:`float` Lower energy state")
      .def_rw("glow",
              &AbsorptionSingleLine::glow,
              ":class:`float` Lower level statistical weight")
      .def_rw("gupp",
              &AbsorptionSingleLine::gupp,
              ":class:`float` Upper level statistical weight")
      .def_rw("A",
              &AbsorptionSingleLine::A,
              ":class:`float` Einstein spontaneous emission coefficient")
      .def_rw("zeeman",
              &AbsorptionSingleLine::zeeman,
              ":class:`pyarts.arts.ZeemanModel` Zeeman model")
      .def_rw("lineshape",
              &AbsorptionSingleLine::lineshape,
              ":class:`pyarts.arts.LineShapeModel` Line shape model")
      .def_rw(
          "localquanta",
          &AbsorptionSingleLine::localquanta,
          ":class:`pyarts.arts.QuantumNumberLocalState` Local quantum numbers")
      .def("__getstate__",
           [](const AbsorptionSingleLine& t) {
             return std::tuple<Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Zeeman::Model,
                               LineShape::Model,
                               Quantum::Number::LocalState>(t.F0,
                                                            t.I0,
                                                            t.E0,
                                                            t.glow,
                                                            t.gupp,
                                                            t.A,
                                                            t.zeeman,
                                                            t.lineshape,
                                                            t.localquanta);
           })
      .def("__setstate__",
           [](AbsorptionSingleLine* a,
              const std::tuple<Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Numeric,
                               Zeeman::Model,
                               LineShape::Model,
                               Quantum::Number::LocalState>& state) {
             new (a) AbsorptionSingleLine(std::get<0>(state),
                                          std::get<1>(state),
                                          std::get<2>(state),
                                          std::get<3>(state),
                                          std::get<4>(state),
                                          std::get<5>(state),
                                          std::get<6>(state),
                                          std::get<7>(state),
                                          std::get<8>(state));
           })
      .doc() = "Single absorption line";

  py::bind_vector<Array<AbsorptionSingleLine>,
                  py::rv_policy::reference_internal>(
      m, "ArrayOfAbsorptionSingleLine")
      .doc() = "List of :class:`~pyarts.arts.AbsorptionSingleLine`";

  py::class_<AbsorptionLines> al(m, "AbsorptionLines");
  workspace_group_interface(al);
  al.def(
        "__init__",
        [](AbsorptionLines* a,
           bool selfbroadening,
           bool bathbroadening,
           AbsorptionCutoffTypeOld cutoff,
           AbsorptionMirroringTypeOld mirroring,
           AbsorptionPopulationTypeOld population,
           AbsorptionNormalizationTypeOld normalization,
           LineShapeTypeOld lineshapetype,
           Numeric T0,
           Numeric cutofffreq,
           Numeric linemixinglimit,
           QuantumIdentifier quantumidentity,
           ArrayOfSpeciesEnum broadeningspecies,
           Array<AbsorptionSingleLine> lines) {
          ARTS_USER_ERROR_IF(not good_enum(quantumidentity.Species()),
                             "Bad quantumidentity, must specify species and "
                             "isotoplogue, e.g., H2O-161")
          ARTS_USER_ERROR_IF(
              quantumidentity.Isotopologue().joker() or
                  Species::is_predefined_model(quantumidentity.Isotopologue()),
              "Bad quantumidentity, must specify isotopologue, e.g., H2O-161")
          ARTS_USER_ERROR_IF(
              std::any_of(lines.begin(),
                          lines.end(),
                          [count = broadeningspecies.size()](auto& line) {
                            return line.lineshape.size() not_eq count;
                          }),
              "Incorrect size of broadeningspecies vs the line shape model")
          ARTS_USER_ERROR_IF(
              broadeningspecies.size() <
                  static_cast<std::size_t>(selfbroadening + bathbroadening),
              "Must have atleast {}"
              " broadening species to support settings",
              (selfbroadening + bathbroadening))
          new (a) AbsorptionLines(selfbroadening,
                                  bathbroadening,
                                  cutoff,
                                  mirroring,
                                  population,
                                  normalization,
                                  lineshapetype,
                                  T0,
                                  cutofffreq,
                                  linemixinglimit,
                                  std::move(quantumidentity),
                                  std::move(broadeningspecies),
                                  std::move(lines));
        },
        "selfbroadening"_a    = false,
        "bathbroadening"_a    = false,
        "cutoff"_a            = AbsorptionCutoffTypeOld::None,
        "mirroring"_a         = AbsorptionMirroringTypeOld::None,
        "population"_a        = AbsorptionPopulationTypeOld::LTE,
        "normalization"_a     = AbsorptionNormalizationTypeOld::None,
        "lineshapetype"_a     = LineShapeTypeOld::DP,
        "T0"_a                = 296,
        "cutofffreq"_a        = -1,
        "linemixinglimit"_a   = -1,
        "quantumidentity"_a   = QuantumIdentifier("H2O-161"),
        "broadeningspecies"_a = ArrayOfSpeciesEnum{},
        "lines"_a             = Array<AbsorptionSingleLine>{},
        "From values")
      .def_rw("selfbroadening",
              &AbsorptionLines::selfbroadening,
              ":class:`bool` Does the line broadening have self broadening?")
      .def_rw("bathbroadening",
              &AbsorptionLines::bathbroadening,
              ":class:`bool` Does the line broadening have bath broadening?")
      .def_rw(
          "cutoff",
          &AbsorptionLines::cutoff,
          ":class:`~pyarts.arts.options.AbsorptionCutoffTypeOld` Cutoff type")
      .def_rw(
          "mirroring",
          &AbsorptionLines::mirroring,
          ":class:`~pyarts.arts.options.AbsorptionMirroringype` Mirroring type")
      .def_rw(
          "population",
          &AbsorptionLines::population,
          ":class:`~pyarts.arts.options.AbsorptionPopulationTypeOld` Line population distribution")
      .def_rw(
          "normalization",
          &AbsorptionLines::normalization,
          ":class:`~pyarts.arts.options.AbsorptionNormalizationTypeOld` Normalization type")
      .def_rw("lineshapetype",
              &AbsorptionLines::lineshapetype,
              ":class:`~pyarts.arts.options.LineShapeTypeOld` Line shape type")
      .def_rw(
          "T0",
          &AbsorptionLines::T0,
          ":class:`float` Reference temperature for all parameters of the lines")
      .def_rw("cutofffreq",
              &AbsorptionLines::cutofffreq,
              ":class:`float` Cutoff frequency")
      .def_rw("linemixinglimit",
              &AbsorptionLines::linemixinglimit,
              ":class:`float` Linemixing limit")
      .def_rw("quantumidentity",
              &AbsorptionLines::quantumidentity,
              ":class:`~pyarts.arts.QuantumIdentifier` Catalog ID")
      .def_rw(
          "broadeningspecies",
          &AbsorptionLines::broadeningspecies,
          ":class:`~pyarts.arts.ArrayOfSpecies` A list of broadening specie")
      .def_rw(
          "lines",
          &AbsorptionLines::lines,
          ":class:`~pyarts.arts.AbsorptionSingleLine` A list of individual lines")
      .def_prop_ro(
          "ok",
          &AbsorptionLines::OK,
          R"(:class:`bool` If False, the catalog cannot be used for any calculations)")
      .def_prop_ro("meta_data",
                   &AbsorptionLines::MetaData,
                   R"(:class:`~pyarts.arts.String` Catalog meta data string)")
      .def(
          "LineShapeOutput",
          [](AbsorptionLines& band,
             Size line,
             Numeric T,
             Numeric P,
             const Vector& VMR) {
            ARTS_USER_ERROR_IF(not band.OK(), "Band in bad shape")
            ARTS_USER_ERROR_IF(not(T > 0) or not(P >= 0),
                               "Bad atmospheric state (T P): ",
                               Vector{T, P})
            ARTS_USER_ERROR_IF(
                static_cast<Size>(VMR.size()) not_eq
                    band.broadeningspecies.size(),
                "Mismatch between VMRs and broadening species.\nVMR: {:B,}"
                "\nSpecies: {:B,}",
                VMR,
                band.broadeningspecies)
            if (line >= band.lines.size())
              throw std::out_of_range("Line index out of range");

            return band.ShapeParameters(line, T, P, VMR);
          },
          "line"_a,
          "T"_a,
          "P"_a,
          "VMR"_a,

          R"--(Computes the line shape paramters for the given atmospheric state

Note that the normalization assumes sum(VMR) is 1 for good results but does not enforce it

Parameters
----------
line : int
    Line index
T : float
    Temperature
P : float
    Pressure
VMR : ~pyarts.arts.Vector
    List of VMR values, must match length of broadening species

Returns
-------
X : ~pyarts.arts.LineShapeOutput
    The computed line shape parameters
)--")
      .def("__getstate__",
           [](const AbsorptionLines& t) {
             return std::tuple<bool,
                               bool,
                               AbsorptionCutoffTypeOld,
                               AbsorptionMirroringTypeOld,
                               AbsorptionPopulationTypeOld,
                               AbsorptionNormalizationTypeOld,
                               LineShapeTypeOld,
                               Numeric,
                               Numeric,
                               Numeric,
                               QuantumIdentifier,
                               ArrayOfSpeciesEnum,
                               Array<AbsorptionSingleLine>>(t.selfbroadening,
                                                            t.bathbroadening,
                                                            t.cutoff,
                                                            t.mirroring,
                                                            t.population,
                                                            t.normalization,
                                                            t.lineshapetype,
                                                            t.T0,
                                                            t.cutofffreq,
                                                            t.linemixinglimit,
                                                            t.quantumidentity,
                                                            t.broadeningspecies,
                                                            t.lines);
           })
      .def("__setstate__",
           [](AbsorptionLines* a,
              const std::tuple<bool,
                               bool,
                               AbsorptionCutoffTypeOld,
                               AbsorptionMirroringTypeOld,
                               AbsorptionPopulationTypeOld,
                               AbsorptionNormalizationTypeOld,
                               LineShapeTypeOld,
                               Numeric,
                               Numeric,
                               Numeric,
                               QuantumIdentifier,
                               ArrayOfSpeciesEnum,
                               Array<AbsorptionSingleLine>>& state) {
             new (a) AbsorptionLines(std::get<0>(state),
                                     std::get<1>(state),
                                     std::get<2>(state),
                                     std::get<3>(state),
                                     std::get<4>(state),
                                     std::get<5>(state),
                                     std::get<6>(state),
                                     std::get<7>(state),
                                     std::get<8>(state),
                                     std::get<9>(state),
                                     std::get<10>(state),
                                     std::get<11>(state),
                                     std::get<12>(state));
           });

  auto a1 =
      py::bind_vector<ArrayOfAbsorptionLines,
                      py::rv_policy::reference_internal>(
          m, "ArrayOfAbsorptionLines")
          .def(
              "fuzzy_find_all",
              [](const ArrayOfAbsorptionLines& a, const QuantumIdentifier& q) {
                return fuzzy_find_all(a, q);
              },
              "q"_a,
              "Find all the indexes that could match the given quantum identifier");
  workspace_group_interface(a1);
  vector_interface(a1);

  auto a2 = py::bind_vector<ArrayOfArrayOfAbsorptionLines,
                            py::rv_policy::reference_internal>(
      m, "ArrayOfArrayOfAbsorptionLines");
  workspace_group_interface(a2);
  vector_interface(a2);
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize spectroscopy\n", e.what()));
}
}  // namespace Python
