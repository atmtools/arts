#include <algorithm>
#include <lineshapemodel.h>
#include <python_interface.h>

#include <zeemandata.h>

#include <memory>
#include <string>
#include <utility>

#include "absorptionlines.h"
#include "cia.h"
#include "debug.h"
#include "isotopologues.h"
#include "physics_funcs.h"
#include "py_macros.h"
#include "quantum_numbers.h"
#include "species_tags.h"

namespace Python {
void py_spectroscopy(py::module_& m) try {
  static_assert(LineShapeModelParameters::N == 4);
  artsclass<LineShapeModelParameters>(m, "LineShapeModelParameters")
      .def(py::init([](LineShapeTemperatureModelOld a,
                       Numeric b,
                       Numeric c,
                       Numeric d,
                       Numeric e) {
             return std::make_shared<LineShapeModelParameters>(a, b, c, d, e);
           }),
           py::arg("type") = LineShapeTemperatureModelOld::None,
           py::arg("X0") = 0,
           py::arg("X1") = 0,
           py::arg("X2") = 0,
           py::arg("X3") = 0, "From values")
      .PythonInterfaceCopyValue(LineShapeModelParameters)
      .def_readwrite("type", &LineShapeModelParameters::type, ":class:`~pyarts.arts.options.LineShapeTemperatureModel` The temperature model")
      .def_readwrite("X0", &LineShapeModelParameters::X0, ":class:`float` 1st coefficient")
      .def_readwrite("X1", &LineShapeModelParameters::X1, ":class:`float` 2nd coefficient")
      .def_readwrite("X2", &LineShapeModelParameters::X2, ":class:`float` 3rd coefficient")
      .def_readwrite("X3", &LineShapeModelParameters::X3, ":class:`float` 4th coefficient")
      .PythonInterfaceBasicRepresentation(LineShapeModelParameters)
      .def(py::pickle(
          [](const LineShapeModelParameters& t) {
            return py::make_tuple(t.type, t.X0, t.X1, t.X2, t.X3);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 5, "Invalid state!")
            return std::make_shared<LineShapeModelParameters>(
                t[0].cast<LineShapeTemperatureModelOld>(),
                t[1].cast<Numeric>(),
                t[2].cast<Numeric>(),
                t[3].cast<Numeric>(),
                t[4].cast<Numeric>());
          })).doc() = "A temperature model calculator";

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
    ARTS_USER_ERROR("Bad enum value ", c);
  };
  artsclass<Zeeman::Polarization>(m, "ZeemanPolarization")
      .def(py::init([]() { return std::make_shared<Zeeman::Polarization>(); }), "Default polarization")
      .def(py::init([ZeemanPolarizationEnumGetter](const std::string& c) {
             return ZeemanPolarizationEnumGetter(c);
           }),
           py::arg("str"), "From :class:`str`")
      .PythonInterfaceCopyValue(Zeeman::Polarization)
      .def("__repr__",
           [ZeemanPolarizationStringGetter](Zeeman::Polarization c) {
             return ZeemanPolarizationStringGetter(c);
           })
      .def("__str__",
           [ZeemanPolarizationStringGetter](Zeeman::Polarization c) {
             return ZeemanPolarizationStringGetter(c);
           })
      .def(py::pickle(
          [ZeemanPolarizationStringGetter](const Zeeman::Polarization& c) {
            return py::make_tuple(ZeemanPolarizationStringGetter(c));
          },
          [ZeemanPolarizationEnumGetter](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return ZeemanPolarizationEnumGetter(t[0].cast<std::string>());
          })).doc() = "Options for ZeemanPolarization";
  py::implicitly_convertible<std::string, Zeeman::Polarization>();

  artsclass<Zeeman::Model>(m, "ZeemanModel")
      .def(py::init([]() { return std::make_shared<Zeeman::Model>(); }), "Empty model")
      .def(py::init([](Numeric a, Numeric b) {
             return std::make_shared<Zeeman::Model>(a, b);
           }),
           py::arg("gu"),
           py::arg("gl"), "From two numeric values")
      .def(py::init([](std::array<Numeric, 2> a) {
             return std::make_shared<Zeeman::Model>(a[0], a[1]);
           }),
           py::arg("gs"), "From list of two values")
      .PythonInterfaceCopyValue(Zeeman::Model)
      .PythonInterfaceBasicRepresentation(Zeeman::Model)
      .def_property("gu", &Zeeman::Model::gu, &Zeeman::Model::gu, ":class:`~float` The upper state g")
      .def_property("gl", &Zeeman::Model::gl, &Zeeman::Model::gl, ":class:`~float` The lower state g")
      .def(py::pickle(
          [](const Zeeman::Model& t) { return py::make_tuple(t.gu(), t.gl()); },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 2, "Invalid state!")
            return std::make_shared<Zeeman::Model>(t[0].cast<Numeric>(),
                                     t[1].cast<Numeric>());
          })).doc() = "A Zeeman model";
  py::implicitly_convertible<std::array<Numeric, 2>, Zeeman::Model>();

  artsclass<LineShape::Output>(m, "LineShapeOutput")
      .def(py::init([](Numeric a,
                       Numeric b,
                       Numeric c,
                       Numeric d,
                       Numeric e,
                       Numeric f,
                       Numeric g,
                       Numeric h,
                       Numeric i) {
             return std::make_shared< LineShape::Output>(a, b, c, d, e, f, g, h, i);
           }),
           py::arg("g0") = 0,
           py::arg("d0") = 0,
           py::arg("g2") = 0,
           py::arg("d2") = 0,
           py::arg("fvc") = 0,
           py::arg("eta") = 0,
           py::arg("y") = 0,
           py::arg("g") = 0,
           py::arg("dv") = 0, "From values")
      .PythonInterfaceCopyValue(LineShape::Output)
      .PythonInterfaceReadWriteData(LineShape::Output, G0, ":class:`float`: Pressure broadening speed-independent")
      .PythonInterfaceReadWriteData(LineShape::Output, D0, ":class:`float`: Pressure f-shifting speed-independent")
      .PythonInterfaceReadWriteData(LineShape::Output, G2, ":class:`float`: Pressure broadening speed-dependent")
      .PythonInterfaceReadWriteData(LineShape::Output, D2, ":class:`float`: Pressure f-shifting speed-dependent")
      .PythonInterfaceReadWriteData(LineShape::Output, ETA, ":class:`float`: Correlation")
      .PythonInterfaceReadWriteData(LineShape::Output, FVC, ":class:`float`: Frequency of velocity-changing collisions")
      .PythonInterfaceReadWriteData(LineShape::Output, Y, ":class:`float`: First order line mixing coefficient")
      .PythonInterfaceReadWriteData(LineShape::Output, G, ":class:`float`: Second order line mixing coefficient")
      .PythonInterfaceReadWriteData(LineShape::Output, DV, ":class:`float`: Second order line mixing f-shifting")
      .PythonInterfaceBasicRepresentation(LineShape::Output)
      .def(py::pickle(
          [](const LineShape::Output& t) {
            return py::make_tuple(
                t.G0, t.D0, t.G2, t.D2, t.FVC, t.ETA, t.Y, t.G, t.DV);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 9, "Invalid state!")
            return std::make_shared< LineShape::Output>(t[0].cast<Numeric>(),
                                         t[1].cast<Numeric>(),
                                         t[2].cast<Numeric>(),
                                         t[3].cast<Numeric>(),
                                         t[4].cast<Numeric>(),
                                         t[5].cast<Numeric>(),
                                         t[6].cast<Numeric>(),
                                         t[7].cast<Numeric>(),
                                         t[8].cast<Numeric>());
          })).doc() = "Derived line shape parameters";

  artsclass<LineShapeSingleSpeciesModel>(m, "LineShapeSingleSpeciesModel")
      .def(py::init([](LineShape::ModelParameters G0,
                       LineShape::ModelParameters D0,
                       LineShape::ModelParameters G2,
                       LineShape::ModelParameters D2,
                       LineShape::ModelParameters FVC,
                       LineShape::ModelParameters ETA,
                       LineShape::ModelParameters Y,
                       LineShape::ModelParameters G,
                       LineShape::ModelParameters DV) {
             return std::make_shared<LineShapeSingleSpeciesModel>(G0, D0, G2, D2, FVC, ETA, Y, G, DV);
           }),
           py::arg("G0") = LineShape::ModelParameters{},
           py::arg("D0") = LineShape::ModelParameters{},
           py::arg("G2") = LineShape::ModelParameters{},
           py::arg("D2") = LineShape::ModelParameters{},
           py::arg("FVC") = LineShape::ModelParameters{},
           py::arg("ETA") = LineShape::ModelParameters{},
           py::arg("Y") = LineShape::ModelParameters{},
           py::arg("G") = LineShape::ModelParameters{},
           py::arg("DV") = LineShape::ModelParameters{}, "From values")
      .PythonInterfaceCopyValue(LineShapeSingleSpeciesModel)
      .PythonInterfaceBasicRepresentation(LineShapeSingleSpeciesModel)
      .def_property(
          "G0",
          [](const LineShapeSingleSpeciesModel& x) { return x.G0(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.G0() = y;
          }, ":class:`~pyarts.arts.LineShapeModelParameters`: Pressure broadening speed-independent")
      .def_property(
          "D0",
          [](const LineShapeSingleSpeciesModel& x) { return x.D0(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.D0() = y;
          }, ":class:`~pyarts.arts.LineShapeModelParameters`: Pressure f-shifting speed-independent")
      .def_property(
          "G2",
          [](const LineShapeSingleSpeciesModel& x) { return x.G2(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.G2() = y;
          }, ":class:`~pyarts.arts.LineShapeModelParameters`: Pressure broadening speed-dependent")
      .def_property(
          "D2",
          [](const LineShapeSingleSpeciesModel& x) { return x.D2(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.D2() = y;
          }, ":class:`~pyarts.arts.LineShapeModelParameters`: Pressure f-shifting speed-dependent")
      .def_property(
          "FVC",
          [](const LineShapeSingleSpeciesModel& x) { return x.FVC(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.FVC() = y;
          }, ":class:`~pyarts.arts.LineShapeModelParameters`: Frequency of velocity-changing collisions")
      .def_property(
          "ETA",
          [](const LineShapeSingleSpeciesModel& x) { return x.ETA(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.ETA() = y;
          }, ":class:`~pyarts.arts.LineShapeModelParameters`: Correlation")
      .def_property(
          "Y",
          [](const LineShapeSingleSpeciesModel& x) { return x.Y(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.Y() = y;
          }, ":class:`~pyarts.arts.LineShapeModelParameters`: First order line mixing coefficient")
      .def_property(
          "G",
          [](const LineShapeSingleSpeciesModel& x) { return x.G(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.G() = y;
          }, ":class:`~pyarts.arts.LineShapeModelParameters`: Second order line mixing coefficient")
      .def_property(
          "DV",
          [](const LineShapeSingleSpeciesModel& x) { return x.DV(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.DV() = y;
          }, ":class:`~pyarts.arts.LineShapeModelParameters`: Second order line mixing f-shifting")
      .def(py::pickle(
          [](const LineShapeSingleSpeciesModel& t) {
            return py::make_tuple(t.G0(),
                                  t.D0(),
                                  t.G2(),
                                  t.D2(),
                                  t.FVC(),
                                  t.ETA(),
                                  t.Y(),
                                  t.G(),
                                  t.DV());
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 9, "Invalid state!")
            return std::make_shared<LineShapeSingleSpeciesModel>(
                t[0].cast<LineShape::ModelParameters>(),
                t[1].cast<LineShape::ModelParameters>(),
                t[2].cast<LineShape::ModelParameters>(),
                t[3].cast<LineShape::ModelParameters>(),
                t[4].cast<LineShape::ModelParameters>(),
                t[5].cast<LineShape::ModelParameters>(),
                t[6].cast<LineShape::ModelParameters>(),
                t[7].cast<LineShape::ModelParameters>(),
                t[8].cast<LineShape::ModelParameters>());
          })).doc() = "Single species line shape model";

  artsclass<LineShapeModel>(m, "LineShapeModel")
      .def(py::init([]() { return std::make_shared<LineShapeModel>(); }), "Empty model")
      .def(py::init([](const std::vector<LineShapeSingleSpeciesModel>& v) {
        return std::make_shared<LineShapeModel>(v);
      }), "From :class:`list`")
      .PythonInterfaceCopyValue(LineShapeModel)
      .PythonInterfaceBasicRepresentation(LineShapeModel)
      .PythonInterfaceIndexItemAccess(LineShapeModel)
      .def_property(
          "data",
          [](const LineShapeModel& x) { return x.Data(); },
          [](LineShapeModel& x, std::vector<LineShapeSingleSpeciesModel> y) {
            x.Data() = std::move(y);
          }, ":class:`list` of single species models")
      .def(py::pickle(
          [](const LineShapeModel& t) { return py::make_tuple(t.Data()); },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return std::make_shared<LineShapeModel>(
                t[0].cast<std::vector<LineShapeSingleSpeciesModel>>());
          })).doc() = "Multi-species line shape model";
  py::implicitly_convertible<std::vector<LineShapeSingleSpeciesModel>,
                             LineShapeModel>();

  artsclass<AbsorptionSingleLine>(m, "AbsorptionSingleLine")
      .def(py::init([](Numeric a,
                       Numeric b,
                       Numeric c,
                       Numeric d,
                       Numeric e,
                       Numeric f,
                       ZeemanModel g,
                       LineShapeModel h,
                       Quantum::Number::LocalState i) {
             return std::make_shared<AbsorptionSingleLine>(
                 a, b, c, d, e, f, g, std::move(h), std::move(i));
           }),
           py::arg("F0") = 0,
           py::arg("I0") = 0,
           py::arg("E0") = 0,
           py::arg("glow") = 0,
           py::arg("gupp") = 0,
           py::arg("A") = 0,
           py::arg("zeeman") = Zeeman::Model(),
           py::arg("lineshape") = LineShape::Model(),
           py::arg("localquanta") = Quantum::Number::LocalState{}, "From values")
      .PythonInterfaceCopyValue(AbsorptionSingleLine)
      .PythonInterfaceBasicRepresentation(AbsorptionSingleLine)
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, F0, ":class:`float` Central frequency")
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, I0, ":class:`float` Reference intensity")
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, E0, ":class:`float` Lower energy state")
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, glow, ":class:`float` Lower level statistical weight")
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, gupp, ":class:`float` Upper level statistical weight")
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, A, ":class:`float` Einstein spontaneous emission coefficient")
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, zeeman, ":class:`pyarts.arts.ZeemanModel` Zeeman model")
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, lineshape, ":class:`pyarts.arts.LineShapeModel` Line shape model")
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, localquanta, ":class:`pyarts.arts.QuantumNumberLocalState` Local quantum numbers")
      .def(py::pickle(
          [](const AbsorptionSingleLine& t) {
            return py::make_tuple(t.F0,
                                  t.I0,
                                  t.E0,
                                  t.glow,
                                  t.gupp,
                                  t.A,
                                  t.zeeman,
                                  t.lineshape,
                                  t.localquanta);
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 9, "Invalid state!")
            return std::make_shared<AbsorptionSingleLine>(
                t[0].cast<Numeric>(),
                t[1].cast<Numeric>(),
                t[2].cast<Numeric>(),
                t[3].cast<Numeric>(),
                t[4].cast<Numeric>(),
                t[5].cast<Numeric>(),
                t[6].cast<Zeeman::Model>(),
                t[7].cast<LineShape::Model>(),
                t[8].cast<Quantum::Number::LocalState>());
          })).doc() = "Single absorption line";

  artsarray<Array<AbsorptionSingleLine>>(m, "ArrayOfAbsorptionSingleLine")
      .doc() = "List of :class:`~pyarts.arts.AbsorptionSingleLine`";

  artsclass<AbsorptionLines>(m, "AbsorptionLines")
      .def(
          py::init([](bool selfbroadening, bool bathbroadening,
                      AbsorptionCutoffTypeOld cutoff,
                      AbsorptionMirroringTypeOld mirroring,
                      AbsorptionPopulationTypeOld population,
                      AbsorptionNormalizationTypeOld normalization,
                      LineShapeTypeOld lineshapetype, Numeric T0,
                      Numeric cutofffreq, Numeric linemixinglimit,
                      QuantumIdentifier quantumidentity,
                      ArrayOfSpeciesEnum broadeningspecies,
                      Array<AbsorptionSingleLine> lines) {
            ARTS_USER_ERROR_IF(not good_enum(quantumidentity.Species()),
                               "Bad quantumidentity, must specify species and "
                               "isotoplogue, e.g., H2O-161")
            ARTS_USER_ERROR_IF(
                quantumidentity.Isotopologue().joker() or
                    Species::is_predefined_model(
                        quantumidentity.Isotopologue()),
                "Bad quantumidentity, must specify isotopologue, e.g., H2O-161")
            ARTS_USER_ERROR_IF(
                std::any_of(lines.begin(), lines.end(),
                            [count = broadeningspecies.size()](auto &line) {
                              return line.lineshape.size() not_eq count;
                            }),
                "Incorrect size of broadeningspecies vs the line shape model")
            ARTS_USER_ERROR_IF(
                broadeningspecies.size() < static_cast<std::size_t>(selfbroadening + bathbroadening),
                "Must have atleast ", (selfbroadening + bathbroadening),
                " broadening species to support settings")
            return std::make_shared<AbsorptionLines>(selfbroadening,
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
          }),
          py::arg("selfbroadening") = false, py::arg("bathbroadening") = false,
          py::arg("cutoff") = AbsorptionCutoffTypeOld::None,
          py::arg("mirroring") = AbsorptionMirroringTypeOld::None,
          py::arg("population") = AbsorptionPopulationTypeOld::LTE,
          py::arg("normalization") = AbsorptionNormalizationTypeOld::None,
          py::arg("lineshapetype") = LineShapeTypeOld::DP, py::arg("T0") = 296,
          py::arg("cutofffreq") = -1, py::arg("linemixinglimit") = -1,
          py::arg("quantumidentity") = QuantumIdentifier("H2O-161"),
          py::arg("broadeningspecies") = ArrayOfSpeciesEnum{},
          py::arg("lines") = Array<AbsorptionSingleLine>{}, "From values")
      .PythonInterfaceCopyValue(AbsorptionLines)
      .PythonInterfaceWorkspaceVariableConversion(AbsorptionLines)
      .def(
          "__str__", [](const AbsorptionLines &x) { return var_string(x); })
      .def(
          "__repr__",
          [](const AbsorptionLines &x) {
            return var_string("'",
                              x.quantumidentity,
                              "'-band of ",
                              x.lines.size(),
                              " lines");
          })
      .PythonInterfaceFileIO(AbsorptionLines)
      .PythonInterfaceReadWriteData(AbsorptionLines, selfbroadening, ":class:`bool` Does the line broadening have self broadening?")
      .PythonInterfaceReadWriteData(AbsorptionLines, bathbroadening, ":class:`bool` Does the line broadening have bath broadening?")
      .PythonInterfaceReadWriteData(AbsorptionLines, cutoff, ":class:`~pyarts.arts.options.AbsorptionCutoffTypeOld` Cutoff type")
      .PythonInterfaceReadWriteData(AbsorptionLines, mirroring, ":class:`~pyarts.arts.options.AbsorptionMirroringype` Mirroring type")
      .PythonInterfaceReadWriteData(AbsorptionLines, population, ":class:`~pyarts.arts.options.AbsorptionPopulationTypeOld` Line population distribution")
      .PythonInterfaceReadWriteData(AbsorptionLines, normalization, ":class:`~pyarts.arts.options.AbsorptionNormalizationTypeOld` Normalization type")
      .PythonInterfaceReadWriteData(AbsorptionLines, lineshapetype, ":class:`~pyarts.arts.options.LineShapeTypeOld` Line shape type")
      .PythonInterfaceReadWriteData(AbsorptionLines, T0, ":class:`float` Reference temperature for all parameters of the lines")
      .PythonInterfaceReadWriteData(AbsorptionLines, cutofffreq, ":class:`float` Cutoff frequency")
      .PythonInterfaceReadWriteData(AbsorptionLines, linemixinglimit, ":class:`float` Linemixing limit")
      .PythonInterfaceReadWriteData(AbsorptionLines, quantumidentity, ":class:`~pyarts.arts.QuantumIdentifier` Catalog ID")
      .PythonInterfaceReadWriteData(AbsorptionLines, broadeningspecies, ":class:`~pyarts.arts.ArrayOfSpecies` A list of broadening specie")
      .PythonInterfaceReadWriteData(AbsorptionLines, lines, ":class:`~pyarts.arts.AbsorptionSingleLine` A list of individual lines")
      .def_property_readonly(
          "ok", &AbsorptionLines::OK,
          py::doc(R"(:class:`bool` If False, the catalog cannot be used for any calculations)"))
      .def_property_readonly("meta_data", &AbsorptionLines::MetaData,
                             py::doc(R"(:class:`~pyarts.arts.String` Catalog meta data string)"))
      .def(
          "LineShapeOutput",
          [](AbsorptionLines &band, Index line, Numeric T, Numeric P,
             const Vector &VMR) {
            ARTS_USER_ERROR_IF(not band.OK(), "Band in bad shape")
            ARTS_USER_ERROR_IF(not(T > 0) or not(P >= 0),
                               "Bad atmospheric state (T P): ",
                               Vector{T, P})
            ARTS_USER_ERROR_IF(
                static_cast<Size>(VMR.size()) not_eq band.broadeningspecies.size(),
                "Mismatch between VMRs and broadening species.\nVMR: ",
                VMR,
                "\nSpecies: ",
                band.broadeningspecies)
            band.lines.at(line);

            return band.ShapeParameters(line, T, P, VMR);
          }, py::arg("line"), py::arg("T"), py::arg("P"), py::arg("VMR"),
          py::doc(
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
)--"))
      .def(py::pickle(
          [](const AbsorptionLines &t) {
            return py::make_tuple(t.selfbroadening,
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
          },
          [](const py::tuple &t) {
            ARTS_USER_ERROR_IF(t.size() != 13, "Invalid state!")
            return std::make_shared<AbsorptionLines>(
                t[0].cast<bool>(),
                t[1].cast<bool>(),
                t[2].cast<AbsorptionCutoffTypeOld>(),
                t[3].cast<AbsorptionMirroringTypeOld>(),
                t[4].cast<AbsorptionPopulationTypeOld>(),
                t[5].cast<AbsorptionNormalizationTypeOld>(),
                t[6].cast<LineShapeTypeOld>(),
                t[7].cast<Numeric>(),
                t[8].cast<Numeric>(),
                t[9].cast<Numeric>(),
                t[10].cast<QuantumIdentifier>(),
                t[11].cast<ArrayOfSpeciesEnum>(),
                t[12].cast<Array<AbsorptionSingleLine>>());
          }))
      .PythonInterfaceWorkspaceDocumentation(AbsorptionLines);

  artsarray<ArrayOfAbsorptionLines>(m, "ArrayOfAbsorptionLines")
      .PythonInterfaceFileIO(ArrayOfAbsorptionLines)
.def("fuzzy_find_all", [](const ArrayOfAbsorptionLines& a, const QuantumIdentifier& q) {
      return fuzzy_find_all(a, q);
    }, py::arg("q"), "Find all the indexes that could match the given quantum identifier")
      .PythonInterfaceWorkspaceDocumentation(ArrayOfAbsorptionLines);

  artsarray<ArrayOfArrayOfAbsorptionLines>(m,
                                                "ArrayOfArrayOfAbsorptionLines")
      .PythonInterfaceFileIO(ArrayOfArrayOfAbsorptionLines)
      .PythonInterfaceWorkspaceDocumentation(ArrayOfAbsorptionLines);
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize spectroscopy\n", e.what()));
}
}  // namespace Python
