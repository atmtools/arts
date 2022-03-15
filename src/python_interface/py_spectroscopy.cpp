#include <py_auto_interface.h>

#include "py_macros.h"

#include <lineshape.h>
#include <lineshapemodel.h>
#include <zeemandata.h>

namespace Python {
void py_spectroscopy(py::module_& m) {
  py::class_<LineShape::TemperatureModel>(m, "LineShapeTemperatureModel")
      .def(py::init<>())
      .def(py::init([](const char* c) {
             return LineShape::toTemperatureModelOrThrow(c);
           }),
           py::arg("str").none(false))
      .PythonInterfaceBasicRepresentation(LineShape::TemperatureModel);
  py::implicitly_convertible<py::str, LineShape::TemperatureModel>();

  static_assert(LineShapeModelParameters::N == 4);
  py::class_<LineShapeModelParameters>(m, "LineShapeModelParameters")
      .def(py::init<LineShape::TemperatureModel,
                    Numeric_,
                    Numeric_,
                    Numeric_,
                    Numeric_>(),
           py::arg("type") = LineShape::TemperatureModel::None,
           py::arg("X0") = 0,
           py::arg("X1") = 0,
           py::arg("X2") = 0,
           py::arg("X3") = 0)
      .def_readwrite("type", &LineShapeModelParameters::type)
      .def_readwrite("X0", &LineShapeModelParameters::X0)
      .def_readwrite("X1", &LineShapeModelParameters::X1)
      .def_readwrite("X2", &LineShapeModelParameters::X2)
      .def_readwrite("X3", &LineShapeModelParameters::X3)
      .PythonInterfaceBasicRepresentation(LineShapeModelParameters);

  py::class_<AbsorptionCutoffType>(m, "AbsorptionCutoffType")
      .def(py::init<>())
      .def(py::init([](const char* c) {
             return Absorption::toCutoffTypeOrThrow(c);
           }),
           py::arg("str").none(false))
      .PythonInterfaceBasicRepresentation(AbsorptionCutoffType);
  py::implicitly_convertible<py::str, AbsorptionCutoffType>();

  py::class_<Zeeman::Polarization>(m, "ZeemanPolarization")
      .def(py::init<>())
      .def(py::init([](String c) {
             if (c == "SigmaMinus") return Zeeman::Polarization::SigmaMinus;
             if (c == "Pi") return Zeeman::Polarization::Pi;
             if (c == "SigmaPlus") return Zeeman::Polarization::SigmaPlus;
             if (c == "None") return Zeeman::Polarization::None;
             ARTS_USER_ERROR("Bad enum value ", c);
           }),
           py::arg("str"))
      .def("__repr__",
           [](Zeeman::Polarization c) {
             if (c == Zeeman::Polarization::SigmaMinus) return "SigmaMinus";
             if (c == Zeeman::Polarization::Pi) return "Pi";
             if (c == Zeeman::Polarization::SigmaPlus) return "SigmaPlus";
             if (c == Zeeman::Polarization::None) return "None";
             ARTS_USER_ERROR("Bad enum state")
           })
      .def("__str__", [](Zeeman::Polarization c) {
        if (c == Zeeman::Polarization::SigmaMinus) return "SigmaMinus";
        if (c == Zeeman::Polarization::Pi) return "Pi";
        if (c == Zeeman::Polarization::SigmaPlus) return "SigmaPlus";
        if (c == Zeeman::Polarization::None) return "None";
        ARTS_USER_ERROR("Bad enum state")
      });
  py::implicitly_convertible<py::str, Zeeman::Polarization>();

  py::class_<AbsorptionMirroringType>(m, "AbsorptionMirroringType")
      .def(py::init<>())
      .def(py::init([](const char* c) {
             return Absorption::toMirroringTypeOrThrow(c);
           }),
           py::arg("str").none(false))
      .PythonInterfaceBasicRepresentation(AbsorptionMirroringType);
  py::implicitly_convertible<py::str, AbsorptionMirroringType>();

  py::class_<AbsorptionPopulationType>(m, "AbsorptionPopulationType")
      .def(py::init<>())
      .def(py::init([](const char* c) {
             return Absorption::toPopulationTypeOrThrow(c);
           }),
           py::arg("str").none(false))
      .PythonInterfaceBasicRepresentation(AbsorptionPopulationType);
  py::implicitly_convertible<py::str, AbsorptionPopulationType>();

  py::class_<AbsorptionNormalizationType>(m, "AbsorptionNormalizationType")
      .def(py::init<>())
      .def(py::init([](const char* c) {
             return Absorption::toNormalizationTypeOrThrow(c);
           }),
           py::arg("str").none(false))
      .PythonInterfaceBasicRepresentation(AbsorptionNormalizationType);
  py::implicitly_convertible<py::str, AbsorptionNormalizationType>();

  py::class_<LineShapeType>(m, "LineShapeType")
      .def(py::init<>())
      .def(py::init([](const char* c) { return LineShape::toTypeOrThrow(c); }),
           py::arg("str").none(false))
      .PythonInterfaceBasicRepresentation(LineShapeType);
  py::implicitly_convertible<py::str, LineShapeType>();

  py::class_<Zeeman::Model>(m, "ZeemanModel")
      .def(py::init<>())
      .def(py::init<Numeric, Numeric>())
      .PythonInterfaceBasicRepresentation(Zeeman::Model)
      .def_property("gu", &Zeeman::Model::gu, &Zeeman::Model::gu)
      .def_property("gl", &Zeeman::Model::gl, &Zeeman::Model::gl);

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
           py::arg("g0") = 0,
           py::arg("d0") = 0,
           py::arg("g2") = 0,
           py::arg("d2") = 0,
           py::arg("fvc") = 0,
           py::arg("eta") = 0,
           py::arg("y") = 0,
           py::arg("g") = 0,
           py::arg("dv") = 0)
      .PythonInterfaceReadWriteData(LineShape::Output, G0)
      .PythonInterfaceReadWriteData(LineShape::Output, D0)
      .PythonInterfaceReadWriteData(LineShape::Output, G2)
      .PythonInterfaceReadWriteData(LineShape::Output, D2)
      .PythonInterfaceReadWriteData(LineShape::Output, ETA)
      .PythonInterfaceReadWriteData(LineShape::Output, FVC)
      .PythonInterfaceReadWriteData(LineShape::Output, Y)
      .PythonInterfaceReadWriteData(LineShape::Output, G)
      .PythonInterfaceReadWriteData(LineShape::Output, DV)
      .PythonInterfaceBasicRepresentation(LineShape::Output);

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
           py::arg("G0") = LineShape::ModelParameters{},
           py::arg("D0") = LineShape::ModelParameters{},
           py::arg("G2") = LineShape::ModelParameters{},
           py::arg("D2") = LineShape::ModelParameters{},
           py::arg("FVC") = LineShape::ModelParameters{},
           py::arg("ETA") = LineShape::ModelParameters{},
           py::arg("Y") = LineShape::ModelParameters{},
           py::arg("G") = LineShape::ModelParameters{},
           py::arg("DV") = LineShape::ModelParameters{})
      .def_property(
          "G0",
          [](const LineShapeSingleSpeciesModel& x) { return x.G0(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.G0() = std::move(y);
          })
      .def_property(
          "D0",
          [](const LineShapeSingleSpeciesModel& x) { return x.D0(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.D0() = std::move(y);
          })
      .def_property(
          "G2",
          [](const LineShapeSingleSpeciesModel& x) { return x.G2(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.G2() = std::move(y);
          })
      .def_property(
          "D2",
          [](const LineShapeSingleSpeciesModel& x) { return x.D2(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.D2() = std::move(y);
          })
      .def_property(
          "FVC",
          [](const LineShapeSingleSpeciesModel& x) { return x.FVC(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.FVC() = std::move(y);
          })
      .def_property(
          "ETA",
          [](const LineShapeSingleSpeciesModel& x) { return x.ETA(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.ETA() = std::move(y);
          })
      .def_property(
          "Y",
          [](const LineShapeSingleSpeciesModel& x) { return x.Y(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.Y() = std::move(y);
          })
      .def_property(
          "G",
          [](const LineShapeSingleSpeciesModel& x) { return x.G(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.G() = std::move(y);
          })
      .def_property(
          "DV",
          [](const LineShapeSingleSpeciesModel& x) { return x.DV(); },
          [](LineShapeSingleSpeciesModel& x, LineShapeModelParameters y) {
            return x.DV() = std::move(y);
          });

  py::class_<LineShapeModel>(m, "LineShapeModel")
      .def(py::init<>())
      .def(py::init<std::vector<LineShapeSingleSpeciesModel>>())
      .PythonInterfaceBasicRepresentation(LineShapeModel)
      .PythonInterfaceIndexItemAccess(LineShapeModel)
      .def_property(
          "data",
          [](const LineShapeModel& x) { return x.Data(); },
          [](LineShapeModel& x, std::vector<LineShapeSingleSpeciesModel> y) {
            x.Data() = std::move(y);
          });
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
           py::arg("F0") = 0,
           py::arg("I0") = 0,
           py::arg("E0") = 0,
           py::arg("glow") = 0,
           py::arg("gupp") = 0,
           py::arg("A") = 0,
           py::arg("zeeman") = Zeeman::Model(),
           py::arg("lineshape") = LineShape::Model(),
           py::arg("localquanta") = Quantum::Number::LocalState{})
      .PythonInterfaceBasicRepresentation(AbsorptionSingleLine)
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, F0)
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, I0)
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, E0)
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, glow)
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, gupp)
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, A)
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, zeeman)
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, lineshape)
      .PythonInterfaceReadWriteData(AbsorptionSingleLine, localquanta);

  py::class_<Array<AbsorptionSingleLine>>(m, "ArrayOfAbsorptionSingleLine")
      .PythonInterfaceBasicRepresentation(Array<AbsorptionSingleLine>)
      .PythonInterfaceArrayDefault(AbsorptionSingleLine)
      .doc() = "The Arts ArrayOfAbsorptionSingleLine class";
  py::implicitly_convertible<std::vector<AbsorptionSingleLine>,
                             Array<AbsorptionSingleLine>>();

  py::class_<AbsorptionLines>(m, "AbsorptionLines")
      .def(py::init<>())
      .def(py::init<bool,
                    bool,
                    AbsorptionCutoffType,
                    AbsorptionMirroringType,
                    AbsorptionPopulationType,
                    AbsorptionNormalizationType,
                    LineShape::Type,
                    Numeric,
                    Numeric,
                    Numeric,
                    QuantumIdentifier,
                    ArrayOfSpecies,
                    Array<AbsorptionSingleLine>>(),
           py::arg("selfbroadening") = false,
           py::arg("bathbroadening") = false,
           py::arg("cutoff") = AbsorptionCutoffType::None,
           py::arg("mirroring") = AbsorptionMirroringType::None,
           py::arg("population") = AbsorptionPopulationType::LTE,
           py::arg("normalization") = AbsorptionNormalizationType::None,
           py::arg("lineshapetype") = LineShape::Type::DP,
           py::arg("T0") = 296,
           py::arg("cutofffreq") = -1,
           py::arg("linemixinglimit") = -1,
           py::arg("quantumidentity") = QuantumIdentifier(),
           py::arg("broadeningspecies") = ArrayOfSpecies{},
           py::arg("lines") = Array<AbsorptionSingleLine>{})
      .PythonInterfaceWorkspaceVariableConversion(AbsorptionLines)
      .PythonInterfaceBasicRepresentation(AbsorptionLines)
      .PythonInterfaceFileIO(AbsorptionLines)
      .PythonInterfaceReadWriteData(AbsorptionLines, selfbroadening)
      .PythonInterfaceReadWriteData(AbsorptionLines, bathbroadening)
      .PythonInterfaceReadWriteData(AbsorptionLines, cutoff)
      .PythonInterfaceReadWriteData(AbsorptionLines, mirroring)
      .PythonInterfaceReadWriteData(AbsorptionLines, population)
      .PythonInterfaceReadWriteData(AbsorptionLines, normalization)
      .PythonInterfaceReadWriteData(AbsorptionLines, lineshapetype)
      .PythonInterfaceReadWriteData(AbsorptionLines, T0)
      .PythonInterfaceReadWriteData(AbsorptionLines, cutofffreq)
      .PythonInterfaceReadWriteData(AbsorptionLines, linemixinglimit)
      .PythonInterfaceReadWriteData(AbsorptionLines, quantumidentity)
      .PythonInterfaceReadWriteData(AbsorptionLines, broadeningspecies)
      .PythonInterfaceReadWriteData(AbsorptionLines, lines)
      .def_property_readonly("ok", &AbsorptionLines::OK, py::doc("If false, the catalog cannot be used for any calculations"))
      .def_property_readonly("meta_data", &AbsorptionLines::MetaData, py::doc("Catalog meta data string"))
      .def("LineShapeOutput",
           [](AbsorptionLines& band,
              Index line,
              Numeric T,
              Numeric P,
              const Vector& VMR) {
             ARTS_USER_ERROR_IF(not band.OK(), "Band in bad shape")
             ARTS_USER_ERROR_IF(not (T > 0) or not (P >= 0),
                 "Bad atmospheric state (T P): ",
                 Vector{T, P})
             ARTS_USER_ERROR_IF(
                 VMR.size() not_eq band.broadeningspecies.nelem(),
                 "Mismatch between VMRs and broadening species.\nVMR: ",
                 VMR,
                 "\nSpecies: ",
                 band.broadeningspecies)
             band.lines.at(line);

             return band.ShapeParameters(line, T, P, VMR);
           }, py::doc("Computes the line shape paramters for the given atmospheric state\n\nNote that the normalization assumes sum(VMR) is 1 for good results but does not enforce it"));

  PythonInterfaceWorkspaceArray(AbsorptionLines);
  PythonInterfaceWorkspaceArray(ArrayOfAbsorptionLines);

  py::class_<LineShape::Calculator>(m, "LineShapeCalculator")
      .def(py::init([](AbsorptionLines& band,
                       Index line,
                       Numeric T,
                       Numeric P,
                       const Vector& VMR,
                       Zeeman::Polarization zeeman,
                       Numeric H,
                       Index iz) {
             ARTS_USER_ERROR_IF(not band.OK(), "Band in bad shape")
             ARTS_USER_ERROR_IF(not (T > 0) or not (P >= 0) or not (H >= 0),
                                "Bad atmospheric state (T P H): ",
                                Vector{T, P, H})
             ARTS_USER_ERROR_IF(
                 VMR.size() not_eq band.broadeningspecies.nelem(),
                 "Mismatch between VMRs and broadening species.\nVMR: ",
                 VMR,
                 "\nSpecies: ",
                 band.broadeningspecies)
             auto F0 = band.lines.at(line).F0;
             auto DC = band.DopplerConstant(T);
             auto mirror = band.mirroring;
             auto type = band.lineshapetype;
             auto X = band.ShapeParameters(line, T, P, VMR);
             auto DZ = band.ZeemanSplitting(line, zeeman, iz) * H;

             if (mirror == AbsorptionMirroringType::Manual)
               return LineShape::Calculator(mirror, type, F0, X, DC, DZ);
             return LineShape::Calculator(type, F0, X, DC, DZ, false);
           }),
           py::arg("band"),
           py::arg("line"),
           py::arg("T"),
           py::arg("P"),
           py::arg("VMR"),
           py::arg("zeeman") = Zeeman::Polarization::None,
           py::arg("H") = 0,
           py::arg("iz") = 0)
      .def("dFdT",
           [](LineShape::Calculator& LS,
              const LineShape::Output& dXdT,
              Numeric T) { return LS.dFdT(dXdT, T); })
      .def("dFdf", [](LineShape::Calculator& LS) { return LS.dFdf(); })
      .def("dFdF0", [](LineShape::Calculator& LS) { return LS.dFdF0(); })
      .def(
          "dFdH",
          [](LineShape::Calculator& LS, Numeric dfdH) { return LS.dFdH(dfdH); })
      .def("dFdVMR",
           [](LineShape::Calculator& LS, const LineShape::Output& dXdVMR) {
             return LS.dFdVMR(dXdVMR);
           })
      .def("dFdFVC",
           [](LineShape::Calculator& LS, Numeric d) { return LS.dFdFVC(d); })
      .def("dFdETA",
           [](LineShape::Calculator& LS, Numeric d) { return LS.dFdETA(d); })
      .def("dFdDV",
           [](LineShape::Calculator& LS, Numeric d) { return LS.dFdDV(d); })
      .def("dFdD0",
           [](LineShape::Calculator& LS, Numeric d) { return LS.dFdD0(d); })
      .def("dFdG0",
           [](LineShape::Calculator& LS, Numeric d) { return LS.dFdG0(d); })
      .def("dFdD2",
           [](LineShape::Calculator& LS, Numeric d) { return LS.dFdD2(d); })
      .def("dFdG2",
           [](LineShape::Calculator& LS, Numeric d) { return LS.dFdG2(d); })
      .def("F", [](LineShape::Calculator& LS) { return LS.F(); })
      .def("F", [](LineShape::Calculator& LS, Numeric f) { return LS(f); })
    .doc()="Class to compute the line shape\n\nNote that the normalization assumes sum(VMR) is 1 for good results but does not enforce it";

  py::class_<SpeciesErrorCorrectedSuddenData>(m,
                                              "SpeciesErrorCorrectedSuddenData")
      .def(py::init<>())
      .PythonInterfaceBasicRepresentation(SpeciesErrorCorrectedSuddenData)
      .def_readwrite("spec", &SpeciesErrorCorrectedSuddenData::spec)
      .def_readwrite("scaling", &SpeciesErrorCorrectedSuddenData::scaling)
      .def_readwrite("beta", &SpeciesErrorCorrectedSuddenData::beta)
      .def_readwrite("lambda", &SpeciesErrorCorrectedSuddenData::lambda)
      .def_readwrite("collisional_distance",
                     &SpeciesErrorCorrectedSuddenData::collisional_distance)
      .def_readwrite("mass", &SpeciesErrorCorrectedSuddenData::mass);

  py::class_<ErrorCorrectedSuddenData>(m, "ErrorCorrectedSuddenData")
      .def(py::init<>())
      .PythonInterfaceBasicRepresentation(ErrorCorrectedSuddenData)
      .def(
          "__getitem__",
          [](ErrorCorrectedSuddenData& x, Species::Species& y)
              -> SpeciesErrorCorrectedSuddenData& { return x[y]; },
          py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](ErrorCorrectedSuddenData& x,
             Species::Species& y,
             SpeciesErrorCorrectedSuddenData& z) { x[y] = z; },
          py::return_value_policy::reference_internal);

  py::class_<MapOfErrorCorrectedSuddenData>(m, "MapOfErrorCorrectedSuddenData")
      .def(py::init<>())
      .PythonInterfaceWorkspaceVariableConversion(MapOfErrorCorrectedSuddenData)
      .PythonInterfaceBasicRepresentation(MapOfErrorCorrectedSuddenData)
      .PythonInterfaceFileIO(MapOfErrorCorrectedSuddenData)
      .def(
          "__getitem__",
          [](MapOfErrorCorrectedSuddenData& x, QuantumIdentifier& y)
              -> ErrorCorrectedSuddenData& { return x[y]; },
          py::return_value_policy::reference_internal)
      .def(
          "__setitem__",
          [](MapOfErrorCorrectedSuddenData& x,
             QuantumIdentifier& y,
             ErrorCorrectedSuddenData& z) { x[y] = z; },
          py::return_value_policy::reference_internal);

  py::class_<EnergyLevelMap>(m, "EnergyLevelMap")
      .def(py::init<>())
      .PythonInterfaceWorkspaceVariableConversion(EnergyLevelMap)
      .PythonInterfaceBasicRepresentation(EnergyLevelMap)
      .PythonInterfaceFileIO(EnergyLevelMap);

  py::class_<CIARecord>(m, "CIARecord")
      .def(py::init<>())
      .PythonInterfaceWorkspaceVariableConversion(CIARecord)
      .PythonInterfaceBasicRepresentation(CIARecord)
      .PythonInterfaceFileIO(CIARecord);

  PythonInterfaceWorkspaceArray(CIARecord);

  py::class_<HitranRelaxationMatrixData>(m, "HitranRelaxationMatrixData")
      .def(py::init<>())
      .PythonInterfaceWorkspaceVariableConversion(HitranRelaxationMatrixData)
      .PythonInterfaceBasicRepresentation(HitranRelaxationMatrixData)
      .PythonInterfaceFileIO(HitranRelaxationMatrixData)
      .def_readwrite("W0pp", &HitranRelaxationMatrixData::W0pp)
      .def_readwrite("B0pp", &HitranRelaxationMatrixData::B0pp)
      .def_readwrite("W0rp", &HitranRelaxationMatrixData::W0rp)
      .def_readwrite("B0rp", &HitranRelaxationMatrixData::B0rp)
      .def_readwrite("W0qp", &HitranRelaxationMatrixData::W0qp)
      .def_readwrite("B0qp", &HitranRelaxationMatrixData::B0qp)
      .def_readwrite("W0pr", &HitranRelaxationMatrixData::W0pr)
      .def_readwrite("B0pr", &HitranRelaxationMatrixData::B0pr)
      .def_readwrite("W0rr", &HitranRelaxationMatrixData::W0rr)
      .def_readwrite("B0rr", &HitranRelaxationMatrixData::B0rr)
      .def_readwrite("W0qr", &HitranRelaxationMatrixData::W0qr)
      .def_readwrite("B0qr", &HitranRelaxationMatrixData::B0qr)
      .def_readwrite("W0pq", &HitranRelaxationMatrixData::W0pq)
      .def_readwrite("B0pq", &HitranRelaxationMatrixData::B0pq)
      .def_readwrite("W0rq", &HitranRelaxationMatrixData::W0rq)
      .def_readwrite("B0rq", &HitranRelaxationMatrixData::B0rq)
      .def_readwrite("W0qq", &HitranRelaxationMatrixData::W0qq)
      .def_readwrite("B0qq", &HitranRelaxationMatrixData::B0qq);

  py::class_<XsecRecord>(m, "XsecRecord")
      .def(py::init<>())
      .PythonInterfaceBasicRepresentation(XsecRecord);

  PythonInterfaceWorkspaceArray(XsecRecord);
}
}  // namespace Python