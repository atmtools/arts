#include <py_auto_interface.h>

#include "py_macros.h"

namespace Python {
void py_spectroscopy(py::module_& m) {
  py::class_<LineShape::TemperatureModel>(m, "LineShapeTemperatureModel")
      .def(py::init<>())
      .def(py::init(
          [](char* c) { return LineShape::toTemperatureModelOrThrow(c); }))
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
           py::arg("X0") = std::numeric_limits<Numeric>::quiet_NaN(),
           py::arg("X1") = std::numeric_limits<Numeric>::quiet_NaN(),
           py::arg("X2") = std::numeric_limits<Numeric>::quiet_NaN(),
           py::arg("X3") = std::numeric_limits<Numeric>::quiet_NaN())
      .def_readwrite("type", &LineShapeModelParameters::type)
      .def_readwrite("X0", &LineShapeModelParameters::X0)
      .def_readwrite("X1", &LineShapeModelParameters::X1)
      .def_readwrite("X2", &LineShapeModelParameters::X2)
      .def_readwrite("X3", &LineShapeModelParameters::X3);

  py::class_<AbsorptionCutoffType>(m, "AbsorptionCutoffType")
      .def(py::init<>())
      .def(py::init([](char* c) { return Absorption::toCutoffTypeOrThrow(c); }))
      .PythonInterfaceBasicRepresentation(AbsorptionCutoffType);
  py::implicitly_convertible<py::str, AbsorptionCutoffType>();

  py::class_<AbsorptionMirroringType>(m, "AbsorptionMirroringType")
      .def(py::init<>())
      .def(py::init(
          [](char* c) { return Absorption::toMirroringTypeOrThrow(c); }))
      .PythonInterfaceBasicRepresentation(AbsorptionMirroringType);
  py::implicitly_convertible<py::str, AbsorptionMirroringType>();

  py::class_<AbsorptionPopulationType>(m, "AbsorptionPopulationType")
      .def(py::init<>())
      .def(py::init(
          [](char* c) { return Absorption::toPopulationTypeOrThrow(c); }))
      .PythonInterfaceBasicRepresentation(AbsorptionPopulationType);
  py::implicitly_convertible<py::str, AbsorptionPopulationType>();

  py::class_<AbsorptionNormalizationType>(m, "AbsorptionNormalizationType")
      .def(py::init<>())
      .def(py::init(
          [](char* c) { return Absorption::toNormalizationTypeOrThrow(c); }))
      .PythonInterfaceBasicRepresentation(AbsorptionNormalizationType);
  py::implicitly_convertible<py::str, AbsorptionNormalizationType>();

  py::class_<LineShapeType>(m, "LineShapeType")
      .def(py::init<>())
      .def(py::init([](char* c) { return LineShape::toTypeOrThrow(c); }))
      .PythonInterfaceBasicRepresentation(LineShapeType);
  py::implicitly_convertible<py::str, LineShapeType>();

  py::class_<Zeeman::Model>(m, "ZeemanModel")
      .def(py::init<>())
      .def(py::init<Numeric, Numeric>())
      .PythonInterfaceBasicRepresentation(Zeeman::Model)
      .def_property("gu", &Zeeman::Model::gu, &Zeeman::Model::gu)
      .def_property("gl", &Zeeman::Model::gl, &Zeeman::Model::gl);

  py::class_<LineShapeSingleSpeciesModel>(m, "LineShapeSingleSpeciesModel")
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
      .PythonInterfaceBasicRepresentation(LineShapeModel)
      .PythonInterfaceIndexItemAccess(LineShapeModel)
      .def_property(
          "data",
          [](const LineShapeModel& x) { return x.Data(); },
          [](LineShapeModel& x, std::vector<LineShapeSingleSpeciesModel> y) {
            x.Data() = std::move(y);
          });

  py::class_<AbsorptionSingleLine>(m, "AbsorptionSingleLine")
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

  py::class_<AbsorptionLines, std::shared_ptr<AbsorptionLines>>(
      m, "AbsorptionLines")
      .def(py::init<>())
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
      .PythonInterfaceReadWriteData(AbsorptionLines, lines);

  PythonInterfaceWorkspaceArray(AbsorptionLines);
  PythonInterfaceWorkspaceArray(ArrayOfAbsorptionLines);

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