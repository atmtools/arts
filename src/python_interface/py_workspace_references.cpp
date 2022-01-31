#include <auto_md.h>
#include <xml_io.h>

#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_workspace_references(py::module_& m) {
  m.def(                                                                  
      "Index",                                                              
      [](std::size_t x) -> Index_& { return *reinterpret_cast<Index_*>(x); }, 
      py::return_value_policy::reference);
  m.def(                                                                  
      "Numeric",                                                              
      [](std::size_t x) -> Numeric_& { return *reinterpret_cast<Numeric_*>(x); }, 
      py::return_value_policy::reference);

PythonInterfaceReferenceFromPointer(AbsorptionLines)
//Agenda
//Any
PythonInterfaceReferenceFromPointer(ArrayOfAbsorptionLines)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfAbsorptionLines)
//ArrayOfAgenda
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfGriddedField1)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfGriddedField2)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfGriddedField3)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfIndex)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfMatrix)
PythonInterfaceReferenceFromPointer(ArrayOfPpath)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfPropagationMatrix)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfRadiationVector)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfScatteringMetaData)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfSingleScatteringData)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfSpeciesTag)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfStokesVector)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfString)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfTensor3)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfTensor6)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfTime)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfTransmissionMatrix)
PythonInterfaceReferenceFromPointer(ArrayOfArrayOfVector)
PythonInterfaceReferenceFromPointer(ArrayOfCIARecord)
PythonInterfaceReferenceFromPointer(ArrayOfGriddedField1)
PythonInterfaceReferenceFromPointer(ArrayOfGriddedField2)
PythonInterfaceReferenceFromPointer(ArrayOfGriddedField3)
PythonInterfaceReferenceFromPointer(ArrayOfGriddedField4)
PythonInterfaceReferenceFromPointer(ArrayOfIndex)
PythonInterfaceReferenceFromPointer(ArrayOfJacobianTarget)
PythonInterfaceReferenceFromPointer(ArrayOfMatrix)
PythonInterfaceReferenceFromPointer(ArrayOfPropagationMatrix)
PythonInterfaceReferenceFromPointer(ArrayOfQuantumIdentifier)
PythonInterfaceReferenceFromPointer(ArrayOfRadiationVector)
PythonInterfaceReferenceFromPointer(ArrayOfRetrievalQuantity)
PythonInterfaceReferenceFromPointer(ArrayOfScatteringMetaData)
PythonInterfaceReferenceFromPointer(ArrayOfSingleScatteringData)
PythonInterfaceReferenceFromPointer(ArrayOfSpeciesTag)
PythonInterfaceReferenceFromPointer(ArrayOfSparse)
PythonInterfaceReferenceFromPointer(ArrayOfStokesVector)
PythonInterfaceReferenceFromPointer(ArrayOfString)
PythonInterfaceReferenceFromPointer(ArrayOfTelsemAtlas)
PythonInterfaceReferenceFromPointer(ArrayOfTensor3)
PythonInterfaceReferenceFromPointer(ArrayOfTensor4)
PythonInterfaceReferenceFromPointer(ArrayOfTensor5)
PythonInterfaceReferenceFromPointer(ArrayOfTensor6)
PythonInterfaceReferenceFromPointer(ArrayOfTensor7)
PythonInterfaceReferenceFromPointer(ArrayOfTime)
PythonInterfaceReferenceFromPointer(ArrayOfTransmissionMatrix)
PythonInterfaceReferenceFromPointer(ArrayOfVector)
PythonInterfaceReferenceFromPointer(ArrayOfXsecRecord)
PythonInterfaceReferenceFromPointer(CIARecord)
PythonInterfaceReferenceFromPointer(CovarianceMatrix)
PythonInterfaceReferenceFromPointer(EnergyLevelMap)
PythonInterfaceReferenceFromPointer(GasAbsLookup)
PythonInterfaceReferenceFromPointer(GridPos)
PythonInterfaceReferenceFromPointer(GriddedField1)
PythonInterfaceReferenceFromPointer(GriddedField2)
PythonInterfaceReferenceFromPointer(GriddedField3)
PythonInterfaceReferenceFromPointer(GriddedField4)
PythonInterfaceReferenceFromPointer(GriddedField5)
PythonInterfaceReferenceFromPointer(GriddedField6)
PythonInterfaceReferenceFromPointer(HitranRelaxationMatrixData)
PythonInterfaceReferenceFromPointer(Index)
PythonInterfaceReferenceFromPointer(JacobianTarget)
PythonInterfaceReferenceFromPointer(MapOfErrorCorrectedSuddenData)
PythonInterfaceReferenceFromPointer(MCAntenna)
PythonInterfaceReferenceFromPointer(Matrix)
PythonInterfaceReferenceFromPointer(Numeric)
PythonInterfaceReferenceFromPointer(Ppath)
PythonInterfaceReferenceFromPointer(PropagationMatrix)
PythonInterfaceReferenceFromPointer(QuantumIdentifier)
PythonInterfaceReferenceFromPointer(RadiationVector)
PythonInterfaceReferenceFromPointer(Rational)
PythonInterfaceReferenceFromPointer(ScatteringMetaData)
PythonInterfaceReferenceFromPointer(SingleScatteringData)
PythonInterfaceReferenceFromPointer(Sparse)
PythonInterfaceReferenceFromPointer(SpeciesIsotopologueRatios)
PythonInterfaceReferenceFromPointer(StokesVector)
PythonInterfaceReferenceFromPointer(String)
PythonInterfaceReferenceFromPointer(TelsemAtlas)
PythonInterfaceReferenceFromPointer(Tensor3)
PythonInterfaceReferenceFromPointer(Tensor4)
PythonInterfaceReferenceFromPointer(Tensor5)
PythonInterfaceReferenceFromPointer(Tensor6)
PythonInterfaceReferenceFromPointer(Tensor7)
PythonInterfaceReferenceFromPointer(Timer)
PythonInterfaceReferenceFromPointer(Time)
PythonInterfaceReferenceFromPointer(TessemNN)
PythonInterfaceReferenceFromPointer(TransmissionMatrix)
PythonInterfaceReferenceFromPointer(Vector)
PythonInterfaceReferenceFromPointer(Verbosity)
}
}  // namespace Python