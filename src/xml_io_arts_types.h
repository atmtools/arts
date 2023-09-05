////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io_arts_types.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains private function declarations and
         template instantiation to handle XML data files.

*/

#ifndef xml_io_arts_types_h
#define xml_io_arts_types_h

#include <cfloat>
#include <stdexcept>

// All workspace groups are known by the agenda class through tokval.h
#include <workspace.h>

// Extras
#include "mc_interp.h"
#include "operators.h"
#include "template_partfun.h"

#define TMPL_XML_READ_WRITE_STREAM(what)                                       \
  void xml_read_from_stream(std::istream &, what &, bifstream *);              \
  void xml_write_to_stream(std::ostream &, const what &, bofstream *,          \
                           const String &);

////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for XML streams
////////////////////////////////////////////////////////////////////////////

//=== Basic Types ==========================================================

TMPL_XML_READ_WRITE_STREAM(PartitionFunctionsData)
TMPL_XML_READ_WRITE_STREAM(JacobianTarget)
TMPL_XML_READ_WRITE_STREAM(Rational)
TMPL_XML_READ_WRITE_STREAM(Time)
TMPL_XML_READ_WRITE_STREAM(VibrationalEnergyLevels)

//=== Compound Types =======================================================

TMPL_XML_READ_WRITE_STREAM(AbsorptionLines)
TMPL_XML_READ_WRITE_STREAM(Agenda)
TMPL_XML_READ_WRITE_STREAM(AtmField)
TMPL_XML_READ_WRITE_STREAM(AtmPoint)
TMPL_XML_READ_WRITE_STREAM(CIARecord)
TMPL_XML_READ_WRITE_STREAM(CovarianceMatrix)
TMPL_XML_READ_WRITE_STREAM(GasAbsLookup)
TMPL_XML_READ_WRITE_STREAM(GriddedField)
TMPL_XML_READ_WRITE_STREAM(GriddedField1)
TMPL_XML_READ_WRITE_STREAM(GriddedField2)
TMPL_XML_READ_WRITE_STREAM(GriddedField3)
TMPL_XML_READ_WRITE_STREAM(GriddedField4)
TMPL_XML_READ_WRITE_STREAM(GriddedField5)
TMPL_XML_READ_WRITE_STREAM(GriddedField6)
TMPL_XML_READ_WRITE_STREAM(GridPos)
TMPL_XML_READ_WRITE_STREAM(HitranRelaxationMatrixData)
TMPL_XML_READ_WRITE_STREAM(SpeciesIsotopologueRatios)
TMPL_XML_READ_WRITE_STREAM(MapOfErrorCorrectedSuddenData)
TMPL_XML_READ_WRITE_STREAM(MCAntenna)
TMPL_XML_READ_WRITE_STREAM(Ppath)
TMPL_XML_READ_WRITE_STREAM(PredefinedModelData)
TMPL_XML_READ_WRITE_STREAM(QuantumIdentifier)
TMPL_XML_READ_WRITE_STREAM(RetrievalQuantity)
TMPL_XML_READ_WRITE_STREAM(ScatteringMetaData)
TMPL_XML_READ_WRITE_STREAM(SLIData2)
TMPL_XML_READ_WRITE_STREAM(SingleScatteringData)
TMPL_XML_READ_WRITE_STREAM(SpeciesTag)
TMPL_XML_READ_WRITE_STREAM(Sun)
TMPL_XML_READ_WRITE_STREAM(SurfaceField)
TMPL_XML_READ_WRITE_STREAM(SurfacePoint)
TMPL_XML_READ_WRITE_STREAM(TelsemAtlas)
TMPL_XML_READ_WRITE_STREAM(TessemNN)
TMPL_XML_READ_WRITE_STREAM(XsecRecord)

//=== Array Types ==========================================================

TMPL_XML_READ_WRITE_STREAM(ArrayOfAbsorptionLines)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfAbsorptionLines)
TMPL_XML_READ_WRITE_STREAM(ArrayOfAgenda)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfArrayOfArrayOfGridPos)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfGriddedField1)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfGriddedField2)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfGriddedField3)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfGridPos)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfArrayOfGridPos)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfScatteringMetaData)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfSingleScatteringData)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfSpeciesTag)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfString)
TMPL_XML_READ_WRITE_STREAM(ArrayOfCIARecord)
TMPL_XML_READ_WRITE_STREAM(ArrayOfGriddedField1)
TMPL_XML_READ_WRITE_STREAM(ArrayOfGriddedField2)
TMPL_XML_READ_WRITE_STREAM(ArrayOfGriddedField3)
TMPL_XML_READ_WRITE_STREAM(ArrayOfGriddedField4)
TMPL_XML_READ_WRITE_STREAM(ArrayOfGridPos)
TMPL_XML_READ_WRITE_STREAM(ArrayOfJacobianTarget)
TMPL_XML_READ_WRITE_STREAM(ArrayOfQuantumIdentifier)
TMPL_XML_READ_WRITE_STREAM(ArrayOfPpath)
TMPL_XML_READ_WRITE_STREAM(ArrayOfRetrievalQuantity)
TMPL_XML_READ_WRITE_STREAM(ArrayOfScatteringMetaData)
TMPL_XML_READ_WRITE_STREAM(ArrayOfSingleScatteringData)
TMPL_XML_READ_WRITE_STREAM(ArrayOfSpeciesTag)
TMPL_XML_READ_WRITE_STREAM(ArrayOfSun)
TMPL_XML_READ_WRITE_STREAM(ArrayOfString)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTelsemAtlas)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTime)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfTime)
TMPL_XML_READ_WRITE_STREAM(ArrayOfXsecRecord)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfIndex)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfString)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfTensor3)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfTensor6)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfAtmPoint)
TMPL_XML_READ_WRITE_STREAM(ArrayOfIndex)
TMPL_XML_READ_WRITE_STREAM(ArrayOfMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfSparse)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTensor3)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTensor4)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTensor5)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTensor6)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTensor7)
TMPL_XML_READ_WRITE_STREAM(ArrayOfVector)


//=== Not storable Types ===================================================

TMPL_XML_READ_WRITE_STREAM(CallbackOperator)
TMPL_XML_READ_WRITE_STREAM(SpectralRadianceProfileOperator)
TMPL_XML_READ_WRITE_STREAM(NumericUnaryOperator)

//=== rtepack types ========================================================

TMPL_XML_READ_WRITE_STREAM(Propmat)
TMPL_XML_READ_WRITE_STREAM(PropmatVector)
TMPL_XML_READ_WRITE_STREAM(PropmatMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfPropmatVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfPropmatMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfPropmatVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfPropmatMatrix)

TMPL_XML_READ_WRITE_STREAM(Stokvec)
TMPL_XML_READ_WRITE_STREAM(StokvecVector)
TMPL_XML_READ_WRITE_STREAM(StokvecMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfStokvecVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfStokvecMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfStokvecVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfStokvecMatrix)

TMPL_XML_READ_WRITE_STREAM(Muelmat)
TMPL_XML_READ_WRITE_STREAM(MuelmatVector)
TMPL_XML_READ_WRITE_STREAM(MuelmatMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfMuelmatVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfMuelmatMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfMuelmatVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfMuelmatMatrix)

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_XML_READ_WRITE_STREAM

#endif /* xml_io_arts_types_h */
