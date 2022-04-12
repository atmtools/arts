/* Copyright (C) 2003-2012 Oliver Lemke <olemke@core-dump.info>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

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
#include "xml_io_general_types.h"
#include "absorption.h"
#include "agenda_class.h"
#include "array.h"
#include "artstime.h"
#include "bifstream.h"
#include "bofstream.h"
#include "cia.h"
#include "covariance_matrix.h"
#include "gas_abs_lookup.h"
#include "gridded_fields.h"
#include "hitran_xsec.h"
#include "jacobian.h"
#include "linemixing_hitran.h"
#include "m_general.h"
#include "matpackII.h"
#include "matpackVII.h"
#include "mc_antenna.h"
#include "mc_interp.h"
#include "messages.h"
#include "optproperties.h"
#include "template_partfun.h"
#include "ppath.h"
#include "propagationmatrix.h"
#include "telsem.h"
#include "tessem.h"
#include "transmissionmatrix.h"
#include "hitran_xsec.h"
#include "absorptionlines.h"
#include "linemixing.h"
#include "callback.h"

#define TMPL_XML_READ_WRITE_STREAM(what)                  \
  void xml_read_from_stream(                              \
      istream &, what &, bifstream *, const Verbosity &); \
  void xml_write_to_stream(ostream &,                     \
                           const what &,                  \
                           bofstream *,                   \
                           const String &,                \
                           const Verbosity &);

////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for XML streams
////////////////////////////////////////////////////////////////////////////

//=== Basic Types ==========================================================

TMPL_XML_READ_WRITE_STREAM(PartitionFunctionsData)
TMPL_XML_READ_WRITE_STREAM(JacobianTarget)
TMPL_XML_READ_WRITE_STREAM(Rational)
TMPL_XML_READ_WRITE_STREAM(Timer)
TMPL_XML_READ_WRITE_STREAM(Time)

//=== Compound Types =======================================================

TMPL_XML_READ_WRITE_STREAM(AbsorptionLines)
TMPL_XML_READ_WRITE_STREAM(Agenda)
TMPL_XML_READ_WRITE_STREAM(CIARecord)
TMPL_XML_READ_WRITE_STREAM(CovarianceMatrix)
TMPL_XML_READ_WRITE_STREAM(EnergyLevelMap)
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
TMPL_XML_READ_WRITE_STREAM(QuantumIdentifier)
TMPL_XML_READ_WRITE_STREAM(RetrievalQuantity)
TMPL_XML_READ_WRITE_STREAM(ScatteringMetaData)
TMPL_XML_READ_WRITE_STREAM(SLIData2)
TMPL_XML_READ_WRITE_STREAM(SingleScatteringData)
TMPL_XML_READ_WRITE_STREAM(SpeciesTag)
TMPL_XML_READ_WRITE_STREAM(TelsemAtlas)
TMPL_XML_READ_WRITE_STREAM(TessemNN)
TMPL_XML_READ_WRITE_STREAM(XsecRecord)
TMPL_XML_READ_WRITE_STREAM(Verbosity)

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
TMPL_XML_READ_WRITE_STREAM(ArrayOfString)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTelsemAtlas)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTime)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfTime)
TMPL_XML_READ_WRITE_STREAM(PropagationMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfPropagationMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfPropagationMatrix)
TMPL_XML_READ_WRITE_STREAM(TransmissionMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTransmissionMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfTransmissionMatrix)
TMPL_XML_READ_WRITE_STREAM(StokesVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfStokesVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfStokesVector)
TMPL_XML_READ_WRITE_STREAM(RadiationVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfRadiationVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfRadiationVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfXsecRecord)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfIndex)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfString)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfTensor3)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfTensor6)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfVector)
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

TMPL_XML_READ_WRITE_STREAM(CallbackFunction)

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_XML_READ_WRITE_STREAM

#endif /* xml_io_arts_types_h */
