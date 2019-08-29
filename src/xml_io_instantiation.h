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
  \file   xml_io_instantiation.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2008-07-31

  \brief This file contains template instantiations to handle XML data files.

*/

#ifndef xml_io_instantiation_h
#define xml_io_instantiation_h

#include <stdexcept>
#include <cfloat>
#include "xml_io.h"
#include "xml_io_types.h"

#define TMPL_XML_READ_WRITE(what) \
  template void xml_read_from_file<what>(const String&, what&, const Verbosity&); \
  template void xml_write_to_file<what>(const String&, const what&, FileType, const Index, const Verbosity&);


////////////////////////////////////////////////////////////////////////////
//   Explicit instantiation of template functions we need
////////////////////////////////////////////////////////////////////////////

//=== Basic Types ==========================================================

TMPL_XML_READ_WRITE(Index)
TMPL_XML_READ_WRITE(Matrix)
TMPL_XML_READ_WRITE(Numeric)
TMPL_XML_READ_WRITE(Rational)
TMPL_XML_READ_WRITE(Sparse)
TMPL_XML_READ_WRITE(String)
TMPL_XML_READ_WRITE(Tensor3)
TMPL_XML_READ_WRITE(Tensor4)
TMPL_XML_READ_WRITE(Tensor5)
TMPL_XML_READ_WRITE(Tensor6)
TMPL_XML_READ_WRITE(Tensor7)
TMPL_XML_READ_WRITE(Timer)
TMPL_XML_READ_WRITE(Vector)
TMPL_XML_READ_WRITE(TransmissionMatrix)
TMPL_XML_READ_WRITE(RadiationVector)

//=== Compound Types =======================================================

TMPL_XML_READ_WRITE(Agenda)
TMPL_XML_READ_WRITE(CIARecord)
TMPL_XML_READ_WRITE(CovarianceMatrix)
TMPL_XML_READ_WRITE(GriddedField1)
TMPL_XML_READ_WRITE(GriddedField2)
TMPL_XML_READ_WRITE(GriddedField3)
TMPL_XML_READ_WRITE(GriddedField4)
TMPL_XML_READ_WRITE(GriddedField5)
TMPL_XML_READ_WRITE(GriddedField6)
TMPL_XML_READ_WRITE(GasAbsLookup)
TMPL_XML_READ_WRITE(GridPos)
TMPL_XML_READ_WRITE(IsotopologueRecord)
TMPL_XML_READ_WRITE(MCAntenna)
TMPL_XML_READ_WRITE(Ppath)
TMPL_XML_READ_WRITE(QuantumIdentifier)
TMPL_XML_READ_WRITE(QuantumNumbers)
TMPL_XML_READ_WRITE(RetrievalQuantity)
TMPL_XML_READ_WRITE(ScatteringMetaData)
TMPL_XML_READ_WRITE(SLIData2)
TMPL_XML_READ_WRITE(SingleScatteringData)
TMPL_XML_READ_WRITE(SpeciesAuxData)
TMPL_XML_READ_WRITE(SpeciesRecord)
TMPL_XML_READ_WRITE(SpeciesTag)
TMPL_XML_READ_WRITE(TelsemAtlas)
TMPL_XML_READ_WRITE(TessemNN)
TMPL_XML_READ_WRITE(XsecRecord)
TMPL_XML_READ_WRITE(Verbosity)

//=== Array Types ==========================================================

TMPL_XML_READ_WRITE(ArrayOfAgenda)
TMPL_XML_READ_WRITE(Array<IsotopologueRecord> )
TMPL_XML_READ_WRITE(Array<SpeciesRecord> )
TMPL_XML_READ_WRITE(ArrayOfArrayOfArrayOfArrayOfGridPos)
TMPL_XML_READ_WRITE(ArrayOfArrayOfGriddedField1)
TMPL_XML_READ_WRITE(ArrayOfArrayOfGriddedField2)
TMPL_XML_READ_WRITE(ArrayOfArrayOfGriddedField3)
TMPL_XML_READ_WRITE(ArrayOfArrayOfGridPos)
TMPL_XML_READ_WRITE(ArrayOfArrayOfIndex)
TMPL_XML_READ_WRITE(ArrayOfArrayOfLineRecord)
TMPL_XML_READ_WRITE(ArrayOfArrayOfMatrix)
TMPL_XML_READ_WRITE(ArrayOfArrayOfScatteringMetaData)
TMPL_XML_READ_WRITE(ArrayOfArrayOfSingleScatteringData)
TMPL_XML_READ_WRITE(ArrayOfArrayOfSpeciesTag)
TMPL_XML_READ_WRITE(ArrayOfArrayOfString)
TMPL_XML_READ_WRITE(ArrayOfArrayOfTensor3)
TMPL_XML_READ_WRITE(ArrayOfArrayOfTensor6)
TMPL_XML_READ_WRITE(ArrayOfArrayOfVector)
TMPL_XML_READ_WRITE(ArrayOfCIARecord)
TMPL_XML_READ_WRITE(ArrayOfGriddedField1)
TMPL_XML_READ_WRITE(ArrayOfGriddedField2)
TMPL_XML_READ_WRITE(ArrayOfGriddedField3)
TMPL_XML_READ_WRITE(ArrayOfGriddedField4)
TMPL_XML_READ_WRITE(ArrayOfGridPos)
TMPL_XML_READ_WRITE(ArrayOfIndex)
TMPL_XML_READ_WRITE(ArrayOfLineRecord)
TMPL_XML_READ_WRITE(ArrayOfLineshapeSpec)
TMPL_XML_READ_WRITE(ArrayOfMatrix)
TMPL_XML_READ_WRITE(ArrayOfSparse)
TMPL_XML_READ_WRITE(ArrayOfPpath)
TMPL_XML_READ_WRITE(ArrayOfQuantumIdentifier)
TMPL_XML_READ_WRITE(ArrayOfRetrievalQuantity)
TMPL_XML_READ_WRITE(ArrayOfScatteringMetaData)
TMPL_XML_READ_WRITE(ArrayOfSingleScatteringData)
TMPL_XML_READ_WRITE(ArrayOfSpeciesTag)
TMPL_XML_READ_WRITE(ArrayOfString)
TMPL_XML_READ_WRITE(ArrayOfTelsemAtlas)
TMPL_XML_READ_WRITE(ArrayOfTensor3)
TMPL_XML_READ_WRITE(ArrayOfTensor4)
TMPL_XML_READ_WRITE(ArrayOfTensor5)
TMPL_XML_READ_WRITE(ArrayOfTensor6)
TMPL_XML_READ_WRITE(ArrayOfTensor7)
TMPL_XML_READ_WRITE(ArrayOfVector)
TMPL_XML_READ_WRITE(ArrayOfTransmissionMatrix)
TMPL_XML_READ_WRITE(ArrayOfArrayOfTransmissionMatrix)
TMPL_XML_READ_WRITE(ArrayOfRadiationVector)
TMPL_XML_READ_WRITE(ArrayOfArrayOfRadiationVector)
TMPL_XML_READ_WRITE(PropagationMatrix)
TMPL_XML_READ_WRITE(ArrayOfPropagationMatrix)
TMPL_XML_READ_WRITE(ArrayOfArrayOfPropagationMatrix)
TMPL_XML_READ_WRITE(StokesVector)
TMPL_XML_READ_WRITE(ArrayOfStokesVector)
TMPL_XML_READ_WRITE(ArrayOfArrayOfStokesVector)
TMPL_XML_READ_WRITE(ArrayOfXsecRecord)

//==========================================================================

// Undefine the macro to avoid it being used anywhere else
#undef TMPL_XML_READ_WRITE

#endif /* xml_io_instantiation_h */

