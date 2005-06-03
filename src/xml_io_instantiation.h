/* Copyright (C) 2003 Oliver Lemke <olemke@uni-bremen.de>

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
  \author Oliver Lemke <olemke@uni-bremen.de>
  \date   2002-11-06

  \brief This file contains private function declarations and
         template instantiation to handle XML data files.

*/

#ifndef xml_io_instantiation_h
#define xml_io_instantiation_h

#include "xml_io.h"
#include <stdexcept>
#include <cfloat>
#include "xml_io_basic_types.h"
#include "xml_io_compound_types.h"
#include "xml_io_array_types.h"
#include "bifstream.h"
#include "bofstream.h"


////////////////////////////////////////////////////////////////////////////
//   Explicit instantiation of template functions we need
////////////////////////////////////////////////////////////////////////////

//=== Basic Types ==========================================================

template void
xml_read_from_file<Index> (const String&, Index&);

template void
xml_read_from_file<Matrix> (const String&, Matrix&);

template void
xml_read_from_file<Numeric> (const String&, Numeric&);

template void
xml_read_from_file<Sparse> (const String&, Sparse&);

template void
xml_read_from_file<String> (const String&, String&);

template void
xml_read_from_file<Tensor3> (const String&, Tensor3&);

template void
xml_read_from_file<Tensor4> (const String&, Tensor4&);

template void
xml_read_from_file<Tensor5> (const String&, Tensor5&);

template void
xml_read_from_file<Tensor6> (const String&, Tensor6&);

template void
xml_read_from_file<Tensor7> (const String&, Tensor7&);

template void
xml_read_from_file<Timer> (const String&, Timer&);

template void
xml_read_from_file<Vector> (const String&, Vector&);

template void
xml_write_to_file<Index> (const String&, const Index&, FileType);

template void
xml_write_to_file<Matrix> (const String&, const Matrix&, FileType);

template void
xml_write_to_file<Numeric> (const String&, const Numeric&, FileType);

template void
xml_write_to_file<Sparse> (const String&, const Sparse&, FileType);

template void
xml_write_to_file<String> (const String&, const String&, FileType);

template void
xml_write_to_file<Tensor3> (const String&, const Tensor3&, FileType);

template void
xml_write_to_file<Tensor4> (const String&, const Tensor4&, FileType);

template void
xml_write_to_file<Tensor5> (const String&, const Tensor5&, FileType);

template void
xml_write_to_file<Tensor6> (const String&, const Tensor6&, FileType);

template void
xml_write_to_file<Tensor7> (const String&, const Tensor7&, FileType);

template void
xml_write_to_file<Timer> (const String&, const Timer&, FileType);

template void
xml_write_to_file<Vector> (const String&, const Vector&, FileType);


//=== Compound Types =======================================================

template void
xml_read_from_file<Agenda> (const String&, Agenda&);

template void
xml_read_from_file<GasAbsLookup> (const String&, GasAbsLookup&);

template void
xml_read_from_file<GriddedField3> (const String&, GriddedField3&);

template void
xml_read_from_file<GridPos> (const String&, GridPos&);

template void
xml_read_from_file<IsotopeRecord> (const String&, IsotopeRecord&);

template void
xml_read_from_file<Ppath> (const String&, Ppath&);

template void
xml_read_from_file<RetrievalQuantity> (const String&, RetrievalQuantity&);

template void
xml_read_from_file<SingleScatteringData> (const String&, SingleScatteringData&);

template void
xml_read_from_file<SLIData2> (const String&, SLIData2&);

template void
xml_read_from_file<SpeciesRecord> (const String&, SpeciesRecord&);

template void
xml_read_from_file<SpeciesTag> (const String&, SpeciesTag&);

template void
xml_write_to_file<Agenda> (const String&, const Agenda&, FileType);

template void
xml_write_to_file<GasAbsLookup> (const String&, const GasAbsLookup&, FileType);

template void
xml_write_to_file<GriddedField3> (const String&, const GriddedField3&,
                                  FileType);

template void
xml_write_to_file<GridPos> (const String&, const GridPos&, FileType);

template void
xml_write_to_file<IsotopeRecord> (const String&, const IsotopeRecord&,
                                  FileType);

template void
xml_write_to_file<Ppath> (const String&, const Ppath&, FileType);

template void
xml_write_to_file<RetrievalQuantity> (const String&,
                                      const RetrievalQuantity&, FileType);

template void
xml_write_to_file<SingleScatteringData> (const String&,
                                         const SingleScatteringData&, FileType);

template void
xml_write_to_file<SLIData2> (const String&,
                             const SLIData2&, FileType);

template void
xml_write_to_file<SpeciesRecord> (const String&, const SpeciesRecord&,
                                  FileType);

template void
xml_write_to_file<SpeciesTag> (const String&, const SpeciesTag&, FileType);


//=== Array Types ==========================================================

template void
xml_read_from_file< Array<IsotopeRecord> > (const String&,
                                            Array<IsotopeRecord>&);

template void
xml_read_from_file< Array<SpeciesRecord> > (const String&,
                                            Array<SpeciesRecord>&);

template void
xml_read_from_file<ArrayOfArrayOfSpeciesTag> (const String&,
                                              ArrayOfArrayOfSpeciesTag&);

template void
xml_read_from_file<ArrayOfSingleScatteringData> (const String&,
                                              ArrayOfSingleScatteringData&);

template void
xml_read_from_file<ArrayOfGriddedField3> (const String&,
                                          ArrayOfGriddedField3&);

template void
xml_read_from_file<ArrayOfArrayOfTensor3> (const String&,
                                           ArrayOfArrayOfTensor3&);

template void
xml_read_from_file<ArrayOfArrayOfTensor6> (const String&,
                                           ArrayOfArrayOfTensor6&);

template void
xml_read_from_file<ArrayOfGridPos> (const String&, ArrayOfGridPos&);

template void
xml_read_from_file<ArrayOfArrayOfGridPos> 
(const String&, ArrayOfArrayOfGridPos&);

template void
xml_read_from_file<ArrayOfArrayOfArrayOfGridPos> 
(const String&, ArrayOfArrayOfArrayOfGridPos&);

template void
xml_read_from_file<ArrayOfArrayOfArrayOfArrayOfGridPos>
(const String&, ArrayOfArrayOfArrayOfArrayOfGridPos&);

template void
xml_read_from_file<ArrayOfIndex> (const String&, ArrayOfIndex&);

template void
xml_read_from_file<ArrayOfArrayOfIndex> (const String&, ArrayOfArrayOfIndex&);

template void
xml_read_from_file<ArrayOfMatrix> (const String&, ArrayOfMatrix&);

template void
xml_read_from_file<ArrayOfArrayOfMatrix> 
(const String&, ArrayOfArrayOfMatrix&);

template void
xml_read_from_file<ArrayOfRetrievalQuantity> (const String&,
                                              ArrayOfRetrievalQuantity&);

template void
xml_read_from_file<ArrayOfSpeciesTag> (const String&, ArrayOfSpeciesTag&);

template void
xml_read_from_file<ArrayOfTensor3> (const String&, ArrayOfTensor3&);

template void
xml_read_from_file<ArrayOfTensor6> (const String&, ArrayOfTensor6&);

template void
xml_read_from_file<ArrayOfTensor7> (const String&, ArrayOfTensor7&);

template void
xml_read_from_file<ArrayOfString> (const String&, ArrayOfString&);

template void
xml_read_from_file<ArrayOfVector> (const String&, ArrayOfVector&);

template void
xml_write_to_file<Array<IsotopeRecord> > (const String&,
                                          const Array<IsotopeRecord>&,
                                          FileType);

template void
xml_write_to_file<Array<SpeciesRecord> > (const String&,
                                          const Array<SpeciesRecord>&,
                                          FileType);

template void
xml_write_to_file<ArrayOfSingleScatteringData> (const String&,
                                            const ArrayOfSingleScatteringData&,
                                                FileType);

template void
xml_write_to_file<ArrayOfGriddedField3> (const String&,
                                            const ArrayOfGriddedField3&,
                                                FileType);

template void
xml_write_to_file<ArrayOfIndex> (const String&, const ArrayOfIndex&, FileType);

template void
xml_write_to_file<ArrayOfArrayOfIndex> (const String&, const ArrayOfArrayOfIndex&, FileType);

template void
xml_write_to_file<ArrayOfMatrix> (const String&, const ArrayOfMatrix&,
                                  FileType);

template void
xml_write_to_file<ArrayOfArrayOfMatrix>
(const String&, const ArrayOfArrayOfMatrix&, FileType);

template void
xml_write_to_file<ArrayOfRetrievalQuantity> (const String&,
                                             const ArrayOfRetrievalQuantity&,
                                             FileType);

template void
xml_write_to_file<ArrayOfTensor3> (const String&, const ArrayOfTensor3&,
                                   FileType);

template void
xml_write_to_file<ArrayOfArrayOfTensor3> (const String&,
                                          const ArrayOfArrayOfTensor3&,
                                          FileType);

template void
xml_write_to_file<ArrayOfTensor6> (const String&, const ArrayOfTensor6&,
                                   FileType);

template void
xml_write_to_file<ArrayOfTensor7> (const String&, const ArrayOfTensor7&,
                                   FileType);

template void
xml_write_to_file<ArrayOfArrayOfTensor6> (const String&,
                                          const ArrayOfArrayOfTensor6&,
                                          FileType);

template void
xml_write_to_file<ArrayOfString> (const String&, const ArrayOfString&,
                                  FileType);

template void
xml_write_to_file<ArrayOfVector> (const String&, const ArrayOfVector&,
                                  FileType);

template void
xml_write_to_file<ArrayOfArrayOfSpeciesTag> (const String&,
                                             const ArrayOfArrayOfSpeciesTag&,
                                             FileType);

template void
xml_write_to_file<ArrayOfGridPos> (const String&,
                                   const ArrayOfGridPos&,
                                   FileType);

template void
xml_write_to_file<ArrayOfArrayOfGridPos> (const String&,
                                   const ArrayOfArrayOfGridPos&,
                                   FileType);

template void
xml_write_to_file<ArrayOfArrayOfArrayOfGridPos> (const String&,
                                                 const ArrayOfArrayOfArrayOfGridPos&,
                                                 FileType);

template void
xml_write_to_file<ArrayOfArrayOfArrayOfArrayOfGridPos> (const String&,
                                   const ArrayOfArrayOfArrayOfArrayOfGridPos&,
                                   FileType);

#endif /* xml_io_instantiation_h */

