/* Copyright (C) 2002 Oliver Lemke <olemke@uni-bremen.de>

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
  \file   xml_io_private.h
  \author Oliver Lemke <olemke@uni-bremen.de>
  \date   2002-11-06

  \brief This file contains private function declarations and
         template instantiation to handle XML data files.

*/

#ifndef xml_io_private_h
#define xml_io_private_h

#include "xml_io.h"
#include <stdexcept>
#include <cfloat>
#include "matpackI.h"
#include "matpackIII.h"
#include "matpackIV.h"
#include "matpackV.h"
#include "matpackVI.h"
#include "matpackVII.h"
#include "array.h"
#include "messages.h"
#include "ppath.h"
#include "agenda_class.h"
#include "absorption.h"
#include "gas_abs_lookup.h"
#include "xml_io_private.h"


////////////////////////////////////////////////////////////////////////////
//   General XML handling routines
////////////////////////////////////////////////////////////////////////////

void
xml_parse_error (const String& str_error);

void
xml_read_header_from_stream (istream& is);

void
xml_read_footer_from_stream (istream& is);

void
xml_write_header_to_stream (ostream& os);

void
xml_write_footer_to_stream (ostream& os);

void
xml_set_stream_precision (ostream &os);


////////////////////////////////////////////////////////////////////////////
//   Overloaded reading/writing routines for XML streams
////////////////////////////////////////////////////////////////////////////

void
xml_read_from_stream (istream& is, ArrayOfArrayOfSpeciesTag& aastag);

void
xml_write_to_stream (ostream& os, const ArrayOfArrayOfSpeciesTag& aastag);

void
xml_read_from_stream (istream& is, ArrayOfArrayOfTensor3& aatensor3);

void
xml_write_to_stream (ostream& os, const ArrayOfArrayOfTensor3& aatensor3);

void
xml_read_from_stream (istream& is, ArrayOfArrayOfTensor6& aatensor6);

void
xml_write_to_stream (ostream& os, const ArrayOfArrayOfTensor6& aatensor6);

void
xml_read_from_stream (istream& is, ArrayOfIndex& aindex);

void
xml_write_to_stream (ostream& os, const ArrayOfIndex& aindex);

void
xml_read_from_stream (istream& is, ArrayOfMatrix& amatrix);

void
xml_write_to_stream (ostream& os, const ArrayOfMatrix& amatrix);

void
xml_read_from_stream (istream& is, ArrayOfSpeciesTag& astag);

void
xml_write_to_stream (ostream& os, const ArrayOfSpeciesTag& astag);

void
xml_read_from_stream (istream& is, ArrayOfTensor3& atensor3);

void
xml_write_to_stream (ostream& os, const ArrayOfTensor3& atensor3);

void
xml_read_from_stream (istream& is, ArrayOfTensor6& atensor6);

void
xml_write_to_stream (ostream& os, const ArrayOfTensor6& atensor6);

void
xml_read_from_stream (istream& is, ArrayOfVector& avector);

void
xml_write_to_stream (ostream& os, const ArrayOfVector& avector);

void
xml_read_from_stream (istream& is, GasAbsLookup& gal);

void
xml_write_to_stream (ostream& os, const GasAbsLookup& gal);

void
xml_read_from_stream (istream& is, GridPos& gp);

void
xml_write_to_stream (ostream& os, const GridPos& gp);

void
xml_read_from_stream (istream& is, Index& index);

void
xml_write_to_stream (ostream& os, const Index& index);

void
xml_read_from_stream (istream& is, Matrix& matrix);

void
xml_write_to_stream (ostream& os, const Matrix& matrix);

void
xml_read_from_stream (istream& is, Numeric& numeric);

void
xml_write_to_stream (ostream& os, const Numeric& numeric);

void
xml_read_from_stream (istream& is, SpeciesTag& stag);

void
xml_write_to_stream (ostream& os, const SpeciesTag& stag);

void
xml_read_from_stream (istream& is, String& str);

void
xml_write_to_stream (ostream& os, const String& str);

void
xml_read_from_stream (istream& is, Tensor3& tensor);

void
xml_write_to_stream (ostream& os, const Tensor3& tensor);

void
xml_read_from_stream (istream& is, Tensor4& tensor);

void
xml_write_to_stream (ostream& os, const Tensor4& tensor);

void
xml_read_from_stream (istream& is, Tensor5& tensor);

void
xml_write_to_stream (ostream& os, const Tensor5& tensor);

void
xml_read_from_stream (istream& is, Tensor6& tensor);

void
xml_write_to_stream (ostream& os, const Tensor6& tensor);

void
xml_read_from_stream (istream& is, Tensor7& tensor);

void
xml_write_to_stream (ostream& os, const Tensor7& tensor);

void
xml_read_from_stream (istream& is, Vector& vector);

void
xml_write_to_stream (ostream& os, const Vector& vector);


////////////////////////////////////////////////////////////////////////////
//   Explicit instantiation of template functions we need
////////////////////////////////////////////////////////////////////////////

template void
xml_read_from_file<Agenda> (const String&, Agenda&);

template void
xml_read_from_file<ArrayOfArrayOfSpeciesTag> (const String&, ArrayOfArrayOfSpeciesTag&);

template void
xml_read_from_file<ArrayOfArrayOfTensor3> (const String&,
                                           ArrayOfArrayOfTensor3&);

template void
xml_read_from_file<ArrayOfArrayOfTensor6> (const String&, 
                                           ArrayOfArrayOfTensor6&);

template void
xml_read_from_file<ArrayOfIndex> (const String&, ArrayOfIndex&);

template void
xml_read_from_file<ArrayOfMatrix> (const String&, ArrayOfMatrix&);

template void
xml_read_from_file<ArrayOfSpeciesTag> (const String&, ArrayOfSpeciesTag&);

template void
xml_read_from_file<ArrayOfTensor3> (const String&, ArrayOfTensor3&);

template void
xml_read_from_file<ArrayOfTensor6> (const String&, ArrayOfTensor6&);

template void
xml_read_from_file<ArrayOfString> (const String&, ArrayOfString&);

template void
xml_read_from_file<ArrayOfVector> (const String&, ArrayOfVector&);

template void
xml_read_from_file<GasAbsLookup> (const String&, GasAbsLookup&);

template void
xml_read_from_file<GridPos> (const String&, GridPos&);

template void
xml_read_from_file<Index> (const String&, Index&);

template void
xml_read_from_file<Matrix> (const String&, Matrix&);

template void
xml_read_from_file<Numeric> (const String&, Numeric&);

template void
xml_read_from_file<Ppath> (const String&, Ppath&);

template void
xml_read_from_file<SpeciesTag> (const String&, SpeciesTag&);

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
xml_read_from_file<Vector> (const String&, Vector&);

template void
xml_write_to_file<Agenda> (const String&, const Agenda&);

template void
xml_write_to_file<ArrayOfIndex> (const String&, const ArrayOfIndex&);

template void
xml_write_to_file<ArrayOfMatrix> (const String&, const ArrayOfMatrix&);

template void
xml_write_to_file<ArrayOfTensor3> (const String&, const ArrayOfTensor3&);

template void
xml_write_to_file<ArrayOfArrayOfTensor3> (const String&,
                                          const ArrayOfArrayOfTensor3&);

template void
xml_write_to_file<ArrayOfTensor6> (const String&, const ArrayOfTensor6&);

template void
xml_write_to_file<ArrayOfArrayOfTensor6> (const String&, 
                                          const ArrayOfArrayOfTensor6&);

template void
xml_write_to_file<ArrayOfString> (const String&, const ArrayOfString&);

template void
xml_write_to_file<ArrayOfVector> (const String&, const ArrayOfVector&);

template void
xml_write_to_file<GasAbsLookup> (const String&, const GasAbsLookup&);

template void
xml_write_to_file<GridPos> (const String&, const GridPos&);

template void
xml_write_to_file<Index> (const String&, const Index&);

template void
xml_write_to_file<Matrix> (const String&, const Matrix&);

template void
xml_write_to_file<Numeric> (const String&, const Numeric&);

template void
xml_write_to_file<Ppath> (const String&, const Ppath&);

template void
xml_write_to_file<String> (const String&, const String&);

template void
xml_write_to_file<Tensor3> (const String&, const Tensor3&);

template void
xml_write_to_file<Tensor4> (const String&, const Tensor4&);

template void
xml_write_to_file<Tensor5> (const String&, const Tensor5&);

template void
xml_write_to_file<Tensor6> (const String&, const Tensor6&);

template void
xml_write_to_file<Tensor7> (const String&, const Tensor7&);

template void
xml_write_to_file<Vector> (const String&, const Vector&);

template void
xml_write_to_file<ArrayOfArrayOfSpeciesTag> (const String&, const ArrayOfArrayOfSpeciesTag&);

#endif  /* xml_io_private_h */
