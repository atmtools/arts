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
#include "optproperties.h"
#include "bifstream.h"
#include "bofstream.h"


////////////////////////////////////////////////////////////////////////////
//   Functions to open and read XML files
////////////////////////////////////////////////////////////////////////////

void
xml_open_output_file (ofstream& file, const String& name);

void
xml_open_input_file (ifstream& file, const String& name);


////////////////////////////////////////////////////////////////////////////
//   XML parser classes
////////////////////////////////////////////////////////////////////////////

//! XML attribute class
/*!
  Holds the name and value of an XML attribute.
*/

class XMLAttribute
{
public:
  String name;                  /*!< Attribute name */
  String value;                 /*!< Attribute value */
};


//! The ARTS XML tag class
/*!
  Handles reading, writing and constructing of XML tags.
*/
class ArtsXMLTag
{
public:

  String&
  get_name () { return (name); }

  void
  check_name (const String& expected_name);

  void
  set_name (const String& new_name) { name = new_name; }

  void
  add_attribute (const String& aname, const String& value);

  void
  add_attribute (const String& aname, const Index& value);

  void
  check_attribute (const String& aname, const String& value);

  void
  get_attribute_value (const String& aname, String& value);

  void
  get_attribute_value (const String& aname, Index& value);

  void
  read_from_stream (istream& is);

  void
  write_to_stream (ostream& os);

private:
  String name;                  /*!< Tag name */
  Array<XMLAttribute> attribs;  /*!< List of attributes */
};


////////////////////////////////////////////////////////////////////////////
//   General XML handling routines
////////////////////////////////////////////////////////////////////////////

void
xml_parse_error (const String& str_error);

void
xml_read_header_from_stream (istream& is, FileType &ftype,
                             NumericType &ntype, EndianType &etype);

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
xml_read_from_stream (istream&, ArrayOfArrayOfSpeciesTag&, bifstream *);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfSpeciesTag&, bofstream *);

void
xml_read_from_stream (istream&, ArrayOfArrayOfTensor3&, bifstream *);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfTensor3&, bofstream *);

void
xml_read_from_stream (istream&, ArrayOfArrayOfTensor6&, bifstream *);

void
xml_write_to_stream (ostream&, const ArrayOfArrayOfTensor6&, bofstream *);

void
xml_read_from_stream (istream&, ArrayOfGridPos&, bifstream *);

void
xml_write_to_stream (ostream&, const ArrayOfGridPos&, bofstream *);

void
xml_read_from_stream (istream&, ArrayOfIndex&, bifstream *);

void
xml_write_to_stream (ostream&, const ArrayOfIndex&, bofstream *);

void
xml_read_from_stream (istream&, ArrayOfMatrix&, bifstream *);

void
xml_write_to_stream (ostream&, const ArrayOfMatrix&, bofstream *);

void
xml_read_from_stream (istream&, ArrayOfSpeciesTag&, bifstream *);

void
xml_write_to_stream (ostream&, const ArrayOfSpeciesTag&, bofstream *);

void
xml_read_from_stream (istream&, ArrayOfTensor3&, bifstream *);

void
xml_write_to_stream (ostream&, const ArrayOfTensor3&, bofstream *);

void
xml_read_from_stream (istream&, ArrayOfTensor6&, bifstream *);

void
xml_write_to_stream (ostream&, const ArrayOfTensor6&, bofstream *);

void
xml_read_from_stream (istream&, ArrayOfVector&, bifstream *);

void
xml_write_to_stream (ostream&, const ArrayOfVector&, bofstream *);

void
xml_read_from_stream (istream&, GasAbsLookup&, bifstream *);

void
xml_write_to_stream (ostream&, const GasAbsLookup&, bofstream *);

void
xml_read_from_stream (istream&, GridPos&, bifstream *);

void
xml_write_to_stream (ostream&, const GridPos&, bofstream *);

void
xml_read_from_stream (istream&, Index&, bifstream *);

void
xml_write_to_stream (ostream&, const Index&, bofstream *);

void
xml_read_from_stream (istream&, Matrix&, bifstream *);

void
xml_write_to_stream (ostream&, const Matrix&, bofstream *);

void
xml_read_from_stream (istream&, Numeric&, bifstream *);

void
xml_write_to_stream (ostream&, const Numeric&, bofstream *);

void
xml_read_from_stream (istream&, SingleScatteringData&,
                      bifstream *);

void
xml_write_to_stream (ostream&, const SingleScatteringData&,
                     bofstream *);

void
xml_read_from_stream (istream&, SpeciesTag&, bifstream *);

void
xml_write_to_stream (ostream&, const SpeciesTag&, bofstream *);

void
xml_read_from_stream (istream&, String&, bifstream *);

void
xml_write_to_stream (ostream&, const String&, bofstream *);

void
xml_read_from_stream (istream&, Tensor3&, bifstream *);

void
xml_write_to_stream (ostream&, const Tensor3&, bofstream *);

void
xml_read_from_stream (istream&, Tensor4&, bifstream *);

void
xml_write_to_stream (ostream&, const Tensor4&, bofstream *);

void
xml_read_from_stream (istream&, Tensor5&, bifstream *);

void
xml_write_to_stream (ostream&, const Tensor5&, bofstream *);

void
xml_read_from_stream (istream&, Tensor6&, bifstream *);

void
xml_write_to_stream (ostream&, const Tensor6&, bofstream *);

void
xml_read_from_stream (istream&, Tensor7&, bifstream *);

void
xml_write_to_stream (ostream&, const Tensor7&, bofstream *);

void
xml_read_from_stream (istream&, Vector&, bifstream *);

void
xml_write_to_stream (ostream&, const Vector&, bofstream *);


////////////////////////////////////////////////////////////////////////////
//   Explicit instantiation of template functions we need
////////////////////////////////////////////////////////////////////////////

template void
xml_read_from_file<Agenda> (const String&, Agenda&);

template void
xml_read_from_file<ArrayOfArrayOfSpeciesTag> (const String&,
                                              ArrayOfArrayOfSpeciesTag&);

template void
xml_read_from_file<ArrayOfArrayOfTensor3> (const String&,
                                           ArrayOfArrayOfTensor3&);

template void
xml_read_from_file<ArrayOfArrayOfTensor6> (const String&, 
                                           ArrayOfArrayOfTensor6&);

template void
xml_read_from_file<ArrayOfGridPos> (const String&, ArrayOfGridPos&);

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
xml_read_from_file<SingleScatteringData> (const String&, SingleScatteringData&);

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
xml_write_to_file<Agenda> (const String&, const Agenda&, FileType);

template void
xml_write_to_file<ArrayOfGridPos> (const String&, const ArrayOfGridPos&,
                                   FileType);

template void
xml_write_to_file<ArrayOfIndex> (const String&, const ArrayOfIndex&, FileType);

template void
xml_write_to_file<ArrayOfMatrix> (const String&, const ArrayOfMatrix&,
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
xml_write_to_file<GasAbsLookup> (const String&, const GasAbsLookup&, FileType);

template void
xml_write_to_file<GridPos> (const String&, const GridPos&, FileType);

template void
xml_write_to_file<Index> (const String&, const Index&, FileType);

template void
xml_write_to_file<Matrix> (const String&, const Matrix&, FileType);

template void
xml_write_to_file<Numeric> (const String&, const Numeric&, FileType);

template void
xml_write_to_file<Ppath> (const String&, const Ppath&, FileType);

template void
xml_write_to_file<SingleScatteringData> (const String&,
                                         const SingleScatteringData&, FileType);

template void
xml_write_to_file<SpeciesTag> (const String&, const SpeciesTag&, FileType);

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
xml_write_to_file<Vector> (const String&, const Vector&, FileType);

template void
xml_write_to_file<ArrayOfArrayOfSpeciesTag> (const String&,
                                             const ArrayOfArrayOfSpeciesTag&,
                                             FileType);

#endif  /* xml_io_private_h */
