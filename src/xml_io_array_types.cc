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
  \file   xml_io_array_types.cc
  \author Oliver Lemke <olemke@uni-bremen.de>
  \date   2003-06-11

  \brief This file contains basic functions to handle XML data files.

*/

#include "arts.h"
#include "xml_io_private.h"
#include "xml_io_basic_types.h"
#include "xml_io_compound_types.h"
#include "xml_io_array_types.h"


////////////////////////////////////////////////////////////////////////////
//   Overloaded functions for reading/writing data from/to XML stream
////////////////////////////////////////////////////////////////////////////

//=== Array<SpeciesRecord> ================================================

//! Reads SpeciesData from XML input stream
/*!
  \param is_xml    XML Input stream
  \param asrecord  SpeciesData return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      Array<SpeciesRecord>& asrecord,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "SpeciesRecord");

  tag.get_attribute_value ("nelem", nelem);
  asrecord.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, asrecord[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading SpeciesData: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes SpeciesData to XML output stream
/*!
  \param os_xml    XML Output stream
  \param asrecord  SpeciesData
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const Array<SpeciesRecord>& asrecord,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "SpeciesData");
  open_tag.add_attribute ("nelem", asrecord.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < asrecord.nelem (); n++)
    {
      xml_write_to_stream (os_xml, asrecord[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfArrayOfSpeciesTag ================================================

//! Reads ArrayOfArrayOfSpeciesTag from XML input stream
/*!
  \param is_xml  XML Input stream
  \param aastag  ArrayOfArrayOfSpeciesTag return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfArrayOfSpeciesTag& aastag,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfSpeciesTag");

  tag.get_attribute_value ("nelem", nelem);
  aastag.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, aastag[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfSpeciesTag: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfSpeciesTag to XML output stream
/*!
  \param os_xml  XML Output stream
  \param aastag  ArrayOfArrayOfSpeciesTag
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfArrayOfSpeciesTag& aastag,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "ArrayOfSpeciesTag");
  open_tag.add_attribute ("nelem", aastag.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aastag.nelem (); n++)
    {
      xml_write_to_stream (os_xml, aastag[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfArrayOfTensor3==================================================

//! Reads ArrayOfArrayOfTensor3 from XML input stream
/*!
  \param is_xml     XML Input stream
  \param aatensor3  ArrayOfArrayOfTensor3 return value
  \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfArrayOfTensor3& aatensor3,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfTensor3");

  tag.get_attribute_value ("nelem", nelem);
  aatensor3.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, aatensor3[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfTensor3: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfTensor3 to XML output stream
/*!
  \param os_xml     XML Output stream
  \param aatensor3  ArrayOfArrayOfTensor3
  \param pbofs      Pointer to binary file stream. NULL for ASCII output.
  \param name       Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfArrayOfTensor3& aatensor3,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "ArrayOfTensor3");
  open_tag.add_attribute ("nelem", aatensor3.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aatensor3.nelem (); n++)
    {
      xml_write_to_stream (os_xml, aatensor3[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfArrayOfTensor6==================================================

//! Reads ArrayOfArrayOfTensor6 from XML input stream
/*!
  \param is_xml     XML Input stream
  \param aatensor6  ArrayOfArrayOfTensor6 return value
  \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfArrayOfTensor6& aatensor6,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfTensor6");

  tag.get_attribute_value ("nelem", nelem);
  aatensor6.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, aatensor6[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfTensor6: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfTensor6 to XML output stream
/*!
  \param os_xml     XML Output stream
  \param aatensor6  ArrayOfArrayOfTensor6
  \param pbofs      Pointer to binary file stream. NULL for ASCII output.
  \param name       Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfArrayOfTensor6& aatensor6,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "ArrayOfTensor6");
  open_tag.add_attribute ("nelem", aatensor6.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aatensor6.nelem (); n++)
    {
      xml_write_to_stream (os_xml, aatensor6[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfGridPos =========================================================

//! Reads ArrayOfGridPos from XML input stream
/*!
  \param is_xml  XML Input stream
  \param agpos   ArrayOfGridPos return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfGridPos& agpos,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "GridPos");

  tag.get_attribute_value ("nelem", nelem);
  agpos.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, agpos[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfGridPos: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfGridPos to XML output stream
/*!
  \param os_xml  XML Output stream
  \param agpos   ArrayOfGridPos
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfGridPos& agpos,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "GridPos");
  open_tag.add_attribute ("nelem", agpos.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < agpos.nelem (); n++)
    {
      xml_write_to_stream (os_xml, agpos[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfIndex ===========================================================

//! Reads ArrayOfIndex from XML input stream
/*!
  \param is_xml  XML Input stream
  \param aindex  ArrayOfIndex return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfIndex& aindex,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Index");

  tag.get_attribute_value ("nelem", nelem);
  aindex.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, aindex[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfIndex: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfIndex to XML output stream
/*!
  \param os_xml  XML Output stream
  \param aindex  ArrayOfIndex
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfIndex& aindex,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "Index");
  open_tag.add_attribute ("nelem", aindex.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aindex.nelem (); n++)
    {
      xml_write_to_stream (os_xml, aindex[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}

//=== ArrayOfMatrix ==========================================================

//! Reads ArrayOfMatrix from XML input stream
/*!
  \param is_xml   XML Input stream
  \param amatrix  ArrayOfMatrix return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfMatrix& amatrix,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Matrix");

  tag.get_attribute_value ("nelem", nelem);
  amatrix.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, amatrix[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfMatrix: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfMatrix to XML output stream
/*!
  \param os_xml   XML Output stream
  \param amatrix  ArrayOfMatrix
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfMatrix& amatrix,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "Matrix");
  open_tag.add_attribute ("nelem", amatrix.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < amatrix.nelem (); n++)
    {
      xml_write_to_stream (os_xml, amatrix[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfSpeciesTag ================================================

//! Reads ArrayOfSpeciesTag from XML input stream
/*!
  \param is_xml  XML Input stream
  \param astag   ArrayOfSpeciesTag return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfSpeciesTag& astag,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "SpeciesTag");

  tag.get_attribute_value ("nelem", nelem);
  astag.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, astag[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfSpeciesTag: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfSpeciesTag to XML output stream
/*!
  \param os_xml  XML Output stream
  \param astag   ArrayOfSpeciesTag
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfSpeciesTag& astag,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "SpeciesTag");
  open_tag.add_attribute ("nelem", astag.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < astag.nelem (); n++)
    {
      xml_write_to_stream (os_xml, astag[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfSingleScatteringData===========================================

//! Reads ArrayOfSingleScatteringData from XML input stream
/*!
  \param is_xml   XML Input stream
  \param assdata  ArrayOfSingleScatteringData return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfSingleScatteringData& assdata,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "SingleScatteringData");

  tag.get_attribute_value ("nelem", nelem);
  assdata.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, assdata[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfSingleScatteringData: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfSpeciesTag to XML output stream
/*!
  \param os_xml   XML Output stream
  \param assdata  ArrayOfSpeciesTag
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfSingleScatteringData& assdata,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "SingleScatteringData");
  open_tag.add_attribute ("nelem", assdata.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < assdata.nelem (); n++)
    {
      xml_write_to_stream (os_xml, assdata[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfTensor3=========================================================

//! Reads ArrayOfTensor3 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param atensor3  ArrayOfTensor3 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfTensor3& atensor3,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Tensor3");

  tag.get_attribute_value ("nelem", nelem);
  atensor3.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, atensor3[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfTensor3: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfTensor3 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param atensor3  ArrayOfTensor3
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfTensor3& atensor3,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "Tensor3");
  open_tag.add_attribute ("nelem", atensor3.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < atensor3.nelem (); n++)
    {
      xml_write_to_stream (os_xml, atensor3[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfTensor6=========================================================

//! Reads ArrayOfTensor6 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param atensor6  ArrayOfTensor6 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfTensor6& atensor6,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Tensor6");

  tag.get_attribute_value ("nelem", nelem);
  atensor6.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, atensor6[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfTensor6: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfTensor6 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param atensor6  ArrayOfTensor6
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfTensor6& atensor6,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "Tensor6");
  open_tag.add_attribute ("nelem", atensor6.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < atensor6.nelem (); n++)
    {
      xml_write_to_stream (os_xml, atensor6[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfString ==========================================================

//! Reads ArrayOfString from XML input stream
/*!
  \param is_xml   XML Input stream
  \param astring  ArrayOfString return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfString& astring,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "String");

  tag.get_attribute_value ("nelem", nelem);
  astring.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, astring[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfString: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfString to XML output stream
/*!
  \param os_xml   XML Output stream
  \param astring  ArrayOfString
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfString& astring,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "String");
  open_tag.add_attribute ("nelem", astring.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < astring.nelem (); n++)
    {
      xml_write_to_stream (os_xml, astring[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}

//=== ArrayOfVector ==========================================================

//! Reads ArrayOfVector from XML input stream
/*!
  \param is_xml   XML Input stream
  \param avector  ArrayOfVector return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfVector& avector,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Vector");

  tag.get_attribute_value ("nelem", nelem);
  avector.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, avector[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfVector: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfVector to XML output stream
/*!
  \param os_xml   XML Output stream
  \param avector  ArrayOfVector
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfVector& avector,
                     bofstream *pbofs,
                     const String &name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "Vector");
  open_tag.add_attribute ("nelem", avector.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < avector.nelem (); n++)
    {
      xml_write_to_stream (os_xml, avector[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}

