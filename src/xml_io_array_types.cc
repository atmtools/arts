/* Copyright (C) 2003-2008 Oliver Lemke <olemke@core-dump.info>

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
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains basic functions to handle XML data files.

*/

#include "arts.h"
#include "xml_io_private.h"
#include "xml_io_types.h"
#include "jacobian.h"


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
  tag.check_attribute ("type", "SpeciesData");

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
                     const String& name)
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
                     const String& name)
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



//=== ArrayOfPpath =========================================================

//! Reads ArrayOfPpath from XML input stream
/*!
  \param is_xml  XML Input stream
  \param appath  ArrayOfPpath return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfPpath& appath,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Ppath");

  tag.get_attribute_value ("nelem", nelem);
  appath.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, appath[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfPpath: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}



//! Writes ArrayOfPpath to XML output stream
/*!
  \param os_xml  XML Output stream
  \param appath   ArrayOfPpath
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfPpath& appath,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "Ppath");
  open_tag.add_attribute ("nelem", appath.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < appath.nelem (); n++)
    {
      xml_write_to_stream (os_xml, appath[n], pbofs);
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
                     const String& name)
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
                     const String& name)
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
                     const String& name)
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

//=== ArrayOfArrayOfGridPos =====================================

//! Reads ArrayOfArrayOfGridPos from XML input stream
/*!
  \param is_xml  XML Input stream
  \param aagpos  ArrayOfArrayOfGridPos return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfArrayOfGridPos& aagpos,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfGridPos");

  tag.get_attribute_value ("nelem", nelem);
  aagpos.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, aagpos[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfGridPos: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfGridPos to XML output stream
/*!
  \param os_xml  XML Output stream
  \param aagpos  ArrayOfArrayOfGridPos
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfArrayOfGridPos& aagpos,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "ArrayOfGridPos");
  open_tag.add_attribute ("nelem", aagpos.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aagpos.nelem (); n++)
    {
      xml_write_to_stream (os_xml, aagpos[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfArrayOfGridPos =====================================

//! Reads ArrayOfArrayOfArrayOfGridPos from XML input stream
/*!
  \param is_xml  XML Input stream
  \param aaagpos ArrayOfArrayOfArrayOfGridPos return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfArrayOfArrayOfGridPos& aaagpos,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfArrayOfGridPos");

  tag.get_attribute_value ("nelem", nelem);
  aaagpos.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, aaagpos[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfArrayOfGridPos: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfArrayOfGridPos to XML output stream
/*!
  \param os_xml  XML Output stream
  \param aaagpos ArrayOfArrayOfArrayOfGridPos
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfArrayOfArrayOfGridPos& aaagpos,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "ArrayOfArrayOfGridPos");
  open_tag.add_attribute ("nelem", aaagpos.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aaagpos.nelem (); n++)
    {
      xml_write_to_stream (os_xml, aaagpos[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfArrayOfArrayOfArrayOfGridPos =====================================

//! Reads ArrayOfArrayOfArrayOfArrayOfGridPos from XML input stream
/*!
  \param is_xml  XML Input stream
  \param aaaagpos ArrayOfArrayOfArrayOfArrayOfGridPos return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfArrayOfArrayOfArrayOfGridPos& aaaagpos,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfArrayOfArrayOfGridPos");

  tag.get_attribute_value ("nelem", nelem);
  aaaagpos.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, aaaagpos[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfArrayOfArrayOfGridPos: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfArrayOfArrayOfGridPos to XML output stream
/*!
  \param os_xml  XML Output stream
  \param aaaagpos   ArrayOfArrayOfArrayOfArrayOfGridPos
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfArrayOfArrayOfArrayOfGridPos& aaaagpos,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "ArrayOfArrayOfArrayOfGridPos");
  open_tag.add_attribute ("nelem", aaaagpos.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aaaagpos.nelem (); n++)
    {
      xml_write_to_stream (os_xml, aaaagpos[n], pbofs);
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
                     const String& name)
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



//=== ArrayOfArrayOfIndex =====================================================

//! Reads ArrayOfArrayOfIndex from XML input stream
/*!
  \param is_xml     XML Input stream
  \param aaindex    ArrayOfArrayOfIndex return value
  \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfArrayOfIndex& aaindex,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfIndex");

  tag.get_attribute_value ("nelem", nelem);
  aaindex.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, aaindex[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfIndex: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfIndex to XML output stream
/*!
  \param os_xml     XML Output stream
  \param aaindex    ArrayOfArrayOfIndex
  \param pbofs      Pointer to binary file stream. NULL for ASCII output.
  \param name       Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfArrayOfIndex& aaindex,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "ArrayOfIndex");
  open_tag.add_attribute ("nelem", aaindex.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aaindex.nelem (); n++)
    {
      xml_write_to_stream (os_xml, aaindex[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}




//=== ArrayOfIsotopeRecord ===================================================

//! Reads ArrayOfIsotopeRecord from XML input stream
/*!
  \param is_xml    XML Input stream
  \param airecord  ArrayOfIsotopeRecord return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      Array<IsotopeRecord>& airecord,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "IsotopeRecord");

  tag.get_attribute_value ("nelem", nelem);
  airecord.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, airecord[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfIsotopeRecord: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfIsotopeRecord to XML output stream
/*!
  \param os_xml    XML Output stream
  \param airecord  ArrayOfIsotopeRecord
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const Array<IsotopeRecord>& airecord,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "IsotopeRecord");
  open_tag.add_attribute ("nelem", airecord.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < airecord.nelem (); n++)
    {
      xml_write_to_stream (os_xml, airecord[n], pbofs);
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
                     const String& name)
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


//=== ArrayOfArrayOfMatrix====================================================

//! Reads ArrayOfArrayOfMatrix from XML input stream
/*!
  \param is_xml     XML Input stream
  \param aamatrix   ArrayOfArrayOfMatrix return value
  \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfArrayOfMatrix& aamatrix,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfMatrix");

  tag.get_attribute_value ("nelem", nelem);
  aamatrix.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, aamatrix[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfMatrix: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfMatrix to XML output stream
/*!
  \param os_xml     XML Output stream
  \param aamatrix   ArrayOfArrayOfMatrix
  \param pbofs      Pointer to binary file stream. NULL for ASCII output.
  \param name       Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfArrayOfMatrix& aamatrix,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "ArrayOfMatrix");
  open_tag.add_attribute ("nelem", aamatrix.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aamatrix.nelem (); n++)
    {
      xml_write_to_stream (os_xml, aamatrix[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfRetrievalQuantity =======================================

//! Reads ArrayOfRetrievalQuantity from XML input stream
/*!
  \param is_xml    XML Input stream
  \param arq       ArrayOfRetrievalQuantity return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfRetrievalQuantity& arq,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "RetrievalQuantity");

  tag.get_attribute_value ("nelem", nelem);
  arq.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, arq[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfRetrievalQuantity: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfRetrivalQuantity to XML output stream
/*!
  \param os_xml    XML Output stream
  \param arq       ArrayOfRetrievalQuantity
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfRetrievalQuantity& arq,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "RetrievalQuantity");
  open_tag.add_attribute ("nelem", arq.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < arq.nelem (); n++)
    {
      xml_write_to_stream (os_xml, arq[n], pbofs);
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
                     const String& name)
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

//! Writes ArrayOfSingleScatteringData to XML output stream
/*!
  \param os_xml   XML Output stream
  \param assdata  ArrayOfSingleScatteringData
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfSingleScatteringData& assdata,
                     bofstream *pbofs,
                     const String& name)
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


//=== ArrayOfGField1 ===========================================

//! Reads ArrayOfGField1 from XML input stream
/*!
  \param is_xml   XML Input stream
  \param agfield  ArrayOfGField1 return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfGField1& agfield,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "GriddedField1");

  tag.get_attribute_value ("nelem", nelem);
  agfield.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, agfield[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfGField1: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfGField1 to XML output stream
/*!
  \param os_xml   XML Output stream
  \param agfield  ArrayOfGField1
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfGField1& agfield,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "GriddedField1");
  open_tag.add_attribute ("nelem", agfield.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < agfield.nelem (); n++)
    {
      xml_write_to_stream (os_xml, agfield[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfGField2 ===========================================

//! Reads ArrayOfGField2 from XML input stream
/*!
  \param is_xml   XML Input stream
  \param agfield  ArrayOfGField2 return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfGField2& agfield,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "GriddedField2");

  tag.get_attribute_value ("nelem", nelem);
  agfield.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, agfield[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfGField2: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfGField2 to XML output stream
/*!
  \param os_xml   XML Output stream
  \param agfield  ArrayOfGField2
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfGField2& agfield,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "GriddedField2");
  open_tag.add_attribute ("nelem", agfield.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < agfield.nelem (); n++)
    {
      xml_write_to_stream (os_xml, agfield[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfGField3 ===========================================

//! Reads ArrayOfGField3 from XML input stream
/*!
  \param is_xml   XML Input stream
  \param agfield  ArrayOfGField3 return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfGField3& agfield,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "GriddedField3");

  tag.get_attribute_value ("nelem", nelem);
  agfield.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, agfield[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfGField3: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfGField3 to XML output stream
/*!
  \param os_xml   XML Output stream
  \param agfield  ArrayOfGField3
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfGField3& agfield,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "GriddedField3");
  open_tag.add_attribute ("nelem", agfield.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < agfield.nelem (); n++)
    {
      xml_write_to_stream (os_xml, agfield[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfArrayOfGField1 ===========================================

//! Reads ArrayOfArrayOfGField1 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param aagfield  ArrayOfArrayOfGField1 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfArrayOfGField1& aagfield,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfGriddedField1");

  tag.get_attribute_value ("nelem", nelem);
  aagfield.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, aagfield[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfGField1: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfGField1 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param aagfield  ArrayOfArrayOfGField1
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfArrayOfGField1& aagfield,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "ArrayGriddedField1");
  open_tag.add_attribute ("nelem", aagfield.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aagfield.nelem (); n++)
    {
      xml_write_to_stream (os_xml, aagfield[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfArrayOfGField3 ===========================================

//! Reads ArrayOfArrayOfGField3 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param aagfield  ArrayOfArrayOfGField3 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfArrayOfGField3& aagfield,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfGriddedField3");

  tag.get_attribute_value ("nelem", nelem);
  aagfield.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, aagfield[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfGField3: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfGField3 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param aagfield  ArrayOfArrayOfGField3
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfArrayOfGField3& aagfield,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "ArrayGriddedField3");
  open_tag.add_attribute ("nelem", aagfield.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aagfield.nelem (); n++)
    {
      xml_write_to_stream (os_xml, aagfield[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfGField4 ===========================================

//! Reads ArrayOfGField4 from XML input stream
/*!
  \param is_xml   XML Input stream
  \param agfield  ArrayOfGField4 return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfGField4& agfield,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "GriddedField4");

  tag.get_attribute_value ("nelem", nelem);
  agfield.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, agfield[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfGField4: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfGField4 to XML output stream
/*!
  \param os_xml   XML Output stream
  \param agfield  ArrayOfGField4
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfGField4& agfield,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "GriddedField4");
  open_tag.add_attribute ("nelem", agfield.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < agfield.nelem (); n++)
    {
      xml_write_to_stream (os_xml, agfield[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfLineRecord ===========================================

//! Reads ArrayOfLineRecord from XML input stream
/*!
  \param is_xml   XML Input stream
  \param alrecord ArrayOfLineRecord return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfLineRecord& alrecord,
                      bifstream * pbifs _U_)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("ArrayOfLineRecord");

  tag.get_attribute_value ("nelem", nelem);

  LineRecord dummy_line_record;
  String version;
  tag.get_attribute_value ("version", version);

  if (version != dummy_line_record.Version())
    {
      ostringstream os;

      if (9 <= version.nelem())
        {
          if ("ARTSCAT" == version.substr (0,7))
            {
              os << "The ARTS line file you are trying contains a version tag\n"
                << "different from the current version.\n"
                << "Tag in file:     " << version << "\n"
                << "Current version: " << dummy_line_record.Version();
              throw runtime_error (os.str());
            }
        }

      os << "The ARTS line file you are trying to read does not contain a valid version tag.\n"
        << "Probably it was created with an older version of ARTS that used different units.";
      throw runtime_error (os.str());
    }

  alrecord.resize (0);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          LineRecord lr;
          if (lr.ReadFromArtsStream (is_xml))
            throw runtime_error ("Cannot read line from file");

          alrecord.push_back (lr);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfLineRecord: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/ArrayOfLineRecord");
}


//! Reads ArrayOfLineRecord from XML input stream within specified frequency
//! range
/*!
  \param is_xml   XML Input stream
  \param alrecord ArrayOfLineRecord return value
  \param fmin     Lowest frequency
  \param fmax     Highest frequency
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfLineRecord& alrecord,
                      const Numeric fmin,
                      const Numeric fmax,
                      bifstream * pbifs _U_)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("ArrayOfLineRecord");

  tag.get_attribute_value ("nelem", nelem);

  LineRecord dummy_line_record;
  String version;
  tag.get_attribute_value ("version", version);

  if (version != dummy_line_record.Version())
    {
      ostringstream os;

      if (9 <= version.nelem())
        {
          if ("ARTSCAT" == version.substr (0,7))
            {
              os << "The ARTS line file you are trying contains a version tag\n"
                << "different from the current version.\n"
                << "Tag in file:     " << version << "\n"
                << "Current version: " << dummy_line_record.Version();
              throw runtime_error (os.str());
            }
        }

      os << "The ARTS line file you are trying to read does not contain a valid version tag.\n"
        << "Probably it was created with an older version of ARTS that used different units.";
      throw runtime_error (os.str());
    }

  alrecord.resize (0);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          LineRecord lr;
          if (lr.ReadFromArtsStream (is_xml))
            throw runtime_error ("Cannot read line from file");

          if ( fmin <= lr.F() && lr.F() <= fmax )
            alrecord.push_back (lr);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfLineRecord: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is_xml);
  tag.check_name ("/ArrayOfLineRecord");
}


//! Writes ArrayOfLineRecord to XML output stream
/*!
  \param os_xml   XML Output stream
  \param alrecord ArrayOfLineRecord
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfLineRecord& alrecord,
                     bofstream * pbofs _U_,
                     const String& name)

{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;
  LineRecord dummy_line_record;

  open_tag.set_name ("ArrayOfLineRecord");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("version", dummy_line_record.Version ());
  open_tag.add_attribute ("nelem", alrecord.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for ( Index n = 0; n < alrecord.nelem(); n++ )
    {
      os_xml << alrecord[n] << "\n";
    }

  close_tag.set_name ("/ArrayOfLineRecord");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
} 


//=== ArrayOfArrayOfLineRecord ===========================================

//! Reads ArrayOfArrayOfLineRecord from XML input stream
/*!
  \param is_xml   XML Input stream
  \param aalrecord ArrayOfArrayOfLineRecord return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfArrayOfLineRecord& aalrecord,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfLineRecord");

  tag.get_attribute_value ("nelem", nelem);
  aalrecord.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, aalrecord[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfLineRecord: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfLineRecord to XML output stream
/*!
  \param os_xml   XML Output stream
  \param aalrecord ArrayOfArrayOfLineRecord
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfArrayOfLineRecord& aalrecord,
                     bofstream *pbofs,
                     const String& name)

{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "ArrayOfLineRecord");
  open_tag.add_attribute ("nelem", aalrecord.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aalrecord.nelem (); n++)
    {
      xml_write_to_stream (os_xml, aalrecord[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
} 


//=== ArrayOfLineshapeSpec ===========================================

//! Reads ArrayOfLineshapeSpec from XML input stream
/*!
  \param is_xml   XML Input stream
  \param alspec   ArrayOfLineshapeSpec return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml _U_,
                      ArrayOfLineshapeSpec& alspec _U_,
                      bifstream * pbifs _U_)
{
  // FIXME OLE: Implement this.
  throw runtime_error ("Boo. Not yet implemented.");
}


//! Writes ArrayOfLineshapeSpec to XML output stream
/*!
  \param os_xml   XML Output stream
  \param alspec   ArrayOfLineshapeSpec
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml _U_,
                     const ArrayOfLineshapeSpec& alspec _U_,
                     bofstream * pbofs _U_,
                     const String&  name _U_)

{
  // FIXME OLE: Implement this.
  throw runtime_error ("Boo. Not yet implemented.");
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
                     const String& name)
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



//=== ArrayOfTensor4=========================================================

//! Reads ArrayOfTensor4 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param atensor4  ArrayOfTensor4 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfTensor4& atensor4,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Tensor4");

  tag.get_attribute_value ("nelem", nelem);
  atensor4.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, atensor4[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfTensor4: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is_xml);
  tag.check_name ("/Array");
}


//! Writes ArrayOfTensor4 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param atensor4  ArrayOfTensor4
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfTensor4& atensor4,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "Tensor4");
  open_tag.add_attribute ("nelem", atensor4.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < atensor4.nelem (); n++)
    {
      xml_write_to_stream (os_xml, atensor4[n], pbofs);
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
                     const String& name)
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


//=== ArrayOfTensor7=========================================================

//! Reads ArrayOfTensor7 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param atensor7  ArrayOfTensor6 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void
xml_read_from_stream (istream& is_xml,
                      ArrayOfTensor7& atensor7,
                      bifstream *pbifs)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Tensor7");

  tag.get_attribute_value ("nelem", nelem);
  atensor7.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is_xml, atensor7[n], pbifs);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfTensor7: "
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
  \param atensor7  ArrayOfTensor7
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void
xml_write_to_stream (ostream& os_xml,
                     const ArrayOfTensor7& atensor7,
                     bofstream *pbofs,
                     const String& name)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");
  if (name.length ())
    open_tag.add_attribute ("name", name);

  open_tag.add_attribute ("type", "Tensor7");
  open_tag.add_attribute ("nelem", atensor7.nelem ());

  open_tag.write_to_stream (os_xml);
  os_xml << '\n';

  for (Index n = 0; n < atensor7.nelem (); n++)
    {
      xml_write_to_stream (os_xml, atensor7[n], pbofs);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os_xml);

  os_xml << '\n';
}


//=== ArrayOfString ==========================================================

//! Parse ArrayOfString from XML input stream
/*!
  \param is_xml   XML Input stream
  \param astring  ArrayOfString return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
  \param tag      XML tag object
*/
void
xml_parse_from_stream (istream& is_xml,
                       ArrayOfString& astring,
                       bifstream *pbifs,
                       ArtsXMLTag& tag)
{
  Index nelem;

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
}


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

  tag.read_from_stream (is_xml);
  tag.check_name ("Array");

  xml_parse_from_stream (is_xml, astring, pbifs, tag);

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
                     const String& name)
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
                     const String& name)
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

