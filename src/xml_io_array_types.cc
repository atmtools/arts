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
  \file   xml_io_array_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains basic functions to handle XML data files.

*/

#include "absorptionlines.h"
#include "arts.h"
#include "debug.h"
#include "propagationmatrix.h"
#include "species_tags.h"
#include "tokval.h"
#include "transmissionmatrix.h"
#include "workspace_ng.h"
#include "xml_io.h"
#include "xml_io_arts_types.h"
#include <stdexcept>
#include <type_traits>

////////////////////////////////////////////////////////////////////////////
//   Overloaded functions for reading/writing data from/to XML stream
////////////////////////////////////////////////////////////////////////////

//=== ArrayOfGridPos =========================================================

//! Reads ArrayOfGridPos from XML input stream
/*!
  \param is_xml  XML Input stream
  \param agpos   ArrayOfGridPos return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfGridPos& agpos,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "GridPos");

  tag.get_attribute_value("nelem", nelem);
  agpos.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, agpos[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfGridPos: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfGridPos to XML output stream
/*!
  \param os_xml  XML Output stream
  \param agpos   ArrayOfGridPos
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfGridPos& agpos,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "GridPos");
  open_tag.add_attribute("nelem", agpos.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < agpos.nelem(); n++)
    xml_write_to_stream(os_xml, agpos[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfRetrievalQuantity =======================================

//! Reads ArrayOfRetrievalQuantity from XML input stream
/*!
  \param is_xml    XML Input stream
  \param arq       ArrayOfRetrievalQuantity return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfRetrievalQuantity& arq,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "RetrievalQuantity");

  tag.get_attribute_value("nelem", nelem);
  arq.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, arq[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfRetrievalQuantity: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfRetrivalQuantity to XML output stream
/*!
  \param os_xml    XML Output stream
  \param arq       ArrayOfRetrievalQuantity
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfRetrievalQuantity& arq,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "RetrievalQuantity");
  open_tag.add_attribute("nelem", arq.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < arq.nelem(); n++)
    xml_write_to_stream(os_xml, arq[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfSpeciesTag ================================================

//! Reads ArrayOfSpeciesTag from XML input stream
/*!
  \param is_xml  XML Input stream
  \param astag   ArrayOfSpeciesTag return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfSpeciesTag& astag,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "SpeciesTag");

  tag.get_attribute_value("nelem", nelem);
  astag.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, astag[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfSpeciesTag: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfSpeciesTag to XML output stream
/*!
  \param os_xml  XML Output stream
  \param astag   ArrayOfSpeciesTag
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfSpeciesTag& astag,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "SpeciesTag");
  open_tag.add_attribute("nelem", astag.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < astag.nelem(); n++)
    xml_write_to_stream(os_xml, astag[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}


//=== ArrayOfSun =========================================================

//! Reads ArrayOfSun from XML input stream
/*!
  \param is_xml  XML Input stream
  \param astar   astar return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfSun& astar,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Sun");

  tag.get_attribute_value("nelem", nelem);
  astar.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++) {
      xml_read_from_stream(is_xml, astar[n], pbifs, verbosity);
    }
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfSun: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfSun to XML output stream
/*!
  \param os_xml  XML Output stream
  \param astar   ArrayOfSun
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfSun& astar,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Sun");
  open_tag.add_attribute("nelem", astar.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < astar.nelem(); n++) {
    xml_write_to_stream(os_xml, astar[n], pbofs, "", verbosity);
  }

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

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
void xml_parse_from_stream(istream& is_xml,
                           ArrayOfString& astring,
                           bifstream* pbifs,
                           XMLTag& tag,
                           const Verbosity& verbosity) {
  Index nelem;

  tag.check_attribute("type", "String");

  tag.get_attribute_value("nelem", nelem);
  astring.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, astring[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfString: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }
}

//=== ArrayOfXsecRecord ======================================================

//! Reads ArrayOfXsecData from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param axd        ArrayOfXsecData return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfXsecRecord& axd,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "XsecRecord");

  tag.get_attribute_value("nelem", nelem);
  axd.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, axd[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfXsecRecord: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfXsecData to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param axd     ArrayOfXsecData
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfXsecRecord& axd,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "XsecRecord");
  open_tag.add_attribute("nelem", axd.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < axd.nelem(); n++)
    xml_write_to_stream(os_xml, axd[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//! Both T and T{}[0] are ARTS groups exposed to the user if this is true
template <typename T>
concept array_of_group = ArtsType<T> and ArtsType<decltype(T{}[0])>;

template <array_of_group T>
void xml_read(istream &is_xml, T &at, bifstream *pbifs,
              const Verbosity &verbosity) try {
  const static String subtype =
      WorkspaceGroupNameValue<std::remove_cvref_t<decltype(T{}[0])>>;

  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", subtype);

  tag.get_attribute_value("nelem", nelem);
  at.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, at[n], pbifs, verbosity);
  } catch (const std::runtime_error &e) {
    ostringstream os;
    os << "Error reading "
       << WorkspaceGroupNameValue<std::remove_cvref_t<T>> << ": "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
} catch (std::runtime_error &e) {
  throw std::runtime_error(
      var_string("Failed reading routine for ",
                 WorkspaceGroupNameValue<std::remove_cvref_t<T>>,
                 "\nError reads:\n", e.what()));
}

template <array_of_group T>
void xml_write(ostream &os_xml, const T &at, bofstream *pbofs,
               const String &name, const Verbosity &verbosity) try {
  const static String subtype =
      WorkspaceGroupNameValue<std::remove_cvref_t<decltype(T{}[0])>>;

  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length())
    open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", subtype);
  open_tag.add_attribute("nelem", at.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < at.nelem(); n++)
    xml_write_to_stream(os_xml, at[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
} catch (std::runtime_error &e) {
  throw std::runtime_error(
      var_string("Failed saving routine for ",
                 WorkspaceGroupNameValue<std::remove_cvref_t<T>>,
                 "\nError reads:\n", e.what()));
}

//! Helper macro for when both Array<T> and T are ARTS groups
#define TMPL_XML_READ_WRITE_STREAM(T)                                          \
  void xml_read_from_stream(istream &is_xml, T &at, bifstream *pbifs,          \
                            const Verbosity &verbosity) {                      \
    xml_read(is_xml, at, pbifs, verbosity);                                    \
  }                                                                            \
                                                                               \
  void xml_write_to_stream(ostream &os_xml, const T &at, bofstream *pbofs,     \
                           const String &name, const Verbosity &verbosity) {   \
    xml_write(os_xml, at, pbofs, name, verbosity);                             \
  }

TMPL_XML_READ_WRITE_STREAM(ArrayOfAbsorptionLines)
TMPL_XML_READ_WRITE_STREAM(ArrayOfAgenda)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfAbsorptionLines)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfArrayOfTransmissionMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfGriddedField1)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfGriddedField2)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfGriddedField3)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfIndex)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfPropagationMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfRadiationVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfScatteringMetaData)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfSingleScatteringData)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfSpeciesTag)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfStokesVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfString)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfTensor3)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfTensor6)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfTime)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfTransmissionMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfArrayOfVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfAtmPoint)
TMPL_XML_READ_WRITE_STREAM(ArrayOfCIARecord)
TMPL_XML_READ_WRITE_STREAM(ArrayOfGriddedField1)
TMPL_XML_READ_WRITE_STREAM(ArrayOfGriddedField2)
TMPL_XML_READ_WRITE_STREAM(ArrayOfGriddedField3)
TMPL_XML_READ_WRITE_STREAM(ArrayOfGriddedField4)
TMPL_XML_READ_WRITE_STREAM(ArrayOfIndex)
TMPL_XML_READ_WRITE_STREAM(ArrayOfJacobianTarget)
TMPL_XML_READ_WRITE_STREAM(ArrayOfMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfPpath)
TMPL_XML_READ_WRITE_STREAM(ArrayOfPropagationMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfQuantumIdentifier)
TMPL_XML_READ_WRITE_STREAM(ArrayOfRadiationVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfScatteringMetaData)
TMPL_XML_READ_WRITE_STREAM(ArrayOfSingleScatteringData)
TMPL_XML_READ_WRITE_STREAM(ArrayOfSparse)
TMPL_XML_READ_WRITE_STREAM(ArrayOfStokesVector)
TMPL_XML_READ_WRITE_STREAM(ArrayOfString)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTelsemAtlas)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTensor3)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTensor4)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTensor5)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTensor6)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTensor7)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTime)
TMPL_XML_READ_WRITE_STREAM(ArrayOfTransmissionMatrix)
TMPL_XML_READ_WRITE_STREAM(ArrayOfVector)
