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
#include "debug.h"
#include "isotopologues.h"
#include "sorted_grid.h"
#include "species_tags.h"
#include <workspace.h>
#include "xml_io.h"
#include "xml_io_array_macro.h"
#include "xml_io_arts_types.h"
#include <stdexcept>
#include <type_traits>

#include <rtepack.h>

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
void xml_read_from_stream(std::istream& is_xml,
                          ArrayOfGridPos& agpos,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "GridPos");

  tag.get_attribute_value("nelem", nelem);
  agpos.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, agpos[n], pbifs);
  } catch (const std::runtime_error& e) {
    std::ostringstream os;
    os << "Error reading ArrayOfGridPos: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw std::runtime_error(os.str());
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
void xml_write_to_stream(std::ostream& os_xml,
                         const ArrayOfGridPos& agpos,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "GridPos");
  open_tag.add_attribute("nelem", static_cast<Index>(agpos.size()));

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (const auto & agpo : agpos)
    xml_write_to_stream(os_xml, agpo, pbofs, "");

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
void xml_read_from_stream(std::istream& is_xml,
                          ArrayOfSpeciesTag& astag,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "SpeciesTag");

  tag.get_attribute_value("nelem", nelem);
  astag.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, astag[n], pbifs);
  } catch (const std::runtime_error& e) {
    std::ostringstream os;
    os << "Error reading ArrayOfSpeciesTag: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw std::runtime_error(os.str());
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
void xml_write_to_stream(std::ostream& os_xml,
                         const ArrayOfSpeciesTag& astag,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "SpeciesTag");
  open_tag.add_attribute("nelem", static_cast<Index>(astag.size()));

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (auto n : astag)
    xml_write_to_stream(os_xml, n, pbofs, "");

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfScatteringSpecies =========================================================

//! Reads ArrayOfScatteringSpecis from XML input stream
/*!
  \param is_xml  XML Input stream
  \param astar   astar return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(std::istream& is_xml,
                          ArrayOfScatteringSpecies& species,
                          bifstream* pbifs) {
  // Placeholder for now.
}

//! Writes ArrayOfScatteringSpecies to XML output stream
/*!
  \param os_xml  XML Output stream
  \param species   ArrayOfScatteringSpecies
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(std::ostream& os_xml,
                         const ArrayOfScatteringSpecies& species,
                         bofstream* pbofs,
                         const String& name) {
}


//=== ArrayOfSun =========================================================

//! Reads ArrayOfSun from XML input stream
/*!
  \param is_xml  XML Input stream
  \param astar   astar return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(std::istream& is_xml,
                          ArrayOfSun& astar,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Sun");

  tag.get_attribute_value("nelem", nelem);
  astar.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++) {
      xml_read_from_stream(is_xml, astar[n], pbifs);
    }
  } catch (const std::runtime_error& e) {
    std::ostringstream os;
    os << "Error reading ArrayOfSun: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw std::runtime_error(os.str());
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
void xml_write_to_stream(std::ostream& os_xml,
                         const ArrayOfSun& astar,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Sun");
  open_tag.add_attribute("nelem", static_cast<Index>(astar.size()));

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (const auto & n : astar) {
    xml_write_to_stream(os_xml, n, pbofs, "");
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
void xml_parse_from_stream(std::istream& is_xml,
                           ArrayOfString& astring,
                           bifstream* pbifs,
                           XMLTag& tag) {
  Index nelem;

  tag.check_attribute("type", "String");

  tag.get_attribute_value("nelem", nelem);
  astring.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, astring[n], pbifs);
  } catch (const std::runtime_error& e) {
    std::ostringstream os;
    os << "Error reading ArrayOfString: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw std::runtime_error(os.str());
  }
}

//=== ArrayOfXsecRecord ======================================================

//! Reads ArrayOfXsecData from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param axd        ArrayOfXsecData return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          ArrayOfXsecRecord& axd,
                          bifstream* pbifs) {
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "XsecRecord");

  tag.get_attribute_value("nelem", nelem);
  axd.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, axd[n], pbifs);
  } catch (const std::runtime_error& e) {
    std::ostringstream os;
    os << "Error reading ArrayOfXsecRecord: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw std::runtime_error(os.str());
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
void xml_write_to_stream(std::ostream& os_xml,
                         const ArrayOfXsecRecord& axd,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "XsecRecord");
  open_tag.add_attribute("nelem", static_cast<Index>(axd.size()));

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (const auto & n : axd)
    xml_write_to_stream(os_xml, n, pbofs, "");

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfAbsorptionLines)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfAgenda)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfAbsorptionLines)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfGriddedField1)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfGriddedField2)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfGriddedField3)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfIndex)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfMatrix)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfScatteringMetaData)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfSingleScatteringData)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfSpeciesTag)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfString)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfTensor3)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfTensor6)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfTime)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfAtmPoint)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfCIARecord)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfGriddedField1)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfGriddedField2)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfGriddedField3)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfGriddedField4)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfIndex)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfMatrix)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfQuantumIdentifier)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfScatteringMetaData)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfSingleScatteringData)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfSparse)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfString)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfTelsemAtlas)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfTensor3)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfTensor4)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfTensor5)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfTensor6)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfTensor7)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfTime)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfSpeciesEnum)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfPropagationPathPoint)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfPropagationPathPoint)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfArrayOfPropagationPathPoint)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfNamedGriddedField2)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfGriddedField1Named)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfAscendingGrid)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfVector2)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfVector3)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfAbsorptionBand)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfSpeciesIsotope)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfSensorObsel)
