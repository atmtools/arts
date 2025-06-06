////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io_basic_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains basic functions to handle XML data files.

*/

#include <algorithm>
#include <utility>
#include <vector>

#include "debug.h"
#include "double_imanip.h"
#include "enums.h"
#include "isotopologues.h"
#include "quantum_numbers.h"
#include "xml_io.h"
#include "xml_io_general_types.h"

////////////////////////////////////////////////////////////////////////////
//   Overloaded functions for reading/writing data from/to XML stream
////////////////////////////////////////////////////////////////////////////

//=== Rational =========================================================

//! Reads Rational from XML input stream
/*!
 * \param is_xml   XML Input stream
 * \param rational  Rational return value
 * \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          Rational& rational,
                          bifstream* pbifs) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("Rational");

  if (pbifs) {
    *pbifs >> rational;
    if (pbifs->fail()) {
      xml_data_parse_error(tag, "");
    }
  } else {
    is_xml >> rational;
    if (is_xml.fail()) {
      xml_data_parse_error(tag, "");
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Rational");
}

//! Writes Rational to XML output stream
/*!
 * \param os_xml   XML Output stream
 * \param rational Rational value
 * \param pbofs    Pointer to binary file stream. NULL for ASCII output.
 * \param name     Optional name attribute
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const Rational& rational,
                         bofstream* pbofs,
                         const String& name) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Rational");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);

  if (pbofs)
    *pbofs << rational;
  else
    std::print(os_xml, " {:IO} ", rational);

  close_tag.set_name("/Rational");
  close_tag.write_to_stream(os_xml);
  std::println(os_xml);
}

//=== Time ================================================================

//! Reads Time from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param t       Time return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          Time& t,
                          bifstream* pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("Time");

  Index version;
  tag.get_attribute_value("version", version);
  ARTS_USER_ERROR_IF(version not_eq 1,
                     "Your version of ARTS can only handle version 1 of Time");

  is_xml >> t;
  if (is_xml.fail()) {
    xml_data_parse_error(tag, "Time is poorly formatted");
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Time");
}

//! Writes Time to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param t       Time
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const Time& t,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Time");
  open_tag.add_attribute("version", t.Version());
  open_tag.write_to_stream(os_xml);

  xml_set_stream_precision(os_xml);

  std::print(os_xml, " {:IO} ", t);

  close_tag.set_name("/Time");
  close_tag.write_to_stream(os_xml);
  std::println(os_xml);
}

//=== SurfacePropertyTag ================================================================

//! Reads SurfacePropertyTag from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param vib     SurfacePropertyTag return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          SurfacePropertyTag& s,
                          bifstream* pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("SurfacePropertyTag");
  xml_read_from_stream(is_xml, s.name, pbifs);
  tag.read_from_stream(is_xml);
  tag.check_name("/SurfacePropertyTag");
}

//! Writes SurfacePropertyTag to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param vib     SurfacePropertyTag
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const SurfacePropertyTag& s,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("SurfacePropertyTag");
  open_tag.write_to_stream(os_xml);
  std::println(os_xml);

  xml_write_to_stream(os_xml, s.name, pbofs, "");

  close_tag.set_name("/SurfacePropertyTag");
  close_tag.write_to_stream(os_xml);
  std::println(os_xml);
}

//=== ScatteringSpeciesProperty ================================================================

//! Reads ScatteringSpeciesProperty from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param vib     ScatteringSpeciesProperty return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& /*is_xml*/,
                          ScatteringSpeciesProperty& /*s*/,
                          bifstream* pbifs [[maybe_unused]]) {}

//! Writes ScatteringSpeciesProperty to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param vib     ScatteringSpeciesProperty
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& /*os_xml*/,
                         const ScatteringSpeciesProperty& /*s*/,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {}

////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////
