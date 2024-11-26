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

#include "absorptionlines.h"
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
    os_xml << ' ' << rational << ' ';

  close_tag.set_name("/Rational");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
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

  os_xml << ' ' << t << ' ';

  close_tag.set_name("/Time");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== AbsorptionLines ================================================================

//! Reads AbsorptionLines from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param al      AbsorptionLines return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          AbsorptionLines& al,
                          bifstream* pbifs) {
  static_assert(
      AbsorptionLines::version == 2,
      "The reading routine expects version 1 of the absorption lines data type to work");

  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("AbsorptionLines");

  Index version;
  if (tag.has_attribute("version")) {
    tag.get_attribute_value("version", version);
  } else {
    version = 0;
  }

  ARTS_USER_ERROR_IF(
      AbsorptionLines::version < version,
      "The version of this catalog is too new.  You need to upgrade ARTS to use it.")

  ARTS_USER_ERROR_IF(
      version < AbsorptionLines::version - 1 or
          (pbifs and version not_eq AbsorptionLines::version),
      "Using descoped version of the catalog; version: {}"
      "\n"
      "We only ever support limited number of versions. Your compilation supports versions {}"
      " and {}"
      " in ascii but only version {}"
      " in binary\n\n"
      "To update from versions, please check-out the following branch, compile, and save your catalog again (as ascii):\n"
      "0 to 1: 3b6565fb93702308c4cdd660ec63c71d63dcaf26\n"
      "1 to 2: Current version\n",
      version,
      AbsorptionLines::version,
      AbsorptionLines::version - 1,
      AbsorptionLines::version)

  // Number of lines
  Index nlines;
  tag.get_attribute_value("nlines", nlines);

  // Identity of the lines (Changes between versions)
  QuantumIdentifier id;
  if (version == 2) {
    String id_str;
    tag.get_attribute_value("id", id_str);
    id = QuantumIdentifier(id_str);
  } else if (version == 1) {
    String spec;
    tag.get_attribute_value("species", spec);

    Index spec_ind = Species::find_species_index(spec);
    ARTS_USER_ERROR_IF(spec_ind < 0, "Bad species index for: {}", spec)
    id.isotopologue_index = spec_ind;
  }

  // Cutoff type
  String s_cutoff;
  tag.get_attribute_value("cutofftype", s_cutoff);
  const AbsorptionCutoffTypeOld cutoff = to<AbsorptionCutoffTypeOld>(s_cutoff);

  // Mirroring type
  String s_mirroring;
  tag.get_attribute_value("mirroringtype", s_mirroring);
  const AbsorptionMirroringTypeOld mirroring =
      to<AbsorptionMirroringTypeOld>(s_mirroring);

  // Line population type
  String s_population;
  tag.get_attribute_value("populationtype", s_population);
  const AbsorptionPopulationTypeOld population =
      to<AbsorptionPopulationTypeOld>(s_population);

  // Normalization type
  String s_normalization;
  tag.get_attribute_value("normalizationtype", s_normalization);
  const AbsorptionNormalizationTypeOld normalization =
      to<AbsorptionNormalizationTypeOld>(s_normalization);

  // Shape type
  String s_lineshapetype;
  tag.get_attribute_value("lineshapetype", s_lineshapetype);
  const LineShapeTypeOld lineshapetype = to<LineShapeTypeOld>(s_lineshapetype);

  /** Reference temperature for all parameters of the lines */
  Numeric T0;
  tag.get_attribute_value("T0", T0);

  /** cutoff frequency */
  Numeric cutofffreq;
  tag.get_attribute_value("cutofffreq", cutofffreq);

  /** linemixing limit */
  Numeric linemixinglimit;
  tag.get_attribute_value("linemixinglimit", linemixinglimit);

  /** List of local quantum numbers, these must be defined and not strings */
  Quantum::Number::LocalState meta_localstate;
  String localquanta_str;
  tag.get_attribute_value("localquanta", localquanta_str);
  const Index nlocal = Quantum::Number::count_items(localquanta_str);
  Array<QuantumNumberType> qn_key;
  for (Index i = 0; i < nlocal; i++)
    qn_key.push_back(
        to<QuantumNumberType>(Quantum::Number::items(localquanta_str, i)));
  ARTS_USER_ERROR_IF(
      std::any_of(
          qn_key.begin(),
          qn_key.end(),
          [](auto& qn) {
            return Quantum::Number::common_value_type(
                       Quantum::Number::common_value_type(qn),
                       Quantum::Number::QuantumNumberValueType::H) not_eq
                   Quantum::Number::QuantumNumberValueType::H;
          }),
      "Quantum number list contains a string type, this is not allowed: [{}]",
      qn_key)
  meta_localstate.set_unsorted_qns(qn_key);

  /** Catalog ID */
  if (version == 1) {
    String uid, lid;
    tag.get_attribute_value("upperglobalquanta", uid);
    tag.get_attribute_value("lowerglobalquanta", lid);
    id.val = Quantum::Number::ValueList(uid, lid);
  }

  /** A list of broadening species */
  ArrayOfSpeciesEnum broadeningspecies;
  bool selfbroadening;
  bool bathbroadening;
  tag.get_attribute_value(
      "broadeningspecies", broadeningspecies, selfbroadening, bathbroadening);
  if (selfbroadening) broadeningspecies.front() = id.Species();

  String temperaturemodes;
  tag.get_attribute_value("temperaturemodes", temperaturemodes);
  auto metamodel = LineShape::MetaData2ModelShape(temperaturemodes);

  al = AbsorptionLines(selfbroadening,
                       bathbroadening,
                       nlines,
                       cutoff,
                       mirroring,
                       population,
                       normalization,
                       lineshapetype,
                       T0,
                       cutofffreq,
                       linemixinglimit,
                       id,
                       broadeningspecies,
                       meta_localstate,
                       metamodel);

  if (pbifs) {
    al.read(*pbifs);
    if (pbifs->fail()) {
      std::ostringstream os;
      os << "AbsorptionLines has wrong dimensions";
      xml_data_parse_error(tag, os.str());
    }
  } else {
    is_xml >> al;
    if (is_xml.fail()) {
      std::ostringstream os;
      os << "AbsorptionLines has wrong dimensions";
      xml_data_parse_error(tag, os.str());
    }
  }

  // Finalize the sorting because we have to
  for (auto& line : al.lines) line.localquanta.val.finalize();

  tag.read_from_stream(is_xml);
  tag.check_name("/AbsorptionLines");
}

//! Writes AbsorptionLines to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param al      AbsorptionLines
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const AbsorptionLines& al,
                         bofstream* pbofs,
                         const String&) {
  ArtsXMLTag open_comment_tag;
  ArtsXMLTag close_comment_tag;
  open_comment_tag.set_name("comment");
  open_comment_tag.write_to_stream(os_xml);
  os_xml << al.MetaData();
  close_comment_tag.set_name("/comment");
  close_comment_tag.write_to_stream(os_xml);
  os_xml << '\n';

  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("AbsorptionLines");
  open_tag.add_attribute("version", al.version);
  open_tag.add_attribute("id", std::format("{}", al.quantumidentity));
  open_tag.add_attribute("nlines", al.NumLines());
  open_tag.add_attribute("cutofftype", String{toString(al.cutoff)});
  open_tag.add_attribute("mirroringtype", String{toString(al.mirroring)});
  open_tag.add_attribute("populationtype", String{toString(al.population)});
  open_tag.add_attribute("normalizationtype",
                         String{toString(al.normalization)});
  open_tag.add_attribute("lineshapetype", String{toString(al.lineshapetype)});
  open_tag.add_attribute("T0", al.T0);
  open_tag.add_attribute("cutofffreq", al.cutofffreq);
  open_tag.add_attribute("linemixinglimit", al.linemixinglimit);

  const String localquanta_str =
      al.NumLines() ? al.lines.front().localquanta.keys() : "";
  open_tag.add_attribute("localquanta", localquanta_str);

  open_tag.add_attribute("broadeningspecies",
                         al.broadeningspecies,
                         al.selfbroadening,
                         al.bathbroadening);
  open_tag.add_attribute("temperaturemodes", al.LineShapeMetaData());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  if (pbofs)
    al.write(*pbofs);
  else
    os_xml << al;

  close_tag.set_name("/AbsorptionLines");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== VibrationalEnergyLevels ================================================================

//! Reads VibrationalEnergyLevels from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param vib     VibrationalEnergyLevels return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          VibrationalEnergyLevels& vib,
                          bifstream* pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("VibrationalEnergyLevels");
  const Index n = [&]() {
    Index x;
    tag.get_attribute_value("nelem", x);
    return x;
  }();
  vib.data.clear();
  vib.data.reserve(n);

  for (Index i = 0; i < n; i++) {
    tag.read_from_stream(is_xml);
    tag.check_name("Data");
    QuantumIdentifier qkey = [&]() {
      String k;
      tag.get_attribute_value("key", k);
      return QuantumIdentifier{k};
    }();
    Numeric data;
    if (pbifs)
      *pbifs >> data;
    else
      is_xml >> data;
    vib[std::move(qkey)] = data;
    tag.read_from_stream(is_xml);
    tag.check_name("/Data");
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/VibrationalEnergyLevels");
}

//! Writes VibrationalEnergyLevels to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param vib     VibrationalEnergyLevels
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const VibrationalEnergyLevels& vib,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("VibrationalEnergyLevels");
  open_tag.add_attribute("nelem", vib.nelem());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (auto& a : vib) {
    ArtsXMLTag open_data_tag;
    open_data_tag.set_name("Data");
    open_data_tag.add_attribute("key", std::format("{}", a.first));
    open_data_tag.write_to_stream(os_xml);
    if (pbofs)
      *pbofs << a.second;
    else
      os_xml << ' ' << a.second << ' ';
    ArtsXMLTag close_data_tag;
    close_data_tag.set_name("/Data");
    close_data_tag.write_to_stream(os_xml);
    os_xml << '\n';
  }

  close_tag.set_name("/VibrationalEnergyLevels");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== SurfaceTypeTag ================================================================

//! Reads SurfaceTypeTag from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param vib     SurfaceTypeTag return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream& is_xml,
                          SurfaceTypeTag& s,
                          bifstream* pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("SurfaceTypeTag");
  xml_read_from_stream(is_xml, s.name, pbifs);
  tag.read_from_stream(is_xml);
  tag.check_name("/SurfaceTypeTag");
}

//! Writes SurfaceTypeTag to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param vib     SurfaceTypeTag
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const SurfaceTypeTag& s,
                         bofstream* pbofs [[maybe_unused]],
                         const String&) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("SurfaceTypeTag");
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_write_to_stream(os_xml, s.name, pbofs, "");

  close_tag.set_name("/SurfaceTypeTag");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
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
  os_xml << '\n';

  xml_write_to_stream(os_xml, s.name, pbofs, "");

  close_tag.set_name("/SurfacePropertyTag");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
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
                          bifstream* pbifs [[maybe_unused]]) {
}

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
                         const String&) {
}

////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////
