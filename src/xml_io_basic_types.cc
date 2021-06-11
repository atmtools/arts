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
  \file   xml_io_basic_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains basic functions to handle XML data files.

*/

#include "arts.h"
#include "xml_io.h"

////////////////////////////////////////////////////////////////////////////
//   Overloaded functions for reading/writing data from/to XML stream
////////////////////////////////////////////////////////////////////////////

//=== JacobianTarget ==================================================================

//! Reads JacobianTarget from XML input stream
/*!
  \param is_xml  XML Input stream
  \param jt      JacobianTarget return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          JacobianTarget& jt,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("JacobianTarget");
  
  // Type information
  String typestr, subtypestr;
  tag.get_attribute_value("Type", typestr);
  tag.get_attribute_value("SubType", subtypestr);
  jt.TargetType(typestr);
  jt.TargetSubType(subtypestr);
  
  /** Catalog ID */
  if (jt.needQuantumIdentity()) {
    SpeciesTag spec;
    tag.get_attribute_value("species", spec);
    
    if (jt == Jacobian::Line::VMR) {
      jt.QuantumIdentity() = QuantumIdentifier(spec.Isotopologue(), Quantum::IdentifierType::All);
    } else if (jt == Jacobian::Line::NLTE) {
      QuantumNumbers levelquanta;
      tag.get_attribute_value("levelquanta", levelquanta);
      jt.QuantumIdentity() = QuantumIdentifier(spec.Isotopologue(), levelquanta);
    } else {
      QuantumNumbers upperglobalquanta, lowerglobalquanta;
      tag.get_attribute_value("upperglobalquanta", upperglobalquanta);
      tag.get_attribute_value("lowerglobalquanta", lowerglobalquanta);
      jt.QuantumIdentity() = QuantumIdentifier(spec.Isotopologue(), upperglobalquanta, lowerglobalquanta);
    }
  }

  if (jt.needArrayOfSpeciesTag()) {
    String key;
    tag.get_attribute_value("species", key);
    jt.SpeciesList() = ArrayOfSpeciesTag(key);
  }

  if (jt.needString()) {
    tag.get_attribute_value("string_key", jt.StringKey());
  }
  
  if (pbifs) {
    *pbifs >> jt.Perturbation();
    if (pbifs->fail()) {
      xml_data_parse_error(tag, "");
    }
  } else {
    is_xml >> jt.Perturbation();
    if (is_xml.fail()) {
      xml_data_parse_error(tag, "");
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/JacobianTarget");
  
  ARTS_USER_ERROR_IF (not jt.TargetSubTypeOK(),
    "Bad input: ", typestr, " or ", subtypestr, '\n', "\tCannot be interpreted as a type or substype...\n")
}

//! Writes JacobianTarget to XML output stream
/*!
  \param os_xml  XML Output stream
  \param jt      JacobianTarget value
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const JacobianTarget& jt,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  os_xml << '\n';
  open_tag.set_name("JacobianTarget");
  if (name.length()) open_tag.add_attribute("name", name);
  
  // Type information
  open_tag.add_attribute("Type", jt.TargetType());
  open_tag.add_attribute("SubType", jt.TargetSubType());
  
  /** Catalog ID */
  if (jt.needQuantumIdentity()) {
    open_tag.add_attribute("species", String(Species::toShortName(jt.QuantumIdentity().Species())));
    
    if (jt == Jacobian::Line::VMR) {
    } else if (jt == Jacobian::Line::NLTE) {
      open_tag.add_attribute("levelquanta", jt.QuantumIdentity().Level().toString());
    } else {
      open_tag.add_attribute("upperglobalquanta", jt.QuantumIdentity().Upper().toString());
      open_tag.add_attribute("lowerglobalquanta", jt.QuantumIdentity().Lower().toString());
    }
  }

  if (jt.needArrayOfSpeciesTag()) {
    open_tag.add_attribute("species", jt.SpeciesList().Name());
  }

  if (jt.needString()) {
    open_tag.add_attribute("string_key", jt.StringKey());
  }
  open_tag.write_to_stream(os_xml);

  if (pbofs)
    *pbofs << jt.Perturbation();
  else
    os_xml << ' ' << jt.Perturbation() << ' ';

  close_tag.set_name("/JacobianTarget");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== Rational =========================================================

//! Reads Rational from XML input stream
/*!
 * \param is_xml   XML Input stream
 * \param rational  Rational return value
 * \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          Rational& rational,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);

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
void xml_write_to_stream(ostream& os_xml,
                         const Rational& rational,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Rational");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.write_to_stream(os_xml);

  if (pbofs)
    *pbofs << rational;
  else
    os_xml << rational;

  close_tag.set_name("/Rational");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== TransmissionMatrix ================================================================

//! Reads TransmissionMatrix from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param tm      TransmissionMatrix return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          TransmissionMatrix& tm,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index stokes_dim, nf;

  tag.read_from_stream(is_xml);
  tag.check_name("TransmissionMatrix");

  tag.get_attribute_value("Stokes", stokes_dim);
  tag.get_attribute_value("Freqs", nf);
  tm = TransmissionMatrix(nf, stokes_dim);
  if (pbifs) {
    *pbifs >> tm;
    if (pbifs->fail()) {
      ostringstream os;
      os << "TransmissionMatrix has wrong dimensions";
      xml_data_parse_error(tag, os.str());
    }
  } else {
    is_xml >> tm;
    if (is_xml.fail()) {
      ostringstream os;
      os << "TransmissionMatrix has wrong dimensions";
      xml_data_parse_error(tag, os.str());
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/TransmissionMatrix");
}

//! Writes RadiationVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param tm      TransmissionMatrix
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const TransmissionMatrix& tm,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("TransmissionMatrix");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("Stokes", tm.StokesDim());
  open_tag.add_attribute("Freqs", tm.Frequencies());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  if (pbofs)
    *pbofs << tm;
  else
    os_xml << tm;

  close_tag.set_name("/TransmissionMatrix");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== RadiationVector ================================================================

//! Reads RadiationVector from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param rv      RadiationVector return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          RadiationVector& rv,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index stokes_dim, nf;

  tag.read_from_stream(is_xml);
  tag.check_name("RadiationVector");

  tag.get_attribute_value("Stokes", stokes_dim);
  tag.get_attribute_value("Freqs", nf);
  rv = RadiationVector(nf, stokes_dim);
  if (pbifs) {
    *pbifs >> rv;
    if (pbifs->fail()) {
      ostringstream os;
      os << "RadiationVector has wrong dimensions";
      xml_data_parse_error(tag, os.str());
    }
  } else {
    is_xml >> rv;
    if (is_xml.fail()) {
      ostringstream os;
      os << "RadiationVector has wrong dimensions";
      xml_data_parse_error(tag, os.str());
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/RadiationVector");
}

//! Writes RadiationVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param rv      RadiationVector
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const RadiationVector& rv,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("RadiationVector");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("Stokes", rv.StokesDim());
  open_tag.add_attribute("Freqs", rv.Frequencies());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  if (pbofs)
    *pbofs << rv;
  else
    os_xml << rv;

  close_tag.set_name("/RadiationVector");
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
void xml_read_from_stream(istream& is_xml,
                          Time& t,
                          bifstream* pbifs [[maybe_unused]],
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  
  tag.read_from_stream(is_xml);
  tag.check_name("Time");
  
  Index version;
  tag.get_attribute_value("version", version);
  ARTS_USER_ERROR_IF (version not_eq 1,
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
void xml_write_to_stream(ostream& os_xml,
                         const Time& t,
                         bofstream* pbofs [[maybe_unused]],
                         const String&,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

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
void xml_read_from_stream(istream& is_xml,
                          AbsorptionLines& al,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  static_assert(AbsorptionLines::version == 1, "The reading routine expects version 1 of the absorption lines data type to work");
  
  ArtsXMLTag tag(verbosity);
  
  tag.read_from_stream(is_xml);
  tag.check_name("AbsorptionLines");
  
  Index version;
  if (tag.has_attribute("version")) {
    tag.get_attribute_value("version", version);
  } else {
    version = 0;
  }
  
  // Number of lines
  Index nlines;
  tag.get_attribute_value("nlines", nlines);
  
  // Species of the lines
  SpeciesTag spec;
  if (version == 1) {
    tag.get_attribute_value("species", spec);
  } else {
    String specname;
    tag.get_attribute_value("species", specname);
    specname = Species::update_isot_name(specname);
    
    spec = SpeciesTag(specname);
  }
  
  // Cutoff type
  String s_cutoff;
  tag.get_attribute_value("cutofftype", s_cutoff);
  const Absorption::CutoffType cutoff = Absorption::toCutoffType(s_cutoff);
  
  // Mirroring type
  String s_mirroring;
  tag.get_attribute_value("mirroringtype", s_mirroring);
  const Absorption::MirroringType mirroring = Absorption::toMirroringType(s_mirroring);
  
  // Line population type
  String s_population;
  tag.get_attribute_value("populationtype", s_population);
  const Absorption::PopulationType population = Absorption::toPopulationType(s_population);
  
  // Normalization type
  String s_normalization;
  tag.get_attribute_value("normalizationtype", s_normalization);
  const Absorption::NormalizationType normalization = Absorption::toNormalizationType(s_normalization);
  
  // Shape type
  String s_lineshapetype;
  tag.get_attribute_value("lineshapetype", s_lineshapetype);
  const LineShape::Type lineshapetype = LineShape::toType(s_lineshapetype);
  
  /** Reference temperature for all parameters of the lines */
  Numeric T0;
  tag.get_attribute_value("T0", T0);
  
  /** cutoff frequency */
  Numeric cutofffreq;
  tag.get_attribute_value("cutofffreq", cutofffreq);
  
  /** linemixing limit */
  Numeric linemixinglimit;
  tag.get_attribute_value("linemixinglimit", linemixinglimit);
  
  /** List of local quantum numbers, these must be defined */
  std::vector<QuantumNumberType> localquanta;
  tag.get_attribute_value("localquanta", localquanta);
  
  /** Catalog ID */  
  QuantumNumbers upperglobalquanta;
  tag.get_attribute_value("upperglobalquanta", upperglobalquanta);
  QuantumNumbers lowerglobalquanta;
  tag.get_attribute_value("lowerglobalquanta", lowerglobalquanta);
  
  /** A list of broadening species */
  ArrayOfSpecies broadeningspecies;
  bool selfbroadening;
  bool bathbroadening;
  tag.get_attribute_value("broadeningspecies", broadeningspecies, selfbroadening, bathbroadening);
  if (selfbroadening) broadeningspecies.front() = spec.Spec();
  
  String temperaturemodes;
  tag.get_attribute_value("temperaturemodes", temperaturemodes);
  auto metamodel = LineShape::MetaData2ModelShape(temperaturemodes);

  al = AbsorptionLines(selfbroadening, bathbroadening,
                       nlines, cutoff, mirroring,
                       population, normalization,
                       lineshapetype, T0, cutofffreq,
                       linemixinglimit, QuantumIdentifier(
                         spec.Isotopologue(), upperglobalquanta, lowerglobalquanta),
                       localquanta, broadeningspecies, metamodel);
  
  if (pbifs) {
    al.read(*pbifs);
     if (pbifs->fail()) {
       ostringstream os;
       os << "AbsorptionLines has wrong dimensions";
       xml_data_parse_error(tag, os.str());
     }
  } else {
    is_xml >> al;
    if (is_xml.fail()) {
      ostringstream os;
      os << "AbsorptionLines has wrong dimensions";
      xml_data_parse_error(tag, os.str());
    }
  }

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
void xml_write_to_stream(ostream& os_xml,
                         const AbsorptionLines& al,
                         bofstream* pbofs,
                         const String&,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_comment_tag(verbosity);
  ArtsXMLTag close_comment_tag(verbosity);
  open_comment_tag.set_name("comment");
  open_comment_tag.write_to_stream(os_xml);
  os_xml << al.MetaData();
  close_comment_tag.set_name("/comment");
  close_comment_tag.write_to_stream(os_xml);
  os_xml << '\n';

  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("AbsorptionLines");
  open_tag.add_attribute("version", al.version);
  open_tag.add_attribute("nlines", al.NumLines());
  open_tag.add_attribute("species", al.SpeciesName());
  open_tag.add_attribute("cutofftype", Absorption::toString(al.Cutoff()));
  open_tag.add_attribute("mirroringtype", Absorption::toString(al.Mirroring()));
  open_tag.add_attribute("populationtype", Absorption::toString(al.Population()));
  open_tag.add_attribute("normalizationtype", Absorption::toString(al.Normalization()));
  open_tag.add_attribute("lineshapetype", LineShape::toString(al.LineShapeType()));
  open_tag.add_attribute("T0", al.T0());
  open_tag.add_attribute("cutofffreq", al.CutoffFreqValue());
  open_tag.add_attribute("linemixinglimit", al.LinemixingLimit());
  open_tag.add_attribute("localquanta", al.LocalQuanta());
  open_tag.add_attribute("upperglobalquanta", al.UpperQuantumNumbers());
  open_tag.add_attribute("lowerglobalquanta", al.LowerQuantumNumbers());
  open_tag.add_attribute("broadeningspecies", al.BroadeningSpecies(), al.Self(), al.Bath());
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

////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

// FIXME: These should be implemented, sooner or later...

void xml_read_from_stream(istream&,
                          Timer&,
                          bifstream* /* pbifs */,
                          const Verbosity&) {
  ARTS_USER_ERROR_IF(true, "Method not implemented!");
}

void xml_write_to_stream(ostream&,
                         const Timer&,
                         bofstream* /* pbofs */,
                         const String& /* name */,
                         const Verbosity&) {
  ARTS_USER_ERROR_IF(true, "Method not implemented!");
}
