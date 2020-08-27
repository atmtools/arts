/* Copyright (C) 2019
 Richard Larsson <larsson@mps.mpg.de>

 
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

/** Contains the absorption namespace
 * @file   m_absorptionlines.cc
 * @author Richard Larsson
 * @date   2019-09-11
 * 
 * @brief  Contains the user interaction with absorption lines
 **/

#include "absorptionlines.h"
#include "auto_md.h"
#include "file.h"
#include "global_data.h"
#include "xml_io_private.h"
#include "m_xml.h"

/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////// Reading old/external functions
/////////////////////////////////////////////////////////////////////////////////////

/** Get a list of quantum numbers from a string
 * 
 * @param[in] qnstr A string such as "J N v1"
 * 
 * @return List of quantum numbers
 */
std::vector<QuantumNumberType> string2vecqn(const String& qnstr)
{
  std::vector<QuantumNumberType> nums(0);
  
  String part;
  if (qnstr not_eq "") {
    std::istringstream str(qnstr);
    while (not str.eof()) {
      str >> part; 
      if (IsValidQuantumNumberName(part)) {
        nums.push_back(string2quantumnumbertype(part));
      } else {
        std::ostringstream os;
        os << "The quantum number key: \"" << part << "\" is invalid.\n";
        throw std::runtime_error(os.str());
      }
    }
  }
  
  return nums;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ReadArrayOfARTSCAT(ArrayOfAbsorptionLines& abs_lines,
                        const String& artscat_file,
                        const Numeric& fmin,
                        const Numeric& fmax,
                        const String& globalquantumnumbers,
                        const String& localquantumnumbers,
                        const String& normalization_option,
                        const String& mirroring_option,
                        const String& population_option,
                        const String& lineshapetype_option,
                        const String& cutoff_option,
                        const Numeric& cutoff_value,
                        const Numeric& linemixinglimit_value,
                        const Verbosity& verbosity)
{
  // Global numbers
  const std::vector<QuantumNumberType> global_nums = string2vecqn(globalquantumnumbers);
  
  // Local numbers
  const std::vector<QuantumNumberType> local_nums = string2vecqn(localquantumnumbers);
  
  CREATE_OUT2;
  
  ArtsXMLTag tag(verbosity);
  Index nelem;
  
  // ARTSCAT data
  shared_ptr<istream> ifs = nullptr;
  xml_find_and_open_input_file(ifs, artscat_file, verbosity);
  istream& is_xml = *ifs;
  auto a = FILE_TYPE_ASCII;
  auto b = NUMERIC_TYPE_DOUBLE;
  auto c = ENDIAN_TYPE_LITTLE;
  xml_read_header_from_stream(is_xml, a,b,c,
                              verbosity);
  
  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  
  Index num_arrays;
  tag.get_attribute_value("nelem", num_arrays);
  
  std::vector<Absorption::SingleLineExternal> v(0);
  
  for (Index i=0; i<num_arrays; i++) {
    tag.read_from_stream(is_xml);
    tag.check_name("ArrayOfLineRecord");
    
    tag.get_attribute_value("nelem", nelem);
    
    String version;
    tag.get_attribute_value("version", version);
    
    Index artscat_version;
    
    if (version == "3") {
      artscat_version = 3;
    } else if (version.substr(0, 8) != "ARTSCAT-") {
      ostringstream os;
      os << "The ARTS line file you are trying to read does not contain a valid version tag.\n"
      << "Probably it was created with an older version of ARTS that used different units.";
      throw runtime_error(os.str());
    } else {
      istringstream is(version.substr(8));
      is >> artscat_version;
    }
    
    if (artscat_version < 3 or artscat_version > 5) {
      ostringstream os;
      os << "Unknown ARTS line file version: " << version;
      throw runtime_error(os.str());
    }
    
    bool go_on = true;
    Index n = 0;
    while (go_on and n<nelem) {
      switch(artscat_version) {
        case 3:
          v.push_back(Absorption::ReadFromArtscat3Stream(is_xml));
          break;
        case 4:
          v.push_back(Absorption::ReadFromArtscat4Stream(is_xml));
          break;
        case 5:
          v.push_back(Absorption::ReadFromArtscat5Stream(is_xml));
          break;
        default:
          throw std::runtime_error("Bad version!");
      }
      
      if (v.back().bad) {
        v.pop_back();
        go_on = false;
      } else if (v.back().line.F0() < fmin) {
        v.pop_back();
      } else if (v.back().line.F0() > fmax) {
        v.pop_back();
        go_on = false;
      }
      
      n++;
    }
    
    tag.read_from_stream(is_xml);
    tag.check_name("/ArrayOfLineRecord");
  }
  
  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
  
  for (auto& x: v)
    x.line.Zeeman() = Zeeman::GetAdvancedModel(x.quantumidentity);
  
  auto x = Absorption::split_list_of_external_lines(v, local_nums, global_nums);
  abs_lines.resize(0);
  abs_lines.reserve(x.size());
  while (x.size()) {
    abs_lines.push_back(x.back());
    abs_lines.back().sort_by_frequency();
    x.pop_back();
  }
  
  abs_linesSetNormalization(abs_lines, normalization_option, verbosity);
  abs_linesSetMirroring(abs_lines, mirroring_option, verbosity);
  abs_linesSetPopulation(abs_lines, population_option, verbosity);
  abs_linesSetLineShapeType(abs_lines, lineshapetype_option, verbosity);
  abs_linesSetCutoff(abs_lines, cutoff_option, cutoff_value, verbosity);
  abs_linesSetLinemixingLimit(abs_lines, linemixinglimit_value, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ReadARTSCAT(ArrayOfAbsorptionLines& abs_lines,
                 const String& artscat_file,
                 const Numeric& fmin,
                 const Numeric& fmax,
                 const String& globalquantumnumbers,
                 const String& localquantumnumbers,
                 const String& normalization_option,
                 const String& mirroring_option,
                 const String& population_option,
                 const String& lineshapetype_option,
                 const String& cutoff_option,
                 const Numeric& cutoff_value,
                 const Numeric& linemixinglimit_value,
                 const Verbosity& verbosity)
{
  // Global numbers
  const std::vector<QuantumNumberType> global_nums = string2vecqn(globalquantumnumbers);
  
  // Local numbers
  const std::vector<QuantumNumberType> local_nums = string2vecqn(localquantumnumbers);
  
  CREATE_OUT2;
  
  ArtsXMLTag tag(verbosity);
  Index nelem;
  
  // ARTSCAT data
  shared_ptr<istream> ifs = nullptr;
  xml_find_and_open_input_file(ifs, artscat_file, verbosity);
  istream& is_xml = *ifs;
  auto a = FILE_TYPE_ASCII;
  auto b = NUMERIC_TYPE_DOUBLE;
  auto c = ENDIAN_TYPE_LITTLE;
  xml_read_header_from_stream(is_xml, a,b,c,
                              verbosity);
  
  tag.read_from_stream(is_xml);
  tag.check_name("ArrayOfLineRecord");
  
  tag.get_attribute_value("nelem", nelem);
  
  String version;
  tag.get_attribute_value("version", version);
  
  Index artscat_version;
  
  if (version == "3") {
    artscat_version = 3;
  } else if (version.substr(0, 8) != "ARTSCAT-") {
    ostringstream os;
    os << "The ARTS line file you are trying to read does not contain a valid version tag.\n"
    << "Probably it was created with an older version of ARTS that used different units.";
    throw runtime_error(os.str());
  } else {
    istringstream is(version.substr(8));
    is >> artscat_version;
  }
  
  if (artscat_version < 3 or artscat_version > 5) {
    ostringstream os;
    os << "Unknown ARTS line file version: " << version;
    throw runtime_error(os.str());
  }
  
  std::vector<Absorption::SingleLineExternal> v(0);
  
  bool go_on = true;
  Index n = 0;
  while (n<nelem) {
    if (go_on) {
      switch(artscat_version) {
        case 3:
          v.push_back(Absorption::ReadFromArtscat3Stream(is_xml));
          break;
        case 4:
          v.push_back(Absorption::ReadFromArtscat4Stream(is_xml));
          break;
        case 5:
          v.push_back(Absorption::ReadFromArtscat5Stream(is_xml));
          break;
        default:
          throw std::runtime_error("Bad version!");
      }
      
      if (v.back().bad) {
        v.pop_back();
        go_on = false;
      } else if (v.back().line.F0() < fmin) {
        v.pop_back();
      } else if (v.back().line.F0() > fmax) {
        v.pop_back();
        go_on = false;
      }
    } else {
      String line;
      getline(is_xml, line);
    }
    
    n++;
  }
  
  tag.read_from_stream(is_xml);
  tag.check_name("/ArrayOfLineRecord");
  
  for (auto& x: v)
    x.line.Zeeman() = Zeeman::GetAdvancedModel(x.quantumidentity);
  
  auto x = Absorption::split_list_of_external_lines(v, local_nums, global_nums);
  abs_lines.resize(0);
  abs_lines.reserve(x.size());
  while (x.size()) {
    abs_lines.push_back(x.back());
    abs_lines.back().sort_by_frequency();
    x.pop_back();
  }
  
  abs_linesSetNormalization(abs_lines, normalization_option, verbosity);
  abs_linesSetMirroring(abs_lines, mirroring_option, verbosity);
  abs_linesSetPopulation(abs_lines, population_option, verbosity);
  abs_linesSetLineShapeType(abs_lines, lineshapetype_option, verbosity);
  abs_linesSetCutoff(abs_lines, cutoff_option, cutoff_value, verbosity);
  abs_linesSetLinemixingLimit(abs_lines, linemixinglimit_value, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ReadSplitARTSCAT(ArrayOfAbsorptionLines& abs_lines,
                      const ArrayOfArrayOfSpeciesTag& abs_species,
                      const String& basename,
                      const Numeric& fmin,
                      const Numeric& fmax,
                      const String& globalquantumnumbers,
                      const String& localquantumnumbers,
                      const Index& ignore_missing,
                      const String& normalization_option,
                      const String& mirroring_option,
                      const String& population_option,
                      const String& lineshapetype_option,
                      const String& cutoff_option,
                      const Numeric& cutoff_value,
                      const Numeric& linemixinglimit_value,
                      const Verbosity& verbosity)
{
  using global_data::species_data;
  
  // Build a set of species indices. Duplicates are ignored.
  std::set<Index> unique_species;
  for (auto asp = abs_species.begin(); asp != abs_species.end(); asp++) {
    for (ArrayOfSpeciesTag::const_iterator sp = asp->begin(); sp != asp->end(); sp++) {
      if (sp->Type() == SpeciesTag::TYPE_PLAIN || sp->Type() == SpeciesTag::TYPE_ZEEMAN) {
        unique_species.insert(sp->Species());
      }
    }
  }
  
  String tmpbasename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/') {
    tmpbasename += '.';
  }
  
  // Read catalogs for each identified species and put them all into
  // abs_lines.
  abs_lines.resize(0);
  for (auto it = unique_species.begin(); it != unique_species.end(); it++) {
    ArrayOfAbsorptionLines more_abs_lines;
    
    try {
      ReadARTSCAT(more_abs_lines,
                  tmpbasename + (species_data[*it].Name()) + ".xml",
                  fmin,
                  fmax,
                  globalquantumnumbers,
                  localquantumnumbers,
                  normalization_option,
                  mirroring_option,
                  population_option,
                  lineshapetype_option,
                  cutoff_option,
                  cutoff_value,
                  linemixinglimit_value,
                  verbosity);
      
      // Either find a line like this in the list of lines or start a new Lines
      for (auto& newband: more_abs_lines) {
        bool found = false;
        for (auto& band: abs_lines) {
          if (band.Match(newband)) {
            for (Index k=0; k<newband.NumLines(); k++) {
              band.AppendSingleLine(newband.Line(k));
              found = true;
            }
          }
        }
        if (not found) {
          abs_lines.push_back(newband);
        }
      }
    } catch (const std::exception& e) {
      if (not ignore_missing) {
        std::ostringstream os;
        os << "Errors in calls by *propmat_clearskyAddZeeman*:\n";
        os << e.what();
        throw std::runtime_error(os.str());
      } else {
        continue;
      }
    }
  }
  
  for (auto& band: abs_lines)
    band.sort_by_frequency();
  
  if (normalization_option != "None")
    abs_linesSetNormalization(abs_lines, normalization_option, verbosity);
  if (mirroring_option != "None")
    abs_linesSetMirroring(abs_lines, mirroring_option, verbosity);
  if (population_option != "None")
    abs_linesSetPopulation(abs_lines, population_option, verbosity);
  if (lineshapetype_option != "None")
    abs_linesSetLineShapeType(abs_lines, lineshapetype_option, verbosity);
  if (cutoff_option != "None")
    abs_linesSetCutoff(abs_lines, cutoff_option, cutoff_value, verbosity);
  abs_linesSetLinemixingLimit(abs_lines, linemixinglimit_value, verbosity);
}

enum class HitranType {
  Pre2004,
  Post2004,
  Online,
};

HitranType string2hitrantype(const String& s) {
  if (s == "Pre2004")
    return HitranType::Pre2004;
  else if (s == "Post2004")
    return HitranType::Post2004;
  else if (s == "Online")
    return HitranType::Online;
  else {
    std::ostringstream os;
    os << "The type \"" << s << "\" is an invalid hitran type\n";
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ReadHITRAN(ArrayOfAbsorptionLines& abs_lines,
                const String& hitran_file,
                const Numeric& fmin,
                const Numeric& fmax,
                const String& globalquantumnumbers,
                const String& localquantumnumbers,
                const String& hitran_type,
                const String& normalization_option,
                const String& mirroring_option,
                const String& population_option,
                const String& lineshapetype_option,
                const String& cutoff_option,
                const Numeric& cutoff_value,
                const Numeric& linemixinglimit_value,
                const Verbosity& verbosity)
{
  // Global numbers
  const std::vector<QuantumNumberType> global_nums = string2vecqn(globalquantumnumbers);
  
  // Local numbers
  const std::vector<QuantumNumberType> local_nums = string2vecqn(localquantumnumbers);
  
  // HITRAN type
  const auto hitran_version = string2hitrantype(hitran_type);
  
  // Hitran data
  ifstream is;
  open_input_file(is, hitran_file);
  
  std::vector<Absorption::SingleLineExternal> v(0);
  
  bool go_on = true;
  while (go_on) {
    switch (hitran_version) {
      case HitranType::Post2004:
        v.push_back(Absorption::ReadFromHitran2004Stream(is));
        break;
      case HitranType::Pre2004:
        v.push_back(Absorption::ReadFromHitran2001Stream(is));
        break;
      case HitranType::Online:
        v.push_back(Absorption::ReadFromHitranOnlineStream(is));
        break;
      default:
        throw std::runtime_error("A bad developer did not throw in time to stop this message.\nThe HitranType enum class has to be fully updated!\n");
    }
    
    if (v.back().bad) {
      v.pop_back();
      go_on = false;
    } else if (v.back().line.F0() < fmin) {
        v.pop_back();
    } else if (v.back().line.F0() > fmax) {
      v.pop_back();
      go_on = false;
    }
  }
  
  for (auto& x: v)
    x.line.Zeeman() = Zeeman::GetAdvancedModel(x.quantumidentity);
  
  auto x = Absorption::split_list_of_external_lines(v, local_nums, global_nums);
  abs_lines.resize(0);
  abs_lines.reserve(x.size());
  while (x.size()) {
    abs_lines.push_back(x.back());
    abs_lines.back().sort_by_frequency();
    x.pop_back();
  }
  
  abs_linesSetNormalization(abs_lines, normalization_option, verbosity);
  abs_linesSetMirroring(abs_lines, mirroring_option, verbosity);
  abs_linesSetPopulation(abs_lines, population_option, verbosity);
  abs_linesSetLineShapeType(abs_lines, lineshapetype_option, verbosity);
  abs_linesSetCutoff(abs_lines, cutoff_option, cutoff_value, verbosity);
  abs_linesSetLinemixingLimit(abs_lines, linemixinglimit_value, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ReadLBLRTM(ArrayOfAbsorptionLines& abs_lines,
                const String& lblrtm_file,
                const Numeric& fmin,
                const Numeric& fmax,
                const String& globalquantumnumbers,
                const String& localquantumnumbers,
                const String& normalization_option,
                const String& mirroring_option,
                const String& population_option,
                const String& lineshapetype_option,
                const String& cutoff_option,
                const Numeric& cutoff_value,
                const Numeric& linemixinglimit_value,
                const Verbosity& verbosity)
{
  // Global numbers
  const std::vector<QuantumNumberType> global_nums = string2vecqn(globalquantumnumbers);
  
  // Local numbers
  const std::vector<QuantumNumberType> local_nums = string2vecqn(localquantumnumbers);
  
  // LBLRTM data
  ifstream is;
  open_input_file(is, lblrtm_file);
  
  std::vector<Absorption::SingleLineExternal> v(0);
  
  bool go_on = true;
  while (go_on) {
    v.push_back(Absorption::ReadFromLBLRTMStream(is));
    
    if (v.back().bad) {
      v.pop_back();
      go_on = false;
    } else if (v.back().line.F0() < fmin) {
      v.pop_back();
    } else if (v.back().line.F0() > fmax) {
      v.pop_back();
      go_on = false;
    }
  }
  
  for (auto& x: v)
    x.line.Zeeman() = Zeeman::GetAdvancedModel(x.quantumidentity);
  
  auto x = Absorption::split_list_of_external_lines(v, local_nums, global_nums);
  abs_lines.resize(0);
  abs_lines.reserve(x.size());
  while (x.size()) {
    abs_lines.push_back(x.back());
    abs_lines.back().sort_by_frequency();
    x.pop_back();
  }
  
  abs_linesSetNormalization(abs_lines, normalization_option, verbosity);
  abs_linesSetMirroring(abs_lines, mirroring_option, verbosity);
  abs_linesSetPopulation(abs_lines, population_option, verbosity);
  abs_linesSetLineShapeType(abs_lines, lineshapetype_option, verbosity);
  abs_linesSetCutoff(abs_lines, cutoff_option, cutoff_value, verbosity);
  abs_linesSetLinemixingLimit(abs_lines, linemixinglimit_value, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ReadMytran2(ArrayOfAbsorptionLines& abs_lines,
                 const String& mytran2_file,
                 const Numeric& fmin,
                 const Numeric& fmax,
                 const String& globalquantumnumbers,
                 const String& localquantumnumbers,
                 const String& normalization_option,
                 const String& mirroring_option,
                 const String& population_option,
                 const String& lineshapetype_option,
                 const String& cutoff_option,
                 const Numeric& cutoff_value,
                 const Numeric& linemixinglimit_value,
                 const Verbosity& verbosity)
{
  // Global numbers
  const std::vector<QuantumNumberType> global_nums = string2vecqn(globalquantumnumbers);
  
  // Local numbers
  const std::vector<QuantumNumberType> local_nums = string2vecqn(localquantumnumbers);
  
  // LBLRTM data
  ifstream is;
  open_input_file(is, mytran2_file);
  
  std::vector<Absorption::SingleLineExternal> v(0);
  
  bool go_on = true;
  while (go_on) {
    v.push_back(Absorption::ReadFromMytran2Stream(is));
    
    if (v.back().bad) {
      v.pop_back();
      go_on = false;
    } else if (v.back().line.F0() < fmin) {
      v.pop_back();
    } else if (v.back().line.F0() > fmax) {
      v.pop_back();
      go_on = false;
    }
  }
  
  for (auto& x: v)
    x.line.Zeeman() = Zeeman::GetAdvancedModel(x.quantumidentity);
  
  auto x = Absorption::split_list_of_external_lines(v, local_nums, global_nums);
  abs_lines.resize(0);
  abs_lines.reserve(x.size());
  while (x.size()) {
    abs_lines.push_back(x.back());
    abs_lines.back().sort_by_frequency();
    x.pop_back();
  }
  
  abs_linesSetNormalization(abs_lines, normalization_option, verbosity);
  abs_linesSetMirroring(abs_lines, mirroring_option, verbosity);
  abs_linesSetPopulation(abs_lines, population_option, verbosity);
  abs_linesSetLineShapeType(abs_lines, lineshapetype_option, verbosity);
  abs_linesSetCutoff(abs_lines, cutoff_option, cutoff_value, verbosity);
  abs_linesSetLinemixingLimit(abs_lines, linemixinglimit_value, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ReadJPL(ArrayOfAbsorptionLines& abs_lines,
             const String& jpl_file,
             const Numeric& fmin,
             const Numeric& fmax,
             const String& globalquantumnumbers,
             const String& localquantumnumbers,
             const String& normalization_option,
             const String& mirroring_option,
             const String& population_option,
             const String& lineshapetype_option,
             const String& cutoff_option,
             const Numeric& cutoff_value,
             const Numeric& linemixinglimit_value,
             const Verbosity& verbosity)
{
  // Global numbers
  const std::vector<QuantumNumberType> global_nums = string2vecqn(globalquantumnumbers);
  
  // Local numbers
  const std::vector<QuantumNumberType> local_nums = string2vecqn(localquantumnumbers);
  
  // LBLRTM data
  ifstream is;
  open_input_file(is, jpl_file);
  
  std::vector<Absorption::SingleLineExternal> v(0);
  
  bool go_on = true;
  while (go_on) {
    v.push_back(Absorption::ReadFromJplStream(is));
    
    if (v.back().bad) {
      v.pop_back();
      go_on = false;
    } else if (v.back().line.F0() < fmin) {
      v.pop_back();
    } else if (v.back().line.F0() > fmax) {
      v.pop_back();
      go_on = false;
    }
  }
  
  for (auto& x: v)
    x.line.Zeeman() = Zeeman::GetAdvancedModel(x.quantumidentity);
  
  auto x = Absorption::split_list_of_external_lines(v, local_nums, global_nums);
  abs_lines.resize(0);
  abs_lines.reserve(x.size());
  while (x.size()) {
    abs_lines.push_back(x.back());
    abs_lines.back().sort_by_frequency();
    x.pop_back();
  }
  
  abs_linesSetNormalization(abs_lines, normalization_option, verbosity);
  abs_linesSetMirroring(abs_lines, mirroring_option, verbosity);
  abs_linesSetPopulation(abs_lines, population_option, verbosity);
  abs_linesSetLineShapeType(abs_lines, lineshapetype_option, verbosity);
  abs_linesSetCutoff(abs_lines, cutoff_option, cutoff_value, verbosity);
  abs_linesSetLinemixingLimit(abs_lines, linemixinglimit_value, verbosity);
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////// IO of AbsorptionLines
/////////////////////////////////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesWriteSplitXML(const String& output_format,
                            const ArrayOfAbsorptionLines& abs_lines,
                            const String& basename,
                            const Verbosity& verbosity)
{
  std::map<String, int> names;

  String true_basename = basename;
  if (not(true_basename.back() == '.' or true_basename.back() == '/'))
    true_basename += '.';

  for (auto& lines : abs_lines) {
    auto name = lines.SpeciesName();
    const String fname = true_basename + name;

    WriteXML(output_format, lines,
             fname + '.' + std::to_string(names[name]++) + ".xml",
             0, "", "", "", verbosity);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesWriteSpeciesSplitXML(const String& output_format,
                                   const ArrayOfAbsorptionLines& abs_lines,
                                   const String& basename,
                                   const Verbosity& verbosity)
{
  // Set the true name of the saving
  String true_basename = basename;
  if (not(true_basename.back() == '.' or true_basename.back() == '/'))
    true_basename += '.';
  
  // Find all species
  ArrayOfString specs(0);
  for (auto& band: abs_lines) {
    auto specname = band.SpeciesName();
    
    bool any = false;
    for (auto& thisname: specs) {
      if (any) break;
      else if (thisname == specname) any = true;
    }
    
    if (not any)
      specs.push_back(specname);
  }
  
  // Make all species into a species tag array
  Index throwaway;
  ArrayOfArrayOfSpeciesTag as;
  abs_speciesSet(as, throwaway, throwaway, specs, verbosity);
  
  // Split lines by species
  ArrayOfArrayOfAbsorptionLines alps;
  abs_lines_per_speciesCreateFromLines(alps, abs_lines, as, verbosity);
  
  // Save the arrays
  for (Index i=0; i<specs.nelem(); i++) {
    auto& name = specs[i];
    auto& lines = alps[i];
    
    WriteXML(output_format, lines,
             true_basename + name + ".xml",
             0, "", "", "", verbosity);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesWriteSplitXML(const String& output_format,
                                        const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                        const String& basename,
                                        const Verbosity& verbosity)
{
  std::map<String, int> names;

  String true_basename = basename;
  if (not(true_basename.back() == '.' or true_basename.back() == '/'))
    true_basename += '.';

  for (auto& spec_lines : abs_lines_per_species) {
    for (auto& lines : spec_lines) {
      auto name = lines.SpeciesName();
      const String fname = true_basename + name;

      WriteXML(output_format, lines,
              fname + '.' + std::to_string(names[name]++) + ".xml",
              0, "", "", "", verbosity);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesWriteSpeciesSplitXML(const String& output_format,
                                               const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                               const String& basename,
                                               const Verbosity& verbosity)
{ 
  // Compact to abs_lines
  ArrayOfAbsorptionLines abs_lines(0);
  for (auto& lines: abs_lines_per_species) {
    for (auto& band: lines) {
      abs_lines.push_back(band);
    }
  }
  
  // Save using the other function
  abs_linesWriteSpeciesSplitXML(output_format, abs_lines, basename, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesReadSplitCatalog(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                           const ArrayOfArrayOfSpeciesTag& abs_species,
                                           const String& basename,
                                           const Verbosity& verbosity)
{
  using global_data::species_data;
  
  // Build a set of species indices. Duplicates are ignored.
  std::set<Index> unique_species;
  for (auto asp = abs_species.begin(); asp != abs_species.end(); asp++) {
    for (ArrayOfSpeciesTag::const_iterator sp = asp->begin(); sp != asp->end(); sp++) {
      if (sp->Type() == SpeciesTag::TYPE_PLAIN || sp->Type() == SpeciesTag::TYPE_ZEEMAN) {
        unique_species.insert(sp->Species());
      }
    }
  }
  
  String tmpbasename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/') {
    tmpbasename += '.';
  }
  
  // Read catalogs for each identified species and put them all into
  // abs_lines
  ArrayOfAbsorptionLines abs_lines(0);
  for (auto it = unique_species.begin(); it != unique_species.end(); it++) {
    for (Index k=0; k<species_data[*it].Isotopologue().nelem(); k++) {
      String filename;
      Index i = 0;
      do {
        filename = tmpbasename + species_data[*it].FullName(k) + '.' + std::to_string(i) + ".xml";
        if (find_xml_file_existence(filename)) {
          abs_lines.push_back(AbsorptionLines());
          xml_read_from_file(filename, abs_lines.back(), verbosity);
          i++;
        } else {
          break;
        }
      } while (true);
    }
  }
  
  abs_lines_per_speciesCreateFromLines(abs_lines_per_species, abs_lines, abs_species, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReadSpeciesSplitCatalog(ArrayOfAbsorptionLines& abs_lines,
                                      const String& basename,
                                      const Index& robust,
                                      const Verbosity& verbosity)
{
  using global_data::species_data;
  
  CREATE_OUT3;
  std::size_t bands_found{0};
  
  String tmpbasename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/') {
    tmpbasename += '.';
  }
  
  // Read catalogs for each identified species and put them all into
  // abs_lines
  abs_lines.resize(0);
  for (auto it = species_data.begin(); it != species_data.end(); it++) {
    for (Index k=0; k<(*it).Isotopologue().nelem(); k++) {
      String filename;
      filename = tmpbasename + (*it).FullName(k) + ".xml";
      if (find_xml_file_existence(filename)) {
        ArrayOfAbsorptionLines speclines;
        xml_read_from_file(filename, speclines, verbosity);
        for (auto& band: speclines) {
          abs_lines.push_back(band);
          bands_found++;
        }
      }
    }
  }
  
  if (not bands_found and not robust)
    throw std::runtime_error("Cannot find any bands in the directory you are reading");
  else
    out3 << "Found " << bands_found << " lines\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesReadSpeciesSplitCatalog(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                  const ArrayOfArrayOfSpeciesTag& abs_species,
                                                  const String& basename,
                                                  const Verbosity& verbosity)
{
  using global_data::species_data;
  
  // Build a set of species indices. Duplicates are ignored.
  std::set<Index> unique_species;
  for (auto asp = abs_species.begin(); asp != abs_species.end(); asp++) {
    for (ArrayOfSpeciesTag::const_iterator sp = asp->begin(); sp != asp->end(); sp++) {
      if (sp->Type() == SpeciesTag::TYPE_PLAIN || sp->Type() == SpeciesTag::TYPE_ZEEMAN) {
        unique_species.insert(sp->Species());
      }
    }
  }
  
  String tmpbasename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/') {
    tmpbasename += '.';
  }
  
  // Read catalogs for each identified species and put them all into
  // abs_lines
  ArrayOfAbsorptionLines abs_lines(0);
  for (auto it = unique_species.begin(); it != unique_species.end(); it++) {
    for (Index k=0; k<species_data[*it].Isotopologue().nelem(); k++) {
      String filename;
      filename = tmpbasename + species_data[*it].FullName(k) + ".xml";
      if (find_xml_file_existence(filename)) {
        ArrayOfAbsorptionLines speclines;
        xml_read_from_file(filename, speclines, verbosity);
        for (auto& band: speclines) {
          abs_lines.push_back(band);
        }
      }
    }
  }
  
  abs_lines_per_speciesCreateFromLines(abs_lines_per_species, abs_lines, abs_species, verbosity);
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// Manipulation of quantum numbers
/////////////////////////////////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetQuantumNumberForMatch(ArrayOfAbsorptionLines& abs_lines,
                                       const String& qn,
                                       const Rational& x,
                                       const QuantumIdentifier& QI,
                                       const Verbosity&)
{
  auto QN = string2quantumnumbertype(qn);
  if (QN == QuantumNumberType::FINAL_ENTRY) {
    ostringstream os;
    os << "Usupported quantum number key: " << qn << '\n';
    throw std::runtime_error(os.str());
  }
  
  for (auto& band: abs_lines) {
    for (Index k=0; k<band.NumLines(); k++) {
      if (Absorption::id_in_line_lower(band, QI, k)) {
        band.LowerQuantumNumber(k, QN) = x;
      } else if (Absorption::id_in_line_upper(band, QI, k)) {
        band.UpperQuantumNumber(k, QN) = x;
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetQuantumNumberForMatch(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                   const String& qn,
                                                   const Rational& x,
                                                   const QuantumIdentifier& QI,
                                                   const Verbosity& v)
{
  for (auto& band: abs_lines_per_species)
    abs_linesSetQuantumNumberForMatch(band, qn, x, QI, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesTruncateGlobalQuantumNumbers(ArrayOfAbsorptionLines& abs_lines,
                                           const Verbosity&)
{
  ArrayOfAbsorptionLines x(0);
  
  for (auto& lines: abs_lines) {
    lines.truncate_global_quantum_numbers();
    
    Index match = -1;
    for (Index ind=0; ind<x.nelem(); ind++) {
      if (x[ind].Match(lines)) {
        match = ind;
        break;
      }
    }
    
    if (match < 0)
      x.push_back(lines);
    else {
      for (auto& line: lines.AllLines())
        x[match].AppendSingleLine(line);
    }
  }
  
  abs_lines = std::move(x);
  for (auto& lines: abs_lines)
    lines.sort_by_frequency();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesRemoveUnusedLocalQuantumNumbers(ArrayOfAbsorptionLines& abs_lines,
                                              const Verbosity&)
{
  for (auto& lines: abs_lines) {
    lines.RemoveUnusedLocalQuantums();
  }
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// Incremental change of abs_lines content
/////////////////////////////////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReplaceWithLines(ArrayOfAbsorptionLines& abs_lines, const ArrayOfAbsorptionLines& replacing_lines, const Verbosity&)
{
  for (auto& rlines: replacing_lines) {
    Index number_of_matching_bands = 0;
    for (auto& tlines: abs_lines) {
      if (tlines.Match(rlines)) {
        number_of_matching_bands++;
        for (auto& rline: rlines.AllLines()) {
          Index number_of_matching_single_lines = 0;
          for (auto& tline: tlines.AllLines()) {
            if (tline.SameQuantumNumbers(rline)) {
              number_of_matching_single_lines++;
              tline = rline;
            }
          }
          
          if (number_of_matching_single_lines not_eq 1) {
            throw std::runtime_error("Error: Did not match to a single single line.  This means the input data has not been understood.  This function needs exactly one match.");
          }
        }
        tlines.sort_by_frequency();
      }
    }
    
    if (number_of_matching_bands not_eq 1) {
      throw std::runtime_error("Error: Did not match to a single set of absorption lines.  This means the input data has not been understood.  This function needs exactly one match.");
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesAppendWithLines(ArrayOfAbsorptionLines& abs_lines, const ArrayOfAbsorptionLines& appending_lines, const Index& safe, const Verbosity&)
{
  if (safe) {
    std::vector<AbsorptionLines> addedlines(0);
    
    for (auto& alines: appending_lines) {
      Index number_of_matching_bands = 0;
      for (auto& tlines: abs_lines) {
        if (tlines.Match(alines)) {
          number_of_matching_bands++;
          for (auto& aline: alines.AllLines()) {
            Index number_of_matching_single_lines = 0;
            for (auto& tline: tlines.AllLines()) {
              if (tline.SameQuantumNumbers(aline)) {
                number_of_matching_single_lines++;
              }
            }
            if (number_of_matching_single_lines not_eq 0) {
              throw std::runtime_error("Error: Did match to a single single line.  This means the input data has not been understood.  This function needs exactly zero matches.");
            }
            tlines.AppendSingleLine(aline);
          }
          tlines.sort_by_frequency();
        }
      }
      
      if (number_of_matching_bands == 0)
        addedlines.push_back(alines);
      else if (number_of_matching_bands not_eq 1) {
        throw std::runtime_error("Error: Did not match to a single set of absorption lines.  This means the input data has not been understood.  This function needs exactly one or zero matches.");
      }
    }
    
    for (auto& lines: addedlines) {
      abs_lines.push_back(std::move(lines));
    }
  } else {
    for (auto& band: appending_lines)
      abs_lines.push_back(band);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesDeleteWithLines(ArrayOfAbsorptionLines& abs_lines, const ArrayOfAbsorptionLines& deleting_lines, const Verbosity&)
{
  for (auto& dlines: deleting_lines) {
    for (auto& tlines: abs_lines) {
      std::vector<Index> hits(0);
      
      if (tlines.Match(dlines)) {
        for (auto& dline: dlines.AllLines()) {
          for (Index i=0; i<tlines.NumLines(); i++) {
            if (tlines.AllLines()[i].SameQuantumNumbers(dline)) {
              hits.push_back(i);
            }
          }
        }
        
        // Sort and test the input
        std::sort(hits.begin(), hits.end());
        auto n = hits.size();
        hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
        if(n not_eq hits.size()) {
          throw std::runtime_error("Removing the same line more than once is not accepted");
        }
        
        // Remove the bad values
        while(not hits.empty()) {
          tlines.RemoveLine(hits.back());
          hits.pop_back();
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesDeleteLinesWithUndefinedLocalQuanta(ArrayOfAbsorptionLines& abs_lines, const Verbosity& verbosity)
{
  CREATE_OUT2;
  Index i = 0;
  
  for (auto& band: abs_lines) {
    std::vector<Index> deleters(0);
    
    for (Index iline=0; iline<band.NumLines(); iline++) {
      if (std::any_of(band.Line(iline).LowerQuantumNumbers().cbegin(),
                      band.Line(iline).LowerQuantumNumbers().cend(),
                      [](auto x) -> bool {return x.isUndefined();}) or
          std::any_of(band.Line(iline).UpperQuantumNumbers().cbegin(),
                      band.Line(iline).UpperQuantumNumbers().cend(),
                      [](auto x) -> bool {return x.isUndefined();})) {
        deleters.push_back(iline);
      }
    }
    
    while (deleters.size()) {
      band.RemoveLine(deleters.back());
      deleters.pop_back();
      i++;
    }
  }
  
  out2 << "Deleted " << i << " lines.\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesDeleteLinesWithBadOrHighChangingJs(ArrayOfAbsorptionLines& abs_lines, const Verbosity& verbosity)
{
  CREATE_OUT2;
  Index i = 0;
  
  for (auto& band: abs_lines) {
    std::vector<Index> deleters(0);
    
    for (Index iline=0; iline<band.NumLines(); iline++) {
      auto Jlo = band.LowerQuantumNumber(iline, QuantumNumberType::J);
      auto Jup = band.UpperQuantumNumber(iline, QuantumNumberType::J);
      
      if (Jlo.isUndefined() or Jup.isUndefined() or 1 < abs(Jup - Jlo)) {
        deleters.push_back(iline);
      }
    }
    
    while (deleters.size()) {
      band.RemoveLine(deleters.back());
      deleters.pop_back();
      i++;
    }
  }
  
  out2 << "Deleted " << i << " lines.\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesDeleteLinesWithQuantumNumberAbove(ArrayOfAbsorptionLines& abs_lines, const String& qn_id, const Index& qn_val, const Verbosity&)
{
  const auto qn = string2quantumnumbertype(qn_id);
  
  for (auto& band: abs_lines) {
    for (Index i=band.NumLines() - 1; i>=0; i--) {
      if (band.UpperQuantumNumber(i, qn) > qn_val or band.LowerQuantumNumber(i, qn) > qn_val) {
        band.RemoveLine(i);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetEmptyBroadeningParametersToEmpty(ArrayOfAbsorptionLines& abs_lines, const Verbosity& /*verbosity*/)
{
  for (auto& band: abs_lines) {
    std::array<bool, LineShape::nVars> var_is_empty;
    
    // Species by species can be empty, so loop each species by themselves
    for (Index ispec=0; ispec<band.NumBroadeners(); ispec++) {
      var_is_empty.fill(true);
      
      // Check if any variable in this band for any line is non-empty
      for (Index iline=0; iline<band.NumLines(); iline++) {
        for (Index ivar=0; ivar < LineShape::nVars; ivar++) {
          if (not LineShape::modelparameterEmpty(band.Line(iline).LineShape().Data()[ispec].Data()[ivar])) {
            var_is_empty[ivar] = false;
          }
        }
      }
      
      // Remove empty variables from the writing.  This will also speed up some calculations
      for (Index iline=0; iline<band.NumLines(); iline++) {
        for (Index ivar=0; ivar < LineShape::nVars; ivar++) {
          if (var_is_empty[ivar]) {
            band.Line(iline).LineShape().Data()[ispec].Data()[ivar].type = LineShape::TemperatureModel::None;
          }
        }
      }
    }
  }
}

void abs_linesKeepBands(ArrayOfAbsorptionLines& abs_lines, const QuantumIdentifier& qid, const Index& ignore_spec, const Index& ignore_isot, const Verbosity&)
{
  // Invalid setting
  if (ignore_spec and not ignore_isot) {
    throw std::runtime_error("Cannot ignore the species and not ignore the isotopologue");
  }
  
  // local ID is a transition even if an energy level is given
  QuantumIdentifier this_id=qid;
  if (QuantumIdentifier::ENERGY_LEVEL == this_id.Type()) {
    this_id.SetTransition();
    this_id.UpperQuantumNumbers() = qid.EnergyLevelQuantumNumbers();
    this_id.LowerQuantumNumbers() = qid.EnergyLevelQuantumNumbers();
  }
  
  for (auto& band: abs_lines) {
    if (ignore_spec)
      this_id.Species() = band.Species();
    if (ignore_isot)
      this_id.Isotopologue() = band.Isotopologue();
    
    const bool in_lines = band.QuantumIdentity().In(this_id);
    while (not in_lines and band.NumLines()) {
      band.RemoveLine(0);
    }
  }
}


void abs_linesCleanupEmpty(ArrayOfAbsorptionLines& abs_lines, const Verbosity&)
{
  for (Index i=abs_lines.nelem() - 1; i>=0; i--) {
    if (abs_lines[i].NumLines() == 0) {
      abs_lines.erase(abs_lines.begin() + i);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Change of style of computations for whole bands
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////
///////////////////////////////// Change Cutoff Frequency
/////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetCutoff(ArrayOfAbsorptionLines& abs_lines,
                        const String& type,
                        const Numeric& x,
                        const Verbosity&) 
{
  auto t = Absorption::string2cutofftype(type);
  for (auto& lines: abs_lines) {
    lines.Cutoff(t);
    lines.CutoffFreqValue(x);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetCutoff(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                    const String& type,
                                    const Numeric& x,
                                    const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetCutoff(abs_lines, type, x, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetCutoffForMatch(
  ArrayOfAbsorptionLines& abs_lines,
  const String& type,
  const Numeric& x,
  const QuantumIdentifier& QI,
  const Verbosity&)
{
  auto t = Absorption::string2cutofftype(type);
  for (auto& band: abs_lines) {
    if (QI.In(band.QuantumIdentity())) {
      band.Cutoff(t);
      band.CutoffFreqValue(x);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetCutoffForMatch(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const String& type,
  const Numeric& x,
  const QuantumIdentifier& QI,
  const Verbosity& verbosity)
{
  for (auto& lines: abs_lines_per_species) {
    abs_linesSetCutoffForMatch(lines, type, x, QI, verbosity);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetCutoffForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const String& type,
  const Numeric& x,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesSetCutoff(abs_lines_per_species[ispec], type, x, verbosity);
    }
  }
}

/////////////////////////////////////////////////////////
////////////////////////////////// Change Mirroring Style
/////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetMirroring(ArrayOfAbsorptionLines& abs_lines,
                           const String& type,
                           const Verbosity&) 
{
  auto t = Absorption::string2mirroringtype(type);
  for (auto& lines: abs_lines)
    lines.Mirroring(t);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetMirroring(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                       const String& type,
                                       const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetMirroring(abs_lines, type, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetMirroringForMatch(ArrayOfAbsorptionLines& abs_lines,
                                   const String& type,
                                   const QuantumIdentifier& QI,
                                   const Verbosity&) 
{
  auto t = Absorption::string2mirroringtype(type);
  for (auto& band: abs_lines) {
    if (QI.In(band.QuantumIdentity())) {
      band.Mirroring(t);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetMirroringForMatch(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                               const String& type,
                                               const QuantumIdentifier& QI,
                                               const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetMirroringForMatch(abs_lines, type, QI, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetMirroringForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const String& type,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesSetMirroring(abs_lines_per_species[ispec], type, verbosity);
    }
  }
}

/////////////////////////////////////////////////////////
///////////////////////////////// Change Population Style
/////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetPopulation(ArrayOfAbsorptionLines& abs_lines,
                            const String& type,
                            const Verbosity&) 
{
  auto t = Absorption::string2populationtype(type);
  for (auto& lines: abs_lines)
    lines.Population(t);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetPopulation(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                        const String& type,
                                        const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetPopulation(abs_lines, type, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetPopulationForMatch(ArrayOfAbsorptionLines& abs_lines,
                                    const String& type,
                                    const QuantumIdentifier& QI,
                                    const Verbosity&) 
{
  auto t = Absorption::string2populationtype(type);
  for (auto& lines: abs_lines)
    if (QI.In(lines.QuantumIdentity()))
      lines.Population(t);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetPopulationForMatch(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                const String& type,
                                                const QuantumIdentifier& QI,
                                                const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetPopulationForMatch(abs_lines, type, QI, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetPopulationForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const String& type,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesSetPopulation(abs_lines_per_species[ispec], type, verbosity);
    }
  }
}

/////////////////////////////////////////////////////////
////////////////////////////// Change Normalization Style
/////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetNormalization(ArrayOfAbsorptionLines& abs_lines,
                           const String& type,
                           const Verbosity&) 
{
  auto t = Absorption::string2normalizationtype(type);
  for (auto& lines: abs_lines)
    lines.Normalization(t);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetNormalization(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                           const String& type,
                                           const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetNormalization(abs_lines, type, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetNormalizationForMatch(ArrayOfAbsorptionLines& abs_lines,
                                       const String& type,
                                       const QuantumIdentifier& QI,
                                       const Verbosity&) 
{
  auto t = Absorption::string2normalizationtype(type);
  for (auto& lines: abs_lines)
    if (QI.In(lines.QuantumIdentity()))
      lines.Normalization(t);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetNormalizationForMatch(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                   const String& type,
                                                   const QuantumIdentifier& QI,
                                                   const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetNormalizationForMatch(abs_lines, type, QI, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetNormalizationForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const String& type,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesSetNormalization(abs_lines_per_species[ispec], type, verbosity);
    }
  }
}

/////////////////////////////////////////////////////////
///////////////////////////////// Change Line Shape Style
/////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetLineShapeType(ArrayOfAbsorptionLines& abs_lines,
                               const String& type,
                               const Verbosity&) 
{
  auto t = LineShape::string2shapetype(type);
  for (auto& lines: abs_lines)
    lines.LineShapeType(t);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetLineShapeType(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                           const String& type,
                                           const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetLineShapeType(abs_lines, type, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetLineShapeTypeForMatch(ArrayOfAbsorptionLines& abs_lines,
                                       const String& type,
                                       const QuantumIdentifier& QI,
                                       const Verbosity&) 
{
  auto t = LineShape::string2shapetype(type);
  for (auto& lines: abs_lines)
    if (QI.In(lines.QuantumIdentity()))
      lines.LineShapeType(t);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetLineShapeTypeForMatch(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                   const String& type,
                                                   const QuantumIdentifier& QI,
                                                   const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetLineShapeTypeForMatch(abs_lines, type, QI, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetLineShapeTypeForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const String& type,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesSetLineShapeType(abs_lines_per_species[ispec], type, verbosity);
    }
  }
}

/////////////////////////////////////////////////////////
//////////////////////////////// Change Line Mixing Limit
/////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetLinemixingLimit(ArrayOfAbsorptionLines& abs_lines,
                                 const Numeric& x,
                                 const Verbosity&) 
{
  for (auto& lines: abs_lines)
    lines.LinemixingLimit(x);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetLinemixingLimit(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                             const Numeric& x,
                                             const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetLinemixingLimit(abs_lines, x, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetLinemixingLimitForMatch(ArrayOfAbsorptionLines& abs_lines,
                                         const Numeric& x,
                                         const QuantumIdentifier& QI,
                                         const Verbosity&) 
{
  for (auto& lines: abs_lines)
    if (QI.In(lines.QuantumIdentity()))
      lines.LinemixingLimit(x);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetLinemixingLimitForMatch(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                     const Numeric& x,
                                                     const QuantumIdentifier& QI,
                                                     const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetLinemixingLimitForMatch(abs_lines, x, QI, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetLinemixingLimitForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const Numeric& x,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesSetLinemixingLimit(abs_lines_per_species[ispec], x, verbosity);
    }
  }
}

/////////////////////////////////////////////////////////
//////////////////////////// Change Reference Temperature
/////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetT0(ArrayOfAbsorptionLines& abs_lines,
                    const Numeric& x,
                    const Verbosity&) 
{
  for (auto& lines: abs_lines)
    lines.T0(x);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetT0(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                const Numeric& x,
                                const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetT0(abs_lines, x, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetT0ForMatch(ArrayOfAbsorptionLines& abs_lines,
                            const Numeric& x,
                            const QuantumIdentifier& QI,
                            const Verbosity&) 
{
  for (auto& lines: abs_lines)
    if (QI.In(lines.QuantumIdentity()))
      lines.T0(x);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetT0ForMatch(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                        const Numeric& x,
                                        const QuantumIdentifier& QI,
                                        const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetT0ForMatch(abs_lines, x, QI, v);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetT0ForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const Numeric& x,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesSetT0(abs_lines_per_species[ispec], x, verbosity);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////// Change values for individual lines
/////////////////////////////////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesChangeBaseParameterForMatchingLines(ArrayOfAbsorptionLines& abs_lines,
                                                  const QuantumIdentifier& QI,
                                                  const String& parameter_name,
                                                  const Numeric& change,
                                                  const Index& relative,
                                                  const Index& loose_matching,
                                                  const Verbosity&)
{
  Index parameter_switch = -1;

  if (parameter_name.nelem() == 0)
    throw std::runtime_error("parameter_name is empty.\n");
  else if (parameter_name == "Central Frequency" or
           parameter_name == "Line Center")
    parameter_switch = 0;
  else if (parameter_name == "Line Strength")
    parameter_switch = 1;
  else if (parameter_name == "Lower State Energy")
    parameter_switch = 4;
  else if (parameter_name == "Einstein Coefficient")
    parameter_switch = 5;
  else if (parameter_name == "Lower Statistical Weight")
    parameter_switch = 6;
  else if (parameter_name == "Upper Statistical Weight")
    parameter_switch = 7;
  else if (parameter_name == "Lower Zeeman Coefficient")
    parameter_switch = 8;
  else if (parameter_name == "Upper Zeeman Coefficient")
    parameter_switch = 9;

  for (auto& band: abs_lines) {
    for (Index k=0; k<band.NumLines(); k++) {
      if (loose_matching ? Absorption::id_in_line(band, QI, k)
                         : Absorption::line_is_id(band, QI, k)) {
        switch (parameter_switch) {
          case 0:  // "Central Frequency":
            if (relative == 0)
              band.F0(k) += change;
            else
              band.F0(k) *= 1.0e0 + change;
            break;
          case 1:  // "Line Strength":
            if (relative == 0)
              band.I0(k) += change;
            else
              band.I0(k) *= 1.0e0 + change;
            break;
          case 4:  // "Lower State Energy":
            if (relative == 0)
              band.E0(k) += change;
            else
              band.E0(k) *= 1.0e0 + change;
            break;
          case 5:  // "Einstein":
            if (relative == 0)
              band.A(k) += change;
            else
              band.A(k) *= 1.0e0 + change;
            break;
          case 6:  // "Lower Statistical Weight":
            if (relative == 0)
              band.g_low(k) += change;
            else
              band.g_low(k) *= 1.0e0 + change;
            break;
          case 7:  // "Upper Statistical Weight":
            if (relative == 0)
              band.g_upp(k) += change;
            else
              band.g_upp(k) *= 1.0e0 + change;
            break;
          case 8:  // "Lower Zeeman Coefficient":
            if (relative == 0)
              band.Line(k).Zeeman().gl() += change;
            else
              band.Line(k).Zeeman().gl() *= 1.0e0 + change;
            break;
          case 9:  // "Upper Zeeman Coefficient":
            if (relative == 0)
              band.Line(k).Zeeman().gu() += change;
            else
              band.Line(k).Zeeman().gu() *= 1.0e0 + change;
            break;
          default: {
            ostringstream os;
            os << "Usupported paramter_name\n"
              << parameter_name
              << "\nSee method description for supported parameter names.\n";
            throw std::runtime_error(os.str());
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesChangeBaseParameterForMatchingLines(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                              const QuantumIdentifier& QI,
                                                              const String& parameter_name,
                                                              const Numeric& change,
                                                              const Index& relative,
                                                              const Index& loose_matching,
                                                              const Verbosity& verbosity)
{
  for (auto& lines: abs_lines_per_species)
    abs_linesChangeBaseParameterForMatchingLines(lines, QI, parameter_name, change, relative, loose_matching, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesChangeBaseParameterForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const QuantumIdentifier& QI,
  const String& parameter_name,
  const Numeric& change,
  const Index& relative,
  const Index& loose_matching,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesChangeBaseParameterForMatchingLines(abs_lines_per_species[ispec], QI, parameter_name, change, relative, loose_matching, verbosity);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetBaseParameterForMatchingLines(ArrayOfAbsorptionLines& abs_lines,
                                               const QuantumIdentifier& QI,
                                               const String& parameter_name,
                                               const Numeric& x,
                                               const Index& loose_matching,
                                               const Verbosity&)
{
  Index parameter_switch = -1;
  
  if (parameter_name.nelem() == 0)
    throw std::runtime_error("parameter_name is empty.\n");
  else if (parameter_name == "Central Frequency" or
    parameter_name == "Line Center")
    parameter_switch = 0;
  else if (parameter_name == "Line Strength")
    parameter_switch = 1;
  else if (parameter_name == "Lower State Energy")
    parameter_switch = 4;
  else if (parameter_name == "Einstein Coefficient")
    parameter_switch = 5;
  else if (parameter_name == "Lower Statistical Weight")
    parameter_switch = 6;
  else if (parameter_name == "Upper Statistical Weight")
    parameter_switch = 7;
  else if (parameter_name == "Lower Zeeman Coefficient")
    parameter_switch = 8;
  else if (parameter_name == "Upper Zeeman Coefficient")
    parameter_switch = 9;
  
  for (auto& band: abs_lines) {
    for (Index k=0; k<band.NumLines(); k++) {
      if (loose_matching ? Absorption::id_in_line(band, QI, k)
        : Absorption::line_is_id(band, QI, k)) {
        switch (parameter_switch) {
          case 0:  // "Central Frequency":
            band.F0(k) = x;
            break;
          case 1:  // "Line Strength":
            band.I0(k) = x;
            break;
          case 4:  // "Lower State Energy":
            band.E0(k) = x;
            break;
          case 5:  // "Einstein":
            band.A(k) = x;
            break;
          case 6:  // "Lower Statistical Weight":
            band.g_low(k) = x;
            break;
          case 7:  // "Upper Statistical Weight":
            band.g_upp(k) = x;
            break;
          case 8:
            band.Line(k).Zeeman().gl() = x;
            break;
          case 9:
            band.Line(k).Zeeman().gu() = x;
            break;
          default: {
            ostringstream os;
            os << "Usupported paramter_name\n"
            << parameter_name
            << "\nSee method description for supported parameter names.\n";
            throw std::runtime_error(os.str());
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetBaseParameterForMatchingLines(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                           const QuantumIdentifier& QI,
                                                           const String& parameter_name,
                                                           const Numeric& change,
                                                           const Index& loose_matching,
                                                           const Verbosity& verbosity)
{
  for (auto& lines: abs_lines_per_species)
    abs_linesSetBaseParameterForMatchingLines(lines, QI, parameter_name, change, loose_matching, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetBaseParameterForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const QuantumIdentifier& QI,
  const String& parameter_name,
  const Numeric& change,
  const Index& loose_matching,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesSetBaseParameterForMatchingLines(abs_lines_per_species[ispec], QI, parameter_name, change, loose_matching, verbosity);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesChangeLineShapeModelParameterForMatchingLines(
    ArrayOfAbsorptionLines& abs_lines,
    const QuantumIdentifier& QI,
    const String& parameter,
    const String& coefficient,
    const String& species,
    const Numeric& x,
    const Index& relative,
    const Verbosity&)
{
  const bool do_self = species == LineShape::self_broadening;
  const bool do_bath = species == LineShape::bath_broadening;
  
  // Set the spec index if possible
  const Index spec = (do_self or do_bath) ? -1 : SpeciesTag(species).Species();

  const LineShape::Variable var = LineShape::string2variable(parameter);
  
  for (auto& band: abs_lines) {
    for (Index k=0; k<band.NumLines(); k++) {
      if (Absorption::id_in_line(band, QI, k)) {
        if (do_self and band.Self()) {
          if (relative) {
            SingleModelParameter(band.Line(k).LineShape().Data().front().Data()[Index(var)], coefficient) *= 1 + x;
          } else {
            SingleModelParameter(band.Line(k).LineShape().Data().front().Data()[Index(var)], coefficient) += x;
          }
        } else if (do_bath and band.Bath()) {
          if (relative) {
            SingleModelParameter(band.Line(k).LineShape().Data().back().Data()[Index(var)], coefficient) *= 1 + x;
          } else {
            SingleModelParameter(band.Line(k).LineShape().Data().back().Data()[Index(var)], coefficient) += x;
          }
        } else {
          for (Index i=band.Self(); i<band.BroadeningSpecies().nelem()-band.Bath(); i++) {
            if (spec == band.BroadeningSpecies()[i].Species()) {
              if (relative) {
                SingleModelParameter(band.Line(k).LineShape().Data()[i].Data()[Index(var)], coefficient) *= 1 + x;
              } else {
                SingleModelParameter(band.Line(k).LineShape().Data()[i].Data()[Index(var)], coefficient) += x;
              }
            }
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesChangeLineShapeModelParameterForMatchingLines(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const QuantumIdentifier& QI,
  const String& parameter,
  const String& coefficient,
  const String& species,
  const Numeric& x,
  const Index& relative,
  const Verbosity& verbosity)
{
  for (auto& lines: abs_lines_per_species)
    abs_linesChangeLineShapeModelParameterForMatchingLines(
      lines, QI, parameter, coefficient, species, x, relative, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesChangeLineShapeModelParameterForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const QuantumIdentifier& QI,
  const String& parameter,
  const String& coefficient,
  const String& species,
  const Numeric& x,
  const Index& relative,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesChangeLineShapeModelParameterForMatchingLines(
        abs_lines_per_species[ispec], QI, parameter, coefficient, species, x, relative, verbosity);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetLineShapeModelParameterForMatchingLines(
    ArrayOfAbsorptionLines& abs_lines,
    const QuantumIdentifier& QI,
    const String& parameter,
    const String& coefficient,
    const String& species,
    const Numeric& new_value,
    const Verbosity&)
{
  const bool do_self = species == LineShape::self_broadening;
  const bool do_bath = species == LineShape::bath_broadening;
  
  // Set the spec index if possible
  const Index spec = (do_self or do_bath) ? -1 : SpeciesTag(species).Species();

  const LineShape::Variable var = LineShape::string2variable(parameter);
  
  for (auto& band: abs_lines) {
    for (Index k=0; k<band.NumLines(); k++) {
      if (Absorption::id_in_line(band, QI, k)) {
        if (do_self and band.Self()) {
          SingleModelParameter(band.Line(k).LineShape().Data().front().Data()[Index(var)], coefficient) = new_value;
        } else if (do_bath and band.Bath()) {
          SingleModelParameter(band.Line(k).LineShape().Data().back().Data()[Index(var)], coefficient) = new_value;
        } else {
          for (Index i=band.Self(); i<band.BroadeningSpecies().nelem()-band.Bath(); i++) {
            if (spec == band.BroadeningSpecies()[i].Species()) {
              SingleModelParameter(band.Line(k).LineShape().Data()[i].Data()[Index(var)], coefficient) = new_value;
            }
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetLineShapeModelParameterForMatchingLines(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const QuantumIdentifier& QI,
  const String& parameter,
  const String& coefficient,
  const String& species,
  const Numeric& new_value,
  const Verbosity& verbosity)
{
  for (auto& lines: abs_lines_per_species)
    abs_linesSetLineShapeModelParameterForMatchingLines(
      lines, QI, parameter, coefficient, species, new_value, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetLineShapeModelParameterForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const QuantumIdentifier& QI,
  const String& parameter,
  const String& coefficient,
  const String& species,
  const Numeric& new_value,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesSetLineShapeModelParameterForMatchingLines(
        abs_lines_per_species[ispec], QI, parameter, coefficient, species, new_value, verbosity);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// Change values for individual levels
/////////////////////////////////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesChangeBaseParameterForMatchingLevel(ArrayOfAbsorptionLines& abs_lines,
                                                  const QuantumIdentifier& QI,
                                                  const String& parameter_name,
                                                  const Numeric& change,
                                                  const Index& relative,
                                                  const Verbosity&)
{
  if (QI.Type() not_eq QuantumIdentifier::ENERGY_LEVEL) {
    std::ostringstream os;
    os << "Bad input.  Must be energy level.  Is: " << QI << '\n';
    throw std::runtime_error(os.str());
  }
  
  Index parameter_switch = -1;
  
  if (parameter_name.nelem() == 0)
    throw std::runtime_error("parameter_name is empty.\n");
  else if (parameter_name == "Statistical Weight")
    parameter_switch = 1;
  else if (parameter_name == "Zeeman Coefficient")
    parameter_switch = 2;
  
  for (auto& band: abs_lines) {
    for (Index k=0; k<band.NumLines(); k++) {
      if (Absorption::id_in_line_lower(band, QI, k)) {
        switch (parameter_switch) {
          case 1:  // "Statistical Weight":
            if (relative == 0)
              band.g_low(k) += change;
            else
              band.g_low(k) *= 1.0e0 + change;
            break;
          case 2:  // "Zeeman Coefficient":
            if (relative == 0)
              band.Line(k).Zeeman().gl() += change;
            else
              band.Line(k).Zeeman().gl() *= 1.0e0 + change;
            break;
          default: {
            ostringstream os;
            os << "Usupported paramter_name\n"
            << parameter_name
            << "\nSee method description for supported parameter names.\n";
            throw std::runtime_error(os.str());
          }
        }
      } else if (Absorption::id_in_line_upper(band, QI, k)) {
        switch (parameter_switch) {
          case 1:  // "Statistical Weight":
            if (relative == 0)
              band.g_upp(k) += change;
            else
              band.g_upp(k) *= 1.0e0 + change;
            break;
          case 2:  // "Zeeman Coefficient":
            if (relative == 0)
              band.Line(k).Zeeman().gu() += change;
            else
              band.Line(k).Zeeman().gu() *= 1.0e0 + change;
            break;
          default: {
            ostringstream os;
            os << "Usupported paramter_name\n"
            << parameter_name
            << "\nSee method description for supported parameter names.\n";
            throw std::runtime_error(os.str());
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesChangeBaseParameterForMatchingLevel(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                              const QuantumIdentifier& QI,
                                                              const String& parameter_name,
                                                              const Numeric& change,
                                                              const Index& relative,
                                                              const Verbosity& verbosity)
{
  for (auto& lines: abs_lines_per_species)
    abs_linesChangeBaseParameterForMatchingLevel(lines, QI, parameter_name, change, relative, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesChangeBaseParameterForMatchingLevels(ArrayOfAbsorptionLines& abs_lines,
                                                   const ArrayOfQuantumIdentifier& QID,
                                                   const String& parameter_name,
                                                   const Vector& change,
                                                   const Index& relative,
                                                   const Verbosity& verbosity)
{
  if (QID.nelem() not_eq change.nelem()) {
    throw std::runtime_error("Mismatch between QID and change input lengths not allowed");
  }
  
  for (Index iq=0; iq<QID.nelem(); iq++)
    abs_linesChangeBaseParameterForMatchingLevel(abs_lines, QID[iq], parameter_name, change[iq], relative, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesChangeBaseParameterForMatchingLevels(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                               const ArrayOfQuantumIdentifier& QID,
                                                               const String& parameter_name,
                                                               const Vector& change,
                                                               const Index& relative,
                                                               const Verbosity& verbosity)
{
  if (QID.nelem() not_eq change.nelem()) {
    throw std::runtime_error("Mismatch between QID and change input lengths not allowed");
  }
  
  for (Index iq=0; iq<QID.nelem(); iq++)
    for (auto& lines: abs_lines_per_species)
      abs_linesChangeBaseParameterForMatchingLevel(lines, QID[iq], parameter_name, change[iq], relative, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetBaseParameterForMatchingLevel(ArrayOfAbsorptionLines& abs_lines,
                                               const QuantumIdentifier& QI,
                                               const String& parameter_name,
                                               const Numeric& x,
                                               const Verbosity&)
{
  if (QI.Type() not_eq QuantumIdentifier::ENERGY_LEVEL) {
    std::ostringstream os;
    os << "Bad input.  Must be energy level.  Is: " << QI << '\n';
    throw std::runtime_error(os.str());
  }
  
  Index parameter_switch = -1;
  
  if (parameter_name.nelem() == 0)
    throw std::runtime_error("parameter_name is empty.\n");
  else if (parameter_name == "Statistical Weight")
    parameter_switch = 1;
  else if (parameter_name == "Zeeman Coefficient")
    parameter_switch = 2;
  
  for (auto& band: abs_lines) {
    for (Index k=0; k<band.NumLines(); k++) {
      if (Absorption::id_in_line_lower(band, QI, k)) {
        switch (parameter_switch) {
          case 1:  // "Statistical Weight":
            band.g_low(k) = x;
            break;
          case 2:  // "Zeeman Coefficient":
            band.Line(k).Zeeman().gl() = x;
            break;
          default: {
            ostringstream os;
            os << "Usupported paramter_name\n"
            << parameter_name
            << "\nSee method description for supported parameter names.\n";
            throw std::runtime_error(os.str());
          }
        }
      } else if (Absorption::id_in_line_upper(band, QI, k)) {
        switch (parameter_switch) {
          case 1:  // "Statistical Weight":
            band.g_upp(k) = x;
            break;
          case 2:  // "Zeeman Coefficient":
            band.Line(k).Zeeman().gu() = x;
            break;
          default: {
            ostringstream os;
            os << "Usupported paramter_name\n"
            << parameter_name
            << "\nSee method description for supported parameter names.\n";
            throw std::runtime_error(os.str());
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetBaseParameterForMatchingLevel(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                           const QuantumIdentifier& QI,
                                                           const String& parameter_name,
                                                           const Numeric& change,
                                                           const Verbosity& verbosity)
{
  for (auto& lines: abs_lines_per_species)
    abs_linesSetBaseParameterForMatchingLevel(lines, QI, parameter_name, change, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetBaseParameterForMatchingLevels(ArrayOfAbsorptionLines& abs_lines,
                                                const ArrayOfQuantumIdentifier& QID,
                                                const String& parameter_name,
                                                const Vector& change,
                                                const Verbosity& verbosity)
{
  if (QID.nelem() not_eq change.nelem()) {
    throw std::runtime_error("Mismatch between QID and change input lengths not allowed");
  }
  
  for (Index iq=0; iq<QID.nelem(); iq++)
    abs_linesSetBaseParameterForMatchingLevel(abs_lines, QID[iq], parameter_name, change[iq], verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetBaseParameterForMatchingLevels(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                            const ArrayOfQuantumIdentifier& QID,
                                                            const String& parameter_name,
                                                            const Vector& change,
                                                            const Verbosity& verbosity)
{
  if (QID.nelem() not_eq change.nelem()) {
    throw std::runtime_error("Mismatch between QID and change input lengths not allowed");
  }
  
  for (Index iq=0; iq<QID.nelem(); iq++)
    for (auto& lines: abs_lines_per_species)
      abs_linesSetBaseParameterForMatchingLevel(lines, QID[iq], parameter_name, change[iq], verbosity);
}

/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Special functions for special problems
/////////////////////////////////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void nlteSetByQuantumIdentifiers(
    Index& nlte_do,
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const EnergyLevelMap& nlte_field,
    const Verbosity&)
{
  nlte_field.ThrowIfNotOK();
  
  if (nlte_field.Data().empty()) {
    nlte_do = 0;
    return;
  } else {
    nlte_do = 1;
  }
  
  const Absorption::PopulationType poptyp = nlte_field.Energies().empty() ? 
        Absorption::PopulationType::ByNLTEPopulationDistribution :
        Absorption::PopulationType::ByNLTEVibrationalTemperatures;

  for (auto& spec_lines: abs_lines_per_species) {
    for (auto& band: spec_lines) {
      for (auto& id: nlte_field.Levels()) {
        if (band.QuantumIdentity().UpperQuantumId().In(id) or
            band.QuantumIdentity().LowerQuantumId().In(id)) {
          for (Index k=0; k<band.NumLines(); k++) {
            if (poptyp==Absorption::PopulationType::ByNLTEPopulationDistribution and
                (not std::isnormal(band.A(k)) or band.A(k) < 0)) {
              std::ostringstream os;
              os << "Error in band deemed for NLTE calculations by population distribution\n"
                 << "some of the lines in the band below have a bad Einstein coefficient:\n"
                 << band << '\n';
              throw std::runtime_error(os.str());
            }
          }
          band.Population(poptyp);
        }
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Destructively manipulate content of whole catalog
/////////////////////////////////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetEmpty(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                   const ArrayOfArrayOfSpeciesTag& abs_species,
                                   const Verbosity&) {
  abs_lines_per_species = ArrayOfArrayOfAbsorptionLines(abs_species.nelem(), ArrayOfAbsorptionLines(0));
}

/* Workspace method: Doxygen documentation will be auto-generated */                               
void abs_linesCompact(ArrayOfAbsorptionLines& abs_lines,
                      const Vector& f_grid,
                      const Verbosity&)
{
  const Numeric fmax = max(f_grid);
  const Numeric fmin = min(f_grid);
  
  for (auto& band: abs_lines) {
    const Numeric fmean = (band.Cutoff() == Absorption::CutoffType::BandFixedFrequency) ? band.F_mean() : 0;
    for (Index k=band.NumLines()-1; k>=0; k--) {
      const Numeric fcut_upp = band.CutoffFreq(k);
      const Numeric fcut_low = band.CutoffFreqMinus(k, fmean);
      
      if (fmax < fcut_low or fmin > fcut_upp) {
        band.RemoveLine(k);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesCompact(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                  const Vector& f_grid,
                                  const Verbosity& verbosity)
{
  for (auto& lines: abs_lines_per_species) {
    abs_linesCompact(lines, f_grid, verbosity);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesRemoveBand(ArrayOfAbsorptionLines& abs_lines,
                                     const QuantumIdentifier& qid,
                                     const Verbosity&)
{
  for (Index i=0; i<abs_lines.nelem(); i++) {
    if (qid.In(abs_lines[i].QuantumIdentity())) {
      abs_lines.erase(abs_lines.begin()+i);
      break;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Manipulation of other ARTS variables based on AbsorptionLines
/////////////////////////////////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromAbsorptionLines(Vector& f_grid,
                               const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                               const Numeric& delta_f_low,
                               const Numeric& delta_f_upp,
                               const Index& num_freqs,
                               const Verbosity&)
{
  const Index n=nelem(abs_lines_per_species);
  
  if (delta_f_low >= delta_f_upp) {
    throw std::runtime_error("The lower frequency delta has to be smaller "
                             "than the upper frequency delta");
  } else if (num_freqs < 1) {
    throw std::runtime_error("Need more than zero frequency points");
  } else if (n < 1) {
    throw std::runtime_error("No lines found.  Error?  Use *VectorSet* "
                             "to resize *f_grid*");
  }
  
  f_grid.resize(n*num_freqs);
  
  Index i=0;
  for (auto& lines: abs_lines_per_species) {
    for (auto& band: lines) {
      for (Index k=0; k<band.NumLines(); k++) {
        if (num_freqs > 1) {
          nlinspace(f_grid[Range(i, num_freqs)], band.F0(k)+delta_f_low, band.F0(k)+delta_f_upp, num_freqs);
        } else {
          f_grid[i] = band.F0(k);
        }
        i += num_freqs;
      }
    }
  }
  
  auto tmp = MapToEigen(f_grid);
  std::sort(tmp.data(), tmp.data()+tmp.size(), [](auto a, auto b){return a < b;});
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// Print meta data about the lines
/////////////////////////////////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesPrintDefinedQuantumNumbers(const ArrayOfAbsorptionLines& abs_lines,
                                         const Verbosity& verbosity)
{
  CREATE_OUT0;
  
  std::map<Index, Index> qns;
  
  for (auto& band: abs_lines) {
    for (Index iline=0; iline<band.NumLines(); iline++) {
      for (Index iqn=0; iqn<Index(QuantumNumberType::FINAL_ENTRY); iqn++) {
        if (band.LowerQuantumNumber(iline, QuantumNumberType(iqn)).isDefined() or
            band.UpperQuantumNumber(iline, QuantumNumberType(iqn)).isDefined()) {
          qns[iqn]++;
        }
      }
    }
  }
  
  for (auto& qn: qns) {
    out0 << quantumnumbertype2string(QuantumNumberType(qn.first)) << ':' << ' ' << qn.second << '\n';
  }
}
