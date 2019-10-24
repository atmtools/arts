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
#include "xml_io_private.h"
#include "auto_md.h"
#include "file.h"


void ReadArrayOfARTSCAT(ArrayOfAbsorptionLines& abs_lines,
                        const String& artscat_file,
                        const Numeric& fmin,
                        const Numeric& fmax,
                        const String& globalquantumnumbers,
                        const String& localquantumnumbers,
                        const Verbosity& verbosity)
{
  // Take care of quantum numbers
  String tmp_string;
  
  // Global numbers
  std::vector<QuantumNumberType> global_nums(0);
  if (globalquantumnumbers not_eq "") {
    std::istringstream global_str(globalquantumnumbers);
    while (not global_str.eof()) {
      global_str >> tmp_string; 
      global_nums.push_back(string2quantumnumbertype(tmp_string));
    }
  }
  
  // Local numbers
  std::vector<QuantumNumberType> local_nums(0);
  if (localquantumnumbers not_eq "") {
    std::istringstream local_str(localquantumnumbers);
    while (not local_str.eof()) {
      local_str >> tmp_string;
      local_nums.push_back(string2quantumnumbertype(tmp_string));
    }
  }
  
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
  for (auto& lines: x) {
    lines.sort_by_frequency();
    abs_lines.push_back(lines);
  }
}


void ReadARTSCAT(ArrayOfAbsorptionLines& abs_lines,
                 const String& artscat_file,
                 const Numeric& fmin,
                 const Numeric& fmax,
                 const String& globalquantumnumbers,
                 const String& localquantumnumbers,
                 const Verbosity& verbosity)
{
  // Take care of quantum numbers
  String tmp_string;
  
  // Global numbers
  std::vector<QuantumNumberType> global_nums(0);
  if (globalquantumnumbers not_eq "") {
    std::istringstream global_str(globalquantumnumbers);
    while (not global_str.eof()) {
      global_str >> tmp_string; 
      global_nums.push_back(string2quantumnumbertype(tmp_string));
    }
  }
  
  // Local numbers
  std::vector<QuantumNumberType> local_nums(0);
  if (localquantumnumbers not_eq "") {
    std::istringstream local_str(localquantumnumbers);
    while (not local_str.eof()) {
      local_str >> tmp_string;
      local_nums.push_back(string2quantumnumbertype(tmp_string));
    }
  }
  
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
  
  for (auto& x: v)
    x.line.Zeeman() = Zeeman::GetAdvancedModel(x.quantumidentity);
  
  auto x = Absorption::split_list_of_external_lines(v, local_nums, global_nums);
  abs_lines.resize(0);
  abs_lines.reserve(x.size());
  for (auto& lines: x)
    abs_lines.push_back(lines);
}

void ReadHITRAN(ArrayOfAbsorptionLines& abs_lines,
                const String& hitran_file,
                const Numeric& fmin,
                const Numeric& fmax,
                const String& globalquantumnumbers,
                const String& localquantumnumbers,
                const String& hitran_type,
                const Verbosity&)
{
  // Take care of quantum numbers
  String tmp_string;
  
  // Global numbers
  std::vector<QuantumNumberType> global_nums(0);
  if (globalquantumnumbers not_eq "") {
    std::istringstream global_str(globalquantumnumbers);
    while (not global_str.eof()) {
      global_str >> tmp_string; 
      global_nums.push_back(string2quantumnumbertype(tmp_string));
    }
  }
  
  // Local numbers
  std::vector<QuantumNumberType> local_nums(0);
  if (localquantumnumbers not_eq "") {
    std::istringstream local_str(localquantumnumbers);
    while (not local_str.eof()) {
      local_str >> tmp_string;
      local_nums.push_back(string2quantumnumbertype(tmp_string));
    }
  }
  
  // HITRAN type
  bool new_type = hitran_type=="Pre2004" ? false : true;
  
  // Hitran data
  ifstream is;
  open_input_file(is, hitran_file);
  
  std::vector<Absorption::SingleLineExternal> v(0);
  
  bool go_on = true;
  while (go_on) {
    if (new_type)
      v.push_back(Absorption::ReadFromHitran2004Stream(is));
    else
      v.push_back(Absorption::ReadFromHitran2001Stream(is));
    
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
  for (auto& lines: x)
    abs_lines.push_back(lines);
}

void ReadLBLRTM(ArrayOfAbsorptionLines& abs_lines,
                const String& lblrtm_file,
                const Numeric& fmin,
                const Numeric& fmax,
                const String& globalquantumnumbers,
                const String& localquantumnumbers,
                const Verbosity&)
{
  // Take care of quantum numbers
  String tmp_string;
  
  // Global numbers
  std::vector<QuantumNumberType> global_nums(0);
  if (globalquantumnumbers not_eq "") {
    std::istringstream global_str(globalquantumnumbers);
    while (not global_str.eof()) {
      global_str >> tmp_string; 
      global_nums.push_back(string2quantumnumbertype(tmp_string));
    }
  }
  
  // Local numbers
  std::vector<QuantumNumberType> local_nums(0);
  if (localquantumnumbers not_eq "") {
    std::istringstream local_str(localquantumnumbers);
    while (not local_str.eof()) {
      local_str >> tmp_string;
      local_nums.push_back(string2quantumnumbertype(tmp_string));
    }
  }
  
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
  for (auto& lines: x)
    abs_lines.push_back(lines);
}

void ReadMytran2(ArrayOfAbsorptionLines& abs_lines,
                 const String& mytran2_file,
                 const Numeric& fmin,
                 const Numeric& fmax,
                 const String& globalquantumnumbers,
                 const String& localquantumnumbers,
                 const Verbosity&)
{
  // Take care of quantum numbers
  String tmp_string;
  
  // Global numbers
  std::vector<QuantumNumberType> global_nums(0);
  if (globalquantumnumbers not_eq "") {
    std::istringstream global_str(globalquantumnumbers);
    while (not global_str.eof()) {
      global_str >> tmp_string; 
      global_nums.push_back(string2quantumnumbertype(tmp_string));
    }
  }
  
  // Local numbers
  std::vector<QuantumNumberType> local_nums(0);
  if (localquantumnumbers not_eq "") {
    std::istringstream local_str(localquantumnumbers);
    while (not local_str.eof()) {
      local_str >> tmp_string;
      local_nums.push_back(string2quantumnumbertype(tmp_string));
    }
  }
  
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
  for (auto& lines: x)
    abs_lines.push_back(lines);
}

void ReadJPL(ArrayOfAbsorptionLines& abs_lines,
             const String& jpl_file,
             const Numeric& fmin,
             const Numeric& fmax,
             const String& globalquantumnumbers,
             const String& localquantumnumbers,
             const Verbosity&)
{
  // Take care of quantum numbers
  String tmp_string;
  
  // Global numbers
  std::vector<QuantumNumberType> global_nums(0);
  if (globalquantumnumbers not_eq "") {
    std::istringstream global_str(globalquantumnumbers);
    while (not global_str.eof()) {
      global_str >> tmp_string; 
      global_nums.push_back(string2quantumnumbertype(tmp_string));
    }
  }
  
  // Local numbers
  std::vector<QuantumNumberType> local_nums(0);
  if (localquantumnumbers not_eq "") {
    std::istringstream local_str(localquantumnumbers);
    while (not local_str.eof()) {
      local_str >> tmp_string;
      local_nums.push_back(string2quantumnumbertype(tmp_string));
    }
  }
  
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
  for (auto& lines: x)
    abs_lines.push_back(lines);
}

void abs_linesWriteSplitXML(const ArrayOfAbsorptionLines& abs_lines,
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

    WriteXML("ascii", lines,
             fname + '.' + std::to_string(names[name]++) + ".xml",
             0, "", "", "", verbosity);
  }
}

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

void abs_linesRemoveUnusedLocalQuantumNumbers(ArrayOfAbsorptionLines& abs_lines,
                                              const Verbosity&)
{
  for (auto& lines: abs_lines) {
    lines.RemoveUnusedLocalQuantums();
  }
}


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


void abs_linesAppendWithLines(ArrayOfAbsorptionLines& abs_lines, const ArrayOfAbsorptionLines& appending_lines, const Verbosity&)
{
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
}


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
        std::unique(hits.begin(), hits.end());
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

void abs_linesSetMirroring(ArrayOfAbsorptionLines& abs_lines,
                           const String& type,
                           const Verbosity&) 
{
  auto t = Absorption::string2mirroringtype(type);
  for (auto& lines: abs_lines)
    lines.Mirroring(t);
}

void abs_linesSetPopulation(ArrayOfAbsorptionLines& abs_lines,
                            const String& type,
                            const Verbosity&) 
{
  auto t = Absorption::string2populationtype(type);
  for (auto& lines: abs_lines)
    lines.Population(t);
}

void abs_linesSetNormalization(ArrayOfAbsorptionLines& abs_lines,
                           const String& type,
                           const Verbosity&) 
{
  auto t = Absorption::string2normalizationtype(type);
  for (auto& lines: abs_lines)
    lines.Normalization(t);
}

void abs_linesSetLineShapeType(ArrayOfAbsorptionLines& abs_lines,
                               const String& type,
                               const Verbosity&) 
{
  auto t = LineShape::string2shapetype(type);
  for (auto& lines: abs_lines)
    lines.LineShapeType(t);
}

void abs_linesSetLinemixingLimit(ArrayOfAbsorptionLines& abs_lines,
                                 const Numeric& x,
                                 const Verbosity&) 
{
  for (auto& lines: abs_lines)
    lines.LinemixingLimit(x);
}

void abs_lines_per_speciesSetCutoff(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                    const String& type,
                                    const Numeric& x,
                                    const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetCutoff(abs_lines, type, x, v);
}

void abs_lines_per_speciesSetMirroring(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                       const String& type,
                                       const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetMirroring(abs_lines, type, v);
}

void abs_lines_per_speciesSetPopulation(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                        const String& type,
                                        const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetPopulation(abs_lines, type, v);
}

void abs_lines_per_speciesSetNormalization(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                           const String& type,
                                           const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetNormalization(abs_lines, type, v);
}

void abs_lines_per_speciesSetLineShapeType(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                           const String& type,
                                           const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetLineShapeType(abs_lines, type, v);
}

void abs_lines_per_speciesSetLinemixingLimit(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                             const Numeric& x,
                                             const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetLinemixingLimit(abs_lines, x, v);
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

void abs_linesSetT0(ArrayOfAbsorptionLines& abs_lines,
                    const Numeric& x,
                    const Verbosity&) 
{
  for (auto& lines: abs_lines)
    lines.T0(x);
}

void abs_lines_per_speciesSetT0(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                const Numeric& x,
                                const Verbosity& v) 
{
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesSetT0(abs_lines, x, v);
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
  else if (parameter_name == "Einstein")
    parameter_switch = 5;
  else if (parameter_name == "Lower Statistical Weight")
    parameter_switch = 6;
  else if (parameter_name == "Upper Statistical Weight")
    parameter_switch = 7;

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
          default: {
            ostringstream os;
            os << "Usupported paramter_name\n"
              << parameter_name
              << "\nSee method description for supported parameter names.\n";
            throw std::runtime_error(os.str());
            break;
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
  else if (parameter_name == "Einstein")
    parameter_switch = 5;
  else if (parameter_name == "Lower Statistical Weight")
    parameter_switch = 6;
  else if (parameter_name == "Upper Statistical Weight")
    parameter_switch = 7;
  
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
          default: {
            ostringstream os;
            os << "Usupported paramter_name\n"
            << parameter_name
            << "\nSee method description for supported parameter names.\n";
            throw std::runtime_error(os.str());
            break;
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
  // Set the spec index to negative 101 if species is self and negative 102 if species is bath and to the species otherwise
  const Index spec = species == LineShape::self_broadening ? -101 :
                     species == LineShape::bath_broadening ? -102 :
                     SpeciesTag(species).Species();

  const LineShape::Variable var = LineShape::string2variable(parameter);
  
  for (auto& band: abs_lines) {
    for (Index k=0; k<band.NumLines(); k++) {
      if (Absorption::id_in_line(band, QI, k)) {
        if (spec == -101 and band.Self()) {
          if (relative) {
            SingleModelParameter(band.Line(k).LineShape().Data().front().Data()[Index(var)], coefficient) *= 1 + x;
          } else {
            SingleModelParameter(band.Line(k).LineShape().Data().front().Data()[Index(var)], coefficient) += x;
          }
        } else if (spec == -102 and band.Bath()) {
          if (relative) {
            SingleModelParameter(band.Line(k).LineShape().Data().back().Data()[Index(var)], coefficient) *= 1 + x;
          } else {
            SingleModelParameter(band.Line(k).LineShape().Data().back().Data()[Index(var)], coefficient) += x;
          }
        } else {
          for (Index i=0; i<band.BroadeningSpecies().nelem(); i++) {
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
  // Set the spec index to negative 101 if species is self and negative 102 if species is bath and to the species otherwise
  const Index spec = species == LineShape::self_broadening ? -101 :
                     species == LineShape::bath_broadening ? -102 :
                     SpeciesTag(species).Species();

  const LineShape::Variable var = LineShape::string2variable(parameter);
  
  for (auto& band: abs_lines) {
    for (Index k=0; k<band.NumLines(); k++) {
      if (Absorption::id_in_line(band, QI, k)) {
        if (spec == -101 and band.Self()) {
          SingleModelParameter(band.Line(k).LineShape().Data().front().Data()[Index(var)], coefficient) = new_value;
        } else if (spec == -102 and band.Bath()) {
          SingleModelParameter(band.Line(k).LineShape().Data().back().Data()[Index(var)], coefficient) = new_value;
        } else {
          for (Index i=0; i<band.BroadeningSpecies().nelem(); i++) {
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

/* Workspace method: Doxygen documentation will be auto-generated */
void nlteSetByQuantumIdentifiers(
    Index& nlte_do,
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const EnergyLevelMap& nlte_field,
    const Verbosity&) {
  nlte_field.ThrowIfNotOK();
  
  if (nlte_field.Data().empty()) {
    nlte_do = 0;
    return;
  } else {
    nlte_do = 1;
  }
  
  const Absorption::PopulationType poptyp = nlte_field.Energies().nelem() == 0 ? 
        Absorption::PopulationType::ByNLTEPopulationDistribution :
        Absorption::PopulationType::ByNLTEVibrationalTemperatures;

  for (auto& id: nlte_field.Levels()) {
    for (auto& spec_lines: abs_lines_per_species) {
      for (auto& band: spec_lines) {
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
