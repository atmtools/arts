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


void ReadARTSCAT(ArrayOfAbsorptionLines& abs_lines,
                 const String& artscat_file,
                 const String& globalquantumnumbers,
                 const String& localquantumnumbers,
                 const Verbosity& verbosity)
{
  // Take care of quantum numbers
  String tmp_string;
  
  // Global numbers
  std::vector<QuantumNumberType> global_nums(0);
  std::istringstream global_str(globalquantumnumbers);
  while (not global_str.eof()) {
    global_str >> tmp_string; 
    global_nums.push_back(string2quantumnumbertype(tmp_string));
  }
  
  // Local numbers
  std::vector<QuantumNumberType> local_nums(0);
  std::istringstream local_str(localquantumnumbers);
  while (not local_str.eof()) {
    local_str >> tmp_string;
    local_nums.push_back(string2quantumnumbertype(tmp_string));
  }
  CREATE_OUT2;
  
  ArtsXMLTag tag(verbosity);
  Index nelem;
  
  // ARTSCAT data
  ifstream is_xml;
  xml_open_input_file(is_xml, artscat_file, verbosity);
  auto a = FILE_TYPE_ASCII;
  auto b = NUMERIC_TYPE_DOUBLE;
  auto c = ENDIAN_TYPE_LITTLE;
  xml_read_header_from_stream(is_xml, a,b,c,
                              verbosity);
  
  tag.read_from_stream(is_xml);
  tag.check_name("ArrayOfLineRecord");
  
  tag.get_attribute_value("nelem", nelem);
  
  LineRecord dummy_line_record;
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
    }
    
    n++;
  }
  
  if (not go_on)
    throw std::runtime_error("Bad file?  Cannot continue reading!  Reached end of file or encountered bad line");
  
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
                const Numeric& fmax,
                const String& globalquantumnumbers,
                const String& localquantumnumbers,
                const Verbosity&)
{
  // Take care of quantum numbers
  String tmp_string;
  
  // Global numbers
  std::vector<QuantumNumberType> global_nums(0);
  std::istringstream global_str(globalquantumnumbers);
  while (not global_str.eof()) {
    global_str >> tmp_string; 
    global_nums.push_back(string2quantumnumbertype(tmp_string));
  }
  
  // Local numbers
  std::vector<QuantumNumberType> local_nums(0);
  std::istringstream local_str(localquantumnumbers);
  while (not local_str.eof()) {
    local_str >> tmp_string;
    local_nums.push_back(string2quantumnumbertype(tmp_string));
  }
  
  // Hitran data
  ifstream is;
  open_input_file(is, hitran_file);
  
  std::vector<Absorption::SingleLineExternal> v(0);
  
  bool go_on = true;
  while (go_on) {
    v.push_back(Absorption::ReadFromHitran2004Stream(is));
    
    if (v.back().bad) {
      v.pop_back();
      go_on = false;
    }
    else if (v.back().line.F0() > fmax) {
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
      for(auto& line: lines.AllLines())
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
  for(auto& lines: abs_lines) {
    lines.RemoveUnusedLocalQuantums();
  }
}
