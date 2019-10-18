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

void ReadLBLRTM(ArrayOfAbsorptionLines& abs_lines,
                const String& lblrtm_file,
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


void abs_linesReplaceWithLines2(ArrayOfAbsorptionLines& abs_lines, const ArrayOfAbsorptionLines& replacing_lines, const Verbosity&)
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


void abs_linesAppendWithLines2(ArrayOfAbsorptionLines& abs_lines, const ArrayOfAbsorptionLines& appending_lines, const Verbosity&)
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


void abs_linesDeleteWithLines2(ArrayOfAbsorptionLines& abs_lines, const ArrayOfAbsorptionLines& deleting_lines, const Verbosity&)
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
void abs_linesChangeBaseParameterForMatchingLines2(ArrayOfAbsorptionLines& abs_lines,
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
    abs_linesChangeBaseParameterForMatchingLines2(lines, QI, parameter_name, change, relative, loose_matching, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetBaseParameterForMatchingLines2(ArrayOfAbsorptionLines& abs_lines,
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
    abs_linesSetBaseParameterForMatchingLines2(lines, QI, parameter_name, change, loose_matching, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesChangeLineShapeModelParameterForMatchingLines2(
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
    abs_linesChangeLineShapeModelParameterForMatchingLines2(
      lines, QI, parameter, coefficient, species, x, relative, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetLineShapeModelParameterForMatchingLines2(
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
    abs_linesSetLineShapeModelParameterForMatchingLines2(
      lines, QI, parameter, coefficient, species, new_value, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void nlteSetByQuantumIdentifiers2(
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const EnergyLevelMap& nlte_field,
    const Verbosity&) {
  nlte_field.ThrowIfNotOK();
  
  const Absorption::PopulationType poptyp = nlte_field.Energies().nelem() == 0 ? 
        Absorption::PopulationType::ByNLTEPopulationDistribution :
        Absorption::PopulationType::ByNLTEVibrationalTemperatures;

  for (auto& id: nlte_field.Levels()) {
    for (auto& spec_lines: abs_lines_per_species) {
      for (auto& band: spec_lines) {
        if (band.QuantumIdentity().UpperQuantumId().In(id) or
            band.QuantumIdentity().LowerQuantumId().In(id))
          band.Population(poptyp);
      }
    }
  }
}
