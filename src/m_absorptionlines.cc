/* Copyright (C) 2019
 Richard Larsson <ric.larsson@gmail.com>

 
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
#include "array.h"
#include "auto_md.h"
#include "debug.h"
#include "enums.h"
#include "file.h"
#include "global_data.h"
#include "m_xml.h"
#include "quantum_numbers.h"
#include <algorithm>

/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////// Basic removal and flattening
/////////////////////////////////////////////////////////////////////////////////////

void abs_linesRemoveEmptyBands(ArrayOfAbsorptionLines& abs_lines,
                               const Verbosity&)
{
  abs_lines.erase(std::remove_if(abs_lines.begin(), abs_lines.end(),
                                 [](auto& x){return x.NumLines() == 0;}),
                  abs_lines.end());
}

void abs_linesFlatten(ArrayOfAbsorptionLines& abs_lines,
                      const Verbosity& verbosity)
{
  const Index n = abs_lines.nelem();
  
  for (Index i=0; i<n; i++) {
    AbsorptionLines& band = abs_lines[i];
    if (band.NumLines()) {
      for (Index j=i+1; j<n; j++) {
        if (band.Match(abs_lines[j]).first) {
          for (auto& line: abs_lines[j].lines) {
            band.AppendSingleLine(line);
          }
          abs_lines[j].lines.clear();
        }
      }
    }
  }
  
  abs_linesRemoveEmptyBands(abs_lines, verbosity);
}

void abs_lines_per_speciesFlatten(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                  const Verbosity& verbosity)
{
  for (auto& abs_lines: abs_lines_per_species) abs_linesFlatten(abs_lines, verbosity);
}

/////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////// Reading old/external functions
/////////////////////////////////////////////////////////////////////////////////////



/** Merge an external line to abs_lines
 * 
 * @param abs_lines As WSV
 * @param sline A single local line
 * @param global_qid The band id of the local line
 */
void merge_external_line(ArrayOfAbsorptionLines& abs_lines,
                         const Absorption::SingleLineExternal& sline,
                         const QuantumIdentifier& global_qid) {
  auto band =
      std::find_if(abs_lines.begin(), abs_lines.end(), [&](const auto& li) {
        return li.MatchWithExternal(sline, global_qid);
      });
  if (band not_eq abs_lines.end()) {
    band->AppendSingleLine(sline.line);
  } else {
    abs_lines.emplace_back(sline.selfbroadening,
                           sline.bathbroadening,
                           sline.cutoff,
                           sline.mirroring,
                           sline.population,
                           sline.normalization,
                           sline.lineshapetype,
                           sline.T0,
                           sline.cutofffreq,
                           sline.linemixinglimit,
                           global_qid,
                           sline.species,
                           Array<AbsorptionSingleLine>{sline.line});
  }
}

constexpr Index merge_local_lines_size = 499;

/** Merge lines to abs_lines
 *
 * @param abs_lines As WSV
 * @param local_lines A local list of lines
 */
void merge_local_lines(ArrayOfAbsorptionLines& abs_lines,
                       const ArrayOfAbsorptionLines& local_lines) {
  for (auto& band : local_lines) {
    if (auto ptr = std::find_first_of(
            abs_lines.begin(),
            abs_lines.end(),
            &band,
            &band,
            [](const auto& abs_band, const auto& local_band) {
              return abs_band.Match(local_band).first;
            });
        ptr == abs_lines.end())
      abs_lines.push_back(band);
    else
      for (auto& line : band.lines) ptr->lines.push_back(line);
  }
}

/** Selects the global quantum numbers
 *
 * @param qns Quantum numbers to select
 * @param qid Identifeier to select from
 * @return QuantumIdentifier of all qns in qid
 */
QuantumIdentifier global_quantumidentifier(const Array<QuantumNumberType>& qns, const QuantumIdentifier& qid) {
  QuantumIdentifier out(qid.Isotopologue());
  for(auto qn: qns) {
    if (qid.val.has(qn)) {
      out.val.set(qid.val[qn]);
    }
  }
  return out;
}

/** Get a list of quantum numbers from a string
 * 
 * @param[in] qnstr A string such as "J N v1"
 * 
 * @return List of quantum numbers
 */
Array<QuantumNumberType> string2vecqn(std::string_view qnstr) {
  using namespace Quantum::Number;

  if (qnstr == "DEFAULT_GLOBAL") {
    return Array<Type>(global_types);
  }
  if (qnstr == "DEFAULT_LOCAL") {
    return Array<Type>(local_types);
  }

  const Index N = count_items(qnstr);
  Array<Type> nums(N);
  for (Index i = 0; i < N; i++) {
    nums[i] = toTypeOrThrow(items(qnstr, i));
  }

  return nums;
}

bool check_local(const Array<QuantumNumberType>& local_state) {
  using namespace Quantum::Number;
  for (auto qn: local_state) if (common_value_type(common_value_type(qn), ValueType::H) not_eq ValueType::H) return false;
  return true;
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
  const Array<QuantumNumberType> global_nums = string2vecqn(globalquantumnumbers);
  
  // Local numbers
  const Array<QuantumNumberType> local_nums = string2vecqn(localquantumnumbers);
  ARTS_USER_ERROR_IF(not check_local(local_nums), "Can only have non-string values in the local state")
  
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

  ArrayOfAbsorptionLines local_bands(0);
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
      ARTS_USER_ERROR (
        "The ARTS line file you are trying to read does not contain a valid version tag.\n"
        "Probably it was created with an older version of ARTS that used different units.")
    } else {
      istringstream is(version.substr(8));
      is >> artscat_version;
    }
    
    ARTS_USER_ERROR_IF (artscat_version < 3 or artscat_version > 5,
      "Unknown ARTS line file version: ", version)
    
    bool go_on = true;
    Index n = 0;
    while (n<nelem) {
      n++;

      if (go_on) {
        Absorption::SingleLineExternal sline;
        switch(artscat_version) {
          case 3:
            sline = Absorption::ReadFromArtscat3Stream(is_xml);
            break;
          case 4:
            sline = Absorption::ReadFromArtscat4Stream(is_xml);
            break;
          case 5:
            sline = Absorption::ReadFromArtscat5Stream(is_xml);
            break;
          default:
            ARTS_ASSERT (false, "Bad version!");
        }
        
        ARTS_USER_ERROR_IF(sline.bad, "Cannot read line ", n)
        if (sline.line.F0 < fmin) continue;
        if (sline.line.F0 > fmax) {go_on = false; continue;}

        sline.line.zeeman = Zeeman::GetAdvancedModel(sline.quantumidentity);

        // Get the global quantum number identifier
        const QuantumIdentifier global_qid = global_quantumidentifier(global_nums, sline.quantumidentity);
      
        // Get local quantum numbers into the line
        for(auto qn: local_nums) {
          if (sline.quantumidentity.val.has(qn)) {
            sline.line.localquanta.val.set(sline.quantumidentity.val[qn]);
          }
        }

        merge_external_line(local_bands, sline, global_qid);
        if (local_bands.nelem() > merge_local_lines_size) {
          merge_local_lines(abs_lines, local_bands);
          local_bands.resize(0);
        }
      } else {
        String line;
        getline(is_xml, line);
      }
    }
    
    tag.read_from_stream(is_xml);
    tag.check_name("/ArrayOfLineRecord");
  }

  merge_local_lines(abs_lines, local_bands);
  
  abs_linesSetNormalization(abs_lines, normalization_option, verbosity);
  abs_linesSetMirroring(abs_lines, mirroring_option, verbosity);
  abs_linesSetPopulation(abs_lines, population_option, verbosity);
  abs_linesSetLineShapeType(abs_lines, lineshapetype_option, verbosity);
  abs_linesSetCutoff(abs_lines, cutoff_option, cutoff_value, verbosity);
  abs_linesSetLinemixingLimit(abs_lines, linemixinglimit_value, verbosity);
  
  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
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
  const Array<QuantumNumberType> global_nums = string2vecqn(globalquantumnumbers);
  
  // Local numbers
  const Array<QuantumNumberType> local_nums = string2vecqn(localquantumnumbers);
  ARTS_USER_ERROR_IF(not check_local(local_nums), "Can only have non-string values in the local state")
  
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
    ARTS_USER_ERROR (
      "The ARTS line file you are trying to read does not contain a valid version tag.\n"
      "Probably it was created with an older version of ARTS that used different units.")
  } else {
    istringstream is(version.substr(8));
    is >> artscat_version;
  }
  
  ARTS_USER_ERROR_IF (artscat_version < 3 or artscat_version > 5,
                      "Unknown ARTS line file version: ", version)
  
  bool go_on = true;
  Index n = 0;
  ArrayOfAbsorptionLines local_bands(0);
  while (n<nelem) {
    n++;

    if (go_on) {
      Absorption::SingleLineExternal sline;
      switch(artscat_version) {
        case 3:
          sline = Absorption::ReadFromArtscat3Stream(is_xml);
          break;
        case 4:
          sline = Absorption::ReadFromArtscat4Stream(is_xml);
          break;
        case 5:
          sline = Absorption::ReadFromArtscat5Stream(is_xml);
          break;
        default:
          ARTS_ASSERT (false, "Bad version!");
      }
        
      ARTS_USER_ERROR_IF(sline.bad, "Cannot read line ", n)
      if (sline.line.F0 < fmin) continue;
      if (sline.line.F0 > fmax) {go_on = false; continue;}
    
      sline.line.zeeman = Zeeman::GetAdvancedModel(sline.quantumidentity);

      // Get the global quantum number identifier
      const QuantumIdentifier global_qid = global_quantumidentifier(global_nums, sline.quantumidentity);
    
      // Get local quantum numbers into the line
      for(auto qn: local_nums) {
        if (sline.quantumidentity.val.has(qn)) {
          sline.line.localquanta.val.set(sline.quantumidentity.val[qn]);
        }
      }

      merge_external_line(local_bands, sline, global_qid);
      if (local_bands.nelem() > merge_local_lines_size) {
        merge_local_lines(abs_lines, local_bands);
        local_bands.resize(0);
      }
    } else {
      String line;
      getline(is_xml, line);
    }
  }

  merge_local_lines(abs_lines, local_bands);
  
  abs_linesSetNormalization(abs_lines, normalization_option, verbosity);
  abs_linesSetMirroring(abs_lines, mirroring_option, verbosity);
  abs_linesSetPopulation(abs_lines, population_option, verbosity);
  abs_linesSetLineShapeType(abs_lines, lineshapetype_option, verbosity);
  abs_linesSetCutoff(abs_lines, cutoff_option, cutoff_value, verbosity);
  abs_linesSetLinemixingLimit(abs_lines, linemixinglimit_value, verbosity);
  
  tag.read_from_stream(is_xml);
  tag.check_name("/ArrayOfLineRecord");
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
  // Build a set of species indices. Duplicates are ignored.
  const std::set<Species::Species> unique_species = lbl_species(abs_species);
  
  String tmpbasename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/') {
    tmpbasename += '.';
  }
  
  // Read catalogs for each identified species and put them all into
  // abs_lines.
  abs_lines.resize(0);
  for (auto& spec: unique_species) {
    ArrayOfAbsorptionLines more_abs_lines;
    
    try {
      ReadARTSCAT(more_abs_lines,
                  tmpbasename + String(Species::toShortName(spec)) + ".xml",
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
          if (band.Match(newband).first) {
            for (Index k=0; k<newband.NumLines(); k++) {
              band.AppendSingleLine(newband.lines[k]);
              found = true;
            }
          }
        }
        if (not found) {
          abs_lines.push_back(newband);
        }
      }
    } catch (const std::exception& e) {
      ARTS_USER_ERROR_IF (not ignore_missing,
        "Errors in calls by *propmat_clearskyAddZeeman*:\n",
        e.what())
      continue;
    }
  }
  
  for (auto& band: abs_lines)
    band.sort_by_frequency();
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
  abs_lines.resize(0);
  
  // Global numbers
  const Array<QuantumNumberType> global_nums = string2vecqn(globalquantumnumbers);
  
  // Local numbers
  const Array<QuantumNumberType> local_nums = string2vecqn(localquantumnumbers);
  ARTS_USER_ERROR_IF(not check_local(local_nums), "Can only have non-string values in the local state")
  
  // HITRAN type
  const Options::HitranType hitran_version = Options::toHitranTypeOrThrow(hitran_type);
  
  // Hitran data
  ifstream is;
  open_input_file(is, hitran_file);

  ArrayOfAbsorptionLines local_bands(0);
  bool go_on = true;
  while (go_on) {

    Absorption::SingleLineExternal sline;
    switch (hitran_version) {
      case Options::HitranType::Post2004:
        sline = Absorption::ReadFromHitran2004Stream(is);
        break;
      case Options::HitranType::Pre2004:
        sline = Absorption::ReadFromHitran2001Stream(is);
        break;
      case Options::HitranType::Online:
        sline = Absorption::ReadFromHitranOnlineStream(is);
        break;
      default:
        ARTS_ASSERT(false, "The HitranType enum class has to be fully updated!\n");
    }
    
    if (sline.bad) {
      if (is.eof())
        break;
      ARTS_USER_ERROR("Cannot read line ", nelem(abs_lines) + nelem(local_bands) + 1);
    }
    if (sline.line.F0 < fmin)
      continue; // Skip this line
    if (sline.line.F0 > fmax)
      break;  // We assume sorted so quit here

    // Set Zeeman if implemented
    sline.line.zeeman = Zeeman::GetAdvancedModel(sline.quantumidentity);

    // Get the global quantum number identifier
    const QuantumIdentifier global_qid = global_quantumidentifier(global_nums, sline.quantumidentity);
    
    // Get local quantum numbers into the line
    for(auto qn: local_nums) {
      if (sline.quantumidentity.val.has(qn)) {
        sline.line.localquanta.val.set(sline.quantumidentity.val[qn]);
      }
    }

    // Either find a line like this in the list of lines or start a new Lines
    merge_external_line(local_bands, sline, global_qid);
    if (local_bands.nelem() > merge_local_lines_size) {
      merge_local_lines(abs_lines, local_bands);
      local_bands.resize(0);
    }
  }

  merge_local_lines(abs_lines, local_bands);

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
  const Array<QuantumNumberType> global_nums = string2vecqn(globalquantumnumbers);
  
  // Local numbers
  const Array<QuantumNumberType> local_nums = string2vecqn(localquantumnumbers);
  ARTS_USER_ERROR_IF(not check_local(local_nums), "Can only have non-string values in the local state")
  
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
    } else if (v.back().line.F0 < fmin) {
      v.pop_back();
    } else if (v.back().line.F0 > fmax) {
      v.pop_back();
      go_on = false;
    }
  }
  
  for (auto& x: v)
    x.line.zeeman = Zeeman::GetAdvancedModel(x.quantumidentity);
  
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
  const Array<QuantumNumberType> global_nums = string2vecqn(globalquantumnumbers);
  
  // Local numbers
  const Array<QuantumNumberType> local_nums = string2vecqn(localquantumnumbers);
  ARTS_USER_ERROR_IF(not check_local(local_nums), "Can only have non-string values in the local state")
  
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
    } else if (v.back().line.F0 < fmin) {
      v.pop_back();
    } else if (v.back().line.F0 > fmax) {
      v.pop_back();
      go_on = false;
    }
  }
  
  for (auto& x: v)
    x.line.zeeman = Zeeman::GetAdvancedModel(x.quantumidentity);
  
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
      if (thisname == specname) any = true;
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
  // Build a set of species indices. Duplicates are ignored.
  const std::set<Species::Species> unique_species = lbl_species(abs_species);
  
  String tmpbasename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/') {
    tmpbasename += '.';
  }
  
  // Read catalogs for each identified species and put them all into
  // abs_lines
  ArrayOfAbsorptionLines abs_lines(0);
  for (auto& spec: unique_species) {
    ArrayOfIsotopeRecord specs = Species::isotopologues(spec);
    for (auto& ir: specs) {
      Index i=0;
      do {
        String filename = tmpbasename + ir.FullName() + '.' + std::to_string(i) + ".xml";
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
  CREATE_OUT3;
  std::size_t bands_found{0};
  
  String tmpbasename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/') {
    tmpbasename += '.';
  }
  
  // Read catalogs for each identified species and put them all into
  // abs_lines
  abs_lines.resize(0);
  for (auto& ir: Species::Isotopologues) {
    String filename = tmpbasename + ir.FullName() + ".xml";
    if (find_xml_file_existence(filename)) {
      ArrayOfAbsorptionLines speclines;
      xml_read_from_file(filename, speclines, verbosity);
      for (auto& band: speclines) {
        abs_lines.push_back(band);
        bands_found++;
      }
    }
  }
  
  ARTS_USER_ERROR_IF (not bands_found and not robust,
                      "Cannot find any bands in the directory you are reading");
  out3 << "Found " << bands_found << " bands\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesReadSpeciesSplitCatalog(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                  const ArrayOfArrayOfSpeciesTag& abs_species,
                                                  const String& basename,
                                                  const Index& robust,
                                                  const Verbosity& verbosity)
{
  CREATE_OUT3;
  std::size_t bands_found{0};
  
  // Build a set of species indices. Duplicates are ignored.
  const std::set<Species::Species> unique_species = lbl_species(abs_species);
  
  String tmpbasename = basename;
  if (basename.length() && basename[basename.length() - 1] != '/') {
    tmpbasename += '.';
  }
  
  // Read catalogs for each identified species and put them all into
  // abs_lines
  ArrayOfAbsorptionLines abs_lines(0);
  for (auto spec: unique_species) {
    auto isots = Species::isotopologues(spec);
    for (auto& isot: isots) {
      String filename = tmpbasename + isot.FullName() + ".xml";
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
  
  ARTS_USER_ERROR_IF (not bands_found and not robust,
                      "Cannot find any bands in the directory you are reading");
  out3 << "Found " << bands_found << " bands\n";
  
  abs_lines_per_speciesCreateFromLines(abs_lines_per_species, abs_lines, abs_species, verbosity);
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
      if (tlines.Match(rlines).first) {
        number_of_matching_bands++;
        for (auto& rline: rlines.lines) {
          Index number_of_matching_single_lines = 0;
          for (auto& tline: tlines.lines) {
            if (tline.localquanta.val == rline.localquanta.val) {
              number_of_matching_single_lines++;
              tline = rline;
            }
          }
          
          ARTS_USER_ERROR_IF (number_of_matching_single_lines not_eq 1,
                              "Error: Did not match to a single single line.  "
                              "This means the input data has not been understood.  "
                              "This function needs exactly one match.");
        }
        tlines.sort_by_frequency();
      }
    }
    
    ARTS_USER_ERROR_IF (number_of_matching_bands not_eq 1,
                        "Error: Did not match to a single set of absorption lines.  "
                        "This means the input data has not been understood.  "
                        "This function needs exactly one match.");
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
        if (tlines.Match(alines).first) {
          number_of_matching_bands++;
          for (auto& aline: alines.lines) {
            Index number_of_matching_single_lines = 0;
            for (auto& tline: tlines.lines) {
              if (tline.localquanta.val == aline.localquanta.val) {
                number_of_matching_single_lines++;
              }
            }
            ARTS_USER_ERROR_IF (number_of_matching_single_lines not_eq 0,
                                "Error: Did match to a single single line.  "
                                "This means the input data has not been understood.  "
                                "This function needs exactly zero matches.");
            tlines.AppendSingleLine(aline);
          }
          tlines.sort_by_frequency();
        }
      }
      
      if (number_of_matching_bands == 0) addedlines.push_back(alines);
      ARTS_USER_ERROR_IF (number_of_matching_bands not_eq 1,
                          "Error: Did not match to a single set of absorption lines.  "
                          "This means the input data has not been understood.  "
                          "This function needs exactly one or zero matches.");
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
void abs_linesDeleteBadF0(ArrayOfAbsorptionLines& abs_lines, const Numeric& f0, const Index& lower, const Verbosity&)
{
  for (auto& lines: abs_lines) {
    std::vector<Index> hits;
    for (Index i=0; i<lines.NumLines(); i++) {
      if (lower and lines.lines[i].F0 < f0)
        hits.push_back(i);
      else if (not lower and lines.lines[i].F0 > f0)
        hits.push_back(i);
    }
    
    // Remove the bad values (sort by descending firs)
    std::sort(hits.begin(), hits.end());
    while(not hits.empty()) {
      lines.RemoveLine(hits.back());
      hits.pop_back();
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesDeleteWithLines(ArrayOfAbsorptionLines& abs_lines, const ArrayOfAbsorptionLines& deleting_lines, const Verbosity&)
{
  for (auto& dlines: deleting_lines) {
    for (auto& tlines: abs_lines) {
      std::vector<Index> hits(0);
      
      if (tlines.Match(dlines).first) {
        for (auto& dline: dlines.lines) {
          for (Index i=0; i<tlines.NumLines(); i++) {
            if (tlines.lines[i].localquanta.val == dline.localquanta.val) {
              hits.push_back(i);
            }
          }
        }
        
        // Sort and test the input
        std::sort(hits.begin(), hits.end());
        auto n = hits.size();
        hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
        ARTS_USER_ERROR_IF(n not_eq hits.size(),
                           "Removing the same line more than once is not accepted");
        
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
          if (not LineShape::modelparameterEmpty(band.lines[iline].lineshape.Data()[ispec].Data()[ivar])) {
            var_is_empty[ivar] = false;
          }
        }
      }
      
      // Remove empty variables from the writing.  This will also speed up some calculations
      for (Index iline=0; iline<band.NumLines(); iline++) {
        for (Index ivar=0; ivar < LineShape::nVars; ivar++) {
          if (var_is_empty[ivar]) {
            band.lines[iline].lineshape.Data()[ispec].Data()[ivar].type = LineShape::TemperatureModel::None;
          }
        }
      }
    }
  }
}

void abs_linesKeepBand(ArrayOfAbsorptionLines& abs_lines, const QuantumIdentifier& qid, const Verbosity&)
{
  for (auto& band: abs_lines) {
    const Quantum::Number::StateMatch lt(qid, band.quantumidentity);
    while (lt not_eq Quantum::Number::StateMatchType::Full and band.NumLines()) {
      band.RemoveLine(0);
    }
  }
}

void CheckUnique(const ArrayOfAbsorptionLines& lines,
                 const Verbosity&) {
  const Index nb = lines.nelem();
  for (Index i=0; i<nb; i++) {
    for (Index j=i+1; j<nb; j++) {
      ARTS_USER_ERROR_IF(
        Quantum::Number::StateMatch(lines[i].quantumidentity, lines[j].quantumidentity) ==
        Quantum::Number::StateMatchType::Full,
        "Not unique, these bands match:\n", lines[i], "\nand\n", lines[j])
    }
  }
}

void abs_linesReplaceBands(ArrayOfAbsorptionLines& abs_lines,
                           const ArrayOfAbsorptionLines& replacing_bands,
                           const Verbosity&) {
  for (auto& replacement: replacing_bands) {
    struct {Index band;} pos{-1};
    
    for (Index i=0; i<abs_lines.nelem(); i++) {
      auto& band = abs_lines[i];
      
      if (Quantum::Number::StateMatch(band.quantumidentity, replacement.quantumidentity) ==
        Quantum::Number::StateMatchType::Full) {
        ARTS_USER_ERROR_IF (pos.band not_eq -1, "Duplicate band matches for replacement line:\n",
                            replacement, "\nThese are for band indexes ", pos.band, " and ", i)
        
        pos.band = i;
      }
    }
    
    ARTS_USER_ERROR_IF(pos.band == -1,
      "There is no match for replacement band:\n", replacement,
      "\nYou need to append the entire band")
    abs_lines[pos.band] = replacement;
  }
}

void abs_linesReplaceLines(ArrayOfAbsorptionLines& abs_lines,
                           const ArrayOfAbsorptionLines& replacing_lines,
                           const Verbosity&) {
  const Index nb = abs_lines.nelem();  // Evaluate first so new bands are not counted
  
  for (auto& replacement: replacing_lines) {
    const Index nl=replacement.NumLines();
    
    struct {Index band; ArrayOfIndex lines;} pos{-1, ArrayOfIndex(nl, -1)};
    
    for (Index i=0; i<nb; i++) {
      auto& band = abs_lines[i];
      
      if (Quantum::Number::StateMatch(band.quantumidentity, replacement.quantumidentity) ==
        Quantum::Number::StateMatchType::Full) {
        ARTS_USER_ERROR_IF (pos.band not_eq -1, "Duplicate band matches for replacement line:\n",
                            replacement, "\nThese are for band indexes ", pos.band, " and ", i, '\n',
                            "ID ", i, ": ", band.quantumidentity, '\n',
                            "ID ", pos.band, ": ", abs_lines[pos.band].quantumidentity, '\n',
                            "ID target: ", replacement.quantumidentity, '\n')
        
        pos.band = i;
        
        for (Index k=0; k<band.NumLines(); k++) {
          for (Index j=0; j<nl; j++) {
            if (replacement.lines[j].localquanta.val == band.lines[k].localquanta.val) {
              // We cannot have multiple entries
              ARTS_USER_ERROR_IF(pos.lines[j] not_eq -1 or std::any_of(pos.lines.begin(), pos.lines.end(), [k](auto& a){return a == k;}),
                                "Found multiple matches of lines in:\n", replacement, "\n\nin mathcing band:\n", band)
              
              pos.lines[j] = k;
            }
          }
        }
      }
    }
    
    ARTS_USER_ERROR_IF(pos.band == -1 or std::any_of(pos.lines.begin(), pos.lines.end(), [](auto& a){return a == -1;}),
      "There is no match for replacement line:\n", replacement,
      "\nYou need to append the entire band")
    
    // Add or change the line catalog
    auto& band = abs_lines[pos.band];
    if (const auto [match, nullable] = band.Match(replacement); nullable) {
      for (Index j=0; j<nl; j++) band.lines[pos.lines[j]] = replacement.lines[j];
      band.MakeLineShapeModelCommon();
    } else {
      // Sort to remove from behind
      std::sort(pos.lines.begin(), pos.lines.end(), std::greater<>());
      for (auto& k: pos.lines) band.RemoveLine(k);
      abs_lines.push_back(replacement);
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
  auto t = Absorption::toCutoffTypeOrThrow(type);
  for (auto& lines: abs_lines) {
    lines.cutoff = t;
    lines.cutofffreq = x;
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
  auto t = Absorption::toCutoffTypeOrThrow(type);
  for (auto& band: abs_lines) {
    const Quantum::Number::StateMatch lt(QI, band.quantumidentity);
    if (lt == Quantum::Number::StateMatchType::Full) {
      band.cutoff = t;
      band.cutofffreq = x;
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
  auto t = Absorption::toMirroringTypeOrThrow(type);
  for (auto& lines: abs_lines)
    lines.mirroring = t;
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
  auto t = Absorption::toMirroringTypeOrThrow(type);
  for (auto& band: abs_lines) {
    const Quantum::Number::StateMatch lt(QI, band.quantumidentity);
    if (lt == Quantum::Number::StateMatchType::Full) {
      band.mirroring = t;
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

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesMakeManualMirroring(ArrayOfAbsorptionLines& abs_lines,
                                  const Verbosity&) 
{
  const ArrayOfAbsorptionLines abs_lines_copy = abs_lines;
  for (AbsorptionLines band: abs_lines_copy) {
    band.mirroring = Absorption::MirroringType::Manual;
    
    //! Don't allow running this function twice
    ARTS_USER_ERROR_IF(
      std::find_if(abs_lines_copy.cbegin(), abs_lines_copy.cend(),
                   [&band](const AbsorptionLines& li){
                     return band.Match(li).first;
                   }) not_eq abs_lines_copy.cend(),
      "Dual bands with same setup is not allowed for mirroring of band:\n",
      band, '\n')
    
    for (auto& line: band.lines) {
      line.F0 *= -1;
    }
    
    abs_lines.emplace_back(std::move(band));
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesMakeManualMirroring(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                              const Verbosity& verbosity) 
{
  for (auto& abs_lines: abs_lines_per_species) abs_linesMakeManualMirroring(abs_lines, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesMakeManualMirroringSpecies(ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                                                     const ArrayOfArrayOfSpeciesTag& abs_species,
                                                     const ArrayOfSpeciesTag& species,
                                                     const Verbosity& verbosity) 
{
  ARTS_USER_ERROR_IF(abs_species.size() not_eq abs_lines_per_species.size(),
    "Mismatch abs_species and abs_lines_per_species sizes [",
    abs_species.size(), " vs ", abs_lines_per_species.size(),
    ", respectively]")
  
  if (auto ind = std::distance(abs_species.cbegin(), std::find(abs_species.cbegin(), abs_species.cend(), species)); ind not_eq abs_species.nelem()) {
    abs_linesMakeManualMirroring(abs_lines_per_species[ind], verbosity);
  } else {
    ARTS_USER_ERROR("Cannot find species: ", species, "\nIn abs_species: [", abs_species, ']')
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
  auto t = Absorption::toPopulationTypeOrThrow(type);
  for (auto& lines: abs_lines)
    lines.population = t;
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
  auto t = Absorption::toPopulationTypeOrThrow(type);
  for (auto& lines: abs_lines) {
    const Quantum::Number::StateMatch lt(QI, lines.quantumidentity);
    if (lt == Quantum::Number::StateMatchType::Full) {
      lines.population = t;
    }
  }
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
  auto t = Absorption::toNormalizationTypeOrThrow(type);
  for (auto& lines: abs_lines)
    lines.normalization = t;
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
  auto t = Absorption::toNormalizationTypeOrThrow(type);
  for (auto& lines: abs_lines) {
    const Quantum::Number::StateMatch lt(QI, lines.quantumidentity);
    if (lt == Quantum::Number::StateMatchType::Full) {
      lines.normalization = t;
    }
  }
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
  auto t = LineShape::toType(type);
  check_enum_error(t, "Cannot understand type: ", type);
  for (auto& lines: abs_lines)
    lines.lineshapetype = t;
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
  auto t = LineShape::toTypeOrThrow(type);
  for (auto& lines: abs_lines) {
    const Quantum::Number::StateMatch lt(QI, lines.quantumidentity);
    if (lt == Quantum::Number::StateMatchType::Full) {
      lines.lineshapetype = t;
    }
  }
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
    lines.linemixinglimit = x;
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
  for (auto& lines: abs_lines) {
    const Quantum::Number::StateMatch lt(QI, lines.quantumidentity);
    if (lt == Quantum::Number::StateMatchType::Full) {
      lines.linemixinglimit = x;
    }
  }
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
    lines.T0 = x;
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
  for (auto& lines: abs_lines) {
    const Quantum::Number::StateMatch lt(QI, lines.quantumidentity);
    if (lt == Quantum::Number::StateMatchType::Full) {
      lines.T0 = x;
    }
  }
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
                                                  const Verbosity&)
{
  Index parameter_switch = -1;

  ARTS_USER_ERROR_IF (parameter_name.nelem() == 0,
                      "parameter_name is empty.\n");
  
  if (parameter_name == "Central Frequency" or
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
      const Quantum::Number::StateMatch lt(QI, band.lines[k].localquanta, band.quantumidentity);
      if (lt == Quantum::Number::StateMatchType::Full) {
        switch (parameter_switch) {
          case 0:  // "Central Frequency":
            if (relative == 0)
              band.lines[k].F0 += change;
            else
              band.lines[k].F0 *= 1.0e0 + change;
            break;
          case 1:  // "Line Strength":
            if (relative == 0)
              band.lines[k].I0 += change;
            else
              band.lines[k].I0 *= 1.0e0 + change;
            break;
          case 4:  // "Lower State Energy":
            if (relative == 0)
              band.lines[k].E0 += change;
            else
              band.lines[k].E0 *= 1.0e0 + change;
            break;
          case 5:  // "Einstein":
            if (relative == 0)
              band.lines[k].A += change;
            else
              band.lines[k].A *= 1.0e0 + change;
            break;
          case 6:  // "Lower Statistical Weight":
            if (relative == 0)
              band.lines[k].glow += change;
            else
              band.lines[k].glow *= 1.0e0 + change;
            break;
          case 7:  // "Upper Statistical Weight":
            if (relative == 0)
              band.lines[k].gupp += change;
            else
              band.lines[k].gupp *= 1.0e0 + change;
            break;
          case 8:  // "Lower Zeeman Coefficient":
            if (relative == 0)
              band.lines[k].zeeman.gl() += change;
            else
              band.lines[k].zeeman.gl() *= 1.0e0 + change;
            break;
          case 9:  // "Upper Zeeman Coefficient":
            if (relative == 0)
              band.lines[k].zeeman.gu() += change;
            else
              band.lines[k].zeeman.gu() *= 1.0e0 + change;
            break;
          default: {
            ARTS_USER_ERROR (
                                "Usupported paramter_name\n", parameter_name,
                                "\nSee method description for supported parameter names.\n")
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
                                                              const Verbosity& verbosity)
{
  for (auto& lines: abs_lines_per_species)
    abs_linesChangeBaseParameterForMatchingLines(lines, QI, parameter_name, change, relative, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesChangeBaseParameterForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const QuantumIdentifier& QI,
  const String& parameter_name,
  const Numeric& change,
  const Index& relative,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesChangeBaseParameterForMatchingLines(abs_lines_per_species[ispec], QI, parameter_name, change, relative, verbosity);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetBaseParameterForMatchingLines(ArrayOfAbsorptionLines& abs_lines,
                                               const QuantumIdentifier& QI,
                                               const String& parameter_name,
                                               const Numeric& x,
                                               const Verbosity&)
{
  Index parameter_switch = -1;
  
  ARTS_USER_ERROR_IF (parameter_name.nelem() == 0,
                      "parameter_name is empty.\n");
  
  if (parameter_name == "Central Frequency" or
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
      const Quantum::Number::StateMatch lt(QI, band.lines[k].localquanta, band.quantumidentity);
      if (lt == Quantum::Number::StateMatchType::Full) {
        switch (parameter_switch) {
          case 0:  // "Central Frequency":
            band.lines[k].F0 = x;
            break;
          case 1:  // "Line Strength":
            band.lines[k].I0 = x;
            break;
          case 4:  // "Lower State Energy":
            band.lines[k].E0 = x;
            break;
          case 5:  // "Einstein":
            band.lines[k].A = x;
            break;
          case 6:  // "Lower Statistical Weight":
            band.lines[k].glow = x;
            break;
          case 7:  // "Upper Statistical Weight":
            band.lines[k].gupp = x;
            break;
          case 8:
            band.lines[k].zeeman.gl() = x;
            break;
          case 9:
            band.lines[k].zeeman.gu() = x;
            break;
          default: {
            ARTS_USER_ERROR (
                                "Usupported paramter_name\n", parameter_name,
                                "\nSee method description for supported parameter names.\n")
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
                                                           const Verbosity& verbosity)
{
  for (auto& lines: abs_lines_per_species)
    abs_linesSetBaseParameterForMatchingLines(lines, QI, parameter_name, change, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetBaseParameterForSpecies(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const ArrayOfArrayOfSpeciesTag& abs_species,
  const QuantumIdentifier& QI,
  const String& parameter_name,
  const Numeric& change,
  const String& species_tag,
  const Verbosity& verbosity)
{
  Index t1, t2;
  ArrayOfArrayOfSpeciesTag target_species;
  abs_speciesSet(target_species, t1, t2, {species_tag}, verbosity);
  for (Index ispec=0; ispec<abs_species.nelem(); ispec++) {
    if (std::equal(abs_species[ispec].begin(), abs_species[ispec].end(), target_species[0].begin())) {
      abs_linesSetBaseParameterForMatchingLines(abs_lines_per_species[ispec], QI, parameter_name, change, verbosity);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetLineShapeModelParametersForMatchingLines(
  ArrayOfAbsorptionLines& abs_lines,
  const QuantumIdentifier& QI,
  const String& parameter,
  const String& species,
  const String& temperaturemodel,
  const Vector& new_values,
  const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  const bool do_self = species == LineShape::self_broadening;
  const bool do_bath = species == LineShape::bath_broadening;
  
  // Set the spec index if possible
  const Species::Species spec = do_self ? Species::Species::FINAL : do_bath ? Species::Species::Bath : SpeciesTag(species).Spec();
  
  const LineShape::Variable var = LineShape::toVariableOrThrow(parameter);
  
  ARTS_USER_ERROR_IF (new_values.nelem() not_eq LineShape::ModelParameters::N,
    "Mismatch between input and expected number of variables\n"
    "\tInput is: ", new_values.nelem(), " long but expects: ", LineShape::ModelParameters::N, " values\n")
  
  LineShape::ModelParameters newdata;
  newdata.type = LineShape::toTemperatureModelOrThrow(temperaturemodel);
  newdata.X0 = new_values[0];
  newdata.X1 = new_values[1];
  newdata.X2 = new_values[2];
  newdata.X3 = new_values[3];
  
  for (auto& band: abs_lines) {
    for (Index k=0; k<band.NumLines(); k++) {
      const Quantum::Number::StateMatch lt(QI, band.lines[k].localquanta, band.quantumidentity);
      if (lt == Quantum::Number::StateMatchType::Full) {
        if (do_self and band.selfbroadening) {
          out3 << "Changing self\n";
          band.lines[k].lineshape.Data().front().Data()[Index(var)] = newdata;
        } else if (do_bath and band.bathbroadening) {
          out3 << "Changing bath\n";
          band.lines[k].lineshape.Data().back().Data()[Index(var)] = newdata;
        } else {
          for (Index i=band.selfbroadening; i<band.broadeningspecies.nelem()-band.bathbroadening; i++) {
            if (spec == band.broadeningspecies[i]) {
              out3 << "Changing species: " << Species::toShortName(spec) << '\n';
              band.lines[k].lineshape.Data()[i].Data()[Index(var)] = newdata;
            }
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetLineShapeModelParametersForMatchingLines(
  ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
  const QuantumIdentifier& QI,
  const String& parameter,
  const String& species,
  const String& temperaturemodel,
  const Vector& new_values,
  const Verbosity& verbosity)
{
  for (auto& lines: abs_lines_per_species)
    abs_linesSetLineShapeModelParametersForMatchingLines(
      lines, QI, parameter, species, temperaturemodel, new_values, verbosity);
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
  Index parameter_switch = -1;
  
  ARTS_USER_ERROR_IF (parameter_name.nelem() == 0,
                      "parameter_name is empty.\n");
  if (parameter_name == "Statistical Weight")
    parameter_switch = 1;
  else if (parameter_name == "Zeeman Coefficient")
    parameter_switch = 2;
  
  for (auto& band: abs_lines) {
    for (Index k=0; k<band.NumLines(); k++) {
      const Quantum::Number::StateMatch lt(QI, band.lines[k].localquanta, band.quantumidentity);
      if (lt == Quantum::Number::StateMatchType::Level and lt.low) {
        switch (parameter_switch) {
          case 1:  // "Statistical Weight":
            if (relative == 0)
              band.lines[k].glow += change;
            else
              band.lines[k].glow *= 1.0e0 + change;
            break;
          case 2:  // "Zeeman Coefficient":
            if (relative == 0)
              band.lines[k].zeeman.gl() += change;
            else
              band.lines[k].zeeman.gl() *= 1.0e0 + change;
            break;
          default: {
            ARTS_USER_ERROR (
                                "Usupported paramter_name\n", parameter_name,
                                "\nSee method description for supported parameter names.\n")
          }
        }
      }
      
      if (lt == Quantum::Number::StateMatchType::Level and lt.upp) {
        switch (parameter_switch) {
          case 1:  // "Statistical Weight":
            if (relative == 0)
              band.lines[k].gupp += change;
            else
              band.lines[k].gupp *= 1.0e0 + change;
            break;
          case 2:  // "Zeeman Coefficient":
            if (relative == 0)
              band.lines[k].zeeman.gu() += change;
            else
              band.lines[k].zeeman.gu() *= 1.0e0 + change;
            break;
          default: {
            ARTS_USER_ERROR (
                                "Usupported paramter_name\n", parameter_name,
                                "\nSee method description for supported parameter names.\n")
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
  ARTS_USER_ERROR_IF (QID.nelem() not_eq change.nelem(),
                      "Mismatch between QID and change input lengths not allowed");
  
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
  ARTS_USER_ERROR_IF (QID.nelem() not_eq change.nelem(),
                      "Mismatch between QID and change input lengths not allowed");
  
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
  Index parameter_switch = -1;
  
  ARTS_USER_ERROR_IF (parameter_name.nelem() == 0,
                      "parameter_name is empty.\n");
  if (parameter_name == "Statistical Weight")
    parameter_switch = 1;
  else if (parameter_name == "Zeeman Coefficient")
    parameter_switch = 2;
  
  for (auto& band: abs_lines) {
    for (Index k=0; k<band.NumLines(); k++) {
      const Quantum::Number::StateMatch lt(QI, band.lines[k].localquanta, band.quantumidentity);
      if (lt == Quantum::Number::StateMatchType::Level and lt.low) {
        switch (parameter_switch) {
          case 1:  // "Statistical Weight":
            band.lines[k].glow = x;
            break;
          case 2:  // "Zeeman Coefficient":
            band.lines[k].zeeman.gl() = x;
            break;
          default: {
            ARTS_USER_ERROR (
                                "Usupported paramter_name\n", parameter_name,
                                "\nSee method description for supported parameter names.\n")
          }
        }
      }
      
      if (lt == Quantum::Number::StateMatchType::Level and lt.upp) {
        switch (parameter_switch) {
          case 1:  // "Statistical Weight":
            band.lines[k].gupp = x;
            break;
          case 2:  // "Zeeman Coefficient":
            band.lines[k].zeeman.gu() = x;
            break;
          default: {
            ARTS_USER_ERROR (
                                "Usupported paramter_name\n", parameter_name,
                                "\nSee method description for supported parameter names.\n")
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
  ARTS_USER_ERROR_IF (QID.nelem() not_eq change.nelem(),
                      "Mismatch between QID and change input lengths not allowed");
  
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
  ARTS_USER_ERROR_IF (QID.nelem() not_eq change.nelem(),
                      "Mismatch between QID and change input lengths not allowed");
  
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
    const Verbosity&) {
  nlte_field.ThrowIfNotOK();

  if (nlte_field.Data().empty()) {
    nlte_do = 0;
    return;
  }
  nlte_do = 1;

  const Absorption::PopulationType poptyp =
      nlte_field.Energies().empty() ? Absorption::PopulationType::NLTE
                                    : Absorption::PopulationType::VibTemps;
  
  for (auto& spec_lines : abs_lines_per_species) {
    for (auto& band : spec_lines) {
      Index low = 0, upp = 0;
      for (auto& id : nlte_field.Levels()) {
        for (auto& line : band.lines) {
          const auto lt =
              poptyp == Absorption::PopulationType::NLTE
                  ? Quantum::Number::StateMatch(
                        id, line.localquanta, band.quantumidentity)
                  : Quantum::Number::StateMatch(id, band.quantumidentity);
          low += lt.low;
          upp += lt.upp;
        }
      }

      ARTS_USER_ERROR_IF(
          not(low == 0 or low == band.NumLines()) or
              not(upp == 0 or upp == band.NumLines()),
          "The band with metadata:\n",
          band.MetaData(),
          "\n\nwith lines:\n",
          band.lines,
          "\n\nThe number of levels don't add upp correctly.\n There are ",
          band.NumLines(),
          " lines but there are ",
          low,
          " lower level matches and ",
          upp,
          " upper level matches.")

      if (upp or low) band.population = poptyp;
    }
  }
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Destructively manipulate content of whole catalog
/////////////////////////////////////////////////////////////////////////////////////

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesSetEmpty(
    ArrayOfArrayOfAbsorptionLines & abs_lines_per_species,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Verbosity&) {
  abs_lines_per_species = ArrayOfArrayOfAbsorptionLines(
      abs_species.nelem(), ArrayOfAbsorptionLines(0));
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesCompact(ArrayOfAbsorptionLines & abs_lines,
                      const Vector& f_grid,
                      const Verbosity&) {
  const Numeric fmax = max(f_grid);
  const Numeric fmin = min(f_grid);

  for (auto& band : abs_lines) {
    for (Index k = band.NumLines() - 1; k >= 0; k--) {
      const Numeric fcut_upp = band.CutoffFreq(k);
      const Numeric fcut_low = band.CutoffFreqMinus(k);

      if (fmax < fcut_low or fmin > fcut_upp) {
        band.RemoveLine(k);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesCompact(
    ArrayOfArrayOfAbsorptionLines & abs_lines_per_species,
    const Vector& f_grid,
    const Verbosity& verbosity) {
  for (auto& lines : abs_lines_per_species) {
    abs_linesCompact(lines, f_grid, verbosity);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesRemoveBand(ArrayOfAbsorptionLines & abs_lines,
                         const QuantumIdentifier& qid,
                         const Verbosity&) {
  for (Index i = 0; i < abs_lines.nelem(); i++) {
    const Quantum::Number::StateMatch lt(qid, abs_lines[i].quantumidentity);
    if (lt == Quantum::Number::StateMatchType::Full) {
      abs_lines.erase(abs_lines.begin() + i);
      break;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////// Manipulation of other ARTS variables based on AbsorptionLines
/////////////////////////////////////////////////////////////////////////////////////

template <class T>
std::vector<T> linspace(
    T s, T e, typename std::vector<T>::size_type count) noexcept {
  std::vector<T> ls(count);

  if (count == 0) return ls;
  if (count == 1) {
    ls.front() = (e + s) / 2;
    return ls;
  }
  const T step = (e - s) / T(count - 1);
  ls.front() = s;
  ls.back() = e;
  for (typename std::vector<T>::size_type i = 1; i < count - 1; ++i)
    ls[i] = s + step * T(i);
  return ls;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromAbsorptionLines(
    Vector & f_grid,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const Numeric& delta_f_low,
    const Numeric& delta_f_upp,
    const Index& num_freqs,
    const Verbosity&) {
  const Index n = nelem(abs_lines_per_species);

  ARTS_USER_ERROR_IF(delta_f_low >= delta_f_upp,
                     "The lower frequency delta has to be smaller "
                     "than the upper frequency delta");
  ARTS_USER_ERROR_IF(num_freqs == 0, "Need more than zero frequency points");
  ARTS_USER_ERROR_IF(n < 1,
                     "No lines found.  Error?  Use *VectorSet* "
                     "to resize *f_grid*");

  std::vector<Numeric> fout(0);
  for (auto& lines : abs_lines_per_species) {
    for (auto& band : lines) {
      for (Index k = 0; k < band.NumLines(); k++) {
        if (num_freqs > 1) {
          auto ftmp =
              linspace<Numeric>(band.lines[k].F0 + delta_f_low,
                                band.lines[k].F0 + delta_f_upp,
                                std::vector<Numeric>::size_type(num_freqs));
          for (auto& f : ftmp) {
            if (f > 0) fout.push_back(f);
          }
        } else {
          fout.push_back(band.lines[k].F0);
        }
      }
    }
  }

  std::sort(fout.begin(), fout.end());
  fout.erase(std::unique(fout.begin(), fout.end()), fout.end());
  f_grid.resize(fout.size());
  for (Index i = 0; i < f_grid.nelem(); i++) f_grid[i] = fout[i];
}

/////////////////////////////////////////////////////////
//////////////////// Remove lines safely from the catalog
/////////////////////////////////////////////////////////

void remove_impl(ArrayOfAbsorptionLines & abs_lines,
                 const ArrayOfSpeciesTag& species,
                 const Numeric lower_frequency,
                 const Numeric upper_frequency,
                 const Numeric lower_intensity,
                 const Index safe,
                 const Index flip_flims,                 
                 const Verbosity& verbosity) {
  const bool care_about_species = species.nelem();

  for (auto& band : abs_lines) {
    if (care_about_species and
        species[0].Isotopologue() not_eq band.Isotopologue())
      continue;

    auto& lines = band.lines;

    if (not safe) {
      std::vector<std::size_t> rem;
      if (flip_flims) {
        for (std::size_t i=lines.size()-1; i<lines.size(); i--) {
          auto& line = lines[i];
          if ((line.F0 >= lower_frequency and
               line.F0 <= upper_frequency) or
              line.I0 < lower_intensity)
            rem.push_back(i);
        }
      } else {      
        for (std::size_t i = lines.size() - 1; i < lines.size(); i--) {
          auto& line = lines[i];
          if (line.F0 < lower_frequency or line.F0 > upper_frequency or
              line.I0 < lower_intensity)
            rem.push_back(i);
        }
      }
      for (auto i : rem) band.RemoveLine(i);
    } else {
      ARTS_USER_ERROR_IF(flip_flims,
                         "Not allowed to combine GINs flip_flims and safe.")      
      const bool all_low = std::all_of(
          lines.begin(), lines.end(), [lower_frequency](auto& line) {
            return line.F0 < lower_frequency;
          });
      const bool all_upp = std::all_of(
          lines.begin(), lines.end(), [upper_frequency](auto& line) {
            return line.F0 > upper_frequency;
          });
      const bool low_int = std::all_of(
          lines.begin(), lines.end(), [lower_intensity](auto& line) {
            return line.I0 < lower_intensity;
          });
      if (all_low or all_upp or low_int) lines.resize(0);
    }
  }

  // Removes empty bands
  abs_linesRemoveEmptyBands(abs_lines, verbosity);
}

void abs_linesRemoveLines(ArrayOfAbsorptionLines & abs_lines,
                          const Numeric& lower_frequency,
                          const Numeric& upper_frequency,
                          const Numeric& lower_intensity,
                          const Index& safe,
                          const Index& flip_flims,
                          const Verbosity& verbosity) {
  remove_impl(abs_lines,
              {},
              lower_frequency,
              upper_frequency,
              lower_intensity,
              safe,
              flip_flims,
              verbosity);
}

void abs_lines_per_speciesRemoveLines(
    ArrayOfArrayOfAbsorptionLines & abs_lines_per_species,
    const Numeric& lower_frequency,
    const Numeric& upper_frequency,
    const Numeric& lower_intensity,
    const Index& safe,
    const Index& flip_flims,
    const Verbosity& verbosity) {
  for (auto& abs_lines : abs_lines_per_species)
    abs_linesRemoveLines(abs_lines,
                         lower_frequency,
                         upper_frequency,
                         lower_intensity,
                         safe,
                         flip_flims,
                         verbosity);
}

void abs_linesRemoveLinesFromSpecies(ArrayOfAbsorptionLines & abs_lines,
                                     const ArrayOfSpeciesTag& species,
                                     const Numeric& lower_frequency,
                                     const Numeric& upper_frequency,
                                     const Numeric& lower_intensity,
                                     const Index& safe,
                                     const Index& flip_flims,
                                     const Verbosity& verbosity) {
  ARTS_USER_ERROR_IF(
      species.nelem() not_eq 1, "Must have a single species, got: ", species)
  ARTS_USER_ERROR_IF(species[0].Isotopologue().joker(),
                     "Cannot give joker species, got: ",
                     species)

  remove_impl(abs_lines,
              species,
              lower_frequency,
              upper_frequency,
              lower_intensity,
              safe,
              flip_flims,
              verbosity);
}

void abs_lines_per_speciesRemoveLinesFromSpecies(
    ArrayOfArrayOfAbsorptionLines & abs_lines_per_species,
    const ArrayOfSpeciesTag& species,
    const Numeric& lower_frequency,
    const Numeric& upper_frequency,
    const Numeric& lower_intensity,
    const Index& safe,
    const Index& flip_flims,
    const Verbosity& verbosity) {
  for (auto& abs_lines : abs_lines_per_species)
    abs_linesRemoveLinesFromSpecies(abs_lines,
                                    species,
                                    lower_frequency,
                                    upper_frequency,
                                    lower_intensity,
                                    safe,
                                    flip_flims,
                                    verbosity);
}

void abs_linesSort(ArrayOfAbsorptionLines & abs_lines,
                   const String& option,
                   const Verbosity& verbosity) {
  const auto opt = Options::toSortingOptionOrThrow(option);

  abs_linesRemoveEmptyBands(abs_lines, verbosity);

  switch (opt) {
    case Options::SortingOption::ByFrequency:
      for (auto& band : abs_lines) band.sort_by_frequency();
      std::sort(abs_lines.begin(), abs_lines.end(), [](auto& a, auto& b) {
        return a.lines[0].F0 <= b.lines[0].F0;
      });
      break;
    case Options::SortingOption::ByEinstein:
      for (auto& band : abs_lines) band.sort_by_einstein();
      std::sort(abs_lines.begin(), abs_lines.end(), [](auto& a, auto& b) {
        return a.lines[0].A <= b.lines[0].A;
      });
      break;
    case Options::SortingOption::FINAL: { /*leave last*/
    }
  }
}
