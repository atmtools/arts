/* Copyright (C) 2015
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

/**
 * @file m_linerecord.cc
 * @author Richard Larsson (larsson (at) mps.mpg.de)
 * @date 2015-03-20
 * 
 * @brief Functions for altering the line catalogs
 */

#include "absorption.h"
#include "auto_md.h"
#include "sorting.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesShiftFrequency(ArrayOfLineRecord& abs_lines,
                             const Numeric& freqeuncy_shift,
                             const Verbosity&) {
  for (auto& line: abs_lines)
    line.setF(line.F() + freqeuncy_shift);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesShiftFrequency(
    ArrayOfArrayOfLineRecord& abs_lines_per_species,
    const Numeric& freqeuncy_shift,
    const Verbosity& verbosity) {
  for (auto& abs_lines: abs_lines_per_species)
    abs_linesShiftFrequency(abs_lines, freqeuncy_shift, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesRelativeLineStrengthShift(
    ArrayOfLineRecord& abs_lines,
    const Numeric& relative_line_strength_shift,
    const Verbosity&) {
  const auto r = 1.0 + relative_line_strength_shift;
  for (auto& line: abs_lines)
    line.setI0(line.I0() * r);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesRelativeLineStrengthShift(
    ArrayOfArrayOfLineRecord& abs_lines_per_species,
    const Numeric& relative_line_strength_shift,
    const Verbosity& verbosity) {
  for(auto& abs_lines: abs_lines_per_species)
    abs_linesRelativeLineStrengthShift(abs_lines, relative_line_strength_shift, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReplaceWithLines(ArrayOfLineRecord& abs_lines,
                               const ArrayOfLineRecord& replacement_lines,
                               const Verbosity&) {
  for (const auto& rline : replacement_lines) {
    auto n = 0;  // Counter for finds, API says a find must be unique

    for (auto& line : abs_lines) {
      if (line.QuantumIdentity() == rline.QuantumIdentity()) {
        line = rline;
        n++;
      }
    }

    if (n > 1) {
      std::ostringstream os;
      os << "\"" << rline.QuantumIdentity()
         << "\" is not a unique identifier wrt input catalog\n";
      throw std::runtime_error(os.str());
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesReplaceParameterWithLinesParameter(
    ArrayOfLineRecord& abs_lines,
    const ArrayOfLineRecord& replacement_lines,
    const String& parameter_name,
    const Verbosity&) {
  Index parameter_switch = -1;

  if (parameter_name.nelem() == 0)
    throw std::runtime_error("parameter_name is empty.\n");
  
  if (replacement_lines.nelem() == 0)
    throw std::runtime_error("replacement_lines is empty.\n");
  
  if (parameter_name == "Central Frequency")
    parameter_switch = 0;
  else if (parameter_name == "Line Strength")
    parameter_switch = 1;
  else if (parameter_name == "Line Shape Model")
    parameter_switch = 2;
  else if (parameter_name == "Lower State Energy")
    parameter_switch = 4;

  for (const auto& rline : replacement_lines) {
    auto n = 0;  // Counter for finds, API says a find must be unique

    for (auto& line : abs_lines) {
      if (line.QuantumIdentity() == rline.QuantumIdentity()) {
        n++;

        switch (parameter_switch) {
          case 0:  //"Central Frequency":
            line.setF(rline.F());
            break;
          case 1:  //"Line Strength":
            line.setI0(rline.I0());
            break;
          case 2:  //"Shape data":
            line.SetLineShapeModel(rline.GetLineShapeModel());
            break;
          case 4:  //"Lower State Energy":
            line.SetElow(rline.Elow());
            break;
          default: {
            ostringstream os;
            os << "Unsupported paramter_name\n"
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
void abs_linesChangeBaseParameterForMatchingLines(ArrayOfLineRecord& abs_lines,
                                                  const QuantumIdentifier& QI,
                                                  const String& parameter_name,
                                                  const Numeric& change,
                                                  const Index& relative,
                                                  const Index& loose_matching,
                                                  const Verbosity&) {
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

  for (auto& line : abs_lines) {
    if (loose_matching ? QI.In(line.QuantumIdentity())
                       : line.QuantumIdentity() == QI) {
      switch (parameter_switch) {
        case 0:  //"Central Frequency":
          if (relative == 0)
            line.setF(line.F() + change);
          else
            line.setF(line.F() * (1.0e0 + change));
          break;
        case 1:  //"Line Strength":
          if (relative == 0)
            line.setI0(line.I0() + change);
          else
            line.setI0(line.I0() * (1.0e0 + change));
          break;
        case 4:  //"Lower State Energy":
          if (relative == 0)
            line.SetElow(line.Elow() + change);
          else
            line.SetElow(line.Elow() * (1.0e0 + change));
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

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetBaseParameterForMatchingLines(ArrayOfLineRecord& abs_lines,
                                               const QuantumIdentifier& QI,
                                               const String& parameter_name,
                                               const Numeric& new_value,
                                               const Index& loose_matching,
                                               const Verbosity&) {
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

  for (auto& line : abs_lines) {
    if (loose_matching ? QI.In(line.QuantumIdentity())
                       : line.QuantumIdentity() == QI) {
      switch (parameter_switch) {
        case 0:  //"Central Frequency":
          line.setF(new_value);
          break;
        case 1:  //"Line Strength":
          line.setI0(new_value);
          break;
        case 4:  //"Lower State Energy":
          line.SetElow(new_value);
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

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesSetLineShapeModelParameterForMatchingLines(
    ArrayOfLineRecord& abs_lines,
    const QuantumIdentifier& QI,
    const String& parameter,
    const String& coefficient,
    const String& species,
    const Numeric& new_value,
    const Verbosity&) {
  bool any = false;
  for (auto& lr : abs_lines) {
    if (QI.In(lr.QuantumIdentity())) {
      lr.SetLineShapeModelParameter(new_value, species, parameter, coefficient);
      if (not any) any = true;
    }
  }

  if (not any)
    throw std::runtime_error(
        "You have no matches.  This is not accepted as a valid use case.  (Is your matching information correct?)\n");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_linesChangeLineShapeModelParameterForMatchingLines(
    ArrayOfLineRecord& abs_lines,
    const QuantumIdentifier& QI,
    const String& parameter,
    const String& coefficient,
    const String& species,
    const Numeric& change,
    const Index& relative,
    const Verbosity&) {
  bool any = false;
  for (auto& lr : abs_lines) {
    if (QI.In(lr.QuantumIdentity())) {
      auto x = lr.GetLineShapeModelParameters(species, parameter);
      auto& old = LineShape::SingleModelParameter(x, coefficient);
      if (relative)
        lr.SetLineShapeModelParameter(
            old * (1 + change), species, parameter, coefficient);
      else
        lr.SetLineShapeModelParameter(
            old + change, species, parameter, coefficient);
      if (not any) any = true;
    }
  }

  if (not any) {
    std::ostringstream os;
    os << "You have no matches.  This is not accepted as a valid use case.  (Is your matching information correct?)\n";
    os << "MATCHING INFORMATION:\t" << QI << '\n';
    throw std::runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void nlteSetByQuantumIdentifiers(
    Index& nlte_do,
    ArrayOfArrayOfLineRecord& abs_lines_per_species,
    const ArrayOfQuantumIdentifier& nlte_quantum_identifiers,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Vector& vibrational_energies,
    const String& population_type,
    const Verbosity&) {
  nlte_do = 1;

  const bool do_ev = vibrational_energies.nelem();

  if (do_ev) {
    if (vibrational_energies.nelem() != nlte_quantum_identifiers.nelem()) {
      ostringstream os;
      os << "Your vibrational energy levels vector is not the same size as\n"
         << "your *nlte_quantum_identifiers* array.  These must be the same\n"
         << "size and the content should match.\n";
      throw std::runtime_error(os.str());
    }
  }

  const LinePopulationType poptyp =
      LinePopulationTypeFromString(population_type);

  // All energies must be positive
  for (Index ii = 0; ii < vibrational_energies.nelem(); ii++) {
    if (vibrational_energies[ii] < 0) {
      ostringstream os;
      os << "Some of your vibrational energy levels are negative.  They should be positive.\n"
         << "Your vibrational levels are:\n"
         << vibrational_energies;
      throw std::runtime_error(os.str());
    }
  }

  for (Index qi = 0; qi < nlte_quantum_identifiers.nelem(); qi++) {
    auto& id = nlte_quantum_identifiers[qi];
    for (Index s = 0; s < abs_lines_per_species.nelem(); s++) {
      if (abs_species[s][0].Species() !=
          nlte_quantum_identifiers[qi].Species()) {
        continue;
      }

      ArrayOfLineRecord& species_lines = abs_lines_per_species[s];
      for (Index i = 0; i < species_lines.nelem(); i++) {
        LineRecord& lr = species_lines[i];
        if (lr.QuantumIdentity().UpperQuantumId().In(id)) {
          if (lr.NLTEUpperIndex() not_eq -1) {
            ostringstream os;
            os << "The linerecord:\n"
               << lr << "\nhad the energy state level of "
               << "this quantum identifier:\n"
               << nlte_quantum_identifiers[qi]
               << "\nset twice by the input quantum identifiers.  All levels must "
               << "point at a unique state. " << qi;
            throw std::runtime_error(os.str());
          }

          lr.SetNLTEUpperIndex(qi);
          if (do_ev) lr.SetEvupp(vibrational_energies[qi]);
          lr.SetLinePopulationType(poptyp);
        }

        if (lr.QuantumIdentity().LowerQuantumId().In(id)) {
          if (lr.NLTELowerIndex() not_eq -1) {
            ostringstream os;
            os << "The linerecord:\n"
               << lr << "\nhad the energy state level of "
               << "this quantum identifier:\n"
               << nlte_quantum_identifiers[qi]
               << "\nset twice by the input quantum identifiers.  All levels must "
               << "point at a unique state. " << qi;
            throw std::runtime_error(os.str());
          }

          lr.SetNLTELowerIndex(qi);
          if (do_ev) lr.SetEvlow(vibrational_energies[qi]);
          lr.SetLinePopulationType(poptyp);
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromabs_linesSet(Vector& f_grid,
                            const ArrayOfLineRecord& abs_lines,
                            const Numeric& half_width,
                            const Index& nr_f_per_line,
                            const Index& line_nr,
                            const Verbosity& verbosity) {
  CREATE_OUT2;

  const Index lines = abs_lines.nelem();

  if (line_nr < 0)
    f_grid.resize(lines * nr_f_per_line);
  else
    f_grid.resize(nr_f_per_line);

  out2 << "  Creating f_grid vector of length " << f_grid.nelem();

  if (lines == 0)
    throw std::runtime_error("You need at least one line to run this code.\n");
  if (line_nr >= lines)
    throw std::runtime_error(
        "You specified a line number that is outside the range of abs_lines.\n");
  if (nr_f_per_line < 1)
    throw std::runtime_error(
        "You need more than 0 frequencies per line to execute this function.\n");

  // Helper variable to ensure that there are no overlaps
  Numeric f_max = 0.0;

  if (line_nr >= 0)  // there is a line, then set frequency for a single line
  {
    if ((abs_lines[line_nr].F() - half_width) < f_max)
      throw std::runtime_error(
          "Frequencies below 0 Hz are not supported by this function.\n");
    VectorNLinSpace(f_grid,
                    nr_f_per_line,
                    abs_lines[line_nr].F() - half_width,
                    abs_lines[line_nr].F() + half_width,
                    verbosity);
  } else  // if there are many lines, then set frequency from the many lines
  {
    Vector tmp;
    for (Index ii = 0; ii < lines; ii++) {
      if ((abs_lines[ii].F() - half_width) < f_max)
        throw std::runtime_error(
            "Frequency overlaps are not supported by this function.\n");
      else
        f_max = abs_lines[ii].F() + half_width;
      VectorNLinSpace(tmp,
                      nr_f_per_line,
                      abs_lines[ii].F() - half_width,
                      abs_lines[ii].F() + half_width,
                      verbosity);
      f_grid[Range(ii * nr_f_per_line, nr_f_per_line)] = tmp;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void f_gridFromabs_lines_per_speciesSetFromSpeciesTag(
    Vector& f_grid,
    const ArrayOfArrayOfLineRecord& abs_lines_per_species,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Numeric& half_width,
    const Index& nr_f_per_line,
    const String& species_tag,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  if (species_tag == "")
    throw std::runtime_error("You need at least one tag in this code.\n");

  ArrayOfSpeciesTag st;
  array_species_tag_from_string(st, species_tag);

  for (Index ii = 0; ii < abs_species.nelem(); ii++) {
    bool test = false;
    if (abs_species[ii].nelem() == st.nelem()) {
      for (Index jj = 0; jj < st.nelem(); jj++) {
        if (st[jj] == abs_species[ii][jj])
          test = true;
        else {
          test = false;
          break;
        }
      }
      if (test) {
        f_gridFromabs_linesSet(f_grid,
                               abs_lines_per_species[ii],
                               half_width,
                               nr_f_per_line,
                               -1,
                               verbosity);
        return;
      }
    }
  }

  throw std::runtime_error("No frequency set for the given species_tag.\n");
}

void abs_linesSetQuantumNumberForAll(ArrayOfLineRecord& abs_lines,
                                     const Index& where,
                                     const String& quantum_number_name,
                                     const Rational& quantum_number_value,
                                     const Verbosity& verbosity) {
  CREATE_OUT3;
  const bool for_lower = where < 1;
  const bool for_upper = where > -1;
  for (auto& line : abs_lines) {
    if (for_lower)
      line.SetQuantumNumberLower(quantum_number_name, quantum_number_value);
    if (for_upper)
      line.SetQuantumNumberUpper(quantum_number_name, quantum_number_value);
  }

  out3 << "Set " << quantum_number_name << " to " << quantum_number_value
       << " at "
       << ((for_lower and for_upper)
               ? "both levels"
               : for_lower ? "the lower level" : "the upper level")
       << " of all " << abs_lines.nelem() << " line(s)\n";
}

void abs_linesSetNormalizationForAll(ArrayOfLineRecord& abs_lines,
                                     const String& option,
                                     const Verbosity&) {
  LineNormalizationType a;
  if (option == "VVH")
    a = LineNormalizationType::VVH;
  else if (option == "VVW")
    a = LineNormalizationType::VVW;
  else if (option == "RosenkranzQuadratic")
    a = LineNormalizationType::RosenkranzQuadratic;
  else if (option == "None")
    a = LineNormalizationType::None;
  else
    throw std::runtime_error(
        "Cannot understand normalization type option, see builtin documentation for details");

  for (LineRecord& line : abs_lines) line.SetLineNormalizationType(a);
}

void abs_linesSetMirroringForAll(ArrayOfLineRecord& abs_lines,
                                 const String& option,
                                 const Verbosity&) {
  MirroringType a;
  if (option == "Lorentz")
    a = MirroringType::Lorentz;
  else if (option == "SameAsLineShape")
    a = MirroringType::SameAsLineShape;
  else if (option == "None")
    a = MirroringType::None;
  else
    throw std::runtime_error(
        "Cannot understand mirroring type option, see builtin documentation for details");

  for (LineRecord& line : abs_lines) line.SetMirroringType(a);
}

void abs_linesCutOffForAll(ArrayOfLineRecord& abs_lines,
                           const Numeric& option,
                           const Verbosity&) {
  if (option < 0 and option not_eq -1)
    throw std::runtime_error("Cannot cutoff frequency");

  for (LineRecord& line : abs_lines) line.SetCutOff(option);
}

void abs_linesSetNlteOffForAll(ArrayOfLineRecord& abs_lines, const Verbosity&) {
  for (LineRecord& line : abs_lines)
    line.SetLinePopulationType(LinePopulationType::ByLTE);
}

void abs_lines_per_speciesSetNlteOffForAll(
    ArrayOfArrayOfLineRecord& abs_lines_per_species,
    const Verbosity& verbosity) {
  for (auto& abs_lines : abs_lines_per_species)
    abs_linesSetNlteOffForAll(abs_lines, verbosity);
}

void abs_lines_per_speciesSetNormalizationForAll(
    ArrayOfArrayOfLineRecord& abs_lines_per_species,
    const String& option,
    const Verbosity& verbosity) {
  for (auto& abs_lines : abs_lines_per_species)
    abs_linesSetNormalizationForAll(abs_lines, option, verbosity);
}

void abs_lines_per_speciesSetMirroringForAll(
    ArrayOfArrayOfLineRecord& abs_lines_per_species,
    const String& option,
    const Verbosity& verbosity) {
  for (auto& abs_lines : abs_lines_per_species)
    abs_linesSetMirroringForAll(abs_lines, option, verbosity);
}

void abs_lines_per_speciesCutOffForAll(
    ArrayOfArrayOfLineRecord& abs_lines_per_species,
    const Numeric& option,
    const Verbosity& verbosity) {
  for (auto& abs_lines : abs_lines_per_species)
    abs_linesCutOffForAll(abs_lines, option, verbosity);
}

void abs_linesFromSplitLines(
    ArrayOfLineRecord& abs_lines,
    const ArrayOfArrayOfLineRecord& abs_lines_per_species,
    const Verbosity&) {
  abs_lines.resize(0);
  for (const auto& lines : abs_lines_per_species)
    for (const auto& line : lines) abs_lines.push_back(line);
}
