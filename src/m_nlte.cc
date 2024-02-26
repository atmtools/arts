/**
 * @file m_nlte.cc
 * @author Richard Larsson
 * @date 2018-03-07
 * 
 * @brief User interface to NLTE variables and functions
 */

#include <workspace.h>

#include "atm.h"
#include "debug.h"
#include "nlte.h"
#include "xml_io.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfQuantumIdentifierFromLines(
    ArrayOfQuantumIdentifier& qid,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const Index& global) {
  // Defined only as output not input so resizing
  qid.resize(0);

  QuantumIdentifier lower, upper;

  // For all lines' upper and lower energy levels
  for (const auto& lines : abs_lines_per_species) {
    for (const auto& band : lines) {
      for (Index k = 0; k < band.NumLines() and (global ? (k == 0) : false);
           k++) {
        if (global) {
          lower = band.quantumidentity.LowerLevel();
          upper = band.quantumidentity.UpperLevel();
        } else {
          auto x = band.QuantumIdentityOfLine(k);
          lower = x.LowerLevel();
          upper = x.UpperLevel();
        }

        if (std::none_of(
                qid.begin(), qid.end(), [&](auto& x) { return x == lower; }))
          qid.push_back(lower);
        if (std::none_of(
                qid.begin(), qid.end(), [&](auto& x) { return x == upper; }))
          qid.push_back(upper);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void atmospheric_fieldRescalePopulationLevels(AtmField& atmospheric_field,
                                              const Numeric& scale) {
  for (auto& nlte : atmospheric_field.nlte()) {
    nlte.second.rescale(scale);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void collision_coefficientsFromSplitFiles(
    ArrayOfArrayOfGriddedField1& collision_coefficients,
    ArrayOfQuantumIdentifier& collision_line_identifiers,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const String& basename) {
  // Standard ARTS basename-assumption
  String tmp_basename = basename, filename;
  if (basename.length() && basename[basename.length() - 1] != '/')
    tmp_basename += ".";

  // Read the identification file
  filename = tmp_basename + "qid.xml";
  xml_read_from_file(filename, collision_line_identifiers);
  check_collision_line_identifiers(collision_line_identifiers);

  // Inner array size has to be this constantly
  const Size n = collision_line_identifiers.size();

  // Set species dimensions and fill the array
  collision_coefficients.resize(abs_species.size());
  for (Size i = 0; i < collision_coefficients.size(); i++) {
    ArrayOfGriddedField1 aogf1;

    // Read the file for a species and check that the size is correct of the array
    filename = tmp_basename +
               String(toString<1>(abs_species[i].Species())) + ".xml";
    xml_read_from_file(filename, aogf1);
    ARTS_USER_ERROR_IF(
        aogf1.size() not_eq n,
        "Mismatch between collision_line_identifiers and some collision_coefficients");
    collision_coefficients[i] = aogf1;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void nlteOff(Index& nlte_do) { nlte_do = 0; }
