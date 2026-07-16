#include "hitran_species.h"

#include <isotopologues.h>

#include <cassert>
#include <map>

//  Contributes the molparam_map
#include "auto_hitran_species_map.h"

namespace Hitran {
namespace {
using HitranMap = std::map<Index, std::map<char, std::pair<Index, Numeric>>>;

/** The latest version of the HITRAN online molparam.txt file as a map
 * 
 * Note that ARTS does not use AFGL notation as HITRAN so several species has
 * to be changed manually when recreating this map.  Please keep the comments
 * of this change around so that it is easy to do this recreation.  Thanks!
 * 
 * To recreate this map, run the pyarts.hitran.gen_latest_molparam_map.  If more
 * species names mismatch between ARTS and HITRAN, do please add a comment to indicate
 * this when you update this variable.  If any isotopologue ratio is changed, or any
 * isotopologue key (the char) is changed, also update methods.cc ReadHITRAN as it
 * depends on this map for versions of Hitran
 *
 * Last Updated: 2021-11-25
 */

using OurHitranMap = std::map<Index, std::map<char, SpeciesIsotope>>;

/** Turns the string-map required at compile time into a species-map
 * to be used as a static runtime map
 */
OurHitranMap to_species_map(const HitranMap& string_map) {
  OurHitranMap species_map;
  for (auto& specs : string_map) {
    assert(specs.second.find('1') not_eq specs.second.cend());
    for (auto& isot : specs.second) {
      assert(isot.second.first >= 0);
      species_map[specs.first][isot.first] = Species::Isotopologues[isot.second.first];
    }
  }
  return species_map;
}

/** Returns the species if possible or throws an error if it cannot be found
 * 
 * @param[in] molnum The hitran molecular number
 * @param[in] isonum The hitran character representing an isotopologue
 * @return An ARTS identifier of the species' as a transition
 */
QuantumIdentifier from_mol_iso(Index molnum, char isonum) {
  // Map of all HITRAN species
  static const OurHitranMap hitmap = to_species_map(molparam_map);

  // Search the map with pointers so that we can throw manually if something is bad

  if (auto species_list = hitmap.find(molnum); species_list not_eq hitmap.cend()) {
    if (auto species_info = species_list->second.find(isonum); species_info not_eq species_list->second.cend()) {
      return QuantumIdentifier{species_info->second};
    }
  }

  ARTS_USER_ERROR(R"(
One of two things have happened:

    1.  The molecular number {0} and isotopologue character '{1}' do not correspond to any valid HITRAN species, or
    2.  The species corresponding to the molecular number {0} and isotopologue character '{1}' has not yet been added to ARTS.

Please confirm that the molecular number {0} and isotopologue character '{1}' are correct HITRAN identifiers for the species you want to use.

If they are, ARTS must be recompiled with information to support the species.  Please see the ARTS developer documentation for instructions on
how to add a new species to ARTS.  If you do add a new species, please consider contributing the required information to the ARTS development
team so that we can add it to the next release of ARTS and avoid this issue for other users in the future.  Thanks!
)",
                  molnum,
                  isonum);
}

SpeciesIsotopologueRatios isotopologue_ratios_impl(const HitranMap& data) {
  SpeciesIsotopologueRatios out;
  for (auto& x : data) {
    for (auto& y : x.second) { out.data[y.second.first] = y.second.second; }
  }
  return out;
}
}  // namespace

QuantumIdentifier id_from_lookup(Index mol, char isochar) try {
  return from_mol_iso(mol, isochar);
} catch (const std::exception& e) {
  ARTS_USER_ERROR("Error occurred while looking up species index {} {}:\n{}", mol, isochar, e.what());
}

Numeric ratio_from_lookup(Index mol, char isochar) { return molparam_map.at(mol).at(isochar).second; }

SpeciesIsotopologueRatios isotopologue_ratios() { return isotopologue_ratios_impl(molparam_map); }
}  // namespace Hitran
