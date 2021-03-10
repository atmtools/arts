#include <map>

#include "hitran_species.h"

#include "abs_species_tags.h"

namespace Hitran {
/** In 2012 the order if isotopologues were changed in HITRAN
 * 
 * This version takes that into account.  Note that the other species
 * are left as they were in the latest_molparam_map map at the time
 * of creating this file (2021-03-10)
 */
const std::map<Index, std::map<char, const char * const>> pre2012co2change_molparam_map{
  {1, {  // H2O
    {'1', "H2O-161"},
    {'2', "H2O-181"},
    {'3', "H2O-171"},
    {'4', "H2O-162"},
    {'5', "H2O-182"},
    {'6', "H2O-172"},
    {'7', "H2O-262"},
  }},
  {2, {  // CO2
    {'1', "CO2-626"},
    {'2', "CO2-636"},
    {'3', "CO2-628"},
    {'4', "CO2-627"},
    {'5', "CO2-638"},
    {'6', "CO2-637"},
    {'7', "CO2-828"},
    {'8', "CO2-728"},  // Modified from HITRAN because ARTS does not use AFGL notation
    {'9', "CO2-838"},  // This is different for HITRAN2008 and earlier cf original map
  }},
  {3, {  // O3
    {'1', "O3-666"},
    {'2', "O3-668"},
    {'3', "O3-686"},
    {'4', "O3-667"},
    {'5', "O3-676"},
  }},
  {4, {  // N2O
    {'1', "N2O-446"},
    {'2', "N2O-456"},
    {'3', "N2O-546"},
    {'4', "N2O-448"},
    {'5', "N2O-447"},
  }},
  {5, {  // CO
    {'1', "CO-26"},
    {'2', "CO-36"},
    {'3', "CO-28"},
    {'4', "CO-27"},
    {'5', "CO-38"},
    {'6', "CO-37"},
  }},
  {6, {  // CH4
    {'1', "CH4-211"},
    {'2', "CH4-311"},
    {'3', "CH4-212"},
    {'4', "CH4-312"},
  }},
  {7, {  // O2
    {'1', "O2-66"},
    {'2', "O2-68"},
    {'3', "O2-67"},
  }},
  {8, {  // NO
    {'1', "NO-46"},
    {'2', "NO-56"},
    {'3', "NO-48"},
  }},
  {9, {  // SO2
    {'1', "SO2-626"},
    {'2', "SO2-646"},
  }},
  {10, {  // NO2
    {'1', "NO2-646"},
    {'2', "NO2-656"},
  }},
  {11, {  // NH3
    {'1', "NH3-4111"},
    {'2', "NH3-5111"},
  }},
  {12, {  // HNO3
    {'1', "HNO3-146"},
    {'2', "HNO3-156"},
  }},
  {13, {  // OH
    {'1', "OH-61"},
    {'2', "OH-81"},
    {'3', "OH-62"},
  }},
  {14, {  // HF
    {'1', "HF-19"},
    {'2', "HF-29"},
  }},
  {15, {  // HCl
    {'1', "HCl-15"},
    {'2', "HCl-17"},
    {'3', "HCl-25"},
    {'4', "HCl-27"},
  }},
  {16, {  // HBr
    {'1', "HBr-19"},
    {'2', "HBr-11"},
    {'3', "HBr-29"},
    {'4', "HBr-21"},
  }},
  {17, {  // HI
    {'1', "HI-17"},
    {'2', "HI-27"},
  }},
  {18, {  // ClO
    {'1', "ClO-56"},
    {'2', "ClO-76"},
  }},
  {19, {  // OCS
    {'1', "OCS-622"},
    {'2', "OCS-624"},
    {'3', "OCS-632"},
    {'4', "OCS-623"},
    {'5', "OCS-822"},
    {'6', "OCS-634"},
  }},
  {20, {  // H2CO
    {'1', "H2CO-1126"},  // Modified from HITRAN because ARTS does not use AFGL notation
    {'2', "H2CO-1136"},  // Modified from HITRAN because ARTS does not use AFGL notation
    {'3', "H2CO-1128"},  // Modified from HITRAN because ARTS does not use AFGL notation
  }},
  {21, {  // HOCl
    {'1', "HOCl-165"},
    {'2', "HOCl-167"},
  }},
  {22, {  // N2
    {'1', "N2-44"},
    {'2', "N2-45"},
  }},
  {23, {  // HCN
    {'1', "HCN-124"},
    {'2', "HCN-134"},
    {'3', "HCN-125"},
  }},
  {24, {  // CH3Cl
    {'1', "CH3Cl-215"},
    {'2', "CH3Cl-217"},
  }},
  {25, {  // H2O2
    {'1', "H2O2-1661"},
  }},
  {26, {  // C2H2
    {'1', "C2H2-1221"},
    {'2', "C2H2-1231"},
    {'3', "C2H2-1222"},
  }},
  {27, {  // C2H6
    {'1', "C2H6-1221"},
    {'2', "C2H6-1231"},
  }},
  {28, {  // PH3
    {'1', "PH3-1111"},
  }},
  {29, {  // COF2
    {'1', "COF2-269"},
    {'2', "COF2-369"},
  }},
  {30, {  // SF6
    {'1', "SF6-29"},
  }},
  {31, {  // H2S
    {'1', "H2S-121"},
    {'2', "H2S-141"},
    {'3', "H2S-131"},
  }},
  {32, {  // HCOOH
    {'1', "HCOOH-1261"},  // Modified from HITRAN because ARTS does not use AFGL notation
  }},
  {33, {  // HO2
    {'1', "HO2-166"},
  }},
  {34, {  // O
    {'1', "O-6"},
  }},
  {35, {  // ClONO2
    {'1', "ClONO2-5646"},
    {'2', "ClONO2-7646"},
  }},
  {36, {  // NO+
    {'1', "NO+-46"},
  }},
  {37, {  // HOBr
    {'1', "HOBr-169"},
    {'2', "HOBr-161"},
  }},
  {38, {  // C2H4
    {'1', "C2H4-221"},
    {'2', "C2H4-231"},
  }},
  {39, {  // CH3OH
    {'1', "CH3OH-2161"},
  }},
  {40, {  // CH3Br
    {'1', "CH3Br-219"},
    {'2', "CH3Br-211"},
  }},
  {41, {  // CH3CN
    {'1', "CH3CN-211124"},  // Modified from HITRAN because ARTS does not use AFGL notation
  }},
  {42, {  // CF4
    {'1', "CF4-29"},
  }},
  {43, {  // C4H2
    {'1', "C4H2-2211"},
  }},
  {44, {  // HC3N
    {'1', "HC3N-12224"},  // Modified from HITRAN because ARTS does not use AFGL notation
  }},
  {45, {  // H2
    {'1', "H2-11"},
    {'2', "H2-12"},
  }},
  {46, {  // CS
    {'1', "CS-22"},
    {'2', "CS-24"},
    {'3', "CS-32"},
    {'4', "CS-23"},
  }},
  {47, {  // SO3
    {'1', "SO3-26"},
  }},
  {48, {  // C2N2
    {'1', "C2N2-4224"},
  }},
  {49, {  // COCl2
    {'1', "COCl2-2655"},
    {'2', "COCl2-2657"},
  }},
  {53, {  // CS2
    {'1', "CS2-222"},
    {'2', "CS2-224"},
    {'3', "CS2-223"},
    {'4', "CS2-232"},
  }},
};

/** The latest version of the HITRAN online molparam.txt file as a map
 * 
 * Note that ARTS does not use AFGL notation as HITRAN so several species has
 * to be changed manually when recreating this map.  Please keep the comments
 * of this change around so that it is easy to do this recreation.  Thanks!
 * 
 * To recreate this map, run the pyarts.hitran.gen_latest_molparam_map.  If more
 * species names mismatch between ARTS and HITRAN, do please add a comment to indicate
 * this when you update this variable
 */
const std::map<Index, std::map<char, const char * const>> latest_molparam_map{
  {1, {  // H2O
    {'1', "H2O-161"},
    {'2', "H2O-181"},
    {'3', "H2O-171"},
    {'4', "H2O-162"},
    {'5', "H2O-182"},
    {'6', "H2O-172"},
    {'7', "H2O-262"},
  }},
  {2, {  // CO2
    {'1', "CO2-626"},
    {'2', "CO2-636"},
    {'3', "CO2-628"},
    {'4', "CO2-627"},
    {'5', "CO2-638"},
    {'6', "CO2-637"},
    {'7', "CO2-828"},
    {'8', "CO2-728"},  // Modified from HITRAN because ARTS does not use AFGL notation
    {'9', "CO2-727"},
    {'0', "CO2-838"},
    {'A', "CO2-837"},
    {'B', "CO2-737"},
  }},
  {3, {  // O3
    {'1', "O3-666"},
    {'2', "O3-668"},
    {'3', "O3-686"},
    {'4', "O3-667"},
    {'5', "O3-676"},
  }},
  {4, {  // N2O
    {'1', "N2O-446"},
    {'2', "N2O-456"},
    {'3', "N2O-546"},
    {'4', "N2O-448"},
    {'5', "N2O-447"},
  }},
  {5, {  // CO
    {'1', "CO-26"},
    {'2', "CO-36"},
    {'3', "CO-28"},
    {'4', "CO-27"},
    {'5', "CO-38"},
    {'6', "CO-37"},
  }},
  {6, {  // CH4
    {'1', "CH4-211"},
    {'2', "CH4-311"},
    {'3', "CH4-212"},
    {'4', "CH4-312"},
  }},
  {7, {  // O2
    {'1', "O2-66"},
    {'2', "O2-68"},
    {'3', "O2-67"},
  }},
  {8, {  // NO
    {'1', "NO-46"},
    {'2', "NO-56"},
    {'3', "NO-48"},
  }},
  {9, {  // SO2
    {'1', "SO2-626"},
    {'2', "SO2-646"},
  }},
  {10, {  // NO2
    {'1', "NO2-646"},
    {'2', "NO2-656"},
  }},
  {11, {  // NH3
    {'1', "NH3-4111"},
    {'2', "NH3-5111"},
  }},
  {12, {  // HNO3
    {'1', "HNO3-146"},
    {'2', "HNO3-156"},
  }},
  {13, {  // OH
    {'1', "OH-61"},
    {'2', "OH-81"},
    {'3', "OH-62"},
  }},
  {14, {  // HF
    {'1', "HF-19"},
    {'2', "HF-29"},
  }},
  {15, {  // HCl
    {'1', "HCl-15"},
    {'2', "HCl-17"},
    {'3', "HCl-25"},
    {'4', "HCl-27"},
  }},
  {16, {  // HBr
    {'1', "HBr-19"},
    {'2', "HBr-11"},
    {'3', "HBr-29"},
    {'4', "HBr-21"},
  }},
  {17, {  // HI
    {'1', "HI-17"},
    {'2', "HI-27"},
  }},
  {18, {  // ClO
    {'1', "ClO-56"},
    {'2', "ClO-76"},
  }},
  {19, {  // OCS
    {'1', "OCS-622"},
    {'2', "OCS-624"},
    {'3', "OCS-632"},
    {'4', "OCS-623"},
    {'5', "OCS-822"},
    {'6', "OCS-634"},
  }},
  {20, {  // H2CO
    {'1', "H2CO-1126"},  // Modified from HITRAN because ARTS does not use AFGL notation
    {'2', "H2CO-1136"},  // Modified from HITRAN because ARTS does not use AFGL notation
    {'3', "H2CO-1128"},  // Modified from HITRAN because ARTS does not use AFGL notation
  }},
  {21, {  // HOCl
    {'1', "HOCl-165"},
    {'2', "HOCl-167"},
  }},
  {22, {  // N2
    {'1', "N2-44"},
    {'2', "N2-45"},
  }},
  {23, {  // HCN
    {'1', "HCN-124"},
    {'2', "HCN-134"},
    {'3', "HCN-125"},
  }},
  {24, {  // CH3Cl
    {'1', "CH3Cl-215"},
    {'2', "CH3Cl-217"},
  }},
  {25, {  // H2O2
    {'1', "H2O2-1661"},
  }},
  {26, {  // C2H2
    {'1', "C2H2-1221"},
    {'2', "C2H2-1231"},
    {'3', "C2H2-1222"},
  }},
  {27, {  // C2H6
    {'1', "C2H6-1221"},
    {'2', "C2H6-1231"},
  }},
  {28, {  // PH3
    {'1', "PH3-1111"},
  }},
  {29, {  // COF2
    {'1', "COF2-269"},
    {'2', "COF2-369"},
  }},
  {30, {  // SF6
    {'1', "SF6-29"},
  }},
  {31, {  // H2S
    {'1', "H2S-121"},
    {'2', "H2S-141"},
    {'3', "H2S-131"},
  }},
  {32, {  // HCOOH
    {'1', "HCOOH-1261"},  // Modified from HITRAN because ARTS does not use AFGL notation
  }},
  {33, {  // HO2
    {'1', "HO2-166"},
  }},
  {34, {  // O
    {'1', "O-6"},
  }},
  {35, {  // ClONO2
    {'1', "ClONO2-5646"},
    {'2', "ClONO2-7646"},
  }},
  {36, {  // NO+
    {'1', "NO+-46"},
  }},
  {37, {  // HOBr
    {'1', "HOBr-169"},
    {'2', "HOBr-161"},
  }},
  {38, {  // C2H4
    {'1', "C2H4-221"},
    {'2', "C2H4-231"},
  }},
  {39, {  // CH3OH
    {'1', "CH3OH-2161"},
  }},
  {40, {  // CH3Br
    {'1', "CH3Br-219"},
    {'2', "CH3Br-211"},
  }},
  {41, {  // CH3CN
    {'1', "CH3CN-211124"},  // Modified from HITRAN because ARTS does not use AFGL notation
  }},
  {42, {  // CF4
    {'1', "CF4-29"},
  }},
  {43, {  // C4H2
    {'1', "C4H2-2211"},
  }},
  {44, {  // HC3N
    {'1', "HC3N-12224"},  // Modified from HITRAN because ARTS does not use AFGL notation
  }},
  {45, {  // H2
    {'1', "H2-11"},
    {'2', "H2-12"},
  }},
  {46, {  // CS
    {'1', "CS-22"},
    {'2', "CS-24"},
    {'3', "CS-32"},
    {'4', "CS-23"},
  }},
  {47, {  // SO3
    {'1', "SO3-26"},
  }},
  {48, {  // C2N2
    {'1', "C2N2-4224"},
  }},
  {49, {  // COCl2
    {'1', "COCl2-2655"},
    {'2', "COCl2-2657"},
  }},
  {53, {  // CS2
    {'1', "CS2-222"},
    {'2', "CS2-224"},
    {'3', "CS2-223"},
    {'4', "CS2-232"},
  }},
};


/** Turns the string-map required at compile time into a species-map
 * to be used as a static runtime map
 */
std::map<Index, std::map<char, SpeciesTag>> to_species_map(const std::map<Index, std::map<char, const char * const>>& string_map) {
  std::map<Index, std::map<char, SpeciesTag>> species_map;
  for (auto& specs: string_map) {
    ARTS_ASSERT(specs.second.find('1') not_eq specs.second.cend(),
                "Must have species '1' in map")
    for (auto& isot: specs.second) {
      species_map[specs.first][isot.first] = SpeciesTag(isot.second);
    }
  }
  return species_map;
}

/** Selects the map and returns it.  Uses templates to avoid weird warnings */
template <Type t>
const std::map<Index, std::map<char, SpeciesTag>> select_hitran_map() {
  static_assert(good_enum(t), "Bad enum encountered, something is amiss");
  
  /*
  // Debug code to uncomment to check that the output is as expected
  auto x = to_species_map(latest_molparam_map);
  for (auto y: x) {
    std::cout << y.first << '\n';
    for (auto z: y.second) {
      std:: cout << z.first << ' ' << z.second << '\n';
    }
  }
  */
  
  if constexpr (t == Type::Newest) {
    return to_species_map(latest_molparam_map);
  } else if constexpr (t == Type::Pre2012CO2Change) {
    return to_species_map(pre2012co2change_molparam_map);
  } else {
    std::terminate();
  }
}

/** Returns the species if possible or throws an error if it cannot be found
 * 
 * @param[in] molnum The hitran molecular number
 * @param[in] isonum The hitran character representing an isotopologue
 * @return An ARTS identifier of the species' as a transition
 */
template <Type t> QuantumIdentifier from_mol_iso(Index molnum, char isonum) {
  static_assert(good_enum(t), "Bad enum encountered, something is amiss. enum is: ");
  
  // Map of all HITRAN species
  static const std::map<Index, std::map<char, SpeciesTag>> hitmap = select_hitran_map<t>();
  
  // Search the map with pointers so that we can throw manually if something is bad
  if (auto species_list = hitmap.find(molnum); species_list not_eq hitmap.cend()) {
    if (auto species_name = species_list -> second.find(isonum); species_name not_eq species_list -> second.cend()) {
      const SpeciesTag& spec = species_name -> second;
      return QuantumIdentifier(QuantumIdentifier::TRANSITION, spec.Species(), spec.Isotopologue());
    } else {
      ARTS_USER_ERROR("Species ", molnum, " has no isotopologue ", isonum,
                      " in ARTS' HITRAN implementation.\n",
                      "(Species is ", hitmap.at(molnum).at('1').SpeciesNameMain(), ")\n"
                      "If you are using a new version of HITRAN that has added the isotopologue, please consider\n"
                      "contacting the ARTS developers so we can append the species to our list and make this work.\n"
                      "Note that you are calling this templated function as the ", t, " template")
    }
  } else {
    ARTS_USER_ERROR("Species ", molnum, " does not exist in ARTS' HITRAN implementation\n"
                    "If you are using a new version of HITRAN that has added the species, please consider\n"
                    "contacting the ARTS developers so we can append the species to our list and make this work.\n"
                    "Note that you are calling this templated function as the ", t, " template");
  }
}

QuantumIdentifier from_lookup(Index mol, char isochar, Type type) {
  switch (type) {
    case Type::Pre2012CO2Change:
      return from_mol_iso<Type::Pre2012CO2Change>(mol, isochar);
    case Type::Newest:
      return from_mol_iso<Type::Newest>(mol, isochar);
    case Type::FINAL: {/* leave last */}
  }
  return {};
}

}  // namespace Hitran
