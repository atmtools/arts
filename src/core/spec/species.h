#pragma once

#include <array.h>

#include <algorithm>
#include <string_view>

#include "enums.h"

namespace Species {
/** named species */
ENUMCLASS(Species,
          unsigned char,
          Bath,  // This must be first, it represents any non-specific gas!
          Water,
          CarbonDioxide,
          Ozone,
          NitrogenOxide,
          CarbonMonoxide,
          Methane,
          Oxygen,
          NitricOxide,
          SulfurDioxide,
          NitrogenDioxide,
          Ammonia,
          NitricAcid,
          Hydroxyl,
          HydrogenFluoride,
          HydrogenChloride,
          HydrogenBromide,
          HydrogenIodide,
          ChlorineMonoxide,
          CarbonylSulfide,
          Formaldehyde,
          HeavyFormaldehyde,
          VeryHeavyFormaldehyde,
          HypochlorousAcid,
          Nitrogen,
          HydrogenCyanide,
          Chloromethane,
          HydrogenPeroxide,
          Acetylene,
          Ethane,
          Phosphine,
          CarbonylFluoride,
          SulfurHexafluoride,
          HydrogenSulfide,
          FormicAcid,
          LeftHeavyFormicAcid,
          RightHeavyFormicAcid,
          Hydroperoxyl,
          OxygenAtom,
          ChlorineNitrate,
          NitricOxideCation,
          ChlorineDioxide,
          BromineMonoxide,
          SulfuricAcid,
          ChlorineMonoxideDimer,
          HypobromousAcid,
          Ethylene,
          Methanol,
          Bromomethane,
          Acetonitrile,
          HeavyAcetonitrile,
          CarbonTetrafluoride,
          Diacetylene,
          Cyanoacetylene,
          Hydrogen,
          CarbonMonosulfide,
          HydrogenIsocyanide,
          SulfurMonoxide,
          Propane,
          HydrogenAtom,
          Helium,
          Argon,
          SulfurTrioxide,
          Cyanogen,
          Phosgene,
          CarbonDisulfide,
          Methyl,
          Cyclopropene,
          Hexafluoroethane,
          Perfluoropropane,
          Perfluorobutane,
          Perfluoropentane,
          Perfluorohexane,
          Perfluorooctane,
          Perfluorocyclobutane,
          CarbonTetrachloride,
          CFC11,
          CFC113,
          CFC114,
          CFC115,
          CFC12,
          Dichloromethane,
          Trichloroethane,
          Trichloromethane,
          Bromochlorodifluoromethane,
          Bromotrifluoromethane,
          Dibromotetrafluoroethane,
          HCFC141b,
          HCFC142b,
          HCFC22,
          HFC125,
          HFC134a,
          HFC143a,
          HFC152a,
          HFC227ea,
          HFC23,
          HFC236fa,
          HFC245fa,
          HFC32,
          HFC365mfc,
          NitrogenTrifluoride,
          SulfurylFluoride,
          HFC4310mee,
          Germane,
          Iodomethane,
          Fluoromethane,
          liquidcloud,
          icecloud,
          rain,
          free_electrons,
          particles)  // Species

struct short_name {
  std::string_view name;
  Species species;

  friend std::ostream& operator<<(std::ostream& os, const short_name& x);
};

enum class short_name_sort_by : bool { name, spec };

constexpr std::array<short_name, static_cast<Size>(Species::FINAL)>
short_name_species(short_name_sort_by Sort) noexcept {
  auto x =
      std::array{short_name{"AIR", Species::Bath},
                 short_name{"H2O", Species::Water},
                 short_name{"CO2", Species::CarbonDioxide},
                 short_name{"O3", Species::Ozone},
                 short_name{"N2O", Species::NitrogenOxide},
                 short_name{"CO", Species::CarbonMonoxide},
                 short_name{"CH4", Species::Methane},
                 short_name{"O2", Species::Oxygen},
                 short_name{"NO", Species::NitricOxide},
                 short_name{"SO2", Species::SulfurDioxide},
                 short_name{"NO2", Species::NitrogenDioxide},
                 short_name{"NH3", Species::Ammonia},
                 short_name{"HNO3", Species::NitricAcid},
                 short_name{"OH", Species::Hydroxyl},
                 short_name{"HF", Species::HydrogenFluoride},
                 short_name{"HCl", Species::HydrogenChloride},
                 short_name{"HBr", Species::HydrogenBromide},
                 short_name{"HI", Species::HydrogenIodide},
                 short_name{"ClO", Species::ChlorineMonoxide},
                 short_name{"OCS", Species::CarbonylSulfide},
                 short_name{"H2CO", Species::Formaldehyde},
                 short_name{"HDCO", Species::HeavyFormaldehyde},
                 short_name{"D2CO", Species::VeryHeavyFormaldehyde},
                 short_name{"HOCl", Species::HypochlorousAcid},
                 short_name{"N2", Species::Nitrogen},
                 short_name{"HCN", Species::HydrogenCyanide},
                 short_name{"CH3Cl", Species::Chloromethane},
                 short_name{"H2O2", Species::HydrogenPeroxide},
                 short_name{"C2H2", Species::Acetylene},
                 short_name{"C2H6", Species::Ethane},
                 short_name{"PH3", Species::Phosphine},
                 short_name{"COF2", Species::CarbonylFluoride},
                 short_name{"SF6", Species::SulfurHexafluoride},
                 short_name{"H2S", Species::HydrogenSulfide},
                 short_name{"HCOOH", Species::FormicAcid},
                 short_name{"DCOOH", Species::LeftHeavyFormicAcid},
                 short_name{"HCOOD", Species::RightHeavyFormicAcid},
                 short_name{"HO2", Species::Hydroperoxyl},
                 short_name{"O", Species::OxygenAtom},
                 short_name{"ClONO2", Species::ChlorineNitrate},
                 short_name{"NO+", Species::NitricOxideCation},
                 short_name{"HOBr", Species::HypobromousAcid},
                 short_name{"C2H4", Species::Ethylene},
                 short_name{"CH3OH", Species::Methanol},
                 short_name{"CH3Br", Species::Bromomethane},
                 short_name{"CH3CN", Species::Acetonitrile},
                 short_name{"CH2DCN", Species::HeavyAcetonitrile},
                 short_name{"CF4", Species::CarbonTetrafluoride},
                 short_name{"C4H2", Species::Diacetylene},
                 short_name{"HC3N", Species::Cyanoacetylene},
                 short_name{"H2", Species::Hydrogen},
                 short_name{"CS", Species::CarbonMonosulfide},
                 short_name{"SO3", Species::SulfurTrioxide},
                 short_name{"C2N2", Species::Cyanogen},
                 short_name{"COCl2", Species::Phosgene},
                 short_name{"SO", Species::SulfurMonoxide},
                 short_name{"CS2", Species::CarbonDisulfide},
                 short_name{"CH3", Species::Methyl},
                 short_name{"C3H4", Species::Cyclopropene},
                 short_name{"H2SO4", Species::SulfuricAcid},
                 short_name{"HNC", Species::HydrogenIsocyanide},
                 short_name{"BrO", Species::BromineMonoxide},
                 short_name{"OClO", Species::ChlorineDioxide},
                 short_name{"C3H8", Species::Propane},
                 short_name{"He", Species::Helium},
                 short_name{"Cl2O2", Species::ChlorineMonoxideDimer},
                 short_name{"H", Species::HydrogenAtom},
                 short_name{"Ar", Species::Argon},
                 short_name{"C2F6", Species::Hexafluoroethane},
                 short_name{"C3F8", Species::Perfluoropropane},
                 short_name{"C4F10", Species::Perfluorobutane},
                 short_name{"C5F12", Species::Perfluoropentane},
                 short_name{"C6F14", Species::Perfluorohexane},
                 short_name{"C8F18", Species::Perfluorooctane},
                 short_name{"cC4F8", Species::Perfluorocyclobutane},
                 short_name{"CCl4", Species::CarbonTetrachloride},
                 short_name{"CFC11", Species::CFC11},
                 short_name{"CFC113", Species::CFC113},
                 short_name{"CFC114", Species::CFC114},
                 short_name{"CFC115", Species::CFC115},
                 short_name{"CFC12", Species::CFC12},
                 short_name{"CH2Cl2", Species::Dichloromethane},
                 short_name{"CH3CCl3", Species::Trichloroethane},
                 short_name{"CHCl3", Species::Trichloromethane},
                 short_name{"Halon1211", Species::Bromochlorodifluoromethane},
                 short_name{"Halon1301", Species::Bromotrifluoromethane},
                 short_name{"Halon2402", Species::Dibromotetrafluoroethane},
                 short_name{"HCFC141b", Species::HCFC141b},
                 short_name{"HCFC142b", Species::HCFC142b},
                 short_name{"HCFC22", Species::HCFC22},
                 short_name{"HFC125", Species::HFC125},
                 short_name{"HFC134a", Species::HFC134a},
                 short_name{"HFC143a", Species::HFC143a},
                 short_name{"HFC152a", Species::HFC152a},
                 short_name{"HFC227ea", Species::HFC227ea},
                 short_name{"HFC23", Species::HFC23},
                 short_name{"HFC236fa", Species::HFC236fa},
                 short_name{"HFC245fa", Species::HFC245fa},
                 short_name{"HFC32", Species::HFC32},
                 short_name{"HFC365mfc", Species::HFC365mfc},
                 short_name{"NF3", Species::NitrogenTrifluoride},
                 short_name{"SO2F2", Species::SulfurylFluoride},
                 short_name{"HFC4310mee", Species::HFC4310mee},
                 short_name{"GeH4", Species::Germane},
                 short_name{"CH3I", Species::Iodomethane},
                 short_name{"CH3F", Species::Fluoromethane},
                 short_name{"liquidcloud", Species::liquidcloud},
                 short_name{"icecloud", Species::icecloud},
                 short_name{"rain", Species::rain},
                 short_name{"free_electrons", Species::free_electrons},
                 short_name{"particles", Species::particles}};


    
  if (Sort == short_name_sort_by::name) {
    std::ranges::sort(x, {}, &short_name::name);
  } else if (Sort == short_name_sort_by::spec) {
    std::ranges::sort(x, {}, &short_name::species);
  }
  return x;
};

static constexpr auto short_names_spec =
    short_name_species(short_name_sort_by::spec);
static constexpr auto short_names_name =
    short_name_species(short_name_sort_by::name);

static_assert (std::ranges::adjacent_find(short_names_spec, {}, &short_name::species) ==
               short_names_spec.end(), "Non-unique species enum in short-name set");

static_assert (std::ranges::adjacent_find(short_names_name, {}, &short_name::name) ==
               short_names_name.end(), "Non-unique species name in short-name set");

constexpr std::string_view toShortName(Species x) noexcept {
  const auto* it =
      std::ranges::lower_bound(short_names_spec, x, {}, &short_name::species);
  if (it != short_names_spec.end() && it->species == x) {
    return it->name;
  }
  return "InvalidSpecies";
}

constexpr Species fromShortName(const std::string_view x) noexcept {
  const auto* it =
      std::ranges::lower_bound(short_names_name, x, {}, &short_name::name);
  if (it != short_names_name.end() && it->name == x) {
    return it->species;
  }
  return Species::FINAL;
}

Species toSpeciesEnumOrThrow(const std::string_view x);
}  // namespace Species

using SpeciesEnum = Species::Species;
using ArrayOfSpecies = Array<SpeciesEnum>;
