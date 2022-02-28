#ifndef species_h
#define species_h

#include "enums.h"
#include "matpack.h"

namespace Species {
/** named species */
ENUMCLASS(Species, unsigned char,
          Bath,
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
          particles
         )  // Species

constexpr std::string_view toShortName(Species x) noexcept {
  switch (x) {
    case Species::Bath:
      return "AIR";
    case Species::Water:
      return "H2O";
    case Species::CarbonDioxide:
      return "CO2";
    case Species::Ozone:
      return "O3";
    case Species::NitrogenOxide:
      return "N2O";
    case Species::CarbonMonoxide:
      return "CO";
    case Species::Methane:
      return "CH4";
    case Species::Oxygen:
      return "O2";
    case Species::NitricOxide:
      return "NO";
    case Species::SulfurDioxide:
      return "SO2";
    case Species::NitrogenDioxide:
      return "NO2";
    case Species::Ammonia:
      return "NH3";
    case Species::NitricAcid:
      return "HNO3";
    case Species::Hydroxyl:
      return "OH";
    case Species::HydrogenFluoride:
      return "HF";
    case Species::HydrogenChloride:
      return "HCl";
    case Species::HydrogenBromide:
      return "HBr";
    case Species::HydrogenIodide:
      return "HI";
    case Species::ChlorineMonoxide:
      return "ClO";
    case Species::CarbonylSulfide:
      return "OCS";
    case Species::Formaldehyde:
      return "H2CO";
    case Species::HeavyFormaldehyde:
      return "HDCO";
    case Species::VeryHeavyFormaldehyde:
      return "D2CO";
    case Species::HypochlorousAcid:
      return "HOCl";
    case Species::Nitrogen:
      return "N2";
    case Species::HydrogenCyanide:
      return "HCN";
    case Species::Chloromethane:
      return "CH3Cl";
    case Species::HydrogenPeroxide:
      return "H2O2";
    case Species::Acetylene:
      return "C2H2";
    case Species::Ethane:
      return "C2H6";
    case Species::Phosphine:
      return "PH3";
    case Species::CarbonylFluoride:
      return "COF2";
    case Species::SulfurHexafluoride:
      return "SF6";
    case Species::HydrogenSulfide:
      return "H2S";
    case Species::FormicAcid:
      return "HCOOH";
    case Species::LeftHeavyFormicAcid:
      return "DCOOH";
    case Species::RightHeavyFormicAcid:
      return "HCOOD";
    case Species::Hydroperoxyl:
      return "HO2";
    case Species::OxygenAtom:
      return "O";
    case Species::ChlorineNitrate:
      return "ClONO2";
    case Species::NitricOxideCation:
      return "NO+";
    case Species::HypobromousAcid:
      return "HOBr";
    case Species::Ethylene:
      return "C2H4";
    case Species::Methanol:
      return "CH3OH";
    case Species::Bromomethane:
      return "CH3Br";
    case Species::Acetonitrile:
      return "CH3CN";
    case Species::HeavyAcetonitrile:
      return "CH2DCN";
    case Species::CarbonTetrafluoride:
      return "CF4";
    case Species::Diacetylene:
      return "C4H2";
    case Species::Cyanoacetylene:
      return "HC3N";
    case Species::Hydrogen:
      return "H2";
    case Species::CarbonMonosulfide:
      return "CS";
    case Species::SulfurTrioxide:
      return "SO3";
    case Species::Cyanogen:
      return "C2N2";
    case Species::Phosgene:
      return "COCl2";
    case Species::SulfurMonoxide:
      return "SO";
    case Species::CarbonDisulfide:
      return "CS2";
    case Species::Methyl:
      return "CH3";
    case Species::Cyclopropene:
      return "C3H4";
    case Species::SulfuricAcid:
      return "H2SO4";
    case Species::HydrogenIsocyanide:
      return "HNC";
    case Species::BromineMonoxide:
      return "BrO";
    case Species::ChlorineDioxide:
      return "OClO";
    case Species::Propane:
      return "C3H8";
    case Species::Helium:
      return "He";
    case Species::ChlorineMonoxideDimer:
      return "Cl2O2";
    case Species::HydrogenAtom:
      return "H";
    case Species::Argon:
      return "Ar";
    case Species::Hexafluoroethane: return "C2F6";
    case Species::Perfluoropropane: return "C3F8";
    case Species::Perfluorobutane: return "C4F10";
    case Species::Perfluoropentane: return "C5F12";
    case Species::Perfluorohexane: return "C6F14";
    case Species::Perfluorooctane: return "C8F18";
    case Species::Perfluorocyclobutane: return "cC4F8";
    case Species::CarbonTetrachloride: return "CCl4";
    case Species::CFC11: return "CFC11";
    case Species::CFC113: return "CFC113";
    case Species::CFC114: return "CFC114";
    case Species::CFC115: return "CFC115";
    case Species::CFC12: return "CFC12";
    case Species::Dichloromethane: return "CH2Cl2";
    case Species::Trichloroethane: return "CH3CCl3";
    case Species::Trichloromethane: return "CHCl3";
    case Species::Bromochlorodifluoromethane: return "Halon1211";
    case Species::Bromotrifluoromethane: return "Halon1301";
    case Species::Dibromotetrafluoroethane: return "Halon2402";
    case Species::HCFC141b: return "HCFC141b";
    case Species::HCFC142b: return "HCFC142b";
    case Species::HCFC22: return "HCFC22";
    case Species::HFC125: return "HFC125";
    case Species::HFC134a: return "HFC134a";
    case Species::HFC143a: return "HFC143a";
    case Species::HFC152a: return "HFC152a";
    case Species::HFC227ea: return "HFC227ea";
    case Species::HFC23: return "HFC23";
    case Species::HFC236fa: return "HFC236fa";
    case Species::HFC245fa: return "HFC245fa";
    case Species::HFC32: return "HFC32";
    case Species::HFC365mfc: return "HFC365mfc";
    case Species::NitrogenTrifluoride: return "NF3";
    case Species::SulfurylFluoride: return "SO2F2";
    case Species::HFC4310mee: return "HFC4310mee";
    case Species::Germane: return "GeH4";
    case Species::Iodomethane: return "CH3I";
    case Species::Fluoromethane: return "CH3F";
    case Species::liquidcloud:
      return "liquidcloud";
    case Species::icecloud:
      return "icecloud";
    case Species::rain:
      return "rain";
    case Species::free_electrons:
      return "free_electrons";
    case Species::particles:
      return "particles";
    case Species::FINAL: { /* Leave last */
    }
  }
  return "InvalidSpecies";
}

constexpr Species fromShortName(const std::string_view x) noexcept {
  if (x == "AIR") 
    return Species::Bath;
  if (x == "H2O")
    return Species::Water;
  if (x == "CO2")
    return Species::CarbonDioxide;
  if (x == "O3")
    return Species::Ozone;
  if (x == "N2O")
    return Species::NitrogenOxide;
  if (x == "CO")
    return Species::CarbonMonoxide;
  if (x == "CH4")
    return Species::Methane;
  if (x == "O2")
    return Species::Oxygen;
  if (x == "NO")
    return Species::NitricOxide;
  if (x == "SO2")
    return Species::SulfurDioxide;
  if (x == "NO2")
    return Species::NitrogenDioxide;
  if (x == "NH3")
    return Species::Ammonia;
  if (x == "HNO3")
    return Species::NitricAcid;
  if (x == "OH")
    return Species::Hydroxyl;
  if (x == "HF")
    return Species::HydrogenFluoride;
  if (x == "HCl")
    return Species::HydrogenChloride;
  if (x == "HBr")
    return Species::HydrogenBromide;
  if (x == "HI")
    return Species::HydrogenIodide;
  if (x == "ClO")
    return Species::ChlorineMonoxide;
  if (x == "OCS")
    return Species::CarbonylSulfide;
  if (x == "H2CO")
    return Species::Formaldehyde;
  if (x == "HDCO")
    return Species::HeavyFormaldehyde;
  if (x == "D2CO")
    return Species::VeryHeavyFormaldehyde;
  if (x == "HOCl")
    return Species::HypochlorousAcid;
  if (x == "N2")
    return Species::Nitrogen;
  if (x == "HCN")
    return Species::HydrogenCyanide;
  if (x == "CH3Cl")
    return Species::Chloromethane;
  if (x == "H2O2")
    return Species::HydrogenPeroxide;
  if (x == "C2H2")
    return Species::Acetylene;
  if (x == "C2H6")
    return Species::Ethane;
  if (x == "PH3")
    return Species::Phosphine;
  if (x == "COF2")
    return Species::CarbonylFluoride;
  if (x == "SF6")
    return Species::SulfurHexafluoride;
  if (x == "H2S")
    return Species::HydrogenSulfide;
  if (x == "HCOOH")
    return Species::FormicAcid;
  if (x == "DCOOH")
    return Species::LeftHeavyFormicAcid;
  if (x == "HCOOD")
    return Species::RightHeavyFormicAcid;
  if (x == "HO2")
    return Species::Hydroperoxyl;
  if (x == "O")
    return Species::OxygenAtom;
  if (x == "ClONO2")
    return Species::ChlorineNitrate;
  if (x == "NO+")
    return Species::NitricOxideCation;
  if (x == "HOBr")
    return Species::HypobromousAcid;
  if (x == "C2H4")
    return Species::Ethylene;
  if (x == "CH3OH")
    return Species::Methanol;
  if (x == "CH3Br")
    return Species::Bromomethane;
  if (x == "CH3CN")
    return Species::Acetonitrile;
  if (x == "CH2DCN")
    return Species::HeavyAcetonitrile;
  if (x == "CF4")
    return Species::CarbonTetrafluoride;
  if (x == "C4H2")
    return Species::Diacetylene;
  if (x == "HC3N")
    return Species::Cyanoacetylene;
  if (x == "H2")
    return Species::Hydrogen;
  if (x == "CS")
    return Species::CarbonMonosulfide;
  if (x == "SO3")
    return Species::SulfurTrioxide;
  if (x == "C2N2")
    return Species::Cyanogen;
  if (x == "COCl2")
    return Species::Phosgene;
  if (x == "SO")
    return Species::SulfurMonoxide;
  if (x == "CS2")
    return Species::CarbonDisulfide;
  if (x == "CH3")
    return Species::Methyl;
  if (x == "C3H4")
    return Species::Cyclopropene;
  if (x == "H2SO4")
    return Species::SulfuricAcid;
  if (x == "HNC")
    return Species::HydrogenIsocyanide;
  if (x == "BrO")
    return Species::BromineMonoxide;
  if (x == "OClO")
    return Species::ChlorineDioxide;
  if (x == "C3H8")
    return Species::Propane;
  if (x == "He")
    return Species::Helium;
  if (x == "Cl2O2")
    return Species::ChlorineMonoxideDimer;
  if (x == "H")
    return Species::HydrogenAtom;
  if (x == "Ar")
    return Species::Argon;
  if (x == "C2F6") return Species::Hexafluoroethane;
  if (x == "C3F8") return Species::Perfluoropropane;
  if (x == "C4F10") return Species::Perfluorobutane;
  if (x == "C5F12") return Species::Perfluoropentane;
  if (x == "C6F14") return Species::Perfluorohexane;
  if (x == "C8F18") return Species::Perfluorooctane;
  if (x == "cC4F8") return Species::Perfluorocyclobutane;
  if (x == "CCl4") return Species::CarbonTetrachloride;
  if (x == "CFC11") return Species::CFC11;
  if (x == "CFC113") return Species::CFC113;
  if (x == "CFC114") return Species::CFC114;
  if (x == "CFC115") return Species::CFC115;
  if (x == "CFC12") return Species::CFC12;
  if (x == "CH2Cl2") return Species::Dichloromethane;
  if (x == "CH3CCl3") return Species::Trichloroethane;
  if (x == "CHCl3") return Species::Trichloromethane;
  if (x == "Halon1211") return Species::Bromochlorodifluoromethane;
  if (x == "Halon1301") return Species::Bromotrifluoromethane;
  if (x == "Halon2402") return Species::Dibromotetrafluoroethane;
  if (x == "HCFC141b") return Species::HCFC141b;
  if (x == "HCFC142b") return Species::HCFC142b;
  if (x == "HCFC22") return Species::HCFC22;
  if (x == "HFC125") return Species::HFC125;
  if (x == "HFC134a") return Species::HFC134a;
  if (x == "HFC143a") return Species::HFC143a;
  if (x == "HFC152a") return Species::HFC152a;
  if (x == "HFC227ea") return Species::HFC227ea;
  if (x == "HFC23") return Species::HFC23;
  if (x == "HFC236fa") return Species::HFC236fa;
  if (x == "HFC245fa") return Species::HFC245fa;
  if (x == "HFC32") return Species::HFC32;
  if (x == "HFC365mfc") return Species::HFC365mfc;
  if (x == "NF3") return Species::NitrogenTrifluoride;
  if (x == "SO2F2") return Species::SulfurylFluoride;
  if (x == "HFC4310mee") return Species::HFC4310mee;
  if (x == "GeH4") return Species::Germane;
  if (x == "CH3I") return Species::Iodomethane;
  if (x == "CH3F") return Species::Fluoromethane;
  if (x == "liquidcloud")
    return Species::liquidcloud;
  if (x == "icecloud")
    return Species::icecloud;
  if (x == "rain")
    return Species::rain;
  if (x == "free_electrons")
    return Species::free_electrons;
  if (x == "particles")
    return Species::particles;
  return Species::FINAL;
}
} // namespace Species

#endif  // species_h
