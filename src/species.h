#ifndef species_h
#define species_h

#include "enums.h"
#include "matpack.h"

namespace Species {
/** named species */
ENUMCLASS(Species, unsigned char, Bath, Water, CarbonDioxide, Ozone,
          NitrogenOxide, CarbonMonoxide, Methane, Oxygen, NitricOxide,
          SulfurDioxide, NitrogenDioxide, Ammonia, NitricAcid, Hydroxyl,
          HydrogenFluoride, HydrogenChloride, HydrogenBromide, HydrogenIodide,
          ChlorineMonoxide, CarbonylSulfide, Formaldehyde, HypochlorousAcid,
          Nitrogen, HydrogenCyanide, MethylChloride, HydrogenPeroxide,
          Acetylene, Ethane, Phosphine, CarbonylFluoride, SulfurHexafluoride,
          HydrogenSulfide, FormicAcid, Hydroperoxyl, OxygenAtom,
          ChlorineNitrate, NitricOxideCation, HypobromousAcid, Ethylene,
          Methanol, MethylBromide, Acetonitrile, CarbonTetrafluoride, Diacetylene,
          Cyanoacetylene, Hydrogen, CarbonMonosulfide, SulfurTrioxide, Cyanogen,
          Phosgene, SulfurMonoxide, CarbonDisulfide, Methyl, Cyclopropene,
          SulfuricAcid, HydrogenIsocyanide, BromineMonoxide, ChlorineDioxide, Propane,
          Helium, ChlorineMonoxideDimer, HydrogenAtom, Argon,
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
          HFC245fa,
          HFC32,
          NitrogenTrifluoride,
          SulfurylFluoride,
          HFC4310mee,
          liquidcloud, icecloud, rain, free_electrons, particles
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
    case Species::HypochlorousAcid:
      return "HOCl";
    case Species::Nitrogen:
      return "N2";
    case Species::HydrogenCyanide:
      return "HCN";
    case Species::MethylChloride:
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
    case Species::MethylBromide:
      return "CH3Br";
    case Species::Acetonitrile:
      return "CH3CN";
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
    case Species::HFC245fa: return "HFC245fa";
    case Species::HFC32: return "HFC32";
    case Species::NitrogenTrifluoride: return "NF3";
    case Species::SulfurylFluoride: return "SO2F2";
    case Species::HFC4310mee: return "HFC4310mee";
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
  return "";
}

constexpr Species fromShortName(const std::string_view x) noexcept {
  if (x == "AIR") 
    return Species::Bath;
  else if (x == "H2O")
    return Species::Water;
  else if (x == "CO2")
    return Species::CarbonDioxide;
  else if (x == "O3")
    return Species::Ozone;
  else if (x == "N2O")
    return Species::NitrogenOxide;
  else if (x == "CO")
    return Species::CarbonMonoxide;
  else if (x == "CH4")
    return Species::Methane;
  else if (x == "O2")
    return Species::Oxygen;
  else if (x == "NO")
    return Species::NitricOxide;
  else if (x == "SO2")
    return Species::SulfurDioxide;
  else if (x == "NO2")
    return Species::NitrogenDioxide;
  else if (x == "NH3")
    return Species::Ammonia;
  else if (x == "HNO3")
    return Species::NitricAcid;
  else if (x == "OH")
    return Species::Hydroxyl;
  else if (x == "HF")
    return Species::HydrogenFluoride;
  else if (x == "HCl")
    return Species::HydrogenChloride;
  else if (x == "HBr")
    return Species::HydrogenBromide;
  else if (x == "HI")
    return Species::HydrogenIodide;
  else if (x == "ClO")
    return Species::ChlorineMonoxide;
  else if (x == "OCS")
    return Species::CarbonylSulfide;
  else if (x == "H2CO")
    return Species::Formaldehyde;
  else if (x == "HOCl")
    return Species::HypochlorousAcid;
  else if (x == "N2")
    return Species::Nitrogen;
  else if (x == "HCN")
    return Species::HydrogenCyanide;
  else if (x == "CH3Cl")
    return Species::MethylChloride;
  else if (x == "H2O2")
    return Species::HydrogenPeroxide;
  else if (x == "C2H2")
    return Species::Acetylene;
  else if (x == "C2H6")
    return Species::Ethane;
  else if (x == "PH3")
    return Species::Phosphine;
  else if (x == "COF2")
    return Species::CarbonylFluoride;
  else if (x == "SF6")
    return Species::SulfurHexafluoride;
  else if (x == "H2S")
    return Species::HydrogenSulfide;
  else if (x == "HCOOH")
    return Species::FormicAcid;
  else if (x == "HO2")
    return Species::Hydroperoxyl;
  else if (x == "O")
    return Species::OxygenAtom;
  else if (x == "ClONO2")
    return Species::ChlorineNitrate;
  else if (x == "NO+")
    return Species::NitricOxideCation;
  else if (x == "HOBr")
    return Species::HypobromousAcid;
  else if (x == "C2H4")
    return Species::Ethylene;
  else if (x == "CH3OH")
    return Species::Methanol;
  else if (x == "CH3Br")
    return Species::MethylBromide;
  else if (x == "CH3CN")
    return Species::Acetonitrile;
  else if (x == "CF4")
    return Species::CarbonTetrafluoride;
  else if (x == "C4H2")
    return Species::Diacetylene;
  else if (x == "HC3N")
    return Species::Cyanoacetylene;
  else if (x == "H2")
    return Species::Hydrogen;
  else if (x == "CS")
    return Species::CarbonMonosulfide;
  else if (x == "SO3")
    return Species::SulfurTrioxide;
  else if (x == "C2N2")
    return Species::Cyanogen;
  else if (x == "COCl2")
    return Species::Phosgene;
  else if (x == "SO")
    return Species::SulfurMonoxide;
  else if (x == "CS2")
    return Species::CarbonDisulfide;
  else if (x == "CH3")
    return Species::Methyl;
  else if (x == "C3H4")
    return Species::Cyclopropene;
  else if (x == "H2SO4")
    return Species::SulfuricAcid;
  else if (x == "HNC")
    return Species::HydrogenIsocyanide;
  else if (x == "BrO")
    return Species::BromineMonoxide;
  else if (x == "OClO")
    return Species::ChlorineDioxide;
  else if (x == "C3H8")
    return Species::Propane;
  else if (x == "He")
    return Species::Helium;
  else if (x == "Cl2O2")
    return Species::ChlorineMonoxideDimer;
  else if (x == "H")
    return Species::HydrogenAtom;
  else if (x == "Ar")
    return Species::Argon;
  else if (x == "C2F6") return Species::Hexafluoroethane;
  else if (x == "C3F8") return Species::Perfluoropropane;
  else if (x == "C4F10") return Species::Perfluorobutane;
  else if (x == "C5F12") return Species::Perfluoropentane;
  else if (x == "C6F14") return Species::Perfluorohexane;
  else if (x == "C8F18") return Species::Perfluorooctane;
  else if (x == "cC4F8") return Species::Perfluorocyclobutane;
  else if (x == "CCl4") return Species::CarbonTetrachloride;
  else if (x == "CFC11") return Species::CFC11;
  else if (x == "CFC113") return Species::CFC113;
  else if (x == "CFC114") return Species::CFC114;
  else if (x == "CFC115") return Species::CFC115;
  else if (x == "CFC12") return Species::CFC12;
  else if (x == "CH2Cl2") return Species::Dichloromethane;
  else if (x == "CH3CCl3") return Species::Trichloroethane;
  else if (x == "CHCl3") return Species::Trichloromethane;
  else if (x == "Halon1211") return Species::Bromochlorodifluoromethane;
  else if (x == "Halon1301") return Species::Bromotrifluoromethane;
  else if (x == "Halon2402") return Species::Dibromotetrafluoroethane;
  else if (x == "HCFC141b") return Species::HCFC141b;
  else if (x == "HCFC142b") return Species::HCFC142b;
  else if (x == "HCFC22") return Species::HCFC22;
  else if (x == "HFC125") return Species::HFC125;
  else if (x == "HFC134a") return Species::HFC134a;
  else if (x == "HFC143a") return Species::HFC143a;
  else if (x == "HFC152a") return Species::HFC152a;
  else if (x == "HFC227ea") return Species::HFC227ea;
  else if (x == "HFC23") return Species::HFC23;
  else if (x == "HFC245fa") return Species::HFC245fa;
  else if (x == "HFC32") return Species::HFC32;
  else if (x == "NF3") return Species::NitrogenTrifluoride;
  else if (x == "SO2F2") return Species::SulfurylFluoride;
  else if (x == "HFC4310mee") return Species::HFC4310mee;
  else if (x == "liquidcloud")
    return Species::liquidcloud;
  else if (x == "icecloud")
    return Species::icecloud;
  else if (x == "rain")
    return Species::rain;
  else if (x == "free_electrons")
    return Species::free_electrons;
  else if (x == "particles")
    return Species::particles;
  else
    return Species::FINAL;
}
}

#endif  // species_h
