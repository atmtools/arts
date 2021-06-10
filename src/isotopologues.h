#ifndef isotopologues_h
#define isotopologues_h

#include <limits>

#include "mystring.h"
#include "nonstd.h"
#include "species.h"

namespace Species {
constexpr std::string_view Joker = "*";

/** Struct containing all information needed about one isotope */
struct IsotopeRecord {
  //! Species type as defined in species.h
  Species spec;
  
  //! A custom name that is unique for this Species type
  std::string_view isotname;
  
  //! The mass of the isotope in units of grams per mol.  It is Nan if not defined
  Numeric mass;
  
  //! The degeneracy of states of the molecule.  It is -1 if not defined.
  Index gi;
  constexpr explicit IsotopeRecord(Species spec_, const std::string_view isotname_=Joker, Numeric mass_=std::numeric_limits<Numeric>::quiet_NaN(), Index gi_=-1) noexcept
  : spec(spec_), isotname(isotname_), mass(mass_), gi(gi_) {}
  constexpr IsotopeRecord() noexcept : IsotopeRecord(Species::FINAL) {}
  friend std::ostream& operator<<(std::ostream& os, const IsotopeRecord& ir) {
    return os << ir.spec << ' ' << ir.isotname << ' ' << ir.mass << ' ' << ir.gi;
  }
  constexpr bool operator==(const IsotopeRecord& that) const noexcept {
    return that.spec == spec and that.isotname == isotname;
  }
  
  //! A comparison with pure named string (this is not an exact comparison)
  constexpr bool operator==(const std::string_view specstr) const noexcept {
    auto lim = specstr.find('-');
    return (lim not_eq specstr.npos) and (fromShortName(specstr.substr(0, lim)) == spec) and (specstr.substr(lim+1) == isotname);
  }
  
  template <typename T> constexpr bool operator!=(T x)const noexcept {return not operator==(x);}
  
  String FullName() const noexcept {return String(toShortName(spec)) + String("-") + String(isotname);}
  constexpr bool joker() const noexcept {return isotname == Joker;}
};

#define deal_with_spec(SPEC) IsotopeRecord(Species::SPEC),

constexpr std::array Isotopologues {
  /** Water species **/
  deal_with_spec(Water)
  IsotopeRecord(fromShortName("H2O"), "161", 18.010565, 1),
  IsotopeRecord(fromShortName("H2O"), "162", 19.016740, 6),
  IsotopeRecord(fromShortName("H2O"), "171", 19.014780, 6),
  IsotopeRecord(fromShortName("H2O"), "172", 20.020956, 36),
  IsotopeRecord(fromShortName("H2O"), "181", 20.014811, 1),
  IsotopeRecord(fromShortName("H2O"), "182", 21.020985, 6),
  IsotopeRecord(fromShortName("H2O"), "262", 20.022915, 1),
  IsotopeRecord(fromShortName("H2O"), "CP98"),
  IsotopeRecord(fromShortName("H2O"), "ContMPM93"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContATM01"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContCKD222"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContCKD24"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContCKD242"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContCKDMT100"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContCKDMT252"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContCKDMT320"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContMaTippingType"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContStandardType"),
  IsotopeRecord(fromShortName("H2O"), "MPM87"),
  IsotopeRecord(fromShortName("H2O"), "MPM89"),
  IsotopeRecord(fromShortName("H2O"), "MPM93"),
  IsotopeRecord(fromShortName("H2O"), "PWR98"),
  IsotopeRecord(fromShortName("H2O"), "SelfContCKD222"),
  IsotopeRecord(fromShortName("H2O"), "SelfContCKD24"),
  IsotopeRecord(fromShortName("H2O"), "SelfContCKD242"),
  IsotopeRecord(fromShortName("H2O"), "SelfContCKDMT100"),
  IsotopeRecord(fromShortName("H2O"), "SelfContCKDMT252"),
  IsotopeRecord(fromShortName("H2O"), "SelfContCKDMT320"),
  IsotopeRecord(fromShortName("H2O"), "SelfContStandardType"),
  /** Water species **/
  
  /** Carbon dioxide species **/
  deal_with_spec(CarbonDioxide)
  IsotopeRecord(fromShortName("CO2"), "626", 43.989830, 1),
  IsotopeRecord(fromShortName("CO2"), "627", 44.994045, 6),
  IsotopeRecord(fromShortName("CO2"), "628", 45.994076, 1),
  IsotopeRecord(fromShortName("CO2"), "636", 44.993185, 2),
  IsotopeRecord(fromShortName("CO2"), "637", 45.997400, 12),
  IsotopeRecord(fromShortName("CO2"), "638", 46.997431, 2),
  IsotopeRecord(fromShortName("CO2"), "727", 45.998262, 1),
  IsotopeRecord(fromShortName("CO2"), "737", 47.001618, 2),
  IsotopeRecord(fromShortName("CO2"), "827", 46.998291, 6),
  IsotopeRecord(fromShortName("CO2"), "828", 47.998322, 1),
  IsotopeRecord(fromShortName("CO2"), "837", 48.001646, 12),
  IsotopeRecord(fromShortName("CO2"), "838", 49.001675, 2),
  IsotopeRecord(fromShortName("CO2"), "CKD241"),
  IsotopeRecord(fromShortName("CO2"), "CKDMT100"),
  IsotopeRecord(fromShortName("CO2"), "CKDMT252"),
  IsotopeRecord(fromShortName("CO2"), "ForeignContHo66"),
  IsotopeRecord(fromShortName("CO2"), "ForeignContPWR93"),
  IsotopeRecord(fromShortName("CO2"), "SelfContHo66"),
  IsotopeRecord(fromShortName("CO2"), "SelfContPWR93"),
  /** Carbon dioxide species **/
  
  /** Ozone species **/
  deal_with_spec(Ozone)
  IsotopeRecord(fromShortName("O3"), "666", 47.984745, 1),
  IsotopeRecord(fromShortName("O3"), "667", 48.988960, 6),
  IsotopeRecord(fromShortName("O3"), "668", 49.988991, 1),
  IsotopeRecord(fromShortName("O3"), "676", 48.988960, 6),
  IsotopeRecord(fromShortName("O3"), "686", 49.988991, 1),
  /** Ozone species **/
  
  /** N2O species **/
  deal_with_spec(NitrogenOxide)
  IsotopeRecord(fromShortName("N2O"), "446", 44.001062, 9),
  IsotopeRecord(fromShortName("N2O"), "447", 45.005278, 54),
  IsotopeRecord(fromShortName("N2O"), "448", 46.005308, 9),
  IsotopeRecord(fromShortName("N2O"), "456", 44.998096, 6),
  IsotopeRecord(fromShortName("N2O"), "546", 44.998096, 6),
  /** N2O species **/
  
  /** CO species **/
  deal_with_spec(CarbonMonoxide)
  IsotopeRecord(fromShortName("CO"), "26", 27.994915, 1),
  IsotopeRecord(fromShortName("CO"), "27", 28.999130, 6),
  IsotopeRecord(fromShortName("CO"), "28", 29.999161, 1),
  IsotopeRecord(fromShortName("CO"), "36", 28.998270, 2),
  IsotopeRecord(fromShortName("CO"), "37", 30.002485, 12),
  IsotopeRecord(fromShortName("CO"), "38", 31.002516, 2),
  /** CO species **/
  
  /** CH4 species **/
  deal_with_spec(Methane)
  IsotopeRecord(fromShortName("CH4"), "211", 16.031300, 1),
  IsotopeRecord(fromShortName("CH4"), "212", 17.037475, 3),
  IsotopeRecord(fromShortName("CH4"), "311", 17.034655, 2),
  IsotopeRecord(fromShortName("CH4"), "312", 18.040830, 6),
  /** CH4 species **/
  
  /** Oxygen species **/
  deal_with_spec(Oxygen)
  IsotopeRecord(fromShortName("O2"), "66", 31.989830, 1),
  IsotopeRecord(fromShortName("O2"), "67", 32.994045, 6),
  IsotopeRecord(fromShortName("O2"), "68", 33.994076, 1),
  IsotopeRecord(fromShortName("O2"), "CIAfunCKDMT100"),
  IsotopeRecord(fromShortName("O2"), "MPM2020"),
  IsotopeRecord(fromShortName("O2"), "MPM85"),
  IsotopeRecord(fromShortName("O2"), "MPM87"),
  IsotopeRecord(fromShortName("O2"), "MPM89"),
  IsotopeRecord(fromShortName("O2"), "MPM92"),
  IsotopeRecord(fromShortName("O2"), "MPM93"),
  IsotopeRecord(fromShortName("O2"), "PWR88"),
  IsotopeRecord(fromShortName("O2"), "PWR93"),
  IsotopeRecord(fromShortName("O2"), "PWR98"),
  IsotopeRecord(fromShortName("O2"), "SelfContMPM93"),
  IsotopeRecord(fromShortName("O2"), "SelfContPWR93"),
  IsotopeRecord(fromShortName("O2"), "SelfContStandardType"),
  IsotopeRecord(fromShortName("O2"), "TRE05"),
  IsotopeRecord(fromShortName("O2"), "v0v0CKDMT100"),
  IsotopeRecord(fromShortName("O2"), "v1v0CKDMT100"),
  IsotopeRecord(fromShortName("O2"), "visCKDMT252"),
  /** Oxygen species **/
  
  /** NO species **/
  deal_with_spec(NitricOxide)
  IsotopeRecord(fromShortName("NO"), "46", 29.997989, 3),
  IsotopeRecord(fromShortName("NO"), "48", 32.002234, 3),
  IsotopeRecord(fromShortName("NO"), "56", 30.995023, 2),
  /** NO species **/
  
  /** SO2 species **/
  deal_with_spec(SulfurDioxide)
  IsotopeRecord(fromShortName("SO2"), "626", 63.961901, 1),
  IsotopeRecord(fromShortName("SO2"), "628", 66),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("SO2"), "636", 65),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("SO2"), "646", 65.957695, 1),
  /** SO2 species **/
  
  /** NO2 species **/
  deal_with_spec(NitrogenDioxide)
  IsotopeRecord(fromShortName("NO2"), "646", 45.992904, 3),
  IsotopeRecord(fromShortName("NO2"), "656", 46.989938, 2),
  /** NO2 species **/
  
  /** NH3 species **/
  deal_with_spec(Ammonia)
  IsotopeRecord(fromShortName("NH3"), "4111", 17.026549, 3),
  IsotopeRecord(fromShortName("NH3"), "4112", 18),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("NH3"), "5111", 18.023583, 2),
  /** NH3 species **/
  
  /** HNO3 species **/
  deal_with_spec(NitricAcid)
  IsotopeRecord(fromShortName("HNO3"), "146", 62.995644, 6),
  IsotopeRecord(fromShortName("HNO3"), "156", 63.992680, 4),
  /** HNO3 species **/
  
  /** OH species **/
  deal_with_spec(Hydroxyl)
  IsotopeRecord(fromShortName("OH"), "61", 17.002740, 2),
  IsotopeRecord(fromShortName("OH"), "62", 18.008915, 3),
  IsotopeRecord(fromShortName("OH"), "81", 19.006986, 2),
  /** OH species **/
  
  /** HF species **/
  deal_with_spec(HydrogenFluoride)
  IsotopeRecord(fromShortName("HF"), "19", 20.006229, 4),
  IsotopeRecord(fromShortName("HF"), "29", 21.012404, 6),
  /** HF species **/
  
  /** HCl species **/
  deal_with_spec(HydrogenChloride)
  IsotopeRecord(fromShortName("HCl"), "15", 35.976678, 8),
  IsotopeRecord(fromShortName("HCl"), "17", 37.973729, 8),
  IsotopeRecord(fromShortName("HCl"), "25", 36.982853, 12),
  IsotopeRecord(fromShortName("HCl"), "27", 38.979904, 12),
  /** HCl species **/
  
  /** HBr species **/
  deal_with_spec(HydrogenBromide)
  IsotopeRecord(fromShortName("HBr"), "11", 81.924115, 8),
  IsotopeRecord(fromShortName("HBr"), "19", 79.926160, 8),
  IsotopeRecord(fromShortName("HBr"), "21", 82.930289, 12),
  IsotopeRecord(fromShortName("HBr"), "29", 80.932336, 12),
  /** HBr species **/
  
  /** HI species **/
  deal_with_spec(HydrogenIodide)
  IsotopeRecord(fromShortName("HI"), "17", 127.912297, 12),
  IsotopeRecord(fromShortName("HI"), "27", 128.918472, 18),
  /** HI species **/
  
  /** ClO species **/
  deal_with_spec(ChlorineMonoxide)
  IsotopeRecord(fromShortName("ClO"), "56", 50.963768, 4),
  IsotopeRecord(fromShortName("ClO"), "76", 52.960819, 4),
  /** ClO species **/
  
  /** OCS species **/
  deal_with_spec(CarbonylSulfide)
  IsotopeRecord(fromShortName("OCS"), "622", 59.966986, 1),
  IsotopeRecord(fromShortName("OCS"), "623", 60.966371, 4),
  IsotopeRecord(fromShortName("OCS"), "624", 61.962780, 1),
  IsotopeRecord(fromShortName("OCS"), "632", 60.970341, 2),
  IsotopeRecord(fromShortName("OCS"), "634", 62.966137, 2),
  IsotopeRecord(fromShortName("OCS"), "822", 61.971231, 1),
  /** OCS species **/
  
  /** H2CO species **/
  deal_with_spec(Formaldehyde)
  IsotopeRecord(fromShortName("H2CO"), "126", 30.010565, 1),
  IsotopeRecord(fromShortName("H2CO"), "128", 32.014811, 1),
  IsotopeRecord(fromShortName("H2CO"), "136", 31.013920, 2),
  /** H2CO species **/
  
  /** HDCO species nb. If the order D matters, rename this to indicate how **/
  deal_with_spec(HeavyFormaldehyde)
  IsotopeRecord(fromShortName("HDCO"), "26", 31),  // FIXME: Better mass and some gj?  What is the AFGL code???
  /** HDCO species **/
  
  /** D2CO species **/
  deal_with_spec(VeryHeavyFormaldehyde)
  IsotopeRecord(fromShortName("D2CO"), "26", 32),  // FIXME: Better mass and some gj?  What is the AFGL code???
  /** D2CO species **/
  
  /** HOCl species **/
  deal_with_spec(HypochlorousAcid)
  IsotopeRecord(fromShortName("HOCl"), "165", 51.971593, 8),
  IsotopeRecord(fromShortName("HOCl"), "167", 53.968644, 8),
  /** HOCl species **/
  
  /** N2 species **/
  deal_with_spec(Nitrogen)
  IsotopeRecord(fromShortName("N2"), "44", 28.006148, 1),
  IsotopeRecord(fromShortName("N2"), "45", 29.003182, 6),
  IsotopeRecord(fromShortName("N2"), "CIAfunCKDMT100"),
  IsotopeRecord(fromShortName("N2"), "CIAfunCKDMT252"),
  IsotopeRecord(fromShortName("N2"), "CIArotCKDMT100"),
  IsotopeRecord(fromShortName("N2"), "CIArotCKDMT252"),
  IsotopeRecord(fromShortName("N2"), "DryContATM01"),
  IsotopeRecord(fromShortName("N2"), "SelfContBorysow"),
  IsotopeRecord(fromShortName("N2"), "SelfContMPM93"),
  IsotopeRecord(fromShortName("N2"), "SelfContPWR93"),
  IsotopeRecord(fromShortName("N2"), "SelfContStandardType"),
  /** N2 species **/
  
  /** HCN species **/
  deal_with_spec(HydrogenCyanide)
  IsotopeRecord(fromShortName("HCN"), "124", 27.010899, 6),
  IsotopeRecord(fromShortName("HCN"), "125", 28.007933, 4),
  IsotopeRecord(fromShortName("HCN"), "134", 28.014254, 12),
  IsotopeRecord(fromShortName("HCN"), "224", 28),  // FIXME: Better mass and some gj?
  /** HCN species **/
  
  /** CH3Cl species **/
  deal_with_spec(Chloromethane)
  IsotopeRecord(fromShortName("CH3Cl"), "215", 49.992328, 4),
  IsotopeRecord(fromShortName("CH3Cl"), "217", 51.989379, 4),
  /** CH3Cl species **/
  
  /** H2O2 species **/
  deal_with_spec(HydrogenPeroxide)
  IsotopeRecord(fromShortName("H2O2"), "1661", 34.005480, 1),
  /** H2O2 species **/
  
  /** C2H2 species **/
  deal_with_spec(Acetylene)
  IsotopeRecord(fromShortName("C2H2"), "1221", 26.015650, 1),
  IsotopeRecord(fromShortName("C2H2"), "1222", 27.021825, 6),
  IsotopeRecord(fromShortName("C2H2"), "1231", 27.019005, 8),
  /** C2H2 species **/
  
  /** C2H6 species **/
  deal_with_spec(Ethane)
  IsotopeRecord(fromShortName("C2H6"), "1221", 30.046950, 1),
  IsotopeRecord(fromShortName("C2H6"), "1231", 31.050305, 2),
  /** C2H6 species **/
  
  /** PH3 species **/
  deal_with_spec(Phosphine)
  IsotopeRecord(fromShortName("PH3"), "1111", 33.997238, 2),
  /** PH3 species **/
  
  /** COF2 species **/
  deal_with_spec(CarbonylFluoride)
  IsotopeRecord(fromShortName("COF2"), "269", 65.991722, 1),
  IsotopeRecord(fromShortName("COF2"), "369", 66.995083, 2),
  /** COF2 species **/
  
  /** SF6 species **/
  deal_with_spec(SulfurHexafluoride)
  IsotopeRecord(fromShortName("SF6"), "29", 145.962492, 1),
  /** SF6 species **/
  
  /** H2S species **/
  deal_with_spec(HydrogenSulfide)
  IsotopeRecord(fromShortName("H2S"), "121", 33.987721, 1),
  IsotopeRecord(fromShortName("H2S"), "122", 35),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("H2S"), "131", 34.987105, 4),
  IsotopeRecord(fromShortName("H2S"), "141", 35.983515, 1),
  /** H2S species **/
  
  /** HCOOH species **/
  deal_with_spec(FormicAcid)
  IsotopeRecord(fromShortName("HCOOH"), "126", 46.005480, 4),
  IsotopeRecord(fromShortName("HCOOH"), "136", 47),  // FIXME: Better mass and some gj?
  /** HCOOH species **/
  
  /** DCOOH species **/
  deal_with_spec(LeftHeavyFormicAcid)
  IsotopeRecord(fromShortName("DCOOH"), "266", 47),  // FIXME: Better mass and some gj?  What is the AFGL code???
  /** DCOOH species **/
  
  /** HCOOD species **/
  deal_with_spec(RightHeavyFormicAcid)
  IsotopeRecord(fromShortName("HCOOD"), "266", 47),  // FIXME: Better mass and some gj?  What is the AFGL code???
  /** HCOOD species **/
  
  /** HO2 species **/
  deal_with_spec(Hydroperoxyl)
  IsotopeRecord(fromShortName("HO2"), "166", 32.997655, 2),
  /** HO2 species **/
  
  /** O species **/
  deal_with_spec(OxygenAtom)
  IsotopeRecord(fromShortName("O"), "6", 15.994915, 1),
  /** O species **/
  
  /** ClONO2 species **/
  deal_with_spec(ChlorineNitrate)
  IsotopeRecord(fromShortName("ClONO2"), "5646", 96.956672, 12),
  IsotopeRecord(fromShortName("ClONO2"), "7646", 98.953723, 12),
  /** ClONO2 species **/
  
  /** NO+ species **/
  deal_with_spec(NitricOxideCation)
  IsotopeRecord(fromShortName("NO+"), "46", 29.997989, 3),
  /** NO+ species **/
  
  /** OClO species **/
  deal_with_spec(ChlorineDioxide)
  IsotopeRecord(fromShortName("OClO"), "656", 67),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("OClO"), "676", 69),  // FIXME: Better mass and some gj?
  /** OClO species **/
  
  /** BrO species **/
  deal_with_spec(BromineMonoxide)
  IsotopeRecord(fromShortName("BrO"), "16", 97),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("BrO"), "96", 95),  // FIXME: Better mass and some gj?
  /** BrO species **/
  
  /** H2SO4 species **/
  deal_with_spec(SulfuricAcid)
  IsotopeRecord(fromShortName("H2SO4"), "126", 98),  // FIXME: Better mass and some gj?
  /** H2SO4 species **/
  
  /** Cl2O2 species **/
  deal_with_spec(ChlorineMonoxideDimer)
  IsotopeRecord(fromShortName("Cl2O2"), "565", 102),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("Cl2O2"), "765", 104),  // FIXME: Better mass and some gj?
  /** Cl2O2 species **/
  
  /** HOBr species **/
  deal_with_spec(HypobromousAcid)
  IsotopeRecord(fromShortName("HOBr"), "161", 97.919027, 8),
  IsotopeRecord(fromShortName("HOBr"), "169", 95.921076, 8),
  /** HOBr species **/
  
  /** C2H4 species **/
  deal_with_spec(Ethylene)
  IsotopeRecord(fromShortName("C2H4"), "221", 28.031300, 1),
  IsotopeRecord(fromShortName("C2H4"), "231", 29.034655, 2),
  /** C2H4 species **/
  
  /** CH3OH species **/
  deal_with_spec(Methanol)
  IsotopeRecord(fromShortName("CH3OH"), "2161", 32.026215, 2),
  /** CH3OH species **/
  
  /** CH3Br species **/
  deal_with_spec(Bromomethane)
  IsotopeRecord(fromShortName("CH3Br"), "211", 95.939764, 4),
  IsotopeRecord(fromShortName("CH3Br"), "219", 93.941811, 4),
  /** CH3Br species **/
  
  /** CH3CN species **/
  deal_with_spec(Acetonitrile)
  IsotopeRecord(fromShortName("CH3CN"), "2124", 41.026549, 3),
  IsotopeRecord(fromShortName("CH3CN"), "2125", 42),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("CH3CN"), "2134", 42),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("CH3CN"), "3124", 42),  // FIXME: Better mass and some gj?
  /** CH3CN species **/
  
  /** CH2DCN species nb. If the order D matters, rename this to indicate how **/
  deal_with_spec(HeavyAcetonitrile)
  IsotopeRecord(fromShortName("CH2DCN"), "224", 42),  // FIXME: Better mass and some gj?  What is the AFGL code???
  /** CH2DCN species **/
  
  /** CF4 species **/
  deal_with_spec(CarbonTetrafluoride)
  IsotopeRecord(fromShortName("CF4"), "29", 87.993616, 1),
  /** CF4 species **/
  
  /** C4H2 species **/
  deal_with_spec(Diacetylene)
  IsotopeRecord(fromShortName("C4H2"), "2211", 50.015650, 1),
  /** C4H2 species **/
  
  /** HC3N species **/
  deal_with_spec(Cyanoacetylene)
  IsotopeRecord(fromShortName("HC3N"), "12224", 51.010899, 6),
  IsotopeRecord(fromShortName("HC3N"), "12225", 52),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HC3N"), "12234", 52),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HC3N"), "12324", 52),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HC3N"), "13224", 52),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HC3N"), "22224", 52),  // FIXME: Better mass and some gj?
  /** HC3N species **/
  
  /** H2 species **/
  deal_with_spec(Hydrogen)
  IsotopeRecord(fromShortName("H2"), "11", 2.015650, 1),
  IsotopeRecord(fromShortName("H2"), "12", 3.021825, 6),
  /** H2 species **/
  
  /** CS species **/
  deal_with_spec(CarbonMonosulfide)
  IsotopeRecord(fromShortName("CS"), "22", 43.971036, 1),
  IsotopeRecord(fromShortName("CS"), "23", 44.970399, 4),
  IsotopeRecord(fromShortName("CS"), "24", 45.966787, 1),
  IsotopeRecord(fromShortName("CS"), "32", 44.974368, 2),
  /** CS species **/
  
  /** HNC species **/
  deal_with_spec(HydrogenIsocyanide)
  IsotopeRecord(fromShortName("HNC"), "142", 27),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HNC"), "143", 28),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HNC"), "152", 28),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HNC"), "242", 28),  // FIXME: Better mass and some gj?
  /** HNC species **/
  
  /** SO species **/
  deal_with_spec(SulfurMonoxide)
  IsotopeRecord(fromShortName("SO"), "26", 48),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("SO"), "28", 50),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("SO"), "46", 50),  // FIXME: Better mass and some gj?
  /** SO species **/
  
  /** C3H8 species **/
  deal_with_spec(Propane)
  IsotopeRecord(fromShortName("C3H8"), "21", 54),  // FIXME: Better mass and some gj?
  /** C3H8 species **/
  
  /** H species **/
  deal_with_spec(HydrogenAtom)
  IsotopeRecord(fromShortName("H"), "1", 1),  // FIXME: Better mass and some gj?
  /** H species **/
  
  /** He species **/
  deal_with_spec(Helium)
  IsotopeRecord(fromShortName("He"), "4", 4),  // FIXME: Better mass and some gj?
  /** He species **/
  
  /** Ar species **/
  deal_with_spec(Argon)
  IsotopeRecord(fromShortName("Ar"), "8", 18),  // FIXME: Better mass and some gj?
  /** Ar species **/
  
  /** SO3 species **/
  deal_with_spec(SulfurTrioxide)
  IsotopeRecord(fromShortName("SO3"), "26", 79.956820, 1),
  /** SO3 species **/
  
  /** C2N2 species **/
  deal_with_spec(Cyanogen)
  IsotopeRecord(fromShortName("C2N2"), "4224", 52.006148, 1),
  /** C2N2 species **/
  
  /** COCl2 species **/
  deal_with_spec(Phosgene)
  IsotopeRecord(fromShortName("COCl2"), "2655", 97.932620, 1),
  IsotopeRecord(fromShortName("COCl2"), "2657", 99.929670, 16),
  /** COCl2 species **/
  
  /** CS2 species **/
  deal_with_spec(CarbonDisulfide)
  IsotopeRecord(fromShortName("CS2"), "222", 75.944140, 1),
  IsotopeRecord(fromShortName("CS2"), "223", 76.943256, 4),
  IsotopeRecord(fromShortName("CS2"), "224", 77.939940, 1),
  IsotopeRecord(fromShortName("CS2"), "232", 76.947495, 2),
  /** CS2 species **/
  
  /** All species need a default joker **/
  deal_with_spec(Methyl)
  deal_with_spec(Cyclopropene)
  deal_with_spec(Hexafluoroethane)
  deal_with_spec(Perfluoropropane)
  deal_with_spec(Perfluorobutane)
  deal_with_spec(Perfluoropentane)
  deal_with_spec(Perfluorohexane)
  deal_with_spec(Perfluorooctane)
  deal_with_spec(Perfluorocyclobutane)
  deal_with_spec(CarbonTetrachloride)
  deal_with_spec(CFC11)
  deal_with_spec(CFC113)
  deal_with_spec(CFC114)
  deal_with_spec(CFC115)
  deal_with_spec(CFC12)
  deal_with_spec(Dichloromethane)
  deal_with_spec(Trichloroethane)
  deal_with_spec(Trichloromethane)
  deal_with_spec(Bromochlorodifluoromethane)
  deal_with_spec(Bromotrifluoromethane)
  deal_with_spec(Dibromotetrafluoroethane)
  deal_with_spec(HCFC141b)
  deal_with_spec(HCFC142b)
  deal_with_spec(HCFC22)
  deal_with_spec(HFC125)
  deal_with_spec(HFC134a)
  deal_with_spec(HFC143a)
  deal_with_spec(HFC152a)
  deal_with_spec(HFC227ea)
  deal_with_spec(HFC23)
  deal_with_spec(HFC245fa)
  deal_with_spec(HFC32)
  deal_with_spec(NitrogenTrifluoride)
  deal_with_spec(SulfurylFluoride)
  deal_with_spec(HFC4310mee)
  deal_with_spec(Germane)
  deal_with_spec(Iodomethane)
  deal_with_spec(Fluoromethane)
  /** All species need a default joker **/
  
  /** Model species **/
  deal_with_spec(liquidcloud)
  IsotopeRecord(Species::liquidcloud, "ELL07"),
  IsotopeRecord(Species::liquidcloud, "MPM93"),
  deal_with_spec(icecloud)
  IsotopeRecord(Species::icecloud, "MPM93"),
  deal_with_spec(rain)
  IsotopeRecord(Species::rain, "MPM93"),
  deal_with_spec(free_electrons)
  deal_with_spec(particles)
  /** Model species **/
};

#undef deal_with_spec

constexpr std::array<std::size_t, std::size_t(Species::FINAL)+1> start_positions() noexcept {
  std::array<bool, std::size_t(Species::FINAL)> found{};
  for (auto& x: found) x = false;
  
  std::array<std::size_t, std::size_t(Species::FINAL)+1> out{};
  for (auto& x: out) x = Isotopologues.size();
  
  for (std::size_t i=0; i<Isotopologues.size(); i++) {
    const std::size_t ind = std::size_t(Isotopologues[i].spec);
    if (not found[ind]) {
      found[ind] = true;
      out[ind] = i;
    }
  }
  return out;
}

constexpr auto IsotopologuesStart = start_positions();

template <Species spec>
constexpr std::size_t count_isotopologues() noexcept {
  return IsotopologuesStart[std::size_t(spec) + 1] - IsotopologuesStart[std::size_t(spec)];
}

template <Species spec>
constexpr std::array<IsotopeRecord, count_isotopologues<spec>()> isotopologues() noexcept {
  static_assert(count_isotopologues<spec>() not_eq 0, "All species must be defined in the Isotopologues!");
  std::array<IsotopeRecord, count_isotopologues<spec>()> isots;
  for (std::size_t i=0; i<count_isotopologues<spec>(); i++) {
    isots[i] = Isotopologues[i + IsotopologuesStart[std::size_t(spec)]];
  }
  return isots;
}

Array<IsotopeRecord> isotopologues(Species spec);

constexpr Index find_species_index(const Species spec,
                                   const std::string_view isot) noexcept {
  for (std::size_t i=IsotopologuesStart[std::size_t(spec)]; i<IsotopologuesStart[std::size_t(spec) + 1]; i++) {
    if (isot == Isotopologues[i].isotname) {
      return i;
    }
  }
  return -1;
}


constexpr Index find_species_index(const IsotopeRecord ir) noexcept {
  return find_species_index(ir.spec, ir.isotname);
}

constexpr Index find_species_index(const std::string_view spec,
                                   const std::string_view isot) noexcept {
  return find_species_index(fromShortName(spec), isot);
}

constexpr const IsotopeRecord& select(Species spec, const std::string_view isotname) noexcept {
  return Isotopologues[find_species_index(spec, isotname)];
}

constexpr const IsotopeRecord& select_joker(Species spec) noexcept {
  return select(spec, Joker);
}

constexpr const IsotopeRecord& select_joker(std::string_view spec) noexcept {
  return select(fromShortName(spec), Joker);
}

String isotopologues_names(Species spec);

constexpr bool is_predefined_model(const IsotopeRecord& ir) noexcept {
  return not (nonstd::isdigit(ir.isotname[0]) or ir.isotname == Joker);
}

String predefined_model_names() noexcept;

constexpr bool same_or_joker(const IsotopeRecord& ir1, const IsotopeRecord& ir2) noexcept {
  if (ir1.spec not_eq ir2.spec) return false;
  if (ir1.joker() or ir2.joker()) return true;
  else return ir1.isotname == ir2.isotname;
}

struct IsotopologueRatios {
  static constexpr Index maxsize = Index(Isotopologues.size());
  std::array<Numeric, maxsize> data;
  
  constexpr IsotopologueRatios() noexcept : data() {
    for (auto& x: data) x = std::numeric_limits<Numeric>::quiet_NaN();
  }
  
  constexpr Numeric operator[](const Index spec_ind) const ARTS_NOEXCEPT {
    ARTS_ASSERT(spec_ind < maxsize and spec_ind >= 0)
    return data[spec_ind];
  }
  
  constexpr Numeric operator[](const IsotopeRecord& ir) const {
    const Index spec_ind = find_species_index(ir);
    ARTS_USER_ERROR_IF(spec_ind >= maxsize and spec_ind < 0,
                       "Cannot understand: ", ir.FullName(), " as a valid species")
    return data[spec_ind];
  }
  
  friend std::ostream& operator<<(std::ostream& os, const IsotopologueRatios& iso_rat) {
    for (size_t i=0; i<iso_rat.maxsize; i++) {
      if (i not_eq 0)
        os << '\n';
      os << Isotopologues[i].FullName() << ' ' << iso_rat.data[i];
    }
    return os;
  }
  
  constexpr bool all_isotopes_have_a_value() const noexcept {
    for (Index i=0; i<maxsize; i++) {
      if (not is_predefined_model(Isotopologues[i]) and not Isotopologues[i].joker() and nonstd::isnan(data[i])) {
        return false;
      }
    }
    return true;
  }
};

constexpr IsotopologueRatios isotopologue_ratiosInitFromBuiltin() {
  IsotopologueRatios isotopologue_ratios;
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("H2O", ISOT)] = VAL
  set_isot_val("161", .997317E+00);
  set_isot_val("181", 1.99983E-03);
  set_isot_val("171", 3.71884E-04);
  set_isot_val("162", 3.10693E-04);
  set_isot_val("182", 6.23003E-07);
  set_isot_val("172", 1.15853E-07);
  set_isot_val("262", 2.41970E-08);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("CO2", ISOT)] = VAL
  set_isot_val("626", .984204E+00);
  set_isot_val("636", 1.10574E-02);
  set_isot_val("628", 3.94707E-03);
  set_isot_val("627", 7.33989E-04);
  set_isot_val("638", 4.43446E-05);
  set_isot_val("637", 8.24623E-06);
  set_isot_val("828", 3.95734E-06);
  set_isot_val("827", 1.47180E-06);
  set_isot_val("727", 1.36847E-07);
  set_isot_val("838", 4.44600E-08);
  set_isot_val("837", 1.65354E-08);
  set_isot_val("737", 1.537500E-08);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("O3", ISOT)] = VAL
  set_isot_val("666", .992901E+00);
  set_isot_val("668", 3.98194E-03);
  set_isot_val("686", 1.99097E-03);
  set_isot_val("667", 7.40475E-04);
  set_isot_val("676", 3.70237E-04);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("N2O", ISOT)] = VAL
  set_isot_val("446", .990333E+00);
  set_isot_val("456", 3.64093E-03);
  set_isot_val("546", 3.64093E-03);
  set_isot_val("448", 1.98582E-03);
  set_isot_val("447", 3.69280E-04);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("CO", ISOT)] = VAL
  set_isot_val("26", .986544E+00);
  set_isot_val("36", 1.10836E-02);
  set_isot_val("28", 1.97822E-03);
  set_isot_val("27", 3.67867E-04);
  set_isot_val("38", 2.22250E-05);
  set_isot_val("37", 4.13292E-06);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("CH4", ISOT)] = VAL
  set_isot_val("211", .988274E+00);
  set_isot_val("311", 1.11031E-02);
  set_isot_val("212", 6.15751E-04);
  set_isot_val("312", 6.91785E-06);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("O2", ISOT)] = VAL
  set_isot_val("66", .995262E+00);
  set_isot_val("68", 3.99141E-03);
  set_isot_val("67", 7.42235E-04);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("NO", ISOT)] = VAL
  set_isot_val("46", .993974E+00);
  set_isot_val("56", 3.65431E-03);
  set_isot_val("48", 1.99312E-03);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("SO2", ISOT)] = VAL
  set_isot_val("626", .945678E+00);
  set_isot_val("646", 4.19503E-02);
  set_isot_val("636", 0.0074989421);
  set_isot_val("628", 0.0020417379);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("NO2", ISOT)] = VAL
  set_isot_val("646", .991616E+00);
  set_isot_val("656", 3.64564E-03);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("NH3", ISOT)] = VAL
  set_isot_val("4111", .995872E+00);
  set_isot_val("5111", 3.66129E-03);
  set_isot_val("4112", 0.00044792294);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HNO3", ISOT)] = VAL
  set_isot_val("146", .989110E+00);
  set_isot_val("156", 3.63600E-03);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("OH", ISOT)] = VAL
  set_isot_val("61", .997473E+00);
  set_isot_val("81", 2.00014E-03);
  set_isot_val("62", 1.55371E-04);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HF", ISOT)] = VAL
  set_isot_val("19", .999844E+00);
  set_isot_val("29", 1.55741E-04);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HCl", ISOT)] = VAL
  set_isot_val("15", .757587E+00);
  set_isot_val("17", .242257E+00);
  set_isot_val("25", 1.18005E-04);
  set_isot_val("27", 3.77350E-05);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HBr", ISOT)] = VAL
  set_isot_val("19", .506781E+00);
  set_isot_val("11", .493063E+00);
  set_isot_val("29", 7.89384E-05);
  set_isot_val("21", 7.68016E-05);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HI", ISOT)] = VAL
  set_isot_val("17", .999844E+00);
  set_isot_val("27", 1.55741E-04);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("ClO", ISOT)] = VAL
  set_isot_val("56", .755908E+00);
  set_isot_val("76", .241720E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("OCS", ISOT)] = VAL
  set_isot_val("622", .937395E+00);
  set_isot_val("624", 4.15828E-02);
  set_isot_val("632", 1.05315E-02);
  set_isot_val("623", 7.39908E-03);
  set_isot_val("822", 1.87967E-03);
  set_isot_val("634", 4.67508E-04);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("H2CO", ISOT)] = VAL
  set_isot_val("126", .986237E+00);
  set_isot_val("136", 1.10802E-02);
  set_isot_val("128", 1.97761E-03);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HDCO", ISOT)] = VAL
  set_isot_val("26", 0.00029578940);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("D2CO", ISOT)] = VAL
  set_isot_val("26", 2.2181076E-08);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HOCl", ISOT)] = VAL
  set_isot_val("165", .755790E+00);
  set_isot_val("167", .241683E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("N2", ISOT)] = VAL
  set_isot_val("44", .992687E+00);
  set_isot_val("45", 7.47809E-03);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HCN", ISOT)] = VAL
  set_isot_val("124", .985114E+00);
  set_isot_val("134", 1.10676E-02);
  set_isot_val("125", 3.62174E-03);
  set_isot_val("224", 0.00014773545);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("CH3Cl", ISOT)] = VAL
  set_isot_val("215", .748937E+00);
  set_isot_val("217", .239491E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("H2O2", ISOT)] = VAL
  set_isot_val("1661", .994952E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("C2H2", ISOT)] = VAL
  set_isot_val("1221", .977599E+00);
  set_isot_val("1231", 2.19663E-02);
  set_isot_val("1222", 3.04550E-04);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("C2H6", ISOT)] = VAL
  set_isot_val("1221", .976990E+00);
  set_isot_val("1231", 2.19526E-02);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("PH3", ISOT)] = VAL
  set_isot_val("1111", .999533E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("COF2", ISOT)] = VAL
  set_isot_val("269", .986544E+00);
  set_isot_val("369", 1.10834E-02);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("SF6", ISOT)] = VAL
  set_isot_val("29", .950180E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("H2S", ISOT)] = VAL
  set_isot_val("121", .949884E+00);
  set_isot_val("141", 4.21369E-02);
  set_isot_val("131", 7.49766E-03);
  set_isot_val("122", 0.00029991625);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HCOOH", ISOT)] = VAL
  set_isot_val("126", .983898E+00);
  set_isot_val("136", 0.010913149);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("DCOOH", ISOT)] = VAL
  set_isot_val("266", 0.00014755369);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HCOOD", ISOT)] = VAL
  set_isot_val("266", 0.00014755369);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HO2", ISOT)] = VAL
  set_isot_val("166", .995107E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("O", ISOT)] = VAL
  set_isot_val("6", .997628E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("ClONO2", ISOT)] = VAL
  set_isot_val("5646", .749570E+00);
  set_isot_val("7646", .239694E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("NO+", ISOT)] = VAL
  set_isot_val("46", .993974E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("OClO", ISOT)] = VAL
  set_isot_val("656", 0.75509223);
  set_isot_val("676", 0.24490632);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("BrO", ISOT)] = VAL
  set_isot_val("96", 0.50582466);
  set_isot_val("16", 0.49431069);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("H2SO4", ISOT)] = VAL
  set_isot_val("126", 0.95060479);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("Cl2O2", ISOT)] = VAL
  set_isot_val("565", 0.57016427);
  set_isot_val("765", 0.36982818);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HOBr", ISOT)] = VAL
  set_isot_val("169", .505579E+00);
  set_isot_val("161", .491894E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("C2H4", ISOT)] = VAL
  set_isot_val("221", .977294E+00);
  set_isot_val("231", .219595E-01);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("CH3OH", ISOT)] = VAL
  set_isot_val("2161", .985930E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("CH3Br", ISOT)] = VAL
  set_isot_val("219", .500995E+00);
  set_isot_val("211", .487433E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("CH3CN", ISOT)] = VAL
  set_isot_val("2124", .973866E+00);
  set_isot_val("3124", .102683e-01);
  set_isot_val("2134", .102683e-01);
  set_isot_val("2125", .347136e-02);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("CH2DCN", ISOT)] = VAL
  set_isot_val("224", .441185e-03);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("CF4", ISOT)] = VAL
  set_isot_val("29", .988890E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HC3N", ISOT)] = VAL
  set_isot_val("12224", .963346E+00);
  set_isot_val("12234", .106852e-01);
  set_isot_val("12324", .106852e-01);
  set_isot_val("13224", .106852e-01);
  set_isot_val("12225", .356272e-02);
  set_isot_val("22224", .144472e-03);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("CS", ISOT)] = VAL
  set_isot_val("22", .939624E+00);
  set_isot_val("24", .416817E-01);
  set_isot_val("32", .105565E-01);
  set_isot_val("23", .741668E-02);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("HNC", ISOT)] = VAL
  set_isot_val("142", .985280e+00);
  set_isot_val("143", .109285e-01);
  set_isot_val("152", .364384e-02);
  set_isot_val("242", .147761e-03);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("SO", ISOT)] = VAL
  set_isot_val("26", .950605e+00);
  set_isot_val("46", .420727e-01);
  set_isot_val("28", .194089e-02);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("C3H8", ISOT)] = VAL
  set_isot_val("21", 9.66290e-01);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("H2", ISOT)] = VAL
  set_isot_val("11", .999688E+00);
  set_isot_val("12", 3.11432E-04);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("H", ISOT)] = VAL
  set_isot_val("1", 1.00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("He", ISOT)] = VAL
  set_isot_val("4", 1.00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("Ar", ISOT)] = VAL
  set_isot_val("8", 1.00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("C4H2", ISOT)] = VAL
  set_isot_val("2211", .955998E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("SO3", ISOT)] = VAL
  set_isot_val("26", .943400E+00);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("CS2", ISOT)] = VAL
  set_isot_val("222", 8.92811E-01);
  set_isot_val("224", 7.92600E-02);
  set_isot_val("223", 1.40940E-02);
  set_isot_val("232", 1.03100E-02);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("C2N2", ISOT)] = VAL
  set_isot_val("4224", 9.70752E-01);
  #undef set_isot_val
  
  #define set_isot_val(ISOT, VAL) isotopologue_ratios.data[find_species_index("COCl2", ISOT)] = VAL
  set_isot_val("2655", 5.66392E-01);
  set_isot_val("2657", 5.66392E-01);
  #undef set_isot_val
  
  return isotopologue_ratios;
}

/*! Computes the mean mass for all defined isotopes of the species with mass and isotopologue ratio
 * 
 * \f[ 
 * m = \frac{ \sum_i r_i m_i }{ \sum_i r_i }
 * \f]
 * 
 * @param[in] spec A species
 * @param[in] ir All isotopologue ratios
 * @return mean mass
 */
constexpr Numeric mean_mass(Species spec, const IsotopologueRatios& ir) noexcept {
  Numeric sum_rm=0;
  Numeric sum_r =0;
  for (std::size_t i=IsotopologuesStart[std::size_t(spec)]; i <IsotopologuesStart[std::size_t(spec) + 1]; i++) {
    if (not nonstd::isnan(Isotopologues[i].mass) and not nonstd::isnan(ir[i])) {
      sum_rm += ir[i] * Isotopologues[i].mass;
      sum_r += ir[i];
    }
  }
  if (sum_r not_eq 0) return sum_rm / sum_r;
  else return 0 * std::numeric_limits<Numeric>::signaling_NaN();
}

/** Updates the name of the isotopologue based on
 * updates of the isotopologues.
 * 
 * This should only be invoked by versioned code as it is
 * not very efficient.
 * 
 * @param[in] old_name A valid isotopologue name in any version of ARTS
 * @return A name that is valid and equivalent in ARTS today (or a copy of old_name)
 */
String update_isot_name(const String& old_name);

constexpr bool all_have_ratio(const Species spec, const IsotopologueRatios& ir) noexcept {
  for (std::size_t i=IsotopologuesStart[std::size_t(spec)]; i<IsotopologuesStart[std::size_t(spec) + 1]; i++) {
    if (not Isotopologues[i].joker() and not is_predefined_model(Isotopologues[i]) and nonstd::isnan(ir[i])) {
      return false;
    }
  }
  return true;
}

std::pair<ArrayOfString, ArrayOfString> names_of_have_and_havenot_ratio(const Species spec, const IsotopologueRatios& ir) noexcept;
}  // Species

using SpeciesIsotopeRecord = Species::IsotopeRecord;

using ArrayOfIsotopeRecord = Array<SpeciesIsotopeRecord>;

using ArrayOfSpecies = Array<Species::Species>;

using SpeciesIsotopologueRatios = Species::IsotopologueRatios;

#endif  // isotopologues_h
