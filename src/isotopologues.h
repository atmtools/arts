#ifndef isotopologues_h
#define isotopologues_h

#include "mystring.h"
#include "nonstd.h"
#include "species.h"

namespace Species {
constexpr std::string_view Joker = "*";

/** Struct containing all information needed about one isotope */
struct IsotopeRecord {
  Species spec;
  std::string_view isotname;
  Numeric mass;
  Index gi;
  constexpr explicit IsotopeRecord(Species spec_, const std::string_view isotname_=Joker, Numeric mass_=std::numeric_limits<Numeric>::quiet_NaN(), Index gi_=-1) noexcept
  : spec(spec_), isotname(isotname_), mass(mass_), gi(gi_) {}
  constexpr IsotopeRecord() noexcept : spec(Species::FINAL), isotname(""), mass(std::numeric_limits<Numeric>::signaling_NaN()), gi(-1) {}
  friend std::ostream& operator<<(std::ostream& os, const IsotopeRecord& ir) {
    return os << ir.spec << ' ' << ir.isotname << ' ' << ir.mass << ' ' << ir.gi;
  }
  constexpr bool operator==(const IsotopeRecord& that) const noexcept {
    return that.spec == spec and that.isotname == isotname;
  }
  
  //! A comparison with pure named string (this is not an exact comparison)
  constexpr bool operator==(const std::string_view specstr) const noexcept {
    auto namepos = specstr.find(toShortName(spec));
    auto isotpos = specstr.find(isotname);
    return isotpos not_eq specstr.npos and namepos not_eq specstr.npos;
  }
  
  template <typename T> constexpr bool operator!=(T x)const noexcept {return not operator==(x);}
  
  String FullName() const noexcept {return String(toShortName(spec)) + String("-") + String(isotname);}
  constexpr bool joker() const noexcept {return isotname == Joker;}
};

#define deal_with_spec(SPEC) IsotopeRecord(Species::SPEC),

constexpr std::array Isotopologues {
  /** Water species **/
  IsotopeRecord(fromShortName("H2O"), "161", 18.010565, 1),
  IsotopeRecord(fromShortName("H2O"), "181", 20.014811, 1),
  IsotopeRecord(fromShortName("H2O"), "171", 19.014780, 6),
  IsotopeRecord(fromShortName("H2O"), "162", 19.016740, 6),
  IsotopeRecord(fromShortName("H2O"), "182", 21.020985, 6),
  IsotopeRecord(fromShortName("H2O"), "172", 20.020956, 36),
  IsotopeRecord(fromShortName("H2O"), "262", 20.022915, 1),
  IsotopeRecord(fromShortName("H2O"), "SelfContStandardType"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContStandardType"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContMaTippingType"),
  IsotopeRecord(fromShortName("H2O"), "ContMPM93"),
  IsotopeRecord(fromShortName("H2O"), "SelfContCKDMT100"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContCKDMT100"),
  IsotopeRecord(fromShortName("H2O"), "SelfContCKDMT252"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContCKDMT252"),
  IsotopeRecord(fromShortName("H2O"), "SelfContCKDMT320"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContCKDMT320"),
  IsotopeRecord(fromShortName("H2O"), "SelfContCKD222"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContCKD222"),
  IsotopeRecord(fromShortName("H2O"), "SelfContCKD242"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContCKD242"),
  IsotopeRecord(fromShortName("H2O"), "SelfContCKD24"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContCKD24"),
  IsotopeRecord(fromShortName("H2O"), "ForeignContATM01"),
  IsotopeRecord(fromShortName("H2O"), "CP98"),
  IsotopeRecord(fromShortName("H2O"), "MPM87"),
  IsotopeRecord(fromShortName("H2O"), "MPM89"),
  IsotopeRecord(fromShortName("H2O"), "MPM93"),
  IsotopeRecord(fromShortName("H2O"), "PWR98"),
  /** Water species **/
  
  /** Carbon dioxide species **/
  IsotopeRecord(fromShortName("CO2"), "626", 43.989830, 1),
  IsotopeRecord(fromShortName("CO2"), "636", 44.993185, 2),
  IsotopeRecord(fromShortName("CO2"), "628", 45.994076, 1),
  IsotopeRecord(fromShortName("CO2"), "627", 44.994045, 6),
  IsotopeRecord(fromShortName("CO2"), "638", 46.997431, 2),
  IsotopeRecord(fromShortName("CO2"), "637", 45.997400, 12),
  IsotopeRecord(fromShortName("CO2"), "828", 47.998322, 1),
  IsotopeRecord(fromShortName("CO2"), "728", 46.998291, 6),
  IsotopeRecord(fromShortName("CO2"), "727", 45.998262, 1),
  IsotopeRecord(fromShortName("CO2"), "838", 49.001675, 2),
  IsotopeRecord(fromShortName("CO2"), "837", 48.001646, 12),
  IsotopeRecord(fromShortName("CO2"), "737", 47.001618, 2),
  IsotopeRecord(fromShortName("CO2"), "CKD241"),
  IsotopeRecord(fromShortName("CO2"), "CKDMT100"),
  IsotopeRecord(fromShortName("CO2"), "CKDMT252"),
  IsotopeRecord(fromShortName("CO2"), "SelfContPWR93"),
  IsotopeRecord(fromShortName("CO2"), "ForeignContPWR93"),
  IsotopeRecord(fromShortName("CO2"), "SelfContHo66"),
  IsotopeRecord(fromShortName("CO2"), "ForeignContHo66"),
  /** Carbon dioxide species **/
  
  /** Ozone species **/
  IsotopeRecord(fromShortName("O3"), "666", 47.984745, 1),
  IsotopeRecord(fromShortName("O3"), "668", 49.988991, 1),
  IsotopeRecord(fromShortName("O3"), "686", 49.988991, 1),
  IsotopeRecord(fromShortName("O3"), "667", 48.988960, 6),
  IsotopeRecord(fromShortName("O3"), "676", 48.988960, 6),
  /** Ozone species **/
  
  /** N2O species **/
  IsotopeRecord(fromShortName("N2O"), "446", 44.001062, 9),
  IsotopeRecord(fromShortName("N2O"), "456", 44.998096, 6),
  IsotopeRecord(fromShortName("N2O"), "546", 44.998096, 6),
  IsotopeRecord(fromShortName("N2O"), "448", 46.005308, 9),
  IsotopeRecord(fromShortName("N2O"), "447", 45.005278, 54),
  /** N2O species **/
  
  /** CO species **/
  IsotopeRecord(fromShortName("CO"), "26", 27.994915, 1),
  IsotopeRecord(fromShortName("CO"), "36", 28.998270, 2),
  IsotopeRecord(fromShortName("CO"), "28", 29.999161, 1),
  IsotopeRecord(fromShortName("CO"), "27", 28.999130, 6),
  IsotopeRecord(fromShortName("CO"), "38", 31.002516, 2),
  IsotopeRecord(fromShortName("CO"), "37", 30.002485, 12),
  /** CO species **/
  
  /** CH4 species **/
  IsotopeRecord(fromShortName("CH4"), "211", 16.031300, 1),
  IsotopeRecord(fromShortName("CH4"), "311", 17.034655, 2),
  IsotopeRecord(fromShortName("CH4"), "212", 17.037475, 3),
  IsotopeRecord(fromShortName("CH4"), "312", 18.040830, 6),
  /** CH4 species **/
  
  /** Oxygen species **/
  IsotopeRecord(fromShortName("O2"), "66", 31.989830, 1),
  IsotopeRecord(fromShortName("O2"), "68", 33.994076, 1),
  IsotopeRecord(fromShortName("O2"), "67", 32.994045, 6),
  IsotopeRecord(fromShortName("O2"), "CIAfunCKDMT100"),
  IsotopeRecord(fromShortName("O2"), "v0v0CKDMT100"),
  IsotopeRecord(fromShortName("O2"), "v1v0CKDMT100"),
  IsotopeRecord(fromShortName("O2"), "visCKDMT252"),
  IsotopeRecord(fromShortName("O2"), "SelfContStandardType"),
  IsotopeRecord(fromShortName("O2"), "SelfContMPM93"),
  IsotopeRecord(fromShortName("O2"), "SelfContPWR93"),
  IsotopeRecord(fromShortName("O2"), "PWR98"),
  IsotopeRecord(fromShortName("O2"), "PWR93"),
  IsotopeRecord(fromShortName("O2"), "PWR88"),
  IsotopeRecord(fromShortName("O2"), "MPM93"),
  IsotopeRecord(fromShortName("O2"), "TRE05"),
  IsotopeRecord(fromShortName("O2"), "MPM92"),
  IsotopeRecord(fromShortName("O2"), "MPM89"),
  IsotopeRecord(fromShortName("O2"), "MPM87"),
  IsotopeRecord(fromShortName("O2"), "MPM85"),
  IsotopeRecord(fromShortName("O2"), "MPM2020"),
  /** Oxygen species **/
  
  /** NO species **/
  IsotopeRecord(fromShortName("NO"), "46", 29.997989, 3),
  IsotopeRecord(fromShortName("NO"), "56", 30.995023, 2),
  IsotopeRecord(fromShortName("NO"), "48", 32.002234, 3),
  /** NO species **/
  
  /** SO2 species **/
  IsotopeRecord(fromShortName("SO2"), "626", 63.961901, 1),
  IsotopeRecord(fromShortName("SO2"), "646", 65.957695, 1),
  IsotopeRecord(fromShortName("SO2"), "636", 65),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("SO2"), "628", 66),  // FIXME: Better mass and some gj?
  /** SO2 species **/
  
  /** NO2 species **/
  IsotopeRecord(fromShortName("NO2"), "646", 45.992904, 3),
  IsotopeRecord(fromShortName("NO2"), "656", 46.989938, 2),
  /** NO2 species **/
  
  /** NH3 species **/
  IsotopeRecord(fromShortName("NH3"), "4111", 17.026549, 3),
  IsotopeRecord(fromShortName("NH3"), "5111", 18.023583, 2),
  IsotopeRecord(fromShortName("NH3"), "4112", 18),  // FIXME: Better mass and some gj?
  /** NH3 species **/
  
  /** HNO3 species **/
  IsotopeRecord(fromShortName("HNO3"), "146", 62.995644, 6),
  IsotopeRecord(fromShortName("HNO3"), "156", 63.992680, 4),
  /** HNO3 species **/
  
  /** OH species **/
  IsotopeRecord(fromShortName("OH"), "61", 17.002740, 2),
  IsotopeRecord(fromShortName("OH"), "81", 19.006986, 2),
  IsotopeRecord(fromShortName("OH"), "62", 18.008915, 3),
  /** OH species **/
  
  /** HF species **/
  IsotopeRecord(fromShortName("HF"), "19", 20.006229, 4),
  IsotopeRecord(fromShortName("HF"), "29", 21.012404, 6),
  /** HF species **/
  
  /** HCl species **/
  IsotopeRecord(fromShortName("HCl"), "15", 35.976678, 8),
  IsotopeRecord(fromShortName("HCl"), "17", 37.973729, 8),
  IsotopeRecord(fromShortName("HCl"), "25", 36.982853, 12),
  IsotopeRecord(fromShortName("HCl"), "27", 38.979904, 12),
  /** HCl species **/
  
  /** HBr species **/
  IsotopeRecord(fromShortName("HBr"), "19", 79.926160, 8),
  IsotopeRecord(fromShortName("HBr"), "11", 81.924115, 8),
  IsotopeRecord(fromShortName("HBr"), "29", 80.932336, 12),
  IsotopeRecord(fromShortName("HBr"), "21", 82.930289, 12),
  /** HBr species **/
  
  /** HI species **/
  IsotopeRecord(fromShortName("HI"), "17", 127.912297, 12),
  IsotopeRecord(fromShortName("HI"), "27", 128.918472, 18),
  /** HI species **/
  
  /** ClO species **/
  IsotopeRecord(fromShortName("ClO"), "56", 50.963768, 4),
  IsotopeRecord(fromShortName("ClO"), "76", 52.960819, 4),
  /** ClO species **/
  
  /** OCS species **/
  IsotopeRecord(fromShortName("OCS"), "622", 59.966986, 1),
  IsotopeRecord(fromShortName("OCS"), "624", 61.962780, 1),
  IsotopeRecord(fromShortName("OCS"), "632", 60.970341, 2),
  IsotopeRecord(fromShortName("OCS"), "623", 60.966371, 4),
  IsotopeRecord(fromShortName("OCS"), "822", 61.971231, 1),
  IsotopeRecord(fromShortName("OCS"), "634", 62.966137, 2),
  /** OCS species **/
  
  /** H2CO species **/
  IsotopeRecord(fromShortName("H2CO"), "1126", 30.010565, 1),
  IsotopeRecord(fromShortName("H2CO"), "1136", 31.013920, 2),
  IsotopeRecord(fromShortName("H2CO"), "1128", 32.014811, 1),
  IsotopeRecord(fromShortName("H2CO"), "1226", 31),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("H2CO"), "2226", 32),  // FIXME: Better mass and some gj?
  /** H2CO species **/
  
  /** HOCl species **/
  IsotopeRecord(fromShortName("HOCl"), "165", 51.971593, 8),
  IsotopeRecord(fromShortName("HOCl"), "167", 53.968644, 8),
  /** HOCl species **/
  
  /** N2 species **/
  IsotopeRecord(fromShortName("N2"), "44", 28.006148, 1),
  IsotopeRecord(fromShortName("N2"), "45", 29.003182, 6),
  IsotopeRecord(fromShortName("N2"), "SelfContMPM93"),
  IsotopeRecord(fromShortName("N2"), "SelfContPWR93"),
  IsotopeRecord(fromShortName("N2"), "SelfContStandardType"),
  IsotopeRecord(fromShortName("N2"), "SelfContBorysow"),
  IsotopeRecord(fromShortName("N2"), "CIArotCKDMT100"),
  IsotopeRecord(fromShortName("N2"), "CIAfunCKDMT100"),
  IsotopeRecord(fromShortName("N2"), "CIArotCKDMT252"),
  IsotopeRecord(fromShortName("N2"), "CIAfunCKDMT252"),
  IsotopeRecord(fromShortName("N2"), "DryContATM01"),
  /** N2 species **/
  
  /** HCN species **/
  IsotopeRecord(fromShortName("HCN"), "124", 27.010899, 6),
  IsotopeRecord(fromShortName("HCN"), "134", 28.014254, 12),
  IsotopeRecord(fromShortName("HCN"), "125", 28.007933, 4),
  IsotopeRecord(fromShortName("HCN"), "224", 28),  // FIXME: Better mass and some gj?
  /** HCN species **/
  
  /** CH3Cl species **/
  IsotopeRecord(fromShortName("CH3Cl"), "215", 49.992328, 4),
  IsotopeRecord(fromShortName("CH3Cl"), "217", 51.989379, 4),
  /** CH3Cl species **/
  
  /** H2O2 species **/
  IsotopeRecord(fromShortName("H2O2"), "1661", 34.005480, 1),
  /** H2O2 species **/
  
  /** C2H2 species **/
  IsotopeRecord(fromShortName("C2H2"), "1221", 26.015650, 1),
  IsotopeRecord(fromShortName("C2H2"), "1231", 27.019005, 8),
  IsotopeRecord(fromShortName("C2H2"), "1222", 27.021825, 6),
  /** C2H2 species **/
  
  /** C2H6 species **/
  IsotopeRecord(fromShortName("C2H6"), "1221", 30.046950, 1),
  IsotopeRecord(fromShortName("C2H6"), "1231", 31.050305, 2),
  /** C2H6 species **/
  
  /** PH3 species **/
  IsotopeRecord(fromShortName("PH3"), "1111", 33.997238, 2),
  /** PH3 species **/
  
  /** COF2 species **/
  IsotopeRecord(fromShortName("COF2"), "269", 65.991722, 1),
  IsotopeRecord(fromShortName("COF2"), "369", 66.995083, 2),
  /** COF2 species **/
  
  /** SF6 species **/
  IsotopeRecord(fromShortName("SF6"), "29", 145.962492, 1),
  /** SF6 species **/
  
  /** H2S species **/
  IsotopeRecord(fromShortName("H2S"), "121", 33.987721, 1),
  IsotopeRecord(fromShortName("H2S"), "141", 35.983515, 1),
  IsotopeRecord(fromShortName("H2S"), "131", 34.987105, 4),
  IsotopeRecord(fromShortName("H2S"), "122", 35),  // FIXME: Better mass and some gj?
  /** H2S species **/
  
  /** HCOOH species **/
  IsotopeRecord(fromShortName("HCOOH"), "1261", 46.005480, 4),
  IsotopeRecord(fromShortName("HCOOH"), "1361", 47),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HCOOH"), "2261", 47),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HCOOH"), "1262", 47),  // FIXME: Better mass and some gj?
  /** HCOOH species **/
  
  /** HO2 species **/
  IsotopeRecord(fromShortName("HO2"), "166", 32.997655, 2),
  /** HO2 species **/
  
  /** O species **/
  IsotopeRecord(fromShortName("O"), "6", 15.994915, 1),
  /** O species **/
  
  /** ClONO2 species **/
  IsotopeRecord(fromShortName("ClONO2"), "5646", 96.956672, 12),
  IsotopeRecord(fromShortName("ClONO2"), "7646", 98.953723, 12),
  /** ClONO2 species **/
  
  /** NO+ species **/
  IsotopeRecord(fromShortName("NO+"), "46", 29.997989, 3),
  /** NO+ species **/
  
  /** OClO species **/
  IsotopeRecord(fromShortName("OClO"), "656", 67),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("OClO"), "676", 69),  // FIXME: Better mass and some gj?
  /** OClO species **/
  
  /** BrO species **/
  IsotopeRecord(fromShortName("BrO"), "96", 95),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("BrO"), "16", 97),  // FIXME: Better mass and some gj?
  /** BrO species **/
  
  /** H2SO4 species **/
  IsotopeRecord(fromShortName("H2SO4"), "126", 98),  // FIXME: Better mass and some gj?
  /** H2SO4 species **/
  
  /** Cl2O2 species **/
  IsotopeRecord(fromShortName("Cl2O2"), "565", 102),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("Cl2O2"), "765", 104),  // FIXME: Better mass and some gj?
  /** Cl2O2 species **/
  
  /** HOBr species **/
  IsotopeRecord(fromShortName("HOBr"), "169", 95.921076, 8),
  IsotopeRecord(fromShortName("HOBr"), "161", 97.919027, 8),
  /** HOBr species **/
  
  /** C2H4 species **/
  IsotopeRecord(fromShortName("C2H4"), "221", 28.031300, 1),
  IsotopeRecord(fromShortName("C2H4"), "231", 29.034655, 2),
  /** C2H4 species **/
  
  /** CH3OH species **/
  IsotopeRecord(fromShortName("CH3OH"), "2161", 32.026215, 2),
  /** CH3OH species **/
  
  /** CH3Br species **/
  IsotopeRecord(fromShortName("CH3Br"), "219", 93.941811, 4),
  IsotopeRecord(fromShortName("CH3Br"), "211", 95.939764, 4),
  /** CH3Br species **/
  
  /** CH3CN species **/
  IsotopeRecord(fromShortName("CH3CN"), "211124", 41.026549, 3),
  IsotopeRecord(fromShortName("CH3CN"), "311124", 42),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("CH3CN"), "211134", 42),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("CH3CN"), "211125", 42),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("CH3CN"), "211224", 42),  // FIXME: Better mass and some gj?
  /** CH3CN species **/
  
  /** CF4 species **/
  IsotopeRecord(fromShortName("CF4"), "29", 87.993616, 1),
  /** CF4 species **/
  
  /** C4H2 species **/
  IsotopeRecord(fromShortName("C4H2"), "2211", 50.015650, 1),
  /** C4H2 species **/
  
  /** HC3N species **/
  IsotopeRecord(fromShortName("HC3N"), "12224", 51.010899, 6),
  IsotopeRecord(fromShortName("HC3N"), "12234", 52),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HC3N"), "12324", 52),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HC3N"), "13224", 52),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HC3N"), "12225", 52),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HC3N"), "22224", 52),  // FIXME: Better mass and some gj?
  /** HC3N species **/
  
  /** H2 species **/
  IsotopeRecord(fromShortName("H2"), "11", 2.015650, 1),
  IsotopeRecord(fromShortName("H2"), "12", 3.021825, 6),
  /** H2 species **/
  
  /** CS species **/
  IsotopeRecord(fromShortName("CS"), "22", 43.971036, 1),
  IsotopeRecord(fromShortName("CS"), "24", 45.966787, 1),
  IsotopeRecord(fromShortName("CS"), "32", 44.974368, 2),
  IsotopeRecord(fromShortName("CS"), "23", 44.970399, 4),
  /** CS species **/
  
  /** HNC species **/
  IsotopeRecord(fromShortName("HNC"), "142", 27),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HNC"), "143", 28),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HNC"), "152", 28),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("HNC"), "242", 28),  // FIXME: Better mass and some gj?
  /** HNC species **/
  
  /** SO species **/
  IsotopeRecord(fromShortName("SO"), "26", 48),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("SO"), "46", 50),  // FIXME: Better mass and some gj?
  IsotopeRecord(fromShortName("SO"), "28", 50),  // FIXME: Better mass and some gj?
  /** SO species **/
  
  /** C3H8 species **/
  IsotopeRecord(fromShortName("C3H8"), "21", 54),  // FIXME: Better mass and some gj?
  /** C3H8 species **/
  
  /** H species **/
  IsotopeRecord(fromShortName("H"), "1", 1),  // FIXME: Better mass and some gj?
  /** H species **/
  
  /** He species **/
  IsotopeRecord(fromShortName("He"), "4", 4),  // FIXME: Better mass and some gj?
  /** He species **/
  
  /** Ar species **/
  IsotopeRecord(fromShortName("Ar"), "8", 18),  // FIXME: Better mass and some gj?
  /** Ar species **/
  
  /** SO3 species **/
  IsotopeRecord(fromShortName("SO3"), "26", 79.956820, 1),
  /** SO3 species **/
  
  /** C2N2 species **/
  IsotopeRecord(fromShortName("C2N2"), "4224", 52.006148, 1),
  /** C2N2 species **/
  
  /** COCl2 species **/
  IsotopeRecord(fromShortName("COCl2"), "2655", 97.932620, 1),
  IsotopeRecord(fromShortName("COCl2"), "2657", 99.929670, 16),
  /** COCl2 species **/
  
  /** CS2 species **/
  IsotopeRecord(fromShortName("CS2"), "222", 75.944140, 1),
  IsotopeRecord(fromShortName("CS2"), "224", 77.939940, 1),
  IsotopeRecord(fromShortName("CS2"), "223", 76.943256, 4),
  IsotopeRecord(fromShortName("CS2"), "232", 76.947495, 2),
  /** CS2 species **/
  
  /** Model species **/
  IsotopeRecord(Species::liquidcloud, "MPM93"),
  IsotopeRecord(Species::liquidcloud, "ELL07"),
  IsotopeRecord(Species::icecloud, "MPM93"),
  IsotopeRecord(Species::rain, "MPM93"),
  /** Model species **/
  
  /** All species need a default joker **/
  deal_with_spec(Water)
  deal_with_spec(CarbonDioxide)
  deal_with_spec(Ozone)
  deal_with_spec(NitrogenOxide)
  deal_with_spec(CarbonMonoxide)
  deal_with_spec(Methane)
  deal_with_spec(Oxygen)
  deal_with_spec(NitricOxide)
  deal_with_spec(SulfurDioxide)
  deal_with_spec(NitrogenDioxide)
  deal_with_spec(Ammonia)
  deal_with_spec(NitricAcid)
  deal_with_spec(Hydroxyl)
  deal_with_spec(HydrogenFluoride)
  deal_with_spec(HydrogenChloride)
  deal_with_spec(HydrogenBromide)
  deal_with_spec(HydrogenIodide)
  deal_with_spec(ChlorineMonoxide)
  deal_with_spec(CarbonylSulfide)
  deal_with_spec(Formaldehyde)
  deal_with_spec(HypochlorousAcid)
  deal_with_spec(Nitrogen)
  deal_with_spec(HydrogenCyanide)
  deal_with_spec(MethylChloride)
  deal_with_spec(HydrogenPeroxide)
  deal_with_spec(Acetylene)
  deal_with_spec(Ethane)
  deal_with_spec(Phosphine)
  deal_with_spec(CarbonylFluoride)
  deal_with_spec(SulfurHexafluoride)
  deal_with_spec(HydrogenSulfide)
  deal_with_spec(FormicAcid)
  deal_with_spec(Hydroperoxyl)
  deal_with_spec(OxygenAtom)
  deal_with_spec(ChlorineNitrate)
  deal_with_spec(NitricOxideCation)
  deal_with_spec(HypobromousAcid)
  deal_with_spec(Ethylene)
  deal_with_spec(Methanol)
  deal_with_spec(MethylBromide)
  deal_with_spec(Acetonitrile)
  deal_with_spec(CarbonTetrafluoride)
  deal_with_spec(Diacetylene)
  deal_with_spec(Cyanoacetylene)
  deal_with_spec(Hydrogen)
  deal_with_spec(CarbonMonosulfide)
  deal_with_spec(SulfurTrioxide)
  deal_with_spec(Cyanogen)
  deal_with_spec(Phosgene)
  deal_with_spec(SulfurMonoxide)
  deal_with_spec(CarbonDisulfide)
  deal_with_spec(Methyl)
  deal_with_spec(Cyclopropene)
  deal_with_spec(SulfuricAcid)
  deal_with_spec(HydrogenIsocyanide)
  deal_with_spec(BromineMonoxide)
  deal_with_spec(ChlorineDioxide)
  deal_with_spec(Propane)
  deal_with_spec(Helium)
  deal_with_spec(ChlorineMonoxideDimer)
  deal_with_spec(HydrogenAtom)
  deal_with_spec(Argon)
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
  deal_with_spec(liquidcloud)
  deal_with_spec(icecloud)
  deal_with_spec(rain)
  deal_with_spec(free_electrons)
  deal_with_spec(particles)
  /** All species need a default joker **/
};

#undef deal_with_spec

template <Species spec>
constexpr std::size_t count_isotopologues() noexcept {
  std::size_t n=0;
  for (auto& x: Isotopologues) if (x.spec == spec) ++n;
  return n;
}

template <Species spec>
constexpr std::array<IsotopeRecord, count_isotopologues<spec>()> isotopologues() noexcept {
  static_assert(count_isotopologues<spec>() not_eq 0, "All species must be defined in the Isotopologues!");
  std::array<IsotopeRecord, count_isotopologues<spec>()> isots;
  std::size_t i=0;
  for (auto& x: Isotopologues) if (x.spec == spec) isots[i++] = x;
  return isots;
}

Array<IsotopeRecord> isotopologues(Species spec);

constexpr Index find_species_index(const IsotopeRecord& ir) noexcept {
  for (std::size_t i=0; i<Isotopologues.size(); i++) {
    if (ir == Isotopologues[i]) return Index(i);
  }
  return -1;
}

constexpr Index find_species_index(const Species spec,
                                   const std::string_view isot) noexcept {
  return find_species_index(IsotopeRecord(spec, isot));
}

constexpr Index find_species_index(const std::string_view spec,
                                   const std::string_view isot) noexcept {
  return find_species_index(fromShortName(spec), isot);
}

constexpr Numeric first_mass(Species spec) noexcept {
  for (auto& x: Isotopologues) {
    if (spec == x.spec) return x.mass;
  }
  return 0 * std::numeric_limits<Numeric>::signaling_NaN();
}

constexpr const IsotopeRecord& select(Species spec, const std::string_view isotname) noexcept {
  return Isotopologues[find_species_index(IsotopeRecord(spec, isotname))];
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
}  // Species

using ArrayOfIsotopeRecord = Array<Species::IsotopeRecord>;

using ArrayOfSpecies = Array<Species::Species>;

#endif  // isotopologues_h
