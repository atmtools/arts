#pragma once

#include <compare.h>
#include <mystring.h>
#include <nonstd.h>
#include <species.h>
#include <xml_io_stream.h>

#include <array>
#include <compare>
#include <limits>
#include <string_view>

namespace Species {
inline constexpr std::string_view Joker = "*";

/** Struct containing all information needed about one isotope */
struct Isotope {
  SpeciesEnum spec{SpeciesEnum::Bath};

  //! A custom name that is unique for this Species type
  std::string_view isotname{Joker};

  //! The mass of the isotope in units of grams per mol.  It is Nan if not defined
  Numeric mass{std::numeric_limits<Numeric>::quiet_NaN()};

  //! The degeneracy of states of the molecule.  It is -1 if not defined.
  Index gi{-1};

  constexpr Isotope() = default;
  constexpr explicit Isotope(
      SpeciesEnum spec_,
      const std::string_view isotname_ = Joker,
      Numeric mass_ = std::numeric_limits<Numeric>::quiet_NaN(),
      Index gi_     = -1)
      : spec(spec_), isotname(isotname_), mass(mass_), gi(gi_) {}

  Isotope(const std::string_view);

  [[nodiscard]] constexpr bool is_joker() const { return isotname == Joker; }

  [[nodiscard]] constexpr bool is_predefined() const {
    return not(nonstd::isdigit(isotname[0]) or is_joker());
  }
  [[nodiscard]] constexpr bool is_normal() const {
    return nonstd::isdigit(isotname[0]) and not is_joker();
  }

  [[nodiscard]] constexpr bool OK() const { return good_enum(spec); }
  [[nodiscard]] String FullName() const;

  friend std::ostream& operator<<(std::ostream& os, const Isotope& ir);

  std::strong_ordering operator<=>(const Isotope&) const;
  bool operator==(const Isotope&) const;
  bool operator!=(const Isotope&) const;
};
}  // namespace Species

template <>
struct std::formatter<Species::Isotope> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Species::Isotope& v,
                              FmtContext& ctx) const {
    const std::string_view quote = tags.quote();
    return tags.format(ctx, quote, v.FullName(), quote);
  }
};

namespace Species {

#define deal_with_spec(SPEC) Isotope(SpeciesEnum::SPEC)

/** A list of all ARTS isotopologues, note how the species enum class input HAS to be sorted */
inline constexpr std::array Isotopologues{
    deal_with_spec(Bath),
    /** Water species **/
    deal_with_spec(Water),
    Isotope("H2O"_spec, "161", 18.010565, 1),
    Isotope("H2O"_spec, "162", 19.016740, 6),
    Isotope("H2O"_spec, "171", 19.014780, 6),
    Isotope("H2O"_spec, "172", 20.020956, 36),
    Isotope("H2O"_spec, "181", 20.014811, 1),
    Isotope("H2O"_spec, "182", 21.020985, 6),
    Isotope("H2O"_spec, "262", 20.022915, 1),
    Isotope("H2O"_spec, "ForeignContCKDMT320"),
    Isotope("H2O"_spec, "ForeignContCKDMT350"),
    Isotope("H2O"_spec, "ForeignContCKDMT400"),
    Isotope("H2O"_spec, "ForeignContStandardType"),
    Isotope("H2O"_spec, "MPM89"),
    Isotope("H2O"_spec, "PWR2021"),
    Isotope("H2O"_spec, "PWR2022"),
    Isotope("H2O"_spec, "PWR98"),
    Isotope("H2O"_spec, "SelfContCKDMT320"),
    Isotope("H2O"_spec, "SelfContCKDMT350"),
    Isotope("H2O"_spec, "SelfContCKDMT400"),
    Isotope("H2O"_spec, "SelfContStandardType"),
    /** Water species **/

    /** Carbon dioxide species **/
    deal_with_spec(CarbonDioxide),
    Isotope("CO2"_spec, "626", 43.989830, 1),
    Isotope("CO2"_spec, "627", 44.994045, 6),
    Isotope("CO2"_spec, "628", 45.994076, 1),
    Isotope("CO2"_spec, "636", 44.993185, 2),
    Isotope("CO2"_spec, "637", 45.997400, 12),
    Isotope("CO2"_spec, "638", 46.997431, 2),
    Isotope("CO2"_spec, "727", 45.998262, 1),
    Isotope("CO2"_spec, "737", 47.001618, 2),
    Isotope("CO2"_spec, "827", 46.998291, 6),
    Isotope("CO2"_spec, "828", 47.998322, 1),
    Isotope("CO2"_spec, "837", 48.001646, 12),
    Isotope("CO2"_spec, "838", 49.001675, 2),
    Isotope("CO2"_spec, "CKDMT252"),
    /** Carbon dioxide species **/

    /** Ozone species **/
    deal_with_spec(Ozone),
    Isotope("O3"_spec, "666", 47.984745, 1),
    Isotope("O3"_spec, "667", 48.988960, 6),
    Isotope("O3"_spec, "668", 49.988991, 1),
    Isotope("O3"_spec, "676", 48.988960, 6),
    Isotope("O3"_spec, "686", 49.988991, 1),
    /** Ozone species **/

    /** N2O species **/
    deal_with_spec(NitrogenOxide),
    Isotope("N2O"_spec, "446", 44.001062, 9),
    Isotope("N2O"_spec, "447", 45.005278, 54),
    Isotope("N2O"_spec, "448", 46.005308, 9),
    Isotope("N2O"_spec, "456", 44.998096, 6),
    Isotope("N2O"_spec, "546", 44.998096, 6),
    /** N2O species **/

    /** CO species **/
    deal_with_spec(CarbonMonoxide),
    Isotope("CO"_spec, "26", 27.994915, 1),
    Isotope("CO"_spec, "27", 28.999130, 6),
    Isotope("CO"_spec, "28", 29.999161, 1),
    Isotope("CO"_spec, "36", 28.998270, 2),
    Isotope("CO"_spec, "37", 30.002485, 12),
    Isotope("CO"_spec, "38", 31.002516, 2),
    /** CO species **/

    /** CH4 species **/
    deal_with_spec(Methane),
    Isotope("CH4"_spec, "211", 16.031300, 1),
    Isotope("CH4"_spec, "212", 17.037475, 3),
    Isotope("CH4"_spec, "311", 17.034655, 2),
    Isotope("CH4"_spec, "312", 18.040830, 6),
    /** CH4 species **/

    /** Oxygen species **/
    deal_with_spec(Oxygen),
    Isotope("O2"_spec, "66", 31.989830, 1),
    Isotope("O2"_spec, "67", 32.994045, 6),
    Isotope("O2"_spec, "68", 33.994076, 1),
    Isotope("O2"_spec, "CIAfunCKDMT100"),
    Isotope("O2"_spec, "MPM2020"),
    Isotope("O2"_spec, "MPM89"),
    Isotope("O2"_spec, "PWR2021"),
    Isotope("O2"_spec, "PWR2022"),
    Isotope("O2"_spec, "PWR98"),
    Isotope("O2"_spec, "SelfContStandardType"),
    Isotope("O2"_spec, "TRE05"),
    Isotope("O2"_spec, "v0v0CKDMT100"),
    Isotope("O2"_spec, "v1v0CKDMT100"),
    Isotope("O2"_spec, "visCKDMT252"),
    /** Oxygen species **/

    /** NO species **/
    deal_with_spec(NitricOxide),
    Isotope("NO"_spec, "46", 29.997989, 3),
    Isotope("NO"_spec, "48", 32.002234, 3),
    Isotope("NO"_spec, "56", 30.995023, 2),
    /** NO species **/

    /** SO2 species **/
    deal_with_spec(SulfurDioxide),
    Isotope("SO2"_spec, "626", 63.961901, 1),
    Isotope("SO2"_spec, "628", 65.966146, 1),
    Isotope("SO2"_spec, "636", 64.961286, 4),
    Isotope("SO2"_spec, "646", 65.957695, 1),
    /** SO2 species **/

    /** NO2 species **/
    deal_with_spec(NitrogenDioxide),
    Isotope("NO2"_spec, "646", 45.992904, 3),
    Isotope("NO2"_spec, "656", 46.989938, 2),
    /** NO2 species **/

    /** NH3 species **/
    deal_with_spec(Ammonia),
    Isotope("NH3"_spec, "4111", 17.026549, 3),
    Isotope("NH3"_spec, "4112", 18),  // FIXME: Better mass and some gj?
    Isotope("NH3"_spec, "5111", 18.023583, 2),
    /** NH3 species **/

    /** HNO3 species **/
    deal_with_spec(NitricAcid),
    Isotope("HNO3"_spec, "146", 62.995644, 6),
    Isotope("HNO3"_spec, "156", 63.992680, 4),
    /** HNO3 species **/

    /** OH species **/
    deal_with_spec(Hydroxyl),
    Isotope("OH"_spec, "61", 17.002740, 2),
    Isotope("OH"_spec, "62", 18.008915, 3),
    Isotope("OH"_spec, "81", 19.006986, 2),
    /** OH species **/

    /** HF species **/
    deal_with_spec(HydrogenFluoride),
    Isotope("HF"_spec, "19", 20.006229, 4),
    Isotope("HF"_spec, "29", 21.012404, 6),
    /** HF species **/

    /** HCl species **/
    deal_with_spec(HydrogenChloride),
    Isotope("HCl"_spec, "15", 35.976678, 8),
    Isotope("HCl"_spec, "17", 37.973729, 8),
    Isotope("HCl"_spec, "25", 36.982853, 12),
    Isotope("HCl"_spec, "27", 38.979904, 12),
    /** HCl species **/

    /** HBr species **/
    deal_with_spec(HydrogenBromide),
    Isotope("HBr"_spec, "11", 81.924115, 8),
    Isotope("HBr"_spec, "19", 79.926160, 8),
    Isotope("HBr"_spec, "21", 82.930289, 12),
    Isotope("HBr"_spec, "29", 80.932336, 12),
    /** HBr species **/

    /** HI species **/
    deal_with_spec(HydrogenIodide),
    Isotope("HI"_spec, "17", 127.912297, 12),
    Isotope("HI"_spec, "27", 128.918472, 18),
    /** HI species **/

    /** ClO species **/
    deal_with_spec(ChlorineMonoxide),
    Isotope("ClO"_spec, "56", 50.963768, 4),
    Isotope("ClO"_spec, "76", 52.960819, 4),
    /** ClO species **/

    /** OCS species **/
    deal_with_spec(CarbonylSulfide),
    Isotope("OCS"_spec, "622", 59.966986, 1),
    Isotope("OCS"_spec, "623", 60.966371, 4),
    Isotope("OCS"_spec, "624", 61.962780, 1),
    Isotope("OCS"_spec, "632", 60.970341, 2),
    Isotope("OCS"_spec, "634", 62.966137, 2),
    Isotope("OCS"_spec, "822", 61.971231, 1),
    /** OCS species **/

    /** H2CO species **/
    deal_with_spec(Formaldehyde),
    Isotope("H2CO"_spec, "126", 30.010565, 1),
    Isotope("H2CO"_spec, "128", 32.014811, 1),
    Isotope("H2CO"_spec, "136", 31.013920, 2),
    /** H2CO species **/

    /** HDCO species nb. If the order D matters, rename this to indicate how **/
    deal_with_spec(HeavyFormaldehyde),
    Isotope("HDCO"_spec,
            "26",
            31),  // FIXME: Better mass and some gj?  What is the AFGL code???
    /** HDCO species **/

    /** D2CO species **/
    deal_with_spec(VeryHeavyFormaldehyde),
    Isotope("D2CO"_spec,
            "26",
            32),  // FIXME: Better mass and some gj?  What is the AFGL code???
    /** D2CO species **/

    /** HOCl species **/
    deal_with_spec(HypochlorousAcid),
    Isotope("HOCl"_spec, "165", 51.971593, 8),
    Isotope("HOCl"_spec, "167", 53.968644, 8),
    /** HOCl species **/

    /** N2 species **/
    deal_with_spec(Nitrogen),
    Isotope("N2"_spec, "44", 28.006148, 1),
    Isotope("N2"_spec, "45", 29.003182, 6),
    Isotope("N2"_spec, "CIAfunCKDMT252"),
    Isotope("N2"_spec, "CIArotCKDMT252"),
    Isotope("N2"_spec, "SelfContMPM93"),
    Isotope("N2"_spec, "SelfContPWR2021"),
    Isotope("N2"_spec, "SelfContStandardType"),
    /** N2 species **/

    /** HCN species **/
    deal_with_spec(HydrogenCyanide),
    Isotope("HCN"_spec, "124", 27.010899, 6),
    Isotope("HCN"_spec, "125", 28.007933, 4),
    Isotope("HCN"_spec, "134", 28.014254, 12),
    Isotope("HCN"_spec, "224", 28),  // FIXME: Better mass and some gj?
    /** HCN species **/

    /** CH3Cl species **/
    deal_with_spec(Chloromethane),
    Isotope("CH3Cl"_spec, "215", 49.992328, 4),
    Isotope("CH3Cl"_spec, "217", 51.989379, 4),
    /** CH3Cl species **/

    /** H2O2 species **/
    deal_with_spec(HydrogenPeroxide),
    Isotope("H2O2"_spec, "1661", 34.005480, 1),
    /** H2O2 species **/

    /** C2H2 species **/
    deal_with_spec(Acetylene),
    Isotope("C2H2"_spec, "1221", 26.015650, 1),
    Isotope("C2H2"_spec, "1222", 27.021825, 6),
    Isotope("C2H2"_spec, "1231", 27.019005, 8),
    /** C2H2 species **/

    /** C2H6 species **/
    deal_with_spec(Ethane),
    Isotope("C2H6"_spec, "1221", 30.046950, 1),
    Isotope("C2H6"_spec, "1231", 31.050305, 2),
    /** C2H6 species **/

    /** PH3 species **/
    deal_with_spec(Phosphine),
    Isotope("PH3"_spec, "1111", 33.997238, 2),
    /** PH3 species **/

    /** COF2 species **/
    deal_with_spec(CarbonylFluoride),
    Isotope("COF2"_spec, "269", 65.991722, 1),
    Isotope("COF2"_spec, "369", 66.995083, 2),
    /** COF2 species **/

    /** SF6 species **/
    deal_with_spec(SulfurHexafluoride),
    Isotope("SF6"_spec, "29", 145.962492, 1),
    /** SF6 species **/

    /** H2S species **/
    deal_with_spec(HydrogenSulfide),
    Isotope("H2S"_spec, "121", 33.987721, 1),
    Isotope("H2S"_spec, "122", 35),  // FIXME: Better mass and some gj?
    Isotope("H2S"_spec, "131", 34.987105, 4),
    Isotope("H2S"_spec, "141", 35.983515, 1),
    /** H2S species **/

    /** HCOOH species **/
    deal_with_spec(FormicAcid),
    Isotope("HCOOH"_spec, "126", 46.005480, 4),
    Isotope("HCOOH"_spec,
            "136",
            47),  // FIXME: Better mass and some gj?
    /** HCOOH species **/

    /** DCOOH species **/
    deal_with_spec(LeftHeavyFormicAcid),
    Isotope("DCOOH"_spec,
            "266",
            47),  // FIXME: Better mass and some gj?  What is the AFGL code???
    /** DCOOH species **/

    /** HCOOD species **/
    deal_with_spec(RightHeavyFormicAcid),
    Isotope("HCOOD"_spec,
            "266",
            47),  // FIXME: Better mass and some gj?  What is the AFGL code???
    /** HCOOD species **/

    /** HO2 species **/
    deal_with_spec(Hydroperoxyl),
    Isotope("HO2"_spec, "166", 32.997655, 2),
    /** HO2 species **/

    /** O species **/
    deal_with_spec(OxygenAtom),
    Isotope("O"_spec, "6", 15.994915, 1),
    /** O species **/

    /** ClONO2 species **/
    deal_with_spec(ChlorineNitrate),
    Isotope("ClONO2"_spec, "5646", 96.956672, 12),
    Isotope("ClONO2"_spec, "7646", 98.953723, 12),
    /** ClONO2 species **/

    /** NO+ species **/
    deal_with_spec(NitricOxideCation),
    Isotope("NO+"_spec, "46", 29.997989, 3),
    /** NO+ species **/

    /** HOBr species **/
    deal_with_spec(HypobromousAcid),
    Isotope("HOBr"_spec, "161", 97.919027, 8),
    Isotope("HOBr"_spec, "169", 95.921076, 8),
    /** HOBr species **/

    /** C2H4 species **/
    deal_with_spec(Ethylene),
    Isotope("C2H4"_spec, "221", 28.031300, 1),
    Isotope("C2H4"_spec, "231", 29.034655, 2),
    /** C2H4 species **/

    /** CH3OH species **/
    deal_with_spec(Methanol),
    Isotope("CH3OH"_spec, "2161", 32.026215, 2),
    /** CH3OH species **/

    /** CH3Br species **/
    deal_with_spec(Bromomethane),
    Isotope("CH3Br"_spec, "211", 95.939764, 4),
    Isotope("CH3Br"_spec, "219", 93.941811, 4),
    /** CH3Br species **/

    /** CH3CN species **/
    deal_with_spec(Acetonitrile),
    Isotope("CH3CN"_spec, "2124", 41.026549, 3),
    Isotope("CH3CN"_spec,
            "2125",
            42),  // FIXME: Better mass and some gj?
    Isotope("CH3CN"_spec,
            "2134",
            42),  // FIXME: Better mass and some gj?
    Isotope("CH3CN"_spec,
            "3124",
            42),  // FIXME: Better mass and some gj?
    /** CH3CN species **/

    /** CH2DCN species nb. If the order D matters, rename this to indicate how **/
    deal_with_spec(HeavyAcetonitrile),
    Isotope("CH2DCN"_spec,
            "224",
            42),  // FIXME: Better mass and some gj?  What is the AFGL code???
    /** CH2DCN species **/

    /** CF4 species **/
    deal_with_spec(CarbonTetrafluoride),
    Isotope("CF4"_spec, "29", 87.993616, 1),
    /** CF4 species **/

    /** C4H2 species **/
    deal_with_spec(Diacetylene),
    Isotope("C4H2"_spec, "2211", 50.015650, 1),
    /** C4H2 species **/

    /** HC3N species **/
    deal_with_spec(Cyanoacetylene),
    Isotope("HC3N"_spec, "12224", 51.010899, 6),
    Isotope("HC3N"_spec,
            "12225",
            52),  // FIXME: Better mass and some gj?
    Isotope("HC3N"_spec,
            "12234",
            52),  // FIXME: Better mass and some gj?
    Isotope("HC3N"_spec,
            "12324",
            52),  // FIXME: Better mass and some gj?
    Isotope("HC3N"_spec,
            "13224",
            52),  // FIXME: Better mass and some gj?
    Isotope("HC3N"_spec,
            "22224",
            52),  // FIXME: Better mass and some gj?
    /** HC3N species **/

    /** H2 species **/
    deal_with_spec(Hydrogen),
    Isotope("H2"_spec, "11", 2.015650, 1),
    Isotope("H2"_spec, "12", 3.021825, 6),
    /** H2 species **/

    /** CS species **/
    deal_with_spec(CarbonMonosulfide),
    Isotope("CS"_spec, "22", 43.971036, 1),
    Isotope("CS"_spec, "23", 44.970399, 4),
    Isotope("CS"_spec, "24", 45.966787, 1),
    Isotope("CS"_spec, "32", 44.974368, 2),
    /** CS species **/

    /** SO3 species **/
    deal_with_spec(SulfurTrioxide),
    Isotope("SO3"_spec, "26", 79.956820, 1),
    /** SO3 species **/

    /** C2N2 species **/
    deal_with_spec(Cyanogen),
    Isotope("C2N2"_spec, "4224", 52.006148, 1),
    /** C2N2 species **/

    /** COCl2 species **/
    deal_with_spec(Phosgene),
    Isotope("COCl2"_spec, "2655", 97.932620, 1),
    Isotope("COCl2"_spec, "2657", 99.929670, 16),
    /** COCl2 species **/

    /** SO species **/
    deal_with_spec(SulfurMonoxide),
    Isotope("SO"_spec, "26", 47.966986, 1),
    Isotope("SO"_spec, "28", 49.971231, 1),
    Isotope("SO"_spec, "46", 49.962782, 1),
    /** SO species **/

    /** CS2 species **/
    deal_with_spec(CarbonDisulfide),
    Isotope("CS2"_spec, "222", 75.944140, 1),
    Isotope("CS2"_spec, "223", 76.943256, 4),
    Isotope("CS2"_spec, "224", 77.939940, 1),
    Isotope("CS2"_spec, "232", 76.947495, 2),
    /** CS2 species **/

    deal_with_spec(Methyl),
    deal_with_spec(Cyclopropene),

    /** H2SO4 species **/
    deal_with_spec(SulfuricAcid),
    Isotope("H2SO4"_spec,
            "126",
            98),  // FIXME: Better mass and some gj?
    /** H2SO4 species **/

    /** HNC species **/
    deal_with_spec(HydrogenIsocyanide),
    Isotope("HNC"_spec, "142", 27),  // FIXME: Better mass and some gj?
    Isotope("HNC"_spec, "143", 28),  // FIXME: Better mass and some gj?
    Isotope("HNC"_spec, "152", 28),  // FIXME: Better mass and some gj?
    Isotope("HNC"_spec, "242", 28),  // FIXME: Better mass and some gj?
    /** HNC species **/

    /** BrO species **/
    deal_with_spec(BromineMonoxide),
    Isotope("BrO"_spec, "16", 97),  // FIXME: Better mass and some gj?
    Isotope("BrO"_spec, "96", 95),  // FIXME: Better mass and some gj?
    /** BrO species **/

    /** OClO species **/
    deal_with_spec(ChlorineDioxide),
    Isotope("OClO"_spec, "656", 67),  // FIXME: Better mass and some gj?
    Isotope("OClO"_spec, "676", 69),  // FIXME: Better mass and some gj?
    /** OClO species **/

    /** C3H8 species **/
    deal_with_spec(Propane),
    Isotope("C3H8"_spec, "21", 54),  // FIXME: Better mass and some gj?
    /** C3H8 species **/

    /** He species **/
    deal_with_spec(Helium),
    Isotope("He"_spec, "4", 4),  // FIXME: Better mass and some gj?
    /** He species **/

    /** Cl2O2 species **/
    deal_with_spec(ChlorineMonoxideDimer),
    Isotope("Cl2O2"_spec,
            "565",
            102),  // FIXME: Better mass and some gj?
    Isotope("Cl2O2"_spec,
            "765",
            104),  // FIXME: Better mass and some gj?
    /** Cl2O2 species **/

    /** H species **/
    deal_with_spec(HydrogenAtom),
    Isotope("H"_spec, "1", 1),  // FIXME: Better mass and some gj?
    /** H species **/

    /** Ar species **/
    deal_with_spec(Argon),
    Isotope("Ar"_spec, "8", 39.948),  // FIXME: Better mass and some gj?
    /** Ar species **/

    deal_with_spec(Hexafluoroethane),
    deal_with_spec(Perfluoropropane),
    deal_with_spec(Perfluorobutane),
    deal_with_spec(Perfluoropentane),
    deal_with_spec(Perfluorohexane),
    deal_with_spec(Perfluorooctane),
    deal_with_spec(Perfluorocyclobutane),
    deal_with_spec(CarbonTetrachloride),
    deal_with_spec(CFC11),
    deal_with_spec(CFC113),
    deal_with_spec(CFC114),
    deal_with_spec(CFC115),
    deal_with_spec(CFC12),
    deal_with_spec(Dichloromethane),
    deal_with_spec(Trichloroethane),
    deal_with_spec(Trichloromethane),
    deal_with_spec(Bromochlorodifluoromethane),
    deal_with_spec(Bromotrifluoromethane),
    deal_with_spec(Dibromotetrafluoroethane),
    deal_with_spec(HCFC141b),
    deal_with_spec(HCFC142b),
    deal_with_spec(HCFC22),
    deal_with_spec(HFC125),
    deal_with_spec(HFC134a),
    deal_with_spec(HFC143a),
    deal_with_spec(HFC152a),
    deal_with_spec(HFC227ea),
    deal_with_spec(HFC23),
    deal_with_spec(HFC236fa),
    deal_with_spec(HFC245fa),
    deal_with_spec(HFC32),
    deal_with_spec(HFC365mfc),

    /** NF3 species **/
    deal_with_spec(NitrogenTrifluoride),
    Isotope("NF3"_spec, "4999", 70.998286, 3),
    /** NF3 species **/

    deal_with_spec(SulfurylFluoride),
    deal_with_spec(HFC4310mee),

    /** GeH4 species **/
    deal_with_spec(Germane),
    Isotope("GeH4"_spec, "011", 73.955550, 1),
    Isotope("GeH4"_spec, "211", 75.953380, 1),
    Isotope("GeH4"_spec, "311", 76.954764, 10),
    Isotope("GeH4"_spec, "411", 77.952479, 1),
    Isotope("GeH4"_spec, "611", 79.952703, 1),
    /** GeH4 species **/

    /** CH3I species **/
    deal_with_spec(Iodomethane),
    Isotope("CH3I"_spec, "217", 141.927947, 6),
    /** CH3I species **/

    /** CH3F species **/
    deal_with_spec(Fluoromethane),
    Isotope("CH3F"_spec, "219", 34.021878, 2),
    /** CH3F species **/

    /** Model species **/
    deal_with_spec(liquidcloud),
    Isotope(SpeciesEnum::liquidcloud, "ELL07"),
    deal_with_spec(icecloud),
    deal_with_spec(rain),
    deal_with_spec(free_electrons),
    deal_with_spec(particles),
    /** Model species **/
};

#undef deal_with_spec

consteval std::array<std::size_t, enumsize::SpeciesEnumSize + 1>
start_positions() {
  std::array<bool, enumsize::SpeciesEnumSize> found{};
  found.fill(false);

  std::array<std::size_t, enumsize::SpeciesEnumSize + 1> start{};
  start.fill(0);

  for (std::size_t i = 0; i < Isotopologues.size(); i++) {
    const auto spec = static_cast<std::size_t>(Isotopologues[i].spec);
    if (not found[spec]) {
      start[spec] = i;
      found[spec] = true;
    }
  }

  start.back() = Isotopologues.size();
  return start;
}

constexpr inline auto IsotopologuesStart = start_positions();

template <SpeciesEnum spec>
consteval std::size_t count_isotopologues()
  requires(IsotopologuesStart[static_cast<std::size_t>(spec)] <=
           IsotopologuesStart[static_cast<std::size_t>(spec) + 1])
{
  return IsotopologuesStart[static_cast<std::size_t>(spec) + 1] -
         IsotopologuesStart[static_cast<std::size_t>(spec)];
}

template <SpeciesEnum spec>
consteval std::array<Isotope, count_isotopologues<spec>()> isotopologues()
  requires(IsotopologuesStart[static_cast<std::size_t>(spec)] <=
           IsotopologuesStart[static_cast<std::size_t>(spec) + 1])
{
  static_assert(count_isotopologues<spec>() not_eq 0,
                "All species must be defined in the Isotopologues!");
  std::array<Isotope, count_isotopologues<spec>()> isots;
  for (std::size_t i = 0; i < count_isotopologues<spec>(); i++) {
    isots[i] = Isotopologues[i + IsotopologuesStart[std::size_t(spec)]];
  }
  return isots;
}

Array<Isotope> isotopologues(SpeciesEnum spec);

constexpr Index find_species_index(const SpeciesEnum spec,
                                   const std::string_view isot) {
  if (good_enum(spec)) {
    for (std::size_t i = IsotopologuesStart[std::size_t(spec)];
         i < IsotopologuesStart[std::size_t(spec) + 1];
         i++) {
      if (isot == Isotopologues[i].isotname) {
        return i;
      }
    }
  }

  throw std::runtime_error("Cannot find species index");
}

constexpr Index find_species_index(const Isotope& ir) {
  return find_species_index(ir.spec, ir.isotname);
}

constexpr Index find_species_index(const std::string_view spec,
                                   const std::string_view isot) {
  return find_species_index(to<SpeciesEnum>(spec), isot);
}

constexpr Index find_species_index(std::string_view s) {
  auto minus = s.find('-');
  return find_species_index(s.substr(0, minus),
                            minus == s.npos ? Joker : s.substr(minus + 1));
}

constexpr const Isotope& select(SpeciesEnum spec,
                                const std::string_view isotname) {
  return Isotopologues[find_species_index(spec, isotname)];
}

constexpr const Isotope& select(const std::string_view spec,
                                const std::string_view isotname) {
  return Isotopologues[find_species_index(spec, isotname)];
}

constexpr const Isotope& select(const std::string_view name) {
  return Isotopologues[find_species_index(name)];
}

String isotopologues_names(SpeciesEnum spec);

struct IsotopologueRatios {
  static constexpr Index maxsize = Index(Isotopologues.size());
  std::array<Numeric, maxsize> data;

  IsotopologueRatios();

  Numeric operator[](const Index spec_ind) const;

  Numeric operator[](const Isotope& ir) const;

  friend std::ostream& operator<<(std::ostream& os,
                                  const IsotopologueRatios& iso_rat);

  [[nodiscard]] bool all_isotopes_have_a_value() const;
};

IsotopologueRatios isotopologue_ratiosInitFromBuiltin();

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
constexpr Numeric mean_mass(SpeciesEnum spec, const IsotopologueRatios& ir) {
  Numeric sum_rm = 0;
  Numeric sum_r  = 0;
  for (std::size_t i = IsotopologuesStart[std::size_t(spec)];
       i < IsotopologuesStart[std::size_t(spec) + 1];
       i++) {
    if (not nonstd::isnan(Isotopologues[i].mass) and not nonstd::isnan(ir[i])) {
      sum_rm += ir[i] * Isotopologues[i].mass;
      sum_r  += ir[i];
    }
  }
  if (sum_r not_eq 0) return sum_rm / sum_r;
  return std::numeric_limits<Numeric>::signaling_NaN();
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

constexpr bool all_have_ratio(const SpeciesEnum spec,
                              const IsotopologueRatios& ir) {
  for (std::size_t i = IsotopologuesStart[std::size_t(spec)];
       i < IsotopologuesStart[std::size_t(spec) + 1];
       i++) {
    if (not Isotopologues[i].is_joker() and
        not Isotopologues[i].is_predefined() and nonstd::isnan(ir[i])) {
      return false;
    }
  }
  return true;
}

std::pair<ArrayOfString, ArrayOfString> names_of_have_and_havenot_ratio(
    const SpeciesEnum spec, const IsotopologueRatios& ir);

std::ostream& operator<<(std::ostream& os, const std::vector<Isotope>& isots);
}  // namespace Species

using SpeciesIsotope = Species::Isotope;

using ArrayOfSpeciesIsotope = Array<SpeciesIsotope>;

using SpeciesIsotopologueRatios = Species::IsotopologueRatios;

consteval SpeciesIsotope operator""_isot(const char* x, std::size_t) {
  return Species::select(x);
}

consteval Index operator""_isot_index(const char* x, std::size_t) {
  return Species::find_species_index(x);
}

template <>
struct std::hash<SpeciesIsotope> {
  std::size_t operator()(const SpeciesIsotope& g) const {
    return static_cast<std::size_t>(find_species_index(g));
  }
};

template <>
struct std::formatter<SpeciesIsotopologueRatios> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SpeciesIsotopologueRatios& v,
                              FmtContext& ctx) const {
    const auto sep = tags.sep();
    tags.format(ctx, Species::Isotopologues[0].FullName(), sep, v.data[0]);
    for (Index i = 1; i < v.maxsize; i++) {
      tags.format(
          ctx, '\n', Species::Isotopologues[i].FullName(), sep, v.data[i]);
    }

    return ctx.out();
  }
};

template <>
struct xml_io_stream<SpeciesIsotope> {
  static constexpr std::string_view type_name = "SpeciesIsotope"sv;

  static void write(std::ostream& os,
                    const SpeciesIsotope& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   SpeciesIsotope& x,
                   bifstream* pbifs = nullptr);
};
