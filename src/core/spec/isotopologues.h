#pragma once

#include <compare.h>
#include <enumsSpeciesEnum.h>
#include <mystring.h>
#include <nonstd.h>

#include <array>
#include <limits>
#include <string_view>
#include <tuple>

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

  [[nodiscard]] constexpr bool joker() const { return isotname == Joker; }
  [[nodiscard]] constexpr bool OK() const { return good_enum(spec); }
  [[nodiscard]] String FullName() const;

  friend std::ostream& operator<<(std::ostream& os, const Isotope& ir);

  constexpr auto operator==(const Isotope& that) const {
    return std::tie(spec, isotname) == std::tie(that.spec, that.isotname);
  }

  constexpr auto operator!=(const Isotope& that) const {
    return std::tie(spec, isotname) != std::tie(that.spec, that.isotname);
  }

  constexpr auto operator<=(const Isotope& that) const {
    return std::tie(spec, isotname) <= std::tie(that.spec, that.isotname);
  }

  constexpr auto operator>=(const Isotope& that) const {
    return std::tie(spec, isotname) >= std::tie(that.spec, that.isotname);
  }

  constexpr auto operator<(const Isotope& that) const {
    return std::tie(spec, isotname) < std::tie(that.spec, that.isotname);
  }

  constexpr auto operator>(const Isotope& that) const {
    return std::tie(spec, isotname) > std::tie(that.spec, that.isotname);
  }
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
    return std::format_to(ctx.out(), "{}{}{}", quote, v.FullName(), quote);
  }
};

namespace Species {

#define deal_with_spec(SPEC) Isotope(SpeciesEnum::SPEC)

/** A list of all ARTS isotopologues, note how the species enum class input HAS to be sorted */
inline constexpr std::array Isotopologues{
    deal_with_spec(Bath),
    /** Water species **/
    deal_with_spec(Water),
    Isotope(to<SpeciesEnum>("H2O"), "161", 18.010565, 1),
    Isotope(to<SpeciesEnum>("H2O"), "162", 19.016740, 6),
    Isotope(to<SpeciesEnum>("H2O"), "171", 19.014780, 6),
    Isotope(to<SpeciesEnum>("H2O"), "172", 20.020956, 36),
    Isotope(to<SpeciesEnum>("H2O"), "181", 20.014811, 1),
    Isotope(to<SpeciesEnum>("H2O"), "182", 21.020985, 6),
    Isotope(to<SpeciesEnum>("H2O"), "262", 20.022915, 1),
    Isotope(to<SpeciesEnum>("H2O"), "ForeignContCKDMT320"),
    Isotope(to<SpeciesEnum>("H2O"), "ForeignContCKDMT350"),
    Isotope(to<SpeciesEnum>("H2O"), "ForeignContCKDMT400"),
    Isotope(to<SpeciesEnum>("H2O"), "ForeignContStandardType"),
    Isotope(to<SpeciesEnum>("H2O"), "MPM89"),
    Isotope(to<SpeciesEnum>("H2O"), "PWR2021"),
    Isotope(to<SpeciesEnum>("H2O"), "PWR2022"),
    Isotope(to<SpeciesEnum>("H2O"), "PWR98"),
    Isotope(to<SpeciesEnum>("H2O"), "SelfContCKDMT320"),
    Isotope(to<SpeciesEnum>("H2O"), "SelfContCKDMT350"),
    Isotope(to<SpeciesEnum>("H2O"), "SelfContCKDMT400"),
    Isotope(to<SpeciesEnum>("H2O"), "SelfContStandardType"),
    /** Water species **/

    /** Carbon dioxide species **/
    deal_with_spec(CarbonDioxide),
    Isotope(to<SpeciesEnum>("CO2"), "626", 43.989830, 1),
    Isotope(to<SpeciesEnum>("CO2"), "627", 44.994045, 6),
    Isotope(to<SpeciesEnum>("CO2"), "628", 45.994076, 1),
    Isotope(to<SpeciesEnum>("CO2"), "636", 44.993185, 2),
    Isotope(to<SpeciesEnum>("CO2"), "637", 45.997400, 12),
    Isotope(to<SpeciesEnum>("CO2"), "638", 46.997431, 2),
    Isotope(to<SpeciesEnum>("CO2"), "727", 45.998262, 1),
    Isotope(to<SpeciesEnum>("CO2"), "737", 47.001618, 2),
    Isotope(to<SpeciesEnum>("CO2"), "827", 46.998291, 6),
    Isotope(to<SpeciesEnum>("CO2"), "828", 47.998322, 1),
    Isotope(to<SpeciesEnum>("CO2"), "837", 48.001646, 12),
    Isotope(to<SpeciesEnum>("CO2"), "838", 49.001675, 2),
    Isotope(to<SpeciesEnum>("CO2"), "CKDMT252"),
    /** Carbon dioxide species **/

    /** Ozone species **/
    deal_with_spec(Ozone),
    Isotope(to<SpeciesEnum>("O3"), "666", 47.984745, 1),
    Isotope(to<SpeciesEnum>("O3"), "667", 48.988960, 6),
    Isotope(to<SpeciesEnum>("O3"), "668", 49.988991, 1),
    Isotope(to<SpeciesEnum>("O3"), "676", 48.988960, 6),
    Isotope(to<SpeciesEnum>("O3"), "686", 49.988991, 1),
    /** Ozone species **/

    /** N2O species **/
    deal_with_spec(NitrogenOxide),
    Isotope(to<SpeciesEnum>("N2O"), "446", 44.001062, 9),
    Isotope(to<SpeciesEnum>("N2O"), "447", 45.005278, 54),
    Isotope(to<SpeciesEnum>("N2O"), "448", 46.005308, 9),
    Isotope(to<SpeciesEnum>("N2O"), "456", 44.998096, 6),
    Isotope(to<SpeciesEnum>("N2O"), "546", 44.998096, 6),
    /** N2O species **/

    /** CO species **/
    deal_with_spec(CarbonMonoxide),
    Isotope(to<SpeciesEnum>("CO"), "26", 27.994915, 1),
    Isotope(to<SpeciesEnum>("CO"), "27", 28.999130, 6),
    Isotope(to<SpeciesEnum>("CO"), "28", 29.999161, 1),
    Isotope(to<SpeciesEnum>("CO"), "36", 28.998270, 2),
    Isotope(to<SpeciesEnum>("CO"), "37", 30.002485, 12),
    Isotope(to<SpeciesEnum>("CO"), "38", 31.002516, 2),
    /** CO species **/

    /** CH4 species **/
    deal_with_spec(Methane),
    Isotope(to<SpeciesEnum>("CH4"), "211", 16.031300, 1),
    Isotope(to<SpeciesEnum>("CH4"), "212", 17.037475, 3),
    Isotope(to<SpeciesEnum>("CH4"), "311", 17.034655, 2),
    Isotope(to<SpeciesEnum>("CH4"), "312", 18.040830, 6),
    /** CH4 species **/

    /** Oxygen species **/
    deal_with_spec(Oxygen),
    Isotope(to<SpeciesEnum>("O2"), "66", 31.989830, 1),
    Isotope(to<SpeciesEnum>("O2"), "67", 32.994045, 6),
    Isotope(to<SpeciesEnum>("O2"), "68", 33.994076, 1),
    Isotope(to<SpeciesEnum>("O2"), "CIAfunCKDMT100"),
    Isotope(to<SpeciesEnum>("O2"), "MPM2020"),
    Isotope(to<SpeciesEnum>("O2"), "MPM89"),
    Isotope(to<SpeciesEnum>("O2"), "PWR2021"),
    Isotope(to<SpeciesEnum>("O2"), "PWR2022"),
    Isotope(to<SpeciesEnum>("O2"), "PWR98"),
    Isotope(to<SpeciesEnum>("O2"), "SelfContStandardType"),
    Isotope(to<SpeciesEnum>("O2"), "TRE05"),
    Isotope(to<SpeciesEnum>("O2"), "v0v0CKDMT100"),
    Isotope(to<SpeciesEnum>("O2"), "v1v0CKDMT100"),
    Isotope(to<SpeciesEnum>("O2"), "visCKDMT252"),
    /** Oxygen species **/

    /** NO species **/
    deal_with_spec(NitricOxide),
    Isotope(to<SpeciesEnum>("NO"), "46", 29.997989, 3),
    Isotope(to<SpeciesEnum>("NO"), "48", 32.002234, 3),
    Isotope(to<SpeciesEnum>("NO"), "56", 30.995023, 2),
    /** NO species **/

    /** SO2 species **/
    deal_with_spec(SulfurDioxide),
    Isotope(to<SpeciesEnum>("SO2"), "626", 63.961901, 1),
    Isotope(to<SpeciesEnum>("SO2"), "628", 65.966146, 1),
    Isotope(to<SpeciesEnum>("SO2"), "636", 64.961286, 4),
    Isotope(to<SpeciesEnum>("SO2"), "646", 65.957695, 1),
    /** SO2 species **/

    /** NO2 species **/
    deal_with_spec(NitrogenDioxide),
    Isotope(to<SpeciesEnum>("NO2"), "646", 45.992904, 3),
    Isotope(to<SpeciesEnum>("NO2"), "656", 46.989938, 2),
    /** NO2 species **/

    /** NH3 species **/
    deal_with_spec(Ammonia),
    Isotope(to<SpeciesEnum>("NH3"), "4111", 17.026549, 3),
    Isotope(
        to<SpeciesEnum>("NH3"), "4112", 18),  // FIXME: Better mass and some gj?
    Isotope(to<SpeciesEnum>("NH3"), "5111", 18.023583, 2),
    /** NH3 species **/

    /** HNO3 species **/
    deal_with_spec(NitricAcid),
    Isotope(to<SpeciesEnum>("HNO3"), "146", 62.995644, 6),
    Isotope(to<SpeciesEnum>("HNO3"), "156", 63.992680, 4),
    /** HNO3 species **/

    /** OH species **/
    deal_with_spec(Hydroxyl),
    Isotope(to<SpeciesEnum>("OH"), "61", 17.002740, 2),
    Isotope(to<SpeciesEnum>("OH"), "62", 18.008915, 3),
    Isotope(to<SpeciesEnum>("OH"), "81", 19.006986, 2),
    /** OH species **/

    /** HF species **/
    deal_with_spec(HydrogenFluoride),
    Isotope(to<SpeciesEnum>("HF"), "19", 20.006229, 4),
    Isotope(to<SpeciesEnum>("HF"), "29", 21.012404, 6),
    /** HF species **/

    /** HCl species **/
    deal_with_spec(HydrogenChloride),
    Isotope(to<SpeciesEnum>("HCl"), "15", 35.976678, 8),
    Isotope(to<SpeciesEnum>("HCl"), "17", 37.973729, 8),
    Isotope(to<SpeciesEnum>("HCl"), "25", 36.982853, 12),
    Isotope(to<SpeciesEnum>("HCl"), "27", 38.979904, 12),
    /** HCl species **/

    /** HBr species **/
    deal_with_spec(HydrogenBromide),
    Isotope(to<SpeciesEnum>("HBr"), "11", 81.924115, 8),
    Isotope(to<SpeciesEnum>("HBr"), "19", 79.926160, 8),
    Isotope(to<SpeciesEnum>("HBr"), "21", 82.930289, 12),
    Isotope(to<SpeciesEnum>("HBr"), "29", 80.932336, 12),
    /** HBr species **/

    /** HI species **/
    deal_with_spec(HydrogenIodide),
    Isotope(to<SpeciesEnum>("HI"), "17", 127.912297, 12),
    Isotope(to<SpeciesEnum>("HI"), "27", 128.918472, 18),
    /** HI species **/

    /** ClO species **/
    deal_with_spec(ChlorineMonoxide),
    Isotope(to<SpeciesEnum>("ClO"), "56", 50.963768, 4),
    Isotope(to<SpeciesEnum>("ClO"), "76", 52.960819, 4),
    /** ClO species **/

    /** OCS species **/
    deal_with_spec(CarbonylSulfide),
    Isotope(to<SpeciesEnum>("OCS"), "622", 59.966986, 1),
    Isotope(to<SpeciesEnum>("OCS"), "623", 60.966371, 4),
    Isotope(to<SpeciesEnum>("OCS"), "624", 61.962780, 1),
    Isotope(to<SpeciesEnum>("OCS"), "632", 60.970341, 2),
    Isotope(to<SpeciesEnum>("OCS"), "634", 62.966137, 2),
    Isotope(to<SpeciesEnum>("OCS"), "822", 61.971231, 1),
    /** OCS species **/

    /** H2CO species **/
    deal_with_spec(Formaldehyde),
    Isotope(to<SpeciesEnum>("H2CO"), "126", 30.010565, 1),
    Isotope(to<SpeciesEnum>("H2CO"), "128", 32.014811, 1),
    Isotope(to<SpeciesEnum>("H2CO"), "136", 31.013920, 2),
    /** H2CO species **/

    /** HDCO species nb. If the order D matters, rename this to indicate how **/
    deal_with_spec(HeavyFormaldehyde),
    Isotope(to<SpeciesEnum>("HDCO"),
            "26",
            31),  // FIXME: Better mass and some gj?  What is the AFGL code???
    /** HDCO species **/

    /** D2CO species **/
    deal_with_spec(VeryHeavyFormaldehyde),
    Isotope(to<SpeciesEnum>("D2CO"),
            "26",
            32),  // FIXME: Better mass and some gj?  What is the AFGL code???
    /** D2CO species **/

    /** HOCl species **/
    deal_with_spec(HypochlorousAcid),
    Isotope(to<SpeciesEnum>("HOCl"), "165", 51.971593, 8),
    Isotope(to<SpeciesEnum>("HOCl"), "167", 53.968644, 8),
    /** HOCl species **/

    /** N2 species **/
    deal_with_spec(Nitrogen),
    Isotope(to<SpeciesEnum>("N2"), "44", 28.006148, 1),
    Isotope(to<SpeciesEnum>("N2"), "45", 29.003182, 6),
    Isotope(to<SpeciesEnum>("N2"), "CIAfunCKDMT252"),
    Isotope(to<SpeciesEnum>("N2"), "CIArotCKDMT252"),
    Isotope(to<SpeciesEnum>("N2"), "SelfContMPM93"),
    Isotope(to<SpeciesEnum>("N2"), "SelfContPWR2021"),
    Isotope(to<SpeciesEnum>("N2"), "SelfContStandardType"),
    /** N2 species **/

    /** HCN species **/
    deal_with_spec(HydrogenCyanide),
    Isotope(to<SpeciesEnum>("HCN"), "124", 27.010899, 6),
    Isotope(to<SpeciesEnum>("HCN"), "125", 28.007933, 4),
    Isotope(to<SpeciesEnum>("HCN"), "134", 28.014254, 12),
    Isotope(
        to<SpeciesEnum>("HCN"), "224", 28),  // FIXME: Better mass and some gj?
    /** HCN species **/

    /** CH3Cl species **/
    deal_with_spec(Chloromethane),
    Isotope(to<SpeciesEnum>("CH3Cl"), "215", 49.992328, 4),
    Isotope(to<SpeciesEnum>("CH3Cl"), "217", 51.989379, 4),
    /** CH3Cl species **/

    /** H2O2 species **/
    deal_with_spec(HydrogenPeroxide),
    Isotope(to<SpeciesEnum>("H2O2"), "1661", 34.005480, 1),
    /** H2O2 species **/

    /** C2H2 species **/
    deal_with_spec(Acetylene),
    Isotope(to<SpeciesEnum>("C2H2"), "1221", 26.015650, 1),
    Isotope(to<SpeciesEnum>("C2H2"), "1222", 27.021825, 6),
    Isotope(to<SpeciesEnum>("C2H2"), "1231", 27.019005, 8),
    /** C2H2 species **/

    /** C2H6 species **/
    deal_with_spec(Ethane),
    Isotope(to<SpeciesEnum>("C2H6"), "1221", 30.046950, 1),
    Isotope(to<SpeciesEnum>("C2H6"), "1231", 31.050305, 2),
    /** C2H6 species **/

    /** PH3 species **/
    deal_with_spec(Phosphine),
    Isotope(to<SpeciesEnum>("PH3"), "1111", 33.997238, 2),
    /** PH3 species **/

    /** COF2 species **/
    deal_with_spec(CarbonylFluoride),
    Isotope(to<SpeciesEnum>("COF2"), "269", 65.991722, 1),
    Isotope(to<SpeciesEnum>("COF2"), "369", 66.995083, 2),
    /** COF2 species **/

    /** SF6 species **/
    deal_with_spec(SulfurHexafluoride),
    Isotope(to<SpeciesEnum>("SF6"), "29", 145.962492, 1),
    /** SF6 species **/

    /** H2S species **/
    deal_with_spec(HydrogenSulfide),
    Isotope(to<SpeciesEnum>("H2S"), "121", 33.987721, 1),
    Isotope(
        to<SpeciesEnum>("H2S"), "122", 35),  // FIXME: Better mass and some gj?
    Isotope(to<SpeciesEnum>("H2S"), "131", 34.987105, 4),
    Isotope(to<SpeciesEnum>("H2S"), "141", 35.983515, 1),
    /** H2S species **/

    /** HCOOH species **/
    deal_with_spec(FormicAcid),
    Isotope(to<SpeciesEnum>("HCOOH"), "126", 46.005480, 4),
    Isotope(to<SpeciesEnum>("HCOOH"),
            "136",
            47),  // FIXME: Better mass and some gj?
    /** HCOOH species **/

    /** DCOOH species **/
    deal_with_spec(LeftHeavyFormicAcid),
    Isotope(to<SpeciesEnum>("DCOOH"),
            "266",
            47),  // FIXME: Better mass and some gj?  What is the AFGL code???
    /** DCOOH species **/

    /** HCOOD species **/
    deal_with_spec(RightHeavyFormicAcid),
    Isotope(to<SpeciesEnum>("HCOOD"),
            "266",
            47),  // FIXME: Better mass and some gj?  What is the AFGL code???
    /** HCOOD species **/

    /** HO2 species **/
    deal_with_spec(Hydroperoxyl),
    Isotope(to<SpeciesEnum>("HO2"), "166", 32.997655, 2),
    /** HO2 species **/

    /** O species **/
    deal_with_spec(OxygenAtom),
    Isotope(to<SpeciesEnum>("O"), "6", 15.994915, 1),
    /** O species **/

    /** ClONO2 species **/
    deal_with_spec(ChlorineNitrate),
    Isotope(to<SpeciesEnum>("ClONO2"), "5646", 96.956672, 12),
    Isotope(to<SpeciesEnum>("ClONO2"), "7646", 98.953723, 12),
    /** ClONO2 species **/

    /** NO+ species **/
    deal_with_spec(NitricOxideCation),
    Isotope(to<SpeciesEnum>("NO+"), "46", 29.997989, 3),
    /** NO+ species **/

    /** HOBr species **/
    deal_with_spec(HypobromousAcid),
    Isotope(to<SpeciesEnum>("HOBr"), "161", 97.919027, 8),
    Isotope(to<SpeciesEnum>("HOBr"), "169", 95.921076, 8),
    /** HOBr species **/

    /** C2H4 species **/
    deal_with_spec(Ethylene),
    Isotope(to<SpeciesEnum>("C2H4"), "221", 28.031300, 1),
    Isotope(to<SpeciesEnum>("C2H4"), "231", 29.034655, 2),
    /** C2H4 species **/

    /** CH3OH species **/
    deal_with_spec(Methanol),
    Isotope(to<SpeciesEnum>("CH3OH"), "2161", 32.026215, 2),
    /** CH3OH species **/

    /** CH3Br species **/
    deal_with_spec(Bromomethane),
    Isotope(to<SpeciesEnum>("CH3Br"), "211", 95.939764, 4),
    Isotope(to<SpeciesEnum>("CH3Br"), "219", 93.941811, 4),
    /** CH3Br species **/

    /** CH3CN species **/
    deal_with_spec(Acetonitrile),
    Isotope(to<SpeciesEnum>("CH3CN"), "2124", 41.026549, 3),
    Isotope(to<SpeciesEnum>("CH3CN"),
            "2125",
            42),  // FIXME: Better mass and some gj?
    Isotope(to<SpeciesEnum>("CH3CN"),
            "2134",
            42),  // FIXME: Better mass and some gj?
    Isotope(to<SpeciesEnum>("CH3CN"),
            "3124",
            42),  // FIXME: Better mass and some gj?
    /** CH3CN species **/

    /** CH2DCN species nb. If the order D matters, rename this to indicate how **/
    deal_with_spec(HeavyAcetonitrile),
    Isotope(to<SpeciesEnum>("CH2DCN"),
            "224",
            42),  // FIXME: Better mass and some gj?  What is the AFGL code???
    /** CH2DCN species **/

    /** CF4 species **/
    deal_with_spec(CarbonTetrafluoride),
    Isotope(to<SpeciesEnum>("CF4"), "29", 87.993616, 1),
    /** CF4 species **/

    /** C4H2 species **/
    deal_with_spec(Diacetylene),
    Isotope(to<SpeciesEnum>("C4H2"), "2211", 50.015650, 1),
    /** C4H2 species **/

    /** HC3N species **/
    deal_with_spec(Cyanoacetylene),
    Isotope(to<SpeciesEnum>("HC3N"), "12224", 51.010899, 6),
    Isotope(to<SpeciesEnum>("HC3N"),
            "12225",
            52),  // FIXME: Better mass and some gj?
    Isotope(to<SpeciesEnum>("HC3N"),
            "12234",
            52),  // FIXME: Better mass and some gj?
    Isotope(to<SpeciesEnum>("HC3N"),
            "12324",
            52),  // FIXME: Better mass and some gj?
    Isotope(to<SpeciesEnum>("HC3N"),
            "13224",
            52),  // FIXME: Better mass and some gj?
    Isotope(to<SpeciesEnum>("HC3N"),
            "22224",
            52),  // FIXME: Better mass and some gj?
    /** HC3N species **/

    /** H2 species **/
    deal_with_spec(Hydrogen),
    Isotope(to<SpeciesEnum>("H2"), "11", 2.015650, 1),
    Isotope(to<SpeciesEnum>("H2"), "12", 3.021825, 6),
    /** H2 species **/

    /** CS species **/
    deal_with_spec(CarbonMonosulfide),
    Isotope(to<SpeciesEnum>("CS"), "22", 43.971036, 1),
    Isotope(to<SpeciesEnum>("CS"), "23", 44.970399, 4),
    Isotope(to<SpeciesEnum>("CS"), "24", 45.966787, 1),
    Isotope(to<SpeciesEnum>("CS"), "32", 44.974368, 2),
    /** CS species **/

    /** SO3 species **/
    deal_with_spec(SulfurTrioxide),
    Isotope(to<SpeciesEnum>("SO3"), "26", 79.956820, 1),
    /** SO3 species **/

    /** C2N2 species **/
    deal_with_spec(Cyanogen),
    Isotope(to<SpeciesEnum>("C2N2"), "4224", 52.006148, 1),
    /** C2N2 species **/

    /** COCl2 species **/
    deal_with_spec(Phosgene),
    Isotope(to<SpeciesEnum>("COCl2"), "2655", 97.932620, 1),
    Isotope(to<SpeciesEnum>("COCl2"), "2657", 99.929670, 16),
    /** COCl2 species **/

    /** SO species **/
    deal_with_spec(SulfurMonoxide),
    Isotope(to<SpeciesEnum>("SO"), "26", 47.966986, 1),
    Isotope(to<SpeciesEnum>("SO"), "28", 49.971231, 1),
    Isotope(to<SpeciesEnum>("SO"), "46", 49.962782, 1),
    /** SO species **/

    /** CS2 species **/
    deal_with_spec(CarbonDisulfide),
    Isotope(to<SpeciesEnum>("CS2"), "222", 75.944140, 1),
    Isotope(to<SpeciesEnum>("CS2"), "223", 76.943256, 4),
    Isotope(to<SpeciesEnum>("CS2"), "224", 77.939940, 1),
    Isotope(to<SpeciesEnum>("CS2"), "232", 76.947495, 2),
    /** CS2 species **/

    deal_with_spec(Methyl),
    deal_with_spec(Cyclopropene),

    /** H2SO4 species **/
    deal_with_spec(SulfuricAcid),
    Isotope(to<SpeciesEnum>("H2SO4"),
            "126",
            98),  // FIXME: Better mass and some gj?
    /** H2SO4 species **/

    /** HNC species **/
    deal_with_spec(HydrogenIsocyanide),
    Isotope(
        to<SpeciesEnum>("HNC"), "142", 27),  // FIXME: Better mass and some gj?
    Isotope(
        to<SpeciesEnum>("HNC"), "143", 28),  // FIXME: Better mass and some gj?
    Isotope(
        to<SpeciesEnum>("HNC"), "152", 28),  // FIXME: Better mass and some gj?
    Isotope(
        to<SpeciesEnum>("HNC"), "242", 28),  // FIXME: Better mass and some gj?
    /** HNC species **/

    /** BrO species **/
    deal_with_spec(BromineMonoxide),
    Isotope(
        to<SpeciesEnum>("BrO"), "16", 97),  // FIXME: Better mass and some gj?
    Isotope(
        to<SpeciesEnum>("BrO"), "96", 95),  // FIXME: Better mass and some gj?
    /** BrO species **/

    /** OClO species **/
    deal_with_spec(ChlorineDioxide),
    Isotope(
        to<SpeciesEnum>("OClO"), "656", 67),  // FIXME: Better mass and some gj?
    Isotope(
        to<SpeciesEnum>("OClO"), "676", 69),  // FIXME: Better mass and some gj?
    /** OClO species **/

    /** C3H8 species **/
    deal_with_spec(Propane),
    Isotope(
        to<SpeciesEnum>("C3H8"), "21", 54),  // FIXME: Better mass and some gj?
    /** C3H8 species **/

    /** He species **/
    deal_with_spec(Helium),
    Isotope(to<SpeciesEnum>("He"), "4", 4),  // FIXME: Better mass and some gj?
    /** He species **/

    /** Cl2O2 species **/
    deal_with_spec(ChlorineMonoxideDimer),
    Isotope(to<SpeciesEnum>("Cl2O2"),
            "565",
            102),  // FIXME: Better mass and some gj?
    Isotope(to<SpeciesEnum>("Cl2O2"),
            "765",
            104),  // FIXME: Better mass and some gj?
    /** Cl2O2 species **/

    /** H species **/
    deal_with_spec(HydrogenAtom),
    Isotope(to<SpeciesEnum>("H"), "1", 1),  // FIXME: Better mass and some gj?
    /** H species **/

    /** Ar species **/
    deal_with_spec(Argon),
    Isotope(
        to<SpeciesEnum>("Ar"), "8", 39.948),  // FIXME: Better mass and some gj?
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
    Isotope(to<SpeciesEnum>("NF3"), "4999", 70.998286, 3),
    /** NF3 species **/

    deal_with_spec(SulfurylFluoride),
    deal_with_spec(HFC4310mee),

    /** GeH4 species **/
    deal_with_spec(Germane),
    Isotope(to<SpeciesEnum>("GeH4"), "011", 73.955550, 1),
    Isotope(to<SpeciesEnum>("GeH4"), "211", 75.953380, 1),
    Isotope(to<SpeciesEnum>("GeH4"), "311", 76.954764, 10),
    Isotope(to<SpeciesEnum>("GeH4"), "411", 77.952479, 1),
    Isotope(to<SpeciesEnum>("GeH4"), "611", 79.952703, 1),
    /** GeH4 species **/

    /** CH3I species **/
    deal_with_spec(Iodomethane),
    Isotope(to<SpeciesEnum>("CH3I"), "217", 141.927947, 6),
    /** CH3I species **/

    /** CH3F species **/
    deal_with_spec(Fluoromethane),
    Isotope(to<SpeciesEnum>("CH3F"), "219", 34.021878, 2),
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

constexpr Index find_species_index(const Isotope ir) {
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

constexpr const Isotope& select_joker(SpeciesEnum spec) {
  return select(spec, Joker);
}

constexpr const Isotope& select_joker(std::string_view spec) {
  return select(to<SpeciesEnum>(spec), Joker);
}

String isotopologues_names(SpeciesEnum spec);

constexpr bool is_predefined_model(const Isotope& ir) {
  return not(nonstd::isdigit(ir.isotname[0]) or ir.isotname == Joker);
}

constexpr bool is_normal_isotopologue(const Isotope& ir) {
  return nonstd::isdigit(ir.isotname[0]) and ir.isotname not_eq Joker;
}

String predefined_model_names();

constexpr bool same_or_joker(const Isotope& ir1, const Isotope& ir2) {
  if (ir1.spec not_eq ir2.spec) return false;
  if (ir1.joker() or ir2.joker()) return true;
  return ir1.isotname == ir2.isotname;
}

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

IsotopologueRatios isotopologue_ratiosInitFromBuiltin() ;

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
    if (not Isotopologues[i].joker() and
        not is_predefined_model(Isotopologues[i]) and nonstd::isnan(ir[i])) {
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
    for (Index i = 0; i < v.maxsize; i++) {
      tags.format(ctx,
                  Species::Isotopologues[i].FullName(),
                  tags.sep(),
                  v.data[i],
                  tags.sep(true));
    }

    return ctx.out();
  }
};
