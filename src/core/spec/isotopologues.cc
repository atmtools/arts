#include "isotopologues.h"

#include <debug.h>
#include <xml_io_base.h>

#include <cassert>
#include <compare>
#include <sstream>

namespace Species {
ArrayOfSpeciesIsotope isotopologues(SpeciesEnum spec) {
#define deal_with_spec(SPEC)                                      \
  case SpeciesEnum::SPEC: {                                       \
    static constexpr auto v = isotopologues<SpeciesEnum::SPEC>(); \
    return {v.begin(), v.end()};                                  \
  } break

  switch (spec) {
    case SpeciesEnum::Bath:
      break;
      deal_with_spec(Water);
      deal_with_spec(CarbonDioxide);
      deal_with_spec(Ozone);
      deal_with_spec(NitrogenOxide);
      deal_with_spec(CarbonMonoxide);
      deal_with_spec(Methane);
      deal_with_spec(Oxygen);
      deal_with_spec(NitricOxide);
      deal_with_spec(SulfurDioxide);
      deal_with_spec(NitrogenDioxide);
      deal_with_spec(Ammonia);
      deal_with_spec(NitricAcid);
      deal_with_spec(Hydroxyl);
      deal_with_spec(HydrogenFluoride);
      deal_with_spec(HydrogenChloride);
      deal_with_spec(HydrogenBromide);
      deal_with_spec(HydrogenIodide);
      deal_with_spec(ChlorineMonoxide);
      deal_with_spec(CarbonylSulfide);
      deal_with_spec(Formaldehyde);
      deal_with_spec(HeavyFormaldehyde);
      deal_with_spec(VeryHeavyFormaldehyde);
      deal_with_spec(HypochlorousAcid);
      deal_with_spec(Nitrogen);
      deal_with_spec(HydrogenCyanide);
      deal_with_spec(Chloromethane);
      deal_with_spec(HydrogenPeroxide);
      deal_with_spec(Acetylene);
      deal_with_spec(Ethane);
      deal_with_spec(Phosphine);
      deal_with_spec(CarbonylFluoride);
      deal_with_spec(SulfurHexafluoride);
      deal_with_spec(HydrogenSulfide);
      deal_with_spec(FormicAcid);
      deal_with_spec(LeftHeavyFormicAcid);
      deal_with_spec(RightHeavyFormicAcid);
      deal_with_spec(Hydroperoxyl);
      deal_with_spec(OxygenAtom);
      deal_with_spec(ChlorineNitrate);
      deal_with_spec(NitricOxideCation);
      deal_with_spec(HypobromousAcid);
      deal_with_spec(Ethylene);
      deal_with_spec(Methanol);
      deal_with_spec(Bromomethane);
      deal_with_spec(Acetonitrile);
      deal_with_spec(HeavyAcetonitrile);
      deal_with_spec(CarbonTetrafluoride);
      deal_with_spec(Diacetylene);
      deal_with_spec(Cyanoacetylene);
      deal_with_spec(Hydrogen);
      deal_with_spec(CarbonMonosulfide);
      deal_with_spec(SulfurTrioxide);
      deal_with_spec(Cyanogen);
      deal_with_spec(Phosgene);
      deal_with_spec(SulfurMonoxide);
      deal_with_spec(CarbonDisulfide);
      deal_with_spec(Methyl);
      deal_with_spec(Cyclopropene);
      deal_with_spec(SulfuricAcid);
      deal_with_spec(HydrogenIsocyanide);
      deal_with_spec(BromineMonoxide);
      deal_with_spec(ChlorineDioxide);
      deal_with_spec(Propane);
      deal_with_spec(Helium);
      deal_with_spec(ChlorineMonoxideDimer);
      deal_with_spec(HydrogenAtom);
      deal_with_spec(Argon);
      deal_with_spec(Hexafluoroethane);
      deal_with_spec(Perfluoropropane);
      deal_with_spec(Perfluorobutane);
      deal_with_spec(Perfluoropentane);
      deal_with_spec(Perfluorohexane);
      deal_with_spec(Perfluorooctane);
      deal_with_spec(Perfluorocyclobutane);
      deal_with_spec(CarbonTetrachloride);
      deal_with_spec(CFC11);
      deal_with_spec(CFC113);
      deal_with_spec(CFC114);
      deal_with_spec(CFC115);
      deal_with_spec(CFC12);
      deal_with_spec(Dichloromethane);
      deal_with_spec(Trichloroethane);
      deal_with_spec(Trichloromethane);
      deal_with_spec(Bromochlorodifluoromethane);
      deal_with_spec(Bromotrifluoromethane);
      deal_with_spec(Dibromotetrafluoroethane);
      deal_with_spec(HCFC141b);
      deal_with_spec(HCFC142b);
      deal_with_spec(HCFC22);
      deal_with_spec(HFC125);
      deal_with_spec(HFC134a);
      deal_with_spec(HFC143a);
      deal_with_spec(HFC152a);
      deal_with_spec(HFC227ea);
      deal_with_spec(HFC23);
      deal_with_spec(HFC236fa);
      deal_with_spec(HFC245fa);
      deal_with_spec(HFC32);
      deal_with_spec(HFC365mfc);
      deal_with_spec(NitrogenTrifluoride);
      deal_with_spec(SulfurylFluoride);
      deal_with_spec(HFC4310mee);
      deal_with_spec(Germane);
      deal_with_spec(Iodomethane);
      deal_with_spec(Fluoromethane);
      deal_with_spec(Arsine);
      deal_with_spec(Benzene);
      deal_with_spec(liquidcloud);
      deal_with_spec(icecloud);
      deal_with_spec(rain);
      deal_with_spec(free_electrons);
      deal_with_spec(particles);
      deal_with_spec(unused);
  }

#undef deal_with_spec

  ARTS_USER_ERROR("Cannot understand: {}", spec)
}

String isotopologues_names(SpeciesEnum spec) {
  auto x = isotopologues(spec);
  std::ostringstream os;
  for (auto& s : x) os << s.FullName() << '\n';
  return os.str();
}

String predefined_model_names() {
  std::ostringstream os;
  for (auto& x : Isotopologues) {
    if (x.is_predefined()) {
      os << x.FullName() << '\n';
    }
  }
  return os.str();
}

String update_isot_name(const String& old_name) {
  if (old_name == "CH3CN-211224") return "CH2DCN-224";
  if (old_name == "CH3CN-211124") return "CH3CN-2124";
  if (old_name == "CH3CN-211125") return "CH3CN-2125";
  if (old_name == "CH3CN-211134") return "CH3CN-2134";
  if (old_name == "CH3CN-311124") return "CH3CN-3124";
  if (old_name == "CO2-728") return "CO2-827";
  if (old_name == "HCOOH-2261") return "DCOOH-266";
  if (old_name == "HCOOH-1262") return "HCOOD-266";
  if (old_name == "HCOOH-1261") return "HCOOH-126";
  if (old_name == "HCOOH-1361") return "HCOOH-136";
  if (old_name == "H2CO-1126") return "H2CO-126";
  if (old_name == "H2CO-1128") return "H2CO-128";
  if (old_name == "H2CO-1136") return "H2CO-136";
  if (old_name == "H2CO-1226") return "HDCO-26";
  if (old_name == "H2CO-2226") return "D2CO-26";
  return old_name;
}

Isotope::Isotope(const std::string_view name) { *this = select(name); }

String Isotope::FullName() const {
  return is_joker() ? String{toString<1>(spec)}
                    : std::format("{}-{}", toString<1>(spec), isotname);
}

std::ostream& operator<<(std::ostream& os, const Isotope& ir) {
  return os << ir.FullName();
}

std::ostream& operator<<(std::ostream& os, const std::vector<Isotope>& isots) {
  for (const auto& i : isots) os << i << ' ';
  return os;
}

IsotopologueRatios::IsotopologueRatios() : data() {
  for (auto& x : data) x = std::numeric_limits<Numeric>::quiet_NaN();
}

Numeric IsotopologueRatios::operator[](const Index spec_ind) const {
  assert(spec_ind < maxsize and spec_ind >= 0);
  return data[spec_ind];
}

Numeric IsotopologueRatios::operator[](const Isotope& ir) const {
  const Index spec_ind = find_species_index(ir);
  ARTS_USER_ERROR_IF(
      spec_ind >= maxsize or spec_ind < 0, "Invalid species {}", ir.FullName())
  return data[spec_ind];
}

bool IsotopologueRatios::all_isotopes_have_a_value() const {
  for (Index i = 0; i < maxsize; i++) {
    if (not Isotopologues[i].is_predefined() and
        not Isotopologues[i].is_joker() and nonstd::isnan(data[i])) {
      return false;
    }
  }
  return true;
}

IsotopologueRatios isotopologue_ratiosInitFromBuiltin() {
  IsotopologueRatios isotopologue_ratios;

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("H2O", ISOT)] = VAL
  set_isot_val("161", .997317E+00);
  set_isot_val("181", 1.99983E-03);
  set_isot_val("171", 3.71884E-04);
  set_isot_val("162", 3.10693E-04);
  set_isot_val("182", 6.23003E-07);
  set_isot_val("172", 1.15853E-07);
  set_isot_val("262", 2.41970E-08);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CO2", ISOT)] = VAL
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
  set_isot_val("737", 1.53750E-09);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("O3", ISOT)] = VAL
  set_isot_val("666", .992901E+00);
  set_isot_val("668", 3.98194E-03);
  set_isot_val("686", 1.99097E-03);
  set_isot_val("667", 7.40475E-04);
  set_isot_val("676", 3.70237E-04);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("N2O", ISOT)] = VAL
  set_isot_val("446", .990333E+00);
  set_isot_val("456", 3.64093E-03);
  set_isot_val("546", 3.64093E-03);
  set_isot_val("448", 1.98582E-03);
  set_isot_val("447", 3.69280E-04);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CO", ISOT)] = VAL
  set_isot_val("26", .986544E+00);
  set_isot_val("36", 1.10836E-02);
  set_isot_val("28", 1.97822E-03);
  set_isot_val("27", 3.67867E-04);
  set_isot_val("38", 2.22250E-05);
  set_isot_val("37", 4.13292E-06);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CH4", ISOT)] = VAL
  set_isot_val("211", .988274E+00);
  set_isot_val("311", 1.11031E-02);
  set_isot_val("212", 6.15751E-04);
  set_isot_val("312", 6.91785E-06);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("O2", ISOT)] = VAL
  set_isot_val("66", .995262E+00);
  set_isot_val("68", 3.99141E-03);
  set_isot_val("67", 7.42235E-04);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("NO", ISOT)] = VAL
  set_isot_val("46", .993974E+00);
  set_isot_val("56", 3.65431E-03);
  set_isot_val("48", 1.99312E-03);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("SO2", ISOT)] = VAL
  set_isot_val("626", .945678E+00);
  set_isot_val("646", 4.19503E-02);
  set_isot_val("636", 0.0074989421);
  set_isot_val("628", 0.0020417379);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("NO2", ISOT)] = VAL
  set_isot_val("646", .991616E+00);
  set_isot_val("656", 3.64564E-03);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("NH3", ISOT)] = VAL
  set_isot_val("4111", .995872E+00);
  set_isot_val("5111", 3.66129E-03);
  set_isot_val("4112", 0.00044792294);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HNO3", ISOT)] = VAL
  set_isot_val("146", .989110E+00);
  set_isot_val("156", 3.63600E-03);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("OH", ISOT)] = VAL
  set_isot_val("61", .997473E+00);
  set_isot_val("81", 2.00014E-03);
  set_isot_val("62", 1.55371E-04);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HF", ISOT)] = VAL
  set_isot_val("19", .999844E+00);
  set_isot_val("29", 1.55741E-04);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HCl", ISOT)] = VAL
  set_isot_val("15", .757587E+00);
  set_isot_val("17", .242257E+00);
  set_isot_val("25", 1.18005E-04);
  set_isot_val("27", 3.77350E-05);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HBr", ISOT)] = VAL
  set_isot_val("19", .506781E+00);
  set_isot_val("11", .493063E+00);
  set_isot_val("29", 7.89384E-05);
  set_isot_val("21", 7.68016E-05);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HI", ISOT)] = VAL
  set_isot_val("17", .999844E+00);
  set_isot_val("27", 1.55741E-04);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("ClO", ISOT)] = VAL
  set_isot_val("56", .755908E+00);
  set_isot_val("76", .241720E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("OCS", ISOT)] = VAL
  set_isot_val("622", .937395E+00);
  set_isot_val("624", 4.15828E-02);
  set_isot_val("632", 1.05315E-02);
  set_isot_val("623", 7.39908E-03);
  set_isot_val("822", 1.87967E-03);
  set_isot_val("634", 4.67508E-04);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("H2CO", ISOT)] = VAL
  set_isot_val("126", .986237E+00);
  set_isot_val("136", 1.10802E-02);
  set_isot_val("128", 1.97761E-03);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HDCO", ISOT)] = VAL
  set_isot_val("26", 0.00029578940);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("D2CO", ISOT)] = VAL
  set_isot_val("26", 2.2181076E-08);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HOCl", ISOT)] = VAL
  set_isot_val("165", .755790E+00);
  set_isot_val("167", .241683E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("N2", ISOT)] = VAL
  set_isot_val("44", .992687E+00);
  set_isot_val("45", 7.47809E-03);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HCN", ISOT)] = VAL
  set_isot_val("124", .985114E+00);
  set_isot_val("134", 1.10676E-02);
  set_isot_val("125", 3.62174E-03);
  set_isot_val("224", 0.00014773545);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CH3Cl", ISOT)] = VAL
  set_isot_val("215", .748937E+00);
  set_isot_val("217", .239491E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("H2O2", ISOT)] = VAL
  set_isot_val("1661", .994952E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("C2H2", ISOT)] = VAL
  set_isot_val("1221", .977599E+00);
  set_isot_val("1231", 2.19663E-02);
  set_isot_val("1222", 3.04550E-04);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("C2H6", ISOT)] = VAL
  set_isot_val("1221", .976990E+00);
  set_isot_val("1231", 2.19526E-02);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("PH3", ISOT)] = VAL
  set_isot_val("1111", .999533E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("COF2", ISOT)] = VAL
  set_isot_val("269", .986544E+00);
  set_isot_val("369", 1.10834E-02);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("SF6", ISOT)] = VAL
  set_isot_val("29", .950180E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("H2S", ISOT)] = VAL
  set_isot_val("121", .949884E+00);
  set_isot_val("141", 4.21369E-02);
  set_isot_val("131", 7.49766E-03);
  set_isot_val("122", 0.00029991625);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HCOOH", ISOT)] = VAL
  set_isot_val("126", .983898E+00);
  set_isot_val("136", 0.010913149);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("DCOOH", ISOT)] = VAL
  set_isot_val("266", 0.00014755369);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HCOOD", ISOT)] = VAL
  set_isot_val("266", 0.00014755369);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HO2", ISOT)] = VAL
  set_isot_val("166", .995107E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("O", ISOT)] = VAL
  set_isot_val("6", .997628E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("ClONO2", ISOT)] = VAL
  set_isot_val("5646", .749570E+00);
  set_isot_val("7646", .239694E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("NO+", ISOT)] = VAL
  set_isot_val("46", .993974E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("OClO", ISOT)] = VAL
  set_isot_val("656", 0.75509223);
  set_isot_val("676", 0.24490632);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("BrO", ISOT)] = VAL
  set_isot_val("96", 0.50582466);
  set_isot_val("16", 0.49431069);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("H2SO4", ISOT)] = VAL
  set_isot_val("126", 0.95060479);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("Cl2O2", ISOT)] = VAL
  set_isot_val("565", 0.57016427);
  set_isot_val("765", 0.36982818);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HOBr", ISOT)] = VAL
  set_isot_val("169", .505579E+00);
  set_isot_val("161", .491894E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("C2H4", ISOT)] = VAL
  set_isot_val("221", .977294E+00);
  set_isot_val("231", .219595E-01);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CH3OH", ISOT)] = VAL
  set_isot_val("2161", .985930E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CH3Br", ISOT)] = VAL
  set_isot_val("219", .500995E+00);
  set_isot_val("211", .487433E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CH3CN", ISOT)] = VAL
  set_isot_val("2124", .973866E+00);
  set_isot_val("3124", .102683e-01);
  set_isot_val("2134", .102683e-01);
  set_isot_val("2125", .347136e-02);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CH2DCN", ISOT)] = VAL
  set_isot_val("224", .441185e-03);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CF4", ISOT)] = VAL
  set_isot_val("29", .988890E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HC3N", ISOT)] = VAL
  set_isot_val("12224", .963346E+00);
  set_isot_val("12234", .106852e-01);
  set_isot_val("12324", .106852e-01);
  set_isot_val("13224", .106852e-01);
  set_isot_val("12225", .356272e-02);
  set_isot_val("22224", .144472e-03);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CS", ISOT)] = VAL
  set_isot_val("22", .939624E+00);
  set_isot_val("24", .416817E-01);
  set_isot_val("32", .105565E-01);
  set_isot_val("23", .741668E-02);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("HNC", ISOT)] = VAL
  set_isot_val("142", .985280e+00);
  set_isot_val("143", .109285e-01);
  set_isot_val("152", .364384e-02);
  set_isot_val("242", .147761e-03);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("SO", ISOT)] = VAL
  set_isot_val("26", .950605e+00);
  set_isot_val("46", .420727e-01);
  set_isot_val("28", .194089e-02);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("C3H8", ISOT)] = VAL
  set_isot_val("21", 9.66290e-01);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("H2", ISOT)] = VAL
  set_isot_val("11", .999688E+00);
  set_isot_val("12", 3.11432E-04);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("H", ISOT)] = VAL
  set_isot_val("1", 1.00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("He", ISOT)] = VAL
  set_isot_val("4", 1.00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("Ar", ISOT)] = VAL
  set_isot_val("8", 1.00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("C4H2", ISOT)] = VAL
  set_isot_val("2211", .955998E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("SO3", ISOT)] = VAL
  set_isot_val("26", .943400E+00);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CS2", ISOT)] = VAL
  set_isot_val("222", 8.92811E-01);
  set_isot_val("224", 7.92600E-02);
  set_isot_val("223", 1.40940E-02);
  set_isot_val("232", 1.03100E-02);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("C2N2", ISOT)] = VAL
  set_isot_val("4224", 9.70752E-01);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("COCl2", ISOT)] = VAL
  set_isot_val("2655", 5.66392E-01);
  set_isot_val("2657", 3.62235E-01);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CH3F", ISOT)] = VAL
  set_isot_val("219", 9.88428E-01);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("GeH4", ISOT)] = VAL
  set_isot_val("411", 3.65172E-01);
  set_isot_val("211", 2.74129E-01);
  set_isot_val("011", 2.05072E-01);
  set_isot_val("311", 7.75517E-01);
  set_isot_val("611", 7.75517E-01);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("CH3I", ISOT)] = VAL
  set_isot_val("217", 9.88428E-01);
#undef set_isot_val

#define set_isot_val(ISOT, VAL) \
  isotopologue_ratios.data[find_species_index("NF3", ISOT)] = VAL
  set_isot_val("4999", 9.96337E-01);
#undef set_isot_val

  assert(isotopologue_ratios.all_isotopes_have_a_value());
  return isotopologue_ratios;
}

std::strong_ordering Isotope::operator<=>(const Isotope& other) const {
  if (std::strong_ordering test = spec <=> other.spec;
      test != std::strong_ordering::equal)
    return test;
  return isotname <=> other.isotname;
}
static_assert(std::strong_ordering::equal == std::strong_ordering::equivalent);

bool Isotope::operator==(const Isotope& other) const {
  return this->operator<=>(other) == std::strong_ordering::equal;
}

bool Isotope::operator!=(const Isotope& other) const {
  return not this->operator==(other);
}
}  // namespace Species


void xml_io_stream<SpeciesIsotope>::write(std::ostream& os,
                                          const SpeciesIsotope& x,
                                          bofstream*,
                                          std::string_view name) {
  XMLTag tag(type_name, "name", name, "isot", x.FullName());
  tag.write_to_stream(os);
  tag.write_to_end_stream(os);
}

void xml_io_stream<SpeciesIsotope>::read(std::istream& is,
                                         SpeciesIsotope& x,
                                         bifstream*) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  String v;
  tag.get_attribute_value("isot", v);
  x = SpeciesIsotope(v);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
