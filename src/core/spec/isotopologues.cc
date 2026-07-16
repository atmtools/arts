#include "isotopologues.h"

#include <debug.h>
#include <xml_io_base.h>

#include <cassert>
#include <compare>
#include <sstream>

namespace Species {
String isotopologues_names(SpeciesEnum spec) {
  auto               x = isotopologues(spec);
  std::ostringstream os;
  for (auto& s : x) os << s.FullName() << '\n';
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
  if (old_name == "HC3N-1224") return "HC3N-12224";
  return old_name;
}

const Isotope& Isotope::from_name(const std::string_view name) { return select(name); }

String Isotope::FullName() const {
  return is_joker() ? String{toString<1>(spec)} : std::format("{}-{}", toString<1>(spec), isotname);
}

std::ostream& operator<<(std::ostream& os, const Isotope& ir) { return os << ir.FullName(); }

std::ostream& operator<<(std::ostream& os, const std::vector<Isotope>& isots) {
  for (const auto& i : isots) os << i << ' ';
  return os;
}

Numeric IsotopologueRatios::operator[](const Index spec_ind) const {
  assert(spec_ind < maxsize and spec_ind >= 0);
  return data[spec_ind];
}

Numeric IsotopologueRatios::operator[](const Isotope& ir) const {
  const Index spec_ind = find_species_index(ir);
  ARTS_USER_ERROR_IF(spec_ind >= maxsize or spec_ind < 0, "Invalid species {}", ir.FullName())
  return data[spec_ind];
}

std::vector<std::string> IsotopologueRatios::valueless_isotopes() const {
  std::vector<std::string> names;

  for (Index i = 0; i < maxsize; i++) {
    if (not Isotopologues[i].is_predefined() and not Isotopologues[i].is_joker() and nonstd::isnan(data[i])) {
      names.push_back(Isotopologues[i].FullName());
    }
  }

  return names;
}

namespace {
consteval IsotopologueRatios from_builtin() {
  IsotopologueRatios isotopologue_ratios{};

  stdr::copy(Isotopologues | std::views::transform(&Isotope::builtin_ratio), isotopologue_ratios.data.begin());

  return isotopologue_ratios;
}

consteval bool all_values() {
  for (const auto& iso : Isotopologues) {
    if (not(iso.is_predefined() or iso.is_joker())) {
      if (nonstd::isnan(iso.builtin_ratio) or nonstd::isnan(iso.mass) or iso.gi == 0) { return false; }
    }
  }
  return true;
}

static_assert(all_values(),
              "Some isotopologues do not have all values defined.  "
              "Please check the source files.  "
              "You must define the missing values.");
}  // namespace

const IsotopologueRatios& isotopologue_ratiosInitFromBuiltin() {
  static constexpr IsotopologueRatios isotopologue_ratios = from_builtin();
  return isotopologue_ratios;
}

std::strong_ordering Isotope::operator<=>(const Isotope& other) const {
  if (std::strong_ordering test = spec <=> other.spec; test != std::strong_ordering::equal) return test;
  return isotname <=> other.isotname;
}
static_assert(std::strong_ordering::equal == std::strong_ordering::equivalent);

bool Isotope::operator==(const Isotope& other) const { return this->operator<=>(other) == std::strong_ordering::equal; }

bool Isotope::operator!=(const Isotope& other) const { return not this->operator==(other); }
}  // namespace Species
