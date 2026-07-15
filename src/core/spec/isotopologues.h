#pragma once

#include <auto_isotopologues.h>
#include <compare.h>
#include <mystring.h>
#include <nonstd.h>
#include <species.h>
#include <xml.h>

#include <array>
#include <limits>
#include <string_view>

namespace Species {
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
  return find_species_index(
      s.substr(0, minus),
      minus == s.npos ? Isotope::Joker : s.substr(minus + 1));
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

  constexpr IsotopologueRatios() = default;

  Numeric operator[](const Index spec_ind) const;

  Numeric operator[](const Isotope& ir) const;

  friend std::ostream& operator<<(std::ostream& os,
                                  const IsotopologueRatios& iso_rat);

  [[nodiscard]] std::vector<std::string> valueless_isotopes() const;
};

const IsotopologueRatios& isotopologue_ratiosInitFromBuiltin();

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

using SpeciesIsotopologueRatios = Species::IsotopologueRatios;

consteval SpeciesIsotope operator""_isot(const char* x, std::size_t) {
  return Species::select(x);
}

consteval Index operator""_isot_index(const char* x, std::size_t) {
  return Species::find_species_index(x);
}

template <>
struct std::hash<SpeciesIsotope> {
  static std::size_t operator()(const SpeciesIsotope& g) {
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
          ctx, "\n"sv, Species::Isotopologues[i].FullName(), sep, v.data[i]);
    }

    return ctx.out();
  }
};
