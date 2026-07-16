#pragma once

#include <compare.h>
#include <mystring.h>
#include <nonstd.h>
#include <species.h>
#include <xml_io_stream.h>

#include <compare>
#include <limits>
#include <string_view>

namespace Species {
/** Struct containing all information needed about one isotope */
struct Isotope {
  static constexpr std::string_view Joker = "*";

  SpeciesEnum spec{SpeciesEnum::Bath};

  //! A custom name that is unique for this Species type
  std::string_view isotname{Joker};

  //! The mass of the isotope in units of grams per mol.  It is Nan if not defined
  Numeric mass{std::numeric_limits<Numeric>::quiet_NaN()};

  //! The builtin ratio of the isotope to the normal isotopologue.  It is Nan if not defined.
  Numeric builtin_ratio{std::numeric_limits<Numeric>::quiet_NaN()};

  //! The degeneracy of states of the molecule.  It is -1 if not defined.
  Index gi{-1};

  static const Isotope& from_name(const std::string_view);

  [[nodiscard]] constexpr bool is_joker() const { return isotname == Joker; }

  [[nodiscard]] constexpr bool is_predefined() const { return not(nonstd::isdigit(isotname[0]) or is_joker()); }
  [[nodiscard]] constexpr bool is_normal() const { return nonstd::isdigit(isotname[0]) and not is_joker(); }

  [[nodiscard]] constexpr bool OK() const { return good_enum(spec); }
  [[nodiscard]] String         FullName() const;

  friend std::ostream& operator<<(std::ostream& os, const Isotope& ir);

  std::strong_ordering operator<=>(const Isotope&) const;
  bool                 operator==(const Isotope&) const;
  bool                 operator!=(const Isotope&) const;
};
}  // namespace Species

using SpeciesIsotope = Species::Isotope;

using ArrayOfSpeciesIsotope = Array<SpeciesIsotope>;

template <> struct std::formatter<SpeciesIsotope> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const SpeciesIsotope& v, FmtContext& ctx) const {
    const std::string_view quote = tags.quote();
    return tags.format(ctx, quote, v.FullName(), quote);
  }
};

template <> struct xml_io_stream<SpeciesIsotope> {
  static constexpr std::string_view type_name = "SpeciesIsotope"sv;

  static void write(std::ostream&         os,
                    const SpeciesIsotope& x,
                    bofstream*            pbofs = nullptr,
                    std::string_view      name  = ""sv);

  static void read(std::istream& is, SpeciesIsotope& x, bifstream* pbifs = nullptr);
};
