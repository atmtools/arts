#pragma once

#include <configtypes.h>
#include <enumsSpeciesTagType.h>
#include <isotopologues.h>
#include <mystring.h>

#include <boost/container_hash/hash.hpp>
#include <set>

namespace Species {
struct Tag {
  //! Molecular species index in Species::Isotopologues
  Index spec_ind{"Ar"_isot_index};

  //! Flag for the type
  SpeciesTagType type{SpeciesTagType::Plain};

  //! 2nd CIA species index.
  /*! Contains the second CIA species that should be used for this tag. */
  SpeciesEnum cia_2nd_species{SpeciesEnum::Bath};

  constexpr Tag() noexcept = default;

  // Documentation is with implementation.
  explicit Tag(std::string_view text);

  Tag(const SpeciesIsotope& isot) noexcept;

  // Documentation is with implementation.
  [[nodiscard]] String Name() const;

  [[nodiscard]] SpeciesIsotope Isotopologue() const noexcept;

  void Isotopologue(const SpeciesIsotope& ir);

  [[nodiscard]] Numeric Mass() const noexcept;

  [[nodiscard]] String FullName() const noexcept;

  [[nodiscard]] SpeciesEnum Spec() const noexcept;

  [[nodiscard]] SpeciesTagType Type() const noexcept;

  [[nodiscard]] constexpr auto operator<=>(const Tag& other) const noexcept =
      default;

  [[nodiscard]] bool is_joker() const;
};
}  // namespace Species

using SpeciesTag               = Species::Tag;
using ArrayOfSpeciesTag        = Array<SpeciesTag>;
using ArrayOfArrayOfSpeciesTag = Array<ArrayOfSpeciesTag>;

//! Struct to test of an ArrayOfArrayOfSpeciesTag contains a Speciestagtype
struct SpeciesTagTypeStatus {
  bool Plain{false}, Predefined{false}, Cia{false}, XsecFit{false};
  SpeciesTagTypeStatus(const ArrayOfSpeciesTag& abs_species);
};

namespace std {
//! Allow SpeciesTag to be used in hashes
template <>
struct hash<SpeciesTag> {
  std::size_t operator()(const SpeciesTag& g) const {
    std::size_t seed = 0;

    boost::hash_combine(seed, std::hash<Index>{}(g.spec_ind));
    boost::hash_combine(seed, std::hash<SpeciesTagType>{}(g.type));
    boost::hash_combine(seed, std::hash<SpeciesEnum>{}(g.cia_2nd_species));

    return seed;
  }
};
}  // namespace std

template <>
struct std::formatter<SpeciesTag> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SpeciesTag& v, FmtContext& ctx) const {
    const std::string_view quote = tags.quote();
    tags.format(ctx, quote, v.Spec());

    switch (v.type) {
      case SpeciesTagType::Plain:
        if (not v.is_joker()) tags.format(ctx, '-', v.Isotopologue().isotname);
        break;
      case SpeciesTagType::Predefined:
        tags.format(ctx, '-', v.Isotopologue().isotname);
        break;
      case SpeciesTagType::Cia:
        tags.format(ctx, "-CIA-", v.cia_2nd_species);
        break;
      case SpeciesTagType::XsecFit: tags.format(ctx, "-XFIT"sv); break;
    }

    tags.format(ctx, quote);
    return ctx.out();
  }
};

template <>
struct xml_io_stream<SpeciesTag> {
  static constexpr std::string_view type_name = "SpeciesTag"sv;

  static void write(std::ostream& os,
                    const SpeciesTag& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, SpeciesTag& x, bifstream* pbifs = nullptr);
};
