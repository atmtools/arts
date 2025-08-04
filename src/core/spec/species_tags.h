#pragma once

#include <configtypes.h>
#include <enumsSpeciesTagType.h>
#include <isotopologues.h>
#include <mystring.h>

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

using SpeciesTag = Species::Tag;

class ArrayOfSpeciesTag final : public Array<SpeciesTag> {
 public:
  ArrayOfSpeciesTag() noexcept : Array<SpeciesTag>() {}
  explicit ArrayOfSpeciesTag(Index n) : Array<SpeciesTag>(n) {}
  ArrayOfSpeciesTag(Index n, const SpeciesTag& fillvalue)
      : Array<SpeciesTag>(n, fillvalue) {}
  ArrayOfSpeciesTag(const ArrayOfSpeciesTag& A) = default;
  ArrayOfSpeciesTag(ArrayOfSpeciesTag&& A) noexcept
      : Array<SpeciesTag>(std::move(A)) {}
  explicit ArrayOfSpeciesTag(std::vector<SpeciesTag> x)
      : Array<SpeciesTag>(std::move(x)) {}

  // Assignment operators:
  ArrayOfSpeciesTag& operator=(SpeciesTag x);

  ArrayOfSpeciesTag& operator=(const ArrayOfSpeciesTag& A);

  ArrayOfSpeciesTag& operator=(ArrayOfSpeciesTag&& A) noexcept;

  [[nodiscard]] bool operator==(const ArrayOfSpeciesTag& x) const;

  ArrayOfSpeciesTag(std::string_view text);

  /*! Returns the species of the first elements, it is not allowed to have an empty list calling this */
  [[nodiscard]] SpeciesEnum Species() const;

  //   /*! Returns the species of the first elements, it is not allowed to have an empty list calling this */
  [[nodiscard]] SpeciesTagType Type() const;

  [[nodiscard]] String Name() const;

  [[nodiscard]] bool Plain() const noexcept;

  [[nodiscard]] bool RequireLines() const noexcept;

  [[nodiscard]] bool FreeElectrons() const noexcept;

  [[nodiscard]] bool Particles() const noexcept;
};

using ArrayOfArrayOfSpeciesTag = Array<ArrayOfSpeciesTag>;

//! Struct to test of an ArrayOfArrayOfSpeciesTag contains a Speciestagtype
struct SpeciesTagTypeStatus {
  bool Plain{false}, Predefined{false}, Cia{false}, XsecFit{false};
  SpeciesTagTypeStatus(const ArrayOfArrayOfSpeciesTag& abs_species);
};

/*! Find the next species of this type inclusively after the start index
 * 
 * @param[in] abs_species As WSV
 * @param[in] spec A Species
 * @param[in] i The starting index in the outermost 
 * @return An index larger or equal to i pointing to the next species, or -1 if there's no next species
 */
Index find_next_species(const ArrayOfArrayOfSpeciesTag& abs_species,
                        SpeciesEnum spec,
                        Index i) noexcept;

/*! Find the first species of this type
 * 
 * @param[in] abs_species As WSV
 * @param[in] spec A Species
 * @return An index larger or equal to 0 pointing to the first species, or -1 if there's no such species
 */
Index find_first_species(const ArrayOfArrayOfSpeciesTag& abs_species,
                         SpeciesEnum spec) noexcept;

/*! Find the first species tag of this type
 * 
 * @param[in] abs_species As WSV
 * @param[in] tag A Species tag
 * @return [i, j] so that abs_species[i][j] is tag, or [-1, -1] if there's no such thing
 */
std::pair<Index, Index> find_first_species_tag(
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const SpeciesTag& tag) noexcept;

/*! Find the first species tag of this type
 * 
 * @param[in] abs_species As WSV
 * @param[in] isot An isotopologue record
 * @return [i, j] so that abs_species[i][j] is the same isotopologue record as isot, or [-1, -1] if there's no such thing
 */
std::pair<Index, Index> find_first_isotologue(
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const SpeciesIsotope& isot) noexcept;

/*! Find species that requires line-by-line calculations
 * 
 * @param abs_species  As WSV
 * @return The set of unique species that requires line-by-line calculations
 */
std::set<SpeciesEnum> lbl_species(const ArrayOfArrayOfSpeciesTag&) noexcept;

namespace Species {
/** Parse a list of species tags into an Array<Tag>
 *
 * This is the function call that ArrayOfSpeciesTag uses
 * to create itself from a string_view, but it also performs
 * vital error checking necessary for ARTS internals
 *
 * @param text A textual representation of comma-separated species tags
 * @return Array<Tag> List of species tags with no constraints
 */
Array<Tag> parse_tags(std::string_view text);
}  // namespace Species

namespace std {
//! Allow SpeciesTag to be used in hashes
template <>
struct hash<SpeciesTag> {
  std::size_t operator()(const SpeciesTag& g) const {
    return std::hash<Index>{}(g.spec_ind) ^
           (std::hash<SpeciesTagType>{}(g.type) << 20) ^
           (std::hash<SpeciesEnum>{}(g.cia_2nd_species) << 24);
  }
};

//! Allow ArrayOfSpeciesTag to be used in hashes
template <>
struct hash<ArrayOfSpeciesTag> {
  std::size_t operator()(const ArrayOfSpeciesTag& g) const {
    std::size_t out = 0;
    std::size_t i   = 1;
    for (auto& a : g) {
      out ^= (std::hash<SpeciesTag>{}(a) << i);
      i++;
    }
    return out;
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
struct std::formatter<ArrayOfSpeciesTag>
    : std::formatter<std::vector<SpeciesTag>> {};

template <>
struct xml_io_stream<SpeciesTag> {
  static constexpr std::string_view type_name = "SpeciesTag"sv;

  static void write(std::ostream& os,
                    const SpeciesTag& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, SpeciesTag& x, bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<ArrayOfSpeciesTag> {
  static constexpr std::string_view type_name = "ArrayOfSpeciesTag"sv;

  static void write(std::ostream& os,
                    const ArrayOfSpeciesTag& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   ArrayOfSpeciesTag& x,
                   bifstream* pbifs = nullptr);
};
