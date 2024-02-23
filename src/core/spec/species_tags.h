#ifndef species_tags_h
#define species_tags_h

#include <enums.h>
#include <isotopologues.h>
#include <matpack.h>
#include <mystring.h>

#include <algorithm>
#include <set>

#include "species.h"

namespace Species {
struct Tag {
  //! Molecular species index in Species::Isotopologues
  Index spec_ind{-1};

  //! Flag for the type
  SpeciesTagType type{SpeciesTagType::Plain};

  //! 2nd CIA species index.
  /*! Contains the second CIA species that should be used for this tag. */
  SpeciesEnum cia_2nd_species{SpeciesEnum::Bath};

  constexpr Tag() noexcept = default;

  // Documentation is with implementation.
  explicit Tag(std::string_view text);

  constexpr Tag(const SpeciesIsotope& isot) noexcept
      : spec_ind(find_species_index(isot)),
        type(is_predefined_model(isot)
                 ? SpeciesTagType::Predefined
                 : SpeciesTagType::Plain) { /* Nothing to be done here. */
  }

  // Documentation is with implementation.
  [[nodiscard]] String Name() const;

  [[nodiscard]] constexpr SpeciesIsotope Isotopologue() const noexcept {
    return spec_ind < 0 ? SpeciesIsotope{} : Isotopologues[spec_ind];
  }

  constexpr void Isotopologue(const SpeciesIsotope& ir) ARTS_NOEXCEPT {
    Index ind = find_species_index(ir);
    ARTS_ASSERT(ind < 0, "Bad species extracted from: ", ir)
    spec_ind = ind;
  }

  [[nodiscard]] constexpr Numeric Mass() const noexcept {
    return Isotopologue().mass;
  }

  [[nodiscard]] Numeric Q(Numeric T) const;

  [[nodiscard]] Numeric dQdT(Numeric T) const;

  [[nodiscard]] String FullName() const noexcept {
    return Isotopologue().FullName();
  }

  [[nodiscard]] constexpr SpeciesEnum Spec() const noexcept {
    return Isotopologue().spec;
  }

  [[nodiscard]] constexpr SpeciesTagType Type() const noexcept { return type; }

  friend std::ostream& operator<<(std::ostream& os, const Tag& ot);

  [[nodiscard]] constexpr auto operator<=>(const Tag& other) const noexcept = default;

  [[nodiscard]] constexpr bool is_joker() const ARTS_NOEXCEPT {
    ARTS_ASSERT(spec_ind >= 0) return Joker == Isotopologue().isotname;
  }
};
}  // namespace Species

using SpeciesTagType = SpeciesTagType;

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
  ArrayOfSpeciesTag& operator=(SpeciesTag x) {
    std::fill(this->begin(), this->end(), x);
    return *this;
  }

  ArrayOfSpeciesTag& operator=(const ArrayOfSpeciesTag& A) {
    this->resize(A.size());
    std::copy(A.begin(), A.end(), this->begin());
    return *this;
  }

  ArrayOfSpeciesTag& operator=(ArrayOfSpeciesTag&& A) noexcept {
    Array<SpeciesTag>::operator=(std::move(A));
    return *this;
  }

  [[nodiscard]] bool operator==(const ArrayOfSpeciesTag& x) const {
    return std::ranges::equal(*this, x);
  }

  ArrayOfSpeciesTag(std::string_view text);

  friend std::ostream& operator<<(std::ostream& os,
                                  const ArrayOfSpeciesTag& ot) {
    bool first = true;
    for (auto& x : ot) {
      if (not first)
        os << ' ';
      else
        first = false;
      os << x;
    }
    return os;
  }

  /*! Returns the species of the first elements, it is not allowed to have an empty list calling this */
  [[nodiscard]] SpeciesEnum Species() const ARTS_NOEXCEPT {
    ARTS_ASSERT(size() not_eq 0,
                "Invalid ArrayOfSpeciesTag without any species")
    return operator[](0).Spec();
  }

  //   /*! Returns the species of the first elements, it is not allowed to have an empty list calling this */
  [[nodiscard]] SpeciesTagType Type() const ARTS_NOEXCEPT {
    ARTS_ASSERT(size() not_eq 0,
                "Invalid ArrayOfSpeciesTag without any species")
    return operator[](0).Type();
  }

  [[nodiscard]] String Name() const;

  [[nodiscard]] bool Plain() const noexcept {
    return std::any_of(cbegin(), cend(), [](auto& spec) {
      return spec.Type() == SpeciesTagType::Plain;
    });
  }

  [[nodiscard]] bool RequireLines() const noexcept { return Plain(); }

  [[nodiscard]] bool FreeElectrons() const noexcept {
    return std::any_of(cbegin(), cend(), [](auto& spec) {
      return spec.Spec() == SpeciesEnum::free_electrons;
    });
  }

  [[nodiscard]] bool Particles() const noexcept {
    return std::any_of(cbegin(), cend(), [](auto& spec) {
      return spec.Spec() == SpeciesEnum::particles;
    });
  }
};

using ArrayOfArrayOfSpeciesTag = Array<ArrayOfSpeciesTag>;

std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfSpeciesTag& a);

//! Struct to test of an ArrayOfArrayOfSpeciesTag contains a Speciestagtype
struct SpeciesTagTypeStatus {
  bool Plain{false}, Predefined{false}, Cia{false}, XsecFit{false};
  SpeciesTagTypeStatus(const ArrayOfArrayOfSpeciesTag& abs_species);
  friend std::ostream& operator<<(std::ostream&, SpeciesTagTypeStatus);
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
std::set<SpeciesEnum> lbl_species(
    const ArrayOfArrayOfSpeciesTag&) noexcept;

namespace Species {
/*! First VMR or 0
 * 
 * @param abs_species  As WSV
 * @param rtp_vmr  As WSV
 * @param spec A species
 */
Numeric first_vmr(const ArrayOfArrayOfSpeciesTag& abs_species,
                  const Vector& rtp_vmr,
                  const SpeciesEnum spec) ARTS_NOEXCEPT;

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

std::ostream& operator<<(std::ostream& os, const Array<SpeciesEnum>& a);
std::ostream& operator<<(std::ostream& os, const Array<Array<SpeciesEnum>>& a);
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
    std::size_t i = 1;
    for (auto& a : g) {
      out ^= (std::hash<SpeciesTag>{}(a) << i);
      i++;
    }
    return out;
  }
};
}  // namespace std

#endif  // species_tags_h
