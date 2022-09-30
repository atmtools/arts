#ifndef species_tags_h
#define species_tags_h

#include <algorithm>
#include <set>

#include "array.h"
#include "enums.h"
#include "isotopologues.h"
#include "matpackI.h"
#include "mystring.h"

namespace Species {
ENUMCLASS(TagType, unsigned char,
    Plain,
    Zeeman,
    Predefined,
    Cia,
    FreeElectrons,
    Particles,
    XsecFit,
    NoLines)
struct Tag {
  //! Molecular species index in Species::Isotopologues
  Index spec_ind{-1};

  //! The lower limit line center frequency in Hz.
  /*! If this is < 0 it means no lower limit. */
  Numeric lower_freq{-1};

  //! The upper line center frequency in Hz.
  /*! If this is < 0 it means no upper limit. */
  Numeric upper_freq{-1};
  
  //! Flag for the type
  TagType type{TagType::FINAL};

  //! 2nd CIA species index.
  /*! Contains the second CIA species that should be used for this tag. */
  Species cia_2nd_species{Species::FINAL};

  //! CIA dataset index.
  /*! A CIA file contains several datasets. This index specifies which one we want. */
  Index cia_dataset_index{-1};

  constexpr Tag() noexcept = default;
  
  // Documentation is with implementation.
  explicit Tag(String def);
  
  constexpr Tag(const IsotopeRecord& isot) noexcept :
  spec_ind(find_species_index(isot)),
  type(is_predefined_model(isot) ? TagType::Predefined : TagType::Plain) { /* Nothing to be done here. */
  }
  
  // Documentation is with implementation.
  [[nodiscard]] String Name() const;
  
  [[nodiscard]] constexpr IsotopeRecord Isotopologue() const noexcept {return spec_ind < 0 ? IsotopeRecord{} : Isotopologues[spec_ind];}
  
  constexpr void Isotopologue(const IsotopeRecord& ir) ARTS_NOEXCEPT {
    Index ind = find_species_index(ir);
    ARTS_ASSERT(ind < 0, "Bad species extracted from: ", ir)
    spec_ind = ind;
  }
  
  [[nodiscard]] constexpr Numeric Mass() const noexcept {return Isotopologue().mass;}
  
  [[nodiscard]] Numeric Q(Numeric T) const;
  
  [[nodiscard]] Numeric dQdT(Numeric T) const;
  
  [[nodiscard]] String FullName() const noexcept {return Isotopologue().FullName();}
  
  [[nodiscard]] constexpr Species Spec() const noexcept {return Isotopologue().spec;}
  
  [[nodiscard]] constexpr TagType Type() const noexcept {return type;}
  
  friend std::ostream& operator<<(std::ostream& os, const Tag& ot) {return os << ot.Name();}
  
  constexpr bool operator==(const Tag& other) const noexcept {
    return  other.spec_ind == spec_ind and
            other.lower_freq == lower_freq and
            other.upper_freq == upper_freq and
            other.type == type and
            other.cia_2nd_species == cia_2nd_species and
            other.cia_dataset_index == cia_dataset_index;
  }
  
  constexpr bool operator!=(const Tag& other) const noexcept {
    return not operator==(other);
  }
  
  [[nodiscard]] constexpr bool is_joker() const ARTS_NOEXCEPT {ARTS_ASSERT(spec_ind >= 0) return Joker == Isotopologue().isotname;}
};
} // namespace Species

using SpeciesTagType = Species::TagType;

using SpeciesTag = Species::Tag;

class ArrayOfSpeciesTag final : public Array<SpeciesTag> {
public:
  ArrayOfSpeciesTag() noexcept : Array<SpeciesTag>() {}
  explicit ArrayOfSpeciesTag(Index n) :  Array<SpeciesTag>(n) {}
  ArrayOfSpeciesTag(Index n, const SpeciesTag& fillvalue) : Array<SpeciesTag>(n, fillvalue) {}
  ArrayOfSpeciesTag(const ArrayOfSpeciesTag& A) = default;
  ArrayOfSpeciesTag(ArrayOfSpeciesTag&& A) noexcept : Array<SpeciesTag>(std::move(A)) {}
  explicit ArrayOfSpeciesTag(std::vector<SpeciesTag> x) : Array<SpeciesTag>(std::move(x)) {}
  
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
  
  ArrayOfSpeciesTag(String names);
  
  friend std::ostream& operator<<(std::ostream& os, const ArrayOfSpeciesTag& ot) {
    bool first = true;
    for (auto& x: ot) {if (not first) os << ' '; else first = false; os << x;}
    return os;
  }
  
  /*! Returns the species of the first elements, it is not allowed to have an empty list calling this */
  [[nodiscard]] Species::Species Species() const ARTS_NOEXCEPT {
    ARTS_ASSERT(size() not_eq 0, "Invalid ArrayOfSpeciesTag without any species")
    return operator[](0).Spec();
  }
  
//   /*! Returns the species of the first elements, it is not allowed to have an empty list calling this */
  [[nodiscard]] Species::TagType Type() const ARTS_NOEXCEPT {
    ARTS_ASSERT(size() not_eq 0, "Invalid ArrayOfSpeciesTag without any species")
    return operator[](0).Type();
  }
  
  [[nodiscard]] String Name() const;
  
  [[nodiscard]] bool Plain() const noexcept {
    return std::any_of(cbegin(), cend(), [](auto& spec){return spec.Type() == Species::TagType::Plain;});
  }
  
  [[nodiscard]] bool Zeeman() const noexcept {
    return std::any_of(cbegin(), cend(), [](auto& spec){return spec.Type() == Species::TagType::Zeeman;});
  }
  
  [[nodiscard]] bool RequireLines() const noexcept {
    return Plain() or Zeeman();
  }
  
  [[nodiscard]] bool FreeElectrons() const noexcept {
    return std::any_of(cbegin(), cend(), [](auto& spec){return spec.Type() == Species::TagType::FreeElectrons;});
  }
  
  [[nodiscard]] bool Particles() const noexcept {
    return std::any_of(cbegin(), cend(), [](auto& spec){return spec.Type() == Species::TagType::Particles;});
  }
};

using ArrayOfArrayOfSpeciesTag = Array<ArrayOfSpeciesTag>;

//! Struct to test of an ArrayOfArrayOfSpeciesTag contains a tagtype
struct SpeciesTagTypeStatus {
  bool Plain{false},
    Zeeman{false},
    Predefined{false},
    Cia{false},
    FreeElectrons{false},
    Particles{false},
    XsecFit{false},
    NoLines{false};
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
Index find_next_species(const ArrayOfArrayOfSpeciesTag& abs_species, Species::Species spec, Index i) noexcept;

/*! Find the first species of this type
 * 
 * @param[in] abs_species As WSV
 * @param[in] spec A Species
 * @return An index larger or equal to 0 pointing to the first species, or -1 if there's no such species
 */
Index find_first_species(const ArrayOfArrayOfSpeciesTag& abs_species, Species::Species spec) noexcept;

/*! Find the first species tag of this type
 * 
 * @param[in] abs_species As WSV
 * @param[in] tag A Species tag
 * @return [i, j] so that abs_species[i][j] is tag, or [-1, -1] if there's no such thing
 */
std::pair<Index, Index> find_first_species_tag(const ArrayOfArrayOfSpeciesTag& abs_species, const SpeciesTag& tag) noexcept;

/*! Find the first species tag of this type
 * 
 * @param[in] abs_species As WSV
 * @param[in] isot An isotopologue record
 * @return [i, j] so that abs_species[i][j] is the same isotopologue record as isot, or [-1, -1] if there's no such thing
 */
std::pair<Index, Index> find_first_isotologue(const ArrayOfArrayOfSpeciesTag& abs_species, const SpeciesIsotopeRecord& isot) noexcept;

/*!
 *  Checks on the correctness of the tags will be performed,
 *  e.g. free_electrons and particles species are only allowed once in
 *  abs_species.
 *  \param tags  Array of Array of SpeciesTag.
 *  \author Oliver Lemke
 *  \date   2013-04-23
 */
void check_abs_species(const ArrayOfArrayOfSpeciesTag& abs_species);

/*! Find species that requires line-by-line calculations
 * 
 * @param abs_species  As WSV
 * @return The set of unique species that requires line-by-line calculations
 */
std::set<Species::Species> lbl_species(const ArrayOfArrayOfSpeciesTag&) noexcept;

namespace Species {
/*! First VMR or 0
 * 
 * @param abs_species  As WSV
 * @param rtp_vmr  As WSV
 * @param spec A species
 */
Numeric first_vmr(const ArrayOfArrayOfSpeciesTag& abs_species,
                  const Vector& rtp_vmr,
                  const Species spec) ARTS_NOEXCEPT;
} // namespace Species

#endif  // species_tags_h
