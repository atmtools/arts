#ifndef species_tags_h
#define species_tags_h

#include <algorithm>

#include "array.h"
#include "mystring.h"
#include "partfun.h"

namespace Species {
ENUMCLASS(TagType, unsigned char,
    Plain,
    Zeeman,
    Predefined,
    Cia,
    FreeElectrons,
    Particles,
    HitranXsec,
    NoLines)

struct Tag {
  //! Molecular species index in Species::Isotopologues
  Index spec_ind;

  //! The lower limit line center frequency in Hz.
  /*! If this is < 0 it means no lower limit. */
  Numeric lower_freq;

  //! The upper line center frequency in Hz.
  /*! If this is < 0 it means no upper limit. */
  Numeric upper_freq;
  
  //! Flag for the type
  TagType type;

  //! 2nd CIA species index.
  /*! Contains the second CIA species that should be used for this tag. */
  Species cia_2nd_species;

  //! CIA dataset index.
  /*! A CIA file contains several datasets. This index specifies which one we want. */
  Index cia_dataset_index;

  constexpr Tag() noexcept
  : spec_ind(-1),
    lower_freq(0.),
    upper_freq(0.),
    type(TagType::FINAL),
    cia_2nd_species(Species::FINAL),
    cia_dataset_index(-1) { /* Nothing to be done here. */
  }
  
  // Documentation is with implementation.
  explicit Tag(String def);
  
  // Documentation is with implementation.
  String Name() const;
  
  constexpr const IsotopeRecord& Isotopologue() const noexcept {return Isotopologues[spec_ind];}
  
  constexpr void Isotopologue(const IsotopeRecord& ir) ARTS_NOEXCEPT {
    Index ind = find_species_index(ir);
    ARTS_ASSERT(ind < 0, "Bad species extracted from: ", ir)
    spec_ind = ind;
  }
  
  constexpr Numeric Mass() const noexcept {return Isotopologue().mass;}
  
  constexpr Numeric Q(Numeric T) const {return PartitionFunctions::Q(T, Isotopologue());}
  
  constexpr Numeric dQdT(Numeric T) const {return PartitionFunctions::dQdT(T, Isotopologue());}
  
  String FullName() const noexcept {return Isotopologue().FullName();}
  
  constexpr Species Spec() const noexcept {return Isotopologue().spec;}
  
  constexpr TagType Type() const noexcept {return type;}
  
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
  
  constexpr bool is_joker() const ARTS_NOEXCEPT {ARTS_ASSERT(spec_ind >= 0) return Joker == Isotopologue().isotname;}
};
}  // Species

using SpeciesTagType = Species::TagType;

using SpeciesTag = Species::Tag;

class ArrayOfSpeciesTag final : public Array<SpeciesTag> {
public:
  ArrayOfSpeciesTag() noexcept : Array<SpeciesTag>() {}
  explicit ArrayOfSpeciesTag(Index n) :  Array<SpeciesTag>(n) {}
  ArrayOfSpeciesTag(Index n, const SpeciesTag& fillvalue) : Array<SpeciesTag>(n, fillvalue) {}
  ArrayOfSpeciesTag(const ArrayOfSpeciesTag& A) : Array<SpeciesTag>(A) {}
  ArrayOfSpeciesTag(ArrayOfSpeciesTag&& A) noexcept : Array<SpeciesTag>(std::move(A)) {}
  
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
  Species::Species Species() const ARTS_NOEXCEPT {
    ARTS_ASSERT(size() not_eq 0, "Invalid ArrayOfSpeciesTag without any species")
    return operator[](0).Spec();
  }
  
  /*! Returns the species of the first elements, it is not allowed to have an empty list calling this */
  Species::TagType Type() const ARTS_NOEXCEPT {
    ARTS_ASSERT(size() not_eq 0, "Invalid ArrayOfSpeciesTag without any species")
    return operator[](0).Type();
  }
  
  String Name() const {
    String out="";
    for (auto& x: *this) out += x.Name();
    return out;
  }
  
  bool Zeeman() const noexcept {
    return std::any_of(cbegin(), cend(), [](auto& spec){return spec.Type() == Species::TagType::Zeeman;});
  }
  
  bool FreeElectrons() const noexcept {
    return std::any_of(cbegin(), cend(), [](auto& spec){return spec.Type() == Species::TagType::FreeElectrons;});
  }
  
  bool Particles() const noexcept {
    return std::any_of(cbegin(), cend(), [](auto& spec){return spec.Type() == Species::TagType::Particles;});
  }
};

using ArrayOfArrayOfSpeciesTag = Array<ArrayOfSpeciesTag>;

Index find_next_species(const ArrayOfArrayOfSpeciesTag&, Species::Species, Index) noexcept;

Index find_first_species(const ArrayOfArrayOfSpeciesTag&, Species::Species) noexcept;

Index find_first_species_tag(const ArrayOfArrayOfSpeciesTag& specs, const SpeciesTag& tag) noexcept;

/*!
 *  Checks on the correctness of the tags will be performed,
 *  e.g. free_electrons and particles species are only allowed once in
 *  abs_species.
 *  \param tags  Array of Array of SpeciesTag.
 *  \author Oliver Lemke
 *  \date   2013-04-23
 */
void check_abs_species(const ArrayOfArrayOfSpeciesTag& abs_species);

#endif  // species_tags_h
