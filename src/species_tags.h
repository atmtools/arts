#ifndef species_tags_h
#define species_tags_h

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
  
  constexpr Numeric Mass() const noexcept {return Isotopologues[spec_ind].mass;}
  
  constexpr Numeric Q(Numeric T) const {return PartitionFunctions::Q(T, Isotopologues[spec_ind]);}
  
  constexpr Numeric dQdT(Numeric T) const {return PartitionFunctions::dQdT(T, Isotopologues[spec_ind]);}
  
  String SpeciesName() const noexcept {return Isotopologues[spec_ind].FullName();}
  
  const IsotopeRecord& Isotopologue() const noexcept {return Isotopologues[spec_ind];}
  
  friend std::ostream& operator<<(std::ostream& os, const Tag& ot) {return os << ot.Name();}
};
}  // Species

struct ArrayOfSpeciesTag2 : Array<Species::Tag> {
  ArrayOfSpeciesTag2(String names);
  friend std::ostream& operator<<(std::ostream& os, const ArrayOfSpeciesTag2& ot) {
    bool first = true;
    for (auto& x: ot) {if(not first) os << ' '; os << x; first = false;}
    return os;
  }
};

using ArrayOfArrayOfSpeciesTag2 = Array<ArrayOfSpeciesTag2>;

#endif  // species_tags_h
