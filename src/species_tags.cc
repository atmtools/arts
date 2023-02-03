#include "species_tags.h"

#include <cfloat>
#include <iterator>
#include <string_view>

#include "debug.h"
#include "partfun.h"
#include "matpack_concepts.h" 

namespace Species {

Numeric Tag::Q(Numeric T) const {
  return PartitionFunctions::Q(T, Isotopologue());
}

Numeric Tag::dQdT(Numeric T) const {
  return PartitionFunctions::dQdT(T, Isotopologue());
}

Tag::Tag(String def) : type(TagType::Plain) {
  // Save input string for error messages:
  String def_original = def;

  // Name of species and isotopologue (aux variables):
  String name, isoname;
  // Aux index:
  Index n;

  // We cannot set a default value for the isotopologue, because the
  // default should be `ALL' and the value for `ALL' depends on the
  // species.

  // Extract the species name:
  n = def.find('-');  // find the '-'
  if (n != def.npos) {
    name = def.substr(0, n);  // Extract before '-'
    def.erase(0, n + 1);      // Remove from def
  } else {
    // n==def.npos means that def does not contain a '-'. In that
    // case we assume that it contains just a species name and
    // nothing else
    name = def;
    def = "";
  }

  // Remove whitespace
  name.trim();

  // Obtain species index from species name.
  // (This will also remove possible whitespace.)
  const Species spec = fromShortName(name);
  ARTS_USER_ERROR_IF(not good_enum(spec),
                     "Bad species name: ",
                     name,
                     " extracted from: ",
                     def_original)

  // Check if species name contains the special tag for
  // Faraday Rotation
  if (spec == Species::free_electrons) {
    constexpr Index ind =
        find_species_index(IsotopeRecord(Species::free_electrons));
    spec_ind = ind;
    type = TagType::FreeElectrons;
    return;
  }

  // Check if species name contains the special tag for
  // Particles
  if (spec == Species::particles) {
    constexpr Index ind = find_species_index(IsotopeRecord(Species::particles));
    spec_ind = ind;
    type = TagType::Particles;
    return;
  }

  // Set "all" species per default by leaving the joker in
  spec_ind = find_species_index(IsotopeRecord(spec));
  ARTS_USER_ERROR_IF(
      spec_ind < 0, "Bad species extracted: ", spec, " from ", def_original)
  if (0 == def.nelem()) {
    return;
  }

  // Extract the isotopologue name/Zeeman flag:
  n = def.find('-');  // find the '-'
  if (n != def.npos) {
    isoname = def.substr(0, n);  // Extract before '-'
    def.erase(0, n + 1);         // Remove from def

    if ("Z" == isoname) {
      type = TagType::Zeeman;
      // Zeeman flag was present, now extract the isotopologue name:
      n = def.find('-');  // find the '-'
      if (n != def.npos) {
        isoname = def.substr(0, n);  // Extract before '-'
        def.erase(0, n + 1);         // Remove from def
      } else {
        // n==def.npos means that def does not contain a '-'. In that
        // case we assume that it contains just the isotopologue name and
        // nothing else.
        isoname = def;
        def = "";
      }
    }

    if ("XFIT" == isoname) {
      type = TagType::XsecFit;
      // Hitran Xsec flag was present, now extract the isotopologue name:
      n = def.find('-');  // find the '-'
      if (n != def.npos) {
        isoname = def.substr(0, n);  // Extract before '-'
        def.erase(0, n + 1);         // Remove from def
      } else {
        // n==def.npos means that def does not contain a '-'. In that
        // case we assume that it contains just the isotopologue name and
        // nothing else.
        isoname = def;
        def = "";
      }
    }
  } else {
    // n==def.npos means that def does not contain a '-'. In that
    // case we assume that it contains just the isotopologue name or
    // Zeeman flag and nothing else.
    isoname = def;
    def = "";
    if ("Z" == isoname) {
      type = TagType::Zeeman;
      // This means that there is nothing else to parse. Apparently
      // the user wants all isotopologues and no frequency limits.
      return;
    }
    if ("XFIT" == isoname) {
      type = TagType::XsecFit;
      // This means that there is nothing else to parse. Apparently
      // the user wants all isotopologues and no frequency limits.
      return;
    }
  }

  if (Joker == isoname) {
    // The user wants all isotopologues. This already set
  } else if ("CIA" == isoname)  // Check for "cia":
  {
    // The user wants this to use the CIA catalog:
    type = TagType::Cia;

    // We have to read in the second species, and the dataset index
    n = def.find('-');  // find the '-'

    ARTS_USER_ERROR_IF(
        n == def.npos,
        "Invalid species tag ",
        def_original,
        ".\n"
        "I am missing a minus sign (and a dataset index after that.)")

    String otherspec = def.substr(0, n);  // Extract before '-'
    def.erase(0, n + 1);                  // Remove from def

    cia_2nd_species = fromShortName(otherspec);
    ARTS_USER_ERROR_IF(not good_enum(cia_2nd_species),
                       "Bad species name: ",
                       otherspec,
                       " extracted from: ",
                       def_original)

    // Convert remaining def to dataset index.

    // Check that everything remaining is just numbers.
    for (Index i = 0; i < def.nelem(); ++i)
      ARTS_USER_ERROR_IF(!isdigit(def[i]),
                         "Invalid species tag ",
                         def_original,
                         ".\n"
                         "The tag should end with a dataset index")

    // Do the conversion from string to index:
    std::istringstream is(def);
    is >> cia_dataset_index;

    def = "";
  } else {
    spec_ind = find_species_index(IsotopeRecord(spec, isoname));

    // Check if we found a matching isotopologue:
    ARTS_USER_ERROR_IF(spec_ind < 0,
                       "Isotopologue ",
                       isoname,
                       " is not a valid isotopologue or "
                       "absorption model for species ",
                       name,
                       ".\n"
                       "Valid options are:\n",
                       isotopologues_names(spec))

    // Check if the found isotopologue represents a predefined model
    // (continuum or full absorption model) and set the type accordingly:
    if (is_predefined_model(Isotopologue()))
      type = TagType::Predefined;
  }

  if (0 == def.nelem()) {
    // This means that there is nothing else to parse. Apparently
    // the user wants no frequency limits.  Frequency defaults are
    // already set, so we can return directly.

    return;
  }

  ARTS_USER_ERROR_IF(def[0] != Joker[0] && !isdigit(def[0]),
                     "Expected frequency limits, but got \"",
                     def,
                     "\"")

  // Look for the two frequency limits:

  // Extract first frequency
  n = def.find('-');  // find the '-'
  if (n != def.npos) {
    // Frequency as a String:
    String fname;
    fname = def.substr(0, n);  // Extract before '-'
    def.erase(0, n + 1);       // Remove from def

    // Check for joker:
    if (Joker == fname) {
      // The default for lower_freq is already -1, meaning `ALL'.
      // So there is nothing to do here.
    } else {
      ARTS_USER_ERROR_IF(!isdigit(fname[0]),
                         "Expected frequency limit, but got \"",
                         fname,
                         "\"")
      // Convert to Numeric:
      char* endptr;
      lower_freq = strtod(fname.c_str(), &endptr);
      ARTS_USER_ERROR_IF(endptr != fname.c_str() + fname.nelem(),
                         "Error parsing frequency limit \"",
                         fname,
                         "\"")
    }
  } else {
    // n==def.npos means that def does not contain a '-'. In this
    // case that is not allowed!
    ARTS_USER_ERROR(
        "You must either specify both frequency limits\n"
        "(at least with jokers), or none.");
  }

  // Now there should only be the upper frequency left in def.
  // Check for joker:
  if (Joker== def) {
    // The default for upper_freq is already -1, meaning `ALL'.
    // So there is nothing to do here.
  } else {
    ARTS_USER_ERROR_IF(
        !isdigit(def[0]), "Expected frequency limit, but got \"", def, "\"")
    // Convert to Numeric:
    char* endptr;
    upper_freq = strtod(def.c_str(), &endptr);
    ARTS_USER_ERROR_IF(endptr != def.c_str() + def.nelem(),
                       "Error parsing frequency limit \"",
                       def,
                       "\"")
  }
}

//! Return the full name of the tag.
/*! 
 * Examples:
 * 
 * \verbatim
 * O3-*-*-*         : All O3 lines
 * O3-nl            : O3, but without any lines
 * O3-666-*-*       : All O3-666 lines
 * O3-*-500e9-501e9 : All O3 lines between 500 and 501 GHz.
 * \endverbatim
 * 
 * \return The tag name as a string.
 */
String Tag::Name() const {
  std::ostringstream os;

  // First the species name:
  os << toShortName(Isotopologue().spec) << "-";

  // Is this a CIA tag?
  if (type == TagType::Cia) {
    os << "CIA-" << toShortName(cia_2nd_species) << "-" << cia_dataset_index;

  } else if (type == TagType::FreeElectrons || type == TagType::Particles) {
    os << toShortName(Isotopologue().spec);
  }
  // Hitran Xsec flag.
  else if (type == TagType::XsecFit) {
    os << "XFIT";
  } else {
    // Zeeman flag.
    if (type == TagType::Zeeman) os << "Z-";

    // Now the isotopologue. Can be a single isotopologue or ALL.
    os << Isotopologue().isotname << '-';

    // Now the frequency limits, if there are any. For this we first
    // need to determine the floating point precision.

    // Determine the precision, depending on whether Numeric is double
    // or float:
    constexpr int precision = std::same_as<Numeric, double> ? DBL_DIG : FLT_DIG;

    if (0 > lower_freq) {
      // lower_freq < 0 means no lower limit.
      os << Joker << '-';
    } else {
      os << std::setprecision(precision);
      os << lower_freq << "-";
    }

    if (0 > upper_freq) {
      // upper_freq < 0 means no upper limit.
      os << Joker;
    } else {
      os << std::setprecision(precision);
      os << upper_freq;
    }
  }
  return os.str();
}
}  // namespace Species

ArrayOfSpeciesTag::ArrayOfSpeciesTag(String these_names) {
  // There can be a comma separated list of tag definitions, so we
  // need to break the String apart at the commas.
  ArrayOfString tag_def;

  bool go_on = true;

  while (go_on) {
    //          Index n = find_first( these_names, ',' );
    Index n = these_names.find(',');
    if (n == these_names.npos)  // Value npos indicates not found.
    {
      // There are no more commas.
      //              cout << "these_names: (" << these_names << ")\n";
      tag_def.push_back(these_names);
      go_on = false;
    } else {
      tag_def.push_back(these_names.substr(0, n));
      these_names.erase(0, n + 1);
    }
  }
  // tag_def now holds the different tag Strings for this group.

  // Set size to zero, in case the method has been called before.
  resize(0);

  for (Index s = 0; s < tag_def.nelem(); ++s) {
    Species::Tag this_tag(tag_def[s]);

    // Safety checks:
    if (s > 0) {
      // Tags inside a group must belong to the same species.
      ARTS_USER_ERROR_IF(
          front().Isotopologue().spec != this_tag.Isotopologue().spec,
          "Tags in a tag group must belong to the same species.");

      // Zeeman tags and plain line by line tags must not be mixed. (Because
      // there can be only one line list per tag group.)
      ARTS_USER_ERROR_IF(((front().type == Species::TagType::Zeeman) &&
                          (this_tag.type == Species::TagType::Plain)) ||
                             ((front().type == Species::TagType::Plain) &&
                              (this_tag.type == Species::TagType::Zeeman)),
                         "Zeeman tags and plain line-by-line tags must "
                         "not be mixed in the same tag group.");
    }

    push_back(this_tag);
  }
}

Index find_next_species(const ArrayOfArrayOfSpeciesTag& specs,
                        Species::Species spec,
                        Index i) noexcept {
  const Index n = specs.nelem();
  for (; i < n; i++)
    if (specs[i].Species() == spec) return i;
  return -1;
}

Index find_first_species(const ArrayOfArrayOfSpeciesTag& specs,
                         Species::Species spec) noexcept {
  return find_next_species(specs, spec, 0);
}

std::pair<Index, Index> find_first_species_tag(
    const ArrayOfArrayOfSpeciesTag& specs, const SpeciesTag& tag) noexcept {
  for (Index i = 0; i < specs.nelem(); i++) {
    if (auto ptr = std::find(specs[i].cbegin(), specs[i].cend(), tag);
        ptr not_eq specs[i].cend())
      return {i, std::distance(specs[i].cbegin(), ptr)};
  }
  return {-1, -1};
}

std::pair<Index, Index> find_first_isotologue(
    const ArrayOfArrayOfSpeciesTag& specs,
    const SpeciesIsotopeRecord& isot) noexcept {
  for (Index i = 0; i < specs.nelem(); i++) {
    if (auto ptr =
            std::find_if(specs[i].cbegin(),
                         specs[i].cend(),
                         [&](auto& tag) { return tag.Isotopologue() == isot; });
        ptr not_eq specs[i].cend())
      return {i, std::distance(specs[i].cbegin(), ptr)};
  }
  return {-1, -1};
}

void check_abs_species(const ArrayOfArrayOfSpeciesTag& abs_species) {
  Index num_free_electrons = 0;
  for (Index i = 0; i < abs_species.nelem(); ++i) {
    bool has_free_electrons = false;
    bool has_particles = false;
    bool has_hitran_xsec = false;
    for (Index s = 0; s < abs_species[i].nelem(); ++s) {
      if (abs_species[i][s].Type() == Species::TagType::FreeElectrons) {
        num_free_electrons++;
        has_free_electrons = true;
      }

      if (abs_species[i][s].Type() == Species::TagType::Particles) {
        has_particles = true;
      }

      if (abs_species[i][s].Type() == Species::TagType::XsecFit) {
        has_hitran_xsec = true;
      }
    }

    ARTS_USER_ERROR_IF(abs_species[i].nelem() > 1 && has_free_electrons,
                       "'free_electrons' must not be combined "
                       "with other tags in the same group.");
    ARTS_USER_ERROR_IF(num_free_electrons > 1,
                       "'free_electrons' must not be defined "
                       "more than once.");

    ARTS_USER_ERROR_IF(abs_species[i].nelem() > 1 && has_particles,
                       "'particles' must not be combined "
                       "with other tags in the same group.");

    ARTS_USER_ERROR_IF(abs_species[i].nelem() > 1 && has_hitran_xsec,
                       "'hitran_xsec' must not be combined "
                       "with other tags in the same group.");
  }
}

String ArrayOfSpeciesTag::Name() const {
  String out = "";
  bool first = true;
  for (auto& x : *this) {
    if (not first) out += ", ";
    out += x.Name();
    first = false;
  }
  return out;
}

std::set<Species::Species> lbl_species(
    const ArrayOfArrayOfSpeciesTag& abs_species) noexcept {
  std::set<Species::Species> unique_species;
  for (auto& specs : abs_species) {
    if (specs.RequireLines()) unique_species.insert(specs.front().Spec());
  }
  return unique_species;
}

Numeric Species::first_vmr(const ArrayOfArrayOfSpeciesTag& abs_species,
                           const Vector& rtp_vmr,
                           const Species spec) ARTS_NOEXCEPT {
  ARTS_ASSERT(abs_species.nelem() == rtp_vmr.nelem())

  auto pos =
      std::find_if(abs_species.begin(),
                   abs_species.end(),
                   [spec](const ArrayOfSpeciesTag& tag_group) {
                     return tag_group.nelem() and tag_group.Species() == spec;
                   });
  return pos == abs_species.end()
             ? 0
             : rtp_vmr[std::distance(abs_species.begin(), pos)];
}

SpeciesTagTypeStatus::SpeciesTagTypeStatus(const ArrayOfArrayOfSpeciesTag& abs_species) {
  for (auto& species_list : abs_species) {
    for (auto& tag : species_list) {
      switch (tag.type) {
        case Species::TagType::Plain:
          Plain = true;
          break;
        case Species::TagType::Zeeman:
          Zeeman = true;
          break;
        case Species::TagType::Predefined:
          Predefined = true;
          break;
        case Species::TagType::Cia:
          Cia = true;
          break;
        case Species::TagType::FreeElectrons:
          FreeElectrons = true;
          break;
        case Species::TagType::Particles:
          Particles = true;
          break;
        case Species::TagType::XsecFit:
          XsecFit = true;
          break;
        case Species::TagType::NoLines:
          NoLines = true;
          break;
        case Species::TagType::FINAL: { /* leave last */
        }
      }
    }
  }
}

std::ostream& operator<<(std::ostream& os, SpeciesTagTypeStatus val) {
  Species::TagType x{Species::TagType::FINAL};
  switch (x) {
    case Species::TagType::FINAL:
      os << "Species tag types:\n";
      [[fallthrough]];
    case Species::TagType::Plain:
      os << "    Plain:            " << val.Plain << '\n';
      [[fallthrough]];
    case Species::TagType::Zeeman:
      os << "    Zeeman:           " << val.Zeeman << '\n';
      [[fallthrough]];
    case Species::TagType::Predefined:
      os << "    Predefined: " << val.Predefined << '\n';
      [[fallthrough]];
    case Species::TagType::Cia:
      os << "    Cia:              " << val.Cia << '\n';
      [[fallthrough]];
    case Species::TagType::FreeElectrons:
      os << "    FreeElectrons:    " << val.FreeElectrons << '\n';
      [[fallthrough]];
    case Species::TagType::Particles:
      os << "    Particles:        " << val.Particles << '\n';
      [[fallthrough]];
    case Species::TagType::XsecFit:
      os << "    XsecFit:       " << val.XsecFit << '\n';
      [[fallthrough]];
    case Species::TagType::NoLines:
      os << "    NoLines:          " << val.NoLines;
  }
  return os;
}
