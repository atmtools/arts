#include "species_tags.h"

#include <algorithm>
#include <cfloat>
#include <fast_float/fast_float.h>
#include <iomanip>
#include <iterator>
#include <string_view>
#include <system_error>

#include "debug.h"
#include "isotopologues.h"
#include "nonstd.h"
#include "partfun.h"
#include "species.h"

namespace Species {
namespace detail {
constexpr void trim(std::string_view &text) {
  while (text.size() and nonstd::isspace(text.front()))
    text.remove_prefix(1);
  while (text.size() and nonstd::isspace(text.back()))
    text.remove_suffix(1);
}

constexpr std::string_view next(std::string_view &text) {
  std::string_view next = text.substr(
      0, text.find_first_of('-', text.size() > 0 and text.front() == '-'));
  text.remove_prefix(std::min(text.size(), next.size()+1));
  trim(next);
  trim(text);
  return next;
}

constexpr std::string_view next_tag(std::string_view &text) {
  std::string_view next = text.substr(0, text.find_first_of(','));
  text.remove_prefix(std::min(text.size(), next.size()+1));
  trim(next);
  trim(text);
  return next;
}

constexpr Species spec(std::string_view part,
                       std::string_view orig [[maybe_unused]]) {
  Species s = fromShortName(part);
  ARTS_USER_ERROR_IF(not good_enum(s), "The species extraction from ",
                     std::quoted(part),
                     " is not good.  "
                     "The original tag reads ",
                     std::quoted(orig))
  return s;
}

constexpr Index isot(Species species, std::string_view isot,
                     std::string_view orig [[maybe_unused]]) {
  Index spec_ind = find_species_index(species, isot);
  ARTS_USER_ERROR_IF(spec_ind < 0, "Bad isotopologue: ", std::quoted(isot),
                     "\nValid options are:\n",
                     isotopologues_names(species),
                     "\nThe original tag reads ", std::quoted(orig))
  return spec_ind;
}

constexpr Numeric freq(std::string_view part,
                       std::string_view orig [[maybe_unused]]) {
  Numeric f;
  if (part.size() == 1 and part.front() == '*')
    f = -1;
  else {
    auto err =
        fast_float::from_chars(part.data(), part.data() + part.size(), f);
    ARTS_USER_ERROR_IF(err.ec == std::errc::invalid_argument,
                       "Invalid argument: ", std::quoted(part),
                       "\nThe original tag reads ", std::quoted(orig));
    ARTS_USER_ERROR_IF(err.ec == std::errc::result_out_of_range,
                       "Out-of-range: ", std::quoted(part),
                       "\nThe original tag reads ", std::quoted(orig));
  }
  return f;
}

constexpr void check(std::string_view text, std::string_view orig) {
  ARTS_USER_ERROR_IF(text.size(), "Parsing error.  The text ",
                     std::quoted(text),
                     " remains to be parsed at end of a complete species tag\n"
                     "The original tag reads ", std::quoted(orig))
}
} // namespace detail

SpeciesTag parse_tag(std::string_view text) {
  using namespace detail;
  trim(text);

  const std::string_view orig = text;
  SpeciesTag tag;

  // The first part is a species --- we do not know what to do with it yet
  const Species species = spec(next(text), orig);

  // Check if species name contains the special tag for
  // Faraday Rotation
  if (species == Species::free_electrons) {
    constexpr Index ind =
        find_species_index(IsotopeRecord(Species::free_electrons));
    tag.spec_ind = ind;
    tag.type = TagType::FreeElectrons;
  }

  // Check if species name contains the special tag for
  // Particles
  if (species == Species::particles) {
    constexpr Index ind = find_species_index(IsotopeRecord(Species::particles));
    tag.spec_ind = ind;
    tag.type = TagType::Particles;
  }
  
  // If there is no text remaining after the previous next(), then we are a
  // wild-tag species. Otherwise we have to process the tag a bit more
  if (text.size() == 0) {
    if (tag.type == TagType::FINAL) {
      tag.spec_ind = isot(species, Joker, orig);
      tag.type = TagType::Plain;
    }
    check(text, orig);
  } else {
    if (const std::string_view tag_key = next(text); tag_key.front() == 'Z') {
      tag.spec_ind = isot(species, next(text), orig);
      tag.type = TagType::Zeeman;
      ARTS_USER_ERROR_IF(
          is_predefined_model(Isotopologues[tag.spec_ind]),
          "Bad Zeeman tag using predefined model in tag: ", std::quoted(orig))
    } else if (tag_key == "CIA") {
      tag.spec_ind = isot(species, Joker, orig);
      tag.cia_2nd_species = spec(next(text), orig);
      tag.type = TagType::Cia;
      check(text, orig);
    } else if (tag_key == "XFIT") {
      tag.spec_ind = isot(species, Joker, orig);
      tag.type = TagType::XsecFit;
      check(text, orig);
    } else {
      tag.spec_ind = isot(species, tag_key, orig);
      tag.type = is_predefined_model(Isotopologues[tag.spec_ind])
                     ? TagType::Predefined
                     : TagType::Plain;
    }
  }

  if (text.size()) {
    tag.lower_freq = freq(next(text), orig);
    tag.upper_freq = freq(next(text), orig);
    ARTS_USER_ERROR_IF(tag.upper_freq >= 0 and tag.lower_freq >= 0 and
                           tag.upper_freq <= tag.lower_freq,
                       "Invalid frequency range [", tag.lower_freq, ", ",
                       tag.upper_freq, "]\nOriginal tag: ", std::quoted(orig))
    check(text, orig);
  }

  return tag;
}

Array<Tag> parse_tags(std::string_view text) {
  using namespace detail;
  trim(text);

  Array<Tag> tags;
  while (text.size()) {
    tags.emplace_back(parse_tag(next_tag(text)));
  }

  return tags;
}

Numeric Tag::Q(Numeric T) const {
  return PartitionFunctions::Q(T, Isotopologue());
}

Numeric Tag::dQdT(Numeric T) const {
  return PartitionFunctions::dQdT(T, Isotopologue());
}

Tag::Tag(std::string_view text) : Tag(parse_tag(text)) {}

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
    os << "CIA-" << toShortName(cia_2nd_species);

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
} // namespace Species

ArrayOfSpeciesTag::ArrayOfSpeciesTag(std::string_view text)
    : ArrayOfSpeciesTag(Species::parse_tags(text)) {
  ARTS_USER_ERROR_IF(
      size() and
          std::any_of(begin(), end(),
                      [front_species = front().Spec()](const SpeciesTag &tag) {
                        return tag.Spec() not_eq front_species;
                      }),
      "All species in a group of tags must be the same\n"
      "Your list of tags have been parsed as: ",
      *this, "\nThe original tags-list read ", std::quoted(text))

  ARTS_USER_ERROR_IF(
      size() and
          std::any_of(
              begin(), end(),
              [front_is_zeeman = front().type == Species::TagType::Zeeman](
                  const SpeciesTag &tag) {
                return front_is_zeeman ? tag.type != Species::TagType::Zeeman
                                       : tag.type == Species::TagType::Zeeman;
              }),
      "All species in a group of tags must be Zeeman-tagged if any one is\n"
      "Your list of tags have been parsed as: ",
      *this, "\nThe original tags-list read ", std::quoted(text)) {}
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
  }
  return os;
}
