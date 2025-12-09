#include "species_tags.h"

#include <debug.h>
#include <isotopologues.h>
#include <nonstd.h>
#include <xml.h>

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <string_view>

namespace Species {
Tag::Tag(const SpeciesIsotope& isot) noexcept
    : spec_ind(find_species_index(isot)),
      type(isot.is_predefined()
               ? SpeciesTagType::Predefined
               : SpeciesTagType::Plain) { /* Nothing to be done here. */ }

SpeciesIsotope Tag::Isotopologue() const noexcept {
  return spec_ind < 0 ? SpeciesIsotope{} : Isotopologues[spec_ind];
}

void Tag::Isotopologue(const SpeciesIsotope& ir) {
  Index ind = find_species_index(ir);
  assert(ind >= 0);
  spec_ind = ind;
}

Numeric Tag::Mass() const noexcept { return Isotopologue().mass; }

SpeciesEnum Tag::Spec() const noexcept { return Isotopologue().spec; }

SpeciesTagType Tag::Type() const noexcept { return type; }

bool Tag::is_joker() const {
  assert(spec_ind >= 0);
  return Joker == Isotopologue().isotname;
}

namespace {
constexpr void trim(std::string_view& text) {
  while (text.size() and nonstd::isspace(text.front())) text.remove_prefix(1);
  while (text.size() and nonstd::isspace(text.back())) text.remove_suffix(1);
}

constexpr std::string_view next(std::string_view& text) {
  std::string_view next = text.substr(
      0, text.find_first_of('-', text.size() > 0 and text.front() == '-'));
  text.remove_prefix(std::min(text.size(), next.size() + 1));
  trim(next);
  trim(text);
  return next;
}

constexpr SpeciesEnum spec(std::string_view part,
                           std::string_view orig [[maybe_unused]]) {
  return to<SpeciesEnum>(part);
}

constexpr Index isot(SpeciesEnum species,
                     std::string_view isot,
                     std::string_view orig [[maybe_unused]]) {
  Index spec_ind = find_species_index(species, isot);
  ARTS_USER_ERROR_IF(spec_ind < 0,
                     "Bad isotopologue: \"{}"
                     "\"\nValid options are:\n{}"
                     "\nThe original tag reads \"{}\"",
                     isot,
                     species,
                     orig)
  return spec_ind;
}

constexpr void check(std::string_view text, std::string_view orig) {
  ARTS_USER_ERROR_IF(
      text.size(),
      "Parsing error.  The text \"{}\""
      "\" remains to be parsed at end of a complete species tag\n"
      "The original tag reads \"{}\"",
      text,
      orig)
}
}  // namespace

String SpeciesTag::FullName() const noexcept {
  return Isotopologue().FullName();
}

namespace {
SpeciesTag parse_tag(std::string_view text) {
  trim(text);

  const std::string_view orig = text;
  SpeciesTag tag;

  // The first part is a species --- we do not know what to do with it yet
  SpeciesEnum species = SpeciesEnum::Bath;  // Default value
  try {
    species = spec(next(text), orig);
  } catch (std::exception& e) {
    throw std::runtime_error(std::format(
        "Error parsing species in tag \"{}\".\n\nValid SpeciesEnum are:\n{:Bq,}\n\nInternal error: {}",
        orig,
        enumtyps::SpeciesEnumTypes,
        e.what()));
  }

  // If there is no text remaining after the previous next(), then we are a
  // wild-tag species. Otherwise we have to process the tag a bit more
  if (text.size() == 0) {
    tag.spec_ind = isot(species, Joker, orig);
    tag.type     = SpeciesTagType::Plain;
    check(text, orig);
  } else {
    if (const std::string_view tag_key = next(text); tag_key == "CIA") {
      tag.spec_ind = isot(species, Joker, orig);
      try {
        tag.cia_2nd_species = spec(next(text), orig);
      } catch (std::exception& e) {
        throw std::runtime_error(std::format(
            "Error parsing second species in tag \"{}\".\n\nValid SpeciesEnum are:\n{:Bq,}\n\nInternal error: {}",
            orig,
            enumtyps::SpeciesEnumTypes,
            e.what()));
      }
      tag.type = SpeciesTagType::Cia;
      check(text, orig);
    } else if (tag_key == "XFIT") {
      tag.spec_ind = isot(species, Joker, orig);
      tag.type     = SpeciesTagType::XsecFit;
      check(text, orig);
    } else {
      try {
        tag.spec_ind = isot(species, tag_key, orig);
      } catch (std::exception& e) {
        throw std::runtime_error(std::format(
            "Error parsing isotopologue in tag \"{}\".\n\nValid isotopologues of species {} are:\n{:Bq,}\n\nInternal error: {}",
            orig,
            species,
            isotopologues(species),
            e.what()));
      }
      tag.type = Isotopologues[tag.spec_ind].is_predefined()
                     ? SpeciesTagType::Predefined
                     : SpeciesTagType::Plain;
    }
  }

  if (text.size()) {
    check(text, orig);
  }

  return tag;
}
}  // namespace

Tag::Tag(std::string_view text) : Tag(parse_tag(text)) {}

/** Return the full name of the tag.
 * 
 * \return The tag name as a string.
 */
String Tag::Name() const { return std::format("{}", *this); }
}  // namespace Species

SpeciesTagTypeStatus::SpeciesTagTypeStatus(
    const ArrayOfSpeciesTag& abs_species) {
  for (auto& tag : abs_species) {
    switch (tag.type) {
      case SpeciesTagType::Plain:      Plain = true; break;
      case SpeciesTagType::Predefined: Predefined = true; break;
      case SpeciesTagType::Cia:        Cia = true; break;
      case SpeciesTagType::XsecFit:    XsecFit = true; break;
    }
  }
}

void xml_io_stream<SpeciesTag>::write(std::ostream& os,
                                      const SpeciesTag& x,
                                      bofstream* pbofs,
                                      std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.Isotopologue(), pbofs);
  xml_write_to_stream(os, x.type, pbofs);
  xml_write_to_stream(os, x.cia_2nd_species, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<SpeciesTag>::read(std::istream& is,
                                     SpeciesTag& x,
                                     bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  SpeciesIsotope isot;
  xml_read_from_stream(is, isot, pbifs);
  xml_read_from_stream(is, x.type, pbifs);
  xml_read_from_stream(is, x.cia_2nd_species, pbifs);
  x.spec_ind = find_species_index(isot);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
