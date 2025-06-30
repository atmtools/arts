#include "optproperties.h"

#include <format>
#include <iterator>

std::string_view PTypeToString(const PType ptype) {
  switch (ptype) {
    case PTYPE_GENERAL:     return "general"sv; break;
    case PTYPE_TOTAL_RND:   return "totally_random"sv; break;
    case PTYPE_AZIMUTH_RND: return "azimuthally_random"sv; break;
    default:
      ARTS_USER_ERROR(
          "Internal error: Cannot map PType enum value {} to String.", ptype)
      break;
  }

  return "unknown PType"sv;
}
PType PType2FromString(const std::string_view ptype_string) {
  PType ptype;
  if (ptype_string == "general"sv)
    ptype = PTYPE_GENERAL;
  else if (ptype_string == "macroscopically_isotropic"sv)
    ptype = PTYPE_TOTAL_RND;
  else if (ptype_string == "horizontally_aligned"sv)
    ptype = PTYPE_AZIMUTH_RND;
  else {
    ARTS_USER_ERROR(
        "Unknown ptype: {}"
        "\n"
        "Valid types are: general, macroscopically_isotropic and "
        "horizontally_aligned.",
        ptype_string)
  }

  return ptype;
}

PType PTypeFromString(const std::string_view ptype_string) {
  PType ptype;
  if (ptype_string == "general"sv)
    ptype = PTYPE_GENERAL;
  else if (ptype_string == "totally_random"sv)
    ptype = PTYPE_TOTAL_RND;
  else if (ptype_string == "azimuthally_random"sv)
    ptype = PTYPE_AZIMUTH_RND;
  else {
    ARTS_USER_ERROR(
        "Unknown ptype: {}"
        "\n"
        "Valid types are: general, totally_random and "
        "azimuthally_random.",
        ptype_string)
  }

  return ptype;
}

std::string std::formatter<SingleScatteringData>::to_string(
    const SingleScatteringData& v) const {
  const std::string_view sep   = tags.sep(true);
  const std::string_view quote = tags.quote();

  std::string out = tags.vformat(v.ptype,
                                 sep,
                                 quote,
                                 v.description,
                                 quote,
                                 sep,
                                 v.f_grid,
                                 sep,
                                 v.T_grid,
                                 sep,
                                 v.za_grid,
                                 sep,
                                 v.aa_grid,
                                 sep,
                                 v.pha_mat_data,
                                 sep,
                                 v.ext_mat_data,
                                 sep,
                                 v.abs_vec_data);

  return tags.bracket ? ("[" + out + "]") : out;
}

void xml_io_stream<SingleScatteringData>::write(std::ostream& os,
                                                const SingleScatteringData& x,
                                                bofstream* pbofs,
                                                std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, Index{x.ptype}, pbofs);
  xml_write_to_stream(os, x.description, pbofs);
  xml_write_to_stream(os, x.f_grid, pbofs);
  xml_write_to_stream(os, x.T_grid, pbofs);
  xml_write_to_stream(os, x.za_grid, pbofs);
  xml_write_to_stream(os, x.aa_grid, pbofs);
  xml_write_to_stream(os, x.pha_mat_data, pbofs);
  xml_write_to_stream(os, x.ext_mat_data, pbofs);
  xml_write_to_stream(os, x.abs_vec_data, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<SingleScatteringData>::read(std::istream& is,
                                               SingleScatteringData& x,
                                               bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  Index VAL;
  xml_read_from_stream(is, VAL, pbifs);
  xml_read_from_stream(is, x.description, pbifs);
  xml_read_from_stream(is, x.f_grid, pbifs);
  xml_read_from_stream(is, x.T_grid, pbifs);
  xml_read_from_stream(is, x.za_grid, pbifs);
  xml_read_from_stream(is, x.aa_grid, pbifs);
  xml_read_from_stream(is, x.pha_mat_data, pbifs);
  xml_read_from_stream(is, x.ext_mat_data, pbifs);
  xml_read_from_stream(is, x.abs_vec_data, pbifs);

  x.ptype = static_cast<PType>(VAL);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<ScatteringMetaData>::write(std::ostream& os,
                                              const ScatteringMetaData& x,
                                              bofstream* pbofs,
                                              std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.description, pbofs);
  xml_write_to_stream(os, x.source, pbofs);
  xml_write_to_stream(os, x.refr_index, pbofs);
  xml_write_to_stream(os, x.mass, pbofs);
  xml_write_to_stream(os, x.diameter_max, pbofs);
  xml_write_to_stream(os, x.diameter_volume_equ, pbofs);
  xml_write_to_stream(os, x.diameter_area_equ_aerodynamical, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<ScatteringMetaData>::read(std::istream& is,
                                             ScatteringMetaData& x,
                                             bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.description, pbifs);
  xml_read_from_stream(is, x.source, pbifs);
  xml_read_from_stream(is, x.refr_index, pbifs);
  xml_read_from_stream(is, x.mass, pbifs);
  xml_read_from_stream(is, x.diameter_max, pbifs);
  xml_read_from_stream(is, x.diameter_volume_equ, pbifs);
  xml_read_from_stream(is, x.diameter_area_equ_aerodynamical, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
