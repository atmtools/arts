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

void chk_scat_data(const SingleScatteringData& scat_data_single) {
  assert(scat_data_single.ptype == PTYPE_GENERAL ||
         scat_data_single.ptype == PTYPE_TOTAL_RND ||
         scat_data_single.ptype == PTYPE_AZIMUTH_RND);

  ARTS_USER_ERROR_IF(scat_data_single.za_grid[0] != 0.,
                     "The first value of the zenith angle grid in the single"
                     " scattering properties data must be 0.")

  ARTS_USER_ERROR_IF(last(scat_data_single.za_grid) != 180.,
                     "The last value of the zenith angle grid in the single"
                     " scattering properties data must be 180.")

  ARTS_USER_ERROR_IF(scat_data_single.ptype == PTYPE_GENERAL &&
                         scat_data_single.aa_grid[0] != -180.,
                     "For ptype = \"general\" the first value"
                     " of the azimuth angle grid in the single scattering"
                     " properties data must be -180.")

  ARTS_USER_ERROR_IF(scat_data_single.ptype == PTYPE_AZIMUTH_RND &&
                         scat_data_single.aa_grid[0] != 0.,
                     "For ptype = \"azimuthally_random\""
                     " the first value"
                     " of the azimuth angle grid in the single scattering"
                     " properties data must be 0.")

  ARTS_USER_ERROR_IF(scat_data_single.ptype != PTYPE_TOTAL_RND &&
                         last(scat_data_single.aa_grid) != 180.,
                     "For ptypes = \"azimuthally_random\" and \"general\""
                     " the last value of the azimuth angle grid in the single"
                     " scattering properties data must be 180.")

  std::ostringstream os_pha_mat;
  os_pha_mat << "pha_mat ";
  std::ostringstream os_ext_mat;
  os_ext_mat << "ext_mat ";
  std::ostringstream os_abs_vec;
  os_abs_vec << "abs_vec ";

  switch (scat_data_single.ptype) {
    case PTYPE_GENERAL:

      ARTS_USER_ERROR_IF(not same_shape<7>({scat_data_single.f_grid.size(),
                                            scat_data_single.T_grid.size(),
                                            scat_data_single.za_grid.size(),
                                            scat_data_single.aa_grid.size(),
                                            scat_data_single.za_grid.size(),
                                            scat_data_single.aa_grid.size(),
                                            16},
                                           scat_data_single.pha_mat_data),
                         R"(The shape of the pha_mat_data is wrong: {}
      
  Expected size: [{}, {}, {}, {}, {}, {}, {}]
  Got shape:     {:B,}
)",
                         os_pha_mat.str(),
                         scat_data_single.f_grid.size(),
                         scat_data_single.T_grid.size(),
                         scat_data_single.za_grid.size(),
                         scat_data_single.aa_grid.size(),
                         scat_data_single.za_grid.size(),
                         scat_data_single.aa_grid.size(),
                         16,
                         scat_data_single.pha_mat_data.shape());

      ARTS_USER_ERROR_IF(not same_shape<5>({scat_data_single.f_grid.size(),
                                            scat_data_single.T_grid.size(),
                                            scat_data_single.za_grid.size(),
                                            scat_data_single.aa_grid.size(),
                                            7},
                                           scat_data_single.ext_mat_data),
                         R"(The shape of the ext_mat_data is wrong: {}
      
  Expected size: [{}, {}, {}, {}, {}]
  Got shape:     {:B,}
)",
                         os_ext_mat.str(),
                         scat_data_single.f_grid.size(),
                         scat_data_single.T_grid.size(),
                         scat_data_single.za_grid.size(),
                         scat_data_single.aa_grid.size(),
                         7,
                         scat_data_single.ext_mat_data.shape());

      ARTS_USER_ERROR_IF(not same_shape<5>({scat_data_single.f_grid.size(),
                                            scat_data_single.T_grid.size(),
                                            scat_data_single.za_grid.size(),
                                            scat_data_single.aa_grid.size(),
                                            4},
                                           scat_data_single.abs_vec_data),
                         R"(The shape of the abs_vec_data is wrong: {}
      
  Expected size: [{}, {}, {}, {}, {}]
  Got shape:     {:B,}
)",
                         os_abs_vec.str(),
                         scat_data_single.f_grid.size(),
                         scat_data_single.T_grid.size(),
                         scat_data_single.za_grid.size(),
                         scat_data_single.aa_grid.size(),
                         4,
                         scat_data_single.abs_vec_data.shape());

      break;

    case PTYPE_TOTAL_RND:

      ARTS_USER_ERROR_IF(not same_shape<7>({scat_data_single.f_grid.size(),
                                            scat_data_single.T_grid.size(),
                                            scat_data_single.za_grid.size(),
                                            1,
                                            1,
                                            1,
                                            6},
                                           scat_data_single.pha_mat_data),
                         R"(The shape of the pha_mat data is wrong: {}
      
  Expected size: [{}, {}, {}, {}, {}, {}, {}]
  Got shape:     {:B,}
)",
                         os_pha_mat.str(),
                         scat_data_single.f_grid.size(),
                         scat_data_single.T_grid.size(),
                         scat_data_single.za_grid.size(),
                         1,
                         1,
                         1,
                         6,
                         scat_data_single.pha_mat_data.shape());

      ARTS_USER_ERROR_IF(not same_shape<5>({scat_data_single.f_grid.size(),
                                            scat_data_single.T_grid.size(),
                                            1,
                                            1,
                                            1},
                                           scat_data_single.ext_mat_data),
                         R"(The shape of the ext_mat_data is wrong: {}
      
  Expected size: [{}, {}, {}, {}, {}]
  Got shape:     {:B,}
)",
                         os_ext_mat.str(),
                         scat_data_single.f_grid.size(),
                         scat_data_single.T_grid.size(),
                         1,
                         1,
                         1,
                         scat_data_single.ext_mat_data.shape());

      ARTS_USER_ERROR_IF(not same_shape<5>({scat_data_single.f_grid.size(),
                                            scat_data_single.T_grid.size(),
                                            1,
                                            1,
                                            1},
                                           scat_data_single.abs_vec_data),
                         R"(The shape of the abs_vec_data is wrong: {}
      
  Expected size: [{}, {}, {}, {}, {}]
  Got shape:     {:B,}
)",
                         os_abs_vec.str(),
                         scat_data_single.f_grid.size(),
                         scat_data_single.T_grid.size(),
                         1,
                         1,
                         1,
                         scat_data_single.abs_vec_data.shape());

      break;

    case PTYPE_AZIMUTH_RND:

      ARTS_USER_ERROR_IF(not same_shape<7>({scat_data_single.f_grid.size(),
                                            scat_data_single.T_grid.size(),
                                            scat_data_single.za_grid.size(),
                                            scat_data_single.aa_grid.size(),
                                            scat_data_single.za_grid.size(),
                                            1,
                                            16},
                                           scat_data_single.pha_mat_data),
                         R"(The shape of the pha_mat_data is wrong: {}
      
  Expected size: [{}, {}, {}, {}, {}, {}, {}]
  Got shape:     {:B,}
)",
                         os_pha_mat.str(),
                         scat_data_single.f_grid.size(),
                         scat_data_single.T_grid.size(),
                         scat_data_single.za_grid.size(),
                         scat_data_single.aa_grid.size(),
                         scat_data_single.za_grid.size(),
                         1,
                         16,
                         scat_data_single.pha_mat_data.shape());

      ARTS_USER_ERROR_IF(not same_shape<5>({scat_data_single.f_grid.size(),
                                            scat_data_single.T_grid.size(),
                                            scat_data_single.za_grid.size(),
                                            1,
                                            3},
                                           scat_data_single.ext_mat_data),
                         R"(The shape of the ext_mat_data data is wrong: {}
      
  Expected size: [{}, {}, {}, {}, {}]
  Got shape:     {:B,}
)",
                         os_ext_mat.str(),
                         scat_data_single.f_grid.size(),
                         scat_data_single.T_grid.size(),
                         scat_data_single.za_grid.size(),
                         1,
                         3,
                         scat_data_single.ext_mat_data.shape());

      ARTS_USER_ERROR_IF(not same_shape<5>({scat_data_single.f_grid.size(),
                                            scat_data_single.T_grid.size(),
                                            scat_data_single.za_grid.size(),
                                            1,
                                            2},
                                           scat_data_single.abs_vec_data),
                         R"(The shape of the abs_vec_data is wrong: {}
      
  Expected size: [{}, {}, {}, {}, {}]
  Got shape:     {:B,}
)",
                         os_abs_vec.str(),
                         scat_data_single.f_grid.size(),
                         scat_data_single.T_grid.size(),
                         scat_data_single.za_grid.size(),
                         1,
                         2,
                         scat_data_single.abs_vec_data.shape());

      break;
  }

  // Here we only check whether the temperature grid is of the unit K, not
  // whether it corresponds to the required values in t_field. The second
  // option is not trivial since here one has to look whether the pnd_field
  // is non-zero for the corresponding temperature. This check is done in the
  // functions where the multiplication with the particle number density is
  // done.
  ARTS_USER_ERROR_IF(
      scat_data_single.T_grid[0] < 0. || last(scat_data_single.T_grid) > 1001.,
      "The temperature values in the single scattering data"
      " are negative or very large. Check whether you use the "
      "right unit [Kelvin].")
}

void ConvertAzimuthallyRandomSingleScatteringData(SingleScatteringData& ssd) {
  // First check that input fulfills requirements on older data formats:
  // 1) Is za_grid symmetric and includes 90deg?
  Index nza = ssd.za_grid.size();
  for (Index i = 0; i < nza / 2; i++) {
    ARTS_USER_ERROR_IF(
        !is_same_within_epsilon(
            180. - ssd.za_grid[nza - 1 - i], ssd.za_grid[i], 2 * std::numeric_limits<double>::epsilon()),
        "Zenith grid of azimuthally_random single scattering data\n"
        "is not symmetric with respect to 90degree.")
  }
  ARTS_USER_ERROR_IF(
      !is_same_within_epsilon(ssd.za_grid[nza / 2], 90., 2 * std::numeric_limits<double>::epsilon()),
      "Zenith grid of azimuthally_random single scattering data\n"
      "does not contain 90 degree grid point.")

  // 2) Are data sizes correct?
  std::ostringstream os_pha_mat;
  os_pha_mat << "pha_mat ";
  std::ostringstream os_ext_mat;
  os_ext_mat << "ext_mat ";
  std::ostringstream os_abs_vec;
  os_abs_vec << "abs_vec ";

  ARTS_USER_ERROR_IF(
      not same_shape<7>({static_cast<Index>(ssd.f_grid.size()),
                         static_cast<Index>(ssd.T_grid.size()),
                         static_cast<Index>(ssd.za_grid.size()),
                         static_cast<Index>(ssd.aa_grid.size()),
                         static_cast<Index>(ssd.za_grid.size()) / 2 + 1,
                         1,
                         16},
                        ssd.pha_mat_data),
      "Error in {0}.\n\tGrid shape [{1}, {2}, {3}, {4}, {5}, {6}, {7}] versus data shape {8:B,}",
      os_pha_mat.str(),
      ssd.f_grid.size(),
      ssd.T_grid.size(),
      ssd.za_grid.size(),
      ssd.aa_grid.size(),
      ssd.za_grid.size() / 2 + 1,
      1,
      16,
      ssd.pha_mat_data.shape());

  ARTS_USER_ERROR_IF(
      not same_shape<5>({static_cast<Index>(ssd.f_grid.size()),
                         static_cast<Index>(ssd.T_grid.size()),
                         static_cast<Index>(ssd.za_grid.size()) / 2 + 1,
                         1,
                         3},
                        ssd.ext_mat_data),
      "Error in {0}.\n\tGrid shape [{1}, {2}, {3}, {4}, {5}] versus data shape {6:B,}",
      os_ext_mat.str(),
      static_cast<Index>(ssd.f_grid.size()),
      static_cast<Index>(ssd.T_grid.size()),
      static_cast<Index>(ssd.za_grid.size()) / 2 + 1,
      1,
      3,
      ssd.ext_mat_data.shape());

  ARTS_USER_ERROR_IF(
      not same_shape<5>({static_cast<Index>(ssd.f_grid.size()),
                         static_cast<Index>(ssd.T_grid.size()),
                         static_cast<Index>(ssd.za_grid.size()) / 2 + 1,
                         1,
                         2},
                        ssd.abs_vec_data),
      "Error in {0}.\n\tGrid shape [{1}, {2}, {3}, {4}, {5}] versus data shape {6:B,}",
      os_abs_vec.str(),
      static_cast<Index>(ssd.f_grid.size()),
      static_cast<Index>(ssd.T_grid.size()),
      static_cast<Index>(ssd.za_grid.size()) / 2 + 1,
      1,
      2,
      ssd.abs_vec_data.shape());

  // Now that we are sure that za_grid is properly symmetric, we just need to
  // copy over the data (ie no interpolation).
  Tensor5 tmpT5 = ssd.abs_vec_data;
  ssd.abs_vec_data.resize(tmpT5.nshelves(),
                          tmpT5.nbooks(),
                          ssd.za_grid.size(),
                          tmpT5.nrows(),
                          tmpT5.ncols());
  ssd.abs_vec_data[joker, joker, Range(0, nza / 2 + 1), joker, joker] = tmpT5;
  for (Index i = 0; i < nza / 2; i++) {
    ssd.abs_vec_data[joker, joker, nza - 1 - i, joker, joker] =
        tmpT5[joker, joker, i, joker, joker];
  }

  tmpT5 = ssd.ext_mat_data;
  ssd.ext_mat_data.resize(tmpT5.nshelves(),
                          tmpT5.nbooks(),
                          ssd.za_grid.size(),
                          tmpT5.nrows(),
                          tmpT5.ncols());
  ssd.ext_mat_data[joker, joker, Range(0, nza / 2 + 1), joker, joker] = tmpT5;
  for (Index i = 0; i < nza / 2; i++) {
    ssd.ext_mat_data[joker, joker, nza - 1 - i, joker, joker] =
        tmpT5[joker, joker, i, joker, joker];
  }

  Tensor7 tmpT7 = ssd.pha_mat_data;
  ssd.pha_mat_data.resize(tmpT7.nlibraries(),
                          tmpT7.nvitrines(),
                          tmpT7.nshelves(),
                          tmpT7.nbooks(),
                          ssd.za_grid.size(),
                          tmpT7.nrows(),
                          tmpT7.ncols());
  ssd.pha_mat_data
      [joker, joker, joker, joker, Range(0, nza / 2 + 1), joker, joker] = tmpT7;

  // scatt. matrix elements 13,23,31,32 and 14,24,41,42 (=elements 2,6,8,9 and
  // 3,7,12,13 in ARTS' flattened format, respectively) change sign.
  tmpT7[joker, joker, joker, joker, joker, joker, Range(2, 2)]  *= -1.;
  tmpT7[joker, joker, joker, joker, joker, joker, Range(6, 4)]  *= -1.;
  tmpT7[joker, joker, joker, joker, joker, joker, Range(12, 2)] *= -1.;

  // For second half of incident polar angles (>90deg), we need to mirror the
  // original data in both incident and scattered polar angle around 90deg "planes".
  for (Index i = 0; i < nza / 2; i++)
    for (Index j = 0; j < nza; j++)
      ssd.pha_mat_data
          [joker, joker, nza - 1 - j, joker, nza - 1 - i, joker, joker] =
          tmpT7[joker, joker, j, joker, i, joker, joker];
}

void xml_io_stream<SingleScatteringData>::write(
    std::ostream& os_xml,
    const SingleScatteringData& ssdata,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag open_tag;
  XMLTag close_tag;

  open_tag.name=("SingleScatteringData");
  if (name.length()) open_tag.add_attribute("name", name);
  open_tag.add_attribute("version", "3");
  open_tag.write_to_stream(os_xml);

  std::println(os_xml);
  xml_write_to_stream(os_xml, String{PTypeToString(ssdata.ptype)}, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.description, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.f_grid, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.T_grid, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.za_grid, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.aa_grid, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.pha_mat_data, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.ext_mat_data, pbofs, "");
  xml_write_to_stream(os_xml, ssdata.abs_vec_data, pbofs, "");

  close_tag.name=("/SingleScatteringData");
  close_tag.write_to_stream(os_xml);
  std::println(os_xml);
}

void xml_io_stream<SingleScatteringData>::read(std::istream& is_xml,
                                               SingleScatteringData& ssdata,
                                               bifstream* pbifs) try {
  XMLTag tag;
  String version;

  tag.read_from_stream(is_xml);
  tag.check_name("SingleScatteringData");
  tag.get_attribute_value("version", version);

  if (version == "3") {
    String ptype_string;
    xml_read_from_stream(is_xml, ptype_string, pbifs);
    ssdata.ptype = PTypeFromString(ptype_string);
  } else if (version == "2") {
    String ptype_string;
    xml_read_from_stream(is_xml, ptype_string, pbifs);
    ssdata.ptype = PType2FromString(ptype_string);
  } else {
    Index ptype;
    xml_read_from_stream(is_xml, ptype, pbifs);
    if (ptype != PTYPE_GENERAL && ptype != PTYPE_TOTAL_RND &&
        ptype != PTYPE_AZIMUTH_RND) {
      std::ostringstream os;
      os << "Ptype value (" << ptype << ") is wrong."
         << "It must be \n"
         << PTYPE_TOTAL_RND << " - totally randomly oriented particles,\n"
         << PTYPE_AZIMUTH_RND
         << " - azimuthally randomly oriented particles, or\n"
         << PTYPE_GENERAL << " - arbitrary oriented particles.\n";
      throw std::runtime_error(os.str());
    }
    ssdata.ptype = PType(ptype);
  }
  xml_read_from_stream(is_xml, ssdata.description, pbifs);
  xml_read_from_stream(is_xml, ssdata.f_grid, pbifs);
  xml_read_from_stream(is_xml, ssdata.T_grid, pbifs);
  xml_read_from_stream(is_xml, ssdata.za_grid, pbifs);
  /* Verify that we have a good coverage for the za grid */
  if ((ssdata.za_grid[0] > 1) ||
      ssdata.za_grid[ssdata.za_grid.size() - 1] < 179) {
    std::ostringstream os;
    os << "Missing data in xml-stream. Expected za_grid: [0, 180]. "
       << "Found za_grid: [" << ssdata.za_grid[0] << ", "
       << ssdata.za_grid[ssdata.za_grid.size() - 1] << "]";
    throw std::runtime_error(os.str());
  }
  xml_read_from_stream(is_xml, ssdata.aa_grid, pbifs);

  xml_read_from_stream(is_xml, ssdata.pha_mat_data, pbifs);
  if (static_cast<Size>(ssdata.pha_mat_data.nlibraries()) !=
      ssdata.f_grid.size()) {
    throw std::runtime_error(
        "Number of frequencies in f_grid and pha_mat_data "
        "not matching!!!");
  }

  xml_read_from_stream(is_xml, ssdata.ext_mat_data, pbifs);
  xml_read_from_stream(is_xml, ssdata.abs_vec_data, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_name("/SingleScatteringData");

  if (version != "3" && ssdata.ptype == PTYPE_AZIMUTH_RND) {
    ConvertAzimuthallyRandomSingleScatteringData(ssdata);
  }

  chk_scat_data(ssdata);
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error reading {}:\n{}", type_name.data(), e.what()));
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
                                             bifstream* pbifs) try {
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
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error reading {}:\n{}", type_name.data(), e.what()));
}
