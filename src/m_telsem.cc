/*!
  \file   m_telsem.cc

  \brief  This file contains functions to read TELSEM atlases.
*/

#include "arts_conversions.h"
#include "check_input.h"
#include "double_imanip.h"
#include "file.h"
#include "matpack_data.h"
#include "mystring.h"
#include "telsem.h"

inline constexpr Numeric EARTH_RADIUS = Constant::earth_radius;
inline constexpr Numeric DEG2RAD = Conversion::deg2rad(1);

/* Workspace method: Doxygen documentation will be auto-generated */
void emissivitiesTelsemAtlasLookup(Vector &emis,
                                   const Numeric &lat,
                                   const Numeric &lon,
                                   const TelsemAtlas &atlas) {
  chk_if_in_range("Latitude input to TELSEM2", lat, -90.0, 90.0);
  chk_if_in_range("Longitude input to TELSEM2", lon, 0.0, 360.0);

  Index cellnumber = atlas.calc_cellnum(lat, lon);
  if (atlas.contains(cellnumber)) {
    emis = atlas[cellnumber];
  } else {
    emis.resize(0);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void atlasReadAscii(TelsemAtlas &atlas,
                    const String &directory,
                    const Index &month,
                    const String &filename_pattern) {
  const Size imonth = filename_pattern.find("@MM@");
  ARTS_USER_ERROR_IF(imonth == String::npos,
                     "Substring '@MM@' not found in filename_pattern for\n"
                     "month number replacement: {}",
                     filename_pattern)

  std::ifstream is;

  std::ostringstream month_ss;
  if (month < 10) {
    month_ss << 0;
  }
  month_ss << month;

  String this_filename = filename_pattern;
  this_filename.replace(imonth, 4, month_ss.str());
  this_filename = directory + '/' + this_filename;

  open_input_file(is, this_filename);
  atlas.read(is);
  atlas.set_month(month);

  String corr_filename = directory + '/' + "correlations";
  std::ifstream corr_is;
  open_input_file(corr_is, corr_filename);
  Tensor3 correlation(10, 7, 7);
  String s;
  for (Index i = 0; i < 10; i++) {
    std::getline(corr_is, s);
    for (Index j = 0; j < 7; j++) {
      for (Index k = 0; k < 7; k++) {
        corr_is >> double_imanip() >> correlation(i, j, k);
        ARTS_USER_ERROR_IF(corr_is.fail(), "Error reading correlation.");
      }
      std::getline(corr_is, s);
    }
  }
  atlas.set_correl(correlation);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void telsem_atlasesReadAscii(ArrayOfTelsemAtlas &telsem_atlases,
                             const String &directory,
                             const String &filename_pattern) {
  const Size imonth = filename_pattern.find("@MM@");
  ARTS_USER_ERROR_IF(imonth == String::npos,
                     "Substring '@MM@' not found in filename_pattern for\n"
                     "month number replacement: {}",
                     filename_pattern)

  telsem_atlases.resize(12);
  for (Index i = 1; i <= 12; i++) {
    std::ifstream is;
    std::ostringstream month;
    if (i < 10) month << 0;
    month << i;
    String this_filename = filename_pattern;
    this_filename.replace(imonth, 4, month.str());
    this_filename = directory + '/' + this_filename;

    open_input_file(is, this_filename);
    telsem_atlases[i - 1].read(is);
    telsem_atlases[i - 1].set_month(i);
  }

  std::ifstream is;
  String corr_filename = directory + '/' + "correlations";
  open_input_file(is, corr_filename);
  Tensor3 correlation(10, 7, 7);
  String s;
  for (Index i = 0; i < 10; i++) {
    std::getline(is, s);
    for (Index j = 0; j < 7; j++) {
      for (Index k = 0; k < 7; k++) {
        is >> double_imanip() >> correlation(i, j, k);
        ARTS_USER_ERROR_IF(is.fail(), "Error reading correlation.");
      }
      std::getline(is, s);
    }
  }

  for (Index i = 0; i < 12; i++) {
    telsem_atlases[i].set_correl(correlation);
  }
}
