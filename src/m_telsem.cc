/*!
  \file   m_telsem.cc

  \brief  This file contains functions to read TELSEM atlases.
*/

#include "arts_conversions.h"
#include "file.h"
#include "geodetic.h"
#include "matpack_data.h"
#include "mystring.h"
#include "rte.h"
#include "telsem.h"
#include "check_input.h"

inline constexpr Numeric EARTH_RADIUS=Constant::earth_radius;
inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);

/* Workspace method: Doxygen documentation will be auto-generated */
void telsemStandalone(Matrix &emis,
                      const Numeric &lat,
                      const Numeric &lon,
                      const Numeric &theta,
                      const Vector &f,
                      const TelsemAtlas &atlas,
                      const Numeric &d_max,
                      const Verbosity &) {
  chk_if_in_range("Latitude input to TELSEM2", lat, -90.0, 90.0);
  chk_if_in_range("Longitude input to TELSEM2", lon, 0.0, 360.0);

  Index cellnumber = atlas.calc_cellnum(lat, lon);
  // Check if cell is in atlas.
  if (!atlas.contains(cellnumber)) {
    ARTS_USER_ERROR_IF (d_max <= 0.0,
          "Given coordinates are not contained in "
          " TELSEM atlas. To enable nearest neighbor"
          "interpolation set *d_max* to a positive "
          "value.");
    cellnumber = atlas.calc_cellnum_nearest_neighbor(lat, lon);
    Numeric lat_nn, lon_nn;
    std::tie(lat_nn, lon_nn) = atlas.get_coordinates(cellnumber);
    Numeric d = sphdist(lat, lon, lat_nn, lon_nn);
    ARTS_USER_ERROR_IF (d > d_max,
      "Distance of nearest neighbor exceeds provided limit (",
      d, " > ", d_max, ").")
  }

  Index class1 = atlas.get_class1(cellnumber);
  Index class2 = atlas.get_class2(cellnumber);
  Vector emis_v = atlas.get_emis_v(cellnumber);
  Vector emis_h = atlas.get_emis_h(cellnumber);

  emis.resize(f.nelem(), 2);
  for (Index i = 0; i < f.nelem(); ++i) {
    std::tie(emis(i, 0), emis(i, 1)) =
        atlas.emis_interp(theta, f[i] * 1e-9, class1, class2, emis_v, emis_h);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void telsemSurfaceTypeLandSea(Index &surface_type,
                              const Index &atmosphere_dim,
                              const Vector &lat_grid,
                              const Vector &lat_true,
                              const Vector &lon_true,
                              const Vector &rtp_pos,
                              const TelsemAtlas &atlas,
                              const Verbosity &) {
  // Checks
  chk_latlon_true(atmosphere_dim, lat_grid, lat_true, lon_true);

  Numeric lat, lon;
  pos2true_latlon(
      lat, lon, atmosphere_dim, lat_grid, lat_true, lon_true, rtp_pos);
  chk_if_in_range("Latitude input to TELSEM2", lat, -90.0, 90.0);
  chk_if_in_range("Longitude input to TELSEM2", lon, 0.0, 360.0);

  Index cellnumber = atlas.calc_cellnum(lat, lon);
  if (atlas.contains(cellnumber)) {
    surface_type = 1;
  } else {
    surface_type = 0;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void telsemAtlasLookup(Vector &emis,
                       const Numeric &lat,
                       const Numeric &lon,
                       const TelsemAtlas &atlas,
                       const Verbosity &) {
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
void telsem_atlasReadAscii(TelsemAtlas &atlas,
                           const String &directory,
                           const Index &month,
                           const String &filename_pattern,
                           const Verbosity &verbosity) {
  CREATE_OUT2;
  const Index imonth = filename_pattern.find("@MM@");
  ARTS_USER_ERROR_IF(imonth == String::npos,
                     "Substring '@MM@' not found in filename_pattern for\n",
                     "month number replacement: ",
                     filename_pattern)

  std::ifstream is;

  ostringstream month_ss;
  if (month < 10) {
    month_ss << 0;
  }
  month_ss << month;

  String this_filename = filename_pattern;
  this_filename.replace(imonth, 4, month_ss.str());
  this_filename = directory + '/' + this_filename;

  out2 << "Reading TELSEM atlas: " << this_filename << '\n';
  open_input_file(is, this_filename);
  atlas.read(is);
  atlas.set_month(month);

  String corr_filename = directory + '/' + "correlations";
  out2 << "Reading correlations: " << corr_filename << '\n';
  std::ifstream corr_is;
  open_input_file(corr_is, corr_filename);
  Tensor3 correlation(10, 7, 7);
  String s;
  for (Index i = 0; i < 10; i++) {
    std::getline(corr_is, s);
    for (Index j = 0; j < 7; j++) {
      for (Index k = 0; k < 7; k++) {
        corr_is >> double_imanip() >> correlation(i, j, k);
        ARTS_USER_ERROR_IF (corr_is.fail(),
                            "Error reading correlation.");
      }
      std::getline(corr_is, s);
    }
  }
  atlas.set_correl(correlation);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void telsem_atlasesReadAscii(ArrayOfTelsemAtlas &telsem_atlases,
                             const String &directory,
                             const String &filename_pattern,
                             const Verbosity &verbosity) {
  CREATE_OUT2;
  const Index imonth = filename_pattern.find("@MM@");
  ARTS_USER_ERROR_IF(imonth == String::npos,
                     "Substring '@MM@' not found in filename_pattern for\n",
                     "month number replacement: ",
                     filename_pattern)

  telsem_atlases.resize(12);
  for (Index i = 1; i <= 12; i++) {
    std::ifstream is;
    ostringstream month;
    if (i < 10) month << 0;
    month << i;
    String this_filename = filename_pattern;
    this_filename.replace(imonth, 4, month.str());
    this_filename = directory + '/' + this_filename;

    out2 << "Reading TELSEM atlas: " << this_filename << '\n';
    open_input_file(is, this_filename);
    telsem_atlases[i - 1].read(is);
    telsem_atlases[i - 1].set_month(i);
  }

  std::ifstream is;
  String corr_filename = directory + '/' + "correlations";
  out2 << "Reading correlations: " << corr_filename << '\n';
  open_input_file(is, corr_filename);
  Tensor3 correlation(10, 7, 7);
  String s;
  for (Index i = 0; i < 10; i++) {
    std::getline(is, s);
    for (Index j = 0; j < 7; j++) {
      for (Index k = 0; k < 7; k++) {
        is >> double_imanip() >> correlation(i, j, k);
        ARTS_USER_ERROR_IF (is.fail(), "Error reading correlation.");
      }
      std::getline(is, s);
    }
  }

  for (Index i = 0; i < 12; i++) {
    telsem_atlases[i].set_correl(correlation);
  }
}
