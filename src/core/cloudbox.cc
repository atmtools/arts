/*!
  \file   cloudbox.cc
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   Thu May  23 10:59:55 2002
  
  \brief  Internal functions for scattering calculations.
*/

#include "cloudbox.h"

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <algorithm>
#include <cmath>

#include "check_input.h"
#include "math_funcs.h"

//! Check particle number density files
/*!
  This function checks, whether the particle number density file
  has the right atmospheric dimension (check for no non-zero pnd values outside
  cloudbox removed. done in pnd_fieldCalcFrompnd_field_raw). 

  \param pnd_field_raw   pnd field data
  \param pnd_field_file  pnd field filename

  \author Claudia Emde
  \date   2005-04-05
*/
void chk_pnd_data(const GriddedField3& pnd_field_raw,
                  const String& pnd_field_file) {
  const Vector& pfr_lat_grid =
      pnd_field_raw.grid<1>();
  const Vector& pfr_lon_grid =
      pnd_field_raw.grid<2>();

  // The consistency of the dimensions is checked in the reading routine.
  // Here we have to check whether the atmospheric dimension is correct and whether
  // the particle number density is 0 on the cloudbox boundary and outside the cloudbox.

    ARTS_USER_ERROR_IF (pfr_lat_grid.size() == 1 || pfr_lon_grid.size() == 1,
         "The atmospheric dimension is 3D but the particle "
         "number density file * ", pnd_field_file,
         " is for a 1D or a 2D atmosphere. \n")
}

//! Check particle number density files (pnd_field_raw)
/*!
  
  \param pnd_field_raw   pnd field raw data (array for all scattering elements)
  \param pnd_field_file  pnd field filename
 
  \author Claudia Emde
  \date   2005-04-05
*/
void chk_pnd_raw_data(const ArrayOfGriddedField3& pnd_field_raw,
                      const String& pnd_field_file) {
  for (Size i = 0; i < pnd_field_raw.size(); i++) {
    chk_pnd_data(pnd_field_raw[i], pnd_field_file);
  }
}

//! chk_pnd_field_raw_only_in_cloudbox
/*! 
    Checks whether the pnd_field is zero outside the cloudbox.
    This is of a higher level than chk_pnd_data because it does
    not require any filename and because it works on all pnd_field_raw
    rather than just one element. Otherwise, it is mostly a new
    implementation of the same functionality.

    \param    dim                The atmospheric dimensionality.
    \param    pnd_field_raw      All pnd_field_raw data.
    \param    p_grid             Pressure grid.
    \param    lat_grid           Latitude grid.
    \param    lon_grid           Longitude grid.
    \param    cloudbox_limits    The edges of the cloudbox.

    \author Gerrit Holl
    \date   2011-03-24
*/
void chk_pnd_field_raw_only_in_cloudbox(
    const Index& dim,
    const ArrayOfGriddedField3& pnd_field_raw,
    ConstVectorView p_grid,
    ConstVectorView lat_grid,
    ConstVectorView lon_grid,
    const ArrayOfIndex& cloudbox_limits) {
  Numeric p, lat, lon, v;
  Index n, p_i, lat_i, lon_i;
  // For any non-zero point, verify we're outside the cloudbox
  for (n = 0; n < static_cast<Index>(pnd_field_raw.size()); n++) {
    for (p_i = 0; p_i < pnd_field_raw[n].data.npages(); p_i++) {
      for (lat_i = 0; lat_i < pnd_field_raw[n].data.nrows(); lat_i++) {
        for (lon_i = 0; lon_i < pnd_field_raw[n].data.ncols(); lon_i++) {
          v = pnd_field_raw[n].data(p_i, lat_i, lon_i);
          if (v != 0) {
            // Verify pressure is between cloudbox limits
            p = pnd_field_raw[n].grid<0>()[p_i];
            //                        if (!((p <= p_grid[cloudbox_limits[0]]) &
            //                              (p >= p_grid[cloudbox_limits[1]]))) {
            ARTS_USER_ERROR_IF ((p <= p_grid[cloudbox_limits[1]]) ||
                ((p >= p_grid[cloudbox_limits[0]]) &&
                 (cloudbox_limits[0] != 0)),
                 "Found non-zero pnd outside cloudbox. "
                 "Cloudbox extends from p=", p_grid[cloudbox_limits[0]],
                 " Pa to p=", p_grid[cloudbox_limits[1]],
                 " Pa, but found pnd=", v, "/m³ at p=", p,
                 " Pa for scattering "
                 "element #", n, ".")
            // Verify latitude is too
            if (dim > 1) {
              lat = pnd_field_raw[n].grid<1>()[lat_i];
              ARTS_USER_ERROR_IF (!((lat > lat_grid[cloudbox_limits[2]]) and
                    (lat < lat_grid[cloudbox_limits[3]])),
                  "Found non-zero pnd outside cloudbox. "
                  "Cloudbox extends from lat=",
                  lat_grid[cloudbox_limits[2]],
                  "° to lat=", lat_grid[cloudbox_limits[3]],
                  "°, but found pnd=", v, "/m³ at lat=", lat,
                  "° for scattering "
                  "element #", n, ".")
            }
            // Etc. for longitude
            if (dim > 2) {
              lon = pnd_field_raw[n].grid<2>()[lon_i];
              ARTS_USER_ERROR_IF (!((lon > lon_grid[cloudbox_limits[4]]) and
                    (lon < lon_grid[cloudbox_limits[5]])),
                  "Found non-zero pnd outside cloudbox. "
                  "Cloudbox extends from lon=",
                  lon_grid[cloudbox_limits[4]],
                  "° to lat=", lon_grid[cloudbox_limits[5]],
                  "°, but found pnd=", v, "/m³ at lon=", lon,
                  "° for scattering "
                  "element #", n, ".")
            }
          }
        }
      }
    }
  }
}

//!  Check validity of scat_species setting
/*!
  This function checks, whether number of elements in each scattering species
  string is ok, and whether the entries for size limits are indeed numbers (or
  '*').

	\param scat_species Array of scattering species tags.
  \param delim        delimiter string of *scat_species* elements.

  \author Jana Mendrok
  \date 2012-10-25

*/
void chk_scat_species(const ArrayOfString& scat_species, const String& delim) {
  ArrayOfString strarr;
  Size size = 2;

  for (Size k = 0; k < scat_species.size(); k++) {
    split(strarr, scat_species[k], delim);
    ARTS_USER_ERROR_IF (strarr.size() < size,
         "Individual strings in scat_species must contain at least ", size,
         " elements,\n"
         "but entry #", k, " contains only the following ",
         strarr.size(), ":\n",
         strarr, "\n")
  }
}

//! Check scattering data general
/*!
  FIXME
  
  \param scat_data Array of single scattering data
  \param scat_meta Array of scattering meta data

  \author Daniel Kreyling
  \date 2010-12-02
*/

void chk_scattering_data(const ArrayOfSingleScatteringData& scat_data,
                         const ArrayOfScatteringMetaData& scat_meta) {
  ARTS_USER_ERROR_IF (scat_data.size() != scat_meta.size(),
      "The number of elements in in current scat_species'  *scat_data* and "
      "*scat_meta* do not match.\n"
      "Each *scat_data* entry must correspond to one entry in *scat_meta*.")
}

//! Check single scattering data
/*!
  This function checks the self consistency of the data by checking the
  dimensions of pha_mat, ext_mat and abs_vec depending on the ptype case.
  It furthermore checks whether the angular grids are defined correctly
  depending on ptype and the sanity of the temperature grid.
  
  \param scat_data_single[in]  Single scattering data of a single scattering element

  \author Claudia Emde
  \date   2005-04-04
*/
void chk_scat_data(const SingleScatteringData& scat_data_single) {
  ARTS_ASSERT(scat_data_single.ptype == PTYPE_GENERAL ||
         scat_data_single.ptype == PTYPE_TOTAL_RND ||
         scat_data_single.ptype == PTYPE_AZIMUTH_RND);

  ARTS_USER_ERROR_IF (scat_data_single.za_grid[0] != 0.,
      "The first value of the zenith angle grid in the single"
      " scattering properties data must be 0.")

  ARTS_USER_ERROR_IF (last(scat_data_single.za_grid) != 180.,
      "The last value of the zenith angle grid in the single"
      " scattering properties data must be 180.")

  ARTS_USER_ERROR_IF (scat_data_single.ptype == PTYPE_GENERAL &&
                      scat_data_single.aa_grid[0] != -180.,
      "For ptype = \"general\" the first value"
      " of the azimuth angle grid in the single scattering"
      " properties data must be -180.")

  ARTS_USER_ERROR_IF (scat_data_single.ptype == PTYPE_AZIMUTH_RND &&
                      scat_data_single.aa_grid[0] != 0.,
      "For ptype = \"azimuthally_random\""
      " the first value"
      " of the azimuth angle grid in the single scattering"
      " properties data must be 0.")

  ARTS_USER_ERROR_IF (scat_data_single.ptype != PTYPE_TOTAL_RND &&
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

      chk_size(os_pha_mat.str(),
               scat_data_single.pha_mat_data,
               scat_data_single.f_grid.size(),
               scat_data_single.T_grid.size(),
               scat_data_single.za_grid.size(),
               scat_data_single.aa_grid.size(),
               scat_data_single.za_grid.size(),
               scat_data_single.aa_grid.size(),
               16);

      chk_size(os_ext_mat.str(),
               scat_data_single.ext_mat_data,
               scat_data_single.f_grid.size(),
               scat_data_single.T_grid.size(),
               scat_data_single.za_grid.size(),
               scat_data_single.aa_grid.size(),
               7);

      chk_size(os_abs_vec.str(),
               scat_data_single.abs_vec_data,
               scat_data_single.f_grid.size(),
               scat_data_single.T_grid.size(),
               scat_data_single.za_grid.size(),
               scat_data_single.aa_grid.size(),
               4);
      break;

    case PTYPE_TOTAL_RND:

      chk_size(os_pha_mat.str(),
               scat_data_single.pha_mat_data,
               scat_data_single.f_grid.size(),
               scat_data_single.T_grid.size(),
               scat_data_single.za_grid.size(),
               1,
               1,
               1,
               6);

      chk_size(os_ext_mat.str(),
               scat_data_single.ext_mat_data,
               scat_data_single.f_grid.size(),
               scat_data_single.T_grid.size(),
               1,
               1,
               1);

      chk_size(os_abs_vec.str(),
               scat_data_single.abs_vec_data,
               scat_data_single.f_grid.size(),
               scat_data_single.T_grid.size(),
               1,
               1,
               1);
      break;

    case PTYPE_AZIMUTH_RND:

      chk_size(os_pha_mat.str(),
               scat_data_single.pha_mat_data,
               scat_data_single.f_grid.size(),
               scat_data_single.T_grid.size(),
               scat_data_single.za_grid.size(),
               scat_data_single.aa_grid.size(),
               scat_data_single.za_grid.size(),
               1,
               16);

      chk_size(os_ext_mat.str(),
               scat_data_single.ext_mat_data,
               scat_data_single.f_grid.size(),
               scat_data_single.T_grid.size(),
               scat_data_single.za_grid.size(),
               1,
               3);

      chk_size(os_abs_vec.str(),
               scat_data_single.abs_vec_data,
               scat_data_single.f_grid.size(),
               scat_data_single.T_grid.size(),
               scat_data_single.za_grid.size(),
               1,
               2);
      break;
  }

  // Here we only check whether the temperature grid is of the unit K, not
  // whether it corresponds to the required values in t_field. The second
  // option is not trivial since here one has to look whether the pnd_field
  // is non-zero for the corresponding temperature. This check is done in the
  // functions where the multiplication with the particle number density is
  // done.
  ARTS_USER_ERROR_IF (scat_data_single.T_grid[0] < 0. ||
                      last(scat_data_single.T_grid) > 1001.,
        "The temperature values in the single scattering data"
        " are negative or very large. Check whether you use the "
        "right unit [Kelvin].")
}

/*! Checks, whether a gridpoint is inside the cloudbox.

    \return true is returned if the point is inside the 
          cloudbox.
          
  \param gp_p  pressure GridPos
  \param gp_lat latitude GridPos
  \param gp_lon longitude GridPos
  \param cloudbox_limits The limits of the cloudbox.
  \param include_boundaries boolean: determines whther or not points on the 
  boundary are considered to be inside the cloudbox.

  \author Claudia Emde (rewritten by Cory Davis 2005-07-03)
  \date 2003-06-06

*/
bool is_gp_inside_cloudbox(const GridPos& gp_p,
                           const GridPos& gp_lat,
                           const GridPos& gp_lon,
                           const ArrayOfIndex& cloudbox_limits,
                           const bool& include_boundaries)

{
  if (include_boundaries) {
    // Pressure dimension
    double ipos = fractional_gp(gp_p);
    if (ipos < double(cloudbox_limits[0]) ||
        ipos > double(cloudbox_limits[1])) {
      return false;
    }

      // Latitude dimension
      ipos = fractional_gp(gp_lat);
      if (ipos < double(cloudbox_limits[2]) ||
          ipos > double(cloudbox_limits[3])) {
        return false;
      }

        // Longitude dimension
        ipos = fractional_gp(gp_lon);
        if (ipos < double(cloudbox_limits[4]) ||
            ipos > double(cloudbox_limits[5])) {
          return false;
        }

    return true;
  } else {
    // Pressure dimension
    double ipos = fractional_gp(gp_p);
    if (ipos <= double(cloudbox_limits[0]) ||
        ipos >= double(cloudbox_limits[1])) {
      return false;
    }

      // Latitude dimension
      ipos = fractional_gp(gp_lat);
      if (ipos <= double(cloudbox_limits[2]) ||
          ipos >= double(cloudbox_limits[3])) {
        return false;
      }

        // Longitude dimension
        ipos = fractional_gp(gp_lon);
        if (ipos <= double(cloudbox_limits[4]) ||
            ipos >= double(cloudbox_limits[5])) {
          return false;
        }
    return true;
  }
}

/*! Derives weights of a bin-type quadrature for arbitrary wide bins.
 *
 * Note: Rectangular and trapezoidal rule essentially give the same weights
 * when considering the same nodes x (not the mid-points between the nodes as
 * rectangular is sometimes using) and limiting the quadrature range by the
 * first and last node).
 * Keyword order=0 calculates rectangular bins (ie bins extend beyond the first
 * and last nodes), order=1 to trapezoidal bins (ie bins end exactly at nodes).
         
 \param w      resulting weights at ordinates.
 \param x      ordinates.
 \param order  order of quadrature (see above for details).
  
  \author Jana Mendrok, Daniel Kreyling
  \date 2017-06-16

*/
void bin_quadweights(Vector& w, const Vector& x, const Index& order) {
  Index nx = x.size();

  ARTS_ASSERT(nx > 1);
  ARTS_ASSERT(is_increasing(x));

  if (order == 0) {
    w[0] = std::min(x[1] - x[0],
               0.5 * (x[1] + x[0]));  // the latter is the half distance
                                      // from x0 to x1 plus the distance
                                      // to 0, ie 0.5(x1-x0)+x0.
    w[nx - 1] = x[nx - 1] - x[nx - 2];
  } else {
    w[0] = 0.5 * (x[1] - x[0]);
    w[nx - 1] = 0.5 * (x[nx - 1] - x[nx - 2]);
  }
  for (Index i = 1; i < nx - 1; i++) {
    w[i] = 0.5 * (x[i + 1] - x[i - 1]);
  }
}

//! Check whether field of a specific scattering species zero everywhere.
/*!
  \return empty_flag        flag whether all field entries are zero
  \param scat_species_field scattering species field (e.g. mass density,
                             mass flux, total number density)
  \param fieldname          name of scattering species field (just for info)
  \param dim                the atmosphere dimension 
  \param p_grid             pressure grid of current atmosphere
  \param lat_grid           latitude grid of current atmosphere
  \param lon_grid           longitude grid of current atmosphere

  \author Daniel Kreyling
  \date   2011-01-27
*/
void chk_scat_species_field(bool& empty_flag,
                            const Tensor3& scat_species_field,
                            const String& fieldname,
                            const Index& dim,
                            const Vector& p_grid,
                            const Vector& lat_grid,
                            const Vector& lon_grid) {
  // check p
  ARTS_USER_ERROR_IF (scat_species_field.npages() != p_grid.size(),
                      "The size of *p_grid* (", p_grid.size(),
                      ") is unequal the number of pages of *", fieldname, "* (",
                      scat_species_field.npages(), ").")

  // check lat
  if (dim >= 2) {
    ARTS_USER_ERROR_IF (scat_species_field.nrows() != lat_grid.size(),
        "The size of *lat_grid* (", lat_grid.size(),
         ") is unequal the number of rows of *", fieldname, "* (",
         scat_species_field.nrows(), ").")
  }

  // check lon
  if (dim == 3) {
    ARTS_USER_ERROR_IF (scat_species_field.ncols() != lon_grid.size(),
        "The size of *lon_grid* (", lon_grid.size(),
        ") is unequal the number of columns of *", fieldname, "* (",
        scat_species_field.ncols(), ").")
  }

  empty_flag = false;
  // set empty_flag to true if a single value of hydromet_field is unequal zero
  for (Index j = 0; j < scat_species_field.npages(); j++) {
    for (Index k = 0; k < scat_species_field.nrows(); k++) {
      for (Index l = 0; l < scat_species_field.ncols(); l++) {
        if (scat_species_field(j, k, l) != 0.0 &&
            !std::isnan(scat_species_field(j, k, l)))
          empty_flag = true;
        //	      if ( scat_species_field(j,k,l) != 0.0 ) empty_flag = true;
      }
    }
  }
}

/*! Parse atm_field_compact fieldname for species type

  \param  species_type  species indentifier (first part of field_name)
  \param  field_name    fieldname of atm_field_compact data entry
  \param  delim         delimiter string of field_name

  \author Jana Mendrok
  \date 2014-11-20

*/
void parse_atmcompact_speciestype(  //WS Output:
    String& species_type,
    // WS Input:
    const String& field_name,
    const String& delim) {
  ArrayOfString strarr;

  // split field_name string at '-' and write to ArrayOfString
  split(strarr, field_name, delim);

  // first entry is species type
  // (i.e. "abs_species" or "scat_species". or "T" or "z", which are ignored.)
  if (strarr.size() > 0 && field_name[0] != '-') {
    species_type = strarr[0];
  } else {
    ARTS_USER_ERROR (
                        "No information on field species type found in '",
                        field_name, "'\n")
  }
}

/*! Parse atm_field_compact fieldname for species name

  \param  species_type  species name (second part of field_name)
  \param  field_name    fieldname of atm_field_compact data entry
  \param  delim         delimiter string of field_name

  \author Jana Mendrok
  \date 2014-11-20

*/
void parse_atmcompact_speciesname(  //WS Output:
    String& species_name,
    // WS Input:
    const String& field_name,
    const String& delim) {
  ArrayOfString strarr;

  // split field_name string at '-' and write to ArrayOfString
  split(strarr, field_name, delim);

  // second entry is species name
  // (e.g. "H2O, "O3" etc. for abs_species or "IWC", "LWC" etc. for scat_species)
  if (strarr.size() > 1) {
    species_name = strarr[1];
  } else {
    ARTS_USER_ERROR (
                        "No information on field species name found in '",
                        field_name, "'\n")
  }
}

/*! Parse atm_field_compact fieldname for type of scat_species field

  \param  scat_type     species name (second part of field_name)
  \param  field_name    fieldname of atm_field_compact data entry
  \param  delim         delimiter string of field_name

  \author Jana Mendrok
  \date 2014-11-20

*/
void parse_atmcompact_scattype(  //WS Output:
    String& scat_type,
    // WS Input:
    const String& field_name,
    const String& delim) {
  ArrayOfString strarr;

  // split field_name string at '-' and write to ArrayOfString
  split(strarr, field_name, delim);

  // third entry is type of scat_species field
  // (e.g. "mass_density", "mass_flux", "number_density")
  if (strarr.size() > 2) {
    scat_type = strarr[2];
  } else {
    ARTS_USER_ERROR (
                        "No information on type of scat_species field found in '",
                        field_name, "'\n")
  }
}

/*! Splitting scat_species string and parse type of scattering species field

  \param  partfield_name name of atmospheric scattering species field
  \param  part_string    scattering species tag from *scat_species*
  \param  delim          delimiter string of *scat_species* elements

  \author Daniel Kreyling
  \date 2011-02-21

*/
void parse_partfield_name(  //WS Output:
    String& partfield_name,
    // WS Input:
    const String& part_string,
    const String& delim) {
  ArrayOfString strarr;

  // split scat_species string at delim and write to ArrayOfString
  split(strarr, part_string, delim);

  //first entry is scattering species field name (e.g. "IWC", "LWC" etc.)
  if (strarr.size() > 0 && part_string[0] != delim[0]) {
    partfield_name = strarr[0];
  } else {
    ARTS_USER_ERROR ("No information on scattering species field name in '",
                        part_string, "'\n")
  }
}

