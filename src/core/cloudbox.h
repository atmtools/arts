/*!
     \file   cloudbox.h
     \author Claudia Emde <claudia.emde@dlr.de>
     \date   Thu May  23 14:34:05 2002
     
     \brief  Internal cloudbox functions.
     
   */

#ifndef cloudbox_h
#define cloudbox_h

#include "array.h"
#include "interpolation.h"
#include <matpack.h>
#include "optproperties.h"
#include "ppath_struct.h"

namespace Cloudbox {
   /** Global constant, minimum distance of cloudbox to lat/lon_grid edges.
    \author Patrick Eriksson, Jana Mendrok
    \date   2016-09-08
*/
  inline constexpr Numeric LAT_LON_MIN = 20;
}  // namespace Cloudbox

void chk_pnd_data(const GriddedField3& pnd_field_raw,
                  const String& pnd_field_file);

void chk_pnd_raw_data(const ArrayOfGriddedField3& pnd_field_raw,
                      const String& pnd_field_file);

void chk_pnd_field_raw_only_in_cloudbox(
    const Index& dim,
    const ArrayOfGriddedField3& pnd_field_raw,
    ConstVectorView p_grid,
    ConstVectorView lat_grid,
    ConstVectorView lon_grid,
    const ArrayOfIndex& cloudbox_limits);

void chk_scat_species(const ArrayOfString& scat_species, const String& delim);

void chk_scattering_data(const ArrayOfSingleScatteringData& scat_data,
                         const ArrayOfScatteringMetaData& scat_meta);

void chk_scat_data(const SingleScatteringData& scat_data);

bool is_gp_inside_cloudbox(const GridPos& gp_p,
                           const GridPos& gp_lat,
                           const GridPos& gp_lon,
                           const ArrayOfIndex& cloudbox_limits,
                           const bool& include_boundaries);

bool is_inside_cloudbox(const Ppath& ppath_step,
                        const ArrayOfIndex& cloudbox_limits,
                        const bool include_boundaries);

void bin_quadweights(Vector& w, const Vector& x, const Index& order = 1);

void chk_scat_species_field(bool& empty_flag,
                            const Tensor3& scat_species_field,
                            const String& fieldname,
                            const Index& dim,
                            const Vector& p_grid,
                            const Vector& lat_grid,
                            const Vector& lon_grid);

void parse_atmcompact_speciestype(String& species_type,
                                  const String& field_name,
                                  const String& delim);

void parse_atmcompact_speciesname(String& species_name,
                                  const String& field_name,
                                  const String& delim);

void parse_atmcompact_scattype(String& scat_type,
                               const String& field_name,
                               const String& delim);

void parse_partfield_name(String& partfield_name,
                          const String& part_string,
                          const String& delim);

#endif  //cloudbox_h
