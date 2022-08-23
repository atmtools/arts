/* Copyright (C) 2002-2012 Claudia Emde <claudia.emde@dlr.de>
                      
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. 
*/

/*!
     \file   cloudbox.h
     \author Claudia Emde <claudia.emde@dlr.de>
     \date   Thu May  23 14:34:05 2002
     
     \brief  Internal cloudbox functions.
     
   */

#ifndef cloudbox_h
#define cloudbox_h

#include "array.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "matpackVII.h"
#include "messages.h"
#include "optproperties.h"
#include "ppath.h"

namespace Cloudbox {
   /** Global constant, minimum distance of cloudbox to lat/lon_grid edges.
    \author Patrick Eriksson, Jana Mendrok
    \date   2016-09-08
*/
  inline constexpr Numeric LAT_LON_MIN = 20;
}  // namespace Cloudbox

void chk_pnd_data(const GriddedField3& pnd_field_raw,
                  const String& pnd_field_file,
                  const Index& atmosphere_dim,
                  const Verbosity& verbosity);

void chk_pnd_raw_data(const ArrayOfGriddedField3& pnd_field_raw,
                      const String& pnd_field_file,
                      const Index& atmosphere_dim,
                      const Verbosity& verbosity);

void chk_pnd_field_raw_only_in_cloudbox(
    const Index& dim,
    const ArrayOfGriddedField3& pnd_field_raw,
    ConstVectorView p_grid,
    ConstVectorView lat_grid,
    ConstVectorView lon_grid,
    const ArrayOfIndex& cloudbox_limits);

void chk_scat_species(const ArrayOfString& scat_species, const String& delim);

void chk_scattering_data(const ArrayOfSingleScatteringData& scat_data,
                         const ArrayOfScatteringMetaData& scat_meta,
                         const Verbosity& verbosity);

void chk_scattering_meta_data(const ScatteringMetaData& scat_meta_single,
                              const String& scat_meta_file,
                              const Verbosity& verbosity);

void chk_scat_data(const SingleScatteringData& scat_data,
                   const Verbosity& verbosity);

bool is_gp_inside_cloudbox(const GridPos& gp_p,
                           const GridPos& gp_lat,
                           const GridPos& gp_lon,
                           const ArrayOfIndex& cloudbox_limits,
                           const bool& include_boundaries,
                           const Index& atmosphere_dim = 3);

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

void find_cloudlimits(Index& lower,
                      Index& upper,
                      const Tensor3& scat_species_field,
                      const Index& atmosphere_dim,
                      const Numeric& cloudbox_margin);

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
