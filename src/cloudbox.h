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

#include "matpackVII.h"
#include "interpolation.h"
#include "optproperties.h"
#include "array.h"
#include "gridded_fields.h"
#include "ppath.h"
#include "messages.h"

void chk_if_pnd_zero_p (const Index& i_p,
                        const GriddedField3& pnd_field_raw,
                        const String& pnd_field_file,
                        const Verbosity& verbosity);

void chk_if_pnd_zero_lat (const Index& i_lat,
                          const GriddedField3& pnd_field_raw,
                          const String& pnd_field_file,
                          const Verbosity& verbosity);

void chk_if_pnd_zero_lon (const Index& i_lon,
                          const GriddedField3& pnd_field_raw,
                          const String& pnd_field_file,
                          const Verbosity& verbosity);

void chk_pnd_data (const GriddedField3& pnd_field_raw,
                   const String& pnd_field_file,
                   const Index& atmosphere_dim,
                   const Verbosity& verbosity);

void chk_pnd_raw_data (const ArrayOfGriddedField3& pnd_field_raw,
                       const String& pnd_field_file,
                       const Index& atmosphere_dim,
                       const Verbosity& verbosity);

void chk_pnd_field_raw_only_in_cloudbox(
        const Index&                 dim,
        const ArrayOfGriddedField3&  pnd_field_raw,
        ConstVectorView              p_grid,
        ConstVectorView              lat_grid,
        ConstVectorView              lon_grid,
        const ArrayOfIndex&          cloudbox_limits);

void chk_scat_species (const ArrayOfString& scat_species,
                       const String& delim);

void chk_scattering_data (const ArrayOfSingleScatteringData& scat_data,
                          const ArrayOfScatteringMetaData& scat_meta,
                          const Verbosity& verbosity);

void chk_scattering_meta_data (const ScatteringMetaData& scat_meta_single,
                               const String& scat_meta_file,
                               const Verbosity& verbosity);

void chk_scat_data_fgrid (const SingleScatteringData& scat_data,
                          ConstVectorView f_grid,
                          const String& infostring);

void chk_scat_data (const SingleScatteringData& scat_data,
                    const Verbosity& verbosity);

bool is_gp_inside_cloudbox(
   const GridPos&      gp_p,
   const GridPos&      gp_lat,
   const GridPos&      gp_lon,
   const ArrayOfIndex& cloudbox_limits,
   const bool&         include_boundaries,
   const Index&        atmosphere_dim=3 );

bool is_inside_cloudbox (const Ppath& ppath_step,
                         const ArrayOfIndex& cloudbox_limits,
                         const bool include_boundaries);


void pnd_fieldMH97 (Tensor4View pnd_field,
                    const Tensor3& IWC_field,
                    const Tensor3& t_field,
                    const ArrayOfIndex& limits,
                    const ArrayOfArrayOfScatteringMetaData& scat_meta,
                    const Index& scat_species,
                    const String& part_string,
                    const String& delim,
                    const Verbosity& verbosity);

void pnd_fieldH11 (Tensor4View pnd_field,
                   const Tensor3& IWC_field,
                   const Tensor3& t_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfArrayOfScatteringMetaData& scat_meta,
                   const Index& scat_species,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity);

void pnd_fieldH13 (Tensor4View pnd_field,
                   const Tensor3& IWC_field,
                   const Tensor3& t_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfArrayOfScatteringMetaData& scat_meta,
                   const Index& scat_species,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity);

void pnd_fieldH13Shape (Tensor4View pnd_field,
                        const Tensor3& IWC_field,
                        const Tensor3& t_field,
                        const ArrayOfIndex& limits,
                        const ArrayOfArrayOfScatteringMetaData& scat_meta,
                        const Index& scat_species,
                        const String& part_string,
                        const String& delim,
                        const Verbosity& verbosity);

void pnd_fieldF07 (Tensor4View pnd_field,
                   const Tensor3& SWC_field,
                   const Tensor3& t_field,
                   const String& regime,
                   const ArrayOfIndex& limits,
                   const ArrayOfArrayOfScatteringMetaData& scat_meta,
                   const Index& scat_species,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity);

void pnd_fieldS2M (Tensor4View pnd_field,
                   const Tensor3& WC_field,
                   const Tensor3& N_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfArrayOfScatteringMetaData& scat_meta,
                   const Index& scat_species,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity);

void pnd_fieldMGD_LWC (Tensor4View pnd_field,
                       const Tensor3& LWC_field,
                       const ArrayOfIndex& limits,
                       const ArrayOfArrayOfScatteringMetaData& scat_meta,
                       const Index& scat_species,
                       const String& part_string,
                       const String& delim,
                       const Verbosity& verbosity);


void pnd_fieldMGD_IWC (Tensor4View pnd_field,
                       const Tensor3& IWC_field,
                       const ArrayOfIndex& limits,
                       const ArrayOfArrayOfScatteringMetaData& scat_meta,
                       const Index& scat_species,
                       const String& part_string,
                       const String& delim,
                       const Verbosity& verbosity);


void pnd_fieldMP48 (Tensor4View pnd_field,
                    const Tensor3& PR_field,
                    const ArrayOfIndex& limits,
                    const ArrayOfArrayOfScatteringMetaData& scat_meta,
                    const Index& scat_species,
                    const String& part_string,
                    const String& delim,
                    const Verbosity& verbosity);

void pnd_fieldW16 (Tensor4View pnd_field,
                    const Tensor3& RWC_field,
                    const ArrayOfIndex& limits,
                    const ArrayOfArrayOfScatteringMetaData& scat_meta,
                    const Index& scat_species,
                    const String& part_string,
                    const String& delim,
                    const Verbosity& verbosity);

void pnd_fieldH98 (Tensor4View pnd_field,
                   const Tensor3& LWC_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfArrayOfScatteringMetaData& scat_meta,
                   const Index& scat_species,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity);

void psd_general_MGD ( Vector& psd,
                   const Vector& diameter,
                   const Numeric& N0,
                   const Numeric& mu,
                   const Numeric& lambda,
                   const Numeric& gamma );

void psd_cloudice_MH97 ( Vector& psd, 
                   const Vector& diameter,
                   const Numeric& iwc,
                   const Numeric& t,
                   const bool noisy );

void psd_rain_W16 ( Vector& psd,
                   const Vector& diameter,
                   const Numeric& rwc );

void psd_snow_F07 ( Vector& psd,
                   const Vector& diameter,
                   const Numeric& swc,
                   const Numeric& t,
                   const Numeric alpha,
                   const Numeric beta,
                   const String& regime );

Numeric IWCtopnd_H11 (const Numeric diameter_mass_equivalent,
                      const Numeric t);

Numeric IWCtopnd_H13 (const Numeric diameter_mass_equivalent,
                      const Numeric t);

Numeric IWCtopnd_H13Shape (const Numeric diameter_mass_equivalent,
                           const Numeric t);

Numeric area_ratioH13 (const Numeric diameter_mass_equivalent,
                       const Numeric t);

Numeric WCtopnd_S2M (const Numeric mass,
                     const Numeric N_tot,
                     const Numeric M,
                     const String psd_type);

Numeric LWCtopnd_MGD_LWC ( const Numeric d, const Numeric m, const Numeric lwc);

Numeric IWCtopnd_MGD_IWC ( const Numeric d, const Numeric m, const Numeric iwc);

Numeric LWCtopnd (const Numeric lwc,
                  const Numeric radius);

Numeric PRtopnd_MP48 (const Numeric R,
                      const Numeric diameter_melted_equivalent);

void bin_quadweights( Vector& w,
                const Vector& x,
                const Index& order=1 );

void bin_integral( Vector& w,
             const Vector& x,
             const Vector& y);

void chk_pndsum (Vector& pnd,
                 const Numeric xwc,
                 const Vector& mass,
                 const Index& p,
                 const Index& lat,
                 const Index& lon,
                 const String& part_type,
                 const Verbosity& verbosity);

void chk_scat_species_field(bool& empty_flag,
                            const Tensor3& scat_species_field, 
                            const String& fieldname,
                            const Index&  dim,	
                            const Vector& p_grid,
                            const Vector& lat_grid,
                            const Vector& lon_grid);

void find_cloudlimits(Index&          lower,
                      Index&          upper,
                      const Tensor3&  scat_species_field,
                      const Index&    atmosphere_dim,
                      const Numeric&  cloudbox_margin);

void parse_atmcompact_speciestype (String& species_type,
                                   const String& field_name,
                                   const String& delim);

void parse_atmcompact_speciesname (String& species_name,
                                   const String& field_name,
                                   const String& delim);

void parse_atmcompact_scattype (String& scat_type,
                                const String& field_name,
                                const String& delim);

void parse_partfield_name (String& partfield_name,
                      const String& part_string,
                      const String& delim);

void parse_psd_param (String& psd_param,
                      const String& part_string,
                      const String& delim);

void parse_psd_options (ArrayOfString& psd_options,
                        const String& part_string,
                        const String& delim);

#endif //cloudbox_h

