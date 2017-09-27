/* Copyright (C) 2011-2017 Jana Mendrok <jana.mendrok@gmail.com>
                      
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
  \file   microphysics.h
  \author Jana Mendrok, Daniel Kreyling, Manfred Brath, Patrick Eriksson
  \date   2017-07-10 
  
  \brief  Internal functions for microphysics calculations (size distributions etc.)
*/


#ifndef microphysics_h
#define microphysics_h

#include "matpackVII.h"
#include "interpolation.h"
#include "optproperties.h"
#include "array.h"
#include "gridded_fields.h"
#include "ppath.h"
#include "messages.h"

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

void pnd_fieldSB06 (Tensor4View pnd_field,
                   const Tensor3& WC_field,
                   const Tensor3& N_field,
                   const ArrayOfIndex& limits,
                   const ArrayOfArrayOfScatteringMetaData& scat_meta,
                   const Index& scat_species,
                   const String& part_string,
                   const String& delim,
                   const Verbosity& verbosity);

void pnd_fieldMY05 (Tensor4View pnd_field,
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

void psd_SB06 (Vector& psd,
              Matrix& dpsd,
              const Vector& mass,
              const Numeric& N_tot,
              const Numeric& WC,
              const String& hydrometeor_type);

void psd_MY05 (Vector& psd,
              Matrix& dpsd,
              const Vector& diameter_max,
              const Numeric N_tot,
              const Numeric WC,
              const String psd_type);

Numeric IWCtopnd_H11 (const Numeric diameter_mass_equivalent,
                      const Numeric t);

Numeric IWCtopnd_H13 (const Numeric diameter_mass_equivalent,
                      const Numeric t);

Numeric IWCtopnd_H13Shape (const Numeric diameter_mass_equivalent,
                           const Numeric t);

Numeric area_ratioH13 (const Numeric diameter_mass_equivalent,
                       const Numeric t);

Numeric LWCtopnd_MGD_LWC ( const Numeric d, const Numeric m, const Numeric lwc);

Numeric IWCtopnd_MGD_IWC ( const Numeric d, const Numeric m, const Numeric iwc);

Numeric LWCtopnd (const Numeric lwc,
                  const Numeric radius);

Numeric PRtopnd_MP48 (const Numeric R,
                      const Numeric diameter_melted_equivalent);

#endif //microphysics_h

