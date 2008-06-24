/* Copyright (C) 2002-2008 Claudia Emde <claudia.emde@dlr.de>
                      
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

void chk_if_pnd_zero_p(
                       const Index& i_p,
                       const GField3& pnd_field_raw,
                       const String& pnd_field_file);

void chk_if_pnd_zero_lat(
                       const Index& i_lat,
                       const GField3& pnd_field_raw,
                       const String& pnd_field_file);

void chk_if_pnd_zero_lon(
                       const Index& i_lon,
                       const GField3& pnd_field_raw,
                       const String& pnd_field_file);

void chk_pnd_data(
                  const GField3& pnd_field_raw,
                  const String& pnd_field_file,
                  const Index& atmosphere_dim,
                  ConstVectorView p_grid,
                  ConstVectorView lat_grid,
                  ConstVectorView lon_grid,
                  const ArrayOfIndex& cloudbox_limits);

void chk_pnd_raw_data(
                      const ArrayOfGField3& pnd_field_raw,
                      const String& pnd_field_file,
                      const Index& atmosphere_dim,
                      ConstVectorView p_grid,
                      ConstVectorView lat_grid,
                      ConstVectorView lon_grid,
                      const ArrayOfIndex& cloudbox_limits);

void chk_single_scattering_data(
                                const SingleScatteringData& scat_data_raw,
                                const String& scat_data_file,
                                ConstVectorView f_grid);

void iy_interp_cloudbox_field(
            Matrix&         iy,
      const Tensor7&        scat_i_p,
      const Tensor7&        scat_i_lat,
      const Tensor7&        scat_i_lon,
      const Tensor4&        doit_i_field1D_spectrum, 
      const GridPos&        rte_gp_p,
      const GridPos&        rte_gp_lat,
      const GridPos&        rte_gp_lon,
      const Vector&         rte_los,
      const Index&          cloudbox_on,
      const ArrayOfIndex&   cloudbox_limits,
      const Index&          atmosphere_dim,
      const Index&          stokes_dim,
      const Vector&         scat_za_grid,
      const Vector&         scat_aa_grid,
      const Vector&         f_grid,
      const String&         interpmeth );

#endif //cloudbox_h

