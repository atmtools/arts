/* Copyright (C) 2003 Claudia Emde <claudia@sat.physik.uni-bremen.de>

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
   USA. */



/*===========================================================================
  ===  File description 
  ===========================================================================*/

/*!
  \file   scatrte.h
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
  \date   2003-06-03
  
  \brief  Radiative transfer in cloudbox.
  
  This file contains functions related to the radiative transfer in the 
  cloudbox.
*/


#ifndef scatrte_h
#define scatrte_h

void cloud_fieldsCalc(// Output:
                        Tensor5View ext_mat_field,
                        Tensor4View abs_vec_field,
                        Index& scat_p_index,
                        Index& scat_lat_index,
                        Index& scat_lon_index, 
                        Tensor3& ext_mat,
                        Matrix& abs_vec,  
                        // Input:
                        const Agenda& spt_calc_agenda,
                        const Agenda& opt_prop_part_agenda,
                        const ArrayOfIndex& cloudbox_limits
                        );

void cloud_ppath_update1D(
                  Tensor6View i_field,
                  // scalar_gas_abs_agenda:
                  Numeric& rte_pressure,
                  Numeric& rte_temperature,
                  Vector& rte_vmr_list,
                  // opt_prop_xxx_agenda:
                  Tensor3& ext_mat,
                  Matrix& abs_vec, 
		  // ground related variables STR
                  Matrix& ground_los,
                  Matrix& ground_emission,
                  Tensor4& ground_refl_coeffs,
                  Vector& rte_los,
		  Vector& rte_pos,
		  GridPos& rte_gp_p,
                  // ppath_step_agenda:
                  Ppath& ppath_step, 
                  const Index& p_index,
                  const Index& scat_za_index,
                  ConstVectorView scat_za_grid,
                  const ArrayOfIndex& cloudbox_limits,
                  ConstTensor6View scat_field,
                  // Calculate scalar gas absorption:
                  const Agenda& scalar_gas_absorption_agenda,
                  ConstTensor4View vmr_field,
                  // Gas absorption:
                  const Agenda& opt_prop_gas_agenda,
                  // Propagation path calculation:
                  const Agenda& ppath_step_agenda,
                  ConstVectorView p_grid,
                  ConstTensor3View z_field,
                  ConstMatrixView r_geoid,
                  // Calculate thermal emission:
                  ConstTensor3View t_field,
                  ConstVectorView f_grid,
                  const Index& f_index,
                  //particle opticla properties
                  ConstTensor5View ext_mat_field,
                  ConstTensor4View abs_vec_field,
		  const Agenda& ground_refl_agenda //STR
                  );

void cloud_ppath_update1D_planeparallel(
                  Tensor6View i_field,
                  // scalar_gas_abs_agenda:
                  Numeric& rte_pressure,
                  Numeric& rte_temperature,
                  Vector& rte_vmr_list,
                  // opt_prop_xxx_agenda:
                  Tensor3& ext_mat,
                  Matrix& abs_vec, 
		  // ground related variables STR
                  Matrix& ground_los,
                  Matrix& ground_emission,
                  Tensor4& ground_refl_coeffs,
                  Vector& rte_los,
		  Vector& rte_pos,
		  GridPos& rte_gp_p,
                  // ppath_step_agenda:
                  Ppath& ppath_step, 
                  const Index& p_index,
                  const Index& scat_za_index,
                  ConstVectorView scat_za_grid,
                  const ArrayOfIndex& cloudbox_limits,
                  ConstTensor6View scat_field,
                  // Calculate scalar gas absorption:
                  const Agenda& scalar_gas_absorption_agenda,
                  ConstTensor4View vmr_field,
                  // Gas absorption:
                  const Agenda& opt_prop_gas_agenda,
                  // Propagation path calculation:
                  const Agenda& ppath_step_agenda,
                  ConstVectorView p_grid,
                  ConstTensor3View z_field,
                  ConstMatrixView r_geoid,
                  // Calculate thermal emission:
                  ConstTensor3View t_field,
                  ConstVectorView f_grid,
                  const Index& f_index,
                  //particle opticla properties
                  ConstTensor5View ext_mat_field,
                  ConstTensor4View abs_vec_field,
		  const Agenda& ground_refl_agenda //STR
                  );


void cloud_ppath_update3D(
			  Tensor6View i_field,
			  // scalar_gas_abs_agenda:
			  Numeric& rte_pressure,
			  Numeric& rte_temperature,
			  Vector& rte_vmr_list,
			  // opt_prop_xxx_agenda:
			  Tensor3& ext_mat,
			  Matrix& abs_vec,  
			  // ppath_step_agenda:
			  Ppath& ppath_step, 
			  const Index& p_index,
                          const Index& lat_index,
                          const Index& lon_index,
			  const Index& scat_za_index,
                          const Index& scat_aa_index,
			  ConstVectorView scat_za_grid,
                          ConstVectorView scat_aa_grid,
			  const ArrayOfIndex& cloudbox_limits,
			  ConstTensor6View scat_field,
			  // Calculate scalar gas absorption:
			  const Agenda& scalar_gas_absorption_agenda,
			  ConstTensor4View vmr_field,
			  // Gas absorption:
			  const Agenda& opt_prop_gas_agenda,
			  // Propagation path calculation:
			  const Agenda& ppath_step_agenda,
			  ConstVectorView p_grid,
                          ConstVectorView lat_grid,
                          ConstVectorView lon_grid,
			  ConstTensor3View z_field,
			  ConstMatrixView r_geoid,
			  // Calculate thermal emission:
			  ConstTensor3View t_field,
			  ConstVectorView f_grid,
			  const Index& f_index,
			  //particle opticla properties
			  ConstTensor5View ext_mat_field,
			  ConstTensor4View abs_vec_field
			  );


void ppath_step_in_cloudbox(Ppath& ppath_step,
                            const Agenda& ppath_step_agenda,
                            const Index& p,
                            const Index& lat, 
                            const Index& lon,
                            ConstTensor3View z_field,
                            ConstMatrixView r_geoid,
                            ConstVectorView scat_za_grid,
                            ConstVectorView aa_grid,
                            const Index& scat_za_index,
                            const Index& scat_aa_index,
                            ConstVectorView lat_grid,
                            ConstVectorView lon_grid);

bool is_inside_cloudbox(const Ppath& ppath_step,
                        const ArrayOfIndex& cloudbox_limits);

#endif //scatrte_h
