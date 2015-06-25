/* Copyright (C) 2015
   Richard Larsson <ric.larsson@gmail.com>

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

#include "absorption.h"


void GetLineScalingData(Numeric& partition_ratio, 
                        Numeric& boltzmann_ratio, 
                        Numeric& abs_nlte_ratio, 
                        Numeric& src_nlte_ratio, 
                        const SpeciesAuxData::AuxType& partition_type,
                        const ArrayOfGriddedField1& partition_data,
                        const Numeric& atm_t,
                        const Numeric& line_t,
                        const Numeric& line_f,
                        const Numeric& line_elow,
                        const bool&    do_nlte,
                        const Numeric& line_evlow,
                        const Numeric& line_evupp,
                        const Index& line_evlow_index,
                        const Index& line_evupp_index,
                        ConstVectorView atm_t_nlte);

void CalculatePartitionFctFromData( Numeric& q_ref, 
                                    Numeric& q_t, 
                                    const Numeric& ref, 
                                    const Numeric& t,
                                    ConstVectorView t_grid, 
                                    ConstVectorView q_grid, 
                                    const Index& interp_order);

void CalculatePartitionFctFromCoeff(Numeric& q_ref, 
                                    Numeric& q_t, 
                                    const Numeric& ref, 
                                    const Numeric& t,
                                    ConstVectorView q_grid);