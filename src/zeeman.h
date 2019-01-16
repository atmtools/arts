/* Copyright (C) 2014
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
   
   
#include "abs_species_tags.h"
#include "physics_funcs.h"
#include "rte.h"
#include "global_data.h"
#include "quantum.h"


void alter_linerecord( LineRecord& new_LR,
                       Numeric& Test_RS,
                       const Numeric& old_LS,
                       const Rational& J_up,
                       const Rational& J_lo,
                       const Rational& M_up,
                       const Rational& M_lo);


void create_Zeeman_linerecordarrays(ArrayOfArrayOfLineRecord& aoaol,
                                    const ArrayOfArrayOfSpeciesTag& abs_species,
                                    const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                    const Verbosity& verbosity);

void zeeman_on_the_fly(ArrayOfPropagationMatrix& propmat_clearsky, 
                       ArrayOfStokesVector& nlte_source,
                       ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                       ArrayOfStokesVector& dnlte_dx_source,
                       ArrayOfStokesVector& nlte_dsource_dx,
                       const ArrayOfArrayOfSpeciesTag& abs_species, 
                       const ArrayOfRetrievalQuantity& flag_partials,
                       const ArrayOfArrayOfLineRecord& zeeman_linerecord_precalc,
                       const SpeciesAuxData& isotopologue_ratios, 
                       const SpeciesAuxData& partition_functions,
                       const ConstVectorView f_grid,
                       const ConstVectorView rtp_vmrs, 
                       const ConstVectorView rtp_nlte, 
                       const ConstVectorView rtp_mag,
                       const ConstVectorView rtp_los,
                       const Numeric& rtp_pressure,
                       const Numeric& rtp_temperature,
                       const Index& manual_zeeman_tag,
                       const Numeric& manual_zeeman_magnetic_field_strength,
                       const Numeric& manual_zeeman_theta,
                       const Numeric& manual_zeeman_eta);
