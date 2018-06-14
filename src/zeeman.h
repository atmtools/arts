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

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PLANCK_CONST;
extern const Numeric BOHR_MAGNETON;
extern const Numeric LANDE_GS;


void xsec_species_line_mixing_wrapper_with_zeeman(  ArrayOfPropagationMatrix& propmat_clearsky, 
                                                    ArrayOfStokesVector& nlte_source,
                                                    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                                                    ArrayOfStokesVector&  dnlte_dx_source,
                                                    ArrayOfStokesVector&  nlte_dsource_dx,
                                                    const ArrayOfArrayOfSpeciesTag& abs_species, 
                                                    const ArrayOfRetrievalQuantity& flag_partials,
                                                    const ArrayOfIndex& flag_partials_position,
                                                    const Index& abs_lineshape_ls, 
                                                    const Index& abs_lineshape_lsn, 
                                                    const Numeric& abs_lineshape_cutoff, 
                                                    const ArrayOfLineRecord& lr,
                                                    const Vector&  planck_BT,
                                                    const Matrix&  dplanck_BT,
                                                    const SpeciesAuxData& isotopologue_ratios, 
                                                    const SpeciesAuxData& partition_functions,
                                                    const Matrix& abs_nlte, 
                                                    const Matrix& abs_vmrs, 
                                                    const Vector& abs_p,
                                                    const Vector& abs_t, 
                                                    const Vector& f_grid,
                                                    const Vector& rtp_mag,
                                                    const Vector& r_path_los,
                                                    const Numeric& lm_p_lim,
                                                    const Numeric& theta, 
                                                    const Numeric& eta, 
                                                    const Numeric& H_mag, 
                                                    const Index& DM, 
                                                    const Index& this_species,
                                                    const Verbosity& verbosity );

void set_magnetic_parameters(Numeric& H_mag,
                             Numeric& eta,
                             Numeric& theta,
                             const Index manual_zeeman_tag,
                             const Numeric& manual_zeeman_eta,
                             const Numeric& manual_zeeman_theta,
                             const Numeric& manual_zeeman_magnetic_field_strength,
                             ConstVectorView rtp_mag,
                             ConstVectorView r_path_los);

void set_magnetic_parameters_derivative(Numeric& dH_du,
                                        Numeric& dH_dv,
                                        Numeric& dH_dw,
                                        Numeric& deta_du,
                                        Numeric& deta_dv,
                                        Numeric& deta_dw,
                                        Numeric& dtheta_du,
                                        Numeric& dtheta_dv,
                                        Numeric& dtheta_dw,
                                        ConstVectorView rtp_mag,
                                        ConstVectorView r_path_los);

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

Index part_mag_strength(const ArrayOfRetrievalQuantity& flag_partials);
Index part_mag_theta(const ArrayOfRetrievalQuantity& flag_partials);

void zeeman_on_the_fly(ArrayOfPropagationMatrix& propmat_clearsky, 
                       ArrayOfStokesVector& nlte_source,
                       ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
                       ArrayOfStokesVector& dnlte_dx_source,
                       ArrayOfStokesVector& nlte_dsource_dx,
                       const ArrayOfArrayOfSpeciesTag& abs_species, 
                       const ArrayOfRetrievalQuantity& flag_partials,
                       const ArrayOfIndex& flag_partials_positions,
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
                       const Numeric& lm_p_lim,
                       const Index& manual_zeeman_tag,
                       const Numeric& manual_zeeman_magnetic_field_strength,
                       const Numeric& manual_zeeman_theta,
                       const Numeric& manual_zeeman_eta,
                       const Verbosity& verbosity);

Numeric zeeman_magnetic_magnitude(const Numeric& u, const Numeric& v, const Numeric& w);
Numeric zeeman_magnetic_dmagnitude_du(const Numeric& u, const Numeric& v, const Numeric& w);
Numeric zeeman_magnetic_dmagnitude_dv(const Numeric& u, const Numeric& v, const Numeric& w);
Numeric zeeman_magnetic_dmagnitude_dw(const Numeric& u, const Numeric& v, const Numeric& w);
Numeric zeeman_magnetic_theta(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a);
Numeric zeeman_magnetic_dtheta_du(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a);
Numeric zeeman_magnetic_dtheta_dv(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a);
Numeric zeeman_magnetic_dtheta_dw(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a);
Numeric zeeman_magnetic_eta(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a);
Numeric zeeman_magnetic_deta_du(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a);
Numeric zeeman_magnetic_deta_dv(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a);
Numeric zeeman_magnetic_deta_dw(const Numeric& u, const Numeric& v, const Numeric& w, const Numeric& z, const Numeric& a);
