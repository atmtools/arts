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

void phase_matrix(MatrixView K, const Numeric& theta, const Numeric& eta, const Index& DM);

void attenuation_matrix(MatrixView K, const Numeric& theta, const Numeric& eta, const Index& DM);

void dphase_matrix_dtheta(MatrixView dK, const Numeric& theta, const Numeric& eta, const Index& DM);
void dphase_matrix_deta(  MatrixView dK, const Numeric& theta, const Numeric& eta, const Index& DM);

void dattenuation_matrix_dtheta(MatrixView dK, const Numeric& theta, const Numeric& eta, const Index& DM);
void dattenuation_matrix_deta(  MatrixView dK, const Numeric& theta, const Numeric& eta, const Index& DM);

Rational relative_strength(const Rational& m, const Rational& j, const Index& dj, const Index& dm);

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
                                                    const Matrix& abs_t_nlte, 
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
