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
   
   

#include <cmath>
#include <stdexcept>
#include <string>
#include "ppath.h"
#include "messages.h"
#include "math_funcs.h"
#include "absorption.h"
#include "abs_species_tags.h"
#include "physics_funcs.h"
#include "matpackIII.h"
#include "rte.h"
#include "rational.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PLANCK_CONST;
extern const Numeric BOHR_MAGNETON;
extern const Numeric LANDE_GS;

void phase_matrix(MatrixView K, const Numeric& theta, const Numeric& eta, const Index& DM);

void attenuation_matrix(MatrixView K, const Numeric& theta, const Numeric& eta, const Index& DM);

Numeric gs_caseb(const Numeric& N, const Numeric& J, const Numeric& S, const Numeric& GS);
Numeric gs_casea(const Numeric& Omega, const Numeric& J, const Numeric& Sigma, const Numeric& GS);

Numeric relative_strength(const Rational& m, const Rational& j, const Index& dj, const Index& dm);

Numeric frequency_change_caseb(const Rational& n, const Rational& m, const Rational& j, const Numeric& S, const Index& DJ, const Index& DM, const Index& DN, const Numeric& H_mag, const Numeric& GS);
Numeric frequency_change_casea(const Rational& omega, const Rational& m, const Rational& j, const Numeric& Sigma, const Index& DJ, const Index& DM, const Index& Domega, const Numeric& H_mag, const Numeric& GS);

void xsec_species_line_mixing_wrapper_with_zeeman(  
        Tensor3View part_abs_mat, 
	Tensor3View part_src_mat,
        const ArrayOfArrayOfSpeciesTag& abs_species, 
        const ArrayOfLineshapeSpec& abs_lineshape, 
        const ArrayOfLineRecord& lr, 
        const Vector&  Zeeman_DF,
        const SpeciesAuxData& isotopologue_ratios, 
        const Matrix& abs_t_nlte, 
        const Matrix& abs_vmrs, 
        const Vector& abs_p,
        const Vector& abs_t, 
        const Vector& f_grid, 
        const Numeric& lm_p_lim,
        const Numeric& theta, 
        const Numeric& eta, 
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

void set_quantum_numbers(Rational& Main,
                         Index& DMain,
                         Rational& J,
                         Index& DJ,
                         Rational& M,
                         Index& DM,
                         Numeric& S,
                         const LineRecord& temp_LR,
                         const Index hund,
                         const SpeciesAuxData& isotopologue_quantum,
                         const Index DO_Main,
                         const Index DO_J,
                         const Index DO_M);

void alter_linerecord(LineRecord& new_LR,
                      Numeric& Test_RS,
                      const LineRecord& old_LR,
                      const Rational& Main,
                      const Rational& M,
                      const Rational& J,
                      const Numeric&  S,
                      const Index&    DJ,
                      const Index&    DM,
                      const Index&    DMain,
                      const Numeric&  H_mag,
                      const Numeric&  GS,
                      Numeric (*frequency_change)(const Rational&, const  Rational&, const Rational&, 
                                    const Numeric&, const Index&, const Index&, 
                                    const Index&, const Numeric&, const Numeric&),
                      const Index& DO_RS,
                      const Index& DO_DF,
                      const Index& DO_QR);


void create_Zeeman_linerecordarrays(ArrayOfArrayOfLineRecord& aoaol,
                                    const ArrayOfArrayOfSpeciesTag& abs_species,
                                    const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                    const SpeciesAuxData& isotopologue_quantum,
                                    const Numeric& H_mag,
                                    const Index&DO_RS,
                                    const Index&DO_DF,
                                    const Index&DO_QR,
                                    const Index&DO_Main,
                                    const Index&DO_J,
                                    const Index&DO_M,
                                    const Verbosity& verbosity);

void set_part_isotopolouge_constants(Index& hund,Numeric& GS,const SpeciesAuxData& isotopologue_quantum,const LineRecord& temp_LR);