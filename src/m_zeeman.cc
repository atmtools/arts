/* Copyright (C) 2012
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

#include "auto_md.h"
#include "zeeman.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddZeeman( Tensor4& propmat_clearsky,
			        Tensor4& propmat_source_clearsky,
				const Vector& f_grid,
				const ArrayOfArrayOfSpeciesTag& abs_species,
				const ArrayOfArrayOfLineRecord& abs_lines_per_species,
				const ArrayOfLineshapeSpec& abs_lineshape,
				const SpeciesAuxData& isotopologue_ratios,
				const SpeciesAuxData& isotopologue_quantum,
				const Numeric& rtp_pressure,
				const Numeric& rtp_temperature,
				const Numeric& lm_p_lim,
                                const Vector& rtp_temperature_nlte,
				const Vector& rtp_vmr,
				const Vector& rtp_mag,
				const Vector& ppath_los,
				const Index& atmosphere_dim,
				const Index& manual_zeeman_tag,
				const Numeric& manual_zeeman_magnetic_field_strength,
				const Numeric& manual_zeeman_theta,
				const Numeric& manual_zeeman_eta,
				const Verbosity& verbosity)
{
  CREATE_OUT3;

  // Check that correct isotopologue ratios are defined for the species
  // we want to calculate
  checkIsotopologueRatios(abs_species, isotopologue_ratios);

  bool do_src = true;
  {// Begin TEST(s)
  if (abs_species.nelem() != abs_lines_per_species.nelem())
      throw std::runtime_error("Dimension of *abs_species* and *abs_lines_per_species* don't match.");
  if( propmat_clearsky.ncols()  != 4 )
      throw std::runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
  if( propmat_clearsky.nrows()  != 4 )
      throw std::runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
  if( propmat_clearsky.npages() != f_grid.nelem() )
      throw std::runtime_error("Frequency dimension of *propmat_clearsky* not equal to length of *f_grid*.");
  if( propmat_clearsky.nbooks() != abs_species.nelem() )
      throw std::runtime_error("Species dimension of *propmat_clearsky* not equal to length of *abs_species*.");
  if( rtp_mag.nelem() != 3 )
      throw std::runtime_error("*rtp_mag* must have length 3.");
  if( atmosphere_dim != 3 )
      throw std::runtime_error("*atmosphere_dim* must be 3.  Zeeman Effect is only implemented for 3D geometry.");
  if( ppath_los.nelem() != 2 )
      throw std::runtime_error("*ppath_los* is not set correctly.");
  if( propmat_source_clearsky.ncols()  != 4 )
    if( propmat_source_clearsky.ncols()  != 0 ) 
      throw std::runtime_error("Zeeman Effect is only implemented for source Stokes dimension 4 or 0.");
  if( propmat_source_clearsky.nrows()  != 4 )
    if( propmat_source_clearsky.nrows()  != 0 ) 
      throw std::runtime_error("Zeeman Effect is only implemented for source Stokes dimension 4 or 0.");
  if( propmat_source_clearsky.npages() != f_grid.nelem() )
    if( propmat_source_clearsky.npages() != 0 ) 
      throw std::runtime_error("Frequency dimension of *propmat_source_clearsky* not equal to length of *f_grid* or 0.");
  if( propmat_source_clearsky.nbooks() != abs_species.nelem() )
  {
    if( propmat_source_clearsky.nbooks() != 0 )
      throw std::runtime_error("Species dimension of *propmat_clearsky* not equal to length of *abs_species* or 0.");
    else
      do_src = false; // Note that test is here beacuse it makes things easier.
  }
  }// End   TEST(s)

  Vector R_path_los;
  mirror_los(R_path_los, ppath_los, atmosphere_dim);

  // Using the standard scalar absorption functions to get physics parameters,
  Vector abs_p, abs_t; Matrix abs_vmrs, abs_t_nlte;
  AbsInputFromRteScalars( abs_p, abs_t, abs_t_nlte, abs_vmrs,            // Output
          rtp_pressure, rtp_temperature, rtp_temperature_nlte, rtp_vmr,  //Input
          verbosity);                                                    // Verbose!

  // FOR LOG:  Loss of speed when mag == 0
  // Set the magnetic parameters...
  Numeric H_mag,eta,theta;
  set_magnetic_parameters(H_mag,
                          eta,
                          theta,
                          manual_zeeman_tag,
                          manual_zeeman_eta,
                          manual_zeeman_theta,
                          manual_zeeman_magnetic_field_strength,
                          rtp_mag,
                          R_path_los);
  
  ArrayOfArrayOfLineRecord aoaol;
  aoaol.reserve(12);
  
  // creating the ArrayOfArrayOfLineRecord
  create_Zeeman_linerecordarrays(aoaol,
                                 abs_species,
                                 abs_lines_per_species,
                                 isotopologue_quantum,
                                 H_mag,
                                 1,//DO_RS,
                                 1,//DO_DF,
                                 0,//DO_QR,
                                 1,//DO_Main,
                                 1,//DO_J,
                                 0,//DO_M,
                                 verbosity);

  // holder variable
  Tensor3 part_abs_mat(f_grid.nelem(), 4, 4);
  Tensor3 part_src_mat(f_grid.nelem(), 4, 4);
  if( !do_src )
      part_src_mat.resize(0,0,0);
  Index zeeman_ind =0; // This is necessary for more than 1 Zeeman species
  
  for(Index II = 0; II<abs_species.nelem(); II++)
  {
    // If the species isn't Zeeman, look at the next species
    if(!is_zeeman(abs_species[II])) continue;
    
    // Add Pi contribution to final propmat_clearsky
    xsec_species_line_mixing_wrapper_with_zeeman( part_abs_mat, part_src_mat, abs_species, abs_lineshape,
                                                  aoaol[zeeman_ind+1], Vector(), isotopologue_ratios, abs_t_nlte,
                                                  abs_vmrs, abs_p, abs_t, f_grid, lm_p_lim, theta, eta, 0, II, verbosity );
    propmat_clearsky(II, joker, joker, joker) += part_abs_mat;
    if( do_src )
      propmat_source_clearsky(II, joker, joker, joker) += part_src_mat;

    // Add Sigma minus contribution to final propmat_clearsky
    xsec_species_line_mixing_wrapper_with_zeeman( part_abs_mat, part_src_mat, abs_species, abs_lineshape,
                                                  aoaol[zeeman_ind+0], Vector(), isotopologue_ratios, abs_t_nlte,  
                                                  abs_vmrs, abs_p, abs_t, f_grid, lm_p_lim, theta, eta, -1, II, verbosity );
    propmat_clearsky(II, joker, joker, joker) += part_abs_mat;
    if( do_src )
      propmat_source_clearsky(II, joker, joker, joker) += part_src_mat;

    // Add Sigma plus contribution to final propmat_clearsky
    xsec_species_line_mixing_wrapper_with_zeeman( part_abs_mat, part_src_mat, abs_species, abs_lineshape,
                                                  aoaol[zeeman_ind+2],Vector(), isotopologue_ratios, abs_t_nlte,
                                                  abs_vmrs, abs_p, abs_t, f_grid, lm_p_lim, theta, eta, 1, II, verbosity );
    propmat_clearsky(II, joker, joker, joker) += part_abs_mat;
    if( do_src )
      propmat_source_clearsky(II, joker, joker, joker) += part_src_mat;
    
    zeeman_ind +=3;
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void zeeman_linerecord_precalcCreateFromLines( ArrayOfArrayOfLineRecord& zeeman_linerecord_precalc,
                                               const ArrayOfArrayOfSpeciesTag& abs_species,
                                               const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                               const SpeciesAuxData& isotopologue_quantum,
                                               const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  zeeman_linerecord_precalc.resize(0);
  zeeman_linerecord_precalc.reserve(24);//will always be multiple of three, default is high

  {// Begin TEST(s)
  if (abs_species.nelem() != abs_lines_per_species.nelem())
      throw std::runtime_error("Dimension of *abs_species* and *abs_lines_per_species* don't match.");
  }// End   TEST(s)
  
  // creating the ArrayOfArrayOfLineRecord
  create_Zeeman_linerecordarrays(zeeman_linerecord_precalc,
                                 abs_species,
                                 abs_lines_per_species,
                                 isotopologue_quantum,
                                 0.,//H_mag
                                 1,//DO_RS,
                                 0,//DO_DF,
                                 1,//DO_QR,
                                 1,//DO_Main,
                                 1,//DO_J,
                                 0,//DO_M,
                                 verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddZeemanFromPreCalc(Tensor4& propmat_clearsky,
					  Tensor4& propmat_source_clearsky,
                                          const ArrayOfArrayOfLineRecord& zeeman_linerecord_precalc,
                                          const Vector& f_grid,
                                          const ArrayOfArrayOfSpeciesTag& abs_species,
                                          const ArrayOfLineshapeSpec& abs_lineshape,
                                          const SpeciesAuxData& isotopologue_ratios,
                                          const SpeciesAuxData& isotopologue_quantum,
                                          const Numeric& rtp_pressure,
                                          const Numeric& rtp_temperature,
                                          const Numeric& lm_p_lim,
                                          const Vector& rtp_temperature_nlte,
                                          const Vector& rtp_vmr,
                                          const Vector& rtp_mag,
                                          const Vector& ppath_los,
                                          const Index& atmosphere_dim,
                                          const Index& manual_zeeman_tag,
                                          const Numeric& manual_zeeman_magnetic_field_strength,
                                          const Numeric& manual_zeeman_theta,
                                          const Numeric& manual_zeeman_eta,
                                          const Verbosity& verbosity)
{
  CREATE_OUT3;

  // Check that correct isotopologue ratios are defined for the species
  // we want to calculate
  checkIsotopologueRatios(abs_species, isotopologue_ratios);

  bool do_src = true;
  {// Begin TEST(s)
  if (abs_species.nelem() == 0)
      throw std::runtime_error("No Zeeman species have been defined.");
  if( propmat_clearsky.ncols()  != 4 )
      throw std::runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
  if( propmat_clearsky.nrows()  != 4 )
      throw std::runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
  if( propmat_clearsky.npages() != f_grid.nelem() )
      throw std::runtime_error("Frequency dimension of *propmat_clearsky* not equal to length of *f_grid*.");
  if( propmat_clearsky.nbooks() != abs_species.nelem() )
      throw std::runtime_error("Species dimension of *propmat_clearsky* not equal to length of *abs_species*.");
  if( rtp_mag.nelem() != 3 )
      throw std::runtime_error("*rtp_mag* must have length 3.");
  if( atmosphere_dim != 3 )
      throw std::runtime_error("*atmosphere_dim* must be 3.  Zeeman Effect is only implemented for 3D geometry.");
  if( ppath_los.nelem() != 2 )
      throw std::runtime_error("*ppath_los* is not set correctly.");
  if( zeeman_linerecord_precalc.nelem() % 3 != 0 )
      throw std::runtime_error("Length of *zeeman_linerecord_precalc* must be multiple of 3 for polarization states.  It is not.");
  if( propmat_source_clearsky.ncols()  != 4 )
    if( propmat_source_clearsky.ncols()  != 0 ) 
      throw std::runtime_error("Zeeman Effect is only implemented for source Stokes dimension 4 or 0.");
  if( propmat_source_clearsky.nrows()  != 4 )
    if( propmat_source_clearsky.nrows()  != 0 ) 
      throw std::runtime_error("Zeeman Effect is only implemented for source Stokes dimension 4 or 0.");
  if( propmat_source_clearsky.npages() != f_grid.nelem() )
    if( propmat_source_clearsky.npages() != 0 ) 
      throw std::runtime_error("Frequency dimension of *propmat_source_clearsky* not equal to length of *f_grid* or 0.");
  if( propmat_source_clearsky.nbooks() != abs_species.nelem() )
  {
    if( propmat_source_clearsky.nbooks() != 0 )
      throw std::runtime_error("Species dimension of *propmat_clearsky* not equal to length of *abs_species* or 0.");
    else
      do_src = false; // Note that test is here beacuse it makes things easier.
  }
  }// End   TEST(s)

  Vector R_path_los;
  mirror_los(R_path_los, ppath_los, atmosphere_dim);

  // Using the standard scalar absorption functions to get physics parameters,
  Vector abs_p, abs_t; Matrix abs_vmrs, abs_t_nlte;
  AbsInputFromRteScalars( abs_p, abs_t, abs_t_nlte, abs_vmrs,           // Output
          rtp_pressure, rtp_temperature, rtp_temperature_nlte, rtp_vmr, //Input
          verbosity);                                                   // Verbose!

  // FOR LOG:  Loss of speed when mag == 0
  // Set the magnetic parameters...
  Numeric H_mag,eta,theta;
  set_magnetic_parameters(H_mag,eta,theta,manual_zeeman_tag,manual_zeeman_eta,
                          manual_zeeman_theta,manual_zeeman_magnetic_field_strength,
                          rtp_mag,R_path_los);
  
  Numeric S;
  Index hund,DM,DJ,DMain;
  Rational Main,J,M;
  Numeric GS;
  Numeric (*frequency_change)(const Rational&, const  Rational&, const Rational&, 
                                  const Numeric&, const Index&, const Index&, 
                                  const Index&, const Numeric&, const Numeric&);
  
  // Store central frequency here
  ArrayOfVector FreqShift(zeeman_linerecord_precalc.nelem());
  
  // Section to fix central line frequency
  for(Index II=0;II<zeeman_linerecord_precalc.nelem();II++)
  {
    FreqShift[II].resize(zeeman_linerecord_precalc[II].nelem());
    for(Index JJ=0;JJ<zeeman_linerecord_precalc[II].nelem();JJ++)
    {
      
      // Set necessary parameters from isotopologue_quantum
      set_part_isotopolouge_constants(hund,GS,isotopologue_quantum,zeeman_linerecord_precalc[II][JJ]);
      
      // Set quantum numbers
      set_quantum_numbers(Main,DMain,J,DJ,M,DM,S,zeeman_linerecord_precalc[II][JJ],hund,isotopologue_quantum,1,1,1);
      
      // Separate setting of the frequency_change function...
      if( hund ==0 )//Case a
          frequency_change=frequency_change_casea;
      else if( hund == 1 )// Case b
          frequency_change=frequency_change_caseb;
      else
      {
          std::ostringstream os;
          os << "There are undefined Hund cases: " << zeeman_linerecord_precalc[II][JJ] << 
          "\nThe case is: "<<hund<<", allowed are (a): "<<0<<" and (b): " << 1<<"\n";
          throw std::runtime_error(os.str());
      }
      
      // Store the frequency shift
    FreqShift[II][JJ] = frequency_change(Main, M, J, S, DJ, DM, DMain,H_mag,GS);
    }
  }
  
  
  
  
  // holder variable
  Tensor3 part_abs_mat(f_grid.nelem(), 4, 4);
  Tensor3 part_src_mat(f_grid.nelem(), 4, 4);
  if( !do_src )
      part_src_mat.resize(0,0,0);
  Index zeeman_ind =0; // This is necessary for more than 1 Zeeman species
  
  for(Index II = 0; II<abs_species.nelem(); II++)
  {
    // If the species isn't Zeeman, look at the next species
    if(!is_zeeman(abs_species[II])) continue;
    
    // Add Pi contribution to final propmat_clearsky
    xsec_species_line_mixing_wrapper_with_zeeman( part_abs_mat, part_src_mat, abs_species, abs_lineshape,
                                                  zeeman_linerecord_precalc[zeeman_ind+1], FreqShift[zeeman_ind+1], 
                                                  isotopologue_ratios, abs_t_nlte, abs_vmrs, abs_p, abs_t, f_grid, lm_p_lim,
                                                  theta, eta, 0, II, verbosity );
    propmat_clearsky(II, joker, joker, joker) += part_abs_mat;
    if( do_src )
      propmat_source_clearsky(II, joker, joker, joker) += part_src_mat;

    // Add Sigma minus contribution to final propmat_clearsky
    xsec_species_line_mixing_wrapper_with_zeeman( part_abs_mat, part_src_mat, abs_species, abs_lineshape,
                                                  zeeman_linerecord_precalc[zeeman_ind+0], FreqShift[zeeman_ind+0], 
                                                  isotopologue_ratios, abs_t_nlte, abs_vmrs, abs_p, abs_t, f_grid, lm_p_lim,
                                                  theta, eta, -1, II, verbosity );
    propmat_clearsky(II, joker, joker, joker) += part_abs_mat;
    if( do_src )
      propmat_source_clearsky(II, joker, joker, joker) += part_src_mat;

    // Add Sigma plus contribution to final propmat_clearsky
    xsec_species_line_mixing_wrapper_with_zeeman( part_abs_mat, part_src_mat, abs_species, abs_lineshape,
                                                  zeeman_linerecord_precalc[zeeman_ind+2], FreqShift[zeeman_ind+2],
                                                  isotopologue_ratios, abs_t_nlte, abs_vmrs, abs_p, abs_t, f_grid, lm_p_lim,
                                                  theta, eta, 1, II, verbosity );
    propmat_clearsky(II, joker, joker, joker) += part_abs_mat;
    if( do_src )
      propmat_source_clearsky(II, joker, joker, joker) += part_src_mat;
    
    zeeman_ind +=3;
  }
    
}