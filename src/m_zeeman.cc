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
#include "global_data.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddZeeman( Tensor4& propmat_clearsky,
                                Tensor3& nlte_source,
                                ArrayOfTensor3& dpropmat_clearsky_dx,
                                ArrayOfMatrix& dnlte_dx_source,
                                ArrayOfMatrix& nlte_dsource_dx,
                                const Vector& f_grid,
                                const ArrayOfArrayOfSpeciesTag& abs_species,
                                const ArrayOfRetrievalQuantity& jacobian_quantities,
                                const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                                const ArrayOfLineshapeSpec& abs_lineshape,
                                const SpeciesAuxData& isotopologue_ratios,
                                const SpeciesAuxData& isotopologue_quantum,
                                const SpeciesAuxData& partition_functions,
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
    ArrayOfArrayOfLineRecord zeeman_linerecord_precalc;
    zeeman_linerecord_precalcCreateFromLines( zeeman_linerecord_precalc,
                                              abs_species,
                                              abs_lines_per_species,
                                              isotopologue_quantum,
                                              verbosity);
    
    propmat_clearskyAddZeemanFromPreCalc(propmat_clearsky,
                                         nlte_source,
                                         dpropmat_clearsky_dx,
                                         dnlte_dx_source,
                                         nlte_dsource_dx,
                                         zeeman_linerecord_precalc,
                                         f_grid,
                                         abs_species,
                                         jacobian_quantities,
                                         abs_lineshape,
                                         isotopologue_ratios,
                                         isotopologue_quantum,
                                         partition_functions,
                                         rtp_pressure,
                                         rtp_temperature,
                                         lm_p_lim,
                                         rtp_temperature_nlte,
                                         rtp_vmr,
                                         rtp_mag,
                                         ppath_los,
                                         atmosphere_dim,
                                         manual_zeeman_tag,
                                         manual_zeeman_magnetic_field_strength,
                                         manual_zeeman_theta,
                                         manual_zeeman_eta,
                                         verbosity);
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
                                 verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddZeemanFromPreCalc(Tensor4& propmat_clearsky,
                                          Tensor3& nlte_source,
                                          ArrayOfTensor3& dpropmat_clearsky_dx,
                                          ArrayOfMatrix& dnlte_dx_source,
                                          ArrayOfMatrix& nlte_dsource_dx,
                                          const ArrayOfArrayOfLineRecord& zeeman_linerecord_precalc,
                                          const Vector& f_grid,
                                          const ArrayOfArrayOfSpeciesTag& abs_species,
                                          const ArrayOfRetrievalQuantity& jacobian_quantities,
                                          const ArrayOfLineshapeSpec& abs_lineshape,
                                          const SpeciesAuxData& isotopologue_ratios,
                                          const SpeciesAuxData& isotopologue_quantum,
                                          const SpeciesAuxData& partition_functions,
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

  bool do_src = !nlte_source.empty();
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
    if(do_src)
    {
        if(nlte_source.npages()!=propmat_clearsky.nbooks())
        throw std::runtime_error("Species dimension of *nlte_source* not equal to length of *abs_species*.");
        if(nlte_source.npages()!=propmat_clearsky.npages())
        throw std::runtime_error("Frequency dimension of *nlte_source* not equal to length of *f_grid*.");
        if(nlte_source.npages()!=propmat_clearsky.nrows())
        throw std::runtime_error("Zeeman Effect is only implemented for Stokes dimension 4.");
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
  
  Numeric GS;
  
  // Store central frequency here
  ArrayOfVector FreqShift(zeeman_linerecord_precalc.nelem());
  
  // Section to fix central line frequency
  for(Index II=0;II<zeeman_linerecord_precalc.nelem();II++)
  {
    FreqShift[II].resize(zeeman_linerecord_precalc[II].nelem());
    for(Index JJ=0;JJ<zeeman_linerecord_precalc[II].nelem();JJ++)
    {
      
      // Set necessary parameters from isotopologue_quantum
      set_GS(GS, isotopologue_quantum, zeeman_linerecord_precalc[II][JJ]);
      
      // Store the frequency shift
      FreqShift[II][JJ] = frequency_change(zeeman_linerecord_precalc[II][JJ], H_mag, GS);
    }
  }
  
  const PropmatPartialsData pps(jacobian_quantities);
  pps.supportsZeemanPrecalc();
  const bool do_freq_jac = pps.do_frequency();
  const bool do_temp_jac = pps.do_temperature();
  
  Vector planck_BT(0);
  if(do_src)
      planck_BT.resize(f_grid.nelem());
  
  Matrix dplanck_BT(0,0);
  
  if(do_freq_jac||do_temp_jac)
      dplanck_BT.resize(2,f_grid.nelem());
  
  for(Index iv=0;iv<planck_BT.nelem();iv++)
  {
      planck_BT[iv] = planck(f_grid[iv],abs_t[0]);
      if(dplanck_BT.ncols()>0)
      {
          dplanck_BT(0,iv) = dplanck_dt(f_grid[iv],abs_t[0]);
          dplanck_BT(1,iv) = dplanck_df(f_grid[iv],abs_t[0]);
      }
  }
  
  Index zeeman_ind =0; // This is necessary for more than 1 Zeeman species
  
  for(Index II = 0; II<abs_species.nelem(); II++)
  {
    const Index ls_index = (1==abs_lineshape.nelem())?0:II;
    
    // If the species isn't Zeeman, look at the next species
    if(!is_zeeman(abs_species[II])) continue;
    
    // Add Pi contribution to final propmat_clearsky
    xsec_species_line_mixing_wrapper_with_zeeman( propmat_clearsky, nlte_source, dpropmat_clearsky_dx, dnlte_dx_source, nlte_dsource_dx, 
                                                  abs_species, pps, 
                                                  abs_lineshape[ls_index].Ind_ls(), abs_lineshape[ls_index].Ind_lsn(), abs_lineshape[ls_index].Cutoff(), 
                                                  zeeman_linerecord_precalc[zeeman_ind+1], FreqShift[zeeman_ind+1], planck_BT, dplanck_BT,
                                                  isotopologue_ratios, partition_functions, abs_t_nlte, abs_vmrs, abs_p, abs_t, f_grid, 
                                                  rtp_mag, R_path_los,lm_p_lim,theta, eta, H_mag, 0, II, verbosity );

    // Add Sigma minus contribution to final propmat_clearsky
    xsec_species_line_mixing_wrapper_with_zeeman( propmat_clearsky, nlte_source, dpropmat_clearsky_dx, dnlte_dx_source, nlte_dsource_dx, 
                                                  abs_species, pps, 
                                                  abs_lineshape[ls_index].Ind_ls(), abs_lineshape[ls_index].Ind_lsn(), abs_lineshape[ls_index].Cutoff(), 
                                                  zeeman_linerecord_precalc[zeeman_ind], FreqShift[zeeman_ind], planck_BT, dplanck_BT,
                                                  isotopologue_ratios, partition_functions, abs_t_nlte, abs_vmrs, abs_p, abs_t, f_grid, 
                                                  rtp_mag, R_path_los,lm_p_lim,theta, eta, H_mag, -1, II, verbosity );

    // Add Sigma plus contribution to final propmat_clearsky
    xsec_species_line_mixing_wrapper_with_zeeman( propmat_clearsky, nlte_source, dpropmat_clearsky_dx, dnlte_dx_source, nlte_dsource_dx, 
                                                  abs_species, pps, 
                                                  abs_lineshape[ls_index].Ind_ls(), abs_lineshape[ls_index].Ind_lsn(), abs_lineshape[ls_index].Cutoff(), 
                                                  zeeman_linerecord_precalc[zeeman_ind+2], FreqShift[zeeman_ind+2], planck_BT, dplanck_BT,
                                                  isotopologue_ratios, partition_functions, abs_t_nlte, abs_vmrs, abs_p, abs_t, f_grid, 
                                                  rtp_mag, R_path_los,lm_p_lim,theta, eta, H_mag, 1, II, verbosity );
    
    // The flat structure reminder for 3-component ArrayOfArrayOfLineRecord
    
    zeeman_ind +=3;
  }
    
}