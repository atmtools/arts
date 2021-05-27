/* Copyright (C) 2013
   Patrick Eriksson <patrick.eriksson@chalmers.se>
                            
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
  === File description 
  ===========================================================================*/

/*!
  \file   m_fos.cc
  \author Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2013-04-29

  \brief  Workspace functions associated with the FOS scattering scheme.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "doit.h"
#include "m_xml.h"
#include "math_funcs.h"
#include "montecarlo.h"
#include "rte.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;

// FOS implemented as an internal function, to allow an recursive algorithm
/*
void fos(
         Workspace&                     ws,
         Matrix&                        iy,
         ArrayOfTensor4&                iy_aux,
         Ppath&                         ppath,
         ArrayOfTensor3&                diy_dx,
   const Index&                         stokes_dim,
   const Vector&                        f_grid,
   const Index&                         atmosphere_dim,
   const Vector&                        p_grid,
   const Tensor3&                       z_field,
   const Tensor3&                       t_field,
   const Tensor4&                       vmr_field,
   const ArrayOfArrayOfSpeciesTag&      abs_species,
   const Tensor3&                       wind_u_field,
   const Tensor3&                       wind_v_field,
   const Tensor3&                       wind_w_field,
   const Tensor3&                       mag_u_field,
   const Tensor3&                       mag_v_field,
   const Tensor3&                       mag_w_field,
   const Index&                         cloudbox_on,
   const ArrayOfIndex&                  cloudbox_limits,
   const Tensor4&                       pnd_field,
   const ArrayOfArrayOfSingleScatteringData&   scat_data,
   const Matrix&                        particle_masses,
   const String&                        iy_unit,
   const ArrayOfString&                 iy_aux_vars,
   const Index&                         jacobian_do,
   const Agenda&                        ppath_agenda,
   const Agenda&                        propmat_clearsky_agenda,
   const Agenda&                        iy_main_agenda,
   const Agenda&                        iy_space_agenda,
   const Agenda&                        iy_surface_agenda,
   const Index&                         iy_agenda_call1,
   const Tensor3&                       iy_transmittance,
   const Vector&                        rte_pos,      
   const Vector&                        rte_los,      
   const Vector&                        rte_pos2, 
   const Numeric&                       rte_alonglos_v,      
   const Numeric&                       ppath_lmax,
   const Numeric&                       ppath_lraytrace,
   const Matrix&                        fos_scatint_angles,
   const Vector&                        fos_iyin_za_angles,
   const Index&                         fos_za_interporder,
   const Index&                         fos_n,
   const Index&                         fos_i,
   const Verbosity&                     verbosity )
{
  // A temporary restriction
  if( atmosphere_dim > 1 )
    throw runtime_error( "FOS is so far only handling 1D atmospheres." );

  ARTS_ASSERT( fos_i >= 0  &&  fos_i <= fos_n );


  // Determine propagation path
  //
  ppath_agendaExecute( ws, ppath, ppath_lmax, ppath_lraytrace,
                       rte_pos, rte_los, rte_pos2, 
                       0, 0, t_field, z_field, vmr_field, f_grid, 
                       ppath_agenda );

  // Some basic sizes
  //
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;

  // The below copied from iyEmission. Activate for-loop when jacobians
  // introduced.

  // Set up variable with index of species where we want abs_per_species.
  // This variable can below be extended in iy_aux part.
  //
  ArrayOfIndex iaps(0);
  //
  //for( Index i=0; i<jac_species_i.nelem(); i++ )
  //  {
  //    if( jac_species_i[i] >= 0 )
  //      { iaps.push_back( jac_species_i[i] ); }
  //  }
  
  //=== iy_aux part ===========================================================
  Index auxPressure    = -1,
        auxTemperature = -1,
        auxAbsSum      = -1,
        auxBackground  = -1,
        auxIy          = -1,
        auxOptDepth    = -1;
  ArrayOfIndex auxAbsSpecies(0), auxAbsIsp(0);
  ArrayOfIndex auxVmrSpecies(0), auxVmrIsp(0);
  ArrayOfIndex auxPartCont(0),   auxPartContI(0);
  ArrayOfIndex auxPartField(0),  auxPartFieldI(0);
  //
  if( !iy_agenda_call1 )
    { iy_aux.resize( 0 ); }
  else
    {
      const Index naux = iy_aux_vars.nelem();
      iy_aux.resize( naux );
      //
      for( Index i=0; i<naux; i++ )
        {
          if( iy_aux_vars[i] == "Pressure" )
            { auxPressure = i;      iy_aux[i].resize( 1, 1, 1, np ); }
          else if( iy_aux_vars[i] == "Temperature" )
            { auxTemperature = i;   iy_aux[i].resize( 1, 1, 1, np ); }
          else if( iy_aux_vars[i].substr(0,13) == "VMR, species " )
            { 
              Index ispecies;
              istringstream is(iy_aux_vars[i].substr(13,2));
              is >> ispecies;
              if( ispecies < 0  ||  ispecies>=abs_species.nelem() )
                {
                  ostringstream os;
                  os << "You have selected VMR of species with index "
                     << ispecies << ".\nThis species does not exist!";
                  throw runtime_error( os.str() );
                }
              auxVmrSpecies.push_back(i);
              auxVmrIsp.push_back(ispecies);
              iy_aux[i].resize( 1, 1, 1, np );               
            }
          else if( iy_aux_vars[i] == "Absorption, summed" )
            { auxAbsSum = i;   iy_aux[i].resize( nf, ns, ns, np ); }
          else if( iy_aux_vars[i].substr(0,20) == "Absorption, species " )
            { 
              Index ispecies;
              istringstream is(iy_aux_vars[i].substr(20,2));
              is >> ispecies;
              if( ispecies < 0  ||  ispecies>=abs_species.nelem() )
                {
                  ostringstream os;
                  os << "You have selected absorption species with index "
                     << ispecies << ".\nThis species does not exist!";
                  throw runtime_error( os.str() );
                }
              auxAbsSpecies.push_back(i);
              const Index ihit = find_first( iaps, ispecies );
              if( ihit >= 0 )
                { auxAbsIsp.push_back( ihit ); }
              else
                { 
                  iaps.push_back(ispecies); 
                  auxAbsIsp.push_back( iaps.nelem()-1 ); 
                }
              iy_aux[i].resize( nf, ns, ns, np );               
            }
          else if( iy_aux_vars[i] == "Radiative background" )
            { auxBackground = i;   iy_aux[i].resize( nf, 1, 1, 1 ); }
          else if( iy_aux_vars[i] == "iy"   &&  auxIy < 0 )
            { auxIy = i;           iy_aux[i].resize( nf, ns, 1, np ); }
          else if( iy_aux_vars[i] == "Optical depth" )
            { auxOptDepth = i;     iy_aux[i].resize( nf, 1, 1, 1 ); }
          else if( iy_aux_vars[i].substr(0,14) == "Mass content, " )
            { 
              Index icont;
              istringstream is(iy_aux_vars[i].substr(14,2));
              is >> icont;
              if( icont < 0  ||  icont>=particle_masses.ncols() )
                {
                  ostringstream os;
                  os << "You have selected particle mass content category with "
                     << "index " << icont << ".\nThis category is not defined!";
                  throw runtime_error( os.str() );
                }
              auxPartCont.push_back(i);
              auxPartContI.push_back(icont);
              iy_aux[i].resize( 1, 1, 1, np );
            }
          else if( iy_aux_vars[i].substr(0,10) == "PND, type " )
            { 
              Index ip;
              istringstream is(iy_aux_vars[i].substr(10,2));
              is >> ip;
              if( ip < 0  ||  ip>=pnd_field.nbooks() )
                {
                  ostringstream os;
                  os << "You have selected particle number density field with "
                     << "index " << ip << ".\nThis field is not defined!";
                  throw runtime_error( os.str() );
                }
              auxPartField.push_back(i);
              auxPartFieldI.push_back(ip);
              iy_aux[i].resize( 1, 1, 1, np );
            }
          else
            {
              ostringstream os;
              os << "In *iy_aux_vars* you have included: \"" << iy_aux_vars[i]
                 << "\"\nThis choice is not recognised.";
              throw runtime_error( os.str() );
            }
        }
    }
  //===========================================================================

  // Get atmospheric and RT quantities for each ppath point/step
  //
  Vector       ppath_p, ppath_t;
  Matrix       ppath_vmr, ppath_pnd, ppath_wind, ppath_mag, ppath_f, ppath_nlte;
  Matrix       ppath_blackrad;
  ArrayOfArrayOfPropagationMatrix      abs_per_species, dummy_dppath_ext_dx;
  ArrayOfPropagationMatrix      ppath_ext, pnd_ext_mat;
  Tensor4 trans_partial, trans_cumulat;
  ArrayOfArrayOfStokesVector dummy_dppath_nlte_dx;
  Tensor3      pnd_abs_vec;
  ArrayOfStokesVector ppath_nlte_source;
  Vector       scalar_tau;
  ArrayOfIndex   clear2cloudy, lte;
  ArrayOfMatrix   dummy_ppath_dpnd_dx;
  ArrayOfTensor4  dummy_dpnd_field_dx;
  const Tensor4   nlte_field_empty(0,0,0,0);
  //
  Array<ArrayOfArrayOfSingleScatteringData> scat_data_single;
  ArrayOfArrayOfIndex                extmat_case;  
  //
  if( np > 1 )
    {
      get_ppath_atmvars(  ppath_p, ppath_t, ppath_nlte, ppath_vmr,
                          ppath_wind, ppath_mag, 
                          ppath, atmosphere_dim, p_grid, t_field,
                          nlte_field_empty, vmr_field,
                          wind_u_field, wind_v_field, wind_w_field,
                          mag_u_field, mag_v_field, mag_w_field );
      
      get_ppath_f( ppath_f, ppath, f_grid,  atmosphere_dim, 
                   rte_alonglos_v, ppath_wind );
      
      get_ppath_pmat( ws, ppath_ext, ppath_nlte_source, lte, abs_per_species, 
                      dummy_dppath_ext_dx, dummy_dppath_nlte_dx,
                      propmat_clearsky_agenda, ArrayOfRetrievalQuantity(0), ppath, 
                      ppath_p, ppath_t, ppath_nlte, ppath_vmr, ppath_f, 
                      ppath_mag, f_grid, stokes_dim, iaps );
      
      for( Index i=0; i<lte.nelem(); i++ )
        {
          if( lte[i] == 0 )
            throw runtime_error( "FOS can so far only handle LTE conditions." );
        }
      
      get_ppath_blackrad( ppath_blackrad, ppath, ppath_t, ppath_f );
      
      if( !cloudbox_on )
        { 
          get_ppath_trans( trans_partial, extmat_case, trans_cumulat, 
                           scalar_tau, ppath, ppath_ext, f_grid, stokes_dim );
        }
      else
        {
          get_ppath_cloudvars( clear2cloudy, ppath_pnd, dummy_ppath_dpnd_dx,
                               ppath, atmosphere_dim, cloudbox_limits,
                               pnd_field, dummy_dpnd_field_dx );
          get_ppath_partopt( pnd_abs_vec, pnd_ext_mat, scat_data_single,
                             clear2cloudy, ppath_pnd, ppath, ppath_t,
                             stokes_dim, ppath_f, atmosphere_dim, scat_data,
                             verbosity );
          get_ppath_trans2( trans_partial, extmat_case, trans_cumulat, 
                            scalar_tau, ppath, ppath_ext, f_grid, stokes_dim, 
                            clear2cloudy, pnd_ext_mat );
        }      
    }
  else // For cases totally outside the atmosphere,
    {  // set zero optical thickness and unit transmission
      scalar_tau.resize( nf );
      scalar_tau = 0;
      trans_cumulat.resize( nf, ns, ns, np );
      for( Index iv=0; iv<nf; iv++ )
        { id_mat( trans_cumulat(iv,joker,joker,np-1) ); }
    }

  // iy_transmittance
  //
  Tensor3 iy_trans_new;
  //
  if( iy_agenda_call1 )
    { iy_trans_new = trans_cumulat(joker,joker,joker,np-1); }
  else
    { 
      iy_transmittance_mult( iy_trans_new, iy_transmittance, 
                            trans_cumulat(joker,joker,joker,np-1) ); }

  // Radiative background
  //
  {
    Agenda iy_cbox_agenda;
    const Index iy_id = 0;
    get_iy_of_background( ws, iy, diy_dx, 
                          iy_trans_new, iy_id,  jacobian_do, ppath, rte_pos2, 
                          atmosphere_dim, t_field, z_field, vmr_field, 
                          cloudbox_on, stokes_dim, f_grid, iy_unit,
                          iy_main_agenda, iy_space_agenda, iy_surface_agenda, 
                          iy_cbox_agenda, verbosity );
  }

  //=== iy_aux part ===========================================================
  // Fill parts of iy_aux that are defined even for np=1.
  // Radiative background
  if( auxBackground >= 0 ) 
    { iy_aux[auxBackground](0,0,0,0) = (Numeric)min( (Index)2,
                                              ppath_what_background(ppath)-1); }
  // Radiance 
  if( auxIy >= 0 ) 
    { iy_aux[auxIy](joker,joker,0,np-1) = iy; }
  // Transmission, total
  if( auxOptDepth >= 0 ) 
    { iy_aux[auxOptDepth](joker,0,0,0) = scalar_tau; }
  //===========================================================================

  
  // Do RT calculations
  //
  if( np > 1 )
    {
      //=== iy_aux part =======================================================
      // iy_aux for point np-1:
      // Pressure
      if( auxPressure >= 0 ) 
        { iy_aux[auxPressure](0,0,0,np-1) = ppath_p[np-1]; }
      // Temperature
      if( auxTemperature >= 0 ) 
        { iy_aux[auxTemperature](0,0,0,np-1) = ppath_t[np-1]; }
      // VMR
      for( Index j=0; j<auxVmrSpecies.nelem(); j++ )
        { iy_aux[auxVmrSpecies[j]](0,0,0,np-1) = ppath_vmr(auxVmrIsp[j],np-1); }
      // Absorption
      if( auxAbsSum >= 0 ) 
        { for( Index iv=0; iv<nf; iv++ ) {
            for( Index is1=0; is1<ns; is1++ ){
              for( Index is2=0; is2<ns; is2++ ){
                iy_aux[auxAbsSum](iv,is1,is2,np-1) = 
                                             ppath_ext[np-1](iv,is1,is2); } } } }
      for( Index j=0; j<auxAbsSpecies.nelem(); j++ )
        { for( Index iv=0; iv<nf; iv++ ) {
            for( Index is1=0; is1<stokes_dim; is1++ ){
              for( Index is2=0; is2<stokes_dim; is2++ ){
                iy_aux[auxAbsSpecies[j]](iv,is1,is2,np-1) = 
                          abs_per_species[np-1][auxAbsIsp[j]](iv,is1,is2); } } } }
      // Particle properties
      if( cloudbox_on  )
        {
          // Particle mass content
          for( Index j=0; j<auxPartCont.nelem(); j++ )
            { iy_aux[auxPartCont[j]](0,0,0,np-1) = ppath_pnd(joker,np-1) *
                                      particle_masses(joker,auxPartContI[j]); }
          // Particle number density
          for( Index j=0; j<auxPartField.nelem(); j++ )
            { iy_aux[auxPartField[j]](0,0,0,np-1) = 
                                            ppath_pnd(auxPartFieldI[j],np-1); }
        }
      // Radiance for this point is handled above
      //=======================================================================


      // Scattering source term at ip (0) and ip+1 (1)
      // (If any particles at ip=np-1, this is ignored. Could happen if
      // particles at surface level.)
      Matrix ssource0(nf,ns,0), ssource1(nf,ns);

      // Dummy variables for non-LTE
      const bool nonlte = false;
      const Matrix sourcebar_dummy(0,0);
      const Tensor3 extbar_dummy(0,0,0);

      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          // Path step average of emission source function: Bbar
          Vector bbar(nf);
          for( Index iv=0; iv<nf; iv++ )  
            { bbar[iv] = 0.5 * ( ppath_blackrad(iv,ip) +
                                 ppath_blackrad(iv,ip+1) ); }

          // Check if any particles to consider
          bool any_particles = clear2cloudy[ip] >= 0 || 
                               clear2cloudy[ip+1] >= 0;

          // -----------------------------------------------------------------
          // i = N (only absorption/emission)
          // -----------------------------------------------------------------
          if( fos_i == fos_n ) 
            {
              // No particle absorption to consider
              if( !any_particles )
                {
                  emission_rtstep( iy, stokes_dim, bbar, extmat_case[ip],
                                   trans_partial(joker,joker,joker,ip),
                                   nonlte, extbar_dummy, sourcebar_dummy );
                }

              else  // We want to include particle absorption, but not
                {   // extinction. trans_partial is then not valid.
                  Tensor3      t(nf,ns,ns);
                  ArrayOfIndex extmat_cas2(nf);
                  //
                  for( Index iv=0; iv<nf; iv++ )
                    {
                      // Particle absorption
                      //
                      Matrix pabs_mat(ns,ns,0);
                      //
                      if( clear2cloudy[ip] >= 0 )
                        { ext_matFromabs_vec( pabs_mat, pnd_abs_vec(iv,joker,
                                          clear2cloudy[ip]), stokes_dim ); } 
                      if( clear2cloudy[ip+1] >= 0 )
                        { ext_matFromabs_vec( pabs_mat, pnd_abs_vec(iv,joker,
                                        clear2cloudy[ip+1]), stokes_dim ); } 

                      // Transmission of step
                      Matrix ext_mat(stokes_dim,stokes_dim);  
                      for( Index is1=0; is1<stokes_dim; is1++ ) {
                        for( Index is2=0; is2<stokes_dim; is2++ ) {
                          ext_mat(is1,is2) = 0.5 * ( pabs_mat(is1,is2) +
                                                  ppath_ext[ip](iv,is1,is2) +
                                                  ppath_ext[ip+1](iv,is1,is2) );
                        } }
                      //
                      extmat_cas2[iv] = 0;
                      ext2trans( t(iv,joker,joker), extmat_cas2[iv], 
                                 ext_mat, ppath.lstep[ip] );
                    }
                                
                  // Perform RT
                  emission_rtstep( iy, stokes_dim, bbar, extmat_cas2, t,
                                   nonlte, extbar_dummy, sourcebar_dummy );
                }
            }
          

          // -----------------------------------------------------------------
          // i < N
          // -----------------------------------------------------------------
          else
            {
              // Shift scattering source term (new 1 is old 0)
              ssource1 = ssource0;

              // Clear-sky path step
              if( !any_particles )
                {
                  // Perform RT
                  emission_rtstep( iy, stokes_dim, bbar, extmat_case[ip],
                                   trans_partial(joker,joker,joker,ip),
                                   nonlte, extbar_dummy, sourcebar_dummy );

                  // Scattering source term at ip is zero:
                  ssource0 = 0;
                }

              // Include scattering
              else
                {
                  // Determine scattering source term at ip
                  if( clear2cloudy[ip] < 0 )
                    { ssource0 = 0; }
                  else
                    {
                      // Present position 
                      // (Note that the Ppath positions (ppath.pos) for 1D have
                      // one column more than expected by most functions. Only 
                      // the first atmosphere_dim values shall be copied.)
                      Vector pos = ppath.pos(ip,Range(0,atmosphere_dim));

                      // Determine incoming radiation
                      //
                      const Index nin = fos_scatint_angles.nrows();
                      Tensor3 Y(nin,nf,ns);
                      {
                        // Do RT
                        const Index  nza = fos_iyin_za_angles.nelem();
                        Tensor3 Y1(nza,nf,ns);
                        //
                        for( Index i=0; i<nza; i++ )
                          {
                            // LOS
                            Vector los( 1, fos_iyin_za_angles[i] );

                            // Call recursively, with fos_i increased
                            // 
                            Matrix           iyl;
                            ArrayOfTensor4   iy_auxl;
                            Ppath            ppathl;
                            ArrayOfTensor3   diy_dxl;
                            //
                            fos( ws, iyl, iy_auxl, ppathl, diy_dxl, 
                                 stokes_dim, f_grid, atmosphere_dim,
                                 p_grid, z_field, t_field, vmr_field, 
                                 abs_species, wind_u_field, wind_v_field, 
                                 wind_w_field, mag_u_field, mag_v_field, 
                                 mag_w_field, cloudbox_on, cloudbox_limits, 
                                 pnd_field,
                                 scat_data, 
                                 particle_masses, iy_unit, iy_aux_vars, 
                                 jacobian_do, ppath_agenda, 
                                 propmat_clearsky_agenda, iy_main_agenda, 
                                 iy_space_agenda, iy_surface_agenda, 0,
                                 iy_trans_new, pos, los, rte_pos2, 
                                 rte_alonglos_v, ppath_lmax, ppath_lraytrace, 
                                 fos_scatint_angles, fos_iyin_za_angles, 
                                 fos_za_interporder, fos_n, fos_i+1,
                                 verbosity );
                            
                            Y1(i,joker,joker) = iyl;
                          }

                        // Zenith angle interpolation of Y
                        //ArrayOfGridPosPoly gp( nin );  FIXME: REMOVE WHEN UNCOMMENTING THIS CODE
                        //gridpos_poly( gp, fos_iyin_za_angles, 
                        //              fos_scatint_angles(joker,0),
                        //              fos_za_interporder );
                        const auto lag = Interpolation::LagrangeVector(fos_scatint_angles(joker, 0), fos_iyin_za_angles, 0.5, false, Interpolation::LagrangeType::Linear);
                        //Matrix itw( nin, fos_za_interporder+1 );
                        //interpweights( itw, gp );
                        const auto itw = interpweights( lag )
                        //
                        for( Index iv=0; iv<nf; iv++ ) 
                          { 
                            for( Index is1=0; is1<stokes_dim; is1++ ) 
                              { 
                                //interp( Y(joker,iv,is1), itw, 
                                //        Y1(joker,iv,is1), gp );
                                reinterp( Y(joker,iv,is1), Y1(joker,iv,is1), itw, lag );
                              }
                          }
                      }

                      // Direction of outgoing scattered radiation (which is
                      // reversed to LOS). Note that this outlos is only used
                      // for extracting scattering properties.
                      Vector outlos;
                      mirror_los( outlos, ppath.los(ip,joker), atmosphere_dim );

                      // Determine phase matrix 
                      Tensor4  P( nin, nf, stokes_dim, stokes_dim );
                      Matrix   P1( stokes_dim, stokes_dim );
                      //
                      for( Index ii=0; ii<nin; ii++ )
                        {
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              pha_mat_singleCalc( P1, outlos[0], outlos[1], 
                                                  fos_scatint_angles(ii,0), 
                                                  fos_scatint_angles(ii,1), 
                                                  scat_data_single[iv], stokes_dim, 
                                                  ppath_pnd(joker,ip),
                                                  ppath_t[ip], verbosity );
                              P(ii,iv,joker,joker) = P1;
                            }
                        }


                      // Scattering source term
                      ssource0 = 0.0;
                      for( Index iv=0; iv<nf; iv++ )
                        { 
                          Vector sp(stokes_dim);
                          for( Index ii=0; ii<nin; ii++ )
                            { 
                              mult( sp, P(ii,iv,joker,joker), 
                                        Y(ii,iv,joker) );
                              ssource0(iv,joker) += sp;
                            }
                        }
                      ssource0 *= 4*PI/(Numeric)nin;
                    }
                 
                  // RT of ppath step 
                  for( Index iv=0; iv<nf; iv++ )
                    {
                      // Calculate average of absorption, extinction etc.
                      Matrix  ext_mat( stokes_dim, stokes_dim );
                      Vector  abs_vec( stokes_dim );
                      Vector  sbar( stokes_dim, 0 );
                      //
                      // Contribution from abs_species
                      for( Index is1=0; is1<stokes_dim; is1++ )
                        { 
                          abs_vec[is1] = 0.5 * (                
                                                    ppath_ext[ip](iv,is1,0) +
                                                    ppath_ext[ip+1](iv,is1,0) );
                          for( Index is2=0; is2<stokes_dim; is2++ )
                            {
                              ext_mat(is1,is2) = 0.5 * (
                                                  ppath_ext[ip](iv,is1,is2) +
                                                  ppath_ext[ip+1](iv,is1,is2) );
                            }
                        }
                      // Particle contribution
                      if( clear2cloudy[ip] >= 0 )
                        {
                          for( Index is1=0; is1<stokes_dim; is1++ )
                            { 
                              sbar[is1]    += 0.5 * ssource0(iv,is1); 
                              abs_vec[is1] += 0.5 * (
                                     pnd_abs_vec(iv,is1,clear2cloudy[ip]) );
                              for( Index is2=0; is2<stokes_dim; is2++ )
                                {
                                  ext_mat(is1,is2) += 0.5 * (
                                    pnd_ext_mat[clear2cloudy[ip]](iv,is1,is2) );
                                }
                            }
                        }
                      if( clear2cloudy[ip+1] >= 0 )
                        {
                          for( Index is1=0; is1<stokes_dim; is1++ )
                            { 
                              sbar[is1]    += 0.5 * ssource1(iv,is1); 
                              abs_vec[is1] += 0.5 * (
                                    pnd_abs_vec(iv,is1,clear2cloudy[ip+1]) );
                              for( Index is2=0; is2<stokes_dim; is2++ )
                                {
                                  ext_mat(is1,is2) += 0.5 * (
                                    pnd_ext_mat[clear2cloudy[ip+1]](iv,is1,is2));
                                }
                            }
                        }

                      // Perform RT
                      //
                      Matrix trans_mat = trans_partial(iv,joker,joker,ip);
                      rte_step_doit( iy(iv,joker), trans_mat, ext_mat, abs_vec,
                                     sbar, ppath.lstep[ip], bbar[iv], true );
                    }
                }
            }

          //=== iy_aux part ===================================================
          // Pressure
          if( auxPressure >= 0 ) 
            { iy_aux[auxPressure](0,0,0,ip) = ppath_p[ip]; }
          // Temperature
          if( auxTemperature >= 0 ) 
            { iy_aux[auxTemperature](0,0,0,ip) = ppath_t[ip]; }
          // VMR
          for( Index j=0; j<auxVmrSpecies.nelem(); j++ )
            { iy_aux[auxVmrSpecies[j]](0,0,0,ip) = ppath_vmr(auxVmrIsp[j],ip);}
          // Absorption
          if( auxAbsSum >= 0 ) 
            { for( Index iv=0; iv<nf; iv++ ) {
                for( Index is1=0; is1<ns; is1++ ){
                  for( Index is2=0; is2<ns; is2++ ){
                    iy_aux[auxAbsSum](iv,is1,is2,ip) = 
                                              ppath_ext[ip](iv,is1,is2); } } } }
          for( Index j=0; j<auxAbsSpecies.nelem(); j++ )
            { for( Index iv=0; iv<nf; iv++ ) {
                for( Index is1=0; is1<stokes_dim; is1++ ){
                  for( Index is2=0; is2<stokes_dim; is2++ ){
                    iy_aux[auxAbsSpecies[j]](iv,is1,is2,ip) = 
                           abs_per_species[ip][auxAbsIsp[j]](iv,is1,is2); } } } }
          // Particle properties
          if( cloudbox_on ) 
            {
              // Particle mass content
              for( Index j=0; j<auxPartCont.nelem(); j++ )
                { iy_aux[auxPartCont[j]](0,0,0,ip) = ppath_pnd(joker,ip) *
                                      particle_masses(joker,auxPartContI[j]); }
              // Particle number density
              for( Index j=0; j<auxPartField.nelem(); j++ )
                { iy_aux[auxPartField[j]](0,0,0,ip) = 
                                              ppath_pnd(auxPartFieldI[j],ip); }
            }
          // Radiance 
          if( auxIy >= 0 ) 
            { iy_aux[auxIy](joker,joker,0,ip) = iy; }
          //===================================================================
        } 
    } // if np>1


  // Unit conversions
  if( iy_agenda_call1 )
    {
      // Determine refractive index to use for the n2 radiance law
      Numeric n = 1.0; // First guess is that sensor is in space
      //
      if( ppath.end_lstep == 0 ) // If true, sensor inside the atmosphere
        { n = ppath.nreal[np-1]; }

      // Polarisation index variable
      ArrayOfIndex i_pol(stokes_dim);
      for( Index is=0; is<stokes_dim; is++ )
        { i_pol[is] = is + 1; }

      // iy
      apply_iy_unit( iy, iy_unit, f_grid, n, i_pol );

      // iy_aux
      for( Index q=0; q<iy_aux.nelem(); q++ )
        {
          if( iy_aux_vars[q] == "iy")
            { 
              for( Index ip=0; ip<ppath.np; ip++ )
                { apply_iy_unit( iy_aux[q](joker,joker,0,ip), iy_unit, f_grid, 
                                 ppath.nreal[ip], i_pol ); }
            }
        }
    }
}
*/

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
/*
void iyFOS(
         Workspace&                   ws,
         Matrix&                      iy,
         ArrayOfTensor4&              iy_aux,
         Ppath&                       ppath,
         ArrayOfTensor3&              diy_dx,
   const Index&                       stokes_dim,
   const Vector&                      f_grid,
   const Index&                       atmosphere_dim,
   const Vector&                      p_grid,
   const Tensor3&                     z_field,
   const Tensor3&                     t_field,
   const Tensor4&                     vmr_field,
   const ArrayOfArrayOfSpeciesTag&    abs_species,
   const Tensor3&                     wind_u_field,
   const Tensor3&                     wind_v_field,
   const Tensor3&                     wind_w_field,
   const Tensor3&                     mag_u_field,
   const Tensor3&                     mag_v_field,
   const Tensor3&                     mag_w_field,
   const Index&                       cloudbox_on,
   const ArrayOfIndex&                cloudbox_limits,
   const Tensor4&                     pnd_field,
   const ArrayOfArrayOfSingleScatteringData& scat_data,
   const Matrix&                      particle_masses,
   const String&                      iy_unit,
   const ArrayOfString&               iy_aux_vars,
   const Index&                       jacobian_do,
   const Agenda&                      ppath_agenda,
   const Agenda&                      propmat_clearsky_agenda,
   const Agenda&                      iy_main_agenda,
   const Agenda&                      iy_space_agenda,
   const Agenda&                      iy_surface_agenda,
   const Index&                       iy_agenda_call1,
   const Tensor3&                     iy_transmittance,
   const Vector&                      rte_pos,      
   const Vector&                      rte_los,      
   const Vector&                      rte_pos2, 
   const Numeric&                     rte_alonglos_v,      
   const Numeric&                     ppath_lmax,
   const Numeric&                     ppath_lraytrace,
   const Matrix&                      fos_scatint_angles,
   const Vector&                      fos_iyin_za_angles,
   const Index&                       fos_za_interporder,
   const Index&                       fos_n,
   const Verbosity&                   verbosity )
{
  // Input checks
  if( jacobian_do )
    throw runtime_error( 
     "This method does not yet provide any jacobians (jacobian_do must be 0)" );
  if( fos_scatint_angles.ncols() != 2 )
    throw runtime_error( "The WSV *fos_scatint_angles* must have two columns." );
  if( min(fos_scatint_angles(joker,0))<0 || 
      max(fos_scatint_angles(joker,0))>180 )
    throw runtime_error( 
          "The zenith angles in *fos_scatint_angles* shall be inside [0,180]." );
  if( min(fos_scatint_angles(joker,1))<-180 || 
      max(fos_scatint_angles(joker,1))>180 )
    throw runtime_error( 
      "The azimuth angles in *fos_scatint_angles* shall be inside [-180,180]." );
  if( min(fos_iyin_za_angles)<0 || max(fos_iyin_za_angles)>180 )
    throw runtime_error( 
          "The zenith angles in *fos_iyin_za_angles* shall be inside [0,180]." );
  if( fos_iyin_za_angles[0] != 0 ) 
    throw runtime_error( "The first value in *fos_iyin_za_angles* must be 0." );
  if( last(fos_iyin_za_angles) != 180 ) 
    throw runtime_error( "The last value in *fos_iyin_za_angles* must be 180." );
  if( fos_za_interporder < 1 )
    throw runtime_error( "The argument *fos_za_interporder* must be >= 1." );
  if( fos_iyin_za_angles.nelem() <= fos_za_interporder )
    throw runtime_error( "The length of *fos_iyin_za_angles* must at least "
                         "be *fos_za_interporder* + 1." );
  if( fos_n < 0 )
    throw runtime_error( "The argument *fos_n* must be >= 0." );

  // Switch to order 0 if not primary call
  // (This happens after surface reflection. If fos_n used (and >=1), new
  // surface relections are created ..., and recursion never ends.)
  Index n = fos_n;
  if( !iy_agenda_call1 )
    { n = 0; }

  fos( ws, iy, iy_aux, ppath, diy_dx, stokes_dim, f_grid, atmosphere_dim,
       p_grid, z_field, t_field, vmr_field, abs_species, wind_u_field, 
       wind_v_field, wind_w_field, mag_u_field, mag_v_field, mag_w_field,
       cloudbox_on, cloudbox_limits, pnd_field,
       scat_data, particle_masses, iy_unit, iy_aux_vars, jacobian_do, 
       ppath_agenda, propmat_clearsky_agenda,
       iy_main_agenda, iy_space_agenda, iy_surface_agenda, iy_agenda_call1,
       iy_transmittance, rte_pos, rte_los, rte_pos2, rte_alonglos_v, 
       ppath_lmax, ppath_lraytrace, fos_scatint_angles, fos_iyin_za_angles, 
       fos_za_interporder, n, 0, verbosity );
}
*/

/* Workspace method: Doxygen documentation will be auto-generated */
void iyHybrid(Workspace& ws,
              Matrix& iy,
              ArrayOfMatrix& iy_aux,
              ArrayOfTensor3& diy_dx,
              Vector& ppvar_p,
              Vector& ppvar_t,
              EnergyLevelMap& ppvar_nlte,
              Matrix& ppvar_vmr,
              Matrix& ppvar_wind,
              Matrix& ppvar_mag,
              Matrix& ppvar_pnd,
              Matrix& ppvar_f,
              Tensor3& ppvar_iy,
              Tensor4& ppvar_trans_cumulat,
              Tensor4& ppvar_trans_partial,
              const Index& iy_id,
              const Index& stokes_dim,
              const Vector& f_grid,
              const Index& atmosphere_dim,
              const Vector& p_grid,
              const Tensor3& t_field,
              const EnergyLevelMap& nlte_field,
              const Tensor4& vmr_field,
              const ArrayOfArrayOfSpeciesTag& abs_species,
              const Tensor3& wind_u_field,
              const Tensor3& wind_v_field,
              const Tensor3& wind_w_field,
              const Tensor3& mag_u_field,
              const Tensor3& mag_v_field,
              const Tensor3& mag_w_field,
              const Index& cloudbox_on,
              const ArrayOfIndex& cloudbox_limits,
              const Tensor4& pnd_field,
              const ArrayOfTensor4& dpnd_field_dx,
              const ArrayOfString& scat_species,
              const ArrayOfArrayOfSingleScatteringData& scat_data,
              const String& iy_unit,
              const ArrayOfString& iy_aux_vars,
              const Index& jacobian_do,
              const ArrayOfRetrievalQuantity& jacobian_quantities,
              const Agenda& propmat_clearsky_agenda,
              const Agenda& water_p_eq_agenda,
              const String& rt_integration_option,               
              const Agenda& iy_main_agenda,
              const Agenda& iy_space_agenda,
              const Agenda& iy_surface_agenda,
              const Agenda& iy_cloudbox_agenda,
              const Index& iy_agenda_call1,
              const Tensor3& iy_transmittance,
              const Ppath& ppath,
              const Vector& rte_pos2,
              const Numeric& rte_alonglos_v,
              const Tensor3& surface_props_data,
              const Tensor7& cloudbox_field,
              const Vector& za_grid,
              const Index& Naa,
              const Index& t_interp_order,
              const Verbosity& verbosity) {
  // If cloudbox off, switch to use clearsky method
  if (!cloudbox_on) {
    Tensor4 dummy;
    iyEmissionStandard(ws,
                       iy,
                       iy_aux,
                       diy_dx,
                       ppvar_p,
                       ppvar_t,
                       ppvar_nlte,
                       ppvar_vmr,
                       ppvar_wind,
                       ppvar_mag,
                       ppvar_f,
                       ppvar_iy,
                       ppvar_trans_cumulat,
                       ppvar_trans_partial,
                       iy_id,
                       stokes_dim,
                       f_grid,
                       atmosphere_dim,
                       p_grid,
                       t_field,
                       nlte_field,
                       vmr_field,
                       abs_species,
                       wind_u_field,
                       wind_v_field,
                       wind_w_field,
                       mag_u_field,
                       mag_v_field,
                       mag_w_field,
                       cloudbox_on,
                       iy_unit,
                       iy_aux_vars,
                       jacobian_do,
                       jacobian_quantities,
                       ppath,
                       rte_pos2,
                       propmat_clearsky_agenda,
                       water_p_eq_agenda,
                       rt_integration_option,
                       iy_main_agenda,
                       iy_space_agenda,
                       iy_surface_agenda,
                       iy_cloudbox_agenda,
                       iy_agenda_call1,
                       iy_transmittance,
                       rte_alonglos_v,
                       surface_props_data,
                       verbosity);
    return;
  }
  //  Init Jacobian quantities?
  const Index j_analytical_do = jacobian_do ? do_analytical_jacobian<1>(jacobian_quantities) : 0;

  // Some basic sizes
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index nq = j_analytical_do ? jacobian_quantities.nelem() : 0;

  // Radiative background index
  const Index rbi = ppath_what_background(ppath);

  // Throw error if unsupported features are requested
  if (atmosphere_dim != 1)
    throw runtime_error(
        "With cloudbox on, this method handles only 1D calculations.");
  if (Naa < 3) throw runtime_error("Naa must be > 2.");
  if (jacobian_do)
    if (dpnd_field_dx.nelem() != jacobian_quantities.nelem())
      throw runtime_error(
          "*dpnd_field_dx* not properly initialized:\n"
          "Number of elements in dpnd_field_dx must be equal number of jacobian"
          " quantities.\n(Note: jacobians have to be defined BEFORE *pnd_field*"
          " is calculated/set.");
  if (rbi < 1 || rbi > 9)
    throw runtime_error(
        "ppath.background is invalid. Check your "
        "calculation of *ppath*?");
  if (rbi == 3 || rbi == 4)
    throw runtime_error(
        "The propagation path ends inside or at boundary of "
        "the cloudbox.\nFor this method, *ppath* must be "
        "calculated in this way:\n   ppathCalc( cloudbox_on = 0 ).");
  // iy_aux_vars checked below
  // Checks of i_field
  if (cloudbox_field.ncols() != stokes_dim)
    throw runtime_error(
        "Obtained *cloudbox_field* number of Stokes elements inconsistent with "
        "*stokes_dim*.");
  if (cloudbox_field.nrows() != 1)
    throw runtime_error(
        "Obtained *cloudbox_field* has wrong number of azimuth angles.");
  if (cloudbox_field.npages() != za_grid.nelem())
    throw runtime_error(
        "Obtained *cloudbox_field* number of zenith angles inconsistent with "
        "*za_grid*.");
  if (cloudbox_field.nbooks() != 1)
    throw runtime_error(
        "Obtained *cloudbox_field* has wrong number of longitude points.");
  if (cloudbox_field.nshelves() != 1)
    throw runtime_error(
        "Obtained *cloudbox_field* has wrong number of latitude points.");
  if (cloudbox_field.nvitrines() != cloudbox_limits[1] - cloudbox_limits[0] + 1)
    throw runtime_error(
        "Obtained *cloudbox_field* number of pressure points inconsistent with "
        "*cloudbox_limits*.");
  if (cloudbox_field.nlibraries() != nf)
    throw runtime_error(
        "Obtained *cloudbox_field* number of frequency points inconsistent with "
        "*f_grid*.");
  
  // Set diy_dpath if we are doing are doing jacobian calculations
  ArrayOfTensor3 diy_dpath = j_analytical_do ? get_standard_diy_dpath(jacobian_quantities, np, nf, ns, false) : ArrayOfTensor3(0);
  
  // Set the species pointers if we are doing jacobian
  const ArrayOfIndex jac_species_i = j_analytical_do ? get_pointers_for_analytical_species(jacobian_quantities, abs_species) : ArrayOfIndex(0);
  
  // Start diy_dx out if we are doing the first run and are doing jacobian calculations
  if (j_analytical_do and iy_agenda_call1) diy_dx = get_standard_starting_diy_dx(jacobian_quantities, np, nf, ns, false);
  
  // Checks that the scattering species are treated correctly if their derivatives are needed (we can here discard the Array)
  if (j_analytical_do and iy_agenda_call1) get_pointers_for_scat_species(jacobian_quantities, scat_species, cloudbox_on);

  // Init iy_aux and fill where possible
  const Index naux = iy_aux_vars.nelem();
  iy_aux.resize(naux);
  //
  for (Index i = 0; i < naux; i++) {
    iy_aux[i].resize(nf, ns);

    if (iy_aux_vars[i] == "Optical depth") { /*pass*/
    }                                        // Filled below
    else if (iy_aux_vars[i] == "Radiative background")
      iy_aux[i] = (Numeric)min((Index)2, rbi - 1);
    else {
      ostringstream os;
      os << "The only allowed strings in *iy_aux_vars* are:\n"
         << "  \"Radiative background\"\n"
         << "  \"Optical depth\"\n"
         << "but you have selected: \"" << iy_aux_vars[i] << "\"";
      throw runtime_error(os.str());
    }
  }

  // Get atmospheric and radiative variables along the propagation path
  ppvar_trans_cumulat.resize(np, nf, ns, ns);
  ppvar_trans_partial.resize(np, nf, ns, ns);
  ppvar_iy.resize(nf, ns, np);

  ArrayOfTransmissionMatrix lyr_tra(np, TransmissionMatrix(nf, ns));
  ArrayOfRadiationVector lvl_rad(np, RadiationVector(nf, ns));
  ArrayOfArrayOfRadiationVector dlvl_rad(
      np, ArrayOfRadiationVector(nq, RadiationVector(nf, ns)));
  ArrayOfRadiationVector src_rad(np, RadiationVector(nf, ns));
  ArrayOfArrayOfRadiationVector dsrc_rad(
      np, ArrayOfRadiationVector(nq, RadiationVector(nf, ns)));

  ArrayOfArrayOfTransmissionMatrix dlyr_tra_above(
      np, ArrayOfTransmissionMatrix(nq, TransmissionMatrix(nf, ns)));
  ArrayOfArrayOfTransmissionMatrix dlyr_tra_below(
      np, ArrayOfTransmissionMatrix(nq, TransmissionMatrix(nf, ns)));

  ArrayOfIndex clear2cloudy;
  //
  if (np == 1 && rbi == 1) {  // i.e. ppath is totally outside the atmosphere:
    ppvar_p.resize(0);
    ppvar_t.resize(0);
    ppvar_vmr.resize(0, 0);
    ppvar_wind.resize(0, 0);
    ppvar_mag.resize(0, 0);
    ppvar_f.resize(0, 0);
    ppvar_trans_cumulat = 0;
    ppvar_trans_partial = 0;
    for (Index iv = 0; iv < nf; iv++) {
      for (Index is = 0; is < ns; is++) {
        ppvar_trans_cumulat(0,iv,is,is) = 1;
        ppvar_trans_partial(0,iv,is,is) = 1;
      }
    }
  } else {
    // Basic atmospheric variables
    get_ppath_atmvars(ppvar_p,
                      ppvar_t,
                      ppvar_nlte,
                      ppvar_vmr,
                      ppvar_wind,
                      ppvar_mag,
                      ppath,
                      atmosphere_dim,
                      p_grid,
                      t_field,
                      nlte_field,
                      vmr_field,
                      wind_u_field,
                      wind_v_field,
                      wind_w_field,
                      mag_u_field,
                      mag_v_field,
                      mag_w_field);

    get_ppath_f(
        ppvar_f, ppath, f_grid, atmosphere_dim, rte_alonglos_v, ppvar_wind);

    // here, the cloudbox is on, ie we don't need to check and branch this here
    // anymore.
    ArrayOfMatrix ppvar_dpnd_dx;
    //
    get_ppath_cloudvars(clear2cloudy,
                        ppvar_pnd,
                        ppvar_dpnd_dx,
                        ppath,
                        atmosphere_dim,
                        cloudbox_limits,
                        pnd_field,
                        dpnd_field_dx);

    // Size radiative variables always used
    Vector B(nf);
    PropagationMatrix K_this(nf, ns), K_past(nf, ns), Kp(nf, ns);
    StokesVector a(nf, ns), S(nf, ns), Sp(nf, ns);
    ArrayOfIndex lte(np);

    // Init variables only used if analytical jacobians done
    Vector dB_dT(0);
    ArrayOfPropagationMatrix dK_this_dx(nq), dK_past_dx(nq), dKp_dx(nq);
    ArrayOfStokesVector da_dx(nq), dS_dx(nq), dSp_dx(nq);

    // HSE variables
    Index temperature_derivative_position = -1;
    bool do_hse = false;

    if (j_analytical_do) {
      dB_dT.resize(nf);
      FOR_ANALYTICAL_JACOBIANS_DO(
          dK_this_dx[iq] = PropagationMatrix(nf, ns);
          dK_past_dx[iq] = PropagationMatrix(nf, ns);
          dKp_dx[iq] = PropagationMatrix(nf, ns);
          da_dx[iq] = StokesVector(nf, ns);
          dS_dx[iq] = StokesVector(nf, ns);
          dSp_dx[iq] = StokesVector(nf, ns);
          if (jacobian_quantities[iq] == Jacobian::Atm::Temperature) {
            temperature_derivative_position = iq;
            do_hse = jacobian_quantities[iq].Subtag() == "HSE on";
          })
    }
    const bool temperature_jacobian =
        j_analytical_do and do_temperature_jacobian(jacobian_quantities);

    // Loop ppath points and determine radiative properties
    for (Index ip = 0; ip < np; ip++) {
      get_stepwise_blackbody_radiation(
          B, dB_dT, ppvar_f(joker, ip), ppvar_t[ip], temperature_jacobian);

      get_stepwise_clearsky_propmat(ws,
                                    K_this,
                                    S,
                                    lte[ip],
                                    dK_this_dx,
                                    dS_dx,
                                    propmat_clearsky_agenda,
                                    jacobian_quantities,
                                    ppvar_f(joker, ip),
                                    ppvar_mag(joker, ip),
                                    ppath.los(ip, joker),
                                    ppvar_nlte[ip],
                                    ppvar_vmr(joker, ip),
                                    ppvar_t[ip],
                                    ppvar_p[ip],
                                    j_analytical_do);

      if (j_analytical_do)
        adapt_stepwise_partial_derivatives(dK_this_dx,
                                           dS_dx,
                                           jacobian_quantities,
                                           ppvar_f(joker, ip),
                                           ppath.los(ip, joker),
                                           ppvar_vmr(joker, ip),
                                           ppvar_t[ip],
                                           ppvar_p[ip],
                                           jac_species_i,
                                           lte[ip],
                                           atmosphere_dim,
                                           j_analytical_do);

      if (clear2cloudy[ip] + 1) {
        get_stepwise_scattersky_propmat(a,
                                        Kp,
                                        da_dx,
                                        dKp_dx,
                                        jacobian_quantities,
                                        ppvar_pnd(joker, Range(ip, 1)),
                                        ppvar_dpnd_dx,
                                        ip,
                                        scat_data,
                                        ppath.los(ip, joker),
                                        ppvar_t[Range(ip, 1)],
                                        atmosphere_dim,
                                        jacobian_do);
        a += K_this;
        K_this += Kp;

        if (j_analytical_do)
          FOR_ANALYTICAL_JACOBIANS_DO(da_dx[iq] += dK_this_dx[iq];
                                      dK_this_dx[iq] += dKp_dx[iq];)

        Vector aa_grid;
        nlinspace(aa_grid, 0, 360, Naa);
        //
        get_stepwise_scattersky_source(Sp,
                                       dSp_dx,
                                       jacobian_quantities,
                                       ppvar_pnd(joker, ip),
                                       ppvar_dpnd_dx,
                                       ip,
                                       scat_data,
                                       cloudbox_field,
                                       za_grid,
                                       aa_grid,
                                       ppath.los(Range(ip, 1), joker),
                                       ppath.gp_p[ip],
                                       ppvar_t[Range(ip, 1)],
                                       atmosphere_dim,
                                       jacobian_do,
                                       t_interp_order);
        S += Sp;

        if (j_analytical_do)
          FOR_ANALYTICAL_JACOBIANS_DO(dS_dx[iq] += dSp_dx[iq];)
      } else {  // no particles present at this level
        a = K_this;
        if (j_analytical_do)
          FOR_ANALYTICAL_JACOBIANS_DO(da_dx[iq] = dK_this_dx[iq];)
      }

      if (ip not_eq 0) {
        const Numeric dr_dT_past =
            do_hse ? ppath.lstep[ip - 1] / (2.0 * ppvar_t[ip - 1]) : 0;
        const Numeric dr_dT_this =
            do_hse ? ppath.lstep[ip - 1] / (2.0 * ppvar_t[ip]) : 0;
        stepwise_transmission(lyr_tra[ip],
                              dlyr_tra_above[ip],
                              dlyr_tra_below[ip],
                              K_past,
                              K_this,
                              dK_past_dx,
                              dK_this_dx,
                              ppath.lstep[ip - 1],
                              dr_dT_past,
                              dr_dT_this,
                              temperature_derivative_position);
      }

      stepwise_source(src_rad[ip],
                      dsrc_rad[ip],
                      K_this,
                      a,
                      S,
                      dK_this_dx,
                      da_dx,
                      dS_dx,
                      B,
                      dB_dT,
                      jacobian_quantities,
                      jacobian_do);

      swap(K_past, K_this);
      swap(dK_past_dx, dK_this_dx);
    }
  }

  const ArrayOfTransmissionMatrix tot_tra =
      cumulative_transmission(lyr_tra, CumulativeTransmission::Forward);

  // iy_transmittance
  Tensor3 iy_trans_new;
  if (iy_agenda_call1)
    iy_trans_new = tot_tra[np - 1];
  else
    iy_transmittance_mult(iy_trans_new, iy_transmittance, tot_tra[np - 1]);

  // Copy transmission to iy_aux
  for (Index i = 0; i < naux; i++)
    if (iy_aux_vars[i] == "Optical depth")
      for (Index iv = 0; iv < nf; iv++)
        iy_aux[i](iv, joker) = -log(ppvar_trans_cumulat(np - 1, iv, 0, 0));

  // Radiative background
  get_iy_of_background(ws,
                       iy,
                       diy_dx,
                       iy_trans_new,
                       iy_id,
                       jacobian_do,
                       jacobian_quantities,
                       ppath,
                       rte_pos2,
                       atmosphere_dim,
                       nlte_field,
                       cloudbox_on,
                       stokes_dim,
                       f_grid,
                       iy_unit,
                       surface_props_data,
                       iy_main_agenda,
                       iy_space_agenda,
                       iy_surface_agenda,
                       iy_cloudbox_agenda,
                       iy_agenda_call1,
                       verbosity);

  lvl_rad[np - 1] = iy;

  // Radiative transfer calculations
  for (Index ip = np - 2; ip >= 0; ip--) {
    lvl_rad[ip] = lvl_rad[ip + 1];
    update_radiation_vector(lvl_rad[ip],
                            dlvl_rad[ip],
                            dlvl_rad[ip + 1],
                            src_rad[ip],
                            src_rad[ip + 1],
                            dsrc_rad[ip],
                            dsrc_rad[ip + 1],
                            lyr_tra[ip + 1],
                            tot_tra[ip],
                            dlyr_tra_above[ip + 1],
                            dlyr_tra_below[ip + 1],
                            PropagationMatrix(),
                            PropagationMatrix(),
                            ArrayOfPropagationMatrix(),
                            ArrayOfPropagationMatrix(),
                            Numeric(),
                            Vector(),
                            Vector(),
                            0,
                            0,
                            RadiativeTransferSolver::Emission);
  }

  // Copy back to ARTS external style
  iy = lvl_rad[0];
  for (Index ip = 0; ip < lvl_rad.nelem(); ip++) {
    ppvar_trans_cumulat(ip, joker, joker, joker) = tot_tra[ip];
    ppvar_trans_partial(ip, joker, joker, joker) = lyr_tra[ip];
    ppvar_iy(joker, joker, ip) = lvl_rad[ip];
    if (j_analytical_do)
      FOR_ANALYTICAL_JACOBIANS_DO(diy_dpath[iq](ip, joker, joker) =
                                      dlvl_rad[ip][iq];);
  }

  // Finalize analytical Jacobians
  if (j_analytical_do)
    rtmethods_jacobian_finalisation(ws,
                                    diy_dx,
                                    diy_dpath,
                                    ns,
                                    nf,
                                    np,
                                    atmosphere_dim,
                                    ppath,
                                    ppvar_p,
                                    ppvar_t,
                                    ppvar_vmr,
                                    iy_agenda_call1,
                                    iy_transmittance,
                                    water_p_eq_agenda,
                                    jacobian_quantities,
                                    jac_species_i);

  // Unit conversions
  if (iy_agenda_call1)
    rtmethods_unit_conversion(iy,
                              diy_dx,
                              ppvar_iy,
                              ns,
                              np,
                              f_grid,
                              ppath,
                              jacobian_quantities,
                              j_analytical_do,
                              iy_unit);
}
