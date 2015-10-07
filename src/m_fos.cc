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
#include "math_funcs.h"
#include "montecarlo.h"
#include "rte.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;


// FOS implemented as an internal function, to allow an recursive algorithm
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
   const Index&                         use_mean_scat_data,
   const ArrayOfArrayOfSingleScatteringData&   scat_data,
   const Matrix&                        particle_masses,
   const String&                        iy_unit,
   const ArrayOfString&                 iy_aux_vars,
   const Index&                         jacobian_do,
   const Agenda&                        ppath_agenda,
   const Agenda&                        blackbody_radiation_agenda,
   const Agenda&                        propmat_clearsky_agenda,
   const Agenda&                        iy_main_agenda,
   const Agenda&                        iy_space_agenda,
   const Agenda&                        iy_surface_agenda,
   const Index&                         iy_agenda_call1,
   const Tensor3&                       iy_transmission,
   const Vector&                        rte_pos,      
   const Vector&                        rte_los,      
   const Vector&                        rte_pos2, 
   const Numeric&                       rte_alonglos_v,      
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

  assert( fos_i >= 0  &&  fos_i <= fos_n );


  // Determine propagation path
  //
  ppath_agendaExecute( ws, ppath, ppath_lraytrace, rte_pos, rte_los, rte_pos2, 
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
  Matrix       ppath_vmr, ppath_pnd, ppath_wind, ppath_mag, ppath_f, ppath_t_nlte;
  Matrix       ppath_blackrad;
  Tensor5      abs_per_species, dummy_dppath_ext_dx;
  Tensor4      ppath_ext, trans_partial, trans_cumulat, pnd_ext_mat;
  Tensor3      pnd_abs_vec, ppath_nlte_source;
  Vector       scalar_tau;
  ArrayOfIndex clear2cloudbox, lte;
  const Tensor4 t_nlte_field_dummy;
  //
  Array<ArrayOfArrayOfSingleScatteringData> scat_data_single;
  ArrayOfArrayOfIndex                extmat_case;  
  //
  if( np > 1 )
    {
      get_ppath_atmvars(  ppath_p, ppath_t, ppath_t_nlte, ppath_vmr,
                          ppath_wind, ppath_mag, 
                          ppath, atmosphere_dim, p_grid, t_field,  t_nlte_field_dummy, vmr_field,
                          wind_u_field, wind_v_field, wind_w_field,
                          mag_u_field, mag_v_field, mag_w_field );
      get_ppath_f(        ppath_f, ppath, f_grid,  atmosphere_dim, 
                          rte_alonglos_v, ppath_wind );
      get_ppath_pmat(     ws, ppath_ext, ppath_nlte_source, lte, abs_per_species, dummy_dppath_ext_dx,
                          propmat_clearsky_agenda, ArrayOfRetrievalQuantity(0), ppath, 
                          ppath_p, ppath_t, ppath_t_nlte, ppath_vmr, ppath_f, 
                          ppath_mag, f_grid, stokes_dim, iaps );
      for( Index i=0; i<lte.nelem(); i++ )
        {
          if( lte[i] == 0 )
            throw runtime_error( "FOS can so far only handle LTE conditions." );
        }
      get_ppath_blackrad( ws, ppath_blackrad, blackbody_radiation_agenda, 
                          ppath, ppath_t, ppath_f );
      if( !cloudbox_on )
        { 
          get_ppath_trans( trans_partial, extmat_case, trans_cumulat, 
                           scalar_tau, ppath, ppath_ext, f_grid, stokes_dim );
        }
      else
        {
          get_ppath_ext(    clear2cloudbox, pnd_abs_vec, pnd_ext_mat, scat_data_single,
                            ppath_pnd, ppath, ppath_t, stokes_dim, ppath_f, 
                            atmosphere_dim, cloudbox_limits, pnd_field, 
                            use_mean_scat_data, scat_data, verbosity );
          get_ppath_trans2( trans_partial, extmat_case, trans_cumulat, 
                            scalar_tau, ppath, ppath_ext, f_grid, stokes_dim, 
                            clear2cloudbox, pnd_ext_mat );
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

  // iy_transmission
  //
  Tensor3 iy_trans_new;
  //
  if( iy_agenda_call1 )
    { iy_trans_new = trans_cumulat(joker,joker,joker,np-1); }
  else
    { 
      iy_transmission_mult( iy_trans_new, iy_transmission, 
                            trans_cumulat(joker,joker,joker,np-1) ); }

  // Radiative background
  //
  {
    Agenda iy_cbox_agenda;
    get_iy_of_background( ws, iy, diy_dx, 
                          iy_trans_new, jacobian_do, ppath, rte_pos2, 
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
                                             ppath_ext(iv,is1,is2,np-1); } } } }
      for( Index j=0; j<auxAbsSpecies.nelem(); j++ )
        { for( Index iv=0; iv<nf; iv++ ) {
            for( Index is1=0; is1<stokes_dim; is1++ ){
              for( Index is2=0; is2<stokes_dim; is2++ ){
                iy_aux[auxAbsSpecies[j]](iv,is1,is2,np-1) = 
                          abs_per_species(auxAbsIsp[j],iv,is1,is2,np-1); } } } }
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

      // Help variables for handling of *use_mean_scat_data*
      Index   nfs=nf, ivf=1;
      if( use_mean_scat_data )
        { nfs = 1;  ivf = 0; }

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
          bool any_particles = clear2cloudbox[ip] >= 0 || 
                               clear2cloudbox[ip+1] >= 0;

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
                      if( clear2cloudbox[ip] >= 0 )
                        { ext_matFromabs_vec( pabs_mat, pnd_abs_vec(iv,joker,
                                          clear2cloudbox[ip]), stokes_dim ); } 
                      if( clear2cloudbox[ip+1] >= 0 )
                        { ext_matFromabs_vec( pabs_mat, pnd_abs_vec(iv,joker,
                                        clear2cloudbox[ip+1]), stokes_dim ); } 

                      // Transmission of step
                      Matrix ext_mat(stokes_dim,stokes_dim);  
                      for( Index is1=0; is1<stokes_dim; is1++ ) {
                        for( Index is2=0; is2<stokes_dim; is2++ ) {
                          ext_mat(is1,is2) = 0.5 * ( pabs_mat(is1,is2) +
                                                  ppath_ext(iv,is1,is2,ip) +
                                                  ppath_ext(iv,is1,is2,ip+1) );
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
                  if( clear2cloudbox[ip] < 0 )
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
                                 pnd_field, use_mean_scat_data, scat_data, 
                                 particle_masses, iy_unit, iy_aux_vars, 
                                 jacobian_do, ppath_agenda, 
                                 blackbody_radiation_agenda, 
                                 propmat_clearsky_agenda, iy_main_agenda, 
                                 iy_space_agenda, iy_surface_agenda, 0,
                                 iy_trans_new, pos, los, rte_pos2, 
                                 rte_alonglos_v, ppath_lraytrace, 
                                 fos_scatint_angles, fos_iyin_za_angles, 
                                 fos_za_interporder, fos_n, fos_i+1,
                                 verbosity );
                            
                            Y1(i,joker,joker) = iyl;
                          }

                        // Zenith angle interpolation of Y
                        ArrayOfGridPosPoly gp( nin );
                        gridpos_poly( gp, fos_iyin_za_angles, 
                                      fos_scatint_angles(joker,0),
                                      fos_za_interporder );
                        Matrix itw( nin, fos_za_interporder+1 );
                        interpweights( itw, gp );
                        //
                        for( Index iv=0; iv<nf; iv++ ) 
                          { 
                            for( Index is1=0; is1<stokes_dim; is1++ ) 
                              { 
                                interp( Y(joker,iv,is1), itw, 
                                        Y1(joker,iv,is1), gp );
                              }
                          }
                      }

                      // Direction of outgoing scattered radiation (which is
                      // reversed to LOS). Note that this outlos is only used
                      // for extracting scattering properties.
                      Vector outlos;
                      mirror_los( outlos, ppath.los(ip,joker), atmosphere_dim );

                      // Determine phase matrix 
                      Tensor4  P( nin, nfs, stokes_dim, stokes_dim );
                      Matrix   P1( stokes_dim, stokes_dim );
                      //
                      for( Index ii=0; ii<nin; ii++ )
                        {
                          for( Index iv=0; iv<nfs; iv++ )
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
                              mult( sp, P(ii,iv*ivf,joker,joker), 
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
                                                    ppath_ext(iv,is1,0,ip) +
                                                    ppath_ext(iv,is1,0,ip+1) );
                          for( Index is2=0; is2<stokes_dim; is2++ )
                            {
                              ext_mat(is1,is2) = 0.5 * (
                                                  ppath_ext(iv,is1,is2,ip) +
                                                  ppath_ext(iv,is1,is2,ip+1) );
                            }
                        }
                      // Particle contribution
                      if( clear2cloudbox[ip] >= 0 )
                        {
                          for( Index is1=0; is1<stokes_dim; is1++ )
                            { 
                              sbar[is1]    += 0.5 * ssource0(iv,is1); 
                              abs_vec[is1] += 0.5 * (
                                     pnd_abs_vec(iv,is1,clear2cloudbox[ip]) );
                              for( Index is2=0; is2<stokes_dim; is2++ )
                                {
                                  ext_mat(is1,is2) += 0.5 * (
                                  pnd_ext_mat(iv,is1,is2,clear2cloudbox[ip]) );
                                }
                            }
                        }
                      if( clear2cloudbox[ip+1] >= 0 )
                        {
                          for( Index is1=0; is1<stokes_dim; is1++ )
                            { 
                              sbar[is1]    += 0.5 * ssource1(iv,is1); 
                              abs_vec[is1] += 0.5 * (
                                    pnd_abs_vec(iv,is1,clear2cloudbox[ip+1]) );
                              for( Index is2=0; is2<stokes_dim; is2++ )
                                {
                                  ext_mat(is1,is2) += 0.5 * (
                                  pnd_ext_mat(iv,is1,is2,clear2cloudbox[ip+1]));
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
                                              ppath_ext(iv,is1,is2,ip); } } } }
          for( Index j=0; j<auxAbsSpecies.nelem(); j++ )
            { for( Index iv=0; iv<nf; iv++ ) {
                for( Index is1=0; is1<stokes_dim; is1++ ){
                  for( Index is2=0; is2<stokes_dim; is2++ ){
                    iy_aux[auxAbsSpecies[j]](iv,is1,is2,ip) = 
                           abs_per_species(auxAbsIsp[j],iv,is1,is2,ip); } } } }
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
      // If any conversion, check that standard form of Planck used
      if( !chk_if_std_blackbody_agenda( ws, blackbody_radiation_agenda ) )
        {
          ostringstream os;
          os << "When any unit conversion is applied, "
             << "*blackbody_radiation_agenda\nmust use "
             << "*blackbody_radiationPlanck* or a corresponding WSM.\nA test "
             << "call of the agenda indicates that this is not the case.";
          throw runtime_error( os.str() );
        }
        
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




/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
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
   const Index&                       use_mean_scat_data,
   const ArrayOfArrayOfSingleScatteringData& scat_data,
   const Matrix&                      particle_masses,
   const String&                      iy_unit,
   const ArrayOfString&               iy_aux_vars,
   const Index&                       jacobian_do,
   const Agenda&                      ppath_agenda,
   const Agenda&                      blackbody_radiation_agenda,
   const Agenda&                      propmat_clearsky_agenda,
   const Agenda&                      iy_main_agenda,
   const Agenda&                      iy_space_agenda,
   const Agenda&                      iy_surface_agenda,
   const Index&                       iy_agenda_call1,
   const Tensor3&                     iy_transmission,
   const Vector&                      rte_pos,      
   const Vector&                      rte_los,      
   const Vector&                      rte_pos2, 
   const Numeric&                     rte_alonglos_v,      
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
       cloudbox_on, cloudbox_limits, pnd_field, use_mean_scat_data,
       scat_data, particle_masses, iy_unit, iy_aux_vars, jacobian_do, 
       ppath_agenda, blackbody_radiation_agenda, propmat_clearsky_agenda,
       iy_main_agenda, iy_space_agenda, iy_surface_agenda, iy_agenda_call1,
       iy_transmission, rte_pos, rte_los, rte_pos2, rte_alonglos_v, 
       ppath_lraytrace, fos_scatint_angles, fos_iyin_za_angles, 
       fos_za_interporder, n, 0, verbosity );
}
