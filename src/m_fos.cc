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
#include "m_xml.h"

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
   const Index&                         scat_data_checked,
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
   const Tensor3&                       iy_transmission,
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

  assert( fos_i >= 0  &&  fos_i <= fos_n );


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
  Matrix       ppath_vmr, ppath_pnd, ppath_wind, ppath_mag, ppath_f, ppath_t_nlte;
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
  const Tensor4   t_nlte_field_empty(0,0,0,0);
  //
  Array<ArrayOfArrayOfSingleScatteringData> scat_data_single;
  ArrayOfArrayOfIndex                extmat_case;  
  //
  if( np > 1 )
    {
      get_ppath_atmvars(  ppath_p, ppath_t, ppath_t_nlte, ppath_vmr,
                          ppath_wind, ppath_mag, 
                          ppath, atmosphere_dim, p_grid, t_field,
                          t_nlte_field_empty, vmr_field,
                          wind_u_field, wind_v_field, wind_w_field,
                          mag_u_field, mag_v_field, mag_w_field );
      
      get_ppath_f( ppath_f, ppath, f_grid,  atmosphere_dim, 
                   rte_alonglos_v, ppath_wind );
      
      get_ppath_pmat( ws, ppath_ext, ppath_nlte_source, lte, abs_per_species, 
                      dummy_dppath_ext_dx, dummy_dppath_nlte_dx,
                      propmat_clearsky_agenda, ArrayOfRetrievalQuantity(0), ppath, 
                      ppath_p, ppath_t, ppath_t_nlte, ppath_vmr, ppath_f, 
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
                             clear2cloudy, ppath_pnd,
                             ppath, ppath_t, stokes_dim, ppath_f, atmosphere_dim,
                             use_mean_scat_data, scat_data, scat_data_checked,
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
                                 pnd_field, use_mean_scat_data,
                                 scat_data, scat_data_checked, 
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
   const Index&                       scat_data_checked,
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
   const Tensor3&                     iy_transmission,
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
       cloudbox_on, cloudbox_limits, pnd_field, use_mean_scat_data,
       scat_data, scat_data_checked, particle_masses, iy_unit, iy_aux_vars, jacobian_do, 
       ppath_agenda, propmat_clearsky_agenda,
       iy_main_agenda, iy_space_agenda, iy_surface_agenda, iy_agenda_call1,
       iy_transmission, rte_pos, rte_los, rte_pos2, rte_alonglos_v, 
       ppath_lmax, ppath_lraytrace, fos_scatint_angles, fos_iyin_za_angles, 
       fos_za_interporder, n, 0, verbosity );
}







/* Workspace method: Doxygen documentation will be auto-generated */
void iyHybrid(
         Workspace&                   ws,
         Matrix&                      iy,
         ArrayOfTensor4&              iy_aux,
         Ppath&                       ppath,
         ArrayOfTensor3&              diy_dx,
   const Index&                       iy_id,
   const Index&                       stokes_dim,
   const Vector&                      f_grid,
   const Index&                       atmosphere_dim,
   const Vector&                      p_grid,
   const Tensor3&                     z_field,
   const Tensor3&                     t_field,
   const Tensor4&                     t_nlte_field,
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
   const ArrayOfTensor4&              dpnd_field_dx,
   const ArrayOfArrayOfSingleScatteringData& scat_data,
   const Index&                       scat_data_checked,
   const Matrix&                      particle_masses,
   const String&                      iy_unit,
   const ArrayOfString&               iy_aux_vars,
   const Index&                       jacobian_do,
   const ArrayOfRetrievalQuantity&    jacobian_quantities,
   const ArrayOfArrayOfIndex&         jacobian_indices,
   const Agenda&                      ppath_agenda,
   const Agenda&                      propmat_clearsky_agenda,
   const Agenda&                      iy_main_agenda,
   const Agenda&                      iy_space_agenda,
   const Agenda&                      iy_surface_agenda,
   const Agenda&                      iy_cloudbox_agenda,
   const Agenda&                      doit_i_field_agenda,
   const Index&                       iy_agenda_call1,
   const Tensor3&                     iy_transmission,
   const Vector&                      rte_pos,      
   const Vector&                      rte_los,      
   const Vector&                      rte_pos2,
   const Numeric&                     rte_alonglos_v,      
   const Numeric&                     ppath_lmax,      
   const Numeric&                     ppath_lraytrace,
   const Index&                       Naa,   
   const String&                      pfct_method,
   const Verbosity&                   verbosity )
{
  // If cloudbox off, switch to use clearsky method
  if( !cloudbox_on )
    {
      iyEmissionStandard( ws, iy, iy_aux, ppath, diy_dx, iy_id, stokes_dim,
                          f_grid, atmosphere_dim, p_grid, z_field, t_field,
                          t_nlte_field, vmr_field, abs_species,
                          wind_u_field, wind_v_field, wind_w_field,
                          mag_u_field, mag_v_field, mag_w_field,
                          cloudbox_on, iy_unit, iy_aux_vars,
                          jacobian_do, jacobian_quantities, jacobian_indices,
                          ppath_agenda, propmat_clearsky_agenda, iy_main_agenda,
                          iy_space_agenda, iy_surface_agenda,
                          iy_cloudbox_agenda, iy_agenda_call1, iy_transmission,
                          rte_pos, rte_los, rte_pos2, rte_alonglos_v,
                          ppath_lmax, ppath_lraytrace, verbosity );
      return;
    }
  

  // Throw error if unsupported features are requested
  if( atmosphere_dim != 1 )
    throw runtime_error( "With cloudbox on, this method handles only "
                         "1D calculations." );
  if( !iy_agenda_call1 )
    throw runtime_error( "With cloudbox on, recursive usage not possible "
                         "(iy_agenda_call1 must be 1)." );
  if( iy_transmission.ncols() )
    throw runtime_error( "*iy_transmission* must be empty." );
  if( cloudbox_limits[0] != 0  ||  cloudbox_limits[1] != p_grid.nelem()-1 )
        throw runtime_error( "The cloudbox must be set to cover the complete "
                             "atmosphere." );
  // for now have that here. when all iy* WSM using scat_data are fixed to new
  // type scat_data, then put check inot (i)yCalc and remove here.
  if( scat_data_checked != 1 )
    throw runtime_error( "The scat_data must be flagged to have "
                         "passed a consistency check (scat_data_checked=1)." );
  
  
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
  const Index ne = pnd_field.nbooks();
  const Index nq = jacobian_quantities.nelem();

  // Obtain i_field
  //
  Tensor7 doit_i_field;
  Vector  scat_za_grid;
  //
  {
    Vector scat_aa_grid;
    //
    doit_i_field_agendaExecute( ws, doit_i_field, scat_za_grid, scat_aa_grid,
                                doit_i_field_agenda );
    if( doit_i_field.ncols() != stokes_dim  )
      throw runtime_error(
        "Obtained *doit_i_field* has wrong number of Stokes elements."  );
    if( doit_i_field.nrows() != 1  )
      throw runtime_error(
        "Obtained *doit_i_field* has wrong number of azimuth angles."  );
    if( doit_i_field.nbooks() != 1  )
      throw runtime_error(
        "Obtained *doit_i_field* has wrong number of longitude points."  );
    if( doit_i_field.nshelves() != 1  )
      throw runtime_error(
        "Obtained *doit_i_field* has wrong number of latitude points."  );
    if( doit_i_field.nvitrines() != cloudbox_limits[1]-cloudbox_limits[0]+1  )
      throw runtime_error(
        "Obtained *doit_i_field* has wrong number of pressure points."  );
    if( doit_i_field.nlibraries() != nf  )
      throw runtime_error(
        "Obtained *doit_i_field* has wrong number of frequency points."  );
  }

  // For a brief description of internal variables used, see
  // iyEmissionStandard. 


  //### jacobian part #########################################################
  // Initialise analytical jacobians (diy_dx and help variables)
  //
  Index           j_analytical_do = 0;
  ArrayOfTensor3  diy_dpath; 
  ArrayOfIndex    jac_species_i(0), jac_scat_i(0), jac_is_t(0), jac_wind_i(0);
  ArrayOfIndex    jac_mag_i(0), jac_other(0), jac_to_integrate(0); 
  // Flags for partial derivatives of propmat
  const PropmatPartialsData ppd(jacobian_quantities);
  //
  if( jacobian_do ) 
    { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }
  //
  if( !j_analytical_do )
    { diy_dx.resize( 0 ); }
  else 
    {
      diy_dpath.resize( nq ); 
      jac_species_i.resize( nq ); 
      jac_scat_i.resize( nq ); 
      jac_is_t.resize( nq ); 
      jac_wind_i.resize( nq ); 
      jac_mag_i.resize( nq ); 
      jac_other.resize(nq);
      jac_to_integrate.resize(nq);
      //
      FOR_ANALYTICAL_JACOBIANS_DO( 
        if( jacobian_quantities[iq].Integration() )
        {
            diy_dpath[iq].resize( 1, nf, ns ); 
            diy_dpath[iq] = 0.0;
        }
        else
        {
            diy_dpath[iq].resize( np, nf, ns ); 
            diy_dpath[iq] = 0.0;
        }
      )

      const ArrayOfString scat_species(0);
      
      get_pointers_for_analytical_jacobians( jac_species_i, jac_scat_i, jac_is_t, 
                                             jac_wind_i, jac_mag_i, jac_to_integrate, 
                                             jacobian_quantities,
                                             abs_species, scat_species );
      FOR_ANALYTICAL_JACOBIANS_DO( 
        jac_other[iq] = ppd.is_this_propmattype(iq)?JAC_IS_OTHER:JAC_IS_NONE; 
        if( jac_to_integrate[iq] == JAC_IS_FLUX )
          throw std::runtime_error("This method can not perform flux calculations.\n");
      )

      if( iy_agenda_call1 )
        {
          diy_dx.resize( nq ); 
          //
          FOR_ANALYTICAL_JACOBIANS_DO( diy_dx[iq].resize( 
            jacobian_indices[iq][1]-jacobian_indices[iq][0]+1, nf, ns ); 
            diy_dx[iq] = 0.0;
          )
        }
    } 
  //###########################################################################


  //=== iy_aux part ===========================================================
  Index auxPressure    = -1,
        auxTemperature = -1,
        auxAbsSum      = -1,
        auxPartExt     = -1,
        auxIy          = -1,
        auxTrans       = -1,
        auxOptDepth    = -1,
        auxFarRotTotal = -1,
        auxFarRotSpeed = -1;
  Index ife = -1;     // Index in abs_per:species matching free electons
  ArrayOfIndex iaps(0);
  ArrayOfIndex auxAbsSpecies(0), auxAbsIsp(0);
  ArrayOfIndex auxVmrSpecies(0), auxVmrIsp(0);
  ArrayOfIndex auxPartCont(0),   auxPartContI(0);
  ArrayOfIndex auxPartField(0),  auxPartFieldI(0);
  //
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
        else if( iy_aux_vars[i] == "Particle extinction, summed" )
          { 
            auxPartExt = i;   
            iy_aux[i].resize( nf, ns, ns, np ); 
            iy_aux[i] = 0;
          }
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
        else if( iy_aux_vars[i] == "iy"   &&  auxIy < 0 )
          { auxIy = i;           iy_aux[i].resize( nf, ns, 1, np ); }
        else if( iy_aux_vars[i] == "Transmission"   &&  auxTrans < 0 )
          { auxTrans = i;        iy_aux[i].resize( nf, ns, ns, np ); }
        else if( iy_aux_vars[i] == "Optical depth" )
          { auxOptDepth = i;     iy_aux[i].resize( nf, 1, 1, 1 ); }
        else if( iy_aux_vars[i] == "Faraday rotation" )
          { auxFarRotTotal = i; iy_aux[i].resize( nf, 1, 1, 1 ); iy_aux[i] = 0; }
        else if( iy_aux_vars[i] == "Faraday speed" )
          { auxFarRotSpeed = i; iy_aux[i].resize( nf, 1, 1, np ); iy_aux[i] = 0; }
        else
          {
            ostringstream os;
            os << "In *iy_aux_vars* you have included: \"" << iy_aux_vars[i]
               << "\"\nThis choice is not recognised.";
            throw runtime_error( os.str() );
          }
      }
    
    // Special stuff to handle Faraday rotation
    if( auxFarRotTotal>=0  ||  auxFarRotSpeed>=0 )  
      {
        if( stokes_dim < 3 )
          throw runtime_error( 
                   "To include Faraday rotation, stokes_dim >= 3 is required." );

        // Determine species index of free electrons
        for( Index sp = 0; sp < abs_species.nelem() && ife < 0; sp++ )
          {
            if (abs_species[sp][0].Type() == SpeciesTag::TYPE_FREE_ELECTRONS)
              { ife = sp; }
          }
        // If not found, then aux values already set to zero
        if( ife < 0 )
          {
            auxFarRotTotal = -1;
            auxFarRotSpeed = -1;
          }
        else
          {
            const Index ihit = find_first( iaps, ife );
            if( ihit >= 0 )
              { ife = ihit; }
            else
              { 
                iaps.push_back(ife); 
                ife = iaps.nelem() - 1; 
              }
          }
      }
  }
  //===========================================================================


  // Get atmospheric and RT quantities for each ppath point/step
  //
  Vector              ppath_p, ppath_t;
  Matrix              ppath_vmr, ppath_pnd, ppath_wind, ppath_mag;
  Matrix              ppath_f, ppath_t_nlte;
  Matrix              ppath_blackrad, dppath_blackrad_dt;
  ArrayOfMatrix       ppath_dpnd_dx;
  Tensor5             dtrans_partial_dx_above, dtrans_partial_dx_below;
  Vector              scalar_tau;
  Tensor3             pnd_abs_vec;
  ArrayOfPropagationMatrix pnd_ext_mat, ppath_ext;
  ArrayOfArrayOfPropagationMatrix dppath_ext_dx;
  Tensor4             ppath_scat_source;
  ArrayOfStokesVector ppath_nlte_source;
  ArrayOfArrayOfStokesVector  dppath_nlte_dx, dppath_nlte_source_dx;
  Tensor4             trans_partial, trans_cumulat;
  ArrayOfArrayOfPropagationMatrix abs_per_species;
  ArrayOfIndex        clear2cloudy;
  ArrayOfIndex        lte;
  ArrayOfArrayOfIndex extmat_case;   
  Array<ArrayOfArrayOfSingleScatteringData> scat_data_single;
  
  if( np > 1 )
    {
      get_ppath_atmvars( ppath_p, ppath_t, ppath_t_nlte, ppath_vmr,
                         ppath_wind, ppath_mag, 
                         ppath, atmosphere_dim, p_grid, t_field, t_nlte_field, 
                         vmr_field, wind_u_field, wind_v_field, wind_w_field,
                         mag_u_field, mag_v_field, mag_w_field );      

      get_ppath_f( ppath_f, ppath, f_grid,  atmosphere_dim, 
                   rte_alonglos_v, ppath_wind );

      get_ppath_pmat_and_tmat( ws, ppath_ext, ppath_nlte_source, lte, abs_per_species,
                               dppath_ext_dx, dppath_nlte_source_dx,
                               trans_partial, dtrans_partial_dx_above,
                               dtrans_partial_dx_below, extmat_case, clear2cloudy,
                               trans_cumulat, scalar_tau, pnd_ext_mat, pnd_abs_vec,
                               ppath_pnd, ppath_dpnd_dx, scat_data_single,
                               propmat_clearsky_agenda, jacobian_quantities,
                               ppd, ppath, ppath_p, ppath_t, ppath_t_nlte,
                               ppath_vmr, ppath_mag, ppath_f, f_grid, 
                               jac_species_i, jac_is_t, jac_wind_i, jac_mag_i,
                               jac_to_integrate, jac_other, iaps,
                               scat_data, scat_data_checked,
                               pnd_field, dpnd_field_dx,
                               cloudbox_limits, 0,
                               atmosphere_dim, stokes_dim,
                               jacobian_do, cloudbox_on, verbosity );
      
      get_ppath_blackrad( ppath_blackrad, ppath, ppath_t, ppath_f );

      get_dppath_blackrad_dt( dppath_blackrad_dt, ppath_t, ppath_f, jac_is_t, 
                              j_analytical_do );

      // for debugging
      //WriteXML( "ascii", doit_i_field, "doit_i_field_iyHybrid.xml", 0,
      //          "doit_i_field", "", "", verbosity );
      if( pfct_method=="interpolate" )
      {
        //cout << "T-interpolating pha_mat\n";
        get_ppath_scat_source( ppath_scat_source,
                               scat_data, scat_data_checked, doit_i_field, scat_za_grid,
                               f_grid, stokes_dim, ppath, ppath_t, ppath_pnd, 
                               j_analytical_do, Naa, verbosity );
      }
      else
      {
        Numeric rtp_temp_id;
        if( pfct_method=="low" )
          rtp_temp_id = -9.;
        else if( pfct_method=="high" )
          rtp_temp_id = -19.;
        else //if( pfct_method=="median" )
          rtp_temp_id = -99.;

        get_ppath_scat_source_fixT( ppath_scat_source,
                               scat_data, scat_data_checked, doit_i_field, scat_za_grid,
                               f_grid, stokes_dim, ppath, ppath_pnd,
                               j_analytical_do, Naa, rtp_temp_id, verbosity );
      }
      // for debugging
      //WriteXML( "ascii", ppath_scat_source, "ppath_scat_source.xml", 0,
      //          "ppath_scat_source", "", "", verbosity );

    }
  else
    {  
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
    { iy_transmission_mult( iy_trans_new, iy_transmission, 
                            trans_cumulat(joker,joker,joker,np-1) ); }


  // Get *iy* at end ppath by interpoling doit_i_field
  {
    Tensor4 i_field = doit_i_field(joker,joker,0,0,joker,0,joker);

    ArrayOfGridPos gp_p(1), gp_za(1);
    gridpos_copy( gp_p[0], ppath.gp_p[np-1] );
    gridpos( gp_za, scat_za_grid, Vector(1,ppath.los(np-1,0)) );

    Tensor3 itw( 1, 1, 4 );
    interpweights( itw, gp_p, gp_za );

    iy.resize( nf, stokes_dim );

    Matrix m(1,1);
    for( Index f=0; f<nf; f++ )
      {
        for( Index s=0; s<stokes_dim; s++ )
          {
            interp( m, itw, doit_i_field(f,joker,0,0,joker,0,s), gp_p, gp_za );
            iy(f,s) = m(0,0);
          }
      }
  }


  //=== iy_aux part ===========================================================
  // Fill parts of iy_aux that are defined even for np=1.
  // Radiance 
  if( auxIy >= 0 ) 
    { iy_aux[auxIy](joker,joker,0,np-1) = iy; }
  if( auxOptDepth >= 0 ) 
    {
      if( np == 1 )
        { iy_aux[auxOptDepth] = 0; }
      else
        { iy_aux[auxOptDepth](joker,0,0,0) = scalar_tau; }
    } 
  if( auxTrans >= 0 ) // Complete tensor filled!
    { 
      if( np == 1 )
        { for( Index iv=0; iv<nf; iv++ ) {
            id_mat( iy_aux[auxTrans](iv,joker,joker,0) ); } }
      else
        { iy_aux[auxTrans] = trans_cumulat; }
    }
  // Faraday rotation, total
  if( auxFarRotTotal >= 0 )
    { for( Index iv=0; iv<nf; iv++ ) {
        iy_aux[auxFarRotTotal](iv,0,0,0) = 0; } }
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
          // Extinction
          if( auxPartExt >= 0  && clear2cloudy[np-1] >= 0 ) 
            { 
              const Index ic = clear2cloudy[np-1];
              for( Index iv=0; iv<nf; iv++ ) {
                for( Index is1=0; is1<ns; is1++ ){
                  for( Index is2=0; is2<ns; is2++ ){
                    iy_aux[auxPartExt](iv,is1,is2,np-1) = 
                                              pnd_ext_mat[ic](iv,is1,is2); } } } 
            } 
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
      // Faraday speed
      if( auxFarRotSpeed >= 0 )
        { for( Index iv=0; iv<nf; iv++ ) {
            iy_aux[auxFarRotSpeed](iv,0,0,np-1) = 0.5 *
                                          abs_per_species[np-1][ife](iv,1,2); } }
      //=======================================================================

      
      //=======================================================================
      // Loop ppath steps
/*
      Tensor3 layer_bulk_scatsource( nf, stokes_dim, np-1, 0. ); // just for debugging
      Tensor3 layer_bulk_emissource( nf, stokes_dim, np-1, 0. ); // just for debugging
      Tensor3 layer_gas_emissource( nf, stokes_dim, np-1, 0. ); // just for debugging
      Tensor3 layer_bulk_ext( nf, stokes_dim, np-1, 0. ); // just for debugging
      Tensor3 layer_gas_ext( nf, stokes_dim, np-1, 0. ); // just for debugging
*/
      Vector bulk_scat_source(2), pabs(2);
      Matrix pext(2,stokes_dim);
      Matrix  sourcebar(nf,stokes_dim);
      Tensor3 extbar(nf,stokes_dim, stokes_dim);
      for( Index ip=np-2; ip>=0; ip-- )
        {

          //cout << "atmlev #" << ip << "\n";

          // Path step average of K: extbar
          // Path step average of source terms: sourcebar

          for( Index iv=0; iv<nf; iv++ )  
          { 
            //cout << "freq #" << iv << "\n";

            for( Index is1=0; is1<stokes_dim; is1++ )  
            {
              // Scattering source term
              //
              bulk_scat_source = 0.;
              for( Index ise=0; ise<ne; ise++ )
              {
                bulk_scat_source[0] += ppath_pnd(ise,ip) *
                                       ppath_scat_source(iv,is1,ise,ip);
                bulk_scat_source[1] += ppath_pnd(ise,ip+1) *
                                       ppath_scat_source(iv,is1,ise,ip+1);
              }

              // Emission and total source term
              //
              pabs = 0.;
              if( clear2cloudy[ip] > -1 )
                pabs[0] = pnd_abs_vec(iv,is1,clear2cloudy[ip]);
              if( clear2cloudy[ip+1] > -1 )
                pabs[1] = pnd_abs_vec(iv,is1,clear2cloudy[ip+1]);

              sourcebar(iv,is1) = 0.5 * (
                                  bulk_scat_source[0] + bulk_scat_source[1] //particle scattering contrib
                                + ( ppath_ext[ip](iv,is1)+pabs[0] ) *
                                  ppath_blackrad(iv,ip)
                                + ( ppath_ext[ip+1](iv,is1)+pabs[1] ) *
                                  ppath_blackrad(iv,ip+1) );

              // Extinction
              pext = 0.;
              if( clear2cloudy[ip] > -1 )
                for( Index is2=0; is2<stokes_dim; is2++ )
                  pext(0,is2) = pnd_ext_mat[clear2cloudy[ip]](iv,is1,is2);
              if( clear2cloudy[ip+1] > -1 )
                for( Index is2=0; is2<stokes_dim; is2++ )
                  pext(1,is2) = pnd_ext_mat[clear2cloudy[ip+1]](iv,is1,is2);
              for( Index is2=0; is2<stokes_dim; is2++ )  
              {
                extbar(iv,is1,is2) = 0.5 * ( 
                                     ppath_ext[ip](iv,is1,is2)
                                   + ppath_ext[ip+1](iv,is1,is2)
                                   + pext(0,is2) + pext(1,is2) );
              }

/*
              // just for debugging
              layer_bulk_scatsource(iv, is1, ip) =
                0.5 * ( bulk_scat_source[0] + bulk_scat_source[1] );
              layer_bulk_emissource(iv, is1, ip) =
                0.5 * ( pabs[0]*ppath_blackrad(iv,ip) + 
                        pabs[1]*ppath_blackrad(iv,ip+1) );
              layer_gas_emissource(iv, is1, ip) =
                0.5 * ( ppath_ext[ip](iv,is1)*ppath_blackrad(iv,ip) + 
                        ppath_ext[ip+1](iv,is1)*ppath_blackrad(iv,ip+1) );
              layer_bulk_ext(iv, is1, ip) =
                0.5 * ( pext(0,0) + pext(1,0) );
              layer_gas_ext(iv, is1, ip) =
                0.5 * ( ppath_ext[ip](iv,is1,0) +
                        ppath_ext[ip+1](iv,is1,0) );
              cout << "part scat source = " << layer_bulk_scatsource(iv, is1, ip) << "\n";
              cout << "part emis source = " << layer_bulk_emissource(iv, is1, ip) << "\n";
              cout << "gas emis source  = " << layer_gas_emissource(iv, is1, ip) << "\n";
              cout << "total source     = " << sourcebar(iv,is1) << "\n";
              cout << "part extinction  = " << layer_bulk_ext(iv, is1, ip) << "\n";
              cout << "gas extinction   = " << layer_gas_ext(iv, is1, ip) << "\n";
              cout << "total extinction = " << extbar(iv,is1,0) << "\n";
              cout << "\n";
*/
            }
          }

          // Jacobian code temporarily removed

          for( Index iv=0; iv<nf; iv++ )
            {
              //cout << "freq #" << iv << "\n";

              // Transmitted part
              Vector part1(stokes_dim);
              mult( part1, trans_partial(iv,joker,joker,ip), iy(iv,joker) );

              // Inverse of extbar
              Matrix extbarinv(stokes_dim,stokes_dim);
              inv( extbarinv, extbar(iv,joker,joker) );

              // Emission and scattering
              Vector part2(stokes_dim);
              Vector tmp_vector(stokes_dim);
              mult( tmp_vector, extbarinv, sourcebar(iv,joker) );
              Matrix tmp_matrix(stokes_dim,stokes_dim);
              id_mat(tmp_matrix); 
              tmp_matrix -= trans_partial(iv,joker,joker,ip);
              mult( part2, tmp_matrix, tmp_vector); 

/*
              cout << "iy0                  = " << iy(iv,0) << "\n";
              cout << "T                    = " << trans_partial(iv,0,0,ip) << "\n";
              cout << "part1 = iy0*T        = " << part1[0] << "\n";
              cout << "Kinv                 = " << extbarinv(0,0) << "\n";
              cout << "S                    = " << sourcebar(iv,0) << "\n";
              cout << "(1-T)                = " << tmp_matrix(0,0) << "\n";
              cout << "part2 = (1-T)*Kinv*S = " << part2[0] << "\n";
*/

              // Sum up
              for( Index i=0; i<stokes_dim; i++ )
                { iy(iv,i) = part1[i] + part2[i]; }

//              cout << "iy1 = part1+part2    = " << iy(iv,0) << "\n";
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
            { iy_aux[auxVmrSpecies[j]](0,0,0,ip) =  ppath_vmr(auxVmrIsp[j],ip);}
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
              // Extinction
              if( auxPartExt >= 0  &&  clear2cloudy[ip] >= 0 ) 
                { 
                  const Index ic = clear2cloudy[ip];
                  for( Index iv=0; iv<nf; iv++ ) {
                    for( Index is1=0; is1<ns; is1++ ){
                      for( Index is2=0; is2<ns; is2++ ){
                        iy_aux[auxPartExt](iv,is1,is2,ip) = 
                                              pnd_ext_mat[ic](iv,is1,is2); } } }
                }
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
          // Faraday rotation, total
          if( auxFarRotTotal >= 0 )
            { for( Index iv=0; iv<nf; iv++ ) {
                iy_aux[auxFarRotTotal](iv,0,0,0) += RAD2DEG * ppath.lstep[ip] *
                                0.25 * ( abs_per_species[ip][ife](iv,1,2) +
                                         abs_per_species[ip+1][ife](iv,1,2)); } }
          // Faraday speed
          if( auxFarRotSpeed >= 0 )
            { for( Index iv=0; iv<nf; iv++ ) {  
                iy_aux[auxFarRotSpeed](iv,0,0,ip) = 0.5 *
                                            abs_per_species[ip][ife](iv,1,2); } }
        } // path point loop
      //=======================================================================

      // for debugging
      //WriteXML( "ascii", layer_bulk_scatsource, "layer_bulk_scatsource.xml", 0,
      //          "layer_bulk_scatsource", "", "", verbosity );


      //### jacobian part #####################################################
      // Map jacobians from ppath to retrieval grids
      // (this operation corresponds to the term Dx_i/Dx)
      if( j_analytical_do )
        { 
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_from_path_to_rgrids( diy_dx[iq], jacobian_quantities[iq], 
                               diy_dpath[iq], atmosphere_dim, ppath, ppath_p );
          )
        }
      //#######################################################################
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
      ArrayOfIndex i_pol(ns);
      for( Index is=0; is<ns; is++ )
        { i_pol[is] = is + 1; }

      // Jacobian part (must be converted to Tb before iy for PlanckBT)
      // 
      if( j_analytical_do )
        {
          FOR_ANALYTICAL_JACOBIANS_DO( apply_iy_unit2( diy_dx[iq], iy, iy_unit,
                                                       f_grid, n, i_pol ); )
        } 

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


/* Workspace method: Doxygen documentation will be auto-generated */
void iyHybrid2(
  Workspace&                                ws,
  Matrix&                                   iy,
  ArrayOfTensor4&                           iy_aux,
  Ppath&                                    ppath,
  ArrayOfTensor3&                           diy_dx,
  const Index&                              iy_id,
  const Index&                              stokes_dim,
  const Vector&                             f_grid,
  const Index&                              atmosphere_dim,
  const Vector&                             p_grid,
  const Tensor3&                            z_field,
  const Tensor3&                            t_field,
  const Tensor4&                            t_nlte_field,
  const Tensor4&                            vmr_field,
  const ArrayOfArrayOfSpeciesTag&           abs_species,
  const Tensor3&                            wind_u_field,
  const Tensor3&                            wind_v_field,
  const Tensor3&                            wind_w_field,
  const Tensor3&                            mag_u_field,
  const Tensor3&                            mag_v_field,
  const Tensor3&                            mag_w_field,
  const Index&                              cloudbox_on,
  const ArrayOfIndex&                       cloudbox_limits,
  const Tensor4&                            pnd_field,
  const ArrayOfTensor4&                     dpnd_field_dx,
  const ArrayOfArrayOfSingleScatteringData& scat_data,
  const Index&                              scat_data_checked,
  const Matrix&                             particle_masses _U_,
  const String&                             iy_unit,
  const ArrayOfString&                      iy_aux_vars,
  const Index&                              jacobian_do,
  const ArrayOfRetrievalQuantity&           jacobian_quantities,
  const ArrayOfArrayOfIndex&                jacobian_indices,
  const Agenda&                             ppath_agenda,
  const Agenda&                             propmat_clearsky_agenda,
  const Agenda&                             iy_main_agenda,
  const Agenda&                             iy_space_agenda,
  const Agenda&                             iy_surface_agenda,
  const Agenda&                             iy_cloudbox_agenda,
  const Agenda&                             doit_i_field_agenda,
  const Index&                              iy_agenda_call1,
  const Tensor3&                            iy_transmission,
  const Vector&                             rte_pos,      
  const Vector&                             rte_los,      
  const Vector&                             rte_pos2,
  const Numeric&                            rte_alonglos_v,      
  const Numeric&                            ppath_lmax,      
  const Numeric&                            ppath_lraytrace,
  const Index&                              Naa,   
  const String&                             pfct_method _U_,
  const Verbosity&                          verbosity)
{
  // If cloudbox off, switch to use clearsky method
  if( !cloudbox_on )
  {
    iyEmissionStandard( ws, iy, iy_aux, ppath, diy_dx, iy_id, stokes_dim,
                        f_grid, atmosphere_dim, p_grid, z_field, t_field,
                        t_nlte_field, vmr_field, abs_species,
                        wind_u_field, wind_v_field, wind_w_field,
                        mag_u_field, mag_v_field, mag_w_field,
                        cloudbox_on, iy_unit, iy_aux_vars,
                        jacobian_do, jacobian_quantities, jacobian_indices,
                        ppath_agenda, propmat_clearsky_agenda, iy_main_agenda,
                        iy_space_agenda, iy_surface_agenda,
                        iy_cloudbox_agenda, iy_agenda_call1, iy_transmission,
                        rte_pos, rte_los, rte_pos2, rte_alonglos_v,
                        ppath_lmax, ppath_lraytrace, verbosity );
    return;
  }
  
  
  // Throw error if unsupported features are requested
  if( atmosphere_dim != 1 )
    throw runtime_error( "With cloudbox on, this method handles only "
    "1D calculations." );
  if( !iy_agenda_call1 )
    throw runtime_error( "With cloudbox on,  recursive usage not possible "
    "(iy_agenda_call1 must be 1)." );
  if( iy_transmission.ncols() )
    throw runtime_error( "*iy_transmission* must be empty." );
  if( cloudbox_limits[0] != 0  ||  cloudbox_limits[1] != p_grid.nelem()-1 )
    throw runtime_error( "The cloudbox must be set to cover the complete "
    "atmosphere." );
  if( Naa < 3 )
    throw runtime_error( "Naa must be > 2." );
  // for now have that here. when all iy* WSM using scat_data are fixed to new
  // type scat_data, then put check inot (i)yCalc and remove here.
  if( scat_data_checked != 1 )
    throw runtime_error( "The scat_data must be flagged to have "
                         "passed a consistency check (scat_data_checked=1)." );
  
  
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
  const Index nq = jacobian_quantities.nelem();
  
  // Obtain i_field
  //
  Tensor7 doit_i_field;
  Vector  scat_za_grid;
  Vector scat_aa_grid;
  //
  {
    //
    doit_i_field_agendaExecute( ws, doit_i_field, scat_za_grid, scat_aa_grid,
                                doit_i_field_agenda );
    if( doit_i_field.ncols() != stokes_dim  )
      throw runtime_error(
        "Obtained *doit_i_field* has wrong number of Stokes elements."  );
    if( doit_i_field.nrows() != 1  )
      throw runtime_error(
        "Obtained *doit_i_field* has wrong number of azimuth angles."  );
    if( doit_i_field.nbooks() != 1  )
      throw runtime_error(
        "Obtained *doit_i_field* has wrong number of longitude points."  );
    if( doit_i_field.nshelves() != 1  )
      throw runtime_error(
        "Obtained *doit_i_field* has wrong number of latitude points."  );
    if( doit_i_field.nvitrines() != cloudbox_limits[1]-cloudbox_limits[0]+1  )
      throw runtime_error(
        "Obtained *doit_i_field* has wrong number of pressure points."  );
    if( doit_i_field.nlibraries() != nf  )
      throw runtime_error(
        "Obtained *doit_i_field* has wrong number of frequency points."  );
  }

  // Reset azimuth grid for scattering source calc later on
  nlinspace(scat_aa_grid, 0, 360, Naa);
  
  // For a brief description of internal variables used, see
  // iyEmissionStandard. 
  
  
  //### jacobian part #########################################################
  // Initialise analytical jacobians (diy_dx and help variables)
  //
  Index           j_analytical_do = 0;
  ArrayOfTensor3  diy_dpath; 
  ArrayOfIndex    jac_species_i(0), jac_scat_i(0), jac_is_t(0), jac_wind_i(0);
  ArrayOfIndex    jac_mag_i(0), jac_other(0), jac_to_integrate(0); 
  // Flags for partial derivatives of propmat
  const PropmatPartialsData ppd(jacobian_quantities);
  //
  if( jacobian_do ) 
  { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }
  //
  if( !j_analytical_do )
  { diy_dx.resize( 0 ); }
  else 
  {
    diy_dpath.resize( nq ); 
    jac_species_i.resize( nq ); 
    jac_scat_i.resize( nq ); 
    jac_is_t.resize( nq ); 
    jac_wind_i.resize( nq ); 
    jac_mag_i.resize( nq ); 
    jac_other.resize(nq);
    jac_to_integrate.resize(nq);
    //
    FOR_ANALYTICAL_JACOBIANS_DO( 
    if( jacobian_quantities[iq].Integration() )
    {
      diy_dpath[iq].resize( 1, nf, ns ); 
      diy_dpath[iq] = 0.0;
    }
    else
    {
      diy_dpath[iq].resize( np, nf, ns ); 
      diy_dpath[iq] = 0.0;
    }
    )
    
    const ArrayOfString scat_species(0);
    
    get_pointers_for_analytical_jacobians( jac_species_i, jac_scat_i, jac_is_t, 
                                           jac_wind_i, jac_mag_i, jac_to_integrate, 
                                           jacobian_quantities,
                                           abs_species, scat_species );
    FOR_ANALYTICAL_JACOBIANS_DO( 
    jac_other[iq] = ppd.is_this_propmattype(iq)?JAC_IS_OTHER:JAC_IS_NONE; 
    if( jac_to_integrate[iq] == JAC_IS_FLUX )
      throw std::runtime_error("This method can not perform flux calculations.\n");
    )
    
    if( iy_agenda_call1 )
    {
      diy_dx.resize( nq ); 
      //
      FOR_ANALYTICAL_JACOBIANS_DO( 
        diy_dx[iq].resize( jacobian_indices[iq][1]-jacobian_indices[iq][0]+1,
                           nf, ns ); 
      diy_dx[iq] = 0.0;
      )
    }
  } 
  //###########################################################################
  
  // Get atmospheric and RT quantities for each ppath point/step
  //
  Vector              ppath_p, ppath_t;
  Matrix              ppath_vmr, ppath_pnd, ppath_wind, ppath_mag;
  Matrix              ppath_f, ppath_t_nlte;
  ArrayOfMatrix       ppath_dpnd_dx;
  Tensor5             dtrans_partial_dx_above, dtrans_partial_dx_below;
  Tensor4             trans_partial, trans_cumulat;
  ArrayOfIndex        clear2cloudy;
  ArrayOfIndex        lte;
  
  if( np > 1 )
  {
    get_ppath_atmvars( ppath_p, ppath_t, ppath_t_nlte, ppath_vmr,
                       ppath_wind, ppath_mag, 
                       ppath, atmosphere_dim, p_grid, t_field, t_nlte_field, 
                       vmr_field, wind_u_field, wind_v_field, wind_w_field,
                       mag_u_field, mag_v_field, mag_w_field );      
    
    get_ppath_f( ppath_f, ppath, f_grid,  atmosphere_dim, 
                 rte_alonglos_v, ppath_wind );
    
    // here, the cloudbox is on, ie we don't need to check and branch this here
    // anymore.
    get_ppath_cloudvars( clear2cloudy, ppath_pnd, ppath_dpnd_dx,
                         ppath, atmosphere_dim, cloudbox_limits,
                         pnd_field, dpnd_field_dx );
    
    PropagationMatrix K_this, K_past, Kp(nf, stokes_dim);
    StokesVector S(nf, stokes_dim), Sp(nf, stokes_dim), a(nf, stokes_dim);
    lte.resize(np);
    ArrayOfPropagationMatrix dK_this_dx(nq), dK_past_dx(nq), dKp_dx(nq);
    ArrayOfStokesVector da_dx(nq), dS_dx(nq), dSp_dx(nq);

    trans_cumulat.resize(np, nf, stokes_dim, stokes_dim);
    trans_partial.resize(np, nf, stokes_dim, stokes_dim);
    dtrans_partial_dx_above.resize(np, nq, nf, stokes_dim, stokes_dim);
    dtrans_partial_dx_below.resize(np, nq, nf, stokes_dim, stokes_dim);
    Vector B(nf), dB_dT(nf);
    Tensor3 J(np, nf, stokes_dim);
    Tensor4 dJ_dx(np, nq, nf, stokes_dim);
    
    for(Index ip = 0; ip < np; ip++)
    {
      get_stepwise_blackbody_radiation(B, dB_dT,
                                       ppath_f(joker, ip), ppath_t[ip],
                                       ppd.do_temperature());
      
      get_stepwise_clearsky_propmat(ws,
                                    K_this,
                                    S,
                                    lte[ip],
                                    dK_this_dx,
                                    dS_dx,
                                    propmat_clearsky_agenda, 
                                    jacobian_quantities,
                                    ppd,
                                    ppath_f(joker, ip),
                                    ppath_mag(joker, ip),
                                    ppath.los(ip, joker),
                                    ppath_t_nlte.nrows()?
                                    ppath_t_nlte(joker, ip):
                                    Vector(0),
                                    ppath_vmr.nrows()?
                                    ppath_vmr(joker, ip):
                                    Vector(0),
                                    ppath_t[ip],
                                    ppath_p[ip],
                                    jac_species_i,
                                    jacobian_do);

      get_stepwise_scattersky_propmat(a, Kp, da_dx, dKp_dx,
                                      clear2cloudy[ip]+1,
                                      ppath_pnd(joker, ip),
                                      ppath_dpnd_dx,
                                      ip,
                                      scat_data,
                                      ppath.los(ip, joker),
                                      ppath_t[ip],
                                      atmosphere_dim,
                                      jacobian_do,
                                      verbosity);

      a += K_this;
      K_this += Kp;
      for( Index iq = 0; iq < nq; iq++ )
      {
        da_dx[iq] += dK_this_dx[iq];
        dK_this_dx[iq] += dKp_dx[iq];
      }
      
      get_stepwise_transmission_matrix(trans_cumulat(ip, joker, joker, joker),
                                       trans_partial(ip, joker, joker, joker),
                                       nq?dtrans_partial_dx_above(ip, joker, joker, 
                                           joker, joker):Tensor4(0,0,0,0),
                                       nq?dtrans_partial_dx_below(ip, joker, joker, 
                                           joker, joker):Tensor4(0,0,0,0),
                                       (ip > 0)?
                                       trans_cumulat(ip-1, joker, joker, joker):
                                       Tensor3(0,0,0),
                                       K_past,
                                       K_this,
                                       dK_past_dx,
                                       dK_this_dx,
                                       (ip > 0)?ppath.lstep[ip-1]:Numeric(1.0),
                                       ip==0);

      get_stepwise_scattersky_source(Sp, dSp_dx,
                                     clear2cloudy[ip]+1,
                                     ppath_pnd(joker, ip),
                                     ppath_dpnd_dx,
                                     ip,
                                     scat_data,
                                     scat_data_checked,
                                     doit_i_field,
                                     scat_za_grid,
                                     scat_aa_grid,
                                     ppath.los(ip, joker),
                                     ppath.gp_p[ip],
                                     ppath_t[ip],
                                     atmosphere_dim,
                                     jacobian_do,
                                     verbosity);
      S += Sp;
      for( Index iq = 0; iq < nq; iq++ )
      {
        dS_dx[iq] += dSp_dx[iq];
      }

      get_stepwise_effective_source(J(ip, joker, joker),
                                    nq?dJ_dx(ip, joker, joker, joker):Tensor3(0,0,0),
                                    K_this,
                                    a,
                                    S,
                                    dK_this_dx,
                                    ArrayOfStokesVector(0),  // da_dx,
                                    dS_dx,
                                    B,
                                    dB_dT,
                                    jacobian_quantities);
      
      swap(K_past, K_this);
      swap(dK_past_dx, dK_this_dx);
      
    }
    
    // Get *iy* at end ppath by interpoling doit_i_field
    {
      Tensor4 i_field = doit_i_field(joker,joker,0,0,joker,0,joker);
      
      ArrayOfGridPos gp_p(1), gp_za(1);
      gridpos_copy( gp_p[0], ppath.gp_p[np-1] );
      gridpos( gp_za, scat_za_grid, Vector(1,ppath.los(np-1,0)) );
      
      Tensor3 itw( 1, 1, 4 );
      interpweights( itw, gp_p, gp_za );
      
      iy.resize( nf, stokes_dim );
      
      Matrix m(1,1);
      for( Index f=0; f<nf; f++ )
      {
        for( Index s=0; s<stokes_dim; s++ )
        {
          interp( m, itw, doit_i_field(f,joker,0,0,joker,0,s), gp_p, gp_za );
          iy(f,s) = m(0,0);
        }
      }
    }
    
    // Radiative transfer solver
    for(Index iv = 0; iv < nf; iv++)
    {
      Matrix ImT(stokes_dim, stokes_dim);
      Vector one(stokes_dim), two(stokes_dim), three(stokes_dim);
      
      for(Index ip = np - 2; ip >= 0; ip--)
      {
        
        MatrixView T = trans_partial(ip+1, iv, joker, joker);
        
        mult(one, T, iy(iv, joker));
        id_mat(ImT);
        ImT -= T;
        three = J(ip, iv, joker);
        three += J(ip+1, iv, joker);
        three *= 0.5;
        
        mult(two, ImT, three);
        iy(iv, joker) = one;
        iy(iv, joker) += two;
      }
    }
  }
  else
  {  
    trans_cumulat.resize( np, nf, ns, ns );
    for( Index iv=0; iv<nf; iv++ )
    { id_mat( trans_cumulat(np-1, iv,joker,joker) ); }
  }
  
  // Unit conversions
  if( iy_agenda_call1 )
  {
    // Determine refractive index to use for the n2 radiance law
    Numeric n = 1.0; // First guess is that sensor is in space
    //
    if( ppath.end_lstep == 0 ) // If true, sensor inside the atmosphere
    { n = ppath.nreal[np-1]; }
    
    // Polarisation index variable
    ArrayOfIndex i_pol(ns);
    for( Index is=0; is<ns; is++ )
    { i_pol[is] = is + 1; }
    
    // Jacobian part (must be converted to Tb before iy for PlanckBT)
    // 
    if( j_analytical_do )
    {
      FOR_ANALYTICAL_JACOBIANS_DO( apply_iy_unit2( diy_dx[iq], iy, iy_unit,
                                                    f_grid, n, i_pol ); )
    } 
    
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
