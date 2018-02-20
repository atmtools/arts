/* Copyright (C) 2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler(at)ltu.se>

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
  ===  File description
  ===========================================================================*/

/*!
  \file   m_cloudradar.cc
  \author Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2010-10-31

  \brief  Workspace functions related to simulation of cloud radars.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "auto_md.h"
#include "complex.h"
#include "logic.h"
#include "messages.h"
#include "montecarlo.h"
#include "refraction.h"
#include "rte.h"
#include "sensor.h"

extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT;
extern const String SCATSPECIES_MAINTAG;



/* Workspace method: Doxygen documentation will be auto-generated */
void iyActiveSingleScat(
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
   const ArrayOfTensor4&              dpnd_field_dx,
   const ArrayOfString&               scat_species,
   const ArrayOfArrayOfSingleScatteringData& scat_data,
   const Matrix&                      particle_masses,
   const String&                      iy_unit,
   const ArrayOfString&               iy_aux_vars,
   const Index&                       jacobian_do,
   const ArrayOfRetrievalQuantity&    jacobian_quantities,
   const Agenda&                      ppath_agenda,
   const Agenda&                      propmat_clearsky_agenda,
   const Agenda&                      iy_transmitter_agenda,
   const Index&                       iy_agenda_call1,
   const Tensor3&                     iy_transmission,
   const Vector&                      rte_pos,      
   const Vector&                      rte_los,      
   const Numeric&                     rte_alonglos_v,      
   const Numeric&                     ppath_lmax,
   const Numeric&                     ppath_lraytrace,
   const Numeric&                     ze_tref,
   const Numeric&                     k2,
   const Verbosity&                   verbosity )
{
  // Throw error if unsupported features are requested
  if( !iy_agenda_call1 )
    throw runtime_error( 
                  "Recursive usage not possible (iy_agenda_call1 must be 1)." );
  if( !iy_transmission.empty() )
    throw runtime_error( "*iy_transmission* must be empty." );

  
  // Determine propagation path
  //
  ppath_agendaExecute( ws, ppath, ppath_lmax, ppath_lraytrace, rte_pos, rte_los,
                       Vector(0), 0, 0, t_field, z_field, vmr_field, 
                       f_grid, ppath_agenda );

  // Some basic sizes
  //
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index nq = jacobian_quantities.nelem();
  const Index nsetot = pnd_field.nbooks();

      
  //###### jacobian part #######################################################
  // Initialise analytical jacobians (diy_dx and help variables)
  //
  Index           j_analytical_do = 0;
  ArrayOfTensor3  diy_dpath; 
  ArrayOfIndex    jac_species_i(0), jac_scat_i(0), jac_is_t(0), jac_wind_i(0);
  ArrayOfIndex    jac_mag_i(0), jac_other(0); 
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
      //
      FOR_ANALYTICAL_JACOBIANS_DO( 
        if( jacobian_quantities[iq].Integration() )
        {
            diy_dpath[iq].resize( 1, nf*np, ns ); 
            diy_dpath[iq] = 0.0;
        }
        else
        {
            diy_dpath[iq].resize( np, nf*np, ns ); 
            diy_dpath[iq] = 0.0;
        }
      )
        
      get_pointers_for_analytical_jacobians( jac_species_i, jac_scat_i, jac_is_t, 
                                             jac_wind_i, jac_mag_i, 
                                             jacobian_quantities,
                                             abs_species, scat_species );
      
      FOR_ANALYTICAL_JACOBIANS_DO( 
      jac_other[iq] = ppd.is_this_propmattype(iq)?Index(JacobianType::Other):Index(JacobianType::None); 
        if( jac_scat_i[iq]+1 )
        {
          if( dpnd_field_dx[iq].empty() )
            throw std::runtime_error("*dpnd_field_dx* not allowed to be empty"
                                     " for scattering Jacobian species.\n");
        }
        // FIXME: should we indeed check for that? remove if it causes issues.
        else
        {
          if( !dpnd_field_dx[iq].empty() )
            throw std::runtime_error("*dpnd_field_dx* must be empty for"
                                     " non-scattering Jacobian species.\n");
        }
      )

      if( iy_agenda_call1 )
        {
          diy_dx.resize( nq ); 
          //
          bool any_affine;
          ArrayOfArrayOfIndex jacobian_indices;
          jac_ranges_indices( jacobian_indices, any_affine,
                              jacobian_quantities, true );
          //
          FOR_ANALYTICAL_JACOBIANS_DO( diy_dx[iq].resize( 
            jacobian_indices[iq][1]-jacobian_indices[iq][0]+1, nf*np, ns ); 
            diy_dx[iq] = 0.0;
          )
        }
    } 
  //###########################################################################
    
  
  //=== iy_aux part ===========================================================
  Index auxPressure     = -1,
        auxTemperature  = -1,
        auxBackScat     = -1,
        auxTrans        = -1,
        auxRoTrTime     = -1;
  ArrayOfIndex auxPartCont(0), auxPartContI(0);
  ArrayOfIndex auxPartField(0), auxPartFieldI(0);
  //
  const Index naux = iy_aux_vars.nelem();
  iy_aux.resize( naux );
  //
  for( Index i=0; i<naux; i++ )
    {
      if( iy_aux_vars[i] == "Pressure" )
        { auxPressure = i;      iy_aux[i].resize( 1, 1, 1, np ); }
      else if( iy_aux_vars[i] == "Temperature" )
        { auxTemperature = i;   iy_aux[i].resize( 1, 1, 1, np ); }
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
      else if( iy_aux_vars[i] == "Backscattering" )
        { auxBackScat = i;  iy_aux[i].resize( nf, ns, 1, np ); iy_aux[i] = 0; }
      else if( iy_aux_vars[i] == "Transmission" )
        { auxTrans = i;  iy_aux[i].resize( nf, ns, ns, np ); }
      else if( iy_aux_vars[i] == "Round-trip time" )
        { auxRoTrTime = i;   iy_aux[i].resize( 1, 1, 1, np ); }
      else
        {
          ostringstream os;
          os << "In *iy_aux_vars* you have included: \"" << iy_aux_vars[i]
             << "\"\nThis choice is not recognised.";
          throw runtime_error( os.str() );
        }
    }
  //===========================================================================


  // Size iy and set to zero
  iy.resize( nf*np, ns );
  iy = 0;

  // If neither cloudbox or aux variables, nothing more to do
  if( !cloudbox_on  &&  !naux )
    { return; }
  
  // Transmitted signal
  //
  Matrix iy0;
  //
  iy_transmitter_agendaExecute( ws, iy0, f_grid, rte_pos, rte_los,
                                                       iy_transmitter_agenda ); 
  if( iy0.ncols() != stokes_dim  ||  iy0.nrows() != nf )
    throw runtime_error( "The size of *iy* returned from "
                                   "*iy_transmitter_agenda* is not correct." );
  for( Index iv=0; iv<nf; iv++ )
    {
      if( iy0(iv,0) != 1 )
        throw runtime_error( "The *iy* returned from *iy_transmitter_agenda* "
                             "must have the value 1 in the first column." );
    }

  // Get atmospheric and RT quantities for each ppath point/step
  //
  ArrayOfIndex        iaps(0);
  //
  Vector              ppath_p, ppath_t;
  Matrix              ppath_vmr, ppath_pnd, ppath_wind, ppath_mag;
  Matrix              ppath_f, ppath_t_nlte;
  ArrayOfArrayOfPropagationMatrix abs_per_species;
  Tensor5             dtrans_partial_dx_above, dtrans_partial_dx_below;
  ArrayOfPropagationMatrix ppath_ext, pnd_ext_mat;
  Tensor3             pnd_abs_vec;
  Tensor4 trans_partial, trans_cumulat;
  Vector              scalar_tau;
  ArrayOfIndex        clear2cloudy;
  ArrayOfArrayOfIndex extmat_case;   
  ArrayOfArrayOfPropagationMatrix dppath_ext_dx;
  ArrayOfArrayOfStokesVector dppath_nlte_dx, dppath_nlte_source_dx;
  ArrayOfStokesVector        ppath_nlte_source;
  ArrayOfIndex        lte;  
  ArrayOfMatrix       ppath_dpnd_dx;
  const Tensor4       t_nlte_field_empty(0,0,0,0);
  Array<ArrayOfArrayOfSingleScatteringData> scat_data_single;
  //
  if( np > 1 )
    {
      get_ppath_atmvars( ppath_p, ppath_t, ppath_t_nlte, ppath_vmr,
                         ppath_wind, ppath_mag, 
                         ppath, atmosphere_dim, p_grid, t_field,
                         t_nlte_field_empty, vmr_field,
                         wind_u_field, wind_v_field, wind_w_field,
                         mag_u_field, mag_v_field, mag_w_field );
      
      get_ppath_f( ppath_f, ppath, f_grid,  atmosphere_dim, 
                   rte_alonglos_v, ppath_wind );

      if( jacobian_do )
        {
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
                               jac_other, iaps,
                               scat_data,
                               pnd_field, dpnd_field_dx,
                               cloudbox_limits,
                               atmosphere_dim, stokes_dim,
                               jacobian_do, cloudbox_on, verbosity );
        }
      else
        {
          get_ppath_pmat( ws, ppath_ext, ppath_nlte_source, lte, 
                      abs_per_species, 
                      dppath_ext_dx, dppath_nlte_dx,
                      propmat_clearsky_agenda, jacobian_quantities, ppath, 
                      ppath_p, ppath_t, ppath_t_nlte, ppath_vmr, ppath_f, ppath_mag,
                      f_grid, stokes_dim, ArrayOfIndex(0) );
      
          // Extract basic scattering data
          get_ppath_cloudvars( clear2cloudy, ppath_pnd, ppath_dpnd_dx,
                       ppath, atmosphere_dim, cloudbox_limits,
                       pnd_field, dpnd_field_dx );
          get_ppath_partopt( pnd_abs_vec, pnd_ext_mat, scat_data_single,
                             clear2cloudy, ppath_pnd, ppath, ppath_t,
                             stokes_dim, ppath_f, atmosphere_dim, scat_data,
                             verbosity );

          get_ppath_trans2( trans_partial, extmat_case, trans_cumulat, scalar_tau, 
                        ppath, ppath_ext, f_grid, stokes_dim, 
                        clear2cloudy, pnd_ext_mat );
        }
    }


  // Transmission for reversed direction
  Tensor3 tr_rev( nf, ns, ns );
  for( Index iv=0; iv<nf; iv++ ) 
    { id_mat( tr_rev(iv,joker,joker) ); }


  // Loop ppath steps
  for( Index ip=0; ip<np; ip++ )
    {
      // Direction of outgoing scattered radiation (which is reverse to LOS).
      Vector los_sca;
      mirror_los( los_sca, ppath.los(ip,joker), atmosphere_dim );

      // Obtain a length-2 vector for incoming direction
      Vector los_inc;
      if( atmosphere_dim == 3 )
        { los_inc = ppath.los(ip,joker); }
      else // Mirror back to get a correct 3D LOS
        { mirror_los( los_inc, los_sca, 3 ); }

      
      // Radar return only possible if inside cloudbox
      if( clear2cloudy[ip] >= 0 )
        {
          for( Index iv=0; iv<nf; iv++ )
            {
              // Obtain scattering matrix data
              //
              Matrix  P(ns,ns);
              Tensor3 Pe(nsetot,ns,ns);
              //
              if( !j_analytical_do )
                {
                  // Here we just need the total matrix
                  pha_mat_singleCalc( P, los_sca[0], los_sca[1],
                                      los_inc[0], los_inc[1],
                                      scat_data_single[iv],
                                      stokes_dim, ppath_pnd(joker,ip),
                                      ppath_t[ip], verbosity );
                }
              else
                {
                  // Here we need the matrix per scattering element. 
                  // We extract data for pnd=1 (but only when either
                  // ppath_pnd or ppath_dpnd_dx is non-zero
                  //
                  Vector pnd_probe(nsetot,0);
                  for( Index i=0; i<nsetot; i++ )
                    {
                      if( ppath_pnd(i,ip) != 0 )
                        {
                          pnd_probe[i] = 1;
                        }
                      else
                        {
                          for( Index iq=0; iq<nq && !pnd_probe[i]; iq++ ) 
                            {
                              if( jac_scat_i[iq] >= 0 )
                                {
                                  if( ppath_dpnd_dx[iq](i,ip) != 0 )
                                    {
                                      pnd_probe[i] = 1;
                                      break;
                                    }
                                }
                            }
                        }
                    }
                  //
                  pha_mat_singleCalcScatElement( Pe, los_sca[0], los_sca[1],
                                                 los_inc[0], los_inc[1],
                                                 scat_data_single[iv],
                                                 stokes_dim, pnd_probe,
                                                 ppath_t[ip], verbosity );

                  // Total scattering matrix
                  P = 0.0;
                  for( Index i=0; i<nsetot; i++ ) {
                    if( ppath_pnd(i,ip) != 0 ) {
                      for( Index is1=0; is1<ns; is1++ ) {
                        for( Index is2=0; is2<ns; is2++ )
                          { P(is1,is2) += ppath_pnd(i,ip) * Pe(i,is1,is2); }
                    } } }
                }
              
              // Combine iy0, double transmission and scattering matrix
              Vector iy1(ns), iy2(ns);
              const Index iout = iv*np + ip;
              mult( iy1, tr_rev(iv,joker,joker), iy0(iv,joker) );
              mult( iy2, P, iy1 );
              mult( iy(iout,joker), trans_cumulat(iv,joker,joker,ip), iy2 );
              
              //=== iy_aux part ===========================================
              // Backscattering
              if( auxBackScat >= 0 ) {
                mult( iy_aux[auxBackScat](iv,joker,0,ip), P, iy0(iv,joker) ); }
              
              // Jacobians
              if( j_analytical_do )
                {
                  for( Index iq=0; iq<nq; iq++ ) 
                    {
                      if( jacobian_quantities[iq].Analytical() )
                        {
                          // Impact on back-scattering.
                          // Only non-zero for scattering species
                          if( jac_scat_i[iq] >= 0 )
                            {
                              // Change of scattering scattering matrix
                              P = 0.0;
                              for( Index i=0; i<nsetot; i++ ) {
                                if( ppath_dpnd_dx[iq](i,ip) != 0 ) {
                                  for( Index is1=0; is1<ns; is1++ ) {
                                    for( Index is2=0; is2<ns; is2++ )
                                      { P(is1,is2) += ppath_dpnd_dx[iq](i,ip) *
                                                                Pe(i,is1,is2); }
                               } } }

                              // Apply transmissions as above
                              Vector iy_tmp(ns);
                              mult( iy_tmp, P, iy1 );
                              mult( diy_dpath[iq](ip,iout,joker),
                                    trans_cumulat(iv,joker,joker,ip), iy_tmp );
                            }
                        }

                      // Impact on attenuation
                      if( ip > 0  &&  jac_species_i[iq] >= 0 )
                        {
                          Vector iy_tmp0(ns), iy_tmp1(ns);
                          for( Index ip2=ip-1; ip2>=0; ip2-- )
                            {
                              mult( iy_tmp0,
                                    dtrans_partial_dx_below(iq,iv,joker,joker,ip2),
                                    iy2 );
                              mult( iy_tmp1,
                                    dtrans_partial_dx_above(iq,iv,joker,joker,ip2+1),
                                    iy2 );
                              for( Index is1=0; is1<ns; is1++ )
                                {
                                  diy_dpath[iq](ip2,iout,is1)   += 2*iy_tmp0[is1];
                                  diy_dpath[iq](ip2+1,iout,is1) += 2*iy_tmp1[is1];
                                }
                            }
                        }
                    }
                }  // j_analytical_do
            }
        }  // if inside cloudbox

      
      // Stuff to do even outside cloudbox:      
      
      //=== iy_aux part ===================================================
      // Pressure
      if( auxPressure >= 0 ) 
        { iy_aux[auxPressure](0,0,0,ip) = ppath_p[ip]; }
      // Temperature
      if( auxTemperature >= 0 ) 
        { iy_aux[auxTemperature](0,0,0,ip) = ppath_t[ip]; }
      // Particle mass content
      for( Index j=0; j<auxPartCont.nelem(); j++ )
        { iy_aux[auxPartCont[j]](0,0,0,ip) = ppath_pnd(joker,ip) *
                                      particle_masses(joker,auxPartContI[j]); }
      // Particle number density
      for( Index j=0; j<auxPartField.nelem(); j++ )
        { iy_aux[auxPartField[j]](0,0,0,ip) = ppath_pnd(auxPartFieldI[j],ip); }
      if( auxRoTrTime >= 0 ) 
        { 
          if( ip == 0 )
            { iy_aux[auxRoTrTime](0,0,0,ip) = 2 * ppath.end_lstep / 
                                                              SPEED_OF_LIGHT; }
          else
            { iy_aux[auxRoTrTime](0,0,0,ip) = iy_aux[auxRoTrTime](0,0,0,ip-1)+
                ppath.lstep[ip-1] * ( ppath.ngroup[ip-1]+ppath.ngroup[ip] ) /
                                                              SPEED_OF_LIGHT; }
        }
      // Transmission
      if( auxTrans >= 0 ) 
        { for( Index iv=0; iv<nf; iv++ ) {
            for( Index is1=0; is1<ns; is1++ ){
              for( Index is2=0; is2<ns; is2++ ){
                iy_aux[auxTrans](iv,is1,is2,ip) = 
                  tr_rev(iv,is1,is2); } } } }
      //===================================================================

      // Update tr_rev
      if( ip<np-1 )
        {
          for( Index iv=0; iv<nf; iv++ )
            {
              Matrix tmp = tr_rev(iv,joker,joker);
              mult( tr_rev(iv,joker,joker), 
                    trans_partial(iv,joker,joker,ip), tmp );
            }
        }
    } // for ip


  //### jacobian part #####################################################
  if( j_analytical_do )
    {
      // Map to retrieval grids
      FOR_ANALYTICAL_JACOBIANS_DO( 
        diy_from_path_to_rgrids( diy_dx[iq], jacobian_quantities[iq], 
                                 diy_dpath[iq], atmosphere_dim, ppath, ppath_p );
      )
    }
  //#######################################################################
  
  
  if( iy_unit == "1" )
    {}
  else if( iy_unit == "Ze" )
    {
      
      // Refractive index for water (if needed)
      Matrix complex_n(0,0);
      if( k2 <= 0 )
        { complex_n_water_liebe93( complex_n, f_grid, ze_tref ); }
          
      // Common conversion factor
      const Numeric a = 4e18/(PI*PI*PI*PI);

      for( Index iv=0; iv<nf; iv++ )
        {
          // Reference dielectric factor.
          Numeric K2;
          if( k2 > 0 )
            { K2 = k2; }
          else
            {
              Complex n( complex_n(iv,0), complex_n(iv,1) );
              Complex n2 = n*n;
              Complex K  = ( n2 - Numeric(1.0) ) / ( n2 + Numeric(2.0) );
              Numeric absK = abs( K );
              K2 = absK * absK;
            }

          // Conversion factor
          Numeric la = SPEED_OF_LIGHT / f_grid[iv];
          Numeric fac = a*la*la*la*la / K2; 
      
          iy(Range(iv*np,np),joker) *= fac;

          // Jacobian part
          // 
          if( j_analytical_do )
            {
              FOR_ANALYTICAL_JACOBIANS_DO( diy_dx[iq] *= fac; )
            }
          
          //=== iy_aux part ===========================================
          // Backscattering
          if( auxBackScat >= 0 ) {
            iy_aux[auxBackScat](iv,joker,joker,joker) *= fac; }

        }
    }
  else
    {
      throw runtime_error( "For this method, *iy_unit* must be set to \"1\" "
                                                                "or \"Ze\"." );
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void iyActiveSingleScat2(
        Workspace&                          ws,
        Matrix&                             iy,
        ArrayOfMatrix&                      iy_aux,
        ArrayOfTensor3&                     diy_dx,
        Vector&                             ppvar_p,
        Vector&                             ppvar_t,
        Matrix&                             ppvar_t_nlte,
        Matrix&                             ppvar_vmr,
        Matrix&                             ppvar_wind,
        Matrix&                             ppvar_mag,         
        Matrix&                             ppvar_pnd,
        Matrix&                             ppvar_f,  
        Tensor3&                            ppvar_iy,  
  const Index&                              stokes_dim,
  const Vector&                             f_grid,
  const Index&                              atmosphere_dim,
  const Vector&                             p_grid,
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
  const ArrayOfString&                      scat_species,
  const ArrayOfArrayOfSingleScatteringData& scat_data,
  const ArrayOfString&                      iy_aux_vars,
  const Index&                              jacobian_do,
  const ArrayOfRetrievalQuantity&           jacobian_quantities,
  const Ppath&                              ppath,
  const Agenda&                             propmat_clearsky_agenda,
  const Agenda&                             iy_transmitter_agenda,
  const Index&                              iy_agenda_call1,
  const Tensor3&                            iy_transmission,
  const Numeric&                            rte_alonglos_v,
  const Verbosity& )
{
  // Some basic sizes
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index ne = pnd_field.nbooks();

  // Radiative background index
  const Index rbi = ppath_what_background( ppath );
  
  // Checks of input
  // Throw error if unsupported features are requested
  if( !iy_agenda_call1 )
    throw runtime_error( 
                  "Recursive usage not possible (iy_agenda_call1 must be 1)" );
  if( !iy_transmission.empty() )
    throw runtime_error( "*iy_transmission* must be empty" );
  if( rbi < 1  ||  rbi > 9 )
    throw runtime_error( "ppath.background is invalid. Check your "
                         "calculation of *ppath*?" );
  if( jacobian_do )
    {
      if( dpnd_field_dx.nelem() != jacobian_quantities.nelem() )
        throw runtime_error(
          "*dpnd_field_dx* not properly initialized:\n"
          "Number of elements in dpnd_field_dx must be equal number of jacobian"
          " quantities.\n(Note: jacobians have to be defined BEFORE *pnd_field*"
          " is calculated/set." );
    }
  // iy_aux_vars checked below

  // Transmitted signal
  //
  Matrix iy0;
  //
  iy_transmitter_agendaExecute( ws, iy0, f_grid, 
                                ppath.pos(np-1,Range(0,atmosphere_dim)),
                                ppath.los(np-1,joker), iy_transmitter_agenda );
  if( iy0.ncols() != ns  ||  iy0.nrows() != nf )
    {
      ostringstream os;
      os << "The size of *iy* returned from *iy_transmitter_agenda* is\n"
         << "not correct:\n"
         << "  expected size = [" << nf << "," << stokes_dim << "]\n"
         << "  size of iy    = [" << iy0.nrows() << "," << iy0.ncols()<< "]\n";
      throw runtime_error( os.str() );      
    }
  for( Index iv=0; iv<nf; iv++ )
    {
      if( iy0(iv,0) != 1 )
        throw runtime_error( "The *iy* returned from *iy_transmitter_agenda* "
                             "must have the value 1 in the first column." );
    }
  
  
  //  Init Jacobian quantities
  Index   j_analytical_do = 0;
  if( jacobian_do ) 
    { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }
  //
  const Index     nq = j_analytical_do ? jacobian_quantities.nelem() : 0;
  ArrayOfTensor3  diy_dpath(nq); 
  ArrayOfIndex    jac_species_i(nq), jac_scat_i(nq), jac_is_t(nq), jac_wind_i(nq);
  ArrayOfIndex    jac_mag_i(nq), jac_other(nq);
  //
  // Flags for partial derivatives of propmat
  const PropmatPartialsData ppd(jacobian_quantities);
  //
  if( j_analytical_do )
    {
      rtmethods_jacobian_init( jac_species_i, jac_scat_i, jac_is_t, jac_wind_i,
                               jac_mag_i, jac_other, diy_dx, diy_dpath,
                               ns, nf, np, nq, abs_species,
                               scat_species, dpnd_field_dx, ppd,
                               jacobian_quantities, iy_agenda_call1 );
    }
  
  // Init iy_aux and fill where possible
  const Index naux = iy_aux_vars.nelem();
  iy_aux.resize( naux );
  //
  for( Index i=0; i<naux; i++ )
    {
      iy_aux[i].resize(nf,ns); 
      
      if( iy_aux_vars[i] == "Radiative background" )
        { iy_aux[i] = (Numeric)min( (Index)2, rbi-1 ); }
      else
        {
          ostringstream os;
          os << "The only allowed strings in *iy_aux_vars* are:\n"
             << "  \"Radiative background\"\n"
             << "but you have selected: \"" << iy_aux_vars[i] << "\"";
          throw runtime_error( os.str() );      
        }
    }

  // Get atmospheric and radiative variables along the propagation path
  //
  Tensor3 J(np,nf,ns);
  Tensor4 trans_cumulat(np,nf,ns,ns), trans_partial(np,nf,ns,ns);
  Tensor5 dtrans_partial_dx_above(np,nq,nf,ns,ns);
  Tensor5 dtrans_partial_dx_below(np,nq,nf,ns,ns);
  Tensor5 Pe( ne, np, nf, ns, ns, 0 ); 
  ArrayOfMatrix ppvar_dpnd_dx(0);
  ArrayOfIndex clear2cloudy;
  //
  if( np == 1  &&  rbi == 1 )  // i.e. ppath is totally outside the atmosphere:
    {
      ppvar_p.resize(0);
      ppvar_t.resize(0);
      ppvar_vmr.resize(0,0);
      ppvar_t_nlte.resize(0,0);
      ppvar_wind.resize(0,0);
      ppvar_mag.resize(0,0);
      ppvar_pnd.resize(0,0);
      ppvar_f.resize(0,0);
      ppvar_iy.resize(0,0,0);
    }
  else
    {
      // ppvar_iy
      ppvar_iy.resize(nf,ns,np);
      ppvar_iy(joker,joker,np-1) = iy;
      
      // Basic atmospheric variables
      get_ppath_atmvars( ppvar_p, ppvar_t, ppvar_t_nlte, ppvar_vmr,
                         ppvar_wind, ppvar_mag, 
                         ppath, atmosphere_dim, p_grid, t_field, t_nlte_field, 
                         vmr_field, wind_u_field, wind_v_field, wind_w_field,
                         mag_u_field, mag_v_field, mag_w_field );
      
      get_ppath_f( ppvar_f, ppath, f_grid,  atmosphere_dim, 
                   rte_alonglos_v, ppvar_wind );

      // pnd_field
      if( cloudbox_on )
        {
          get_ppath_cloudvars( clear2cloudy, ppvar_pnd, ppvar_dpnd_dx,
                               ppath, atmosphere_dim, cloudbox_limits,
                               pnd_field, dpnd_field_dx );
        }
      else
        {
          clear2cloudy.resize(np);
          for( Index ip=0; ip<np; ip++ )
            { clear2cloudy[ip] = -1; }
        }
      
      // Size radiative variables always used
      PropagationMatrix K_this(nf,ns), K_past(nf,ns), Kp(nf,ns);
      StokesVector a(nf,ns), S(nf,ns), Sp(nf,ns);
      ArrayOfIndex lte(np);

      // Init variables only used if analytical jacobians done
      Vector dB_dT(0);
      ArrayOfPropagationMatrix dK_this_dx(nq), dK_past_dx(nq), dKp_dx(nq);
      ArrayOfStokesVector da_dx(nq), dS_dx(nq), dSp_dx(nq);
      //
      if( j_analytical_do )
        {
          dB_dT.resize(nf);
          FOR_ANALYTICAL_JACOBIANS_DO
            (
              dK_this_dx[iq] = PropagationMatrix(nf,ns);
              dK_past_dx[iq] = PropagationMatrix(nf,ns);
              dKp_dx[iq]     = PropagationMatrix(nf,ns);
              dS_dx[iq]      = StokesVector(nf,ns);
              dSp_dx[iq]     = StokesVector(nf,ns);
            )
        }

      // Loop ppath points and determine radiative properties
      for( Index ip=0; ip<np; ip++ )
        { 
          get_stepwise_clearsky_propmat( ws,
                                         K_this,
                                         S,
                                         lte[ip],
                                         dK_this_dx,
                                         dS_dx,
                                         propmat_clearsky_agenda,
                                         jacobian_quantities,
                                         ppd,
                                         ppvar_f(joker,ip),
                                         ppvar_mag(joker,ip),
                                         ppath.los(ip,joker),
                                         ppvar_t_nlte(joker,ip),
                                         ppvar_vmr(joker,ip),
                                         ppvar_t[ip],
                                         ppvar_p[ip],
                                         jac_species_i,
                                         j_analytical_do );

          if( j_analytical_do )
            {
              adapt_stepwise_partial_derivatives( dK_this_dx,
                                                  dS_dx,
                                                  jacobian_quantities,
                                                  ppvar_f(joker,ip),
                                                  ppath.los(ip,joker),
                                                  ppvar_vmr(joker,ip),
                                                  ppvar_t[ip],
                                                  ppvar_p[ip],
                                                  jac_species_i,
                                                  jac_wind_i,
                                                  lte[ip],
                                                  atmosphere_dim,
                                                  j_analytical_do );
            }

          if( clear2cloudy[ip]+1 )
            {
              get_stepwise_scattersky_propmat( a,
                                               Kp,
                                               da_dx,
                                               dKp_dx,
                                               jacobian_quantities,
                                               ppvar_pnd(joker,Range(ip,1)),
                                               ppvar_dpnd_dx,
                                               ip,
                                               scat_data,
                                               ppath.los(ip,joker),
                                               ppvar_t[Range(ip,1)],
                                               atmosphere_dim,
                                               jacobian_do );
              K_this += Kp;
              
              if( j_analytical_do )
                {
                  FOR_ANALYTICAL_JACOBIANS_DO
                    (
                      dK_this_dx[iq] += dKp_dx[iq];
                    )
                 }

              // Get back-scattering per particle, where relevant
              {
                // Direction of outgoing scattered radiation (which is reverse to LOS).
                Vector los_sca;
                mirror_los( los_sca, ppath.los(ip,joker), atmosphere_dim );

                // Obtain a length-2 vector for incoming direction
                Vector los_inc;
                if( atmosphere_dim == 3 )
                  { los_inc = ppath.los(ip,joker); }
                else // Mirror back to get a correct 3D LOS
                  { mirror_los( los_inc, los_sca, 3 ); }

                // Fill vector with ones where we need back-scattering
                Vector pnd_unit(ne,0);
                for( Index i=0; i<ne; i++ )
                  {
                    if( ppvar_pnd(i,ip) != 0 )
                        { pnd_unit[i] = 1; }
                    else if( j_analytical_do )
                      {
                        for( Index iq=0; iq<nq && !pnd_unit[i]; iq++ ) 
                          {
                            if( jac_scat_i[iq] >= 0 )
                              {
                                if( ppvar_dpnd_dx[iq](i,ip) != 0 )
                                  {
                                    pnd_unit[i] = 1;
                                    break;
                                  }
                              }
                          }
                      }
                  }

                // Exctract back-scattering
                /*
                pha_mat_singleCalcScatElement( Pe, los_sca[0], los_sca[1],
                                               los_inc[0], los_inc[1],
                                               scat_data_single[iv],
                                               stokes_dim, pnd_probe,
                                               ppath_t[ip], verbosity );
                */
              }
              
            }  // clear2cloudy
      
          get_stepwise_transmission_matrix(
                                 trans_cumulat(ip,joker,joker,joker),
                                 trans_partial(ip,joker,joker,joker),
                                 dtrans_partial_dx_above(ip,joker,joker,joker,joker),
                                 dtrans_partial_dx_below(ip,joker,joker,joker,joker),
                                 (ip>0)?
                                   trans_cumulat(ip-1,joker,joker,joker):
                                   Tensor3(0,0,0),
                                 K_past,
                                 K_this,
                                 dK_past_dx,
                                 dK_this_dx,
                                 (ip>0)?
                                   ppath.lstep[ip-1]:
                                   Numeric(1.0),
                                 ip==0 );
      
          swap( K_past, K_this );
          swap( dK_past_dx, dK_this_dx );
        }
    }

  // Transmission for reversed direction
  Tensor3 tr_rev( nf, ns, ns );
  for( Index iv=0; iv<nf; iv++ ) 
    { id_mat( tr_rev(iv,joker,joker) ); }
  
  // Radiative transfer calculations
  for( Index ip=0; ip<np; ip++ )
    {
      // Radar return only possible if inside cloudbox
      if( clear2cloudy[ip] >= 0 )
        {
          for( Index iv=0; iv<nf; iv++ )
            {
              // Calculate bulk back-scattering
              Matrix  P(ns,ns,0);
              for( Index i=0; i<ne; i++ ) {
                if( ppvar_pnd(i,ip) != 0 ) {
                  for( Index is1=0; is1<ns; is1++ ) {
                    for( Index is2=0; is2<ns; is2++ )
                      { P(is1,is2) += ppvar_pnd(i,ip) * Pe(i,ip,iv,is1,is2); }
                  } } }

              // Combine iy0, double transmission and scattering matrix
              Vector iy1(ns), iy2(ns);
              const Index iout = iv*np + ip;
              mult( iy1, tr_rev(iv,joker,joker), iy0(iv,joker) );
              mult( iy2, P, iy1 );
              mult( iy(iout,joker), trans_cumulat(ip,iv,joker,joker), iy2 );

              // Jacobians
              if( j_analytical_do )
                {
                  // All this is to cover impact on back-scattering of
                  // scattering species 
                  for( Index iq=0; iq<nq; iq++ ) 
                    {
                      if( jacobian_quantities[iq].Analytical() )
                        {
                          if( jac_scat_i[iq] >= 0 )
                            {
                              // Change of scattering scattering matrix
                              P = 0.0;
                              for( Index i=0; i<ne; i++ ) {
                                if( ppvar_dpnd_dx[iq](i,ip) != 0 ) {
                                  for( Index is1=0; is1<ns; is1++ ) {
                                    for( Index is2=0; is2<ns; is2++ )
                                      { P(is1,is2) += ppvar_dpnd_dx[iq](i,ip) *
                                          Pe(i,ip,iv,is1,is2); }
                               } } }

                              // Apply transmissions as above
                              Vector iy_tmp(ns);
                              mult( iy_tmp, P, iy1 );
                              mult( diy_dpath[iq](ip,iout,joker),
                                    trans_cumulat(ip,iv,joker,joker), iy_tmp );
                            }
                        }
                    }
                }  // j_analytical_do
            } // for iv
        } // if cloudy

      // Update tr_rev
      if( ip<np-1 )
        {
          for( Index iv=0; iv<nf; iv++ )
            {          
              Matrix tmp = tr_rev(iv,joker,joker);
              mult( tr_rev(iv,joker,joker), trans_partial(ip+1,iv,joker,joker),
                    tmp );
            }
        }
    } // for ip

  
  // Finalize analytical Jacobians
  if( j_analytical_do )
    {
      rtmethods_jacobian_finalisation( diy_dx, diy_dpath,
                                       ns, nf, np, atmosphere_dim, ppath,
                                       ppvar_p, ppvar_t, ppvar_vmr,
                                       iy_agenda_call1, iy_transmission,
                                       jacobian_quantities, jac_species_i,
                                       jac_is_t );
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void yActive(
         Workspace&                ws,
         Vector&                   y,
         Vector&                   y_f,
         ArrayOfIndex&             y_pol,
         Matrix&                   y_pos,
         Matrix&                   y_los,
         ArrayOfVector&            y_aux,
         Matrix&                   y_geo,
         Matrix&                   jacobian,
   const Index&                    atmgeom_checked,
   const Index&                    atmfields_checked,
   const String&                   iy_unit,   
   const ArrayOfString&            iy_aux_vars,
   const Index&                    stokes_dim,
   const Vector&                   f_grid,
   const Index&                    atmosphere_dim,
   const Tensor3&                  t_field,
   const Tensor3&                  z_field,
   const Tensor4&                  vmr_field,
   const Index&                    cloudbox_on,
   const Index&                    cloudbox_checked,
   const Matrix&                   sensor_pos,
   const Matrix&                   sensor_los,
   const Index&                    sensor_checked,
   const Index&                    jacobian_do,     
   const ArrayOfRetrievalQuantity& jacobian_quantities,
   const Agenda&                   iy_main_agenda,
   const Agenda&                   geo_pos_agenda,
   const ArrayOfArrayOfIndex&      instrument_pol_array,
   const Vector&                   range_bins,
   const Verbosity& )
{
  // Important sizes
  const Index npos  = sensor_pos.nrows();
  const Index nbins = range_bins.nelem() - 1;
  const Index nf    = f_grid.nelem();
  const Index ns    = stokes_dim;
  const Index naux  = iy_aux_vars.nelem();

  //---------------------------------------------------------------------------
  // Input checks
  //---------------------------------------------------------------------------

  // Basics
  //
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  //
  if ( f_grid.empty() )
    { throw runtime_error ( "The frequency grid is empty." ); }
  chk_if_increasing ( "f_grid", f_grid );
  if( f_grid[0] <= 0) 
    { throw runtime_error( "All frequencies in *f_grid* must be > 0." ); }
  //
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  if( atmfields_checked != 1 )
    throw runtime_error( "The atmospheric fields must be flagged to have "
                         "passed a consistency check (atmfields_checked=1)." );
  if( atmgeom_checked != 1 )
    throw runtime_error( "The atmospheric geometry must be flagged to have "
                         "passed a consistency check (atmgeom_checked=1)." );
  if( cloudbox_checked != 1 )
    throw runtime_error( "The cloudbox must be flagged to have "
                         "passed a consistency check (cloudbox_checked=1)." );
  if( sensor_checked != 1 )
    throw runtime_error( "The sensor variables must be flagged to have "
                         "passed a consistency check (sensor_checked=1)." );

  // Method specific variables
  bool is_z = max(range_bins) > 1;
  if( !is_increasing( range_bins ) )
    throw runtime_error( "The vector *range_bins* must contain strictly "
                         "increasing values." );
  if( !is_z && min(range_bins) < 0 )
    throw runtime_error( "The vector *range_bins* is not allowed to contain "
                         "negative times." );
  if( instrument_pol_array.nelem() != nf )
    throw runtime_error( "The main length of *instrument_pol_array* must match "
                         "the number of frequencies." );

  
  //---------------------------------------------------------------------------
  // Various  initializations
  //---------------------------------------------------------------------------

  // Conversion from Stokes to instrument_pol
  ArrayOfVector   s2p;
  stokes2pol( s2p, 0.5 );
  ArrayOfIndex npolcum(nf+1); npolcum[0]=0;
  for( Index i=0; i<nf; i++ )
    { 
      npolcum[i+1] = npolcum[i] + instrument_pol_array[i].nelem(); 
      for( Index j=0; j<instrument_pol_array[i].nelem(); j++ )
        {
          if( s2p[instrument_pol_array[i][j]-1].nelem() > stokes_dim )
            throw runtime_error( "Your definition of *instrument_pol_array* " 
                                 "requires a higher value for *stokes_dim*." );
        }
    }

  // Resize and init *y* and *y_XXX*
  //
  const Index ntot = npos * npolcum[nf] * nbins;
  y.resize( ntot );
  y = NAN;
  y_f.resize( ntot );
  y_pol.resize( ntot );
  y_pos.resize( ntot, sensor_pos.ncols() );
  y_los.resize( ntot, sensor_los.ncols() );
  y_geo.resize( ntot, atmosphere_dim );
  y_geo = NAN;   // Will be replaced if relavant data are provided (*geo_pos*)

  // y_aux
  y_aux.resize( naux );
  for( Index i=0; i<naux; i++ )
    { 
      y_aux[i].resize( ntot ); 
      y_aux[i] = NAN; 
    }

  // Jacobian variables
  //
  Index j_analytical_do = 0;
  Index njq             = 0;
  ArrayOfArrayOfIndex jacobian_indices;
  //
  if( jacobian_do )
    {
      bool any_affine;
      jac_ranges_indices( jacobian_indices, any_affine,
                          jacobian_quantities, true );
      
      jacobian.resize( ntot, 
                       jacobian_indices[jacobian_indices.nelem()-1][1]+1 );
      jacobian = 0;
      njq      = jacobian_quantities.nelem();
      //
      FOR_ANALYTICAL_JACOBIANS_DO(
        j_analytical_do  = 1; 
      )
    }
  else
    { jacobian.resize( 0, 0 ); }
  

  //---------------------------------------------------------------------------
  // The calculations
  //---------------------------------------------------------------------------

  // Loop positions
  for( Index p=0; p<npos; p++ )
    {
      // RT part
      Tensor3        iy_transmission(0,0,0);
      ArrayOfTensor3 diy_dx;
      Vector         rte_pos2(0);
      Matrix         iy;
      Ppath          ppath;
      ArrayOfTensor4 iy_aux;
      const Index    iy_id = (Index)1e6*p;
      //
      iy_main_agendaExecute( ws, iy, iy_aux, ppath, diy_dx, 
                             1, iy_unit, iy_transmission, iy_aux_vars,
                             iy_id, cloudbox_on, jacobian_do, t_field,
                             z_field, vmr_field, f_grid, 
                             sensor_pos(p,joker), sensor_los(p,joker), 
                             rte_pos2, iy_main_agenda );

      // Check if path and size OK
      const Index np = ppath.np;
      if( np == 1 )
        throw runtime_error( "A path consisting of a single point found. "
                             "This is not allowed." );
      const bool isupward = ppath.pos(1,0) > ppath.pos(0,0);
      if( !( ( isupward & is_increasing( ppath.pos(joker,0) ) ) ||
             ( !isupward & is_decreasing( ppath.pos(joker,0) ) ) ) )
        throw runtime_error( "A path with strictly increasing or decreasing "
                             "altitudes must be obtained (e.g. limb "
                             "measurements are not allowed)." );
      if( iy.nrows() != nf*np )
        throw runtime_error( "The size of *iy* returned from *iy_main_agenda* "
                             "is not correct (for this method)." );

      // Range of ppath, in altitude or time
      Vector range(np);
      if( is_z )
        { range = ppath.pos(joker,0); }
      else
        { // Calculate round-trip time
          range[0] = 2 * ppath.end_lstep / SPEED_OF_LIGHT;
          for( Index i=1; i<np; i++ )
            { range[i] = range[i-1] + ppath.lstep[i-1] * 
                      ( ppath.ngroup[i-1]+ppath.ngroup[i] ) / SPEED_OF_LIGHT; }
        }
      const Numeric range_end1 = min( range[0], range[np-1] );
      const Numeric range_end2 = max( range[0], range[np-1] );

      // Help variables to calculate jacobians
      ArrayOfTensor3 dI(njq);
      ArrayOfMatrix drefl(njq);
      if( j_analytical_do )
        {
          FOR_ANALYTICAL_JACOBIANS_DO(
            dI[iq].resize( jacobian_indices[iq][1]-jacobian_indices[iq][0]+1, np, ns ); 
            drefl[iq].resize( dI[iq].npages(), np ); 
          )
        }
      
      // Loop radar bins
      for( Index b=0; b<nbins; b++ )
        {
          if( ! ( range_bins[b] >= range_end2  ||     // Otherwise bin totally outside
                  range_bins[b+1] <= range_end1 ) )   // range of ppath
            {
              // Bin limits
              Numeric blim1 = max( range_bins[b], range_end1 );
              Numeric blim2 = min( range_bins[b+1], range_end2 );

              // Determine weight vector to obtain mean inside bin
              Vector hbin(np);
              integration_bin_by_vecmult( hbin, range, blim1, blim2 );
              // The function above handles integration over the bin, while we
              // want the average, so divide weights with bin width
              hbin /= (blim2-blim1);

              for( Index iv=0; iv<nf; iv++ )
                {
                  // Pick out part of iy for frequency
                  //
                  Matrix I = iy( Range(iv*np,np), joker );
                  //
                  if( j_analytical_do )
                    {
                      FOR_ANALYTICAL_JACOBIANS_DO(
                        dI[iq] = diy_dx[iq]( joker, Range(iv*np,np), joker );
                      )
                    }

                  
                  for( Index ip=0; ip<instrument_pol_array[iv].nelem(); ip++ )
                    {
                      // Extract reflectivity for received polarisation
                      Vector refl( np, 0.0 );
                      if( j_analytical_do )
                        { FOR_ANALYTICAL_JACOBIANS_DO( drefl[iq] = 0.0; ) }
                      Vector w = s2p[instrument_pol_array[iv][ip]-1];
                      for( Index i=0; i<w.nelem(); i++ )     // Note that w can
                        {                                    // be shorter than
                          for( Index j=0; j<np; j++ )        // stokes_dim (and
                            {                                // dot product can 
                              refl[j] += I(j,i) * w[i];      // not be used)
                              if( j_analytical_do )
                                {
                                  FOR_ANALYTICAL_JACOBIANS_DO(
                                    for( Index k=0; k<drefl[iq].nrows(); k++ )
                                      { drefl[iq](k,j) += dI[iq](k,j,i) * w[i]; }
                                  )
                                }
                            }
                        }  
                                                             
                                                             
                      // Apply weight vector to get final value
                      //
                      Index iout = nbins * ( p*npolcum[nf] + 
                                             npolcum[iv] + ip ) + b;
                      y[iout] = hbin * refl;
                      //
                      if( j_analytical_do )
                        {
                          FOR_ANALYTICAL_JACOBIANS_DO(
                            for( Index k=0; k<drefl[iq].nrows(); k++ )
                              { jacobian(iout,jacobian_indices[iq][0]+k) =
                                                   hbin * drefl[iq](k,joker); }
                          )
                        }
                      
                      // Same for aux variables
                      for( Index a=0; a<naux; a++ )
                        {
                          if( iy_aux_vars[a] == "Backscattering" )
                            {
                              // I2 matches I, but is transposed
                              Matrix I2 = iy_aux[a](iv,joker,0,joker);
                              refl = 0;
                              for( Index i=0; i<w.nelem(); i++ )
                                { for( Index j=0; j<np; j++ )
                                    { refl[j] += I2(i,j) * w[i]; } }
                              y_aux[a][iout] = hbin * refl;
                            }                          
                          else if( iy_aux[a].nbooks()==1 && iy_aux[a].npages()==1 &&
                              iy_aux[a].nrows() ==1 && iy_aux[a].ncols() ==np )
                            {
                              y_aux[a][iout] = hbin * iy_aux[a](0,0,0,joker);
                            }
                          else
                            {
                              ostringstream os;
                              os << "The auxiliary variable " << iy_aux_vars[a]
                                 << "is not handled by this WSM (but OK when "
                                 << "using *iyCalc*).";
                              throw runtime_error( os.str() );
                            }
                        }
                    }
                }
            }
        }

      // Other aux variables
      //
      Vector geo_pos;
      geo_pos_agendaExecute( ws, geo_pos, ppath, geo_pos_agenda );
      if( geo_pos.nelem() &&  geo_pos.nelem() != atmosphere_dim )
        {
          throw runtime_error( "Wrong size of *geo_pos* obtained from "
                               "*geo_pos_agenda*.\nThe length of *geo_pos* must "
                               "be zero or equal to *atmosphere_dim*." );
        }
      //
      for( Index b=0; b<nbins; b++ )
        {
          for( Index iv=0; iv<nf; iv++ )
            {
              for( Index ip=0; ip<instrument_pol_array[iv].nelem(); ip++ )
                {
                  const Index iout = nbins * ( p*npolcum[nf] + 
                                               npolcum[iv] + ip ) + b;
                  y_f[iout]         = f_grid[iv];
                  y_pol[iout]       = instrument_pol_array[iv][ip];
                  y_pos(iout,joker) = sensor_pos(p,joker);
                  y_los(iout,joker) = sensor_los(p,joker);
                  if( geo_pos.nelem() )
                    { y_geo(iout,joker) = geo_pos; }
                }
            }
        }
    }
}
