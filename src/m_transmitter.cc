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

/**
  @file   m_transmitter.cc
  @author Patrick Eriksson <patrick.eriksson@chalmers.se>
  @date   2012-10-31

  @brief  Workspace functions related to transmitters and radiative transfer
  for transmitted signals.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "auto_md.h"
#include "matpack_complex.h"
#include "geodetic.h"
#include "jacobian.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "rte.h"
#include "sensor.h"

extern const Numeric DEG2RAD;
extern const Numeric PI;
extern const Numeric RAD2DEG;
extern const Numeric SPEED_OF_LIGHT;

/* Workspace method: Doxygen documentation will be auto-generated */
/* Has not been updated since v2.2
void iyRadioLink(
         Workspace&                   ws,
         Matrix&                      iy,
         ArrayOfTensor4&              iy_aux,
         Ppath&                       ppath,
         ArrayOfTensor3&              diy_dx,
   const Index&                       stokes_dim,
   const Vector&                      f_grid,
   const Index&                       atmosphere_dim,
   const Vector&                      p_grid,
   const Vector&                      lat_grid,
   const Vector&                      lon_grid,
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
   const Vector&                      refellipsoid,
   const Matrix&                      z_surface,
   const Index&                       cloudbox_on,
   const ArrayOfIndex&                cloudbox_limits,
   const Tensor4&                     pnd_field,
   const ArrayOfArrayOfSingleScatteringData& scat_data,
   const Matrix&                      particle_masses,
   const ArrayOfString&               iy_aux_vars,
   const Index&                       jacobian_do,
   const Agenda&                      ppath_agenda,
   const Agenda&                      ppath_step_agenda,
   const Agenda&                      propmat_clearsky_agenda,
   const Agenda&                      iy_transmitter_agenda,
   const Index&                       iy_agenda_call1,
   const Tensor3&                     iy_transmittance,
   const Vector&                      rte_pos,      
   const Vector&                      rte_los,      
   const Vector&                      rte_pos2,      
   const Numeric&                     rte_alonglos_v,      
   const Numeric&                     ppath_lmax,
   const Numeric&                     ppath_lraytrace,
   const Index&                       defocus_method,
   const Numeric&                     defocus_shift,
   const Verbosity&                   verbosity )
{
  // Throw error if unsupported features are requested
  if( !iy_agenda_call1 )
    throw runtime_error( 
                  "Recursive usage not possible (iy_agenda_call1 must be 1)" );
  if( !iy_transmittance.empty() )
    throw runtime_error( "*iy_transmittance* must be empty" );
  if( jacobian_do )
    throw runtime_error( "This method does not provide any jacobians and "
                         "*jacobian_do* must be 0." );
  if( defocus_method < 1 || defocus_method > 2 )
    throw runtime_error( "Allowed choices for *defocus_method* is 1 and 2." );
  diy_dx.resize(0);


  //- Determine propagation path
  ppath_agendaExecute( ws, ppath, ppath_lmax, ppath_lraytrace, rte_pos, rte_los, 
                       rte_pos2, cloudbox_on, 0, t_field, z_field, vmr_field, 
                       f_grid, ppath_agenda );

  //- Set np to zero if ground intersection
  const Index radback = ppath_what_background(ppath);
  // np should already be 1 fon non-OK cases, but for extra safety ...
  if( radback == 0  || radback == 2 )
    { ppath.np = 1; }

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
  Index auxPressure        = -1,
        auxTemperature     = -1,
        auxAbsSum          = -1,
        auxPartExt         = -1,
        auxImpactParam     = -1,
        auxFreeSpaceLoss   = -1,
        auxFreeSpaceAtte   = -1,
        auxAtmosphericLoss = -1,
        auxDefocusingLoss  = -1,
        auxFarRotTotal     = -1,
        auxFarRotSpeed     = -1,
        auxExtraPathDelay  = -1,
        auxBendingAngle    = -1;
  ArrayOfIndex auxAbsSpecies(0), auxAbsIsp(0);
  ArrayOfIndex auxVmrSpecies(0), auxVmrIsp(0);
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
      else if( iy_aux_vars[i] == "Impact parameter" )
        { auxImpactParam = i;       iy_aux[i].resize( 1, 1, 1, 1 ); }
      else if( iy_aux_vars[i] == "Free space loss" )
        { auxFreeSpaceLoss = i;     iy_aux[i].resize( 1, 1, 1, 1 ); }
      else if( iy_aux_vars[i] == "Free space attenuation" )
        { auxFreeSpaceAtte = i;     iy_aux[i].resize( 1, 1, 1, np ); }
      else if( iy_aux_vars[i] == "Atmospheric loss" )
        { auxAtmosphericLoss = i;   iy_aux[i].resize( nf, 1, 1, 1 ); } 
      else if( iy_aux_vars[i] == "Defocusing loss" )
        { auxDefocusingLoss = i;    iy_aux[i].resize( 1, 1, 1, 1 ); }
      else if( iy_aux_vars[i] == "Faraday rotation" )
        { auxFarRotTotal = i; iy_aux[i].resize( nf, 1, 1, 1 ); iy_aux[i] = 0; }
      else if( iy_aux_vars[i] == "Faraday speed" )
        { auxFarRotSpeed = i; iy_aux[i].resize( nf, 1, 1, np ); iy_aux[i] = 0; }
      else if( iy_aux_vars[i] == "Extra path delay" )
        { auxExtraPathDelay = i; iy_aux[i].resize( 1, 1, 1, 1 ); }
      else if( iy_aux_vars[i] == "Bending angle" )
        { auxBendingAngle = i;      iy_aux[i].resize( 1, 1, 1, 1 ); } 
      else
        {
          ostringstream os;
          os << "In *iy_aux_vars* you have included: \"" << iy_aux_vars[i]
             << "\"\nThis choice is not recognised.";
          throw runtime_error( os.str() );
        }
    }
  
  // Special stuff to handle Faraday rotation
  Index ife = -1;  // When ready, ife refers abs_per_species
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
  //===========================================================================


  // Handle cases whn no link was establsihed
  if( radback == 0  || radback == 2 )
    {
      Numeric fillvalue = 0;
      if( radback == 0 )
        { fillvalue = NAN; }
      //
      iy.resize( nf, stokes_dim );
      iy = fillvalue;
      //
      for( Index i=0; i<naux; i++ )
        { iy_aux[i] = fillvalue; }
      //
      return;
    }


  // Transmitted signal
  //
  iy_transmitter_agendaExecute( ws, iy, f_grid, 
                                ppath.pos(np-1,Range(0,atmosphere_dim)),
                                ppath.los(np-1,joker), iy_transmitter_agenda );
  if( iy.ncols() != stokes_dim  ||  iy.nrows() != nf )
    { throw runtime_error( "The size of *iy* returned from "
                                 "*iy_transmitter_agenda* is not correct." ); }


  // Get atmospheric and attenuation quantities for each ppath point/step
  //
  Vector       ppath_p, ppath_t;
  Matrix       ppath_vmr, ppath_pnd, ppath_mag, ppath_wind, ppath_f, ppath_t_nlte;
  ArrayOfArrayOfPropagationMatrix      abs_per_species, dummy_dppath_ext_dx;
  ArrayOfPropagationMatrix      ppath_ext, pnd_ext_mat;
  Tensor4 trans_partial, trans_cumulat;
  ArrayOfArrayOfStokesVector dummy_dppath_nlte_dx;
  ArrayOfStokesVector      dummy_ppath_nlte_source;
  Vector       scalar_tau;
  ArrayOfIndex clear2cloudy, dummy_lte;
  ArrayOfMatrix        dummy_ppath_dpnd_dx;
  ArrayOfTensor4       dummy_dpnd_field_dx;
  const Tensor4 t_nlte_field_empty(0,0,0,0);
  //
  if( np > 1 )
    {
      get_ppath_atmvars( ppath_p, ppath_t, ppath_t_nlte, ppath_vmr,
                         ppath_wind, ppath_mag, 
                         ppath, atmosphere_dim, p_grid, t_field, t_nlte_field_empty,
                         vmr_field, wind_u_field, wind_v_field, wind_w_field,
                         mag_u_field, mag_v_field, mag_w_field );
      
      get_ppath_f( ppath_f, ppath, f_grid,  atmosphere_dim, 
                   rte_alonglos_v, ppath_wind );
      
      get_ppath_pmat( ws, ppath_ext, dummy_ppath_nlte_source, dummy_lte, abs_per_species, 
                      dummy_dppath_ext_dx, dummy_dppath_nlte_dx,
                      propmat_clearsky_agenda, ArrayOfRetrievalQuantity(0), ppath, 
                      ppath_p, ppath_t, ppath_t_nlte, ppath_vmr, ppath_f, 
                      ppath_mag, f_grid, stokes_dim, iaps );
      
      if( !cloudbox_on )
        { 
          ArrayOfArrayOfIndex  extmat_case;          
          get_ppath_trans( trans_partial, extmat_case, trans_cumulat,
                           scalar_tau, ppath, ppath_ext, f_grid, stokes_dim );
        }
      else
        {
          // Extract basic scattering data
          Array<ArrayOfArrayOfSingleScatteringData> scat_data_single;
          ArrayOfArrayOfIndex                       extmat_case;          
          Tensor3                                   pnd_abs_vec;

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

  // Ppath length variables
  //
  Numeric lbg;  // Bent geometrical length of ray path
  Numeric lba;  // Bent apparent length of ray path
  //
  lbg = ppath.end_lstep;
  lba = lbg;

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
        { iy_aux[auxVmrSpecies[j]](0,0,0,np-1) = ppath_vmr(auxVmrIsp[j],np-1);}
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
      // Free space
      if( auxFreeSpaceAtte >= 0 )
        { iy_aux[auxFreeSpaceAtte](joker,0,0,np-1) = 2/lbg; }
      // Faraday speed
      if( auxFarRotSpeed >= 0 )
        { for( Index iv=0; iv<nf; iv++ ) {
            iy_aux[auxFarRotSpeed](iv,0,0,np-1) = 0.5 *
                                          abs_per_species[np-1][ife](iv,1,2); } }
      //=======================================================================

      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          // Lengths
          lbg += ppath.lstep[ip];
          lba += ppath.lstep[ip] * (ppath.ngroup[ip]+ppath.ngroup[ip+1]) / 2.0;

          // Atmospheric loss of path step + Faraday rotation
          if( stokes_dim == 1 )
            {
              for( Index iv=0; iv<nf; iv++ )  
                { iy(iv,0) = iy(iv,0) * trans_partial(iv,0,0,ip); }
            }
          else
            {
              for( Index iv=0; iv<nf; iv++ )  
                {
                  // Unpolarised:
                  if( is_diagonal( trans_partial(iv,joker,joker,ip) ) )
                    {
                      for( Index is=0; is<ns; is++ )
                        { iy(iv,is) = iy(iv,is) * trans_partial(iv,is,is,ip); }
                    }
                  // The general case:
                  else
                    {
                      Vector t1(ns);
                      mult( t1, trans_partial(iv,joker,joker,ip), iy(iv,joker));
                      iy(iv,joker) = t1;
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
          // Free space loss
          if( auxFreeSpaceAtte >= 0 )
            { iy_aux[auxFreeSpaceAtte](joker,0,0,ip) = 2/lbg; }
          // Faraday rotation, total
          if( auxFarRotTotal >= 0 )
            { for( Index iv=0; iv<nf; iv++ ) {
                iy_aux[auxFarRotTotal](iv,0,0,0) += RAD2DEG * ppath.lstep[ip] *
                  0.25 * ( abs_per_species[ip][ife](iv,1,2)+
                           abs_per_species[ip+1][ife](iv,1,2)); } }
          // Faraday speed
          if( auxFarRotSpeed >= 0 )
            { for( Index iv=0; iv<nf; iv++ ) {  
                iy_aux[auxFarRotSpeed](iv,0,0,ip) = 0.5 *
                                            abs_per_species[ip][ife](iv,1,2); } }
          //===================================================================
        }


      //=== iy_aux part =======================================================
      if( auxAtmosphericLoss >= 0 )
        { iy_aux[auxAtmosphericLoss](joker,0,0,0) = iy(joker,0); }      
      if( auxImpactParam >= 0 )
        { 
          ARTS_ASSERT( ppath.constant >= 0 );
          iy_aux[auxImpactParam](joker,0,0,0) = ppath.constant; 
        }
      //=======================================================================


      // Remaing length of ppath
      lbg += ppath.start_lstep;
      lba += ppath.start_lstep;


      // Determine total free space loss
      Numeric fspl = 1 / ( 4 * PI * lbg*lbg ); 
      //
      if( auxFreeSpaceLoss >= 0 )
        { iy_aux[auxFreeSpaceLoss] = fspl; }


      // Determine defocusing loss
      Numeric dfl = 1;
      if( defocus_method == 1 )
        {
          defocusing_general( ws, dfl, ppath_step_agenda, atmosphere_dim, 
                              p_grid, lat_grid, lon_grid, t_field, z_field, 
                              vmr_field, f_grid, refellipsoid, 
                              z_surface, ppath, ppath_lmax, ppath_lraytrace,
                              defocus_shift, verbosity );
        }
      else if( defocus_method == 2 )
        { defocusing_sat2sat( ws, dfl, ppath_step_agenda, atmosphere_dim, 
                              p_grid, lat_grid, lon_grid, t_field, z_field, 
                              vmr_field, f_grid, refellipsoid, 
                              z_surface, ppath, ppath_lmax, ppath_lraytrace, 
                              defocus_shift, verbosity ); 
        }
      if( auxDefocusingLoss >= 0 )
        { iy_aux[auxDefocusingLoss] = dfl; }



      // Include free space and defocusing losses
      iy *= fspl*dfl;


      // Extra path delay
      if( auxExtraPathDelay >= 0 )
        {
          // Radius of rte_pos and rte_pos2
          const Numeric r1 = ppath.end_pos[0] +
                             pos2refell_r( atmosphere_dim, refellipsoid, 
                                           lat_grid, lon_grid, ppath.end_pos );
          const Numeric r2 = ppath.start_pos[0] +
                             pos2refell_r( atmosphere_dim, refellipsoid, 
                                         lat_grid, lon_grid, ppath.start_pos );

          // Geometrical distance between start and end point
          Numeric lgd ;
          if( atmosphere_dim <= 2 )
            { distance2D( lgd, r1, ppath.end_pos[1], r2, ppath.start_pos[1] ); }
          else 
            { distance3D( lgd, r1, ppath.end_pos[1],   ppath.end_pos[2],
                               r2, ppath.start_pos[1], ppath.start_pos[2] ); }
          //
          iy_aux[auxExtraPathDelay] = ( lba - lgd ) / SPEED_OF_LIGHT;
        }


      // Bending angle
      if( auxBendingAngle >= 0 )
        { 
          Numeric ba = -999;
          bending_angle1d( ba, ppath ); 
          //
          iy_aux[auxBendingAngle] = ba;
        }
    }
}
*/

/* Workspace method: Doxygen documentation will be auto-generated */
void iyTransmissionStandard(Workspace& ws,
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
                            const ArrayOfString& iy_aux_vars,
                            const Index& jacobian_do,
                            const ArrayOfRetrievalQuantity& jacobian_quantities,
                            const Ppath& ppath,
                            const Agenda& propmat_clearsky_agenda,
                            const Agenda& water_p_eq_agenda,
                            const Agenda& iy_transmitter_agenda,
                            const Index& iy_agenda_call1,
                            const Tensor3& iy_transmittance,
                            const Numeric& rte_alonglos_v,
                            const Verbosity&) {
  //  Init Jacobian quantities?
  const Index j_analytical_do = jacobian_do ? do_analytical_jacobian<1>(jacobian_quantities) : 0;
    
  // Some basic sizes
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index nq = j_analytical_do ? jacobian_quantities.nelem() : 0;

  // Radiative background index
  const Index rbi = ppath_what_background(ppath);

  // Checks of input
  // Throw error if unsupported features are requested
  if (!iy_agenda_call1)
    throw runtime_error(
        "Recursive usage not possible (iy_agenda_call1 must be 1)");
  if (!iy_transmittance.empty())
    throw runtime_error("*iy_transmittance* must be empty");
  if (rbi < 1 || rbi > 9)
    throw runtime_error(
        "ppath.background is invalid. Check your "
        "calculation of *ppath*?");
  if (jacobian_do) {
    if (dpnd_field_dx.nelem() != jacobian_quantities.nelem())
      throw runtime_error(
          "*dpnd_field_dx* not properly initialized:\n"
          "Number of elements in dpnd_field_dx must be equal number of jacobian"
          " quantities.\n(Note: jacobians have to be defined BEFORE *pnd_field*"
          " is calculated/set.");
  }
  // iy_aux_vars checked below

  // Transmitted signal
  iy_transmitter_agendaExecute(ws,
                               iy,
                               f_grid,
                               ppath.pos(np - 1, Range(0, atmosphere_dim)),
                               ppath.los(np - 1, joker),
                               iy_transmitter_agenda);
  if (iy.ncols() != ns || iy.nrows() != nf) {
    ostringstream os;
    os << "The size of *iy* returned from *iy_transmitter_agenda* is\n"
       << "not correct:\n"
       << "  expected size = [" << nf << "," << stokes_dim << "]\n"
       << "  size of iy    = [" << iy.nrows() << "," << iy.ncols() << "]\n";
    throw runtime_error(os.str());
  }
  
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
  Index auxOptDepth = -1;
  //
  for (Index i = 0; i < naux; i++) {
    iy_aux[i].resize(nf, ns);
    iy_aux[i] = 0;

    if (iy_aux_vars[i] == "Radiative background")
      iy_aux[i](joker, 0) = (Numeric)min((Index)2, rbi - 1);
    else if (iy_aux_vars[i] == "Optical depth")
      auxOptDepth = i;
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
  //
  ppvar_trans_cumulat.resize(np, nf, ns, ns);
  ppvar_trans_partial.resize(np, nf, ns, ns);

  ArrayOfRadiationVector lvl_rad(np, RadiationVector(nf, ns));
  ArrayOfArrayOfRadiationVector dlvl_rad(
      np, ArrayOfRadiationVector(nq, RadiationVector(nf, ns)));

  ArrayOfTransmissionMatrix lyr_tra(np, TransmissionMatrix(nf, ns));
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
    ppvar_pnd.resize(0, 0);
    ppvar_f.resize(0, 0);
    ppvar_iy.resize(0, 0, 0);
    ppvar_trans_cumulat = 0;
    ppvar_trans_partial = 0;
    for (Index iv = 0; iv < nf; iv++) {
      for (Index is = 0; is < ns; is++) {
        ppvar_trans_cumulat(0,iv,is,is) = 1;
        ppvar_trans_partial(0,iv,is,is) = 1;        
      }
    }
    
  } else {
    // ppvar_iy
    ppvar_iy.resize(nf, ns, np);
    ppvar_iy(joker, joker, np - 1) = iy;

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

    // pnd_field
    ArrayOfMatrix ppvar_dpnd_dx(0);
    //
    if (cloudbox_on)
      get_ppath_cloudvars(clear2cloudy,
                          ppvar_pnd,
                          ppvar_dpnd_dx,
                          ppath,
                          atmosphere_dim,
                          cloudbox_limits,
                          pnd_field,
                          dpnd_field_dx);
    else {
      clear2cloudy.resize(np);
      for (Index ip = 0; ip < np; ip++) clear2cloudy[ip] = -1;
    }

    // Size radiative variables always used
    PropagationMatrix K_this(nf, ns), K_past(nf, ns), Kp(nf, ns);
    StokesVector a(nf, ns), S(nf, ns), Sp(nf, ns);
    ArrayOfIndex lte(np);

    // Init variables only used if analytical jacobians done
    Vector dB_dT(0);
    ArrayOfPropagationMatrix dK_this_dx(nq), dK_past_dx(nq), dKp_dx(nq);
    ArrayOfStokesVector da_dx(nq), dS_dx(nq), dSp_dx(nq);
    //
    Index temperature_derivative_position = -1;
    bool do_hse = false;

    if (j_analytical_do) {
      dB_dT.resize(nf);
      FOR_ANALYTICAL_JACOBIANS_DO(dK_this_dx[iq] = PropagationMatrix(nf, ns);
                                  dK_past_dx[iq] = PropagationMatrix(nf, ns);
                                  dKp_dx[iq] = PropagationMatrix(nf, ns);
                                  da_dx[iq] = StokesVector(nf, ns);
                                  dS_dx[iq] = StokesVector(nf, ns);
                                  dSp_dx[iq] = StokesVector(nf, ns);
                                  if (jacobian_quantities[iq] == Jacobian::Atm::Temperature) {
                                    temperature_derivative_position = iq;
                                    do_hse = jacobian_quantities[iq].Subtag() ==
                                             "HSE on";
                                  })
    }

    // Loop ppath points and determine radiative properties
    for (Index ip = 0; ip < np; ip++) {
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

      if (j_analytical_do) {
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
      }

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
        K_this += Kp;

        if (j_analytical_do) {
          FOR_ANALYTICAL_JACOBIANS_DO(dK_this_dx[iq] += dKp_dx[iq];)
        }
      }

      // Transmission
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

      swap(K_past, K_this);
      swap(dK_past_dx, dK_this_dx);
    }
  }

  const ArrayOfTransmissionMatrix tot_tra =
      cumulative_transmission(lyr_tra, CumulativeTransmission::Forward);

  // iy_aux: Optical depth
  if (auxOptDepth >= 0) {
    for (Index iv = 0; iv < nf; iv++)
      iy_aux[auxOptDepth](iv, 0) = -std::log(tot_tra[np - 1](iv, 0, 0));
  }

  lvl_rad[np - 1] = iy;

  // Radiative transfer calculations
  for (Index ip = np - 2; ip >= 0; ip--) {
    lvl_rad[ip] = lvl_rad[ip + 1];
    update_radiation_vector(lvl_rad[ip],
                            dlvl_rad[ip],
                            dlvl_rad[ip + 1],
                            RadiationVector(),
                            RadiationVector(),
                            ArrayOfRadiationVector(),
                            ArrayOfRadiationVector(),
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
                            RadiativeTransferSolver::Transmission);
  }

  // Copy back to ARTS external style
  iy = lvl_rad[0];
  //
  if (np > 1) {
    for (Index ip = 0; ip < lvl_rad.nelem(); ip++) {
      ppvar_trans_cumulat(ip, joker, joker, joker) = tot_tra[ip];
      ppvar_trans_partial(ip, joker, joker, joker) = lyr_tra[ip];
      ppvar_iy(joker, joker, ip) = lvl_rad[ip];
      if (j_analytical_do)
        FOR_ANALYTICAL_JACOBIANS_DO(diy_dpath[iq](ip, joker, joker) =
                                        dlvl_rad[ip][iq];);
    }
  }

  // Finalize analytical Jacobians
  if (j_analytical_do) {
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
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iy_transmitterMultiplePol(Matrix& iy,
                               const Index& stokes_dim,
                               const Vector& f_grid,
                               const ArrayOfIndex& instrument_pol,
                               const Verbosity&) {
  const Index nf = f_grid.nelem();

  if (instrument_pol.nelem() != nf)
    throw runtime_error(
        "The length of *f_grid* and the number of elements "
        "in *instrument_pol* must be equal.");

  iy.resize(nf, stokes_dim);

  for (Index i = 0; i < nf; i++) {
    stokes2pol(iy(i, joker), stokes_dim, instrument_pol[i], 1);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iy_transmitterSinglePol(Matrix& iy,
                             const Index& stokes_dim,
                             const Vector& f_grid,
                             const ArrayOfIndex& instrument_pol,
                             const Verbosity&) {
  const Index nf = f_grid.nelem();

  if (instrument_pol.nelem() != 1)
    throw runtime_error(
        "The number of elements in *instrument_pol* must be 1.");

  iy.resize(nf, stokes_dim);

  stokes2pol(iy(0, joker), stokes_dim, instrument_pol[0], 1);

  for (Index i = 1; i < nf; i++) {
    iy(i, joker) = iy(0, joker);
  }
}
