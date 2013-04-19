/* Copyright (C) 2012
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
  \date   2010-04-06 

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
#include "montecarlo.h"
#include "rte.h"
#include "special_interp.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PI;





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
   const ArrayOfSingleScatteringData& scat_data_raw,
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
   const Verbosity&                   verbosity )
{
  // Throw error if unsupported features are requested
  if( jacobian_do )
    throw runtime_error( "This method does not provide any jacobians and "
                         "*jacobian_do* must be 0." );
  diy_dx.resize(0);


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
              auxAbsIsp.push_back(ispecies);
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
  Matrix       ppath_vmr, ppath_pnd, ppath_wind, ppath_mag, ppath_f;
  Matrix       ppath_blackrad;
  Tensor5      ppath_abs;
  Tensor4      trans_partial, trans_cumulat, pnd_ext_mat;
  Tensor3      pnd_abs_vec;
  Vector       scalar_tau;
  ArrayOfIndex clear2cloudbox;
  //
  if( np > 1 )
    {
      get_ppath_atmvars(  ppath_p, ppath_t, ppath_vmr, ppath_wind, ppath_mag, 
                          ppath, atmosphere_dim, p_grid, t_field, 
                          vmr_field, wind_u_field, wind_v_field, wind_w_field,
                          mag_u_field, mag_v_field, mag_w_field );      
      get_ppath_f(        ppath_f, ppath, f_grid,  atmosphere_dim, 
                          rte_alonglos_v, ppath_wind );
      get_ppath_abs(      ws, ppath_abs, propmat_clearsky_agenda, ppath, 
                          ppath_p, ppath_t, ppath_vmr, ppath_f, 
                          ppath_mag, f_grid, stokes_dim );
      get_ppath_blackrad( ws, ppath_blackrad, blackbody_radiation_agenda, 
                          ppath, ppath_t, ppath_f );
      if( !cloudbox_on )
        { 
          get_ppath_trans( trans_partial, trans_cumulat, scalar_tau,
                           ppath, ppath_abs, f_grid, stokes_dim );
        }
      else
        {
          Array<ArrayOfSingleScatteringData> scat_data;
          //
          get_ppath_ext(    clear2cloudbox, pnd_abs_vec, pnd_ext_mat, scat_data,
                            ppath_pnd, ppath, ppath_t, stokes_dim, ppath_f, 
                            atmosphere_dim, cloudbox_limits, pnd_field, 
                            use_mean_scat_data, scat_data_raw, verbosity );
          get_ppath_trans2( trans_partial, trans_cumulat, scalar_tau,
                            ppath, ppath_abs, f_grid, stokes_dim, 
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
    { iy_transmission_mult( iy_trans_new, iy_transmission, 
                            trans_cumulat(joker,joker,joker,np-1) ); }

  // Radiative background
  //
  {
    Agenda iy_cb_agenda;  // This OK ???
    get_iy_of_background( ws, iy, diy_dx, 
                          iy_trans_new, jacobian_do, ppath, rte_pos2, 
                          atmosphere_dim, t_field, z_field, vmr_field, 
                          cloudbox_on, stokes_dim, f_grid, iy_main_agenda, 
                          iy_space_agenda, iy_surface_agenda, iy_cb_agenda,
                          verbosity );
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
                                ppath_abs(joker,iv,is1,is2,np-1).sum(); } } } } 
      for( Index j=0; j<auxAbsSpecies.nelem(); j++ )
        { for( Index iv=0; iv<nf; iv++ ) {
            for( Index is1=0; is1<stokes_dim; is1++ ){
              for( Index is2=0; is2<stokes_dim; is2++ ){
                iy_aux[auxAbsSpecies[j]](iv,is1,is2,np-1) = 
                               ppath_abs(auxAbsIsp[j],iv,is1,is2,np-1); } } } }
      // Particle properties
      if( cloudbox_on  )
        {
          // Particle mass content
          for( Index j=0; j<auxPartCont.nelem(); j++ )
            { iy_aux[auxPartCont[j]](0,0,0,np-1) = ppath_pnd(joker,np-1) *
                                      particle_masses(joker,auxPartContI[j]); }
          // Particle field
          for( Index j=0; j<auxPartField.nelem(); j++ )
            { iy_aux[auxPartField[j]](0,0,0,np-1) = 
                                            ppath_pnd(auxPartFieldI[j],np-1); }
        }
      // Radiance for this point is handled above
      //=======================================================================


      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          // Path step average of B: Bbar
          //
          Vector bbar(nf);
          //
          for( Index iv=0; iv<nf; iv++ )  
            { bbar[iv] = 0.5 * ( ppath_blackrad(iv,ip) +
                                 ppath_blackrad(iv,ip+1) ); }

          // Order zero case (ignore scattering, but include particle absorpt.
          if( 0 == 0 ) 
            {
              // Transmission of step
              Tensor3 t = trans_partial(joker,joker,joker,ip);
              
              // Add particle absorption (if any)
              Vector pabs(nf,0);
              bool   any_pabs = false;
              if( clear2cloudbox[ip] >= 0 )
                {
                  any_pabs = true;
                  for( Index iv=0; iv<nf; iv++ )  
                    { pabs[iv] += pnd_abs_vec(iv,0,clear2cloudbox[ip]); }
                }
              if( clear2cloudbox[ip+1] >= 0 )
                {
                  any_pabs = true;
                  for( Index iv=0; iv<nf; iv++ )  
                    { pabs[iv] += pnd_abs_vec(iv,0,clear2cloudbox[ip+1]); }
                }
              if( any_pabs )
                {
                  for( Index iv=0; iv<nf; iv++ )  
                    { t(iv,0,0) *= exp( -0.5 * ppath.lstep[ip] * pabs[iv] ); }
                }
              
              // Perform RT
              emission_rtstep( iy, stokes_dim, bbar, t );
            }
          
          // Order >=1, include scattering
          else
            {
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
                                  ppath_abs(joker,iv,is1,is2,ip).sum(); } } } } 
          for( Index j=0; j<auxAbsSpecies.nelem(); j++ )
            { for( Index iv=0; iv<nf; iv++ ) {
                for( Index is1=0; is1<stokes_dim; is1++ ){
                  for( Index is2=0; is2<stokes_dim; is2++ ){
                    iy_aux[auxAbsSpecies[j]](iv,is1,is2,ip) = 
                                 ppath_abs(auxAbsIsp[j],iv,is1,is2,ip); } } } }
          // Particle properties
          if( cloudbox_on ) 
            {
              // Particle mass content
              for( Index j=0; j<auxPartCont.nelem(); j++ )
                { iy_aux[auxPartCont[j]](0,0,0,ip) = ppath_pnd(joker,ip) *
                                      particle_masses(joker,auxPartContI[j]); }
              // Particle field
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
















































































/* Workspace method: Doxygen documentation will be auto-generated */
/*
void fos_yStandard(Workspace&          ws,
                   Tensor3&            fos_y,
                   Matrix&             iy_aux,
                   ArrayOfTensor3&     diy_dx,
                   const Vector&       rte_pos,
                   const Index&        atmosphere_dim,
                   const Vector&       p_grid,
                   const Vector&       lat_grid,
                   const Tensor3&      z_field,
                   const Tensor3&      t_field,
                   const Tensor4&      vmr_field,
                   const Tensor3&      edensity_field,
                   const Index&        cloudbox_on,
                   const ArrayOfIndex& cloudbox_limits,
                   const Index&        stokes_dim,
                   const Vector&       f_grid,
                   const Agenda&       ppath_agenda,
                   const Agenda&       blackbody_radiation_agenda,
                   const Agenda&       abs_scalar_gas_agenda,
                   const Agenda&       iy_main_agenda,
                   const Tensor3&      iy_transmission,
                   const Tensor4&      pnd_field,
                   const ArrayOfSingleScatteringData&   scat_data_raw,
                   const Agenda&       opt_prop_gas_agenda,
                   const Agenda&       fos_y_agenda,
                   const Matrix&       fos_angles,
                   const Index&        use_mean_scat_data,
                   const Index&        fos_n,
                   const Index&        fos_i,
                   const Verbosity&    verbosity)
{
  // Angles inside these ranges are considered to be equal for 1D and 2D
  const Numeric dza = 0.01;
  const Numeric daa = 1;

  const Index jacobian_do = 0;

  const Index   nfosa = fos_angles.nrows();
        Matrix  tmp;

  fos_y.resize(nfosa,f_grid.nelem(),stokes_dim);

  if( fos_i == fos_n-1 )
    {
      if( atmosphere_dim == 1 )
        { 
          for( Index ia=0; ia<nfosa; ia++ )
            { 
              // To use already calculated data we demand that difference
              // in za is < dza
              Index ihit = -1;
              //
              for( Index it=ia-1; it>=0 && ihit<0; it-- )
                {
                   if( fabs( fos_angles(ia,0) - fos_angles(it,0) ) < dza )
                     { ihit = it; }
                } 
              
              if( ihit >= 0 )
                { 
                  fos_y(ia,joker,joker) = fos_y(ihit,joker,joker); 
                }
              else
                {
                  iy_main_agendaExecute( ws, tmp, 
                                             iy_aux, diy_dx, 0, iy_transmission,
                                             rte_pos, fos_angles(ia,Range(0,1)),
                                             0, jacobian_do, t_field, z_field, 
                                             vmr_field, -1,
                                             iy_main_agenda );
                  fos_y(ia,joker,joker) = tmp;
                }
            }
          
        }
      else if( atmosphere_dim == 2 )
        { 
          Vector rte_los(1);

          for( Index ia=0; ia<nfosa; ia++ )
            { 
              // To use already calculated data we demand that difference
              // in za is < dza, and aa is mirrored with respect to the
              // orbit plane (aa(it) = -aa(ia)).
              Index ihit = -1;
              //
              for( Index it=ia-1; it>=0 && ihit<0; it-- )
                {
                   if( fabs( fos_angles(ia,0) - fos_angles(it,0) ) < dza )
                     { 
                       if( fabs( fos_angles(ia,1) + fos_angles(it,1) ) < daa )
                         { ihit = it; }
                     }
                } 
              
              if( ihit >= 0 )
                { 
                  fos_y(ia,joker,joker) = fos_y(ihit,joker,joker); 
                }
              else
                {
                  // LOS
                  if( fabs(fos_angles(ia,1)) <= 90 )
                    { rte_los[0] = fos_angles(ia,0); }
                  else
                    { rte_los[0] = -fos_angles(ia,0); }

                  // Create stretched latitude grid
                  //
                  Vector lat_stretched( lat_grid.nelem() );
                  //
                  // No strect needed for zenith, nadir and aa= 0 or +-180
                  if( fos_angles(ia,0) > 0  &&  fabs(fos_angles(ia,0)) < 180 &&
                      fos_angles(ia,1) != 0  && fabs(fos_angles(ia,1)) < 180 )
                    {
                      // Stretch factor (a max of 100 is applied)
                      const Numeric stretch = max( 100.0, 
                                     1.0/fabs(cos(DEG2RAD*fos_angles(ia,1))) );
                      const Numeric lat0 = rte_pos[1];
                      for( Index i=0; i<lat_grid.nelem(); i++ )
                        { 
                          lat_stretched[i] = lat0 + stretch*(lat_grid[i]-lat0); 
                        }
                    }
                  else
                    { lat_stretched = lat_grid; }
                  
                  throw runtime_error( "2D FOS can not be used presently "
                       "as the \"lat-stretch\" approach not can be handled!" );

                  iy_main_agendaExecute( ws, tmp, 
                                             iy_aux, diy_dx, 0, iy_transmission,
                                             rte_pos, rte_los, 0, jacobian_do, 
                                             t_field, z_field, vmr_field, -1, 
                                             iy_main_agenda );
                  fos_y(ia,joker,joker) = tmp;
                }
            }
          
        }
      else if( atmosphere_dim == 3 )
        {
          for( Index ia=0; ia<nfosa; ia++ )
            { 
              iy_main_agendaExecute( ws, tmp, 
                                         iy_aux, diy_dx, 0, iy_transmission, 
                                         rte_pos, fos_angles(ia,Range(0,2)),
                                         0, jacobian_do, t_field, z_field, 
                                         vmr_field, -1, iy_main_agenda );
              fos_y(ia,joker,joker) = tmp;
            }
        }
    }
  else
    {
      // The azimuth information is lost for 1D and 2D. Then not possible to
      // handle the latitude stretching for 2D. In fact, the latitude grid
      // should be compressed for many angle combinations!
      if( atmosphere_dim == 2 )
        throw runtime_error( 
             "For atmosphere_dim = 2, only single scattering can be used." );

      Index nlos = 1 + (atmosphere_dim==3);

      for( Index ia=0; ia<nfosa; ia++ )
        { 
          iyFOS( ws, tmp, iy_aux, diy_dx,
                 iy_transmission, rte_pos, fos_angles(ia,Range(0,nlos)), 
                 jacobian_do, atmosphere_dim, p_grid, 
                 z_field, t_field, vmr_field, edensity_field,
                 cloudbox_on, cloudbox_limits, stokes_dim, f_grid, 
                 ppath_agenda, blackbody_radiation_agenda, 
                 abs_scalar_gas_agenda, iy_main_agenda, 
                 pnd_field, scat_data_raw, opt_prop_gas_agenda, fos_y_agenda, 
                 fos_angles, use_mean_scat_data, fos_n, fos_i+1, verbosity);

          fos_y(ia,joker,joker) = tmp;
        }
    }
}
*/


/* Workspace method: Doxygen documentation will be auto-generated */
/*
void iyFOS(Workspace&          ws,
           Matrix&             iy,
           Matrix&             iy_aux,
           ArrayOfTensor3&     diy_dx,
           const Tensor3&      iy_transmission,
           const Vector&       rte_pos,
           const Vector&       rte_los,
           const Index&        jacobian_do,
           const Index&        atmosphere_dim,
           const Vector&       p_grid,
           const Tensor3&      z_field,
           const Tensor3&      t_field,
           const Tensor4&      vmr_field,
           const Tensor3&      edensity_field,
           const Index&        cloudbox_on,
           const ArrayOfIndex& cloudbox_limits,
           const Index&        stokes_dim,
           const Vector&       f_grid,
           const Agenda&       ppath_agenda,
           const Agenda&       blackbody_radiation_agenda,
           const Agenda&       abs_scalar_gas_agenda,
           const Agenda&       iy_main_agenda,
           const Tensor4&      pnd_field,
           const ArrayOfSingleScatteringData&   scat_data_raw,
           const Agenda&       opt_prop_gas_agenda,
           const Agenda&       fos_y_agenda,
           const Matrix&       fos_angles,
           const Index&        use_mean_scat_data,
           const Index&        fos_n,
           const Index&        fos_i,
           const Verbosity&    verbosity)
{
  // Input checks
  if( jacobian_do )
    throw runtime_error( 
     "This method does not yet provide any jacobians (jacobian_do must be 0)" );
  if( !cloudbox_on )
    throw runtime_error( "The cloudbox must be defined to use this method." );
  if( fos_angles.ncols() != 3 )
    throw runtime_error( "The WSV *fos_angles* must have three columns." );
  if( max(fos_angles) <= PI )
    throw runtime_error( 
                  "The WSV *fos_angles* shall be in degrees (not radians)." );
  if( min(fos_angles(joker,0))<0 || max(fos_angles(joker,0))>180 )
    throw runtime_error( 
                 "The zenith angles in *fos_angles* shall be inside [0,180]." );
  if( min(fos_angles(joker,1))<-180 || max(fos_angles(joker,1))>180 )
    throw runtime_error( 
             "The azimuth angles in *fos_angles* shall be inside [-180,180]." );
  if( fos_angles(joker,2).sum() < 6 || fos_angles(joker,2).sum() > 20 )
    throw runtime_error( 
     "The sum of integration weights in *fos_angles* shall be inside [2,20]." );
  if( fos_n < 0 )
    throw runtime_error( "The WSV *fos_n* must be >= 0." );
  if( fos_i < 0 )
    throw runtime_error( "The WSV *fos_i* must be >= 0." );

  // Determine ppath through the cloudbox
  //
  Ppath  ppath;
  //
  ppath_agendaExecute( ws, ppath, rte_pos, rte_los, cloudbox_on, 1, -1,
                       t_field, z_field, vmr_field, edensity_field, -1,
                       ppath_agenda );

  // Check radiative background
  const Index bkgr = ppath_what_background( ppath );
  if( bkgr == 2 )
    throw runtime_error( "Observations where (unscattered) propagation path "
                         "hits the surface inside the cloudbox are not yet "
                         "handled by this method." );
  assert( bkgr == 3 );

  // Get atmospheric and RT quantities for each ppath point/step (inside box)
  // 
  // If np = 1, we only need to determine the radiative background
  //
  // "atmvars"
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w;
  Matrix    ppath_vmr, ppath_pnd;
  // "rtvars"
  Matrix    ppath_emission, ppath_tau;
  Tensor3   wind_field_dummy(0,0,0), iy_trans_new;
  Tensor3   ppath_asp_abs_vec, ppath_pnd_abs_vec, total_transmission;
  Tensor4   ppath_asp_ext_mat, ppath_pnd_ext_mat, ppath_transmission;
  Array<ArrayOfSingleScatteringData>  scat_data;
  //
  const Index np  = ppath.np;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      get_ppath_atmvars( ppath_p, ppath_t, ppath_vmr, 
                         ppath_wind_u, ppath_wind_v, ppath_wind_w,
                         ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                         wind_field_dummy, wind_field_dummy, wind_field_dummy );

      // Particle number densities
      get_ppath_pnd( ppath_pnd, 
                     ppath, atmosphere_dim, cloudbox_limits, pnd_field );

      // Absorption and optical thickness for each step
      get_ppath_cloudrtvars( ws, ppath_asp_abs_vec, ppath_asp_ext_mat,
                            ppath_pnd_abs_vec, ppath_pnd_ext_mat, ppath_transmission, 
                            total_transmission, ppath_emission, scat_data,
                            abs_scalar_gas_agenda, blackbody_radiation_agenda, opt_prop_gas_agenda,
                            ppath, ppath_p, ppath_t, ppath_vmr, ppath_wind_u,
                            ppath_wind_v, ppath_wind_w, ppath_pnd, use_mean_scat_data,
                            scat_data_raw, stokes_dim, f_grid, atmosphere_dim, 1,
                            verbosity);
    }
  else // Just in case, should not happen
    { assert( 0 ); }

  // iy_transmission
  //
  iy_transmission_mult( iy_trans_new, iy_transmission, total_transmission );

  // Get iy for unscattered direction
  //
  // Note that the Ppath positions (ppath.pos) for 1D have one column more
  // than expected by most functions. Only the first atmosphere_dim values
  // shall be copied.
  //
  Vector  rte_pos2;
  {
    Vector  rte_los2;
    rte_pos2 = ppath.pos(ppath.np-1,Range(0,atmosphere_dim));
    rte_los2 = ppath.los(ppath.np-1,joker);
    //
    iy_main_agendaExecute( ws, iy, 
                               iy_aux, diy_dx, 0, iy_trans_new,
                               rte_pos2, rte_los2, 0, jacobian_do, t_field, 
                               z_field, vmr_field, -1, iy_main_agenda );
  }

  // RT for part inside cloudbox
  //
  if( np > 1 )
    {
      // General variables
      const Index   nf    = f_grid.nelem();
      const Index   nfosa = fos_angles.nrows();   // Number of FOS angles
      Matrix  s1(nf,stokes_dim);                  // Scattering source term
      Matrix  s2(nf,stokes_dim,0.0);

      // Help variables for handling of *use_mean_scat_data*
      Index   nfs, ivf;
      //
      if( use_mean_scat_data )
        { nfs = 1;  ivf = 0; }
      else
        { nfs = nf; ivf = 1; }

      // Loop ppath steps
      for( Index ip=np-1; ip>=0; ip-- )
        {              
          // Update scattering source term (new 1 is old 2)
          s1 = s2;

          // Scattering source term (is zero if no particles)
          if( max(ppath_pnd(joker,ip)) < 1e-3 )
            {
              s2 = 0.0;
            }
          else
            {
              // Determine incoming radiation (here Y, WSV equals fos_y)
              Tensor3  Y;
              fos_y_agendaExecute( ws, Y, rte_pos2, fos_angles, fos_n, fos_i, 
                                                                fos_y_agenda );


              // Direction of outgoing scattered radiation (which is reversed
              // to LOS). Note that this rte_los2 is only used for extracting
              // scattering properties.
              Vector rte_los2;
              mirror_los( rte_los2, ppath.los(ip,joker), atmosphere_dim );

              // Determine phase matrix for fov_angles
              Tensor4  P( nfosa, nfs, stokes_dim, stokes_dim );
              Matrix   P1( stokes_dim, stokes_dim );
              //
              for( Index ia=0; ia<nfosa; ia++ )
                { 
                  for( Index iv=0; iv<nfs; iv++ )
                    {
                      pha_mat_singleCalc( P1, rte_los2[0], rte_los2[1],
                                          fos_angles(ia,0), fos_angles(ia,1),
                                          scat_data[iv], stokes_dim, 
                                          ppath_pnd(joker,ip), ppath_t[ip],
                                          verbosity );
                      P(ia,iv,joker,joker) = P1;
                    }
                }

              // Scattering source term
              s2 = 0.0;
              for( Index iv=0; iv<nf; iv++ )
                { 
                  Vector sp(stokes_dim);
                  for( Index ia=0; ia<nfosa; ia++ )
                    { 
                      mult( sp, P(ia,iv*ivf,joker,joker), Y(ia,iv,joker) );
                      sp           *= fos_angles(ia,2);
                      s2(iv,joker) += sp;
                    }
                }
            }

          // RT of ppath step (nothing to do when at upper point)
          if( ip < np-1 )
            {
              // Loop frequencies
              for( Index iv=0; iv<nf; iv++ )
                {
                  // Calculate average of absorption, extinction etc.
                  Matrix  ext_mat( stokes_dim, stokes_dim,0 );
                  Vector  abs_vec(stokes_dim,0.0);
                  Vector  s(stokes_dim,0.0);
                  Numeric b = 0.5 * ( ppath_emission(iv,ip) + 
                                      ppath_emission(iv,ip+1) );
                  for( Index is1=0; is1<stokes_dim; is1++ )
                    { 
                      s[is1]  = 0.5 * ( s1(iv,is1) + s2(iv,is1) );
                      abs_vec[is1] = 0.5 * ( 
                                          ppath_asp_abs_vec(iv,is1,ip+1) +
                                          ppath_asp_abs_vec(iv,is1,ip)   +
                                          ppath_pnd_abs_vec(iv,is1,ip+1) +
                                          ppath_pnd_abs_vec(iv,is1,ip) );
                      for( Index is2=0; is2<stokes_dim; is2++ )
                        {
                          ext_mat(is1,is2) = 0.5 * (
                                          ppath_asp_ext_mat(iv,is1,is2,ip+1) +
                                          ppath_asp_ext_mat(iv,is1,is2,ip)   +
                                          ppath_pnd_ext_mat(iv,is1,is2,ip+1) +
                                          ppath_pnd_ext_mat(iv,is1,is2,ip) );
                        }
                    }

                  // RT for step
                  Matrix trans_mat(stokes_dim,stokes_dim); 
                  //
                  rte_step_doit( iy(iv,joker), trans_mat, ext_mat, abs_vec,
                                                       s, ppath.lstep[ip], b );
                }
            }
        } 
    }
}
*/
