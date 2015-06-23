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



/* Workspace method: Doxygen documentation will be auto-generated */
void iyCloudRadar(
         Workspace&                   ws,
         Matrix&                      iy,
         ArrayOfTensor4&              iy_aux,
         Ppath&                       ppath,
   const Index&                       stokes_dim,
   const Vector&                      f_grid,
   const Index&                       atmosphere_dim,
   const Vector&                      p_grid,
   const Tensor3&                     z_field,
   const Tensor3&                     t_field,
   const Tensor4&                     vmr_field,
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
   const Agenda&                      propmat_clearsky_agenda,
   const Agenda&                      iy_transmitter_agenda,
   const Index&                       iy_agenda_call1,
   const Tensor3&                     iy_transmission,
   const Vector&                      rte_pos,      
   const Vector&                      rte_los,      
   const Numeric&                     rte_alonglos_v,      
   const Numeric&                     ppath_lraytrace,
   const Numeric&                     ze_tref,
   const Verbosity&                   verbosity )
{
  // Throw error if unsupported features are requested
  if( !iy_agenda_call1 )
    throw runtime_error( 
                  "Recursive usage not possible (iy_agenda_call1 must be 1)" );
  if( iy_transmission.ncols() )
    throw runtime_error( "*iy_transmission* must be empty" );
  if( !cloudbox_on )
    throw runtime_error( 
                    "The cloudbox must be activated (cloudbox_on must be 1)" );
  if( jacobian_do )
    throw runtime_error( 
        "This method does not provide any jacobians (jacobian_do must be 0)" );

  // Determine propagation path
  //
  ppath_agendaExecute( ws, ppath, ppath_lraytrace, rte_pos, rte_los, Vector(0),
                       0, 0, t_field, z_field, vmr_field, 
                       f_grid, ppath_agenda );

  // Some basic sizes
  //
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;


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
  Vector    ppath_p, ppath_t;
  Matrix    ppath_vmr, ppath_pnd, ppath_wind, ppath_mag, ppath_f, ppath_t_nlte;
  Tensor3   dummy_ppath_nlte_source;
  Tensor4   ppath_ext, trans_partial, trans_cumulat, pnd_ext_mat;
  Tensor5   dummy_abs_per_species;
  Vector    scalar_tau;
  ArrayOfIndex clear2cloudbox, dummy_lte;
  Array<ArrayOfArrayOfSingleScatteringData> scat_data_single;
  const Tensor4 t_nlte_field_dummy;
  //
  if( np > 1 )
    {
      get_ppath_atmvars( ppath_p, ppath_t, ppath_t_nlte, ppath_vmr,
                         ppath_wind, ppath_mag, 
                         ppath, atmosphere_dim, p_grid, t_field, t_nlte_field_dummy, vmr_field,
                         wind_u_field, wind_v_field, wind_w_field,
                         mag_u_field, mag_v_field, mag_w_field );
      get_ppath_f(       ppath_f, ppath, f_grid,  atmosphere_dim, 
                         rte_alonglos_v, ppath_wind );
      get_ppath_pmat(    ws, ppath_ext, dummy_ppath_nlte_source, dummy_lte, 
                         dummy_abs_per_species, 
                         propmat_clearsky_agenda, ppath, 
                         ppath_p, ppath_t, ppath_t_nlte, ppath_vmr, ppath_f, ppath_mag,
                         f_grid, stokes_dim, ArrayOfIndex(0) );
      if( !cloudbox_on )
        { 
          ArrayOfArrayOfIndex  extmat_case;
          get_ppath_trans( trans_partial, extmat_case, trans_cumulat, 
                           scalar_tau, ppath, ppath_ext, f_grid, stokes_dim );
        }
      else
        {
          // Extract basic scattering data
          ArrayOfArrayOfIndex  extmat_case;
          Tensor3              pnd_abs_vec;
          //
          get_ppath_ext(    clear2cloudbox, pnd_abs_vec, pnd_ext_mat, scat_data_single,
                            ppath_pnd, ppath, ppath_t, stokes_dim, ppath_f, 
                            atmosphere_dim, cloudbox_limits, pnd_field, 
                            use_mean_scat_data, scat_data, verbosity );
          get_ppath_trans2( trans_partial, extmat_case, trans_cumulat, scalar_tau, 
                            ppath, ppath_ext, f_grid, stokes_dim, 
                            clear2cloudbox, pnd_ext_mat );
        }
    }


  // Size iy and set to zero (as not filled if pnd=0)
  iy.resize( nf*np, ns );
  iy = 0;

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
      
      if( clear2cloudbox[ip] >= 0 )
        {
          for( Index iv=0; iv<nf; iv++ )
            {
              // Obtain scattering matrix
              Matrix P(ns,ns);
              pha_mat_singleCalc( P, los_sca[0], los_sca[1], los_inc[0], 
                                  los_inc[1], scat_data_single[iv], stokes_dim, 
                                  ppath_pnd(joker,ip), ppath_t[ip], verbosity );

              // Combine iy0, double transmission and scattering matrix
              Vector iy1(stokes_dim), iy2(stokes_dim);
              mult( iy1, tr_rev(iv,joker,joker), iy0(iv,joker) );
              mult( iy2, P,                      iy1 );
              mult( iy(iv*np+ip,joker), trans_cumulat(iv,joker,joker,ip), iy2 );

              //=== iy_aux part ===========================================
              // Backscattering
              if( auxBackScat >= 0 ) {
                mult( iy_aux[auxBackScat](iv,joker,0,ip), P, iy0(iv,joker) ); }
              // Transmission
              if( auxTrans >= 0 ) 
                { for( Index is1=0; is1<ns; is1++ ){
                    for( Index is2=0; is2<ns; is2++ ){
                      iy_aux[auxTrans](iv,is1,is2,ip) = 
                                                      tr_rev(iv,is1,is2); } } }

              // Update tr_rev
              if( ip<np-1 )
                {
                  Matrix tmp = tr_rev(iv,joker,joker);
                  mult( tr_rev(iv,joker,joker), 
                        trans_partial(iv,joker,joker,ip), tmp );
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
      //===================================================================
    } 

  if( iy_unit == "1" )
    {}
  else if( iy_unit == "Ze" )
    {
      // Get refractive index for water
      Matrix complex_n;
      complex_n_water_liebe93( complex_n, f_grid, ze_tref );

      // Common conversion factor
      const Numeric a = 4e18/(PI*PI*PI*PI);

      for( Index iv=0; iv<nf; iv++ )
        {
          // Calculate the dielectric factor
          Complex n( complex_n(iv,0), complex_n(iv,1) );
          Complex n2 = n*n;
          Complex K  = ( n2 - Numeric(1.0) ) / ( n2 + Numeric(2.0) );
          Numeric absK = abs( K );
          Numeric la = SPEED_OF_LIGHT / f_grid[iv];
          Numeric fac = a*la*la*la*la / ( absK * absK ); // Alalalala!
      
          iy(Range(iv*np,np),joker) *= fac;

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
void yCloudRadar(
         Workspace&             ws,
         Vector&                y,
         ArrayOfVector&         y_aux,
   const Index&                 atmfields_checked,
   const Index&                 atmgeom_checked,
   const String&                iy_unit,   
   const ArrayOfString&         iy_aux_vars,
   const Index&                 stokes_dim,
   const Vector&                f_grid,
   const Tensor3&               t_field,
   const Tensor3&               z_field,
   const Tensor4&               vmr_field,
   const Index&                 cloudbox_on,
   const Index&                 cloudbox_checked,
   const Matrix&                sensor_pos,
   const Matrix&                sensor_los,
   const Index&                 sensor_checked,
   const Agenda&                iy_main_agenda,
   const ArrayOfArrayOfIndex&   sensor_pol_array,
   const Vector&                range_bins,
   const Verbosity& )
{
  // Important sizes
  const Index npos    = sensor_pos.nrows();
  const Index nbins   = range_bins.nelem() - 1;
  const Index nf      = f_grid.nelem();
  const Index naux    = iy_aux_vars.nelem();

  //---------------------------------------------------------------------------
  // Input checks
  //---------------------------------------------------------------------------

  // Basics
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
  if( sensor_pol_array.nelem() != nf )
    throw runtime_error( "The main length of *sensor_pol_array* must match "
                         "the number of frequencies." );


  //---------------------------------------------------------------------------
  // The calculations
  //---------------------------------------------------------------------------

  // Conversion from Stokes to sensor_pol
  ArrayOfVector   s2p;
  stokes2pol( s2p, 0.5 );

  ArrayOfIndex npolcum(nf+1); npolcum[0]=0;
  for( Index i=0; i<nf; i++ )
    { 
      npolcum[i+1] = npolcum[i] + sensor_pol_array[i].nelem(); 
      for( Index j=0; j<sensor_pol_array[i].nelem(); j++ )
        {
          if( s2p[sensor_pol_array[i][j]-1].nelem() > stokes_dim )
            throw runtime_error( "Your definition of *sensor_pol_array* " 
                                 "requires a higher value for *stokes_dim*." );
        }
    }

  // Size output arguments, and set to 0
  const Index ntot = npos * npolcum[nf] * nbins;
  y.resize( ntot );
  y = 0;
  //
  y_aux.resize( naux );
  for( Index i=0; i<naux; i++ )
    { 
      y_aux[i].resize( ntot ); 
      y_aux[i] = 0; 
    }


  // Are range_bins z or t?

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
      //
      iy_main_agendaExecute( ws, iy, iy_aux, ppath, diy_dx, 
                             1, iy_unit, iy_transmission, 
                             iy_aux_vars, cloudbox_on, 0, t_field, z_field, 
                             vmr_field, f_grid, 
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
      if( is_z)
        { range = ppath.pos(joker,0); }
      else
        { // Calculate round-trip time
          range[0] = 2 * ppath.end_lstep / SPEED_OF_LIGHT;
          for( Index i=1; i<np; i++ )
            { range[i] = range[i-1] + ppath.lstep[i-1] * 
                      ( ppath.ngroup[i-1]+ppath.ngroup[i] ) / SPEED_OF_LIGHT; }
        }

      // Loop radar bins
      for( Index b=0; b<nbins; b++ )
        {
          // Get effective limits for bin
          Numeric ztlow  = max( range_bins[b], min(range) );
          Numeric zthigh = min( range_bins[b+1], max(range) );

          if( ztlow < zthigh )   // Otherwise bin outside range of ppath
            {
              // Get ppath altitudes inside bin + edges
              Index n_in = 0;
              ArrayOfIndex i_in(np);
              for( Index i=0; i<np; i++ )
                { if( range[i] > ztlow  &&  range[i] < zthigh )
                    { i_in[n_in] = i; n_in += 1; } }
              Vector zt(n_in+2);
              zt[0] = ztlow;
              for( Index i=0; i<n_in; i++ )
                { zt[i+1] = range[i_in[i]]; }
              zt[n_in+1]  = zthigh;
              n_in      += 2;

              // Height of layer, for reflectivity and aux, respectively
              Numeric dl1 = range_bins[b+1] - range_bins[b];
              Numeric dl2 = zt[n_in-1] - zt[0];

              // Convert to interpolation weights
              ArrayOfGridPos gp( n_in );
              gridpos( gp, range, zt );
              Matrix itw( n_in, 2 );
              interpweights( itw, gp );

              for( Index iv=0; iv<nf; iv++ )
                {
                  // Pick out part of iy for frequency
                  Matrix I = iy( Range(iv*np,np), joker );

                  for( Index ip=0; ip<sensor_pol_array[iv].nelem(); ip++ )
                    {
                      // Extract reflectivity for recieved polarisation
                      Vector refl( np, 0 );
                      Vector w = s2p[sensor_pol_array[iv][ip]-1];
                      for( Index i=0; i<w.nelem(); i++ )     // Note that w can
                        { for( Index j=0; j<np; j++ )        // be shorter than
                            { refl[j] += I(j,i) * w[i]; } }  // stokes_dim (and
                                                             // dot product can 
                      // Interpolate and sum up                 not be used)
                      Vector rv( n_in );
                      interp( rv, itw, refl, gp );
                      Index iout = nbins * ( p*npolcum[nf] + 
                                               npolcum[iv] + ip ) + b;
                      for( Index i=0; i<n_in-1; i++ )
                        { y[iout] += (rv[i]+rv[i+1])*(zt[i+1]-zt[i])/(2*dl1); }

                      // Same for aux variables
                      for( Index a=0; a<naux; a++ )
                        { 
                          if( iy_aux[a].npages()==1 && iy_aux[a].nbooks()==1 &&
                              iy_aux[a].nrows() ==1 && iy_aux[a].ncols() ==np )
                            {
                              interp( rv, itw, iy_aux[a](0,0,0,joker), gp );
                              for( Index i=0; i<n_in-1; i++ )
                                { y_aux[a][iout] += 
                                    (rv[i]+rv[i+1])*(zt[i+1]-zt[i])/(2*dl2); }
                            }
                          else if( iy_aux_vars[a] == "Backscattering" )
                            {
                              // I2 matches I, but is transposed
                              Matrix I2 = iy_aux[a](iv,joker,0,joker);
                              refl = 0;
                              for( Index i=0; i<w.nelem(); i++ )
                                { for( Index j=0; j<np; j++ )
                                    { refl[j] += I2(i,j) * w[i]; } }
                              interp( rv, itw, refl, gp );
                              for( Index i=0; i<n_in-1; i++ )
                                { y_aux[a][iout] += (rv[i]+rv[i+1])*
                                                    (zt[i+1]-zt[i])/(2*dl1); }
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
    }
}



