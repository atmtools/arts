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
#include "rte.h"
#include "sensor.h"

extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT;



/* Workspace method: Doxygen documentation will be auto-generated */
void iy_transmitterCloudRadar(
        Matrix&        iy,
  const Index&         stokes_dim,
  const Vector&        f_grid,
  const ArrayOfIndex&  sensor_pol,
  const Verbosity&  )
{
  const Index nf = f_grid.nelem();
  
  if( sensor_pol.nelem() != nf )
    throw runtime_error( "The length of *f_grid* and the number of elements "
                         "in *sensor_pol* must be equal." );

  iy.resize( nf, stokes_dim );
  iy = 0;

  ArrayOfVector   s2p;
  stokes2pol( s2p, 1 );

  for( Index i=0; i<nf; i++ )
    {
      for( Index j=0; j<s2p[sensor_pol[i]-1].nelem(); j++ )
        {
          iy(i,j) = s2p[sensor_pol[i]-1][j];
        }
    }
}



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
   const Tensor3&                     edensity_field,
   const Index&                       cloudbox_on,
   const ArrayOfIndex&                cloudbox_limits,
   const Tensor4&                     pnd_field,
   const ArrayOfSingleScatteringData& scat_data_raw,
   const Matrix&                      particle_masses,
   const String&                      iy_unit,
   const ArrayOfString&               iy_aux_vars,
   const Index&                       jacobian_do,
   const Agenda&                      ppath_agenda,
   const Agenda&                      abs_mat_per_species_agenda,
   const Agenda&                      iy_transmitter_agenda,
   const Index&                       iy_agenda_call1,
   const Tensor3&                     iy_transmission,
   const Vector&                      rte_pos,      
   const Vector&                      rte_los,      
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
  ppath_agendaExecute( ws, ppath, rte_pos, rte_los, Vector(0), 0, 0,
                       t_field, z_field, vmr_field, edensity_field, f_grid, 
                       ppath_agenda );

  // Some basic sizes
  //
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;


  //=== iy_aux part ===========================================================
  Index auxPressure     = -1,
        auxTemperature  = -1,
        auxBackScat     = -1;
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
      else if( iy_aux_vars[i] == "Backscattering" )
        { auxBackScat = i;  iy_aux[i].resize( nf, ns, 1, np ); iy_aux[i] = 0; }
      else if( iy_aux_vars[i].substr(0,18) == "Particle content, " )
        { 
          Index icont;
          istringstream is(iy_aux_vars[i].substr(18,2));
          is >> icont;
          if( icont < 0  ||  icont>=particle_masses.ncols() )
            {
              ostringstream os;
              os << "You have selected particle content category with index "
                 << icont << ".\nThis category is not defined!";
              throw runtime_error( os.str() );
            }
          auxPartCont.push_back(i);
          auxPartContI.push_back(icont);
          iy_aux[i].resize( 1, 1, 1, np );
        }
      else if( iy_aux_vars[i].substr(0,16) == "Particle field, " )
        { 
          Index ip;
          istringstream is(iy_aux_vars[i].substr(16,2));
          is >> ip;
          if( ip < 0  ||  ip>=pnd_field.nbooks() )
            {
              ostringstream os;
              os << "You have selected particle field with index "
                 << ip << ".\nThis field is not defined!";
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
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w,
                              ppath_mag_u,  ppath_mag_v,  ppath_mag_w;
  Matrix    ppath_vmr, ppath_pnd;
  Tensor5   ppath_abs;
  Tensor4   trans_partial, trans_cumulat, pnd_ext_mat;
  Vector    scalar_tau;
  ArrayOfIndex clear2cloudbox;
  Array<ArrayOfSingleScatteringData> scat_data;
  //
  if( np > 1 )
    {
      get_ppath_atmvars(   ppath_p, ppath_t, ppath_vmr,
                           ppath_wind_u, ppath_wind_v, ppath_wind_w,
                           ppath_mag_u,  ppath_mag_v,  ppath_mag_w,
                           ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                           wind_u_field, wind_v_field, wind_w_field ,
                           mag_u_field, mag_v_field, mag_w_field );      
      get_ppath_abs(       ws, ppath_abs, abs_mat_per_species_agenda, ppath, 
                           ppath_p, ppath_t, ppath_vmr, 
                           ppath_wind_u, ppath_wind_v, ppath_wind_w, 
                           ppath_mag_u, ppath_mag_v, ppath_mag_w, 
                           f_grid, stokes_dim, atmosphere_dim );
      if( !cloudbox_on )
        { get_ppath_trans( trans_partial, trans_cumulat, scalar_tau, 
                           ppath, ppath_abs, f_grid, stokes_dim );
        }
      else
        {
          // Extract basic scattering data
          Index        use_mean_scat_data = 0;
          Tensor3      pnd_abs_vec;
          //
          get_ppath_ext( clear2cloudbox, pnd_abs_vec, pnd_ext_mat, scat_data, 
                         ppath_pnd, ppath, ppath_t, stokes_dim, f_grid, 
                         atmosphere_dim, cloudbox_limits, pnd_field, 
                         use_mean_scat_data, scat_data_raw, verbosity );
          // Transmission
          get_ppath_trans2( trans_partial, trans_cumulat, scalar_tau, 
                            ppath, ppath_abs, f_grid, stokes_dim,
                            clear2cloudbox, pnd_ext_mat );
        }
    }


  // Do RT calculations
  //
  iy.resize( nf*np, ns );
  iy = 0;
  //
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
                                  los_inc[1], scat_data[iv], stokes_dim, 
                                  ppath_pnd(joker,ip), ppath_t[ip], verbosity );

              // Combine iy0, double ransmission and scattering matrix
              Vector iy1(stokes_dim), iy2(stokes_dim);
              mult( iy1, trans_cumulat(iv,joker,joker,ip), iy0(iv,joker) );
              mult( iy2, P,                                iy1 );
              mult( iy(iv*np+ip,joker), 
                         trans_cumulat(iv,joker,joker,ip), iy2 );

              //=== iy_aux part ===========================================
              // Backscattering
              if( auxBackScat >= 0 ) {
                mult( iy_aux[auxBackScat](iv,joker,0,ip), P, iy0(iv,joker) ); }
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
      // Particle field
      for( Index j=0; j<auxPartField.nelem(); j++ )
        { iy_aux[auxPartField[j]](0,0,0,ip) = ppath_pnd(auxPartFieldI[j],ip); }
      //===================================================================
    } 

  if( iy_unit == "1" )
    {}
  else if( iy_unit == "Ze" )
    {
      // Get refractive index for water
      Matrix complex_n;
      complex_nWaterLiebe93( complex_n, f_grid, ze_tref, verbosity );

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
   const Index&                 basics_checked,
   const Index&                 atmosphere_dim,
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
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );

  // Sensor position and LOS.
  if( sensor_pos.ncols() != atmosphere_dim )
    throw runtime_error( "The number of columns of sensor_pos must be "
                         "equal to the atmospheric dimensionality." );
  if( atmosphere_dim <= 2  &&  sensor_los.ncols() != 1 )
    throw runtime_error( "For 1D and 2D, sensor_los shall have one column." );
  if( atmosphere_dim == 3  &&  sensor_los.ncols() != 2 )
    throw runtime_error( "For 3D, sensor_los shall have two columns." );
  if( sensor_los.nrows() != npos )
    {
      ostringstream os;
      os << "The number of rows of sensor_pos and sensor_los must be "
         << "identical, but sensor_pos has " << npos << " rows,\n"
         << "while sensor_los has " << sensor_los.nrows() << " rows.";
      throw runtime_error( os.str() );
    }
  if( max( sensor_los(joker,0) ) > 180 )
    throw runtime_error( 
     "First column of *sensor_los* is not allowed to have values above 180." );
  if( atmosphere_dim == 2 )
    {
      if( min( sensor_los(joker,0) ) < -180 )
          throw runtime_error( "For atmosphere_dim = 2, first column of "
                    "*sensor_los* is not allowed to have values below -180." );
    }     
  else
    {
      if( min( sensor_los(joker,0)  ) < 0 )
          throw runtime_error( "For atmosphere_dim != 2, first column of "
                       "*sensor_los* is not allowed to have values below 0." );
    }    
  if( atmosphere_dim == 3  &&  max( sensor_los(joker,1) ) > 180 )
    throw runtime_error( 
    "Second column of *sensor_los* is not allowed to have values above 180." );

  // Method specific variables
  if( !is_increasing( range_bins ) )
    throw runtime_error( "The vector *range_bins* must contain strictly "
                         "increasing values." );
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
  y.resize( npos*npolcum[nf]*nbins );
  y = 0;
  y_aux.resize( naux );
  for( Index i=0; i<naux; i++ )
    { 
      y_aux[i].resize( npos*nbins ); 
      y_aux[i] = 0; 
    }


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
      iy_main_agendaExecute( ws, iy, iy_aux, ppath, diy_dx, 1, iy_transmission, 
                             iy_aux_vars, cloudbox_on, 0, t_field, z_field, 
                             vmr_field, f_grid, 
                             sensor_pos(p,joker), sensor_los(p,joker), 
                             rte_pos2, iy_main_agenda );

      // Check if path and size OK
      if( ppath.np == 1 )
        throw runtime_error( "A path consisting of a single point found. "
                             "This is not allowed." );
      const bool isupward = ppath.pos(1,0) > ppath.pos(0,0);
      if( !( ( isupward & is_increasing( ppath.pos(joker,0) ) ) ||
             ( !isupward & is_decreasing( ppath.pos(joker,0) ) ) ) )
        throw runtime_error( "A path with strictly increasing or decreasing "
                             "altitudes must be obtained (e.g. limb "
                             "measurements are not allowed)." );
      if( iy.nrows() != nf*ppath.np )
        throw runtime_error( "The size of *iy* returned from *iy_main_agenda* "
                             "is not correct (for this method)." );

      // Loop radar bins
      for( Index b=0; b<nbins; b++ )
        {
          // Get effective limits for bin
          Numeric zlow  = max( range_bins[b], min(ppath.pos(joker,0)) );
          Numeric zhigh = min( range_bins[b+1], max(ppath.pos(joker,0)) );

          if( zlow < zhigh )   // Otherwise bin outside range of ppath
            {
              // Get ppath altitudes inside bin + edges
              Index n_in = 0;
              ArrayOfIndex i_in(ppath.np);
              for( Index i=0; i<ppath.np; i++ )
                { if( ppath.pos(i,0) > zlow  &&  ppath.pos(i,0) < zhigh )
                    { i_in[n_in] = i; n_in += 1; } }
              Vector z(n_in+2);
              z[0] = zlow;
              for( Index i=0; i<n_in; i++ )
                { z[i+1] = ppath.pos(i_in[i],0); }
              z[n_in+1]  = zhigh;
              n_in      += 2;

              // Height of layer, for reflectivity and aux, respectively
              Numeric dl1 = range_bins[b+1] - range_bins[b];
              Numeric dl2 = z[n_in-1] - z[0];

              // Convert to interpolation weights
              ArrayOfGridPos gp( n_in );
              gridpos( gp, ppath.pos(joker,0), z );
              Matrix itw( n_in, 2 );
              interpweights( itw, gp );

              for( Index iv=0; iv<nf; iv++ )
                {
                  // Pick out part of iy for frequency
                  Matrix I = iy( Range(iv*ppath.np,ppath.np), joker );

                  for( Index ip=0; ip<sensor_pol_array[iv].nelem(); ip++ )
                    {
                      // Extract reflectivity for recieved polarisation
                      Vector refl( ppath.np, 0 );
                      Vector w = s2p[sensor_pol_array[iv][ip]-1];
                      for( Index i=0; i<w.nelem(); i++ )
                        { for( Index j=0; j<ppath.np; j++ )
                            { refl[j] += I(j,i) * w[i]; } }

                      // Interpolate and sum up
                      Vector rv( n_in );
                      interp( rv, itw, refl, gp );
                      Index iout = nbins * ( p*npolcum[nf] + 
                                               npolcum[iv] + ip ) + b;
                      for( Index i=0; i<n_in-1; i++ )
                        { y[iout] += (rv[i]+rv[i+1])*(z[i+1]-z[i])/(2*dl1); }

                      // Same for aux variables (if first freq and pol)
                      if( iv==0  && ip==0 )
                        {
                          for( Index a=0; a<naux; a++ )
                            { 
                              if( iy_aux[a].nrows()>1 || iy_aux[a].npages()>1 ||
                                  iy_aux[a].nbooks()>1 )
                                throw runtime_error( "Auxilary variables must "
                                                   "have size 1 for frequency "
                                                   "and Stokes dimensions.");
                              interp( rv, itw, iy_aux[a](0,0,0,joker), gp );
                              iout = nbins*p + b;
                              for( Index i=0; i<n_in-1; i++ )
                                { y_aux[a][iout] += 
                                       (rv[i]+rv[i+1])*(z[i+1]-z[i])/(2*dl2); }
                            }
                        }
                    }
                }
            }
        }
    }
}



