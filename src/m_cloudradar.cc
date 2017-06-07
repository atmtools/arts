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
void pnd_size_gridFromScatMeta(
          Vector&                             pnd_size_grid,
    const ArrayOfArrayOfScatteringMetaData&   scat_meta,
    const Index&                              scat_index,
    const String&                             unit,
    const Verbosity& )
{
  // Sizes
  const Index nss = scat_meta.nelem();

  // Checks
  if( nss == 0 )
    throw runtime_error( "*scat_meta* is empty!" );
  if( nss < scat_index+1 )
    {
      ostringstream os;
      os << "Selected scattering species index is " << scat_index << " but this "
         << "is not allowed since *scat_meta* has only" << scat_meta.nelem()
         << " elements.";
      throw runtime_error(os.str());
    }

  // Create size grid
  //
  const Index nse = scat_meta[scat_index].nelem();
  //
  pnd_size_grid.resize( nse );
  //
  if( unit == "dveq" )
    {
      for ( Index i=0; i<nse; i++ )
        { pnd_size_grid[i] = scat_meta[scat_index][i].diameter_volume_equ; }
    }
  else if( unit == "dmax" )
    {
      for ( Index i=0; i<nse; i++ )
        { pnd_size_grid[i] = scat_meta[scat_index][i].diameter_max; }
    }
  else if( unit == "mass" )
    {
      for ( Index i=0; i<nse; i++ )
        { pnd_size_grid[i] = scat_meta[scat_index][i].mass; }
    }
  else if( unit == "area" )
    {
      for ( Index i=0; i<nse; i++ )
        { pnd_size_grid[i] = scat_meta[scat_index][i].diameter_area_equ_aerodynamical; }
    }
  else
    {
      ostringstream os;
      os << "You have selected the unit: " << unit 
         << "while accepted choices are: \"dveq\", \"dmax\", \"mass\" and \"area\"";
      throw runtime_error(os.str());
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void pndFromPsdBasic(
         Matrix&    pnd_data,
         Tensor3&   dpnd_data_dx,
   const Vector&    pnd_size_grid,
   const Matrix&    psd_data,
   const Vector&    psd_size_grid,
   const Tensor3&   dpsd_data_dx,
   const Verbosity& )
{
  // This code is just for testing. Shall be replaced with code using
  // piece-wise linear representation when mapping psd to pnd.
  
  // Some sizes 
  const Index np  = psd_data.nrows();
  const Index ng  = psd_size_grid.nelem();
        Index ndx = 0;
  const bool  do_dx = !dpsd_data_dx.empty();

  // Checks
  if( ng < 2 )
    throw runtime_error( "The method requires that length of *psd_size_grid* is >= 2." );
  if( ng != pnd_size_grid.nelem() )
    throw runtime_error( "The method requires that *psd_size_grid* and "
                         "*pnd_size_grid* have same length." );
  for( Index i=0; i<ng; i++ )
    {
      if( psd_size_grid[i] != pnd_size_grid[i] )
        throw runtime_error( "The method requires that *psd_size_grid* and "
                             "*pnd_size_grid* are identical." );
    }
  if( psd_data.ncols() != ng )
    throw runtime_error( "Number of columns in *psd_data* and length of "
                         "*psd_size_grid* must match." );
  if( do_dx )
    {
      if( dpsd_data_dx.ncols() != ng )
        throw runtime_error( "Number of columns in *dpsd_data_dx* and length of "
                             "*psd_size_grid* must match." );
      ndx = dpsd_data_dx.npages();
      dpnd_data_dx.resize( ndx, np, ng );
    }
  else
    { dpnd_data_dx.resize( 0, 0, 0 ); }

  if( !is_increasing( psd_size_grid ) )
    throw runtime_error( "*psd_size_grid* must be strictly increasing." );    


  // Also, is my scheme to define bin sizes the same as you used?
  // Asking as I get deviating results, and I suspect it could come from this
  // part 

  // dpnd_dx is sized above    
  pnd_data.resize( np, ng );

  
  // Calculate
  Numeric binsize;
  for ( Index i=0; i<ng; i++ )
    {
      // This bin is twice the half-distance to point 1, but could be limited
      // by 0 in lower end.
      if( i == 0 )
        {
          const Numeric dd = ( psd_size_grid[1] - psd_size_grid[0] ) / 2.0;
          binsize = dd + min( dd, psd_size_grid[0]/ 2.0 );
        }
      // This bin is twice the half-distance to closest point      
      else if( i == ng-1 )
        { binsize = 2.0 * ( psd_size_grid[i] - psd_size_grid[i-1] ) / 2.0; }
      // This bin is the sum of the two half-distances      
      else
        { binsize = ( psd_size_grid[i+1] - psd_size_grid[i-1] ) / 2.0; }

      for( Index ip=0; ip<np; ip++ )
        { pnd_data(ip,i) = binsize * psd_data(ip,i); }

      if( do_dx )
        {
          for( Index ip=0; ip<np; ip++ )
            {
              for ( Index ix=0; ix<ndx; ix++ )
                { dpnd_data_dx(ix,ip,i) = binsize * dpsd_data_dx(ix,ip,i); }
            }
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void psdMH97 (
          Matrix&                             psd_data,
          Tensor3&                            dpsd_data_dx,
    const Vector&                             psd_size_grid,
    const Matrix&                             pnd_agenda_input,
    const ArrayOfString&                      pnd_agenda_input_names,
    const ArrayOfString&                      dpnd_data_dx_names,
    const Numeric&                            t_min, 
    const Numeric&                            t_max, 
    const Index&                              picky, 
    const Index&                              noisy,
    const Verbosity&)
{
  // Some sizes 
  const Index nin = pnd_agenda_input_names.nelem();
  const Index ndx = dpnd_data_dx_names.nelem();
  const Index np  = pnd_agenda_input.nrows();
  const Index nsi = psd_size_grid.nelem();
  
  // Checks
  if( pnd_agenda_input.ncols() != nin )
    throw runtime_error( "Length of *pnd_agenda_input_names* and number of "
                         "columns in *pnd_agenda_input* must be equal." );
  if( pnd_agenda_input.ncols() != 2 )
    throw runtime_error( "*pnd_input* must have two columns (IWC and Temperature)" );
  if( pnd_agenda_input_names[0] != "IWC" )
    throw runtime_error( "Element 0 of *pnd_agenda_input_names* must be \"IWC\"." );
  if( pnd_agenda_input_names[1] != "Temperature" )
    throw runtime_error( "Element 1 of *pnd_agenda_input_names* must be "
                         "\"Temperature\"." );
  if( ndx )
    {
      if( ndx != 1 )
        throw runtime_error( "*dpnd_data_dx_names* must have length 0 or 1." );
      if( dpnd_data_dx_names[0] != "IWC" )
        throw runtime_error( "With MH97, the only valid option for "
                             "*dpnd_data_dx_names* is: \"IWC\"." );
    }        
  if( noisy  &&   ndx )
    throw runtime_error( "Jacobian calculations and \"noisy\" can not be "
                         "combined." );

  
  // Init psd_data and dpsd_data_dx with zeros
  psd_data.resize( np, nsi );
  psd_data = 0.0;
  if( ndx )
    {
      dpsd_data_dx.resize( 1, np, nsi );   // IWC only possible retrieval quantity
      dpsd_data_dx = 0.0;
    }
  else
    { dpsd_data_dx.resize( 0, 0, 0  ); }  

  for( Index ip=0; ip<np; ip++ )
    {
      
      // Extract the input variables
      Numeric iwc = pnd_agenda_input(ip,0);
      Numeric   t = pnd_agenda_input(ip,1);

      // Negative iwc?
      Numeric psd_weight = 1.0;
      if( iwc < 0 )
        {
          psd_weight = -1.0;
          iwc       *= -1.0;
        }
      
      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: ["
                 << t_min << "," << t_max << "]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }

      // PSD assumed to be constant outside [200,273.15]
      // Shall these temperature limits be GIN?
      if( t < 200 )
        { t = 200; }
      else if( t > 273.15 )
        { t = 273.15; }
  
      // To Patrick & Oliver: Shall the core psd function (psd_cloudice_MH97)
      // handle even more input (i.e. Vectors of IWC and T) for calculation
      // overhead avoidance?

      Vector psd_1p(nsi);
      // Calculate PSD
      if( iwc != 0 )
        {
          psd_cloudice_MH97 ( psd_1p, psd_size_grid, iwc, t, noisy );
          for ( Index i=0; i<nsi; i++ )
            { psd_data(ip,i) = psd_weight * psd_1p[i]; }
        }

      // Calculate derivative with respect to IWC
      if( ndx )
        {
          // Obtain derivative by perturbation of 0.1%, but not less than 0.1 mg/m3.
          // Note that the last value becomes the perturbation for IWC=0.
          const Numeric diwc = max( 0.001*iwc, 1e-7 );
          const Numeric iwcp = iwc + diwc;
          psd_cloudice_MH97 ( psd_1p, psd_size_grid, iwcp, t, noisy );
          for ( Index i=0; i<nsi; i++ )
            { dpsd_data_dx(0,ip,i) = ( psd_1p[i] - psd_data(ip,i) ) / diwc; }
        }   
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldCalcFromParticleBulkProps(
         Workspace&                   ws,
         Tensor4&                     pnd_field,
         ArrayOfTensor4&              dpnd_field_dx,
   const Index&                       atmosphere_dim,
   const Vector&                      p_grid,
   const Vector&                      lat_grid,
   const Vector&                      lon_grid,
   const Tensor3&                     t_field,
   const Index&                       cloudbox_on,
   const ArrayOfIndex&                cloudbox_limits,
   const ArrayOfString&               scat_species,
   const ArrayOfArrayOfSingleScatteringData& scat_data,
   const ArrayOfArrayOfScatteringMetaData&   scat_meta,
   const Tensor4&                     particle_bulkprop_field,
   const ArrayOfString&               particle_bulkprop_names,
   const ArrayOfAgenda&               pnd_agenda_array,
   const ArrayOfArrayOfString&        pnd_agenda_array_input_names,
   const Index&                       jacobian_do,
   const ArrayOfRetrievalQuantity&    jacobian_quantities,
   const Verbosity&)
{
  // Number of scattering species
  const Index nss = scat_data.nelem();

  // Checks (not totally complete, but should cover most mistakes)
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_field( "t_field", t_field, atmosphere_dim, 
                 p_grid, lat_grid, lon_grid );
  chk_atm_field( "particle_bulkprop_field", particle_bulkprop_field,
                 atmosphere_dim, particle_bulkprop_names.nelem(),
                 p_grid, lat_grid, lon_grid );
  // Further checks of *particle_bulkprop_field* below
  if( !cloudbox_on )
    throw runtime_error( "*cloudbox_on* must be true to use this method." );
  if( cloudbox_limits.nelem() != 2*atmosphere_dim )
    throw runtime_error( "Length of *cloudbox_limits* incorrect with respect "
                         "to *atmosphere_dim*." );
  if( cloudbox_limits[1]<=cloudbox_limits[0] || cloudbox_limits[0]<0 ||
                                           cloudbox_limits[1]>=p_grid.nelem() )
    throw runtime_error( "Invalid data in pressure part of *cloudbox_limits*." );
  if( atmosphere_dim > 1 )
    {
      if( cloudbox_limits[3]<=cloudbox_limits[2] || cloudbox_limits[2]<0 ||
                                           cloudbox_limits[3]>=lat_grid.nelem() )
        throw runtime_error( "Invalid data in latitude part of *cloudbox_limits*." );
      if( atmosphere_dim > 2 )
        {
          if( cloudbox_limits[5]<=cloudbox_limits[4] || cloudbox_limits[4]<0 ||
                                           cloudbox_limits[5]>=lon_grid.nelem() )
            throw runtime_error( "Invalid data in longitude part of *cloudbox_limits*." );
        }
    }
  if( nss < 1 )
    throw runtime_error( "*scat_data* is empty!." );
  if( scat_species.nelem() != nss )
    throw runtime_error( "*scat_data* and *scat_species* are inconsistent in size." );
  if( scat_meta.nelem() != nss )
    throw runtime_error( "*scat_data* and *scat_meta* are inconsistent in size." );
  if( pnd_agenda_array.nelem() != nss )
    throw runtime_error( "*scat_data* and *pnd_agenda_array* are inconsistent "
                         "in size." );
  if( pnd_agenda_array_input_names.nelem() != nss )
    throw runtime_error( "*scat_data* and *pnd_agenda_array_input_names* are "
                         "inconsistent in size." );
  // Further checks of scat_data vs. scat_meta below  

  
  // Effective lengths of cloudbox
  const Index np = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  const Index ip_offset = cloudbox_limits[0];
  Index nlat = 1;
  Index ilat_offset = 0;
  if( atmosphere_dim > 1 )
    {
      nlat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
      ilat_offset = cloudbox_limits[2];
    }
  Index nlon = 1;
  Index ilon_offset = 0;
  if( atmosphere_dim > 2 )
    {
      nlat = cloudbox_limits[5] - cloudbox_limits[4] + 1;
      ilon_offset = cloudbox_limits[4];
    }

  // Check that *particle_bulkprop_field* contains zeros outside and at
  // cloudbox boundaries
  const String estring = "*particle_bulkprop_field* can only contain non-zero "
    "values inside the cloudbox.";
  // Pressure end ranges
  for( Index ilon=0; ilon<nlon; ilon++ ) {
    for( Index ilat=0; ilat<nlat; ilat++ ) {
      for( Index ip=0; ip<=cloudbox_limits[0]; ip++ ) {
        if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
          throw runtime_error( estring ); } 
      for( Index ip=cloudbox_limits[1]; ip<p_grid.nelem(); ip++ ) {
        if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
          throw runtime_error( estring ); } } }
  if( atmosphere_dim > 1 )
    {
      // Latitude end ranges
      for( Index ilon=0; ilon<nlon; ilon++ ) {
        for( Index ip=cloudbox_limits[0]+1; ip<cloudbox_limits[1]-1; ip++ ) {
          for( Index ilat=0; ilat<=cloudbox_limits[2]; ilat++ ) {
            if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
              throw runtime_error( estring ); }
          for( Index ilat=cloudbox_limits[3]; ilat<lat_grid.nelem(); ilat++ ) {
            if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
              throw runtime_error( estring ); } } }
      if( atmosphere_dim > 2 )
        {
          // Longitude end ranges
          for( Index ip=cloudbox_limits[0]+1; ip<cloudbox_limits[1]-1; ip++ ) {
            for( Index ilat=cloudbox_limits[2]+1; ilat<cloudbox_limits[3]-1; ilat++ ) {
              for( Index ilon=0; ilon<=cloudbox_limits[4]; ilon++ ) {
                if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
                  throw runtime_error( estring ); }
              for( Index ilon=cloudbox_limits[5]; ilon<lon_grid.nelem(); ilon++ ) {
                if( max(particle_bulkprop_field(joker,ip,ilat,ilon)) > 0 )
                  throw runtime_error( estring ); } } }
        }
    }
  
  // Cumulative number of scattering elements
  ArrayOfIndex ncumse(nss+1);
  ncumse[0] = 0;
  for( Index i=0; i<nss; i++ )
    {
      if( scat_data[i].nelem() != scat_meta[i].nelem() )
        throw runtime_error( "*scat_data* and *scat_meta* have inconsistent sizes." );
      ncumse[i+1] = ncumse[i+1] + scat_data[i].nelem();
    }

  // Allocate output variables
  //
  pnd_field.resize( ncumse[nss], np, nlat, nlon );
  pnd_field = 0.0;  // To set all end values to zero
  //
  // Help variables for partial derivatives
  Index                nq = 0;
  ArrayOfArrayOfIndex  scatspecies_to_jq;
  //
  if( !jacobian_do )
    { dpnd_field_dx.resize(0); }
  else
    {
      nq = jacobian_quantities.nelem();
      dpnd_field_dx.resize( nq );
      scatspecies_to_jq.resize( nss );
      //
      for( Index iq=0; iq<nq; iq++ )
        {
          if( jacobian_quantities[iq].MainTag() == SCATSPECIES_MAINTAG  )
            {
              const Index ihit = find_first( scat_species,
                                             jacobian_quantities[iq].Subtag() );
              if( ihit < 0 )
                {
                  ostringstream os;
                  os << "Jacobian quantity with index " << iq << " refers to\n"
                     << "  " << jacobian_quantities[iq].Subtag()
                     << "\nbut this species could not be found in *scat_species*.";
                  throw runtime_error(os.str());
                }
              scatspecies_to_jq[ihit].push_back( iq );
              dpnd_field_dx[iq].resize( ncumse[nss], np, nlat, nlon );
              dpnd_field_dx[iq] = 0.0;  // To set all end values to zero
            }
        }
    }

  // Extract data from pnd-agenda array
  for( Index is=0; is<nss; is++ )
    {
      // Index range with respect to pnd_field
      Range se_range( ncumse[is], ncumse[is+1]-ncumse[is] );

      // Determine how pnd_agenda_array_input_names are related to input fields
      //
      const Index nin = pnd_agenda_array_input_names[is].nelem();
      ArrayOfIndex i_pbulkprop(nin);
      //
      for( Index i=0; i<nin; i++ )
        {
          // We flag temperature with -100
          if( pnd_agenda_array_input_names[is][i] == "Temperature" )
            { i_pbulkprop[i] = -100; }
          else
            {
              i_pbulkprop[i] = find_first( particle_bulkprop_names,
                                           pnd_agenda_array_input_names[is][i] );
              if( i_pbulkprop[i] < 0 )
                {
                  ostringstream os;
                  os << "Pnd-agenda with index " << is << " is set to require \""
                     << pnd_agenda_array_input_names[is][i] << "\",\nbut this quantity "
                     << "could not found in *particle_bulkprop_names*.\n"
                     << "(Note that temperature must be written as \"Temperature\")";
                  throw runtime_error(os.str());
                }
            }
        }

      // Set *dpnd_data_dx_names*
      //
      Index ndx = 0;
      ArrayOfString dpnd_data_dx_names( 0 );
      //
      if( jacobian_do )
        {
          ndx = scatspecies_to_jq[is].nelem();
          dpnd_data_dx_names.resize( ndx );
          for( Index ix=0; ix<ndx; ix++ )
            { dpnd_data_dx_names[ix] =
                jacobian_quantities[scatspecies_to_jq[is][ix]].SubSubtag(); }
        }
      
      // Loop lat/lon positions and call *pnd_agenda*
      for( Index ilon=0; ilon<nlon; ilon++ )
        { 
          for( Index ilat=0; ilat<nlat; ilat++ )
            {
              // Note that we don't need any calculations for end points

              // Here we consider this for lat and lon
              if( ( nlat > 1  &&  ( ilat == 0  ||  ilat == nlat-1 ) )  || 
                  ( nlon > 1  &&  ( ilon == 0  ||  ilon == nlon-1 ) ) )
                { continue; }
                  
              // Pressure handled here, by not including end points in loops
              Matrix pnd_agenda_input( np-2, nin );

              
              for( Index i=0; i<nin; i++ )
                {
                  if( i_pbulkprop[i] == -100 )
                    {
                      for( Index ip=1; ip<np-1; ip++ )
                        { pnd_agenda_input(ip-1 ,i) = t_field(
                                                         ip_offset   + ip,
                                                         ilat_offset + ilat,
                                                         ilon_offset + ilon ); }
                    }
                  else
                    {
                      for( Index ip=1; ip<np-1; ip++ )
                        { pnd_agenda_input(ip-1,i) = particle_bulkprop_field(
                                                         i_pbulkprop[i],
                                                         ip_offset   + ip,
                                                         ilat_offset + ilat,
                                                         ilon_offset + ilon ); }
                    }
                }
              
              // Call pnd-agenda array
              Matrix pnd_data;
              Tensor3 dpnd_data_dx;
              //
              pnd_agenda_arrayExecute( ws, pnd_data, dpnd_data_dx, is,
                                       pnd_agenda_input,
                                       pnd_agenda_array_input_names[is],
                                       dpnd_data_dx_names, pnd_agenda_array );

              // Copy to output variables
              for( Index ip=1; ip<np-1; ip++ )
                { pnd_field(se_range,ip,ilat,ilon) = pnd_data(ip-1,joker); }
              for( Index ix=0; ix<ndx; ix++ )
                {
                  for( Index ip=1; ip<np-1; ip++ )
                    { dpnd_field_dx[scatspecies_to_jq[is][ix]]
                            (se_range,ip,ilat,ilon) = dpnd_data_dx(ix,ip-1,joker); }
                }
            }
        }
    }
}
   
                           


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
   const ArrayOfArrayOfIndex&         jacobian_indices,
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
  if( !cloudbox_on )
    throw runtime_error( 
                    "The cloudbox must be activated (cloudbox_on must be 1)." );

  
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
  Tensor4   ppath_ext, trans_partial, trans_cumulat, pnd_ext_mat, dummy_dppath_nlte_dx;
  Tensor5   dummy_abs_per_species, dummy_dppath_ext_dx;
  Vector    scalar_tau;
  ArrayOfIndex   clear2cloudbox, dummy_lte;
  ArrayOfMatrix  ppath_dpnd_dx;
  Array<ArrayOfArrayOfSingleScatteringData> scat_data_single;
  const Tensor4 t_nlte_field_empty(0,0,0,0);
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
      
      get_ppath_pmat( ws, ppath_ext, dummy_ppath_nlte_source, dummy_lte, 
                      dummy_abs_per_species, 
                      dummy_dppath_ext_dx, dummy_dppath_nlte_dx,
                      propmat_clearsky_agenda, ArrayOfRetrievalQuantity(0), ppath, 
                      ppath_p, ppath_t, ppath_t_nlte, ppath_vmr, ppath_f, ppath_mag,
                      f_grid, stokes_dim, ArrayOfIndex(0) );
      
      // Extract basic scattering data
      ArrayOfArrayOfIndex  extmat_case;
      Tensor3              pnd_abs_vec;
      //
      get_ppath_ext( clear2cloudbox, pnd_abs_vec, pnd_ext_mat, scat_data_single,
                     ppath_pnd, ppath_dpnd_dx, ppath, ppath_t, stokes_dim, ppath_f, 
                     atmosphere_dim, cloudbox_limits, pnd_field, dpnd_field_dx,
                     0, scat_data, verbosity );
      
      get_ppath_trans2( trans_partial, extmat_case, trans_cumulat, scalar_tau, 
                        ppath, ppath_ext, f_grid, stokes_dim, 
                        clear2cloudbox, pnd_ext_mat );
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

      
      // Radar return only possible if inside cloudbox
      if( clear2cloudbox[ip] >= 0 )
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
                  // We extract data for pnd=1.
                  //
                  Vector pnd_probe(nsetot); pnd_probe = 1;
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
                          // Scattering species
                          if( jac_scat_i[iq] >= 0 )
                            {
                              // Change of scattering scattering matrix
                              P = 0.0;
                              for( Index i=0; i<nsetot; i++ ) {
                                if( ppath_pnd(i,ip) != 0 ) {
                                  for( Index is1=0; is1<ns; is1++ ) {
                                    for( Index is2=0; is2<ns; is2++ )
                                      { P(is1,is2) += ppath_dpnd_dx[iq](i,ip) *
                                                                Pe(i,is1,is2); }
                               } } }

                              // Apply transmissions as above
                              mult( iy2, P, iy1 );
                              mult( diy_dpath[iq](ip,iv*np+ip,joker),
                                    trans_cumulat(iv,joker,joker,ip), iy2 );
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
   const ArrayOfArrayOfIndex&      jacobian_indices,
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
  //
  if( jacobian_do )
    {
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
