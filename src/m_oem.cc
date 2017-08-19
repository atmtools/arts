/* Copyright (C) 2015
   Patrick Eriksson <patrick.eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>

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
  \file   m_oem.cc
  \author Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2015-09-08

  \brief  Workspace functions related to making OEM inversions.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <iterator>
#include <stdexcept>
#include <string>
#include <sstream>
#include "arts.h"
#include "arts_omp.h"
#include "array.h"
#include "auto_md.h"
#include "math_funcs.h"
#include "physics_funcs.h"

extern const String ABSSPECIES_MAINTAG;
extern const String TEMPERATURE_MAINTAG;
extern const String POINTING_MAINTAG;
extern const String POINTING_SUBTAG_A;
extern const String POLYFIT_MAINTAG;
extern const String SCATSPECIES_MAINTAG;
extern const String SINEFIT_MAINTAG;



/*===========================================================================
  === Help functions for grid handling
  ===========================================================================*/


//! Determines grid positions for regridding of atmospheric fields to retrieval
//  grids
/*!
  The grid positions arrays are sized inside the function. gp_lat is given
  length 0 for atmosphere_dim=1 etc.

  This regridding uses extpolfac=0.

  \param[out] gp_p                 Pressure grid positions.
  \param[out] gp_lat               Latitude grid positions.
  \param[out] gp_lon               Longitude grid positions.
  \param[in]  rq                   Retrieval quantity structure.
  \param[in]  atmosphere_dim       As the WSV with same name.
  \param[in]  p_grid               As the WSV with same name.
  \param[in]  lat_grid             As the WSV with same name.
  \param[in]  lon_grid             As the WSV with same name.

  \author Patrick Eriksson
  \date   2015-09-09
*/
void get_gp_atmgrids_to_rq(
         ArrayOfGridPos&      gp_p,
         ArrayOfGridPos&      gp_lat,
         ArrayOfGridPos&      gp_lon,
   const RetrievalQuantity&   rq,
   const Index&               atmosphere_dim,
   const Vector&              p_grid,
   const Vector&              lat_grid,
   const Vector&              lon_grid )
{
  gp_p.resize( rq.Grids()[0].nelem() );
  p2gridpos( gp_p, p_grid, rq.Grids()[0], 0 );  
  //
  if( atmosphere_dim >= 2 )
    {
      gp_lat.resize( rq.Grids()[1].nelem() );
      gridpos( gp_lat, lat_grid, rq.Grids()[1], 0 );  
    }
  else
    { gp_lat.resize(0); }
  //
  if( atmosphere_dim >= 3 )
    {
      gp_lon.resize( rq.Grids()[2].nelem() );
      gridpos( gp_lon, lon_grid, rq.Grids()[2], 0 );  
    }
  else
    { gp_lon.resize(0); }
}



//! Determines grid positions for regridding of atmospheric fields to retrieval
//  grids 
/*!
  The grid positions arrays are sized inside the function. gp_lat is given
  length 0 for atmosphere_dim=1 etc.

  This regridding uses extpolfac=Inf (where Inf is a very large value).

  Note that the length output arguments (n_p etc.) are for the retrieval grids
  (not the length of grid positions arrays). n-Lat is set to 1 for
  atmosphere_dim=1 etc.

  \param[out] gp_p                 Pressure grid positions.
  \param[out] gp_lat               Latitude grid positions.
  \param[out] gp_lon               Longitude grid positions.
  \param[out] n_p                  Length of retrieval pressure grid.
  \param[out] n_lat                Length of retrieval lataitude grid.
  \param[out] n_lon                Length of retrieval longitude grid.
  \param[in]  rq                   Retrieval quantity structure.
  \param[in]  atmosphere_dim       As the WSV with same name.
  \param[in]  p_grid               As the WSV with same name.
  \param[in]  lat_grid             As the WSV with same name.
  \param[in]  lon_grid             As the WSV with same name.

  \author Patrick Eriksson 
  \date   2015-09-09
*/
void get_gp_rq_to_atmgrids( 
         ArrayOfGridPos&      gp_p,
         ArrayOfGridPos&      gp_lat,
         ArrayOfGridPos&      gp_lon,
         Index&               n_p,
         Index&               n_lat,
         Index&               n_lon,
   const RetrievalQuantity&   rq,
   const Index&               atmosphere_dim,
   const Vector&              p_grid,
   const Vector&              lat_grid,
   const Vector&              lon_grid )
{
  // We want here an extrapolation to infinity -> 
  //                                        extremly high extrapolation factor
  const Numeric inf_proxy = 1.0e99;

  gp_p.resize( p_grid.nelem() );
  n_p = rq.Grids()[0].nelem();
  if( n_p > 1 )
    { 
      p2gridpos( gp_p, rq.Grids()[0], p_grid, inf_proxy ); 
      jacobian_type_extrapol( gp_p );
    }
  else
    { gp4length1grid( gp_p ); }        

  if( atmosphere_dim >= 2 )
    {
      gp_lat.resize( lat_grid.nelem() );
      n_lat = rq.Grids()[1].nelem();
      if( n_lat > 1 )
        { 
          gridpos( gp_lat, rq.Grids()[1], lat_grid, inf_proxy );  
          jacobian_type_extrapol( gp_lat );
        }
      else
        { gp4length1grid( gp_lat ); }        
    }
  else
    { 
      gp_lat.resize(0); 
      n_lat = 1;
    }
  //
  if( atmosphere_dim >= 3 )
    {
      gp_lon.resize( lon_grid.nelem() );
      n_lon = rq.Grids()[2].nelem();
      if( n_lon > 1 )
        { 
          gridpos( gp_lon, rq.Grids()[2], lon_grid, inf_proxy ); 
          jacobian_type_extrapol( gp_lon );
        }
      else
        { gp4length1grid( gp_lon ); }        
    }
  else
    { 
      gp_lon.resize(0); 
      n_lon = 1;
    }
}




/*===========================================================================
  === Workspace methods associated with OEM
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void xaStandard( 
         Vector&                     xa,
   const Index&                      atmfields_checked,
   const Index&                      atmgeom_checked,
   const ArrayOfRetrievalQuantity&   jq,
   const ArrayOfArrayOfIndex&        ji,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const Index&                      cloudbox_on,
   const Index&                      cloudbox_checked,
   const Tensor4&                    particle_bulkprop_field,
   const ArrayOfString&              particle_bulkprop_names,         
   const Verbosity&)
{
  // Basics
  //
  if( atmfields_checked != 1 )
    throw runtime_error( "The atmospheric fields must be flagged to have "
                         "passed a consistency check (atmfields_checked=1)." );
  if( atmgeom_checked != 1 )
    throw runtime_error( "The atmospheric geometry must be flagged to have "
                         "passed a consistency check (atmgeom_checked=1)." );
  if( cloudbox_checked != 1 )
    throw runtime_error( "The cloudbox must be flagged to have "
                         "passed a consistency check (cloudbox_checked=1)." );

  // Sizes
  const Index nq = jq.nelem();
  //
  xa.resize( ji[nq-1][1]+1 );

  // Loop retrieval quantities and fill *xa*
  for( Index q=0; q<jq.nelem(); q++ )
    {
      // Index range of this retrieval quantity
      const Index np = ji[q][1] - ji[q][0] + 1;
      Range ind( ji[q][0], np );

      
      // Atmospheric temperatures
      if( jq[q].MainTag() == TEMPERATURE_MAINTAG )
        {
          // Here we need to interpolate *t_field*
          ArrayOfGridPos gp_p, gp_lat, gp_lon;
          get_gp_atmgrids_to_rq( gp_p, gp_lat, gp_lon, jq[q], atmosphere_dim,
                                 p_grid, lat_grid, lon_grid );
          Tensor3 t_x(gp_p.nelem(),gp_lat.nelem(),gp_lon.nelem());
          regrid_atmfield_by_gp( t_x, atmosphere_dim, t_field, 
                                 gp_p, gp_lat, gp_lon );
          flat( xa[ind], t_x );
        }

      
      // Abs species
      else if( jq[q].MainTag() == ABSSPECIES_MAINTAG )
        {
          // Index position of species
          ArrayOfSpeciesTag  atag;
          array_species_tag_from_string( atag, jq[q].Subtag() );
          const Index isp = chk_contains( "abs_species", abs_species, atag );

          if( jq[q].Mode() == "rel" )
            {
              // This one is simple, just a vector of ones
              xa[ind] = 1; 
            }
          else if( jq[q].Mode() == "logrel" )
            {
              // Also simple, just a vector of zeros
              xa[ind] = 0; 
            }
          else if( jq[q].Mode() == "vmr" )
            {
              // Here we need to interpolate *vmr_field*
              ArrayOfGridPos gp_p, gp_lat, gp_lon;
              get_gp_atmgrids_to_rq( gp_p, gp_lat, gp_lon, jq[q], atmosphere_dim,
                                     p_grid, lat_grid, lon_grid );
              Tensor3 vmr_x(gp_p.nelem(),gp_lat.nelem(),gp_lon.nelem());
              regrid_atmfield_by_gp( vmr_x, atmosphere_dim, 
                                     vmr_field(isp,joker,joker,joker), 
                                     gp_p, gp_lat, gp_lon );
              flat( xa[ind], vmr_x );
            }
          else if( jq[q].Mode() == "nd" )
            {
              // Here we need to interpolate both *vmr_field* and *t_field*
              ArrayOfGridPos gp_p, gp_lat, gp_lon;
              get_gp_atmgrids_to_rq( gp_p, gp_lat, gp_lon, jq[q], atmosphere_dim,
                                     p_grid, lat_grid, lon_grid );
              Tensor3 vmr_x(gp_p.nelem(),gp_lat.nelem(),gp_lon.nelem());
              Tensor3 t_x(gp_p.nelem(),gp_lat.nelem(),gp_lon.nelem());
              regrid_atmfield_by_gp( vmr_x, atmosphere_dim, 
                                     vmr_field(isp,joker,joker,joker), 
                                     gp_p, gp_lat, gp_lon );
              regrid_atmfield_by_gp( t_x, atmosphere_dim,  t_field,
                                     gp_p, gp_lat, gp_lon );
              // Calculate number density for species (vmr*nd_tot)
              Index i = 0;
              for( Index i3=0; i3<=vmr_x.ncols(); i3++ )
                { for( Index i2=0; i2<=vmr_x.nrows(); i2++ )
                    { for( Index i1=0; i1<=vmr_x.npages(); i1++ )
                        { 
                          xa[ji[q][0]+i] = vmr_x(i1,i2,i3) *
                            number_density( jq[q].Grids()[0][i1], t_x(i1,i2,i3) );
                          i += 1;
                }   }   }
            }
          else
            { assert(0); }
        }

      
      // Scattering species
      else if( jq[q].MainTag() == SCATSPECIES_MAINTAG )
        {
          if( cloudbox_on )
            {
              if( particle_bulkprop_field.empty() )
                {
                  throw runtime_error( "One jacobian quantity belongs to the "
                    "scattering species category, but *particle_bulkprop_field* "
                    "is empty." );
                }
              if( particle_bulkprop_field.nbooks() != particle_bulkprop_names.nelem() )
                {
                  throw runtime_error( "Mismatch in size between "
                    "*particle_bulkprop_field* and *particle_bulkprop_field*." );
                }

              const Index isp = find_first( particle_bulkprop_names,
                                        jq[q].SubSubtag() );
              if( isp < 0 )
                {
                  ostringstream os;
                  os << "Jacobian quantity with index " << q << " covers a "
                     << "scattering species, and the field quantity is set to \""
                     << jq[q].SubSubtag() << "\", but this quantity "
                     << "could not found in *particle_bulkprop_names*.";
                  throw runtime_error(os.str());
                }
              
              ArrayOfGridPos gp_p, gp_lat, gp_lon;
              get_gp_atmgrids_to_rq( gp_p, gp_lat, gp_lon, jq[q], atmosphere_dim,
                                     p_grid, lat_grid, lon_grid );
              Tensor3 pbp_x(gp_p.nelem(),gp_lat.nelem(),gp_lon.nelem());
              regrid_atmfield_by_gp( pbp_x, atmosphere_dim, 
                                     particle_bulkprop_field(isp,joker,joker,joker), 
                                     gp_p, gp_lat, gp_lon );
              flat( xa[ind], pbp_x );
            }
          else
            { xa[ind] = 0; }
        }

      // All variables having zero as a priori
      // ----------------------------------------------------------------------------
      else if( jq[q].MainTag() == POINTING_MAINTAG ||
               jq[q].MainTag() == POLYFIT_MAINTAG  ||  
               jq[q].MainTag() == SINEFIT_MAINTAG )
        {
          xa[ind] = 0;
        }

      else
        {
          ostringstream os;
          os << "Found a retrieval quantity that is not yet handled by\n"
             << "internal retrievals: " << jq[q].MainTag() << endl;
          throw runtime_error(os.str());
        }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void x2artsStandard(
         Vector&                     y_baseline,
         Tensor4&                    vmr_field,
         Tensor3&                    t_field,
         Tensor4&                    particle_bulkprop_field,
         Matrix&                     sensor_los,
   const Index&                      atmfields_checked,
   const Index&                      atmgeom_checked,
   const Matrix&                     jacobian,
   const ArrayOfRetrievalQuantity&   jq,
   const ArrayOfArrayOfIndex&        ji,
   const Vector&                     x,
   const Vector&                     xa,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const Index&                      cloudbox_on,
   const Index&                      cloudbox_checked,
   const ArrayOfString&              particle_bulkprop_names,         
   const Vector&                     sensor_time,
   const Verbosity&)
{
  // Basics
  //
  if( atmfields_checked != 1 )
    throw runtime_error( "The atmospheric fields must be flagged to have "
                         "passed a consistency check (atmfields_checked=1)." );
  if( atmgeom_checked != 1 )
    throw runtime_error( "The atmospheric geometry must be flagged to have "
                         "passed a consistency check (atmgeom_checked=1)." );
  if( cloudbox_checked != 1 )
    throw runtime_error( "The cloudbox must be flagged to have "
                         "passed a consistency check (cloudbox_checked=1)." );

  // Main sizes
  const Index nq = jq.nelem();

  // Check input
  if( x.nelem() != ji[nq-1][1]+1 )
    throw runtime_error( "Length of *x* does not match information in "
                         "*jacobian_quantities*.");

  // Note that when this method is called, vmr_field and other output variables
  // have original values, i.e. matching the a priori state.

  // Flag indicating that y_baseline is not set
  bool yb_set = false;

  // Loop retrieval quantities and fill *xa*
  for( Index q=0; q<nq; q++ )
    {
      // Index range of this retrieval quantity
      const Index np = ji[q][1] - ji[q][0] + 1;
      Range ind( ji[q][0], np );

      
      // Atmospheric temperatures
      // ----------------------------------------------------------------------------
      if( jq[q].MainTag() == TEMPERATURE_MAINTAG )
        {
          // Determine grid positions for interpolation from retrieval grids back
          // to atmospheric grids
          ArrayOfGridPos gp_p, gp_lat, gp_lon;
          Index          n_p, n_lat, n_lon;
          get_gp_rq_to_atmgrids( gp_p, gp_lat, gp_lon, n_p, n_lat, n_lon,
                                 jq[q], atmosphere_dim, p_grid, lat_grid, lon_grid );

          // Map values in x back to t_field
          Tensor3 t_x( n_p, n_lat, n_lon );
          reshape( t_x, x[ind] );
          regrid_atmfield_by_gp( t_field, atmosphere_dim, t_x,
                                 gp_p, gp_lat, gp_lon );
        }

      
      // Abs species
      // ----------------------------------------------------------------------------
      else if( jq[q].MainTag() == ABSSPECIES_MAINTAG )
        {
          // Index position of species
          ArrayOfSpeciesTag  atag;
          array_species_tag_from_string( atag, jq[q].Subtag() );
          const Index isp = chk_contains( "abs_species", abs_species, atag );

          // Determine grid positions for interpolation from retrieval grids back
          // to atmospheric grids
          ArrayOfGridPos gp_p, gp_lat, gp_lon;
          Index          n_p, n_lat, n_lon;
          get_gp_rq_to_atmgrids( gp_p, gp_lat, gp_lon, n_p, n_lat, n_lon,
                                 jq[q], atmosphere_dim, p_grid, lat_grid, lon_grid );
          //
          if( jq[q].Mode() == "rel"  ||  jq[q].Mode() == "logrel" )
            {
              // Find multiplicate factor for elements in vmr_field
              Tensor3 fac_x( n_p, n_lat, n_lon );
              reshape( fac_x, x[ind] ); 
              // Take exp(x) if logrel
              if( jq[q].Mode() == "logrel" )
                { transform( fac_x, exp, fac_x ); }
              Tensor3 factor( vmr_field.npages(), vmr_field.nrows(),
                                                  vmr_field.ncols() );
              regrid_atmfield_by_gp( factor, atmosphere_dim, fac_x,
                                     gp_p, gp_lat, gp_lon );
              for( Index i3=0; i3<vmr_field.ncols(); i3++ )
                { for( Index i2=0; i2<vmr_field.nrows(); i2++ )
                    { for( Index i1=0; i1<vmr_field.npages(); i1++ )
                        { 
                            vmr_field(isp,i1,i2,i3) *= factor(i1,i2,i3); 
                }   }   }
            }
          else if( jq[q].Mode() == "vmr" )
            {
              // Here we just need to map back state x
              Tensor3 vmr_x( n_p, n_lat, n_lon );
              reshape( vmr_x, x[ind] ); 
              Tensor3 vmr( vmr_field.npages(), vmr_field.nrows(),
                                               vmr_field.ncols() );
              regrid_atmfield_by_gp( vmr, atmosphere_dim, vmr_x,
                                     gp_p, gp_lat, gp_lon );
              vmr_field(isp,joker,joker,joker) = vmr;
            }
          else if( jq[q].Mode() == "nd" )
            {
              Tensor3 nd_x( n_p, n_lat, n_lon );
              reshape( nd_x, x[ind] ); 
              Tensor3 nd( vmr_field.npages(), vmr_field.nrows(),
                                              vmr_field.ncols() );
              regrid_atmfield_by_gp( nd, atmosphere_dim, nd_x,
                                     gp_p, gp_lat, gp_lon );
              // Calculate vmr for species (=nd/nd_tot)
              for( Index i3=0; i3<vmr_field.ncols(); i3++ )
                { for( Index i2=0; i2<vmr_field.nrows(); i2++ )
                    { for( Index i1=0; i1<vmr_field.npages(); i1++ )
                        { 
                          vmr_field(isp,i1,i2,i3) = nd(i1,i2,i3) /
                            number_density( p_grid[i1], t_field(i1,i2,i3) );
                }   }   }
            }
          else
            { assert(0); }
        }

      
      // Scattering species
      // ----------------------------------------------------------------------------
      else if( jq[q].MainTag() == SCATSPECIES_MAINTAG )
        {
          // If no cloudbox, we assume that there is nothing to do
          if( cloudbox_on )
            {
              if( particle_bulkprop_field.empty() )
                {
                  throw runtime_error( "One jacobian quantity belongs to the "
                    "scattering species category, but *particle_bulkprop_field* "
                    "is empty." );
                }
              if( particle_bulkprop_field.nbooks() != particle_bulkprop_names.nelem() )
                {
                  throw runtime_error( "Mismatch in size between "
                    "*particle_bulkprop_field* and *particle_bulkprop_field*." );
                }
              
              const Index isp = find_first( particle_bulkprop_names,
                                        jq[q].SubSubtag() );
              if( isp < 0 )
                {
                  ostringstream os;
                  os << "Jacobian quantity with index " << q << " covers a "
                     << "scattering species, and the field quantity is set to \""
                     << jq[q].SubSubtag() << "\", but this quantity "
                     << "could not found in *particle_bulkprop_names*.";
                  throw runtime_error(os.str());
                }
          
              // Determine grid positions for interpolation from retrieval grids back
              // to atmospheric grids
              ArrayOfGridPos gp_p, gp_lat, gp_lon;
              Index          n_p, n_lat, n_lon;
              get_gp_rq_to_atmgrids( gp_p, gp_lat, gp_lon, n_p, n_lat, n_lon,
                                     jq[q], atmosphere_dim, p_grid, lat_grid, lon_grid );
              // Map x to particle_bulkprop_field
              Tensor3 pbfield_x( n_p, n_lat, n_lon );
              reshape( pbfield_x, x[ind] ); 
              Tensor3 pbfield( particle_bulkprop_field.npages(),
                               particle_bulkprop_field.nrows(),
                               particle_bulkprop_field.ncols() );
              regrid_atmfield_by_gp( pbfield, atmosphere_dim, pbfield_x,
                                     gp_p, gp_lat, gp_lon );
              particle_bulkprop_field(isp,joker,joker,joker) = pbfield;
            }
        }

      
      // Pointing off-set
      // ----------------------------------------------------------------------------
      else if( jq[q].MainTag() == POINTING_MAINTAG )
        {
          if( jq[q].Subtag() != POINTING_SUBTAG_A )
            {
              ostringstream os;
              os << "Only pointing off-sets treated by *jacobianAddPointingZa* "
                 << "are so far handled.";
              throw runtime_error(os.str());
            }
          // Handle pointing "jitter" seperately
          if( jq[q].Grids()[0][0] == -1 ) 
            {                          
              if( sensor_los.nrows() != np )
                throw runtime_error( 
                     "Mismatch between pointing jacobian and *sensor_los* found." );
              // Simply add retrieved off-set(s) to za column of *sensor_los*
              for( Index i=0; i<np; i++ )
                { sensor_los(i,0) += x[ji[q][0]+i]; }
            }                                
          // Polynomial representation
          else
            {
              if( sensor_los.nrows() != sensor_time.nelem() )
                throw runtime_error( 
                     "Sizes of *sensor_los* and *sensor_time* do not match." );
              Vector w;
              for( Index c=0; c<np; c++ )
                {
                  polynomial_basis_func( w, sensor_time, c );
                  for( Index i=0; i<w.nelem(); i++ )
                    {  sensor_los(i,0) += w[i] * x[ji[q][0]+c]; }
                }
            }
        }

      
      // Baseline fit: polynomial or sinusoidal
      // ----------------------------------------------------------------------------
      else if( jq[q].MainTag() == POLYFIT_MAINTAG  ||  
               jq[q].MainTag() == SINEFIT_MAINTAG )
        {
          // As the baseline is set using *jacobian*, there must exist a
          // calculated Jacobian, but this is not these case for the first
          // iteration. This should in general be no problem as the a priori
          // for baseline variables should throughout be zero. But this must
          // anyhow be checked:
          if( jacobian.empty() )
            {
              if( min(xa[ind]) != 0  ||  max(xa[ind]) != 0 )
                throw runtime_error(
                   "If any value in *x* that matches a baseline variable "
                   "deviates from zero, *jacobian* must be set." );
            }
          else
            {
              if( jacobian.ncols() != ji[nq-1][1]+1 ) 
                throw runtime_error( "Number of columns in *jacobian* is "
                                     "inconsistent with *jacobian_indices*.");
              if( yb_set )
                {
                    Vector bl( y_baseline.nelem() );
                    mult( bl, jacobian(joker,ind), x[ind] );
                    y_baseline += bl;
                }
                else
                {
                    yb_set = true;
                    y_baseline.resize( jacobian.nrows() );
                    mult(y_baseline, jacobian(joker,ind), x[ind]);
                }
            }
        }

      
      // Or we have to throw an error
      // ----------------------------------------------------------------------------
      else
        {
          ostringstream os;
          os << "Found a retrieval quantity that is not yet handled by\n"
             << "ARTS internal retrievals: " << jq[q].MainTag() << endl;
          throw runtime_error(os.str());
        }
    }

  
  // *y_baseline* not yet set?
  if( !yb_set )
    {
      y_baseline.resize(1);
      y_baseline[0] = 0;
    }
}



/*===========================================================================
  === OEM itself (with wrappers and tempate definitions)
  ===========================================================================*/

// Include only if compiling with C++11.
#ifdef OEM_SUPPORT

// isnan macro breaks invlib.
#ifdef isnan
#undef isnan
#define ISNAN_SET
#endif

#include "oem.h"

#ifdef ISNAN_SET
#define isnan std::isnan
#endif

#include "agenda_wrapper.h"

//
// Check input OEM input arguments.
//
void OEM_checks(
          Workspace&                  ws,
          Vector&                     x,
          Vector&                     yf,
          Matrix&                     jacobian,
    const Agenda&                     inversion_iterate_agenda,
    const Vector&                     xa,
    const CovarianceMatrix&           covmat_sx,
    const CovarianceMatrix&           covmat_se,
    const Index&                      jacobian_do,
    const ArrayOfRetrievalQuantity&   jacobian_quantities,
    const ArrayOfArrayOfIndex&        jacobian_indices,
    const String&                     method,
    const Vector&                     x_norm,
    const Index&                      max_iter,
    const Numeric&                    stop_dx,
    const Vector&                     ml_ga_settings,
    const Index&                      clear_matrices,
    const Index&                      display_progress)
{
  const Index nq = jacobian_quantities.nelem();
  const Index n  = xa.nelem();
  const Index m  = covmat_se.nrows();

  if((x.nelem() != n) && (x.nelem() !=0))
    throw runtime_error( "The length of *x* must be either the same as *xa* or 0." );
  if( covmat_sx.ncols() != covmat_sx.nrows() )
    throw runtime_error( "*covmat_sx* must be a square matrix." );
  if( covmat_sx.ncols() != n )
    throw runtime_error( "Inconsistency in size between *x* and *covmat_sx*." );
  if((yf.nelem() != m) && (yf.nelem() != 0))
    throw runtime_error( "The length of *yf* must be either the same as *y* or 0." );
  if( covmat_se.ncols() != covmat_se.nrows() )
    throw runtime_error( "*covmat_se* must be a square matrix." );
  if( covmat_se.ncols() != m )
    throw runtime_error( "Inconsistency in size between *y* and *covmat_se*." );
  if( !jacobian_do )
    throw runtime_error( "Jacobian calculations must be turned on (but jacobian_do=0)." );
  if((jacobian.nrows() != m) && (!jacobian.empty()))
    throw runtime_error( "The number of rows of the jacobian must be either the number of elements in *y* or 0." );
  if((jacobian.ncols() != n) && (!jacobian.empty()))
      throw runtime_error( "The number of cols of the jacobian must be either the number of elements in *xa* or 0." );
  if( jacobian_indices.nelem() != nq )
    throw runtime_error( "Different number of elements in *jacobian_quantities* "
                          "and *jacobian_indices*." );
  if( nq  &&  jacobian_indices[nq-1][1]+1 != n )
    throw runtime_error( "Size of *covmat_sx* do not agree with Jacobian " 
                          "information (*jacobian_indices*)." );


  // Check GINs
  if( !( method == "li"    || method == "gn"    ||
         method == "ml"    || method == "lm"    ||
         method == "li_cg" || method == "gn_cg" ||
         method == "lm_cg" || method == "ml_cg" ) )
  {
    throw runtime_error( "Valid options for *method* are \"nl\", \"gn\" and " 
                         "\"ml\" or \"lm\"." );
  }

  if( !( x_norm.nelem() == 0  ||  x_norm.nelem() == n ) )
  {
    throw runtime_error( "The vector *x_norm* must have length 0 or match "
                         "*covmat_sx*." );
  }

  if( x_norm.nelem() > 0  &&  min( x_norm ) <= 0 )
  {
    throw runtime_error( "All values in *x_norm* must be > 0." );
  }

  if( max_iter <= 0 )
  {
    throw runtime_error( "The argument *max_iter* must be > 0." );
  }

  if( stop_dx <= 0 )
  {
    throw runtime_error( "The argument *stop_dx* must be > 0." );
  }

  if( method == "ml" )
  {
      if( ml_ga_settings.nelem() != 6 )
      {
          throw runtime_error( "When using \"ml\", *ml_ga_setings* must be a "
                             "vector of length 6." );
      }
      if( min(ml_ga_settings) < 0 )
      {
          throw runtime_error( "The vector *ml_ga_setings* can not contain any "
                               "negative value." );
      }
  }

  if( clear_matrices < 0  ||  clear_matrices > 1 )
    throw runtime_error( "Valid options for *clear_matrices* are 0 and 1." );
  if( display_progress < 0  ||  display_progress > 1 )
    throw runtime_error( "Valid options for *display_progress* are 0 and 1." );

  // If necessary compute yf and jacobian.
  if(x.nelem() == 0) {
    x = xa;
    inversion_iterate_agendaExecute( ws, yf, jacobian, xa, 1,
                                     inversion_iterate_agenda );
  }
  if((yf.nelem() == 0) || (jacobian.empty())) {
    inversion_iterate_agendaExecute( ws, yf, jacobian, x, 1,
                                     inversion_iterate_agenda );
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void OEM(
         Workspace&                  ws,
         Vector&                     x,
         Vector&                     yf,
         Matrix&                     jacobian,
         Matrix&                     dxdy,
         Vector&                     oem_diagnostics,
         Vector&                     ml_ga_history,
         ArrayOfString&              errors,
   const Vector&                     xa,
   const CovarianceMatrix&           covmat_sx,
   const Vector&                     y,
   const CovarianceMatrix&           covmat_se,
   const Index&                      jacobian_do,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const Agenda&                     inversion_iterate_agenda,
   const String&                     method,
   const Numeric&                    max_start_cost,
   const Vector&                     x_norm,
   const Index&                      max_iter,
   const Numeric&                    stop_dx,
   const Vector&                     ml_ga_settings,
   const Index&                      clear_matrices,
   const Index&                      display_progress,
   const Verbosity& )
{
    // Main sizes
    const Index n = covmat_sx.nrows();
    const Index m = y.nelem();

    // Checks
    covmat_sx.compute_inverse();
    covmat_se.compute_inverse();

    OEM_checks(ws, x, yf, jacobian, inversion_iterate_agenda, xa, covmat_sx,
               covmat_se, jacobian_do, jacobian_quantities, jacobian_indices,
               method, x_norm, max_iter, stop_dx, ml_ga_settings,
               clear_matrices, display_progress);

    // Size diagnostic output and init with NaNs
    oem_diagnostics.resize( 5 );
    oem_diagnostics = NAN;
    //
    if( method == "ml" || method == "lm"
        || method == "ml_cg" || method == "lm_cg" )
    {
        ml_ga_history.resize(max_iter + 1);
        ml_ga_history = NAN;
    }
    else
    {
        ml_ga_history.resize(0);
    }

    // Start value of cost function
    Numeric cost_start = NAN;
    if( method == "ml" || method == "lm" || display_progress || max_start_cost > 0 )
    {
        Vector dy  = y; dy -= yf;
        Vector sdy = y; mult(sdy, covmat_se, dy);
        Vector dx  = x; dx -= xa;
        Vector sdx = x; mult(sdx, covmat_sx, dx);
        cost_start = dx * sdx + dy * sdy;
    }
    oem_diagnostics[1] = cost_start;

    // Check for start vector and precomputed yf, jacobian
    if (x.nelem() != n) {
        yf.resize(0);
        jacobian.resize(0,0);
    }

    // Handle cases with too large start cost
    if( max_start_cost > 0  &&  cost_start > max_start_cost )  
    {
        // Flag no inversion in oem_diagnostics, and let x to be undefined 
        oem_diagnostics[0] = 99;
        //
        if( display_progress )
        {
            cout << "\n   No OEM inversion, too high start cost:\n" 
                 << "        Set limit : " << max_start_cost << endl
                 << "      Found value : " << cost_start << endl << endl;
        }
    }
    // Otherwise do inversion
    else {
      bool apply_norm = false;
      OEMMatrix T{};
      if (x_norm.nelem() == n) {
          T.resize(n,n);
          T *= 0.0;
          T.diagonal() = x_norm;
          for (Index i = 0; i < n; i++) {
              T(i,i) = x_norm[i];
          }
          apply_norm = true;
      }

        OEMCovarianceMatrix SeInv(covmat_se), SaInv(covmat_sx);
        OEMVector xa_oem(xa), y_oem(y), x_oem(x);
        AgendaWrapper aw(&ws, (unsigned int) m, (unsigned int) n,
                         jacobian, yf, &inversion_iterate_agenda);
        OEM_STANDARD<AgendaWrapper> oem(aw, xa_oem, SaInv, SeInv);
        int oem_verbosity = static_cast<int>(display_progress);

        int return_code = 0;

        try
        {
            if (method == "li")
            {
                Normed<> s(T, apply_norm);
                GN gn(stop_dx, 1, s); // Linear case, only one step.
                return_code = oem.compute<GN, ArtsLog>(
                    x_oem, y_oem, gn, oem_verbosity,
                    ml_ga_history, true);
            }
            else if (method == "li_cg")
            {
                Normed<CG> cg(T, apply_norm, 1e-12, oem_verbosity);
                GN_CG gn(stop_dx, 1, cg); // Linear case, only one step.
                return_code = oem.compute<GN_CG, ArtsLog>(
                    x_oem, y_oem, gn, oem_verbosity,
                    ml_ga_history, true);
            }
            else if (method == "gn")
            {
                Normed<> s(T, apply_norm);
                GN gn(stop_dx, (unsigned int) max_iter, s); // Linear case, only one step.
                return_code = oem.compute<GN, ArtsLog>(
                    x_oem, y_oem, gn, oem_verbosity,
                    ml_ga_history);
            }
            else if (method == "gn_cg")
            {
                Normed<CG> cg(T, apply_norm, 1e-12, oem_verbosity);
                GN_CG gn(stop_dx, (unsigned int) max_iter, cg);
                return_code = oem.compute<GN_CG, ArtsLog>(
                    x_oem, y_oem, gn, oem_verbosity,
                    ml_ga_history);
            }
            else if ( (method == "lm") || (method == "ml") )
            {
                Normed<> s(T, apply_norm);
                LM_S lm(SaInv, s);

                lm.set_tolerance(stop_dx);
                lm.set_maximum_iterations((unsigned int) max_iter);
                lm.set_lambda(ml_ga_settings[0]);
                lm.set_lambda_decrease(ml_ga_settings[1]);
                lm.set_lambda_increase(ml_ga_settings[2]);
                lm.set_lambda_maximum(ml_ga_settings[3]);
                lm.set_lambda_threshold(ml_ga_settings[4]);
                lm.set_lambda_constraint(ml_ga_settings[5]);

                return_code = oem.compute<LM_S, ArtsLog>(
                    x_oem, y_oem, lm, oem_verbosity,
                    ml_ga_history);
            }
            else if ( (method == "lm_cg") || (method == "ml_cg") )
            {
                Normed<CG> cg(T, apply_norm, 1e-12, oem_verbosity);
                LM_CG_S lm(SaInv, cg);

                lm.set_maximum_iterations((unsigned int) max_iter);
                lm.set_lambda(ml_ga_settings[0]);
                lm.set_lambda_decrease(ml_ga_settings[1]);
                lm.set_lambda_increase(ml_ga_settings[2]);
                lm.set_lambda_threshold(ml_ga_settings[3]);
                lm.set_lambda_maximum(ml_ga_settings[4]);

                return_code = oem.compute<LM_CG_S, ArtsLog>(
                    x_oem, y_oem, lm, oem_verbosity,
                    ml_ga_history);
            }

            oem_diagnostics[0] = static_cast<Index>(return_code);
            oem_diagnostics[2] = oem.cost / static_cast<Numeric>(m);
            oem_diagnostics[3] = oem.cost_y / static_cast<Numeric>(m);
            oem_diagnostics[4] = static_cast<Numeric>(oem.iterations);
        }
        catch (const std::exception & e)
        {
            oem_diagnostics[0]  = 99;
            oem_diagnostics[2] = oem.cost;
            oem_diagnostics[3] = oem.cost_y;
            oem_diagnostics[4] = static_cast<Numeric>(oem.iterations);
            x_oem *= NAN;
            std::vector<std::string> sv = handle_nested_exception(e);
            for (auto & s : sv)
            {

                std::stringstream ss{s};
                std::string t{};
                while (std::getline(ss, t))
                {
                    errors.push_back(t.c_str());
                }
            }
        }
        catch(...)
        {
            throw;
        }

        x  = x_oem;
        yf = aw.yi;

        // Shall empty jacobian and dxdy be returned?
        if(clear_matrices)
        {
            jacobian.resize(0,0);
            dxdy.resize(0,0);
        }
        else if (oem_diagnostics[0] == 0)
        {
            dxdy.resize(n, m);
            Matrix tmp1(n,m), tmp2(n,n), tmp3(n,n);
            mult_inv(tmp1, transpose(jacobian), covmat_se);
            mult(tmp2, tmp1, jacobian);
            add_inv(tmp2, covmat_sx);
            inv(tmp3, tmp2);
            mult(dxdy, tmp3, tmp1);
        }
    }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void covmat_soCalc(
    Matrix& covmat_so,
    const Matrix& dxdy,
    const CovarianceMatrix& covmat_se,
    const Verbosity& /*v*/)
{
    Index n(dxdy.nrows()), m(dxdy.ncols());
    Matrix tmp1(m,n);

    if ((m == 0) || (n == 0)) {
        throw runtime_error("The gain matrix *dxdy* is required to compute the observation error covariance matrix.");
    }
    if ((covmat_se.nrows() != m) || (covmat_se.ncols() != m)) {
        throw runtime_error("The covariance matrix covmat_se has invalid dimensions.");
    }

    covmat_so.resize(n,n);
    mult(tmp1, covmat_se, transpose(dxdy));
    mult(covmat_so, dxdy, tmp1);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void covmat_ssCalc(
    Matrix& covmat_ss,
    const Matrix& avk,
    const CovarianceMatrix& covmat_sx,
    const Verbosity& /*v*/)
{
    Index n(avk.ncols());
    Matrix tmp1(n,n), tmp2(n,n);

    if (n == 0) {
      throw runtime_error("The averaging kernel matrix *dxdy* is required to compute the smoothing error covariance matrix.");
    }
    if ((covmat_sx.nrows() != n) || (covmat_sx.ncols() != n)) {
      throw runtime_error("The covariance matrix *covmat_sx* invalid dimensions.");
    }

    covmat_ss.resize(n,n);

    // Sign doesn't matter since we're dealing with a quadratic form.
    id_mat(tmp1);
    tmp1 -= avk;

    mult(tmp2, covmat_sx, tmp1);
    mult(covmat_ss, tmp1, tmp2);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void avkCalc(
    Matrix& avk,
    const Matrix& dxdy,
    const Matrix& jacobian,
    const Verbosity& /*v*/)
{
    Index m(jacobian.nrows()), n(jacobian.ncols());
    if ((m == 0) || (n == 0)) {
        throw runtime_error("A Jacobian is required to compute the averaging kernel matrix.");
    }
    if ((dxdy.nrows() != n) || (dxdy.ncols() != m)) {
        std::cout << jacobian.nrows() << " | " << jacobian.ncols() << std::endl;
        std::cout << dxdy.nrows() << " | " << dxdy.ncols() << std::endl;
        throw runtime_error("The gain matrix has invalid dimensions.");
    }

    avk.resize(n,n);
    mult(avk, dxdy, jacobian);
}

#else

void covmat_soCalc(
          Matrix& /* covmat_so */,
    const Matrix& /* dxdy */,
    const CovarianceMatrix& /* covmat_se*/,
    const Verbosity& /*v*/)
{
  throw runtime_error("WSM is not available because ARTS was compiled without "
                      "OEM support.");
}

void covmat_ssCalc(
          Matrix& /*covmat_ss*/,
    const Matrix& /*avk*/,
    const CovarianceMatrix& /*covmat_sx*/,
    const Verbosity& /*v*/)
{
  throw runtime_error("WSM is not available because ARTS was compiled without "
                          "OEM support.");
}

void avkCalc(
    Matrix& /* avk */,
    const Matrix& /* dxdy */,
    const Matrix& /* jacobian */,
    const Verbosity& /*v*/)
{
  throw runtime_error("WSM is not available because ARTS was compiled without "
                      "OEM support.");
}


void OEM(Workspace&,
         Vector&,
         Vector&,
         Matrix&,
         Matrix&,
         Vector&,
         Vector&,
         ArrayOfString&,
         const Vector&,
         const CovarianceMatrix&,
         const Vector&,
         const CovarianceMatrix&,
         const Index&,
         const ArrayOfRetrievalQuantity&,
         const ArrayOfArrayOfIndex&,
         const Agenda&,
         const String&,
         const Numeric&,
         const Vector&,
         const Index&,
         const Numeric&,
         const Vector&,
         const Index&,
         const Index&,
         const Verbosity&)
{
  throw runtime_error("WSM is not available because ARTS was compiled without "
                      "OEM support.");
}

#endif // OEM_SUPPORT

#if defined(OEM_SUPPORT) && defined (ENABLE_MPI)

#include "oem_mpi.h"
#include "agenda_wrapper_mpi.h"

//
// Performs manipulations of workspace variables necessary for distributed
// retrievals with MPI:
//
//   - Splits up sensor positions evenly over processes
//   - Splits up inverse covariance matrices.
//
void MPI_Initialize(Matrix&    sensor_los,
                    Matrix&    sensor_pos,
                    Vector&    sensor_time)
{
    int initialized;

    MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(nullptr, nullptr);
    }

    int rank, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int nmblock = (int) sensor_pos.nrows();
    int mblock_range = nmblock / nprocs;
    int mblock_start = mblock_range * rank;
    int remainder = nmblock % std::max(mblock_range, nprocs);

    //
    // Split up sensor positions.
    //

    if (rank < remainder)
    {
        mblock_range += 1;
        mblock_start += rank;
    }
    else
    {
        mblock_start += remainder;
    }

    if (nmblock > 0)
    {
        Range range = Range(mblock_start, mblock_range);

        Matrix tmp_m = sensor_los(range, joker);
        sensor_los = tmp_m;

        tmp_m = sensor_pos(range, joker);
        sensor_pos = tmp_m;

        Vector tmp_v = sensor_time[range];
        sensor_time = tmp_v;
    }
    else
    {
        sensor_los.resize(0,0);
        sensor_pos.resize(0,0);
        sensor_time.resize(0);
    }

}

void OEM_MPI(
    Workspace&                  ws,
    Vector&                     x,
    Vector&                     yf,
    Matrix&                     jacobian,
    Matrix&                     dxdy,
    Vector&                     oem_diagnostics,
    Vector&                     ml_ga_history,
    Matrix&                     sensor_los,
    Matrix&                     sensor_pos,
    Vector&                     sensor_time,
    CovarianceMatrix&           covmat_sx,
    CovarianceMatrix&           covmat_se,
    const Vector&                     xa,
    const Vector&                     y,
    const Index&                      jacobian_do,
    const ArrayOfRetrievalQuantity&   jacobian_quantities,
    const ArrayOfArrayOfIndex&        jacobian_indices,
    const Agenda&                     inversion_iterate_agenda,
    const String&                     method,
    const Numeric&                    max_start_cost,
    const Vector&                     x_norm,
    const Index&                      max_iter,
    const Numeric&                    stop_dx,
    const Vector&                     ml_ga_settings,
    const Index&                      clear_matrices,
    const Index&                      display_progress,
    const Verbosity&                  /*v*/ )
{
    // Main sizes
    const Index n = covmat_sx.nrows();
    const Index m = y.nelem();

    // Check WSVs
    OEM_checks(ws, x, yf, jacobian, inversion_iterate_agenda, xa, covmat_sx,
               covmat_se, jacobian_do, jacobian_quantities, jacobian_indices,
               method, x_norm, max_iter, stop_dx, ml_ga_settings,
               clear_matrices, display_progress);

    // Calculate spectrum and Jacobian for a priori state
    // Jacobian is also input to the agenda, and to flag this is this first
    // call, this WSV must be set to be empty.
    jacobian.resize(0,0);

    // Initialize MPI environment.
    MPI_Initialize(sensor_los, sensor_pos, sensor_time);

    // Setup distributed matrices.
    MPICovarianceMatrix SeInvMPI(covmat_se);
    MPICovarianceMatrix SaInvMPI(covmat_sx);

    // Create temporary MPI vector from local results and use conversion to
    // standard vector to broadcast results to all processes.
    OEMVector tmp;
    inversion_iterate_agendaExecute( ws, tmp, jacobian, xa, 1,
                                     inversion_iterate_agenda );
    yf = MPIVector(tmp);

    // Size diagnostic output and init with NaNs
    oem_diagnostics.resize( 5 );
    oem_diagnostics = NAN;
    //
    if( method == "ml"  || method == "lm" )
    {
        ml_ga_history.resize( max_iter );
        ml_ga_history = NAN;
    }
    else
    {
        ml_ga_history.resize( 0 );
    }

    // Start value of cost function. Covariance matrices are already distributed
    // over processes, so we need to use invlib matrix algebra.
    Numeric cost_start = NAN;
    if( method == "ml" || method == "lm" || display_progress ||  
        max_start_cost > 0 )
    {
        OEMVector dy = y;
        dy -= yf;
        cost_start = dot(dy, SeInvMPI * dy);
    }
    oem_diagnostics[1] = cost_start;


    // Handle cases with too large start cost
    if( max_start_cost > 0  &&  cost_start > max_start_cost )  
    {
        // Flag no inversion in oem_diagnostics, and let x to be undefined 
        oem_diagnostics[0] = 99;
        //
        if( display_progress )
        {
            cout << "\n   No OEM inversion, too high start cost:\n" 
                 << "        Set limit : " << max_start_cost << endl
                 << "      Found value : " << cost_start << endl << endl;
        }
    }

    // Otherwise do inversion
    else
    {
        // Size remaining output arguments
        x.resize( n );
        dxdy.resize( n, m );

        OEMVector        xa_oem(xa), y_oem(y), x_oem;
        AgendaWrapperMPI aw(&ws, &inversion_iterate_agenda, m, n);

        OEM_PS_PS_MPI<AgendaWrapperMPI> oem(aw, xa_oem, SaInvMPI, SeInvMPI);

        // Call selected method
        int return_value = 99;

        if (method == "li")
        {
            CG cg(1e-12, 0);
            GN_CG gn(stop_dx, (unsigned int) max_iter, cg);
            return_value =  oem.compute<GN_CG, invlib::MPILog>(
                x_oem, y_oem, gn,
                2 * (int) display_progress);
        }
        else if (method == "gn")
        {
            CG cg(1e-12, 0);
            GN_CG gn(stop_dx, (unsigned int) max_iter, cg);
            return_value =  oem.compute<GN_CG, invlib::MPILog>(
                x_oem, y_oem, gn,
                2 * (int) display_progress);
        }
        else if ( (method == "lm") || (method == "ml") )
        {
            CG cg(1e-12, 0);
            LM_CG_S_MPI lm(SaInvMPI, cg);

            lm.set_tolerance(stop_dx);
            lm.set_maximum_iterations((unsigned int) max_iter);
            lm.set_lambda(ml_ga_settings[0]);
            lm.set_lambda_decrease(ml_ga_settings[1]);
            lm.set_lambda_increase(ml_ga_settings[2]);
            lm.set_lambda_threshold(ml_ga_settings[3]);
            lm.set_lambda_maximum(ml_ga_settings[4]);

            return_value = oem.compute<LM_CG_S_MPI, invlib::MPILog>(
                x_oem, y_oem, lm,
                2 * (int) display_progress);
        }

        oem_diagnostics[0] = return_value;
        oem_diagnostics[2] = oem.cost;
        oem_diagnostics[3] = oem.cost_y;
        oem_diagnostics[4] = static_cast<Numeric>(oem.iterations);

        x = x_oem;
        // Shall empty jacobian and dxdy be returned?
        if( clear_matrices && (oem_diagnostics[0]))
        {
            jacobian.resize(0,0);
            dxdy.resize(0,0);
        }
    }
    MPI_Finalize();
}

#else

void OEM_MPI(
    Workspace&,
    Vector&,
    Vector&,
    Matrix&,
    Matrix&,
    Vector&,
    Vector&,
    Matrix&,
    Matrix&,
    Vector&,
    CovarianceMatrix&,
    CovarianceMatrix&,
    const Vector&,
    const Vector&,
    const Index&,
    const ArrayOfRetrievalQuantity&,
    const ArrayOfArrayOfIndex&,
    const Agenda&,
    const String&,
    const Numeric&,
    const Vector&,
    const Index&,
    const Numeric&,
    const Vector&,
    const Index&,
    const Index&,
    const Verbosity&)
{
    throw runtime_error("You have to compile ARTS with OEM support "
                        " and enable MPI to use OEM_MPI.");
}

#endif // OEM_SUPPORT && ENABLE_MPI
