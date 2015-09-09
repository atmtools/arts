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
#include <stdexcept>
#include "arts.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "math_funcs.h"
#include "physics_funcs.h"

extern const String ABSSPECIES_MAINTAG;


/*===========================================================================
  === Help functions 
  ===========================================================================*/


//! Determines grid positions for regridding of atmospheric fields
/*!
  The grid positions arrays are sized inside the function. gp_lat is given
  length 0 for atmosphere_dim=1 etc.

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
void get_gp_for_jq_grids( 
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



//! Sets up xa for inversion methods.
/*!
  The function analyses jacobian_quantities and jacobian_indices to creat xa.

  \param[out] xa                   As the WSV with same name.
  \param[in]  jq                   Matches the WSV *jacobian_quantities*.
  \param[in]  ji                   Matches the WSV *jacobian_indices*.
  \param[in]  atmosphere_dim       As the WSV with same name.
  \param[in]  p_grid               As the WSV with same name.
  \param[in]  lat_grid             As the WSV with same name.
  \param[in]  lon_grid             As the WSV with same name.
  \param[in]  t_field              As the WSV with same name.
  \param[in]  vmr_field            As the WSV with same name.
  \param[in]  abs_species          As the WSV with same name.

  \author Patrick Eriksson 
  \date   2015-09-09
*/
void setup_xa( 
         Vector&                     xa,
   const ArrayOfRetrievalQuantity&   jq,
   const ArrayOfArrayOfIndex&        ji,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const ArrayOfArrayOfSpeciesTag&   abs_species )
{
  const Index nq = jq.nelem();

  xa.resize( ji[nq-1][1]+1 );

  for( Index q=0; q<jq.nelem(); q++ )
    {
      // Index range of this retrieval quantity
      const Index np = ji[q][1] - ji[q][0] + 1;
      Range ind( ji[q][0], np );

      // Index position of species
      ArrayOfSpeciesTag  atag;
      array_species_tag_from_string( atag, jq[q].Subtag() );
      const Index isp = chk_contains( "abs_species", abs_species, atag );

      // Abs species
      if( jq[q].MainTag() == ABSSPECIES_MAINTAG )
        {
          if( jq[q].Mode() == "rel" )
            {
              // This one is simple, just a vector of ones
              xa[ind] = 1; 
            }
          else if( jq[q].Mode() == "vmr" )
            {
              // Here we need to interpolate *vmr_field*
              ArrayOfGridPos gp_p, gp_lat, gp_lon;
              get_gp_for_jq_grids( gp_p, gp_lat, gp_lon, jq[q], atmosphere_dim,
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
              get_gp_for_jq_grids( gp_p, gp_lat, gp_lon, jq[q], atmosphere_dim,
                                   p_grid, lat_grid, lon_grid );
              Tensor3 vmr_x(gp_p.nelem(),gp_lat.nelem(),gp_lon.nelem());
              Tensor3 t_x(gp_p.nelem(),gp_lat.nelem(),gp_lon.nelem());
              regrid_atmfield_by_gp( vmr_x, atmosphere_dim, 
                                     vmr_field(isp,joker,joker,joker), 
                                     gp_p, gp_lat, gp_lon );
              regrid_atmfield_by_gp( t_x, atmosphere_dim,  t_field,
                                     gp_p, gp_lat, gp_lon );
              // Calculate number density
              Index i = 0;
              for( Index i3=0; i3<=vmr_x.ncols(); i3++ )
                {
                  for( Index i2=0; i2<=vmr_x.nrows(); i2++ )
                    {
                      for( Index i1=0; i1<=vmr_x.npages(); i1++ )
                        {
                          xa[ji[q][0]+i] = vmr_x(i1,i2,i3) *
                            number_density( jq[q].Grids()[0][i1], t_x(i1,i2,i3) );
                          i += 1;
                        }
                    }
                }
            }
          else
            { assert(0); }
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



/*===========================================================================
  === Workspace methods 
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void oem(
         Vector&                     x,
         Vector&                     xa,
         Vector&                     yf,
   const Vector&                     y,
   const Matrix&                     Sx,
   const Matrix&                     So,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const String&                     method,
   const Numeric&                    start_ga,
   const Verbosity&)
{
  // Main sizes
  const Index n = Sx.nrows();
  const Index m = y.nelem();
  const Index nq = jacobian_quantities.nelem();

  // Check input
  if( Sx.ncols() != n )
    throw runtime_error( "*covmat_sx* must be a square matrix." );
  if( So.ncols() != So.nrows() )
    throw runtime_error( "*covmat_so* must be a square matrix." );
  if( So.ncols() != m )
    throw runtime_error( "Inconsistency in size between *y* and *covmat_so*." );
  if( jacobian_indices.nelem() != nq )
    throw runtime_error( "Different number of elements in *jacobian_quantities* "
                         "and *jacobian_indices*." );
  if( jacobian_indices[nq-1][1]+1 != n )
    throw runtime_error( "Length of *x* and last value in *jacobian_indices* ." 
                         "do not agree." );
  if( !( method == "li"  ||  method == "gn"  ||  method == "ml" ) )  
    throw runtime_error( "Valid options for *method* are \"nl\", \"gn\" and " 
                         "\"ml\"." );
  // Special checks for ML
  if( method == "ml" )  
    {
      if( start_ga < 0 )
        throw runtime_error( "*start_ga must be >= 0." );
    }

  // Temporary limitations
  if( atmosphere_dim > 1 )
    throw runtime_error( "Only 1D is handled so far." );
    

  // Create xa
  setup_xa( xa, jacobian_quantities, jacobian_indices, atmosphere_dim,
            p_grid, lat_grid, lon_grid, t_field, vmr_field, abs_species );

  // So far dummy x and yf
  x = xa;
  yf.resize(m);
  yf = 0;
}
