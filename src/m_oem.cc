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
#include "oem.h"

extern const String ABSSPECIES_MAINTAG;
extern const String TEMPERATURE_MAINTAG;


/*===========================================================================
  === Help functions 
  ===========================================================================*/

//! Wrapper class for forward model.
/*!
  Wrapper class for the inversion_iterate_agendaExecute function to implement
  the forward model interface used by the non-linear oem function in oem.cc.
  The object is constructed with the pointers to the variables used as arguments
  for the function and then simply forwards the calls made to
  ForwardModel::evaluate() and ForwardModel::evaluate_jacobian() to
  inversion_iterate_agendaExecute.

 */
class AgendaWrapper : public ForwardModel
{
    Workspace *ws;
    Matrix *jacobian;
    const ArrayOfRetrievalQuantity *jacobian_quantities;
    const ArrayOfArrayOfIndex *jacobian_indices;
    const Tensor4 *vmr_field;
    const Tensor3 *t_field;
    const Agenda *inversion_iterate_agenda;
public:
    AgendaWrapper( Workspace *ws_,
                   Matrix *jacobian_,
                   const ArrayOfRetrievalQuantity *jacobian_quantities_,
                   const ArrayOfArrayOfIndex *jacobian_indices_,
                   const Tensor4 *vmr_field_,
                   const Tensor3 *t_field_,
                   const Agenda* inversion_iterate_agenda_) :
        ws( ws_ ), jacobian( jacobian_),
        jacobian_quantities( jacobian_quantities_ ),
        jacobian_indices( jacobian_indices_ ), vmr_field( vmr_field_ ),
        t_field( t_field_ ), inversion_iterate_agenda( inversion_iterate_agenda_ )
        {}

//! Evaluate forward model and compute Jacobian.
/*!

  Forwards the call to evaluate_jacobian() and evaluate() that is made by
  Gauss-Newton and Levenberg-Marquardt OEM methods using the variables pointed
  to by the pointers provided to the constructor as arguments.

  \param[out] y The measurement vector y = K(x) for the current state vector x
  as computed by the forward model.
  \param[out] J The Jacobian Ki=d/dx(K(x)) of the forward model.
  \param[in] x The current state vector x.
*/
    void evaluate_jacobian( VectorView yi,
                            MatrixView Ki,
                            ConstVectorView xi )
        {
            inversion_iterate_agendaExecute( *ws,
                                             dynamic_cast<Vector&>( yi ),
                                             dynamic_cast<Matrix&>( Ki ),
                                             1, *jacobian_quantities,
                                             *jacobian_indices, xi,
                                             *vmr_field, *t_field,
                                             *inversion_iterate_agenda );
        }

//! Evaluate forward model.
/*!

  Forwards the call to evaluate that is made by Gauss-Newton and
  Levenberg-Marquardt OEM methods to the function pointers provided.

  \param[out] y The measurement vector y = K(x) for the current state vector x.
  \param[in] x The current state vector x.
*/
    void evaluate( VectorView yi,
                   ConstVectorView xi )
        {
            inversion_iterate_agendaExecute( *ws,
                                             dynamic_cast<Vector&>( yi ),
                                             *jacobian, 0,
                                             *jacobian_quantities,
                                             *jacobian_indices, xi,
                                             *vmr_field, *t_field,
                                             *inversion_iterate_agenda );
        }
};

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



//! Creates xa for inversion methods.
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

      // Abs species
      if( jq[q].MainTag() == TEMPERATURE_MAINTAG )
        {
          // Here we need to interpolate *vmr_field*
          ArrayOfGridPos gp_p, gp_lat, gp_lon;
          get_gp_atmgrids_to_rq( gp_p, gp_lat, gp_lon, jq[q], atmosphere_dim,
                                 p_grid, lat_grid, lon_grid );
          Tensor3 t_x(gp_p.nelem(),gp_lat.nelem(),gp_lon.nelem());
          regrid_atmfield_by_gp( t_x, atmosphere_dim, t_field(joker,joker,joker), 
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
void x2arts_std(
         Tensor4&                    vmr_field,
         Tensor3&                    t_field,
   const ArrayOfRetrievalQuantity&   jq,
   const ArrayOfArrayOfIndex&        ji,
   const Vector&                     x,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const Verbosity&)
{
  // Main sizes
  const Index nq = jq.nelem();

  // Check input
  if( x.nelem() != ji[nq-1][1]+1 ) 
    throw runtime_error( "Length of *x* does not match information in "
                         "*jacobian_quantities*.");

  // Note that when this method is called, vmr_field and other output variables
  // have original values, i.e. matching the a priori state.

  // Loop retrieval quantities and fill *xa*
  for( Index q=0; q<jq.nelem(); q++ )
    {
      // Index range of this retrieval quantity
      const Index np = ji[q][1] - ji[q][0] + 1;
      Range ind( ji[q][0], np );

      // Abs species
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
          Tensor3 t( t_field.npages(), t_field.nrows(), t_field.ncols() );
          regrid_atmfield_by_gp( t, atmosphere_dim, t_x,
                                 gp_p, gp_lat, gp_lon );
          t_field = t;
        }

      // Abs species
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
          if( jq[q].Mode() == "rel" )
            {
              // Find multiplicate factor for elements in vmr_field
              Tensor3 fac_x( n_p, n_lat, n_lon );
              reshape( fac_x, x[ind] ); 
              Tensor3 factor( vmr_field.npages(), vmr_field.nrows(),
                                                  vmr_field.ncols() );
              regrid_atmfield_by_gp( factor, atmosphere_dim, fac_x,
                                     gp_p, gp_lat, gp_lon );
              for( Index i3=0; i3<=vmr_field.ncols(); i3++ )
                { for( Index i2=0; i2<=vmr_field.nrows(); i2++ )
                    { for( Index i1=0; i1<=vmr_field.npages(); i1++ )
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
              for( Index i3=0; i3<=vmr_field.ncols(); i3++ )
                { for( Index i2=0; i2<=vmr_field.nrows(); i2++ )
                    { for( Index i1=0; i1<=vmr_field.npages(); i1++ )
                        { 
                          vmr_field(isp,i1,i2,i3) = nd(i1,i2,i3) /
                            number_density( p_grid[i1], t_field(i1,i2,i3) );
                }   }   }
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



/* Workspace method: Doxygen documentation will be auto-generated */
void oem(
         Workspace&                  ws,
         Vector&                     x,
         Vector&                     xa,
         Vector&                     yf,
         Matrix&                     jacobian,
         Matrix&                     dxdy,
   const Vector&                     y,
   const Matrix&                     Sx,
   const Matrix&                     So,
   const Index&                      jacobian_do,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const Agenda&                     inversion_iterate_agenda,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const String&                     method,
   const Index&                      yf_linear,
   const Numeric&                    start_ga,
   const Index&                      clear_matrices,
   const Verbosity& )
{
  // Main sizes
  const Index n = Sx.nrows();
  const Index m = y.nelem();
  const Index nq = jacobian_quantities.nelem();

  // Check input
  if( !jacobian_do )
    throw runtime_error( "Jacobian calculations must be turned on (but jacobian_do=0)." );
  if( !nq )
    throw runtime_error( "Jacobian quantities are empty, no inversion to do!." );
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
    throw runtime_error( "Size of *Sx* do not agree with Jacobian information " 
                         "(*jacobian_indices*)." );
  if( !( method == "li"  ||  method == "gn"  ||  method == "ml" ) )  
    throw runtime_error( "Valid options for *method* are \"nl\", \"gn\" and " 
                         "\"ml\"." );
  // Special checks for ML
  if( method == "ml" )
    {
      if( start_ga < 0 )
        throw runtime_error( "*start_ga must be >= 0." );
    }
  //--- End of checks ---------------------------------------------------------------




  // Create xa and init x
  setup_xa( xa, jacobian_quantities, jacobian_indices, atmosphere_dim,
            p_grid, lat_grid, lon_grid, t_field, vmr_field, abs_species );
  //
  x = xa;

  // Calculate spectrum and Jacobian for a priori state
  inversion_iterate_agendaExecute( ws, yf, jacobian, 1, jacobian_quantities,
                                   jacobian_indices, x, vmr_field, t_field, 
                                   inversion_iterate_agenda );

  AgendaWrapper aw( &ws, &jacobian, &jacobian_quantities,
                    &jacobian_indices, &vmr_field, &t_field,
                    &inversion_iterate_agenda );

  if (method == "li")
  {
      oem_linear( x, y, xa, jacobian, So, Sx, false );

      if( yf_linear )
      {
        inversion_iterate_agendaExecute( ws, yf, jacobian, 1, jacobian_quantities,
                                         jacobian_indices, x, vmr_field, t_field, 
                                         inversion_iterate_agenda );
      }
  }

  else if (method == "gn")
  {
      oem_gauss_newton( x, y, xa, aw, So, Sx, 10e-5, 1000 );
  }

  else if (method == "lm")
  {
      Numeric gamma_max = 100.0;
      Numeric gamma_scale_dec = 2.0;
      Numeric gamma_scale_inc = 3.0;
      Numeric gamma_threshold = 1.0;
      oem_levenberg_marquardt( x, y, xa, aw, So, Sx, 10e-5, 1000,
                               start_ga, gamma_scale_dec, gamma_scale_inc,
                               gamma_max, gamma_threshold );
  }

  // So far we just create a dummy dxdy, matching zero measurement response
  dxdy.resize(n,m); 
  dxdy = 0;


  // Shall empty jacobian and dxdy be returned
  if( clear_matrices )
    {
      jacobian.resize(0,0); 
      dxdy.resize(0,0); 
    }
}
