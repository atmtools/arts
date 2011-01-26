/* Copyright (C) 2002-2010
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
  \file   m_rte.cc
  \author Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2002-05-11 

  \brief  Workspace functions for solving clear sky radiative transfer.

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
#include "check_input.h"
#include "jacobian.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "montecarlo.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"

extern const String ABSSPECIES_MAINTAG;
extern const String TEMPERATURE_MAINTAG;

// Number of auxilary variables, in y_aux
const Index N_AUX = 3;


// A maqcro to loop analytical jacobian quantities
#define FOR_ANALYTICAL_JACOBIANS_DO(what_to_do) \
  for( Index iq=0; iq<jacobian_quantities.nelem(); iq++ ) \
    { \
      if( jacobian_quantities[iq].Analytical() ) \
        { what_to_do } \
    } 




/*===========================================================================
  === Special sub-functions for analytical jacobians (in alphabetical order)
  === (Main sub-functions get_iy_of_background and iyb_calc placed separately)
  ===========================================================================*/


//! diy_from_path_to_rgrids
/*!
    Maps jacobian data for points along the propagation path, to
    jacobian retrieval grid data.

    \param   diy_dx              Out: Jacobians for selected retrieval grids.
    \param   jacobian_quantity   As the WSV.
    \param   diy_dpath           Jacobians along the propagation path.
    \param   atmosphere_dim      As the WSV.
    \param   ppath               As the WSV.
    \param   ppath_p             The pressure at each ppath point.

    \author Patrick Eriksson 
    \date   2009-10-08
*/
// A small help function, to make the code below cleaner
void from_dpath_to_dx(
         MatrixView   diy_dx,
    ConstMatrixView   diy_dq,
    const Numeric&    w )
{
  for( Index irow=0; irow<diy_dx.nrows(); irow++ )
    { 
      for( Index icol=0; icol<diy_dx.ncols(); icol++ )
        { diy_dx(irow,icol) += w * diy_dq(irow,icol); }
    }
}
//
void diy_from_path_to_rgrids(
         Tensor3View          diy_dx,
   const RetrievalQuantity&   jacobian_quantity,
   ConstTensor3View           diy_dpath,
   const Index&               atmosphere_dim,
   const Ppath&               ppath,
   ConstVectorView            ppath_p )
{
  // We want here an extrapolation to infinity -> 
  //                                        extremly high extrapolation factor
  const Numeric   extpolfac = 1.0e99;

  // Retrieval grid of interest
  Vector r_grid;

  if( ppath.np > 1 )  // Otherwise nothing to do here
    {
      // Pressure
      r_grid = jacobian_quantity.Grids()[0];
      Index            nr1 = r_grid.nelem();
      ArrayOfGridPos   gp_p(ppath.np);
      p2gridpos( gp_p, r_grid, ppath_p, extpolfac );

      // Latitude
      Index            nr2 = 1;
      ArrayOfGridPos   gp_lat;
      if( atmosphere_dim > 1 )
        {
          gp_lat.resize(ppath.np);
          r_grid = jacobian_quantity.Grids()[1];
          nr2    = r_grid.nelem();
          gridpos( gp_lat, r_grid, ppath.pos(joker,1), extpolfac );
        }

      // Longitude
      Index            nr3 = 1;
      ArrayOfGridPos   gp_lon;
      if( atmosphere_dim > 2 )
        {
          gp_lon.resize(ppath.np);
          r_grid = jacobian_quantity.Grids()[2];
          nr3    = r_grid.nelem();
          gridpos( gp_lon, r_grid, ppath.pos(joker,2), extpolfac );
        }
      
      //- 1D
      if( atmosphere_dim == 1 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              if( gp_p[ip].fd[0] < 1 )
                {
                  from_dpath_to_dx( diy_dx(gp_p[ip].idx,joker,joker),
                                    diy_dpath(ip,joker,joker), gp_p[ip].fd[1] );
                }
              if( gp_p[ip].fd[0] > 0 )
                {
                  from_dpath_to_dx( diy_dx(gp_p[ip].idx+1,joker,joker),
                                    diy_dpath(ip,joker,joker), gp_p[ip].fd[0] );
                }
            }
        }

      //- 2D
      else if( atmosphere_dim == 2 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              Index   ix = nr1*gp_lat[ip].idx + gp_p[ip].idx;
              // Low lat, low p
              from_dpath_to_dx( diy_dx(ix,joker,joker),
                                diy_dpath(ip,joker,joker), 
                                gp_lat[ip].fd[1]*gp_p[ip].fd[1] );
              // Low lat, high p
              from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                                gp_lat[ip].fd[1]*gp_p[ip].fd[0] );
              // High lat, low p
              from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                                gp_lat[ip].fd[0]*gp_p[ip].fd[1] );
              // High lat, high p
              from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                                gp_lat[ip].fd[0]*gp_p[ip].fd[0] );
            }
        }

      //- 3D
      else if( atmosphere_dim == 3 )
        {
          for( Index ip=0; ip<ppath.np; ip++ )
            {
              Index   ix = nr2*nr1*gp_lon[ip].idx +
                           nr1*gp_lat[ip].idx + gp_p[ip].idx;
              // Low lon, low lat, low p
              from_dpath_to_dx( diy_dx(ix,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[1]*gp_p[ip].fd[1]);
              // Low lon, low lat, high p
              from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[1]*gp_p[ip].fd[0]);
              // Low lon, high lat, low p
              from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[0]*gp_p[ip].fd[1]);
              // Low lon, high lat, high p
              from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[1]*gp_lat[ip].fd[0]*gp_p[ip].fd[0]);

              // Increase *ix* (to be valid for high lon level)
              ix += nr2*nr1;

              // High lon, low lat, low p
              from_dpath_to_dx( diy_dx(ix,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[1]*gp_p[ip].fd[1]);
              // High lon, low lat, high p
              from_dpath_to_dx( diy_dx(ix+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[1]*gp_p[ip].fd[0]);
              // High lon, high lat, low p
              from_dpath_to_dx( diy_dx(ix+nr1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[0]*gp_p[ip].fd[1]);
              // High lon, high lat, high p
              from_dpath_to_dx( diy_dx(ix+nr1+1,joker,joker),
                                diy_dpath(ip,joker,joker), 
                             gp_lon[ip].fd[0]*gp_lat[ip].fd[0]*gp_p[ip].fd[0]);
            }
        }
    }
}



//! Help function for analytical jacobian calculations
/*!
    The function determines which terms in jacobian_quantities that are 
    analytical absorption species and temperature jacobians. 

    *abs_species_i* and *is_t* shall be sized to have the same length
    as *jacobian_quantities*. For analytical absorption species
    jacobians, *abs_species_i* is set to the matching index in
    *abs_species*. Otherwise, to -1. For analytical temperature
    jacobians, *is_t* is set to 1. Otherwise to 0.

    \param   abs_species_i         Out: Matching index in abs_species 
    \param   is_t                  Out: Flag for: Is a temperature jacobian?
    \param   jacobian_quantities   As the WSV.
    \param   abs_species           As the WSV.


    \author Patrick Eriksson 
    \date   2009-10-07
*/
void get_pointers_for_analytical_jacobians( 
        ArrayOfIndex&               abs_species_i, 
        ArrayOfIndex&               is_t,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfSpeciesTag&   abs_species )
{

  FOR_ANALYTICAL_JACOBIANS_DO( 
    //
    if( jacobian_quantities[iq].MainTag() == TEMPERATURE_MAINTAG )
      { is_t[iq] = 1; }
    else
      { is_t[iq] = 0; }
    //
    if( jacobian_quantities[iq].MainTag() == ABSSPECIES_MAINTAG )
      {
        ArrayOfSpeciesTag  atag;
        array_species_tag_from_string( atag, jacobian_quantities[iq].Subtag() );
        abs_species_i[iq] = chk_contains( "abs_species", abs_species, atag );
      }
    else
      { abs_species_i[iq] = -1; }
  )
}









/*===========================================================================
  === Main help functions (get_iy_of_background and iyb_calc)
  ===========================================================================*/


//! get_iy_of_background
/*!
    Determines iy of the "background" of a propgation path.

    The task is to determine *iy* and related variables for the
    background, or to continue the raditiave calculations
    "backwards". The details here depends on the method selected for
    the agendas.

    Each background is handled by an agenda. Several of these agandes
    can involve recursive calls of *iy_clearsky_agenda*. It is also
    allowed to input *iy_clearsky_basic_agenda* instead of
    *iy_clearsky_agenda*.

    \param   ws                    Out: The workspace
    \param   iy                    Out: As the WSV.
    \param   iy_error              Out: As the WSV.
    \param   iy_error_type         Out: As the WSV.
    \param   iy_aux                Out: As the WSV.
    \param   diy_aux               Out: As the WSV.
    \param   iy_transmission       As the WSV.
    \param   iy_aux_do             As the WSV.
    \param   jacobian_do           As the WSV.
    \param   ppath                 As the WSV.
    \param   atmosphere_dim        As the WSV.
    \param   p_grid                As the WSV.
    \param   lat_grid              As the WSV.
    \param   lon_grid              As the WSV.
    \param   t_field               As the WSV.
    \param   vmr_field             As the WSV.
    \param   cloudbox_on           As the WSV.
    \param   stokes_dim            As the WSV.
    \param   f_grid                As the WSV.
    \param   iy_clearsky_agenda    As the WSV or iy_clearsky_basic_agenda.
    \param   iy_space_agenda       As the WSV.
    \param   surface_prop_agenda   As the WSV.
    \param   iy_cloudbox_agenda    As the WSV.

    \author Patrick Eriksson 
    \date   2009-10-08
*/
void get_iy_of_background(
        Workspace&        ws,
        Matrix&           iy,
        Matrix&           iy_error,
        Index&            iy_error_type,
        Tensor3&          iy_aux,
        ArrayOfTensor3&   diy_dx,
  ConstTensor3View        iy_transmission,
  const Index&            iy_aux_do,
  const Index&            jacobian_do,
  const Ppath&            ppath,
  const Index&            atmosphere_dim,
  ConstVectorView         p_grid,
  ConstVectorView         lat_grid,
  ConstVectorView         lon_grid,
  ConstTensor3View        t_field,
  ConstTensor4View        vmr_field,
  const Matrix&           r_geoid,
  const Matrix&           z_surface,
  const Index&            cloudbox_on,
  const Index&            stokes_dim,
  ConstVectorView         f_grid,
  const Agenda&           iy_clearsky_agenda,
  const Agenda&           iy_space_agenda,
  const Agenda&           surface_prop_agenda,
  const Agenda&           iy_cloudbox_agenda )
{
  // Some sizes
  const Index nf      = f_grid.nelem();
  const Index np      = ppath.np;

  // Set rte_pos, rte_gp_XXX and rte_los to match the last point in ppath.
  // The agendas below use different combinations of these variables.
  //
  // Note that the Ppath positions (ppath.pos) for 1D have one column more
  // than expected by most functions. Only the first atmosphere_dim values
  // shall be copied.
  //
  Vector rte_pos;
  Vector rte_los;
  GridPos rte_gp_p;
  GridPos rte_gp_lat;
  GridPos rte_gp_lon;
  rte_pos.resize( atmosphere_dim );
  rte_pos = ppath.pos(np-1,Range(0,atmosphere_dim));
  rte_los.resize( ppath.los.ncols() );
  rte_los = ppath.los(np-1,joker);
  gridpos_copy( rte_gp_p, ppath.gp_p[np-1] );
  if( atmosphere_dim > 1 )
    { gridpos_copy( rte_gp_lat, ppath.gp_lat[np-1] ); }
  if( atmosphere_dim > 2 )
    { gridpos_copy( rte_gp_lon, ppath.gp_lon[np-1] ); }

  out3 << "Radiative background: " << ppath.background << "\n";


  // Handle the different background cases
  //
  switch( ppath_what_background( ppath ) )
    {

    case 1:   //--- Space ---------------------------------------------------- 
      {
        iy_space_agendaExecute( ws, iy, rte_pos, rte_los, iy_space_agenda );
        
        if( iy.ncols() != stokes_dim  ||  iy.nrows() != nf )
          {
            ostringstream os;
            os << "expected size = [" << stokes_dim << "," << nf << "]\n"
               << "size of iy    = [" << iy.nrows() << "," << iy.ncols()<< "]\n"
               << "The size of *iy* returned from *iy_space_agenda* is "
               << "not correct.";
            throw runtime_error( os.str() );      
          }
      }
      break;

    case 2:   //--- The surface -----------------------------------------------
      {
        // Call *surface_prop_agenda*
        //
        Matrix    surface_los;
        Tensor4   surface_rmatrix;
        Matrix    surface_emission;
        //
        surface_prop_agendaExecute( ws, surface_emission, surface_los, 
                                    surface_rmatrix, rte_pos, rte_los, 
                                    rte_gp_p, rte_gp_lat, rte_gp_lon,
                                    surface_prop_agenda );

        // Check output of *surface_prop_agenda*
        //
        const Index   nlos = surface_los.nrows();
        //
        if( nlos )   // if 0, blackbody ground and some checks are not needed
          {
            if( surface_los.ncols() != rte_los.nelem() )
              throw runtime_error( 
                        "Number of columns in *surface_los* is not correct." );
            if( nlos != surface_rmatrix.nbooks() )
              throw runtime_error( 
                  "Mismatch in size of *surface_los* and *surface_rmatrix*." );
            if( surface_rmatrix.npages() != nf )
              throw runtime_error( 
                       "Mismatch in size of *surface_rmatrix* and *f_grid*." );
            if( surface_rmatrix.nrows() != stokes_dim  ||  
                surface_rmatrix.ncols() != stokes_dim ) 
              throw runtime_error( 
              "Mismatch between size of *surface_rmatrix* and *stokes_dim*." );
          }
        if( surface_emission.ncols() != stokes_dim )
          throw runtime_error( 
             "Mismatch between size of *surface_emission* and *stokes_dim*." );
        if( surface_emission.nrows() != nf )
          throw runtime_error( 
                       "Mismatch in size of *surface_emission* and f_grid*." );
        //---------------------------------------------------------------------

        // Variable to hold down-welling radiation
        Tensor3   I( nlos, nf, stokes_dim );
 
        // Loop *surface_los*-es. If no such LOS, we are ready.
        if( nlos > 0 )
          {
            for( Index ilos=0; ilos<nlos; ilos++ )
              {
                Vector los = surface_los(ilos,joker);

                // Include surface reflection matrix in *iy_transmission*
                // If iy_transmission is empty, this is interpreted as the
                // variable is not needed.
                //
                Tensor3 iy_trans_new;
                //
                if( iy_transmission.npages() )
                  {
                    iy_transmission_mult(  iy_trans_new, iy_transmission, 
                                     surface_rmatrix(ilos,joker,joker,joker) );
                  }

                // Calculate angular tilt of the surface
                // (sign(los[0]) to handle negative za for 2D)
                const Numeric atilt = surfacetilt( atmosphere_dim, lat_grid, 
                                                 lon_grid, r_geoid, z_surface,
                                                 rte_gp_lat, rte_gp_lat, los );
                const Numeric zamax = 90 - sign(los[0])*atilt;

                // I considered to add a check that surface_los is above the
                // horizon, but that would force e.g. surface_specular_los to
                // actually calculate the surface tilt, which causes
                // unnecessary calculation overhead. The angles are now moved
                // to be just above the horizon, which should be acceptable.

                // Make sure that we have some margin to the "horizon"
                // (otherwise numerical problems can create an infinite loop)
                if( atmosphere_dim == 2  && los[0]<0 ) //2D with za<0
                  { los[0] = max( -zamax+0.1, los[0] ); }
                else
                  { los[0] = min( zamax-0.1, los[0] ); }

                // Calculate downwelling radiation for LOS ilos 
                //
                // The variable iy_clearsky_agenda can in fact be 
                // iy_clearsky_BASIC_agenda
                //
                if( iy_clearsky_agenda.name() == "iy_clearsky_basic_agenda" )
                  {
                    iy_clearsky_basic_agendaExecute( ws, iy, rte_pos, los,
                                              cloudbox_on, iy_clearsky_agenda);
                  }
                else
                  {
                    iy_clearsky_agendaExecute( ws, iy, iy_error, iy_error_type,
                            iy_aux, diy_dx, 0, rte_pos, los,
                            iy_trans_new, cloudbox_on, jacobian_do, iy_aux_do, 
                            f_grid, p_grid, lat_grid, lon_grid, t_field, 
                            vmr_field, iy_clearsky_agenda );
                  }

                I(ilos,joker,joker) = iy;
              }
          }

        // Add up
        surface_calc( iy, I, surface_los, surface_rmatrix, surface_emission );
      }
      break;


    case 3:   //--- Cloudbox surface or interior ------------------------------
    case 4:
      {
        // Pass a local copy of ppath to the cloudbox agenda
        Ppath   ppath_local;
        ppath_init_structure( ppath_local, atmosphere_dim, ppath.np );
        ppath_copy( ppath_local, ppath );

        // The cloudbox agenda is not yet handling iy_aux and iy_transmission.
        //
        iy_cloudbox_agendaExecute( ws, iy, ppath_local,
                                   rte_pos, rte_los, rte_gp_p,
                                   rte_gp_lat, rte_gp_lon,
                                   iy_cloudbox_agenda );

        if( iy.nrows() != nf  ||  iy.ncols() != stokes_dim )
          {
            out1 << "expected size = [" << nf << "," << stokes_dim << "]\n";
            out1 << "iy size = [" << iy.nrows() << "," << iy.ncols()<< "]\n";
            throw runtime_error( "The size of *iy* returned from "
                                      "*iy_cloudbox_agenda* is not correct." );
          }
      }
      break;

    default:  //--- ????? ----------------------------------------------------
      // Are we here, the coding is wrong somewhere
      assert( false );
    }
}



//! iyb_calc
/*!
    Calculation of pencil beam monochromatic spectra for 1 measurement block.

    All in- and output variables as the WSV with the same name.

    \author Patrick Eriksson 
    \date   2009-10-16
*/
void iyb_calc(
        Workspace&                  ws,
        Vector&                     iyb,
        Vector&                     iyb_error,
        Index&                      iy_error_type,
        Matrix&                     iyb_aux,
        ArrayOfMatrix&              diyb_dx,
  const Index&                      imblock,
  const Index&                      atmosphere_dim,
  ConstVectorView                   p_grid,
  ConstVectorView                   lat_grid,
  ConstVectorView                   lon_grid,
  ConstTensor3View                  t_field,
  ConstTensor4View                  vmr_field,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  ConstVectorView                   f_grid,
  ConstMatrixView                   sensor_pos,
  ConstMatrixView                   sensor_los,
  ConstVectorView                   mblock_za_grid,
  ConstVectorView                   mblock_aa_grid,
  const Index&                      antenna_dim,
  const Agenda&                     iy_clearsky_agenda,
  const Index&                      iy_aux_do,
  const String&                     y_unit,
  const Index&                      j_analytical_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices )
{
  // Sizes
  const Index   nf      = f_grid.nelem();
  const Index   nza     = mblock_za_grid.nelem();
        Index   naa     = mblock_aa_grid.nelem();   
  if( antenna_dim == 1 )  
    { naa = 1; }
  const Index   niyb    = nf * nza * naa * stokes_dim;

  // Set up size of containers for data of 1 measurement block.
  //
  // We do not know the number of aux variables. The parallelisation makes it
  // unsafe to do the allocation of iyb_aux inside the for-loops. We must use 
  // a hard-coded maximum number.
  //
  iyb.resize( niyb );
  iyb_error.resize( niyb );
  if( iy_aux_do )
    iyb_aux.resize( niyb, N_AUX );
  else
    iyb_aux.resize( 0, 0 );
  //
  if( j_analytical_do )
    {
      diyb_dx.resize( jacobian_indices.nelem() );
      FOR_ANALYTICAL_JACOBIANS_DO(
        diyb_dx[iq].resize( niyb, jacobian_indices[iq][1] -
                                  jacobian_indices[iq][0] + 1 );
      )
    }
  else
    { diyb_dx.resize( 0 ); }

  // Polarisation idex variable
  // If some kind of modified Stokes vector will be introduced, this must be
  // made to a WSV
  ArrayOfIndex i_pol(stokes_dim);
  for( Index is=0; is<stokes_dim; is++ )
    { i_pol[is] = is + 1; }

  // We have to make a local copy of the Workspace and the agendas because
  // only non-reference types can be declared firstprivate in OpenMP
  Workspace l_ws (ws);
  Agenda l_iy_clearsky_agenda (iy_clearsky_agenda);
  
  // Start of actual calculations
/*#pragma omp parallel for                                          \
  if(!arts_omp_in_parallel() && nza>1)                            \
  default(none)                                                   \
  firstprivate(l_ws, l_iy_clearsky_agenda)                        \
  shared(sensor_los, mblock_za_grid, mblock_aa_grid, vmr_field,   \
         t_field, lon_grid, lat_grid, p_grid, f_grid, sensor_pos, \
         joker, naa) */
#pragma omp parallel for                                          \
  if(!arts_omp_in_parallel() && nza>1)                            \
  firstprivate(l_ws, l_iy_clearsky_agenda)
  for( Index iza=0; iza<nza; iza++ )
    {
      // The try block here is necessary to correctly handle
      // exceptions inside the parallel region. 
      try
        {
          for( Index iaa=0; iaa<naa; iaa++ )
            {
              //--- LOS of interest
              //
              Vector los( sensor_los.ncols() );
              //
              los     = sensor_los( imblock, joker );
              los[0] += mblock_za_grid[iza];

              // Handle za/aa_grid "out-of-bounds" and mapping effects
              //
              if( antenna_dim == 2 )
                { map_daa( los[0], los[1], los[0], los[1], 
                                                    mblock_aa_grid[iaa] ); }
              else if( atmosphere_dim == 1  && abs(los[0]-90) > 90 )
                { if( los[0] < 0 )          { los[0] = -los[0]; }
                  else if( los[0] > 180 )   { los[0] = 360 - los[0]; } }
              else if( atmosphere_dim == 2  && abs(los[0]) > 180 )
                { los[0] = -sign(los[0])*360 + los[0]; }
              else if( atmosphere_dim == 3  &&  abs(los[0]-90) > 90 )
                { map_daa( los[0], los[1], los[0], los[1], 0 ); }


              // Calculate iy and associated variables
              //
              Matrix         iy, iy_error;
              Tensor3        iy_aux, iy_transmission;
              ArrayOfTensor3 diy_dx; 
              //
              iy_clearsky_agendaExecute( l_ws, iy, iy_error, iy_error_type, 
                             iy_aux, diy_dx, 1, sensor_pos(imblock,joker), los, 
                             iy_transmission, cloudbox_on, j_analytical_do, 
                             iy_aux_do, f_grid, p_grid, lat_grid, lon_grid, 
                             t_field, vmr_field, l_iy_clearsky_agenda );

              // Start row in iyb etc. for present LOS
              //
              const Index row0 = ( iza*naa + iaa ) * nf * stokes_dim;

              // Jacobian part (must be converted to Tb before iy for PlanckBT)
              // 
              if( j_analytical_do )
                {
                  FOR_ANALYTICAL_JACOBIANS_DO(
                    //
                    apply_y_unit2( diy_dx[iq], iy, y_unit, f_grid, i_pol );
                    //
                    for( Index ip=0; ip<jacobian_indices[iq][1] -
                                        jacobian_indices[iq][0]+1; ip++ )
                      {
                        for( Index is=0; is<stokes_dim; is++ )
                          { 
                            diyb_dx[iq](Range(row0+is,nf,stokes_dim),ip)=
                                                     diy_dx[iq](ip,joker,is); 
                          }
                      }                              
                  )
                }

              // iy       : unit conversion and copy to iyb
              // iy_error : unit conversion and copy to iyb_error
              // iy_aux   : copy to iyb_aux
              //
              if( iy_error_type > 0 )
                { apply_y_unit2( iy_error, iy, y_unit, f_grid, i_pol ); }
              //
              apply_y_unit( iy, y_unit, f_grid, i_pol );
              //
              for( Index is=0; is<stokes_dim; is++ )
                { 
                  iyb[Range(row0+is,nf,stokes_dim)] = iy(joker,is); 
                  if( iy_error_type > 0 )
                    { iyb_error[Range(row0+is,nf,stokes_dim)] = 
                                                            iy_error(joker,is); 
                    }
                  //
                  if( iy_aux_do )
                    {
                      for( Index iaux=0; iaux<N_AUX; iaux++ )
                        { 
                          iyb_aux(Range(row0+is,nf,stokes_dim),iaux) = 
                                                         iy_aux(iaux,joker,is);
                        }
                    }
                }
            }  // End aa loop
        }  // End try

      catch (runtime_error e)
        {
          exit_or_rethrow(e.what());
        }
    }  // End za loop
}










/*===========================================================================
  === The workspace methods (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void iyBeerLambertStandardClearsky(
        Workspace&                  ws,
        Matrix&                     iy,
        Matrix&                     iy_error,
        Index&                      iy_error_type,
        Tensor3&                    iy_aux,
        ArrayOfTensor3&             diy_dx,
  const Index&                      iy_agenda_call1,
  const Vector&                     rte_pos,      
  const Vector&                     rte_los,      
  const Index&                      iy_aux_do,
  const Index&                      jacobian_do,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Tensor3&                    z_field,
  const Tensor3&                    t_field,
  const Tensor4&                    vmr_field,
  const Matrix&                     r_geoid,
  const Matrix&                     z_surface,
  const Index&                      cloudbox_on,
  const ArrayOfIndex&               cloudbox_limits,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const ArrayOfArrayOfSpeciesTag&   abs_species,
  const Agenda&                     ppath_step_agenda,
  const Agenda&                     abs_scalar_agenda,
  const Agenda&                     iy_clearsky_agenda,
  const Agenda&                     iy_space_agenda,
  const Agenda&                     surface_prop_agenda,
  const Agenda&                     iy_cloudbox_agenda,
  const Tensor3&                    iy_transmission,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices )
{
  // A minimal amount of checks, for efficiency reasons
  assert( rte_pos.nelem() == atmosphere_dim );
  assert( ( atmosphere_dim < 3   &&  rte_los.nelem() == 1 )  ||
          ( atmosphere_dim == 3  &&  rte_los.nelem() == 2 ) );

  // Jacobian code exists, but has not been tested at all
  if( jacobian_do )
    throw runtime_error( "Method not yet tested with jacobians." );      

  // Some sizes
  //
  const Index   nf  = f_grid.nelem();

  // Determine if there are any analytical jacobians to handle
  //
  Index j_analytical_do = 0;
  //
  if( jacobian_do ) { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }

  // If primary call, initilise iy_aux and diy_dx
  //
  if( iy_agenda_call1 )
    {
      // Error are considered to be zero
      iy_error.resize( 0, 0 );
      iy_error_type = 0;
      //
      if( iy_aux_do ) { iy_aux.resize( N_AUX, nf, stokes_dim ); iy_aux = -999; }
      else            { iy_aux.resize( 0, 0, 0 ); }
      //
      if( j_analytical_do ) 
        {
          diy_dx.resize( jacobian_indices.nelem() ); 
          //
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_dx[iq].resize( jacobian_indices[iq][1] - 
                               jacobian_indices[iq][0] + 1, nf, stokes_dim ); 
            diy_dx[iq] = 0.0;
          ) 
        }
      else { diy_dx.resize( 0 ); }
    }


  //- Determine propagation path
  //
  Ppath  ppath;
  //
  ppath_calc( ws, ppath, ppath_step_agenda, atmosphere_dim, p_grid, 
              lat_grid, lon_grid, z_field, r_geoid, z_surface,
              cloudbox_on, cloudbox_limits, rte_pos, rte_los, 1 );
  //
  const Index   i_background = ppath_what_background( ppath );
  assert( i_background > 0 );

  // Get quantities required for each ppath point
  //
  const Index     np  = ppath.np;
        Vector    ppath_p, ppath_t, total_tau;
        Matrix    ppath_vmr, mdummy, ppath_tau;
        Tensor3   ppath_abs_scalar, iy_trans_new;
        Agenda    adummy;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      //
      get_ptvmr_for_ppath( ppath_p, ppath_t, ppath_vmr, ppath, atmosphere_dim, 
                           p_grid, t_field, vmr_field );

      // Get emission, absorption and optical thickness for each step
      //
      get_step_vars_for_standardRT( ws, ppath_abs_scalar, mdummy, ppath_tau, 
                                    total_tau, abs_scalar_agenda, adummy,
                                    ppath, ppath_p, ppath_t, ppath_vmr, nf, 0 );
    }
  else // To handle cases inside the cloudbox, or totally outside the atmosphere
    {
      total_tau.resize( nf );
      total_tau = 0;
    }

  // Handle iy_transmission (the variable not needed when background is
  // space or no analytical jacobians will be calculated)
  //
  if( i_background > 1 && j_analytical_do )
    {
      if( iy_agenda_call1 )
        { iy_transmission_for_scalar_tau( iy_trans_new, stokes_dim, 
                                                                 total_tau ); }
      else
        { iy_transmission_mult_scalar_tau( iy_trans_new, iy_transmission, 
                                                                 total_tau ); }
    }

  // Radiative background
  //
  get_iy_of_background( 
                ws, iy, iy_error, iy_error_type, iy_aux, diy_dx, iy_trans_new,  
                iy_aux_do, jacobian_do, ppath, atmosphere_dim, p_grid, 
                lat_grid, lon_grid, t_field, vmr_field, r_geoid, z_surface, 
                cloudbox_on, stokes_dim, f_grid, iy_clearsky_agenda, 
                iy_space_agenda, surface_prop_agenda, iy_cloudbox_agenda );
  
  // Do RT calculations
  //
  if( np > 1 )
    {
      // Create container for the derivatives with respect to changes
      // at each ppath point. And find matching abs_species-index and 
      // "temperature flag" (for analytical jacobians).
      //
      ArrayOfTensor3  diy_dpath; 
      ArrayOfIndex    abs_species_i, is_t; 
      //
      const Numeric   dt = 0.1;
            Tensor3   ppath_as2;
      //
      if( j_analytical_do )
        { 
          diy_dpath.resize( jacobian_indices.nelem() ); 
          abs_species_i.resize( jacobian_indices.nelem() ); 
          is_t.resize( jacobian_indices.nelem() ); 
          //
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_dpath[iq].resize( np, nf, stokes_dim ); 
            diy_dpath[iq] = 0.0;
          )
          get_pointers_for_analytical_jacobians( abs_species_i, is_t, 
                                            jacobian_quantities, abs_species );
          //
          // Determine if temperature is among the analytical jac. quantities.
          // If yes, get emission and absorption for disturbed temperature
          Index do_t=0;
          for( Index iq=0; iq<is_t.nelem() && do_t==0; iq++ )
            { if( is_t[iq] ) { do_t = 1; } }
          if( do_t )
            {
              Matrix mdummy1, mdummy2; 
              Vector vdummy, t2 = ppath_t;
              t2 += dt;
              get_step_vars_for_standardRT( ws, ppath_as2, mdummy1, mdummy2, 
                                   vdummy, abs_scalar_agenda, adummy,
                                   ppath, ppath_p, t2, ppath_vmr, nf, 0 );
            }
        }

      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          // Loop frequencies
          for( Index iv=0; iv<nf; iv++ )
            { 
              // Jacobian quantities
              if( j_analytical_do )
                {
                  const Numeric A = ppath.l_step[ip] * exp( -total_tau[iv] );
                  
                  for( Index iq=0; iq<jacobian_quantities.nelem(); iq++ ) 
                    {
                      if( jacobian_quantities[iq].Analytical() )
                        {
                          // Absorbing species
                          const Index isp = abs_species_i[iq]; 
                          if( isp >= 0 )
                            {
                              // Scaling factors to handle retrieval unit
                              Numeric unitscf1, unitscf2;
                              vmrunitscf( unitscf1, 
                                          jacobian_quantities[iq].Mode(), 
                                          ppath_vmr(isp,ip), ppath_p[ip], 
                                          ppath_t[ip] );
                              vmrunitscf( unitscf2, 
                                          jacobian_quantities[iq].Mode(), 
                                          ppath_vmr(isp,ip+1), ppath_p[ip+1], 
                                          ppath_t[ip+1] );
                              //
                              // All Stokes components equal
                              for( Index is=0; is<stokes_dim; is++ )
                                { 
                                  const Numeric wi = -0.5 * A * iy(iv,is);
                                  diy_dpath[iq](ip  ,iv,is) += unitscf1 * wi *
                                                 ppath_abs_scalar(iv,isp,ip);
                                  diy_dpath[iq](ip+1,iv,is) += unitscf2 * wi *
                                                 ppath_abs_scalar(iv,isp,ip+1);
                                }
                            }

                          // Temperature
                          else if( is_t[iq] )
                            {
                              const Numeric dkdt1 =
                                ( ppath_as2(iv,joker,ip).sum()
                                - ppath_abs_scalar(iv,joker,ip).sum() ) / dt;
                              const Numeric dkdt2 =
                                ( ppath_as2(iv,joker,ip+1).sum()
                                - ppath_abs_scalar(iv,joker,ip+1).sum() ) / dt;
                              //
                              // All Stokes components equal
                              for( Index is=0; is<stokes_dim; is++ )
                                { 
                                  const Numeric wi = -0.5 * A * iy(iv,is);
                                  diy_dpath[iq](ip  ,iv,is) += wi * dkdt1;
                                  diy_dpath[iq](ip+1,iv,is) += wi * dkdt2;
                                }
                            }
                        }
                    }
                  
                  // Update total_tau (not used for spectrum calculation)
                  total_tau[iv] -= ppath_tau(iv,ip);
                }

              // Spectrum
              //
              const Numeric step_tr = exp( -ppath_tau(iv,ip) );
              //
              for( Index is=0; is<stokes_dim; is++ )
                { iy(iv,is) *= step_tr; }
            }
        } 

      // Map jacobians from ppath to retrieval grids
      if( j_analytical_do )
        { 
          // Weight with iy_transmission
          if( !iy_agenda_call1 )
            {
              Matrix X, Y(stokes_dim,diy_dpath[0].npages()); 
              //
              FOR_ANALYTICAL_JACOBIANS_DO( 
                for( Index iv=0; iv<nf; iv++ )
                  { 
                    X = transpose( diy_dpath[iq](joker,iv,joker) );
                    mult( Y, iy_transmission(iv,joker,joker), X );
                    diy_dpath[iq](joker,iv,joker) = transpose( Y );
                  }
               )
            }

          // Map to retrieval grids
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_from_path_to_rgrids( diy_dx[iq], jacobian_quantities[iq], 
                               diy_dpath[iq], atmosphere_dim, ppath, ppath_p );
          )
        }
    }

  // Set iy_aux 
  //
  if( iy_aux_do ) 
    {
      // Set background if primary call
      if( iy_agenda_call1 )
        { iy_aux( 0, joker, joker ) = (Numeric)i_background; }
      // 
      // Flag intersection with cloudbox
      if( i_background >= 3 )
        { iy_aux( 1, joker, joker ) = 1.0; }
    }    
}



/* Workspace method: Doxygen documentation will be auto-generated */
void iyBeerLambertStandardCloudbox(
        Workspace&                     ws,
        Matrix&                        iy,
        Matrix&                        iy_error,
        Index&                         iy_error_type,
        Tensor3&                       iy_aux,
        ArrayOfTensor3&                diy_dx,
  const Vector&                        rte_pos,
  const Vector&                        rte_los,
  const Index&                         iy_aux_do,
  const Index&                         jacobian_do,
  const Index&                         atmosphere_dim,
  const Vector&                        p_grid,
  const Vector&                        lat_grid,
  const Vector&                        lon_grid,
  const Tensor3&                       z_field,
  const Tensor3&                       t_field,
  const Tensor4&                       vmr_field,
  const Matrix&                        r_geoid,
  const Matrix&                        z_surface,
  const Index&                         cloudbox_on,
  const ArrayOfIndex&                  cloudbox_limits,
  const Index&                         stokes_dim,
  const Vector&                        f_grid,
  const Agenda&                        ppath_step_agenda,
  const Agenda&                        abs_scalar_gas_agenda,
  const Agenda&                        iy_clearsky_agenda,
  const Tensor4&                       pnd_field,
  const ArrayOfSingleScatteringData&   scat_data_raw,
  const Agenda&                        opt_prop_gas_agenda )
{
  // Input checks
  if( !cloudbox_on )
    throw runtime_error( "The cloudbox must be defined to use this method." );

  // Determine ppath through the cloudbox
  //
  Ppath  ppath;
  //
  ppath_calc( ws, ppath, ppath_step_agenda, atmosphere_dim, p_grid, 
              lat_grid, lon_grid, z_field, r_geoid, z_surface,
              cloudbox_on, cloudbox_limits, rte_pos, rte_los, 0 );

  // Check radiative background
  const Index bkgr = ppath_what_background( ppath );
  if( bkgr == 2 )
    throw runtime_error( "It is not possible to have an intersection with "
                         "the surface inside the cloudbox." );
  assert( bkgr == 3 );

  // Get iy for unscattered direction
  //
  // Note that the Ppath positions (ppath.pos) for 1D have one column more
  // than expected by most functions. Only the first atmosphere_dim values
  // shall be copied.
  //
  Vector rte_pos2, rte_los2;
  Tensor3 dummy;
  rte_pos2 = ppath.pos(ppath.np-1,Range(0,atmosphere_dim));
  rte_los2 = ppath.los(ppath.np-1,joker);
  //
  //
  iy_clearsky_agendaExecute( ws, iy, iy_error, iy_error_type,
                             iy_aux, diy_dx, 0, rte_pos2, rte_los2,
                             dummy, 0, jacobian_do, iy_aux_do, 
                             f_grid, p_grid, lat_grid, lon_grid, t_field, 
                             vmr_field, iy_clearsky_agenda );

  // RT for part inside cloudbox
  //
  const Index   np  = ppath.np;
  //
  if( np > 1 )
    {
      // Some sizes
      const Index   nf  = f_grid.nelem();
      const Index   nabs = vmr_field.nbooks();
      const Index   npar = pnd_field.nbooks();

      // Containers for atmospheric data along the LOS
      //
      Vector   ppath_t(np), ppath_logp(np);
      Matrix   ppath_vmr(nabs,np);
      Matrix   ppath_pnd(npar,np);

      // Determine log of the pressure
      {
        Vector   ppath_p( np );
        Matrix   itw_p(np,2);
        interpweights( itw_p, ppath.gp_p );
        itw2p( ppath_p, p_grid, ppath.gp_p, itw_p );
        transform( ppath_logp, log, ppath_p  );
      }
      
      // Temperature and VMR
      {
        Matrix   itw_field;
        //
        interp_atmfield_gp2itw( itw_field, atmosphere_dim, 
                                ppath.gp_p, ppath.gp_lat, ppath.gp_lon );
        interp_atmfield_by_itw( ppath_t,  atmosphere_dim, t_field, ppath.gp_p, 
                                ppath.gp_lat, ppath.gp_lon, itw_field );
        for( Index is=0; is<nabs; is++ )
          {
            interp_atmfield_by_itw( ppath_vmr(is,joker), atmosphere_dim, 
                                    vmr_field(is,joker,joker,joker), 
                                    ppath.gp_p, ppath.gp_lat, ppath.gp_lon, 
                                    itw_field );
          }
      }

      // PND
      {
        Matrix   itw_field;
        ArrayOfGridPos gp_p   = ppath.gp_p;
        ArrayOfGridPos gp_lat = ppath.gp_lat;
        ArrayOfGridPos gp_lon = ppath.gp_lon;
        //
        interp_cloudfield_gp2itw( itw_field, gp_p, gp_lat, gp_lon, 
                                  atmosphere_dim, cloudbox_limits );
        for( Index ip=0; ip<npar; ip++ )
          {
            interp_atmfield_by_itw( ppath_pnd(ip,joker), atmosphere_dim,
                                      pnd_field(ip,joker,joker,joker), 
                                      gp_p, gp_lat, gp_lon, itw_field );
          }
      }

      // Scattering properties
      //
      Array<ArrayOfSingleScatteringData>   scat_data( nf );
      //
      for( Index iv=0; iv<nf; iv++ )
        { scat_data_monoCalc( scat_data[iv], scat_data_raw, f_grid, iv ); }

      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          // Mean pressure, temperature, VMR and PND
          const Numeric   p_mean = exp( 0.5*(ppath_logp[ip+1]+ppath_logp[ip] ));
          const Numeric   t_mean = 0.5*( ppath_t[ip+1] + ppath_t[ip] );
          Vector    vmr_mean( nabs ), pnd_mean( npar );
          for( Index ia=0; ia<nabs; ia++ )
            { vmr_mean[ia] = 0.5 * ( ppath_vmr(ia,ip+1) + ppath_vmr(ia,ip) ); }
          for( Index ia=0; ia<npar; ia++ )
            { pnd_mean[ia] = 0.5 * ( ppath_pnd(ia,ip+1) + ppath_pnd(ia,ip) ); }

          // Gas absorption and blackbody radiation
          Matrix   abs_scalar_gas;
          abs_scalar_gas_agendaExecute( ws, abs_scalar_gas, -1, p_mean, t_mean, 
                                        vmr_mean, abs_scalar_gas_agenda );

          // Clearsky step
          if( max(pnd_mean) < 1e-3 )
            {
              // Loop frequencies
              for( Index iv=0; iv<nf; iv++ )
                { 
                  const Numeric step_tr = exp( -ppath.l_step[ip] * 
                                               abs_scalar_gas(iv,joker).sum() );
                  //
                  for( Index is=0; is<stokes_dim; is++ )
                    { iy(iv,is) *= step_tr; }
                }
            }

          
          // RT step with scattering
          else
            {
              // Mean position of step
              for( Index ia=0; ia<atmosphere_dim; ia++ )
                { rte_pos2[ia] = 0.5*( ppath.pos(ip+1,ia)+ppath.pos(ip,ia) ); }

              // Direction of outgoing scattered radiation
              // Note that rte_los2 is only used for extracting scattering
              // properties. 
              // A viewing direction of aa=0 is assumed for 1D. This
              // corresponds to positive za for 2D.
              //
              Numeric za_mean = 0.5*( ppath.los(ip+1,0)+ppath.los(ip,0) );
              //
              rte_los2.resize(2);
              //
              if( atmosphere_dim == 1 )
                { 
                  rte_los2[0] = 180 - za_mean; 
                  rte_los2[1] = 180; 
                }
              else if( atmosphere_dim == 2 )
                {
                  rte_los2[0] = 180 - fabs( za_mean ); 
                  if( za_mean >= 0 )
                    { rte_los2[1] = 180; }
                  else
                    { rte_los2[1] = 0; }
                }
              else if( atmosphere_dim == 3 )
                { 
                  rte_los2[0] = 180 - za_mean; 
                  Numeric aa_mean = 0.5*( ppath.los(ip+1,1)+ppath.los(ip,1) );
                  rte_los2[1] = aa_mean + 180; 
                  if( rte_los2[1] > 180 )
                    { rte_los2[1] -= 360; }
                }

              // Extinction matrix and absorption vector for gases
              //
              Tensor3  ext_mat_gas;
              Matrix   abs_vec_gas;
              //
              opt_prop_gas_agendaExecute( ws, ext_mat_gas, abs_vec_gas, -1,
                                         abs_scalar_gas, opt_prop_gas_agenda );

              // Loop frequencies
              //
              Matrix ext_mat( stokes_dim, stokes_dim );
              Vector abs_vec( stokes_dim );
              Matrix trans_mat(stokes_dim,stokes_dim); 
              //
              for( Index iv=0; iv<nf; iv++ )
                { 
                  // Particle extinction
                  opt_propCalc( ext_mat, abs_vec, rte_los2[0], 
                                rte_los2[1], scat_data[iv], stokes_dim, 
                                    pnd_mean, t_mean );

                  // Add gas part
                  ext_mat += ext_mat_gas(iv,joker,joker);

                  // RT for step
                  //
                  trans_step_std( iy(iv,joker), trans_mat, 
                                  ext_mat, ppath.l_step[ip] );
                }
            }
        } 
    }

}



/* Workspace method: Doxygen documentation will be auto-generated */
void iyEmissionStandardClearsky(
        Workspace&                  ws,
        Matrix&                     iy,
        Matrix&                     iy_error,
        Index&                      iy_error_type,
        Tensor3&                    iy_aux,
        ArrayOfTensor3&             diy_dx,
  const Index&                      iy_agenda_call1,
  const Vector&                     rte_pos,      
  const Vector&                     rte_los,      
  const Index&                      iy_aux_do,
  const Index&                      jacobian_do,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Tensor3&                    z_field,
  const Tensor3&                    t_field,
  const Tensor4&                    vmr_field,
  const Matrix&                     r_geoid,
  const Matrix&                     z_surface,
  const Index&                      cloudbox_on,
  const ArrayOfIndex&               cloudbox_limits,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const ArrayOfArrayOfSpeciesTag&   abs_species,
  const Agenda&                     ppath_step_agenda,
  const Agenda&                     emission_agenda,
  const Agenda&                     abs_scalar_agenda,
  const Agenda&                     iy_clearsky_agenda,
  const Agenda&                     iy_space_agenda,
  const Agenda&                     surface_prop_agenda,
  const Agenda&                     iy_cloudbox_agenda,
  const Tensor3&                    iy_transmission,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices )
{
  // A minimal amount of checks, for efficiency reasons
  assert( rte_pos.nelem() == atmosphere_dim );
  assert( ( atmosphere_dim < 3   &&  rte_los.nelem() == 1 )  ||
          ( atmosphere_dim == 3  &&  rte_los.nelem() == 2 ) );

  // Some sizes
  //
  const Index   nf  = f_grid.nelem();

  // Determine if there are any analytical jacobians to handle
  //
  Index j_analytical_do = 0;
  //
  if( jacobian_do ) { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }

  // If primary call, initilise iy_aux and diy_dx
  //
  if( iy_agenda_call1 )
    {
      // Error are considered to be zero
      iy_error.resize( 0, 0 );
      iy_error_type = 0;
      //
      if( iy_aux_do ) { iy_aux.resize( N_AUX, nf, stokes_dim ); iy_aux = 0.0; }
      else            { iy_aux.resize( 0, 0, 0 ); }
      //
      if( j_analytical_do ) 
        {
          diy_dx.resize( jacobian_indices.nelem() ); 
          //
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_dx[iq].resize( jacobian_indices[iq][1] - 
                               jacobian_indices[iq][0] + 1, nf, stokes_dim ); 
            diy_dx[iq] = 0.0;
          ) 
        }
      else { diy_dx.resize( 0 ); }
    }

  //- Determine propagation path
  //
  Ppath  ppath;
  //
  ppath_calc( ws, ppath, ppath_step_agenda, atmosphere_dim, p_grid, 
              lat_grid, lon_grid, z_field, r_geoid, z_surface,
              cloudbox_on, cloudbox_limits, rte_pos, rte_los, 1 );
  //
  const Index   i_background = ppath_what_background( ppath );
  assert( i_background > 0 );

  // Get quantities required for each ppath point
  //
  const Index     np  = ppath.np;
        Vector    ppath_p, ppath_t, total_tau;
        Matrix    ppath_vmr, ppath_emission, ppath_tau;
        Tensor3   ppath_abs_scalar, iy_trans_new;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      //
      get_ptvmr_for_ppath( ppath_p, ppath_t, ppath_vmr, ppath, atmosphere_dim, 
                           p_grid, t_field, vmr_field );

      // Get emission, absorption, optical thickness for each step, and total
      // optical thickness
      //
      get_step_vars_for_standardRT( ws, ppath_abs_scalar, ppath_emission, 
                                    ppath_tau, total_tau, abs_scalar_agenda, 
                                    emission_agenda, ppath, ppath_p, ppath_t, 
                                    ppath_vmr, nf, 1 );
    }
  else // To handle cases inside the cloudbox, or totally outside the atmosphere
    { 
      total_tau.resize( nf );
      total_tau = 0;
    }

  // Handle iy_transmission (the variable not needed when background is
  // space or no analytical jacobians will be calculated, but always with
  // iy_aux_do set to true)
  //
  if( ( i_background > 1 && j_analytical_do )  ||  iy_aux_do )
    {
      if( iy_agenda_call1 )
        { iy_transmission_for_scalar_tau( iy_trans_new, stokes_dim, 
                                                                 total_tau ); }
      else
        { iy_transmission_mult_scalar_tau( iy_trans_new, iy_transmission, 
                                                                 total_tau ); }
    }

  // Radiative background
  //
  get_iy_of_background( 
                ws, iy, iy_error, iy_error_type, iy_aux, diy_dx, iy_trans_new,  
                iy_aux_do, jacobian_do, ppath, atmosphere_dim, p_grid, 
                lat_grid, lon_grid, t_field, vmr_field, r_geoid, z_surface, 
                cloudbox_on, stokes_dim, f_grid, iy_clearsky_agenda, 
                iy_space_agenda, surface_prop_agenda, iy_cloudbox_agenda );
  
  // Do RT calculations
  //
  if( np > 1 )
    {
      // Create container for the derivatives with respect to changes
      // at each ppath point. And find matching abs_species-index and 
      // "temperature flag" (for analytical jacobians).
      //
      ArrayOfTensor3  diy_dpath; 
      ArrayOfIndex    abs_species_i, is_t; 
      //
      const Numeric   dt = 0.1;
            Matrix    ppath_e2;
            Tensor3   ppath_as2;
      //
      if( j_analytical_do )
        { 
          diy_dpath.resize( jacobian_indices.nelem() ); 
          abs_species_i.resize( jacobian_indices.nelem() ); 
          is_t.resize( jacobian_indices.nelem() ); 
          //
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_dpath[iq].resize( np, nf, stokes_dim ); 
            diy_dpath[iq] = 0.0;
          )
          get_pointers_for_analytical_jacobians( abs_species_i, is_t, 
                                            jacobian_quantities, abs_species );
          //
          // Determine if temperature is among the analytical jac. quantities.
          // If yes, get emission and absorption for disturbed temperature
          Index do_t=0;
          for( Index iq=0; iq<is_t.nelem() && do_t==0; iq++ )
            { if( is_t[iq] ) { do_t = 1; } }
          if( do_t )
            {
              Matrix mdummy; 
              Vector vdummy, t2 = ppath_t;
              t2 += dt;
              get_step_vars_for_standardRT( ws, ppath_as2, ppath_e2, mdummy, 
                                   vdummy, abs_scalar_agenda, emission_agenda,
                                   ppath, ppath_p, t2, ppath_vmr, nf, 1 );
            }
        }

      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          // Loop frequencies
          for( Index iv=0; iv<nf; iv++ )
            { 
              // Average of emission source term
              Numeric esource = 0.5 * ( ppath_emission(iv,ip) + 
                                       ppath_emission(iv,ip+1) );

              // Jacobian quantities
              if( j_analytical_do )
                {
                  const Numeric A = ppath.l_step[ip] * exp( -total_tau[iv] );
                  const Numeric B = 0.5 * A * (esource-iy(iv,0));
                  
                  for( Index iq=0; iq<jacobian_quantities.nelem(); iq++ ) 
                    {
                      if( jacobian_quantities[iq].Analytical() )
                        {
                          // Absorbing species
                          const Index isp = abs_species_i[iq]; 
                          if( isp >= 0 )
                            {
                              // Scaling factors to handle retrieval unit
                              Numeric unitscf1, unitscf2;
                              vmrunitscf( unitscf1, 
                                          jacobian_quantities[iq].Mode(), 
                                          ppath_vmr(isp,ip), ppath_p[ip], 
                                          ppath_t[ip] );
                              vmrunitscf( unitscf2, 
                                          jacobian_quantities[iq].Mode(), 
                                          ppath_vmr(isp,ip+1), ppath_p[ip+1], 
                                          ppath_t[ip+1] );
                              //
                              // Stokes component 1
                              diy_dpath[iq](ip  ,iv,0) += unitscf1 * B *
                                                 ppath_abs_scalar(iv,isp,ip);
                              diy_dpath[iq](ip+1,iv,0) += unitscf2 * B * 
                                                 ppath_abs_scalar(iv,isp,ip+1);
                              //
                              // Higher stokes components
                              for( Index is=1; is<stokes_dim; is++ )
                                { 
                                  const Numeric wi = -0.5 * A * iy(iv,is); 
                                  diy_dpath[iq](ip  ,iv,is) += unitscf1 * wi *
                                                 ppath_abs_scalar(iv,isp,ip);
                                  diy_dpath[iq](ip+1,iv,is) += unitscf2 * wi *
                                                 ppath_abs_scalar(iv,isp,ip+1);
                                }
                            }

                          // Temperature
                          else if( is_t[iq] )
                            {
                              // Stokes component 1
                              //
                              // The terms associated with tau:
                              const Numeric dkdt1 =
                                ( ppath_as2(iv,joker,ip).sum()
                                - ppath_abs_scalar(iv,joker,ip).sum() ) / dt;
                              const Numeric dkdt2 =
                                ( ppath_as2(iv,joker,ip+1).sum()
                                - ppath_abs_scalar(iv,joker,ip+1).sum() ) / dt;
                              diy_dpath[iq](ip  ,iv,0) += B * dkdt1;
                              diy_dpath[iq](ip+1,iv,0) += B * dkdt2;
                              //
                              // The terms associated with B-bar:
                              const Numeric w0 = 0.5 * 
                                      exp(-(total_tau[iv]-ppath_tau(iv,ip))) *
                                              ( 1.0 - exp(-ppath_tau(iv,ip)) );
                              diy_dpath[iq](ip  ,iv,0) += w0 *
                                              ( ppath_e2(iv,ip) -
                                                ppath_emission(iv,ip) ) / dt;
                              diy_dpath[iq](ip+1,iv,0) += w0 * 
                                              ( ppath_e2(iv,ip+1) -
                                                ppath_emission(iv,ip+1) ) / dt;

                              // Higher stokes components
                              for( Index is=1; is<stokes_dim; is++ )
                                { 
                                  const Numeric wi = -0.5 * A * iy(iv,is);
                                  diy_dpath[iq](ip  ,iv,is) += wi * dkdt1;
                                  diy_dpath[iq](ip+1,iv,is) += wi * dkdt2;
                                }
                            }
                        }
                    }
                  
                  // Update total_tau (not used for spectrum calculation)
                  total_tau[iv] -= ppath_tau(iv,ip);
                }

              // Spectrum
              //
              const Numeric step_tr = exp( -ppath_tau(iv,ip) );
              //
              iy(iv,0)  = iy(iv,0) * step_tr + esource * (1-step_tr); 
              //
              for( Index is=1; is<stokes_dim; is++ )
                { iy(iv,is) *= step_tr; }
            }
        } 

      // Map jacobians from ppath to retrieval grids
      if( j_analytical_do )
        { 
          // Weight with iy_transmission
          if( !iy_agenda_call1 )
            {
              Matrix X, Y(stokes_dim,diy_dpath[0].npages()); 
              //
              FOR_ANALYTICAL_JACOBIANS_DO( 
                for( Index iv=0; iv<nf; iv++ )
                  { 
                    X = transpose( diy_dpath[iq](joker,iv,joker) );
                    mult( Y, iy_transmission(iv,joker,joker), X );
                    diy_dpath[iq](joker,iv,joker) = transpose( Y );
                  }
               )
            }

          // Map to retrieval grids
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_from_path_to_rgrids( diy_dx[iq], jacobian_quantities[iq], 
                               diy_dpath[iq], atmosphere_dim, ppath, ppath_p );
          )
        }
    }

  // Set iy_aux 
  //
  if( iy_aux_do ) 
    {
      // Fill transmission if background is TOA or surface
      if( i_background <= 2 )
        {
          for( Index iv=0; iv<nf; iv++ )
            { 
              for( Index is=0; is<stokes_dim; is++ )
                { iy_aux( 0, iv, is ) = iy_trans_new( iv, is, is ); }
            }
        }

      // Set background if primary call
      if( iy_agenda_call1 )
        { iy_aux( 1, joker, joker ) = (Numeric)i_background; }
      // 
      // Flag intersection with cloudbox
      if( i_background >= 3 )
        { iy_aux( 2, joker, joker ) = 1.0; }
    }    
}



/* Workspace method: Doxygen documentation will be auto-generated */
void iyEmissionStandardClearskyBasic(
        Workspace&                  ws,
        Matrix&                     iy,
  const Vector&                     rte_pos,      
  const Vector&                     rte_los,      
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Tensor3&                    z_field,
  const Tensor3&                    t_field,
  const Tensor4&                    vmr_field,
  const Matrix&                     r_geoid,
  const Matrix&                     z_surface,
  const Index&                      cloudbox_on,
  const ArrayOfIndex&               cloudbox_limits,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const Agenda&                     ppath_step_agenda,
  const Agenda&                     emission_agenda,
  const Agenda&                     abs_scalar_agenda,
  const Agenda&                     iy_clearsky_basic_agenda,
  const Agenda&                     iy_space_agenda,
  const Agenda&                     surface_prop_agenda,
  const Agenda&                     iy_cloudbox_agenda )
{
  // A minimal amount of checks, for efficiency reasons
  assert( rte_pos.nelem() == atmosphere_dim );
  assert( ( atmosphere_dim < 3   &&  rte_los.nelem() == 1 )  ||
          ( atmosphere_dim == 3  &&  rte_los.nelem() == 2 ) );


  //- Determine propagation path
  //
  Ppath  ppath;
  //
  ppath_calc( ws, ppath, ppath_step_agenda, atmosphere_dim, p_grid, 
              lat_grid, lon_grid, z_field, r_geoid, z_surface,
              cloudbox_on, cloudbox_limits, rte_pos, rte_los, 1 );


  // Get quantities required for each ppath point
  //
  // If np = 1, we only need to determine the radiative background
  //
  const Index     nf  = f_grid.nelem();
  const Index     np  = ppath.np;
        Vector    ppath_p, ppath_t, total_tau;
        Matrix    ppath_vmr, ppath_emission, ppath_tau;
        Tensor3   ppath_abs_scalar;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      //
      get_ptvmr_for_ppath( ppath_p, ppath_t, ppath_vmr, ppath, atmosphere_dim, 
                           p_grid, t_field, vmr_field );

      // Get emission, absorption and optical thickness for each step
      //
      get_step_vars_for_standardRT( ws, ppath_abs_scalar, ppath_emission, 
                                    ppath_tau, total_tau, 
                                    abs_scalar_agenda, emission_agenda,
                                    ppath, ppath_p, ppath_t, ppath_vmr, nf, 1 );
    }


  // Radiative background
  //
  Matrix          dummy1;
  Index           dummy2;
  Tensor3         dummy3, dummy5;
  ArrayOfTensor3  dummy4;
  //  
  get_iy_of_background( ws, iy, dummy1, dummy2, dummy3, dummy4, dummy5,  
                    0, 0, ppath, atmosphere_dim, p_grid, lat_grid, lon_grid,
                    t_field, vmr_field, r_geoid, z_surface, cloudbox_on, 
                    stokes_dim, f_grid, iy_clearsky_basic_agenda, 
                    iy_space_agenda, surface_prop_agenda, iy_cloudbox_agenda );


  // Do RT calculations
  //
  if( np > 1 )
    {
      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          // Loop frequencies
          for( Index iv=0; iv<nf; iv++ )
            { 

              const Numeric step_tr = exp( -ppath_tau(iv,ip) );
              //
              iy(iv,0)  = iy(iv,0) * step_tr + 0.5 * ( ppath_emission(iv,ip) + 
                                       ppath_emission(iv,ip+1) ) * (1-step_tr); 
              //
              for( Index is=1; is<stokes_dim; is++ )
                { iy(iv,is) *= step_tr; }
            }
        } 
    }
}


 


/* Workspace method: Doxygen documentation will be auto-generated */
void iyMC(
        Workspace&                  ws,
        Matrix&                     iy,
        Matrix&                     iy_error,
        Index&                      iy_error_type,
        Tensor3&                    iy_aux,
        ArrayOfTensor3&             diy_dx,
  const Index&                      iy_agenda_call1,
  const Index&                      atm_checked,
  const Vector&                     rte_pos,      
  const Vector&                     rte_los,      
  const Index&                      iy_aux_do,
  const Index&                      jacobian_do,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Tensor3&                    z_field,
  const Tensor3&                    t_field,
  const Tensor4&                    vmr_field,
  const Matrix&                     r_geoid,
  const Matrix&                     z_surface,
  const Index&                      cloudbox_on,
  const ArrayOfIndex&               cloudbox_limits,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const ArrayOfSingleScatteringData&   scat_data_raw,
  const Agenda&                     iy_space_agenda,
  const Agenda&                     surface_prop_agenda,
  const Agenda&                     abs_scalar_gas_agenda, 
  const Agenda&                     opt_prop_gas_agenda,
  const Tensor4&                    pnd_field,
  const Tensor3&                    iy_transmission,
  const String&                     y_unit,
  const Numeric&                    mc_std_err,
  const Index&                      mc_max_time,
  const Index&                      mc_max_iter,
  const Index&                      mc_z_field_is_1D )
{
  // Throw error if unsupported features are requested
  if( atmosphere_dim != 3 )
    throw runtime_error( 
                "Only 3D atmospheres are allowed (atmosphere_dim must be 3)" );
  if( !cloudbox_on )
    throw runtime_error( 
                    "The cloudbox must be activated (cloudbox_on must be 1)" );
  if( iy_aux_do )
    throw runtime_error( 
      "This method does not provide any auxilary data (iy_aux_do must be 0)" );
  if( jacobian_do )
    throw runtime_error( 
        "This method does not provide any jacobians (jacobian_do must be 0)" );
  if( !iy_agenda_call1 )
    throw runtime_error( 
                  "Recursive usage not possible (iy_agenda_call1 must be 1)" );
  if( iy_transmission.ncols() )
    throw runtime_error( "*iy_transmission* must be empty" );

  // Size output variables
  //
  const Index   nf  = f_grid.nelem();
  //
  iy.resize( nf, stokes_dim );
  iy_error.resize( nf, stokes_dim );
  iy_error_type = 1;
  //
  iy_aux.resize( 0, 0, 0 );
  diy_dx.resize(0);

  // Some MC variables are only local here
  Tensor3  mc_points;
  Index    mc_iteration_count, mc_seed;
  //
  MCAntenna mc_antenna;
  mc_antenna.set_pencil_beam();

  // Pos and los must be matrices 
  Matrix pos(1,3), los(1,2);
  //
  pos(0,joker) = rte_pos;
  los(0,joker) = rte_los;

  for( Index f_index=0; f_index<nf; f_index++ )
    {
      ArrayOfSingleScatteringData   scat_data_mono;
      
      scat_data_monoCalc( scat_data_mono, scat_data_raw, 
                                          f_grid, f_index );

      // Seed reset for each loop. If not done, the errors 
      // appear to be highly correlated.
      MCSetSeedFromTime( mc_seed );

      Vector y, mc_error;
                  
      MCGeneral( ws, y, mc_iteration_count, mc_error, mc_points, 
                 mc_antenna, f_grid, f_index, pos, los, stokes_dim, 
                 iy_space_agenda, surface_prop_agenda, opt_prop_gas_agenda,
                 abs_scalar_gas_agenda, p_grid, lat_grid, lon_grid, z_field, 
                 r_geoid, z_surface, t_field, vmr_field, 
                 cloudbox_on, cloudbox_limits, 
                 pnd_field, scat_data_mono, atm_checked, mc_seed, y_unit, 
                 mc_std_err, mc_max_time, mc_max_iter, mc_z_field_is_1D ); 

      assert( y.nelem() == stokes_dim );

      // Data returned can not be in Tb
      if ( y_unit == "RJBT" )
        { 
          const Numeric scfac = invrayjean( 1, f_grid[f_index] );
          y /= scfac;
          mc_error /= scfac;
        }
     
      iy(f_index,joker) = y;
      iy_error(f_index,joker) = mc_error;
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void yCalc(
        Workspace&                  ws,
        Vector&                     y,
        Vector&                     y_f,
        ArrayOfIndex&               y_pol,
        Matrix&                     y_pos,
        Matrix&                     y_los,
        Vector&                     y_error,
        Matrix&                     y_aux,
        Matrix&                     jacobian,
  const Index&                      atm_checked,
  const Index&                      atmosphere_dim,
  const Vector&                     p_grid,
  const Vector&                     lat_grid,
  const Vector&                     lon_grid,
  const Tensor3&                    t_field,
  const Tensor4&                    vmr_field,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  const Vector&                     f_grid,
  const Matrix&                     sensor_pos,
  const Matrix&                     sensor_los,
  const Vector&                     mblock_za_grid,
  const Vector&                     mblock_aa_grid,
  const Index&                      antenna_dim,
  const Sparse&                     sensor_response,
  const Vector&                     sensor_response_f,
  const ArrayOfIndex&               sensor_response_pol,
  const Vector&                     sensor_response_za,
  const Vector&                     sensor_response_aa,
  const Agenda&                     iy_clearsky_agenda,
  const Index&                      iy_aux_do,
  const String&                     y_unit,
  const Agenda&                     jacobian_agenda,
  const Index&                      jacobian_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices )
{
  // Some sizes
  const Index   nf      = f_grid.nelem();
  const Index   nza     = mblock_za_grid.nelem();
        Index   naa     = mblock_aa_grid.nelem();   
  if( antenna_dim == 1 )  
    { naa = 1; }
  const Index   n1y     = sensor_response.nrows();
  const Index   nmblock = sensor_pos.nrows();
  const Index   niyb    = nf * nza * naa * stokes_dim;


  //---------------------------------------------------------------------------
  // Input checks
  //---------------------------------------------------------------------------

  // Atmosphere OK?
  //
  if( !atm_checked )
    throw runtime_error( "The atmosphere must be flagged to have passed a "
                         "consistency check (atm_checked=1)." );

  // Basic dimensionalities
  //
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );
  chk_if_in_range( "antenna_dim", antenna_dim, 1, 2 );

  // Sensor position and LOS.
  //
  if( sensor_pos.ncols() != atmosphere_dim )
    throw runtime_error( "The number of columns of sensor_pos must be "
                         "equal to the atmospheric dimensionality." );
  if( atmosphere_dim <= 2  &&  sensor_los.ncols() != 1 )
    throw runtime_error( "For 1D and 2D, sensor_los shall have one column." );
  if( atmosphere_dim == 3  &&  sensor_los.ncols() != 2 )
    throw runtime_error( "For 3D, sensor_los shall have two columns." );
  if( sensor_los.nrows() != nmblock )
    {
      ostringstream os;
      os << "The number of rows of sensor_pos and sensor_los must be "
         << "identical, but sensor_pos has " << nmblock << " rows,\n"
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

  // Antenna
  //
  if( nza == 0 )
    throw runtime_error( "The measurement block zenith angle grid is empty." );
  chk_if_increasing( "mblock_za_grid", mblock_za_grid );
  //
  if( antenna_dim == 1 )
    {
      if( mblock_aa_grid.nelem() != 0 )
        throw runtime_error( 
          "For antenna_dim = 1, the azimuthal angle grid must be empty." );
    }
  else
    {
      if( atmosphere_dim < 3 )
        throw runtime_error( "2D antennas (antenna_dim=2) can only be "
                                                 "used with 3D atmospheres." );
      if( mblock_aa_grid.nelem() == 0 )
        throw runtime_error(
                      "The measurement block azimuthal angle grid is empty." );
      chk_if_increasing( "mblock_aa_grid", mblock_aa_grid );
    }

  // Sensor
  //
  if( sensor_response.ncols() != niyb ) 
    {
      ostringstream os;
      os << "The *sensor_response* matrix does not have the right size,\n"
         << "either the method *sensor_responseInit* has not been run or some\n"
         << "of the other sensor response methods has not been correctly\n"
         << "configured.";
      throw runtime_error( os.str() );
    }

  // Sensor aux variables
  //
  if( n1y != sensor_response_f.nelem()  || n1y != sensor_response_pol.nelem() ||
      n1y != sensor_response_za.nelem() || n1y != sensor_response_za.nelem() )
    {
      ostringstream os;
      os << "Sensor auxiliary variables do not have the correct size.\n"
         << "The following variables should all have same size:\n"
         << "length(y) for one block      : " << n1y << "\n"
         << "sensor_response_f.nelem()    : " << sensor_response_f.nelem()
         << "\nsensor_response_pol.nelem(): " << sensor_response_pol.nelem()
         << "\nsensor_response_za.nelem() : " << sensor_response_za.nelem() 
         << "\nsensor_response_aa.nelem() : " << sensor_response_za.nelem() 
         << "\n";
      throw runtime_error( os.str() );
    }


  //---------------------------------------------------------------------------
  // Allocations and resizing
  //---------------------------------------------------------------------------

  // Resize *y* and *y_XXX*
  //
  y.resize( nmblock*n1y );
  y_f.resize( nmblock*n1y );
  y_pol.resize( nmblock*n1y );
  y_pos.resize( nmblock*n1y, sensor_pos.ncols() );
  y_los.resize( nmblock*n1y, sensor_los.ncols() );
  y_error.resize( nmblock*n1y );
  if( iy_aux_do )
    y_aux.resize( nmblock*n1y, N_AUX );
  else 
    y_aux.resize( 0, 0 );

  // Jacobian variables
  //
  Index  j_analytical_do = 0;
  //
  if( jacobian_do )
    {
      jacobian.resize( nmblock*n1y, 
                           jacobian_indices[jacobian_indices.nelem()-1][1]+1 );
      jacobian = 0;
      //
      FOR_ANALYTICAL_JACOBIANS_DO(
        j_analytical_do  = 1; 
      )
    }


  //---------------------------------------------------------------------------
  // The calculations
  //---------------------------------------------------------------------------

  // We have to make a local copy of the Workspace and the agendas because
  // only non-reference types can be declared firstprivate in OpenMP
  Workspace l_ws (ws);
  Agenda l_jacobian_agenda (jacobian_agenda);
  Agenda l_iy_clearsky_agenda (iy_clearsky_agenda);

/*#pragma omp parallel for                                           \
  if(!arts_omp_in_parallel() && nmblock>1 && nmblock>=nza)        \
    default(none)                                                   \
    firstprivate(l_ws, l_jacobian_agenda, l_iy_clearsky_agenda)     \
    shared(j_analytical_do, sensor_los, mblock_za_grid, mblock_aa_grid, \
           vmr_field, t_field, lon_grid, lat_grid, p_grid, f_grid,      \
           sensor_pos, joker, naa)*/
#pragma omp parallel for                                          \
  if(!arts_omp_in_parallel() && nmblock>1 && nmblock>=nza)        \
  firstprivate(l_ws, l_jacobian_agenda, l_iy_clearsky_agenda)
  for( Index imblock=0; imblock<nmblock; imblock++ )
    {
      // Calculate monochromatic pencil beam data for 1 measurement block
      //
      Vector          iyb, iyb_error, yb(n1y);
      Matrix          iyb_aux;
      Index           iy_error_type;
      ArrayOfMatrix   diyb_dx;      
      //
      iyb_calc( l_ws, iyb, iyb_error, iy_error_type, iyb_aux, diyb_dx, 
                imblock, atmosphere_dim, p_grid, lat_grid, lon_grid, t_field, 
                vmr_field, cloudbox_on, stokes_dim, f_grid, sensor_pos, 
                sensor_los, mblock_za_grid, mblock_aa_grid, antenna_dim, 
                l_iy_clearsky_agenda, iy_aux_do, y_unit, j_analytical_do, 
                jacobian_quantities, jacobian_indices );


      // Apply sensor response matrix on iyb, and put into y
      //
      const Range rowind = get_rowindex_for_mblock( sensor_response, imblock );
      const Index row0   = rowind.get_start();
      //
      mult( yb, sensor_response, iyb );
      //
      y[rowind] = yb;  // *yb* also used below, as input to jacobian_agenda

      // Apply sensor response matrix on iyb_error, and put into y_error
      //
      if( iy_error_type == 0 )
        { y_error[rowind] = 0; }
      else if( iy_error_type == 1 )
        { 
          for( Index i=0; i<n1y; i++ )
            {
              const Index ithis = row0+i;
              y_error[ithis] = 0;
              for( Index icol=0; icol<niyb; icol++ )
                { 
                  y_error[ithis] += pow( 
                       sensor_response(i,icol)*iyb_error[icol], (Numeric)2.0 ); 
                }
              y_error[ithis] = sqrt( y_error[ithis] );              
            }
        }
      else if( iy_error_type == 2 )
        { mult( y_error[rowind], sensor_response, iyb ); }
      else
        { assert( 0 ); }

      // Apply sensor response matrix on iyb_aux, and put into y_aux
      //
      if( iy_aux_do )
        { mult( y_aux(rowind,joker), sensor_response, iyb_aux ); }

      // Fill information variables
      //
      for( Index i=0; i<n1y; i++ )
        { 
          y_f[row0+i]         = sensor_response_f[i];
          y_pol[row0+i]       = sensor_response_pol[i]; 
          y_pos(row0+i,joker) = sensor_pos(imblock,joker);
          y_los(row0+i,0)     = sensor_los(imblock,0) + sensor_response_za[i];
          if( sensor_response_aa.nelem() )
            { 
              y_los(row0+i,1) = sensor_los(imblock,0) + sensor_response_aa[i]; 
            }
        }

      // Apply sensor response matrix on diyb_dx, and put into jacobian
      // (that is, analytical jacobian part)
      //
      if( j_analytical_do )
        {
          FOR_ANALYTICAL_JACOBIANS_DO(
            mult( jacobian(rowind, Range(jacobian_indices[iq][0],
                          jacobian_indices[iq][1]-jacobian_indices[iq][0]+1)),
                                                sensor_response, diyb_dx[iq] );
          )
        }

      // Rest of *jacobian*: run jacobian_agenda (can be empty)
      //
      if( jacobian_do  &&  jacobian_agenda.nelem() > 0 )
        { 
          jacobian_agendaExecute( l_ws, jacobian, imblock, iyb, yb, 
                                                            l_jacobian_agenda );
        }
    }  // End mblock loop
}





/* Workspace method: Doxygen documentation will be auto-generated */
void y_unitApply(
        Vector&         y,
        Vector&         y_error,
        Matrix&         jacobian,
  const Vector&         y_f,
  const ArrayOfIndex&   y_pol,
  const String&         y_unit )
{
  if( y_unit == "1" )
    { throw runtime_error(
        "No need to use this method with *y_unit* = \"1\"." ); }

  if( max(y) > 1e-3 )
    {
      ostringstream os;
      os << "The spectrum vector *y* is required to have original radiance\n"
         << "unit, but this seems not to be the case. This as a value above\n"
         << "1e-3 is found in *y*.";
      throw runtime_error( os.str() );      
    }

  // Are jacobian and y_error set?
  //
  const Index ny = y.nelem();
  //
  const bool do_j = jacobian.nrows() == ny;
  const bool do_e = y_error.nelem() == ny;

  // Some jacobian quantities can not be handled
  if( do_j  &&  max(jacobian) > 1e-3 )
    {
      ostringstream os;
      os << "The method can not be used with jacobian quantities that are not\n"
         << "obtained through radiative transfer calculations. One example on\n"
         << "quantity that can not be handled is *jacobianAddPolyfit*.\n"
         << "The maximum value of *jacobian* indicates that one or several\n"
         << "such jacobian quantities are included.";
      throw runtime_error( os.str() );      
    }

  // All conversions could be handled as Planck-BT, but before everything is
  // tested properly and it is clear which version that is most efficient, both
  // versions are kept

  // Planck-Tb 
  //-------------------------------------------------------------------------- 
  if( y_unit == "PlanckBT" )
    {
      // Hard to use telescoping here as the data are sorted differently in y
      // and jacobian, than what is expected apply_y_unit. Copy to temporary
      // variables instead.

      // Handle the elements in "frequency chunks"

      Index i0 = 0;           // Index of first element for present chunk
      //
      while( i0 < ny )
        { 
          // Find number values for this chunk
          Index n = 1;
          //
          while( i0+n < ny &&  y_f[i0] == y_f[i0+n] ) 
            { n++; }                              

          Matrix yv(1,n);  
          ArrayOfIndex i_pol(n);
          bool any_quv = false;
          //
          for( Index i=0; i<n; i++ )
            { 
              const Index ix=i0+i;   
              yv(0,i) = y[ix];   
              i_pol[i] = y_pol[ix]; 
              if( i_pol[i] > 1  &&  i_pol[i] < 5 )
                { any_quv = true; }
            }

          // Index of elements to convert
          Range ii( i0, n );

          if( do_e || do_j )
            {
              if( any_quv  &&  i_pol[0] != 1 )
                {
                  ostringstream os;
                  os << "The conversion to PlanckBT, of the Jacobian and "
                     << "errors for Q, U and V, requires that I (first Stokes "
                     << "element) is at hand and that the data are sorted in "
                     << "such way that I comes first for each frequency.";
                  throw runtime_error( os.str() );      
                }

              // y_error
              if( do_e )
                {
                  Tensor3 e(1,1,n);
                  e(0,0,joker) = y_error[ii];
                  apply_y_unit2( e, yv, y_unit, y_f[i0], i_pol ); 
                  y_error[ii] = e(0,0,joker);
                }

              // Jacobian
              if( do_j )
                {
                  Tensor3 J(jacobian.ncols(),1,n);
                  J(joker,0,joker) = transpose( jacobian(ii,joker) );
                  apply_y_unit2( J, yv, y_unit, y_f[i0], i_pol ); 
                  jacobian(ii,joker) = transpose( J(joker,0,joker) );
                }
            }

          // y (must be done last)
          apply_y_unit( yv, y_unit, y_f[i0], i_pol );
          y[ii] = yv(0,joker);

          i0 += n;
        }
    }


  // Other conversions
  //-------------------------------------------------------------------------- 
  else
    {
      // Here we take each element of y separately. 

      Matrix yv(1,1);
      ArrayOfIndex i_pol(1);

      for( Index i=0; i<ny; i++ )
        {
          yv(0,0)  = y[i];      // To avoid repeated telescoping
          i_pol[0] = y_pol[i];

          // y_error
          if( do_e )
            {
              apply_y_unit2( MatrixView(VectorView( y_error[i] )), yv, 
                             y_unit, y_f[i], i_pol ); 
            }

          // Jacobian
          if( do_j )
            {
              apply_y_unit2( MatrixView( jacobian(i,joker) ), yv, 
                             y_unit, y_f[i], i_pol ); 
            }

          // y (must be done last)
          apply_y_unit( yv, y_unit, y_f[i], i_pol );
          y[i] = yv(0,0);
        }
    }
}
