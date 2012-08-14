/* Copyright (C) 2002-2012
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
#include "geodetic.h"
#include "jacobian.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "montecarlo.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"

extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT;
extern const String ABSSPECIES_MAINTAG;
extern const String TEMPERATURE_MAINTAG;





/*===========================================================================
  === Help sub-functions to handle analytical jacobians (in alphabetical order)
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
// Kept as a local function as long as just used here.
// We trust here gridpos, and "extrapolation points" are identified simply
// by a fd outside [0,1].
void add_extrap( ArrayOfGridPos&   gp )
{
  for( Index i=0; i<gp.nelem(); i++ )
    { 
      if( gp[i].fd[0] < 0 ) 
        {  
          gp[i].fd[0] = 0; 
          gp[i].fd[1] = 1; 
        }
      else if( gp[i].fd[0] > 1 ) 
        {  
          gp[i].fd[0] = 1; 
          gp[i].fd[1] = 0; 
        }
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
      add_extrap( gp_p );

      // Latitude
      Index            nr2 = 1;
      ArrayOfGridPos   gp_lat;
      if( atmosphere_dim > 1 )
        {
          gp_lat.resize(ppath.np);
          r_grid = jacobian_quantity.Grids()[1];
          nr2    = r_grid.nelem();
          gridpos( gp_lat, r_grid, ppath.pos(joker,1), extpolfac );
          add_extrap( gp_lat );
        }

      // Longitude
      ArrayOfGridPos   gp_lon;
      if( atmosphere_dim > 2 )
        {
          gp_lon.resize(ppath.np);
          r_grid = jacobian_quantity.Grids()[2];
          gridpos( gp_lon, r_grid, ppath.pos(joker,2), extpolfac );
          add_extrap( gp_lon );
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
  === Workspace methods 
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void iyApplyYunit(
         Matrix&         iy,
   ArrayOfTensor3&       iy_aux,
   const Index&          stokes_dim,
   const Vector&         f_grid,
   const Index&          jacobian_do,
   const ArrayOfString&  iy_aux_vars,
   const String&         y_unit,
   const Verbosity&)
{
  if( jacobian_do )
    throw runtime_error( "This method does not handle conversion of jacobian "
                     "quantities and should not be used with jacobian_do=1." );

  if( y_unit == "1" )
    { throw runtime_error(
        "No need to use this method with *y_unit* = \"1\"." ); }

  // Polarisation index variable
  ArrayOfIndex i_pol(stokes_dim);
  for( Index is=0; is<stokes_dim; is++ )
    { i_pol[is] = is + 1; }

  apply_y_unit( iy, y_unit, f_grid, i_pol );
  
  for( Index i=0; i<iy_aux_vars.nelem(); i++ )
    {
      if( iy_aux_vars[i] == "Radiance" )
        {
          for( Index j=0; j<iy_aux[i].ncols(); j++ )
            { apply_y_unit( iy_aux[i](joker,joker,j), y_unit, f_grid, i_pol ); }
        }
    }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void iyCalc(
         Workspace&   ws,
         Matrix&      iy,
   ArrayOfTensor3&    iy_aux,
         Ppath&       ppath,
   ArrayOfTensor3&    diy_dx,
   const Index&       basics_checked,
   const ArrayOfString&  iy_aux_vars,
   const Tensor3&     t_field,
   const Tensor3&     z_field,
   const Tensor4&     vmr_field,
   const Index&       cloudbox_on,
   const Index&       cloudbox_checked,
   const Vector&      rte_pos,
   const Vector&      rte_los,
   const Index&       jacobian_do,
   const Index&       mblock_index,
   const Agenda&      iy_main_agenda,
   const Verbosity& )
{
  // Basics and cloudbox OK?
  //
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );

  // iy_transmission is just input and can be left empty for first call
  Tensor3   iy_transmission(0,0,0);


  iy_main_agendaExecute( ws, iy, iy_aux, ppath, diy_dx, 1, iy_transmission, 
                             iy_aux_vars, cloudbox_on, jacobian_do, t_field, 
                             z_field, vmr_field, mblock_index, rte_pos, rte_los,
                             iy_main_agenda );
}





/* Workspace method: Doxygen documentation will be auto-generated */
void iyEmissionStandard(
         Workspace&                  ws,
         Matrix&                     iy,
         ArrayOfTensor3&             iy_aux,
         Ppath&                      ppath,
         ArrayOfTensor3&             diy_dx,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Tensor3&                    z_field,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const Tensor3&                    wind_u_field,
   const Tensor3&                    wind_v_field,
   const Tensor3&                    wind_w_field,
   const Tensor3&                    mag_u_field,
   const Tensor3&                    mag_v_field,
   const Tensor3&                    mag_w_field,
   const Tensor3&                    edensity_field,
   const Index&                      cloudbox_on,
   const ArrayOfString&              iy_aux_vars,
   const Index&                      jacobian_do,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const Agenda&                     ppath_agenda,
   const Agenda&                     blackbody_radiation_agenda,
   const Agenda&                     abs_mat_per_species_agenda,
   const Agenda&                     iy_main_agenda,
   const Agenda&                     iy_space_agenda,
   const Agenda&                     iy_surface_agenda,
   const Agenda&                     iy_cloudbox_agenda,
   const Index&                      iy_agenda_call1,
   const Tensor3&                    iy_transmission,
   const Index&                      mblock_index,
   const Vector&                     rte_pos,      
   const Vector&                     rte_los,      
   const Verbosity&                  verbosity )
{
  // Determine propagation path
  //
  ppath_agendaExecute( ws, ppath, rte_pos, rte_los, cloudbox_on, 0, 
                       mblock_index, t_field, z_field, vmr_field, 
                       edensity_field, -1, ppath_agenda );

  // Some basic sizes
  //
  const Index nf = f_grid.nelem();
  const Index np = ppath.np;
  const Index nq = jacobian_quantities.nelem();
  
  // Initialise iy_aux --------------------------------------------------------
  //
  Index aux_temperature = -1,
        aux_abs_sum     = -1,
        aux_background  = -1,
        aux_radiance    = -1,
        aux_trans_total = -1;
  ArrayOfIndex aux_abs_species(0), aux_abs_isp(0);
  ArrayOfIndex aux_vmr_species(0), aux_vmr_isp(0);
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
          if( iy_aux_vars[i] == "Temperature" )
            { aux_temperature = i;    iy_aux[i].resize( 1, 1, np ); }
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
              for( Index j=0; j<aux_vmr_isp.nelem(); j++ )
                {
                  if( ispecies ==  aux_vmr_isp[j] )
                    {
                      ostringstream os;
                      os << "You have selected VMR of species with index "
                         << ispecies << " twice!";
                      throw runtime_error( os.str() );
                    }
                }
              aux_vmr_species.push_back(i);
              aux_vmr_isp.push_back(ispecies);
              iy_aux[i].resize( 1, 1, np );               
            }
          else if( iy_aux_vars[i] == "Absorption, summed" )
            { aux_abs_sum = i;       iy_aux[i].resize( nf, stokes_dim, np ); }
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
              for( Index j=0; j<aux_abs_isp.nelem(); j++ )
                {
                  if( ispecies ==  aux_abs_isp[j] )
                    {
                      ostringstream os;
                      os << "You have selected absorption species with index "
                         << ispecies << " twice!";
                      throw runtime_error( os.str() );
                    }
                }
              aux_abs_species.push_back(i);
              aux_abs_isp.push_back(ispecies);
              iy_aux[i].resize( nf, stokes_dim, np );               
            }
          else if( iy_aux_vars[i] == "Radiative background" )
            { aux_background = i;     iy_aux[i].resize( 1, 1, 1 ); }
          else if( iy_aux_vars[i] == "Radiance" )
            { aux_radiance = i;       iy_aux[i].resize( nf, stokes_dim, np ); }
          else if( iy_aux_vars[i] == "Transmission" )
            { aux_trans_total = i;    iy_aux[i].resize( nf, 1, 1 ); }
          else
            {
              ostringstream os;
              os << "In *iy_aux_vars* you have included: \"" << iy_aux_vars[i]
                 << "\"\nThis choice is not recognised.";
              throw runtime_error( os.str() );
            }
        }
    }
  //---------------------------------------------------------------------------

  // Initialise analytical jacobians (diy_dx)
  //
  Index j_analytical_do = 0;
  //
  if( jacobian_do ) { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }
  //
  if( !j_analytical_do )
    { diy_dx.resize( 0 ); }
  else if( iy_agenda_call1 )
    {
      diy_dx.resize( nq ); 
      //
      FOR_ANALYTICAL_JACOBIANS_DO( 
        diy_dx[iq].resize( jacobian_indices[iq][1]-jacobian_indices[iq][0]+1,
                           nf, stokes_dim ); 
        diy_dx[iq] = 0.0;
      )
    }

  // Get atmospheric and RT quantities for each ppath point/step
  //
  // If np = 1, we only need to determine the radiative background
  //
  // "atmvars"
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w,
                               ppath_mag_u,  ppath_mag_v,  ppath_mag_w;
  Matrix    ppath_vmr;
  // "rtvars"
  Tensor3    total_tau;
  Matrix    ppath_emission;
  Tensor4   ppath_tau;
  Tensor5   ppath_abs;
  Tensor3   iy_trans_new;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      get_ppath_atmvars( ppath_p, ppath_t, ppath_vmr,
                         ppath_wind_u, ppath_wind_v, ppath_wind_w,
                         ppath_mag_u,  ppath_mag_v,  ppath_mag_w,
                         ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                         wind_u_field, wind_v_field, wind_w_field ,
                         mag_u_field, mag_v_field, mag_w_field );
      
      // Get emission, absorption and optical thickness for each step
      get_ppath_rtvars( ws, ppath_abs, ppath_tau, total_tau,
                        ppath_emission,  abs_mat_per_species_agenda, blackbody_radiation_agenda,
                        ppath, ppath_p, ppath_t, ppath_vmr, ppath_wind_u,
                        ppath_wind_v, ppath_wind_w, ppath_mag_u, ppath_mag_v,
                        ppath_mag_w, 
                        f_grid, -1, stokes_dim, atmosphere_dim, 1 );

    }
  else // To handle cases inside the cloudbox, or totally outside the atmosphere
    { 
      total_tau.resize( nf, stokes_dim, stokes_dim );
      total_tau = 0;
    }

  // iy_transmission
  //
  if( iy_agenda_call1 )
    { iy_transmission_for_tensor3_tau( iy_trans_new, stokes_dim, total_tau ); }
  else
    { iy_transmission_mult_tensor3_tau( iy_trans_new, iy_transmission, total_tau ); }

  // Radiative background
  //
  get_iy_of_background( ws, iy, diy_dx, 
                        iy_trans_new, jacobian_do, ppath, atmosphere_dim, 
                        t_field, z_field, vmr_field, cloudbox_on, 
                        stokes_dim, f_grid, iy_main_agenda, iy_space_agenda,
                        iy_surface_agenda, iy_cloudbox_agenda, verbosity);

  // Fill some parts of iy_aux
  // Radiative background
  if( aux_background >= 0 ) 
    { iy_aux[aux_background](0,0,0) = (Numeric)min((Index)2,
                                                   ppath_what_background(ppath)-1); }
  // Radiance 
  if( aux_radiance >= 0 ) 
    { iy_aux[aux_radiance](joker,joker,np-1) = iy; }
  // Transmission, total
  if( aux_trans_total >= 0 ) 
    { iy_aux[aux_trans_total](joker,0,0) = iy_trans_new(joker,0,0); }


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
            Tensor5   ppath_as2_tensor5;
      //
      if( j_analytical_do )
        { 
          diy_dpath.resize( nq ); 
          abs_species_i.resize( nq ); 
          is_t.resize( nq ); 
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
          Index do_t = 0;
          for( Index iq=0; iq<is_t.nelem() && do_t==0; iq++ )
            { if( is_t[iq] ) { do_t = 1; } }
          if( do_t )
            {
              Tensor4 tau_dummy;
              Tensor3 total_tau_dummy;
              Vector t2 = ppath_t;
              t2 += dt;
              get_ppath_rtvars( ws, ppath_as2_tensor5, tau_dummy, total_tau_dummy,
                                ppath_e2, abs_mat_per_species_agenda, blackbody_radiation_agenda,
                                ppath, ppath_p, t2, ppath_vmr, ppath_wind_u,
                                ppath_wind_v, ppath_wind_w, ppath_mag_u, ppath_mag_v,
                                ppath_mag_w, 
                                f_grid, -1, stokes_dim, atmosphere_dim, 1 );
            }
        }

      Tensor3View ppath_abs_scalar = ppath_abs(joker, joker, 0, 0, joker);
      // iy_aux for point np-1:
      // Temperature
      if( aux_temperature >= 0 ) 
        { iy_aux[aux_temperature](0,0,np-1) = ppath_t[np-1]; }
      // VMR
      for( Index j=0; j<aux_vmr_species.nelem(); j++ )
        { iy_aux[aux_vmr_species[j]](0,0,np-1) = 
                                              ppath_vmr(aux_vmr_isp[j],np-1); }
      // Absorption (change to store "absorption vector")
      if( aux_abs_sum >= 0 ) 
        { for( Index iv=0; iv<nf; iv++ ) {
            iy_aux[aux_abs_sum](iv,0,np-1) = 
                                     ppath_abs_scalar(joker, iv,np-1).sum(); } }
      for( Index j=0; j<aux_abs_species.nelem(); j++ )
        { for( Index iv=0; iv<nf; iv++ ) {
            iy_aux[aux_abs_species[j]](iv,0,np-1) = 
                                  ppath_abs_scalar(aux_abs_isp[j],iv,np-1); } }
      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          Vector esource(nf,0.);  // B-bar

          // Spectrum at end of ppath step )
          for( Index iv=0; iv<nf; iv++ )  
            {
                esource[iv] = 0.5 * ( ppath_emission(iv,ip) +
                                    ppath_emission(iv,ip+1) );
                // compare rte_step_std ... this function is not used to
                // avoid calling the LU-functions, which are unstable in
                // this case
                if( stokes_dim == 1 )
                    iy(iv,0)  = iy(iv,0) * exp(-ppath_tau(iv,0,0,ip)) +
                            esource[iv] * (1 - exp(-ppath_tau(iv,0,0,ip)));
                else
                {
                    Matrix trans_mat(stokes_dim, stokes_dim), temp = ppath_tau(iv,joker,joker,ip);
                        temp *= -1;

                    Vector planck(stokes_dim, 0.), term1(stokes_dim), term2(stokes_dim);
                        planck[0] = esource[iv]; //K⁻¹(a*B+S) = [B,0,0,0] when a from K and S = 0

                    matrix_exp( trans_mat, temp, 10);

                    id_mat(temp);
                        temp -= trans_mat;
                    mult(term1, trans_mat, iy(iv,joker) );
                    mult(term2, temp, planck);

                    iy(iv,joker)  = term1;
                    iy(iv,joker) += term2;
                }
            }

          //Debugging variables
          // Jacobian quantities
          if( j_analytical_do )
            {
              VectorView total_tau_vector = total_tau(joker, 0, 0);
              MatrixView ppath_tau_matrix = ppath_tau(joker,0,0,joker);
              // Common terms introduced for efficiency and clarity
              Vector X(nf), Y(nf);   // See AUG
              //
              for( Index iv=0; iv<nf; iv++ )
                {
                  X[iv] = 0.5 * ppath.lstep[ip] * exp( -total_tau_vector[iv] );
                  Y[iv] = X[iv] * ( esource[iv] - iy(iv,0) );
                }

              // Loop quantities
              for( Index iq=0; iq<nq; iq++ ) 
                {
                  if( jacobian_quantities[iq].Analytical() )
                    {
                      // Absorbing species
                      const Index isp = abs_species_i[iq]; 
                      if( isp >= 0 )
                        {
                          // Scaling factors to handle retrieval unit
                          // (gives the cross-section to apply)
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
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // Stokes component 1
                              diy_dpath[iq](ip  ,iv,0) += Y[iv] * unitscf1 *
                                                 ppath_abs_scalar(isp,iv,ip);
                              diy_dpath[iq](ip+1,iv,0) += Y[iv] * unitscf2 * 
                                                 ppath_abs_scalar(isp,iv,ip+1);
                              // Higher stokes components
                              for( Index is=1; is<stokes_dim; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is); 
                                  diy_dpath[iq](ip  ,iv,is) += Z * unitscf1 *
                                                 ppath_abs_scalar(isp,iv,ip);
                                  diy_dpath[iq](ip+1,iv,is) += Z * unitscf2 *
                                                 ppath_abs_scalar(isp,iv,ip+1);
                                }
                            }
                        }

                      // Temperature
                      else if( is_t[iq] )
                        {
                          Tensor3View ppath_as2 = ppath_as2_tensor5(joker, joker, 0, 0, joker);
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // The terms associated with Dtau/Dk:
                              const Numeric k1 = 
                                         ppath_abs_scalar(joker,iv,ip  ).sum();
                              const Numeric k2 = 
                                         ppath_abs_scalar(joker,iv,ip+1).sum();
                              const Numeric dkdt1 =
                                  ( ppath_as2(joker,iv,ip  ).sum() - k1 ) / dt;
                              const Numeric dkdt2 =
                                  ( ppath_as2(joker,iv,ip+1).sum() - k2 ) / dt;
                              // Stokes 1:
                              diy_dpath[iq](ip  ,iv,0) += Y[iv] * dkdt1;
                              diy_dpath[iq](ip+1,iv,0) += Y[iv] * dkdt2;
                              // Higher Stokes
                              for( Index is=1; is<stokes_dim; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is);
                                  diy_dpath[iq](ip  ,iv,is) += Z * dkdt1;
                                  diy_dpath[iq](ip+1,iv,is) += Z * dkdt2;
                                }
                              //
                              // The terms associated with B-bar:
                              const Numeric V = 0.5 * 
                                      exp(-(total_tau_vector[iv]-ppath_tau_matrix(iv,ip))) *
                                              ( 1.0 - exp(-ppath_tau_matrix(iv,ip)) );
                              diy_dpath[iq](ip  ,iv,0) += V *
                                              ( ppath_e2(iv,ip) -
                                                ppath_emission(iv,ip) ) / dt;
                              diy_dpath[iq](ip+1,iv,0) += V * 
                                              ( ppath_e2(iv,ip+1) -
                                                ppath_emission(iv,ip+1) ) / dt;
                              // Zero for higher Stokes
                              //
                              // The terms associated with Delta-s:
                              if( jacobian_quantities[iq].Subtag() == "HSE on" )
                                {
                                  // Stokes 1:
                                  const Numeric kbar = 0.5 * ( k1 + k2 );
                                  diy_dpath[iq](ip  ,iv,0) += Y[iv] * kbar /
                                                                  ppath_t[ip];
                                  diy_dpath[iq](ip+1,iv,0) += Y[iv] * kbar /
                                                                 ppath_t[ip+1];
                                  // Higher Stokes
                                  for( Index is=1; is<stokes_dim; is++ )
                                    { 
                                      const Numeric Z = -X[iv] * iy(iv,is);
                                      diy_dpath[iq](ip  ,iv,is) += Z * kbar /
                                                                  ppath_t[ip];
                                      diy_dpath[iq](ip+1,iv,is) += Z * kbar /
                                                                 ppath_t[ip+1];
                                    }
                                } //hse
                            }  // frequency
                        }  // if is_t
                    } // if analytical
                } // for iq
              
              // Update total_tau (not used for spectrum calculation)
              for( Index iv=0; iv<nf; iv++ )
                { total_tau_vector[iv] -= ppath_tau_matrix(iv,ip); }
            }

          // --- iy_aux: ------------------------------------------------------
          // Temperature
// FIXME: Richard: iy calculated directly at l. 1147
/*          if( aux_temperature >= 0 ) 
            { iy_aux[aux_temperature](0,0,ip) = ppath_t[ip]; }
          // VMR
          for( Index j=0; j<aux_vmr_species.nelem(); j++ )
            { iy_aux[aux_vmr_species[j]](0,0,ip) = 
                                                ppath_vmr(aux_vmr_isp[j],ip); }
          // Absorption (change to store "absorption vector")
          if( aux_abs_sum >= 0 ) 
            { for( Index iv=0; iv<nf; iv++ ) {
              iy_aux[aux_abs_sum](iv,0,ip) = 
                                       ppath_abs_scalar(iv,joker,ip).sum(); } }
          for( Index j=0; j<aux_abs_species.nelem(); j++ )
            { for( Index iv=0; iv<nf; iv++ ) {
              iy_aux[aux_abs_species[j]](iv,0,ip) = 
                                    ppath_abs_scalar(iv,aux_abs_isp[j],ip); } }
          // Radiance 
          if( aux_radiance >= 0 ) 
            { iy_aux[aux_radiance](joker,joker,ip) = iy_new; }
          // --- iy_aux: ------------------------------------------------------

          // Update iy
          iy = iy_new;
*/        } 

      // Map jacobians from ppath to retrieval grids
      // (this operation corresponds to the term Dx_i/Dx)
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
    } // if np>1
}



/* Workspace method: Doxygen documentation will be auto-generated */
/*
void iyEmissionStandardClearsky(
         Workspace&                  ws,
         Matrix&                     iy,
         ArrayOfTensor3&             iy_aux,
         ArrayOfTensor3&             diy_dx,
   const Index&                      iy_agenda_call1,
   const Tensor3&                    iy_transmission,
   const Vector&                     rte_pos,      
   const Vector&                     rte_los,      
   const Index&                      jacobian_do,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Tensor3&                    z_field,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const Tensor3&                    wind_u_field,
   const Tensor3&                    wind_v_field,
   const Tensor3&                    wind_w_field,
   const Tensor3&                    edensity_field,
   const Index&                      cloudbox_on,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const Index&                      mblock_index,
   const Agenda&                     ppath_agenda,
   const Agenda&                     blackbody_radiation_agenda,
   const Agenda&                     abs_scalar_gas_agenda,
   const Agenda&                     iy_main_agenda,
   const Agenda&                     iy_space_agenda,
   const Agenda&                     iy_surface_agenda,
   const Agenda&                     iy_cloudbox_agenda,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const Verbosity&                  verbosity )
{
  iy_aux.resize(0);

  // The method can in principle be used "stand-alone", but for efficiency
  // reasons we skip all checks needed to handle such usage. Those checks are
  // done by yCalc.

  // iy_aux and diy_dx are initiated by iyb_calc
  // to be empty/zero. Accordingly, if any radiative background wants
  // to flag an error, iy_error must be resized. No subsequent scaling of the
  // error is made. Thus the iy_transmission should be considered when adding
  // a term to the error. A setting of iy_aux must also include a resizing. 

  // The WSV *iy_transmission* holds the transmission between the sensor
  // and the end of the propagation path.

  // Sizes
  //
  const Index nf = f_grid.nelem();
  const Index nq = jacobian_quantities.nelem();

  // Determine if there are any analytical jacobians to handle, and if primary
  // call, resize diy_dx.
  //
  Index j_analytical_do = 0;
  //
  if( jacobian_do ) { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }
  //
  if( iy_agenda_call1 && j_analytical_do )
    {
      diy_dx.resize( nq ); 
      //
      FOR_ANALYTICAL_JACOBIANS_DO( 
        diy_dx[iq].resize( jacobian_indices[iq][1] - jacobian_indices[iq][0] +
                           1, nf, stokes_dim ); 
        diy_dx[iq] = 0.0;
      )
    }

  //- Determine propagation path
  //
  Ppath  ppath;
  //
  ppath_agendaExecute( ws, ppath, rte_pos, rte_los, cloudbox_on, 0, 
                       mblock_index, t_field, z_field, vmr_field, 
                       edensity_field, -1, ppath_agenda );

  // Get atmospheric and RT quantities for each ppath point/step
  //
  // If np = 1, we only need to determine the radiative background
  //
  // "atmvars"
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w,
                               ppath_mag_u,  ppath_mag_v,  ppath_mag_w;
  Matrix    ppath_vmr;
  // "rtvars"
  Tensor3   total_tau;
  Matrix    ppath_emission;
  Tensor4   ppath_tau;
  Tensor5   ppath_abs;
  Tensor3   iy_trans_new;
  //
  const Index np  = ppath.np;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      get_ppath_atmvars( ppath_p, ppath_t, ppath_vmr,
                         ppath_wind_u, ppath_wind_v, ppath_wind_w,
                         ppath_mag_u,  ppath_mag_v,  ppath_mag_w,
                         ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                         wind_u_field, wind_v_field, wind_w_field ,
                         mag_u_field, mag_v_field, mag_w_field );

      // Get emission, absorption and optical thickness for each step
      get_ppath_rtvars( ws, ppath_abs, ppath_tau, total_tau,
                        ppath_emission, abs_mat_per_species_agenda, blackbody_radiation_agenda,
                        ppath, ppath_p, ppath_t, ppath_vmr, ppath_wind_u,
                        ppath_wind_v, ppath_wind_w, ppath_mag_u, ppath_mag_v,
                        ppath_mag_w, 
                        f_grid, -1, stokes_dim, atmosphere_dim, 1 );
    }
  else // To handle cases inside the cloudbox, or totally outside the atmosphere
    { 
      total_tau.resize( nf );
      total_tau = 0;
    }

  // iy_transmission
  //
  if( iy_agenda_call1 )
    { iy_transmission_for_scalar_tau( iy_trans_new, stokes_dim, total_tau ); }
  else
    { iy_transmission_mult_scalar_tau( iy_trans_new, iy_transmission, 
                                                                total_tau ); }

  // Radiative background
  //
  get_iy_of_background( ws, iy, diy_dx, 
                        iy_trans_new, jacobian_do, ppath, atmosphere_dim, 
                        t_field, z_field, vmr_field, cloudbox_on, 
                        stokes_dim, f_grid, iy_main_agenda, iy_space_agenda,
                        iy_surface_agenda, iy_cloudbox_agenda, verbosity);
  
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
          diy_dpath.resize( nq ); 
          abs_species_i.resize( nq ); 
          is_t.resize( nq ); 
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
          Index do_t = 0;
          for( Index iq=0; iq<is_t.nelem() && do_t==0; iq++ )
            { if( is_t[iq] ) { do_t = 1; } }
          if( do_t )
            {
              Matrix tau_dummy; 
              Vector total_tau_dummy, t2 = ppath_t;
              t2 += dt;
              get_ppath_rtvars( ws, ppath_as2, tau_dummy, total_tau_dummy, 
                                ppath_e2, abs_scalar_gas_agenda, 
                                blackbody_radiation_agenda, ppath, ppath_p, 
                                t2, ppath_vmr, 
                                ppath_wind_u, ppath_wind_v, ppath_wind_w, 
                                f_grid, -1, atmosphere_dim, 1 );
            }
        }

      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          Vector esource(nf);            // B-bar
          Matrix iy_new(nf,stokes_dim);  // Intensity at end of step

          // Spectrum at end of ppath step 
          // (iy updated to iy_new at end of ip-loop)
          for( Index iv=0; iv<nf; iv++ )  
            {

                if( stokes_dim == 1 )
                    iy(iv,0)  = iy(iv,0) * exp(-ppath_tau(iv,0,0,ip)) +
                    0.5 * ( ppath_emission(iv,ip) +  ppath_emission(iv,ip+1) ) *
                    (1 - exp(-ppath_tau(iv,0,0,ip)));
                else
                {
                    Matrix trans_mat(stokes_dim, stokes_dim), temp = ppath_tau(iv,joker,joker,ip);
                        temp *= -1;

                    Vector planck(stokes_dim, 0.), term1(stokes_dim), term2(stokes_dim);
                        planck[0] = 0.5 * ( ppath_emission(iv,ip) +
                                    ppath_emission(iv,ip+1) ); //K⁻¹(a*B+S) = [B,0,0,0] when a from K and S = 0

                    matrix_exp( trans_mat, temp, 10);

                    id_mat(temp);
                        temp -= trans_mat;

                    mult(term1, trans_mat, iy(iv,joker) );
                    mult(term2, temp, planck);

                    iy(iv,joker)  = term1;
                    iy(iv,joker) += term2;
                }

            }

          // Jacobian quantities
          if( j_analytical_do )
            {
              // Common terms introduced for efficiency and clarity
              Vector X(nf), Y(nf);   // See AUG
              //
              for( Index iv=0; iv<nf; iv++ )
                {
                  X[iv] = 0.5 * ppath.lstep[ip] * exp( -total_tau[iv] );
                  Y[iv] = X[iv] * ( esource[iv] - iy(iv,0) );
                }

              // Loop quantities
              for( Index iq=0; iq<nq; iq++ ) 
                {
                  if( jacobian_quantities[iq].Analytical() )
                    {
                      // Absorbing species
                      const Index isp = abs_species_i[iq]; 
                      if( isp >= 0 )
                        {
                          // Scaling factors to handle retrieval unit
                          // (gives the cross-section to apply)
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
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // Stokes component 1
                              diy_dpath[iq](ip  ,iv,0) += Y[iv] * unitscf1 *
                                                 ppath_abs_scalar(iv,isp,ip);
                              diy_dpath[iq](ip+1,iv,0) += Y[iv] * unitscf2 * 
                                                 ppath_abs_scalar(iv,isp,ip+1);
                              // Higher stokes components
                              for( Index is=1; is<stokes_dim; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is); 
                                  diy_dpath[iq](ip  ,iv,is) += Z * unitscf1 *
                                                 ppath_abs_scalar(iv,isp,ip);
                                  diy_dpath[iq](ip+1,iv,is) += Z * unitscf2 *
                                                 ppath_abs_scalar(iv,isp,ip+1);
                                }
                            }
                        }

                      // Temperature
                      else if( is_t[iq] )
                        {
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // The terms associated with Dtau/Dk:
                              const Numeric k1 = 
                                         ppath_abs_scalar(iv,joker,ip  ).sum();
                              const Numeric k2 = 
                                         ppath_abs_scalar(iv,joker,ip+1).sum();
                              const Numeric dkdt1 =
                                  ( ppath_as2(iv,joker,ip  ).sum() - k1 ) / dt;
                              const Numeric dkdt2 =
                                  ( ppath_as2(iv,joker,ip+1).sum() - k2 ) / dt;
                              // Stokes 1:
                              diy_dpath[iq](ip  ,iv,0) += Y[iv] * dkdt1;
                              diy_dpath[iq](ip+1,iv,0) += Y[iv] * dkdt2;
                              // Higher Stokes
                              for( Index is=1; is<stokes_dim; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is);
                                  diy_dpath[iq](ip  ,iv,is) += Z * dkdt1;
                                  diy_dpath[iq](ip+1,iv,is) += Z * dkdt2;
                                }
                              //
                              // The terms associated with B-bar:
                              const Numeric V = 0.5 * 
                                      exp(-(total_tau[iv]-ppath_tau(iv,ip))) *
                                              ( 1.0 - exp(-ppath_tau(iv,ip)) );
                              diy_dpath[iq](ip  ,iv,0) += V *
                                              ( ppath_e2(iv,ip) -
                                                ppath_emission(iv,ip) ) / dt;
                              diy_dpath[iq](ip+1,iv,0) += V * 
                                              ( ppath_e2(iv,ip+1) -
                                                ppath_emission(iv,ip+1) ) / dt;
                              // Zero for higher Stokes
                              //
                              // The terms associated with Delta-s:
                              if( jacobian_quantities[iq].Subtag() == "HSE on" )
                                {
                                  // Stokes 1:
                                  const Numeric kbar = 0.5 * ( k1 + k2 );
                                  diy_dpath[iq](ip  ,iv,0) += Y[iv] * kbar /
                                                                  ppath_t[ip];
                                  diy_dpath[iq](ip+1,iv,0) += Y[iv] * kbar /
                                                                 ppath_t[ip+1];
                                  // Higher Stokes
                                  for( Index is=1; is<stokes_dim; is++ )
                                    { 
                                      const Numeric Z = -X[iv] * iy(iv,is);
                                      diy_dpath[iq](ip  ,iv,is) += Z * kbar /
                                                                  ppath_t[ip];
                                      diy_dpath[iq](ip+1,iv,is) += Z * kbar /
                                                                 ppath_t[ip+1];
                                    }
                                } //hse
                            }  // frequency
                        }  // if is_t
                    } // if analytical
                } // for iq
              
              // Update total_tau (not used for spectrum calculation)
              for( Index iv=0; iv<nf; iv++ )
                { total_tau[iv] -= ppath_tau(iv,ip); }
            }

          // Update iy
          iy = iy_new;
        } 

      // Map jacobians from ppath to retrieval grids
      // (this operation corresponds to the term Dx_i/Dx)
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
    } // if np>1
}
*/




/* Workspace method: Doxygen documentation will be auto-generated */
void iyMC(
         Workspace&                  ws,
         Matrix&                     iy,
         ArrayOfTensor3&             iy_aux,
         ArrayOfTensor3&             diy_dx,
   const Index&                      iy_agenda_call1,
   const Tensor3&                    iy_transmission,
   const Vector&                     rte_pos,      
   const Vector&                     rte_los,      
   const ArrayOfString&              iy_aux_vars,
   const Index&                      jacobian_do,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
   const Tensor3&                    z_field,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const Vector&                     refellipsoid,
   const Matrix&                     z_surface,
   const Index&                      cloudbox_on,
   const ArrayOfIndex&               cloudbox_limits,
   const Index&                      cloudbox_checked,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const ArrayOfSingleScatteringData&   scat_data_raw,
   const Agenda&                     iy_space_agenda,
   const Agenda&                     surface_rtprop_agenda,
   const Agenda&                     abs_mat_per_species_agenda, 
   const Tensor4&                    pnd_field,
   const String&                     y_unit,
   const Numeric&                    mc_std_err,
   const Index&                      mc_max_time,
   const Index&                      mc_max_iter,
   const Verbosity&                  verbosity)
{
  iy_aux.resize(iy_aux_vars.nelem()); // !!!

  // Throw error if unsupported features are requested
  if( atmosphere_dim != 3 )
    throw runtime_error( 
                "Only 3D atmospheres are allowed (atmosphere_dim must be 3)" );
  if( !cloudbox_on )
    throw runtime_error( 
                    "The cloudbox must be activated (cloudbox_on must be 1)" );
  if( jacobian_do )
    throw runtime_error( 
        "This method does not provide any jacobians (jacobian_do must be 0)" );
  if( !iy_agenda_call1 )
    throw runtime_error( 
                  "Recursive usage not possible (iy_agenda_call1 must be 1)" );
  if( iy_transmission.ncols() )
    throw runtime_error( "*iy_transmission* must be empty" );

  // See initial comments of iyEmissionStandardClearsky

  // Size output variables
  //
  const Index   nf  = f_grid.nelem();
  //
  iy.resize( nf, stokes_dim );
  //
  iy_aux.resize(0);
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
                                          f_grid, f_index, verbosity );

      // Seed reset for each loop. If not done, the errors 
      // appear to be highly correlated.
      MCSetSeedFromTime( mc_seed, verbosity );

      Vector y, mc_error;
                  
      MCGeneral( ws, y, mc_iteration_count, mc_error, mc_points, mc_antenna, 
                 f_grid, f_index, pos, los, stokes_dim, atmosphere_dim,
                 iy_space_agenda, surface_rtprop_agenda,
                 abs_mat_per_species_agenda, p_grid, lat_grid, lon_grid, z_field, 
                 refellipsoid, z_surface, t_field, vmr_field, 
                 cloudbox_on, cloudbox_limits, 
                 pnd_field, scat_data_mono, 1, cloudbox_checked,
                 mc_seed, y_unit, mc_std_err, mc_max_time, mc_max_iter,
                 verbosity); 
                 // GH 2011-06-17, mc_z_field_is_1D);

      assert( y.nelem() == stokes_dim );

      // Data returned can not be in Tb
      if ( y_unit == "RJBT" )
        { 
          const Numeric scfac = invrayjean( 1, f_grid[f_index] );
          y /= scfac;
          mc_error /= scfac;
        }
     
      iy(f_index,joker) = y;
      //iy_error(f_index,joker) = mc_error;
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void iyRadioLink(
         Workspace&            ws,
         Matrix&               iy,
         ArrayOfTensor3&       iy_aux,
         Ppath&                ppath,
         ArrayOfTensor3&       diy_dx,
   const Index&                iy_agenda_call1,
   const Tensor3&              ,//iy_transmission
   const Vector&               rte_pos,      
   const ArrayOfString&        iy_aux_vars,
   const Index&                jacobian_do,
   const Index&                atmosphere_dim,
   const Vector&               p_grid,
   const Vector&               lat_grid,
   const Vector&               lon_grid,
   const Tensor3&              z_field,
   const Tensor3&              t_field,
   const Tensor4&              vmr_field,
   const Tensor3&              wind_u_field,
   const Tensor3&              wind_v_field,
   const Tensor3&              wind_w_field,
   const Tensor3&              mag_u_field,
   const Tensor3&              mag_v_field,
   const Tensor3&              mag_w_field,
   const Tensor3&              edensity_field,
   const Vector&               refellipsoid,
   const Matrix&               z_surface,
   const Index&                cloudbox_on,
   const Index&                stokes_dim,
   const Vector&               f_grid,
   const Index&                dispersion_do,
   const Index&                mblock_index,
   const Agenda&               ppath_agenda,
   const Agenda&               ppath_step_agenda,
   const Agenda&               abs_mat_per_species_agenda,
   const Agenda&               iy_transmitter_agenda,
   const Verbosity&            verbosity )
{
  // Method used correctly? 
  if( !iy_agenda_call1 )
    throw runtime_error( 
                 "Recursive usage not possible (iy_agenda_call1 must be 1)." );
  //
  if( jacobian_do )
    throw runtime_error( "This method does not yet provide any jacobians and "
                         "*jacobian_do* must be 0." );
  diy_dx.resize(0);

  // Some basic sizes
  //
  Index nf = f_grid.nelem();

  // Initialise iy_aux --------------------------------------------------------
  //
  Index aux_FreeSpaceLoss = -1,
        aux_AtmosphericLoss = -1,
        aux_DefocusingLoss = -1,
        aux_TotalLoss = -1,
        aux_ExtraPathDelay = -1,
        aux_BendingAngle = -1;
  //
  const Index naux = iy_aux_vars.nelem();
  iy_aux.resize( naux );
  //
  for( Index i=0; i<naux; i++ )
    {
      if( iy_aux_vars[i] == "Free space loss" )
        { aux_FreeSpaceLoss = i;     iy_aux[i].resize( nf, 1, 1 ); }
      else if( iy_aux_vars[i] == "Atmospheric loss" )
        { aux_AtmosphericLoss = i;   iy_aux[i].resize( nf, 1, 1 ); }
      else if( iy_aux_vars[i] == "Defocusing loss" )
        { aux_DefocusingLoss = i;    iy_aux[i].resize( nf, 1, 1 ); }
      else if( iy_aux_vars[i] == "Total loss" )
        { aux_TotalLoss = i;         iy_aux[i].resize( nf, 1, 1 ); }
      else if( iy_aux_vars[i] == "Extra path delay" )
        { aux_ExtraPathDelay = i;    iy_aux[i].resize( nf, 1, 1 ); }
      else if( iy_aux_vars[i] == "Bending angle" )
        { aux_BendingAngle = i;      iy_aux[i].resize( nf, 1, 1 ); }      
    }
  //---------------------------------------------------------------------------

  // Variables to handle dispersion
  Index nloops, f_index; 
  Range fr(0,0);
  if( dispersion_do )
    { 
      nloops  = nf; 
      nf      = 1;
    }  
  else
    {
      nloops  = nf; 
      fr      = Range( 0, nf );
      f_index = -1;
    }

  // Frequency loop  
  for( Index i=0; i<nloops; i++ )
    {
      if( dispersion_do )
        { 
          f_index = i; 
          fr      = Range( i, 1 );
        }

      //- Determine propagation path
      Vector rte_los(0);  // Dummy value
      //
      ppath_agendaExecute( ws, ppath, rte_pos, rte_los, cloudbox_on, 0,
                           mblock_index, t_field, z_field, vmr_field,
                           edensity_field, f_index, ppath_agenda );
      if( ppath_what_background(ppath) > 2 )
        { throw runtime_error( "Radiative background not set to \"space\" by "
                      "*ppath_agenda*. Is correct WSM used in the agenda?" ); }

     // if( AtmosphericLoss )
          // Get atmospheric and RT quantities for each ppath point/step
          // "atmvars":
          Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w,
                                       ppath_mag_u,  ppath_mag_v,  ppath_mag_w;
          Matrix    ppath_vmr, emission_dummy;
          // "rtvars":
          Tensor3   total_tau, iy_trans_new;
          Tensor4   ppath_tau;
          Tensor5   ppath_abs;
          Agenda    agenda_dummy;
      if( ppath.np > 1 )
        {
            // Get pressure, temperature and VMRs
            get_ppath_atmvars( ppath_p, ppath_t, ppath_vmr,
                        ppath_wind_u, ppath_wind_v, ppath_wind_w,
                        ppath_mag_u,  ppath_mag_v,  ppath_mag_w,
                        ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                        wind_u_field, wind_v_field, wind_w_field,
                        mag_u_field, mag_v_field, mag_w_field );
            // Absorption and optical thickness for each step
            get_ppath_rtvars( ws, ppath_abs, ppath_tau, total_tau,
                            emission_dummy,  abs_mat_per_species_agenda, agenda_dummy,
                            ppath, ppath_p, ppath_t, ppath_vmr, ppath_wind_u,
                            ppath_wind_v, ppath_wind_w, ppath_mag_u, ppath_mag_v,
                            ppath_mag_w,
                            f_grid, f_index, stokes_dim, atmosphere_dim, 0 );
        }

      // Transmitted signal
      //
      if( dispersion_do )
        { 
          Matrix iy1;  // Single frequency version of iy
          iy_transmitter_agendaExecute( ws, iy1, rte_pos, rte_los, 
                                        Vector(1,f_grid[f_index]), 
                                        iy_transmitter_agenda ); 
              
          if( iy1.ncols() != stokes_dim  ||  iy1.nrows() != 1 )
            { throw runtime_error( "The size of *iy* returned from "
                                 "*iy_transmitter_agenda* is not correct." ); }
              iy(f_index,joker) = iy1(0,joker);
        }
      else
        { 
          iy_transmitter_agendaExecute( ws, iy, rte_pos, rte_los, f_grid,
                                        iy_transmitter_agenda ); 
          if( iy.ncols() != stokes_dim  ||  iy.nrows() != nf )
            { throw runtime_error( "The size of *iy* returned from "
                                 "*iy_transmitter_agenda* is not correct." ); }
        }

      // Different ppath variables
      //
      Numeric lbg;  // Bent geometrical length of ray path
      Numeric lba;  // Bent apparent length of ray path
      //
      lbg = ppath.start_lstep + ppath.end_lstep;
      lba = lbg;

      // Loop path points
      if( ppath.np > 1 )
        {
          // Loop ppath steps
          for( Index ip=ppath.np-2; ip>=0; ip-- )
            {
              lbg += ppath.lstep[ip];
              lba += ppath.lstep[ip] * ( ppath.ngroup[ip] + 
                                         ppath.ngroup[ip+1] ) / 2.0;
              // Include atmospheric loss
              if( dispersion_do )
                {
                    Matrix step_tr(stokes_dim, stokes_dim);
                    Matrix temp = ppath_tau(0,joker,joker,ip);
                    temp *= -1;
                    matrix_exp(step_tr, temp, 10);

                    Vector iy_temp = iy(f_index,joker);
                    mult(iy(f_index,joker), step_tr, iy_temp);

                    /*for( Index is=0; is<stokes_dim; is++ )
                    { iy(f_index,is) *= step_tr; }*/
                }
              else
                {
                  // Loop frequencies
                  for( Index iv=0; iv<nf; iv++ )
                    {
                        Matrix step_tr(stokes_dim, stokes_dim);
                        Matrix temp = ppath_tau(iv,joker,joker,ip);
                        temp *= -1;
                        matrix_exp(step_tr, temp, 10);

                        Vector iy_temp = iy(iv,joker);
                        mult(iy(iv,joker), step_tr, iy_temp);
                        /*
                        const Numeric step_tr = exp( -ppath_tau(iv,ip) );
                        for( Index is=0; is<stokes_dim; is++ )
                        { iy(iv,is) *= step_tr; }*/
                    }
                }
                
            } 
        }

      // Determine free space loss
      Numeric fspl = 1 / ( 4 * PI * lbg*lbg ); 

      // Determine defocusing loss
      Numeric dfl = 1;
      if( 0 )
        { defocusing_sat2sat( ws, dfl, ppath_step_agenda, atmosphere_dim, 
                              p_grid, lat_grid, lon_grid, t_field, z_field, 
                              vmr_field, edensity_field, f_index, 
                              refellipsoid, z_surface, ppath, verbosity ); 
        }
      else
        {
          defocusing_general( ws, dfl, ppath_step_agenda, atmosphere_dim, 
                              p_grid, lat_grid, lon_grid, t_field, z_field, 
                              vmr_field, edensity_field, f_index, 
                              refellipsoid, z_surface, ppath, verbosity ); 
        }

      // Include loss the two last loss terms
      if( dispersion_do )
        { for( Index is=0; is<stokes_dim; is++ )
            { iy(f_index,is) *= fspl*dfl; } }
      else
        { iy *= fspl*dfl; }

      // Fill parts of iy_aux
      if( aux_FreeSpaceLoss >= 0 )
        { iy_aux[aux_FreeSpaceLoss](fr,0,0) = fspl; }
      if( aux_DefocusingLoss >= 0 )
        { iy_aux[aux_DefocusingLoss](fr,0,0) = dfl; }
      if( aux_AtmosphericLoss >= 0 )
        { iy_aux[aux_DefocusingLoss](fr,0,0) = dfl; }


      /*
      // Extra path delay
      Numeric epd = -999;
      if( ExtraPathDelay )
        {
          // Radius of rte_pos and rte_pos2
          const Numeric r1 = ppath.end_pos[0] +
                             pos2refell_r( atmosphere_dim, refellipsoid, 
                                           lat_grid, lon_grid, ppath.end_pos );
          const Numeric r2 = ppath.start_pos[0] +
                             pos2refell_r( atmosphere_dim, refellipsoid, 
                                         lat_grid, lon_grid, ppath.start_pos );
          // Geomtrical distance between start and end point
          Numeric lgd ;
          if( atmosphere_dim <= 2 )
            { distance2D( lgd, r1, ppath.end_pos[1], r2, ppath.start_pos[1] ); }
          else 
            { distance3D( lgd, r1, ppath.end_pos[1],   ppath.end_pos[2],
                               r2, ppath.start_pos[1], ppath.start_pos[2] ); }
          //
          epd = ( lba - lgd ) / SPEED_OF_LIGHT;
        }

      // Bending angle
      Numeric ba = -999;
      if( BendingAngle )
        { bending_angle1d( ba, ppath ); }
      */

      // Atmospheric loss as aux is a special case
      /*
      if( iy_aux_var_old == "AtmosphericLoss" )
        { iy_aux(fr,joker) = iy(fr,joker); }
      */

    }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void iyTransmissionStandard(
         Workspace&                  ws,
         Matrix&                     iy,
         ArrayOfTensor3&             iy_aux,
         Ppath&                      ppath,
         ArrayOfTensor3&             diy_dx,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Tensor3&                    z_field,
   const Tensor3&                    t_field,
   const Tensor4&                    vmr_field,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const Tensor3&                    wind_u_field,
   const Tensor3&                    wind_v_field,
   const Tensor3&                    wind_w_field,
   const Tensor3&                    mag_u_field,
   const Tensor3&                    mag_v_field,
   const Tensor3&                    mag_w_field,
   const Tensor3&                    edensity_field,
   const Index&                      cloudbox_on,
   const ArrayOfIndex&               cloudbox_limits,
   const ArrayOfString&              iy_aux_vars,
   const Index&                      jacobian_do,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const Agenda&                     ppath_agenda,
   const Agenda&                     abs_mat_per_species_agenda,
   const Agenda&                     iy_transmitter_agenda,
   const Index&                      iy_agenda_call1,
   const Tensor3&                 ,//iy_transmission,
   const Index&                      mblock_index,
   const Vector&                     rte_pos,      
   const Vector&                     rte_los,      
   const Verbosity&  )
{
  if( !iy_agenda_call1 )
    throw runtime_error( 
      "This method can not be used recursively (iy_agenda_call1 must be 1)." );

  if( cloudbox_on  || cloudbox_limits.nelem() )
    throw runtime_error( "Cloudbox not yet handled." );

  //- Determine propagation path
  //
  ppath_agendaExecute( ws, ppath, rte_pos, rte_los, 0, 0,
                       mblock_index, t_field, z_field, vmr_field, 
                       edensity_field, -1, ppath_agenda );

  // Some basic sizes
  //
  const Index nf = f_grid.nelem();
  const Index np = ppath.np;
  const Index nq = jacobian_quantities.nelem();

  // Initialise iy_aux --------------------------------------------------------
  //
  Index aux_temperature = -1,
        aux_abs_sum     = -1,
        aux_radiance    = -1;
  ArrayOfIndex aux_abs_species(0), aux_abs_isp(0);
  ArrayOfIndex aux_vmr_species(0), aux_vmr_isp(0);
  //
  const Index naux = iy_aux_vars.nelem();
  iy_aux.resize( naux );
  //
  for( Index i=0; i<naux; i++ )
    {
      if( iy_aux_vars[i] == "Temperature" )
        { aux_temperature = i;    iy_aux[i].resize( 1, 1, np ); }
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
          for( Index j=0; j<aux_vmr_isp.nelem(); j++ )
            {
              if( ispecies ==  aux_vmr_isp[j] )
                {
                  ostringstream os;
                  os << "You have selected VMR of species with index "
                     << ispecies << " twice!";
                  throw runtime_error( os.str() );
                }
            }
          aux_vmr_species.push_back(i);
          aux_vmr_isp.push_back(ispecies);
          iy_aux[i].resize( 1, 1, np );               
        }
      else if( iy_aux_vars[i] == "Absorption, summed" )
        { aux_abs_sum = i;       iy_aux[i].resize( nf, stokes_dim, np ); }
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
          for( Index j=0; j<aux_abs_isp.nelem(); j++ )
            {
              if( ispecies ==  aux_abs_isp[j] )
                {
                  ostringstream os;
                  os << "You have selected absorption species with index "
                     << ispecies << " twice!";
                  throw runtime_error( os.str() );
                }
            }
          aux_abs_species.push_back(i);
          aux_abs_isp.push_back(ispecies);
          iy_aux[i].resize( nf, stokes_dim, np );               
        }
      else if( iy_aux_vars[i] == "Radiance" )
        { aux_radiance = i;       iy_aux[i].resize( nf, stokes_dim, np ); }
      else
        {
          ostringstream os;
          os << "In *iy_aux_vars* you have included: \"" << iy_aux_vars[i]
             << "\"\nThis choice is not recognised.";
          throw runtime_error( os.str() );
        }
    }
  //---------------------------------------------------------------------------

  // Initialise analytical jacobians (diy_dx)
  //
  Index j_analytical_do = 0;
  //
  if( jacobian_do ) { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }
  //
  if( j_analytical_do )
    {
      diy_dx.resize( jacobian_indices.nelem() ); 
      //
      FOR_ANALYTICAL_JACOBIANS_DO( 
        diy_dx[iq].resize( jacobian_indices[iq][1] - jacobian_indices[iq][0] + 
                           1, nf, stokes_dim ); 
        diy_dx[iq] = 0.0;
      ) 
    }

  // Get atmospheric and RT quantities for each ppath point/step
  //
  // If np = 1, we only need to determine the radiative background
  //
  // "atmvars"
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w,
                               ppath_mag_u,  ppath_mag_v,  ppath_mag_w;
  Matrix    ppath_vmr;
  // "rtvars"
  Tensor3   total_tau;
  Matrix    emission_dummy;
  Tensor4   ppath_tau;
  Tensor5   ppath_abs;
  Tensor3   iy_trans_new;
  Agenda    agenda_dummy;
  //
  if( np > 1 )
    {
      // Get pressure, temperature and VMRs
      get_ppath_atmvars( ppath_p, ppath_t, ppath_vmr,
                         ppath_wind_u, ppath_wind_v, ppath_wind_w,
                         ppath_mag_u,  ppath_mag_v,  ppath_mag_w,
                         ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                         wind_u_field, wind_v_field, wind_w_field,
                         mag_u_field, mag_v_field, mag_w_field );

      // Absorption and optical thickness for each step
      get_ppath_rtvars( ws, ppath_abs, ppath_tau, total_tau,
                   emission_dummy,  abs_mat_per_species_agenda, agenda_dummy,
                   ppath, ppath_p, ppath_t, ppath_vmr, ppath_wind_u,
                   ppath_wind_v, ppath_wind_w, ppath_mag_u, ppath_mag_v,
                   ppath_mag_w,
                   f_grid, -1, stokes_dim, atmosphere_dim, 0 );
    }
    Tensor3View ppath_abs_scalar = ppath_abs(joker,joker,0,0,joker);
    // Get transmitted signal
    //
    iy_transmitter_agendaExecute( ws, iy, ppath.pos(np-1,Range(0,atmosphere_dim)),
                            ppath.los(np-1,joker), f_grid, iy_transmitter_agenda );

    // Fill some parts of iy_aux
    // Radiance
    if( aux_radiance >= 0 )
        { iy_aux[aux_radiance](joker,joker,np-1) = iy; }


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
            Tensor5   ppath_as2_tensor5;
      //
      if( j_analytical_do )
        { 
          diy_dpath.resize( nq ); 
          abs_species_i.resize( nq ); 
          is_t.resize( nq ); 
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
              Tensor4 tau_dummy;
              Tensor3 total_tau_dummy;
              Vector t2 = ppath_t;
              t2 += dt;
              get_ppath_rtvars( ws, ppath_as2_tensor5, tau_dummy, total_tau_dummy,
                                emission_dummy,  abs_mat_per_species_agenda, agenda_dummy,
                                ppath, ppath_p, t2, ppath_vmr, ppath_wind_u,
                                ppath_wind_v, ppath_wind_w, ppath_mag_u, ppath_mag_v,
                                ppath_mag_w,
                                f_grid, -1, stokes_dim, atmosphere_dim, 0 );
            }
        }
      
      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {

          // Jacobian quantities
          if( j_analytical_do )
            {
              // Common terms introduced for efficiency and clarity
              Vector X(nf);   // See AUG
              //
              for( Index iv=0; iv<nf; iv++ )
                {
                  X[iv] = 0.5 * ppath.lstep[ip] * exp( -total_tau(iv,0,0) );
                }

              // Loop quantities
              // Loop quantities
              for( Index iq=0; iq<nq; iq++ ) 
                {
                  if( jacobian_quantities[iq].Analytical() )
                    {
                      // Absorbing species
                      const Index isp = abs_species_i[iq]; 
                      if( isp >= 0 )
                        {
                          // Scaling factors to handle retrieval unit
                          // (gives the cross-section to apply)
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
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // All Stokes components equal
                              for( Index is=0; is<stokes_dim; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is);
                                  diy_dpath[iq](ip  ,iv,is) += Z * unitscf1 *
                                                 ppath_abs_scalar(isp,iv,ip);
                                  diy_dpath[iq](ip+1,iv,is) += Z * unitscf2 *
                                                 ppath_abs_scalar(isp,iv,ip+1);
                                }
                            }
                        }
                       // Temperature
                      else if( is_t[iq] )
                        {
                          Tensor3View ppath_as2 = ppath_as2_tensor5(joker, joker, 0, 0, joker);
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // The terms associated with Dtau/Dk:
                              const Numeric k1 = 
                                         ppath_abs_scalar(joker,iv,ip  ).sum();
                              const Numeric k2 = 
                                         ppath_abs_scalar(joker,iv,ip+1).sum();
                              const Numeric dkdt1 =
                                  ( ppath_as2(joker,iv,ip  ).sum() - k1 ) / dt;
                              const Numeric dkdt2 =
                                  ( ppath_as2(joker,iv,ip+1).sum() - k2 ) / dt;
                              for( Index is=0; is<stokes_dim; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is);
                                  diy_dpath[iq](ip, iv, is) += Z * dkdt1;
                                  diy_dpath[iq](ip+1, iv, is) += Z * dkdt2;
                                }
                              //
                              // The terms associated with Delta-s:
                              if( jacobian_quantities[iq].Subtag() == "HSE on" )
                                {
                                  const Numeric kbar = 0.5 * ( k1 + k2 );
                                  for( Index is=0; is<stokes_dim; is++ )
                                    { 
                                      const Numeric Z = -X[iv] * iy(iv,is);
                                      diy_dpath[iq](ip  ,iv,is) += Z * kbar /
                                                                  ppath_t[ip];
                                      diy_dpath[iq](ip+1,iv,is) += Z * kbar /
                                                                 ppath_t[ip+1];
                                    }
                                } //hse
                            }  // frequency
                        }  // if is_t
                    } // if analytical
                } // for iq
              
              // Update total_tau (not used for spectrum calculation)
              for( Index iv=0; iv<nf; iv++ )
                { total_tau(iv,joker,joker) -= ppath_tau(iv,joker,joker,ip); }
            }

          // Spectrum
          //
          Matrix temp_ppath_tau, temp(stokes_dim, stokes_dim);
          Vector temp_iy;
          for( Index iv=0; iv<nf; iv++ )
            {
              //
                temp_ppath_tau  = ppath_tau(iv, joker, joker, ip);
                temp_ppath_tau *= -1;
                temp_iy         = iy(iv, joker);
                matrix_exp(temp, temp_ppath_tau, 10);
                mult(iy(iv,joker), temp , temp_iy);
            } 

          // --- iy_aux: ------------------------------------------------------
          // Temperature
          if( aux_temperature >= 0 ) 
            { iy_aux[aux_temperature](0,0,ip) = ppath_t[ip]; }
          // VMR
          for( Index j=0; j<aux_vmr_species.nelem(); j++ )
            { iy_aux[aux_vmr_species[j]](0,0,ip) = 
                                                ppath_vmr(aux_vmr_isp[j],ip); }
          // Absorption (change to store "absorption vector")
          if( aux_abs_sum >= 0 ) 
            { for( Index iv=0; iv<nf; iv++ ) {
              iy_aux[aux_abs_sum](iv,0,ip) = 
                                       ppath_abs_scalar(iv,joker,ip).sum(); } }
          for( Index j=0; j<aux_abs_species.nelem(); j++ )
            { for( Index iv=0; iv<nf; iv++ ) {
              iy_aux[aux_abs_species[j]](iv,0,ip) = 
                                    ppath_abs_scalar(iv,aux_abs_isp[j],ip); } }
          // Radiance 
          if( aux_radiance >= 0 ) 
            { iy_aux[aux_radiance](joker,joker,ip) = iy; }
          // --- iy_aux: ------------------------------------------------------

        } 

      // Map jacobians from ppath to retrieval grids
      // (this operation corresponds to the term Dx_i/Dx)
      if( j_analytical_do )
        { 
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_from_path_to_rgrids( diy_dx[iq], jacobian_quantities[iq], 
                               diy_dpath[iq], atmosphere_dim, ppath, ppath_p );
          )
        }
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
/*
void iyTransmissionStandardCloudbox(
         Workspace&                     ws,
         Matrix&                        iy,
         ArrayOfTensor3&                diy_dx,
   const Tensor3&                       iy_transmission,
   const Vector&                        rte_pos,
   const Vector&                        rte_los,
   const Index&                         jacobian_do,
   const Index&                         atmosphere_dim,
   const Vector&                        p_grid,
   const Tensor3&                       z_field,
   const Tensor3&                       t_field,
   const Tensor4&                       vmr_field,
   const Tensor3&                       edensity_field,
   const Index&                         cloudbox_on,
   const ArrayOfIndex&                  cloudbox_limits,
   const Index&                         stokes_dim,
   const Vector&                        f_grid,
   const Agenda&                        ppath_agenda,
   const Agenda&                        abs_scalar_gas_agenda,
   const Agenda&                        iy_main_agenda,
   const Tensor4&                       pnd_field,
   const Index&                         use_mean_scat_data,
   const ArrayOfSingleScatteringData&   scat_data_raw,
   const Agenda&                        opt_prop_gas_agenda,
   const Verbosity&                     verbosity)
{
  // Input checks
  if( !cloudbox_on )
    throw runtime_error( "The cloudbox must be defined to use this method." );
  if( jacobian_do )
    throw runtime_error( 
       "This method does not provide any jacobians (jacobian_do must be 0)." );

  // Sizes
  //
  const Index nf = f_grid.nelem();

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
  Matrix    emission_dummy, ppath_tau;
  Tensor3   wind_field_dummy(0,0,0), iy_trans_new;
  Tensor3   ppath_asp_abs_vec, ppath_pnd_abs_vec, total_transmission;
  Tensor4   ppath_asp_ext_mat, ppath_pnd_ext_mat, ppath_transmission;
  Agenda    agenda_dummy;
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
      get_ppath_pnd( ppath_pnd, ppath, atmosphere_dim, cloudbox_limits, 
                     pnd_field );

      // Absorption and optical thickness for each step
      get_ppath_cloudrtvars( ws, ppath_asp_abs_vec, ppath_asp_ext_mat,  
                             ppath_pnd_abs_vec, ppath_pnd_ext_mat, 
                             ppath_transmission, total_transmission, 
                             emission_dummy, scat_data, abs_scalar_gas_agenda, 
                             agenda_dummy, opt_prop_gas_agenda,
                             ppath, ppath_p, ppath_t, ppath_vmr, ppath_wind_u,  
                             ppath_wind_v, ppath_wind_w, ppath_pnd, 
                             use_mean_scat_data, scat_data_raw, stokes_dim, 
                             f_grid, atmosphere_dim, 0, verbosity );

      // *scat_data* not used. Free the memory.
      scat_data.resize( 0 );
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
  {
    Vector  rte_pos2, rte_los2;
    rte_pos2 = ppath.pos(ppath.np-1,Range(0,atmosphere_dim));
    rte_los2 = ppath.los(ppath.np-1,joker);
    //
    iy_main_agendaExecute( ws, iy, iy_aux, diy_dx, 0, iy_trans_new,
                               0, jacobian_do, t_field, z_field, vmr_field, 
                               -1, rte_pos2, rte_los2, iy_main_agenda );
  }

  // Without jacobians the complete RT calculations are just multiplication
  // with *total_transmission*
  for( Index iv=0; iv<nf; iv++ )
    {
      const Vector tmp = iy(iv,joker);
      mult( iy(iv,joker), total_transmission(iv,joker,joker), tmp );
    }
}
*/





/* Workspace method: Doxygen documentation will be auto-generated */
void yCalc(
         Workspace&                  ws,
         Vector&                     y,
         Vector&                     y_f,
         ArrayOfIndex&               y_pol,
         Matrix&                     y_pos,
         Matrix&                     y_los,
         ArrayOfVector&              y_aux,
         Matrix&                     jacobian,
   const Index&                      basics_checked,
   const Index&                      atmosphere_dim,
   const Tensor3&                    t_field,
   const Tensor3&                    z_field,
   const Tensor4&                    vmr_field,
   const Index&                      cloudbox_on,
   const Index&                      cloudbox_checked,
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
   const Agenda&                     iy_main_agenda,
   const String&                     y_unit,
   const Agenda&                     jacobian_agenda,
   const Index&                      jacobian_do,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const ArrayOfString&              iy_aux_vars,
   const Verbosity&                  verbosity )
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

  // Basics and cloudbox OK?
  //
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );

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
  chk_if_in_range( "antenna_dim", antenna_dim, 1, 2 );
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

  // For y_aux we don't know the number of quantities, and we need to 
  // store all output
  ArrayOfArrayOfVector  iyb_aux_array( nmblock );

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
  else
    { jacobian.resize( 0, 0 ); }


  //---------------------------------------------------------------------------
  // The calculations
  //---------------------------------------------------------------------------

  // We have to make a local copy of the Workspace and the agendas because
  // only non-reference types can be declared firstprivate in OpenMP
  Workspace l_ws (ws);
  Agenda l_jacobian_agenda (jacobian_agenda);
  Agenda l_iy_main_agenda (iy_main_agenda);

/*#pragma omp parallel for                                           \
  if(!arts_omp_in_parallel() && nmblock>1 && nmblock>=nza)        \
    default(none)                                                   \
    firstprivate(l_ws, l_jacobian_agenda, l_iy_main_agenda)     \
    shared(j_analytical_do, sensor_los, mblock_za_grid, mblock_aa_grid, \
           vmr_field, t_field, lon_grid, lat_grid, p_grid, f_grid,      \
           sensor_pos, joker, naa)*/
#pragma omp parallel for                                          \
  if(!arts_omp_in_parallel() && nmblock>1 && nmblock>=nza)        \
  firstprivate(l_ws, l_jacobian_agenda, l_iy_main_agenda)
  for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
    {
      // Calculate monochromatic pencil beam data for 1 measurement block
      //
      Vector          iyb, iyb_error, yb(n1y);
      ArrayOfMatrix   diyb_dx;      
      //
      iyb_calc( l_ws, iyb, iyb_aux_array[mblock_index], diyb_dx, 
                mblock_index, atmosphere_dim, t_field, 
                z_field, vmr_field, cloudbox_on, stokes_dim, f_grid, sensor_pos,
                sensor_los, mblock_za_grid, mblock_aa_grid, antenna_dim, 
                l_iy_main_agenda, y_unit, j_analytical_do, 
                jacobian_quantities, jacobian_indices, iy_aux_vars, verbosity );


      // Apply sensor response matrix on iyb, and put into y
      //
      const Range rowind = get_rowindex_for_mblock( sensor_response, 
                                                                mblock_index );
      const Index row0   = rowind.get_start();
      //
      mult( yb, sensor_response, iyb );
      //
      y[rowind] = yb;  // *yb* also used below, as input to jacobian_agenda

      // Fill information variables
      //
      for( Index i=0; i<n1y; i++ )
        { 
          y_f[row0+i]         = sensor_response_f[i];
          y_pol[row0+i]       = sensor_response_pol[i]; 
          y_pos(row0+i,joker) = sensor_pos(mblock_index,joker);
          y_los(row0+i,0)     = sensor_los(mblock_index,0) + 
                                sensor_response_za[i];
          if( sensor_response_aa.nelem() )
            { 
              y_los(row0+i,1) = sensor_los(mblock_index,0) + 
                                sensor_response_aa[i]; 
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

      // Rest of *jacobian*
      //
      if( jacobian_do )
        { 
          jacobian_agendaExecute( l_ws, jacobian, mblock_index, iyb, yb, 
                                                            l_jacobian_agenda );
        }
    }  // End mblock loop

  // Compile y_aux
  //
  const Index nq = iyb_aux_array[0].nelem();
  y_aux.resize( nq );
  //
  for( Index q=0; q<nq; q++ )
    {
      y_aux[q].resize( nmblock*n1y );
      //
      for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
        {
          const Range rowind = get_rowindex_for_mblock( sensor_response, 
                                                                mblock_index );
          mult( y_aux[q][rowind], sensor_response, 
                                              iyb_aux_array[mblock_index][q] );
        }
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void yApplyYunit(
         Vector&         y,
         Matrix&         jacobian,
   const Vector&         y_f,
   const ArrayOfIndex&   y_pol,
   const String&         y_unit,
   const Verbosity&)
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

  // Is jacobian set?
  //
  const Index ny = y.nelem();
  //
  const bool do_j = jacobian.nrows() == ny;

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
          // Find number of values for this chunk
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

          if( do_j )
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
          yv(0,0)  = y[i];    
          i_pol[0] = y_pol[i];

          // Jacobian
          if( do_j )
            { apply_y_unit2( MatrixView( jacobian(i,joker) ), yv, 
                                                     y_unit, y_f[i], i_pol ); }

          // y (must be done last)
          apply_y_unit( yv, y_unit, y_f[i], i_pol );
          y[i] = yv(0,0);
        }
    }
}
