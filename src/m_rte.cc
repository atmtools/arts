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
void iyApplyUnit(
         Matrix&         iy,
   ArrayOfTensor4&       iy_aux,
   const Index&          stokes_dim,
   const Vector&         f_grid,
   const ArrayOfString&  iy_aux_vars,
   const String&         iy_unit,
   const Verbosity&)
{
  if( iy_unit == "1" )
    throw runtime_error( "No need to use this method with *iy_unit* = \"1\"." );

  if( max(iy(joker,0)) > 1e-3 )
    {
      ostringstream os;
      os << "The spectrum matrix *iy* is required to have original radiance\n"
         << "unit, but this seems not to be the case. This as a value above\n"
         << "1e-3 is found in *iy*.";
      throw runtime_error( os.str() );      
    }

  // Polarisation index variable
  ArrayOfIndex i_pol(stokes_dim);
  for( Index is=0; is<stokes_dim; is++ )
    { i_pol[is] = is + 1; }

  apply_iy_unit( iy, iy_unit, f_grid, 1, i_pol );
  
  for( Index i=0; i<iy_aux_vars.nelem(); i++ )
    {
      if( iy_aux_vars[i] == "iy"  ||  iy_aux_vars[i] == "Error"  || 
          iy_aux_vars[i] == "Error (uncorrelated)" )
        {
          if( iy_aux[i].nrows() > 1 )
            throw runtime_error( "Data marked as \"iy\" or \"Error\" "
                                 "have incorrect size." );
          for( Index j=0; j<iy_aux[i].ncols(); j++ )
            { apply_iy_unit( iy_aux[i](joker,joker,0,j), iy_unit, f_grid, 1, 
                                                                      i_pol ); }
        }
    }
}




/* Workspace method: Doxygen documentation will be auto-generated */
void iyCalc(
         Workspace&        ws,
         Matrix&           iy,
         ArrayOfTensor4&   iy_aux,
         Ppath&            ppath,
   const Index&            basics_checked,
   const ArrayOfString&    iy_aux_vars,
   const Vector&           f_grid,
   const Tensor3&          t_field,
   const Tensor3&          z_field,
   const Tensor4&          vmr_field,
   const Index&            cloudbox_on,
   const Index&            cloudbox_checked,
   const Vector&           rte_pos,
   const Vector&           rte_los,
   const Vector&           rte_pos2,
   const Agenda&           iy_main_agenda,
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

  ArrayOfTensor3 diy_dx;

  iy_main_agendaExecute( ws, iy, iy_aux, ppath, diy_dx, 1, iy_transmission, 
                         iy_aux_vars, cloudbox_on, 0, t_field, 
                         z_field, vmr_field, f_grid, rte_pos, rte_los, rte_pos2,
                         iy_main_agenda );
}





/* Workspace method: Doxygen documentation will be auto-generated */
void iyEmissionStandard(
         Workspace&                  ws,
         Matrix&                     iy,
         ArrayOfTensor4&             iy_aux,
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
   const String&                     iy_unit,
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
   const Vector&                     rte_pos,      
   const Vector&                     rte_los,      
   const Vector&                     rte_pos2,      
   const Verbosity&                  verbosity )
{
  // Determine propagation path
  //
  ppath_agendaExecute( ws, ppath, rte_pos, rte_los, rte_pos2, cloudbox_on, 0, 
                       t_field, z_field, vmr_field, edensity_field, f_grid, 
                       ppath_agenda );

  // Some basic sizes
  //
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index nq = jacobian_quantities.nelem();
  
  //=== iy_aux part ===========================================================
  Index auxPressure    = -1,
        auxTemperature = -1,
        auxAbsSum      = -1,
        auxBackground  = -1,
        auxIy          = -1,
        auxOptDepth    = -1;
  ArrayOfIndex auxAbsSpecies(0), auxAbsIsp(0);
  ArrayOfIndex auxVmrSpecies(0), auxVmrIsp(0);
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
            { iy_aux[i].resize( 0, 0, 0, 0 ); }
          else if( iy_aux_vars[i].substr(0,10) == "PND, type " )
            { iy_aux[i].resize( 0, 0, 0, 0 ); }
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


  //### jacobian part #########################################################
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
                           nf, ns ); 
        diy_dx[iq] = 0.0;
      )
    }
  //###########################################################################


  // Get atmospheric and attenuation quantities for each ppath point/step
  //
  // "atmvars"
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w,
                              ppath_mag_u,  ppath_mag_v,  ppath_mag_w;
  Matrix    ppath_vmr;
  // Attenuation vars
  Tensor5   ppath_abs;
  Tensor4   trans_partial, trans_cumulat;
  Matrix    ppath_blackrad;
  Vector    scalar_tau;
  //
  if( np > 1 )
    {
      get_ppath_atmvars(  ppath_p, ppath_t, ppath_vmr,
                          ppath_wind_u, ppath_wind_v, ppath_wind_w,
                          ppath_mag_u,  ppath_mag_v,  ppath_mag_w,
                          ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                          wind_u_field, wind_v_field, wind_w_field ,
                          mag_u_field, mag_v_field, mag_w_field );      
      get_ppath_abs(      ws, ppath_abs, abs_mat_per_species_agenda, ppath, 
                          ppath_p, ppath_t, ppath_vmr, 
                          ppath_wind_u, ppath_wind_v, ppath_wind_w, 
                          ppath_mag_u, ppath_mag_v, ppath_mag_w, 
                          f_grid, stokes_dim, atmosphere_dim );
      get_ppath_trans(    trans_partial, trans_cumulat,scalar_tau, 
                          ppath, ppath_abs, f_grid, stokes_dim );
      get_ppath_blackrad( ws, ppath_blackrad, blackbody_radiation_agenda, 
                          ppath, ppath_t, f_grid );
    }
  else // For cases inside the cloudbox, or totally outside the atmosphere,
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
  get_iy_of_background( ws, iy, diy_dx, 
                        iy_trans_new, jacobian_do, ppath, rte_pos2, 
                        atmosphere_dim, t_field, z_field, vmr_field, 
                        cloudbox_on, stokes_dim, f_grid, iy_main_agenda, 
                        iy_space_agenda, iy_surface_agenda, iy_cloudbox_agenda,
                        verbosity );


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
      //### jacobian part #####################################################
      // Create container for the derivatives with respect to changes
      // at each ppath point. And find matching abs_species-index and 
      // "temperature flag" (for analytical jacobians).
      //
      ArrayOfTensor3  diy_dpath; 
      ArrayOfIndex    abs_species_i, is_t; 
      //
      const Numeric   dt = 0.1;
            Tensor5   ppath_at2;
            Matrix    ppath_bt2;
      //
      if( j_analytical_do )
        { 
          // So far no polarised absorption handled for jacobians
          for( Index iv=0; iv<nf; iv++ ) {
            if( !is_diagonal( trans_cumulat(iv,joker,joker,np-1) ) )
                    throw runtime_error( "The combination of polarised "
                           "absorption and jacobians is not yet handled." ); }
          //------------------------------------------------------------------
          diy_dpath.resize( nq ); 
          abs_species_i.resize( nq ); 
          is_t.resize( nq ); 
          //
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_dpath[iq].resize( np, nf, ns ); 
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
              Vector t2 = ppath_t;
              t2 += dt;
              get_ppath_abs(      ws, ppath_at2, abs_mat_per_species_agenda, 
                                  ppath, ppath_p, t2, ppath_vmr, 
                                  ppath_wind_u, ppath_wind_v, ppath_wind_w, 
                                  ppath_mag_u, ppath_mag_v, ppath_mag_w, 
                                  f_grid, stokes_dim, atmosphere_dim );
              get_ppath_blackrad( ws, ppath_bt2, blackbody_radiation_agenda,
                                  ppath, t2, f_grid );
            }
        }
      //#######################################################################

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

          //### jacobian part #################################################
          if( j_analytical_do )
            {
              // Common terms introduced for efficiency and clarity
              Vector X(nf), Y(nf);   // See AUG
              //
              for( Index iv=0; iv<nf; iv++ )
                {
                  X[iv] = 0.5 * ppath.lstep[ip] * trans_cumulat(iv,0,0,ip+1);
                  Y[iv] = X[iv] * ( bbar[iv] - iy(iv,0) );
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
                                                 ppath_abs(isp,iv,0,0,ip);
                              diy_dpath[iq](ip+1,iv,0) += Y[iv] * unitscf2 * 
                                                 ppath_abs(isp,iv,0,0,ip+1);
                              // Higher stokes components
                              for( Index is=1; is<ns; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is); 
                                  diy_dpath[iq](ip  ,iv,is) += Z * unitscf1 *
                                                 ppath_abs(isp,iv,0,0,ip);
                                  diy_dpath[iq](ip+1,iv,is) += Z * unitscf2 *
                                                 ppath_abs(isp,iv,0,0,ip+1);
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
                                            ppath_abs(joker,iv,0,0,ip  ).sum();
                              const Numeric k2 = 
                                            ppath_abs(joker,iv,0,0,ip+1).sum();
                              const Numeric dkdt1 = 
                               ( ppath_at2(joker,iv,0,0,ip  ).sum() - k1 ) / dt;
                              const Numeric dkdt2 =
                               ( ppath_at2(joker,iv,0,0,ip+1).sum() - k2 ) / dt;
                              // Stokes 1:
                              diy_dpath[iq](ip  ,iv,0) += Y[iv] * dkdt1;
                              diy_dpath[iq](ip+1,iv,0) += Y[iv] * dkdt2;
                              // Higher Stokes
                              for( Index is=1; is<ns; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is);
                                  diy_dpath[iq](ip  ,iv,is) += Z * dkdt1;
                                  diy_dpath[iq](ip+1,iv,is) += Z * dkdt2;
                                }
                              //
                              // The terms associated with B-bar:
                              const Numeric V = 0.5 * trans_cumulat(iv,0,0,ip) *
                                              ( 1.0 - trans_partial(iv,0,0,ip));
                              diy_dpath[iq](ip  ,iv,0) += V *
                                              ( ppath_bt2(iv,ip) -
                                                ppath_blackrad(iv,ip) ) / dt;
                              diy_dpath[iq](ip+1,iv,0) += V * 
                                              ( ppath_bt2(iv,ip+1) -
                                                ppath_blackrad(iv,ip+1) ) / dt;
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
                                  for( Index is=1; is<ns; is++ )
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
            }
          //###################################################################

          // Spectrum at end of ppath step 
          if( stokes_dim == 1 )
            {
              for( Index iv=0; iv<nf; iv++ )  
                { iy(iv,0) = iy(iv,0) * trans_partial(iv,0,0,ip) +
                       bbar[iv] * ( 1 - trans_partial(iv,0,0,ip) ); }
            }
          else
            {
              for( Index iv=0; iv<nf; iv++ )  
                {
                  // Unpolarised absorption:
                  if( is_diagonal( trans_partial(iv,joker,joker,ip) ) )
                    {
                      iy(iv,0) = iy(iv,0) * trans_partial(iv,0,0,ip) +
                           bbar[iv] * ( 1 - trans_partial(iv,0,0,ip) );
                      for( Index is=1; is<ns; is++ )
                        { iy(iv,is) = iy(iv,is) * trans_partial(iv,is,is,ip); }
                    }
                  // The general case:
                  else
                    {
                      // Transmitted term
                      Vector tt(ns);
                      mult( tt, trans_partial(iv,joker,joker,ip), iy(iv,joker));
                      // Add emission, first Stokes element
                      iy(iv,0) = tt[0] + bbar[iv]*(1-trans_partial(iv,0,0,ip));
                      // Remaining Stokes elements
                      for( Index i=1; i<ns; i++ )
                        { iy(iv,i) = tt[i] - bbar[iv]*trans_partial(iv,i,0,ip);}
                      
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
                                  ppath_abs(joker,iv,is1,is2,ip).sum(); } } } } 
          for( Index j=0; j<auxAbsSpecies.nelem(); j++ )
            { for( Index iv=0; iv<nf; iv++ ) {
                for( Index is1=0; is1<stokes_dim; is1++ ){
                  for( Index is2=0; is2<stokes_dim; is2++ ){
                    iy_aux[auxAbsSpecies[j]](iv,is1,is2,ip) = 
                                 ppath_abs(auxAbsIsp[j],iv,is1,is2,ip); } } } }
          // Radiance 
          if( auxIy >= 0 ) 
            { iy_aux[auxIy](joker,joker,0,ip) = iy; }
          //===================================================================
        } 

      //### jacobian part #####################################################
      // Map jacobians from ppath to retrieval grids
      // (this operation corresponds to the term Dx_i/Dx)
      if( j_analytical_do )
        { 
          // Weight with iy_transmission
          if( !iy_agenda_call1 )
            {
              Matrix X, Y(ns,diy_dpath[0].npages()); 
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
      //#######################################################################
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

      // Jacobian part (must be converted to Tb before iy for PlanckBT)
      // 
      if( j_analytical_do )
        {
          FOR_ANALYTICAL_JACOBIANS_DO( apply_iy_unit2( diy_dx[iq], iy, iy_unit,
                                                       f_grid, n, i_pol ); )
        } 

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
void iyLoopFrequencies(
         Workspace&        ws,
         Matrix&           iy,
         ArrayOfTensor4&   iy_aux,
         Ppath&            ppath,
         ArrayOfTensor3&   diy_dx,
   const ArrayOfString&    iy_aux_vars,
   const Index&            stokes_dim,
   const Vector&           f_grid,
   const Tensor3&          t_field,
   const Tensor3&          z_field,
   const Tensor4&          vmr_field,
   const Index&            cloudbox_on,
   const Index&            iy_agenda_call1,
   const Tensor3&          iy_transmission,
   const Vector&           rte_pos,
   const Vector&           rte_los,
   const Vector&           rte_pos2,
   const Index&            jacobian_do,
   const Agenda&           iy_sub_agenda,
   const Verbosity& )
{
  // Throw error if unsupported features are requested
  if( !iy_agenda_call1 )
    throw runtime_error( 
                  "Recursive usage not possible (iy_agenda_call1 must be 1)" );
  if( iy_transmission.ncols() )
    throw runtime_error( "*iy_transmission* must be empty" );

  const Index nf = f_grid.nelem();

  for( Index i=0; i<nf; i++ )
    {
      // Variables for 1 frequency
      Matrix         iy1;
      ArrayOfTensor4 iy_aux1; 
      ArrayOfTensor3 diy_dx1;
      
      iy_sub_agendaExecute( ws, iy1, iy_aux1, ppath, diy_dx1, 
                            1, iy_transmission, iy_aux_vars, cloudbox_on, 
                            jacobian_do, t_field, z_field, vmr_field, 
                            Vector(1,f_grid[i]),
                            rte_pos, rte_los, rte_pos2, iy_sub_agenda );

      // After first frequency, give output its size
      if( i == 0 )
        {
          iy.resize( nf, stokes_dim );
          //
          iy_aux.resize( iy_aux1.nelem() );
          for( Index q=0; q<iy_aux1.nelem(); q++ )
            {
              if( iy_aux1[q].ncols() > 1 )
                throw runtime_error( "When using this method, *iy_aux_vars* "
                        "is not allowed to include along-the-path variables." );
              iy_aux[q].resize(nf,iy_aux1[q].npages(),iy_aux1[q].nrows(),1);
            }
          //
          diy_dx.resize( diy_dx1.nelem() );
          for( Index q=0; q<diy_dx1.nelem(); q++ )
            { diy_dx[q].resize( diy_dx1[q].npages(), nf, stokes_dim ); }
        }

      // Copy to output variables
      iy(i,joker) = iy1(0,joker);
      for( Index q=0; q<iy_aux1.nelem(); q++ )
        { iy_aux[q](i,joker,joker,0) = iy_aux1[q](0,joker,joker,0); }
      for( Index q=0; q<diy_dx1.nelem(); q++ )
        { diy_dx[q](joker,i,joker) = diy_dx1[q](joker,0,joker); }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void iyMC(
         Workspace&                  ws,
         Matrix&                     iy,
         ArrayOfTensor4&             iy_aux,
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
   const Tensor3&                    edensity_field,
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
   const Agenda&                     ppath_step_agenda, 
   const Tensor4&                    pnd_field,
   const String&                     iy_unit,
   const Numeric&                    mc_std_err,
   const Index&                      mc_max_time,
   const Index&                      mc_max_iter,
   const Index&                      mc_min_iter,
   const Verbosity&                  verbosity)
{
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


  // Size output variables
  //
  const Index   nf = f_grid.nelem();
  //
  iy.resize( nf, stokes_dim );
  diy_dx.resize(0);
  //
  //=== iy_aux part ===========================================================
  Index auxError = -1;
  {
    const Index naux = iy_aux_vars.nelem();
    iy_aux.resize( naux );
    //
    for( Index i=0; i<naux; i++ )
      {
        if( iy_aux_vars[i] == "Error (uncorrelated)" )
          { auxError = i;      iy_aux[i].resize( nf, stokes_dim, 1, 1 ); }
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

  // Some MC variables are only local here
  //
  MCAntenna mc_antenna;
  mc_antenna.set_pencil_beam();

  // Pos and los must be matrices 
  Matrix pos(1,3), los(1,2);
  //
  pos(0,joker) = rte_pos;
  los(0,joker) = rte_los;

  Workspace l_ws (ws);
  Agenda l_ppath_step_agenda (ppath_step_agenda);
  Agenda l_iy_space_agenda (iy_space_agenda);
  Agenda l_abs_mat_per_species_agenda (abs_mat_per_species_agenda);
  Agenda l_surface_rtprop_agenda (surface_rtprop_agenda);

  String fail_msg;
  bool failed = false;

#pragma omp parallel for                                          \
  if(!arts_omp_in_parallel()) \
  firstprivate(l_ws, l_ppath_step_agenda, l_iy_space_agenda, l_abs_mat_per_species_agenda, l_surface_rtprop_agenda)
  for( Index f_index=0; f_index<nf; f_index++ )
    {
      if (failed) continue;

      try {
        ArrayOfSingleScatteringData   scat_data_mono;

        scat_data_monoCalc( scat_data_mono, scat_data_raw,
                            f_grid, f_index, verbosity );

        // Seed reset for each loop. If not done, the errors
        // appear to be highly correlated.
        Index    mc_seed;
        MCSetSeedFromTime( mc_seed, verbosity );

        Vector y, mc_error;
        Index    mc_iteration_count;
        Tensor3  mc_points;

        MCGeneral( l_ws, y, mc_iteration_count, mc_error, mc_points, mc_antenna,
                   f_grid, f_index, pos, los, stokes_dim, atmosphere_dim,
                   l_ppath_step_agenda, l_iy_space_agenda, 
                   l_surface_rtprop_agenda, l_abs_mat_per_species_agenda, 
                   p_grid, lat_grid, lon_grid, z_field, 
                   refellipsoid, z_surface, t_field, vmr_field,
                   edensity_field, cloudbox_on, cloudbox_limits,
                   pnd_field, scat_data_mono, 1, cloudbox_checked,
                   mc_seed, iy_unit, mc_std_err, mc_max_time, mc_max_iter,
                   mc_min_iter, verbosity);
          //cout << "Error: "      << mc_error << endl;
          //cout << "N photons: " << mc_iteration_count << endl;
        assert( y.nelem() == stokes_dim );

        iy(f_index,joker) = y;
          
        if( auxError >= 0 ) 
          { iy_aux[auxError](f_index,joker,0,0) = mc_error; }
      } catch (runtime_error e) {
        ostringstream os;
        os << "Error for f_index = " << f_index << " (" << f_grid[f_index] 
           << ")" << endl << e.what();
#pragma omp critical (iyMC_fail)
          { failed = true; fail_msg = os.str(); }
          continue;
      }
    }

  if (failed)
    throw runtime_error(fail_msg);
}





/* Workspace method: Doxygen documentation will be auto-generated */
void iyRadioLink(
         Workspace&                  ws,
         Matrix&                     iy,
         ArrayOfTensor4&             iy_aux,
         Ppath&                      ppath,
         ArrayOfTensor3&             diy_dx,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const Index&                      atmosphere_dim,
   const Vector&                     p_grid,
   const Vector&                     lat_grid,
   const Vector&                     lon_grid,
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
   const Vector&                     refellipsoid,
   const Matrix&                     z_surface,
   const Index&                      cloudbox_on,
   const ArrayOfIndex&               cloudbox_limits,
   const ArrayOfString&              iy_aux_vars,
   const Index&                      jacobian_do,
   const Agenda&                     ppath_agenda,
   const Agenda&                     ppath_step_agenda,
   const Agenda&                     abs_mat_per_species_agenda,
   const Agenda&                     iy_transmitter_agenda,
   const Index&                      iy_agenda_call1,
   const Tensor3&                    iy_transmission,
   const Vector&                     rte_pos,      
   const Vector&                     rte_pos2,      
   const Verbosity&                  verbosity )
{
  // Throw error if unsupported features are requested
  if( cloudbox_on  || cloudbox_limits.nelem() )
    throw runtime_error( "Cloudbox not yet handled." );
  if( !iy_agenda_call1 )
    throw runtime_error( 
                  "Recursive usage not possible (iy_agenda_call1 must be 1)" );
  if( iy_transmission.ncols() )
    throw runtime_error( "*iy_transmission* must be empty" );
  if( jacobian_do )
    throw runtime_error( "This method does not yet provide any jacobians and "
                         "*jacobian_do* must be 0." );
  diy_dx.resize(0);


  //- Determine propagation path
  ppath_agendaExecute( ws, ppath, rte_pos, Vector(0), rte_pos2, cloudbox_on, 0,
                       t_field, z_field, vmr_field, edensity_field, f_grid, 
                       ppath_agenda );
  if( ppath_what_background(ppath) != 9 )
    { throw runtime_error( "Radiative background not set to \"transmitter\" by"
                     " *ppath_agenda*. Is correct WSM used in the agenda?" ); }

  // Some basic sizes
  //
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;

  //=== iy_aux part ===========================================================
  Index auxPressure        = -1,
        auxTemperature     = -1,
        auxAbsSum          = -1,
        auxFreeSpaceLoss   = -1,
        auxFreeSpaceAtte   = -1,
        auxAtmosphericLoss = -1,
        auxDefocusingLoss  = -1,
        auxDefocusingAtte  = -1,
        auxExtraPathDelay  = -1,
        auxBendingAngle    = -1;
  ArrayOfIndex auxAbsSpecies(0), auxAbsIsp(0);
  ArrayOfIndex auxVmrSpecies(0), auxVmrIsp(0);
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
      else if( iy_aux_vars[i] == "Free space loss" )
        { auxFreeSpaceLoss = i;     iy_aux[i].resize( nf, 1, 1, 1 ); }
      else if( iy_aux_vars[i] == "Free space attenuation" )
        { auxFreeSpaceAtte = i;     iy_aux[i].resize( nf, 1, 1, np ); }
      else if( iy_aux_vars[i] == "Atmospheric loss" )
        { auxAtmosphericLoss = i;   iy_aux[i].resize( nf, 1, 1, 1 ); } 
      else if( iy_aux_vars[i] == "Defocusing loss" )
        { auxDefocusingLoss = i;    iy_aux[i].resize( nf, 1, 1, 1 ); }
      else if( iy_aux_vars[i] == "Defocusing attenuation" )
        { auxDefocusingAtte = i;    iy_aux[i].resize( nf, 1, 1, np ); }
      else if( iy_aux_vars[i] == "Extra path delay" )
        { auxExtraPathDelay = i;    iy_aux[i].resize( nf, 1, 1, 1 ); }
      else if( iy_aux_vars[i] == "Bending angle" )
        { auxBendingAngle = i;      iy_aux[i].resize( nf, 1, 1, 1 ); } 
      else
        {
          ostringstream os;
          os << "In *iy_aux_vars* you have included: \"" << iy_aux_vars[i]
             << "\"\nThis choice is not recognised.";
          throw runtime_error( os.str() );
        }
    }
  //===========================================================================


  // Get atmospheric and attenuation quantities for each ppath point/step
  //
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w,
                              ppath_mag_u,  ppath_mag_v,  ppath_mag_w;
  Matrix    ppath_vmr;
  Tensor5   ppath_abs;
  Tensor4   trans_partial, trans_cumulat;
  Vector    scalar_tau;
  //
  if( np > 1 )
    {
      get_ppath_atmvars( ppath_p, ppath_t, ppath_vmr,
                         ppath_wind_u, ppath_wind_v, ppath_wind_w,
                         ppath_mag_u,  ppath_mag_v,  ppath_mag_w,
                         ppath, atmosphere_dim, p_grid, t_field, vmr_field,
                         wind_u_field, wind_v_field, wind_w_field ,
                         mag_u_field, mag_v_field, mag_w_field );      
      get_ppath_abs(     ws, ppath_abs, abs_mat_per_species_agenda, ppath, 
                         ppath_p, ppath_t, ppath_vmr, 
                         ppath_wind_u, ppath_wind_v, ppath_wind_w, 
                         ppath_mag_u, ppath_mag_v, ppath_mag_w, 
                         f_grid, stokes_dim, atmosphere_dim );
      get_ppath_trans(   trans_partial, trans_cumulat,scalar_tau, 
                         ppath, ppath_abs, f_grid, stokes_dim );

    }

  // Transmitted signal
  //
  iy_transmitter_agendaExecute( ws, iy, f_grid, rte_pos, Vector(0),
                                iy_transmitter_agenda ); 
  if( iy.ncols() != stokes_dim  ||  iy.nrows() != nf )
    { throw runtime_error( "The size of *iy* returned from "
                                 "*iy_transmitter_agenda* is not correct." ); }

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
      // Free space
      if( auxFreeSpaceAtte >= 0 )
        { iy_aux[auxFreeSpaceAtte](joker,0,0,np-1) = 2/lbg; }
      //=======================================================================

      // Loop ppath steps
      for( Index ip=np-2; ip>=0; ip-- )
        {
          // Lengths
          lbg += ppath.lstep[ip];
          lba += ppath.lstep[ip] * (ppath.ngroup[ip]+ppath.ngroup[ip+1]) / 2.0;
          
          // Atmospheric loss of path step
          for( Index iv=0; iv<nf; iv++ )
            {
              Vector iy_temp = iy(iv,joker);
              mult( iy(iv,joker), trans_partial(iv,joker,joker,ip), iy_temp );
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
          // Free space loss
          if( auxFreeSpaceAtte >= 0 )
            { iy_aux[auxFreeSpaceAtte](joker,0,0,ip) = 2/lbg; }
          //===================================================================
        }

      // Remaing length of ppath
      lbg += ppath.start_lstep;
      lba += ppath.start_lstep;

      // Determine total free space loss
      Numeric fspl = 1 / ( 4 * PI * lbg*lbg ); 

      // Determine defocusing loss
      Numeric dfl = 1;
      if( 0 )
        { defocusing_sat2sat( ws, dfl, ppath_step_agenda, atmosphere_dim, 
                              p_grid, lat_grid, lon_grid, t_field, z_field, 
                              vmr_field, edensity_field, -1, 
                              refellipsoid, z_surface, ppath, verbosity ); 
        }
      else
        {
          defocusing_general( ws, dfl, ppath_step_agenda, atmosphere_dim, 
                              p_grid, lat_grid, lon_grid, t_field, z_field, 
                              vmr_field, edensity_field, -1, 
                              refellipsoid, z_surface, ppath, verbosity ); 
        }

      //=== iy_aux part =======================================================
      if( auxAtmosphericLoss >= 0 )
        { iy_aux[auxAtmosphericLoss](joker,0,0,0) = iy(joker,0); }
      if( auxFreeSpaceLoss >= 0 )
        { iy_aux[auxFreeSpaceLoss] = fspl; }
      if( auxDefocusingLoss >= 0 )
        { iy_aux[auxDefocusingLoss] = dfl; }
      if( auxDefocusingAtte >= 0 )
        { iy_aux[auxDefocusingAtte] = -999; }  // So far just a dummy value
      //=======================================================================

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
          // Geomtrical distance between start and end point
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



/* Workspace method: Doxygen documentation will be auto-generated */
void iyTransmissionStandard(
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
   const Tensor3&                     edensity_field,
   const Index&                       cloudbox_on,
   const ArrayOfIndex&                cloudbox_limits,
   const Tensor4&                     pnd_field,
   const Index&                       use_mean_scat_data,
   const ArrayOfSingleScatteringData& scat_data_raw,
   const Matrix&                      particle_masses,
   const ArrayOfString&               iy_aux_vars,
   const Index&                       jacobian_do,
   const ArrayOfRetrievalQuantity&    jacobian_quantities,
   const ArrayOfArrayOfIndex&         jacobian_indices,
   const Agenda&                      ppath_agenda,
   const Agenda&                      abs_mat_per_species_agenda,
   const Agenda&                      iy_transmitter_agenda,
   const Index&                       iy_agenda_call1,
   const Tensor3&                     iy_transmission,
   const Vector&                      rte_pos,      
   const Vector&                      rte_los,      
   const Vector&                      rte_pos2,      
   const Verbosity&                   verbosity )
{
  // Throw error if unsupported features are requested
  if( !iy_agenda_call1 )
    throw runtime_error( 
                  "Recursive usage not possible (iy_agenda_call1 must be 1)" );
  if( iy_transmission.ncols() )
    throw runtime_error( "*iy_transmission* must be empty" );


  // Determine propagation path
  //
  ppath_agendaExecute( ws, ppath, rte_pos, rte_los, rte_pos2, 0, 0,
                       t_field, z_field, vmr_field, edensity_field, f_grid, 
                       ppath_agenda );

  // Some basic sizes
  //
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index nq = jacobian_quantities.nelem();

  // Get transmitted signal
  //
  iy_transmitter_agendaExecute( ws, iy, 
                                f_grid, ppath.pos(np-1,Range(0,atmosphere_dim)),
                                ppath.los(np-1,joker), iy_transmitter_agenda );
  //
  if( iy.ncols() != stokes_dim  ||  iy.nrows() != nf )
    {
      ostringstream os;
      os << "The size of *iy* returned from *iy_transmitter_agdna* is\n"
         << "not correct:\n"
         << "  expected size = [" << nf << "," << stokes_dim << "]\n"
         << "  size of iy    = [" << iy.nrows() << "," << iy.ncols()<< "]\n";
      throw runtime_error( os.str() );      
    }


  //=== iy_aux part ===========================================================
  Index auxPressure    = -1,
        auxTemperature = -1,
        auxAbsSum      = -1,
        auxPartExt     = -1,
        auxIy          = -1,
        auxTrans       = -1,
        auxOptDepth    = -1;
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
          auxAbsIsp.push_back(ispecies);
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
      else if( iy_aux_vars[i] == "iy"   &&  auxIy < 0 )
        { auxIy = i;           iy_aux[i].resize( nf, ns, 1, np ); }
      else if( iy_aux_vars[i] == "Transmission"   &&  auxTrans < 0 )
        { auxTrans = i;        iy_aux[i].resize( nf, ns, ns, np ); }
      else if( iy_aux_vars[i] == "Optical depth" )
        { auxOptDepth = i;     iy_aux[i].resize( nf, 1, 1, 1 ); }
      else
        {
          ostringstream os;
          os << "In *iy_aux_vars* you have included: \"" << iy_aux_vars[i]
             << "\"\nThis choice is not recognised.";
          throw runtime_error( os.str() );
        }
    }
  //===========================================================================


  //### jacobian part #########################################################
  // Initialise analytical jacobians (diy_dx)
  //
  Index j_analytical_do = 0;
  //
  if( jacobian_do ) { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }
  //
  if( j_analytical_do )
    {
      if( cloudbox_on )
        throw runtime_error( "The combination of active ckloudbox and "
                             "analytical jacibians is not yet handled." );
      diy_dx.resize( jacobian_indices.nelem() ); 
      //
      FOR_ANALYTICAL_JACOBIANS_DO( 
        diy_dx[iq].resize( jacobian_indices[iq][1] - jacobian_indices[iq][0] + 
                           1, nf, stokes_dim ); 
        diy_dx[iq] = 0.0;
      ) 
    }
  //###########################################################################


  // Get atmospheric and RT quantities for each ppath point/step
  //
  Vector    ppath_p, ppath_t, ppath_wind_u, ppath_wind_v, ppath_wind_w,
                              ppath_mag_u,  ppath_mag_v,  ppath_mag_w;
  Matrix    ppath_vmr, ppath_pnd;
  Tensor5   ppath_abs;
  Tensor4   trans_partial, trans_cumulat, pnd_ext_mat;
  Vector    scalar_tau;
  ArrayOfIndex clear2cloudbox;
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
          Tensor3      pnd_abs_vec;
          Array<ArrayOfSingleScatteringData> scat_data;
          //
          get_ppath_ext( clear2cloudbox, pnd_abs_vec, pnd_ext_mat, scat_data, 
                         ppath_pnd, ppath, ppath_t, stokes_dim, f_grid, 
                         atmosphere_dim, cloudbox_limits, pnd_field, 
                         use_mean_scat_data, scat_data_raw, verbosity );
          get_ppath_trans2( trans_partial, trans_cumulat, scalar_tau, 
                            ppath, ppath_abs, f_grid, stokes_dim,
                            clear2cloudbox, pnd_ext_mat );
        }
    }


  //=== iy_aux part ===========================================================
  // Fill parts of iy_aux that are defined even for np=1.
  // Radiance 
  if( auxIy >= 0 ) 
    { iy_aux[auxIy](joker,joker,0,np-1) = iy; }
  if( auxOptDepth >= 0 ) 
    {
      if( np == 1 )
        { iy_aux[auxOptDepth] = 0; }
      else
        { iy_aux[auxOptDepth](joker,0,0,0) = scalar_tau; }
    } 
  if( auxTrans >= 0 ) // Complete tensor filled!
    { 
      if( np == 1 )
        { for( Index iv=0; iv<nf; iv++ ) {
            id_mat( iy_aux[auxTrans](iv,joker,joker,0) ); } }
      else
        { iy_aux[auxTrans] = trans_cumulat; }
    }
  //===========================================================================


  // Do RT calculations
  //
  if( np > 1 )
    {
      //### jacobian part #####################################################
      // Create container for the derivatives with respect to changes
      // at each ppath point. And find matching abs_species-index and 
      // "temperature flag" (for analytical jacobians).
      //
      ArrayOfTensor3  diy_dpath; 
      ArrayOfIndex    abs_species_i, is_t; 
      //
      const Numeric   dt = 0.1;
            Tensor5   ppath_at2;
      //
      if( j_analytical_do )
        { 
          // So far no polarised absorption handled for jacobians
          for( Index iv=0; iv<nf; iv++ ) {
            if( !is_diagonal( trans_cumulat(iv,joker,joker,np-1) ) )
                    throw runtime_error( "The combination of polarised "
                           "absorption and jacobians is not yet handled." ); }
          //------------------------------------------------------------------
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
              Vector t2 = ppath_t;
              t2 += dt;
              get_ppath_abs( ws, ppath_at2, abs_mat_per_species_agenda, 
                             ppath, ppath_p, t2, ppath_vmr, 
                             ppath_wind_u, ppath_wind_v, ppath_wind_w, 
                             ppath_mag_u, ppath_mag_v, ppath_mag_w, 
                             f_grid, stokes_dim, atmosphere_dim );
            }
        }
      //#######################################################################

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
          // Extinction
          if( auxPartExt >= 0  && clear2cloudbox[np-1] >= 0 ) 
            { 
              const Index ic = clear2cloudbox[np-1];
              for( Index iv=0; iv<nf; iv++ ) {
                for( Index is1=0; is1<ns; is1++ ){
                  for( Index is2=0; is2<ns; is2++ ){
                    iy_aux[auxPartExt](iv,is1,is2,np-1) = 
                                              pnd_ext_mat(iv,is1,is2,ic); } } } 
            } 
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
          //### jacobian part #################################################
          if( j_analytical_do )
            {
              // Common terms introduced for efficiency and clarity
              Vector X(nf);   // See AUG
              //
              for( Index iv=0; iv<nf; iv++ )
                { X[iv] = 0.5 * ppath.lstep[ip] * trans_cumulat(iv,0,0,ip+1); }

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
                                                    ppath_abs(isp,iv,0,0,ip);
                                  diy_dpath[iq](ip+1,iv,is) += Z * unitscf2 *
                                                    ppath_abs(isp,iv,0,0,ip+1);
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
                                             ppath_abs(joker,iv,0,0,ip  ).sum();
                              const Numeric k2 = 
                                             ppath_abs(joker,iv,0,0,ip+1).sum();
                              const Numeric dkdt1 =
                               ( ppath_at2(joker,iv,0,0,ip  ).sum() - k1 ) / dt;
                              const Numeric dkdt2 =
                               ( ppath_at2(joker,iv,0,0,ip+1).sum() - k2 ) / dt;
                              for( Index is=0; is<ns; is++ )
                                { 
                                  const Numeric Z = -X[iv] * iy(iv,is);
                                  diy_dpath[iq](ip  , iv, is) += Z * dkdt1;
                                  diy_dpath[iq](ip+1, iv, is) += Z * dkdt2;
                                }
                              //
                              // The terms associated with Delta-s:
                              if( jacobian_quantities[iq].Subtag() == "HSE on" )
                                {
                                  const Numeric kbar = 0.5 * ( k1 + k2 );
                                  for( Index is=0; is<ns; is++ )
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
            }
          //###################################################################

          // Spectrum at end of ppath step 
          if( stokes_dim == 1 )
            {
              for( Index iv=0; iv<nf; iv++ )  
                { iy(iv,0) = iy(iv,0) * trans_partial(iv,0,0,ip); }
            }
          else
            {
              for( Index iv=0; iv<nf; iv++ )  
                {
                  // Unpolarised absorption:
                  if( is_diagonal( trans_partial(iv,joker,joker,ip) ) )
                    {
                      for( Index is=0; is<ns; is++ )
                        { iy(iv,is) = iy(iv,is) * trans_partial(iv,is,is,ip); }
                    }
                  // The general case:
                  else
                    {
                      // Transmitted term
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
              // Extinction
              if( auxPartExt >= 0  &&  clear2cloudbox[ip] >= 0 ) 
                { 
                  const Index ic = clear2cloudbox[ip];
                  for( Index iv=0; iv<nf; iv++ ) {
                    for( Index is1=0; is1<ns; is1++ ){
                      for( Index is2=0; is2<ns; is2++ ){
                        iy_aux[auxPartExt](iv,is1,is2,ip) = 
                                              pnd_ext_mat(iv,is1,is2,ic); } } }
                }
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

      //### jacobian part #####################################################
      // Map jacobians from ppath to retrieval grids
      // (this operation corresponds to the term Dx_i/Dx)
      if( j_analytical_do )
        { 
          FOR_ANALYTICAL_JACOBIANS_DO( 
            diy_from_path_to_rgrids( diy_dx[iq], jacobian_quantities[iq], 
                               diy_dpath[iq], atmosphere_dim, ppath, ppath_p );
          )
        }
      //#######################################################################
    } // if np>1
}


/* Workspace method: Doxygen documentation will be auto-generated */
void iy_auxFillParticleVariables(
         ArrayOfTensor4&       iy_aux,
   const Index&                basics_checked,
   const Index&                cloudbox_checked,
   const Index&                atmosphere_dim,
   const Index&                cloudbox_on,
   const ArrayOfIndex&         cloudbox_limits,
   const Tensor4&              pnd_field,
   const Matrix&               particle_masses,
   const Ppath&                ppath,
   const ArrayOfString&        iy_aux_vars,
   const Verbosity& )
{
  // Some sizes
  const Index np = ppath.np; 
  const Index naux = iy_aux_vars.nelem();

  // Input checks
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_on )
    throw runtime_error( 
                    "The cloudbox must be activated (cloudbox_on must be 1)" );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );
  if( iy_aux.nelem() != naux )
    throw runtime_error( "*iy_aux_vars* and *iy_aux* must have the same array "
                         "length. (You can not call this WSM before the main "
                         "iy-WSM.)" );

  // Analayse iy_aux_vars
  ArrayOfIndex auxPartCont(0), auxPartContI(0);
  ArrayOfIndex auxPartField(0), auxPartFieldI(0);
  //
  for( Index i=0; i<naux; i++ )
    {
      if( iy_aux_vars[i].substr(0,14) == "Mass content, " )
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
    }

  if( auxPartCont.nelem() + auxPartField.nelem() > 0 )
    {
      // PND along the ppath
      Matrix ppath_pnd( pnd_field.nbooks(), np, 0 );
      //
      for( Index ip=0; ip<np; ip++ )
        {
          Matrix itw( 1, Index(pow(2.0,Numeric(atmosphere_dim))) );

          ArrayOfGridPos gpc_p(1), gpc_lat(1), gpc_lon(1);
          GridPos gp_lat, gp_lon;
          if( atmosphere_dim >= 2 ) { gridpos_copy( gp_lat, ppath.gp_lat[ip] );}
          if( atmosphere_dim == 3 ) { gridpos_copy( gp_lon, ppath.gp_lon[ip] );}
          if( is_gp_inside_cloudbox( ppath.gp_p[ip], gp_lat, gp_lon, 
                                     cloudbox_limits, true, atmosphere_dim ) )
            { 
              interp_cloudfield_gp2itw( itw(0,joker), 
                                        gpc_p[0], gpc_lat[0], gpc_lon[0], 
                                        ppath.gp_p[ip], gp_lat, gp_lon,
                                        atmosphere_dim, cloudbox_limits );
              for( Index i=0; i<pnd_field.nbooks(); i++ )
                {
                  interp_atmfield_by_itw( ppath_pnd(i,ip), atmosphere_dim,
                                          pnd_field(i,joker,joker,joker), 
                                          gpc_p, gpc_lat, gpc_lon, itw );
                }
            }
        }
      
      // Loop ppath steps
      for( Index ip=0; ip<np; ip++ )
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
         ArrayOfVector&              y_aux,
         Matrix&                     jacobian,
   const Index&                      basics_checked,
   const Index&                      atmosphere_dim,
   const Tensor3&                    t_field,
   const Tensor3&                    z_field,
   const Tensor4&                    vmr_field,
   const Index&                      cloudbox_on,
   const Index&                      cloudbox_checked,
   const Index&                      sensor_checked,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const Matrix&                     sensor_pos,
   const Matrix&                     sensor_los,
   const Matrix&                     transmitter_pos,
   const Vector&                     mblock_za_grid,
   const Vector&                     mblock_aa_grid,
   const Index&                      antenna_dim,
   const Sparse&                     sensor_response,
   const Vector&                     sensor_response_f,
   const ArrayOfIndex&               sensor_response_pol,
   const Vector&                     sensor_response_za,
   const Vector&                     sensor_response_aa,
   const Agenda&                     iy_main_agenda,
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

  // Basics, cloudbox, and sensor OK?
  //
  if( !basics_checked )
    throw runtime_error( "The atmosphere and basic control variables must be "
            "flagged to have passed a consistency check (basics_checked=1)." );
  if( !cloudbox_checked )
    throw runtime_error( "The cloudbox must be flagged to have passed a "
                         "consistency check (cloudbox_checked=1)." );
  if( !sensor_checked )
    throw runtime_error( "The sensor variables must be flagged to have passed"
                         "a consistency check (sensor_checked=1)." );

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

  String fail_msg;
  bool failed = false;

#pragma omp parallel for                                          \
  if(!arts_omp_in_parallel() && \
     nmblock>=arts_omp_get_max_threads() && \
     nmblock>=nza)        \
  firstprivate(l_ws, l_jacobian_agenda, l_iy_main_agenda)
  for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
    {
      // Skip remaining iterations if an error occurred
      if (failed) continue;

      try
        {
          // Calculate monochromatic pencil beam data for 1 measurement block
          //
          Vector          iyb, iyb_error, yb(n1y);
          ArrayOfMatrix   diyb_dx;
          //
          iyb_calc(l_ws, iyb, iyb_aux_array[mblock_index], diyb_dx,
                   mblock_index, atmosphere_dim, t_field, z_field, vmr_field,
                   cloudbox_on, stokes_dim, f_grid, sensor_pos, sensor_los,
                   transmitter_pos, mblock_za_grid, mblock_aa_grid, antenna_dim,
                   l_iy_main_agenda, j_analytical_do, jacobian_quantities, 
                   jacobian_indices, iy_aux_vars, verbosity);


          // Apply sensor response matrix on iyb, and put into y
          //
          const Range rowind = get_rowindex_for_mblock(sensor_response,
                                                       mblock_index);
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
                  mult(jacobian(rowind,
                                Range(jacobian_indices[iq][0],
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
        }
        catch (runtime_error e)
        {
#pragma omp critical (yCalc_fail)
            { fail_msg = e.what(); failed = true; }
        }
    }  // End mblock loop

  // Rethrow exception if a runtime error occurred in the mblock loop
  if (failed) throw runtime_error(fail_msg);

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
          const Index row0   = rowind.get_start();

          // The sensor response must be applied in a special way for
          // uncorrelated errors. Schematically: sqrt( H.^2 * y.^2 )
          if( iy_aux_vars[q] == "Error (uncorrelated)" )
            {
              for( Index i=0; i<n1y; i++ )
                {
                  const Index row = row0+i;
                  y_aux[q][row] = 0;
                  for( Index j=0; j<niyb; j++ )
                    { y_aux[q][row] += pow( sensor_response(i,j) * 
                            iyb_aux_array[mblock_index][q][j], (Numeric)2.0 ); }
                  y_aux[q][row] = sqrt( y_aux[q][row] );              
                }
            }
          else
            { mult( y_aux[q][rowind], sensor_response,
                                      iyb_aux_array[mblock_index][q] ); }
        }
    }
}





/* Workspace method: Doxygen documentation will be auto-generated */
void yApplyUnit(
         Vector&         y,
         Matrix&         jacobian,
   const Vector&         y_f,
   const ArrayOfIndex&   y_pol,
   const String&         iy_unit,
   const Verbosity&)
{
  if( iy_unit == "1" )
    { throw runtime_error(
        "No need to use this method with *iy_unit* = \"1\"." ); }

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
  if( iy_unit == "PlanckBT" )
    {
      // Hard to use telescoping here as the data are sorted differently in y
      // and jacobian, than what is expected apply_iy_unit. Copy to temporary
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
                  apply_iy_unit2( J, yv, iy_unit, y_f[i0], 1, i_pol ); 
                  jacobian(ii,joker) = transpose( J(joker,0,joker) );
                }
            }

          // y (must be done last)
          apply_iy_unit( yv, iy_unit, y_f[i0], 1, i_pol );
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
            { apply_iy_unit2( MatrixView( jacobian(i,joker) ), yv, 
                              iy_unit, y_f[i], 1, i_pol ); }

          // y (must be done last)
          apply_iy_unit( yv, iy_unit, y_f[i], 1, i_pol );
          y[i] = yv(0,0);
        }
    }
}
