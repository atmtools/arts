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
extern const String WIND_MAINTAG;



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
   const Index&            atmfields_checked,
   const Index&            atmgeom_checked,
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
   const String&           iy_unit,  
   const Agenda&           iy_main_agenda,
   const Verbosity& )
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


  // iy_transmission is just input and can be left empty for first call
  Tensor3   iy_transmission(0,0,0);

  ArrayOfTensor3 diy_dx;

  iy_main_agendaExecute( ws, iy, iy_aux, ppath, diy_dx, 
                         1, iy_unit, iy_transmission, iy_aux_vars, 
                         cloudbox_on, 0, t_field, 
                         z_field, vmr_field,
                         f_grid, rte_pos, rte_los, rte_pos2,
                         iy_main_agenda );
  
  // Don't allow NaNs (should suffice to check first stokes element)
  for( Index i=0; i<iy.nrows(); i++ )
    { 
      if( isnan(iy(i,0) ) )
        throw runtime_error( "One or several NaNs found in *iy*." );
    }
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
   const Tensor4&                    t_nlte_field,
   const Tensor4&                    vmr_field,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const Tensor3&                    wind_u_field,
   const Tensor3&                    wind_v_field,
   const Tensor3&                    wind_w_field,
   const Tensor3&                    mag_u_field,
   const Tensor3&                    mag_v_field,
   const Tensor3&                    mag_w_field,
   const Index&                      cloudbox_on,
   const String&                     iy_unit,
   const ArrayOfString&              iy_aux_vars,
   const Index&                      jacobian_do,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const Agenda&                     ppath_agenda,
   const Agenda&                     blackbody_radiation_agenda,
   const Agenda&                     propmat_clearsky_agenda,
   const Agenda&                     iy_main_agenda,
   const Agenda&                     iy_space_agenda,
   const Agenda&                     iy_surface_agenda,
   const Agenda&                     iy_cloudbox_agenda,
   const Index&                      iy_agenda_call1,
   const Tensor3&                    iy_transmission,
   const Vector&                     rte_pos,      
   const Vector&                     rte_los,      
   const Vector&                     rte_pos2, 
   const Numeric&                    rte_alonglos_v,      
   const Numeric&                    ppath_lraytrace,
   const Verbosity&                  verbosity )
{
  // Determine propagation path
  //
  ppath_agendaExecute( ws, ppath, ppath_lraytrace, rte_pos, rte_los, rte_pos2, 
                       cloudbox_on, 0, t_field, z_field, vmr_field, f_grid, 
                       ppath_agenda );
  //
  if( !iy_agenda_call1 && ppath_what_background( ppath ) == 2  )
    throw runtime_error( "A secondary propagation path starting at the "
                         "surface and is going directly into the surface "
                         "is found. This is not allowed." );

  // Some basic sizes
  //
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index nq = jacobian_quantities.nelem();
  
  // Helper variable
  const PropmatPartialsData ppd(jacobian_quantities);

  //### jacobian part #########################################################
  // Initialise analytical jacobians (diy_dx and help variables)
  //
  Index j_analytical_do = 0;
  ArrayOfTensor3  diy_dpath; 
  ArrayOfIndex    jac_species_i(0), jac_is_t(0), jac_wind_i(0), jac_mag_i(0), jac_other(0); 
  //
  if( jacobian_do ) { FOR_ANALYTICAL_JACOBIANS_DO( j_analytical_do = 1; ) }
  //
  if( !j_analytical_do )
    { diy_dx.resize( 0 ); }
  else 
    {
      diy_dpath.resize( nq ); 
      jac_species_i.resize( nq ); 
      jac_is_t.resize( nq ); 
      jac_wind_i.resize( nq );  
      jac_mag_i.resize( nq ); 
     //
      FOR_ANALYTICAL_JACOBIANS_DO( 
        diy_dpath[iq].resize( np, nf, ns ); 
        diy_dpath[iq] = 0.0;
      )
      get_pointers_for_analytical_jacobians( jac_species_i, jac_is_t, 
                                             jac_wind_i, jac_mag_i,
                                             jacobian_quantities, abs_species );
      
      jac_other.resize(jac_is_t.nelem());
      FOR_ANALYTICAL_JACOBIANS_DO( jac_other[iq] = ppd.is_this_linetype(iq)?JAC_IS_OTHER:JAC_IS_NONE; )
      
      if( iy_agenda_call1 )
        {
          diy_dx.resize( nq ); 
          //
          FOR_ANALYTICAL_JACOBIANS_DO( diy_dx[iq].resize( 
                  jacobian_indices[iq][1]-jacobian_indices[iq][0]+1, nf, ns ); 
            diy_dx[iq] = 0.0;
           )
        }
    } 
  //###########################################################################

  
  // Set up variable with index of species where we want abs_per_species.
  // This variable can below be extended in iy_aux part.
  //
  ArrayOfIndex iaps(0);
  //
  for( Index i=0; i<jac_species_i.nelem(); i++ )
    {
      if( jac_species_i[i] >= 0 )
        { iaps.push_back( jac_species_i[i] ); }
    }

  
  //=== iy_aux part ===========================================================
  Index auxPressure    = -1,
        auxTemperature = -1,
        auxAbsSum      = -1,
        auxBackground  = -1,
        auxIy          = -1,
        auxTrans       = -1,
        auxOptDepth    = -1;
  ArrayOfIndex auxAbsSpecies(0), auxAbsIsp(0);
  ArrayOfIndex auxVmrSpecies(0), auxVmrIsp(0);
  //
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
            const Index ihit = find_first( iaps, ispecies );
            if( ihit >= 0 )
              { auxAbsIsp.push_back( ihit ); }
            else
              { 
                iaps.push_back(ispecies); 
                auxAbsIsp.push_back( iaps.nelem()-1 ); 
              }
            iy_aux[i].resize( nf, ns, ns, np );               
          }
        else if( iy_aux_vars[i] == "Radiative background" )
          { auxBackground = i;   iy_aux[i].resize( nf, 1, 1, 1 ); }
        else if( iy_aux_vars[i] == "iy"   &&  auxIy < 0 )
          { auxIy = i;           iy_aux[i].resize( nf, ns, 1, np ); }
        else if( iy_aux_vars[i] == "Transmission"   &&  auxTrans < 0 )
          { auxTrans = i;        iy_aux[i].resize( nf, ns, ns, np ); }
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


  // Get atmospheric and attenuation quantities for each ppath point/step
  //
  // "atmvars"
  Vector    ppath_p, ppath_t;
  Matrix    ppath_vmr, ppath_wind, ppath_mag, ppath_f, ppath_t_nlte;
  // Attenuation vars
  Tensor4   ppath_ext;
  Tensor5   abs_per_species, dppath_ext_dx, dtrans_partial_dx_above, dtrans_partial_dx_below;
  Tensor4   trans_partial, trans_cumulat;
  Tensor3   ppath_nlte_source;
  Matrix    ppath_blackrad;
  Vector    scalar_tau;
  ArrayOfIndex   lte;
  ArrayOfArrayOfIndex  extmat_case;
  //
  if( np > 1 )
    {
      get_ppath_atmvars(  ppath_p, ppath_t, ppath_t_nlte, ppath_vmr,
                          ppath_wind, ppath_mag, 
                          ppath, atmosphere_dim, p_grid, t_field, t_nlte_field, vmr_field,
                          wind_u_field, wind_v_field, wind_w_field,
                          mag_u_field, mag_v_field, mag_w_field );      
      get_ppath_f(        ppath_f, ppath, f_grid,  atmosphere_dim, 
                          rte_alonglos_v, ppath_wind );
      get_ppath_pmat(     ws, ppath_ext, ppath_nlte_source, lte, abs_per_species, 
                          dppath_ext_dx,
                          propmat_clearsky_agenda, jacobian_quantities, ppath, 
                          ppath_p, ppath_t, ppath_t_nlte, ppath_vmr, ppath_f, 
                          ppath_mag, f_grid, stokes_dim, iaps );
      
      // Perhaps put all in a function "adapt_dppath_ext_dx"?
      for(Index iq=0;(iq<nq)&&jacobian_do;iq++)
        {
            if( jac_wind_i[iq] == JAC_IS_WIND_U_FROM_PROPMAT ||
                jac_wind_i[iq] == JAC_IS_WIND_V_FROM_PROPMAT ||
                jac_wind_i[iq] == JAC_IS_WIND_W_FROM_PROPMAT ||
                jac_wind_i[iq] == JAC_IS_WIND_ABS_FROM_PROPMAT )
            {
                Index component=-1;
                if(jac_wind_i[iq] == JAC_IS_WIND_ABS_FROM_PROPMAT)
                    component = 0;
                else if(jac_wind_i[iq] == JAC_IS_WIND_U_FROM_PROPMAT)
                    component=1;
                else if(jac_wind_i[iq] == JAC_IS_WIND_V_FROM_PROPMAT)
                    component=2;
                else if(jac_wind_i[iq] == JAC_IS_WIND_W_FROM_PROPMAT)
                    component=3;
                else 
                    throw std::runtime_error("To developer:  You have changed wind jacobians" 
                    "in an incompatible manners with some other code.");
                
                const Index this_ppd_q = ppd.this_jq_index(iq);
                
                Matrix dWdx;
                get_ppath_f_partials(dWdx, component, ppath, f_grid,  atmosphere_dim, ppath_wind );
                
                for(Index is1=0;is1<stokes_dim;is1++)
                    for(Index is2=0;is2<stokes_dim;is2++)
                        dppath_ext_dx(this_ppd_q,joker,is1,is2,joker)*=dWdx;
            }
        }
      
      if(dppath_ext_dx.empty())
        get_ppath_trans(    trans_partial, extmat_case, trans_cumulat, 
                            scalar_tau, ppath, ppath_ext, f_grid, stokes_dim );
      else
        get_ppath_trans_and_dppath_trans_dx(
            trans_partial, dtrans_partial_dx_above, dtrans_partial_dx_below, 
            extmat_case, trans_cumulat, scalar_tau, ppath, ppath_ext, 
            dppath_ext_dx, jacobian_quantities, f_grid, stokes_dim );
      get_ppath_blackrad( ws, ppath_blackrad, blackbody_radiation_agenda, 
                          ppath, ppath_t, ppath_f );
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
                        cloudbox_on, stokes_dim, f_grid, iy_unit,
                        iy_main_agenda, iy_space_agenda, iy_surface_agenda, 
                        iy_cloudbox_agenda, verbosity );


  //=== iy_aux part ===========================================================
  // Fill parts of iy_aux that are defined even for np=1.
  // Radiative background
  if( auxBackground >= 0 ) 
    { iy_aux[auxBackground](joker,0,0,0) = (Numeric)min( (Index)2,
                                              ppath_what_background(ppath)-1); }
  // Radiance 
  if( auxIy >= 0 ) 
    { iy_aux[auxIy](joker,joker,0,np-1) = iy; }
  // Transmission variables
  if( auxTrans >= 0 ) 
    { 
      if( np == 1 )
        { for( Index iv=0; iv<nf; iv++ ) {
            id_mat( iy_aux[auxTrans](iv,joker,joker,0) ); } }
      else
        { iy_aux[auxTrans] = trans_cumulat; }
    }
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
      //
      const Numeric   dt = 0.1;     // Temperature disturbance, K
      const Numeric   dw = 5;       // Wind disturbance, m/s
      const Numeric   dm = 0.1e-6;  // Magnetic field disturbance, T
            Tensor4   ppath_ext_dt;
            ArrayOfTensor4 ppath_ext_dw(3), ppath_ext_dm(3);
            Matrix    ppath_blackrad_dt;
      //
      if( j_analytical_do )
        { 
          // Determine if temperature is among the analytical jac. quantities.
          // If yes, get emission and absorption for disturbed temperature
          // Same for wind and mag. field, but disturb only absorption
          for( Index iq=0; iq<jac_is_t.nelem(); iq++ )
            { 
              // Temperatures
              if( jac_is_t[iq] ) 
                { 
                  Tensor5 dummy_abs_per_species,dummy_dppath_ext_dx;
                  Tensor3 dummy_ppath_nlte_source;
                  ArrayOfIndex dummy_lte;
                  Vector t2 = ppath_t;   t2 += dt;
                  if(jac_is_t[iq]==JAC_IS_T_SEMI_ANALYTIC)
                    get_ppath_pmat( ws, ppath_ext_dt, dummy_ppath_nlte_source,
                                    dummy_lte, dummy_abs_per_species, dummy_dppath_ext_dx,
                                    propmat_clearsky_agenda, ArrayOfRetrievalQuantity(0), ppath, ppath_p,
                                    t2, ppath_t_nlte, ppath_vmr, ppath_f,
                                    ppath_mag, f_grid, 
                                    stokes_dim, ArrayOfIndex(0) );
                  // Use this method or dplanck_dt ?
                  get_ppath_blackrad( ws, ppath_blackrad_dt, blackbody_radiation_agenda,
                                      ppath, t2, ppath_f );
                }
              // Winds
              else if( jac_wind_i[iq] == JAC_IS_WIND_U_SEMI_ANALYTIC ||
                       jac_wind_i[iq] == JAC_IS_WIND_V_SEMI_ANALYTIC ||
                       jac_wind_i[iq] == JAC_IS_WIND_W_SEMI_ANALYTIC )
              {
                  Index this_field=-1;
                  if(jac_wind_i[iq] == JAC_IS_WIND_U_SEMI_ANALYTIC)
                      this_field=0;
                  else if(jac_wind_i[iq] == JAC_IS_WIND_V_SEMI_ANALYTIC)
                      this_field=1;
                  else if(jac_wind_i[iq] == JAC_IS_WIND_W_SEMI_ANALYTIC)
                      this_field=2;
                  else 
                      throw std::runtime_error("To developer:  You have changed wind jacobians" 
                      "in an incompatible manners with some other code.");
                  
                  Tensor5 dummy_abs_per_species,dummy_dppath_ext_dx;
                  Tensor3 dummy_ppath_nlte_source;
                  ArrayOfIndex dummy_lte;
                  Matrix f2, w2 = ppath_wind;   w2(this_field,joker) += dw;
                  get_ppath_f(    f2, ppath, f_grid,  atmosphere_dim, 
                                  rte_alonglos_v, w2 );
                  get_ppath_pmat( ws, ppath_ext_dw[this_field], 
                                  dummy_ppath_nlte_source,
                                  dummy_lte, dummy_abs_per_species, dummy_dppath_ext_dx,
                                  propmat_clearsky_agenda, ArrayOfRetrievalQuantity(0), ppath, ppath_p, 
                                  ppath_t, ppath_t_nlte, ppath_vmr, f2, 
                                  ppath_mag, f_grid,
                                  stokes_dim, ArrayOfIndex(0) );
                }
              // Magnetic field
              else if( jac_mag_i[iq] == JAC_IS_MAG_U_SEMI_ANALYTIC ||
                       jac_mag_i[iq] == JAC_IS_MAG_V_SEMI_ANALYTIC ||
                       jac_mag_i[iq] == JAC_IS_MAG_W_SEMI_ANALYTIC )
                {
                    Index this_field=-1;
                    if(jac_wind_i[iq] == JAC_IS_MAG_U_SEMI_ANALYTIC)
                        this_field=0;
                    else if(jac_wind_i[iq] == JAC_IS_MAG_V_SEMI_ANALYTIC)
                        this_field=1;
                    else if(jac_wind_i[iq] == JAC_IS_MAG_W_SEMI_ANALYTIC)
                        this_field=2;
                    else 
                        throw std::runtime_error("To developer:  You have changed magnetic jacobians" 
                        "in an incompatible manners with some other code.");
                    
                  Tensor5 dummy_abs_per_species, dummy_dppath_ext_dx;
                  Tensor3 dummy_ppath_nlte_source;
                  ArrayOfIndex dummy_lte;
                  Matrix m2 = ppath_mag;   m2(this_field,joker) += dm;
                  get_ppath_pmat( ws, ppath_ext_dm[this_field], 
                                  dummy_ppath_nlte_source,
                                  dummy_lte, dummy_abs_per_species, dummy_dppath_ext_dx,
                                  propmat_clearsky_agenda, ArrayOfRetrievalQuantity(0), ppath, ppath_p, 
                                  ppath_t, ppath_t_nlte, ppath_vmr, ppath_f, 
                                  m2, f_grid,
                                  stokes_dim, ArrayOfIndex(0) );
                }
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
                                             ppath_ext(iv,is1,is2,np-1); } } } }
      for( Index j=0; j<auxAbsSpecies.nelem(); j++ )
        { for( Index iv=0; iv<nf; iv++ ) {
            for( Index is1=0; is1<ns; is1++ ){
              for( Index is2=0; is2<ns; is2++ ){
                iy_aux[auxAbsSpecies[j]](iv,is1,is2,np-1) = 
                          abs_per_species(auxAbsIsp[j],iv,is1,is2,np-1); } } } }
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

          // Extra variables for non-LTE
          //
          const bool nonlte = lte[ip]==0 || lte[ip+1]==0; 
          //
          Matrix sourcebar(0,0);
          Tensor3 extbar(0,0,0);
          //
          if( nonlte )
            {
              if( jacobian_do )
                throw runtime_error( "iyEmissionStandard can so far not handle "
                                     "non-LTE together with Jacobians." );
              
              sourcebar.resize( nf, stokes_dim );
              extbar.resize( nf, stokes_dim, stokes_dim );
              for( Index iv=0; iv<nf; iv++ )  
                { 
                  for( Index is1=0; is1<stokes_dim; is1++ )  
                    {
                      sourcebar(iv,is1) = 0.5 * ( ppath_nlte_source(iv,is1,ip) + 
                                                  ppath_nlte_source(iv,is1,ip+1) );
                      for( Index is2=0; is2<stokes_dim; is2++ )  
                        { extbar(iv,is1,is2) = 0.5 * ( 
                                               ppath_ext(iv,is1,is2,ip) + 
                                               ppath_ext(iv,is1,is2,ip+1) ); }
                    }
                }
              
            }

          //### jacobian part #################################################
          if( j_analytical_do )
            {
              // Calculate (si-bi)
              Matrix sibi(nf,ns);
              for( Index iv=0; iv<nf; iv++ )
                {
                  sibi(iv,0) = iy(iv,0) - bbar[iv];
                  for( Index is=1; is<ns; is++ )
                    { sibi(iv,is) = iy(iv,is); }
                }

              // Loop quantities
              //
              Index iiaps = -1;   // Index with respect to abs_per_species.
              //                     Jacobian species stored first and in order
              for( Index iq=0; iq<nq; iq++ ) 
                {
                  if( jacobian_quantities[iq].Analytical() )
                    {
                      //- Absorbing species -----------------------------------
                      if( jac_species_i[iq] >= 0 )
                        {
                          bool from_propmat = jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG;
                          const Index this_ppd_q = from_propmat?ppd.this_jq_index(iq):-1;
                          
                          // Index with respect to abs_per_species and ppath_vmr
                          iiaps += 1;
                          const Index isp = jac_species_i[iq];

                          // Scaling factors to handle retrieval unit
                          Numeric unitscf1, unitscf2;
                          if(!from_propmat)
                          {
                            vmrunitscf( unitscf1, 
                                        jacobian_quantities[iq].Mode(), 
                                        ppath_vmr(isp,ip), ppath_p[ip], 
                                        ppath_t[ip] );
                            vmrunitscf( unitscf2, 
                                        jacobian_quantities[iq].Mode(), 
                                        ppath_vmr(isp,ip+1), ppath_p[ip+1], 
                                        ppath_t[ip+1] );
                          }
                          /*else implicitly "vmr" is the mode*/
                          

                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // Diagonal transmission matrix
                              if( extmat_case[ip][iv] == 1 )
                                {
                                  const Numeric dkdn1 = (from_propmat?
                                    dppath_ext_dx(this_ppd_q,iv,0,0,ip):
                                    unitscf1 * abs_per_species(iiaps,iv,0,0,ip  ));
                                  const Numeric dkdn2 = (from_propmat?
                                    dppath_ext_dx(this_ppd_q,iv,0,0,ip+1):
                                    unitscf2 * abs_per_species(iiaps,iv,0,0,ip+1));
                                    
                                  const Numeric x = -0.5 * ppath.lstep[ip] * 
                                                    trans_cumulat(iv,0,0,ip+1);
                                  const Numeric y = x * sibi(iv,0);
                                  // Stokes component 1
                                  diy_dpath[iq](ip  ,iv,0) += y * dkdn1;
                                  diy_dpath[iq](ip+1,iv,0) += y * dkdn2;
                                  // Higher stokes components
                                  for( Index is=1; is<ns; is++ )
                                    { 
                                      const Numeric z = x * iy(iv,is); 
                                      diy_dpath[iq](ip  ,iv,is) += z * dkdn1;
                                      diy_dpath[iq](ip+1,iv,is) += z * dkdn2;
                                    }
                                }

                              // General case
                              else
                                { 
                                  // Size of disturbance, a relative number
                                  const Numeric dd = 1e-3;
                                  Matrix ext_mat(ns,ns), dtdx(ns,ns);
                                  
                                  // Disturb for ip
                                  if(!from_propmat)
                                  {
                                    //
                                    for( Index is1=0; is1<ns; is1++ ) {
                                        for( Index is2=0; is2<ns; is2++ ) {
                                        ext_mat(is1,is2) = 0.5 * ( dd *
                                        abs_per_species(iiaps,iv,is1,is2,ip  ) +
                                                    ppath_ext(iv,is1,is2,ip) +
                                                    ppath_ext(iv,is1,is2,ip+1) );
                                        } }
                                    ext2trans( dtdx, extmat_case[ip][iv], 
                                                ext_mat, ppath.lstep[ip] ); 
                                    for( Index is1=0; is1<ns; is1++ ) {
                                        for( Index is2=0; is2<ns; is2++ ) {
                                        dtdx(is1,is2) = (unitscf1/dd) * 
                                                ( dtdx(is1,is2) -
                                                    trans_partial(iv,is1,is2,ip) );
                                        } }
                                  }
                                  else
                                    dtdx = dtrans_partial_dx_below(this_ppd_q,iv,joker,joker,ip);
                                
                                  Vector x(ns), y(ns);
                                  mult( x, dtdx, sibi(iv,joker) );
                                  mult( y, trans_cumulat(iv,joker,joker,ip), x );
                                  diy_dpath[iq](ip,iv,joker) += y;
                                  //
                                  
                                  // Disturb for ip+1
                                  if(!from_propmat)
                                  {
                                      //
                                  for( Index is1=0; is1<ns; is1++ ) {
                                    for( Index is2=0; is2<ns; is2++ ) {
                                      ext_mat(is1,is2) = 0.5 * ( dd *
                                       abs_per_species(iiaps,iv,is1,is2,ip+1) +
                                                  ppath_ext(iv,is1,is2,ip) +
                                                  ppath_ext(iv,is1,is2,ip+1) );
                                    } }
                                  ext2trans( dtdx, extmat_case[ip][iv], 
                                             ext_mat, ppath.lstep[ip] ); 
                                  for( Index is1=0; is1<ns; is1++ ) {
                                    for( Index is2=0; is2<ns; is2++ ) {
                                      dtdx(is1,is2) = (unitscf2/dd) * 
                                              ( dtdx(is1,is2) -
                                                trans_partial(iv,is1,is2,ip) );
                                    } }
                                  }
                                  else
                                      dtdx = dtrans_partial_dx_above(this_ppd_q,iv,joker,joker,ip);
                                  
                                  mult( x, dtdx, sibi(iv,joker) );
                                  mult( y, trans_cumulat(iv,joker,joker,ip), x );
                                  diy_dpath[iq](ip+1,iv,joker) += y;
                                }
                            }
                        }

                      //- Winds and magnetic field -----------------------------------
                      else if( jac_wind_i[iq] || jac_mag_i[iq] || jac_other[iq] )
                        {
                            bool from_propmat = jacobian_quantities[iq].SubSubtag() == PROPMAT_SUBSUBTAG;
                            const Index this_ppd_q = from_propmat?ppd.this_jq_index(iq):-1;
                            
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // Create pointer to disturbed extinction to use
                              // Wind v-component is first guess.
                              Tensor4* ppath_ext2 = &ppath_ext_dw[1];
                              Numeric dd = dw;
                              
                              if(!from_propmat)
                              {
                                if( jac_wind_i[iq] == JAC_IS_WIND_U_SEMI_ANALYTIC )
                                    { ppath_ext2 = &ppath_ext_dw[0]; }
                                else if( jac_wind_i[iq] == JAC_IS_WIND_W_SEMI_ANALYTIC )
                                    { ppath_ext2 = &ppath_ext_dw[2]; }
                                else if( jac_mag_i[iq] == JAC_IS_MAG_U_SEMI_ANALYTIC )
                                    { ppath_ext2 = &ppath_ext_dm[0]; dd = dm; }
                                else if( jac_mag_i[iq] == JAC_IS_MAG_V_SEMI_ANALYTIC )
                                    { ppath_ext2 = &ppath_ext_dm[1]; dd = dm; }
                                else if( jac_mag_i[iq] == JAC_IS_MAG_W_SEMI_ANALYTIC )
                                    { ppath_ext2 = &ppath_ext_dm[2]; dd = dm; }
                              }
                              // Diagonal transmission matrix
                              if( extmat_case[ip][iv] == 1 )
                                {
                                    
                                  const Numeric dkdx1 = from_propmat?dppath_ext_dx(this_ppd_q,iv,0,0,ip):
                                                   (1/dd) * ( 
                                                   (*ppath_ext2)(iv,0,0,ip  ) -
                                                      ppath_ext(iv,0,0,ip  ) );
                                  const Numeric dkdx2 = from_propmat?dppath_ext_dx(this_ppd_q,iv,0,0,ip+1):
                                                   (1/dd) * ( 
                                                   (*ppath_ext2)(iv,0,0,ip+1) -
                                                      ppath_ext(iv,0,0,ip+1) );
                                  const Numeric x = -0.5 * ppath.lstep[ip] * 
                                                 trans_cumulat(iv,0,0,ip+1);
                                  const Numeric y = x * sibi(iv,0);
                                  // Stokes component 1
                                  diy_dpath[iq](ip  ,iv,0) += y * dkdx1;
                                  diy_dpath[iq](ip+1,iv,0) += y * dkdx2;
                                  // Higher stokes components
                                  for( Index is=1; is<ns; is++ )
                                    { 
                                      const Numeric z = x * iy(iv,is); 
                                      diy_dpath[iq](ip  ,iv,is) += z * dkdx1;
                                      diy_dpath[iq](ip+1,iv,is) += z * dkdx2;
                                    }
                                }

                              // General case
                              else
                                { 
                                  Matrix ext_mat(ns,ns), dtdx(ns,ns);
                                  //
                                  
                                  if(!from_propmat)
                                  {
                                    // Disturb for ip
                                    for( Index is1=0; is1<ns; is1++ ) {
                                        for( Index is2=0; is2<ns; is2++ ) {
                                        ext_mat(is1,is2) = 0.5 * (
                                                (*ppath_ext2)(iv,is1,is2,ip  ) +
                                                    ppath_ext(iv,is1,is2,ip+1) );
                                        } }
                                    ext2trans( dtdx, extmat_case[ip][iv], 
                                                ext_mat, ppath.lstep[ip] ); 
                                    for( Index is1=0; is1<ns; is1++ ) {
                                        for( Index is2=0; is2<ns; is2++ ) {
                                        dtdx(is1,is2) = (1/dd) * 
                                                ( dtdx(is1,is2) -
                                                    trans_partial(iv,is1,is2,ip) );
                                        } }
                                  }
                                  else
                                      dtdx = dtrans_partial_dx_below(this_ppd_q,iv,joker,joker,ip);
                                      
                                  Vector x(ns), y(ns);
                                  mult( x, dtdx, sibi(iv,joker) );
                                  mult( y, trans_cumulat(iv,joker,joker,ip), x );
                                  diy_dpath[iq](ip,iv,joker) += y;
                                  //
                                  
                                  if(!from_propmat)
                                  {
                                    // Disturb for ip+1
                                    for( Index is1=0; is1<ns; is1++ ) {
                                        for( Index is2=0; is2<ns; is2++ ) {
                                        ext_mat(is1,is2) = 0.5 * (
                                                    ppath_ext(iv,is1,is2,ip  ) +
                                                (*ppath_ext2)(iv,is1,is2,ip+1) );
                                        } }
                                    ext2trans( dtdx, extmat_case[ip][iv], 
                                                ext_mat, ppath.lstep[ip] ); 
                                    for( Index is1=0; is1<ns; is1++ ) {
                                        for( Index is2=0; is2<ns; is2++ ) {
                                        dtdx(is1,is2) = (1/dd) * 
                                                ( dtdx(is1,is2) -
                                                    trans_partial(iv,is1,is2,ip) );
                                        } }
                                  }
                                  else
                                      dtdx = dtrans_partial_dx_above(this_ppd_q,iv,joker,joker,ip);
                                      
                                  mult( x, dtdx, sibi(iv,joker) );
                                  mult( y, trans_cumulat(iv,joker,joker,ip), x );
                                  diy_dpath[iq](ip+1,iv,joker) += y;
                                }
                            }
                        }
                      
                      //- Temperature -----------------------------------------
                      else if( jac_is_t[iq] )
                        {
                          const bool from_propmat = jac_is_t[iq]==JAC_IS_T_FROM_PROPMAT;
                          const Index this_ppd_q = from_propmat?ppd.this_jq_index(iq):-1;
                          for( Index iv=0; iv<nf; iv++ )
                            {
                              // Diagonal transmission matrix
                              if( extmat_case[ip][iv] == 1 )
                                {
                                  const Numeric dkdt1 = from_propmat?dppath_ext_dx(this_ppd_q,iv,0,0,ip):
                                                      1/dt * (
                                                      ppath_ext_dt(iv,0,0,ip  ) -
                                                      ppath_ext(iv,0,0,ip  ) );
                                  const Numeric dkdt2 = from_propmat?dppath_ext_dx(this_ppd_q,iv,0,0,ip+1):
                                                      1/dt * (
                                                      ppath_ext_dt(iv,0,0,ip+1) - 
                                                      ppath_ext(iv,0,0,ip+1) );
                                  const Numeric x = -0.5 * ppath.lstep[ip] * 
                                                 trans_cumulat(iv,0,0,ip+1);
                                  const Numeric y = x * sibi(iv,0);
                                  // Stokes 1:
                                  diy_dpath[iq](ip  ,iv,0) += y * dkdt1;
                                  diy_dpath[iq](ip+1,iv,0) += y * dkdt2;
                                  // Higher Stokes
                                  for( Index is=1; is<ns; is++ )
                                    { 
                                      const Numeric z = x * iy(iv,is);
                                      diy_dpath[iq](ip  ,iv,is) += z * dkdt1;
                                      diy_dpath[iq](ip+1,iv,is) += z * dkdt2;
                                    }
                                  //
                                  // The terms associated with B-bar:
                                  const Numeric v = 0.5 * 
                                                trans_cumulat(iv,0,0,ip) *
                                             ( 1.0 - trans_partial(iv,0,0,ip));
                                  diy_dpath[iq](ip  ,iv,0) += v/dt * (
                                                     ppath_blackrad_dt(iv,ip) -
                                                     ppath_blackrad(iv,ip) );
                                  diy_dpath[iq](ip+1,iv,0) += v/dt * (
                                                     ppath_blackrad_dt(iv,ip+1) -
                                                     ppath_blackrad(iv,ip+1) );
                                  // Zero for higher Stokes
                                  //
                                  // The terms associated with Delta-s:
                                  if( jacobian_quantities[iq].Subtag() == 
                                                                     "HSE on" )
                                    {
                                      // Stokes 1:
                                      const Numeric kbar = 0.5 * ( 
                                                      ppath_ext(iv,0,0,ip  ) +
                                                      ppath_ext(iv,0,0,ip+1) );
                                      diy_dpath[iq](ip  ,iv,0) += y * kbar /
                                                                  ppath_t[ip];
                                      diy_dpath[iq](ip+1,iv,0) += y * kbar /
                                                                 ppath_t[ip+1];
                                      // Higher Stokes
                                      for( Index is=1; is<ns; is++ )
                                        { 
                                          const Numeric z = x * iy(iv,is);
                                          diy_dpath[iq](ip  ,iv,is) += 
                                                      z * kbar / ppath_t[ip];
                                          diy_dpath[iq](ip+1,iv,is) += 
                                                      z * kbar / ppath_t[ip+1];
                                        }
                                    } //hse
                                }
                              // General case
                              else
                                { 
                                  Matrix ext_mat(ns,ns), dtdx(ns,ns);
                                  if(!from_propmat)
                                  {
                                    // Disturb for ip
                                    
                                    //
                                    for( Index is1=0; is1<ns; is1++ ) {
                                        for( Index is2=0; is2<ns; is2++ ) {
                                        ext_mat(is1,is2) = 0.5 * (
                                                    ppath_ext_dt(iv,is1,is2,ip  ) +
                                                    ppath_ext(iv,is1,is2,ip+1) );
                                        } }
                                    ext2trans( dtdx, extmat_case[ip][iv], 
                                                ext_mat, ppath.lstep[ip] ); 
                                    for( Index is1=0; is1<ns; is1++ ) {
                                        for( Index is2=0; is2<ns; is2++ ) {
                                        dtdx(is1,is2) = (1/dt) * 
                                                ( dtdx(is1,is2) -
                                                    trans_partial(iv,is1,is2,ip) );
                                        } }
                                  }
                                  else
                                      dtdx = dtrans_partial_dx_below(this_ppd_q,iv,joker,joker,ip);
                                  
                                  Vector x(ns), y(ns);
                                  mult( x, dtdx, sibi(iv,joker) );
                                  mult( y, trans_cumulat(iv,joker,joker,ip), 
                                                                           x );
                                  diy_dpath[iq](ip,iv,joker) += y;
                                  
                                  if(!from_propmat)
                                  {
                                    //
                                    // Disturb for ip+1
                                    for( Index is1=0; is1<ns; is1++ ) {
                                        for( Index is2=0; is2<ns; is2++ ) {
                                        ext_mat(is1,is2) = 0.5 * (
                                                    ppath_ext(iv,is1,is2,ip  ) +
                                                    ppath_ext_dt(iv,is1,is2,ip+1) );
                                        } }
                                    ext2trans( dtdx, extmat_case[ip][iv], 
                                                ext_mat, ppath.lstep[ip] ); 
                                    for( Index is1=0; is1<ns; is1++ ) {
                                        for( Index is2=0; is2<ns; is2++ ) {
                                        dtdx(is1,is2) = (1/dt) * 
                                                ( dtdx(is1,is2) -
                                                    trans_partial(iv,is1,is2,ip) );
                                        } }
                                  }
                                  else
                                    dtdx = dtrans_partial_dx_above(this_ppd_q,iv,joker,joker,ip);
                                  
                                  mult( x, dtdx, sibi(iv,joker) );
                                  mult( y, trans_cumulat(iv,joker,joker,ip), 
                                                                           x );
                                  diy_dpath[iq](ip+1,iv,joker) += y; 
                                  //
                                  // The terms associated with B-bar:
                                  const Numeric v = 
                                             ( 1.0 - trans_partial(iv,0,0,ip));
                                  Numeric dbdt = 0.5/dt * ( ppath_blackrad_dt(iv,ip) -
                                                       ppath_blackrad(iv,ip) );
                                  x[0] = v * dbdt;
                                  for( Index is=1; is<ns; is++ ) 
                                    { x[is] = -trans_partial(iv,is,0,ip)*dbdt; }
                                  mult( y, trans_cumulat(iv,joker,joker,ip), 
                                                                           x );
                                  diy_dpath[iq](ip,iv,joker) += y; 
                                  // Some for ip+1
                                  dbdt = 0.5/dt * ( ppath_blackrad_dt(iv,ip+1) -
                                                    ppath_blackrad(iv,ip+1) );
                                  x[0] = v * dbdt;
                                  for( Index is=1; is<ns; is++ ) 
                                    { x[is] = -trans_partial(iv,is,0,ip)*dbdt; }
                                  mult( y, trans_cumulat(iv,joker,joker,ip), x );
                                  diy_dpath[iq](ip+1,iv,joker) += y; 
                                  //
                                  // The terms associated with Delta-s:
                                  if( jacobian_quantities[iq].Subtag() == 
                                                                     "HSE on" )
                                    {
                                      for( Index is1=0; is1<ns; is1++ ) {
                                        for( Index is2=0; is2<ns; is2++ ) {
                                          ext_mat(is1,is2) = 0.5 * (
                                                   ppath_ext(iv,is1,is2,ip  ) +
                                                   ppath_ext(iv,is1,is2,ip+1) );
                                        } }
                                      // dl for disturbed tbar
                                      const Numeric tbar = 0.5 * (
                                                  ppath_t[ip] + ppath_t[ip+1] );
                                      const Numeric dl = ppath.lstep[ip] *
                                                            ( 1 + dt/tbar );
                                      ext2trans( dtdx, extmat_case[ip][iv], 
                                                 ext_mat, dl ); 
                                      for( Index is1=0; is1<ns; is1++ ) {
                                        for( Index is2=0; is2<ns; is2++ ) {
                                          dtdx(is1,is2) = (1/dt) * 
                                            ( dtdx(is1,is2) -
                                                trans_partial(iv,is1,is2,ip) );
                                        } }
                                      mult( x, dtdx, sibi(iv,joker) );
                                      mult( y, trans_cumulat(iv,joker,joker,ip),
                                                                            x );
                                      // Contribution shared between the two
                                      // points  and is proportional to 1/t
                                      // See also AUG.
                                      for( Index is=0; is<ns; is++ ) 
                                        {
                                          diy_dpath[iq](ip  ,iv,is) += y[is] *
                                                    0.5 * tbar / ppath_t[ip];
                                          diy_dpath[iq](ip+1,iv,is) += y[is] *
                                                    0.5 * tbar / ppath_t[ip+1];
                                        }
                                    } // HSE
                                } // General case
                            } // Frequency loop 
                        } // Temperature
                    } // if this analytical
                } // for iq
            } // if any analytical
          //###################################################################


          // Spectrum at end of ppath step
          emission_rtstep( iy, stokes_dim, bbar, extmat_case[ip],
                           trans_partial(joker,joker,joker,ip),
                           nonlte, extbar, sourcebar );


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
                                               ppath_ext(iv,is1,is2,ip); } } } }
          for( Index j=0; j<auxAbsSpecies.nelem(); j++ )
            { for( Index iv=0; iv<nf; iv++ ) {
                for( Index is1=0; is1<ns; is1++ ){
                  for( Index is2=0; is2<ns; is2++ ){
                    iy_aux[auxAbsSpecies[j]](iv,is1,is2,ip) = 
                           abs_per_species(auxAbsIsp[j],iv,is1,is2,ip); } } } }
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
      ArrayOfIndex i_pol(ns);
      for( Index is=0; is<ns; is++ )
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
   const String&           iy_unit,  
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
                            1, iy_unit, iy_transmission, iy_aux_vars, 
                            cloudbox_on, jacobian_do, t_field, z_field, 
                            vmr_field, Vector(1,f_grid[i]),
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
   const Vector&                     refellipsoid,
   const Matrix&                     z_surface,
   const Index&                      cloudbox_on,
   const ArrayOfIndex&               cloudbox_limits,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const ArrayOfArrayOfSingleScatteringData&   scat_data,
   const Agenda&                     iy_space_agenda,
   const Agenda&                     surface_rtprop_agenda,
   const Agenda&                     propmat_clearsky_agenda, 
   const Agenda&                     ppath_step_agenda, 
   const Numeric&                    ppath_lraytrace, 
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
  Agenda l_propmat_clearsky_agenda (propmat_clearsky_agenda);
  Agenda l_surface_rtprop_agenda (surface_rtprop_agenda);

  String fail_msg;
  bool failed = false;

  if (nf)
#pragma omp parallel for                   \
  if (!arts_omp_in_parallel() && nf > 1)   \
  firstprivate(l_ws, l_ppath_step_agenda, l_iy_space_agenda, \
               l_propmat_clearsky_agenda, l_surface_rtprop_agenda)
  for( Index f_index=0; f_index<nf; f_index++ )
    {
      if (failed) continue;

      try {
        ArrayOfArrayOfSingleScatteringData   scat_data_mono;

        scat_data_monoCalc( scat_data_mono, scat_data,
                            f_grid, f_index, verbosity );

        // Seed reset for each loop. If not done, the errors
        // appear to be highly correlated.
        Index    mc_seed;
        MCSetSeedFromTime( mc_seed, verbosity );

        Vector       y, mc_error;
        Index        mc_iteration_count;
        Tensor3      mc_points;
        ArrayOfIndex mc_scat_order, mc_source_domain;

        MCGeneral( l_ws, y, mc_iteration_count, mc_error, mc_points, 
                   mc_scat_order, mc_source_domain, mc_antenna,
                   f_grid, f_index, pos, los, stokes_dim, atmosphere_dim,
                   l_ppath_step_agenda, ppath_lraytrace, l_iy_space_agenda, 
                   l_surface_rtprop_agenda, l_propmat_clearsky_agenda, 
                   p_grid, lat_grid, lon_grid, z_field, 
                   refellipsoid, z_surface, t_field, vmr_field,
                   cloudbox_on, cloudbox_limits,
                   pnd_field, scat_data_mono, 1, 1, 1, iy_unit,
                   mc_seed, mc_std_err, mc_max_time, mc_max_iter,
                   mc_min_iter, 1, verbosity);
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
void iyReplaceFromAux(
         Matrix&            iy,
   const ArrayOfTensor4&    iy_aux,
   const ArrayOfString&     iy_aux_vars,
   const Index&             jacobian_do,
   const String&            aux_var,
   const Verbosity& )
{
  if( iy_aux.nelem() != iy_aux_vars.nelem() )
    throw runtime_error( "*iy_aux* and *iy_aux_vars* must have the same "
                         "number of elements." );
  
  if( jacobian_do )
    throw runtime_error( "This method can not provide any jacobians and "
                         "*jacobian_do* must be 0." );

  bool ready = false;

  for( Index i=0; i<iy_aux.nelem() && !ready; i++ )
    {
      if( iy_aux_vars[i] == aux_var )
        {
          if( iy_aux[i].nrows() > 1  ||  iy_aux[i].ncols() > 1 )
            {
              throw runtime_error( "If an auxiliary variable shall be inserted "
                     "in *iy*, its row and page dimensions must have size 1." );
            }
          if( iy_aux[i].nbooks() != iy.nrows() )
            {
              throw runtime_error( "If an auxiliary variable shall be inserted "
                                   "in *iy*, its frequency dimension must match"
                                   "the length of existing *iy*." );
            }

          iy = 0;
          
          for( Index iv=0; iv<iy.nrows(); iv++ ) 
            {
              for( Index is=0; is<iy_aux[i].npages(); is++ )
                {
                  iy(iv,is) = iy_aux[i](iv,is,0,0);
                }
            }
          
          ready = true;
        }
    }

  if( !ready )
    throw runtime_error( "The selected auxiliary variable to insert in *iy* "
                         "is either not defined at all or is not set." );

}





/* Workspace method: Doxygen documentation will be auto-generated */
void iy_auxFillParticleVariables(
         ArrayOfTensor4&       iy_aux,
   const Index&                atmfields_checked,
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
  if( atmfields_checked != 1 )
    throw runtime_error( "The atmospheric fields must be flagged to have "
                         "passed a consistency check (atmfields_checked=1)." );
  if( cloudbox_checked != 1 )
    throw runtime_error( "The cloudbox must be flagged to have "
                         "passed a consistency check (cloudbox_checked=1)." );
  if( !cloudbox_on )
    throw runtime_error( 
                    "The cloudbox must be activated (cloudbox_on must be 1)" );
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
          // Particle number density
          for( Index j=0; j<auxPartField.nelem(); j++ )
            { iy_aux[auxPartField[j]](0,0,0,ip) = 
                                              ppath_pnd(auxPartFieldI[j],ip); }
        }  
    }
}


void yCalc_mblock_loop_body(
         bool&                       failed,
         String&                     fail_msg,
         ArrayOfArrayOfVector&       iyb_aux_array,
         Workspace&                  ws,
         Vector&                     y,
         Vector&                     y_f,
         ArrayOfIndex&               y_pol,
         Matrix&                     y_pos,
         Matrix&                     y_los,
         Matrix&                     y_geo,
         Matrix&                     jacobian,
   const Index&                      atmosphere_dim,
   const Tensor3&                    t_field,
   const Tensor3&                    z_field,
   const Tensor4&                    vmr_field,
   const Index&                      cloudbox_on,
   const Index&                      stokes_dim,
   const Vector&                     f_grid,
   const Matrix&                     sensor_pos,
   const Matrix&                     sensor_los,
   const Matrix&                     transmitter_pos,
   const Matrix&                     mblock_dlos_grid,
   const Sparse&                     sensor_response,
   const Vector&                     sensor_response_f,
   const ArrayOfIndex&               sensor_response_pol,
   const Matrix&                     sensor_response_dlos,
   const String&                     iy_unit,   
   const Agenda&                     iy_main_agenda,
   const Agenda&                     geo_pos_agenda,
   const Agenda&                     jacobian_agenda,
   const Index&                      jacobian_do,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const ArrayOfString&              iy_aux_vars,
   const Verbosity&                  verbosity,
   const Index&                      mblock_index,
   const Index&                      n1y,
   const Index&                      j_analytical_do)
{
    try
    {
        // Calculate monochromatic pencil beam data for 1 measurement block
        //
        Vector          iyb, iyb_error, yb(n1y);
        ArrayOfMatrix   diyb_dx;
        Matrix          geo_pos_matrix;
        //
        iyb_calc(ws, iyb, iyb_aux_array[mblock_index], diyb_dx, geo_pos_matrix,
                 mblock_index, atmosphere_dim, t_field, z_field, vmr_field,
                 cloudbox_on, stokes_dim, f_grid, sensor_pos, sensor_los,
                 transmitter_pos, mblock_dlos_grid, 
                 iy_unit, iy_main_agenda, geo_pos_agenda,
                 j_analytical_do, jacobian_quantities,
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

        // Fill information variables. And search for NaNs in *y*.
        //
        for( Index i=0; i<n1y; i++ )
          {
            const Index ii = row0 + i; 
            if( isnan(y[ii] ) )
              throw runtime_error( "One or several NaNs found in *y*." );
            y_f[ii]          = sensor_response_f[i];
            y_pol[ii]        = sensor_response_pol[i];
            y_pos(ii,joker)  = sensor_pos(mblock_index,joker);
            y_los(ii,joker)  = sensor_los(mblock_index,joker);
            y_los(ii,0)     += sensor_response_dlos(i,0);
            if( sensor_response_dlos.ncols() > 1 )
              { y_los(ii,1) += sensor_response_dlos(i,1); }
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
            jacobian_agendaExecute( ws, jacobian, mblock_index, iyb, yb,
                                    jacobian_agenda );
        }


        // Handle geo-positioning
        if( geo_pos_matrix(0,0) > -99900 )  // No data are flagged with -99999
          {
            // Find bore sigtht direction be proping sensor_response
            const Index   nf   = f_grid.nelem();
            const Index   nlos = mblock_dlos_grid.nrows();
            const Index   niyb = nf * nlos * stokes_dim;
            ArrayOfIndex i_of_max( n1y );
            Vector max_contr( n1y, -99999 );
            for( Index ilos=0; ilos<nlos; ilos++ )
              {
                Vector itry( niyb, 0 );
                itry[Range(ilos*nf*stokes_dim,nf,stokes_dim)] = 1;
                Vector ytry( n1y );
                mult( ytry, sensor_response, itry );
                for( Index i=0; i<n1y; i++ )
                  {
                    if( ytry[i] > max_contr[i] )
                      {
                        max_contr[i] = ytry[i];
                        i_of_max[i]  = ilos;
                      }
                  }
              }

            // Extract geo_pos_matrix for found bore-sights
            for( Index i=0; i<n1y; i++ )
              { y_geo(row0+i,joker) = geo_pos_matrix(i_of_max[i],joker); }
          }
    }

    catch (runtime_error e)
    {
#pragma omp critical (yCalc_fail)
        { fail_msg = e.what(); failed = true; }
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
         Matrix&                     y_geo,
         Matrix&                     jacobian,
   const Index&                      atmfields_checked,
   const Index&                      atmgeom_checked,
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
   const Matrix&                     mblock_dlos_grid,
   const Sparse&                     sensor_response,
   const Vector&                     sensor_response_f,
   const ArrayOfIndex&               sensor_response_pol,
   const Matrix&                     sensor_response_dlos,
   const String&                     iy_unit,   
   const Agenda&                     iy_main_agenda,
   const Agenda&                     geo_pos_agenda,
   const Agenda&                     jacobian_agenda,
   const Index&                      jacobian_do,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfArrayOfIndex&        jacobian_indices,
   const ArrayOfString&              iy_aux_vars,
   const Verbosity&                  verbosity )
{
  CREATE_OUT3;

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


  // Some sizes
  const Index   nf      = f_grid.nelem();
  const Index   nlos    = mblock_dlos_grid.nrows();
  const Index   n1y     = sensor_response.nrows();
  const Index   nmblock = sensor_pos.nrows();
  const Index   niyb    = nf * nlos * stokes_dim;


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
  y_geo.resize( nmblock*n1y, atmosphere_dim );
  y_geo = -99999;   // Will be replaced if relavant data are provided (*geo_pos*)

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

  String fail_msg;
  bool failed = false;

  if (nmblock >= arts_omp_get_max_threads()
      || (nf <= nmblock && nmblock >= nlos))
  {
      out3 << "  Parallelizing mblock loop (" << nmblock << " iterations)\n";

      // We have to make a local copy of the Workspace and the agendas because
      // only non-reference types can be declared firstprivate in OpenMP
      Workspace l_ws (ws);
      Agenda l_jacobian_agenda (jacobian_agenda);
      Agenda l_iy_main_agenda (iy_main_agenda);
      Agenda l_geo_pos_agenda (geo_pos_agenda);

#pragma omp parallel for                         \
  firstprivate(l_ws, l_jacobian_agenda, l_iy_main_agenda, l_geo_pos_agenda)
      for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
      {
          // Skip remaining iterations if an error occurred
          if (failed) continue;

          yCalc_mblock_loop_body( failed, fail_msg, iyb_aux_array, l_ws,
                                  y, y_f, y_pol, y_pos, y_los, y_geo, jacobian,
                                  atmosphere_dim, t_field, z_field,
                                  vmr_field,
                                  cloudbox_on, stokes_dim, f_grid, 
                                  sensor_pos, sensor_los, transmitter_pos,
                                  mblock_dlos_grid, sensor_response,
                                  sensor_response_f, sensor_response_pol,
                                  sensor_response_dlos, iy_unit,
                                  l_iy_main_agenda, l_geo_pos_agenda,
                                  l_jacobian_agenda,
                                  jacobian_do, jacobian_quantities,
                                  jacobian_indices, iy_aux_vars, verbosity,
                                  mblock_index, n1y, j_analytical_do );          
      }  // End mblock loop
  }
  else
  {
     out3 << "  Not parallelizing mblock loop (" << nmblock << " iterations)\n";

     for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
      {
          // Skip remaining iterations if an error occurred
          if (failed) continue;

          yCalc_mblock_loop_body( failed, fail_msg, iyb_aux_array, ws,
                                  y, y_f, y_pol, y_pos, y_los, y_geo, jacobian,
                                  atmosphere_dim, t_field, z_field,
                                  vmr_field,
                                  cloudbox_on, stokes_dim, f_grid, 
                                  sensor_pos, sensor_los, transmitter_pos,
                                  mblock_dlos_grid, sensor_response,
                                  sensor_response_f, sensor_response_pol,
                                  sensor_response_dlos, iy_unit,
                                  iy_main_agenda, geo_pos_agenda,
                                  jacobian_agenda,
                                  jacobian_do, jacobian_quantities,
                                  jacobian_indices, iy_aux_vars, verbosity,
                                  mblock_index, n1y, j_analytical_do ); 
      }  // End mblock loop
  }

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
void yCalcAppend(
         Workspace&                  ws,
         Vector&                     y,
         Vector&                     y_f,
         ArrayOfIndex&               y_pol,
         Matrix&                     y_pos,
         Matrix&                     y_los,
         ArrayOfVector&              y_aux,
         Matrix&                     y_geo,
         Matrix&                     jacobian,
         ArrayOfRetrievalQuantity&   jacobian_quantities,
         ArrayOfArrayOfIndex&        jacobian_indices,
   const Index&                      atmfields_checked,
   const Index&                      atmgeom_checked,
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
   const Matrix&                     mblock_dlos_grid,
   const Sparse&                     sensor_response,
   const Vector&                     sensor_response_f,
   const ArrayOfIndex&               sensor_response_pol,
   const Matrix&                     sensor_response_dlos,
   const String&                     iy_unit,   
   const Agenda&                     iy_main_agenda,
   const Agenda&                     geo_pos_agenda,
   const Agenda&                     jacobian_agenda,
   const Index&                      jacobian_do,
   const ArrayOfString&              iy_aux_vars,
   const ArrayOfRetrievalQuantity&   jacobian_quantities1,
   const ArrayOfArrayOfIndex&        jacobian_indices1,
   const Index&                      append_instrument_wfs,
   const Verbosity&                  verbosity )
{
  // Check consistency of data representing first measurement
  const Index n1   = y.nelem();
        Index nrq1 = 0; 
  if( y.empty() )
    throw runtime_error( "Input *y* is empty. Use *yCalc*" );
  if( y_f.nelem() != n1 )
    throw runtime_error( "Lengths of input *y* and *y_f* are inconsistent." );
  if( y_pol.nelem() != n1 )
    throw runtime_error( "Lengths of input *y* and *y_pol* are inconsistent." );
  if( y_pos.nrows() != n1 )
    throw runtime_error( "Sizes of input *y* and *y_pos* are inconsistent." );
  if( y_los.nrows() != n1 )
    throw runtime_error( "Sizes of input *y* and *y_los* are inconsistent." );
  if( y_geo.nrows() != n1 )
    throw runtime_error( "Sizes of input *y* and *y_geo* are inconsistent." );
  if( jacobian_do )
    {
      nrq1 = jacobian_quantities1.nelem();
      if( jacobian.nrows() != n1 )
        throw runtime_error( "Sizes of *y* and *jacobian* are inconsistent." );
      if( jacobian_indices1.nelem() != nrq1 )
        throw runtime_error( "Lengths of *jacobian_quantities_copy* and "
                             "*jacobian_indices_copy* are inconsistent." );
      if( jacobian.ncols() != jacobian_indices1[nrq1-1][1]+1 )
        throw runtime_error( "Size of input *jacobian* and max value in " 
                             "*jacobian_indices_copy* are inconsistent." );
    }

  // Calculate new measurement
  //
  Vector        y2, y_f2;
  Matrix        y_pos2, y_los2, y_geo2, jacobian2;
  ArrayOfIndex  y_pol2;
  ArrayOfVector y_aux2;
  //
  yCalc( ws, y2, y_f2, y_pol2, y_pos2, y_los2, y_aux2, y_geo2, jacobian2,
         atmfields_checked, atmgeom_checked, atmosphere_dim, t_field,
         z_field, vmr_field, cloudbox_on, cloudbox_checked, sensor_checked,
         stokes_dim, f_grid, sensor_pos, sensor_los, transmitter_pos,
         mblock_dlos_grid, sensor_response,
         sensor_response_f, sensor_response_pol, sensor_response_dlos, 
         iy_unit, iy_main_agenda, geo_pos_agenda, 
         jacobian_agenda, jacobian_do, jacobian_quantities, jacobian_indices, 
         iy_aux_vars, verbosity );

  // Consistency checks
  if( y_pos.ncols() != y_pos2.ncols() )
    throw runtime_error( 
          "Different number of columns in *y_pos* between the measurements." );
  if( y_los.ncols() != y_los2.ncols() )
    throw runtime_error( 
          "Different number of columns in *y_los* between the measurements." );
  

  // y and y_XXX
  //
  const Index         n2 = y2.nelem(); 
  //
  {
    // Make copy of old measurement
    const Vector        y1=y, y_f1=y_f;
    const Matrix        y_pos1=y_pos, y_los1=y_los, y_geo1=y_geo;
    const ArrayOfIndex  y_pol1=y_pol;
    const ArrayOfVector y_aux1=y_aux;
    //
    y.resize( n1+n2 );
    y[Range(0,n1)] = y1;   y[Range(n1,n2)] = y2; 
    //
    y_f.resize( n1+n2 );
    y_f[Range(0,n1)] = y_f1;   y_f[Range(n1,n2)] = y_f2; 
    //
    y_pos.resize( n1+n2, y_pos1.ncols() );
    y_pos(Range(0,n1),joker) = y_pos1;   y_pos(Range(n1,n2),joker) = y_pos2; 
    //
    y_los.resize( n1+n2, y_los1.ncols() );
    y_los(Range(0,n1),joker) = y_los1;   y_los(Range(n1,n2),joker) = y_los2; 
    //
    y_geo.resize( n1+n2, y_geo1.ncols() );
    y_geo(Range(0,n1),joker) = y_geo1;   y_geo(Range(n1,n2),joker) = y_geo2; 
    //
    y_pol.resize( n1+n2 );
    for( Index i=0; i<n1; i++ )
      { y_pol[i] = y_pol1[i]; }
    for( Index i=0; i<n2; i++ )
      { y_pol[n1+i] = y_pol2[i]; }

    // y_aux
    const Index na1 = y_aux1.nelem();
    const Index na2 = y_aux2.nelem();
    const Index na  = max(na1,na2);
    //
    y_aux.resize( na );
    //
    for( Index a=0; a<na; a++ )
      {
        y_aux[a].resize( n1+n2 );        
        if( a < na1 )
          { y_aux[a][Range(0,n1)] = y_aux1[a]; }
        else
          { y_aux[a][Range(0,n1)] = 0; }
        if( a < na2 )
          { y_aux[a][Range(n1,n2)] = y_aux2[a]; }
        else
          { y_aux[a][Range(n1,n2)] = 0; }
      }
  }

  // Jacobian and friends
  if( jacobian_do )
    {
      // Put in old jacobian_quantities and jacobian_indices as first part in
      // new version of these variables
      ArrayOfRetrievalQuantity  jacobian_quantities2 = jacobian_quantities;
      ArrayOfArrayOfIndex       jacobian_indices2    = jacobian_indices;
      //
      jacobian_quantities = jacobian_quantities1;
      jacobian_indices    = jacobian_indices1;

      // Loop new jacobian_quantities to determine how new jacobian data shall
      // be inserted
      //
      const Index    nrq2 = jacobian_quantities2.nelem();
      ArrayOfIndex   map_table(nrq2);
      //
      for( Index q2=0; q2<nrq2; q2++ )
        {

          Index pos = -1;

          // Compare to old quantities, to determine if append shall be 
          // considered. Some special checks performed here, grids checked later
          if( jacobian_quantities2[q2].MainTag() == ABSSPECIES_MAINTAG   ||
              jacobian_quantities2[q2].MainTag() == TEMPERATURE_MAINTAG  ||
              jacobian_quantities2[q2].MainTag() == WIND_MAINTAG         ||
              append_instrument_wfs )
            {
              for( Index q1=0; q1<nrq1; q1++ && pos < 0 )
                {
                  if( jacobian_quantities2[q2].MainTag() ==
                      jacobian_quantities1[q1].MainTag() )
                    {
                      // Absorption species
                      if( jacobian_quantities2[q2].MainTag() == 
                                                           ABSSPECIES_MAINTAG )
                        {
                          if( jacobian_quantities2[q2].Subtag() ==
                              jacobian_quantities1[q1].Subtag() )
                            {
                              if( jacobian_quantities2[q2].Mode() ==
                                  jacobian_quantities1[q1].Mode() )
                                { pos = q1; }
                              else
                                {
                                  ostringstream os;
                                  os << "Jacobians for " 
                                     << jacobian_quantities2[q2].MainTag()<<"/" 
                                     << jacobian_quantities2[q2].Subtag() 
                                     << " shall be appended.\nThis requires "
                                     << "that the same retrieval unit is used "
                                     << "but it seems that this requirement is "
                                     << "not met.";
                                  throw runtime_error(os.str());
                                }
                            }
                        }
                      // Temperature
                      else if( jacobian_quantities2[q2].MainTag() == 
                                                          TEMPERATURE_MAINTAG )
                        {
                          if( jacobian_quantities2[q2].Subtag() ==
                              jacobian_quantities1[q1].Subtag() )
                            { pos = q1; }
                          else
                            {
                              ostringstream os;
                              os << "Jacobians for " 
                                 << jacobian_quantities2[q2].MainTag()<<"/" 
                                 << jacobian_quantities2[q2].Subtag() 
                                 << " shall be appended.\nThis requires "
                                 << "that HSE is either ON or OFF for both "
                                 << "parts but it seems that this requirement "
                                 << "is not met.";
                              throw runtime_error(os.str());
                            }
                        }
                      // Other
                      else if( jacobian_quantities2[q2].Subtag() ==
                               jacobian_quantities1[q1].Subtag() )
                        { pos = q1; }
                    }
                }
            }


          // New quantity
          if( pos < 0 )
            {
              map_table[q2] = jacobian_quantities.nelem();
              jacobian_quantities.push_back( jacobian_quantities2[q2] );
              ArrayOfIndex indices(2);
              indices[0] = jacobian_indices[jacobian_indices.nelem()-1][1]+1;
              indices[1] = indices[0] + jacobian_indices2[q2][1] -
                                        jacobian_indices2[q2][0];
              jacobian_indices.push_back( indices );
            }
          // Existing quantity
          else
            {
              map_table[q2] = pos;
              // Check if grids are equal
              ArrayOfVector grids1 = jacobian_quantities1[pos].Grids();
              ArrayOfVector grids2 = jacobian_quantities2[q2].Grids();
              bool any_wrong = false;
              if( grids1.nelem() != grids2.nelem() )
                { any_wrong = true; }
              else
                {
                  for( Index g=0; g<grids1.nelem(); g++ )
                    {
                      if( grids1[g].nelem() != grids2[g].nelem() )
                        { any_wrong = true; }
                      else
                        {
                          for( Index e=0; e<grids1[g].nelem(); e++ )
                            {
                              const Numeric v1 = grids1[g][e];
                              const Numeric v2 = grids2[g][e];
                              if( ( v1 == 0  &&  abs(v2) > 1e-9 ) || 
                                                 abs(v1-v2)/v1 > 1e-6 )
                                { any_wrong = true; }
                            }   
                        }   
                    }   
                }
              if( any_wrong )
                {
                  ostringstream os;
                  os << "Jacobians for " 
                     << jacobian_quantities2[q2].MainTag()
                     << "/" << jacobian_quantities2[q2].Subtag() 
                     << " shall be appended.\nThis requires that the "
                     << "same grids are used for both measurements,\nbut "
                     << "it seems that this requirement is not met.";
                  throw runtime_error(os.str());
                }
            }
        }

      // Create and fill *jacobian*
      //
      const Index  nrq = jacobian_quantities.nelem();
      const Matrix jacobian1 = jacobian;
      //
      jacobian.resize( n1+n2, jacobian_indices[nrq-1][1]+1 );
      jacobian = 0;
      //
      // Put in old part in top-left corner
      jacobian(Range(0,n1),Range(0,jacobian_indices1[nrq1-1][1]+1)) = jacobian1;
      // New parts
      for( Index q2=0; q2<nrq2; q2++ )
        {
          jacobian(Range(n1,n2),Range(jacobian_indices[map_table[q2]][0], 
                                      jacobian_indices[map_table[q2]][1]-
                                      jacobian_indices[map_table[q2]][0]+1) ) =
                jacobian2(joker,Range(jacobian_indices2[q2][0], 
                                      jacobian_indices2[q2][1]-
                                      jacobian_indices2[q2][0]+1) );
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
