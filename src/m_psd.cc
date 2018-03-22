/* Copyright (C) 2017
   Patrick Eriksson <patrick.eriksson@chalmers.se>
   Jana Mendrok     <jana.mendrok@gmail.com>
   Manfred Brath    <manfred.brath@uni-hamburg.de>
                         
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
  \file   m_psd.cc
  \author Patrick Eriksson, Jana Mendrok, Manfred Brath
  \date   2017-11-05

  \brief  Workspace functions related to particle size distributions.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h. */



/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <cmath>
#include <stdexcept>
#include <cstdlib>

#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "lin_alg.h"
#include "math_funcs.h"
#include "messages.h"
#include "psd.h"
#include "physics_funcs.h"




// ------------------------------------------------------
// Macros to avoid duplication of code
// ------------------------------------------------------

#define START_OF_PSD_METHODS() \
  const Index nin = pnd_agenda_input_names.nelem(); \
  const Index ndx = dpnd_data_dx_names.nelem(); \
  const Index np  = pnd_agenda_input.nrows(); \
  const Index nsi = psd_size_grid.nelem(); \
  ArrayOfIndex dx2in(ndx); \
   \
  if( pnd_agenda_input.ncols() != nin ) \
    throw runtime_error( "Length of *pnd_agenda_input_names* and number of " \
                         "columns in *pnd_agenda_input* must be equal." ); \
  if( ndx ) \
    { \
      if( ndx > nin ) \
        throw runtime_error( "The length of *dpnd_data_dx_names* can not " \
                             "exceed the one of *pnd_agenda_input_names*." ); \
      for( Index i=0; i<ndx; i++ ) \
        { \
          dx2in[i] = find_first( pnd_agenda_input_names, dpnd_data_dx_names[i] ); \
          if( dx2in[i] < 0 ) \
            { \
              ostringstream os; \
              os << "dpnd_data_dx_names[" << i << "] is " << dpnd_data_dx_names[i] \
                 << "\nThis string could not be found in *pnd_agenda_input_names*.";\
              throw std::runtime_error(os.str()); \
            } \
        } \
    } \
   \
  psd_data.resize( np, nsi ); \
  psd_data = 0.0; \
  if( ndx ) \
    { \
      dpsd_data_dx.resize( ndx, np, nsi ); \
      dpsd_data_dx = 0.0; \
    } \
  else \
    { dpsd_data_dx.resize( 0, 0, 0  ); }  

// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------





/*===========================================================================
  === PSDs of Mono type
  ===========================================================================*/

/* Common code for mono-type PSDs */
void psd_mono_common (
          Matrix&           psd_data,
          Tensor3&          dpsd_data_dx,
    const String&           type,
    const Vector&           pnd_agenda_input_t,
    const Matrix&           pnd_agenda_input,
    const ArrayOfString&    pnd_agenda_input_names,
    const ArrayOfString&    dpnd_data_dx_names,
    const ArrayOfArrayOfScatteringMetaData&   scat_meta,
    const Index&            species_index,          
    const Numeric&          t_min, 
    const Numeric&          t_max, 
    const Index&            picky, 
    const Verbosity&)
{
  // Standard checcks
  const Vector psd_size_grid(1,0);    // As this WSV is not input for thse WSM
  START_OF_PSD_METHODS();

  // Extra checks for this PSD
  const Index nss = scat_meta.nelem();
  if( nss == 0 )
    throw runtime_error( "*scat_meta* is empty!" );
  if( nss < species_index+1 )
    {
      ostringstream os;
      os << "Selected scattering species index is " << species_index << " but this "
         << "is not allowed since *scat_meta* has only " << nss << " elements.";
      throw runtime_error(os.str());
    }
  if( scat_meta[species_index].nelem() != 1 )
    {
      ostringstream os;
      os << "This method only works with scattering species consisting of a\n" 
         << "single element, but your data do not match this demand.\n"
         << "Selected scattering species index is " << species_index << ".\n"
         << "This species has " << scat_meta[species_index].nelem() << " elements.";
      throw runtime_error(os.str());
    }
  //
  if( pnd_agenda_input.ncols() != 1 ) 
    throw runtime_error( "*pnd_agenda_input* must have one column." ); 
  if( nsi != 1 )
    throw runtime_error( "This method demands that length of "
                         "*psd_size_grid* is 1." );                     
  
  // Extract particle mass
  Numeric pmass = 0;
  if( type == "mass" )
    { pmass = scat_meta[species_index][0].mass; }

  for( Index ip=0; ip<np; ip++ )
    {
      // Extract the input variables
      Numeric   x = pnd_agenda_input(ip,0);
      Numeric   t = pnd_agenda_input_t[ip];

      // No calc needed if n==0 and no jacobians requested.
      if( (x==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: [ max(0.,"
                 << t_min << "), " << t_max << " ]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }
      
      // Set PSD
      //
      if( type == "ntot" )
        {
          psd_data(ip,0) = x;
          //
          if( ndx )
            { dpsd_data_dx(0,ip,0) = 1; }
        }
      else if( type == "mass" )
        {
          psd_data(ip,0) = x/pmass;
          //
          if( ndx )
            { dpsd_data_dx(0,ip,0) = 1/pmass; }   
        }
      else
        { assert(0); }
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void psdMono (
          Matrix&           psd_data,
          Tensor3&          dpsd_data_dx,
    const Vector&           pnd_agenda_input_t,
    const Matrix&           pnd_agenda_input,
    const ArrayOfString&    pnd_agenda_input_names,
    const ArrayOfString&    dpnd_data_dx_names,
    const ArrayOfArrayOfScatteringMetaData&   scat_meta,
    const Index&            species_index,          
    const Numeric&          t_min, 
    const Numeric&          t_max, 
    const Index&            picky, 
    const Verbosity&        verbosity )
{
  psd_mono_common( psd_data, dpsd_data_dx,
                   "ntot",
                   pnd_agenda_input_t, pnd_agenda_input,
                   pnd_agenda_input_names, dpnd_data_dx_names,
                   scat_meta,  species_index,
                   t_min, t_max,  picky, verbosity );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void psdMonoMass (
          Matrix&           psd_data,
          Tensor3&          dpsd_data_dx,
    const Vector&           pnd_agenda_input_t,
    const Matrix&           pnd_agenda_input,
    const ArrayOfString&    pnd_agenda_input_names,
    const ArrayOfString&    dpnd_data_dx_names,
    const ArrayOfArrayOfScatteringMetaData&   scat_meta,
    const Index&            species_index,          
    const Numeric&          t_min, 
    const Numeric&          t_max, 
    const Index&            picky, 
    const Verbosity&        verbosity )
{
  psd_mono_common( psd_data, dpsd_data_dx,
                   "mass",
                   pnd_agenda_input_t, pnd_agenda_input,
                   pnd_agenda_input_names, dpnd_data_dx_names,
                   scat_meta,  species_index,
                   t_min, t_max,  picky, verbosity );
}





/*===========================================================================
  === PSDs of MGD type
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void psdMgd(
          Matrix&          psd_data,
          Tensor3&         dpsd_data_dx,
    const Vector&          psd_size_grid,
    const Vector&          pnd_agenda_input_t,
    const Matrix&          pnd_agenda_input,
    const ArrayOfString&   pnd_agenda_input_names,
    const ArrayOfString&   dpnd_data_dx_names,
    const Numeric&         n0, 
    const Numeric&         mu, 
    const Numeric&         la, 
    const Numeric&         ga, 
    const Numeric&         t_min, 
    const Numeric&         t_max, 
    const Index&           picky, 
    const Verbosity& )
{
  // Standard checks
  START_OF_PSD_METHODS();
  
  // Additional (basic) checks
  if( nin > 4 )
    throw runtime_error( "The number of columns in *pnd_agenda_input* must "
                         "be 0, 1, 2, 3 or 4." );
  
  // Check fixed parameters
  const Index n0_fixed = (Index) !( isnan(n0) );
  const Index mu_fixed = (Index) !( isnan(mu) );
  const Index la_fixed = (Index) !( isnan(la) );
  const Index ga_fixed = (Index) !( isnan(ga) );
  //
  if( nin + n0_fixed + mu_fixed + la_fixed + ga_fixed != 4 )
    throw runtime_error( "This PSD has four free parameters. This means that "
                         "the number\nof columns in *pnd_agenda_input* and the "
                         "number of numerics\n(i.e. non-NaN) and among "
                         "the GIN arguments n0, mu, la and\nga must add up to "
                         "four. And this was found not to be the case." );

  // Create vectors to hold the four MGD and the "extra" parameters 
  Vector mgd_pars(4);
  ArrayOfIndex mgd_i_pai = {-1,-1,-1,-1}; // Position in pnd_agenda_input
  {
    Index nhit = 0;
    if( n0_fixed ) { mgd_pars[0]=n0; } else { mgd_i_pai[0]=nhit++; } 
    if( mu_fixed ) { mgd_pars[1]=mu; } else { mgd_i_pai[1]=nhit++; } 
    if( la_fixed ) { mgd_pars[2]=la; } else { mgd_i_pai[2]=nhit++; } 
    if( ga_fixed ) { mgd_pars[3]=ga; } else { mgd_i_pai[3]=nhit++; } 
  }

  // Determine what derivatives to do and their positions
  ArrayOfIndex mgd_do_jac = {0,0,0,0,}; 
  ArrayOfIndex mgd_i_jac = {-1,-1,-1,-1}; // Position among jacobian quantities
  //
  for( Index i=0; i<ndx; i++ )
    {
      for( Index j=0; j<4; j++ )
        {
          if( dx2in[i] == mgd_i_pai[j] )
            {
              mgd_do_jac[j] = 1;
              mgd_i_jac[j]  = i;
              break;
            }
        }
    }

  // Loop input data and calculate PSDs
  for( Index ip=0; ip<np; ip++ )
    {
      // Extract MGD parameters
      for( Index i=0; i<4; i++ )
        {
          if( mgd_i_pai[i] >= 0 )
            { mgd_pars[i] = pnd_agenda_input(ip,mgd_i_pai[i]); }
        }
      Numeric t = pnd_agenda_input_t[ip];
      
      // No calc needed if n0==0 and no jacobians requested.
      if( (mgd_pars[0]==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: [ max(0.,"
                 << t_min << "), " << t_max << " ]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }

      // Check that la and ga are OK
      if( mgd_pars[2] <= 0 )
        throw runtime_error( "Bad MGD parameter detected: la <= 0" );
      if( mgd_pars[3] <= 0 )
        throw runtime_error( "Bad MGD parameter detected: ga <= 0" );
      
      // Calculate PSD and derivatives
      Matrix  jac_data(4,nsi);
      //
      mgd_with_derivatives( psd_data(ip,joker), jac_data, psd_size_grid,
           mgd_pars[0], mgd_pars[1], mgd_pars[2], mgd_pars[3],
           mgd_do_jac[0], mgd_do_jac[1], mgd_do_jac[2], mgd_do_jac[3] );
      //
      for( Index i=0; i<4; i++ )
        {
          if( mgd_do_jac[i] )
            { dpsd_data_dx(mgd_i_jac[i],ip,joker) = jac_data(i,joker); }
        } 
    }  
}



/* Workspace method: Doxygen documentation will be auto-generated */
void psdMgdMass(
          Matrix&          psd_data,
          Tensor3&         dpsd_data_dx,
    const Vector&          psd_size_grid,
    const Vector&          pnd_agenda_input_t,
    const Matrix&          pnd_agenda_input,
    const ArrayOfString&   pnd_agenda_input_names,
    const ArrayOfString&   dpnd_data_dx_names,
    const Numeric&         scat_species_a, 
    const Numeric&         scat_species_b, 
    const Numeric&         n0, 
    const Numeric&         mu, 
    const Numeric&         la, 
    const Numeric&         ga, 
    const Numeric&         t_min, 
    const Numeric&         t_max, 
    const Index&           picky, 
    const Verbosity&)
{
  // Standard checks
  START_OF_PSD_METHODS();
  
  // Additional (basic) checks
  if( nin < 1 || nin > 4 )
    throw runtime_error( "The number of columns in *pnd_agenda_input* must "
                         "be 1, 2, 3 or 4." );
  if( scat_species_a <= 0 )
    throw runtime_error( "*scat_species_a* should be > 0." );
  if( scat_species_b <= 0  ||  scat_species_b >= 5 )
    throw runtime_error( "*scat_species_b* should be > 0 and < 5." );
  
  // Check and determine dependent and fixed parameters
  const Index n0_depend = (Index) n0 == -999;
  const Index mu_depend = (Index) mu == -999;
  const Index la_depend = (Index) la == -999;
  const Index ga_depend = (Index) ga == -999;  
  //
  if( n0_depend + mu_depend + la_depend + ga_depend != 1 )
    throw runtime_error( "One (but only one) of n0, mu, la and ga must be NaN, "
                         "to flag that this parameter is the one dependent of "
                         "mass content." );
  if( mu_depend  ||  ga_depend )
    throw runtime_error( "Sorry, mu and la are not yet allowed to be the "
                         "dependent parameter." );    
  //
  const Index n0_fixed = (Index) !( n0_depend  ||  isnan(n0) );
  const Index mu_fixed = (Index) !( mu_depend  ||  isnan(mu) );
  const Index la_fixed = (Index) !( la_depend  ||  isnan(la) );
  const Index ga_fixed = (Index) !( ga_depend  ||  isnan(ga) );
  //
  if( nin + n0_fixed + mu_fixed + la_fixed + ga_fixed != 4 )
    throw runtime_error( "This PSD has four free parameters. This means that "
                         "the number\nof columns in *pnd_agenda_input* and the "
                         "number of numerics\n(i.e. not -999 or NaN) and among "
                         "the GIN arguments n0, mu, la and\nga must add up to "
                         "four. And this was found not to be the case." );

  // Create vectors to hold the four MGD and the "extra" parameters 
  Vector mgd_pars(4), ext_pars(1);
  ArrayOfIndex mgd_i_pai = {-1,-1,-1,-1}; // Position in pnd_agenda_input
  const ArrayOfIndex ext_i_pai = {0};     // Position in pnd_agenda_input
  {
    Index nhit = 1;  // As mass always occupies first position
    if( n0_fixed ) { mgd_pars[0]=n0; } else if( !n0_depend ) { mgd_i_pai[0]=nhit++; } 
    if( mu_fixed ) { mgd_pars[1]=mu; } else if( !mu_depend ) { mgd_i_pai[1]=nhit++; } 
    if( la_fixed ) { mgd_pars[2]=la; } else if( !la_depend ) { mgd_i_pai[2]=nhit++; } 
    if( ga_fixed ) { mgd_pars[3]=ga; } else if( !ga_depend ) { mgd_i_pai[3]=nhit++; } 
  }

  // Determine what derivatives to do and their positions
  ArrayOfIndex mgd_do_jac = {0,0,0,0,}; 
  ArrayOfIndex ext_do_jac = {0};        
  ArrayOfIndex mgd_i_jac = {-1,-1,-1,-1}; // Position among jacobian quantities
  ArrayOfIndex ext_i_jac = {-1};          // Position among jacobian quantities
  //
  for( Index i=0; i<ndx; i++ )
    {
      if( dx2in[i] == 0 )  // That is,  mass is a derivative
        {
          ext_do_jac[0] = 1;
          ext_i_jac[0]  = i;
        }
      else  // Otherwise, either n0, mu, la or ga
        {
          for( Index j=0; j<4; j++ )
            {
              if( dx2in[i] == mgd_i_pai[j] )
                {
                  mgd_do_jac[j] = 1;
                  mgd_i_jac[j]  = i;
                  break;
                }
            }
        }
    }

  // Loop input data and calculate PSDs
  for( Index ip=0; ip<np; ip++ )
    {
      // Extract mass
      ext_pars[0] = pnd_agenda_input(ip,ext_i_pai[0]);
      // Extract core MGD parameters
      for( Index i=0; i<4; i++ )
        {
          if( mgd_i_pai[i] >= 0 )
            { mgd_pars[i] = pnd_agenda_input(ip,mgd_i_pai[i]); }
        }
      Numeric t = pnd_agenda_input_t[ip];
      
      // No calc needed if mass==0 and no jacobians requested.
      if( (ext_pars[0]==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: [ max(0.,"
                 << t_min << "), " << t_max << " ]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }

      
      // Derive the dependent parameter
      //
      Numeric mub1 = 0, eterm = 0, scfac = 0;
      //
      if( n0_depend )
        {
          mub1  = mgd_pars[1] + scat_species_b + 1;
          eterm = mub1 / mgd_pars[3];
          scfac = ( mgd_pars[3] * pow( mgd_pars[2], eterm ) ) /
            ( scat_species_a * tgamma(eterm) );
          mgd_pars[0] = scfac * ext_pars[0] ;
        }
      else if( la_depend )
        {
          if( ext_pars[0] <= 0 )
            throw runtime_error( "The mass content must be > 0 when la is "
                                 "the dependent parameter." );
          mub1  = mgd_pars[1] + scat_species_b + 1;
          eterm = mub1 / mgd_pars[3];
          scfac = mgd_pars[3] / ( scat_species_a * mgd_pars[0] * tgamma(eterm) );
          scfac = pow( scfac, -1/eterm );
          mgd_pars[2] = scfac * pow( ext_pars[0], -1/eterm );
        }
      else 
        { assert(0); }

      // Now when all four MGD parameters are set, check that they were OK from
      // start, or became OK if set
      if( mub1 <= 0 )
        throw runtime_error( "Bad MGD parameter detected: mu + b + 1 <= 0" );
      if( mgd_pars[2] <= 0 )
        throw runtime_error( "Bad MGD parameter detected: la <= 0" );
      if( mgd_pars[3] <= 0 )
        throw runtime_error( "Bad MGD parameter detected: ga <= 0" );
      
      // Calculate PSS
      Matrix  jac_data(4,nsi);    
      mgd_with_derivatives( psd_data(ip,joker), jac_data, psd_size_grid,
           mgd_pars[0], mgd_pars[1], mgd_pars[2], mgd_pars[3],
           (bool)mgd_do_jac[0] || n0_depend,
           (bool)mgd_do_jac[1] || mu_depend,
           (bool)mgd_do_jac[2] || la_depend,
           (bool)mgd_do_jac[3] || ga_depend );

      // Derivative with respect to mass
      if( ext_do_jac[0] )
        {
          if( n0_depend )
            {
              dpsd_data_dx(ext_i_jac[0],ip,joker) = jac_data(0,joker);
              dpsd_data_dx(ext_i_jac[0],ip,joker) *= scfac; 
            }
          else if( la_depend )
            {
              dpsd_data_dx(ext_i_jac[0],ip,joker) = jac_data(2,joker);
              dpsd_data_dx(ext_i_jac[0],ip,joker) *= scfac * (-1/eterm) *
                pow( ext_pars[0], -(1/eterm+1) ); 
            }
          else 
            { assert(0); }
        }

      // Derivatives for non-dependent native parameters
      for( Index i=0; i<4; i++ )
        {
          if( mgd_do_jac[i] )
            { dpsd_data_dx(mgd_i_jac[i],ip,joker) = jac_data(i,joker); }
        } 
    }  
}



/* Common code for MGD PSDs, taking mass and something as non-native input */
void psd_mgd_mass_and_something(
          Matrix&          psd_data,
          Tensor3&         dpsd_data_dx,
    const String&          something,
    const Vector&          psd_size_grid,
    const Vector&          pnd_agenda_input_t,
    const Matrix&          pnd_agenda_input,
    const ArrayOfString&   pnd_agenda_input_names,
    const ArrayOfString&   dpnd_data_dx_names,
    const Numeric&         scat_species_a, 
    const Numeric&         scat_species_b, 
    const Numeric&         n0, 
    const Numeric&         mu, 
    const Numeric&         la, 
    const Numeric&         ga, 
    const Numeric&         t_min, 
    const Numeric&         t_max, 
    const Index&           picky, 
    const Verbosity&)
{
  // Standard checks
  START_OF_PSD_METHODS();
  
  // Additional (basic) checks
  if( nin < 1 || nin > 4 )
    throw runtime_error( "The number of columns in *pnd_agenda_input* must "
                         "be 2, 3 or 4." );
  if( scat_species_a <= 0 )
    throw runtime_error( "*scat_species_a* should be > 0." );
  if( scat_species_b <= 0  ||  scat_species_b >= 5 )
    throw runtime_error( "*scat_species_b* should be > 0 and < 5." );
  
  // Check and determine dependent and fixed parameters
  const Index n0_depend = (Index) n0 == -999;
  const Index mu_depend = (Index) mu == -999;
  const Index la_depend = (Index) la == -999;
  const Index ga_depend = (Index) ga == -999;  
  //
  if( n0_depend + mu_depend + la_depend + ga_depend != 2 )
    throw runtime_error( "Two (but only two) of n0, mu, la and ga must be NaN, "
                         "to flag that these parameters are the ones dependent of "
                         "mass content and mean particle size." );
  if( mu_depend  ||  ga_depend )
    throw runtime_error( "Sorry, mu and la are not yet allowed to be a "
                         "dependent parameter." );    
  //
  const Index n0_fixed = (Index) !( n0_depend  ||  isnan(n0) );
  const Index mu_fixed = (Index) !( mu_depend  ||  isnan(mu) );
  const Index la_fixed = (Index) !( la_depend  ||  isnan(la) );
  const Index ga_fixed = (Index) !( ga_depend  ||  isnan(ga) );
  //
  if( nin + n0_fixed + mu_fixed + la_fixed + ga_fixed != 4 )
    throw runtime_error( "This PSD has four free parameters. This means that "
                         "the number\nof columns in *pnd_agenda_input* and the "
                         "number of numerics\n(i.e. not -999 or NaN) and among "
                         "the GIN arguments n0, mu, la and\nga must add up to "
                         "four. And this was found not to be the case." );

  // Create vectors to hold the four MGD and the "extra" parameters 
  Vector mgd_pars(4), ext_pars(2);
  ArrayOfIndex mgd_i_pai = {-1,-1,-1,-1}; // Position in pnd_agenda_input
  const ArrayOfIndex ext_i_pai = {0,1};   // Position in pnd_agenda_input
  {
    Index nhit = 2;  // As mass and Dm always occupy first position
    if( n0_fixed ) { mgd_pars[0]=n0; } else if( !n0_depend ) { mgd_i_pai[0]=nhit++; } 
    if( mu_fixed ) { mgd_pars[1]=mu; } else if( !mu_depend ) { mgd_i_pai[1]=nhit++; } 
    if( la_fixed ) { mgd_pars[2]=la; } else if( !la_depend ) { mgd_i_pai[2]=nhit++; } 
    if( ga_fixed ) { mgd_pars[3]=ga; } else if( !ga_depend ) { mgd_i_pai[3]=nhit++; } 
  }

  // Determine what derivatives to do and their positions
  ArrayOfIndex mgd_do_jac = {0,0,0,0,}; 
  ArrayOfIndex ext_do_jac = {0,0};        
  ArrayOfIndex mgd_i_jac = {-1,-1,-1,-1}; // Position among jacobian quantities
  ArrayOfIndex ext_i_jac = {-1,-1};       // Position among jacobian quantities
  //
  for( Index i=0; i<ndx; i++ )
    {
      if( dx2in[i] == 0 )  // That is, mass is a derivative
        {
          ext_do_jac[0] = 1;
          ext_i_jac[0]  = i;
        }
      else if( dx2in[i] == 1 )  // That is, "something" is a derivative
        {
          ext_do_jac[1] = 1;
          ext_i_jac[1]  = i;
        }
      else  // Otherwise, either n0, mu, la or ga
        {
          for( Index j=0; j<4; j++ )
            {
              if( dx2in[i] == mgd_i_pai[j] )
                {
                  mgd_do_jac[j] = 1;
                  mgd_i_jac[j]  = i;
                  break;
                }
            }
        }
    }

  // Loop input data and calculate PSDs
  for( Index ip=0; ip<np; ip++ )
    {
      // Extract mass
      ext_pars[0] = pnd_agenda_input(ip,ext_i_pai[0]);
      ext_pars[1] = pnd_agenda_input(ip,ext_i_pai[1]);
      if( ext_pars[1] <= 0 )
        {
          ostringstream os;
          os << "Negative " << something << "found.\nThis is not allowed.";
          throw std::runtime_error(os.str());
        }
      // Extract core MGD parameters
      for( Index i=0; i<4; i++ )
        {
          if( mgd_i_pai[i] >= 0 )
            { mgd_pars[i] = pnd_agenda_input(ip,mgd_i_pai[i]); }
        }
      Numeric t = pnd_agenda_input_t[ip];
      
      // No calc needed if mass==0 and no jacobians requested.
      if( (ext_pars[0]==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: [ max(0.,"
                 << t_min << "), " << t_max << " ]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }

      // Derive the dependent parameters (see ATD)
      //
      Numeric mu1 = 0, mub1 = 0, eterm = 0, gterm = 0;
      Numeric scfac1 = 0, scfac2 = 0, gab = 0;
      //
      // *** Mean size ***
      if ( something == "mean size" )
        {
          if( n0_depend  &&  la_depend )
            {
              mub1   = mgd_pars[1] + scat_species_b + 1;
              if( mub1 <= 0 )
                throw runtime_error( "Bad MGD parameter detected: mu + b + 1 <= 0" );
              eterm  = mub1 / mgd_pars[3];          
              // Start by deriving la 
              scfac2 = pow( eterm, mgd_pars[3] ); 
              mgd_pars[2] = scfac2 * pow( ext_pars[1], -mgd_pars[3] );
              // We can now derive n0 
              gterm = tgamma( eterm );
              scfac1 = ( mgd_pars[3] * pow( mgd_pars[2], eterm ) ) /
                ( scat_species_a * gterm );
              mgd_pars[0] = scfac1 * ext_pars[0];
            }
          else 
            { assert(0); }
        }
      
      // *** Median size ***
      else if ( something == "median size" )
        {
          if( n0_depend  &&  la_depend )
            {          
              mub1   = mgd_pars[1] + scat_species_b + 1;
              if( mub1 <= 0 )
                throw runtime_error( "Bad MGD parameter detected: mu + b + 1 <= 0" );
              eterm  = mub1 / mgd_pars[3];          
              // Start by deriving la 
              scfac2 = ( mgd_pars[1] + 1 + scat_species_b - 0.327*mgd_pars[3] ) /
                mgd_pars[3];  
              mgd_pars[2] = scfac2 * pow( ext_pars[1], -mgd_pars[3] );
              // We can now derive n0 
              gterm = tgamma( eterm );
              scfac1 = ( mgd_pars[3] * pow( mgd_pars[2], eterm ) ) /
                ( scat_species_a * gterm );
              mgd_pars[0] = scfac1 * ext_pars[0];
            }
          else 
            { assert(0); }
        }
      
      // *** Mean particle size ***
      else if ( something == "mean particle mass" )
        {
          if( n0_depend  &&  la_depend )
            {
              mu1    = mgd_pars[1] + 1;
              if( mu1 <= 0 )
                throw runtime_error( "Bad MGD parameter detected: mu + 1 <= 0" );
              eterm  = ( mgd_pars[1] + scat_species_b + 1 ) / mgd_pars[3];          
              gterm  = tgamma( eterm );
              // Start by deriving la
              gab    = mgd_pars[3] / scat_species_b;
              scfac2 = pow( scat_species_a * gterm /
                            tgamma( mu1/mgd_pars[3] ), gab ); 
              mgd_pars[2] = scfac2 * pow( ext_pars[1], -gab );
              // We can now derive n0 
              scfac1 = ( mgd_pars[3] * pow( mgd_pars[2], eterm ) ) /
                ( scat_species_a * gterm );
              mgd_pars[0] = scfac1 * ext_pars[0];
            }
          else 
            { assert(0); }
        }
      
      else if( something == "Ntot" )
        {
          if( n0_depend  &&  la_depend )
            {
              mu1    = mgd_pars[1] + 1;
              if( mu1 <= 0 )
                throw runtime_error( "Bad MGD parameter detected: mu + 1 <= 0" );
              eterm  = ( mgd_pars[1] + scat_species_b + 1 ) / mgd_pars[3];          
              gterm  = tgamma( eterm );
              // Start by deriving la
              gab    = mgd_pars[3] / scat_species_b;
              scfac2 = pow( scat_species_a * gterm /
                            tgamma( mu1/mgd_pars[3] ), gab ); 
              mgd_pars[2] = scfac2 * pow( ext_pars[1]/ext_pars[0], gab );
              // We can now derive n0 
              scfac1 = ( mgd_pars[3] * pow( mgd_pars[2], eterm ) ) /
                ( scat_species_a * gterm );
              mgd_pars[0] = scfac1 * ext_pars[0];
            }
          else 
            { assert(0); }
        }
      
      // String something not recognised
      else
        { assert(0); }
      
      // Now when all four MGD parameters are set, check that la and ga are OK
      if( mgd_pars[2] <= 0 )
        throw runtime_error( "Bad MGD parameter detected: la <= 0" );
      if( mgd_pars[3] <= 0 )
        throw runtime_error( "Bad MGD parameter detected: ga <= 0" );
      
      // Calculate PSD and derivatives
      Matrix  jac_data(4,nsi);    
      mgd_with_derivatives( psd_data(ip,joker), jac_data, psd_size_grid,
           mgd_pars[0], mgd_pars[1], mgd_pars[2], mgd_pars[3],
           (bool)mgd_do_jac[0] || n0_depend,
           (bool)mgd_do_jac[1] || mu_depend,
           (bool)mgd_do_jac[2] || la_depend,
           (bool)mgd_do_jac[3] || ga_depend );
      
      // Derivatives for mass and something
      if( ext_do_jac[0] | ext_do_jac[1] )
        {
          // *** Mean size ***
          if ( something == "mean size" )
            {
              if( n0_depend  &&  la_depend )
                {
                  // Derivative with respect to mass
                  if( ext_do_jac[0] )
                    {
                      dpsd_data_dx(ext_i_jac[0],ip,joker) = jac_data(0,joker);
                      dpsd_data_dx(ext_i_jac[0],ip,joker) *= scfac1;
                    }
                  // Derivative with respect to mean size
                  if( ext_do_jac[1] )
                    {
                      // 1. Term associated with n0
                      // Calculated as dpsd/dn0 * dn0/dla * dla/dXm
                      dpsd_data_dx(ext_i_jac[1],ip,joker) = jac_data(0,joker);
                      dpsd_data_dx(ext_i_jac[1],ip,joker) *= ext_pars[0] *
                        mgd_pars[3] * eterm * pow(mgd_pars[2],eterm-1) /
                        ( scat_species_a * gterm );
                      // 2. Term associated with la
                      // Calculated as dpsd/dla * dla/dXm
                      dpsd_data_dx(ext_i_jac[1],ip,joker) += jac_data(2,joker);
                      // Apply dla/dXm to sum
                      dpsd_data_dx(ext_i_jac[1],ip,joker) *= -mgd_pars[3] * scfac2 *
                        pow( ext_pars[1], -(mgd_pars[3]+1) );
                    }
                }
              else 
                { assert(0); }              
            }
          
          // *** Median size ***
          else if ( something == "median size" )
            {
              if( n0_depend  &&  la_depend )
                {
                  // Derivative with respect to mass
                  if( ext_do_jac[0] )
                    {
                      dpsd_data_dx(ext_i_jac[0],ip,joker) = jac_data(0,joker);
                      dpsd_data_dx(ext_i_jac[0],ip,joker) *= scfac1;
                    }
                  // Derivative with respect to median size
                  if( ext_do_jac[1] )
                    {
                      // 1. Term associated with n0
                      // Calculated as dpsd/dn0 * dn0/dla * dla/dXm
                      dpsd_data_dx(ext_i_jac[1],ip,joker) = jac_data(0,joker);
                      dpsd_data_dx(ext_i_jac[1],ip,joker) *= ext_pars[0] *
                        mgd_pars[3] * eterm * pow(mgd_pars[2],eterm-1) /
                        ( scat_species_a * gterm );
                      // 2. Term associated with la
                      // Calculated as dpsd/dla * dla/dXm
                      dpsd_data_dx(ext_i_jac[1],ip,joker) += jac_data(2,joker);
                      // Apply dla/dXm to sum
                      dpsd_data_dx(ext_i_jac[1],ip,joker) *= -mgd_pars[3] * scfac2 *
                        pow( ext_pars[1], -(mgd_pars[3]+1) );
                    }
                }
              else 
                { assert(0); }
            }
          
          // *** Mean particle size ***
          else if ( something == "mean particle mass" )
            {
              if( n0_depend  &&  la_depend )
                {
                  // Derivative with respect to mass
                  if( ext_do_jac[0] )
                    {
                      dpsd_data_dx(ext_i_jac[0],ip,joker) = jac_data(0,joker);
                      dpsd_data_dx(ext_i_jac[0],ip,joker) *= scfac1;
                    }
                  // Derivative with respect to mean particle size
                  if( ext_do_jac[1] )
                    {
                      // 1. Term associated with n0
                      // Calculated as dpsd/dn0 * dn0/dla * dla/dMm
                      dpsd_data_dx(ext_i_jac[1],ip,joker) = jac_data(0,joker);
                      dpsd_data_dx(ext_i_jac[1],ip,joker) *= ext_pars[0] *
                        mgd_pars[3] * eterm * pow(mgd_pars[2],eterm-1) /
                        ( scat_species_a * gterm );
                      // 2. Term associated with la
                      // Calculated as dpsd/dla * dla/dMm
                      dpsd_data_dx(ext_i_jac[1],ip,joker) += jac_data(2,joker);
                      // Apply dla/dMm to sum
                      dpsd_data_dx(ext_i_jac[1],ip,joker) *= scfac2 *
                        ( -mgd_pars[3] / scat_species_b ) * pow( ext_pars[1], -(gab+1) );
                    }
                  else 
                    { assert(0); }
                }
            }
          
          else if( something == "Ntot" )
            {
              if( n0_depend  &&  la_depend )
                {
                  // Term part of both derivatives
                  const Numeric dn0dla = ext_pars[0] * mgd_pars[3] * eterm *
                    pow(mgd_pars[2],eterm-1) / ( scat_species_a * gterm );
                  // Derivative with respect to mass
                  if( ext_do_jac[0] )
                    {
                      // Repeated term
                      const Numeric dladw = scfac2 * pow( ext_pars[1], gab ) *
                        ( -mgd_pars[3] / scat_species_b ) *
                        pow( ext_pars[0], -(gab+1) );
                      // 1. Term associated with n0
                      dpsd_data_dx(ext_i_jac[0],ip,joker) = jac_data(0,joker);
                      dpsd_data_dx(ext_i_jac[0],ip,joker) *= scfac1 + dn0dla*dladw;
                      // 2. Term associated with la
                      Vector term2 = jac_data(2,joker);
                      term2 *= dladw;
                      // Sum up
                      dpsd_data_dx(ext_i_jac[0],ip,joker) += term2;                      
                    }
                  // Derivative with respect to Ntot
                  if( ext_do_jac[1] )
                    {
                      // 1. Term associated with n0
                      // Calculated as dpsd/dn0 * dn0/dla * dla/dNtot
                      dpsd_data_dx(ext_i_jac[1],ip,joker) = jac_data(0,joker);
                      dpsd_data_dx(ext_i_jac[1],ip,joker) *= dn0dla;
                      // 2. Term associated with la
                      // Calculated as dpsd/dla * dla/dNtot
                      dpsd_data_dx(ext_i_jac[1],ip,joker) += jac_data(2,joker);
                      // Apply dla/dNtot to sum
                      dpsd_data_dx(ext_i_jac[1],ip,joker) *= scfac2 *
                        pow( ext_pars[0], -gab ) *
                        ( mgd_pars[3] / scat_species_b ) *
                        pow( ext_pars[1], gab-1 );
                    }
                }
              else 
                { assert(0); }
            }

          // String something not recognised
          else
            { assert(0); }
        }

      // Derivatives for non-dependent native parameters
      for( Index i=0; i<4; i++ )
        {
          if( mgd_do_jac[i] )
            { dpsd_data_dx(mgd_i_jac[i],ip,joker) = jac_data(i,joker); }
        } 
    }  
}



/* Workspace method: Doxygen documentation will be auto-generated */
void psdMgdMassNtot(
          Matrix&          psd_data,
          Tensor3&         dpsd_data_dx,
    const Vector&          psd_size_grid,
    const Vector&          pnd_agenda_input_t,
    const Matrix&          pnd_agenda_input,
    const ArrayOfString&   pnd_agenda_input_names,
    const ArrayOfString&   dpnd_data_dx_names,
    const Numeric&         scat_species_a, 
    const Numeric&         scat_species_b, 
    const Numeric&         n0, 
    const Numeric&         mu, 
    const Numeric&         la, 
    const Numeric&         ga, 
    const Numeric&         t_min, 
    const Numeric&         t_max, 
    const Index&           picky, 
    const Verbosity&       verbosity )
{
  psd_mgd_mass_and_something( psd_data, dpsd_data_dx,
                              "Ntot",
                              psd_size_grid, pnd_agenda_input_t,
                              pnd_agenda_input, pnd_agenda_input_names,
                              dpnd_data_dx_names, scat_species_a, scat_species_b,
                              n0, mu, la, ga, t_min, t_max, picky, verbosity );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void psdMgdMassMeanParticleMass(
          Matrix&          psd_data,
          Tensor3&         dpsd_data_dx,
    const Vector&          psd_size_grid,
    const Vector&          pnd_agenda_input_t,
    const Matrix&          pnd_agenda_input,
    const ArrayOfString&   pnd_agenda_input_names,
    const ArrayOfString&   dpnd_data_dx_names,
    const Numeric&         scat_species_a, 
    const Numeric&         scat_species_b, 
    const Numeric&         n0, 
    const Numeric&         mu, 
    const Numeric&         la, 
    const Numeric&         ga, 
    const Numeric&         t_min, 
    const Numeric&         t_max, 
    const Index&           picky, 
    const Verbosity&       verbosity )
{
  psd_mgd_mass_and_something( psd_data, dpsd_data_dx,
                              "mean particle mass",
                              psd_size_grid, pnd_agenda_input_t,
                              pnd_agenda_input, pnd_agenda_input_names,
                              dpnd_data_dx_names, scat_species_a, scat_species_b,
                              n0, mu, la, ga, t_min, t_max, picky, verbosity );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void psdMgdMassXmean(
          Matrix&          psd_data,
          Tensor3&         dpsd_data_dx,
    const Vector&          psd_size_grid,
    const Vector&          pnd_agenda_input_t,
    const Matrix&          pnd_agenda_input,
    const ArrayOfString&   pnd_agenda_input_names,
    const ArrayOfString&   dpnd_data_dx_names,
    const Numeric&         scat_species_a, 
    const Numeric&         scat_species_b, 
    const Numeric&         n0, 
    const Numeric&         mu, 
    const Numeric&         la, 
    const Numeric&         ga, 
    const Numeric&         t_min, 
    const Numeric&         t_max, 
    const Index&           picky, 
    const Verbosity&       verbosity )
{
  psd_mgd_mass_and_something( psd_data, dpsd_data_dx,
                              "mean size",
                              psd_size_grid, pnd_agenda_input_t,
                              pnd_agenda_input, pnd_agenda_input_names,
                              dpnd_data_dx_names, scat_species_a, scat_species_b,
                              n0, mu, la, ga, t_min, t_max, picky, verbosity );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void psdMgdMassXmedian(
          Matrix&          psd_data,
          Tensor3&         dpsd_data_dx,
    const Vector&          psd_size_grid,
    const Vector&          pnd_agenda_input_t,
    const Matrix&          pnd_agenda_input,
    const ArrayOfString&   pnd_agenda_input_names,
    const ArrayOfString&   dpnd_data_dx_names,
    const Numeric&         scat_species_a, 
    const Numeric&         scat_species_b, 
    const Numeric&         n0, 
    const Numeric&         mu, 
    const Numeric&         la, 
    const Numeric&         ga, 
    const Numeric&         t_min, 
    const Numeric&         t_max, 
    const Index&           picky, 
    const Verbosity&       verbosity )
{
  psd_mgd_mass_and_something( psd_data, dpsd_data_dx,
                              "median size",
                              psd_size_grid, pnd_agenda_input_t,
                              pnd_agenda_input, pnd_agenda_input_names,
                              dpnd_data_dx_names, scat_species_a, scat_species_b,
                              n0, mu, la, ga, t_min, t_max, picky, verbosity );
}

extern const Numeric PI;

Numeric dm_from_iwc_n0(Numeric iwc, Numeric n0, Numeric rho) {
    return pow(256.0 * iwc / PI / rho / n0, 0.25);
}

Numeric n0_from_iwc_dm(Numeric iwc, Numeric dm, Numeric rho) {
    return 256.0 * iwc / PI / rho / pow(dm, 4.0);
}

Numeric n0_from_t(Numeric t) {
    return exp(-0.076586 * (t - 273.15) + 17.948);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void psdD14(
    Matrix&                psd_data,
    Tensor3&               dpsd_data_dx,
    const Vector&          psd_size_grid,
    const Vector&          pnd_agenda_input_t,
    const Matrix&          pnd_agenda_input,
    const ArrayOfString&   pnd_agenda_input_names,
    const ArrayOfString&   dpnd_data_dx_names,
    const Numeric&         iwc,
    const Numeric&         n0,
    const Numeric&         dm,
    const Numeric&         rho,
    const Numeric&         alpha,
    const Numeric&         beta,
    const Numeric&         t_min,
    const Numeric&         t_max,
    const Index&           picky,
    const Verbosity& )
{
    // Standard checks
    START_OF_PSD_METHODS();

    // Additional (basic) checks
    if( nin > 2 )
        throw runtime_error( "The number of columns in *pnd_agenda_input* must "
                             "be 0, 1 or 2" );


    // Check and determine dependent and fixed parameters
    const bool n0_depend = (Index) n0 == -999;
    const bool dm_depend = (Index) dm == -999;

    // Check fixed parameters
    const bool iwc_fixed = !(isnan(iwc));
    const bool n0_fixed  = !(isnan(n0)) && !n0_depend;
    const bool dm_fixed  = !(isnan(dm)) && !dm_depend;

    if (!((nin + iwc_fixed + n0_fixed + dm_fixed == 2)
          || (nin + iwc_fixed + n0_fixed + dm_fixed == 1))) {
        throw runtime_error("This PSD can have one or two independent parameters, that is \n"
                            "the sum of the number of rows in pnd_agenda_input and\n"
                            "non-NAN, non-dependent values in iwc, n0, dm must be equal to\n"
                            "one or two.");
    }

    ArrayOfIndex i_pai = {-1,-1,-1}; // Position in pnd_agenda_input

    Index nhit = 0;

    if ((n0_depend || dm_depend) && (!iwc_fixed)) {
        i_pai[0] = nhit++;
    }

    if ((!n0_depend) && (!n0_fixed)) {
        i_pai[1] = nhit++;
    }

    if ((!dm_depend) && (!dm_fixed)) {
        i_pai[2] = nhit++;
    }

    // Determine what derivatives to do and their positions
    ArrayOfIndex do_jac = {0, 0, 0};
    ArrayOfIndex i_jac  = {-1, -1, -1}; // Position among jacobian quantities

    for (Index i=0; i < ndx; ++i) {
        for (Index j = 0; j < 3; ++j) {
            if (dx2in[i] == i_pai[j] ) {
                do_jac[j] = 1;
                i_jac[j]  = i;
                break;
            }
        }
    }

    if (psd_size_grid[0] < std::numeric_limits<Numeric>::epsilon()) {
        if (psd_size_grid.nelem() < 2) {
            throw std::runtime_error("psd_size_grid has only one element which is 0. This is not allowed.");
        }
    }

    Numeric iwc_p(0.0), n0_p(0.0), dm_p(0.0);
    // Loop input data and calculate PSDs
    for(Index ip=0; ip<np; ip++) {

        Numeric t = pnd_agenda_input_t[ip];

        // Extract MGD parameters
        if (i_pai[0] >= 0) {
            iwc_p = pnd_agenda_input(ip, i_pai[0]);
        }
        if (i_pai[1] >= 0) {
            n0_p = pnd_agenda_input(ip, i_pai[1]);
        }
        if (i_pai[2] >= 0) {
            dm_p = pnd_agenda_input(ip, i_pai[2]);
        }

        if (n0_depend && dm_depend) {
            n0_p = n0_from_t(t);
            dm_p = dm_from_iwc_n0(iwc_p, n0_p, rho);
        } else if (n0_depend) {
            n0_p = n0_from_iwc_dm(iwc_p, dm_p, rho);
        } else if (dm_depend) {
            dm_p = dm_from_iwc_n0(iwc_p, n0_p, rho);
        }

        // Outside of [t_min,tmax]?
        if ((t < t_min)  ||  (t > t_max)) {
            if(picky) {
                ostringstream os;
                os << "Method called with a temperature of " << t << " K.\n"
                   << "This is outside the specified allowed range: [ max(0.,"
                   << t_min << "), " << t_max << " ]";
                throw runtime_error(os.str());
            } else {
                continue;
            }
        }

        // Calculate PSD and derivatives
        Matrix  jac_data(1, nsi);
        Vector  x_grid(psd_size_grid);
        x_grid *= 1.0 / dm_p;

        if (x_grid[0] < std::numeric_limits<Numeric>::epsilon()) {
            x_grid[0] = 0.1 * psd_size_grid[1];
        }

        delanoe_shape_with_derivative(psd_data(ip,joker), jac_data,
                                      x_grid, alpha, beta);
        psd_data(ip, joker) *= n0_p;
        jac_data(0, joker)  *= n0_p;

        Vector dndx  = jac_data(0, joker);
        Vector dndn0 = psd_data(ip, joker);
        dndn0 *= (1.0 / n0_p);

        Vector dxddm = x_grid;
        dxddm *= (-1.0 / dm_p);
        Numeric dn0diwc = n0_p / iwc_p;

        if( do_jac[0] ) {

            dpsd_data_dx(i_jac[0], ip, joker) = 0.0;

            if (dm_depend) {

                Numeric ddmdiwc  = 0.25 * dm_p / iwc_p;
                Vector dndiwc = dxddm;
                dndiwc *= dndx;
                dndiwc *= ddmdiwc;
                dpsd_data_dx(i_jac[0],ip,joker) += dndiwc;
            } else if (n0_depend) {
                Vector dndiwc = dndn0;
                dndiwc *= dn0diwc;
                dpsd_data_dx(i_jac[0],ip,joker) += dndiwc;
            }
        }

        if( do_jac[1] ) {
            dpsd_data_dx(i_jac[1], ip, joker) = dndn0;
            if (dm_depend) {
                Vector dndn02 = dndx;
                dndn02 *= dxddm;
                Numeric ddmdn0 = - 0.25 / n0_p * dm_p;
                dndn02 *= ddmdn0;
                dpsd_data_dx(i_jac[1], ip, joker) += dndn02;
            }
        }

        if (do_jac[2]) {
            dpsd_data_dx(i_jac[2], ip, joker) = dxddm;
            dpsd_data_dx(i_jac[2], ip, joker) *= dndx;
            if (n0_depend) {
                Vector dndn02 = dndn0;
                Numeric dn0ddm = - 4.0 * n0_p / dm_p;
                dndn02 *= dn0ddm;
                dpsd_data_dx(i_jac[2], ip, joker) += dndn02;
            }
        }
    }
}

/*===========================================================================
  === Input: IWC and T
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void psdF07 (
          Matrix&          psd_data,
          Tensor3&         dpsd_data_dx,
    const Vector&          psd_size_grid,
    const Vector&          pnd_agenda_input_t,
    const Matrix&          pnd_agenda_input,
    const ArrayOfString&   pnd_agenda_input_names,
    const ArrayOfString&   dpnd_data_dx_names,
    const Numeric&         scat_species_a, 
    const Numeric&         scat_species_b, 
    const String&          regime,
    const Numeric&         t_min,
    const Numeric&         t_max,
    const Numeric&         t_min_psd,
    const Numeric&         t_max_psd,
    const Numeric&         b_min, 
    const Numeric&         b_max, 
    const Index&           picky, 
    const Verbosity&)
{
  // Standard checcks
  START_OF_PSD_METHODS();

  // Additional (basic) checks
  if( pnd_agenda_input.ncols() != 1 ) 
    throw runtime_error( "*pnd_agenda_input* must have one column." ); 
  if( regime!="TR" && regime!="ML" )
    throw runtime_error( "regime must either be \"TR\" or \"ML\"." );
  if( scat_species_a <= 0 )
    throw runtime_error( "*scat_species_a* should be > 0." );
  if( scat_species_b < b_min  ||  scat_species_b > b_max )
    {
      ostringstream os;
      os << "Method called with a mass-dimension-relation exponent b of "
         << scat_species_b << ".\n"
         << "This is outside the specified allowed range: ["
         << b_min << "," << b_max << "]";
      throw runtime_error(os.str());
    }

  
  for( Index ip=0; ip<np; ip++ )
    {
      // Extract the input variables
      Numeric swc = pnd_agenda_input(ip,0);
      Numeric   t = pnd_agenda_input_t[ip];

      // NaN can be generated for extremly small SWC (such as 1e-78)
      // As a solution set a limit, and consider all below as zero
      if( abs(swc) < 1e-15 )
        { swc = 0.0; }
      
      // No calc needed if swc==0 and no jacobians requested.
      if( (swc==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: [ max(0.,"
                 << t_min << "), " << t_max << " ]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }

      // PSD assumed to be constant outside [*t_min_psd*,*t_max_psd*]
      if( t < t_min_psd )
        { t = t_min_psd; }
      else if( t > t_max_psd )
        { t = t_max_psd; }
    
      // Negative swc?
      Numeric psd_weight = 1.0;
      if( swc < 0 )
        {
          psd_weight = -1.0;
          swc       *= -1.0;
        }
      
      // Calculate PSD
      Vector psd_1p(nsi);
      if( swc != 0 )
        {
          psd_snow_F07 ( psd_1p, psd_size_grid, swc, t, scat_species_a,
                         scat_species_b, regime );
          for ( Index i=0; i<nsi; i++ )
            { psd_data(ip,i) = psd_weight * psd_1p[i]; }
        }

      // Calculate derivative with respect to IWC
      if( ndx )
        {
          //const Numeric dswc = max( 0.001*swc, 1e-7 );
          const Numeric dswc = 1e-9;
          const Numeric swcp = swc + dswc;
          psd_snow_F07 ( psd_1p, psd_size_grid, swcp, t, scat_species_a,
                         scat_species_b, regime );
          for ( Index i=0; i<nsi; i++ )
            { dpsd_data_dx(0,ip,i) = ( psd_1p[i] - psd_data(ip,i) ) / dswc; }
        }   
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void psdMH97 (
          Matrix&           psd_data,
          Tensor3&          dpsd_data_dx,
    const Vector&           psd_size_grid,
    const Vector&           pnd_agenda_input_t,
    const Matrix&           pnd_agenda_input,
    const ArrayOfString&    pnd_agenda_input_names,
    const ArrayOfString&    dpnd_data_dx_names,
    const Numeric&          scat_species_a, 
    const Numeric&          scat_species_b, 
    const Numeric&          t_min, 
    const Numeric&          t_max, 
    const Numeric&          t_min_psd, 
    const Numeric&          t_max_psd, 
    const Index&            picky, 
    const Index&            noisy,
    const Verbosity&)
{
  // Standard checcks
  START_OF_PSD_METHODS();

  // Extra checks for this PSD
  if( pnd_agenda_input.ncols() != 1 ) 
    throw runtime_error( "*pnd_agenda_input* must have one column." ); 
  if( noisy  &&   ndx )
    throw runtime_error( "Jacobian calculations and \"noisy\" can not be "
                         "combined." );
  if( scat_species_b < 2.9  ||  scat_species_b > 3.1 )
    {
      ostringstream os;
      os << "This PSD treats pure ice, using Dveq as size grid.\n"
         << "This means that *scat_species_b* should be close to 3,\n"
         << "but it is outside of the tolerated range of [2.9,3.1].\n"
         << "Your value of *scat_species_b* is: " << scat_species_b;
      throw runtime_error(os.str());
    }
  if( scat_species_a < 460  ||  scat_species_a > 500 )
    {
      ostringstream os;
      os << "This PSD treats pure ice, using Dveq as size grid.\n"
         << "This means that *scat_species_a* should be close to 480,\n"
         << "but it is outside of the tolerated range of [460,500].\n"
         << "Your value of *scat_species_a* is: " << scat_species_a;
      throw runtime_error(os.str());
    }

  
  for( Index ip=0; ip<np; ip++ )
    {
      // Extract the input variables
      Numeric iwc = pnd_agenda_input(ip,0);
      Numeric   t = pnd_agenda_input_t[ip];

      // No calc needed if iwc==0 and no jacobians requested.
      if( (iwc==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: [ max(0.,"
                 << t_min << "), " << t_max << " ]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }

      // PSD assumed to be constant outside [*t_min_psd*,*t_max_psd*]
      if( t < t_min_psd )
        { t = t_min_psd; }
      else if( t > t_max_psd )
        { t = t_max_psd; }
  
      // Negative iwc?
      Numeric psd_weight = 1.0;
      if( iwc < 0 )
        {
          psd_weight = -1.0;
          iwc       *= -1.0;
        }
      
      // Calculate PSD
      Vector psd_1p(nsi);
      if( iwc != 0 )
        {
          psd_cloudice_MH97 ( psd_1p, psd_size_grid, iwc, t, noisy );
          for ( Index i=0; i<nsi; i++ )
            { psd_data(ip,i) = psd_weight * psd_1p[i]; }
        }

      // Calculate derivative with respect to IWC
      if( ndx )
        {
          //const Numeric diwc = max( 0.001*iwc, 1e-9 );
          const Numeric diwc = 1e-9;
          const Numeric iwcp = iwc + diwc;
          psd_cloudice_MH97 ( psd_1p, psd_size_grid, iwcp, t, noisy );
          for ( Index i=0; i<nsi; i++ )
            { dpsd_data_dx(0,ip,i) = ( psd_1p[i] - psd_data(ip,i) ) / diwc; }
        }   
    }
}




/*===========================================================================
  === Input: RWC
  ===========================================================================*/

/* Common code for PSDs taking just RWC as input */
void psd_rwc_common(
          Matrix&           psd_data,
          Tensor3&          dpsd_data_dx,
    const String&           psd_name,
    const Vector&           psd_size_grid,
    const Vector&           pnd_agenda_input_t,
    const Matrix&           pnd_agenda_input,
    const ArrayOfString&    pnd_agenda_input_names,
    const ArrayOfString&    dpnd_data_dx_names,
    const Numeric&          scat_species_a, 
    const Numeric&          scat_species_b, 
    const Numeric&          t_min, 
    const Numeric&          t_max, 
    const Index&            picky, 
    const Verbosity&)
{
  // Standard checks
  START_OF_PSD_METHODS();

  // Extra checks for this PSD
  if( pnd_agenda_input.ncols() != 1 )
    throw runtime_error( "*pnd_agenda_input* must have one column." );
  if( scat_species_b < 2.9  ||  scat_species_b > 3.1 )
    {
      ostringstream os;
      os << "This PSD treats rain, using Dveq as size grid.\n"
         << "This means that *scat_species_b* should be close to 3,\n"
         << "but it is outside of the tolerated range of [2.9,3.1].\n"
         << "Your value of *scat_species_b* is: " << scat_species_b;
      throw runtime_error(os.str());
    }
  if( scat_species_a < 500  ||  scat_species_a > 540 )
    {
      ostringstream os;
      os << "This PSD treats rain, using Dveq as size grid.\n"
         << "This means that *scat_species_a* should be close to 520,\n"
         << "but it is outside of the tolerated range of [500,540].\n"
         << "Your value of *scat_species_a* is: " << scat_species_a;
      throw runtime_error(os.str());
    }

  
  for( Index ip=0; ip<np; ip++ )
    {
      // Extract the input variables
      Numeric rwc = pnd_agenda_input(ip,0);
      Numeric   t = pnd_agenda_input_t[ip];

      // No calc needed if swc==0 and no jacobians requested.
      if( (rwc==0.) && (!ndx) )
        { continue; }   // If here, we are ready with this point!

      // Outside of [t_min,tmax]?
      if( t < t_min  ||  t > t_max )
        {
          if( picky )
            {
              ostringstream os;
              os << "Method called with a temperature of " << t << " K.\n"
                 << "This is outside the specified allowed range: [ max(0.,"
                 << t_min << "), " << t_max << " ]";
              throw runtime_error(os.str());
            }
          else  
            { continue; }   // If here, we are ready with this point!
        }
      
      // Negative rwc?
      Numeric psd_weight = 1.0;
      if( rwc < 0 )
        {
          psd_weight = -1.0;
          rwc       *= -1.0;
        }

      // Calculate PSD
      Vector psd_1p(nsi);
      if( rwc != 0 )
        {
          if( psd_name == "A12" )
            { psd_rain_A12( psd_1p, psd_size_grid, rwc ); }
          else if( psd_name == "W16" )
            { psd_rain_W16( psd_1p, psd_size_grid, rwc ); }
          else
            {assert(0); }
          //
          for ( Index i=0; i<nsi; i++ )
            { psd_data(ip,i) = psd_weight * psd_1p[i]; }
        }

      // Calculate derivative with respect to RWC
      if( ndx )
        {
          const Numeric drwc = 1e-9;
          const Numeric rwcp = rwc + drwc;
          if( psd_name == "A12" )
            { psd_rain_A12( psd_1p, psd_size_grid, rwcp ); }
          else if( psd_name == "W16" )
            { psd_rain_W16 ( psd_1p, psd_size_grid, rwcp ); }
          else
            {assert(0); }
          //
          for ( Index i=0; i<nsi; i++ )
            { dpsd_data_dx(0,ip,i) = ( psd_1p[i] - psd_data(ip,i) ) / drwc; }
        }   
    }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void psdA12(
          Matrix&           psd_data,
          Tensor3&          dpsd_data_dx,
    const Vector&           psd_size_grid,
    const Vector&           pnd_agenda_input_t,
    const Matrix&           pnd_agenda_input,
    const ArrayOfString&    pnd_agenda_input_names,
    const ArrayOfString&    dpnd_data_dx_names,
    const Numeric&          scat_species_a, 
    const Numeric&          scat_species_b, 
    const Numeric&          t_min, 
    const Numeric&          t_max, 
    const Index&            picky, 
    const Verbosity&        verbosity )
{
  psd_rwc_common( psd_data, dpsd_data_dx,
                  "A12",
                  psd_size_grid, pnd_agenda_input_t, pnd_agenda_input,
                  pnd_agenda_input_names, dpnd_data_dx_names,
                  scat_species_a,  scat_species_b,  t_min,  t_max,
                  picky, verbosity );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void psdW16(
          Matrix&           psd_data,
          Tensor3&          dpsd_data_dx,
    const Vector&           psd_size_grid,
    const Vector&           pnd_agenda_input_t,
    const Matrix&           pnd_agenda_input,
    const ArrayOfString&    pnd_agenda_input_names,
    const ArrayOfString&    dpnd_data_dx_names,
    const Numeric&          scat_species_a, 
    const Numeric&          scat_species_b, 
    const Numeric&          t_min, 
    const Numeric&          t_max, 
    const Index&            picky, 
    const Verbosity&        verbosity )
{
  psd_rwc_common( psd_data, dpsd_data_dx,
                  "W16",
                  psd_size_grid, pnd_agenda_input_t, pnd_agenda_input,
                  pnd_agenda_input_names, dpnd_data_dx_names,
                  scat_species_a,  scat_species_b,  t_min,  t_max,
                  picky, verbosity );
}

