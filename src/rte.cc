/* Copyright (C) 2002-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
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
  \file   rte.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2002-05-29

  \brief  Functions to solve radiative transfer tasks.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "auto_md.h"
#include "check_input.h"
#include "continua.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "montecarlo.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "refraction.h"
#include "rte.h"
#include "special_interp.h"

extern const String SURFACE_MAINTAG;
extern const String PROPMAT_SUBSUBTAG;
extern const String TEMPERATURE_MAINTAG;
extern const String SCATSPECIES_MAINTAG;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric TEMP_0_C;



/*===========================================================================
  === The functions in alphabetical order
  ===========================================================================*/


//! adjust_los
/*!
    Ensures that the zenith and azimuth angles of a line-of-sight vector are
    inside defined ranges.

    This function should not be used blindly, just when you know that the
    out-of-bounds values are obtained by an OK operation. As when making a
    disturbance calculation where e.g. the zenith angle is shifted with a small
    value. This function then handles the case when the original zenith angle
    is 0 or 180 and the disturbance then moves the angle outside the defined
    range. 

    \param   los              In/Out: LOS vector, defined as e.g. rte_los.
    \param   atmosphere_dim   As the WSV.

    \author Patrick Eriksson 
    \date   2012-04-11
*/
void adjust_los( 
         VectorView   los, 
   const Index &      atmosphere_dim )
{
  if( atmosphere_dim == 1 )
    {
           if( los[0] <   0 ) { los[0] = -los[0];    }
      else if( los[0] > 180 ) { los[0] = 360-los[0]; }
    }
  else if( atmosphere_dim == 2 )
    {
           if( los[0] < -180 ) { los[0] = los[0] + 360; }
      else if( los[0] >  180 ) { los[0] = los[0] - 360; }
    }
  else 
    {
      // If any of the angles out-of-bounds, use cart2zaaa to resolve 
      if( abs(los[0]-90) > 90  ||  abs(los[1]) > 180 )
        {
          Numeric dx, dy, dz;
          zaaa2cart( dx, dy, dz, los[0], los[1] );
          cart2zaaa( los[0], los[1], dx, dy, dz );
        }        
    }
}



//! apply_iy_unit
/*!
    Performs conversion from radiance to other units, as well as applies
    refractive index to fulfill the n2-law of radiance.

    Use *apply_iy_unit2* for conversion of jacobian data.

    \param   iy       In/Out: Tensor3 with data to be converted, where 
                      column dimension corresponds to Stokes dimensionality
                      and row dimension corresponds to frequency.
    \param   iy_unit  As the WSV.
    \param   f_grid   As the WSV.
    \param   n        Refractive index at the observation position.
    \param   i_pol    Polarisation indexes. See documentation of y_pol.

    \author Patrick Eriksson 
    \date   2010-04-07
*/
void apply_iy_unit( 
            MatrixView   iy, 
         const String&   iy_unit, 
       ConstVectorView   f_grid,
   const Numeric&        n,
   const ArrayOfIndex&   i_pol )
{
  // The code is largely identical between the two apply_iy_unit functions.
  // If any change here, remember to update the other function.

  const Index nf = iy.nrows();
  const Index ns = iy.ncols();

  assert( f_grid.nelem() == nf );
  assert( i_pol.nelem() == ns );

  if( iy_unit == "1" )
    {
      if( n != 1 )
        { iy *= (n*n); }
    }

  else if( iy_unit == "RJBT" )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          const Numeric scfac = invrayjean( 1, f_grid[iv] );
          for( Index is=0; is<ns; is++ )
            {
              if( i_pol[is] < 5 )           // Stokes components
                { iy(iv,is) *= scfac; }
              else                          // Measuement single pols
                { iy(iv,is) *= 2*scfac; }
            }
        }
    }

  else if( iy_unit == "PlanckBT" )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          for( Index is=ns-1; is>=0; is-- ) // Order must here be reversed
            {
              if( i_pol[is] == 1 )
                { iy(iv,is) = invplanck( iy(iv,is), f_grid[iv] ); }
              else if( i_pol[is] < 5 )
                { 
                  assert( i_pol[0] == 1 );
                  iy(iv,is) = 
                    invplanck( 0.5*(iy(iv,0)+iy(iv,is)), f_grid[iv] ) -
                    invplanck( 0.5*(iy(iv,0)-iy(iv,is)), f_grid[iv] );
                }
              else
                { iy(iv,is) = invplanck( 2*iy(iv,is), f_grid[iv] ); }
            }
        }
    }
  
  else if ( iy_unit == "W/(m^2 m sr)" )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          const Numeric scfac = n*n * f_grid[iv] * (f_grid[iv]/SPEED_OF_LIGHT);
          for( Index is=0; is<ns; is++ )
            { iy(iv,is) *= scfac; }
        }
    }
  
  else if ( iy_unit == "W/(m^2 m-1 sr)" )
    {
      iy *= ( n * n * SPEED_OF_LIGHT );
    }

  else
    {
      ostringstream os;
      os << "Unknown option: iy_unit = \"" << iy_unit << "\"\n" 
         << "Recognised choices are: \"1\", \"RJBT\", \"PlanckBT\""
         << "\"W/(m^2 m sr)\" and \"W/(m^2 m-1 sr)\""; 
      
      throw runtime_error( os.str() );      
    }
}



//! apply_iy_unit2
/*!
    Largely as *apply_iy_unit* but operates on jacobian data.

    The associated spectrum data *iy* must be in radiance. That is, the
    spectrum can only be converted to Tb after the jacobian data. 

    \param   J        In/Out: Tensor3 with data to be converted, where 
                      column dimension corresponds to Stokes dimensionality
                      and row dimension corresponds to frequency.
    \param   iy       Associated radiance data.
    \param   iy_unit  As the WSV.
    \param   f_grid   As the WSV.
    \param   n        Refractive index at the observation position.
    \param   i_pol    Polarisation indexes. See documentation of y_pol.

    \author Patrick Eriksson 
    \date   2010-04-10
*/
void apply_iy_unit2( 
   Tensor3View           J,
   ConstMatrixView       iy, 
   const String&         iy_unit, 
   ConstVectorView       f_grid,
   const Numeric&        n,
   const ArrayOfIndex&   i_pol )
{
  // The code is largely identical between the two apply_iy_unit functions.
  // If any change here, remember to update the other function.

  const Index nf = iy.nrows();
  const Index ns = iy.ncols();
  const Index np = J.npages();

  assert( J.nrows() == nf );
  assert( J.ncols() == ns );
  assert( f_grid.nelem() == nf );
  assert( i_pol.nelem() == ns );
  
  if( iy_unit == "1" )
    {
      if( n != 1 )
        { J *= (n*n); }
    }

  else if( iy_unit == "RJBT" )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          const Numeric scfac = invrayjean( 1, f_grid[iv] );
          for( Index is=0; is<ns; is++ )
            {
              if( i_pol[is] < 5 )           // Stokes componenets
                {
                  for( Index ip=0; ip<np; ip++ )
                    { J(ip,iv,is) *= scfac; }
                }
              else                          // Measuement single pols
                {
                  for( Index ip=0; ip<np; ip++ )
                    { J(ip,iv,is) *= 2*scfac; }
                }
            }
        }
    }

  else if( iy_unit == "PlanckBT" )
    {
      for( Index iv=0; iv<f_grid.nelem(); iv++ )
        {
          for( Index is=ns-1; is>=0; is-- )
            {
              Numeric scfac = 1;
              if( i_pol[is] == 1 )
                { scfac = dinvplanckdI( iy(iv,is), f_grid[iv] ); }
              else if( i_pol[is] < 5 )
                {
                  assert( i_pol[0] == 1 );
                  scfac = 
                    dinvplanckdI( 0.5*(iy(iv,0)+iy(iv,is)), f_grid[iv] ) +
                    dinvplanckdI( 0.5*(iy(iv,0)-iy(iv,is)), f_grid[iv] );
                }
              else
                { scfac = dinvplanckdI( 2*iy(iv,is), f_grid[iv] ); }
              //
              for( Index ip=0; ip<np; ip++ )
                { J(ip,iv,is) *= scfac; }
            }
        }
    }

  else if ( iy_unit == "W/(m^2 m sr)" )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          const Numeric scfac = n*n * f_grid[iv] * (f_grid[iv]/SPEED_OF_LIGHT);
          for( Index ip=0; ip<np; ip++ )
            {
              for( Index is=0; is<ns; is++ )
                { J(ip,iv,is) *= scfac; }
            }
        }
    }
  
  else if ( iy_unit == "W/(m^2 m-1 sr)" )
    {
      J *= ( n *n * SPEED_OF_LIGHT );
    }
  
  else
    {
      ostringstream os;
      os << "Unknown option: iy_unit = \"" << iy_unit << "\"\n" 
         << "Recognised choices are: \"1\", \"RJBT\", \"PlanckBT\""
         << "\"W/(m^2 m sr)\" and \"W/(m^2 m-1 sr)\""; 
      
      throw runtime_error( os.str() );      
    }  
}



//! ze_cfac
/*!
   Calculates factor to convert back-scattering to Ze

   The vector *fac* shall be sized to match f_grid, before calling the
   function.

   If k2 <= 0, the K" factor is calculated. Otherwise the input k2 is applied
   as "hard-coded".

   \param   fac      Vector with factors.
   \param   f_grid   As the WSV.
   \param   z_tref   Reference temperature for conversion to Ze.
   \param   k2       Reference dielectric factor.

   \author Patrick Eriksson
   \date   2002-05-20
*/
void ze_cfac(
         Vector&    fac,
   const Vector&    f_grid,
   const Numeric&   ze_tref,
   const Numeric&   k2 )
{
  const Index nf = f_grid.nelem();
  
  assert( fac.nelem() == nf );
  
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
      if( k2 >= 0 )
        { K2 = k2; }
      else
        {
          Complex n( complex_n(iv,0), complex_n(iv,1) );
          Complex n2 = n*n;
          Complex K  = ( n2 - Numeric(1.0) ) / ( n2 + Numeric(2.0) );
          Numeric absK = abs( K );
          K2 = absK * absK;
        }
      
      // Wavelength
      Numeric la = SPEED_OF_LIGHT / f_grid[iv];

      fac[iv] = a*la*la*la*la / K2; 
    }
}



//! bending_angle1d
/*!
    Calculates the bending angle for a 1D atmosphere.

    The expression used assumes a 1D atmosphere, that allows the bending angle
    to be calculated by start and end LOS. This is an approximation for 2D and
    3D, but a very small one and the function should in general be OK also for
    2D and 3D.

    The expression is taken from Kursinski et al., The GPS radio occultation
    technique, TAO, 2000.

    \return   alpha   Bending angle
    \param    ppath   Propagation path.

    \author Patrick Eriksson 
    \date   2012-04-05
*/
void bending_angle1d( 
        Numeric&   alpha,
  const Ppath&     ppath )
{
  Numeric theta;
  if( ppath.dim < 3 )
    { theta = abs( ppath.start_pos[1] - ppath.end_pos[1] ); }
  else
    { theta = sphdist( ppath.start_pos[1], ppath.start_pos[2],
                       ppath.end_pos[1], ppath.end_pos[2] ); }

  // Eq 17 in Kursinski et al., TAO, 2000:
  alpha = ppath.start_los[0] - ppath.end_los[0] + theta;

  // This as
  // phi_r = 180 - ppath.end_los[0]
  // phi_t = ppath.start_los[0]
}



//! defocusing_general_sub
/*!
    Just to avoid duplicatuion of code in *defocusing_general*.
   
    rte_los is mainly an input, but is also returned "adjusted" (with zenith
    and azimuth angles inside defined ranges) 
 
    \param    pos                 Out: Position of ppath at optical distance lo0
    \param    rte_los             In/out: Direction for transmitted signal 
                                  (disturbed from nominal value)
    \param    rte_pos             Out: Position of transmitter.
    \param    background          Out: Raditaive background of ppath.
    \param    lo0                 Optical path length between transmitter 
                                  and receiver.
    \param    ppath_step_agenda   As the WSV with the same name.
    \param    atmosphere_dim      As the WSV with the same name.
    \param    p_grid              As the WSV with the same name.
    \param    lat_grid            As the WSV with the same name.
    \param    lon_grid            As the WSV with the same name.
    \param    t_field             As the WSV with the same name.
    \param    z_field             As the WSV with the same name.
    \param    vmr_field           As the WSV with the same name.
    \param    f_grid              As the WSV with the same name.
    \param    refellipsoid        As the WSV with the same name.
    \param    z_surface           As the WSV with the same name.
    \param    verbosity           As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2012-04-11
*/
void defocusing_general_sub( 
        Workspace&   ws,
        Vector&      pos,
        Vector&      rte_los,
        Index&       background,
  ConstVectorView    rte_pos,
  const Numeric&     lo0,
  const Agenda&      ppath_step_agenda,
  const Numeric&     ppath_lmax,
  const Numeric&     ppath_lraytrace,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstVectorView    lat_grid,
  ConstVectorView    lon_grid,
  ConstTensor3View   t_field,
  ConstTensor3View   z_field,
  ConstTensor4View   vmr_field,
  ConstVectorView    f_grid,
  ConstVectorView    refellipsoid,
  ConstMatrixView    z_surface,
  const Verbosity&   verbosity )
{
  // Special treatment of 1D around zenith/nadir
  // (zenith angles outside [0,180] are changed by *adjust_los*)
  bool invert_lat = false;
  if( atmosphere_dim == 1  &&  ( rte_los[0] < 0 || rte_los[0] > 180 ) )
    { invert_lat = true; }

  // Handle cases where angles have moved out-of-bounds due to disturbance
  adjust_los( rte_los, atmosphere_dim );

  // Calculate the ppath for disturbed rte_los
  Ppath ppx;
  //
  ppath_calc( ws, ppx, ppath_step_agenda, atmosphere_dim, p_grid, lat_grid,
              lon_grid, t_field, z_field, vmr_field, 
              f_grid, refellipsoid, z_surface, 0, ArrayOfIndex(0), 
              rte_pos, rte_los, ppath_lmax, ppath_lraytrace, 0, verbosity );
  //
  background = ppath_what_background( ppx );

  // Calcualte cumulative optical path for ppx
  Vector lox( ppx.np );
  Index ilast = ppx.np-1;
  lox[0] = ppx.end_lstep;
  for( Index i=1; i<=ilast; i++ )
    { lox[i] = lox[i-1] + ppx.lstep[i-1] * ( ppx.nreal[i-1] + 
                                             ppx.nreal[i] ) / 2.0; }

  pos.resize( max( Index(2), atmosphere_dim ) );

  // Reciever at a longer distance (most likely out in space):
  if( lox[ilast] < lo0 )
    {
      const Numeric dl = lo0 - lox[ilast];
      if( atmosphere_dim < 3 )
        {
          Numeric x, z, dx, dz;
          poslos2cart( x, z, dx, dz, ppx.r[ilast], ppx.pos(ilast,1), 
                       ppx.los(ilast,0) );
          cart2pol( pos[0], pos[1], x+dl*dx, z+dl*dz, ppx.pos(ilast,1), 
                    ppx.los(ilast,0) );
        }
      else
        {
          Numeric x, y, z, dx, dy, dz;
          poslos2cart( x, y, z, dx, dy, dz, ppx.r[ilast], ppx.pos(ilast,1), 
                       ppx.pos(ilast,2), ppx.los(ilast,0), ppx.los(ilast,1) );
          cart2sph( pos[0], pos[1], pos[2], x+dl*dx, y+dl*dy, z+dl*dz, 
                    ppx.pos(ilast,1), ppx.pos(ilast,2), 
                    ppx.los(ilast,0), ppx.los(ilast,1) );
        }
    }

  // Interpolate to lo0
  else
    { 
      GridPos   gp;
      Vector    itw(2);
      gridpos( gp, lox, lo0 );
      interpweights( itw, gp );
      //
      pos[0] = interp( itw, ppx.r, gp );
      pos[1] = interp( itw, ppx.pos(joker,1), gp );
      if( atmosphere_dim == 3 )
        { pos[2] = interp( itw, ppx.pos(joker,2), gp ); }
    } 

  if( invert_lat )
    { pos[1] = -pos[1]; }
}


//! defocusing_general
/*!
    Defocusing for arbitrary geometry (zenith angle part only)

    Estimates the defocusing loss factor by calculating two paths with zenith
    angle off-sets. The distance between the two path at the optical path
    length between the transmitter and the receiver, divided with the
    corresponding distance for free space propagation, gives the defocusing
    loss. 

    The azimuth (gain) factor is not calculated. The path calculations are here
    done starting from the transmitter, which is the reversed direction
    compared to the ordinary path calculations starting at the receiver.
    
    \return   dlf                 Defocusing loss factor (1 for no loss)
    \param    ppath_step_agenda   As the WSV with the same name.
    \param    atmosphere_dim      As the WSV with the same name.
    \param    p_grid              As the WSV with the same name.
    \param    lat_grid            As the WSV with the same name.
    \param    lon_grid            As the WSV with the same name.
    \param    t_field             As the WSV with the same name.
    \param    z_field             As the WSV with the same name.
    \param    vmr_field           As the WSV with the same name.
    \param    f_grid              As the WSV with the same name.
    \param    refellipsoid        As the WSV with the same name.
    \param    z_surface           As the WSV with the same name.
    \param    ppath               As the WSV with the same name.
    \param    ppath_lmax          As the WSV with the same name.
    \param    ppath_lraytrace     As the WSV with the same name.
    \param    dza                 Size of angular shift to apply.
    \param    verbosity           As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2012-04-11
*/
void defocusing_general( 
        Workspace&   ws,
        Numeric&     dlf,
  const Agenda&      ppath_step_agenda,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstVectorView    lat_grid,
  ConstVectorView    lon_grid,
  ConstTensor3View   t_field,
  ConstTensor3View   z_field,
  ConstTensor4View   vmr_field,
  ConstVectorView    f_grid,
  ConstVectorView    refellipsoid,
  ConstMatrixView    z_surface,
  const Ppath&       ppath,
  const Numeric&     ppath_lmax,
  const Numeric&     ppath_lraytrace,
  const Numeric&     dza,
  const Verbosity&   verbosity )
{
  // Optical and physical path between transmitter and reciver
  Numeric lo = ppath.start_lstep + ppath.end_lstep;
  Numeric lp = lo;
  for( Index i=0; i<=ppath.np-2; i++ )
    { lp += ppath.lstep[i];
      lo += ppath.lstep[i] * ( ppath.nreal[i] + ppath.nreal[i+1] ) / 2.0; 
    }
  // Extract rte_pos and rte_los
  const Vector rte_pos = ppath.start_pos[Range(0,atmosphere_dim)];
  //
  Vector rte_los0(max(Index(1),atmosphere_dim-1)), rte_los;
  mirror_los( rte_los, ppath.start_los, atmosphere_dim );
  rte_los0 = rte_los[Range(0,max(Index(1),atmosphere_dim-1))];

  // A new ppath with positive zenith angle off-set
  //
  Vector  pos1;
  Index   backg1;
  //
  rte_los     = rte_los0;
  rte_los[0] += dza;
  //
  defocusing_general_sub( ws, pos1, rte_los, backg1, rte_pos, lo, 
                          ppath_step_agenda, ppath_lmax, ppath_lraytrace,
                          atmosphere_dim, p_grid, lat_grid, lon_grid, t_field, z_field, 
                          vmr_field, f_grid, refellipsoid, 
                          z_surface, verbosity );

  // Same thing with negative zenit angle off-set
  Vector  pos2;
  Index   backg2;
  //
  rte_los     = rte_los0;  // Use rte_los0 as rte_los can have been "adjusted"
  rte_los[0] -= dza;
  //
  defocusing_general_sub( ws, pos2, rte_los, backg2, rte_pos, lo, 
                          ppath_step_agenda, ppath_lmax, ppath_lraytrace,
                          atmosphere_dim, p_grid, lat_grid, lon_grid, t_field, z_field, 
                          vmr_field, f_grid, refellipsoid, 
                          z_surface, verbosity );

  // Calculate distance between pos1 and 2, and derive the loss factor
  // All appears OK:
  if( backg1 == backg2 )
    {
      Numeric l12;
      if( atmosphere_dim < 3 )
        { distance2D( l12, pos1[0], pos1[1], pos2[0], pos2[1] ); }
      else
        { distance3D( l12, pos1[0], pos1[1], pos1[2], 
                           pos2[0], pos2[1], pos2[2] ); }
      //
      dlf = lp*2*DEG2RAD*dza /  l12;
    }
  // If different backgrounds, then only use the second calculation
  else
    {
      Numeric l12;
      if( atmosphere_dim == 1 )
        { 
          const Numeric r = refellipsoid[0];
          distance2D( l12, r+ppath.end_pos[0], 0, pos2[0], pos2[1] ); 
        }
      else if( atmosphere_dim == 2 )
        { 
          const Numeric r = refell2r( refellipsoid, ppath.end_pos[1] );
          distance2D( l12, r+ppath.end_pos[0], ppath.end_pos[1], 
                                                        pos2[0], pos2[1] ); 
        }
      else
        { 
          const Numeric r = refell2r( refellipsoid, ppath.end_pos[1] );
          distance3D( l12, r+ppath.end_pos[0], ppath.end_pos[1], 
                             ppath.end_pos[2], pos2[0], pos2[1], pos2[2] ); 
        }
      //
      dlf = lp*DEG2RAD*dza /  l12;
    }
}



//! defocusing_sat2sat
/*!
    Calculates defocusing for limb measurements between two satellites.

    The expressions used assume a 1D atmosphere, and can only be applied on
    limb sounding geometry. The function works for 2D and 3D and should give 
    OK estimates. Both the zenith angle (loss) and azimuth angle (gain) terms
    are considered.

    The expressions is taken from Kursinski et al., The GPS radio occultation
    technique, TAO, 2000.

    \return   dlf                 Defocusing loss factor (1 for no loss)
    \param    ppath_step_agenda   As the WSV with the same name.
    \param    atmosphere_dim      As the WSV with the same name.
    \param    p_grid              As the WSV with the same name.
    \param    lat_grid            As the WSV with the same name.
    \param    lon_grid            As the WSV with the same name.
    \param    t_field             As the WSV with the same name.
    \param    z_field             As the WSV with the same name.
    \param    vmr_field           As the WSV with the same name.
    \param    f_grid              As the WSV with the same name.
    \param    refellipsoid        As the WSV with the same name.
    \param    z_surface           As the WSV with the same name.
    \param    ppath               As the WSV with the same name.
    \param    ppath_lmax          As the WSV with the same name.
    \param    ppath_lraytrace     As the WSV with the same name.
    \param    dza                 Size of angular shift to apply.
    \param    verbosity           As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2012-04-11
*/
void defocusing_sat2sat( 
        Workspace&   ws,
        Numeric&     dlf,
  const Agenda&      ppath_step_agenda,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstVectorView    lat_grid,
  ConstVectorView    lon_grid,
  ConstTensor3View   t_field,
  ConstTensor3View   z_field,
  ConstTensor4View   vmr_field,
  ConstVectorView    f_grid,
  ConstVectorView    refellipsoid,
  ConstMatrixView    z_surface,
  const Ppath&       ppath,
  const Numeric&     ppath_lmax,
  const Numeric&     ppath_lraytrace,
  const Numeric&     dza,
  const Verbosity&   verbosity )
{
  if( ppath.end_los[0] < 90  ||  ppath.start_los[0] > 90  )
     throw runtime_error( "The function *defocusing_sat2sat* can only be used "
                         "for limb sounding geometry." );

  // Index of tangent point
  Index it;
  find_tanpoint( it, ppath );
  assert( it >= 0 );

  // Length between tangent point and transmitter/reciver
  Numeric lt = ppath.start_lstep, lr = ppath.end_lstep;
  for( Index i=it; i<=ppath.np-2; i++ )
    { lt += ppath.lstep[i]; }
  for( Index i=0; i<it; i++ )
    { lr += ppath.lstep[i]; }

  // Bending angle and impact parameter for centre ray
  Numeric alpha0, a0;
  bending_angle1d( alpha0, ppath );
  alpha0 *= DEG2RAD;
  a0      = ppath.constant; 

  // Azimuth loss term (Eq 18.5 in Kursinski et al.)
  const Numeric lf = lr*lt / (lr+lt);
  const Numeric alt = 1 / ( 1 - alpha0*lf / refellipsoid[0] );

  // Calculate two new ppaths to get dalpha/da
  Numeric   alpha1, a1, alpha2, a2, dada;
  Ppath     ppt;
  Vector    rte_pos = ppath.end_pos[Range(0,atmosphere_dim)];
  Vector    rte_los = ppath.end_los;
  //
  rte_los[0] -= dza;
  adjust_los( rte_los, atmosphere_dim );
  ppath_calc( ws, ppt, ppath_step_agenda, atmosphere_dim, p_grid, lat_grid,
              lon_grid, t_field, z_field, vmr_field, 
              f_grid, refellipsoid, z_surface, 0, ArrayOfIndex(0), 
              rte_pos, rte_los, ppath_lmax, ppath_lraytrace, 0, verbosity );
  bending_angle1d( alpha2, ppt );
  alpha2 *= DEG2RAD;
  a2      = ppt.constant; 
  //
  rte_los[0] += 2*dza;
  adjust_los( rte_los, atmosphere_dim );
  ppath_calc( ws, ppt, ppath_step_agenda, atmosphere_dim, p_grid, lat_grid,
              lon_grid, t_field, z_field, vmr_field, 
              f_grid, refellipsoid, z_surface, 0, ArrayOfIndex(0), 
              rte_pos, rte_los, ppath_lmax, ppath_lraytrace, 0, verbosity );
  // This path can hit the surface. And we need to check if ppt is OK.
  // (remember this function only deals with sat-to-sat links and OK 
  // background here is be space) 
  // Otherwise use the centre ray as the second one.
  if( ppath_what_background(ppt) == 1 )
    {
      bending_angle1d( alpha1, ppt );
      alpha1 *= DEG2RAD;
      a1      = ppt.constant; 
      dada    = (alpha2-alpha1) / (a2-a1); 
    }
  else
    {
      dada    = (alpha2-alpha0) / (a2-a0);       
    }

  // Zenith loss term (Eq 18 in Kursinski et al.)
  const Numeric zlt = 1 / ( 1 - dada*lf );

  // Total defocusing loss
  dlf = zlt * alt;
}



//! dotprod_with_los
/*!
    Calculates the dot product between a field and a LOS
    The latter three gives the derivatives with respect to the LOS

    The line-of-sight shall be given as in the ppath structure (i.e. the
    viewing direction), but the dot product is calculated for the photon
    direction. The field is specified by its three components.

    The returned value can be written as |f|*cos(theta), where |f| is the field
    strength, and theta the angle between the field and photon vectors.

    \return                    The result of the dot product
    \param   los               Pppath line-of-sight.
    \param   u                 U-component of field.
    \param   v                 V-component of field.
    \param   w                 W-component of field.
    \param   atmosphere_dim    As the WSV.

    \author Patrick Eriksson 
    \date   2012-12-12
*/
Numeric dotprod_with_los(
  ConstVectorView   los, 
  const Numeric&    u,
  const Numeric&    v,
  const Numeric&    w,
  const Index&      atmosphere_dim )
{
  // Strength of field
  const Numeric f = sqrt( u*u + v*v + w*w );

  // Zenith and azimuth angle for field (in radians) 
  const Numeric za_f = acos( w/f );
  const Numeric aa_f = atan2( u, v );

  // Zenith and azimuth angle for photon direction (in radians)
  Vector los_p;
  mirror_los( los_p, los, atmosphere_dim );
  const Numeric za_p = DEG2RAD * los_p[0];
  const Numeric aa_p = DEG2RAD * los_p[1];
  
  return f * ( cos(za_f) * cos(za_p) +
               sin(za_f) * sin(za_p) * cos(aa_f-aa_p) );
}   


void ext_mat_case(Index& icase, ConstMatrixView ext_mat, const Index stokes_dim)
{
  if( icase == 0 )
  { 
    icase = 1;  // Start guess is diagonal
    
    //--- Scalar case ----------------------------------------------------------
    if( stokes_dim == 1 )
    {}
    
    //--- Vector RT ------------------------------------------------------------
    else
    {
      // Check symmetries and analyse structure of exp_mat:
      assert( ext_mat(1,1) == ext_mat(0,0) );
      assert( ext_mat(1,0) == ext_mat(0,1) );
      
      if( ext_mat(1,0) != 0 )
      { icase = 2; }
      
      if( stokes_dim >= 3 )
      {     
        assert( ext_mat(2,2) == ext_mat(0,0) );
        assert( ext_mat(2,1) == -ext_mat(1,2) );
        assert( ext_mat(2,0) == ext_mat(0,2) );
        
        if( ext_mat(2,0) != 0  ||  ext_mat(2,1) != 0 )
        { icase = 3; }
        
        if( stokes_dim > 3 )  
        {
          assert( ext_mat(3,3) == ext_mat(0,0) );
          assert( ext_mat(3,2) == -ext_mat(2,3) );
          assert( ext_mat(3,1) == -ext_mat(1,3) );
          assert( ext_mat(3,0) == ext_mat(0,3) ); 
          
          if( icase < 3 )  // if icase==3, already at most complex case
          {
            if( ext_mat(3,0) != 0  ||  ext_mat(3,1) != 0 )
            { icase = 3; }
            else if( ext_mat(3,2) != 0 )
            { icase = 2; }
          }
        }
      }
    }
  }
}

//! 
/*!
    Converts an extinction matrix to a transmission matrix

    The function performs the calculations differently depending on the
    conditions, to improve the speed. There are three cases: <br>
       1. Scalar RT and/or the matrix ext_mat_av is diagonal. <br>
       2. Special expression for "azimuthally_random" case. <br>
       3. The total general case.

    If the structure of *ext_mat* is known, *icase* can be set to "case index"
    (1, 2 or 3) and some time is saved. This includes that no asserts are
    performed on *ext_mat*.

    Otherwise, *icase* must be set to 0. *ext_mat* is then analysed and *icase*
    is set by the function and is returned.

    trans_mat must be sized before calling the function.

    \param   trans_mat          Input/Output: Transmission matrix of slab.
    \param   icase              Input/Output: Index giving ext_mat case.
    \param   ext_mat            Input: Averaged extinction matrix.
    \param   lstep              Input: The length of the RTE step.

    \author Patrick Eriksson (based on earlier version started by Claudia)
    \date   2013-05-17 
*/
void ext2trans(
         MatrixView   trans_mat,
         Index&       icase,
   ConstMatrixView    ext_mat,
   const Numeric&     lstep )
{
  const Index stokes_dim = ext_mat.ncols();

  assert( ext_mat.nrows()==stokes_dim );
  assert( trans_mat.nrows()==stokes_dim && trans_mat.ncols()==stokes_dim );

  // Theoretically ext_mat(0,0) >= 0, but to demand this can cause problems for
  // iterative retrievals, and the assert is skipped. Negative should be a
  // result of negative vmr, and this issue is checked in basics_checkedCalc.
  //assert( ext_mat(0,0) >= 0 );     

  assert( icase>=0 && icase<=3 );
  assert( !is_singular( ext_mat ) );
  assert( lstep >= 0 );
  
  // Analyse ext_mat?
  ext_mat_case(icase, ext_mat, stokes_dim);

  // Calculation options:
  if( icase == 1 )
    {
      trans_mat = 0;
      trans_mat(0,0) = exp( -ext_mat(0,0) * lstep );
      for( Index i=1; i<stokes_dim; i++ )
        { trans_mat(i,i) = trans_mat(0,0); }
    }
      
  else if( icase == 2  &&  stokes_dim < 3 )
    {
      // Expressions below are found in "Polarization in Spectral Lines" by
      // Landi Degl'Innocenti and Landolfi (2004).
      const Numeric tI = exp( -ext_mat(0,0) * lstep );
      const Numeric HQ = ext_mat(0,1) * lstep;
      trans_mat(0,0) = tI * cosh( HQ );
      trans_mat(1,1) = trans_mat(0,0);
      trans_mat(1,0) = -tI * sinh( HQ );
      trans_mat(0,1) = trans_mat(1,0);
      /* Does not work for stokes_dim==3, and commnted out 180502 by PE: 
      if( stokes_dim >= 3 )
        {
          trans_mat(2,0) = 0;
          trans_mat(2,1) = 0;
          trans_mat(0,2) = 0;
          trans_mat(1,2) = 0;
          const Numeric RQ = ext_mat(2,3) * lstep;
          trans_mat(2,2) = tI * cos( RQ );
          if( stokes_dim > 3 )
            {
              trans_mat(3,0) = 0;
              trans_mat(3,1) = 0;
              trans_mat(0,3) = 0;
              trans_mat(1,3) = 0;
              trans_mat(3,3) = trans_mat(2,2);
              trans_mat(3,2) = tI * sin( RQ );
              trans_mat(2,3) = -trans_mat(3,2); 
            }
        }
      */
    }
  else
    {
      Matrix ext_mat_ds = ext_mat;
      ext_mat_ds *= -lstep; 
      //         
      Index q = 10;  // index for the precision of the matrix exp function
      //
      switch( stokes_dim )
      {
        case 4:
          cayley_hamilton_fitted_method_4x4_propmat_to_transmat__eigen( trans_mat, ext_mat_ds );
          break;
        default :
          matrix_exp( trans_mat, ext_mat_ds, q );
      }
    }
}


//! get_iy
/*!
    Basic call of *iy_main_agenda*.

    This function is an interface to *iy_main_agenda* that can be used when
    only *iy* is of interest. That is, jacobian and auxilary parts are
    deactivated/ignored.

    \param   ws                    Out: The workspace
    \param   iy                    Out: As the WSV.
    \param   t_field               As the WSV.
    \param   z_field               As the WSV.
    \param   vmr_field             As the WSV.
    \param   cloudbox_on           As the WSV.
    \param   rte_pos               As the WSV.
    \param   rte_los               As the WSV.
    \param   iy_unit               As the WSV.
    \param   iy_main_agenda        As the WSV.

    \author Patrick Eriksson 
    \date   2012-08-08
*/
void get_iy(
         Workspace&   ws,
         Matrix&      iy,
   ConstTensor3View   t_field,
   ConstTensor3View   z_field,
   ConstTensor4View   vmr_field,
   ConstTensor4View   nlte_field,
   const Index&       cloudbox_on,
   ConstVectorView    f_grid,
   ConstVectorView    rte_pos,
   ConstVectorView    rte_los,
   ConstVectorView    rte_pos2,
   const String&      iy_unit,
   const Agenda&      iy_main_agenda )
{
  ArrayOfTensor3    diy_dx;
  ArrayOfMatrix     iy_aux;
  Ppath             ppath;
  Tensor3           iy_transmission(0,0,0);
  const Index       iy_agenda_call1 = 1;
  const ArrayOfString iy_aux_vars(0);
  const Index       iy_id = 0;
  const Index       jacobian_do = 0;

  iy_main_agendaExecute( ws, iy, iy_aux, ppath, diy_dx, 
                         iy_agenda_call1, iy_unit, iy_transmission, iy_aux_vars,
                         iy_id, cloudbox_on, jacobian_do, t_field, z_field,
                         vmr_field, nlte_field, f_grid, rte_pos, rte_los, rte_pos2,
                         iy_main_agenda );
}




//! get_iy_of_background
/*!
    Determines iy of the "background" of a propgation path.

    The task is to determine *iy* and related variables for the
    background, or to continue the raditiave calculations
    "backwards". The details here depends on the method selected for
    the agendas.

    Each background is handled by an agenda. Several of these agandes
    can involve recursive calls of *iy_main_agenda*. 

    \param   ws                    Out: The workspace
    \param   iy                    Out: As the WSV.
    \param   diy_dx                Out: As the WSV.
    \param   iy_transmission       As the WSV.
    \param   jacobian_do           As the WSV.
    \param   ppath                 As the WSV.
    \param   atmosphere_dim        As the WSV.
    \param   t_field               As the WSV.
    \param   z_field               As the WSV.
    \param   vmr_field             As the WSV.
    \param   cloudbox_on           As the WSV.
    \param   stokes_dim            As the WSV.
    \param   f_grid                As the WSV.
    \param   iy_unit               As the WSV.    
    \param   surface_props_data    As the WSV.    
    \param   iy_main_agenda        As the WSV.
    \param   iy_space_agenda       As the WSV.
    \param   iy_surface_agenda     As the WSV.
    \param   iy_cloudbox_agenda    As the WSV.

    \author Patrick Eriksson 
    \date   2009-10-08
*/
void get_iy_of_background(
        Workspace&        ws,
        Matrix&           iy,
        ArrayOfTensor3&   diy_dx,
  ConstTensor3View        iy_transmission,
  const Index&            iy_id, 
  const Index&            jacobian_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const Ppath&            ppath,
  ConstVectorView         rte_pos2,
  const Index&            atmosphere_dim,
  ConstTensor3View        t_field,
  ConstTensor3View        z_field,
  ConstTensor4View        vmr_field,
  const Index&            cloudbox_on,
  const Index&            stokes_dim,
  ConstVectorView         f_grid,
  const String&           iy_unit,  
  const Tensor3&          surface_props_data,
  const Agenda&           iy_main_agenda,
  const Agenda&           iy_space_agenda,
  const Agenda&           iy_surface_agenda,
  const Agenda&           iy_cloudbox_agenda,
  const Index&            iy_agenda_call1,
  const Verbosity&        verbosity)
{
  CREATE_OUT3;
  
  // Some sizes
  const Index nf = f_grid.nelem();
  const Index np = ppath.np;

  // Set rtp_pos and rtp_los to match the last point in ppath.
  //
  // Note that the Ppath positions (ppath.pos) for 1D have one column more
  // than expected by most functions. Only the first atmosphere_dim values
  // shall be copied.
  //
  Vector rtp_pos, rtp_los;
  rtp_pos.resize( atmosphere_dim );
  rtp_pos = ppath.pos(np-1,Range(0,atmosphere_dim));
  rtp_los.resize( ppath.los.ncols() );
  rtp_los = ppath.los(np-1,joker);

  out3 << "Radiative background: " << ppath.background << "\n";


  // Handle the different background cases
  //
  String agenda_name;
  // 
  switch( ppath_what_background( ppath ) )
    {

    case 1:   //--- Space ---------------------------------------------------- 
      {
        agenda_name = "iy_space_agenda";
        chk_not_empty( agenda_name, iy_space_agenda );
        iy_space_agendaExecute( ws, iy, f_grid, rtp_pos, rtp_los, 
                                iy_space_agenda );
      }
      break;

    case 2:   //--- The surface -----------------------------------------------
      {
        agenda_name = "iy_surface_agenda";
        chk_not_empty( agenda_name, iy_surface_agenda );
        //
        const Index los_id = iy_id % (Index)1000;
        Index iy_id_new = iy_id + (Index)9*los_id; 
        //
        // Surface jacobian stuff:
        ArrayOfString  dsurface_names(0);
        if( jacobian_do && iy_agenda_call1 )
          {
            for( Index i=0; i<jacobian_quantities.nelem(); i++ )
              {
                if( jacobian_quantities[i].MainTag() == SURFACE_MAINTAG )
                  { dsurface_names.push_back( jacobian_quantities[i].Subtag() ); }
              }
          }
        ArrayOfTensor4 dsurface_rmatrix_dx( dsurface_names.nelem() );
        ArrayOfMatrix  dsurface_emission_dx( dsurface_names.nelem() );
        //
        iy_surface_agendaExecute( ws, iy, diy_dx, dsurface_rmatrix_dx,
                                  dsurface_emission_dx,
                                  iy_unit, iy_transmission, iy_id_new, cloudbox_on,
                                  jacobian_do, t_field, z_field, vmr_field,
                                  f_grid, iy_main_agenda, rtp_pos, rtp_los, 
                                  rte_pos2, surface_props_data, dsurface_names,
                                  iy_surface_agenda );
      }
      break;

    case 3:   //--- Cloudbox boundary or interior ------------------------------
    case 4:
      {
        agenda_name = "iy_cloudbox_agenda";
        chk_not_empty( agenda_name, iy_cloudbox_agenda );
        iy_cloudbox_agendaExecute( ws, iy, f_grid, rtp_pos, rtp_los, 
                                   iy_cloudbox_agenda );
      }
      break;

    default:  //--- ????? ----------------------------------------------------
      // Are we here, the coding is wrong somewhere
      assert( false );
    }

  if( iy.ncols() != stokes_dim  ||  iy.nrows() != nf )
    {
      ostringstream os;
      os << "The size of *iy* returned from *" << agenda_name << "* is\n"
         << "not correct:\n"
         << "  expected size = [" << nf << "," << stokes_dim << "]\n"
         << "  size of iy    = [" << iy.nrows() << "," << iy.ncols()<< "]\n";
      throw runtime_error( os.str() );      
    }
}



//! get_ppath_atmvars
/*!
    Determines pressure, temperature, VMR, winds and magnetic field for each
    propgataion path point.

    The output variables are sized inside the function. For VMR the
    dimensions are [ species, propagation path point ].

    \param   ppath_p           Out: Pressure for each ppath point.
    \param   ppath_t           Out: Temperature for each ppath point.
    \param   ppath_vmr         Out: VMR values for each ppath point.
    \param   ppath_wind        Out: Wind vector for each ppath point.
    \param   ppath_mag         Out: Mag. field vector for each ppath point.
    \param   ppath             As the WSV.
    \param   atmosphere_dim    As the WSV.
    \param   p_grid            As the WSV.
    \param   lat_grid          As the WSV.
    \param   lon_grid          As the WSV.
    \param   t_field           As the WSV.
    \param   vmr_field         As the WSV.
    \param   wind_u_field      As the WSV.
    \param   wind_v_field      As the WSV.
    \param   wind_w_field      As the WSV.
    \param   mag_u_field       As the WSV.
    \param   mag_v_field       As the WSV.
    \param   mag_w_field       As the WSV.

    \author Patrick Eriksson 
    \date   2009-10-05
*/
void get_ppath_atmvars( 
        Vector&      ppath_p, 
        Vector&      ppath_t, 
        Matrix&      ppath_nlte,
        Matrix&      ppath_vmr, 
        Matrix&      ppath_wind, 
        Matrix&      ppath_mag,
  const Ppath&       ppath,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstTensor3View   t_field,
  ConstTensor4View   t_nlte_field,
  ConstTensor4View   vmr_field,
  ConstTensor3View   wind_u_field,
  ConstTensor3View   wind_v_field,
  ConstTensor3View   wind_w_field,
  ConstTensor3View   mag_u_field,
  ConstTensor3View   mag_v_field,
  ConstTensor3View   mag_w_field )
{
  const Index   np  = ppath.np;
  // Pressure:
  ppath_p.resize(np);
  Matrix itw_p(np,2);
  interpweights( itw_p, ppath.gp_p );      
  itw2p( ppath_p, p_grid, ppath.gp_p, itw_p );
  
  // Temperature:
  ppath_t.resize(np);
  Matrix   itw_field;
  interp_atmfield_gp2itw( itw_field, atmosphere_dim, 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon );
  interp_atmfield_by_itw( ppath_t,  atmosphere_dim, t_field, 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );

  // VMR fields:
  const Index ns = vmr_field.nbooks();
  ppath_vmr.resize(ns,np);
  for( Index is=0; is<ns; is++ )
    {
      interp_atmfield_by_itw( ppath_vmr(is, joker), atmosphere_dim,
                              vmr_field( is, joker, joker, joker ), 
                              ppath.gp_p, ppath.gp_lat, ppath.gp_lon, 
                              itw_field );
    }
    
  // NLTE temperatures
  const Index nnlte = t_nlte_field.nbooks();
  ppath_nlte.resize(nnlte,np);
  for( Index itnlte=0; itnlte<nnlte;itnlte++ )
  {
    interp_atmfield_by_itw( ppath_nlte(itnlte, joker),  atmosphere_dim, 
                            t_nlte_field( itnlte, joker, joker, joker ), 
                            ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
  }

  // Winds:
  ppath_wind.resize(3,np);
  ppath_wind = 0;
  //
  if( wind_u_field.npages() > 0 ) 
    { 
      interp_atmfield_by_itw( ppath_wind(0,joker), atmosphere_dim, wind_u_field,
                            ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }
  if( wind_v_field.npages() > 0 ) 
    { 
      interp_atmfield_by_itw( ppath_wind(1,joker), atmosphere_dim, wind_v_field,
                            ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }
  if( wind_w_field.npages() > 0 ) 
    { 
      interp_atmfield_by_itw( ppath_wind(2,joker), atmosphere_dim, wind_w_field,
                            ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }

  // Magnetic field:
  ppath_mag.resize(3,np);
  ppath_mag = 0;
  //
  if( mag_u_field.npages() > 0 )
    {
      interp_atmfield_by_itw( ppath_mag(0,joker), atmosphere_dim, mag_u_field,
                            ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }
  if( mag_v_field.npages() > 0 )
    {
      interp_atmfield_by_itw( ppath_mag(1,joker), atmosphere_dim, mag_v_field,
                            ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }
  if( mag_w_field.npages() > 0 )
    {
      interp_atmfield_by_itw( ppath_mag(2,joker), atmosphere_dim, mag_w_field,
                            ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }
}



//! get_ppath_cloudvars
/*!
    Determines the particle fields along a propagation path.

    \param   clear2cloudy        Out: Mapping of index. See code for details. 
    \param   ppath_pnd           Out: The particle number density for each
                                      path point (also outside cloudbox).
    \param   ppath_dpnd_dx       Out: dpnd_field_dx for each path point
                                      (also outside cloudbox).
    \param   ppath               As the WSV.    
    \param   cloubox_limits      As the WSV.    
    \param   pnd_field           As the WSV.    
    \param   dpnd_field_dx       As the WSV.    

    \author Jana Mendrok, Patrick Eriksson 
    \date   2017-09-18
*/
void get_ppath_cloudvars( 
        ArrayOfIndex&                  clear2cloudy,
        Matrix&                        ppath_pnd,
        ArrayOfMatrix&                 ppath_dpnd_dx,
  const Ppath&                         ppath,
  const Index&                         atmosphere_dim,
  const ArrayOfIndex&                  cloudbox_limits,
  const Tensor4&                       pnd_field,
  const ArrayOfTensor4&                dpnd_field_dx )
{
  const Index np = ppath.np;

  // Pnd along the ppath
  ppath_pnd.resize( pnd_field.nbooks(), np );
  ppath_pnd = 0;
  ppath_dpnd_dx.resize( dpnd_field_dx.nelem() );

  bool any_dpnd = false;
  for( Index iq=0; iq<dpnd_field_dx.nelem(); iq++ )
    {
      if( dpnd_field_dx[iq].empty() )
        { ppath_dpnd_dx[iq].resize( 0, 0 ); }
      else
        {
          any_dpnd = true;
          ppath_dpnd_dx[iq].resize( pnd_field.nbooks(), np );
        }
    }

  // A variable that can map from ppath to particle containers.
  // If outside cloudbox or all (d)pnd=0, this variable holds -1.
  clear2cloudy.resize( np );

  // Determine ppath_pnd and ppath_dpnd_dx
  Index nin = 0;
  for( Index ip=0; ip<np; ip++ ) // PPath point
    {
      Matrix itw( 1, Index(pow(2.0,Numeric(atmosphere_dim))) );

      ArrayOfGridPos gpc_p(1), gpc_lat(1), gpc_lon(1);
      GridPos gp_lat, gp_lon;
      if( atmosphere_dim >= 2 ) { gridpos_copy( gp_lat, ppath.gp_lat[ip] ); } 
      if( atmosphere_dim == 3 ) { gridpos_copy( gp_lon, ppath.gp_lon[ip] ); }

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
          bool any_ppath_dpnd = false;
          if( any_dpnd )
            {
              for( Index iq=0; iq<dpnd_field_dx.nelem(); iq++ ) // Jacobian parameter
                {
                  if( !dpnd_field_dx[iq].empty() )
                    {
                      for( Index i=0; i<pnd_field.nbooks(); i++ ) // Scattering element
                        { interp_atmfield_by_itw( ppath_dpnd_dx[iq](i,ip),
                                                  atmosphere_dim,
                                                  dpnd_field_dx[iq](i,joker,joker,joker), 
                                                  gpc_p, gpc_lat, gpc_lon, itw ); }
                      if( max(ppath_dpnd_dx[iq](joker,ip)) > 0. ||
                          min(ppath_dpnd_dx[iq](joker,ip)) < 0. )
                        any_ppath_dpnd = true;
                    }
                }
            }
          if( max(ppath_pnd(joker,ip)) > 0. || min(ppath_pnd(joker,ip)) < 0. ||
              any_ppath_dpnd )
            { clear2cloudy[ip] = nin;   nin++; }
          else
            { clear2cloudy[ip] = -1; }
        }
      else
        { clear2cloudy[ip] = -1; }
    }
}



//! get_ppath_f
/*!
    Determines the Doppler shifted frequencies along the propagation path.

    ppath_doppler [ nf + np ]

    \param   ppath_f          Out: Doppler shifted f_grid
    \param   ppath            Propagation path.
    \param   f_grid           Original f_grid.
    \param   atmosphere_dim   As the WSV.
    \param   rte_alonglos_v   As the WSV.
    \param   ppath_wind       See get_ppath_atmvars.

    \author Patrick Eriksson 
    \date   2013-02-21
*/
void get_ppath_f( 
        Matrix&    ppath_f,
  const Ppath&     ppath,
  ConstVectorView  f_grid, 
  const Index&     atmosphere_dim,
  const Numeric&   rte_alonglos_v,
  ConstMatrixView  ppath_wind )
{
  // Sizes
  const Index   nf = f_grid.nelem();
  const Index   np = ppath.np;

  ppath_f.resize(nf,np);

  // Doppler relevant velocity
  //
  for( Index ip=0; ip<np; ip++ )
    {
      // Start by adding rte_alonglos_v (most likely sensor effects)
      Numeric v_doppler = rte_alonglos_v;

      // Include wind
      if( ppath_wind(1,ip) != 0  ||  ppath_wind(0,ip) != 0  ||  
                                     ppath_wind(2,ip) != 0  )
        {
          // The dot product below is valid for the photon direction. Winds
          // along this direction gives a positive contribution.
          v_doppler += dotprod_with_los( ppath.los(ip,joker), ppath_wind(0,ip),
                          ppath_wind(1,ip), ppath_wind(2,ip), atmosphere_dim );
        }

      // Determine frequency grid
      if( v_doppler == 0 )
        { ppath_f(joker,ip) = f_grid; }
      else
        { 
          // Positive v_doppler means that sensor measures lower rest
          // frequencies
          const Numeric a = 1 - v_doppler / SPEED_OF_LIGHT;
          for( Index iv=0; iv<nf; iv++ )
            { ppath_f(iv,ip) = a * f_grid[iv]; }
        }
    }
}


//! get_rowindex_for_mblock
/*!
    Returns the "range" of *y* corresponding to a measurement block

    \return  The range.
    \param   sensor_response    As the WSV.
    \param   mblock_index            Index of the measurement block.

    \author Patrick Eriksson 
    \date   2009-10-16
*/
Range get_rowindex_for_mblock( 
  const Sparse&   sensor_response, 
  const Index&    mblock_index )
{
  const Index   n1y = sensor_response.nrows();
  return Range( n1y*mblock_index, n1y );
}


void iyb_calc_body(
        bool&                       failed,
        String&                     fail_msg,
        ArrayOfArrayOfMatrix&       iy_aux_array,
        Workspace&                  ws,
        Ppath&                      ppath,
        Vector&                     iyb,
        ArrayOfMatrix&              diyb_dx,
  const Index&                      mblock_index,
  const Index&                      atmosphere_dim,
  ConstTensor3View                  t_field,
  ConstTensor3View                  z_field,
  ConstTensor4View                  vmr_field,
  ConstTensor4View                  nlte_field,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  ConstVectorView                   f_grid,
  ConstMatrixView                   sensor_pos,
  ConstMatrixView                   sensor_los,
  ConstMatrixView                   transmitter_pos,
  ConstMatrixView                   mblock_dlos_grid,
  const String&                     iy_unit,  
  const Agenda&                     iy_main_agenda,
  const Index&                      j_analytical_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices,
  const ArrayOfString&              iy_aux_vars,
  const Index&                      ilos,
  const Index&                      nf )
{
    // The try block here is necessary to correctly handle
    // exceptions inside the parallel region.
    try
    {
      //--- LOS of interest
      //
      Vector los( sensor_los.ncols() );
      //
      los     = sensor_los( mblock_index, joker );
      los[0] += mblock_dlos_grid(ilos,0);
      if( mblock_dlos_grid.ncols() > 1 )
        { map_daa( los[0], los[1], los[0], los[1],
                   mblock_dlos_grid(ilos,1) ); }
      else
        { adjust_los( los, atmosphere_dim ); }
      //--- rtp_pos 1 and 2
      //
      Vector rtp_pos, rtp_pos2(0);
      //
      rtp_pos = sensor_pos( mblock_index, joker );
      if( transmitter_pos.nrows() )
        { rtp_pos2 = transmitter_pos( mblock_index, joker ); }

      // Calculate iy and associated variables
      //
      Matrix         iy;
      ArrayOfTensor3 diy_dx;
      Tensor3        iy_transmission(0,0,0);
      const Index    iy_agenda_call1 = 1;
      const Index    iy_id = (Index)1e6*(mblock_index+1) + (Index)1e3*(ilos+1);
      //
      iy_main_agendaExecute(ws, iy, iy_aux_array[ilos], ppath, diy_dx, 
                            iy_agenda_call1, iy_unit, iy_transmission, iy_aux_vars,
                            iy_id, cloudbox_on, j_analytical_do, t_field, z_field,
                            vmr_field, nlte_field, f_grid, rtp_pos, los,
                            rtp_pos2, iy_main_agenda );

      // Start row in iyb etc. for present LOS
      //
      const Index row0 = ilos * nf * stokes_dim;

      // Jacobian part
      //
      if( j_analytical_do )
        {
          FOR_ANALYTICAL_JACOBIANS_DO2
            (
                for( Index ip=0;
                           ip<jacobian_indices[iq][1] - jacobian_indices[iq][0]+1;
                           ip++ )
                  {
                    for( Index is=0; is<stokes_dim; is++ )
                      {
                        diyb_dx[iq](Range(row0+is,nf,stokes_dim),ip)=
                          diy_dx[iq](ip,joker,is);
                      }
                  }
            )
        }

      // iy : copy to iyb
      for( Index is=0; is<stokes_dim; is++ )
        { iyb[Range(row0+is,nf,stokes_dim)] = iy(joker,is); }

    }  // End try

    catch (const std::runtime_error &e)
    {
#pragma omp critical (iyb_calc_fail)
        { fail_msg = e.what(); failed = true; }
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
        ArrayOfVector&              iyb_aux,
        ArrayOfMatrix&              diyb_dx,
        Matrix&                     geo_pos_matrix,
  const Index&                      mblock_index,
  const Index&                      atmosphere_dim,
  ConstTensor3View                  t_field,
  ConstTensor3View                  z_field,
  ConstTensor4View                  vmr_field,
  ConstTensor4View                  nlte_field,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  ConstVectorView                   f_grid,
  ConstMatrixView                   sensor_pos,
  ConstMatrixView                   sensor_los,
  ConstMatrixView                   transmitter_pos,
  ConstMatrixView                   mblock_dlos_grid,
  const String&                     iy_unit,  
  const Agenda&                     iy_main_agenda,
  const Agenda&                     geo_pos_agenda,
  const Index&                      j_analytical_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices,
  const ArrayOfString&              iy_aux_vars,
  const Verbosity&                  verbosity)
{
  CREATE_OUT3;

  // Sizes
  const Index   nf   = f_grid.nelem();
  const Index   nlos = mblock_dlos_grid.nrows();
  const Index   niyb = nf * nlos * stokes_dim;
  // Set up size of containers for data of 1 measurement block.
  // (can not be made below due to parallalisation)
  iyb.resize( niyb );
  //
  if( j_analytical_do )
    {
      diyb_dx.resize( jacobian_indices.nelem() );
      FOR_ANALYTICAL_JACOBIANS_DO2(
        diyb_dx[iq].resize( niyb, jacobian_indices[iq][1] -
                                  jacobian_indices[iq][0] + 1 );
      )
    }
  else
    { diyb_dx.resize( 0 ); }
  // Assume that geo_pos_agenda returns empty geo_pos.
  geo_pos_matrix.resize( nlos, 5 );
  geo_pos_matrix = NAN;

  // For iy_aux we don't know the number of quantities, and we have to store
  // all outout
  ArrayOfArrayOfMatrix  iy_aux_array( nlos );

  // We have to make a local copy of the Workspace and the agendas because
  // only non-reference types can be declared firstprivate in OpenMP
  Workspace l_ws (ws);
  Agenda l_iy_main_agenda (iy_main_agenda);
  Agenda l_geo_pos_agenda (geo_pos_agenda);

  String fail_msg;
  bool failed = false;
  if (nlos >= arts_omp_get_max_threads() || nlos*10 >= nf)
    {
      out3 << "  Parallelizing los loop (" << nlos << " iterations, "
      << nf << " frequencies)\n";

      // Start of actual calculations
#pragma omp parallel for                   \
if (!arts_omp_in_parallel()) \
firstprivate(l_ws, l_iy_main_agenda, l_geo_pos_agenda)
      for( Index ilos=0; ilos<nlos; ilos++ )
        {
          // Skip remaining iterations if an error occurred
          if (failed) continue;
          
          Ppath ppath;
          iyb_calc_body( failed, fail_msg, iy_aux_array, l_ws,
                         ppath, iyb, diyb_dx,
                         mblock_index, atmosphere_dim, t_field, z_field,
                         vmr_field, nlte_field, cloudbox_on, stokes_dim, f_grid,
                         sensor_pos, sensor_los, transmitter_pos,
                         mblock_dlos_grid, iy_unit, 
                         l_iy_main_agenda, j_analytical_do, 
                         jacobian_quantities, jacobian_indices,
                         iy_aux_vars, ilos, nf );
          
          // Skip remaining iterations if an error occurred
          if (failed) continue;

          // Note that this code is found in two places inside the function
          Vector geo_pos;
          try
          {
              geo_pos_agendaExecute( l_ws, geo_pos, ppath, l_geo_pos_agenda );
              if( geo_pos.nelem() )
                {
                  if( geo_pos.nelem() != 5 )
                      throw runtime_error(
                       "Wrong size of *geo_pos* obtained from *geo_pos_agenda*.\n"
                       "The length of *geo_pos* must be zero or five." );

                  geo_pos_matrix(ilos,joker) = geo_pos;
                }
          }
          catch (const std::runtime_error &e)
          {
#pragma omp critical (iyb_calc_fail)
              { fail_msg = e.what(); failed = true; }
          }
        }
    }
  else
    {
      out3 << "  Not parallelizing los loop (" << nlos << " iterations, "
      << nf << " frequencies)\n";

      for( Index ilos=0; ilos<nlos; ilos++ )
        {
          // Skip remaining iterations if an error occurred
          if (failed) continue;

          Ppath ppath;
          iyb_calc_body( failed, fail_msg, iy_aux_array, l_ws, 
                         ppath, iyb, diyb_dx,
                         mblock_index, atmosphere_dim, t_field, z_field,
                         vmr_field, nlte_field, cloudbox_on, stokes_dim, f_grid,
                         sensor_pos, sensor_los, transmitter_pos,
                         mblock_dlos_grid, iy_unit, 
                         l_iy_main_agenda, j_analytical_do, 
                         jacobian_quantities, jacobian_indices,
                         iy_aux_vars, ilos, nf );

          // Skip remaining iterations if an error occurred
          if (failed) continue;

          // Note that this code is found in two places inside the function
          Vector geo_pos;
          try
          {
              geo_pos_agendaExecute( l_ws, geo_pos, ppath, l_geo_pos_agenda );
              if( geo_pos.nelem() )
                {
                  if( geo_pos.nelem() != 5 )
                      throw runtime_error(
                       "Wrong size of *geo_pos* obtained from *geo_pos_agenda*.\n"
                       "The length of *geo_pos* must be zero or five." );

                  geo_pos_matrix(ilos,joker) = geo_pos;
                }
          }
          catch (const std::runtime_error &e)
          {
#pragma omp critical (iyb_calc_fail)
              { fail_msg = e.what(); failed = true; }
          }
        }
    }

  if( failed )
    throw runtime_error("Run-time error in function: iyb_calc\n" + fail_msg);

  // Compile iyb_aux
  //
  const Index nq = iy_aux_array[0].nelem();
  iyb_aux.resize( nq );
  //
  for( Index q=0; q<nq; q++ )
    {
      iyb_aux[q].resize( niyb );
      //
      for( Index ilos=0; ilos<nlos; ilos++ )
        {
          const Index row0 = ilos * nf * stokes_dim;
          for( Index iv=0; iv<nf; iv++ )
            { 
              const Index row1 = row0 + iv*stokes_dim;
              const Index i1 = min( iv, iy_aux_array[ilos][q].nrows()-1 );
              for( Index is=0; is<stokes_dim; is++ )
                { 
                  Index i2 = min( is, iy_aux_array[ilos][q].ncols()-1 );
                  iyb_aux[q][row1+is] = iy_aux_array[ilos][q](i1,i2);
                }
            }
        }
    }
}



//! iy_transmission_mult
/*!
    Multiplicates iy_transmission with (vector) transmissions.

    That is, a multiplication of *iy_transmission* with another
    variable having same structure and holding transmission values.

    The "new path" is assumed to be further away from the sensor than 
    the propagtion path already included in iy_transmission. That is,
    the operation can be written as:
    
       Ttotal = Told * Tnew

    where Told is the transmission corresponding to *iy_transmission*
    and Tnew corresponds to *tau*.

    *iy_trans_new* is sized by the function.

    \param   iy_trans_total    Out: Updated version of *iy_transmission*
    \param   iy_trans_old      A variable matching *iy_transmission*.
    \param   iy_trans_new      A variable matching *iy_transmission*.

    \author Patrick Eriksson 
    \date   2009-10-06
*/
void iy_transmission_mult( 
       Tensor3&      iy_trans_total,
  ConstTensor3View   iy_trans_old,
  ConstTensor3View   iy_trans_new )
{
  const Index nf = iy_trans_old.npages();
  const Index ns = iy_trans_old.ncols();

  assert( ns == iy_trans_old.nrows() );
  assert( nf == iy_trans_new.npages() );
  assert( ns == iy_trans_new.nrows() );
  assert( ns == iy_trans_new.ncols() );

  iy_trans_total.resize( nf, ns, ns );

  for( Index iv=0; iv<nf; iv++ )
    {
      mult( iy_trans_total(iv,joker,joker), iy_trans_old(iv,joker,joker),
                                            iy_trans_new(iv,joker,joker) );
    } 
}



//! iy_transmission_mult
/*!
    Multiplicates iy_transmission with iy-variable.

    the operation can be written as:
    
       iy_new = T * iy_old

    where T is the transmission corresponding to *iy_transmission*
    and iy_old is a variable matching iy.

    *iy_new* is sized by the function.

    \param   iy_new        Out: Updated version of iy 
    \param   iy_trans      A variable matching *iy_transmission*.
    \param   iy_old        A variable matching *iy*.

    \author Patrick Eriksson 
    \date   2018-04-10
*/
void iy_transmission_mult( 
       Matrix&       iy_new,
  ConstTensor3View   iy_trans,
  ConstMatrixView    iy_old )
{
  const Index nf = iy_trans.npages();
  const Index ns = iy_trans.ncols();

  assert( ns == iy_trans.nrows() );
  assert( nf == iy_old.nrows() );
  assert( ns == iy_old.ncols() );

  iy_new.resize( nf, ns );

  for( Index iv=0; iv<nf; iv++ )
    { mult( iy_new(iv,joker), iy_trans(iv,joker,joker), iy_old(iv,joker) ); } 
}



//! mirror_los
/*!
    Determines the backward direction for a given line-of-sight.

    This function can be used to get the LOS to apply for extracting single
    scattering properties, if the propagation path LOS is given.

    A viewing direction of aa=0 is assumed for 1D. This corresponds to 
    positive za for 2D.

    \param   los_mirrored      Out: The line-of-sight for reversed direction.
    \param   los               A line-of-sight
    \param   atmosphere_dim    As the WSV.

    \author Patrick Eriksson 
    \date   2011-07-15
*/
void mirror_los(
        Vector&     los_mirrored,
  ConstVectorView   los, 
  const Index&      atmosphere_dim )
{
  los_mirrored.resize(2);
  //
  if( atmosphere_dim == 1 )
    { 
      los_mirrored[0] = 180 - los[0]; 
      los_mirrored[1] = 180; 
    }
  else if( atmosphere_dim == 2 )
    {
      los_mirrored[0] = 180 - fabs( los[0] ); 
      if( los[0] >= 0 )
        { los_mirrored[1] = 180; }
      else
        { los_mirrored[1] = 0; }
    }
  else if( atmosphere_dim == 3 )
    { 
      los_mirrored[0] = 180 - los[0]; 
      los_mirrored[1] = los[1] + 180; 
      if( los_mirrored[1] > 180 )
        { los_mirrored[1] -= 360; }
    }
}    



//! pos2true_latlon
/*!
    Determines the true alt and lon for an "ARTS position"

    The function disentangles if the geographical position shall be taken from
    lat_grid and lon_grid, or lat_true and lon_true.

    \param   lat              Out: True latitude.
    \param   lon              Out: True longitude.
    \param   atmosphere_dim   As the WSV.
    \param   lat_grid         As the WSV.
    \param   lat_true         As the WSV.
    \param   lon_true         As the WSV.
    \param   pos              A position, as defined for rt calculations.

    \author Patrick Eriksson 
    \date   2011-07-15
*/
void pos2true_latlon( 
          Numeric&     lat,
          Numeric&     lon,
    const Index&       atmosphere_dim,
    ConstVectorView    lat_grid,
    ConstVectorView    lat_true,
    ConstVectorView    lon_true,
    ConstVectorView    pos )
{
  assert( pos.nelem() == atmosphere_dim );

  if( atmosphere_dim == 1 )
    {
      assert( lat_true.nelem() == 1 );
      assert( lon_true.nelem() == 1 );
      //
      lat = lat_true[0];
      lon = lon_true[0];
    }

  else if( atmosphere_dim == 2 )
    {
      assert( lat_true.nelem() == lat_grid.nelem() );
      assert( lon_true.nelem() == lat_grid.nelem() );
      GridPos   gp;
      Vector    itw(2);
      gridpos( gp, lat_grid, pos[1] );
      interpweights( itw, gp );
      lat = interp( itw, lat_true, gp );
      lon = interp( itw, lon_true, gp );
    }

  else 
    {
      lat = pos[1];
      lon = pos[2];
    }
}


//! get_stepwise_clearsky_propmat
/*!
 *  Gets the clearsky propgation matrix and NLTE contributions
 * 
 *  \param K                Out: Level propagation matrix
 *  \param S                Out: NLTE source vector for level
 *  \param lte              Out: Index indicating if there is any NLTE source term
 *  \param dK_dx            Out: Unadopted propagation matrix derivatives of level
 *  \param dS_dx            Out: Unadopted NLTE source derivatives of level
...
 * 
 *  \author Richard Larsson 
 *  \adapted from non-stepwise function
 *  \date   2017-09-21
 */
void get_stepwise_clearsky_propmat(Workspace& ws,
                                   PropagationMatrix& K,
                                   StokesVector& S,
                                   Index& lte,
                                   ArrayOfPropagationMatrix& dK_dx,
                                   ArrayOfStokesVector& dS_dx,
                                   const Agenda& propmat_clearsky_agenda,
                                   const ArrayOfRetrievalQuantity& jacobian_quantities,
                                   ConstVectorView ppath_f_grid,
                                   ConstVectorView ppath_magnetic_field,
                                   ConstVectorView ppath_line_of_sight,
                                   ConstVectorView ppath_nlte_temperatures,
                                   ConstVectorView ppath_vmrs,
                                   const Numeric& ppath_temperature,
                                   const Numeric& ppath_pressure,
                                   const ArrayOfIndex& jacobian_species,
                                   const bool& jacobian_do)
{
  // All relevant quantities are extracted first
  const Index nq = jacobian_quantities.nelem();
  
  // Local variables inside Agenda
  ArrayOfPropagationMatrix propmat_clearsky, dpropmat_clearsky_dx;
  ArrayOfStokesVector nlte_source, dnlte_dx_source, nlte_dx_dsource_dx;
  
  // Perform the propagation matrix computations
  propmat_clearsky_agendaExecute( 
                ws, propmat_clearsky, nlte_source, 
                dpropmat_clearsky_dx, dnlte_dx_source, nlte_dx_dsource_dx,
                jacobian_quantities, 
                ppath_f_grid, ppath_magnetic_field, ppath_line_of_sight,
                ppath_pressure, ppath_temperature, 
                ppath_nlte_temperatures, ppath_vmrs, 
                propmat_clearsky_agenda);
  
  // We only now know how large the propagation matrix will be!
  const Index npmat = propmat_clearsky.nelem();
  
  // If therea re no NLTE elements, then set the LTE flag
  lte = nlte_source.nelem()?0:1; 
  
  // Sum the propagation matrix
  K = propmat_clearsky[0];
  for(Index i = 1; i < npmat; i++)
    K += propmat_clearsky[i];
  
  // Set NLTE source term if applicable
  if(not lte)
  {
    // Add all source terms up
    S = nlte_source[0];
    for(Index i = 1; i < npmat; i++)
      S += nlte_source[i];
  }
  else
  {
    S.SetZero();
  }
  
  // Set the partial derivatives
  if(jacobian_do)
  {
    for(Index i = 0; i < nq; i++)
    {
      if(jacobian_quantities[i].MainTag() == SCATSPECIES_MAINTAG)
      {
        dK_dx[i].SetZero();
        dS_dx[i].SetZero();
      }
      else if(jacobian_quantities[i].SubSubtag() == PROPMAT_SUBSUBTAG) 
      {
        // Find position of index in ppd
        const Index j = equivlent_propmattype_index(jacobian_quantities, i);
        
        dK_dx[i] = dpropmat_clearsky_dx[j];
        if(lte)
        {
          dS_dx[i].SetZero();
        }
        else
        {
          dS_dx[i] = dnlte_dx_source[j];
          dS_dx[i] += nlte_dx_dsource_dx[j];
          
          // TEST:  Old routine applied unit conversion on only the last
          // part of the equation. Was this correct or wrong?  If correct,
          // this is an issue.  Otherwise, the old version was incorrect.  
          // Have to setup perturbation test-case to study which is most 
          // reasonable...
        }
      }
      else if(jacobian_species[i] > -1)  // Did not compute values in Agenda
      {
        dK_dx[i] = propmat_clearsky[jacobian_species[i]];
        
        // We cannot know the NLTE jacobian if this method was used
        // because that information is thrown away. It is still faster 
        // to retain this method since it requires less computations 
        // when we do not need NLTE, which is most of the time...
        if(not lte)
        {
          ostringstream os;
          
          os << "We do not yet support species" <<
                " tag and NLTE Jacobians.\n";
          throw std::runtime_error(os.str());
        }
        dS_dx[i].SetZero();
      }
    }
  }
}


//! adapt_stepwise_partial_derivatives
/*!
 *  Adapts clearsky partial derivatives for the following fields:
 * 
 *   Wind
 *   VMR
 * 
 *  Adaptation means changing unit by user input
 * 
 *  \param dK_dx            In/Out: In is unadopted out is adopted propagation matrix derivatives
 *  \param dS_dx            In/Out: In is unadopted extra source and out is adopted extra source derivatives
...
 * 
 *  \author Richard Larsson 
 *  \adapted from non-stepwise function
 *  \date   2017-09-21
 */
void adapt_stepwise_partial_derivatives(ArrayOfPropagationMatrix& dK_dx,
                                        ArrayOfStokesVector& dS_dx,
                                        const ArrayOfRetrievalQuantity& jacobian_quantities,
                                        ConstVectorView ppath_f_grid,
                                        ConstVectorView ppath_line_of_sight,
                                        ConstVectorView ppath_vmrs,
                                        const Numeric& ppath_temperature,
                                        const Numeric& ppath_pressure,
                                        const ArrayOfIndex& jacobian_species,
                                        const ArrayOfIndex& jacobian_wind,
                                        const Index& lte,
                                        const Index& atmosphere_dim,
                                        const bool& jacobian_do)
{
  if(not jacobian_do)
    return;
  
  // All relevant quantities are extracted first
  const Index nq = jacobian_quantities.nelem();
  
  // Computational temporary vector
  Vector a;
  
  for(Index i = 0; i < nq; i++)
  {
    if(not jacobian_quantities[i].Analytical())
      continue;
    
    Index component;
    
    switch(jacobian_wind[i])
    {
      case Index(JacobianType::AbsWind):
        component = 0;
        break;
      case Index(JacobianType::WindFieldU):
        component = 1;
        break;
      case Index(JacobianType::WindFieldV):
        component = 2;
        break;
      case Index(JacobianType::WindFieldW):
        component = 3;
        break;
      default:
        component = -1;  // This could be frequency derivative...
    }
    if(component not_eq -1)
    {
      get_stepwise_f_partials(a, component,
                              ppath_line_of_sight, 
                              ppath_f_grid,  atmosphere_dim);
      
      // Apply conversion to K-matrix partial derivative
      dK_dx[i] *= a;
      
      // Apply conversion to source vector partial derivative
      if(not lte)
        dS_dx[i] *= a;
    }
    else if(jacobian_species[i] > -1)
    {
      const bool from_propmat = jacobian_quantities[i].SubSubtag() == PROPMAT_SUBSUBTAG;
      const Index& isp = jacobian_species[i];
      
      // Computational factor
      Numeric factor;
        
      // Scaling factors to handle retrieval unit
      if(not from_propmat)
      {
        //vmrunitscf(factor, jacobian_quantities[i].Mode(), 
        vmrunitscf(factor, "vmr", 
                   ppath_vmrs[isp], ppath_pressure, 
                   ppath_temperature);
      }
      else 
      {
        //dxdvmrscf(factor, jacobian_quantities[i].Mode(), 
        dxdvmrscf(factor, "vmr", 
                  ppath_vmrs[isp], ppath_pressure, 
                  ppath_temperature);
      }
      // Apply conversion to K-matrix partial derivative
      dK_dx[i] *= factor;
      
      // Apply conversion to source vector partial derivative
      if(not lte)
        dS_dx[i] *= factor;
    }
    else 
    {
      // All other partial derivatives should already be adapted...
    }
  }
}


Numeric guesswork_HSE_derivative(Numeric h, Numeric r, Numeric T) { return 2 * h * std::abs(h/r) / T; /*std::abs to keep sign since one level should be negative and the other positive*/ }


//! get_stepwise_transmission_matrix
/*!
 *  Computes layer transmission matrix and cumulative transmission
 * 
 *  \param cumulative_transmission Out: cumulation of transmission for jacobian computations
 *  \param T                       Out: Layer transmission
 *  \param dT_close_dx             Out: Layer transmission derivative due to closest level
 *  \param dT_far_dx               Out: Layer transmission derivative due to furthest level
...
 * 
 *  \author Richard Larsson 
 *  \adapted from non-stepwise function
 *  \date   2017-09-21
 */
void get_stepwise_transmission_matrix(Tensor3View cumulative_transmission,
                                      Tensor3View T,
                                      Tensor4View dT_close_dx,
                                      Tensor4View dT_far_dx,
                                      ConstTensor3View cumulative_transmission_close,
                                      const PropagationMatrix& K_close,
                                      const PropagationMatrix& K_far,
                                      const ArrayOfPropagationMatrix& dK_close_dx,
                                      const ArrayOfPropagationMatrix& dK_far_dx,
                                      const Numeric& ppath_distance,
                                      const bool& first_level,
                                      const Numeric& dppath_distance_dT_HSE_guesswork_close,
                                      const Numeric& dppath_distance_dT_HSE_guesswork_far,
                                      const Index& temperature_derivative_position_if_hse_is_active)
{
  // Frequency counter
  const Index nf = K_close.NumberOfFrequencies();
  const Index stokes_dim = T.ncols();
  
  if(first_level) {
    if(stokes_dim>1)
      for(Index iv = 0; iv < nf; iv++)
        id_mat(cumulative_transmission(iv, joker, joker));
    else 
      cumulative_transmission = 1;
  }
  else {
    // Compute the transmission of the layer between close and far
    if(not dK_close_dx.nelem())
      compute_transmission_matrix(T, ppath_distance, K_close, K_far);
    else
      compute_transmission_matrix_and_derivative(T, 
                                                dT_close_dx, dT_far_dx, 
                                                ppath_distance, 
                                                K_close, K_far, 
                                                dK_close_dx, dK_far_dx,
                                                dppath_distance_dT_HSE_guesswork_close,
                                                dppath_distance_dT_HSE_guesswork_far,
                                                temperature_derivative_position_if_hse_is_active);
    
    // Cumulate transmission
    if(stokes_dim>1)
      for(Index iv = 0; iv < nf; iv++)
        mult(cumulative_transmission(iv, joker, joker), 
            cumulative_transmission_close(iv, joker, joker),
            T(iv, joker, joker));
    else
    {
      cumulative_transmission = cumulative_transmission_close;
      cumulative_transmission  *= T;
    }
  }
}


void get_stepwise_blackbody_radiation(VectorView B,
                                      VectorView dB_dT,
                                      ConstVectorView ppath_f_grid,
                                      const Numeric& ppath_temperature,
                                      const bool& do_temperature_derivative)
{
  const Index nf = ppath_f_grid.nelem();
  
  for(Index i = 0; i < nf; i++)
    B[i] = planck(ppath_f_grid[i], ppath_temperature);
  
  if(do_temperature_derivative)
    for(Index i = 0; i < nf; i++)
      dB_dT[i] = dplanck_dt(ppath_f_grid[i], ppath_temperature);
}


//! get_stepwise_effective_source
/*!
 *  Computes
 * 
 *  J = K^-1 (a B + S)
 * 
 *  and
 * 
 *  dJ = - K^-1 dK/dx K^-1 (a B + S) + K^-1 (da B + a dB + dS)
 * 
 *  Assumes zeroes for the a and K if nothing is happening but checks all other variables
 * 
 *  \param J                   Out: Source term, frequency times stokes dimension
 *  \param dJ_dx               Out: Source term derivative, quantities times frequency times stokes dimension
 *  \param K                   In: Propagation matrix of level --- all contributions
 *  \param a                   In: Absorption vector of level --- all contributions
 *  \param S                   In: Source terms other than absorption times planck --- all contributions
 *  \param dK_dx               In: Propagation matrix derivatives of level --- all contributions
 *  \param da_dx               In: Absorption vector derivatives of level --- all contributions
 *  \param dS_dx               In: Source terms derivatives other than absorption times planck --- all contributions
 *  \param B                   In: Planck function in Stokes vector form
 *  \param dB_dT               In: Planck function derivative wrt temperatures in Stokes vector form
 *  \param jacobian_quantities In: As wsv
 * 
 *  \author Richard Larsson 
 *  \adapted from non-stepwise function
 *  \date   2017-09-21
 */
void get_stepwise_effective_source(MatrixView J,
                                   Tensor3View dJ_dx,
                                   const PropagationMatrix& K,
                                   const StokesVector& a,
                                   const StokesVector& S,
                                   const ArrayOfPropagationMatrix& dK_dx,
                                   const ArrayOfStokesVector& da_dx,
                                   const ArrayOfStokesVector& dS_dx,
                                   ConstVectorView B,
                                   ConstVectorView dB_dT,
                                   const ArrayOfRetrievalQuantity& jacobian_quantities,
                                   const bool& jacobian_do)
{
  const Index nq = jacobian_quantities.nelem();
  const Index nf = K.NumberOfFrequencies();
  const Index ns = K.StokesDimensions();
  
  dJ_dx = 0;
  for(Index i1 = 0; i1 < nf; i1++)
  {
    Matrix invK(ns, ns);
    Vector j(ns);
    
    if(K.IsRotational(i1)) {
      dJ_dx(joker, i1, joker) = 0;
      J(i1, joker) = 0;
      continue;
    }
    
    K.MatrixInverseAtPosition(invK, i1);
    
    // Set a B to j
    j = a.VectorAtPosition(i1);;
    j *= B[i1];
    
    // Add S to j
    if(not S.IsEmpty())
      j += S.VectorAtPosition(i1);;
    
    // Compute J = K^-1 (a B + S)
    mult(J(i1, joker), invK, j);
    
    // Compute dJ = K^-1((da B + a dB + dS) - dK K^-1(a B + S))
    if( jacobian_do )
      //FOR_ANALYTICAL_JACOBIANS_DO
      //(
      for(Index iq = 0; iq < nq; iq++)
      {
        if( jacobian_quantities[iq].Analytical() )
        {
          Matrix dk(ns, ns), tmp_matrix(ns, ns);
          Vector dj(ns, 0), tmp(ns);
      
          // Control parameters for special jacobians that are computed elsewhere
          //const bool has_dk = (not dK_dx[iq].IsEmpty());   // currently always
          //const bool has_ds = (not dS_dx[iq].IsEmpty());   // evaluate as true
          const bool has_dt = (jacobian_quantities[iq].MainTag() == TEMPERATURE_MAINTAG);
      
          // Sets the -K^-1 dK/dx K^-1 (a B + S) term
          //if(has_dk)
          //{
            dK_dx[iq].MatrixAtPosition(dk, i1);
            mult(tmp, dk, J(i1, joker));
        
            dj = da_dx[iq].VectorAtPosition(i1);
            dj *= B[i1];
        
            dj -= tmp;
        
            // Adds a dB to dj
            if(has_dt)
            {
              tmp = a.VectorAtPosition(i1);
              tmp *= dB_dT[i1];
              dj += tmp;
            }
        
            // Adds dS to dj
            //if(has_ds)
              dj += dS_dx[iq].VectorAtPosition(i1);
        
            mult(dJ_dx(iq, i1, joker), invK, dj);
          //}
        }
        // don't need that anymore now that we zero dJ_dx at the very beginning
        //else
        //{
        //  dJ_dx(iq, i1, joker) = 0;
        //}
      }
      //)
  }
  
  if(nq)
    dJ_dx *= 0.5;
}


//! sum_stepwise_scalar_tau_and_extmat_case
/*!
 *  Sums the scala tau for iy_aux...
 * 
 *  \param scalar_tau                  Out: as in iy_aux
...
 * 
 *  \author Richard Larsson 
 *  \adapted from non-stepwise function
 *  \date   2017-09-21
 */
void sum_stepwise_scalar_tau_and_extmat_case(VectorView scalar_tau,
                                             ArrayOfIndex& extmat_case,
                                             const PropagationMatrix& upper_level,
                                             const PropagationMatrix& lower_level,
                                             const Numeric& distance)
{
  Index nf = scalar_tau.nelem();
  
  for(Index i = 0; i < nf; i++)
  {
    scalar_tau[i] += 0.5 * (upper_level.Kjj()[i] + lower_level.Kjj()[i]) * distance;
    
    if(upper_level.StokesDimensions() > 1 or lower_level.StokesDimensions() > 1)
      extmat_case[i] = 3;
    else
      extmat_case[i] = 1;
  }
}


//! get_stepwise_frequency_grid
/*!
 *  Inverse of get_stepwise_f_partials
 * 
 *  Computes practical frequency grid due to wind for propmat_clearsky_agenda
 * 
...
 * 
 *  \author Richard Larsson 
 *  \adapted from non-stepwise function
 *  \date   2017-09-21
 */
void get_stepwise_frequency_grid(VectorView ppath_f_grid,
                                 ConstVectorView f_grid,
                                 ConstVectorView ppath_wind,
                                 ConstVectorView ppath_line_of_sight,
                                 const Numeric& rte_alonglos_v,
                                 const Index& atmosphere_dim)
{
  Numeric v_doppler = rte_alonglos_v;
  
  if(ppath_wind[0] not_eq 0 or ppath_wind[1] not_eq 0 or ppath_wind[2] not_eq 0)
    v_doppler += dotprod_with_los(ppath_line_of_sight, 
                                  ppath_wind[0], ppath_wind[1], ppath_wind[2],
                                  atmosphere_dim);
  ppath_f_grid = f_grid;
  
  if(v_doppler not_eq 0)
   ppath_f_grid *= 1 - v_doppler / SPEED_OF_LIGHT;
}


//! get_stepwise_f_partials
/*!
 *  Computes the ratio that a partial derivative with regards to frequency
 *  relates to the wind of come component
 * 
 *  \param   f_partial             Out: The frequency vector to multiply each frequency 
 *                                      grid point in clearsky propmat derivatives
 *  \param   component             In: The wind component
 *  \param   line_of_sight         In: The line of sight vector as in workspace rtp_los
 *  \param   f_grid                In: The computational frequency grid
 *  \param   atmosphere_dim        In: atmospheric diemension as wsv
 ...
 * 
 *  \author Richard Larsson 
 *  \adapted from non-stepwise function
 *  \date   2017-09-21
 */
void get_stepwise_f_partials(Vector& f_partials,
                             const Index& component,
                             ConstVectorView& line_of_sight,
                             ConstVectorView f_grid, 
                             const Index& atmosphere_dim)
{
  // component 0 means total speed
  // component 1 means u speed
  // component 2 means v speed
  // component 3 means w speed
  
  // Sizes
  const Index   nf = f_grid.nelem();
  
  // Doppler relevant velocity
  //
  // initialize
  Numeric dv_doppler_dx = 0.0;
    
  switch( component )
  {
    case 0:// this is total and is already initialized to avoid compiler warnings
        dv_doppler_dx = 1.0;
        break;
    case 1:// this is the u-component
        dv_doppler_dx = (dotprod_with_los(line_of_sight, 1, 0, 0, atmosphere_dim));
        break;
    case 2:// this is v-component
        dv_doppler_dx = (dotprod_with_los(line_of_sight, 0, 1, 0, atmosphere_dim));
        break;
    case 3:// this is w-component
        dv_doppler_dx = (dotprod_with_los(line_of_sight, 0, 0, 1, atmosphere_dim));
        break;
    default:
        throw std::runtime_error("This being seen means that there is a development bug in interactions with get_ppath_df_dW.\n");
        break;
  }
    
  // Determine frequency grid
  if( dv_doppler_dx == 0.0 )
  {
    f_partials.resize(nf);
    f_partials = 0.0;
  }
  else
  { 
    f_partials = f_grid;
    f_partials *= - dv_doppler_dx / SPEED_OF_LIGHT;
  }
}


//! get_stepwise_scattersky_propmat
/*!
 *  Computes the contribution by scattering elements towards the absorption 
 *  and emission from a level.
 * 
 *  \param   ap                    Out: The scattering absorption term
 *  \param   Kp                    Out: The scattering propagation matrix term
 *  \param   dap_dx                Out: The scattering absorption term deriative
 *  \param   dKp_dx                Out: The scattering propagation matrix term deriative
 ...
 * 
 *  \author Jana Mendrok, Richard Larsson 
 *  \adapted from non-stepwise function
 *  \date   2017-09-21
 */
void get_stepwise_scattersky_propmat(StokesVector& ap,
                                     PropagationMatrix& Kp,
                                     ArrayOfStokesVector& dap_dx,
                                     ArrayOfPropagationMatrix& dKp_dx,
                                     const ArrayOfRetrievalQuantity& jacobian_quantities,
                                     ConstMatrixView ppath_1p_pnd,       // the ppath_pnd at this ppath point
                                     const ArrayOfMatrix& ppath_dpnd_dx, // the full ppath_dpnd_dx, ie all ppath points
                                     const Index ppath_1p_id,
                                     const ArrayOfArrayOfSingleScatteringData& scat_data,
                                     ConstVectorView ppath_line_of_sight,
                                     ConstVectorView ppath_temperature,
                                     const Index& atmosphere_dim,
                                     const bool& jacobian_do)
{
  const Index nf = Kp.NumberOfFrequencies(),
              stokes_dim = Kp.StokesDimensions();

  //StokesVector da_aux(nf, stokes_dim);
  //PropagationMatrix dK_aux(nf, stokes_dim);

  ArrayOfArrayOfSingleScatteringData scat_data_mono;

  // Direction of outgoing scattered radiation (which is reversed to
  // LOS). Only used for extracting scattering properties.
  Vector dir;
  mirror_los(dir, ppath_line_of_sight, atmosphere_dim);
  Matrix dir_array(1,2,0.);
  dir_array(0,joker) = dir;

  ArrayOfArrayOfTensor5 ext_mat_Nse;
  ArrayOfArrayOfTensor4 abs_vec_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;
  ArrayOfTensor5 ext_mat_ssbulk;
  ArrayOfTensor4 abs_vec_ssbulk;
  ArrayOfIndex ptype_ssbulk;
  Tensor5 ext_mat_bulk;
  Tensor4 abs_vec_bulk;
  Index ptype_bulk;

  // get per-scat-elem data here. and fold with pnd.
  // keep per-scat-elem data to fold with dpnd_dx further down in
  // analyt-jac-loop.
  opt_prop_NScatElems( ext_mat_Nse, abs_vec_Nse, ptypes_Nse, t_ok,
                       scat_data, stokes_dim, ppath_temperature, dir_array, -1 );

  opt_prop_ScatSpecBulk( ext_mat_ssbulk, abs_vec_ssbulk, ptype_ssbulk,
                         ext_mat_Nse, abs_vec_Nse, ptypes_Nse,
                         ppath_1p_pnd, t_ok );
  opt_prop_Bulk( ext_mat_bulk, abs_vec_bulk, ptype_bulk,
                 ext_mat_ssbulk, abs_vec_ssbulk, ptype_ssbulk );

  const Index nf_ssd = abs_vec_bulk.nbooks(); // number of freqs in extracted
                                              // optprops. if 1, we need to
                                              // duplicate the ext/abs output.

  for( Index iv = 0; iv < nf; iv++ )
  {
    if( nf_ssd>1 )
    {
      ap.SetAtPosition(abs_vec_bulk(iv,0,0,joker), iv);
      Kp.SetAtPosition(ext_mat_bulk(iv,0,0,joker,joker), iv);
    }
    else
    {
      ap.SetAtPosition(abs_vec_bulk(0,0,0,joker), iv);
      Kp.SetAtPosition(ext_mat_bulk(0,0,0,joker,joker), iv);
    }
  }

  if( jacobian_do )
    FOR_ANALYTICAL_JACOBIANS_DO
    (
      if( ppath_dpnd_dx[iq].empty() )
      {
        dap_dx[iq].SetZero();
        dKp_dx[iq].SetZero();
      }
      else
      {
        // check, whether we have any non-zero ppath_dpnd_dx in this
        // pnd-affecting x? might speed up things a little bit.
        opt_prop_ScatSpecBulk( ext_mat_ssbulk, abs_vec_ssbulk, ptype_ssbulk,
                               ext_mat_Nse, abs_vec_Nse, ptypes_Nse,
                               ppath_dpnd_dx[iq](joker,Range(ppath_1p_id,1)),
                               t_ok );
        opt_prop_Bulk( ext_mat_bulk, abs_vec_bulk, ptype_bulk,
                       ext_mat_ssbulk, abs_vec_ssbulk, ptype_ssbulk );
        for( Index iv = 0; iv < nf; iv++ )
        {
          if( nf_ssd>1 )
          {
            dap_dx[iq].SetAtPosition(abs_vec_bulk(iv,0,0,joker), iv);
            dKp_dx[iq].SetAtPosition(ext_mat_bulk(iv,0,0,joker,joker), iv);
          }
          else
          {
            dap_dx[iq].SetAtPosition(abs_vec_bulk(0,0,0,joker), iv);
            dKp_dx[iq].SetAtPosition(ext_mat_bulk(0,0,0,joker,joker), iv);
          }
        }
      }
    )
}



//! get_stepwise_scattersky_source
/*!
 *  Calculates the stepwise scattering source terms.
 *  Uses new, unified phase matrix extraction scheme.
 * 
 *  \param   Sp                    Out: The scattering source term
 *  \param   dSp_dx                Out: The derivative of the scattering source term
 ...
 * 
 *  \author Jana Mendrok 
 *  \adapted from non-stepwise function
 *  \date   2018-03-29
 */
void get_stepwise_scattersky_source(StokesVector& Sp,
                                    ArrayOfStokesVector& dSp_dx,
                                    const ArrayOfRetrievalQuantity& jacobian_quantities,
                                    ConstVectorView ppath_1p_pnd,       // the ppath_pnd at this ppath point
                                    const ArrayOfMatrix& ppath_dpnd_dx, // the full ppath_dpnd_dx, ie all ppath points
                                    const Index ppath_1p_id,
                                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                                    ConstTensor7View doit_i_field,
                                    ConstVectorView scat_za_grid,
                                    ConstVectorView scat_aa_grid,
                                    ConstMatrixView ppath_line_of_sight,
                                    const GridPos& ppath_pressure,
                                    const Vector& temperature,
                                    const Index& atmosphere_dim,
                                    const bool& jacobian_do,
                                    const Index& t_interp_order)
{
  if( atmosphere_dim != 1 )
    throw runtime_error( "This function handles so far only 1D atmospheres." );

  const Index nf = Sp.NumberOfFrequencies();
  const Index stokes_dim = Sp.StokesDimensions();
  const Index ne = ppath_1p_pnd.nelem();
  assert( TotalNumberOfElements(scat_data) == ne );
  const Index nza = scat_za_grid.nelem();
  const Index naa = scat_aa_grid.nelem();
  const Index nq = jacobian_do ? jacobian_quantities.nelem() : 0;

  // interpolate incident field to this ppath point (no need to do this
  // separately per scatelem)
  GridPos gp_p;
  gridpos_copy( gp_p, ppath_pressure );
  Vector itw_p(2);
  interpweights( itw_p, gp_p );
  Tensor3 inc_field(nf, nza, stokes_dim, 0.);
  for( Index iv = 0; iv < nf; iv++ )
    {  
      for( Index iza = 0; iza < nza; iza++ )
        {
          for( Index i = 0; i < stokes_dim; i++ )
            {
              inc_field(iv, iza, i) =
                interp( itw_p, doit_i_field(iv,joker,0,0,iza,0,i), gp_p );
            }
        }
    }
  
  // create matrix of incident directions (flat representation of the
  // scat_za_grid * scat_aa_grid matrix)
  Matrix idir(nza*naa,2);
  Index ia=0;
  for( Index iza=0; iza<nza; iza++ )
    {
      for( Index iaa=0; iaa<naa; iaa++ )
        {
          idir(ia,0) = scat_za_grid[iza];
          idir(ia,1) = scat_aa_grid[iaa];
          ia++;
        }
    }
  
  // setting prop (aka scattered) direction
  Matrix pdir(1,2);
  //if( ppath_line_of_sight.ncols()==2 )
  //  pdir(0,joker) = ppath_line_of_sight;
  //else // 1D case only (currently the only handled case). azimuth not defined.
  {
    pdir(0,0) = ppath_line_of_sight(0,0);
    pdir(0,1) = 0.;
  }

  // some further variables needed for pha_mat extraction
  Index nf_ssd = scat_data[0][0].pha_mat_data.nlibraries();
  Index duplicate_freqs = ((nf==nf_ssd)?0:1);
  Tensor6 pha_mat_1se(nf_ssd,1,1,nza*naa,stokes_dim,stokes_dim);
  Vector t_ok(1);
  Index ptype;
  Tensor3 scat_source_1se(ne, nf, stokes_dim, 0.);

  Index ise_flat = 0;
  for( Index i_ss = 0; i_ss<scat_data.nelem(); i_ss++ )
    {
      for( Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++ )
        {
          // determine whether we have some valid pnd for this
          // scatelem (in pnd or dpnd)
          Index val_pnd = 0;
          if( ppath_1p_pnd[ise_flat] != 0 )
            { val_pnd = 1; }
          else if( jacobian_do )
            {
              for( Index iq=0; (!val_pnd) && (iq<nq); iq++ )
                {
                  if( jacobian_quantities[iq].Analytical() &&
                      !ppath_dpnd_dx[iq].empty() &&
                      ppath_dpnd_dx[iq](ise_flat,ppath_1p_id) != 0 )
                        { val_pnd = 1; }
                }
            }

          if( val_pnd )
            {
              pha_mat_1ScatElem( pha_mat_1se, ptype, t_ok,
                                 scat_data[i_ss][i_se],
                                 temperature, pdir, idir,
                                 0, t_interp_order );
              if( t_ok[0] == 0 )
                {
                  ostringstream os;
                  os << "Interpolation error for (flat-array) scattering "
                     << "element #" << ise_flat << "\n"
                     << "at location/temperature point #" << ppath_1p_id << "\n";
                  throw runtime_error( os.str() );
                }
              
              Index this_iv = 0;
              for( Index iv = 0; iv < nf; iv++ )
                {
                  if( !duplicate_freqs )
                    { this_iv = iv; }
                  Tensor3 product_fields(nza, naa, stokes_dim, 0.);

                  ia = 0;
                  for( Index iza = 0; iza < nza; iza++ )
                    {
                      for( Index iaa = 0; iaa < naa; iaa++ )
                        {
                          for ( Index i = 0; i < stokes_dim; i++)
                            {
                              for ( Index j = 0; j < stokes_dim; j++ )
                                {
                                  product_fields(iza, iaa, i) +=
                                    pha_mat_1se(this_iv, 0, 0, ia, i, j) *
                                    inc_field(iv, iza, j);
                                }
                            }
                          ia++;
                        }
                    }
          
                  for ( Index i = 0; i < stokes_dim; i++ )
                    {
                      scat_source_1se( ise_flat, iv, i) = AngIntegrate_trapezoid( 
                        product_fields(joker, joker, i), scat_za_grid,
                                                         scat_aa_grid );
                    }
                }  // for iv
            } // if val_pnd

          ise_flat++;

        } // for i_se
    } // for i_ss

  for( Index iv = 0; iv < nf; iv++ )
    {
      Vector scat_source(stokes_dim, 0.);
      for( ise_flat = 0; ise_flat < ne; ise_flat++ )
        {
          for( Index i=0; i<stokes_dim; i ++ )
            {
              scat_source[i] += scat_source_1se(ise_flat,iv,i) * ppath_1p_pnd[ise_flat];
            }
        }

      Sp.SetAtPosition(scat_source, iv);

      if( jacobian_do )
        {
          FOR_ANALYTICAL_JACOBIANS_DO(
            if( ppath_dpnd_dx[iq].empty() )
              { dSp_dx[iq].SetZero(); }
            else
              {
                scat_source = 0.;
                for( ise_flat = 0; ise_flat < ne; ise_flat++ )
                  {
                    for( Index i=0; i<stokes_dim; i ++ )
                      {
                        scat_source[i] += scat_source_1se(ise_flat,iv,i) *
                          ppath_dpnd_dx[iq](ise_flat,ppath_1p_id);
                        dSp_dx[iq].SetAtPosition(scat_source,iv);
                      }
                  }
              }
           )
        }
    } // for iv
}



//! rtmethods_jacobian_init
/*!
    This function fixes the initial steps around Jacobian calculations, to be
    done inside radiative transfer WSMs.

    See iyEmissonStandard for usage example.

    \author Patrick Eriksson 
    \date   2017-11-20
*/
void rtmethods_jacobian_init(
         ArrayOfIndex&               jac_species_i,
         ArrayOfIndex&               jac_scat_i,
         ArrayOfIndex&               jac_is_t,
         ArrayOfIndex&               jac_wind_i,
         ArrayOfIndex&               jac_mag_i,
         ArrayOfIndex&               jac_other,
         ArrayOfTensor3&             diy_dx,
         ArrayOfTensor3&             diy_dpath,         
   const Index&                      ns,
   const Index&                      nf,
   const Index&                      np,
   const Index&                      nq,
   const ArrayOfArrayOfSpeciesTag&   abs_species,
   const Index&                      cloudbox_on,
   const ArrayOfString&              scat_species,         
   const ArrayOfTensor4&             dpnd_field_dx,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,   
   const Index&                      iy_agenda_call1,
   const bool                        is_active )
{
  const Index nn = is_active ? nf*np : nf;
  
  FOR_ANALYTICAL_JACOBIANS_DO( 
    diy_dpath[iq].resize( np, nn, ns ); 
    diy_dpath[iq] = 0.0;
  )

  get_pointers_for_analytical_jacobians( jac_species_i, jac_scat_i, jac_is_t, 
                                         jac_wind_i, jac_mag_i, 
                                         jacobian_quantities, abs_species,
                                         cloudbox_on, scat_species );
  
  FOR_ANALYTICAL_JACOBIANS_DO( 
    jac_other[iq] = (jacobian_quantities[iq].PropMatType() == JacPropMatType::NotPropagationMatrixType) ? Index(JacobianType::Other) : Index(JacobianType::None); 
      
    if( jac_scat_i[iq]+1 )
      {
        if( dpnd_field_dx[iq].empty() )
          throw runtime_error( "*dpnd_field_dx* not allowed to be empty for "
                               "scattering Jacobian species." );
      }
    // FIXME: should we indeed check for that? remove if it causes issues.
    else
      {
        if( !dpnd_field_dx[iq].empty() )
          throw runtime_error( "*dpnd_field_dx* must be empty for "
                               "non-scattering Jacobian species." );
      }
  )

  if( iy_agenda_call1 )
    {
      diy_dx.resize(nq); 
      //
      bool any_affine;
      ArrayOfArrayOfIndex jacobian_indices;
      jac_ranges_indices( jacobian_indices, any_affine,
                          jacobian_quantities, true );
      //
      FOR_ANALYTICAL_JACOBIANS_DO2( 
        diy_dx[iq].resize( jacobian_indices[iq][1]-jacobian_indices[iq][0]+1,
                           nn, ns );
        diy_dx[iq] = 0.0;
      )
    }
}



//! rtmethods_jacobian_finalisation
/*!
    This function fixes the last steps to made on the Jacobian in some
    radiative transfer WSMs. The method applies iy_transmission, maps from
    ppath to the retrieval grids and applies non-standard Jacobian units.

    See iyEmissonStandard for usage example.

    \author Patrick Eriksson 
    \date   2017-11-19
*/
void rtmethods_jacobian_finalisation(
         Workspace&                  ws,
         ArrayOfTensor3&             diy_dx,
         ArrayOfTensor3&             diy_dpath,  
   const Index&                      ns,
   const Index&                      nf,
   const Index&                      np,
   const Index&                      atmosphere_dim,
   const Ppath&                      ppath,
   const Vector&                     ppvar_p,
   const Vector&                     ppvar_t,
   const Matrix&                     ppvar_vmr,
   const Index&                      iy_agenda_call1,         
   const Tensor3&                    iy_transmission,
   const Agenda&                     water_p_eq_agenda,   
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const ArrayOfIndex                jac_species_i,
   const ArrayOfIndex                jac_is_t)
{
  // Weight with iy_transmission
  if( !iy_agenda_call1 )
    {
      Matrix X, Y; 
      //
      FOR_ANALYTICAL_JACOBIANS_DO
        ( 
          Y.resize(ns,diy_dpath[iq].npages());
          for( Index iv=0; iv<nf; iv++ )
            { 
              X = transpose( diy_dpath[iq](joker,iv,joker) );
              mult( Y, iy_transmission(iv,joker,joker), X );
              diy_dpath[iq](joker,iv,joker) = transpose( Y );
            }
        )
    }

  
  // Handle abs species retrieval units, both internally and impact on T-jacobian
  //
  Tensor3 water_p_eq(0,0,0);
  //
  // Conversion for abs species itself
  for( Index iq=0; iq<jacobian_quantities.nelem(); iq++ )
    {
      // Let x be VMR, and z the selected retrieval unit.
      // We have then that diy/dz = diy/dx * dx/dz
      //
      if( jacobian_quantities[iq].Analytical()  &&  jac_species_i[iq] >= 0 )
        {
          if( jacobian_quantities[iq].Mode() == "vmr" )
            {}
          
          else if( jacobian_quantities[iq].Mode() == "rel" )
            {
              // Here x = vmr*z
              for( Index ip=0; ip<np; ip++ )
                {
                  diy_dpath[iq](ip,joker,joker) *=
                    ppvar_vmr(jac_species_i[iq],ip);
                }
            }
          
          else if( jacobian_quantities[iq].Mode() == "nd" )
            {
              // Here x = z/nd_tot
              for( Index ip=0; ip<np; ip++ )
                {
                  diy_dpath[iq](ip,joker,joker) /=
                    number_density( ppvar_p[ip], ppvar_t[ip] );
                }
            }

          else if( jacobian_quantities[iq].Mode() == "rh" )
            {
              // Here x = (p_sat/p) * z
              Tensor3 t_data(ppvar_t.nelem(),1,1); t_data(joker,0,0) = ppvar_t;
              water_p_eq_agendaExecute( ws, water_p_eq, t_data,
                                        water_p_eq_agenda);
              for( Index ip=0; ip<np; ip++ )
                { diy_dpath[iq](ip,joker,joker) *=
                    water_p_eq(ip,0,0) / ppvar_p[ip]; }
            }

          else if( jacobian_quantities[iq].Mode() == "q" )
            {
              // Here we use the approximation of x = z/0.622              
              diy_dpath[iq](joker,joker,joker) /= 0.622;
            }
          
          else
            { assert(0); }
        }
    }
  
  // Correction of temperature Jacobian
  for( Index iq=0; iq<jacobian_quantities.nelem(); iq++ )
    {
      // Let a be unit for abs species, and iy = f(T,a(T))
      // We have then that diy/dT = df/dT + df/da*da/dT
      // diy_dpath holds already df/dT. Remains is to add
      // df/da*da/dT for which abs species having da/dT != 0
      // This is only true for "nd" and "rh"
      //
      if( jac_is_t[iq] != Index(JacobianType::None) )
        {
          // Loop abs species, again
          for( Index ia=0; ia<jacobian_quantities.nelem(); ia++ )
            {
              if( jacobian_quantities[iq].Analytical() && jac_species_i[ia] >= 0 )
                {
                  if( jacobian_quantities[ia].Mode() == "nd" )
                    {
                      for( Index ip=0; ip<np; ip++ )
                        {
                          Matrix ddterm = diy_dpath[ia](ip,joker,joker);
                          ddterm *= ppvar_vmr(jac_species_i[ia],ip) * (
                             number_density(ppvar_p[ip],ppvar_t[ip]+1) -
                             number_density(ppvar_p[ip],ppvar_t[ip]) );
                          diy_dpath[iq](ip,joker,joker) += ddterm;
                        }
                    }
                  else if( jacobian_quantities[ia].Mode() == "rh" )
                    {
                      Tensor3 t_data(ppvar_t.nelem(),1,1);
                      t_data(joker,0,0) = ppvar_t;
                      // Calculate water sat. pressure if not already done
                      if( water_p_eq.npages() == 0 )
                        { water_p_eq_agendaExecute( ws, water_p_eq, t_data,
                                                    water_p_eq_agenda); }
                      // Sat.pressure for +1K
                      Tensor3 water_p_eq1K;
                      t_data(joker,0,0) += 1;
                      water_p_eq_agendaExecute( ws, water_p_eq1K, t_data,
                                                water_p_eq_agenda);
                      
                      for( Index ip=0; ip<np; ip++ )
                        {
                          const Numeric p_eq = water_p_eq(ip,0,0);
                          const Numeric p_eq1K = water_p_eq1K(ip,0,0);
                          Matrix ddterm = diy_dpath[ia](ip,joker,joker);
                          ddterm *= ppvar_vmr(jac_species_i[ia],ip) *
                            (ppvar_p[ip] / pow(p_eq,2.0) ) * ( p_eq1K - p_eq );
                          diy_dpath[iq](ip,joker,joker) += ddterm;
                        }
                    }
                }
            }
        }
    } 

  
  // Map to retrieval grids
  FOR_ANALYTICAL_JACOBIANS_DO
    ( 
      diy_from_path_to_rgrids( diy_dx[iq], jacobian_quantities[iq], 
                               diy_dpath[iq], atmosphere_dim, ppath, ppvar_p );
    )
}



//! rtmethods_unit_conversion
/*!
    This function handles the unit conversion to be done at the end of some
    radiative transfer WSMs. The method hanldes both *iy* and analytical parts
    of the Jacobian.

    See iyEmissonStandard for usage example.

    \author Patrick Eriksson 
    \date   2017-11-19
*/
void rtmethods_unit_conversion(
         Matrix&                     iy,
         ArrayOfTensor3&             diy_dx,
         Tensor3&                    ppvar_iy,  
   const Index&                      ns,
   const Index&                      np,
   const Vector&                     f_grid,         
   const Ppath&                      ppath,
   const ArrayOfRetrievalQuantity&   jacobian_quantities,
   const Index&                      j_analytical_do,
   const String&                     iy_unit )  
{
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
      FOR_ANALYTICAL_JACOBIANS_DO2( apply_iy_unit2( diy_dx[iq], iy, iy_unit,
                                                    f_grid, n, i_pol ); )
    } 

  // iy
  apply_iy_unit( iy, iy_unit, f_grid, n, i_pol );

  // ppvar_iy 
  for( Index ip=0; ip<ppath.np; ip++ )
    { apply_iy_unit( ppvar_iy(joker,joker,ip), iy_unit, f_grid,
                     ppath.nreal[ip], i_pol ); }
}






