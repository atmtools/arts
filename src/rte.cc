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
#include "geodetic.h"
#include "logic.h"
#include "math_funcs.h"
#include "montecarlo.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"
#include "lin_alg.h"

extern const Numeric BOLTZMAN_CONST;
extern const Numeric PLANCK_CONST;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;

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



//! apply_y_unit
/*!
    Performs conversion from radiance to other units.

    Use *apply_y_unit2* for conversion of jacobian data.

    \param   iy       In/Out: Tensor3 with data to be converted, where 
                      column dimension corresponds to Stokes dimensionality
                      and row dimension corresponds to frequency.
    \param   y_unit   As the WSV.
    \param   f_grid   As the WSV.
    \param   i_pol    Polarisation indexes. See documentation of y_pol.

    \author Patrick Eriksson 
    \date   2010-04-07
*/
void apply_y_unit( 
            MatrixView   iy, 
         const String&   y_unit, 
       ConstVectorView   f_grid,
   const ArrayOfIndex&   i_pol )
{
  // The code is largely identical between the two apply_y_unit functions.
  // If any change here, remember to update the other function.

  const Index nf = iy.nrows();
  const Index ns = iy.ncols();

  assert( f_grid.nelem() == nf );
  assert( i_pol.nelem() == ns );

  if( y_unit == "1" )
    {}

  else if( y_unit == "RJBT" )
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

  else if( y_unit == "PlanckBT" )
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
  
  else if ( y_unit == "W/(m^2 m sr)" )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          const Numeric scfac = f_grid[iv] * ( f_grid[iv] / SPEED_OF_LIGHT );
          for( Index is=0; is<ns; is++ )
            { iy(iv,is) *= scfac; }
        }
    }
  
  else if ( y_unit == "W/(m^2 m-1 sr)" )
    {
      iy *= SPEED_OF_LIGHT;
    }

  else
    {
      ostringstream os;
      os << "Unknown option: y_unit = \"" << y_unit << "\"\n" 
         << "Recognised choices are: \"1\", \"RJBT\", \"PlanckBT\""
         << "\"W/(m^2 m sr)\" and \"W/(m^2 m-1 sr)\""; 
      
      throw runtime_error( os.str() );      
    }
}



//! apply_y_unit2
/*!
    Largely as *apply_y_unit* but operates on jacobian data.

    The associated spectrum data *iy* must be in radiance. That is, the
    spectrum can only be converted to Tb after the jacobian data. 

    \param   J        In/Out: Tensor3 with data to be converted, where 
                      column dimension corresponds to Stokes dimensionality
                      and row dimension corresponds to frequency.
    \param   iy       Associated radiance data.
    \param   y_unit   As the WSV.
    \param   f_grid   As the WSV.

    \author Patrick Eriksson 
    \date   2010-04-10
*/
void apply_y_unit2( 
   Tensor3View           J,
   ConstMatrixView       iy, 
   const String&         y_unit, 
   ConstVectorView       f_grid,
   const ArrayOfIndex&   i_pol )
{
  // The code is largely identical between the two apply_y_unit functions.
  // If any change here, remember to update the other function.

  const Index nf = iy.nrows();
  const Index ns = iy.ncols();
  const Index np = J.npages();

  assert( J.nrows() == nf );
  assert( J.ncols() == ns );
  assert( f_grid.nelem() == nf );
  assert( i_pol.nelem() == ns );

  if( y_unit == "1" )
    {}

  else if( y_unit == "RJBT" )
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

  else if( y_unit == "PlanckBT" )
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

  else if ( y_unit == "W/(m^2 m sr)" )
    {
      for( Index iv=0; iv<nf; iv++ )
        {
          const Numeric scfac = f_grid[iv] * ( f_grid[iv] / SPEED_OF_LIGHT );
          for( Index ip=0; ip<np; ip++ )
            {
              for( Index is=0; is<ns; is++ )
                { J(ip,iv,is) *= scfac; }
            }
        }
    }
  
  else if ( y_unit == "W/(m^2 m-1 sr)" )
    {
      J *= SPEED_OF_LIGHT;
    }
  
  else
    {
      ostringstream os;
      os << "Unknown option: y_unit = \"" << y_unit << "\"\n" 
         << "Recognised choices are: \"1\", \"RJBT\", \"PlanckBT\""
         << "\"W/(m^2 m sr)\" and \"W/(m^2 m-1 sr)\""; 
      
      throw runtime_error( os.str() );      
    }  
}



//! bending_angle1d
/*!
    Calculates the bending angle for a 1D atmosphere.

    The expression used assumes a 1D atmosphere, that allows the bendng angle
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



// Size of disturbance for calculating defocusing.
// Used of both methods below, where +-dza applied.
// 1e-3 corresponds to roughly 100 m difference in 
// tangent altitude for limb sounding cases. 
const Numeric dza = 1e-3;


//! defocusing_general_sub
/*!
    Just to avoid duplicatuion of code in *defocusing_general*.
   
    rte_los is mainly an input, but is also returned "adjusted" (with zenith
    and azimuth angles inside defined ranges) 
 
    \param    pos                 Out: Position of ppath at optical distance lo0
    \param    rte_los             In/out: Direction for transmitted signal 
                                  (disturbed from nominal value)
    \param    rte_pos             Position of transmitter.
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
    \param    edensity_field      As the WSV with the same name.
    \param    f_index             As the WSV with the same name.
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
  ConstVectorView    rte_pos,
  const Numeric&     lo0,
  const Agenda&      ppath_step_agenda,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstVectorView    lat_grid,
  ConstVectorView    lon_grid,
  ConstTensor3View   t_field,
  ConstTensor3View   z_field,
  ConstTensor4View   vmr_field,
  ConstTensor3View   edensity_field,
  const Index&       f_index,
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
              lon_grid, t_field, z_field, vmr_field, edensity_field,
              f_index, refellipsoid, z_surface, 0, ArrayOfIndex(0), 
              rte_pos, rte_los, 0, verbosity );

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
    \param    edensity_field      As the WSV with the same name.
    \param    f_index             As the WSV with the same name.
    \param    refellipsoid        As the WSV with the same name.
    \param    z_surface           As the WSV with the same name.
    \param    ppath               As the WSV with the same name.
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
  ConstTensor3View   edensity_field,
  const Index&       f_index,
  ConstVectorView    refellipsoid,
  ConstMatrixView    z_surface,
  const Ppath&       ppath,
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

  Vector rte_los0(max(Index(1),atmosphere_dim-1)), rte_los;
  mirror_los( rte_los, ppath.start_los, atmosphere_dim );
  rte_los0 = rte_los[Range(0,max(Index(1),atmosphere_dim-1))];

  // A new ppath with positive zenith angle off-set
  //
  Vector  pos1;
  rte_los     = rte_los0;
  rte_los[0] += dza;
  //
  defocusing_general_sub( ws, pos1, rte_los, rte_pos, lo, ppath_step_agenda, 
                          atmosphere_dim, p_grid, lat_grid, lon_grid, t_field, 
                          z_field, vmr_field, edensity_field, f_index, 
                          refellipsoid, z_surface, verbosity );

  // Same thing with negative zenit angle off-set
  Vector  pos2;
  rte_los     = rte_los0;  // Use rte_los0 as rte_los can have been "adjusted"
  rte_los[0] -= dza;
  //
  defocusing_general_sub( ws, pos2, rte_los, rte_pos, lo, ppath_step_agenda, 
                          atmosphere_dim, p_grid, lat_grid, lon_grid, t_field, 
                          z_field, vmr_field, edensity_field, f_index, 
                          refellipsoid, z_surface, verbosity );

  // Calculate distance between pos1 and 2, and derive the loss factor
  Numeric l12;
  if( atmosphere_dim < 3 )
    { distance2D( l12, pos1[0], pos1[1], pos2[0], pos2[1] ); }
  else
    { distance3D( l12, pos1[0], pos1[1], pos1[2], pos2[0], pos2[1], pos2[2] ); }
  //
  dlf = lp*2*DEG2RAD*dza /  l12;
}



//! defocusing_sat2sat
/*!
    Calculates defocusing for limb measurements between two satellites.

    The expressions used assume a 1D atmosphere, and can only be applied on
    limb sounding geoemtry. The function works for 2D and 3D and should give 
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
    \param    edensity_field      As the WSV with the same name.
    \param    f_index             As the WSV with the same name.
    \param    refellipsoid        As the WSV with the same name.
    \param    z_surface           As the WSV with the same name.
    \param    ppath               As the WSV with the same name.
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
  ConstTensor3View   edensity_field,
  const Index&       f_index,
  ConstVectorView    refellipsoid,
  ConstMatrixView    z_surface,
  const Ppath&       ppath,
  const Verbosity&   verbosity )
{
  if( ppath.start_lstep == 0  ||  ppath.end_lstep == 0 )
    throw runtime_error( "The function *defocusing_sat2sat* can only be used "
                         "for satellite-to-satellite cases." );

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

  // Bending angle
  Numeric alpha;
  bending_angle1d( alpha, ppath );
  alpha *= DEG2RAD;

  // Azimuth loss term (Eq 18.5 in Kursinski et al.)
  const Numeric lf = lr*lt / (lr+lt);
  const Numeric alt = 1 / ( 1 - alpha*lf / refellipsoid[0] );

  // Calculate two new ppaths to get dalpha/da
  Numeric   alpha1, a1, alpha2, a2;
  Ppath     ppt;
  Vector    rte_pos = ppath.end_pos[Range(0,atmosphere_dim)];
  Vector    rte_los = ppath.end_los;
  //
  rte_los[0] += dza;
  adjust_los( rte_los, atmosphere_dim );
  ppath_calc( ws, ppt, ppath_step_agenda, atmosphere_dim, p_grid, lat_grid,
              lon_grid, t_field, z_field, vmr_field, edensity_field,
              f_index, refellipsoid, z_surface, 0, ArrayOfIndex(0), 
              rte_pos, rte_los, 0, verbosity );
  bending_angle1d( alpha1, ppt );
  alpha1 *= DEG2RAD;
  find_tanpoint( it, ppt );
  a1 = ppt.pos(it,0);
  //
  rte_los[0] -= 2*dza;
  adjust_los( rte_los, atmosphere_dim );
  ppath_calc( ws, ppt, ppath_step_agenda, atmosphere_dim, p_grid, lat_grid,
              lon_grid, t_field, z_field, vmr_field, edensity_field,
              f_index, refellipsoid, z_surface, 0, ArrayOfIndex(0), 
              rte_pos, rte_los, 0, verbosity );
  bending_angle1d( alpha2, ppt );
  alpha2 *= DEG2RAD;
  find_tanpoint( it, ppt );
  a2 = ppt.pos(it,0);
  //
  const Numeric dada = (alpha2-alpha1) / (a2-a1); 

  // Zenith loss term (Eq 18 in Kursinski et al.)
  const Numeric zlt = 1 / ( 1 - dada*lf );

  // Total defocusing loss
  dlf = zlt * alt;
}



//! ext2trans
/*!
    Converts an extinction matrix to a transmission matrix

    The function performs the calculations differently depending on the
    conditions to improve the speed. There are three cases: <br>
       1. Scalar absorption. <br>
       2. The matrix ext_mat_av is diagonal. <br>
       3. The total general case.

    \param   trans_mat          Input/Output: Transmission matrix of slab.
    \param   ext_mat            Input: Averaged extinction matrix.
    \param   lstep             Input: The length of the RTE step.

    \author Claudia Emde and Patrick Eriksson, 
    \date   2010-10-15
*/
void ext2trans(//Output and Input:
              MatrixView trans_mat,
              //Input
              ConstMatrixView ext_mat_av,
              const Numeric& lstep )
{
  //Stokes dimension:
  Index stokes_dim = ext_mat_av.nrows();

  //Check inputs:
  assert( is_size(trans_mat, stokes_dim, stokes_dim) );
  assert( is_size(ext_mat_av, stokes_dim, stokes_dim) );
  assert( lstep >= 0 );
  assert( !is_singular( ext_mat_av ) );

  // Any changes here should also be implemented in rte_step_std.

  //--- Scalar case: ---------------------------------------------------------
  if( stokes_dim == 1 )
    {
      trans_mat(0,0) = exp(-ext_mat_av(0,0) * lstep);
    }


  //--- Vector case: ---------------------------------------------------------

  //- Unpolarised
  else if( is_diagonal(ext_mat_av) )
    {
      const Numeric tv = exp(-ext_mat_av(0,0) * lstep);

      trans_mat = 0;

      for( Index i=0; i<stokes_dim; i++ )
        {
          trans_mat(i,i)  = tv;
        }
    }


  //- General case
  else
    {
      Matrix ext_mat_ds = ext_mat_av;
      ext_mat_ds *= -lstep; 

      Index q = 10;  // index for the precision of the matrix exp function

      matrix_exp( trans_mat, ext_mat_ds, q );
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
    \param   iy_main_agenda    As the WSV.

    \author Patrick Eriksson 
    \date   2012-08-08
*/
void get_iy(
         Workspace&   ws,
         Matrix&      iy,
   ConstTensor3View   t_field,
   ConstTensor3View   z_field,
   ConstTensor4View   vmr_field,
   const Index&       cloudbox_on,
   ConstVectorView    rte_pos,
   ConstVectorView    rte_los,
   const Agenda&      iy_main_agenda )
{
  ArrayOfTensor3    diy_dx;
  ArrayOfTensor4    iy_aux;
  Ppath             ppath;
  Tensor3           iy_transmission(0,0,0);

  iy_main_agendaExecute( ws, iy, iy_aux, ppath, diy_dx, 1, iy_transmission, 
                             ArrayOfString(0), cloudbox_on, 0, t_field, z_field,
                             vmr_field, -1, rte_pos, rte_los, 
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
    can involve recursive calls of *iy_main_agenda*. It is also
    allowed to input *iy_clearsky_basic_agenda* instead of
    *iy_main_agenda*.

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
    \param   iy_main_agenda    As the WSV.
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
  const Index&            jacobian_do,
  const Ppath&            ppath,
  const Index&            atmosphere_dim,
  ConstTensor3View        t_field,
  ConstTensor3View        z_field,
  ConstTensor4View        vmr_field,
  const Index&            cloudbox_on,
  const Index&            stokes_dim,
  ConstVectorView         f_grid,
  const Agenda&           iy_main_agenda,
  const Agenda&           iy_space_agenda,
  const Agenda&           iy_surface_agenda,
  const Agenda&           iy_cloudbox_agenda,
  const Verbosity&        verbosity)
{
  CREATE_OUT3;
  
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
  Vector rte_pos, rte_los;
  rte_pos.resize( atmosphere_dim );
  rte_pos = ppath.pos(np-1,Range(0,atmosphere_dim));
  rte_los.resize( ppath.los.ncols() );
  rte_los = ppath.los(np-1,joker);

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
        iy_space_agendaExecute( ws, iy, rte_pos, rte_los, f_grid,
                                iy_space_agenda );
      }
      break;

    case 2:   //--- The surface -----------------------------------------------
      {
        agenda_name = "iy_surface_agenda";
        iy_surface_agendaExecute( ws, iy, diy_dx, iy_transmission, cloudbox_on,
                                  jacobian_do, t_field, z_field, vmr_field,
                                  iy_main_agenda, rte_pos, rte_los, 
                                  iy_surface_agenda );
      }
      break;

    case 3:   //--- Cloudbox boundary or interior ------------------------------
    case 4:
      {
        agenda_name = "iy_cloudbox_agenda";
        iy_cloudbox_agendaExecute( ws, iy, rte_pos, rte_los, 
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
    Determines pressure, temperature, VMR and winds for each propgataion path
    point.

    The output variables are sized inside the function. For VMR the
    dimensions are [ species, propagation path point ].

    \param   ppath_p           Out: Pressure for each ppath point.
    \param   ppath_t           Out: Temperature for each ppath point.
    \param   ppath_vmr         Out: VMR values for each ppath point.
    \param   ppath_wind_u      Out: U-wind for each ppath point.
    \param   ppath_wind_v      Out: V-wind for each ppath point.
    \param   ppath_wind_w      Out: W-wind for each ppath point.
    \param   ppath_mag_u       Out: U-mag for each ppath point.
    \param   ppath_mag_v       Out: V-mag for each ppath point.
    \param   ppath_mag_w       Out: W-mag for each ppath point.
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
        Matrix&      ppath_vmr, 
        Vector&      ppath_wind_u, 
        Vector&      ppath_wind_v, 
        Vector&      ppath_wind_w,
        Vector&      ppath_mag_u,
        Vector&      ppath_mag_v,
        Vector&      ppath_mag_w,
  const Ppath&       ppath,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstTensor3View   t_field,
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
                              vmr_field( is, joker, joker,  joker ), 
                              ppath.gp_p, ppath.gp_lat, ppath.gp_lon, 
                              itw_field );
    }

  // Winds:
  //
  ppath_wind_w.resize(np);
  if( wind_w_field.npages() > 0 ) 
    { 
      interp_atmfield_by_itw( ppath_wind_w,  atmosphere_dim, wind_w_field, 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }
  else
    { ppath_wind_w = 0; }
  //
  ppath_wind_v.resize(np);
  if( wind_v_field.npages() > 0 ) 
    { 
      interp_atmfield_by_itw( ppath_wind_v,  atmosphere_dim, wind_v_field, 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }
  else
    { ppath_wind_v = 0; }
  //
  ppath_wind_u.resize(np);
  if( atmosphere_dim > 2 )
    {
      if( wind_u_field.npages() > 0 ) 
        { 
          interp_atmfield_by_itw( ppath_wind_u,  atmosphere_dim, wind_u_field, 
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
        }
      else
        { ppath_wind_u = 0; }
    }
  else
    { ppath_wind_u = 0; }

  // Magnetic field:
  //
  ppath_mag_w.resize(np);
  if( mag_w_field.npages() > 0 )
    {
      interp_atmfield_by_itw( ppath_mag_w,  atmosphere_dim, mag_w_field,
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }
  else
    { ppath_mag_w = 0; }
  //
  ppath_mag_v.resize(np);
  if( mag_v_field.npages() > 0 )
    {
      interp_atmfield_by_itw( ppath_mag_v,  atmosphere_dim, mag_v_field,
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
    }
  else
    { ppath_mag_v = 0; }
  //
  ppath_mag_u.resize(np);
  if( atmosphere_dim > 2 )
    {
      if( mag_u_field.npages() > 0 )
        {
          interp_atmfield_by_itw( ppath_mag_u,  atmosphere_dim, mag_u_field,
                          ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
        }
      else
        { ppath_mag_u = 0; }
    }
  else
    { ppath_mag_u = 0; }
}


//! get_ppath_pnd
/*!
    Determines particle number densities for each propgataion path
    point.

    The output variable is sized inside the function as
    [ particle type, propagation path point ].

    \param   ppath_pnd         Out: PND for each ppath point.
    \param   ppath             As the WSV.
    \param   atmosphere_dim    As the WSV.
    \param   pnd_field         As the WSV.

    \author Patrick Eriksson 
    \date   2011-07-14
*/
void get_ppath_pnd( 
        Matrix&         ppath_pnd, 
  const Ppath&          ppath,
  const Index&          atmosphere_dim,
  const ArrayOfIndex&   cloudbox_limits,
  ConstTensor4View      pnd_field )
{
  ppath_pnd.resize( pnd_field.nbooks(), ppath.np );

  Matrix   itw_field;
  ArrayOfGridPos   gpc_p, gpc_lat, gpc_lon;
  //
  interp_cloudfield_gp2itw( itw_field, gpc_p, gpc_lat, gpc_lon, 
                            ppath.gp_p, ppath.gp_lat, ppath.gp_lon, 
                            atmosphere_dim, cloudbox_limits );
  for( Index ip=0; ip<pnd_field.nbooks(); ip++ )
    {
      interp_atmfield_by_itw( ppath_pnd(ip,joker), atmosphere_dim,
                              pnd_field(ip,joker,joker,joker), 
                              gpc_p, gpc_lat, gpc_lon, itw_field );
    }
}


//! get_ppath_rtvars
/*!
    Determines variables for each step of "standard" RT integration

    The output variables are sized inside the function. The dimension order is 
       [ absorption species, frequency, ppath point ]

    \param   ws                Out: The workspace
    \param   ppath_abs_mat     Out: Absorption coefficients at each ppath point
    \param   ppath_tau         Out: Optical thickness of each ppath step 
    \param   total_tau         Out: Total optical thickness of path
    \param   ppath_emission    Out: Emission source term at each ppath point 
    \param   abs_mat_per_species_agenda As the WSV.    
    \param   blackbody_radiation_agenda   As the WSV.    
    \param   ppath_p           Pressure for each ppath point.
    \param   ppath_t           Temperature for each ppath point.
    \param   ppath_vmr         VMR values for each ppath point.
    \param   ppath_wind_u      U-wind for each ppath point.
    \param   ppath_wind_v      V-wind for each ppath point.
    \param   ppath_wind_w      W-wind for each ppath point.
    \param   ppath_mag_u       Out: U-mag for each ppath point.
    \param   ppath_mag_v       Out: V-mag for each ppath point.
    \param   ppath_mag_w       Out: W-mag for each ppath point.
    \param   f_grid            As the WSV.    
    \param   f_index           As the WSV.
    \param   atmosphere_dim    As the WSV.    
    \param   emission_do       Flag for calculation of emission. Should be
                               set to 0 for pure transmission calculations.

    \author Patrick Eriksson 
    \date   2009-10-06
*/
void get_ppath_rtvars( 
        Workspace&   ws,
        Tensor5&     ppath_abs,
        Tensor4&     ppath_tau,
        Tensor3&     total_tau,
        Matrix&      ppath_emission,
  const Agenda&      abs_mat_per_species_agenda,
  const Agenda&      blackbody_radiation_agenda,
  const Ppath&       ppath,
  ConstVectorView    ppath_p, 
  ConstVectorView    ppath_t, 
  ConstMatrixView    ppath_vmr, 
  ConstVectorView    ppath_wind_u, 
  ConstVectorView    ppath_wind_v, 
  ConstVectorView    ppath_wind_w,
  ConstVectorView    ppath_mag_u,
  ConstVectorView    ppath_mag_v,
  ConstVectorView    ppath_mag_w,
  ConstVectorView    f_grid, 
  const Index&       f_index, 
  const Index&       stokes_dim,
  const Index&       atmosphere_dim,
  const Index&       emission_do )
{
  // Sizes
  const Index   np   = ppath.np;
  const Index   nabs = ppath_vmr.nrows();

  // Frequencies
  Numeric fd;  // Frequency to apply for Doppler shift
  Index   nf;
  if( f_index < 0 )
    { 
      nf = f_grid.nelem(); 
      fd = ( f_grid[0] + f_grid[nf-1] ) / 2.0;
    }
  else
    { 
      nf = 1; 
      fd = f_grid[f_index];
    }
    
  // Init variables
  ppath_abs.resize( nabs, nf, stokes_dim, stokes_dim, np   );
  ppath_tau.resize(       nf, stokes_dim, stokes_dim, np-1 );
  total_tau.resize(       nf, stokes_dim, stokes_dim       );
  total_tau = 0;
  if( emission_do )
    { ppath_emission.resize( nf, np ); }
  else
    { ppath_emission.resize( 0, 0 ); }
    
  for( Index ip=0; ip<np; ip++ )
    {
        // Doppler shift
        //
        Numeric rte_doppler = 0;
        //
        if( ppath_wind_v[ip]!=0 || ppath_wind_u[ip]!=0 || ppath_wind_w[ip]!=0 )
            {
            // za and aa of winds and LOS
            const Numeric v = sqrt( ppath_wind_u[ip]*ppath_wind_u[ip] +
                                    ppath_wind_v[ip]*ppath_wind_v[ip] +
                                    ppath_wind_w[ip]*ppath_wind_w[ip] );

            Numeric aa_w = 0, aa_p = 0; // Azimuth for wind and ppath-los, in RAD
            if( atmosphere_dim < 3 )
                {
                if( ppath_wind_v[ip]<0 )
                    { aa_w = PI; }  // Negative v-winds for 1 and 2D
                if( atmosphere_dim == 2  &&  ppath.los(ip,0) < 0 )
                    { aa_p = PI; }  // Negative sensor za for 2D
                }
            else //3D
                {
                aa_w = atan2( ppath_wind_u[ip], ppath_wind_v[ip] );
                aa_p = DEG2RAD * ppath.los(ip,1);
                }
            const Numeric za_w = acos( ppath_wind_w[ip] / v );
            const Numeric za_p = DEG2RAD * fabs( ppath.los(ip,0) );
            //
            // Actual shift
            const Numeric costerm = cos(za_w) * cos(za_p) +
                                    sin(za_w) * sin(za_p)*cos(aa_w-aa_p);

            rte_doppler = -v * fd * costerm / SPEED_OF_LIGHT;
            }
        
        // Magnetic field
        //
        Vector rte_mag(1,-1.);
         if( ppath_mag_v[ip]!=0 || ppath_mag_u[ip]!=0 || ppath_mag_w[ip]!=0 )
         {
        rte_mag.resize(3); //FIXME: email Patrick about u,v,w order
            rte_mag[0] = ppath_mag_u[ip];
            rte_mag[1] = ppath_mag_v[ip];
            rte_mag[2] = ppath_mag_w[ip];
         }

        // We must use temporary variables as the agenda input must be
        // free to be resized
        Vector   evector;
        Tensor4  sgtensor4;
        //
        
        abs_mat_per_species_agendaExecute(ws, sgtensor4, f_index, rte_doppler, rte_mag,
                                            ppath_p[ip], ppath_t[ip], ppath_vmr(joker,ip),
                                            abs_mat_per_species_agenda );

        ppath_abs(joker,joker,joker,joker,ip) = sgtensor4;

        //
        if( emission_do )
            {
            if( f_index < 0 )
                { blackbody_radiation_agendaExecute( ws, evector, ppath_t[ip],
                                            f_grid, blackbody_radiation_agenda ); }
            else
                { blackbody_radiation_agendaExecute( ws, evector, ppath_t[ip],
                        Vector(1,f_grid[f_index]), blackbody_radiation_agenda ); }
            ppath_emission(joker,ip) = evector;
            }
    
        // Partial and total tau
        //
        if( ip > 0 ){
            for( Index is1=0; is1<stokes_dim; is1++ ){
                for( Index is2=0; is2<stokes_dim; is2++ ){
                    for( Index iv=0; iv<nf; iv++ ){
                        ppath_tau(iv,is1,is2,ip-1) = 0.5 * ppath.lstep[ip-1] * (
                                                    ppath_abs(joker, iv, is1, is2, ip-1).sum() +
                                                    ppath_abs(joker, iv, is1, is2, ip).sum() );
                        total_tau(iv,is1,is2) += ppath_tau(iv, is1, is2, ip-1);
                    }
                }
            }
        }
    }
}



//! get_ppath_cloudrtvars
/*!
    Determines variables for each step of "standard" RT integration

    See code for details about dimensions etc.

    \param   ws                  Out: The workspace
    \param   ppath_asp_abs_vec   Out: Absorption vectors for absorption species
    \param   ppath_asp_ext_vec   Out: Extinction matrices for absorption species
    \param   ppath_pnd_abs_vec   Out: Absorption vectors for particles
    \param   ppath_pnd_ext_vec   Out: Extinction matrices for particles
    \param   ppath_transmission  Out: Transmissions of each ppath step 
    \param   total_transmission  Out: Total transmission of path
    \param   ppath_emission      Out: Emission source term at each ppath point 
    \param   scat_data           Out: Extracted scattering data. Length of
                                      array affected by *use_mean_scat_data*.
    \param   abs_mat_per_species_agenda As the WSV.
    \param   blackbody_radiation_agenda     As the WSV
    \param   ppath_p             Pressure for each ppath point.
    \param   ppath_t             Temperature for each ppath point.
    \param   ppath_vmr           VMR values for each ppath point.
    \param   ppath_wind_u        U-wind for each ppath point.
    \param   ppath_wind_v        V-wind for each ppath point.
    \param   ppath_wind_w        W-wind for each ppath point.
    \param   ppath_wind_u        U-mag for each ppath point.
    \param   ppath_wind_v        V-mag for each ppath point.
    \param   ppath_wind_w        W-mag for each ppath point.
    \param   ppath_pnd           Particle number densities for each ppath point.
    \param   use_mean_scat_data  As the WSV.    
    \param   scat_data_raw       As the WSV.    
    \param   stokes_dim          As the WSV.    
    \param   f_grid              As the WSV.    
    \param   atmosphere_dim      As the WSV.    
    \param   emission_do         Flag for calculation of emission. Should be
                                 set to 0 for pure transmission calculations.

    \author Patrick Eriksson 
    \date   2011-07-14
*/
void get_ppath_cloudrtvars( 
        Workspace&                     ws,
        Tensor3&                       ppath_asp_abs_vec, 
        Tensor4&                       ppath_asp_ext_mat, 
        Tensor3&                       ppath_pnd_abs_vec, 
        Tensor4&                       ppath_pnd_ext_mat, 
        Tensor4&                       ppath_transmission,
        Tensor3&                       total_transmission,
        Matrix&                        ppath_emission, 
  Array<ArrayOfSingleScatteringData>&  scat_data,
  const Agenda&                        abs_mat_per_species_agenda,
  const Agenda&                        blackbody_radiation_agenda,
  const Ppath&                         ppath,
  ConstVectorView                      ppath_p, 
  ConstVectorView                      ppath_t, 
  ConstMatrixView                      ppath_vmr, 
  ConstVectorView                      ppath_wind_u, 
  ConstVectorView                      ppath_wind_v, 
  ConstVectorView                      ppath_wind_w,
  ConstVectorView                      ppath_mag_u,
  ConstVectorView                      ppath_mag_v,
  ConstVectorView                      ppath_mag_w,
  ConstMatrixView                      ppath_pnd, 
  const Index&                         use_mean_scat_data,
  const ArrayOfSingleScatteringData&   scat_data_raw,
  const Index&                         stokes_dim,
  ConstVectorView                      f_grid, 
  const Index&                         atmosphere_dim,
  const Index&                         emission_do,
  const Verbosity&                     verbosity)
{
  // Sizes
  const Index   np   = ppath.np;
  const Index   nf   = f_grid.nelem();

  // Init variables
  ppath_asp_abs_vec.resize( nf, stokes_dim, np );
  ppath_asp_ext_mat.resize( nf, stokes_dim, stokes_dim, np );
  ppath_pnd_abs_vec.resize( nf, stokes_dim, np );
  ppath_pnd_ext_mat.resize( nf, stokes_dim, stokes_dim, np );
  ppath_transmission.resize( nf, stokes_dim, stokes_dim, np-1 );
  total_transmission.resize( nf, stokes_dim, stokes_dim );
  //
  if( emission_do )
    { ppath_emission.resize( nf, np ); }
  else
    { ppath_emission.resize( 0, 0 ); }

  // Mean of extreme frequencies
  const Numeric f0 = (f_grid[0]+f_grid[nf-1])/2.0;


  // Particle single scattering properties (are independent of position)
  //
  if( use_mean_scat_data )
    {
      scat_data.resize( 1 );
      scat_data_monoCalc( scat_data[0], scat_data_raw, Vector(1,f0), 0, 
                          verbosity );
    }
  else
    {
      scat_data.resize( nf );
      for( Index iv=0; iv<nf; iv++ )
        { scat_data_monoCalc( scat_data[iv], scat_data_raw, f_grid, iv, 
                              verbosity ); 
        }
    }


  // Loop ppath points
  //
  for( Index ip=0; ip<np; ip++ )
    {
      // Doppler shift 
      // 
      assert( ppath_wind_u[ip] == 0 );
      assert( ppath_wind_v[ip] == 0 );
      assert( ppath_wind_w[ip] == 0 );
      //
      const Numeric rte_doppler = 0;

      // Magnetic field
      //
      assert( ppath_mag_u[ip] == 0 );
      assert( ppath_mag_v[ip] == 0 );
      assert( ppath_mag_w[ip] == 0 );
      //
      const Vector rte_mag(1,-1.);

      // Emission source term
      //
      if( emission_do )
        {
          Vector   evector;   // Agenda must be free to resize
          blackbody_radiation_agendaExecute( ws, evector, ppath_t[ip], f_grid,
                                  blackbody_radiation_agenda );
          ppath_emission(joker,ip) = evector;
        }

      // Absorption species properties 
        {
          Tensor4  asgtensor4;   // Agendas must be free to resize
          Matrix   abs_vec( nf, stokes_dim, 0 );
          Tensor3  ext_mat( nf, stokes_dim, stokes_dim, 0 );
          //
          abs_mat_per_species_agendaExecute( ws, asgtensor4, -1, rte_doppler,
                                             rte_mag,
                                             ppath_p[ip], ppath_t[ip], 
                                             ppath_vmr(joker,ip),
                                             abs_mat_per_species_agenda );

        opt_prop_add_abs_mat_per_species(ext_mat, abs_vec, asgtensor4, -1);
         
          ppath_asp_ext_mat(joker,joker,joker,ip) = ext_mat;
          ppath_asp_abs_vec(joker,joker,ip)       = abs_vec;
        }

      // Particle properties 
        {
          // Direction of outgoing scattered radiation (which is reversed to
          // LOS). Note that rte_los2 is only used for extracting scattering
          // properties.
          Vector rte_los2;
          mirror_los( rte_los2, ppath.los(ip,joker), atmosphere_dim );

          // Extinction and absorption
          if( use_mean_scat_data )
            {
              Vector   abs_vec( stokes_dim, 0 );
              Matrix   ext_mat( stokes_dim, stokes_dim, 0 );
              opt_propCalc( ext_mat, abs_vec,
                            rte_los2[0], rte_los2[1], scat_data[0], 
                            stokes_dim, ppath_pnd(joker,ip), ppath_t[ip],
                            verbosity);
              for( Index iv=0; iv<nf; iv++ )
                { 
                  ppath_pnd_ext_mat(iv,joker,joker,ip) = ext_mat;
                  ppath_pnd_abs_vec(iv,joker,ip)       = abs_vec;
                }
            }
          else
            {
              for( Index iv=0; iv<nf; iv++ )
                { 
                  Vector   abs_vec( stokes_dim, 0 );
                  Matrix   ext_mat( stokes_dim, stokes_dim, 0 );
                  opt_propCalc( ext_mat, abs_vec,
                                rte_los2[0], rte_los2[1], scat_data[iv], 
                                stokes_dim, ppath_pnd(joker,ip), ppath_t[ip],
                                verbosity );
                  ppath_pnd_ext_mat(iv,joker,joker,ip) = ext_mat;
                  ppath_pnd_abs_vec(iv,joker,ip)       = abs_vec;
                }
            }
        }
       
      // Transmission
      //
      if( ip == 0 )   // Set total_transmission to identity matrices
        {
          for( Index iv=0; iv<nf; iv++ )
            { id_mat( total_transmission(iv,joker,joker) ); } 
        }
      else
        {
          for( Index iv=0; iv<nf; iv++ )
            { 
              // Average extinction matrix
              Matrix  ext_mat_av( stokes_dim, stokes_dim,0 );
              for( Index is1=0; is1<stokes_dim; is1++ )
                { 
                  for( Index is2=0; is2<stokes_dim; is2++ )
                    {
                      ext_mat_av(is1,is2) = 0.5 * (
                                          ppath_asp_ext_mat(iv,is1,is2,ip-1) +
                                          ppath_asp_ext_mat(iv,is1,is2,ip)   +
                                          ppath_pnd_ext_mat(iv,is1,is2,ip-1) +
                                          ppath_pnd_ext_mat(iv,is1,is2,ip) );
                    }
                }
              // Transmission
              ext2trans( ppath_transmission(iv,joker,joker,ip-1),
                         ext_mat_av, ppath.lstep[ip-1] );  
              const Matrix tmp = total_transmission(iv,joker,joker);
              mult( total_transmission(iv,joker,joker), tmp,
                    ppath_transmission(iv,joker,joker,ip-1) );
            }
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
  const Index&                      mblock_index,
  const Index&                      atmosphere_dim,
  ConstTensor3View                  t_field,
  ConstTensor3View                  z_field,
  ConstTensor4View                  vmr_field,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  ConstVectorView                   f_grid,
  ConstMatrixView                   sensor_pos,
  ConstMatrixView                   sensor_los,
  ConstVectorView                   mblock_za_grid,
  ConstVectorView                   mblock_aa_grid,
  const Index&                      antenna_dim,
  const Agenda&                     iy_main_agenda,
  const String&                     y_unit,
  const Index&                      j_analytical_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices,
  const ArrayOfString&              iy_aux_vars,
  const Verbosity&                  verbosity )
{
  // Sizes
  const Index   nf   = f_grid.nelem();
  const Index   nza  = mblock_za_grid.nelem();
        Index   naa  = mblock_aa_grid.nelem();   
  if( antenna_dim == 1 )  
    { naa = 1; }
  const Index   niyb = nf * nza * naa * stokes_dim;
  // Set up size of containers for data of 1 measurement block.
  // (can not be made below due to parallalisation)
  iyb.resize( niyb );
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

  // For iy_aux we don't know the number of quantities, and we have to store
  // all outout
  ArrayOfArrayOfTensor4  iy_aux_array( nza*naa );

  // Polarisation index variable
  ArrayOfIndex i_pol(stokes_dim);
  for( Index is=0; is<stokes_dim; is++ )
    { i_pol[is] = is + 1; }

  // We have to make a local copy of the Workspace and the agendas because
  // only non-reference types can be declared firstprivate in OpenMP
  Workspace l_ws (ws);
  Agenda l_iy_main_agenda (iy_main_agenda);
  
  // Start of actual calculations
/*#pragma omp parallel for                                        \
  if(!arts_omp_in_parallel() && nza>1)                            \
  default(none)                                                   \
  firstprivate(l_ws, l_iy_main_agenda)                        \
  shared(sensor_los, mblock_za_grid, mblock_aa_grid, vmr_field,   \
         t_field, f_grid, sensor_pos, \
         joker, naa) */
#pragma omp parallel for                                          \
  if(!arts_omp_in_parallel() && nza>1)                            \
  firstprivate(l_ws, l_iy_main_agenda)
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
              los     = sensor_los( mblock_index, joker );
              los[0] += mblock_za_grid[iza];
              //
              if( antenna_dim == 2 )  // map_daa handles also "adjustment"
                { map_daa( los[0], los[1], los[0], los[1], 
                                                       mblock_aa_grid[iaa] ); }
              else 
                { adjust_los( los, atmosphere_dim ); }

              // Calculate iy and associated variables
              //
              Matrix         iy;
              ArrayOfTensor3 diy_dx;
              Ppath          ppath;
              Tensor3        iy_transmission(0,0,0);
              Index          iang = iza*naa + iaa;
              //
              iy_main_agendaExecute( l_ws, iy, iy_aux_array[iang], ppath,
                                         diy_dx, 1, iy_transmission, 
                                         iy_aux_vars, cloudbox_on, 
                                         j_analytical_do, t_field, z_field, 
                                         vmr_field, mblock_index, 
                                         sensor_pos(mblock_index,joker), los, 
                                         l_iy_main_agenda );

              // Check that aux data can be handled
              for( Index q=0; q<iy_aux_array[iang].nelem(); q++ )
                {
                  if( iy_aux_array[iang][q].ncols() > 1  ||  
                      iy_aux_array[iang][q].nrows() > 1 )
                    { throw runtime_error( "For calculations using yCalc, "
                       "*iy_aux_vars* can not include varaibles of "
                       "along-the-path variables or extinction matrix type."); }
                }              

              // Start row in iyb etc. for present LOS
              //
              const Index row0 = iang * nf * stokes_dim;

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

              // iy : unit conversion and copy to iyb
              //
              apply_y_unit( iy, y_unit, f_grid, i_pol );
              //
              for( Index is=0; is<stokes_dim; is++ )
                { iyb[Range(row0+is,nf,stokes_dim)] = iy(joker,is); }

            }  // End aa loop
        }  // End try

      catch (runtime_error e)
        {
          CREATE_OUT0;
          exit_or_rethrow(e.what(), out0);
        }
    }  // End za loop


  // Compile iyb_aux
  //
  const Index nq = iy_aux_array[0].nelem();
  iyb_aux.resize( nq );
  //
  for( Index q=0; q<nq; q++ )
    {
      iyb_aux[q].resize( niyb );
      //
      for( Index iza=0; iza<nza; iza++ )
        {
          for( Index iaa=0; iaa<naa; iaa++ )
            {
              const Index iang = iza*naa + iaa;
              const Index row0 = iang * nf * stokes_dim;
              for( Index iv=0; iv<nf; iv++ )
                { 
                  const Index row1 = row0 + iv*stokes_dim;
                  const Index i1 = min( iv, iy_aux_array[iang][q].nbooks()-1 );
                  for( Index is=0; is<stokes_dim; is++ )
                    { 
                      Index i2 = min( is, iy_aux_array[iang][q].npages()-1 );
                      iyb_aux[q][row1+is] = iy_aux_array[iang][q](i1,i2,0,0);
                    }
                }
            }
        }
    }
}



//! iy_transmission_for_scalar_tau
/*!
    Sets *iy_transmission* to match scalar optical thicknesses.

    *iy_transmission* is sized by the function.

    \param   iy_transmission   Out: As the WSV.
    \param   stokes_dim        As the WSV.
    \param   tau               Vector with the optical thickness for each 
                               frequency.

    \author Patrick Eriksson 
    \date   2009-10-06
*/
void iy_transmission_for_tensor3_tau(
       Tensor3&     iy_transmission,
  const Index&      stokes_dim,
  ConstTensor3View  tau )
{
  Matrix temp_tau;
  iy_transmission.resize( tau.npages(), stokes_dim, stokes_dim );
  iy_transmission = 0;
  
  for( Index iv=0; iv<tau.npages(); iv++ )
    { 
        temp_tau  = tau(iv, joker, joker);
        temp_tau *= -1;
        matrix_exp(iy_transmission(iv, joker, joker), temp_tau, 10 ); //Note that 10 is just taken from another place this was called.
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

    \param   iy_trans_new      Out: Updated version of *iy_transmission*
    \param   iy_transmission   As the WSV.
    \param   trans             A variable matching *iy_transmission.

    \author Patrick Eriksson 
    \date   2009-10-06
*/
void iy_transmission_mult( 
       Tensor3&      iy_trans_new,
  ConstTensor3View   iy_transmission,
  ConstTensor3View   trans )
{
  const Index nf = iy_transmission.npages();
  const Index ns = iy_transmission.ncols();

  assert( ns == iy_transmission.nrows() );
  assert( nf == trans.npages() );
  assert( ns == trans.nrows() );
  assert( ns == trans.ncols() );

  iy_trans_new.resize( nf, ns, ns );

  for( Index iv=0; iv<nf; iv++ )
    {
      mult( iy_trans_new(iv,joker,joker), iy_transmission(iv,joker,joker),
                                                    trans(iv,joker,joker) );
    } 
}



//! iy_transmission_mult_scalar_tau
/*!
    Multiplicates iy_transmission with scalar optical thicknesses.

    The functions can incorporate the transmission of a clear-sky
    path. That is, the transmission can be described by a single value
    The transmission of this path is gives as the optical depth for
    each frequency.

    The "new path" is assumed to be further away from the sensor than 
    the propagtion path already included in iy_transmission. That is,
    the operation can be written as:
    
       Ttotal = Told * Tnew

    where Told is the transmission corresponding to *iy_transmission*
    and Tnew corresponds to *tau*.

    *iy_trans_new* is sized by the function.

    \param   iy_trans_new      Out: Updated version of *iy_transmission*
    \param   iy_transmission   As the WSV.
    \param   tau               Vector with the optical thickness for each 
                               frequency.

    \author Patrick Eriksson 
    \date   2009-10-06
*/
void iy_transmission_mult_tensor3_tau(
       Tensor3&      iy_trans_new,
  ConstTensor3View   iy_transmission,
  ConstTensor3View   tau )
{
  const Index nf = iy_transmission.npages();
  Matrix temp_trans(tau.nrows(), tau.ncols(), 0.0), temp_tau;
  assert( iy_transmission.ncols() == iy_transmission.nrows() );
  assert( nf == tau.npages() );

  iy_trans_new = iy_transmission;

  for( Index iv=0; iv<nf; iv++ )
  {
      temp_tau  = tau(iv, joker, joker);
      temp_tau *= -1;
      matrix_exp(temp_trans, temp_tau, 10 ); //Note that 10 is just taken from one other place this was called.
      iy_trans_new(iv,joker,joker) *= temp_trans;
  }
  
}



//! los3d
/*!
    Converts any LOS vector to the implied 3D LOS vector.

    The output argument, *los3d*, is a vector with length 2, with azimuth angle
    set and zenith angle always >= 0. 

    \param   los3d             Out: The line-of-sight in 3D
    \param   los               A line-of-sight
    \param   atmosphere_dim    As the WSV.

    \author Patrick Eriksson 
    \date   2012-07-10
*/
void los3d(
        Vector&     los3d,
  ConstVectorView   los, 
  const Index&      atmosphere_dim )
{
  los3d.resize(2);
  //
  los3d[0] = abs( los[0] ); 
  //
  if( atmosphere_dim == 1 )
    { los3d[1] = 0; }
  else if( atmosphere_dim == 2 )
    {
      if( los[0] >= 0 )
        { los3d[1] = 0; }
      else
        { los3d[1] = 180; }
    }
  else if( atmosphere_dim == 3 )
    { los3d[1] = los[1]; }
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


//! rte_step_std
/*!
    Solves monochromatic VRTE for an atmospheric slab with constant 
    conditions.

    The function can be used for clearsky and cloudbox calculations.

    The function is best explained by considering a homogenous layer. That is,
    the physical conditions inside the layer are constant. In reality they
    are not constant, so in practical all coefficients have to be averaged 
    before calling this function. 
    Total extinction and absorption inside the layer are described by
    *ext_mat_av* and *abs_vec_av* respectively,
    the blackbdody radiation of the layer is given by *rte_planck_value*
    and the propagation path length through the layer is *lstep*.

    There is an additional scattering source term in the 
    VRTE, the scattering integral term. For this function a constant
    scattering term is assumed. The radiative transfer step is only a part 
    the iterative solution of the scattering problem, for more 
    information consider AUG. In the clearsky case this variable has to be
    set to 0.

    When calling the function, the vector *stokes_vec* shall contain the
    Stokes vector for the incoming radiation. The function returns this
    vector, then containing the outgoing radiation on the other side of the 
    layer.

    The function performs the calculations differently depending on the
    conditions to improve the speed. There are three cases: <br>
       1. Scalar absorption (stokes_dim = 1). <br>
       2. The matrix ext_mat_gas is diagonal (unpolarised absorption). <br>
       3. The total general case.

    \param   stokes_vec         Input/Output: A Stokes vector.
    \param   trans_mat          Input/Output: Transmission matrix of slab.
    \param   ext_mat_av         Input: Averaged extinction matrix.
    \param   abs_vec_av         Input: Averaged absorption vector.
    \param   sca_vec_av         Input: averaged scattering vector.
    \param   lstep             Input: The length of the RTE step.
    \param   rte_planck_value   Input: Blackbody radiation.

    \author Claudia Emde and Patrick Eriksson, 
    \date   2002-11-22
*/
void rte_step_std(//Output and Input:
              VectorView stokes_vec,
              MatrixView trans_mat,
              //Input
              ConstMatrixView ext_mat_av,
              ConstVectorView abs_vec_av,
              ConstVectorView sca_vec_av,
              const Numeric& lstep,
              const Numeric& rte_planck_value )
{
  //Stokes dimension:
  Index stokes_dim = stokes_vec.nelem();

  //Check inputs:
  assert(is_size(trans_mat, stokes_dim, stokes_dim));
  assert(is_size(ext_mat_av, stokes_dim, stokes_dim));
  assert(is_size(abs_vec_av, stokes_dim));
  assert(is_size(sca_vec_av, stokes_dim));
  assert( rte_planck_value >= 0 );
  assert( lstep >= 0 );
  assert (!is_singular( ext_mat_av ));

  // Any changes here associated with the extinction matrix should also be
  // implemented in ext2mat.

  // Check, if only the first component of abs_vec is non-zero:
  bool unpol_abs_vec = true;

  for (Index i = 1; i < stokes_dim; i++)
    if (abs_vec_av[i] != 0)
      unpol_abs_vec = false;

  bool unpol_sca_vec = true;

  for (Index i = 1; i < stokes_dim; i++)
    if (sca_vec_av[i] != 0)
      unpol_sca_vec = false;


  //--- Scalar case: ---------------------------------------------------------
  if( stokes_dim == 1 )
    {
      trans_mat(0,0) = exp(-ext_mat_av(0,0) * lstep);
      stokes_vec[0]  = stokes_vec[0] * trans_mat(0,0) +
                       (abs_vec_av[0] * rte_planck_value + sca_vec_av[0]) /
        ext_mat_av(0,0) * (1 - trans_mat(0,0));
    }


  //--- Vector case: ---------------------------------------------------------

  // We have here two cases, diagonal or non-diagonal ext_mat_gas
  // For diagonal ext_mat_gas, we expect abs_vec_gas to only have a
  // non-zero value in position 1.

  //- Unpolarised
  else if( is_diagonal(ext_mat_av) && unpol_abs_vec && unpol_sca_vec )
    {
      trans_mat      = 0;
      trans_mat(0,0) = exp(-ext_mat_av(0,0) * lstep);

      // Stokes dim 1
      //   assert( ext_mat_av(0,0) == abs_vec_av[0] );
      //   Numeric transm = exp( -lstep * abs_vec_av[0] );
      stokes_vec[0] = stokes_vec[0] * trans_mat(0,0) +
                      (abs_vec_av[0] * rte_planck_value + sca_vec_av[0]) /
                      ext_mat_av(0,0) * (1 - trans_mat(0,0));

      // Stokes dims > 1
      for( Index i=1; i<stokes_dim; i++ )
        {
          //      assert( abs_vec_av[i] == 0.);
          trans_mat(i,i) = trans_mat(0,0);
          stokes_vec[i]  = stokes_vec[i] * trans_mat(i,i) +
                       sca_vec_av[i] / ext_mat_av(i,i)  * (1 - trans_mat(i,i));
        }
    }


  //- General case
  else
    {
      //Initialize internal variables:

      // Matrix LU used for LU decompostion and as dummy variable:
      Matrix LU(stokes_dim, stokes_dim);
      ArrayOfIndex indx(stokes_dim); // index for pivoting information
      Vector b(stokes_dim); // dummy variable
      Vector x(stokes_dim); // solution vector for K^(-1)*b
      Matrix I(stokes_dim, stokes_dim);

      Vector B_abs_vec(stokes_dim);
      B_abs_vec = abs_vec_av;
      B_abs_vec *= rte_planck_value;
      
      for (Index i=0; i<stokes_dim; i++)
        b[i] = B_abs_vec[i] + sca_vec_av[i];  // b = abs_vec * B + sca_vec

      // solve K^(-1)*b = x
      ludcmp(LU, indx, ext_mat_av);
      lubacksub(x, LU, b, indx);


      Matrix ext_mat_ds(stokes_dim, stokes_dim);
      ext_mat_ds = ext_mat_av;
      ext_mat_ds *= -lstep; // ext_mat_ds = -ext_mat*ds

      Index q = 10;  // index for the precision of the matrix exp function
      //Matrix exp_ext_mat(stokes_dim, stokes_dim);
      //matrix_exp(exp_ext_mat, ext_mat_ds, q);
      matrix_exp( trans_mat, ext_mat_ds, q);

      Vector term1(stokes_dim);
      Vector term2(stokes_dim);

      id_mat(I);
      for(Index i=0; i<stokes_dim; i++)
        {
          for(Index j=0; j<stokes_dim; j++)
            LU(i,j) = I(i,j) - trans_mat(i,j); // take LU as dummy variable
          // LU(i,j) = I(i,j) - exp_ext_mat(i,j); // take LU as dummy variable
        }
      mult(term2, LU, x); // term2: second term of the solution of the RTE with
                          //fixed scattered field

      // term1: first term of solution of the RTE with fixed scattered field
      //mult(term1, exp_ext_mat, stokes_vec);
      mult( term1, trans_mat, stokes_vec );

      for (Index i=0; i<stokes_dim; i++)
        stokes_vec[i] = term1[i] + term2[i];  // Compute the new Stokes Vector
    }
}



//! surface_calc
/*!
    Weights together downwelling radiation and surface emission.

    *iy* must have correct size when function is called.

    \param   iy                 In/Out: Radiation matrix, amtching 
                                        the WSV with the same name.
    \param   I                  Input: Downwelling radiation, with dimensions
                                       matching: 
                                       (surface_los, f_grid, stokes_dim)
    \param   surface_los        Input: As the WSV with the same name.
    \param   surface_rmatrix    Input: As the WSV with the same name.
    \param   surface_emission   Input: As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2005-04-07
*/
void surface_calc(
              Matrix&         iy,
        ConstTensor3View      I,
        ConstMatrixView       surface_los,
        ConstTensor4View      surface_rmatrix,
        ConstMatrixView       surface_emission )
{
  // Some sizes
  const Index   nf         = I.nrows();
  const Index   stokes_dim = I.ncols();
  const Index   nlos       = surface_los.nrows();

  iy = surface_emission;
  
  // Loop *surface_los*-es. If no such LOS, we are ready.
  if( nlos > 0 )
    {
      for( Index ilos=0; ilos<nlos; ilos++ )
        {
          Vector rtmp(stokes_dim);  // Reflected Stokes vector for 1 frequency

          for( Index iv=0; iv<nf; iv++ )
            {
          mult( rtmp, surface_rmatrix(ilos,iv,joker,joker), I(ilos,iv,joker) );
          iy(iv,joker) += rtmp;
            }
        }
    }
}



//! trans_step_std
/*!
    Solves monochromatic (vector) Beer-Lambert for one step

    The function can be used for clearsky and cloudbox calculations.

    (As rte_step_std, but without emission and scattering source term)

    When calling the function, the vector *stokes_vec* shall contain the
    Stokes vector for the incoming radiation. The function returns this
    vector, then containing the outgoing radiation on the other side of the 
    layer.

    Transmission calculated by *ext2trans*.

    \param   stokes_vec         Input/Output: A Stokes vector.
    \param   trans_mat          Input/Output: Transmission matrix of slab.
    \param   ext_mat            Input: Averaged extinction matrix.
    \param   lstep             Input: The length of the RTE step.

    \author Claudia Emde and Patrick Eriksson, 
    \date   2010-10-15
*/
void trans_step_std(//Output and Input:
              VectorView stokes_vec,
              MatrixView trans_mat,
              //Input
              ConstMatrixView ext_mat_av,
              const Numeric& lstep )
{
  // Checks made in *ext2trans*.
  ext2trans( trans_mat, ext_mat_av, lstep );  

  Vector tmp = stokes_vec;

  mult( stokes_vec, trans_mat, tmp );
}



//! vectorfield2los
/*!
    Calculates the size and direction of a vector field defined as u, v and w
    components.

    \param   l      Size/magnitude of the vector.
    \param   los    Out: The direction, as a LOS vector
    \param   u      Zonal component of the vector field
    \param   v      N-S component of the vector field
    \param   w      Vertical component of the vector field

    \author Patrick Eriksson 
    \date   2012-07-10
*/
void vectorfield2los(
        Numeric&    l,
        Vector&     los,
  const Numeric&    u,
  const Numeric&    v,
  const Numeric&    w )
{
  l= sqrt( u*u + v*v + w*w );
  //
  los.resize(2);
  //
  los[0] = acos( w / l );
  los[1] = atan2( u, v );   
}    



