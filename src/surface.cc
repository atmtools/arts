/* Copyright (C) 2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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
   \file   surface.cc
   \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   \date   2012-02-06 

   This file contains internal functions associated with the surface.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include "auto_md.h"
#include "check_input.h"          
#include "complex.h"          
#include "geodetic.h"          
#include "matpackI.h"
#include "math_funcs.h"
#include "physics_funcs.h"
#include "workspace_ng.h"
#include "surface.h"





/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! surface_specular_los
/*!
    Calculates the specular direction, including the effect of surface 
    topograpghy

    \param   specular_los       Out: LOS corresponding to specular direction
    \param   rte_pos            Input: As the WSV with the same name.
    \param   rte_los            Input: As the WSV with the same name.
    \param   atmosphere_dim     Input: As the WSV with the same name.
    \param   lat_grid           Input: As the WSV with the same name.
    \param   lon_grid           Input: As the WSV with the same name.
    \param   refellipsoid       Input: As the WSV with the same name.
    \param   z_surface          Input: As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2012-07-18
*/
void surface_specular_los(
         Vector&     specular_los,
   ConstVectorView   rte_pos,
   ConstVectorView   rte_los,
   const Index&      atmosphere_dim,
   ConstVectorView   lat_grid,
   ConstVectorView   lon_grid,
   ConstVectorView   refellipsoid,
   ConstMatrixView   z_surface )
{
  specular_los.resize( max( Index(1), atmosphere_dim-1 ) );

  if( atmosphere_dim == 1 )
    { specular_los[0] = 180 - rte_los[0]; }

  else if( atmosphere_dim == 2 )
    { 
      chk_interpolation_grids( "Latitude interpolation", lat_grid, rte_pos[1] );
      GridPos gp_lat;
      gridpos( gp_lat, lat_grid, rte_pos[1] );
      Numeric c1;         // Radial slope of the surface
      plevel_slope_2d( c1, lat_grid, refellipsoid, z_surface(joker,0), 
                                                          gp_lat, rte_los[0] );
      Vector itw(2); interpweights( itw, gp_lat );
      const Numeric rv_surface = refell2d( refellipsoid, lat_grid, gp_lat )
                                + interp( itw, z_surface(joker,0), gp_lat );
      specular_los[0] = sign( rte_los[0] ) * 180 - rte_los[0] - 
                                           2*plevel_angletilt( rv_surface, c1 );
    }

  else if( atmosphere_dim == 3 )
    { 
      // Calculate surface normal in South-North direction
      chk_interpolation_grids( "Latitude interpolation", lat_grid, rte_pos[1] );
      chk_interpolation_grids( "Longitude interpolation", lon_grid, rte_pos[2]);
      GridPos gp_lat, gp_lon;
      gridpos( gp_lat, lat_grid, rte_pos[1] );
      gridpos( gp_lon, lon_grid, rte_pos[2] );
      Numeric c1, c2;
      plevel_slope_3d( c1, c2, lat_grid, lon_grid, refellipsoid, z_surface, 
                       gp_lat, gp_lon, 0 );
      Vector itw(4); interpweights( itw, gp_lat, gp_lon );
      const Numeric rv_surface = refell2d( refellipsoid, lat_grid, gp_lat )
                                 + interp( itw, z_surface, gp_lat, gp_lon );
      const Numeric zaSN = 90 - plevel_angletilt( rv_surface, c1 );
      // The same for East-West
      plevel_slope_3d( c1, c2, lat_grid, lon_grid, refellipsoid, z_surface, 
                       gp_lat, gp_lon, 90 );
      const Numeric zaEW = 90 - plevel_angletilt( rv_surface, c1 );
      // Convert to Cartesian, and determine normal by cross-product
      Vector tangentSN(3), tangentEW(3), normal(3);
      zaaa2cart( tangentSN[0], tangentSN[1], tangentSN[2], zaSN, 0 );
      zaaa2cart( tangentEW[0], tangentEW[1], tangentEW[2], zaEW, 90 );
      cross3( normal, tangentSN, tangentEW );
      // Convert rte_los to cartesian and flip direction
      Vector di(3);
      zaaa2cart( di[0], di[1], di[2], rte_los[0], rte_los[1] );
      di *= -1;
      // Specular direction is 2(dn*di)dn-di, where dn is the normal vector
      Vector speccart(3);      
      const Numeric fac = 2 * (normal * di);
      for( Index i=0; i<3; i++ )
        { speccart[i] = fac*normal[i] - di[i]; }
      cart2zaaa( specular_los[0], specular_los[1], speccart[0], speccart[1], 
                                                                speccart[2] );
   }
}


//! surface_specular_R_and_b
/*!
    Sets up the surface reflection matrix and emission vector for
    the case of specular reflection.

    The function handles only one frequency at the time.

    See further the surface chapter in the user guide.

    \param   surface_rmatrix   Out: As the WSV with the same name, but slice
                                    for one direction and one frequency.
    \param   surface_emission  Out: As the WSV with the same name, but slice
                                    for one direction and one frequency.
    \param   Rv                In: Complex amplitude relection coefficient
                                   for vertical polarisation.
    \param   Rh                In: Complex amplitude relection coefficient
                                   for horisontal polarisation.
    \param   f                 In: Frequency (a scalar).
    \param   stokes_dim        In: As the WSV with the same name.
    \param   blackbody_radiation_agenda In: As the WSV with the same name.
    \param   surface_skin_t    In: As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2004-09-24
*/
void surface_specular_R_and_b(
              Workspace&   ws,
              MatrixView   surface_rmatrix,
              VectorView   surface_emission,
        const Complex&     Rv,
        const Complex&     Rh,
        const Numeric&     f,
        const Index&       stokes_dim,
        const Numeric&     surface_skin_t,
        const Agenda&      blackbody_radiation_agenda )
{
  assert( surface_rmatrix.nrows() == stokes_dim );
  assert( surface_rmatrix.ncols() == stokes_dim );
  assert( surface_emission.nelem() == stokes_dim );

  // Expressions are derived in the surface chapter in the user guide

  surface_rmatrix   = 0.0;
  surface_emission  = 0.0;

  Vector   B;
  blackbody_radiation_agendaExecute( ws, B, surface_skin_t, Vector(1,f), 
                                     blackbody_radiation_agenda ); 

  const Numeric   rv    = pow( abs(Rv), 2.0 );
  const Numeric   rh    = pow( abs(Rh), 2.0 );
  const Numeric   rmean = ( rv + rh ) / 2;

  surface_rmatrix(0,0) = rmean;
  surface_emission[0]  = B[0] * ( 1 - rmean );

  if( stokes_dim > 1 )
    {
      const Numeric   rdiff = ( rv - rh ) / 2;

      surface_rmatrix(1,0) = rdiff;
      surface_rmatrix(0,1) = rdiff;
      surface_rmatrix(1,1) = rmean;
      surface_emission[1]  = -B[0] * rdiff;

        if( stokes_dim > 2 )
          {
            const Complex   a     = Rh * conj(Rv);
            const Complex   b     = Rv * conj(Rh);
            const Numeric   c     = real( a + b ) / 2.0;

            surface_rmatrix(2,2) = c;
      
            if( stokes_dim > 3 )
              {
                const Numeric   d     = imag( a - b ) / 2.0;

                surface_rmatrix(3,2) = d;
                surface_rmatrix(2,3) = d;
                surface_rmatrix(3,3) = c;
              }
          }
    }
}
