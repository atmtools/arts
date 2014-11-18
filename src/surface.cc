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


//! calc_incang
/*!
    Calculates the incidence angle for a flat surface, based on rte_los and
    specular_los.

    \param   incang             Return: Incidence angle.
    \param   rte_los            Input: As the WSV with the same name.
    \param   specular_los       Input: As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2012-11-15
*/
Numeric calc_incang(
   ConstVectorView   rte_los,
   ConstVectorView   specular_los )
{
  return ( 180-abs(rte_los[0]) + abs(specular_los[0]) ) / 2;
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
            const Complex   a = Rh * conj(Rv);
            const Complex   b = Rv * conj(Rh);
            const Numeric   c = real( a + b ) / 2.0;

            surface_rmatrix(2,2) = c;
      
            if( stokes_dim > 3 )
              {
                const Numeric   d = imag( a - b ) / 2.0;

                surface_rmatrix(2,3) = d;
                surface_rmatrix(3,2) = -d;
                surface_rmatrix(3,3) = c;
              }
          }
    }
}
