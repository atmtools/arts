/*-----------------------------------------------------------------------
FILE:      atm_funcs.cc

INCLUDES:  Functions releated to atmospheric physics or geometry.

FUNCTIONS: ztan_geom

HISTORY:   10.04.00 Started by Patrick Eriksson.
-----------------------------------------------------------------------*/

#include "arts.h"
#include "vecmat.h"
#include "messages.h"          
#include "math_funcs.h"          

extern const Numeric EARTH_RADIUS;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;



//==========================================================================
//=== Tangent altitudes.
//==========================================================================

/** Calculates the geometrical tangent altitude.
    That is, refraction is neglected.

    @return        the tangent altitude
    @param view    the angle between zenith and the LOS
    @param z_plat  the platform altitude

    @author Patrick Eriksson 08.04.2000 */
Numeric ztan_geom(
        const Numeric&     view,
        const Numeric&     z_plat )
{
  Numeric  z_tan;
  if ( view >= 90 )   
    z_tan = (EARTH_RADIUS+z_plat)*sin(DEG2RAD*view) - EARTH_RADIUS; 
  else
    z_tan = 9.9999e6;
  return z_tan;
}



//==========================================================================
//=== Conversion between pressures and vertical altitudes.
//==========================================================================

/** Converts an altitude vector to pressures.
    The log. of the pressures are interpolated linearly, 
    i.e. (Matlab notation):

      p = exp(interp1(z0,log(p0),z,'linear'))

    @return        the pressures at z
    @param z0      original altitude grid
    @param p0      original pressure grid
    @param z       new altitude grid

    @author Patrick Eriksson 10.04.2000 */
VECTOR z2p(
        const VECTOR&     z0,
        const VECTOR&     p0,
        const VECTOR&     z )
{
  return exp( interp_lin(z0,log(p0),z) );
}
