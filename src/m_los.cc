/*-----------------------------------------------------------------------
FILE:      m_los.cc

INCLUDES:  This file includes functions for a 1D atmosphere to determine:
            1. the line of sight (LOS)
            2. variables along LOS needed for the radiative transfer 
               equation (RTE)

           !PE! Check this text when ready 
         
           LOS is assumed to be symmetric around either the tangent
           point, or the point of ground reflection. Only points/data
           are stored for half of the LOS.  The first point of LOS is
           placed at the tangent point, the ground or the platform,
           i.e. data are stored with increasing altitude.  The LOS is
           defined in equal long geometrical steps along the path.
           The LOS is defined by a structere of type Los (see
           los_1d.h).  See also the function rte in physical_func.cc.

           The functions to calculate LOS follow to large extent the
           algorithms developed for Odin (i.e. the Skuld forward
           model). The equation numbers below refer to the description
           of Skuld:

           On simulating passive observations of the middle atmosphere
           in the range 1-1000 GHz P.Eriksson and F. Merino Chalmers
           University, 1997.

           also found in (Paper H)

           Microwave radiometric observations of the middle
           atmosphere: Simulations and inversions Patrick Eriksson,
           Ph.D. thesis, Chalmers University, 1999.

FUNCTIONS: 

HISTORY:   27.06.99 Created by Patrick Eriksson.
-----------------------------------------------------------------------*/

#include "arts.h"
#include "vecmat.h"
#include "messages.h"          
#include "workspace.h"          
#include "math_funcs.h"          
#include "atm_funcs.h"          

extern const Numeric EARTH_RADIUS;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;



//==========================================================================
//=== Calculation help functions
//==========================================================================

// 
// A help function to calculate LOS for upward observations
// No refraction.
//
VECTOR upward_geom(
        const Numeric&      z_plat,
        const Numeric&      view,
        const Numeric&      l_step,
        const Numeric&      z_abs_max )
{
  if ( view > 90 )
    throw logic_error("Upward function used for viewing angle > 90 deg."); 

  Numeric  a, b, llim;
  VECTOR   l;             // distance from sensor

  a     = EARTH_RADIUS + z_abs_max;
  b     = (EARTH_RADIUS + z_plat)*sin(DEG2RAD*view);
  llim  = sqrt(a*a-b*b) - (EARTH_RADIUS+z_plat)*cos(DEG2RAD*view) ;
  if ( llim < l_step )  
    return VECTOR(0);               // If only one point, return empty vector
  l = linspace( 0, llim, l_step );
  b     = EARTH_RADIUS + z_plat;  
  return ( sqrt(b*b+l*l+(2.0*b*cos(DEG2RAD*view))*l) - EARTH_RADIUS );
}


// 
// A function to calculate LOS for observations from space
// No refraction.
//
VECTOR space_geom(
         const Numeric&    z_tan,
         const Numeric&    l_step,
         const Numeric&    z_abs_max,
         const Numeric&    z_ground )
{
  Numeric   a, b, llim;
  VECTOR    l;              // length from the tangent point

  // If LOS outside the atmosphere, return empty vector
  if ( z_tan >= z_abs_max )
    return VECTOR(0);

  // Only through the atmosphere
  else if ( z_tan >= z_ground )
  {
    a    = EARTH_RADIUS + z_abs_max;
    b    = EARTH_RADIUS + z_tan;
    llim = sqrt( a*a - b*b );        
    if ( llim >= 2*l_step )           // There must be at least 2 points
    {
      linspace( l, 0, llim, l_step );
      return ( sqrt(b*b+l*l) - EARTH_RADIUS ); 
    }
    else
      return VECTOR(0);               // If only 1 point, return empty vector
  }   

  // Intersection with the ground
  else
  {
    // Determine the "viewing angle" at ground level and call upward function 
    Numeric view = RAD2DEG*asin((EARTH_RADIUS+z_tan)/(EARTH_RADIUS+z_ground));
    return upward_geom( z_ground, view, l_step, z_abs_max );    
  }
}



//==========================================================================
//=== Sub-functions to losGeneral
//==========================================================================

//=== Observation from space, limb or nadir
void los_no_refr_space(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     view,
              const Numeric&    l_step,
              const Numeric&    z_abs_max,
              const Numeric&    z_ground )
{
  Numeric z_tan;

  // Get the number of viewing angles
  int n  =  view.size();

  // Set up the l_step, ground, start and stop vectors
  los.p.newsize(n);
  los.l_step.newsize(n);
  los.ground.newsize(n);
  los.start.newsize(n);
  los.stop.newsize(n);

  // The step length is here throughout l_step
  los.l_step = l_step;

  // Loop the viewing angles
  for ( int i=1; i<=n; i++ )
  { 
    z_tan = ztan_geom(view(i),z_plat);
    los.p(i) = space_geom( z_tan, l_step, z_abs_max, z_ground );
    los.start(i)  = los.p(i).size();
    los.stop(i)   = los.start(i);
    if ( z_tan >= z_ground )
      los.ground(i) = 0;
    else
      los.ground(i) = 1;
  }
}


//=== Observation from point inside the atmosphere, upward or downward
void los_no_refr_inside(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     view,
              const Numeric&    l_step,
              const Numeric&    z_abs_max,
              const Numeric&    z_ground )
{
  Numeric z_tan, a, b, c, l1;

  // Get the number of viewing angles
  int n  =  view.size();

  // Set up the l_step, ground, start and stop vectors
  los.p.newsize(n);
  los.l_step.newsize(n);
  los.ground.newsize(n);
  los.start.newsize(n);
  los.stop.newsize(n);

  // Loop the viewing angles
  for ( int i=1; i<=n; i++ )
  { 
    // Upward
    if ( view(i) <= 90 )
    {
      los.p(i) = upward_geom( z_plat, view(i), l_step, z_abs_max );    
      los.l_step(i) = l_step;
      los.start(i)  = los.p(i).size();
      los.stop(i)   = 1;
      los.ground(i) = 0;
    }

    // Downward
    else
    {
      z_tan = ztan_geom(view(i),z_plat);

      // Only through the atmosphere
      if ( z_tan >= z_ground)
      {
	a      = EARTH_RADIUS + z_plat;
	b      = EARTH_RADIUS + z_tan;
	l1     = sqrt(a*a-b*b);
	los.stop(i)   = (int) ceil( l1 / l_step + 1.0 );  
	los.l_step(i) = l1 / ( (Numeric)los.stop(i) - 1.0 );
	los.p(i) = space_geom( z_tan, los.l_step(i), z_abs_max, z_ground );
        los.start(i)  = los.p(i).size();
        los.ground(i) = 0;
      }   
    
      // Intersection with the ground
      else
      {
        a      = EARTH_RADIUS + z_ground;
	b      = EARTH_RADIUS + z_plat;
	c      = EARTH_RADIUS + z_tan;
	l1     = sqrt(b*b-c*c) - sqrt(a*a-c*c);
	if ( l1 < l_step )  
          throw runtime_error("The distance between platform and ground is to small for downward looking");
	los.stop(i)   = (int) ceil( l1 / l_step + 1.0 );
	los.l_step(i) = l1 / ( (double)los.stop(i) - 1.0 );
	los.p(i) = upward_geom( z_ground, RAD2DEG*asin(c/a), los.l_step(i), 
                                                                  z_abs_max);
        los.start(i)  = los.p(i).size();
        los.ground(i) = 1;        
      }
    }
  }
}



//==========================================================================
//=== The methods
//==========================================================================

/** To be written !PE!

    @param xx
    @param xx

    @author Patrick Eriksson */
void losGeneral(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     view,
              const Numeric&    l_step,
              const VECTOR&     p_abs,
              const VECTOR&     z_abs,
              const int&    refr,
              const Numeric&    l_step_refr,
              const Numeric&    z_ground,
              const VECTOR&     e_ground )

{     
  // Some checks                                                      
  if ( z_ground < z_abs(1) )
    throw runtime_error(
      "There is a gap between the ground and the lowest absorption altitudes");
  if ( z_plat < z_ground )
    throw runtime_error("Your platform altitude is below the ground");
  if ( z_plat < z_abs(1) )  
    throw runtime_error(
      "The platform cannot be below the lowest absorption altitude");

  // Get highest absorption altitude
  Numeric z_abs_max  =  max(z_abs);

  // The functions below calculates the vertical altitudes along LOS, 
  // stored temporarily in los.p.
  // The variables start, stop and ground are given their values,
  // assuming that the ground has an emission < 1.

  // Without refraction
  if ( !refr)
  {
    if ( z_plat >= z_abs_max )
      los_no_refr_space(los,z_plat,view,l_step,z_abs_max,z_ground);
    else
      los_no_refr_inside(los,z_plat,view,l_step,z_abs_max,z_ground);
  } 

  // With refraction
  else
  {
  }

  // Convert altitudes to pressures.
  for ( size_t i=1; i<=view.size(); i++ )
    los.p(i) = z2p( z_abs, p_abs, los.p(i) );

  // If all ground emissions are 1 (>0.999), i.e. the ground acts as a blackbody,
  // the start value is set to 1 if there is a ground reflection.
  if ( any(los.ground) )
  {
    size_t i;
    for ( i=1; (i<=e_ground.size()) && (e_ground(i)>0.999); i++ )
      {}
    if ( i > e_ground.size() )
    {
      for ( i=1; i<=view.size(); i++ )
      {
        if ( los.ground(i) )
	{
          los.start(i) = 1;
          // !PE! remove unnecessary part of los.p
        }
      }  
    }
  }
  out3 << los.p(1);
  out3 << los.start(1) << "\n";
  out3 << los.stop(1) << "\n";
  out3 << los.ground(1) << "\n";
  out3 << los.l_step(1) << "\n";
}
