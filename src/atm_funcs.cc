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
extern const Numeric PLANCK_CONST;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric BOLTZMAN_CONST;



//==========================================================================
//=== Physical functions.
//==========================================================================

void planck (
              MATRIX&     B, 
        const VECTOR&     f,
        const VECTOR&     t )
{
  static const Numeric a = 2.0*PLANCK_CONST/(SPEED_OF_LIGHT*SPEED_OF_LIGHT);
  static const Numeric b = PLANCK_CONST/BOLTZMAN_CONST;
  const size_t  n_f  = f.dim();
  const size_t  n_t  = t.dim();
  size_t        i_f, i_t;
  Numeric       c, d;
  B.newsize(n_f,n_t);
  
  for ( i_f=1; i_f<=n_f; i_f++ )
  {
    c = a * f(i_f)*f(i_f)*f(i_f);
    d = b * f(i_f);
    for ( i_t=1; i_t<=n_t; i_t++ )
      B(i_f,i_t) = c / (exp(d/t(i_t)) - 1.0);
  }
}

void planck (
             VECTOR&    B,
       const VECTOR&    f,
       const Numeric&   t )
{
  static const Numeric a = 2.0*PLANCK_CONST/(SPEED_OF_LIGHT*SPEED_OF_LIGHT);
  static const Numeric b = PLANCK_CONST/BOLTZMAN_CONST;
  
  B = ediv( a*emult(f,emult(f,f)), exp(f*(b/t))-1.0 ) ;
}



//==========================================================================
//=== Tangent altitudes.
//==========================================================================

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
//=== Conversion and interpolation of pressure and altitude grids.
//==========================================================================

void z2p(
              VECTOR&     p,
        const VECTOR&     z0,
        const VECTOR&     p0,
        const VECTOR&     z )
{
  p = exp( interp_lin( z0, log(p0), z ) );
}


void interpp(
              VECTOR&     x, 
        const VECTOR&     p0,
        const VECTOR&     x0,
        const VECTOR&     p )
{
  x = interp_lin( log(p0), x0, log(p) );
}


//==========================================================================
//=== RTE_ITERATE
//==========================================================================

void rte_iterate (
             VECTOR&   y,
       const int&      start_index,
       const int&      stop_index,
       const MATRIX&   Tr,
       const MATRIX&   S,
       const VECTOR&   y_space,
       const int&      ground,
       const VECTOR&   e_ground,
       const VECTOR&   y_ground )
{
  const int   n_f = Tr.dim(1);               // number of frequencies
        int   i_f, i_z;                      // indicies
        int   i_start, i_step, i_break;      // loop variables

  // If LOS starts at ground, init. with Y_GROUND  
  if ( start_index == 1 )
    y = y_ground;

  // If LOS starts in space, init with Y_SPACE
  else
    y = y_space;

  // Check if LOS inside the atmosphere (if START_INDEX=0, the result is Y_SPACE)
  if ( start_index > 0 )
  {
    // Set up index and make first loop. Note that the transmission is 
    // given between the points
    //
    if ( start_index == 1 )
    {
      i_start = 1;
      i_step  = 1;
      i_break = stop_index;
    }
    else 
    {       
      i_start = start_index - 1;         
      i_step  = -1;
      i_break = 0;       
    }
    for ( i_z=i_start; i_z!=i_break; i_z+=i_step ) 
    {
      for ( i_f=1; i_f<=n_f; i_f++ )    
        y(i_f) = y(i_f)*Tr(i_f,i_z) + S(i_f,i_z) * ( 1.0-Tr(i_f,i_z) );
    }
 
    // For limb sounding, also loop upwards. Include ground reflection
    // if necessary.
    //
    if ( (i_step == -1) && (stop_index > 1) )
    {
      if ( ground )               // If intersection with the ground, include 
      {                           // ground reflection
        for ( i_f=1; i_f<=n_f; i_f++ )    
          y(i_f) = y(i_f) * (1.0-e_ground(i_f)) + y_ground(i_f) * e_ground(i_f);
      }
      for ( i_z=1; i_z<stop_index; i_z++ ) // Loop upwards
      {
        for ( i_f=1; i_f<=n_f; i_f++ )    
          y(i_f) = y(i_f)*Tr(i_f,i_z) + S(i_f,i_z) * ( 1.0-Tr(i_f,i_z) );
      }
    } // second part
  } // if any values
}
