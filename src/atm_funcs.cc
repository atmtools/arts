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

// PLANCK (matrix version)
//
// Patrick Eriksson 08.04.00

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



// PLANCK (vector version)
//
// Patrick Eriksson 08.04.00

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

// ZTAN_GEOM
//
// Patrick Eriksson 08.04.00

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
//=== Core functions for RTE and BL 
//==========================================================================

// RTE_ITERATE
//
// Patrick Eriksson 15.06.00

void rte_iterate (
             VECTOR&   y,
       const int&      start_index,
       const int&      stop_index,
       const MATRIX&   Tr,
       const MATRIX&   S,
       const size_t    n_f )
{
        size_t   i_f;        // frequency index
           int   i_z;        // LOS index
           int   i_step;     // step order, -1 or 1

  if ( start_index >= stop_index )
    i_step = -1;
  else
    i_step = 1;

  for ( i_z=start_index; i_z!=(stop_index+i_step); i_z+=i_step ) 
  {
    for ( i_f=1; i_f<=n_f; i_f++ )    
      y(i_f) = y(i_f)*Tr(i_f,i_z) + S(i_f,i_z) * ( 1.0-Tr(i_f,i_z) );
  }
}



// RTE
//
// Patrick Eriksson 22.05.00

void rte (
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
        int   i_f;                           // frequency index
        int   i_break;                       // break index for looping
        int   i_start;                       // variable for second loop

  // If LOS starts at ground, init. with Y_GROUND  
  if ( start_index == 1 )
    y = y_ground;

  // If LOS starts in space, init with Y_SPACE
  else
    y = y_space;

  // Check if LOS inside the atmosphere (if START_INDEX=0, Y=Y_SPACE)
  if ( start_index > 0 )
  {
    // Determine break index for looping
    if ( ground > 0 )
      i_break = ground;
    else
      i_break = 1;       

    // Make first loop
    rte_iterate( y, start_index-1, i_break, Tr, S, n_f );

    // We are now at the sensor, the ground or the tangent point
    // We are ready only if we are at the sensor.
    // If at sensor, we have that STOP_INDEX=1 and not GROUND
    if ( !(stop_index==1 && !ground) )
    {
      // Set most common values for I_START and I_BREAK
      i_start = 1;
      i_break = stop_index - 1;
      
      // If at the ground, include ground reflection. 
      // The loop can continue both downwards or upwards
      if ( ground )
      {            
        for ( i_f=1; i_f<=n_f; i_f++ )    
          y(i_f) = y(i_f)*(1.0-e_ground(i_f)) + y_ground(i_f)*e_ground(i_f);
        
        if ( ground != 1 )  // 2D case, loop downwards
	{
         i_start = ground - 1;
         i_break = 1;
        }
      }

      // Make second loop
      rte_iterate( y, i_start, i_break, Tr, S, n_f );

    } // second part
  } // if any values
}



// BL_ITERATE
//
// Patrick Eriksson 15.06.00

void bl_iterate (
             VECTOR&   y,
       const int&      start_index,
       const int&      stop_index,
       const MATRIX&   Tr,
       const size_t    n_f )
{
        size_t   i_f;        // frequency index
           int   i_z;        // LOS index
           int   i_step;     // step order, -1 or 1

  if ( start_index >= stop_index )
    i_step = -1;
  else
    i_step = 1;

  for ( i_z=start_index; i_z!=(stop_index+i_step); i_z+=i_step ) 
  {
    for ( i_f=1; i_f<=n_f; i_f++ )    
      y(i_f) *= Tr(i_f,i_z);
  }
}



// BL
//
// Patrick Eriksson 22.05.00

void bl (
             VECTOR&   y,
       const int&      start_index,
       const int&      stop_index,
       const MATRIX&   Tr,
       const int&      ground,
       const VECTOR&   e_ground )
{
  const int   nf = Tr.dim(1);          // number of frequencies
  //        int   j;                       // LOS index   
        int   iy;                      // frequency index

  // Init Y
  y.newsize(nf);
  y = 1;

  // Loop steps passed twice
  if ( stop_index > 1 )
  {
    bl_iterate( y, 1, stop_index-1, Tr, nf );
    y = emult( y, y );  
  }

  // Loop remaining steps
  if ( start_index != stop_index )
  {
    bl_iterate( y, stop_index, start_index-1, Tr, nf );
  }

  // Include effect of ground reflection
  if ( ground )
  {
    for ( iy=1; iy<=nf; iy++ )    
      y(iy) *= ( 1.0 - e_ground(iy) );
  }
}



//==========================================================================
//=== Conversion and interpolation of pressure and altitude grids.
//==========================================================================

// Z2P
//
// Patrick Eriksson 10.04.00

void z2p(
              VECTOR&     p,
        const VECTOR&     z0,
        const VECTOR&     p0,
        const VECTOR&     z )
{
  p = exp( interp_lin( z0, log(p0), z ) );
}



// INTERPP (vector version)
//
// Patrick Eriksson 12.04.00

void interpp(
              VECTOR&     x, 
        const VECTOR&     p0,
        const VECTOR&     x0,
        const VECTOR&     p )
{
  interp_lin( x, log(p0), x0, log(p) );
}



// INTERPP (matrix version)
//
// Patrick Eriksson 13.06.00

void interpp(
              MATRIX&  A,
        const VECTOR&  p0, 
        const MATRIX&  A0, 
        const VECTOR&  p )
{
  interp_lin_row( A, log(p0), A0, log(p) );
}


