/*-----------------------------------------------------------------------
FILE:      m_los.cc

INCLUDES:  This file includes functions for a 1D atmosphere to determine:
            1. the line of sight (LOS)
            2. variables along LOS needed for the radiative transfer 
               equation (RTE)
            3. Integration of the RTE

           The LOS is defined by a structure of type Los, defined in
           workspace.h. The first point of LOS (index 1) is here
           placed at the tangent point, the ground or the platform,
           that is, data are stored with increasing altitudes. When
           applicable, LOS is assumed to be symmetric around either
           the tangent point, or the point of ground reflection.
           Points/data are for such cases stored only for half of the
           LOS. This in order to save memory and computational time.
           The user defined step length is applied besides for
           downward observations inside the atmosphere where the step
           length is adjusted to get an integer number of steps
           between the sensor and the tangent aor ground point. The
           new step length is made smaller than the user defined one.

           The algorithms used here are described in the ARTS user guide.

FUNCTIONS: !!

HISTORY:   27.06.99 Created by Patrick Eriksson.
-----------------------------------------------------------------------*/


#include "arts.h"
#include "vecmat.h"
#include "messages.h"          
#include "workspace.h"          
#include "math_funcs.h"          
#include "atm_funcs.h"          
#include "los.h"


extern const Numeric EARTH_RADIUS;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric COSMIC_BG_TEMP;
extern const Numeric SUN_TEMP;



//==========================================================================
//=== Calculation help functions 
//==========================================================================

// 
// A help function to calculate LOS for upward observations.
// No refraction.
//
void upward_geom(
              VECTOR&       z,
        const Numeric&      z_plat,
        const Numeric&      view,
        const Numeric&      l_step,
        const Numeric&      z_abs_max )
{
  Numeric  a, b;          // temporary values
  Numeric  llim;          // distance to atmospheric limit
  VECTOR   l;             // distance from sensor

  if ( view > 90 )
    throw logic_error("Upward function used for viewing angle > 90 deg."); 

  a     = EARTH_RADIUS + z_abs_max;
  b     = (EARTH_RADIUS + z_plat)*sin(DEG2RAD*view);
  llim  = sqrt(a*a-b*b) - (EARTH_RADIUS+z_plat)*cos(DEG2RAD*view) ;
  if ( llim < l_step )  
    z.newsize(0);         // If only one point, return empty vector
  else
  {
    l = linspace( 0, llim, l_step );
    b = EARTH_RADIUS + z_plat;  
    z = sqrt(b*b+emult(l,l)+(2.0*b*cos(DEG2RAD*view))*l) - EARTH_RADIUS;
  }
}



// 
// A function to calculate LOS for observations from space.
// No refraction.
//
void space_geom(
               VECTOR&     z,
         const Numeric&    z_tan,
         const Numeric&    l_step,
         const Numeric&    z_abs_max,
         const Numeric&    z_ground )
{
  Numeric  a, b;          // temporary values
  Numeric  llim;          // distance to atmospheric limit
  VECTOR   l;             // length from the tangent point

  // If LOS outside the atmosphere, return empty vector
  if ( z_tan >= z_abs_max )
    z.newsize(0);

  // Only through the atmosphere
  else if ( z_tan >= z_ground )
  {
    a    = EARTH_RADIUS + z_abs_max;
    b    = EARTH_RADIUS + z_tan;
    llim = sqrt( a*a - b*b );        
    if ( llim >= 2*l_step )       // There must be at least 2 points
    {
      linspace( l, 0, llim, l_step );
      z = sqrt(b*b+emult(l,l)) - EARTH_RADIUS; 
    }
    else
      z.newsize(0);               // If only 1 point, return empty vector
  }   

  // Intersection with the ground
  else
  {
    // Determine the "viewing angle" at ground level and call upward function 
    Numeric view = RAD2DEG*asin((EARTH_RADIUS+z_tan)/(EARTH_RADIUS+z_ground));
    upward_geom( z, z_ground, view, l_step, z_abs_max );    
  }
}



//
//=== Sub-functions to losBasic =============================================
//
// Observation from space, limb or nadir. 
// No refraction.
//
void los_no_refr_space(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     view,
              const Numeric&    l_step,
              const Numeric&    z_abs_max,
              const Numeric&    z_ground )
{
  Numeric z_tan;           // the tangent altitude
  int     n = view.dim();  // the number of viewing angles

  // The step length is here throughout l_step
  los.l_step = l_step;

  // Loop the viewing angles
  for ( int i=1; i<=n; i++ )
  { 
    z_tan = ztan_geom( view(i), z_plat);
    space_geom( los.p(i), z_tan, l_step, z_abs_max, z_ground );
    los.start(i)  = los.p(i).dim();
    los.stop(i)   = los.start(i);
    if ( z_tan >= z_ground )
      los.ground(i) = 0;          // no ground intersection
    else
      los.ground(i) = 1;          // ground at index 1
  }
}

//
// Observation from point inside the atmosphere, upward or downward.
// No refraction.
//
void los_no_refr_inside(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     view,
              const Numeric&    l_step,
              const Numeric&    z_abs_max,
              const Numeric&    z_ground )
{
  Numeric z_tan, a, b, c, l1;  // see below
  int     n = view.dim();      // the number of viewing angles

  // Loop the viewing angles
  for ( int i=1; i<=n; i++ )
  { 
    // Upward
    if ( view(i) <= 90 )
    {
      upward_geom( los.p(i), z_plat, view(i), l_step, z_abs_max );    
      los.l_step(i) = l_step;           // l_step not changed
      los.start(i)  = los.p(i).dim();   // length of LOS
      los.stop(i)   = 1;                // stop index here always 1
      los.ground(i) = 0;                // no ground intersection
    }

    // Downward
    else
    {
      // Calculate the tangent altitude
      z_tan = ztan_geom(view(i),z_plat);

      // Only through the atmosphere
      if ( z_tan >= z_ground )
      {
	a      = EARTH_RADIUS + z_plat;   // help variable
	b      = EARTH_RADIUS + z_tan;    // help variable
	l1     = sqrt(a*a-b*b);           // distance platform-tangent point
        // Adjust l_step downwards to get an integer number of steps
	los.stop(i)   = (int) ceil( l1 / l_step + 1.0 );  
	los.l_step(i) = l1 / ( (Numeric)los.stop(i) - 1.0 );
	space_geom( los.p(i), z_tan, los.l_step(i), z_abs_max, z_ground );
        los.start(i)  = los.p(i).dim();
        los.ground(i) = 0;                // no gound intersection
      }   
    
      // Intersection with the ground
      else
      {
        a      = EARTH_RADIUS + z_ground;
	b      = EARTH_RADIUS + z_plat;
	c      = EARTH_RADIUS + z_tan;
	l1     = sqrt(b*b-c*c) - sqrt(a*a-c*c); // distance platform-ground
        // Adjust l_step downwards to get an integer number of steps
	los.stop(i)   = 1 + (int) ceil( l1 / l_step );
	los.l_step(i) = l1 / ( (double)los.stop(i) - 1.0 );
	upward_geom(los.p(i),z_ground,RAD2DEG*asin(c/a),los.l_step(i),z_abs_max);
        los.start(i)  = los.p(i).dim();
        los.ground(i) = 1;                // ground at index 1
      }
    }
  }
}
//=== End losGeneral =======================================================



//
// rte_iterate
//
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
        int   i_break;                       // break index for looping
        int   i_start, i_step;               // variables for second loop

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
      i_break = ground - 1;
    else
      i_break = 0;       

    // Make first loop
    for ( i_z=(start_index-1); i_z!=i_break; i_z-=1 ) 
    {
      for ( i_f=1; i_f<=n_f; i_f++ )    
        y(i_f) = y(i_f) * Tr(i_f,i_z) + S(i_f,i_z) * ( 1.0-Tr(i_f,i_z) );
    }

    // We are now at the sensor, the ground or the tangent point
    // We are ready only if we are at the sensor.
    // If at sensor, we have that STOP_INDEX=1 and not GROUND
    if ( !(stop_index==1 && !ground) )
    {
      // If at the ground, include ground reflection. 
      // The loop can continue both downwards or upwards
      if ( ground )
      {            
        for ( i_f=1; i_f<=n_f; i_f++ )    
          y(i_f) = y(i_f)*(1.0-e_ground(i_f)) + y_ground(i_f)*e_ground(i_f);
        
        if ( ground == 1 )  // 1D case, loop upwards
        {
          i_start = 1;
          i_step  = 1;
          i_break = stop_index;
	} 
        else                // 2D case, loop downwards
	{
         i_start = ground - 1;
         i_step  = -1;
         i_break = 0;
        }
      }

      // We are at a tangent point, loop upwards
      else
      {
          i_start = 1;
          i_step  = 1;
          i_break = stop_index;        
      }

      // Make second loop
      for ( i_z=i_start; i_z!=i_break; i_z+=i_step ) 
      {
        for ( i_f=1; i_f<=n_f; i_f++ )    
          y(i_f) = y(i_f)*Tr(i_f,i_z) + S(i_f,i_z) * ( 1.0-Tr(i_f,i_z) );
      }
    } // second part
  } // if any values
}



//==========================================================================
//=== The methods
//==========================================================================

void los1d(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     view,
              const Numeric&    l_step,
              const VECTOR&     p_abs,
              const VECTOR&     z_abs,
              const int&        refr,
              const Numeric&    l_step_refr,
              const VECTOR&     refr_index,
              const Numeric&    z_ground )
{     
  size_t n = view.dim();  // number of viewing angles

  // Some checks                                                      
  if ( z_ground < z_abs(1) )
    throw runtime_error(
      "There is a gap between the ground and the lowest absorption altitude");
  if ( z_plat < z_ground )
    throw runtime_error("Your platform altitude is below the ground");
  if ( z_plat < z_abs(1) )  
    throw runtime_error(
      "The platform cannot be below the lowest absorption altitude");

  // Reallocate the l_step, ground, start and stop vectors
  los.p.newsize(n);
  los.l_step.newsize(n);
  los.ground.newsize(n);
  los.start.newsize(n);
  los.stop.newsize(n);

  // Get highest absorption altitude
  Numeric z_abs_max  =  last(z_abs);

  // Without refraction
  if ( refr == 0 )
  {
    // The functions below calculates the vertical altitudes along LOS, 
    // stored temporarily in los.p.
    if ( z_plat >= z_abs_max )
      los_no_refr_space( los, z_plat, view, l_step, z_abs_max, z_ground );
    else
      los_no_refr_inside( los, z_plat, view, l_step, z_abs_max, z_ground );

    // Convert altitudes to pressures.
    for ( size_t i=1; i<=n; i++ )
      z2p( los.p(i), z_abs, p_abs, los.p(i) );
  } 

  // With refraction
  else if ( refr == 1 )
  {
  }

  else
    throw runtime_error("The refraction flag can only be 0 or 1.");
}



void los1dNoRefraction(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     view,
              const Numeric&    l_step,
              const VECTOR&     p_abs,
              const VECTOR&     z_abs,
              const Numeric&    z_ground )
{
  los1d( los, z_plat, view, l_step, p_abs, z_abs, 0, 0, VECTOR(0) ,z_ground); 
}



void los1dUpward(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     view,
              const Numeric&    l_step,
              const VECTOR&     p_abs,
              const VECTOR&     z_abs )
{
  if ( max(view) > 90 )
    throw runtime_error("At least one viewing angle > 90 degrees, that is, not upwards.");

  los1d( los, z_plat, view, l_step, p_abs, z_abs, 0, 0, VECTOR(0) ,z_plat); 
}


void source1d(
                    ARRAYofMATRIX&   source,
              const Los&             los,   
              const VECTOR&          p_abs,
              const VECTOR&          t_abs,
              const VECTOR&          f_abs )
{     
        VECTOR   t1, t2;             // temperatures along the LOS
  const size_t   n=los.start.dim();  // the number of viewing angles  
        size_t   j, n1;              // n1 is the number of pressure points
 
  // Resize the source array
  source.newsize(n);

  // Loop the viewing angles and:
  //  1. interpolate the temperature
  //  2. take the mean of neighbouring values
  //  3. calculate the Planck function
  for (size_t i=1; i<=n; i++ ) 
  {
    if ( los.p(i).size() > 0 )
    {
      interpp( t1, p_abs, t_abs, los.p(i) );
      n1 = t1.dim();
      t2.newsize(n1-1);
      for ( j=1; j<n1; j++ )
        t2(j) = ( t1(j) + t1(j+1) ) / 2.0;
      planck( source(i), f_abs, t2 );
    }
  }  
}



void trans1d(
                    ARRAYofMATRIX&   trans,
              const Los&             los,   
              const VECTOR&          p_abs,
              const MATRIX&          abs )
{    
  // Some variables
  const size_t   n = los.start.dim();// the number of viewing angles
  const size_t   nf = abs.dim(1);    // the number of frequencies
        size_t   np;                 // the number of pressure points
        size_t   row, col;           // counters
        MATRIX   abs2 ;              // matrix to store interpolated abs values
       Numeric   w;                  // = -l_step/2
 
  // Resize the transmission array
  trans.newsize(n);

  // Loop the viewing angles and:
  //  1. interpolate the absorption
  //  2. calculate the transmission using the mean absorption between points
  for (size_t i=1; i<=n; i++ ) 
  {
    np = los.p(i).dim();
    if ( np > 0 )
    {
      interp_lin_row( abs2, p_abs, abs, los.p(i) );
      trans(i).newsize( nf, np-1 );
      w  =  -0.5*los.l_step(i);
      for ( row=1; row<=nf; row++ )
      {
        for ( col=1; col<np; col++ )
          trans(i)(row,col) = exp( w * ( abs2(row,col)+abs2(row,col+1)) );
      }
    }
  }    
}



void y_spaceStd(
                    VECTOR&   y_space,
              const VECTOR&   f,
              const int&      nr )
{
  if ( nr == 0 )
  {
    y_space.newsize(f.dim());
    y_space = 0.0;
  }
  else if ( nr == 1 )
    planck( y_space, f, COSMIC_BG_TEMP );
  else if ( nr == 2 )
    planck( y_space, f, SUN_TEMP );
  else
    throw runtime_error("Possible choices for Y_SPACE are 0 - 2.");
}



void y_spacePlanck(
                    VECTOR&   y_space,
              const VECTOR&   f,
              const Numeric&  t )
{
  if ( t > 0 )
    planck( y_space, f, t );
  else
    throw runtime_error("The temperature must be > 0.");
}



void yRte (
                    VECTOR&          y,
              const Los&             los,   
              const VECTOR&          f_abs,
              const VECTOR&          y_space,
              const ARRAYofMATRIX&   source,
              const ARRAYofMATRIX&   trans,
              const VECTOR&          e_ground,
              const Numeric&         t_ground )
{
  // Some variables
  const size_t   n=los.start.dim();  // Number of viewing angles 
  const size_t   nf=f_abs.size();    // Number of frequencies 
        VECTOR   y_tmp;              // Temporary storage for spectra
        size_t   iy0=0;              // Reference index for output vector
        size_t   iy;                 // Frequency index

  // Resize y
  y.newsize(nf*n);
        
  // Set up vector for ground blackbody radiation if any ground intersection
  // Check also if the ground emission vector has the correct length
  VECTOR   y_ground; 
  if ( any(los.ground) )  
  {
    planck( y_ground, f_abs, t_ground );
    if ( e_ground.dim() != nf )
      throw runtime_error("The frequency and ground emission vectors have different lengths.");
  }

  // Loop viewing angles
  for ( size_t i=1; i<=n; i++ )
  {
    // Iteration is done in seperate function    
    rte_iterate( y_tmp, los.start(i), los.stop(i), trans(i), 
                 source(i), y_space, los.ground(i), e_ground, y_ground);

    // Move values to output vector
    for ( iy=1; iy<=nf; iy++ )    
      y(iy0+iy) = y_tmp(iy);
    iy0 += nf;                    // update iy0
  }
}


void yRteNoGround (
                    VECTOR&          y,
              const Los&             los,   
              const VECTOR&          f_abs,
              const VECTOR&          y_space,
              const ARRAYofMATRIX&   source,
              const ARRAYofMATRIX&   trans )
{
  if ( any(los.ground) )  
    throw runtime_error("There is at least one intersection with the ground and this function cannot be used.");

  yRte( y, los, f_abs, y_space, source, trans, 0.0, 0.0 );
}


void yBl (
                    VECTOR&          y,
              const Los&             los,   
              const ARRAYofMATRIX&   trans,
              const VECTOR&          e_ground )
{
  // Some variables
  const size_t   n=los.start.dim();    // Number of viewing angles 
  const size_t   nf=trans(1).dim(1);   // Number of frequencies 
        size_t   iy0=0;                // Reference index for output vector
        size_t   iy;                   // Frequency index
        int      j;                    // LOS index

  // Resize y and set to 1
  y.newsize(nf*n);
  y = 1.0;

  // Check if the ground emission vector has the correct length
  if ( any(los.ground) )  
  {
    if ( e_ground.dim() != nf )
      throw runtime_error("The frequency and ground emission vectors have different lengths.");
  }
        
  // Loop viewing angles
  for ( size_t i=1; i<=n; i++ )
  {
    // Loop steps passed twice
    for ( j=1; j<los.stop(i); j++ )
    {
      for ( iy=1; iy<=nf; iy++ )    
        y(iy0+iy) *= trans(i)(iy,j)*trans(i)(iy,j);
    }

    // Loop remaining steps
    for ( j=los.stop(i); j<los.start(i); j++ )
    {
      for ( iy=1; iy<=nf; iy++ )    
        y(iy0+iy) *= trans(i)(iy,j);
    }

    // Include effect of ground reflection
    if ( los.ground(i) )
    {
      for ( iy=1; iy<=nf; iy++ )    
        y(iy0+iy) *= ( 1.0 - e_ground(iy) );
    }

    // Update iy0
    iy0 += nf;   
  }
}


void yBlNoGround (
                    VECTOR&          y,
              const Los&             los,   
              const ARRAYofMATRIX&   trans )
{
  if ( any(los.ground) )  
    throw runtime_error("There is at least one intersection with the ground and this function cannot be used.");

  yBl( y, los, trans, VECTOR(0) );
}
