/* Copyright (C) 2000 Patrick Eriksson <patrick@rss.chalmers.se>

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


/////////////////////////////////////////////////////////////////////////////
//
// This file contains functions to determine the line of sight (LOS)
// to calculate variables along the LOS and to solve the radiative transfer
// along the LOS.
// Functions in this file assumes LTE and no scattering
// The LOS is defined by a structure of type Los, defined in los.h.
//
/////////////////////////////////////////////////////////////////////////////


#include "arts.h"
#include "vecmat.h"
#include "messages.h"          
#include "wsv.h"          
#include "math_funcs.h"          
#include "atm_funcs.h"          
#include "los.h"


extern const Numeric EARTH_RADIUS;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric COSMIC_BG_TEMP;
extern const Numeric SUN_TEMP;



//==========================================================================
//=== LOS help functions 
//==========================================================================

// UPWARD_GEOM
// 
// A help function to calculate LOS for upward observations.
// No refraction.
//
// Patrick Eriksson 27.06.99

void upward_geom(
              VECTOR&       z,
              Numeric&      l_step,
        const Numeric&      z_plat,
        const Numeric&      za,
        const Numeric&      z_abs_max )
{
  Numeric  a, b;          // temporary values
  Numeric  llim;          // distance to atmospheric limit
  VECTOR   l;             // distance from sensor

  if ( za > 90 )
    throw logic_error("Upward function used for zenith angle > 90 deg."); 

  a     = EARTH_RADIUS + z_abs_max;
  b     = (EARTH_RADIUS + z_plat)*sin(DEG2RAD*za);
  llim  = sqrt(a*a-b*b) - (EARTH_RADIUS+z_plat)*cos(DEG2RAD*za) ;

  if ( llim < l_step )   // Handle the rare case that llim < l_step
    l_step = llim*0.999; // *0.999 to avoid problem in interpolations

  l = linspace( 0, llim, l_step );
  b = EARTH_RADIUS + z_plat;  
  z = sqrt(b*b+emult(l,l)+(2.0*b*cos(DEG2RAD*za))*l) - EARTH_RADIUS;
}



// SPACE_GEOM
// 
// A function to calculate LOS for observations from space.
// No refraction.
//
// Patrick Eriksson 27.06.99

void space_geom(
               VECTOR&     z,
               Numeric&    l_step,
         const Numeric&    z_tan,
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

    if ( llim < l_step )   // Handle the rare case that llim < l_step
      l_step = llim*0.999; // *0.999 to avoid problem in interpolations


    linspace( l, 0, llim, l_step );
    z = sqrt(b*b+emult(l,l)) - EARTH_RADIUS; 
  }   

  // Intersection with the ground
  else
  {
    // Determine the "zenith angle" at ground level and call upward function 
    Numeric za = RAD2DEG*asin((EARTH_RADIUS+z_tan)/(EARTH_RADIUS+z_ground));
    upward_geom( z, l_step, z_ground, za, z_abs_max );    
  }
}



//==========================================================================
//=== Sub-functions to losBasic 
//==========================================================================

// LOS_NO_REFR_SPACE
//
// Observation from space, limb or nadir. 
// No refraction.
//
// Patrick Eriksson 27.06.99

void los_no_refr_space(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     za,
              const Numeric&    l_step,
              const Numeric&    z_abs_max,
              const Numeric&    z_ground )
{
  Numeric z_tan;           // the tangent altitude
  int     n = za.dim();    // the number of zenith angles

  // Set all step lengths to the user defined value as a first guess.
  // Note that the step length can be changed in SPACE_GEOM.
  los.l_step = l_step;

  // Loop the zenith angles
  for ( int i=1; i<=n; i++ )
  { 
    z_tan = ztan_geom( za(i), z_plat);
    space_geom( los.p(i), los.l_step(i), z_tan, z_abs_max, z_ground );
    los.start(i)  = los.p(i).dim();
    los.stop(i)   = los.start(i);
    if ( z_tan >= z_ground )
      los.ground(i) = 0;          // no ground intersection
    else
      los.ground(i) = 1;          // ground at index 1
  }
}



// LOS_NO_REFR_INSIDE
//
// Observation from point inside the atmosphere, upward or downward.
// No refraction.
//
// Patrick Eriksson 27.06.99

void los_no_refr_inside(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     za,
              const Numeric&    l_step,
              const Numeric&    z_abs_max,
              const Numeric&    z_ground )
{
  Numeric z_tan, a, b, c, l1;  // see below
  int     n = za.dim();        // the number of zenith angles

  // Set all step lengths to the user defined value as a first guess.
  // Note that the step length can be changed in UPWARD_GEOM.
  los.l_step = l_step;

  // Loop the zenith angles
  for ( int i=1; i<=n; i++ )
  { 
    // Upward
    if ( za(i) <= 90 )
    {
      upward_geom( los.p(i), los.l_step(i), z_plat, za(i), z_abs_max );    
      los.start(i)  = los.p(i).dim();   // length of LOS
      los.stop(i)   = 1;                // stop index here always 1
      los.ground(i) = 0;                // no ground intersection
    }

    // Downward
    else
    {
      // Calculate the tangent altitude
      z_tan = ztan_geom(za(i),z_plat);

      // Only through the atmosphere
      if ( z_tan >= z_ground )
      {
	a      = EARTH_RADIUS + z_plat;   // help variable
	b      = EARTH_RADIUS + z_tan;    // help variable
	l1     = sqrt(a*a-b*b);           // distance platform-tangent point
        // Adjust l_step downwards to get an integer number of steps
	los.stop(i)   = (int) ceil( l1 / l_step + 1.0 );  
	los.l_step(i) = l1 / ( (Numeric)los.stop(i) - 1.0 );
	space_geom( los.p(i), los.l_step(i), z_tan, z_abs_max, z_ground );
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
	upward_geom( los.p(i), los.l_step(i), z_ground, RAD2DEG*asin(c/a),
                                                                   z_abs_max);
        los.start(i)  = los.p(i).dim();
        los.ground(i) = 1;                // ground at index 1
      }
    }
  }
}



//==========================================================================
//=== Workspace methods
//==========================================================================

// losCalc
//
// Patrick Eriksson 22.05.00

void losCalc(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     za,
              const Numeric&    l_step,
              const VECTOR&     p_abs,
              const VECTOR&     z_abs,
              const int&        refr,
              const Numeric&    l_step_refr,
              const VECTOR&     refr_index,
              const Numeric&    z_ground )
{     
  size_t n = za.dim();  // number of zenith angles

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

  // Print messages
  if ( refr == 0 )
    out2 << "  Calculating line of sights WITHOUT refraction.\n";
  else if ( refr == 1 )
    out2 << "  Calculating line of sights WITH refraction.\n";
  else
    throw runtime_error("The refraction flag can only be 0 or 1.");
  out3 << "     z_plat: " << z_plat/1e3 << " km\n";
  if ( n == 1 )
    out3 << "         za: " << za(1) << " degs.\n";
  else
  {
    out3 << "     min za: " << min(za) << " degs.\n";
    out3 << "     max za: " << max(za) << " degs.\n";
  }

  // Without refraction
  if ( refr == 0 )
  {
    // The functions below calculates the vertical altitudes along LOS, 
    // stored temporarily in los.p.
    if ( z_plat >= z_abs_max )
      los_no_refr_space( los, z_plat, za, l_step, z_abs_max, z_ground );
    else
      los_no_refr_inside( los, z_plat, za, l_step, z_abs_max, z_ground );

    // Convert altitudes to pressures.
    for ( size_t i=1; i<=n; i++ )
      z2p( los.p(i), z_abs, p_abs, los.p(i) );
  } 

  // With refraction
  else
  {
  }  
}



// losNoRefraction
//
// Patrick Eriksson 07.06.00

void losNoRefraction(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     za,
              const Numeric&    l_step,
              const VECTOR&     p_abs,
              const VECTOR&     z_abs,
              const Numeric&    z_ground )
{
  losCalc( los, z_plat, za, l_step, p_abs, z_abs, 0, 0, VECTOR(0) ,z_ground); 
}



// losUpward
//
// Patrick Eriksson 07.06.00

void losUpward(
                    Los&        los,
              const Numeric&    z_plat,
              const VECTOR&     za,
              const Numeric&    l_step,
              const VECTOR&     p_abs,
              const VECTOR&     z_abs )
{
  if ( max(za) > 90 )
    throw runtime_error("At least one zenith angle > 90 degrees, that is, not upwards. Use losCalc or losNoRefraction.");

  losCalc( los, z_plat, za, l_step, p_abs, z_abs, 0, 0, VECTOR(0) ,z_plat); 
}



// sourceCalc
//
// Patrick Eriksson 07.06.00

void sourceCalc(
                    ARRAYofMATRIX&   source,
              const Los&             los,   
              const VECTOR&          p_abs,
              const VECTOR&          t_abs,
              const VECTOR&          f_mono )
{     
        VECTOR   tlos;                 // temperatures along the LOS
  const size_t   nza=los.start.dim();  // the number of zenith angles  
  const size_t   nf=f_mono.dim();      // the number of frequencies
        size_t   nlos;                 // the number of pressure points
        MATRIX   b;                    // the Planck function for TLOS  
        size_t   iv, ilos;             // frequency and LOS point index

  out2 << "  Calculating the source function for LTE and no scattering.\n";
 
  // Resize the source array
  source.newsize(nza);

  // Loop the zenith angles and:
  //  1. interpolate the temperature
  //  2. calculate the Planck function for the interpolated temperatures
  //  3. take the mean of neighbouring Planck values
  out3 << "    Zenith angle nr:\n      ";
  for (size_t i=1; i<=nza; i++ ) 
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();

    if ( los.p(i).size() > 0 )
    {
      interpp( tlos, p_abs, t_abs, los.p(i) );
      nlos = tlos.dim();
      planck( b, f_mono, tlos );
      source(i).newsize(nf,nlos-1);
      for ( ilos=1; ilos<nlos; ilos++ )
      {
        for ( iv=1; iv<=nf; iv++ )
          source(i)(iv,ilos) = ( b(iv,ilos) + b(iv,ilos+1) ) / 2.0;
      }
    }
  }  
  out3 << "\n";
}



// transCalc
//
// Patrick Eriksson 07.06.00

void transCalc(
                    ARRAYofMATRIX&   trans,
              const Los&             los,   
              const VECTOR&          p_abs,
              const MATRIX&          abs )
{    
  // Some variables
  const size_t   n = los.start.dim();// the number of zenith angles
  const size_t   nf = abs.dim(1);    // the number of frequencies
        size_t   np;                 // the number of pressure points
        size_t   row, col;           // counters
        MATRIX   abs2 ;              // matrix to store interpolated abs values
       Numeric   w;                  // = -l_step/2

  out2 << "  Calculating transmissions WITHOUT scattering.\n";
 
  // Resize the transmission array
  trans.newsize(n);

  // Loop the zenith angles and:
  //  1. interpolate the absorption
  //  2. calculate the transmission using the mean absorption between points
  out3 << "    Zenith angle nr:\n      ";
  for (size_t i=1; i<=n; i++ ) 
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
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
  out3 << "\n";
}



// y_spaceStd
//
// Patrick Eriksson 07.06.00

void y_spaceStd(
                    VECTOR&   y_space,
              const VECTOR&   f,
              const int&      nr )
{
  if ( nr == 0 )
  {
    y_space.newsize(f.dim());
    y_space = 0.0;
    out2 << "  Setting y_space to zero.\n";
  }
  else if ( nr == 1 )
  {
    planck( y_space, f, COSMIC_BG_TEMP );
    out2 << "  Setting y_space to cosmic background radiation.\n";
  }
  else if ( nr == 2 )
  {
    planck( y_space, f, SUN_TEMP );
    out2 << "  Setting y_space to blackbody radiation corresponding to the Sun temperature\n";
  }
  else
    throw runtime_error("Possible choices for Y_SPACE are 0 - 2.");

}



// y_spacePlanck
//
// Patrick Eriksson 07.06.00

void y_spacePlanck(
                    VECTOR&   y_space,
              const VECTOR&   f,
              const Numeric&  t )
{
  if ( t > 0 )
  {
    planck( y_space, f, t );
    out2<<"  Setting y_space to blackbody radiation for "<<t<<" K.\n";
  }
  else
    throw runtime_error("The temperature must be > 0.");
}



// yRte
//
// Patrick Eriksson 07.06.00

void yRte (
                    VECTOR&          y,
              const Los&             los,   
              const VECTOR&          f_mono,
              const VECTOR&          y_space,
              const ARRAYofMATRIX&   source,
              const ARRAYofMATRIX&   trans,
              const VECTOR&          e_ground,
              const Numeric&         t_ground )
{
  // Some variables
  const size_t   n=los.start.dim();  // Number of zenith angles 
  const size_t   nf=f_mono.size();   // Number of frequencies 
        VECTOR   y_tmp;              // Temporary storage for spectra
        size_t   iy0=0;              // Reference index for output vector
        size_t   iy;                 // Frequency index

  out2 << "  Integrating the radiative transfer eq. WITHOUT scattering.\n";

  // Resize y
  y.newsize(nf*n);
        
  // Set up vector for ground blackbody radiation if any ground intersection
  // Check also if the ground emission vector has the correct length
  VECTOR   y_ground; 
  if ( any(los.ground) )  
  {
    out2 << "  There are intersections with the ground.\n";
    planck( y_ground, f_mono, t_ground );
    if ( e_ground.dim() != nf )
      throw runtime_error("The frequency and ground emission vectors have different lengths.");
  }

  // Loop zenith angles
  out3 << "    Zenith angle nr:\n      ";
  for ( size_t i=1; i<=n; i++ )
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    // Iteration is done in seperate function    
    rte( y_tmp, los.start(i), los.stop(i), trans(i), 
                 source(i), y_space, los.ground(i), e_ground, y_ground);

    // Move values to output vector
    for ( iy=1; iy<=nf; iy++ )    
      y(iy0+iy) = y_tmp(iy);
    iy0 += nf;                    // update iy0
  }
  out3 << "\n";
}



// yRteNoGround
//
// Patrick Eriksson 07.06.00

void yRteNoGround (
                    VECTOR&          y,
              const Los&             los,   
              const VECTOR&          f_mono,
              const VECTOR&          y_space,
              const ARRAYofMATRIX&   source,
              const ARRAYofMATRIX&   trans )
{
  if ( any(los.ground) )  
    throw runtime_error("There is at least one intersection with the ground and this function cannot be used.");

  yRte( y, los, f_mono, y_space, source, trans, 0.0, 0.0 );
}



// yBl
//
// Patrick Eriksson 07.06.00

void yBl (
                    VECTOR&          y,
              const Los&             los,   
              const ARRAYofMATRIX&   trans,
              const VECTOR&          e_ground )
{
  // Some variables
  const size_t   n=los.start.dim();    // Number of zenith angles 
  const size_t   nf=trans(1).dim(1);   // Number of frequencies 
        size_t   iy0=0;                // Reference index for output vector
        size_t   iy;                   // Frequency index
        VECTOR   y_tmp;              // Temporary storage for spectra

  out2 << "  Calculating total transmission WITHOUT scattering.\n";

  // Resize y and set to 1
  y.newsize(nf*n);
  y = 1.0;

  // Check if the ground emission vector has the correct length
  if ( any(los.ground) )  
  {
    out2 << "  There are intersections with the ground.\n";
    if ( e_ground.dim() != nf )
      throw runtime_error("The frequency and ground emission vectors have different lengths.");
  }
        
  // Loop zenith angles
  out3 << "    Zenith angle nr:\n      ";
  for ( size_t i=1; i<=n; i++ )
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    // Iteration is done in seperate function    
    bl( y_tmp, los.start(i), los.stop(i), trans(i), los.ground(i), e_ground );

    // Move values to output vector
    for ( iy=1; iy<=nf; iy++ )    
      y(iy0+iy) = y_tmp(iy);
    iy0 += nf;                    // update iy0
  }
  out3 << "\n";
}



// yBlNoGround
//
// Patrick Eriksson 07.06.00

void yBlNoGround (
                    VECTOR&          y,
              const Los&             los,   
              const ARRAYofMATRIX&   trans )
{
  if ( any(los.ground) )  
    throw runtime_error("There is at least one intersection with the ground and this function cannot be used.");

  yBl( y, los, trans, VECTOR(0) );
}


