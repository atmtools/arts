#include <math.h>
#include "arts.h"
#include "matpackI.h"
#include "messages.h"          
#include "auto_wsv.h"          

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric EARTH_RADIUS;

struct LOS {
  Index     dim;
  Index     np;
  Index     i_start;
  Index     i_stop;
  Vector    p;
  Vector    z;
  Vector    ip_p;
  Vector    lat;
  Vector    ip_lat;
  Vector    lon;
  Vector    ip_lon;
  Vector    l_step;
  Index     background;
  Index     ground;
  Index     i_ground;
};


/** Returns the first value of a vector.

    \return       The first value of x.
    \param    x   A vector.

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric first2( ConstVectorView x )
{
  return x[0]; 
}

/** Returns the last value of a vector.

    \return      The last value of x.
    \param   x   A vector.

    \author Patrick Eriksson 
    \date   2000-06-27
*/
Numeric last2( ConstVectorView x )
{
  return x[x.nelem()-1]; 
}


void los_1d_geom (
	          Numeric&    z_tan,
                  Index&      nz,
		  Vector&     z,
		  Vector&     ip_z,
		  Vector&     lat,
		  Vector&     l_step,
 		  Index&      ground,
	    const Numeric&    z_ground,
	    ConstVectorView   z_abs,
	    const Numeric&    r_geoid,
	    const Numeric&    l_max,
	    const Numeric&    z_plat,
	    const Numeric&    za )
{
  // Asserts
  assert( z_ground >= z_abs[0] );
  assert( z_ground < last2(z_abs) );
  assert( z_abs.nelem() > 1 );
  assert( l_max > 0 );
  assert( z_plat >= z_ground );
  assert( za >= 0 );
  assert( za <= 180 );

  // Guess a value for all index return arguments.
  z_tan   = 9999e3;
  nz      = 0;
  ground  = 0;  

  // Get highest absorption altitude and length of z_abs
  const Index     n_zabs = z_abs.nelem();
  const Numeric   z_max  = z_abs[n_zabs-1];

  // Determine the lowest point of the LOS, z1, and the zenith angle of the LOS
  // at this point, za1. The tangent altitude, ground flag and some other 
  // variables are set on the same time.
  // The latitude distance from the sensor and z1 is denoted as lon0.
  // The case with downward observation from inside the atmosphere nees special
  // treatment, handled by the flag do_down.
  //
  Numeric z1, za1, lat0;
  Index   do_down = 0;
  //
  if ( za <= 90 )   // Upward observation (no tangent point)
    {
      z1   = z_plat;
      za1  = za;
      lat0 = 0;
    }
  else              // Downward observation (limb sounding)
    {
      z_tan = (r_geoid+z_plat) * sin(DEG2RAD*za) - r_geoid; 
      if ( z_tan >= z_ground )   // No intersection with the ground
	{
	  z1   = z_tan;
	  za1  = -90;
          lat0 = za - 90.0;
	}
      else                       // Intersection with ground
	{
	  ground = 1;
	  z1     = z_ground;
	  za1    = -RAD2DEG * asin( (r_geoid+z_tan) / (r_geoid+z_ground) );
          lat0   = za - za1 - 180.0;
	}
      if ( z_plat < z_max )
	do_down = 1;
    }

  // The return vectors are set to be empty if z1 >= z_max
  if ( z1 >= z_max )
    {
      z.resize(0);
      ip_z.resize(0);
      l_step.resize(0);
      lat.resize(0);
    }

  // Do anything only if z1 is inside the atmosphere
  else
    {
      // Create vectors for special points.
      // Special points are z1 and z_plat (if do_down).
      // This vector shall start with z1 and end with a dummy value > z_max.
      const Index    n_special = 1 + do_down;  // The dummy value is ignored
            Vector   z_special(n_special+1); 
      //
      z_special[0] = z1;
      if ( do_down )
	z_special[1] = z_plat;
      z_special[do_down+1] = z_max*2;

      // Determine index of first z_abs above z1, i_above.
      Index   i_above;
      for ( i_above=0; z_abs[i_above]<=z1; i_above++ ) {}

      // Create a vector containing z_special and z_abs levels above z1 (zs).
      // The altitudes shall be sorted and there shall be no duplicates.
      // It is assumed that the first shall be taken from z_special.
      Index    n_zs = n_zabs - i_above + n_special; 
      Vector   zs(n_zs);
      Index    i1, i2=1;
      //
      zs[0] = z_special[0];
      n_zs  = 1;            // n_zs counts now the number of values moved to zs
      for ( i1=i_above; i1<n_zabs; i1++  )
	{
	  // Check if values from z_special shall be copied
	  for ( ; z_special[i2] <= z_abs[i1]; i2++ )
	    {
	      if ( zs[n_zs] != z_special[i2] )
		{
		  zs[n_zs] = z_special[i2];
		  n_zs++;
		}
	    }
	  // Copy next z_abs
	  if ( zs[n_zs] != z_abs[i1] )
	    {
	      zs[n_zs] = z_abs[i1];
	      n_zs++;
	    }
	}

      // Calculate the length along the LOS from z1 (ls) and the number of
      // LOS steps needed to reach next altitude in zs (ns).
      // The geometric calculations can go wrong if float is used, and to
      // avoid this, the variables a-e are hard-coded to be double.
      Vector         ls(n_zs);
      ArrayOfIndex   ns(n_zs-1);               
      double         a, b, c, d, e;
      Index          n_sum = 0;       // n_sum is sum(ns)
      //
      // Handle first point seperately
      zs[0] = zs[0]; 
      ls[0] = 0; 
      //
      // Loop zs
      a  = cos( DEG2RAD * za1 );
      b  = sin( DEG2RAD * za1 );
      c  = ( r_geoid + z1 ) * b;
      c *= c;
      d  = ( r_geoid + z1 ) * a;
      for ( i1=1; i1<n_zs; i1++ )
	{
	  e        = r_geoid + zs[i1];
	  ls[i1]   = sqrt( e*e - c ) - d;
          ns[i1-1] = Index( ceil( (ls[i1]-ls[i1-1])/l_max ) );
          n_sum   += ns[i1-1];
	}
      n_sum += 1;   // To account for the point at the last z_abs level

      // The length of the z and lat vectors is n_sum
      nz = n_sum;
      
      // Create the return vectors.
      z.resize(nz);
      ip_z.resize(nz);
      l_step.resize(nz-1);
      lat.resize(nz);
      //
      Numeric   dl, l;              // Different lengts along the LOS
      Numeric   zv;                 // Next altitude
      c = r_geoid + z1;  
      d = c * c;
      n_sum = 0;                 // Re-use this variable
      //
      for ( i1=0; i1<n_zs-1; i1++ )
	{
          dl = (ls[i1+1] - ls[i1]) / ns[i1];
	  for ( i2=0; i2<ns[i1]; i2++ )
	    {
	      l = ls[i1] + i2*dl;
	      if ( i2 == 0 )
		zv = zs[i1];
	      else
		zv = sqrt( d + l*l + 2 * c * l * a ) - r_geoid;
	      z[n_sum+i2]      = zv;
	      l_step[n_sum+i2] = dl;
	      ip_z[n_sum+i2]   = i_above - 1 + ( zv - z_abs[i_above-1]) /
                                           (z_abs[i_above] - z_abs[i_above-1]);
	      lat[n_sum+i2]    = lat0 + RAD2DEG*asin( l*b / (r_geoid+zv) ); 
	    }

          // Increase n_sum with the points done
          n_sum += ns[i1];

          // Check if i_above shall be increased
	  if ( zs[i1] >= z_abs[i_above] )
	    i_above++;
	}
      
      // Put in uppermost z_abs level that is not covered above
      z[n_sum]    = z_abs[n_zabs-1];
      ip_z[n_sum] = n_zabs-1;
      lat[n_sum]  = lat0 + RAD2DEG*asin( ls[n_zs-1]*b / (r_geoid+z[n_sum]) ); 
    }
}



void los_1d (
	          LOS&        los,
                  Numeric&    z_tan,
	    const Numeric&    z_ground,
	    ConstVectorView   z_abs,
	    const Numeric&    r_geoid,
	    const Numeric&    l_max,
	    const Numeric&    z_plat,
	    const Numeric&    za,
            const Index&      refr_on,
            const Index&      blackbody_ground,
	    const Index&      scattering_on )
{
  // Check input 
  // za is inside 0-180

  // LOS variables that always are the same
  los.dim = 1;
  los.ip_lat.resize(0);
  los.lon.resize(0);
  los.ip_lon.resize(0);
  los.i_ground   = 0;

  // Set background to CBGR
  los.background = 0;

  // Do stuff that differ between with and without refraction
  if ( refr_on )
    throw runtime_error(
                  "1D LOS calculations with refraction not yet implemented" );
  else
    los_1d_geom ( z_tan, los.np, los.z, los.ip_p, los.lat, los.l_step, 
               los.ground, z_ground, z_abs, EARTH_RADIUS, l_max, z_plat, za );

  // Set i_start and i_stop assuming blackbody_ground and scattering are 0
  if ( za <= 90 )
    los.i_stop = 0;
  else
    los.i_stop = los.np - 1;
  los.i_start = los.np - 1;

  // Downward observation from inside the atmosphere needs special treatment
  if ( za>90 && z_plat<last2(z_abs) )
      for ( los.i_stop=0; los.z[los.i_stop]!=z_plat; los.i_stop++ ) {}

  // Ignore part of the LOS before ground reflection if blackbody ground
  // The start index is then always 0. Note that only i_start is changed and
  // the vectors are not shorten.
  if ( los.ground && blackbody_ground )
    {
      los.background = 1;
      los.i_start    = 0;
      los.np         = los.i_stop;
    }

  // Without scattering
  if ( scattering_on )
    throw runtime_error(
                   "1D LOS calculations with scattering not yet implemented" );

  // Convert altitudes to pressures
}



void test_new_los()
{
  LOS los;
 
  Numeric z_tan, z_ground, l_max, z_plat, za;

  z_ground = 0;
  l_max    = 100e3;
  z_plat   = 5e3;
  za       = 90;
  
  Vector z_abs( 0, 11, 1e3 );

  los_1d( los, z_tan, z_ground, z_abs, EARTH_RADIUS, l_max, z_plat, za,
                                                                     0, 1, 0 );

  cout << "z = \n" << los.z << "\n";
  //cout << "ip_p = \n" << los.ip_p << "\n";
  cout << "lat = \n" << los.lat << "\n";
  cout << "l_step = \n" << los.l_step << "\n";
  cout << "nz      = " << los.np << "\n";
  cout << "i_start = " << los.i_start << "\n";
  cout << "i_stop  = " << los.i_stop << "\n";
  cout << "bground = " << los.background << "\n";
  cout << "ground  = " << los.ground << "\n";
  cout << "i_ground= " << los.i_ground << "\n";
  cout << "z_tan   = " << z_tan/1e3 << " km\n";
}
