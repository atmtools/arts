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


void los_1d_geom (
	          Numeric&   z_tan,
                  Index&     nz,
		  Vector&    z,
		  Vector&    ip_z,
		  Vector&    lat,
		  Vector&    l_step,
 		  Index&     ground,
	    const Numeric&   z_ground,
	    const Vector&    z_abs,
	    const Numeric&   r_geoid,
	    const Numeric&   l_max,
	    const Numeric&   z_plat,
	    const Numeric&   za )
{
  // Asserts
  assert( z_ground >= z_abs[0] );
  assert( l_max > 0 );

  // Give all index return arguments a value.
  z_tan   = 9999e3;
  nz      = 0;
  ground  = 0;  

  // Determine the lowest point of the LOS, z1, and the zenith angle of the LOS
  // at this point, za1. The tangent altitude and ground flag are set on the 
  // same time.
  //
  Numeric z1, za1;
  //
  if ( za < 90 )   // Upward observation (no tangent point)
    {
      z1    = z_plat;
      za1   = za;
    }
  else              // Downward observation (limb sounding)
    {
      z_tan = (r_geoid+z_plat) * sin(DEG2RAD*za) - r_geoid; 
      if ( z_tan >= z_ground )   // No intersection with the ground
	{
	  z1  = z_tan;
	  za1 = 90;
	}
      else                       // Intersection with ground
	{
	  ground = 1;
	  z1     = z_ground;
	  za1    = RAD2DEG * asin( (r_geoid+z_tan) / (r_geoid+z_ground) );
	}
    }

  // Get highest absorption altitude
  const Index     lz    = z_abs.nelem();
  const Numeric   z_max = z_abs[lz-1];

  cout << "z1  = " << z1 << "\n";
  cout << "za1 = " << za1 << "\n";

  // Do anything only if z1 is inside the atmosphere
  if ( z1 < z_max )
    {
      // Determine index of first z_abs above z1, i1.
      Index   i1;
      for( i1=0; z_abs[i1]<=z1; i1++ ) {}
  cout << "i1  = " << i1 << "\n";

      // Create a vector containing z1 and the z_abs levels above z1 (zs),
      // and calculate the length along the LOS from z1 (ls) and the number of
      // LOS steps needed to reach next point of z1 (ns).
      // There can be numerical problems if a,b and c are float.
      const Index   n = lz - i1 + 1;
            Vector  zs(n), ls(n);
      ArrayOfIndex  ns(n-1);
            double  a;
            double  b = (r_geoid+z1) * sin(DEG2RAD*za1);
            double  c = (r_geoid+z1) * cos(DEG2RAD*za1);
            Index   n_sum = 0, i, j;
      //
      zs[0] = z1; 
      ls[0] = 0; 
      //
      for ( i=1; i<n; i++ )
	{
          zs[i]   = z_abs[i1+i-1];
	  a       = r_geoid + zs[i];
	  ls[i]   = sqrt( a*a - b*b ) - c;
          ns[i-1] = Index( ceil( (ls[i]-ls[i-1])/l_max ) );
          n_sum  += ns[i-1];
	}
      n_sum += 1;   // To account for the point at the last z_abs level

      // The length of the z and lat vectors is n_sum
      nz = n_sum;
  cout << "ns = \n" << ns << "\n";
      
      // Create the return vectors, beside lat
      z.resize( nz );
      ip_z.resize( nz );
      l_step.resize( nz-1 );
      //
      Numeric l, l1, l2;
      b = r_geoid + z1;  
      a = b * b;
      n_sum = 0;    // Re-use this variable
      //
      for ( i=0; i<n-1; i++ )
	{
          // Altitude and index position for next value in zs
          // The index position for i=0 is adjusted below.
          z[n_sum]      = zs[i];
          ip_z[n_sum]   = i1 + i - 1;

          // Only one step to next z_abs level
          if ( ns[i] == 1 )
	    l_step[n_sum] = ls[i+1] - ls[i];

          // More than one step to next z_abs level
	  else
	    {
              l = (ls[i+1] - ls[i]) / ns[i];
	      l_step[n_sum] = l;
              for ( j=1; j<ns[i]; j++ )
		{
                  l1   = ls[i] + j*l;
                  l2   = l1 * l1;
		  l_step[n_sum+j] = l;
		  z[n_sum+j]      = sqrt( a + l2 + 2*b*l1*cos(DEG2RAD*za1) ) - 
                                                                       r_geoid;
		  ip_z[n_sum+j]   = ip_z[n_sum] + 
                      (z[n_sum+j]-z_abs[i1+i-1]) / (z_abs[i1+i]-z_abs[i1+i-1]);
		}
	    }

          // Increase n_sum with the points done
          n_sum += ns[i];
	}
      
      // Adjust the index position for the z1 point
      ip_z[0]   = i1 - 1 + (z1-z_abs[i1-1]) / (z_abs[i1]-z_abs[i1-1]);

      // Put in uppermost z_abs level that is not covered above
      z[n_sum]    = z_abs[lz-1];
      ip_z[n_sum] = lz-1;

      // Calculate latitudes corresponding to z
      lat.resize( nz );
    }
}


void test_new_los()
{
  LOS los;
 
  Numeric z_tan, z_ground, l_max, z_plat, za;

  z_ground = 0;
  l_max    = 20e3;
  z_plat   = 20e3;
  za       = 45;
  
  Vector z_abs( 0, 11, 1e3 );

  los_1d_geom ( z_tan, los.np, los.z, los.ip_p, 
                los.lat, los.l_step, los.ground, z_ground, z_abs, EARTH_RADIUS,
 		                                           l_max, z_plat, za );

  cout << "z = \n" << los.z << "\n";
  cout << "ip_p = \n" << los.ip_p << "\n";
  cout << "l_step = \n" << los.l_step << "\n";
  //  cout << "lat = \n" << los.lat << "\n";
  cout << "z_tan = " << z_tan << "\n";;
}
