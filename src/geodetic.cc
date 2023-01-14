/*===========================================================================
  === File description 
  ===========================================================================*/

/**
    @file    geodetic.cc
    @author  Patrick Eriksson <Patrick.Eriksson@chalmers.se>
    @date    2021-07-29 

    @brief   This file contains the declaration of functions associated
             with the reference ellipsoid, conversion between
             coordinate systems and similar stuff.
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts_conversions.h"
#include "debug.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "math_funcs.h"

#include "arts_constexpr_math.h"
using Math::pow2;

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);


/*===========================================================================
  === The functions, in alphabetical order
  ===========================================================================*/

void ecef2geocentric(VectorView pos,
                     ConstVectorView ecef) {
  pos[0] = sqrt(ecef[0]*ecef[0] + ecef[1]*ecef[1] + ecef[2]*ecef[2]);
  pos[1] = RAD2DEG * asin(ecef[2] / pos[0]);
  pos[2] = RAD2DEG * atan2(ecef[1], ecef[0]);
}


void ecef2geocentric_los(VectorView pos,
                         VectorView los,
                         ConstVectorView ecef,
                         ConstVectorView decef) {
  ecef2geocentric(pos, ecef);

  const Numeric latrad = DEG2RAD * pos[1];
  const Numeric lonrad = DEG2RAD * pos[2];
  const Numeric coslat = cos(latrad);
  const Numeric sinlat = sin(latrad);
  const Numeric coslon = cos(lonrad);
  const Numeric sinlon = sin(lonrad);

  const Numeric dr = coslat*coslon*decef[0] +
                     sinlat*decef[2] +
                     coslat*sinlon*decef[1];
  const Numeric dlat = -sinlat*coslon/pos[0]*decef[0] +
                       coslat/pos[0]*decef[2] -
                       sinlat*sinlon/pos[0]*decef[1];
  const Numeric dlon = -sinlon/coslat/pos[0]*decef[0] +
                       coslon/coslat/pos[0]*decef[1];

  los[0] = acos( dr );  // Conersion to deg below

  // Azimuth needs special treatment at poles
  if (fabs(pos[1]) > POLELATZZZ) {
    los[1] = RAD2DEG * atan2(decef[1], decef[0]);
  } else {
    los[1] = RAD2DEG * acos( pos[0] * dlat / sin( los[0] ) );
    // Corrections of azimuth
    if (std::isnan(los[1])) {
      if (dlat >= 0)
        los[1] = 0;
      else
        los[1] = 180;
    } else if (dlon < 0) {
      los[1] = -los[1];
    }
  }
  los[0] *= RAD2DEG;
}


void ecef2geodetic(VectorView pos,
                   ConstVectorView ecef,
                   ConstVectorView refellipsoid) {
  // Use geocentric function if geoid is spherical
  if (is_ellipsoid_spherical(refellipsoid)) {
      ecef2geocentric(pos, ecef);
      pos[0] -= refellipsoid[0];
    
  // The general algorithm not stable for lat=+-90. Catch these cases
  } else if (ecef[0]==0 && ecef[1]==0) {
    pos[0] = fabs(ecef[2]) - refellipsoid[1];
    pos[1] = ecef[2] >= 0 ? 90 : -90;
    pos[2] = 0;

  // General algorithm
  } else {
    pos[2] = RAD2DEG * atan2(ecef[1],ecef[0]);

    const Numeric sq = sqrt(ecef[0]*ecef[0] + ecef[1]*ecef[1]);
    Numeric B0 = atan2(ecef[2], sq);
    Numeric B = B0-1, N;
    const Numeric e2 = 1 -
      (refellipsoid[1]*refellipsoid[1])/(refellipsoid[0]*refellipsoid[0]);
    // 1e-15 seems to give a accuracy of better than 2 cm
    while (fabs(B-B0)>1e-15) {
      N = refellipsoid[0] / sqrt(1-e2*sin(B0)*sin(B0));
      pos[0] = sq / cos(B0) - N;
      B = B0;
      B0 = atan((ecef[2]/sq) * 1/(1-e2*N/(N+pos[0])));
    }
    pos[1] = RAD2DEG * B;
  }
}


void ecef2geodetic_los(VectorView pos,
                       VectorView los,
                       ConstVectorView ecef,
                       ConstVectorView decef,
                       ConstVectorView refellipsoid) {
  ecef2geodetic(pos, ecef, refellipsoid);
  
  const Numeric latrad = DEG2RAD * pos[1];
  const Numeric lonrad = DEG2RAD * pos[2];
  const Numeric coslat = cos(latrad);
  const Numeric sinlat = sin(latrad);
  const Numeric coslon = cos(lonrad);
  const Numeric sinlon = sin(lonrad);

  // See
  // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
  Vector enu(3);
  enu[0] =        -sinlon*decef[0] +        coslon*decef[1];  
  enu[1] = -sinlat*coslon*decef[0] - sinlat*sinlon*decef[1] + coslat*decef[2];  
  enu[2] =  coslat*coslon*decef[0] + coslat*sinlon*decef[1] + sinlat*decef[2];  

  enu2los(los, enu);

  // Azimuth at poles needs special treatment
  if (fabs(pos[1]) > POLELATZZZ) 
    los[1] = RAD2DEG * atan2(decef[1], decef[0]);
}


void ecef_at_distance(VectorView ecef,
                      ConstVectorView ecef0,
                      ConstVectorView decef,
                      const Numeric l) {
  ecef[0] = ecef0[0] + l * decef[0];
  ecef[1] = ecef0[1] + l * decef[1];
  ecef[2] = ecef0[2] + l * decef[2];
}


Numeric ecef_distance(ConstVectorView ecef1,
                      ConstVectorView ecef2) {
  const Numeric dx = ecef2[0] - ecef1[0]; 
  const Numeric dy = ecef2[1] - ecef1[1]; 
  const Numeric dz = ecef2[2] - ecef1[2];
  return sqrt(dx*dx + dy*dy + dz*dz);
}


void ecef_vector_distance(VectorView ecef,
                          ConstVectorView ecef0,
                          ConstVectorView ecef1) {
  ecef[0] = ecef1[0] - ecef0[0]; 
  ecef[1] = ecef1[1] - ecef0[1]; 
  ecef[2] = ecef1[2] - ecef0[2]; 
}


void enu2los(VectorView los,
             ConstVectorView enu) {
  // los[0] came out as Nan for a case as enu[2] was just below -1
  // So let's be safe and normalise enu[2], and get a cheap assert for free
  const Numeric twonorm = norm2(enu);
  ARTS_ASSERT(fabs(twonorm - 1.0) < 1e-6);
  los[0] = RAD2DEG * acos( enu[2]/twonorm );
  los[1] = RAD2DEG * atan2( enu[0], enu[1] );
}


void geocentric2ecef(VectorView ecef,
                     ConstVectorView pos) {
  const Numeric latrad = DEG2RAD * pos[1];
  const Numeric lonrad = DEG2RAD * pos[2];
  ecef[0] = pos[0] * cos(latrad);  // Common term for x and z
  ecef[1] = ecef[0] * sin(lonrad);
  ecef[0] = ecef[0] * cos(lonrad);
  ecef[2] = pos[0] * sin(latrad);
}


void geocentric_los2ecef(VectorView ecef,
                         VectorView decef,
                         ConstVectorView pos,
                         ConstVectorView los) {
  // lat = +-90
  // For lat = +- 90 the azimuth angle gives the longitude along which the
  // LOS goes
  if (fabs(pos[1]) > POLELATZZZ) {
    const Numeric s = sign(pos[1]);
    const Numeric zarad = DEG2RAD * los[0];
    const Numeric aarad = DEG2RAD * los[1];
    ecef[1] = 0;
    ecef[1] = 0;
    ecef[2] = s * pos[0];
    decef[2] = s * cos(zarad);
    decef[0] = sin(zarad);
    decef[1] = decef[0] * sin(aarad);
    decef[0] = decef[0] * cos(aarad);
  }

  else {
    const Numeric latrad = DEG2RAD * pos[1];
    const Numeric lonrad = DEG2RAD * pos[2];
    const Numeric zarad = DEG2RAD * los[0];
    const Numeric aarad = DEG2RAD * los[1];

    const Numeric coslat = cos(latrad);
    const Numeric sinlat = sin(latrad);
    const Numeric coslon = cos(lonrad);
    const Numeric sinlon = sin(lonrad);
    const Numeric cosza = cos(zarad);
    const Numeric sinza = sin(zarad);
    const Numeric cosaa = cos(aarad);
    const Numeric sinaa = sin(aarad);

    // This part as sph2cart but uses local variables
    ecef[0] = pos[0] * coslat;  // Common term for x and y
    ecef[1] = ecef[0] * sinlon;
    ecef[0] = ecef[0] * coslon;
    ecef[2] = pos[0] * sinlat;

    const Numeric dr = cosza;
    const Numeric dlat = sinza * cosaa;  // r-term cancel out below
    const Numeric dlon = sinza * sinaa / coslat;

    decef[0] = coslat * coslon * dr - sinlat * coslon * dlat - coslat * sinlon * dlon;
    decef[2] = sinlat * dr + coslat * dlat;
    decef[1] = coslat * sinlon * dr - sinlat * sinlon * dlat + coslat * coslon * dlon;
  }
}


void approx_geometrical_tangent_point(VectorView ecef_tan,
                                      ConstVectorView ecef,
                                      ConstVectorView decef,
                                      ConstVectorView refellipsoid) {
  // Spherical case (length simply obtained by dot product)
  if (is_ellipsoid_spherical(refellipsoid)) {
    ecef_at_distance(ecef_tan, ecef, decef, -(decef * ecef));
    
  // General case
  } else {
    // The algorithm used for non-spherical cases is derived by Nick Lloyd at
    // University of Saskatchewan, Canada (nick.lloyd@usask.ca), and is part of
    // the operational code for both OSIRIS and SMR on-board- the Odin
    // satellite.

    // It seems that there is some numerical inaccuracy if the observation is
    // done from above one of the poles (lat = +-90deg)
    
    const Numeric a2 = refellipsoid[0] * refellipsoid[0];
    const Numeric b2 = refellipsoid[1] * refellipsoid[1];
    Vector yunit(3), zunit(3);

    cross3(zunit, decef, ecef);
    zunit /= norm2(zunit);
    cross3(yunit, zunit, decef);
    yunit /= norm2(yunit);
  
    const Numeric yr = ecef * yunit;
    const Numeric xr = ecef * decef;
    const Numeric B = 2.0*((decef[0]*yunit[0] + decef[1]*yunit[1])/a2 +
                           (decef[2]*yunit[2]) / b2);
    Numeric xx;
    if (B == 0.0) {
      xx = 0.0;
    } else {
      const Numeric A = (decef[0]*decef[0] + decef[1]*decef[1])/a2 + decef[2]*decef[2]/b2;
      const Numeric C = (yunit[0]*yunit[0] + yunit[1]*yunit[1])/a2 + yunit[2]*yunit[2]/b2;
      const Numeric K = -2.0*A/B;
      const Numeric factor = 1.0/(A+(B+C*K)*K);
      xx = sqrt(factor);
      const Numeric yy = K*ecef[0];
      const Numeric dist1 = (xr-xx)*(xr-xx) + (yr-yy)*(yr-yy);
      const Numeric dist2 = (xr+xx)*(xr+xx) + (yr+yy)*(yr+yy);
      if (dist1 > dist2)
        xx = -xx;
    }

    ecef_tan[0] = decef[0]*xx + yunit[0]*yr;
    ecef_tan[1] = decef[1]*xx + yunit[1]*yr;
    ecef_tan[2] = decef[2]*xx + yunit[2]*yr;
  }
}


void geodetic2ecef(VectorView ecef,
                   ConstVectorView pos,
                   ConstVectorView refellipsoid ) {
  // Use geocentric function if geoid is spherical
  if (is_ellipsoid_spherical(refellipsoid))
    {
      Vector posc = pos;
      posc[0] += refellipsoid[0];
      geocentric2ecef(ecef, posc);
    }
  else
    {
      // See https://en.wikipedia.org/wiki/
      // Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates
      const Numeric latrad = DEG2RAD * pos[1];
      const Numeric lonrad = DEG2RAD * pos[2];
      const Numeric sinlat = sin(latrad);
      const Numeric coslat = cos(latrad);
      const Numeric a2 = refellipsoid[0]*refellipsoid[0];
      const Numeric b2 = refellipsoid[1]*refellipsoid[1];
      const Numeric N = a2 / sqrt(a2*coslat*coslat+b2*sinlat*sinlat);
      const Numeric nhcos = (N + pos[0]) * coslat;
      ecef[0] = nhcos * cos(lonrad);
      ecef[1] = nhcos * sin(lonrad);
      ecef[2] = ((b2/a2)*N + pos[0]) * sinlat;
    }
}


void geodetic_los2ecef(VectorView ecef,
                       VectorView decef,
                       ConstVectorView pos,
                       ConstVectorView los,
                       ConstVectorView refellipsoid) {
  // lat = +-90
  // For lat = +- 90 the azimuth angle gives the longitude along which the
  // LOS goes
  // At the poles, no difference between geocentric and geodetic zenith
  if (fabs(pos[1]) > POLELATZZZ) {
    const Numeric s = sign(pos[1]);
    const Numeric zarad = DEG2RAD * los[0];
    const Numeric aarad = DEG2RAD * los[1];
    ecef[0] = 0;
    ecef[1] = 0;
    ecef[2] = s * (pos[0] + refellipsoid[1]);
    decef[2] = s * cos(zarad);
    decef[0] = sin(zarad);
    decef[1] = decef[0] * sin(aarad);
    decef[0] = decef[0] * cos(aarad);
  }

  else {      
    const Numeric latrad = DEG2RAD * pos[1];
    const Numeric lonrad = DEG2RAD * pos[2];
    const Numeric coslat = cos(latrad);
    const Numeric sinlat = sin(latrad);
    const Numeric coslon = cos(lonrad);
    const Numeric sinlon = sin(lonrad);

    geodetic2ecef(ecef, pos, refellipsoid);

    Vector enu(3);
    los2enu(enu, los);  

    // See
    // https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
    decef[0] = -sinlon*enu[0] - sinlat*coslon*enu[1] + coslat*coslon*enu[2];  
    decef[1] =  coslon*enu[0] - sinlat*sinlon*enu[1] + coslat*sinlon*enu[2];  
    decef[2] =                  coslat*       enu[1] + sinlat*       enu[2];  
  }
}


Numeric intersection_altitude(ConstVectorView ecef,
                              ConstVectorView decef,
                              ConstVectorView refellipsoid,
                              const Numeric& altitude,
                              const Numeric& l_min) {  
  Numeric l;
  Vector ellipsoid = refellipsoid;
  ellipsoid += altitude;

  // Code taken from Atmlab's ellipsoid_intersection
  
  // Spherical case
  if (is_ellipsoid_spherical(ellipsoid)) {
    const Numeric p  = ecef[0]*decef[0] + ecef[1]*decef[1] + ecef[2]*decef[2];
    const Numeric pp = p*p;
    const Numeric q = ecef[0]*ecef[0] + ecef[1]*ecef[1] + ecef[2]*ecef[2] -
      ellipsoid[0]*ellipsoid[0];
    if (q>pp)
      l = l_min-1.0;
    else {
      const Numeric sq = sqrt(pp - q);
      l = min_geq(-p-sq, -p+sq, l_min);
    }
  }

  // Ellipsoid case
  else {
    // Based on https://medium.com/@stephenhartzell/
    // satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6
    const Numeric a = ellipsoid[0];
    const Numeric b = ellipsoid[0];
    const Numeric c = ellipsoid[1];
    const Numeric a2 = a*a;
    const Numeric b2 = b*b;
    const Numeric c2 = c*c;
    const Numeric x2 = ecef[0]*ecef[0];
    const Numeric y2 = ecef[1]*ecef[1];
    const Numeric z2 = ecef[2]*ecef[2];
    const Numeric dx2 = decef[0]*decef[0];
    const Numeric dy2 = decef[1]*decef[1];
    const Numeric dz2 = decef[2]*decef[2];
    const Numeric rad =
      a2*b2*dz2 + a2*c2*dy2 - a2*dy2*z2 + 2*a2*decef[1]*decef[2]*ecef[1]*ecef[2] -
      a2*dz2*y2 + b2*c2*dx2 - b2*dx2*z2 + 2*b2*decef[0]*decef[2]*ecef[0]*ecef[2] -
      b2*dz2*x2 - c2*dx2*y2 + 2*c2*decef[0]*decef[1]*ecef[0]*ecef[1] - c2*dy2*x2;
    if (rad<0)
      l = -1.0;
    else {
      const Numeric val =
        -a2*b2*decef[2]*ecef[2] - a2*c2*decef[1]*ecef[1] - b2*c2*decef[0]*ecef[0];
      const Numeric mag = a2*b2*dz2 + a2*c2*dy2 + b2*c2*dx2;
      const Numeric abc = a*b*c*sqrt(rad);
      l = min_geq((val-abc)/mag, (val+abc)/mag, l_min);
    }
  }
  return l;
}


Numeric intersection_latitude(ConstVectorView ecef,
                              ConstVectorView decef,
                              ConstVectorView pos,
                              ConstVectorView los,
                              ConstVectorView refellipsoid,
                              const Numeric& lat) {
  // If already at lat, l=0
  if (pos[1] == lat) {
    return 0;
  // No solution if outside of latitude cone and looking away or
  // looking at zenith or nadir
  } else if ((lat<=0 && pos[1]>lat && fabs(los[1])<=90) ||  // lat on SH
             (lat>=0 && pos[1]<lat && fabs(los[1])>=90) ||  // lat on NH
             (los[0]==0) || (los[0]==180) ) {              // zenith/nadir
    return -1;
  // Solution simple if lat = 0 (-z/dz)
  } else if (lat==0) {
    if (decef[2]==0)
      return -1;
    else
      return -ecef[2]/decef[2];
  // Special treatment of lat=90
  } else if (lat>POLELATZZZ) {
    if (los[1] != 0)
      return -1;
    else {
      if (fabs(decef[0]) > fabs(decef[1]))
        return -ecef[0]/decef[0];
      else
        return -ecef[1]/decef[1];
    }
  // Special treatment of lat=-90
  } else if (-lat>POLELATZZZ) {
    if (fabs(los[1]) != 180)
      return -1;
    else {
      if (fabs(decef[0]) > fabs(decef[1]))
        return -ecef[0]/decef[0];
      else
        return -ecef[1]/decef[1];
    }
  } else {
    // Algorithm based on:
    // lousodrome.net/blog/light/2017/01/03/intersection-of-a-ray-and-a-cone/
    // C: Position of tip of cone
    Vector C(3,0.0);
    if (is_ellipsoid_spherical(refellipsoid)) {
      C[2] = 0;  // At (0,0,0) for spherical planet
    } else {     // At (0,0,z) for ellipsoidal planet
      // Calculate normal to ellipsoid at (0,lat,0) and use it to determine z
      // of cone tip. The distance from (0,lat,0) to z-axis is l=x/dx, so
      // z of cone tip is z-l*dz.
      Vector pos2(3,0.0), los2(2,0.0), ecefn(3), n(3);
      pos2[1] = lat;
      geodetic_los2ecef(ecefn, n, pos2, los2, refellipsoid);
      const Numeric l2axis = ecefn[0]/n[0];
      C[2] = ecefn[2] - l2axis*n[2];
    }
    // V: Vector describing centre of cone
    Vector V(3);
    V[0] = 0; V[1] = 0; V[2] = lat>0 ? 1 : -1;
    // Angle term (cos(lat)^2)
    Numeric costerm = cos(DEG2RAD*(90-fabs(lat)));
    costerm *= costerm;
    // Rename to follow nomenclature on web page
    ConstVectorView D = decef;
    // Vector from C to O
    Vector CO(3);
    ecef_vector_distance(CO, C, ecef);
    // Dot products repeated
    const Numeric DVdot = D * V;
    const Numeric COVdot = CO * V;
    // The a, b, c and delta terms
    const Numeric a = DVdot*DVdot - costerm;
    const Numeric b = 2*(DVdot*COVdot-(D*CO)*costerm);
    const Numeric c = COVdot*COVdot - (CO*CO)*costerm;
    const Numeric d = b*b - 4*a*c;
    //
    if (d<0) {
      return -1;
    } else if (d==0) {
      return -b/(2*a);
    } else {
      const Numeric sqrtd = sqrt(d);
      const Numeric aa = 2*a;
      Numeric l1 = (-b-sqrtd)/aa;
      Numeric l2 = (-b+sqrtd)/aa;
      // Check that crossing is not with -lat
      if (l1>0) {
        Vector P(3);
        ecef_at_distance(P, ecef, decef, l1);
        Vector PC(3);
        ecef_vector_distance(PC, C, P);
        if (PC*V <=0)
          l1 = -1;
      }
      if (l2>0) {
        Vector P(3);
        ecef_at_distance(P, ecef, decef, l2);
        Vector PC(3);
        ecef_vector_distance(PC, C, P);
        if (PC*V <=0)
          l2 = -1;
      }      
      if (l1>0 && l2>0)
        return l1<l2 ? l1 : l2;   // min of l1 and l2 is returned
      else
        return l1<l2 ? l2 : l1;   // max of l1 and l2 is returned
    }
  }
}


Numeric intersection_longitude(ConstVectorView ecef,
                               ConstVectorView decef,
                               ConstVectorView pos,
                               ConstVectorView los,
                               const Numeric& lon) {
  // If already at lon or on a pole, l=0
  if ((pos[2]== lon) || fabs(pos[1])>POLELATZZZ) {
    return 0;
  // No solution if looking straight N or S, or azimuth in wrong direction 
  } else if ((los[1]==0) || (fabs(los[1])==180) ||
             (pos[2]>lon && los[1]>0) || (pos[2]<lon && los[1]<0)) {
    return -1;
  } else {
    const Numeric tanlon = tan(DEG2RAD*lon);
    return (ecef[1]-ecef[0]*tanlon) / (decef[0]*tanlon-decef[1]);
  }
}


bool is_ellipsoid_spherical(ConstVectorView ellipsoid) {
  if (fabs(ellipsoid[0]-ellipsoid[1]) < ellipsoid_radii_threshold)
    return true;
  else
    return false;
}


bool is_lon_in_range(const Numeric& lon,
                     const Numeric& lon_min,
                     const Numeric& lon_max) {
  // All lons are asserted in *move_lon_to_range*
  const Numeric lon_adjusted = move_lon_to_range(lon, lon_min, lon_max);
  if (lon_adjusted < lon_min || lon_adjusted > lon_max)
    return false;
  else
    return true;
}


void los2enu(VectorView enu,
             ConstVectorView los) {
  const Numeric zarad = DEG2RAD * los[0];
  const Numeric aarad = DEG2RAD * los[1];
  const Numeric st = sin(zarad);
  enu[0] = st * sin(aarad);
  enu[1] = st * cos(aarad);
  enu[2] = cos(zarad);
}


Numeric move_lon_to_range(const Numeric& lon,
                          const Numeric& lon_min,
                          const Numeric& lon_max) {
  ARTS_ASSERT(lon_max >= lon_min);
  ARTS_ASSERT(lon_max-lon_min <= 360);
  return lon_max <= 180 ? shift_lon_to_pm180(lon) : shift_lon_to_0to360(lon);
}


void poslos_at_distance(VectorView pos,
                        VectorView los,
                        ConstVectorView ecef,
                        ConstVectorView decef,
                        ConstVectorView refellipsoid,
                        const Numeric l) {
  Vector ecef_new(3);
  ecef_at_distance(ecef_new, ecef, decef, l);
  ecef2geodetic_los(pos, los, ecef_new, decef, refellipsoid);
}


void pos_at_distance(VectorView pos,
                     ConstVectorView ecef,
                     ConstVectorView decef,
                     ConstVectorView refellipsoid,
                     const Numeric l) {
  Vector ecef_new(3);
  ecef_at_distance(ecef_new, ecef, decef, l);
  ecef2geodetic(pos, ecef_new, refellipsoid);
}


Numeric prime_vertical_radius(ConstVectorView refellipsoid,
                              const Numeric& lat) {
  return refellipsoid[0] * (refellipsoid[0] /
    sqrt(pow2(refellipsoid[0] * cos(DEG2RAD * lat)) +
         pow2(refellipsoid[1] * sin(DEG2RAD * lat))));
}



Numeric shift_lon_to_pm180(const Numeric& lon) {
  ARTS_ASSERT (lon >= -180 && lon <= 360);
  return lon <= 180 ? lon : lon-360;
}


Numeric shift_lon_to_0to360(const Numeric& lon) {
  ARTS_ASSERT (lon >= -180 && lon <= 360);
  return lon >= 0 ? lon : lon+360;
}
