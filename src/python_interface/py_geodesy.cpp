#include <geodetic.h>
#include <nanobind/stl/pair.h>

#include "python_interface.h"

namespace Python {
void py_geodesy(py::module_& main_module) {
  auto geo  = main_module.def_submodule("geodetic");
  geo.doc() = "Contains various helper functions to deal with geodesy";

  geo.def("ecef2geocentric",
          &ecef2geocentric,
          "ecef"_a,
          R"(Conversion from ECEF to geocentric coordinates

Parameters
----------
ecef : ~pyarts3.arts.Vector3
  Position (x,y,z)

Returns
-------
pos : ~pyarts3.arts.Vector3
  Geocentric position (r,lat,lon)
)");

  geo.def("ecef2geocentric_los",
          &ecef2geocentric_los,
          "ecef"_a,
          "decef"_a,
          R"(Conversion from ECEF to geocentric coordinates including a LOS

Parameters
----------
ecef : ~pyarts3.arts.Vector3
  Position (x,y,z)
decef : ~pyarts3.arts.Vector3
  Normalised direction vector (dx,dy,dz)

Returns
-------
pos : ~pyarts3.arts.Vector3
  Geocentric position (r,lat,lon)
los : ~pyarts3.arts.Vector2
  Geocentric line-of-sight (za,aa)
)");

  geo.def("ecef2geodetic",
          &ecef2geodetic,
          "ecef"_a,
          "ell"_a,
          R"(Conversion from ECEF to geodetic coordinates

Parameters
----------
ecef : ~pyarts3.arts.Vector3
  Position (x,y,z)
ell : ~pyarts3.arts.Vector2
  Ellipsoid (a,b)

Returns
-------
pos : ~pyarts3.arts.Vector3
  Geodetic position (h,lat,lon)
)");

  geo.def("ecef2geodetic_los",
          &ecef2geodetic_los,
          "ecef"_a,
          "decef"_a,
          "ell"_a,
          R"(Conversion from ECEF to geodetic coordinates

Parameters
----------
ecef : ~pyarts3.arts.Vector3
  Position (x,y,z)
decef : ~pyarts3.arts.Vector3
  Normalised direction vector (dx,dy,dz)
ell : ~pyarts3.arts.Vector2
  Ellipsoid (a,b)

Returns
-------
pos : ~pyarts3.arts.Vector3
  Geodetic position (h,lat,lon)
los : ~pyarts3.arts.Vector2
  Geodetic line-of-sight (za,aa)
)");

  geo.def("enu2los",
          &enu2los,
          "enu"_a,
          R"(Converts ENU unit vector vector to local zenith and azimuth

This function and the sister function los2enu handles transformation of
line-of-sights, from and to ENU (east-north-up). The ENU vector is
normalised to have length 1.

Parameters
----------
enu : ~pyarts3.arts.Vector3
  Pointing (de,dn,du)

Returns
-------
los : ~pyarts3.arts.Vector2
  Geodetic line-of-sight (za,aa)
)");

  geo.def("geocentric2ecef",
          &geocentric2ecef,
          "pos"_a,
          R"(Conversion from geocentric to ECEF coordinates

Parameters
----------
pos : ~pyarts3.arts.Vector3
  Geodetic position (h,lat,lon)

Returns
-------
ecef : ~pyarts3.arts.Vector3
  Position (x,y,z)
)");

  geo.def("geocentric_los2ecef",
          &geocentric_los2ecef,
          "pos"_a,
          "los"_a,
          R"(Conversion from geocentric to ECEF coordinates including a LOS

Parameters
----------
pos : ~pyarts3.arts.Vector3
  Geodetic position (h,lat,lon)
los : ~pyarts3.arts.Vector2
  Geodetic line-of-sight (za,aa)

Returns
-------
ecef : ~pyarts3.arts.Vector3
  Position (x,y,z)
decef : ~pyarts3.arts.Vector3
  Normalised direction vector (dx,dy,dz)
)");

  geo.def("geodetic2ecef",
          &geodetic2ecef,
          "pos"_a,
          "ell"_a,
          R"(Conversion from geodetic to ECEF coordinates.

Parameters
----------
pos : ~pyarts3.arts.Vector3
  Geodetic position (h,lat,lon)
ell : ~pyarts3.arts.Vector2
  Ellipsoid (a,b)

Returns
-------
ecef : ~pyarts3.arts.Vector3
  Position (x,y,z)
)");

  geo.def("geodetic_los2ecef",
          &geodetic_los2ecef,
          "pos"_a,
          "los"_a,
          "ell"_a,
          R"(Conversion from geodetic to ECEF coordinates including a LOS.

Parameters
----------
pos : ~pyarts3.arts.Vector3
  Geodetic position (h,lat,lon)
los : ~pyarts3.arts.Vector2
  Geodetic line-of-sight (za,aa)
ell : ~pyarts3.arts.Vector2
  Ellipsoid (a,b)

Returns
-------
ecef : ~pyarts3.arts.Vector3
  Position (x,y,z)
decef : ~pyarts3.arts.Vector3
  Normalised direction vector (dx,dy,dz)
)");

  geo.def("approx_geometrical_tangent_point",
          &approx_geometrical_tangent_point,
          "ecef"_a,
          "decef"_a,
          "ell"_a,
          R"(Calculates the geometrical tangent point, approximately

The tangent point is defined as the lowest altitude above the reference
ellipsiod (not the lowest radius) of the propagation path.

The tangent point can be below the surface. If the zenith angle is abobe 90
deg, the tangent point is behind the sensor.

After testing it was found that the algorithm is not fully numerically
stable and the method is "flagged" as approximative. It was found that the
calculated tangent point could match zenith angles of 90.3 deg (it should be
90 deg exactly).

Parameters
----------
ecef : ~pyarts3.arts.Vector3
  Start position (x,y,z)
decef : ~pyarts3.arts.Vector3
  Normalised direction vector (dx,dy,dz)
ell : ~pyarts3.arts.Vector2
  Ellipsoid (a,b)

Returns
-------
ecef : ~pyarts3.arts.Vector3
  Tangent position (x,y,z)
)");

  geo.def(
      "intersection_altitude",
      &intersection_altitude,
      "ecef"_a,
      "decef"_a,
      "ell"_a,
      "alt"_a,
      "min_l"_a = 0.0,
      R"(Finds the distance to the intersection between an ECEF line and an ellipsoid

A negative distance is returned if there is no intersection. 

If multiple postive solutions, the smallest distance >= l_min is returned. 

Note that there is no check if *pos* is below *altitude*, and the solution 
can match a position on the other side of the planet.

Parameters
----------
ecef : ~pyarts3.arts.Vector3
  Start position (x,y,z)
decef : ~pyarts3.arts.Vector3
  Normalised direction vector (dx,dy,dz)
ell : ~pyarts3.arts.Vector2
  Ellipsoid (a,b)
alt : ~pyarts3.arts.Numeric
  Altitude of the interserciont
min_l : ~pyarts3.arts.Numeric
  Minimum distance (in case you are in an intersection point, this allows you to move)

Returns
-------
l : ~pyarts3.arts.Numeric
  The intersection distance
)");

  geo.def(
      "intersection_latitude",
      &intersection_latitude,
      "ecef"_a,
      "decef"_a,
      "pos"_a,
      "los"_a,
      "ell"_a,
      "lat"_a,
      R"(Finds the distance to the intersection between an ECEF line and a latitude

A negative distance is returned if there is no intersection.

Parameters
----------
ecef : ~pyarts3.arts.Vector3
  Start position (x,y,z)
decef : ~pyarts3.arts.Vector3
  Normalised direction vector (dx,dy,dz)
pos : ~pyarts3.arts.Vector3
  Geodetic position (h,lat,lon)
los : ~pyarts3.arts.Vector2
  Geodetic line-of-sight (za,aa)
ell : ~pyarts3.arts.Vector2
  Ellipsoid (a,b)
lat : ~pyarts3.arts.Numeric
  Latitude of inteserction

Returns
-------
l : ~pyarts3.arts.Numeric
  The intersection distance
)");

  geo.def(
      "intersection_longitude",
      &intersection_longitude,
      "ecef"_a,
      "decef"_a,
      "pos"_a,
      "los"_a,
      "lon"_a,
      R"(Finds the distance to the intersection between an ECEF line and a longitude

Looks only for local solutions. The azimuth angle in *los* must be in the
direction of the target longitude. Otherwise it is considered that there is
no intersection.

A negative distance is returned if there is no intersection.

Parameters
----------
ecef : ~pyarts3.arts.Vector3
  Start position (x,y,z)
decef : ~pyarts3.arts.Vector3
  Normalised direction vector (dx,dy,dz)
pos : ~pyarts3.arts.Vector3
  Geodetic position (h,lat,lon)
los : ~pyarts3.arts.Vector2
  Geodetic line-of-sight (za,aa)
lon : ~pyarts3.arts.Numeric
  Longitude of inteserction

Returns
-------
l : ~pyarts3.arts.Numeric
  The intersection distance
)");

  geo.def("los2enu",
          &los2enu,
          "los"_a,
          R"(Converts ENU unit vector vector to local zenith and azimuth

This function and the sister function los2enu handles transformation of
line-of-sights, from and to ENU (east-north-up). The ENU vector is
normalised to have length 1.

Parameters
----------
los : ~pyarts3.arts.Vector2
  Geodetic line-of-sight (za,aa)

Returns
-------
ebu : ~pyarts3.arts.Vector3
  Pointing (de,dn,du)
)");

  geo.def("poslos_at_distance",
          &poslos_at_distance,
          "ecef"_a,
          "decef"_a,
          "lon"_a,
          "l"_a,
          R"(Geodetic position and line-of-sight at a given distance

Parameters
----------
ecef : ~pyarts3.arts.Vector3
  Start position (x,y,z)
decef : ~pyarts3.arts.Vector3
  Normalised direction vector (dx,dy,dz)
ell : ~pyarts3.arts.Vector2
  Ellipsoid (a,b)
l : ~pyarts3.arts.Numeric
  The distance

Returns
-------
pos : ~pyarts3.arts.Vector3
  Geodetic position (h,lat,lon)
los : ~pyarts3.arts.Vector2
  Geodetic line-of-sight (za,aa)
)");

  geo.def("pos_at_distance",
          &pos_at_distance,
          "ecef"_a,
          "decef"_a,
          "lon"_a,
          "l"_a,
          R"(Geodetic position a given distance

Parameters
----------
ecef : ~pyarts3.arts.Vector3
  Start position (x,y,z)
decef : ~pyarts3.arts.Vector3
  Normalised direction vector (dx,dy,dz)
ell : ~pyarts3.arts.Vector2
  Ellipsoid (a,b)
l : ~pyarts3.arts.Numeric
  The distance

Returns
-------
pos : ~pyarts3.arts.Vector3
  Geodetic position (h,lat,lon)
)");

  geo.def("prime_vertical_radius",
          &prime_vertical_radius,
          "ell"_a,
          "lay"_a,
          R"(The prime vertical radius

This radius is the distance to the polar axis from the reference
ellipsoid surface in the local nadir direction.

See further https://en.wikipedia.org/wiki/Earth_radius#Prime_vertical

Parameters
----------
ell : ~pyarts3.arts.Vector2
  Ellipsoid (a,b)
lat : ~pyarts3.arts.Numeric
  The latitude

Returns
-------
RADIUS : ~pyarts3.arts.Numeric
  The radius
)");

  geo.def("geodetic2geocentric",
          &geodetic2geocentric,
          "pos"_a,
          "ell"_a,
          R"(Converts geodetic coordinates to geocentric coordinates

This radius is the distance to the polar axis from the reference
ellipsoid surface in the local nadir direction.

See further https://en.wikipedia.org/wiki/Earth_radius#Prime_vertical

Parameters
----------
pos : ~pyarts3.arts.Vector3
  Geodetic position (h,lat,lon)
ell : ~pyarts3.arts.Vector2
  Ellipsoid (a,b)

Returns
-------
pos : ~pyarts3.arts.Vector3
  Geocentric position (r,lat,lon)
)");
}

}  // namespace Python
