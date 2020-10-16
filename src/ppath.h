/* Copyright (C) 2002-2012 Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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

/**
 * @file   ppath.h
 * @author Patrick Eriksson <patrick.eriksson@chalmers.se>
 * @date   2002-05-02
 * 
 * @brief  Propagation path structure and functions.
 * 
 * This file contains the definition of the Ppath structure and the
 * functions in ppath.cc that are of interest elsewhere.
 */

#ifndef ppath_h
#define ppath_h

#include "agenda_class.h"
#include "array.h"
#include "arts.h"
#include "interpolation.h"
#include "matpackI.h"
#include "mystring.h"

/*===========================================================================
  === The Ppath structure
  ===========================================================================*/

/** The structure to describe a propagation path and releated quantities.
 * 
 *  The fields of the structure are described more in detail inside the ARTS
 *  user guide (AUG).
 */
struct Ppath {
  /** Atmospheric dimensionality */
  Index dim;
  /** Number of points describing the ppath */
  Index np;
  /** The propagation path constant (only used for 1D) */
  Numeric constant;
  /** Radiative background */
  String background;
  /** Start position */
  Vector start_pos;
  /** Start line-of-sight */
  Vector start_los;
  /** Length between sensor and atmospheric boundary */
  Numeric start_lstep;
  /** The distance between start pos and the last position in pos */
  Matrix pos;
  /** Line-of-sight at each ppath point */
  Matrix los;
  /** Radius of each ppath point */
  Vector r;
  /** The length between ppath points */
  Vector lstep;
  /** End position */
  Vector end_pos;
  /** End line-of-sight */
  Vector end_los;
  /** The distance between end pos and the first position in pos */
  Numeric end_lstep;
  /** The real part of the refractive index at each path position */
  Vector nreal;
  /** The group index of refraction */  
  Vector ngroup;
  /** Index position with respect to the pressure grid */
  ArrayOfGridPos gp_p;
  /** Index position with respect to the latitude grid */
  ArrayOfGridPos gp_lat;
  /** Index position with respect to the longitude grid */
  ArrayOfGridPos gp_lon;
};

/** An array of propagation paths. */
typedef Array<Ppath> ArrayOfPpath;

/** Size of north and south poles
 * 
 * Latitudes with an absolute value > POLELAT are considered to be on
 * the south or north pole for 3D.
 */
const Numeric POLELAT = 90 - 1e-8;

/** Width of zenith and nadir directions
 *
 * This variable defines how much zenith and azimuth angles can
 * deviate from 0, 90 and 180 degrees, but still be treated to be 0,
 * 90 or 180.  For example, an azimuth angle of 180-0.999*ANGTOL will
 * be treated as a strictly southward observation.  However, the
 * angles are not allowed to go outside their defined range.  This
 * means, for example, that values above 180 are never allowed.
 */
const Numeric ANGTOL = 1e-6;

/*===========================================================================
  === Functions from ppath.cc
  ===========================================================================*/

/** Adds up zenith and azimuth angles

   Adds (dza,daa) to (za0,aa0), assuming that a unit changes in za and aa are
   equal where (dza,daa)=(0,0).

   @param[out]   za     End zenith angle
   @param[out]   aa     End azimuth angle
   @param[in]    za0    Start zenith angle
   @param[in]    aa0    Start azimuth angle
   @param[in]    dza    Change in zenith angle
   @param[in]    daa    Change in azimuth angle

   @author Patrick Eriksson
   @date   2018-12-19
 */
void add_za_aa(Numeric& za,
               Numeric& aa,
               const Numeric& za0,
               const Numeric& aa0,
               const Numeric& dza,
               const Numeric& daa);

/** Converts a cartesian directional vector to zenith and azimuth

   This function and the sister function cart2zaaa handle
   transformation of line-of-sights. This in contrast to the sph/poslos
   functions that handles positions, or combinations of positions and
   line-of-sight.

   The cartesian coordinate system used for these two functions can 
   be defined as
    z : za = 0
    x : za=90, aa=0
    y : za=90, aa=90

   @param[out]   za    LOS zenith angle at observation position.
   @param[out]   aa    LOS azimuth angle at observation position.
   @param[in]    dx    x-part of LOS unit vector.
   @param[in]    dy    y-part of LOS unit vector.
   @param[in]    dz    z-part of LOS unit vector.

   @author Patrick Eriksson
   @date   2009-10-02
 */
void cart2zaaa(Numeric& za,
               Numeric& aa,
               const Numeric& dx,
               const Numeric& dy,
               const Numeric& dz);

/** Converts zenith and azimuth angles to a cartesian unit vector.

   This function and the sister function cart2zaaa handles
   transformation of line-of-sights. This in contrast to the sph/poslos
   functions that handles positions, or combinations of positions and
   line-of-sight.

   The cartesian coordinate system used for these two functions can 
   be defined as
    z : za = 0
    x : za=90, aa=0
    y : za=90, aa=90

   @param[out]   dx    x-part of LOS unit vector.
   @param[out]   dy    y-part of LOS unit vector.
   @param[out]   dz    z-part of LOS unit vector.
   @param[in]    za    LOS zenith angle at observation position.
   @param[in]    aa    LOS azimuth angle at observation position.

   @author Patrick Eriksson
   @date   2009-10-02
 */
void zaaa2cart(Numeric& dx,
               Numeric& dy,
               Numeric& dz,
               const Numeric& za,
               const Numeric& aa);

/** Converts ENU unit vector vector to zenith and azimuth

   This function and the sister function enu2zaaa handles transformation of
   line-of-sights, from and to ENU (east-north-up). The ENU vector is
   normalised to have length 1.

   @param[out]   za    LOS zenith angle at observation position.
   @param[out]   aa    LOS azimuth angle at observation position.
   @param[in]    de    e-part of LOS unit vector.
   @param[in]    dn    n-part of LOS unit vector.
   @param[in]    du    u-part of LOS unit vector.

   @author Patrick Eriksson
   @date   2020-09-17
 */
void enu2zaaa(Numeric& za,
              Numeric& aa,
              const Numeric& de,
              const Numeric& dn,
              const Numeric& du);

/** Converts zenith and azimuth angles to ENU unit vector.

   This function and the sister function enu2zaaa handles transformation of
   line-of-sights, from and to ENU (east-north-up). The ENU vector is
   normalised to have length 1.

   @param[out]   de    e-part of LOS unit vector.
   @param[out]   dn    n-part of LOS unit vector.
   @param[out]   du    u-part of LOS unit vector.
   @param[in]    za    LOS zenith angle at observation position.
   @param[in]    aa    LOS azimuth angle at observation position.

   @author Patrick Eriksson
   @date   2020-09-17
 */
void zaaa2enu(Numeric& de,
              Numeric& dn,
              Numeric& du,
              const Numeric& za,
              const Numeric& aa);

/** Takes the difference of zenith and azimuth angles

   Takes the difference between a set of angles (za,aa) and a reference
   direction (za0,aa0). That is, this function is the "inverse" of *add_za_aa*.

   @param[out]    dza   Change in zenith angle
   @param[out]    daa   Change in azimuth angle
   @param[in]     za0   Zenith angle of reference direction
   @param[in]     aa0   Azimuth angle  of reference direction
   @param[in]     za    Zenith angle of second direction
   @param[in]     aa    Azimuth angle  of second direction

   @author Patrick Eriksson
   @date   2018-12-19
 */
void diff_za_aa(Numeric& dza,
                Numeric& daa,
                const Numeric& za0,
                const Numeric& aa0,
                const Numeric& za,
                const Numeric& aa);

/** Throws an error if ppath altitudes not are strictly increasing or decreasing.

   @param[in]   ppath   Propagation path structure.

   @author Patrick Eriksson
   @date   2018-03-07
 */
void error_if_limb_ppath(const Ppath& ppath);

/** Identifies the tangent point of a propagation path

   The tangent points is defined as the point with the lowest altitude. 

   The index of the tangent point is determined. If no tangent point is found,
   the index is set to -1. 

   @param[out]   it      Index of tangent point
   @param[in]    ppath   Propagation path structure.

   @author Patrick Eriksson
   @date   2012-04-07
 */
void find_tanpoint(Index& it, const Ppath& ppath);

/** Determines ppath position just below an altitude

   Returns -1 if the selected altitude is not passed.

   @param[in]    ppath   Propagation path structure.
   @param[in]    alt     Altitude   

   @return   Index of found altitude

   @author Richard Larsson
   @date   20??-??-??
 */
Index first_pos_before_altitude(const Ppath& p, const Numeric& alt);

/** Calculates the propagation path constant for pure geometrical calculations.

   Both positive and negative zenith angles are handled.

   @param[in]   r      Radius of the sensor position.
   @param[in]   za     Zenith angle of the sensor line-of-sight.

   @return   Path constant.

   @author Patrick Eriksson
   @date   2002-05-17
 */
Numeric geometrical_ppc(const Numeric& r, const Numeric& za);

/** Calculates the latitude for a given zenith angle along a geometrical
   propagation path.

   Positive and negative zenith angles are handled. A positive zenith angle
   means a movement towards higher latitudes.

   @param[in]   za0    The zenith angle of the starting point.
   @param[in]   lat0   The latitude of the starting point.
   @param[in]   za     The zenith angle of the second point.

   @return   The latitude of the second point.

   @author Patrick Eriksson
   @date   2002-05-17
 */
Numeric geompath_lat_at_za(const Numeric& za0,
                           const Numeric& lat0,
                           const Numeric& za);

/** Calculates the length from the tangent point for the given radius.

   The tangent point is either real or imaginary depending on the zenith
   angle of the sensor. See geometrical_tangent_radius.

   @param[in]   ppc    Propagation path constant.
   @param[in]   r      Radius of the point of concern.

   @return  Length along the path from the tangent point. Always >= 0.

   @author Patrick Eriksson
   @date   2002-05-20
 */
Numeric geompath_l_at_r(const Numeric& ppc, const Numeric& r);

/** Calculates the radius for a distance from the tangent point.

   The tangent point is either real or imaginary depending on the zenith
   angle of the sensor. See geometrical_tangent_radius.

   @param[in]   ppc    Propagation path constant.
   @param[in]   l      Length from the tangent point (positive or negative).

   @return         Radius. 

   @author Patrick Eriksson
   @date   2002-05-20
 */
Numeric geompath_r_at_l(const Numeric& ppc, const Numeric& l);

/** Calculates the zenith angle for a given radius along a geometrical 
   propagation path.

   For downlooking cases, the two points must be on the same side of 
   the tangent point.

   Both positive and negative zenith angles are handled.

   @param[in]   ppc    Propagation path constant.
   @param[in]   a_za   A zenith angle along the path on the same side of the 
                       tangent point as the point of interest.  
   @param[in]   r      Radius of the point of interest.

   @return   Zenith angle at the point of interest.

   @author Patrick Eriksson
   @date   2002-05-17
 */
Numeric geompath_za_at_r(const Numeric& ppc,
                         const Numeric& a_za,
                         const Numeric& r);

/** Determines if a line-of-sight is downwards compared to the angular tilt
   of the surface or a pressure level.

   For example, this function can be used to determine if the line-of-sight
   goes into the surface for a starting point exactly on the surface radius.
  
   As the radius of the surface and pressure levels varies as a function of
   latitude, it is not clear if a zenith angle of 90 is above or below e.g.
   the surface.
 
   @param[in]   za     Zenith angle of line-of-sight.
   @param[in]   tilt   Angular tilt of the surface or the pressure level (as
                       returned by plevel_angletilt)

   @return   Boolean that is true if LOS is downwards.

   @author Patrick Eriksson
   @date   2002-06-03
 */
bool is_los_downwards(const Numeric& za, const Numeric& tilt);

/** Calculates the angular tilt of the surface or a pressure level.

   Note that the tilt value is a local value. The tilt for a constant
   slope value, is different for different radii.

   @param[in]    r    The radius for the level at the point of interest.
   @param[in]    c1   The radial slope, as returned by e.g. plevel_slope_2d.

   @return   The angular tilt.

   @author Patrick Eriksson
   @date   2002-06-03
 */
Numeric plevel_angletilt(const Numeric& r, const Numeric& c);

/** Calculates the radial slope of the surface or a pressure level for 2D.

   The radial slope is here the derivative of the radius with respect to the
   latitude. The unit is accordingly m/degree. 

   Note that the radius is defined to change linearly between grid points, and
   the slope is constant between to points of the latitude grid. The radius can
   inside the grid range be expressed as r = r0(lat0) + c1*(lat-lat0) .

   Note also that the slope is always calculated with respect to increasing
   latitudes, independently of the zenith angle. The zenith angle is
   only used to determine which grid range that is of interest when the
   position is exactly on top of a grid point. 

   @param[out]  c1             The radial slope [m/degree]
   @param[in]   lat_grid       The latitude grid.
   @param[in]   refellipsoid   As the WSV with the same name.
   @param[in]   z_surf         Geometrical altitude of the surface, or the pressure
                               level of interest, for the latitide dimension
   @param[in]   gp             Latitude grid position for the position of interest
   @param[in]   za             LOS zenith angle.

   @author Patrick Eriksson
   @date   2002-06-03
 */
void plevel_slope_2d(Numeric& c1,
                     ConstVectorView lat_grid,
                     ConstVectorView refellipsoid,
                     ConstVectorView z_surf,
                     const GridPos& gp,
                     const Numeric& za);

/** Calculates the radial slope of the surface or a pressure level for 3D.

   For 2D the radius can be expressed as r = r0 + c1*dalpha, where alpha
   is the latitude. The radius is here for 3D expressed as a second order
   polynomial: r = r0 + c1*dalpha + c2*dalpha^2, where alpha is the angular
   change (in degrees) along the great circle along the given azimuth angle.

   @param[out]  c1     See above. Unit is m/degree.
   @param[out]  c2     See above. Unit is m/degree^2.
   @param[in]   lat1   Lower latitude of grid cell.
   @param[in]   lat3   Upper latitude of grid cell.
   @param[in]   lon5   Lower longitude of grid cell.
   @param[in]   lon6   Upper longitude of grid cell.
   @param[in]   r15    Radius at crossing of *lat1* and *lon5*.
   @param[in]   r35    Radius at crossing of *lat3* and *lon5*.
   @param[in]   r36    Radius at crossing of *lat3* and *lon6*.
   @param[in]   r16    Radius at crossing of *lat1* and *lon6*.
   @param[in]   lat    Latitude for which slope shall be determined.
   @param[in]   lon    Longitude for which slope shall be determined.
   @param[in]   aa     Azimuth angle for which slope shall be determined.

   @author Patrick Eriksson
   @date   2002-12-30
 */
void plevel_slope_3d(Numeric& c1,
                     Numeric& c2,
                     ConstVectorView lat_grid,
                     ConstVectorView lon_grid,
                     ConstVectorView refellipsoid,
                     ConstMatrixView z_surf,
                     const GridPos& gp_lat,
                     const GridPos& gp_lon,
                     const Numeric& aa);

/** This is the core for the WSM ppathStepByStep.

   This function takes mainly the same input as ppathStepByStep (that 
   is, those input arguments are the WSV with the same name).

   @param[in] ws                 Current Workspace
   @param[in] ppath              Output: A Ppath structure
   @param[in] ppath_step_agenda  As the WSM with the same name.
   @param[in] atmosphere_dim     The atmospheric dimensionality.
   @param[in] p_grid             The pressure grid.
   @param[in] lat_grid           The latitude grid.
   @param[in] lon_grid           The longitude grid.
   @param[in] z_field            As the WSM with the same name.
   @param[in] f_grid             As the WSM with the same name.
   @param[in] refellipsoid       As the WSM with the same name.
   @param[in] z_surface          Surface altitude.
   @param[in] cloudbox_on        Flag to activate the cloud box.
   @param[in] cloudbox_limits    Index limits of the cloud box.
   @param[in] rte_pos            The position of the sensor.
   @param[in] rte_los            The line-of-sight of the sensor.
   @param[in] ppath_lmax         As the WSM with the same name.
   @param[in] ppath_lraytrace    As the WSM with the same name.
   @param[in] ppath_inside_cloudbox_do  As the WSM with the same name.

   @author Patrick Eriksson
   @date   2003-01-08
 */
void ppath_calc(Workspace& ws,
                Ppath& ppath,
                const Agenda& ppath_step_agenda,
                const Index& atmosphere_dim,
                const Vector& p_grid,
                const Vector& lat_grid,
                const Vector& lon_grid,
                const Tensor3& z_field,
                const Vector& f_grid,
                const Vector& refellipsoid,
                const Matrix& z_surface,
                const Index& cloudbox_on,
                const ArrayOfIndex& cloudbox_limits,
                const Vector& rte_pos,
                const Vector& rte_los,
                const Numeric& ppath_lmax,
                const Numeric& ppath_lraytrace,
                const bool& ppath_inside_cloudbox_do,
                const Verbosity& verbosity);

/** Copy the content in ppath2 to ppath1.

   The ppath1 structure must be allocated before calling the function. The
   structure can be allocated to hold more points than found in ppath2.
   The data of ppath2 is placed in the first positions of ppath1.

   @param[in]   ppath1    Output: Ppath structure.
   @param[in]   ppath2    The ppath structure to be copied.
   @param[in]   ncopy     Number of points in ppath2 to copy. If set to negative,
                      the number is set to ppath2.np. 

   @author Patrick Eriksson
   @date   2002-07-03
 */
void ppath_copy(Ppath& ppath1, const Ppath& ppath2, const Index& ncopy);

/** Initiates a Ppath structure to hold the given number of points.

   The background field is set to background case 0. The constant field is set
   to -1. The refraction field is set to 0.

   The length of the lstep field is set to np-1.

   @param[in]   ppath            Output: A Ppath structure.
   @param[in]   atmosphere_dim   The atmospheric dimensionality.
   @param[in]   np               Number of points of the path.

   @author Patrick Eriksson
   @date   2002-05-17
 */
void ppath_init_structure(Ppath& ppath,
                          const Index& atmosphere_dim,
                          const Index& np);

/** Sets the background field of a Ppath structure.

   The different background cases have a number coding to simplify a possible
   change of the strings and checking of the what case that is valid.

   The case numbers are:                    <br>
      0. Unvalid.                           <br>
      1. Space.                             <br>
      2. The surface.                       <br>
      3. The cloud box boundary.            <br>
      4. The interior of the cloud box.       

   @param[in]   ppath            Output: A Ppath structure.
   @param[in]   case_nr          Case number (see above)

   @author Patrick Eriksson
   @date   2002-05-17
 */
void ppath_set_background(Ppath& ppath, const Index& case_nr);

/** Initiates a Ppath structure for calculation of a path with *ppath_step*.

   The function performs two main tasks. As mentioned above, it initiates
   a Ppath structure (a), but it also checks that the end point of the path is
   at an allowed location (b).

   (a): The Ppath structure is set to hold the position and LOS of the last
   point of the path inside the atmosphere. This point is either the
   sensor position, or the point where the path leaves the model atmosphere.
   If the path is totally outside the atmosphere, no point is put into the
   structure. If the (practical) end and start points are identical, such
   as when the sensor is inside the cloud box, the background field is set.

   (b): If it is found that the end point of the path is at an illegal position
   a detailed error message is given. Not allowed cases are: <br>  
      1. The sensor is placed below surface level. <br> 
      2. For 2D and 3D, the path leaves the model atmosphere at a latitude or
         longitude end face. <br> 
      3. For 2D and 3D, the path is totally outside the atmosphere and the 
         latitude and longitude of the tangent point is outside the range of
         the corresponding grids. 

   All input variables are identical with the WSV with the same name.
   The output variable is here called ppath for simplicity, but is in
   fact *ppath_step*.

   @param[in]   ppath             Output: A Ppath structure.
   @param[in]   atmosphere_dim    The atmospheric dimensionality.
   @param[in]   p_grid            The pressure grid.
   @param[in]   lat_grid          The latitude grid.
   @param[in]   lon_grid          The longitude grid.
   @param[in]   z_field           The field of geometrical altitudes.
   @param[in]   refellipsoid      As the WSV with the same name.
   @param[in]   z_surface         Surface altitude.
   @param[in]   cloudbox_on       Flag to activate the cloud box.
   @param[in]   cloudbox_limits   Index limits of the cloud box.
   @param[in]   ppath_inside_cloudbox_do   As the WSV with the same name.
   @param[in]   rte_pos           The position of the sensor.
   @param[in]   rte_los           The line-of-sight of the sensor.

   @author Patrick Eriksson
   @date   2002-05-17
 */
void ppath_start_stepping(Ppath& ppath,
                          const Index& atmosphere_dim,
                          ConstVectorView p_grid,
                          ConstVectorView lat_grid,
                          ConstVectorView lon_grid,
                          ConstTensor3View z_field,
                          ConstVectorView refellipsoid,
                          ConstMatrixView z_surface,
                          const Index& cloudbox_on,
                          const ArrayOfIndex& cloudbox_limits,
                          const bool& outside_cloudbox,
                          ConstVectorView rte_pos,
                          ConstVectorView rte_los,
                          const Verbosity& verbosity);

/** Calculates 1D geometrical propagation path steps.

   This is the core function to determine 1D propagation path steps by pure
   geometrical calculations. Path points are included for crossings with the
   grids, tangent points and points of intersection with the surface. In
   addition, points are included in the propgation path to ensure that the
   distance along the path between the points does not exceed the selected
   maximum length (lmax). If lmax is <= 0, this means that no length criterion
   shall be applied.

   Note that the input variables are here compressed to only hold data for
   a 1D atmosphere. For example, z_field is z_field(:,0,0).

   For more information read the chapter on propagation paths in AUG.

   @param[in]   ppath             Output: A Ppath structure.
   @param[in]   z_field           Geometrical altitudes corresponding to p_grid.
   @param[in]   refellipsoid      As the WSV with the same name.
   @param[in]   z_surface         Surface altitude.
   @param[in]   lmax              Maximum allowed length between the path points.

   @author Patrick Eriksson
   @date   2002-05-20
 */
void ppath_step_geom_1d(Ppath& ppath,
                        ConstVectorView z_field,
                        ConstVectorView refellipsoid,
                        const Numeric& z_surface,
                        const Numeric& lmax);

/** Calculates 2D geometrical propagation path steps.

   Works as the same function for 1D, despite that some input arguments are
   of different type.

   @param[in]   ppath             Output: A Ppath structure.
   @param[in]   lat_grid          Latitude grid.
   @param[in]   z_field           Geometrical altitudes
   @param[in]   refellipsoid      As the WSV with the same name.
   @param[in]   z_surface         Surface altitudes.
   @param[in]   lmax              Maximum allowed length between the path points.

   @author Patrick Eriksson
   @date   2002-07-03
 */
void ppath_step_geom_2d(Ppath& ppath,
                        ConstVectorView lat_grid,
                        ConstMatrixView z_field,
                        ConstVectorView refellipsoid,
                        ConstVectorView z_surface,
                        const Numeric& lmax);

/** Calculates 3D geometrical propagation path steps.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   @param[in]   ppath             Output: A Ppath structure.
   @param[in]   lat_grid          Latitude grid.
   @param[in]   lon_grid          Longitude grid.
   @param[in]   z_field           Geometrical altitudes
   @param[in]   refellipsoid      As the WSV with the same name.
   @param[in]   z_surface         Surface altitudes.
   @param[in]   lmax              Maximum allowed length between the path points.

   @author Patrick Eriksson
   @date   2002-12-30
 */
void ppath_step_geom_3d(Ppath& ppath,
                        ConstVectorView lat_grid,
                        ConstVectorView lon_grid,
                        ConstTensor3View z_field,
                        ConstVectorView refellipsoid,
                        ConstMatrixView z_surface,
                        const Numeric& lmax);

/** Calculates 1D propagation path steps including effects of refraction.

   This function works as the function *ppath_step_geom_1d* but considers
   also refraction. The upper length of the ray tracing steps is set by
   the argument *lraytrace*. This argument controls only the internal
   calculations. The maximum distance between the path points is still
   determined by *lmax*.

   @param[in,out]   ws            Current Workspace
   @param[out]  ppath             A Ppath structure.
   @param[in]   p_grid            Pressure grid.
   @param[in]   z_field           As the WSV with the same name.
   @param[in]   t_field           As the WSV with the same name.
   @param[in]   vmr_field         As the WSV with the same name.
   @param[in]   f_grid            As the WSV with the same name.
   @param[in]   refellipsoid      As the WSV with the same name.
   @param[in]   z_surface         Surface altitude (1D).
   @param[in]   lmax              Maximum allowed length between the path points.
   @param[in]   refr_index_air_agenda The WSV with the same name.
   @param[in]   rtrace_method     String giving which ray tracing method to use.
                              See the function for options.
   @param[in]   lraytrace         Maximum allowed length for ray tracing steps.

   @author Patrick Eriksson
   @date   2002-11-26
 */
void ppath_step_refr_1d(Workspace& ws,
                        Ppath& ppath,
                        ConstVectorView p_grid,
                        ConstTensor3View z_field,
                        ConstTensor3View t_field,
                        ConstTensor4View vmr_field,
                        ConstVectorView f_grid,
                        ConstVectorView refellipsoid,
                        const Numeric& z_surface,
                        const Numeric& lmax,
                        const Agenda& refr_index_agenda,
                        const String& rtrace_method,
                        const Numeric& lraytrace);

/** Calculates 2D propagation path steps, with refraction, using a simple
   and fast ray tracing scheme.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   @param[in,out]   ws            Current Workspace
   @param[out]  ppath             A Ppath structure.
   @param[in]   p_grid            Pressure grid.
   @param[in]   lat_grid          Latitude grid.
   @param[in]   z_field           As the WSV with the same name.
   @param[in]   t_field           As the WSV with the same name.
   @param[in]   vmr_field         As the WSV with the same name.
   @param[in]   f_grid            As the WSV with the same name.
   @param[in]   refellipsoid      As the WSV with the same name.
   @param[in]   z_surface         Surface altitudes.
   @param[in]   lmax              Maximum allowed length between the path points.
   @param[in]   refr_index_air_agenda The WSV with the same name.
   @param[in]   rtrace_method     String giving which ray tracing method to use.
                              See the function for options.
   @param[in]   lraytrace         Maximum allowed length for ray tracing steps.

   @author Patrick Eriksson
   @date   2002-12-02
 */
void ppath_step_refr_2d(Workspace& ws,
                        Ppath& ppath,
                        ConstVectorView p_grid,
                        ConstVectorView lat_grid,
                        ConstTensor3View z_field,
                        ConstTensor3View t_field,
                        ConstTensor4View vmr_field,
                        ConstVectorView f_grid,
                        ConstVectorView refellipsoid,
                        ConstVectorView z_surface,
                        const Numeric& lmax,
                        const Agenda& refr_index_agenda,
                        const String& rtrace_method,
                        const Numeric& lraytrace);

/** Calculates 3D propagation path steps, with refraction, using a simple
   and fast ray tracing scheme.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   @param[in,out]   ws            Current Workspace
   @param[out]  ppath             A Ppath structure.
   @param[in]   p_grid            Pressure grid.
   @param[in]   lat_grid          Latitude grid.
   @param[in]   lon_grid          Longitude grid.
   @param[in]   z_field           Geometrical altitudes.
   @param[in]   t_field           Atmospheric temperatures.
   @param[in]   vmr_field         VMR values.
   @param[in]   f_grid            As the WSV with the same name.
   @param[in]   refellipsoid      As the WSV with the same name.
   @param[in]   z_surface         Surface altitudes.
   @param[in]   lmax              Maximum allowed length between the path points.
   @param[in]   refr_index_air_agenda The WSV with the same name.
   @param[in]   rtrace_method     String giving which ray tracing method to use.
                              See the function for options.
   @param[in]   lraytrace         Maximum allowed length for ray tracing steps.

   @author Patrick Eriksson
   @date   2003-01-08
 */
void ppath_step_refr_3d(Workspace& ws,
                        Ppath& ppath,
                        ConstVectorView p_grid,
                        ConstVectorView lat_grid,
                        ConstVectorView lon_grid,
                        ConstTensor3View z_field,
                        ConstTensor3View t_field,
                        ConstTensor4View vmr_field,
                        ConstVectorView f_grid,
                        ConstVectorView refellipsoid,
                        ConstMatrixView z_surface,
                        const Numeric& lmax,
                        const Agenda& refr_index_agenda,
                        const String& rtrace_method,
                        const Numeric& lraytrace);

/** Returns the case number for the radiative background.

   See further the function *ppath_set_background*.

   @param[in]   ppath   A Ppath structure.

   @return   The case number.

   @author Patrick Eriksson
   @date   2002-05-17
 */
Index ppath_what_background(const Ppath& ppath);

/** Resolves which longitude angle that shall be used.

   Longitudes are allowed to vary between -360 and 360 degress, while the
   inverse trigonomtric functions returns values between -180 and 180.
   This function determines if the longitude shall be shifted -360 or
   +360 to fit the longitudes set by the user.
   
   The argument *lon* as input is a value calculated by some inverse
   trigonometric function. The arguments *lon5* and *lon6* are the
   lower and upper limit for the probable range for *lon*. The longitude
   *lon* will be shifted with -360 or +360 degrees if lon is significantly
   outside *lon5* and *lon6*. No error is given if it is not possible to
   obtain a value between *lon5* and *lon6*. 

   @param[in,out]   lon    Longitude, possible shifted when returned.
   @param[in]       lon5   Lower limit of probable range for lon.
   @param[in]       lon6   Upper limit of probable range for lon

   @author Patrick Eriksson
   @date   2003-01-05
 */
void resolve_lon(Numeric& lon, const Numeric& lon5, const Numeric& lon6);

/** Converts zenith and azimuth angles to a cartesian unit vector.

   This function and the sister function cart2zaaa handles
   transformation of line-of-sights. This in contrast to the sph/poslos
   functions that handles positions, or combinations of positions and
   line-of-sight.

   The cartesian coordinate system used for these two functions can 
   be defined as
    z : za = 0
    x : za=90, aa=0
    y : za=90, aa=90

   @param[out]  dx    x-part of LOS unit vector.
   @param[out]  dy    y-part of LOS unit vector.
   @param[out]  dz    z-part of LOS unit vector.
   @param[in]   za    LOS zenith angle at observation position.
   @param[in]   aa    LOS azimuth angle at observation position.

   @author Patrick Eriksson
   @date   2009-10-02
 */
void zaaa2cart(Numeric& dx,
               Numeric& dy,
               Numeric& dz,
               const Numeric& za,
               const Numeric& aa);

#endif  // ppath_h
