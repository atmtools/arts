/**
  @file   special_interp.h
  @author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  @date   2002-11-14   
  
  @brief  Header file for special_interp.cc.
 */

#ifndef special_interp_h
#define special_interp_h

#include "gridded_fields.h"
#include "interpolation.h"
#include "jacobian.h"

/*===========================================================================
  === Interpolation functions for atmospheric grids, fields and surfaces
  ===========================================================================*/

/** Converts atmospheric grid positions to weights for interpolation of an
    atmospheric field.

    The function is intended for "blue" interpolation, that is, interpolation
    for a set of positions. 

    The output matrix for interpolation weights are resized inside the
    function.

    The input atmospheric grids are checked to be consistent.

    @param[out]  itw                Interpolation weights.
    @param[in]   atmosphere_dim     As the WSV with the same name.
    @param[in]   gp_p               Pressure grid positions.
    @param[in]   gp_lat             Latitude grid positions.
    @param[in]   gp_lon             Longitude grid positions.

    @author Patrick Eriksson 
    @date   2002-11-13
 */
void interp_atmfield_gp2itw(Matrix& itw,
                            const Index& atmosphere_dim,
                            const ArrayOfGridPos& gp_p,
                            const ArrayOfGridPos& gp_lat,
                            const ArrayOfGridPos& gp_lon);

/** Interpolates an atmospheric field with pre-calculated weights by
    interp_atmfield_gp2itw.

    The function performs the interpolation for a number of positions. The
    return variable (x) is accordingly a vector. The vector must be set to
    have the same length as the grid position arrays before calling the
    function. 

    The input atmospheric field is checked to be consistent with the 
    *atmosphere_dim*, *p_grid*, *lat_grid* and *lon_grid*. The length of
    the grid position arrays are asserted to be the identical, or for 
    dimensions not used, that the length is zero.

    @param[out]  x                  Values obtained by the interpolation.
    @param[in]   atmosphere_dim     As the WSV with the same name.
    @param[in]   x_field            The atmospheric field to be interpolated.
    @param[in]   gp_p               Pressure grid positions.
    @param[in]   gp_lat             Latitude grid positions.
    @param[in]   gp_lon             Longitude grid positions.
    @param[in]   itw                Interpolation weights from 
                                    interp_atmfield_gp2itw.

    @author Patrick Eriksson 
    @date   2002-11-13
 */
void interp_atmfield_by_itw(VectorView x,
                            const Index& atmosphere_dim,
                            ConstTensor3View x_field,
                            const ArrayOfGridPos& gp_p,
                            const ArrayOfGridPos& gp_lat,
                            const ArrayOfGridPos& gp_lon,
                            ConstMatrixView itw);

/** Interpolates an atmospheric field given the grid positions.

    The function performs the interpolation for a number of positions. The
    return variable (x) is accordingly a vector. The vector must be set to
    have the same length as the grid position arrays before calling the
    function. 

    The input atmospheric field is checked to be consistent with the 
    *atmosphere_dim*, *p_grid*, *lat_grid* and *lon_grid*. The length of
    the grid position arrays are asserted to be the identical, or for 
    dimensions not used, that the length is zero.

    There is also a return version of this function.

    @param[out]  x                  Values obtained by the interpolation.
    @param[in]   atmosphere_dim     As the WSV with the same name.
    @param[in]   x_field            The atmospheric field to be interpolated.
    @param[in]   gp_p               Pressure grid positions.
    @param[in]   gp_lat             Latitude grid positions.
    @param[in]   gp_lon             Longitude grid positions.

    @author Patrick Eriksson 
    @date   2002-11-13
 */
void interp_atmfield_by_gp(VectorView x,
                           const Index& atmosphere_dim,
                           ConstTensor3View x_field,
                           const ArrayOfGridPos& gp_p,
                           const ArrayOfGridPos& gp_lat,
                           const ArrayOfGridPos& gp_lon);

/** Interpolates an atmospheric field given the grid positions.

    The function performs the interpolation for a number of positions. The
    return variable (x) is accordingly a vector. The vector must be set to
    have the same length as the grid position arrays before calling the
    function. 

    The input atmospheric field is checked to be consistent with the 
    *atmosphere_dim*, *p_grid*, *lat_grid* and *lon_grid*. The length of
    the grid position arrays are asserted to be the identical, or for 
    dimensions not used, that the length is zero.

    There is also a vector version of this function.

    @param[in]   atmosphere_dim     As the WSV with the same name.
    @param[in]   x_field            The atmospheric field to be interpolated.
    @param[in]   gp_p               Pressure grid positions.
    @param[in]   gp_lat             Latitude grid positions.
    @param[in]   gp_lon             Longitude grid positions.

    @return   Value obtained by the interpolation.

    @author Patrick Eriksson 
    @date   2002-11-13
 */
Numeric interp_atmfield_by_gp(const Index& atmosphere_dim,
                              ConstTensor3View x_field,
                              const GridPos& gp_p = {0, {0, 1}},
                              const GridPos& gp_lat = {0, {0, 1}},
                              const GridPos& gp_lon = {0, {0, 1}});

/** Converts atmospheric a grid position to weights for interpolation of a
    field defined ONLY inside the cloudbox.

    That is, as interp_atmfield_gp2itw, but for cloudbox only variables.

    The input grid position shall be with respect to total grids. If grid
    positions already refer to grid parts inside the cloudbox, you can use 
    interp_atmfield_gp2itw.

    The output grid positions are created by the function, to match the
    cloudbox field, and can be used for later calls of e.g.
    interp_atmfield_by_itw

    @param[out]  itw                Interpolation weights. Vector must be 
                                    given correct size before call of function.
    @param[in]   gp_p_out           Output: Pressure cloudbox grid position.
    @param[in]   gp_lat_out         Output: Latitude cloudbox grid position.
    @param[in]   gp_lon_out         Output: Longitude cloudbox grid position.
    @param[in]   gp_p_in            Pressure grid position.
    @param[in]   gp_lat_in          Latitude grid position.
    @param[in]   gp_lon_in          Longitude grid position.
    @param[in]   atmosphere_dim     As the WSV with the same name.
    @param[in]   cloudbox_limits    As the WSV with the same name.

    @author Patrick Eriksson 
    @date   2010-02-12
 */
void interp_cloudfield_gp2itw(VectorView itw,
                              GridPos& gp_p_out,
                              GridPos& gp_lat_out,
                              GridPos& gp_lon_out,
                              const GridPos& gp_p_in,
                              const GridPos& gp_lat_in,
                              const GridPos& gp_lon_in,
                              const Index& atmosphere_dim,
                              const ArrayOfIndex& cloudbox_limits);

/** Converts atmospheric grid positions to weights for interpolation of a
    surface-type variable.

    The function is intended for "blue" interpolation, that is, interpolation
    for a set of positions. 

    The output matrix for interpolation weights are resized inside the
    function.

    The input atmospheric grids are checked to be consistent.

    @param[out]  itw                Interpolation weights.
    @param[in]   atmosphere_dim     As the WSV with the same name.
    @param[in]   gp_lat             Latitude grid positions.
    @param[in]   gp_lon             Longitude grid positions.

    @author Patrick Eriksson 
    @date   2002-11-13
 */
void interp_atmsurface_gp2itw(Matrix& itw,
                              const Index& atmosphere_dim,
                              const ArrayOfGridPos& gp_lat,
                              const ArrayOfGridPos& gp_lon);

/** Interpolates a surface-type variable with pre-calculated weights by
    interp_atmsurface_gp2itw.

    The function performs the interpolation for a number of positions. The
    return variable (x) is accordingly a vector. The vector must be set to
    have the same length as the grid position arrays before calling the
    function. 

    The input surface-type variable is checked to be consistent with the 
    *atmosphere_dim*, *lat_grid* and *lon_grid*. The length of
    the grid position arrays are asserted to be the identical, or for 
    dimensions not used, that the length is zero.

    @param[out]  x                  Values obtained by the interpolation.
    @param[in]   atmosphere_dim     As the WSV with the same name.
    @param[in]   x_surface          The atmospheric field to be interpolated.
    @param[in]   gp_lat             Latitude grid positions.
    @param[in]   gp_lon             Longitude grid positions.
    @param[in]   itw                Interpolation weights from 
                                    interp_atmsurface_gp2itw.

    @author Patrick Eriksson 
    @date   2002-11-13
 */
void interp_atmsurface_by_itw(VectorView x,
                              const Index& atmosphere_dim,
                              ConstMatrixView x_surface,
                              const ArrayOfGridPos& gp_lat,
                              const ArrayOfGridPos& gp_lon,
                              ConstMatrixView itw);

/** Interpolates a surface-type variable given the grid positions.

    The function performs the interpolation for a number of positions. The
    return variable (x) is accordingly a vector. The vector must be set to
    have the same length as the grid position arrays before calling the
    function. 

    The input surface-type variable is checked to be consistent with the 
    *atmosphere_dim*, *lat_grid* and *lon_grid*. The length of
    the grid position arrays are asserted to be the identical, or for 
    dimensions not used, that the length is zero.

    There is also a return version of this function.

    @param[out]  x                  Values obtained by the interpolation.
    @param[in]   atmosphere_dim     As the WSV with the same name.
    @param[in]   x_surface          The atmospheric field to be interpolated.
    @param[in]   gp_lat             Latitude grid positions.
    @param[in]   gp_lon             Longitude grid positions.

    @author Patrick Eriksson 
    @date   2002-11-13
 */
void interp_atmsurface_by_gp(VectorView x,
                             const Index& atmosphere_dim,
                             ConstMatrixView x_field,
                             const ArrayOfGridPos& gp_lat,
                             const ArrayOfGridPos& gp_lon);

/** Interpolates a surface-type variable given the grid positions.

    The function performs the interpolation for a number of positions. The
    return variable (x) is accordingly a vector. The vector must be set to
    have the same length as the grid position arrays before calling the
    function. 

    The input surface-type variable is checked to be consistent with the 
    *atmosphere_dim*, *lat_grid* and *lon_grid*. The length of
    the grid position arrays are asserted to be the identical, or for 
    dimensions not used, that the length is zero.

    There is also a vecor version of this function.

    @param[in]   atmosphere_dim     As the WSV with the same name.
    @param[in]   x_surface          The atmospheric field to be interpolated.
    @param[in]   gp_lat             Latitude grid positions.
    @param[in]   gp_lon             Longitude grid positions.

    @return  Values obtained by the interpolation.

    @author Patrick Eriksson 
    @date   2002-11-13
 */
Numeric interp_atmsurface_by_gp(const Index& atmosphere_dim,
                                ConstMatrixView x_field,
                                const GridPos& gp_lat,
                                const GridPos& gp_lon);

/** Regrids an atmospheric field, for precalculated grid positions

  The function adopts automatically to *atmosphere_dim*. Grid positions not
  used are ignored, i.e. gp_lat is ignored for atmosphere_dim=1 etc.

  @param[out]     field_new        Field after interpolation.
  @param[in][in]  atmosphere_dim   As the WSV with same name.
  @param[in][in]  field_old        Field to be interpolated.
  @param[in][in]  gp_p             Pressure grid positions.
  @param[in][in]  gp_lat           Latitude grid positions.
  @param[in][in]  gp_lon           Longitude grid positions.

  @author Patrick Eriksson 
  @date   2015-09-09
 */
void regrid_atmfield_by_gp(Tensor3& field_new,
                           const Index& atmosphere_dim,
                           ConstTensor3View field_old,
                           const ArrayOfGridPos& gp_p,
                           const ArrayOfGridPos& gp_lat,
                           const ArrayOfGridPos& gp_lon);

/** Regrids an atmospheric surface, for precalculated grid positions

  The function adopts automatically to *atmosphere_dim*. Grid positions not
  used are ignored, i.e. gp_lat is ignored for atmosphere_dim=1 etc.

  @param[out]     field_new        Field after interpolation.
  @param[in][in]  atmosphere_dim   As the WSV with same name.
  @param[in][in]  field_old        Field to be interpolated.
  @param[in][in]  gp_lat           Latitude grid positions.
  @param[in][in]  gp_lon           Longitude grid positions.

  @author Patrick Eriksson 
  @date   2018-04-12
 */
void regrid_atmsurf_by_gp(Matrix& field_new,
                          const Index& atmosphere_dim,
                          ConstMatrixView field_old,
                          const ArrayOfGridPos& gp_lat,
                          const ArrayOfGridPos& gp_lon);

/** Determines grid positions for regridding of atmospheric fields to retrieval
 *  grids
 *
 * The grid positions arrays are sized inside the function. gp_lat is given
 * length 0 for atmosphere_dim=1 etc.
 *
 * This regridding uses extpolfac=0.
 *
 * @param[out] gp_p                 Pressure grid positions.
 * @param[out] gp_lat               Latitude grid positions.
 * @param[out] gp_lon               Longitude grid positions.
 * @param[in]  rq                   Retrieval quantity structure.
 * @param[in]  atmosphere_dim       As the WSV with same name.
 * @param[in]  p_grid               As the WSV with same name.
 * @param[in]  lat_grid             As the WSV with same name.
 * @param[in]  lon_grid             As the WSV with same name.
 *
 * @author Patrick Eriksson
 * @date   2015-09-09
 */
void get_gp_atmgrids_to_rq(ArrayOfGridPos& gp_p,
                           ArrayOfGridPos& gp_lat,
                           ArrayOfGridPos& gp_lon,
                           const RetrievalQuantity& rq,
                           const Index& atmosphere_dim,
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid);

/** Determines grid positions for regridding of atmospheric surfaces to retrieval
 *  grids
 *
 * The grid positions arrays are sized inside the function. gp_lat is given
 * length 0 for atmosphere_dim=1 etc.
 *
 * This regridding uses extpolfac=0.
 *
 * @param[out] gp_lat               Latitude grid positions.
 * @param[out] gp_lon               Longitude grid positions.
 * @param[in]  rq                   Retrieval quantity structure.
 * @param[in]  atmosphere_dim       As the WSV with same name.
 * @param[in]  lat_grid             As the WSV with same name.
 * @param[in]  lon_grid             As the WSV with same name.
 *
 * @author Patrick Eriksson
 * @date   2018-04-12
 */
void get_gp_atmsurf_to_rq(ArrayOfGridPos& gp_lat,
                          ArrayOfGridPos& gp_lon,
                          const RetrievalQuantity& rq,
                          const Index& atmosphere_dim,
                          const Vector& lat_grid,
                          const Vector& lon_grid);


/** Determines grid positions for regridding of atmospheric fields to retrieval
 *  grids
 *
 * The grid positions arrays are sized inside the function. gp_lat is given
 * length 0 for atmosphere_dim=1 etc.
 *
 * This regridding uses extpolfac=Inf (where Inf is a very large value).
 *
 * Note that the length output arguments (n_p etc.) are for the retrieval grids
 * (not the length of grid positions arrays). n-Lat is set to 1 for
 * atmosphere_dim=1 etc.
 *
 * @param[out] gp_p                 Pressure grid positions.
 * @param[out] gp_lat               Latitude grid positions.
 * @param[out] gp_lon               Longitude grid positions.
 * @param[out] n_p                  Length of retrieval pressure grid.
 * @param[out] n_lat                Length of retrieval lataitude grid.
 * @param[out] n_lon                Length of retrieval longitude grid.
 * @param[in]  rq                   Retrieval quantity structure.
 * @param[in]  atmosphere_dim       As the WSV with same name.
 * @param[in]  p_grid               As the WSV with same name.
 * @param[in]  lat_grid             As the WSV with same name.
 * @param[in]  lon_grid             As the WSV with same name.
 *
 * @author Patrick Eriksson
 * @date   2015-09-09
 */
void get_gp_rq_to_atmgrids(ArrayOfGridPos& gp_p,
                           ArrayOfGridPos& gp_lat,
                           ArrayOfGridPos& gp_lon,
                           Index& n_p,
                           Index& n_lat,
                           Index& n_lon,
                           const ArrayOfVector& ret_grids,
                           const Index& atmosphere_dim,
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid);

/** Determines grid positions for regridding of atmospheric surfaces to retrieval
 *  grids
 *
 * The grid positions arrays are sized inside the function. gp_lat is given
 * length 0 for atmosphere_dim=1 etc.
 *
 * This regridding uses extpolfac=Inf (where Inf is a very large value).
 *
 * Note that the length output arguments (n_p etc.) are for the retrieval grids
 * (not the length of grid positions arrays). n-Lat is set to 1 for
 * atmosphere_dim=1 etc.
 *
 * @param[out] gp_lat               Latitude grid positions.
 * @param[out] gp_lon               Longitude grid positions.
 * @param[out] n_lat                Length of retrieval lataitude grid.
 * @param[out] n_lon                Length of retrieval longitude grid.
 * @param[in]  rq                   Retrieval quantity structure.
 * @param[in]  atmosphere_dim       As the WSV with same name.
 * @param[in]  lat_grid             As the WSV with same name.
 * @param[in]  lon_grid             As the WSV with same name.
 *
 * @author Patrick Eriksson
 * @date   2018-04-12
 */
void get_gp_rq_to_atmgrids(ArrayOfGridPos& gp_lat,
                           ArrayOfGridPos& gp_lon,
                           Index& n_lat,
                           Index& n_lon,
                           const ArrayOfVector& ret_grids,
                           const Index& atmosphere_dim,
                           const Vector& lat_grid,
                           const Vector& lon_grid);


/** Regridding of atmospheric field OEM-type
 *
 * In this regridding infinite extrapolation (by closest neighbour) is allowed,
 * including the case of that the grid length in original field can be 1.
 *
 * @param[out] field_new        New field.
 * @param[in]  atmosphere_dim   Atmospheric dimensionality.
 * @param[in]  field_old        Original field.
 * @param[in]  gp_p             Pressure grid positions.
 * @param[in]  gp_lat           Latitude grid positions.
 * @param[in]  gp_lon           Longitude grid positions.
 * 
 * @author Patrick Eriksson
 * @date   2018-04-12
 */
void regrid_atmfield_by_gp_oem(Tensor3& field_new,
                               const Index& atmosphere_dim,
                               ConstTensor3View field_old,
                               const ArrayOfGridPos& gp_p,
                               const ArrayOfGridPos& gp_lat,
                               const ArrayOfGridPos& gp_lon);

/** Regridding of surface field OEM-type
 *
 * In this regridding infinite extrapolation (by closest neighbour) is allowed,
 * including the case of that the grid length in original field can be 1.
 *
 * @param[out] field_new        New field.
 * @param[in]  atmosphere_dim   Atmospheric dimensionality.
 * @param[in]  field_old        Original field.
 * @param[in]  gp_lat           Latitude grid positions.
 * @param[in]  gp_lon           Longitude grid positions.
 * 
 * @author Patrick Eriksson
 * @date   2018-04-12
 */
void regrid_atmsurf_by_gp_oem(Matrix& field_new,
                              const Index& atmosphere_dim,
                              ConstMatrixView field_old,
                              const ArrayOfGridPos& gp_lat,
                              const ArrayOfGridPos& gp_lon);

/** Converts interpolation weights to pressures.

    The function takes interpolation weights calculated with respect to the 
    vertical dimension, and determines the corresponding pressures. This 
    function can be used when a geometrical altitude is known and the
    pressure for that altitude shall be determined. The interpolation weights
    are then calculated using the geometrical altitudes for the pressure
    levels for the position of concern.

    This can be seen as a 1D "blue" interpolation. That means that the number
    of columns of itw shall be 2.

    @param[out]  p_values   Found pressure values.
    @param[in]   p_grid     As the WSV with the same name.
    @param[in]   gp         Altitude grid positions.
    @param[in]   itw        Interpolation weights

    @author Patrick Eriksson 
    @date   2002-11-13
 */
void itw2p(VectorView p_values,
           ConstVectorView p_grid,
           const ArrayOfGridPos& gp,
           ConstMatrixView itw);

/** Calculates grid positions for pressure values.

   This function works as *gridpos*, but is adapted to handle
   pressure grids. The ARTS defintions result in that pressures shall
   not be interpolated directly, it is the log of the pressure that
   shall be interpolated. This means that if some values shall be
   interpolated to some given pressures, the grid positions shall be
   calculated with this function. The interpolation can then be
   performed as usual.

   @param[out]  gp          Grid position Array.
   @param[in]   old_pgrid   The original pressure grid.
   @param[in]   new_pgrid   The new pressure grid.
   @param[in]   extpolfac   Extrapolation factor. Default value is 0.5,
                            which means that extrapolation of half of the
                            last grid distance is allowed.
                            You don't have to specify this.

   @author Patrick Eriksson
   @date   2003-01-20

   @author Stefan Buehler
   @date   2008-03-03
*/
void p2gridpos(ArrayOfGridPos& gp,
               ConstVectorView old_pgrid,
               ConstVectorView new_pgrid,
               const Numeric& extpolfac = 0.5);

/** Converts a geographical position (rte_pos) to grid positions for p, 
   lat and lon. (field version)

   The function converts the altitude, latitude and longitude in *rte_pos* to
   matching grid positions. The conversion is straightforwatd for latitude and
   longitude. The altitude shall be converted pressure grid position which
   requires an interpolation of z_field.

   Handles 1D, 2D and 3D (gp_lat and gp_lon untouched if not used).

   Note that the function performs several checks of runtime error type.

   @param[out]  gp_p        Pressure grid position.
   @param[out]  gp_lat      Latitude grid position.
   @param[out]  gp_lon      Longitude grid position.
   @param[in]   atmosphere_dim  As the WSV with the same name.
   @param[in]   p_grid      As the WSV with the same name.
   @param[in]   lat_grid    As the WSV with the same name.
   @param[in]   lon_grid    As the WSV with the same name.
   @param[in]   z_field     As the WSV with the same name.
   @param[in]   rte_pos     As the WSV with the same name.

   @author Patrick Eriksson
   @date   2012-06-22
 */
void rte_pos2gridpos(GridPos& gp_p,
                     GridPos& gp_lat,
                     GridPos& gp_lon,
                     const Index& atmosphere_dim,
                     ConstVectorView p_grid,
                     ConstVectorView lat_grid,
                     ConstVectorView lon_grid,
                     ConstTensor3View z_field,
                     ConstVectorView rte_pos);

/** Converts a geographical position (rte_pos) to grid positions for lat and
   lon. (surface version) 

   The function converts latitude and longitude in *rte_pos* to matching grid
   positions. Handles 1D, 2D and 3D (gp_lat and gp_lon untouched if not used).

   Note that the function performs several checks of runtime error type.

   @param[in]   gp_lat      Output: Latitude grid position.
   @param[in]   gp_lon      Output: Longitude grid position.
   @param[in]   atmosphere_dim  As the WSV with the same name.
   @param[in]   lat_grid    As the WSV with the same name.
   @param[in]   lon_grid    As the WSV with the same name.
   @param[in]   rte_pos     As the WSV with the same name.

   @author Patrick Eriksson
   @date   2018-04-01
 */
void rte_pos2gridpos(GridPos& gp_lat,
                     GridPos& gp_lon,
                     const Index& atmosphere_dim,
                     ConstVectorView lat_grid,
                     ConstVectorView lon_grid,
                     ConstVectorView rte_pos);

/** Returns the geomtrical altitudes of *p_grid* for one latitude.

    The latitude is specified by its grid position, in an
    ArrayOfGridPos of length 1. The altitude field (*z_field*) is then
    interpolated to that latitude.

    @param[out]   z         Found altitudes.
    @param[in]   p_grid     As the WSV with the same name.
    @param[in]   lat_grid   As the WSV with the same name.
    @param[in]   z_field    The pressure and latitude part of the WSV with 
                        the same name (that is, the first column).
    @param[in]   gp_lat     Latitude grid position.

    @author Patrick Eriksson 
    @date   2002-11-18
 */
void z_at_lat_2d(VectorView z,
                 ConstVectorView p_grid,
                 ConstVectorView lat_grid,
                 ConstMatrixView z_field,
                 const GridPos& gp_lat);

/** Returns the geomtrical altitudes of *p_grid* for one latitude and
    one longitude.

    The latitude and longitude are specified by their grid position,
    in an ArrayOfGridPos of length 1. The altitude field (*z_field*)
    is then interpolated to that latitude and longitude.

    @param[in]   z          Out: Found altitudes.
    @param[in]   p_grid     As the WSV with the same name.
    @param[in]   lat_grid   As the WSV with the same name.
    @param[in]   lon_grid   As the WSV with the same name.
    @param[in]   z_field    As the WSV with the same name.
    @param[in]   gp_lat     Latitude grid positions.
    @param[in]   gp_lon     Longitude grid positions.

    @author Patrick Eriksson 
    @date   2002-12-31
 */
void z_at_latlon(VectorView z,
                 ConstVectorView p_grid,
                 ConstVectorView lat_grid,
                 ConstVectorView lon_grid,
                 ConstTensor3View z_field,
                 const GridPos& gp_lat,
                 const GridPos& gp_lon);

/** General function for interpolating data of complex n type.

    See documentation of comples_refr_index for format of complex_n-.

    @param[out]  n_real      Real part [nf,nt]
    @param[out]  n_imag      Imaginary part [nf,nt]
    @param[in]   complex_n   Complex refracton index data.
    @param[in]   varname     The name of complex_n to use in error message.
    @param[in]   f_grid      Output frequency grid [nf]
    @param[in]   t_grid      Output temperature grid [nt]

    @author Patrick Eriksson 
    @date   2013-08-16
 */
void complex_n_interp(MatrixView n_real,
                      MatrixView n_imag,
                      const GriddedField3& complex_n,
                      const String& varname,
                      ConstVectorView f_grid,
                      ConstVectorView t_grid);

#endif  // special_interp_h
