/**
  @file   special_interp.cc
  @author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  @date   2002-11-14   
  
  @brief  Interpolation routines for special purposes.
  
  This file contains functions connected to interpolation of
  non-general character. The total general interpolation routines are
  found in interpolation.cc. 

  - Point(s) interpolation of atmospheric/surface fields: 

  These interpolation functions interpolate a atmospheric/surface fields
  to a set of points, such as the points of a propagation path. That is,
  "blue" interpolation. The functions assume that the grid positions
  are at hand. If several atmospheric fields shall be interpolated,
  the functions to use are *interp_atmfield_gp2itw* and
  *interp_atmfield_by_itw*, where the first function is called once
  and the second is called for each field to be interpolated. If only
  one field shall be interpolated, the function
  *interp_atmfield_by_gp* is a shortcut for calling both functions
  above. There exist an identical set of functions for interpolating
  surface-type variables, with names where _atmfield_ is replaced with
  _atmsurface_. 

  Possible surface-type variables are *z_surface* and one page
  of *z_field*.

  - Regridding of atmospheric/surface fields: 

  These functions interpolate from one set of atmospheric grids to a new 
  set of grids.

  - Conversion of geometric altitudes to pressure values: 

  To convert a geometric altitude to a pressure results in an
  interpolation of the pressure grid, or more exactly the log of the
  *p_grid*. Such functions are placed in this file for that reason.
  These functions have names ending with "2p", for example itw2p.

  - Interpolation of Gridded fields of special types:
 */

#include <cmath>

#include "check_input.h"
#include "special_interp.h"

/*===========================================================================
  === Point interpolation functions for atmospheric grids and fields
  ===========================================================================*/

void interp_atmfield_gp2itw(Matrix& itw,
                            const ArrayOfGridPos& gp_p,
                            const ArrayOfGridPos& gp_lat,
                            const ArrayOfGridPos& gp_lon) {
  const Size n = gp_p.size();
  ARTS_ASSERT(gp_lat.size() == n);
  ARTS_ASSERT(gp_lon.size() == n);
  itw.resize(n, 8);
  interpweights(itw, gp_p, gp_lat, gp_lon);
}

void interp_atmfield_by_itw(VectorView x,
                            ConstTensor3View x_field,
                            const ArrayOfGridPos& gp_p,
                            const ArrayOfGridPos& gp_lat,
                            const ArrayOfGridPos& gp_lon,
                            ConstMatrixView itw) {
  ARTS_ASSERT(static_cast<Size>(x.size()) == gp_p.size());
  ARTS_ASSERT(itw.ncols() == 8);
  interp(x, itw, x_field, gp_p, gp_lat, gp_lon);
}

void interp_atmfield_by_gp(VectorView x,
                           ConstTensor3View x_field,
                           const ArrayOfGridPos& gp_p,
                           const ArrayOfGridPos& gp_lat,
                           const ArrayOfGridPos& gp_lon) {
  Matrix itw;

  interp_atmfield_gp2itw(itw, gp_p, gp_lat, gp_lon);

  interp_atmfield_by_itw(x, x_field, gp_p, gp_lat, gp_lon, itw);
}

Numeric interp_atmfield_by_gp(ConstTensor3View x_field,
                              const GridPos& gp_p,
                              const GridPos& gp_lat,
                              const GridPos& gp_lon) {
  ArrayOfGridPos agp_p(1), agp_lat(0), agp_lon(0);

  gridpos_copy(agp_p[0], gp_p);

    agp_lat.resize(1);
    gridpos_copy(agp_lat[0], gp_lat);
  

    agp_lon.resize(1);
    gridpos_copy(agp_lon[0], gp_lon);
  

  Vector x(1);

  interp_atmfield_by_gp(x, x_field, agp_p, agp_lat, agp_lon);

  return x[0];
}

void interp_cloudfield_gp2itw(VectorView itw,
                              GridPos& gp_p_out,
                              GridPos& gp_lat_out,
                              GridPos& gp_lon_out,
                              const GridPos& gp_p_in,
                              const GridPos& gp_lat_in,
                              const GridPos& gp_lon_in,
                              const ArrayOfIndex& cloudbox_limits) {
  // Shift grid positions to cloud box grids
    gridpos_copy(gp_p_out, gp_p_in);
    gridpos_copy(gp_lat_out, gp_lat_in);
    gridpos_copy(gp_lon_out, gp_lon_in);
    gp_p_out.idx -= cloudbox_limits[0];
    gp_lat_out.idx -= cloudbox_limits[2];
    gp_lon_out.idx -= cloudbox_limits[4];
    gridpos_upperend_check(gp_p_out, cloudbox_limits[1] - cloudbox_limits[0]);
    gridpos_upperend_check(gp_lat_out, cloudbox_limits[3] - cloudbox_limits[2]);
    gridpos_upperend_check(gp_lon_out, cloudbox_limits[5] - cloudbox_limits[4]);
    ARTS_ASSERT(itw.size() == 8);
    interpweights(itw, gp_p_out, gp_lat_out, gp_lon_out);
}

/** Converts atmospheric grid positions to weights for interpolation of a
    surface-type variable.

    The function is intended for "blue" interpolation, that is, interpolation
    for a set of positions. 

    The output matrix for interpolation weights are resized inside the
    function.

    The input atmospheric grids are checked to be consistent.

    \param[out]  itw                Interpolation weights.
    \param[in]   gp_lat             Latitude grid positions.
    \param[in]   gp_lon             Longitude grid positions.

    @author Patrick Eriksson 
    @date   2002-11-13
 */
void interp_atmsurface_gp2itw(Matrix& itw,
                              const ArrayOfGridPos& gp_lat,
                              const ArrayOfGridPos& gp_lon) {
    const Size n = gp_lat.size();
    ARTS_ASSERT(n == gp_lon.size());
    itw.resize(n, 4);
    interpweights(itw, gp_lat, gp_lon);
}

void interp_atmsurface_by_itw(VectorView x,
                              ConstMatrixView x_surface,
                              const ArrayOfGridPos& gp_lat,
                              const ArrayOfGridPos& gp_lon,
                              ConstMatrixView itw) {
    ARTS_ASSERT(static_cast<Size>(x.size()) == gp_lat.size());
    ARTS_ASSERT(itw.ncols() == 4);
    interp(x, itw, x_surface, gp_lat, gp_lon);
}

void interp_atmsurface_by_gp(VectorView x,
                             ConstMatrixView x_surface,
                             const ArrayOfGridPos& gp_lat,
                             const ArrayOfGridPos& gp_lon) {
  Matrix itw;

  interp_atmsurface_gp2itw(itw, gp_lat, gp_lon);

  interp_atmsurface_by_itw(x, x_surface, gp_lat, gp_lon, itw);
}

Numeric interp_atmsurface_by_gp(ConstMatrixView x_surface,
                                const GridPos& gp_lat,
                                const GridPos& gp_lon) {
  ArrayOfGridPos agp_lat(0), agp_lon(0);

    agp_lat.resize(1);
    gridpos_copy(agp_lat[0], gp_lat);

    agp_lon.resize(1);
    gridpos_copy(agp_lon[0], gp_lon);

  Vector x(1);

  interp_atmsurface_by_gp(x, x_surface, agp_lat, agp_lon);

  return x[0];
}

/*===========================================================================
  === Regridding
  ===========================================================================*/

void regrid_atmfield_by_gp(Tensor3& field_new,
                           ConstTensor3View field_old,
                           const ArrayOfGridPos& gp_p,
                           const ArrayOfGridPos& gp_lat,
                           const ArrayOfGridPos& gp_lon) {
  const Size n1 = gp_p.size();

    const Size n2 = gp_lat.size();
    const Size n3 = gp_lon.size();
    field_new.resize(n1, n2, n3);
    Tensor4 itw(n1, n2, n3, 8);
    interpweights(itw, gp_p, gp_lat, gp_lon);
    interp(field_new, itw, field_old, gp_p, gp_lat, gp_lon);
}

void regrid_atmsurf_by_gp(Matrix& field_new,
                          ConstMatrixView field_old,
                          const ArrayOfGridPos& gp_lat,
                          const ArrayOfGridPos& gp_lon) {
    const Size n1 = gp_lat.size();
    const Size n2 = gp_lon.size();
    field_new.resize(n1, n2);
    Tensor3 itw(n1, n2, 4);
    interpweights(itw, gp_lat, gp_lon);
    interp(field_new, itw, field_old, gp_lat, gp_lon);
}

void get_gp_atmgrids_to_rq(ArrayOfGridPos& gp_p,
                           ArrayOfGridPos& gp_lat,
                           ArrayOfGridPos& gp_lon,
                           const RetrievalQuantity& rq,
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid) {
  gp_p.resize(rq.Grids()[0].size());
  p2gridpos(gp_p, p_grid, rq.Grids()[0], 0);
  //
    gp_lat.resize(rq.Grids()[1].size());
    gridpos(gp_lat, lat_grid, rq.Grids()[1], 0);
  //
    gp_lon.resize(rq.Grids()[2].size());
    gridpos(gp_lon, lon_grid, rq.Grids()[2], 0);
}

void get_gp_atmsurf_to_rq(ArrayOfGridPos& gp_lat,
                          ArrayOfGridPos& gp_lon,
                          const RetrievalQuantity& rq,
                          const Vector& lat_grid,
                          const Vector& lon_grid) {
    gp_lat.resize(rq.Grids()[0].size());
    gridpos(gp_lat, lat_grid, rq.Grids()[0], 0);
  //
    gp_lon.resize(rq.Grids()[1].size());
    gridpos(gp_lon, lon_grid, rq.Grids()[1], 0);
}

void get_gp_rq_to_atmgrids(ArrayOfGridPos& gp_p,
                           ArrayOfGridPos& gp_lat,
                           ArrayOfGridPos& gp_lon,
                           Index& n_p,
                           Index& n_lat,
                           Index& n_lon,
                           const ArrayOfVector& ret_grids,
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid) {
  // We want here an extrapolation to infinity ->
  //                                        extremly high extrapolation factor
  const Numeric inf_proxy = 1.0e99;

  gp_p.resize(p_grid.size());
  n_p = ret_grids[0].size();
  if (n_p > 1) {
    p2gridpos(gp_p, ret_grids[0], p_grid, inf_proxy);
    jacobian_type_extrapol(gp_p);
  } else {
    gp4length1grid(gp_p);
  }

    gp_lat.resize(lat_grid.size());
    n_lat = ret_grids[1].size();
    if (n_lat > 1) {
      gridpos(gp_lat, ret_grids[1], lat_grid, inf_proxy);
      jacobian_type_extrapol(gp_lat);
    } else {
      gp4length1grid(gp_lat);
    }
  //
    gp_lon.resize(lon_grid.size());
    n_lon = ret_grids[2].size();
    if (n_lon > 1) {
      gridpos(gp_lon, ret_grids[2], lon_grid, inf_proxy);
      jacobian_type_extrapol(gp_lon);
    } else {
      gp4length1grid(gp_lon);
    }
}

void get_gp_rq_to_atmgrids(ArrayOfGridPos& gp_lat,
                           ArrayOfGridPos& gp_lon,
                           Index& n_lat,
                           Index& n_lon,
                           const ArrayOfVector& ret_grids,
                           const Vector& lat_grid,
                           const Vector& lon_grid) {
  // We want here an extrapolation to infinity ->
  //                                        extremly high extrapolation factor
  const Numeric inf_proxy = 1.0e99;

    gp_lat.resize(lat_grid.size());
    n_lat = ret_grids[0].size();
    if (n_lat > 1) {
      gridpos(gp_lat, ret_grids[0], lat_grid, inf_proxy);
      jacobian_type_extrapol(gp_lat);
    } else {
      gp4length1grid(gp_lat);
    }
  //
    gp_lon.resize(lon_grid.size());
    n_lon = ret_grids[1].size();
    if (n_lon > 1) {
      gridpos(gp_lon, ret_grids[1], lon_grid, inf_proxy);
      jacobian_type_extrapol(gp_lon);
    } else {
      gp4length1grid(gp_lon);
    }
}

void regrid_atmfield_by_gp_oem(Tensor3& field_new,
                               ConstTensor3View field_old,
                               const ArrayOfGridPos& gp_p,
                               const ArrayOfGridPos& gp_lat,
                               const ArrayOfGridPos& gp_lon) {
  const Size n1 = gp_p.size();

  const bool np_is1 = field_old.npages() == 1 ? true : false;
  const bool nlat_is1 =
      field_old.nrows() == 1 ? true : false;
  const bool nlon_is1 =
      field_old.ncols() == 1 ? true : false;

  // If no length 1, we can use standard function
  if (!np_is1 && !nlat_is1 && !nlon_is1) {
    regrid_atmfield_by_gp(
        field_new, field_old, gp_p, gp_lat, gp_lon);
  } else {
      const Size n2 = gp_lat.size();
      const Size n3 = gp_lon.size();
      field_new.resize(n1, n2, n3);
      //
      if (np_is1 && nlat_is1 && nlon_is1)  // 1: No interpolation at all
      {
        field_new(joker, joker, joker) = field_old(0, 0, 0);
      }

      else if (np_is1)  // No pressure interpolation --------------
      {
        if (nlat_is1)  // 2: Just longitude interpolation
        {
          Matrix itw(n3, 2);
          interpweights(itw, gp_lon);
          Vector tmp(n3);
          interp(tmp, itw, field_old(0, 0, joker), gp_lon);
          for (Size p = 0; p < n1; p++) {
            ARTS_ASSERT(gp_p[p].fd[0] < 1e-6);
            for (Size lat = 0; lat < n2; lat++) {
              ARTS_ASSERT(gp_lat[lat].fd[0] < 1e-6);
              field_new(p, lat, joker) = tmp;
            }
          }
        } else if (nlon_is1)  // 3: Just latitude interpolation
        {
          Matrix itw(n2, 2);
          interpweights(itw, gp_lat);
          Vector tmp(n2);
          interp(tmp, itw, field_old(0, joker, 0), gp_lat);
          for (Size p = 0; p < n1; p++) {
            ARTS_ASSERT(gp_p[p].fd[0] < 1e-6);
            for (Size lon = 0; lon < n3; lon++) {
              ARTS_ASSERT(gp_lon[lon].fd[0] < 1e-6);
              field_new(p, joker, lon) = tmp;
            }
          }
        } else  // 4: Both lat and lon interpolation
        {
          Tensor3 itw(n2, n3, 4);
          interpweights(itw, gp_lat, gp_lon);
          Matrix tmp(n2, n3);
          interp(tmp, itw, field_old(0, joker, joker), gp_lat, gp_lon);
          for (Size p = 0; p < n1; p++) {
            ARTS_ASSERT(gp_p[p].fd[0] < 1e-6);
            field_new(p, joker, joker) = tmp;
          }
        }
      }

      else  // Pressure interpolation --------------
      {
        if (nlat_is1 && nlon_is1)  // 5: Just pressure interpolatiom
        {
          Matrix itw(n1, 2);
          interpweights(itw, gp_p);
          Vector tmp(n1);
          interp(tmp, itw, field_old(joker, 0, 0), gp_p);
          for (Size lat = 0; lat < n2; lat++) {
            ARTS_ASSERT(gp_lat[lat].fd[0] < 1e-6);
            for (Size lon = 0; lon < n3; lon++) {
              ARTS_ASSERT(gp_lon[lon].fd[0] < 1e-6);
              field_new(joker, lat, lon) = tmp;
            }
          }
        } else if (nlat_is1)  // 6: Both p and lon interpolation
        {
          Tensor3 itw(n1, n3, 4);
          interpweights(itw, gp_p, gp_lon);
          Matrix tmp(n1, n3);
          interp(tmp, itw, field_old(joker, 0, joker), gp_p, gp_lon);
          for (Size lat = 0; lat < n2; lat++) {
            ARTS_ASSERT(gp_lat[lat].fd[0] < 1e-6);
            field_new(joker, lat, joker) = tmp;
          }
        } else  // 7: Both p and lat interpolation
        {
          Tensor3 itw(n1, n2, 4);
          interpweights(itw, gp_p, gp_lat);
          Matrix tmp(n1, n2);
          interp(tmp, itw, field_old(joker, joker, 0), gp_p, gp_lat);
          for (Size lon = 0; lon < n3; lon++) {
            ARTS_ASSERT(gp_lon[lon].fd[0] < 1e-6);
            field_new(joker, joker, lon) = tmp;
          }
        }
    }
  }
}

/* So far just a temporary test */
void regrid_atmsurf_by_gp_oem(Matrix& field_new,
                              ConstMatrixView field_old,
                              const ArrayOfGridPos& gp_lat,
                              const ArrayOfGridPos& gp_lon) {
  // As 1D is so simple, let's do it here and not go to standard function
    const bool nlat_is1 = field_old.nrows() == 1 ? true : false;
    const bool nlon_is1 =
       field_old.ncols() == 1 ? true : false;

    // If no length 1, we can use standard function
    if (!nlat_is1 && !nlon_is1) {
      regrid_atmsurf_by_gp(
          field_new, field_old, gp_lat, gp_lon);
    } else {
        const Index n1 = gp_lat.size();
        const Index n2 = gp_lon.size();
        field_new.resize(n1, n2);
        //
        if (nlat_is1 && nlon_is1)  // 1: No interpolation at all
        {
          field_new(joker, joker) = field_old(0, 0);
        } else if (nlon_is1)  // 2: Just latitude interpolation
        {
          Matrix itw(n1, 2);
          interpweights(itw, gp_lat);
          Vector tmp(n1);
          interp(tmp, itw, field_old(joker, 0), gp_lat);
          for (Index lon = 0; lon < n2; lon++) {
            ARTS_ASSERT(gp_lon[lon].fd[0] < 1e-6);
            field_new(joker, lon) = tmp;
          }
        } else  // 2: Just longitude interpolation
        {
          Matrix itw(n2, 2);
          interpweights(itw, gp_lon);
          Vector tmp(n2);
          interp(tmp, itw, field_old(0, joker), gp_lon);
          for (Index lat = 0; lat < n1; lat++) {
            ARTS_ASSERT(gp_lat[lat].fd[0] < 1e-6);
            field_new(lat, joker) = tmp;
          }
        }
      }
    }

/*===========================================================================
  === Conversion altitudes / pressure
  ===========================================================================*/

void itw2p(VectorView p_values,
           ConstVectorView p_grid,
           const ArrayOfGridPos& gp,
           ConstMatrixView itw) {
  ARTS_ASSERT(itw.ncols() == 2);
  ARTS_ASSERT(p_values.size() == itw.nrows());

  // Local variable to store log of the pressure grid:
  Vector logpgrid(p_grid.size());

  transform(logpgrid, log, p_grid);

  interp(p_values, itw, logpgrid, gp);

  transform(p_values, exp, p_values);
}

/** Calculates grid positions for pressure values.

   This function works as *gridpos*, but is adapted to handle
   pressure grids. The ARTS defintions result in that pressures shall
   not be interpolated directly, it is the log of the pressure that
   shall be interpolated. This means that if some values shall be
   interpolated to some given pressures, the grid positions shall be
   calculated with this function. The interpolation can then be
   performed as usual.

   \param[out]  gp          Grid position Array.
   \param[in]   old_pgrid   The original pressure grid.
   \param[in]   new_pgrid   The new pressure grid.
   \param[in]   extpolfac   Extrapolation factor. Default value is 0.5,
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
               const Numeric& extpolfac) {
  // Local variable to store log of the pressure grids
  Vector logold(old_pgrid.size());
  Vector lognew(new_pgrid.size());

  transform(logold, log, old_pgrid);
  transform(lognew, log, new_pgrid);

  gridpos(gp, logold, lognew, extpolfac);
}

void rte_pos2gridpos(GridPos& gp_p,
                     GridPos& gp_lat,
                     GridPos& gp_lon,
                     ConstVectorView p_grid,
                     ConstVectorView lat_grid,
                     ConstVectorView lon_grid,
                     ConstTensor3View z_field,
                     ConstVectorView rte_pos) {
  chk_rte_pos(rte_pos);

    // Determine z at lat/lon (z_grid) by blue interpolation
    const Index np = p_grid.size();
    Vector z_grid(np);
    ArrayOfGridPos agp_z, agp_lat(np);
    //
    gridpos_1to1(agp_z, p_grid);
    //
    chk_interpolation_grids("Latitude interpolation", lat_grid, rte_pos[1]);
    gridpos(gp_lat, lat_grid, rte_pos[1]);

      chk_interpolation_grids("Longitude interpolation", lon_grid, rte_pos[2]);
      gridpos(gp_lon, lon_grid, rte_pos[2]);
      ArrayOfGridPos agp_lon(np);
      for (Index i = 0; i < np; i++) {
        agp_lat[i] = gp_lat;
        agp_lon[i] = gp_lon;
      }
      Matrix itw(np, 8);
      interpweights(itw, agp_z, agp_lat, agp_lon);
      interp(z_grid, itw, z_field, agp_z, agp_lat, agp_lon);

    // And use z_grid to get gp_p (gp_al and gp_lon determined above)
    chk_interpolation_grids("Altitude interpolation", z_grid, rte_pos[0]);
    gridpos(gp_p, z_grid, rte_pos[0]);
}

void rte_pos2gridpos(GridPos& gp_lat,
                     GridPos& gp_lon,
                     ConstVectorView lat_grid,
                     ConstVectorView lon_grid,
                     ConstVectorView rte_pos) {
  chk_rte_pos(rte_pos);

    chk_interpolation_grids("Latitude interpolation", lat_grid, rte_pos[1]);
    gridpos(gp_lat, lat_grid, rte_pos[1]);

      chk_interpolation_grids("Longitude interpolation", lon_grid, rte_pos[2]);
      gridpos(gp_lon, lon_grid, rte_pos[2]);
}

void z_at_lat_2d(VectorView z,
                 ConstVectorView p_grid,
// FIXME only used in assertion
#ifndef NDEBUG
                 ConstVectorView lat_grid,
#else
                 ConstVectorView,
#endif
                 ConstMatrixView z_field,
                 const GridPos& gp_lat) {
  const Index np = p_grid.size();

  ARTS_ASSERT(z.size() == np);
  ARTS_ASSERT(z_field.nrows() == np);
  ARTS_ASSERT(z_field.ncols() == lat_grid.size());

  Matrix z_matrix(np, 1);
  ArrayOfGridPos gp_z(np), agp_lat(1);
  Tensor3 itw(np, 1, 4);

  gridpos_copy(agp_lat[0], gp_lat);
  gridpos(gp_z, p_grid, p_grid);
  interpweights(itw, gp_z, agp_lat);
  interp(z_matrix, itw, z_field, gp_z, agp_lat);

  z = z_matrix(Range(joker), 0);
}

void z_at_latlon(VectorView z,
                 ConstVectorView p_grid,
//FIXME only used in assertion
#ifndef NDEBUG
                 ConstVectorView lat_grid,
                 ConstVectorView lon_grid,
#else
                 ConstVectorView,
                 ConstVectorView,
#endif
                 ConstTensor3View z_field,
                 const GridPos& gp_lat,
                 const GridPos& gp_lon) {
  const Index np = p_grid.size();

  ARTS_ASSERT(z.size() == np);
  ARTS_ASSERT(z_field.npages() == np);
  ARTS_ASSERT(z_field.nrows() == lat_grid.size());
  ARTS_ASSERT(z_field.ncols() == lon_grid.size());

  Tensor3 z_tensor(np, 1, 1);
  ArrayOfGridPos agp_z(np), agp_lat(1), agp_lon(1);
  Tensor4 itw(np, 1, 1, 8);

  gridpos_copy(agp_lat[0], gp_lat);
  gridpos_copy(agp_lon[0], gp_lon);
  gridpos(agp_z, p_grid, p_grid);
  interpweights(itw, agp_z, agp_lat, agp_lon);

  interp(z_tensor, itw, z_field, agp_z, agp_lat, agp_lon);

  z = z_tensor(Range(joker), 0, 0);
}

/*===========================================================================
  === Interpolation of GriddedFields
  ===========================================================================*/

void complex_n_interp(MatrixView n_real,
                      MatrixView n_imag,
                      const ComplexGriddedField2& complex_n,
                      const String&,
                      ConstVectorView f_grid,
                      ConstVectorView t_grid) {
  // Set expected order of grids
  // Check of complex_n
  //
  ARTS_USER_ERROR_IF(not complex_n.check(), "Bad complex_n");
  //
  ARTS_USER_ERROR_IF(complex_n.gridname<0>() != "Frequency",
                     "Bad gridname for frequency");
  ARTS_USER_ERROR_IF(complex_n.gridname<1>() != "Temperature",
                     "Bad gridname for temperature");

  // Frequency and temperature grid sizes
  const Index nf_in = complex_n.data.nrows();
  const Index nt_in = complex_n.data.ncols();
  const Index nf_out = f_grid.size();
  const Index nt_out = t_grid.size();

  //Assert size of output variables
  ARTS_ASSERT(n_real.nrows() == nf_out && n_real.ncols() == nt_out);
  ARTS_ASSERT(n_imag.nrows() == nf_out && n_imag.ncols() == nt_out);

  const Vector& f_grid_in = complex_n.grid<0>();
  const Vector& t_grid_in = complex_n.grid<1>();

  // Expand/interpolate in frequency dimension
  //
  Matrix nrf(nf_out, nt_in), nif(nf_out, nt_in);
  //
  if (nf_in == 1) {
    for (Index i = 0; i < nf_out; i++) {
      nrf(i, joker) = complex_n.data(0, joker).real();
      nif(i, joker) = complex_n.data(0, joker).imag();
    }
  } else {
    chk_interpolation_grids("Frequency interpolation", f_grid_in, f_grid);
    //
    ArrayOfGridPos gp(nf_out);
    Matrix itw(nf_out, 2);
    gridpos(gp, f_grid_in, f_grid);
    interpweights(itw, gp);
    for (Index i = 0; i < nt_in; i++) {
      interp(nrf(joker, i), itw, complex_n.data(joker, i).real(), gp);
      interp(nif(joker, i), itw, complex_n.data(joker, i).imag(), gp);
    }
  }

  // Expand/interpolate in temperature dimension
  //
  if (nt_in == 1) {
    for (Index i = 0; i < nt_out; i++) {
      n_real(joker, i) = nrf(joker, 0);
      n_imag(joker, i) = nif(joker, 0);
    }
  } else {
    chk_interpolation_grids("Temperature interpolation", t_grid_in, t_grid);
    //
    ArrayOfGridPos gp(nt_out);
    Matrix itw(nt_out, 2);
    gridpos(gp, t_grid_in, t_grid);
    interpweights(itw, gp);
    for (Index i = 0; i < nf_out; i++) {
      interp(n_real(i, joker), itw, nrf(i, joker), gp);
      interp(n_imag(i, joker), itw, nif(i, joker), gp);
    }
  }
}
