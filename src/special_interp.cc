/* Copyright (C) 2002-2012 Stefan Buehler <sbuehler@ltu.se>

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
#include <iostream>
#include <stdexcept>
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "special_interp.h"

/*===========================================================================
  === Point interpolation functions for atmospheric grids and fields
  ===========================================================================*/

void interp_atmfield_gp2itw(Matrix& itw,
                            const Index& atmosphere_dim,
                            const ArrayOfGridPos& gp_p,
                            const ArrayOfGridPos& gp_lat,
                            const ArrayOfGridPos& gp_lon) {
  const Index n = gp_p.nelem();

  if (atmosphere_dim == 1) {
    itw.resize(n, 2);
    interpweights(itw, gp_p);
  }

  else if (atmosphere_dim == 2) {
    assert(gp_lat.nelem() == n);
    itw.resize(n, 4);
    interpweights(itw, gp_p, gp_lat);
  }

  else if (atmosphere_dim == 3) {
    assert(gp_lat.nelem() == n);
    assert(gp_lon.nelem() == n);
    itw.resize(n, 8);
    interpweights(itw, gp_p, gp_lat, gp_lon);
  }
}

void interp_atmfield_by_itw(VectorView x,
                            const Index& atmosphere_dim,
                            ConstTensor3View x_field,
                            const ArrayOfGridPos& gp_p,
                            const ArrayOfGridPos& gp_lat,
                            const ArrayOfGridPos& gp_lon,
                            ConstMatrixView itw) {
  assert(x.nelem() == gp_p.nelem());

  if (atmosphere_dim == 1) {
    assert(itw.ncols() == 2);
    interp(x, itw, x_field(Range(joker), 0, 0), gp_p);
  }

  else if (atmosphere_dim == 2) {
    assert(itw.ncols() == 4);
    interp(x, itw, x_field(Range(joker), Range(joker), 0), gp_p, gp_lat);
  }

  else if (atmosphere_dim == 3) {
    assert(itw.ncols() == 8);
    interp(x, itw, x_field, gp_p, gp_lat, gp_lon);
  }
}

void interp_atmfield_by_gp(VectorView x,
                           const Index& atmosphere_dim,
                           ConstTensor3View x_field,
                           const ArrayOfGridPos& gp_p,
                           const ArrayOfGridPos& gp_lat,
                           const ArrayOfGridPos& gp_lon) {
  Matrix itw;

  interp_atmfield_gp2itw(itw, atmosphere_dim, gp_p, gp_lat, gp_lon);

  interp_atmfield_by_itw(x, atmosphere_dim, x_field, gp_p, gp_lat, gp_lon, itw);
}

Numeric interp_atmfield_by_gp(const Index& atmosphere_dim,
                              ConstTensor3View x_field,
                              const GridPos& gp_p,
                              const GridPos& gp_lat,
                              const GridPos& gp_lon) {
  ArrayOfGridPos agp_p(1), agp_lat(0), agp_lon(0);

  gridpos_copy(agp_p[0], gp_p);

  if (atmosphere_dim > 1) {
    agp_lat.resize(1);
    gridpos_copy(agp_lat[0], gp_lat);
  }

  if (atmosphere_dim > 2) {
    agp_lon.resize(1);
    gridpos_copy(agp_lon[0], gp_lon);
  }

  Vector x(1);

  interp_atmfield_by_gp(x, atmosphere_dim, x_field, agp_p, agp_lat, agp_lon);

  return x[0];
}

void interp_cloudfield_gp2itw(VectorView itw,
                              GridPos& gp_p_out,
                              GridPos& gp_lat_out,
                              GridPos& gp_lon_out,
                              const GridPos& gp_p_in,
                              const GridPos& gp_lat_in,
                              const GridPos& gp_lon_in,
                              const Index& atmosphere_dim,
                              const ArrayOfIndex& cloudbox_limits) {
  // Shift grid positions to cloud box grids
  if (atmosphere_dim == 1) {
    gridpos_copy(gp_p_out, gp_p_in);
    gp_p_out.idx -= cloudbox_limits[0];
    gridpos_upperend_check(gp_p_out, cloudbox_limits[1] - cloudbox_limits[0]);
    assert(itw.nelem() == 2);
    interpweights(itw, gp_p_out);
  } else if (atmosphere_dim == 2) {
    gridpos_copy(gp_p_out, gp_p_in);
    gridpos_copy(gp_lat_out, gp_lat_in);
    gp_p_out.idx -= cloudbox_limits[0];
    gp_lat_out.idx -= cloudbox_limits[2];
    gridpos_upperend_check(gp_p_out, cloudbox_limits[1] - cloudbox_limits[0]);
    gridpos_upperend_check(gp_lat_out, cloudbox_limits[3] - cloudbox_limits[2]);
    assert(itw.nelem() == 4);
    interpweights(itw, gp_p_out, gp_lat_out);
  } else {
    gridpos_copy(gp_p_out, gp_p_in);
    gridpos_copy(gp_lat_out, gp_lat_in);
    gridpos_copy(gp_lon_out, gp_lon_in);
    gp_p_out.idx -= cloudbox_limits[0];
    gp_lat_out.idx -= cloudbox_limits[2];
    gp_lon_out.idx -= cloudbox_limits[4];
    gridpos_upperend_check(gp_p_out, cloudbox_limits[1] - cloudbox_limits[0]);
    gridpos_upperend_check(gp_lat_out, cloudbox_limits[3] - cloudbox_limits[2]);
    gridpos_upperend_check(gp_lon_out, cloudbox_limits[5] - cloudbox_limits[4]);
    assert(itw.nelem() == 8);
    interpweights(itw, gp_p_out, gp_lat_out, gp_lon_out);
  }
}

/** Converts atmospheric grid positions to weights for interpolation of a
    surface-type variable.

    The function is intended for "blue" interpolation, that is, interpolation
    for a set of positions. 

    The output matrix for interpolation weights are resized inside the
    function.

    The input atmospheric grids are checked to be consistent.

    \param[out]  itw                Interpolation weights.
    \param[in]   atmosphere_dim     As the WSV with the same name.
    \param[in]   gp_lat             Latitude grid positions.
    \param[in]   gp_lon             Longitude grid positions.

    @author Patrick Eriksson 
    @date   2002-11-13
 */
void interp_atmsurface_gp2itw(Matrix& itw,
                              const Index& atmosphere_dim,
                              const ArrayOfGridPos& gp_lat,
                              const ArrayOfGridPos& gp_lon) {
  if (atmosphere_dim == 1) {
    itw.resize(1, 1);
    itw = 1;
  }

  else if (atmosphere_dim == 2) {
    const Index n = gp_lat.nelem();
    itw.resize(n, 2);
    interpweights(itw, gp_lat);
  }

  else if (atmosphere_dim == 3) {
    const Index n = gp_lat.nelem();
    assert(n == gp_lon.nelem());
    itw.resize(n, 4);
    interpweights(itw, gp_lat, gp_lon);
  }
}

void interp_atmsurface_by_itw(VectorView x,
                              const Index& atmosphere_dim,
                              ConstMatrixView x_surface,
                              const ArrayOfGridPos& gp_lat,
                              const ArrayOfGridPos& gp_lon,
                              ConstMatrixView itw) {
  if (atmosphere_dim == 1) {
    assert(itw.ncols() == 1);
    x = x_surface(0, 0);
  }

  else if (atmosphere_dim == 2) {
    assert(x.nelem() == gp_lat.nelem());
    assert(itw.ncols() == 2);
    interp(x, itw, x_surface(Range(joker), 0), gp_lat);
  }

  else if (atmosphere_dim == 3) {
    assert(x.nelem() == gp_lat.nelem());
    assert(itw.ncols() == 4);
    interp(x, itw, x_surface, gp_lat, gp_lon);
  }
}

void interp_atmsurface_by_gp(VectorView x,
                             const Index& atmosphere_dim,
                             ConstMatrixView x_surface,
                             const ArrayOfGridPos& gp_lat,
                             const ArrayOfGridPos& gp_lon) {
  Matrix itw;

  interp_atmsurface_gp2itw(itw, atmosphere_dim, gp_lat, gp_lon);

  interp_atmsurface_by_itw(x, atmosphere_dim, x_surface, gp_lat, gp_lon, itw);
}

Numeric interp_atmsurface_by_gp(const Index& atmosphere_dim,
                                ConstMatrixView x_surface,
                                const GridPos& gp_lat,
                                const GridPos& gp_lon) {
  ArrayOfGridPos agp_lat(0), agp_lon(0);

  if (atmosphere_dim > 1) {
    agp_lat.resize(1);
    gridpos_copy(agp_lat[0], gp_lat);
  }

  if (atmosphere_dim > 2) {
    agp_lon.resize(1);
    gridpos_copy(agp_lon[0], gp_lon);
  }

  Vector x(1);

  interp_atmsurface_by_gp(x, atmosphere_dim, x_surface, agp_lat, agp_lon);

  return x[0];
}

/*===========================================================================
  === Regridding
  ===========================================================================*/

void regrid_atmfield_by_gp(Tensor3& field_new,
                           const Index& atmosphere_dim,
                           ConstTensor3View field_old,
                           const ArrayOfGridPos& gp_p,
                           const ArrayOfGridPos& gp_lat,
                           const ArrayOfGridPos& gp_lon) {
  const Index n1 = gp_p.nelem();

  if (atmosphere_dim == 1) {
    field_new.resize(n1, 1, 1);
    Matrix itw(n1, 2);
    interpweights(itw, gp_p);
    interp(field_new(joker, 0, 0), itw, field_old(joker, 0, 0), gp_p);
  } else if (atmosphere_dim == 2) {
    const Index n2 = gp_lat.nelem();
    field_new.resize(n1, n2, 1);
    Tensor3 itw(n1, n2, 4);
    interpweights(itw, gp_p, gp_lat);
    interp(field_new(joker, joker, 0),
           itw,
           field_old(joker, joker, 0),
           gp_p,
           gp_lat);
  } else if (atmosphere_dim == 3) {
    const Index n2 = gp_lat.nelem();
    const Index n3 = gp_lon.nelem();
    field_new.resize(n1, n2, n3);
    Tensor4 itw(n1, n2, n3, 8);
    interpweights(itw, gp_p, gp_lat, gp_lon);
    interp(field_new, itw, field_old, gp_p, gp_lat, gp_lon);
  }
}

void regrid_atmsurf_by_gp(Matrix& field_new,
                          const Index& atmosphere_dim,
                          ConstMatrixView field_old,
                          const ArrayOfGridPos& gp_lat,
                          const ArrayOfGridPos& gp_lon) {
  if (atmosphere_dim == 1) {
    field_new = field_old;
  } else if (atmosphere_dim == 2) {
    const Index n1 = gp_lat.nelem();
    field_new.resize(n1, 1);
    Matrix itw(n1, 2);
    interpweights(itw, gp_lat);
    interp(field_new(joker, 0), itw, field_old(joker, 0), gp_lat);
  } else if (atmosphere_dim == 3) {
    const Index n1 = gp_lat.nelem();
    const Index n2 = gp_lon.nelem();
    field_new.resize(n1, n2);
    Tensor3 itw(n1, n2, 4);
    interpweights(itw, gp_lat, gp_lon);
    interp(field_new, itw, field_old, gp_lat, gp_lon);
  }
}

void get_gp_atmgrids_to_rq(ArrayOfGridPos& gp_p,
                           ArrayOfGridPos& gp_lat,
                           ArrayOfGridPos& gp_lon,
                           const RetrievalQuantity& rq,
                           const Index& atmosphere_dim,
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid) {
  gp_p.resize(rq.Grids()[0].nelem());
  p2gridpos(gp_p, p_grid, rq.Grids()[0], 0);
  //
  if (atmosphere_dim >= 2) {
    gp_lat.resize(rq.Grids()[1].nelem());
    gridpos(gp_lat, lat_grid, rq.Grids()[1], 0);
  } else {
    gp_lat.resize(0);
  }
  //
  if (atmosphere_dim >= 3) {
    gp_lon.resize(rq.Grids()[2].nelem());
    gridpos(gp_lon, lon_grid, rq.Grids()[2], 0);
  } else {
    gp_lon.resize(0);
  }
}

void get_gp_atmsurf_to_rq(ArrayOfGridPos& gp_lat,
                          ArrayOfGridPos& gp_lon,
                          const RetrievalQuantity& rq,
                          const Index& atmosphere_dim,
                          const Vector& lat_grid,
                          const Vector& lon_grid) {
  if (atmosphere_dim >= 2) {
    gp_lat.resize(rq.Grids()[0].nelem());
    gridpos(gp_lat, lat_grid, rq.Grids()[0], 0);
  } else {
    gp_lat.resize(0);
  }
  //
  if (atmosphere_dim >= 3) {
    gp_lon.resize(rq.Grids()[1].nelem());
    gridpos(gp_lon, lon_grid, rq.Grids()[1], 0);
  } else {
    gp_lon.resize(0);
  }
}

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
                           const Vector& lon_grid) {
  // We want here an extrapolation to infinity ->
  //                                        extremly high extrapolation factor
  const Numeric inf_proxy = 1.0e99;

  gp_p.resize(p_grid.nelem());
  n_p = ret_grids[0].nelem();
  if (n_p > 1) {
    p2gridpos(gp_p, ret_grids[0], p_grid, inf_proxy);
    jacobian_type_extrapol(gp_p);
  } else {
    gp4length1grid(gp_p);
  }

  if (atmosphere_dim >= 2) {
    gp_lat.resize(lat_grid.nelem());
    n_lat = ret_grids[1].nelem();
    if (n_lat > 1) {
      gridpos(gp_lat, ret_grids[1], lat_grid, inf_proxy);
      jacobian_type_extrapol(gp_lat);
    } else {
      gp4length1grid(gp_lat);
    }
  } else {
    gp_lat.resize(0);
    n_lat = 1;
  }
  //
  if (atmosphere_dim >= 3) {
    gp_lon.resize(lon_grid.nelem());
    n_lon = ret_grids[2].nelem();
    if (n_lon > 1) {
      gridpos(gp_lon, ret_grids[2], lon_grid, inf_proxy);
      jacobian_type_extrapol(gp_lon);
    } else {
      gp4length1grid(gp_lon);
    }
  } else {
    gp_lon.resize(0);
    n_lon = 1;
  }
}

void get_gp_rq_to_atmgrids(ArrayOfGridPos& gp_lat,
                           ArrayOfGridPos& gp_lon,
                           Index& n_lat,
                           Index& n_lon,
                           const ArrayOfVector& ret_grids,
                           const Index& atmosphere_dim,
                           const Vector& lat_grid,
                           const Vector& lon_grid) {
  // We want here an extrapolation to infinity ->
  //                                        extremly high extrapolation factor
  const Numeric inf_proxy = 1.0e99;

  if (atmosphere_dim >= 2) {
    gp_lat.resize(lat_grid.nelem());
    n_lat = ret_grids[0].nelem();
    if (n_lat > 1) {
      gridpos(gp_lat, ret_grids[0], lat_grid, inf_proxy);
      jacobian_type_extrapol(gp_lat);
    } else {
      gp4length1grid(gp_lat);
    }
  } else {
    gp_lat.resize(0);
    n_lat = 1;
  }
  //
  if (atmosphere_dim >= 3) {
    gp_lon.resize(lon_grid.nelem());
    n_lon = ret_grids[1].nelem();
    if (n_lon > 1) {
      gridpos(gp_lon, ret_grids[1], lon_grid, inf_proxy);
      jacobian_type_extrapol(gp_lon);
    } else {
      gp4length1grid(gp_lon);
    }
  } else {
    gp_lon.resize(0);
    n_lon = 1;
  }
}

void regrid_atmfield_by_gp_oem(Tensor3& field_new,
                               const Index& atmosphere_dim,
                               ConstTensor3View field_old,
                               const ArrayOfGridPos& gp_p,
                               const ArrayOfGridPos& gp_lat,
                               const ArrayOfGridPos& gp_lon) {
  const Index n1 = gp_p.nelem();

  const bool np_is1 = field_old.npages() == 1 ? true : false;
  const bool nlat_is1 =
      atmosphere_dim > 1 && field_old.nrows() == 1 ? true : false;
  const bool nlon_is1 =
      atmosphere_dim > 2 && field_old.ncols() == 1 ? true : false;

  // If no length 1, we can use standard function
  if (!np_is1 && !nlat_is1 && !nlon_is1) {
    regrid_atmfield_by_gp(
        field_new, atmosphere_dim, field_old, gp_p, gp_lat, gp_lon);
  } else {
    //--- 1D (1 possibilities left) -------------------------------------------
    if (atmosphere_dim == 1) {  // 1: No interpolation at all
      field_new.resize(n1, 1, 1);
      field_new(joker, 0, 0) = field_old(0, 0, 0);
    }

    //--- 2D (3 possibilities left) -------------------------------------------
    else if (atmosphere_dim == 2) {
      const Index n2 = gp_lat.nelem();
      field_new.resize(n1, n2, 1);
      //
      if (np_is1 && nlat_is1)  // 1: No interpolation at all
      {
        // Here we need no interpolation at all
        field_new(joker, joker, 0) = field_old(0, 0, 0);
      } else if (np_is1)  // 2: Latitude interpolation
      {
        Matrix itw(n2, 2);
        interpweights(itw, gp_lat);
        Vector tmp(n2);
        interp(tmp, itw, field_old(0, joker, 0), gp_lat);
        for (Index p = 0; p < n1; p++) {
          assert(gp_p[p].fd[0] < 1e-6);
          field_new(p, joker, 0) = tmp;
        }
      } else  // 3: Pressure interpolation
      {
        Matrix itw(n1, 2);
        interpweights(itw, gp_p);
        Vector tmp(n1);
        interp(tmp, itw, field_old(joker, 0, 0), gp_p);
        for (Index lat = 0; lat < n2; lat++) {
          assert(gp_lat[lat].fd[0] < 1e-6);
          field_new(joker, lat, 0) = tmp;
        }
      }
    }

    //--- 3D (7 possibilities left) -------------------------------------------
    else if (atmosphere_dim == 3) {
      const Index n2 = gp_lat.nelem();
      const Index n3 = gp_lon.nelem();
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
          for (Index p = 0; p < n1; p++) {
            assert(gp_p[p].fd[0] < 1e-6);
            for (Index lat = 0; lat < n2; lat++) {
              assert(gp_lat[lat].fd[0] < 1e-6);
              field_new(p, lat, joker) = tmp;
            }
          }
        } else if (nlon_is1)  // 3: Just latitude interpolation
        {
          Matrix itw(n2, 2);
          interpweights(itw, gp_lat);
          Vector tmp(n2);
          interp(tmp, itw, field_old(0, joker, 0), gp_lat);
          for (Index p = 0; p < n1; p++) {
            assert(gp_p[p].fd[0] < 1e-6);
            for (Index lon = 0; lon < n3; lon++) {
              assert(gp_lon[lon].fd[0] < 1e-6);
              field_new(p, joker, lon) = tmp;
            }
          }
        } else  // 4: Both lat and lon interpolation
        {
          Tensor3 itw(n2, n3, 4);
          interpweights(itw, gp_lat, gp_lon);
          Matrix tmp(n2, n3);
          interp(tmp, itw, field_old(0, joker, joker), gp_lat, gp_lon);
          for (Index p = 0; p < n1; p++) {
            assert(gp_p[p].fd[0] < 1e-6);
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
          for (Index lat = 0; lat < n2; lat++) {
            assert(gp_lat[lat].fd[0] < 1e-6);
            for (Index lon = 0; lon < n3; lon++) {
              assert(gp_lon[lon].fd[0] < 1e-6);
              field_new(joker, lat, lon) = tmp;
            }
          }
        } else if (nlat_is1)  // 6: Both p and lon interpolation
        {
          Tensor3 itw(n1, n3, 4);
          interpweights(itw, gp_p, gp_lon);
          Matrix tmp(n1, n3);
          interp(tmp, itw, field_old(joker, 0, joker), gp_p, gp_lon);
          for (Index lat = 0; lat < n2; lat++) {
            assert(gp_lat[lat].fd[0] < 1e-6);
            field_new(joker, lat, joker) = tmp;
          }
        } else  // 7: Both p and lat interpolation
        {
          Tensor3 itw(n1, n2, 4);
          interpweights(itw, gp_p, gp_lat);
          Matrix tmp(n1, n2);
          interp(tmp, itw, field_old(joker, joker, 0), gp_p, gp_lat);
          for (Index lon = 0; lon < n3; lon++) {
            assert(gp_lon[lon].fd[0] < 1e-6);
            field_new(joker, joker, lon) = tmp;
          }
        }
      }
    }
  }
}

/* So far just a temporary test */
void regrid_atmsurf_by_gp_oem(Matrix& field_new,
                              const Index& atmosphere_dim,
                              ConstMatrixView field_old,
                              const ArrayOfGridPos& gp_lat,
                              const ArrayOfGridPos& gp_lon) {
  // As 1D is so simple, let's do it here and not go to standard function
  if (atmosphere_dim == 1) {
    field_new = field_old;
  } else {
    const bool nlat_is1 = field_old.nrows() == 1 ? true : false;
    const bool nlon_is1 =
        atmosphere_dim > 2 && field_old.ncols() == 1 ? true : false;

    // If no length 1, we can use standard function
    if (!nlat_is1 && !nlon_is1) {
      regrid_atmsurf_by_gp(
          field_new, atmosphere_dim, field_old, gp_lat, gp_lon);
    } else {
      if (atmosphere_dim == 2) {  // 1: No interpolation at all
        const Index n1 = gp_lat.nelem();
        field_new.resize(n1, 1);
        field_new(joker, 0) = field_old(0, 0);
      } else {
        const Index n1 = gp_lat.nelem();
        const Index n2 = gp_lon.nelem();
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
            assert(gp_lon[lon].fd[0] < 1e-6);
            field_new(joker, lon) = tmp;
          }
        } else  // 2: Just longitude interpolation
        {
          Matrix itw(n2, 2);
          interpweights(itw, gp_lon);
          Vector tmp(n2);
          interp(tmp, itw, field_old(0, joker), gp_lon);
          for (Index lat = 0; lat < n1; lat++) {
            assert(gp_lat[lat].fd[0] < 1e-6);
            field_new(lat, joker) = tmp;
          }
        }
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
  assert(itw.ncols() == 2);
  assert(p_values.nelem() == itw.nrows());

  // Local variable to store log of the pressure grid:
  Vector logpgrid(p_grid.nelem());

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
  Vector logold(old_pgrid.nelem());
  Vector lognew(new_pgrid.nelem());

  transform(logold, log, old_pgrid);
  transform(lognew, log, new_pgrid);

  gridpos(gp, logold, lognew, extpolfac);
}

void rte_pos2gridpos(GridPos& gp_p,
                     GridPos& gp_lat,
                     GridPos& gp_lon,
                     const Index& atmosphere_dim,
                     ConstVectorView p_grid,
                     ConstVectorView lat_grid,
                     ConstVectorView lon_grid,
                     ConstTensor3View z_field,
                     ConstVectorView rte_pos) {
  chk_rte_pos(atmosphere_dim, rte_pos);

  if (atmosphere_dim == 1) {
    chk_interpolation_grids(
        "Altitude interpolation", z_field(joker, 0, 0), rte_pos[0]);
    gridpos(gp_p, z_field(joker, 0, 0), rte_pos[0]);
  } else {
    // Determine z at lat/lon (z_grid) by blue interpolation
    const Index np = p_grid.nelem();
    Vector z_grid(np);
    ArrayOfGridPos agp_z, agp_lat(np);
    //
    gridpos_1to1(agp_z, p_grid);
    //
    chk_interpolation_grids("Latitude interpolation", lat_grid, rte_pos[1]);
    gridpos(gp_lat, lat_grid, rte_pos[1]);

    if (atmosphere_dim == 2) {
      for (Index i = 0; i < np; i++) {
        agp_lat[i] = gp_lat;
      }
      Matrix itw(np, 4);
      interpweights(itw, agp_z, agp_lat);
      interp(z_grid, itw, z_field(joker, joker, 0), agp_z, agp_lat);
    } else {
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
    }

    // And use z_grid to get gp_p (gp_al and gp_lon determined above)
    chk_interpolation_grids("Altitude interpolation", z_grid, rte_pos[0]);
    gridpos(gp_p, z_grid, rte_pos[0]);
  }
}

void rte_pos2gridpos(GridPos& gp_lat,
                     GridPos& gp_lon,
                     const Index& atmosphere_dim,
                     ConstVectorView lat_grid,
                     ConstVectorView lon_grid,
                     ConstVectorView rte_pos) {
  chk_rte_pos(atmosphere_dim, rte_pos);

  if (atmosphere_dim == 1) {
  } else {
    chk_interpolation_grids("Latitude interpolation", lat_grid, rte_pos[1]);
    gridpos(gp_lat, lat_grid, rte_pos[1]);

    if (atmosphere_dim == 3) {
      chk_interpolation_grids("Longitude interpolation", lon_grid, rte_pos[2]);
      gridpos(gp_lon, lon_grid, rte_pos[2]);
    }
  }
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
  const Index np = p_grid.nelem();

  assert(z.nelem() == np);
  assert(z_field.nrows() == np);
  assert(z_field.ncols() == lat_grid.nelem());

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
  const Index np = p_grid.nelem();

  assert(z.nelem() == np);
  assert(z_field.npages() == np);
  assert(z_field.nrows() == lat_grid.nelem());
  assert(z_field.ncols() == lon_grid.nelem());

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
                      const GriddedField3& complex_n,
                      const String& varname,
                      ConstVectorView f_grid,
                      ConstVectorView t_grid) {
  // Set expected order of grids
  Index gfield_fID = 0;
  Index gfield_tID = 1;
  Index gfield_compID = 2;

  // Check of complex_n
  //
  complex_n.checksize_strict();
  //
  chk_griddedfield_gridname(complex_n, gfield_fID, "Frequency");
  chk_griddedfield_gridname(complex_n, gfield_tID, "Temperature");
  chk_griddedfield_gridname(complex_n, gfield_compID, "Complex");
  //
  if (complex_n.data.ncols() != 2) {
    ostringstream os;
    os << "The data in *" << varname
       << "* must have exactly two pages. One page "
       << "each\nfor the real and imaginary part of the complex refractive index.";
  }

  // Frequency and temperature grid sizes
  const Index nf_in = complex_n.data.npages();
  const Index nt_in = complex_n.data.nrows();
  const Index nf_out = f_grid.nelem();
  const Index nt_out = t_grid.nelem();

  //Assert size of output variables
  assert(n_real.nrows() == nf_out && n_real.ncols() == nt_out);
  assert(n_imag.nrows() == nf_out && n_imag.ncols() == nt_out);

  const Vector& f_grid_in = complex_n.get_numeric_grid(gfield_fID);
  const Vector& t_grid_in = complex_n.get_numeric_grid(gfield_tID);

  // Expand/interpolate in frequency dimension
  //
  Matrix nrf(nf_out, nt_in), nif(nf_out, nt_in);
  //
  if (nf_in == 1) {
    for (Index i = 0; i < nf_out; i++) {
      nrf(i, joker) = complex_n.data(0, joker, 0);
      nif(i, joker) = complex_n.data(0, joker, 1);
    }
  } else {
    chk_interpolation_grids("Frequency interpolation", f_grid_in, f_grid);
    //
    ArrayOfGridPos gp(nf_out);
    Matrix itw(nf_out, 2);
    gridpos(gp, f_grid_in, f_grid);
    interpweights(itw, gp);
    for (Index i = 0; i < nt_in; i++) {
      interp(nrf(joker, i), itw, complex_n.data(joker, i, 0), gp);
      interp(nif(joker, i), itw, complex_n.data(joker, i, 1), gp);
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
