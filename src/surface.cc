/* Copyright (C) 2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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

/*===========================================================================
  === File description 
  ===========================================================================*/

/**
   @file   surface.cc
   @author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   @date   2012-02-06 

   This file contains internal functions associated with the surface.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "surface.h"
#include <cmath>
#include "auto_md.h"
#include "check_input.h"
#include "complex.h"
#include "geodetic.h"
#include "math_funcs.h"
#include "matpackI.h"
#include "physics_funcs.h"
#include "workspace_ng.h"

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

Numeric calc_incang(ConstVectorView rte_los, ConstVectorView specular_los) {
  return (180 - abs(rte_los[0]) + abs(specular_los[0])) / 2;
}

Index index_of_zsurface(const Numeric& z_surface,
                        ConstVectorView z_profile) {
  Index ip = 0;
  while (z_surface >= z_profile[ip+1]) {
    ip++;
  }
  return ip;
}

void surface_calc(Matrix& iy,
                  ConstTensor3View I,
                  ConstMatrixView surface_los,
                  ConstTensor4View surface_rmatrix,
                  ConstMatrixView surface_emission) {
  // Some sizes
  const Index nf = I.nrows();
  const Index stokes_dim = I.ncols();
  const Index nlos = surface_los.nrows();

  iy = surface_emission;

  // Loop *surface_los*-es. If no such LOS, we are ready.
  if (nlos > 0) {
    for (Index ilos = 0; ilos < nlos; ilos++) {
      Vector rtmp(stokes_dim);  // Reflected Stokes vector for 1 frequency

      for (Index iv = 0; iv < nf; iv++) {
        mult(rtmp, surface_rmatrix(ilos, iv, joker, joker), I(ilos, iv, joker));
        iy(iv, joker) += rtmp;
      }
    }
  }
}

void surface_specular_R_and_b(MatrixView surface_rmatrix,
                              VectorView surface_emission,
                              const Complex& Rv,
                              const Complex& Rh,
                              const Numeric& f,
                              const Index& stokes_dim,
                              const Numeric& surface_skin_t) {
  ARTS_ASSERT(surface_rmatrix.nrows() == stokes_dim);
  ARTS_ASSERT(surface_rmatrix.ncols() == stokes_dim);
  ARTS_ASSERT(surface_emission.nelem() == stokes_dim);

  // Expressions are derived in the surface chapter in the user guide

  surface_rmatrix = 0.0;
  surface_emission = 0.0;

  Numeric B = planck(f, surface_skin_t);

  const Numeric rv = pow(abs(Rv), 2.0);
  const Numeric rh = pow(abs(Rh), 2.0);
  const Numeric rmean = (rv + rh) / 2;

  surface_rmatrix(0, 0) = rmean;
  surface_emission[0] = B * (1 - rmean);

  if (stokes_dim > 1) {
    const Numeric rdiff = (rv - rh) / 2;

    surface_rmatrix(1, 0) = rdiff;
    surface_rmatrix(0, 1) = rdiff;
    surface_rmatrix(1, 1) = rmean;
    surface_emission[1] = -B * rdiff;

    if (stokes_dim > 2) {
      const Complex a = Rh * conj(Rv);
      const Complex b = Rv * conj(Rh);
      const Numeric c = real(a + b) / 2.0;

      surface_rmatrix(2, 2) = c;

      if (stokes_dim > 3) {
        const Numeric d = imag(a - b) / 2.0;

        surface_rmatrix(2, 3) = d;
        surface_rmatrix(3, 2) = -d;
        surface_rmatrix(3, 3) = c;
      }
    }
  }
}

void surface_props_check(const Index& atmosphere_dim,
                         const Vector& lat_grid,
                         const Vector& lon_grid,
                         const Tensor3& surface_props_data,
                         const ArrayOfString& surface_props_names) {
  // Check sizes
  ARTS_USER_ERROR_IF (surface_props_data.npages() != surface_props_names.nelem(),
        "The number of pages in *surface_props_data* and "
        "length of *surface_props_names* differ.");
  // If no surface properties, then we are ready
  if (surface_props_names.nelem() == 0) {
    return;
  }
  ARTS_USER_ERROR_IF (surface_props_data.nrows() !=
      (atmosphere_dim == 1 ? 1 : lat_grid.nelem()),
                      "Row-size of *surface_props_data* not as expected.");
  ARTS_USER_ERROR_IF (surface_props_data.ncols() !=
      (atmosphere_dim <= 2 ? 1 : lon_grid.nelem()),
                      "Column-size of *surface_props_data* not as expected.");

  for (Index i = 0; i < surface_props_names.nelem(); i++) {
    ARTS_USER_ERROR_IF (surface_props_names[i].nelem() == 0,
      "Element ", i, " (0-based) of *surface_props_names* is empty.")
    for (Index j = i + 1; j < surface_props_names.nelem(); j++) {
      ARTS_USER_ERROR_IF (surface_props_names[j] == surface_props_names[i],
        "Two surface properties with same name found!\n"
        "This found for these two properties\n"
        "   index: ", i, '\n',
        "   index: ", j, '\n',
        "    name: ", surface_props_names[i])
    }
  }
}

void surface_props_interp(Vector& v,
                          const String& vname,
                          const Index& atmosphere_dim,
                          const ArrayOfGridPos& gp_lat,
                          const ArrayOfGridPos& gp_lon,
                          const Matrix& itw,
                          const Tensor3& surface_props_data,
                          const ArrayOfString& surface_props_names) {
  ARTS_ASSERT(v.nelem() == 1);
  ARTS_ASSERT(surface_props_data.npages() == surface_props_names.nelem());

  for (Index i = 0; i < surface_props_names.nelem(); i++) {
    if (surface_props_names[i] == vname) {
      interp_atmsurface_by_itw(v,
                               atmosphere_dim,
                               surface_props_data(i, joker, joker),
                               gp_lat,
                               gp_lon,
                               itw);
      return;
    }
  }

  ARTS_USER_ERROR_IF (true,
                      "The following property was requested\n"
                      "   ", vname, '\n',
                      "but it could not be found in *surface_props_names*.")
}

void dsurface_check(const ArrayOfString& surface_props_names,
                    const ArrayOfString& dsurface_names,
                    const ArrayOfTensor4 dsurface_rmatrix_dx,
                    const ArrayOfMatrix& dsurface_emission_dx) {
  const Index nq = dsurface_names.nelem();

  ARTS_USER_ERROR_IF (dsurface_rmatrix_dx.nelem() != nq,
        "The lengths of *dsurface_names* and *dsurface_rmatrix_dx* differ.");
  ARTS_USER_ERROR_IF (dsurface_emission_dx.nelem() != nq,
        "The lengths of *dsurface_names* and *dsurface_emission_dx* differ.");

  for (Index i = 0; i < nq; i++) {
    bool found = false;
    for (Index j = 0; j < surface_props_names.nelem() && !found; j++) {
      if (dsurface_names[i] == surface_props_names[j]) {
        found = true;
      }
    }
    ARTS_USER_ERROR_IF (!found,
        "String ", i, " (0-based) of *dsurface_names* is \"",
        dsurface_names[i], "\"\n"
        "but this string could not be found in *surface_props_names*.\n"
        "This is likely due to incorrect choice of quantity when\n"
        " calling *jacobianAddSurfaceQuantity*.")
  }
}
