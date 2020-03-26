/* Copyright (C) 2006-2012 Claudia Emde <claudia.emde@dlr.de>

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
   USA. 
*/

/*!
  \file   m_disort.cc
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   2006-02-06
  
  \brief  This file contains functions to use the multiple scattering 
  program DISORT.
  
  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include "array.h"
#include "arts.h"
#include "auto_md.h"
#include "disort.h"
#include "m_general.h"
#include "math_funcs.h"
#include "messages.h"
#include "wsv_aux.h"
#include "xml_io.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void DisortCalc(Workspace& ws,
                // WS Output:
                Tensor7& cloudbox_field,
                // WS Input
                const Index& atmfields_checked,
                const Index& atmgeom_checked,
                const Index& scat_data_checked,
                const Index& cloudbox_checked,
                const Index& cloudbox_on,
                const ArrayOfIndex& cloudbox_limits,
                const Agenda& propmat_clearsky_agenda,
                const Index& atmosphere_dim,
                const Tensor4& pnd_field,
                const Tensor3& t_field,
                const Tensor3& z_field,
                const Tensor4& vmr_field,
                const Vector& p_grid,
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const Vector& f_grid,
                const Vector& za_grid,
                const Index& stokes_dim,
                const Matrix& z_surface,
                const Numeric& surface_skin_t,
                const Vector& surface_scalar_reflectivity,
                const Index& nstreams,
                const String& pfct_method,
                const Index& Npfct,
                const Index& cdisort_quiet,
                const Verbosity& verbosity) {
  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) {
    CREATE_OUT0;
    out0 << "  Cloudbox is off, DISORT calculation will be skipped.\n";
    return;
  }

  // FIXME: At the moment, combining scattering elements stored on different
  //  scattering angle grids is only possible for pfct_method 'interpolate'.

  check_disort_input(cloudbox_on,
                     atmfields_checked,
                     atmgeom_checked,
                     cloudbox_checked,
                     scat_data_checked,
                     atmosphere_dim,
                     stokes_dim,
                     cloudbox_limits,
                     scat_data,
                     za_grid,
                     nstreams,
                     pfct_method);

  init_ifield(
      cloudbox_field, f_grid, cloudbox_limits, za_grid.nelem(), stokes_dim);

  Vector albedo(f_grid.nelem(), 0.);
  Numeric btemp;

  get_disortsurf_props(
      albedo, btemp, f_grid, surface_skin_t, surface_scalar_reflectivity);

  run_cdisort(ws,
              cloudbox_field,
              f_grid,
              p_grid,
              z_field(joker, 0, 0),
              z_surface(0, 0),
              t_field(joker, 0, 0),
              vmr_field(joker, joker, 0, 0),
              pnd_field(joker, joker, 0, 0),
              scat_data,
              propmat_clearsky_agenda,
              cloudbox_limits,
              btemp,
              albedo,
              za_grid,
              nstreams,
              Npfct,
              cdisort_quiet,
              verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DisortCalcWithARTSSurface(
    Workspace& ws,
    // WS Output:
    Tensor7& cloudbox_field,
    // WS Input
    const Index& atmfields_checked,
    const Index& atmgeom_checked,
    const Index& scat_data_checked,
    const Index& cloudbox_checked,
    const Index& cloudbox_on,
    const ArrayOfIndex& cloudbox_limits,
    const Agenda& propmat_clearsky_agenda,
    const Agenda& surface_rtprop_agenda,
    const Index& atmosphere_dim,
    const Tensor4& pnd_field,
    const Tensor3& t_field,
    const Tensor3& z_field,
    const Tensor4& vmr_field,
    const Vector& p_grid,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Vector& f_grid,
    const Vector& za_grid,
    const Index& stokes_dim,
    const Index& nstreams,
    const String& pfct_method,
    const Index& Npfct,
    const Index& cdisort_quiet,
    const Verbosity& verbosity) {
  if (!cloudbox_on) {
    CREATE_OUT0;
    out0 << "  Cloudbox is off, DISORT calculation will be skipped.\n";
    return;
  }

  // FIXME: so far surface is implictly assumed at lowest atmospheric level.
  // That should be fixed (using z_surface and allowing other altitudes) at some
  // point.

  // FIXME: At the moment, combining scattering elements stored on different
  //  scattering angle grids is only possible for pfct_method 'interpolate'.

  check_disort_input(cloudbox_on,
                     atmfields_checked,
                     atmgeom_checked,
                     cloudbox_checked,
                     scat_data_checked,
                     atmosphere_dim,
                     stokes_dim,
                     cloudbox_limits,
                     scat_data,
                     za_grid,
                     nstreams,
                     pfct_method);

  init_ifield(
      cloudbox_field, f_grid, cloudbox_limits, za_grid.nelem(), stokes_dim);

  Vector albedo(f_grid.nelem(), 0.);
  Numeric btemp;

  // for now, surface at lowest atm level. later use z_surface or the like
  // for that.
  // at the moment this is only required for groundtype "A", but
  const Numeric surf_altitude = z_field(0, 0, 0);
  //const Numeric surf_altitude = z_surface(0,0);

  surf_albedoCalc(ws,
                  albedo,
                  btemp,
                  surface_rtprop_agenda,
                  f_grid,
                  za_grid,
                  surf_altitude,
                  verbosity);

  run_cdisort(ws,
              cloudbox_field,
              f_grid,
              p_grid,
              z_field(joker, 0, 0),
              surf_altitude,
              t_field(joker, 0, 0),
              vmr_field(joker, joker, 0, 0),
              pnd_field(joker, joker, 0, 0),
              scat_data,
              propmat_clearsky_agenda,
              cloudbox_limits,
              btemp,
              albedo,
              za_grid,
              nstreams,
              Npfct,
              cdisort_quiet,
              verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DisortCalcClearsky(Workspace& ws,
                        Tensor7& spectral_radiance_field,
                        const Index& atmfields_checked,
                        const Index& atmgeom_checked,
                        const Agenda& propmat_clearsky_agenda,
                        const Index& atmosphere_dim,
                        const Tensor3& t_field,
                        const Tensor3& z_field,
                        const Tensor4& vmr_field,
                        const Vector& p_grid,
                        const Vector& f_grid,
                        const Vector& za_grid,
                        const Index& stokes_dim,
                        const Matrix& z_surface,
                        const Numeric& surface_skin_t,
                        const Vector& surface_scalar_reflectivity,
                        const Index& nstreams,
                        const Index& cdisort_quiet,
                        const Verbosity& verbosity) {
  if (atmosphere_dim != 1)
    throw runtime_error(
        "For running DISORT, atmospheric dimensionality "
        "must be 1.\n");

  // Set cloudbox to cover complete atmosphere
  Index cloudbox_on;
  ArrayOfIndex cloudbox_limits;
  const Index cloudbox_checked = 1;
  //
  cloudboxSetFullAtm(cloudbox_on,
                     cloudbox_limits,
                     atmosphere_dim,
                     p_grid,
                     Vector(0),
                     Vector(0),
                     verbosity);

  // Create data matching no particles
  Tensor4 pnd_field;
  ArrayOfTensor4 dpnd_field_dx;
  ArrayOfArrayOfSingleScatteringData scat_data;
  const Index scat_data_checked = 1;
  //
  pnd_fieldZero(pnd_field,
                dpnd_field_dx,
                scat_data,
                atmosphere_dim,
                f_grid,
                cloudbox_limits,
                ArrayOfRetrievalQuantity(0),
                verbosity);

  // Call standard DISORT method
  DisortCalc(ws,
             spectral_radiance_field,
             atmfields_checked,
             atmgeom_checked,
             scat_data_checked,
             cloudbox_checked,
             cloudbox_on,
             cloudbox_limits,
             propmat_clearsky_agenda,
             atmosphere_dim,
             pnd_field,
             t_field,
             z_field,
             vmr_field,
             p_grid,
             scat_data,
             f_grid,
             za_grid,
             stokes_dim,
             z_surface,
             surface_skin_t,
             surface_scalar_reflectivity,
             nstreams,
             "median",
             181,
             cdisort_quiet,
             verbosity);
}
