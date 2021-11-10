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
                const Index& Npfct,
                const Index& only_tro,
                const Index& cdisort_quiet,
                const Verbosity& verbosity) {
  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) {
    CREATE_OUT0;
    out0 << "  Cloudbox is off, DISORT calculation will be skipped.\n";
    return;
  }

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
                     nstreams);

  init_ifield(
      cloudbox_field, f_grid, cloudbox_limits, za_grid.nelem(), 1, stokes_dim);

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
              only_tro,
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
    const Matrix& z_surface,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Vector& f_grid,
    const Vector& za_grid,
    const Index& stokes_dim,
    const Index& nstreams,
    const Index& Npfct,
    const Index& only_tro,
    const Index& cdisort_quiet,
    const Numeric& inc_angle,
    const Verbosity& verbosity) {
  if (!cloudbox_on) {
    CREATE_OUT0;
    out0 << "  Cloudbox is off, DISORT calculation will be skipped.\n";
    return;
  }

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
                     nstreams);

  init_ifield(
      cloudbox_field, f_grid, cloudbox_limits, za_grid.nelem(), 1, stokes_dim);

  Vector albedo(f_grid.nelem(), 0.);
  Numeric btemp;

  if (inc_angle<0 || inc_angle>90) {
    surf_albedoCalc(ws,
                    albedo,
                    btemp,
                    surface_rtprop_agenda,
                    f_grid,
                    za_grid,
                    z_surface(0,0),
                    verbosity);
  } else {
    surf_albedoCalcSingleAngle(ws,
                               albedo,
                               btemp,
                               surface_rtprop_agenda,
                               f_grid,
                               z_surface(0,0),
                               inc_angle);
  }
  
  run_cdisort(ws,
              cloudbox_field,
              f_grid,
              p_grid,
              z_field(joker, 0, 0),
              z_surface(0,0),
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
              only_tro,
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
                     0,
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
             181,
             0,
             cdisort_quiet,
             verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DisortCalcStar(Workspace& ws,
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
                    const Agenda& gas_scattering_agenda,
                    const Index& atmosphere_dim,
                    const Tensor4& pnd_field,
                    const Tensor3& t_field,
                    const Tensor3& z_field,
                    const Tensor4& vmr_field,
                    const Vector& p_grid,
                    const Vector& lat_true,
                    const Vector& lon_true,
                    const Vector& refellipsoid,
                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                    const ArrayOfStar& stars,
                    const Vector& f_grid,
                    const Vector& za_grid,
                    const Vector& aa_grid,
                    const Index& stokes_dim,
                    const Matrix& z_surface,
                    const Numeric& surface_skin_t,
                    const Vector& surface_scalar_reflectivity,
                    const Index& gas_scattering_do,
                    const Index& gas_scattering_output_type,
                    const Index& star_do,
                    const Index& nstreams,
                    const Index& Npfct,
                    const Index& cdisort_quiet,
                    const Verbosity& verbosity) {
  // Don't do anything if there's no cloudbox defined.
  if (!cloudbox_on) {
    CREATE_OUT0;
    out0 << "  Cloudbox is off, DISORT calculation will be skipped.\n";
    return;
  }

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
                     nstreams);

  //Check for number of stars
  ARTS_USER_ERROR_IF(stars.nelem() > 1,
                     "The simulation setup contains ",
                     stars.nelem(),
                     " stars. \n"
                     "Disort can handle only one star.")

  //allocate Varibale for direct (star) source
  Vector star_rte_los;
  Vector star_pos(3);
  Vector cloudboxtop_pos(3);
  Index star_on = star_do;
  Vector lon_grid{lon_true[0] - 0.1, lon_true[0] + 0.1};
  Vector lat_grid{lat_true[0] - 0.1, lat_true[0] + 0.1};

  //Position of star
  star_pos = {stars[0].distance, stars[0].latitude, stars[0].longitude};

  // Position of top of cloudbox
  cloudboxtop_pos = {
      z_field(cloudbox_limits[1], 0, 0), lat_true[0], lon_true[0]};

  // calculate local position of sun at top of cloudbox
  rte_losGeometricFromRtePosToRtePos2(star_rte_los,
                                      3,
                                      lat_grid,
                                      lon_grid,
                                      refellipsoid,
                                      cloudboxtop_pos,
                                      star_pos,
                                      verbosity);

  //FIXME: IF we want to be correct and include refraction, we must calculate the
  // local position of sun via ppathFromRtePos2. The question is, is this needed.

  // Check if sun is above horizon, if not switch it off
  if (star_rte_los[0] >= 90) {
    star_on = 0;

    //TODO: Add warning message that star is switched off because it is below horizon

  }

  //FIXME: Should we add a warning for low sun position near the horizon?

  Vector albedo(f_grid.nelem(), 0.);
  Numeric btemp;

  get_disortsurf_props(
      albedo, btemp, f_grid, surface_skin_t, surface_scalar_reflectivity);

  if (star_on) {
    init_ifield(cloudbox_field,
                f_grid,
                cloudbox_limits,
                za_grid.nelem(),
                aa_grid.nelem(),
                stokes_dim);
  } else {
    init_ifield(cloudbox_field,
                f_grid,
                cloudbox_limits,
                za_grid.nelem(),
                1,
                stokes_dim);
  }

  run_cdisort_star(ws,
                   cloudbox_field,
                   f_grid,
                   p_grid,
                   z_field(joker, 0, 0),
                   z_surface(0, 0),
                   t_field(joker, 0, 0),
                   vmr_field(joker, joker, 0, 0),
                   pnd_field(joker, joker, 0, 0),
                   scat_data,
                   stars,
                   propmat_clearsky_agenda,
                   gas_scattering_agenda,
                   cloudbox_limits,
                   btemp,
                   albedo,
                   za_grid,
                   aa_grid,
                   star_rte_los,
                   gas_scattering_do,
                   gas_scattering_output_type,
                   star_on,
                   nstreams,
                   Npfct,
                   cdisort_quiet,
                   verbosity);
}