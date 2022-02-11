/* Copyright (C) 2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>

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

/*!
  \file   m_surface.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2008-09-17

  \brief  Workspace functions associated wih the surface and its properties.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include "array.h"
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "matpack_complex.h"
#include "fastem.h"
#include "geodetic.h"
#include "interpolation.h"
#include "interpolation_lagrange.h"
#include "math_funcs.h"
#include "messages.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"
#include "surface.h"
#include "tessem.h"

extern Numeric EARTH_RADIUS;
extern const Numeric DEG2RAD;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void FastemStandAlone(Matrix& emissivity,
                      Matrix& reflectivity,
                      const Vector& f_grid,
                      const Numeric& surface_skin_t,
                      const Numeric& za,
                      const Numeric& salinity,
                      const Numeric& wind_speed,
                      const Numeric& rel_aa,
                      const Vector& transmittance,
                      const Index& fastem_version,
                      const Verbosity&) {
  const Index nf = f_grid.nelem();

  chk_if_in_range("zenith angle", za, 90, 180);
  chk_if_in_range_exclude("surface skin temperature", surface_skin_t, 260, 373);
  chk_if_in_range_exclude_high("salinity", salinity, 0, 1);
  chk_if_in_range_exclude_high("wind speed", wind_speed, 0, 100);
  chk_if_in_range("azimuth angle", rel_aa, -180, 180);
  chk_vector_length("transmittance", "f_grid", transmittance, f_grid);
  if (fastem_version < 3 || fastem_version > 6)
    throw std::runtime_error(
        "Invalid fastem version: 3 <= fastem_version <= 6");

  emissivity.resize(nf, 4);
  reflectivity.resize(nf, 4);

  const Numeric t = max(surface_skin_t, Numeric(270));

  for (Index i = 0; i < nf; i++) {
    if (f_grid[i] > 250e9)
      throw std::runtime_error("Only frequency <= 250 GHz are allowed");
    chk_if_in_range("transmittance", transmittance[i], 0, 1);

    Vector e, r;

    fastem(e,
           r,
           f_grid[i],
           za,
           t,
           salinity,
           wind_speed,
           transmittance[i],
           rel_aa,
           fastem_version);

    emissivity(i, joker) = e;
    reflectivity(i, joker) = r;
  }

  // FASTEM does not work close to the horizon (at least v6). Make sure values
  // are inside [0,1]. Then seems best to make sure that e+r=1.
  for (Index i = 0; i < nf; i++) {
    for (Index s = 0; s < 2; s++) {
      if (emissivity(i, s) > 1) {
        emissivity(i, s) = 1;
        reflectivity(i, s) = 0;
      }
      if (emissivity(i, s) < 0) {
        emissivity(i, s) = 0;
        reflectivity(i, s) = 1;
      }
      if (reflectivity(i, s) > 1) {
        emissivity(i, s) = 0;
        reflectivity(i, s) = 1;
      }
      if (reflectivity(i, s) < 0) {
        emissivity(i, s) = 1;
        reflectivity(i, s) = 0;
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void InterpGriddedField2ToPosition(Numeric& outvalue,
                                   const Index& atmosphere_dim,
                                   const Vector& lat_grid,
                                   const Vector& lat_true,
                                   const Vector& lon_true,
                                   const Vector& rtp_pos,
                                   const GriddedField2& gfield2,
                                   const Verbosity&) {
  // Set expected order of grids
  Index gfield_latID = 0;
  Index gfield_lonID = 1;

  // Basic checks and sizes
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_latlon_true(atmosphere_dim, lat_grid, lat_true, lon_true);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  gfield2.checksize_strict();
  //
  chk_griddedfield_gridname(gfield2, gfield_latID, "Latitude");
  chk_griddedfield_gridname(gfield2, gfield_lonID, "Longitude");
  //
  const Index nlat = gfield2.data.nrows();
  const Index nlon = gfield2.data.ncols();
  //
  if (nlat < 2 || nlon < 2) {
    ostringstream os;
    os << "The data in *gfield2* must span a geographical region. That is,\n"
       << "the latitude and longitude grids must have a length >= 2.";
    throw runtime_error(os.str());
  }

  const Vector& GFlat = gfield2.get_numeric_grid(gfield_latID);
  const Vector& GFlon = gfield2.get_numeric_grid(gfield_lonID);

  // Determine true geographical position
  Vector lat(1), lon(1);
  pos2true_latlon(
      lat[0], lon[0], atmosphere_dim, lat_grid, lat_true, lon_true, rtp_pos);

  // Ensure correct coverage of lon grid
  Vector lon_shifted;
  lon_shiftgrid(lon_shifted, GFlon, lon[0]);

  // Check if lat/lon we need are actually covered
  chk_if_in_range("rtp_pos.lat", lat[0], GFlat[0], GFlat[nlat - 1]);
  chk_if_in_range("rtp_pos.lon", lon[0], lon_shifted[0], lon_shifted[nlon - 1]);

  // Interpolate in lat and lon
  //
  GridPos gp_lat, gp_lon;
  gridpos(gp_lat, GFlat, lat[0]);
  gridpos(gp_lon, lon_shifted, lon[0]);
  Vector itw(4);
  interpweights(itw, gp_lat, gp_lon);
  outvalue = interp(itw, gfield2.data, gp_lat, gp_lon);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void InterpSurfaceFieldToPosition(Numeric& outvalue,
                                  const Index& atmosphere_dim,
                                  const Vector& lat_grid,
                                  const Vector& lon_grid,
                                  const Vector& rtp_pos,
                                  const Matrix& z_surface,
                                  const Matrix& field,
                                  const Verbosity& verbosity) {
  // Input checks (dummy p_grid)
  chk_atm_grids(atmosphere_dim, Vector(2, 2, -1), lat_grid, lon_grid);
  chk_atm_surface(
      "input argument *field*", field, atmosphere_dim, lat_grid, lon_grid);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  //
  const Numeric zmax = max(z_surface);
  const Numeric zmin = min(z_surface);
  const Numeric dzok = 1;
  if (rtp_pos[0] < zmin - dzok || rtp_pos[0] > zmax + dzok) {
    ostringstream os;
    os << "The given position does not match *z_surface*.\nThe altitude in "
       << "*rtp_pos* is " << rtp_pos[0] / 1e3 << " km.\n"
       << "The altitude range covered by *z_surface* is [" << zmin / 1e3 << ","
       << zmax / 1e3 << "] km.\n"
       << "One possible mistake is to mix up *rtp_pos* and *rte_los*.";
    throw runtime_error(os.str());
  }

  if (atmosphere_dim == 1) {
    outvalue = field(0, 0);
  } else {
    chk_interpolation_grids("Latitude interpolation", lat_grid, rtp_pos[1]);
    GridPos gp_lat, gp_lon;
    gridpos(gp_lat, lat_grid, rtp_pos[1]);
    if (atmosphere_dim == 3) {
      chk_interpolation_grids("Longitude interpolation", lon_grid, rtp_pos[2]);
      gridpos(gp_lon, lon_grid, rtp_pos[2]);
    }
    //
    outvalue = interp_atmsurface_by_gp(atmosphere_dim, field, gp_lat, gp_lon);
  }

  // Interpolate
  CREATE_OUT3;
  out3 << "    Result = " << outvalue << "\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iySurfaceCallAgendaX(Workspace& ws,
                          Matrix& iy,
                          ArrayOfTensor3& diy_dx,
                          const String& iy_unit,
                          const Tensor3& iy_transmittance,
                          const Index& iy_id,
                          const Index& cloudbox_on,
                          const Index& jacobian_do,
                          const Vector& f_grid,
                          const Agenda& iy_main_agenda,
                          const Vector& rtp_pos,
                          const Vector& rtp_los,
                          const Vector& rte_pos2,
                          const ArrayOfAgenda& iy_surface_agenda_array,
                          const ArrayOfIndex& surface_types,
                          const Vector& surface_types_aux,
                          const Vector& surface_types_weights,
                          const Verbosity&) {
  
  const Index ntypes = surface_types.nelem();
  if (ntypes < 1)
    throw runtime_error("*surface_types* is empty!");

  // Loop surface types and sum up
  //
  Numeric wtot = 0;
  //
  for (Index t=0; t<ntypes; t++ ) {

    if (surface_types[t] < 0)
      throw runtime_error(
                  "No element in *surface_types* is allowed to be negative.");
    if (surface_types[t] >= iy_surface_agenda_array.nelem()) {
      ostringstream os;
      os << "*iy_surface_agenda_array* has only "
         << iy_surface_agenda_array.nelem()
         << " elements,\n while you have selected element " << surface_types[t];
      throw runtime_error(os.str());
    }

    Matrix iy1;
    ArrayOfTensor3 diy_dx1;
    
    iy_surface_agenda_arrayExecute(ws,
                                   iy1,
                                   diy_dx1,
                                   surface_types[t],
                                   iy_unit,
                                   iy_transmittance,
                                   iy_id,
                                   cloudbox_on,
                                   jacobian_do,
                                   iy_main_agenda,
                                   f_grid,
                                   rtp_pos,
                                   rtp_los,
                                   rte_pos2,
                                   surface_types_aux[t],
                                   iy_surface_agenda_array);

    iy1 *= surface_types_weights[t];
    for (Index i=0; diy_dx1.nelem(); i++ )
      diy_dx1[i] *= surface_types_weights[t];
    wtot += surface_types_weights[t];;
      
    if (t==0) {
      iy     = iy1;
      diy_dx = diy_dx1;
    } else {
      iy     += iy1;
      for (Index i=0; diy_dx1.nelem(); i++ )
        diy_dx[i] += diy_dx1[i];      
    }
  }
  if (abs(wtot-1)>1e-4)
    throw runtime_error("Sum of *surface_types_weights* deviates from 1.");  
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iySurfaceFastem(Workspace& ws,
                     Matrix& iy,
                     ArrayOfTensor3& diy_dx,
                     const Tensor3& iy_transmittance,
                     const Index& iy_id,
                     const Index& jacobian_do,
                     const Index& atmosphere_dim,
                     const EnergyLevelMap& nlte_field,
                     const Index& cloudbox_on,
                     const Index& stokes_dim,
                     const Vector& f_grid,
                     const Vector& rtp_pos,
                     const Vector& rtp_los,
                     const Vector& rte_pos2,
                     const String& iy_unit,
                     const Agenda& iy_main_agenda,
                     const Numeric& surface_skin_t,
                     const Numeric& salinity,
                     const Numeric& wind_speed,
                     const Numeric& wind_direction,
                     const Index& fastem_version,
                     const Verbosity& verbosity) {
  // Most obvious input checks are performed in specular_losCalc and surfaceFastem

  // Obtain radiance and transmission for specular direction

  // Determine specular direction
  Vector specular_los, surface_normal;
  specular_losCalcNoTopography(specular_los,
                               surface_normal,
                               rtp_pos,
                               rtp_los,
                               atmosphere_dim,
                               verbosity);

  // Use iy_aux to get optical depth for downwelling radiation.
  ArrayOfString iy_aux_vars(1);
  iy_aux_vars[0] = "Optical depth";

  // Calculate iy for downwelling radiation
  // Note that iy_transmittance used here lacks surface R. Fixed below.
  //
  const Index nf = f_grid.nelem();
  Vector transmittance(nf);
  Vector geo_pos;
  ArrayOfMatrix iy_aux;
  Ppath ppath;
  //
  iy_main_agendaExecute(ws,
                        iy,
                        iy_aux,
                        ppath,
                        diy_dx,
                        geo_pos,
                        0,
                        iy_transmittance,
                        iy_aux_vars,
                        iy_id,
                        iy_unit,
                        cloudbox_on,
                        jacobian_do,
                        f_grid,
                        nlte_field,
                        rtp_pos,
                        specular_los,
                        rte_pos2,
                        iy_main_agenda);

  // Convert tau to transmissions
  for (Index i = 0; i < nf; i++) {
    transmittance[i] = exp(-iy_aux[0](i, 0));
  }

  // Call Fastem by surface_RTprop version
  //
  Matrix surface_los;
  Tensor4 surface_rmatrix;
  Matrix surface_emission;
  //
  surfaceFastem(surface_los,
                surface_rmatrix,
                surface_emission,
                atmosphere_dim,
                stokes_dim,
                f_grid,
                rtp_pos,
                rtp_los,
                surface_skin_t,
                salinity,
                wind_speed,
                wind_direction,
                transmittance,
                fastem_version,
                verbosity);

  // Add up
  //
  Tensor3 I(1, nf, stokes_dim);
  I(0, joker, joker) = iy;
  Matrix sensor_los_dummy(1, 1, 0);
  //
  surface_calc(iy, I, sensor_los_dummy, surface_rmatrix, surface_emission);

  // Adjust diy_dx, if necessary.
  // For vector cases this is a slight approximation, as the order of the
  // different transmission and reflectivities matters.
  if (iy_transmittance.npages()) {
    for (Index q = 0; q < diy_dx.nelem(); q++) {
      for (Index p = 0; p < diy_dx[q].npages(); p++) {
        for (Index i = 0; i < nf; i++) {
          Vector x = diy_dx[q](p, i, joker);
          mult(diy_dx[q](p, i, joker), surface_rmatrix(0, i, joker, joker), x);
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iySurfaceRtpropAgenda(Workspace& ws,
                           Matrix& iy,
                           ArrayOfTensor3& diy_dx,
                           Numeric& surface_skin_t,
                           Matrix& surface_los,
                           Tensor4& surface_rmatrix,
                           Matrix& surface_emission,
                           const Tensor3& iy_transmittance,
                           const Index& iy_id,
                           const Index& jacobian_do,
                           const Index& star_do,
                           const Index& atmosphere_dim,
                           const EnergyLevelMap& nlte_field,
                           const Index& cloudbox_on,
                           const Index& stokes_dim,
                           const Vector& f_grid,
                           const Vector& rtp_pos,
                           const Vector& rtp_los,
                           const Vector& rte_pos2,
                           const String& iy_unit,
                           const Agenda& iy_main_agenda,
                           const Agenda& surface_rtprop_agenda,
                           const Verbosity&) {
  // Input checks
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);

  // Call *surface_rtprop_agenda*
  surface_rtprop_agendaExecute(ws,
                               surface_skin_t,
                               surface_emission,
                               surface_los,
                               surface_rmatrix,
                               f_grid,
                               rtp_pos,
                               rtp_los,
                               surface_rtprop_agenda);

  // Check output of *surface_rtprop_agenda*
  const Index nlos = surface_los.nrows();
  const Index nf = f_grid.nelem();
  //
  if (nlos)  // if 0, blackbody ground and not all checks are needed
  {
    if (surface_los.ncols() != rtp_los.nelem())
      throw runtime_error("Number of columns in *surface_los* is not correct.");
    if (nlos != surface_rmatrix.nbooks())
      throw runtime_error(
          "Mismatch in size of *surface_los* and *surface_rmatrix*.");
    if (surface_rmatrix.npages() != nf)
      throw runtime_error(
          "Mismatch in size of *surface_rmatrix* and *f_grid*.");
    if (surface_rmatrix.nrows() != stokes_dim ||
        surface_rmatrix.ncols() != stokes_dim)
      throw runtime_error(
          "Mismatch between size of *surface_rmatrix* and *stokes_dim*.");
  }
  if (surface_emission.ncols() != stokes_dim)
    throw runtime_error(
        "Mismatch between size of *surface_emission* and *stokes_dim*.");
  if (surface_emission.nrows() != nf)
    throw runtime_error("Mismatch in size of *surface_emission* and f_grid*.");

  // Variable to hold down-welling radiation
  Tensor3 I(nlos, nf, stokes_dim);

  ArrayOfString iy_aux_var(0);
  if (star_do) {
    iy_aux_var.resize(1);
    iy_aux_var[0] = "Direct radiation";
  }

  // Loop *surface_los*-es. If no such LOS, we are ready.
  if (nlos > 0) {
    for (Index ilos = 0; ilos < nlos; ilos++) {
      Vector los = surface_los(ilos, joker);

      // Include surface reflection matrix in *iy_transmittance*
      // If iy_transmittance is empty, this is interpreted as the
      // variable is not needed.
      //
      Tensor3 iy_trans_new;
      //
      if (iy_transmittance.npages()) {
        iy_transmittance_mult(iy_trans_new,
                             iy_transmittance,
                             surface_rmatrix(ilos, joker, joker, joker));
      }

      // Calculate downwelling radiation for LOS ilos
      //
      {
        ArrayOfMatrix iy_aux;
        Ppath ppath;
        Vector geo_pos;
        Index iy_id_new = iy_id + ilos + 1;
        iy_main_agendaExecute(ws,
                              iy,
                              iy_aux,
                              ppath,
                              diy_dx,
                              geo_pos,
                              0,
                              iy_trans_new,
                              iy_aux_var,
                              iy_id_new,
                              iy_unit,
                              cloudbox_on,
                              jacobian_do,
                              f_grid,
                              nlte_field,
                              rtp_pos,
                              los,
                              rte_pos2,
                              iy_main_agenda);

        //For the case that a star is present and the los is towards a star, we
        //subtract the direct radiation, so that only the diffuse radiation is considered here.
        //If star is within los iy_aux[0] is the incoming and attenuated direct (star)
        //radiation at the surface. Otherwise it is zero.
        if (star_do){
          iy-=iy_aux[0];
        }

      }

      if (iy.ncols() != stokes_dim || iy.nrows() != nf) {
        ostringstream os;
        os << "The size of *iy* returned from *" << iy_main_agenda.name()
           << "* is\n"
           << "not correct:\n"
           << "  expected size = [" << nf << "," << stokes_dim << "]\n"
           << "  size of iy    = [" << iy.nrows() << "," << iy.ncols() << "]\n";
        throw runtime_error(os.str());
      }

      I(ilos, joker, joker) = iy;
    }
  }

  // Add up
  surface_calc(iy, I, surface_los, surface_rmatrix, surface_emission);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iySurfaceRtpropCalc(Workspace& ws,
                         Matrix& iy,
                         ArrayOfTensor3& diy_dx,
                         const Matrix& surface_los,
                         const Tensor4& surface_rmatrix,
                         const Matrix& surface_emission,
                         const ArrayOfString& dsurface_names,
                         const ArrayOfTensor4& dsurface_rmatrix_dx,
                         const ArrayOfMatrix& dsurface_emission_dx,
                         const Tensor3& iy_transmittance,
                         const Index& iy_id,
                         const Index& jacobian_do,
                         const ArrayOfRetrievalQuantity& jacobian_quantities,
                         const Index& atmosphere_dim,
                         const EnergyLevelMap& nlte_field,                         
                         const Index& cloudbox_on,
                         const Index& stokes_dim,
                         const Vector& f_grid,
                         const Vector& rtp_pos,
                         const Vector& rtp_los,
                         const Vector& rte_pos2,
                         const String& iy_unit,
                         const Agenda& iy_main_agenda,
                         const Verbosity&) {
  // Input checks
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);

  // Check provided surface rtprop variables
  const Index nlos = surface_los.nrows();
  const Index nf = f_grid.nelem();
  //
  if (nlos)  // if 0, blackbody ground and not all checks are needed
  {
    if (surface_los.ncols() != rtp_los.nelem())
      throw runtime_error("Number of columns in *surface_los* is not correct.");
    if (nlos != surface_rmatrix.nbooks())
      throw runtime_error(
          "Mismatch in size of *surface_los* and *surface_rmatrix*.");
    if (surface_rmatrix.npages() != nf)
      throw runtime_error(
          "Mismatch in size of *surface_rmatrix* and *f_grid*.");
    if (surface_rmatrix.nrows() != stokes_dim ||
        surface_rmatrix.ncols() != stokes_dim)
      throw runtime_error(
          "Mismatch between size of *surface_rmatrix* and *stokes_dim*.");
  }
  if (surface_emission.ncols() != stokes_dim)
    throw runtime_error(
        "Mismatch between size of *surface_emission* and *stokes_dim*.");
  if (surface_emission.nrows() != nf)
    throw runtime_error("Mismatch in size of *surface_emission* and f_grid*.");

  // Variable to hold down-welling radiation
  Tensor3 I(nlos, nf, stokes_dim);

  // Loop *surface_los*-es.
  if (nlos > 0) {
    for (Index ilos = 0; ilos < nlos; ilos++) {
      Vector los = surface_los(ilos, joker);

      // Include surface reflection matrix in *iy_transmittance*
      // If iy_transmittance is empty, this is interpreted as the
      // variable is not needed.
      //
      Tensor3 iy_trans_new;
      //
      if (iy_transmittance.npages()) {
        iy_transmittance_mult(iy_trans_new,
                             iy_transmittance,
                             surface_rmatrix(ilos, joker, joker, joker));
      }

      // Calculate downwelling radiation for LOS ilos
      //
      {
        ArrayOfMatrix iy_aux;
        Ppath ppath;
        Vector geo_pos;
        iy_main_agendaExecute(ws,
                              iy,
                              iy_aux,
                              ppath,
                              diy_dx,
                              geo_pos,
                              0,
                              iy_trans_new,
                              ArrayOfString(0),
                              iy_id,
                              iy_unit,
                              cloudbox_on,
                              jacobian_do,
                              f_grid,
                              nlte_field,
                              rtp_pos,
                              los,
                              rte_pos2,
                              iy_main_agenda);
      }

      if (iy.ncols() != stokes_dim || iy.nrows() != nf) {
        ostringstream os;
        os << "The size of *iy* returned from *" << iy_main_agenda.name()
           << "* is\n"
           << "not correct:\n"
           << "  expected size = [" << nf << "," << stokes_dim << "]\n"
           << "  size of iy    = [" << iy.nrows() << "," << iy.ncols() << "]\n";
        throw runtime_error(os.str());
      }

      I(ilos, joker, joker) = iy;
    }
  }

  // Add up
  surface_calc(iy, I, surface_los, surface_rmatrix, surface_emission);

  // Surface Jacobians
  if (jacobian_do && dsurface_names.nelem()) {
    // Loop dsurface_names
    for (Index i = 0; i < dsurface_names.nelem(); i++) {
      // Error if derivatives not calculated
      // Or should we accept this?
      if (dsurface_emission_dx[i].empty() || dsurface_rmatrix_dx[i].empty()) {
        ostringstream os;
        os << "The derivatives for surface quantity: " << dsurface_names[i]
           << "\nwere not calculated by *iy_surface_agenda*.\n"
           << "That is, *dsurface_emission_dx* and/or *dsurface_rmatrix_dx*\n"
           << "are empty.";
        throw runtime_error(os.str());
      } else {
        // Find index among jacobian quantities
        Index ihit = -1;
        for (Index j = 0; j < jacobian_quantities.nelem() && ihit < 0; j++) {
          if (dsurface_names[i] == jacobian_quantities[j].Subtag()) {
            ihit = j;
          }
        }
        ARTS_ASSERT(ihit >= 0);
        // Derivative, as observed at the surface
        Matrix diy_dpos0, diy_dpos;
        surface_calc(diy_dpos0,
                     I,
                     surface_los,
                     dsurface_rmatrix_dx[i],
                     dsurface_emission_dx[i]);
        // Weight with transmission to sensor
        iy_transmittance_mult(diy_dpos, iy_transmittance, diy_dpos0);
        // Put into diy_dx
        diy_from_pos_to_rgrids(diy_dx[ihit],
                               jacobian_quantities[ihit],
                               diy_dpos,
                               atmosphere_dim,
                               rtp_pos);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void specular_losCalcNoTopography(Vector& specular_los,
                                  Vector& surface_normal,
                                  const Vector& rtp_pos,
                                  const Vector& rtp_los,
                                  const Index& atmosphere_dim,
                                  const Verbosity&) {
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);

  surface_normal.resize(max(Index(1), atmosphere_dim - 1));
  specular_los.resize(max(Index(1), atmosphere_dim - 1));

  if (atmosphere_dim == 1) {
    surface_normal[0] = 0;
    if (rtp_los[0] < 90) {
      throw runtime_error(
          "Invalid zenith angle. The zenith angle corresponds "
          "to observe the surface from below.");
    }
    specular_los[0] = 180 - rtp_los[0];
  }

  else if (atmosphere_dim == 2) {
    specular_los[0] = sign(rtp_los[0]) * 180 - rtp_los[0];
    surface_normal[0] = 0;
  }

  else if (atmosphere_dim == 3) {
    specular_los[0] = 180 - rtp_los[0];
    specular_los[1] = rtp_los[1];
    surface_normal[0] = 0;
    surface_normal[1] = 0;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void specular_losCalc(Vector& specular_los,
                      Vector& surface_normal,
                      const Vector& rtp_pos,
                      const Vector& rtp_los,
                      const Index& atmosphere_dim,
                      const Vector& lat_grid,
                      const Vector& lon_grid,
                      const Vector& refellipsoid,
                      const Matrix& z_surface,
                      const Index& ignore_surface_slope,
                      const Verbosity& verbosity) {
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);
  chk_if_in_range("ignore_surface_slope", ignore_surface_slope, 0, 1);

  // Use special function if there is no slope, or it is ignored
  if (atmosphere_dim == 1 || ignore_surface_slope) {
    specular_losCalcNoTopography(specular_los,
                                 surface_normal,
                                 rtp_pos,
                                 rtp_los,
                                 atmosphere_dim,
                                 verbosity);
    return;
  }

  surface_normal.resize(max(Index(1), atmosphere_dim - 1));
  specular_los.resize(max(Index(1), atmosphere_dim - 1));

  if (atmosphere_dim == 2) {
    chk_interpolation_grids("Latitude interpolation", lat_grid, rtp_pos[1]);
    GridPos gp_lat;
    gridpos(gp_lat, lat_grid, rtp_pos[1]);
    Numeric c1;  // Radial slope of the surface
    plevel_slope_2d(
        c1, lat_grid, refellipsoid, z_surface(joker, 0), gp_lat, rtp_los[0]);
    Vector itw(2);
    interpweights(itw, gp_lat);
    const Numeric rv_surface = refell2d(refellipsoid, lat_grid, gp_lat) +
                               interp(itw, z_surface(joker, 0), gp_lat);
    surface_normal[0] = -plevel_angletilt(rv_surface, c1);
    if (abs(rtp_los[0] - surface_normal[0]) < 90) {
      throw runtime_error(
          "Invalid zenith angle. The zenith angle corresponds "
          "to observe the surface from below.");
    }
    specular_los[0] =
        sign(rtp_los[0]) * 180 - rtp_los[0] + 2 * surface_normal[0];
  }

  else {
    // Calculate surface normal in South-North direction
    chk_interpolation_grids("Latitude interpolation", lat_grid, rtp_pos[1]);
    chk_interpolation_grids("Longitude interpolation", lon_grid, rtp_pos[2]);
    GridPos gp_lat, gp_lon;
    gridpos(gp_lat, lat_grid, rtp_pos[1]);
    gridpos(gp_lon, lon_grid, rtp_pos[2]);
    Numeric c1, c2;
    plevel_slope_3d(
        c1, c2, lat_grid, lon_grid, refellipsoid, z_surface, gp_lat, gp_lon, 0);
    Vector itw(4);
    interpweights(itw, gp_lat, gp_lon);
    const Numeric rv_surface = refell2d(refellipsoid, lat_grid, gp_lat) +
                               interp(itw, z_surface, gp_lat, gp_lon);
    const Numeric zaSN = 90 - plevel_angletilt(rv_surface, c1);
    // The same for East-West
    plevel_slope_3d(c1,
                    c2,
                    lat_grid,
                    lon_grid,
                    refellipsoid,
                    z_surface,
                    gp_lat,
                    gp_lon,
                    90);
    const Numeric zaEW = 90 - plevel_angletilt(rv_surface, c1);
    // Convert to Cartesian, and determine normal by cross-product
    Vector tangentSN(3), tangentEW(3), normal(3);
    zaaa2cart(tangentSN[0], tangentSN[1], tangentSN[2], zaSN, 0);
    zaaa2cart(tangentEW[0], tangentEW[1], tangentEW[2], zaEW, 90);
    cross3(normal, tangentSN, tangentEW);
    // Convert rtp_los to cartesian and flip direction
    Vector di(3);
    zaaa2cart(di[0], di[1], di[2], rtp_los[0], rtp_los[1]);
    di *= -1;
    // Set LOS normal vector
    cart2zaaa(
        surface_normal[0], surface_normal[1], normal[0], normal[1], normal[2]);
    if (abs(rtp_los[0] - surface_normal[0]) < 90) {
      throw runtime_error(
          "Invalid zenith angle. The zenith angle corresponds "
          "to observe the surface from below.");
    }
    // Specular direction is 2(dn*di)dn-di, where dn is the normal vector
    Vector speccart(3);
    const Numeric fac = 2 * (normal * di);
    for (Index i = 0; i < 3; i++) {
      speccart[i] = fac * normal[i] - di[i];
    }
    cart2zaaa(specular_los[0],
              specular_los[1],
              speccart[0],
              speccart[1],
              speccart[2]);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceBlackbody(Matrix& surface_los,
                      Tensor4& surface_rmatrix,
                      Matrix& surface_emission,
                      const Index& atmosphere_dim,
                      const Vector& f_grid,
                      const Index& stokes_dim,
                      const Vector& rtp_pos,
                      const Vector& rtp_los,
                      const Numeric& surface_skin_t,
                      const Verbosity& verbosity) {
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);
  chk_not_negative("surface_skin_t", surface_skin_t);

  CREATE_OUT2;
  out2 << "  Sets variables to model a blackbody surface with a temperature "
       << " of " << surface_skin_t << " K.\n";

  surface_los.resize(0, 0);
  surface_rmatrix.resize(0, 0, 0, 0);

  const Index nf = f_grid.nelem();

  Vector b(nf);
  planck(b, f_grid, surface_skin_t);

  surface_emission.resize(nf, stokes_dim);
  surface_emission = 0.0;

  for (Index iv = 0; iv < nf; iv++) {
    surface_emission(iv, 0) = b[iv];
    for (Index is = 1; is < stokes_dim; is++) {
      surface_emission(iv, is) = 0;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFastem(Matrix& surface_los,
                   Tensor4& surface_rmatrix,
                   Matrix& surface_emission,
                   const Index& atmosphere_dim,
                   const Index& stokes_dim,
                   const Vector& f_grid,
                   const Vector& rtp_pos,
                   const Vector& rtp_los,
                   const Numeric& surface_skin_t,
                   const Numeric& salinity,
                   const Numeric& wind_speed,
                   const Numeric& wind_direction,
                   const Vector& transmittance,
                   const Index& fastem_version,
                   const Verbosity& verbosity) {
  // Input checks
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);
  chk_if_in_range("wind_direction", wind_direction, -180, 180);

  const Index nf = f_grid.nelem();

  // Determine specular direction
  Vector specular_los, surface_normal;
  specular_losCalcNoTopography(specular_los,
                               surface_normal,
                               rtp_pos,
                               rtp_los,
                               atmosphere_dim,
                               verbosity);

  // Determine relative azimuth
  //
  // According to email from James Hocking, UkMet:
  // All azimuth angles are defined as being measured clockwise from north
  // (i.e. if the projection onto the Earth's surface of the view path lies due
  // north the azimuth angle is 0 and if it lies due east the azimuth angle is
  // 90 degrees). The relative azimuth is the wind direction (azimuth) angle
  // minus the satellite azimuth angle.
  Numeric rel_azimuth = wind_direction;  // Always true for 1D
  if (atmosphere_dim == 2 && rtp_los[0] < 0) {
    rel_azimuth -= 180;
    resolve_lon(rel_azimuth, -180, 180);
  } else if (atmosphere_dim == 3) {
    rel_azimuth -= rtp_los[1];
    resolve_lon(rel_azimuth, -180, 180);
  }

  // Call FASTEM
  Matrix emissivity, reflectivity;
  FastemStandAlone(emissivity,
                   reflectivity,
                   f_grid,
                   surface_skin_t,
                   abs(rtp_los[0]),
                   salinity,
                   wind_speed,
                   rel_azimuth,
                   transmittance,
                   fastem_version,
                   verbosity);

  // Set surface_los
  surface_los.resize(1, specular_los.nelem());
  surface_los(0, joker) = specular_los;

  // Surface emission
  //
  Vector b(nf);
  planck(b, f_grid, surface_skin_t);
  //
  surface_emission.resize(nf, stokes_dim);
  for (Index i = 0; i < nf; i++) {
    // I
    surface_emission(i, 0) = b[i] * 0.5 * (emissivity(i, 0) + emissivity(i, 1));
    // Q
    if (stokes_dim >= 2) {
      surface_emission(i, 1) =
          b[i] * 0.5 * (emissivity(i, 0) - emissivity(i, 1));
    }
    // U and V
    for (Index j = 2; j < stokes_dim; j++) {
      surface_emission(i, j) = b[i] * emissivity(i, j);
    }
  }

  // Surface reflectivity matrix
  //
  surface_rmatrix.resize(1, nf, stokes_dim, stokes_dim);
  surface_rmatrix = 0.0;
  for (Index iv = 0; iv < nf; iv++) {
    surface_rmatrix(0, iv, 0, 0) =
        0.5 * (reflectivity(iv, 0) + reflectivity(iv, 1));
    if (stokes_dim >= 2) {
      surface_rmatrix(0, iv, 0, 1) =
          0.5 * (reflectivity(iv, 0) - reflectivity(iv, 1));
      ;
      surface_rmatrix(0, iv, 1, 0) = surface_rmatrix(0, iv, 0, 1);
      surface_rmatrix(0, iv, 1, 1) = surface_rmatrix(0, iv, 0, 0);

      for (Index i = 2; i < stokes_dim; i++) {
        surface_rmatrix(0, iv, i, i) = surface_rmatrix(0, iv, 0, 0);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceMapToLinearPolarisation(Matrix& surface_emission,
                                    Tensor4& surface_rmatrix,
                                    const Index& stokes_dim,
                                    const Numeric& pol_angle,
                                    const Verbosity& ) {
  ARTS_USER_ERROR_IF (stokes_dim != 1,
        "You should only use this method where the main calculations are "
        "done with *stokes_dim* set to 1.");
  
  Index local_stokes = surface_emission.ncols();
  ARTS_USER_ERROR_IF (local_stokes < 2,
        "This method requires that the input surface proporties match a Stokes "
        "dimension of 2-4. Incoming *surface_emission* matches stokes_dim=1.");

  const Index nf = surface_emission.nrows();
  const Index nlos = surface_rmatrix.nbooks();
  Matrix se = surface_emission;
  Tensor4 sr = surface_rmatrix;
  surface_emission.resize(nf, 1);
  surface_rmatrix.resize(nlos, nf, 1, 1);
  
  if (local_stokes == 2) {
    const Numeric alpha = DEG2RAD * pol_angle; 
    const Numeric c2 = pow( cos(alpha), 2 );
    const Numeric s2 = pow( sin(alpha), 2 );
    for (Index f=0; f<nf; ++f) {
      const Numeric bv = se(f,0) + se(f,1); // Note that we have to multiply with 2
      const Numeric bh = se(f,0) - se(f,1); // as we place a single pol as I
      surface_emission(f,0) = c2*bv + s2*bh;
      for (Index l=0; l<nlos; ++l) {
        const Numeric rv = sr(l,f,0,0) + sr(l,f,1,0); 
        const Numeric rh = sr(l,f,0,0) - sr(l,f,1,0);
        surface_rmatrix(l,f,0,0) = c2*rv + s2*rh;
      }
    }
  } else {
    Matrix Cm, Cs, Lp, Lm;
    mueller_stokes2modif(Cm, local_stokes);
    mueller_modif2stokes(Cs, local_stokes);
    mueller_rotation(Lp, local_stokes, pol_angle);
    mueller_rotation(Lm, local_stokes, -pol_angle);
    Matrix Mleft(local_stokes,local_stokes), Mright(local_stokes,local_stokes);
    mult(Mleft, Cm, Lp);
    mult(Mright, Lm, Cs);
    //
    Vector Vr(local_stokes);
    Matrix Tmp(local_stokes,local_stokes), Mr(local_stokes,local_stokes);
    //
    for (Index f=0; f<nf; ++f) {
      mult(Vr, Mleft, se(f,joker) );        // Note that we have to multiply with 2
      surface_emission(f,0) = 2.0 * Vr[0];  // as we place a single pol as I
      for (Index l=0; l<nlos; ++l) {
        mult(Tmp, sr(l,f,joker,joker), Mright);
        mult(Mr, Mleft, Tmp);
        surface_rmatrix(l,f,0,0) = Tmp(0,0);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceTelsem(Matrix& surface_los,
                   Tensor4& surface_rmatrix,
                   Matrix& surface_emission,
                   const Index& atmosphere_dim,
                   const Index& stokes_dim,
                   const Vector& f_grid,
                   const Vector& lat_grid,
                   const Vector& lat_true,
                   const Vector& lon_true,
                   const Vector& rtp_pos,
                   const Vector& rtp_los,
                   const Vector& specular_los,
                   const Numeric& surface_skin_t,
                   const TelsemAtlas& atlas,
                   const Numeric& r_min,
                   const Numeric& r_max,
                   const Numeric& d_max,
                   const Verbosity& verbosity) {
  // Input checks
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);
  chk_rte_los(atmosphere_dim, specular_los);
  chk_if_in_range_exclude(
      "surface skin temperature", surface_skin_t, 190.0, 373.0);

  const Index nf = f_grid.nelem();
  Matrix surface_rv_rh(nf, 2);

  // Lookup cell in atlas.

  Numeric lat, lon;
  pos2true_latlon(
      lat, lon, atmosphere_dim, lat_grid, lat_true, lon_true, rtp_pos);
  chk_if_in_range("Latitude input to TELSEM2", lat, -90.0, 90.0);

  lon = fmod(lon, 360.0);
  if (lon < 0) {
    lon += 360.0;
  }
  chk_if_in_range("Longitude input to TELSEM2", lon, 0.0, 360.0);

  Index cellnumber = atlas.calc_cellnum(lat, lon);
  // Check if cell is in atlas.
  if (!atlas.contains(cellnumber)) {
    if (d_max <= 0.0) {
      throw std::runtime_error(
          "Given coordinates are not contained in "
          " TELSEM atlas. To enable nearest neighbour"
          "interpolation set *d_max* to a positive "
          "value.");
    } else {
      cellnumber = atlas.calc_cellnum_nearest_neighbor(lat, lon);
      Numeric lat_nn, lon_nn;
      std::tie(lat_nn, lon_nn) = atlas.get_coordinates(cellnumber);
      Numeric d = sphdist(lat, lon, lat_nn, lon_nn);
      if (d > d_max) {
        std::ostringstream out{};
        out << "Distance of nearest neighbour exceeds provided limit (";
        out << d << " > " << d_max << ").";
        throw std::runtime_error(out.str());
      }
    }
  }

  Index class1 = atlas.get_class1(cellnumber);
  Index class2 = atlas.get_class2(cellnumber);

  Vector emis_h = atlas.get_emis_h(cellnumber);
  Vector emis_v = atlas.get_emis_v(cellnumber);
  Numeric e_v, e_h;

  // Incidence angle
  const Numeric theta = calc_incang(rtp_los, specular_los);
  ARTS_ASSERT(theta <= 90);

  for (Index i = 0; i < nf; ++i) {
    if (f_grid[i] < 5e9)
      throw std::runtime_error("Only frequency >= 5 GHz are allowed");
    if (f_grid[i] > 900e9)
      throw std::runtime_error("Only frequency <= 900 GHz are allowed");

    // For frequencies 700 <= f <= 900 GHz use 700 GHz to avoid extrapolation
    // outside of TELSEM specifications.
    Numeric f = std::min(f_grid[i], 700e9) * 1e-9;
    std::tie(e_v, e_h) =
        atlas.emis_interp(theta, f, class1, class2, emis_v, emis_h);

    surface_rv_rh(i, 0) = min(max(1.0 - e_v, r_min), r_max);
    surface_rv_rh(i, 1) = min(max(1.0 - e_h, r_min), r_max);
  }

  surfaceFlatRvRh(surface_los,
                  surface_rmatrix,
                  surface_emission,
                  f_grid,
                  stokes_dim,
                  atmosphere_dim,
                  rtp_pos,
                  rtp_los,
                  specular_los,
                  surface_skin_t,
                  surface_rv_rh,
                  verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceTessem(Matrix& surface_los,
                   Tensor4& surface_rmatrix,
                   Matrix& surface_emission,
                   const Index& atmosphere_dim,
                   const Index& stokes_dim,
                   const Vector& f_grid,
                   const Vector& rtp_pos,
                   const Vector& rtp_los,
                   const Numeric& surface_skin_t,
                   const TessemNN& net_h,
                   const TessemNN& net_v,
                   const Numeric& salinity,
                   const Numeric& wind_speed,
                   const Verbosity& verbosity) {
  // Input checks
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);
  chk_if_in_range_exclude(
      "surface skin temperature", surface_skin_t, 260.0, 373.0);
  chk_if_in_range_exclude_high("salinity", salinity, 0, 1);
  chk_if_in_range_exclude_high("wind speed", wind_speed, 0, 100);

  // Determine specular direction
  Vector specular_los, surface_normal;
  specular_losCalcNoTopography(specular_los,
                               surface_normal,
                               rtp_pos,
                               rtp_los,
                               atmosphere_dim,
                               verbosity);

  // TESSEM in and out
  //
  Vector out(2);
  VectorView e_h = out[Range(0, 1)];
  VectorView e_v = out[Range(1, 1)];
  //
  Vector in(5);
  in[1] = 180.0 - abs(rtp_los[0]);
  in[2] = wind_speed;
  in[3] = surface_skin_t;
  in[4] = salinity;

  // Get Rv and Rh
  //
  const Index nf = f_grid.nelem();
  Matrix surface_rv_rh(nf, 2);
  //
  for (Index i = 0; i < nf; ++i) {
    if (f_grid[i] < 5e9)
      throw std::runtime_error("Only frequency >= 5 GHz are allowed");
    if (f_grid[i] > 900e9)
      throw std::runtime_error("Only frequency <= 900 GHz are allowed");

    in[0] = f_grid[i];

    tessem_prop_nn(e_h, net_h, in);
    tessem_prop_nn(e_v, net_v, in);

    surface_rv_rh(i, 0) = min(max(1 - e_v[0], (Numeric)0), (Numeric)1);
    surface_rv_rh(i, 1) = min(max(1 - e_h[0], (Numeric)0), (Numeric)1);
  }

  surfaceFlatRvRh(surface_los,
                  surface_rmatrix,
                  surface_emission,
                  f_grid,
                  stokes_dim,
                  atmosphere_dim,
                  rtp_pos,
                  rtp_los,
                  specular_los,
                  surface_skin_t,
                  surface_rv_rh,
                  verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatRefractiveIndex(Matrix& surface_los,
                                Tensor4& surface_rmatrix,
                                Matrix& surface_emission,
                                const Vector& f_grid,
                                const Index& stokes_dim,
                                const Index& atmosphere_dim,
                                const Vector& rtp_pos,
                                const Vector& rtp_los,
                                const Vector& specular_los,
                                const Numeric& surface_skin_t,
                                const GriddedField3& surface_complex_refr_index,
                                const Verbosity& verbosity) {
  CREATE_OUT2;
  CREATE_OUT3;

  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);
  chk_rte_los(atmosphere_dim, specular_los);
  chk_not_negative("surface_skin_t", surface_skin_t);

  // Interpolate *surface_complex_refr_index*
  //
  const Index nf = f_grid.nelem();
  //
  Matrix n_real(nf, 1), n_imag(nf, 1);
  //
  complex_n_interp(n_real,
                   n_imag,
                   surface_complex_refr_index,
                   "surface_complex_refr_index",
                   f_grid,
                   Vector(1, surface_skin_t));

  out2 << "  Sets variables to model a flat surface\n";
  out3 << "     surface temperature: " << surface_skin_t << " K.\n";

  surface_los.resize(1, specular_los.nelem());
  surface_los(0, joker) = specular_los;

  surface_emission.resize(nf, stokes_dim);
  surface_rmatrix.resize(1, nf, stokes_dim, stokes_dim);

  // Incidence angle
  const Numeric incang = calc_incang(rtp_los, specular_los);
  ARTS_ASSERT(incang <= 90);

  // Complex (amplitude) reflection coefficients
  Complex Rv, Rh;

  for (Index iv = 0; iv < nf; iv++) {
    // Set n2 (refractive index of surface medium)
    Complex n2(n_real(iv, 0), n_imag(iv, 0));

    // Amplitude reflection coefficients
    fresnel(Rv, Rh, Numeric(1.0), n2, incang);

    // Fill reflection matrix and emission vector
    surface_specular_R_and_b(surface_rmatrix(0, iv, joker, joker),
                             surface_emission(iv, joker),
                             Rv,
                             Rh,
                             f_grid[iv],
                             stokes_dim,
                             surface_skin_t);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatReflectivity(Matrix& surface_los,
                             Tensor4& surface_rmatrix,
                             Matrix& surface_emission,
                             const Vector& f_grid,
                             const Index& stokes_dim,
                             const Index& atmosphere_dim,
                             const Vector& rtp_pos,
                             const Vector& rtp_los,
                             const Vector& specular_los,
                             const Numeric& surface_skin_t,
                             const Tensor3& surface_reflectivity,
                             const Verbosity& verbosity) {
  CREATE_OUT2;

  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);
  chk_rte_los(atmosphere_dim, specular_los);
  chk_not_negative("surface_skin_t", surface_skin_t);

  const Index nf = f_grid.nelem();

  if (surface_reflectivity.nrows() != stokes_dim &&
      surface_reflectivity.ncols() != stokes_dim) {
    ostringstream os;
    os << "The number of rows and columnss in *surface_reflectivity* must\n"
       << "match *stokes_dim*."
       << "\n stokes_dim : " << stokes_dim
       << "\n number of rows in *surface_reflectivity* : "
       << surface_reflectivity.nrows()
       << "\n number of columns in *surface_reflectivity* : "
       << surface_reflectivity.ncols() << "\n";
    throw runtime_error(os.str());
  }

  if (surface_reflectivity.npages() != nf &&
      surface_reflectivity.npages() != 1) {
    ostringstream os;
    os << "The number of pages in *surface_reflectivity* should\n"
       << "match length of *f_grid* or be 1."
       << "\n length of *f_grid* : " << nf
       << "\n dimension of *surface_reflectivity* : "
       << surface_reflectivity.npages() << "\n";
    throw runtime_error(os.str());
  }

  out2 << "  Sets variables to model a flat surface\n";

  surface_los.resize(1, specular_los.nelem());
  surface_los(0, joker) = specular_los;

  surface_emission.resize(nf, stokes_dim);
  surface_rmatrix.resize(1, nf, stokes_dim, stokes_dim);

  Matrix R, IR(stokes_dim, stokes_dim);

  Vector b(nf);
  planck(b, f_grid, surface_skin_t);

  Vector B(stokes_dim, 0);

  for (Index iv = 0; iv < nf; iv++) {
    if (iv == 0 || surface_reflectivity.npages() > 1) {
      R = surface_reflectivity(iv, joker, joker);
      for (Index i = 0; i < stokes_dim; i++) {
        for (Index j = 0; j < stokes_dim; j++) {
          if (i == j) {
            IR(i, j) = 1 - R(i, j);
          } else {
            IR(i, j) = -R(i, j);
          }
        }
      }
    }

    surface_rmatrix(0, iv, joker, joker) = R;

    B[0] = b[iv];
    mult(surface_emission(iv, joker), IR, B);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatRvRh(Matrix& surface_los,
                     Tensor4& surface_rmatrix,
                     Matrix& surface_emission,
                     const Vector& f_grid,
                     const Index& stokes_dim,
                     const Index& atmosphere_dim,
                     const Vector& rtp_pos,
                     const Vector& rtp_los,
                     const Vector& specular_los,
                     const Numeric& surface_skin_t,
                     const Matrix& surface_rv_rh,
                     const Verbosity&) {
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);
  chk_rte_los(atmosphere_dim, specular_los);
  chk_not_negative("surface_skin_t", surface_skin_t);

  const Index nf = f_grid.nelem();

  if (surface_rv_rh.ncols() != 2) {
    ostringstream os;
    os << "The number of columns in *surface_rv_rh* must be two,\n"
       << "but the actual number of columns is " << surface_rv_rh.ncols()
       << "\n";
    throw runtime_error(os.str());
  }

  if (surface_rv_rh.nrows() != nf && surface_rv_rh.nrows() != 1) {
    ostringstream os;
    os << "The number of rows in *surface_rv_rh* should\n"
       << "match length of *f_grid* or be 1."
       << "\n length of *f_grid* : " << nf
       << "\n rows in *surface_rv_rh* : " << surface_rv_rh.nrows() << "\n";
    throw runtime_error(os.str());
  }

  if (min(surface_rv_rh) < 0 || max(surface_rv_rh) > 1) {
    throw runtime_error("All values in *surface_rv_rh* must be inside [0,1].");
  }

  surface_los.resize(1, specular_los.nelem());
  surface_los(0, joker) = specular_los;

  surface_emission.resize(nf, stokes_dim);
  surface_rmatrix.resize(1, nf, stokes_dim, stokes_dim);

  surface_emission = 0;
  surface_rmatrix = 0;

  Vector b(nf);
  planck(b, f_grid, surface_skin_t);

  Numeric rmean = 0.0, rdiff = 0.0;

  for (Index iv = 0; iv < nf; iv++) {
    if (iv == 0 || surface_rv_rh.nrows() > 1) {
      rmean = 0.5 * (surface_rv_rh(iv, 0) + surface_rv_rh(iv, 1));
      rdiff = 0.5 * (surface_rv_rh(iv, 0) - surface_rv_rh(iv, 1));
    }

    surface_emission(iv, 0) = (1.0 - rmean) * b[iv];
    surface_rmatrix(0, iv, 0, 0) = rmean;

    if (stokes_dim > 1) {
      surface_emission(iv, 1) = -rdiff * b[iv];

      surface_rmatrix(0, iv, 0, 1) = rdiff;
      surface_rmatrix(0, iv, 1, 0) = rdiff;
      surface_rmatrix(0, iv, 1, 1) = rmean;

      for (Index i = 2; i < stokes_dim; i++) {
        surface_rmatrix(0, iv, i, i) = rmean;
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceFlatScalarReflectivity(Matrix& surface_los,
                                   Tensor4& surface_rmatrix,
                                   Matrix& surface_emission,
                                   const Vector& f_grid,
                                   const Index& stokes_dim,
                                   const Index& atmosphere_dim,
                                   const Vector& rtp_pos,
                                   const Vector& rtp_los,
                                   const Vector& specular_los,
                                   const Numeric& surface_skin_t,
                                   const Vector& surface_scalar_reflectivity,
                                   const Verbosity&) {
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);
  chk_rte_los(atmosphere_dim, specular_los);
  chk_not_negative("surface_skin_t", surface_skin_t);

  const Index nf = f_grid.nelem();

  if (surface_scalar_reflectivity.nelem() != nf &&
      surface_scalar_reflectivity.nelem() != 1) {
    ostringstream os;
    os << "The number of elements in *surface_scalar_reflectivity* should\n"
       << "match length of *f_grid* or be 1."
       << "\n length of *f_grid* : " << nf
       << "\n length of *surface_scalar_reflectivity* : "
       << surface_scalar_reflectivity.nelem() << "\n";
    throw runtime_error(os.str());
  }

  if (min(surface_scalar_reflectivity) < 0 ||
      max(surface_scalar_reflectivity) > 1) {
    throw runtime_error(
        "All values in *surface_scalar_reflectivity* must be inside [0,1].");
  }

  surface_los.resize(1, specular_los.nelem());
  surface_los(0, joker) = specular_los;

  surface_emission.resize(nf, stokes_dim);
  surface_rmatrix.resize(1, nf, stokes_dim, stokes_dim);

  surface_emission = 0;
  surface_rmatrix = 0;

  Vector b(nf);
  planck(b, f_grid, surface_skin_t);

  Numeric r = 0.0;

  for (Index iv = 0; iv < nf; iv++) {
    if (iv == 0 || surface_scalar_reflectivity.nelem() > 1) {
      r = surface_scalar_reflectivity[iv];
    }

    surface_emission(iv, 0) = (1.0 - r) * b[iv];
    surface_rmatrix(0, iv, 0, 0) = r;
    for (Index i = 1; i < stokes_dim; i++) {
      surface_rmatrix(0, iv, i, i) = r;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceLambertianSimple(Matrix& surface_los,
                             Tensor4& surface_rmatrix,
                             Matrix& surface_emission,
                             const Vector& f_grid,
                             const Index& stokes_dim,
                             const Index& atmosphere_dim,
                             const Vector& rtp_pos,
                             const Vector& rtp_los,
                             const Vector& surface_normal,
                             const Numeric& surface_skin_t,
                             const Vector& surface_scalar_reflectivity,
                             const Index& lambertian_nza,
                             const Numeric& za_pos,
                             const Verbosity&) {
  const Index nf = f_grid.nelem();

  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);
  chk_not_negative("surface_skin_t", surface_skin_t);
  chk_if_in_range("za_pos", za_pos, 0, 1);

  if (surface_scalar_reflectivity.nelem() != nf &&
      surface_scalar_reflectivity.nelem() != 1) {
    ostringstream os;
    os << "The number of elements in *surface_scalar_reflectivity* should\n"
       << "match length of *f_grid* or be 1."
       << "\n length of *f_grid* : " << nf
       << "\n length of *surface_scalar_reflectivity* : "
       << surface_scalar_reflectivity.nelem() << "\n";
    throw runtime_error(os.str());
  }

  if (min(surface_scalar_reflectivity) < 0 ||
      max(surface_scalar_reflectivity) > 1) {
    throw runtime_error(
        "All values in *surface_scalar_reflectivity* must be inside [0,1].");
  }

  // Allocate and init everything to zero
  //
  surface_los.resize(lambertian_nza, rtp_los.nelem());
  surface_rmatrix.resize(lambertian_nza, nf, stokes_dim, stokes_dim);
  surface_emission.resize(nf, stokes_dim);
  //
  surface_los = 0.0;
  surface_rmatrix = 0.0;
  surface_emission = 0.0;

  // Help variables
  //
  const Numeric dza = (90.0 - abs(surface_normal[0])) / (Numeric)lambertian_nza;
  const Vector za_lims(0.0, lambertian_nza + 1, dza);

  // surface_los
  for (Index ip = 0; ip < lambertian_nza; ip++) {
    surface_los(ip, 0) = za_lims[ip] + za_pos * dza;
    if (atmosphere_dim == 2) {
      if (rtp_los[0] < 0) {
        surface_los(ip, 0) *= -1.0;
      }

    } else if (atmosphere_dim == 3) {
      surface_los(ip, 1) = rtp_los[1];
    }
  }

  Vector b(nf);
  planck(b, f_grid, surface_skin_t);

  // Loop frequencies and set remaining values
  //
  Numeric r = 0.0;
  //
  for (Index iv = 0; iv < nf; iv++) {
    // Get reflectivity
    if (iv == 0 || surface_scalar_reflectivity.nelem() > 1) {
      r = surface_scalar_reflectivity[iv];
    }

    // surface_rmatrix:
    // Only element (0,0) is set to be non-zero. This follows VDISORT
    // that refers to: K. L. Coulson, Polarization and Intensity of Light in
    // the Atmosphere (1989), page 229
    // (Thanks to Michael Kahnert for providing this information!)
    // Update: Is the above for a later edition? We have not found a copy of
    // that edition. In a 1988 version of the book, the relevant page seems
    // to be 232.
    for (Index ip = 0; ip < lambertian_nza; ip++) {
      const Numeric w =
          r * 0.5 *
          (cos(2 * DEG2RAD * za_lims[ip]) - cos(2 * DEG2RAD * za_lims[ip + 1]));
      surface_rmatrix(ip, iv, 0, 0) = w;
    }

    // surface_emission
    surface_emission(iv, 0) = (1 - r) * b[iv];
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceSemiSpecularBy3beams(Workspace& ws,
                                 Numeric& surface_skin_t,
                                 Matrix& surface_los,
                                 Tensor4& surface_rmatrix,
                                 Matrix& surface_emission,
                                 const Index& atmosphere_dim,
                                 const Vector& f_grid,
                                 const Vector& rtp_pos,
                                 const Vector& rtp_los,
                                 const Agenda& surface_rtprop_sub_agenda,
                                 const Numeric& specular_factor,
                                 const Numeric& dza,
                                 const Verbosity&) {
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);

  // Checks of GIN variables
  if (specular_factor > 1 || specular_factor < 1.0 / 3.0)
    throw runtime_error("The valid range for *specular_factor* is [1/3,1].");
  if (dza > 45 || dza <= 0)
    throw runtime_error("The valid range for *dza* is ]0,45].");

  // Obtain data for specular direction
  //
  Matrix los1, emission1;
  Tensor4 rmatrix1;
  //
  surface_rtprop_sub_agendaExecute(ws,
                                   surface_skin_t,
                                   emission1,
                                   los1,
                                   rmatrix1,
                                   f_grid,
                                   rtp_pos,
                                   rtp_los,
                                   surface_rtprop_sub_agenda);
  if (los1.nrows() != 1)
    throw runtime_error(
        "*surface_rtprop_sub_agenda* must return data "
        "describing a specular surface.");

  // Standard number of beams. Set to 2 if try/catch below fails
  Index nbeams = 3;

  // Test if some lower zenith angle works.
  // It will fail if a higher za results in looking at the surface from below.
  //
  Matrix los2, emission2;
  Tensor4 rmatrix2;
  //
  Numeric skin_t_dummy;
  Numeric dza_try = dza;
  bool failed = true;
  while (failed && dza_try > 0) {
    try {
      Vector los_new = rtp_los;
      los_new[0] -=
          sign(rtp_los[0]) * dza_try;  // Sign to also handle 2D negative za
      adjust_los(los_new, atmosphere_dim);
      surface_rtprop_sub_agendaExecute(ws,
                                       skin_t_dummy,
                                       emission2,
                                       los2,
                                       rmatrix2,
                                       f_grid,
                                       rtp_pos,
                                       los_new,
                                       surface_rtprop_sub_agenda);
      failed = false;
    } catch (const runtime_error& e) {
      dza_try -= 1.0;
    }
  }
  if (failed) {
    nbeams = 2;
  }

  // Allocate output WSVs
  //
  surface_emission.resize(emission1.nrows(), emission1.ncols());
  surface_emission = 0;
  surface_los.resize(nbeams, los1.ncols());
  surface_rmatrix.resize(
      nbeams, rmatrix1.npages(), rmatrix1.nrows(), rmatrix1.ncols());

  // Put in specular direction at index 1
  //
  Numeric w;
  if (nbeams == 3) {
    w = specular_factor;
  } else {
    w = specular_factor + (1.0 - specular_factor) / 2.0;
  }
  //
  surface_los(1, joker) = los1(0, joker);
  for (Index p = 0; p < rmatrix1.npages(); p++) {
    for (Index r = 0; r < rmatrix1.nrows(); r++) {
      surface_emission(p, r) += w * emission1(p, r);
      for (Index c = 0; c < rmatrix1.ncols(); c++) {
        surface_rmatrix(1, p, r, c) = w * rmatrix1(0, p, r, c);
      }
    }
  }

  // Put in lower za as index 2, if worked
  //
  w = (1.0 - specular_factor) / 2.0;
  //
  if (nbeams == 3) {
    surface_los(2, joker) = los2(0, joker);
    for (Index p = 0; p < rmatrix2.npages(); p++) {
      for (Index r = 0; r < rmatrix2.nrows(); r++) {
        surface_emission(p, r) += w * emission2(p, r);
        for (Index c = 0; c < rmatrix1.ncols(); c++) {
          surface_rmatrix(2, p, r, c) = w * rmatrix2(0, p, r, c);
        }
      }
    }
  }

  // Do higher za and put in as index 0 (reusing variables for beam 2)
  //
  Vector los_new = rtp_los;
  los_new[0] += sign(rtp_los[0]) * dza;  // Sign to also handle 2D negative za
  adjust_los(los_new, atmosphere_dim);
  surface_rtprop_sub_agendaExecute(ws,
                                   skin_t_dummy,
                                   emission2,
                                   los2,
                                   rmatrix2,
                                   f_grid,
                                   rtp_pos,
                                   los_new,
                                   surface_rtprop_sub_agenda);
  //
  surface_los(0, joker) = los2(0, joker);
  for (Index p = 0; p < rmatrix2.npages(); p++) {
    for (Index r = 0; r < rmatrix2.nrows(); r++) {
      surface_emission(p, r) += w * emission2(p, r);
      for (Index c = 0; c < rmatrix1.ncols(); c++) {
        surface_rmatrix(0, p, r, c) = w * rmatrix2(0, p, r, c);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surfaceSplitSpecularTo3beams(Matrix& surface_los,
                                  Tensor4& surface_rmatrix,
                                  const Index& atmosphere_dim,
                                  const Vector& rtp_pos,
                                  const Vector& rtp_los,
                                  const Numeric& specular_factor,
                                  const Numeric& dza,
                                  const Verbosity&) {
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);

  // Check that input surface data are of specular type
  if (surface_los.nrows() != 1)
    throw runtime_error(
        "Input surface data must be of specular type. That is, "
        "*surface_los* must contain a single direction.");
  if (surface_rmatrix.nbooks() != 1)
    throw runtime_error(
        "*surface_rmatrix* describes a different number of "
        "directions than *surface_los*.");

  // Checks of GIN variables
  if (specular_factor > 1 || specular_factor < 1.0 / 3.0)
    throw runtime_error("The valid range for *specular_factor* is [1/3,1].");
  if (dza > 45 || dza <= 0)
    throw runtime_error("The valid range for *dza* is ]0,45].");

  // Make copies of input data
  const Matrix los1 = surface_los;
  const Tensor4 rmatrix1 = surface_rmatrix;

  // Use abs(za) in all expressions below, to also handle 2D

  // Calculate highest possible za for downwelling radiation, with 1 degree
  // margin to the surface. This can be derived from za in rtp_los and surface_los.
  // (The directions in surface_los are not allowed to point into the surface)
  const Numeric za_max = 89 + (180 - abs(los1(0, 0)) - abs(rtp_los[0])) / 2.0;

  // Number of downwelling beams
  Index nbeams = 3;
  if (abs(los1(0, 0)) > za_max) {
    nbeams = 2;
  }

  // New los-s
  //
  surface_los.resize(nbeams, los1.ncols());
  //
  for (Index r = 0; r < nbeams; r++) {
    surface_los(r, 0) = ((Numeric)r - 1.0) * dza + abs(los1(0, 0));
    if (r == 2 && surface_los(r, 0) > za_max) {
      surface_los(r, 0) = za_max;
    }
    for (Index c = 1; c < los1.ncols(); c++) {
      surface_los(r, c) = los1(0, c);
    }
  }

  // New rmatrix
  //
  surface_rmatrix.resize(
      nbeams, rmatrix1.npages(), rmatrix1.nrows(), rmatrix1.ncols());
  //
  for (Index b = 0; b < nbeams; b++) {
    Numeric w;
    if (b == 1 && nbeams == 3)  // Specular direction with nbeams==3
    {
      w = specular_factor;
    } else if (b == 1)  // Specular direction with nbeams==2
    {
      w = specular_factor + (1.0 - specular_factor) / 2.0;
    } else  // Side directions
    {
      w = (1.0 - specular_factor) / 2.0;
    }

    for (Index p = 0; p < rmatrix1.npages(); p++) {
      for (Index r = 0; r < rmatrix1.nrows(); r++) {
        for (Index c = 0; c < rmatrix1.ncols(); c++) {
          surface_rmatrix(b, p, r, c) = w * rmatrix1(0, p, r, c);
        }
      }
    }
  }

  // Handle sign of za
  if (atmosphere_dim == 1) {
    // We only need to make sure that first direction has positive za
    surface_los(0, 0) = abs(surface_los(0, 0));
  } else if (atmosphere_dim == 2) {
    // Change sign if specular direction has za < 0
    if (los1(0, 0) < 0) {
      for (Index r = 0; r < rmatrix1.nrows(); r++) {
        surface_los(r, 0) = -surface_los(r, 0);
      }
    }
  } else if (atmosphere_dim == 1) {
    // We only need to make sure that first direction has positive za
    if (surface_los(0, 0) < 0) {
      surface_los(0, 0) = -surface_los(0, 0);
      surface_los(0, 1) += 180;
      if (surface_los(0, 1) > 180) {
        surface_los(0, 1) -= 360;
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_complex_refr_indexFromGriddedField5(
    GriddedField3& surface_complex_refr_index,
    const Index& atmosphere_dim,
    const Vector& lat_grid,
    const Vector& lat_true,
    const Vector& lon_true,
    const Vector& rtp_pos,
    const GriddedField5& complex_n_field,
    const Verbosity&) {
  // Set expected order of grids
  Index gfield_fID = 0;
  Index gfield_tID = 1;
  Index gfield_compID = 2;
  Index gfield_latID = 3;
  Index gfield_lonID = 4;

  // Basic checks and sizes
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_latlon_true(atmosphere_dim, lat_grid, lat_true, lon_true);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  complex_n_field.checksize_strict();
  //
  chk_griddedfield_gridname(complex_n_field, gfield_fID, "Frequency");
  chk_griddedfield_gridname(complex_n_field, gfield_tID, "Temperature");
  chk_griddedfield_gridname(complex_n_field, gfield_compID, "Complex");
  chk_griddedfield_gridname(complex_n_field, gfield_latID, "Latitude");
  chk_griddedfield_gridname(complex_n_field, gfield_lonID, "Longitude");
  //
  const Index nf = complex_n_field.data.nshelves();
  const Index nt = complex_n_field.data.nbooks();
  const Index nn = complex_n_field.data.npages();
  const Index nlat = complex_n_field.data.nrows();
  const Index nlon = complex_n_field.data.ncols();
  //
  if (nlat < 2 || nlon < 2) {
    ostringstream os;
    os << "The data in *complex_refr_index_field* must span a geographical "
       << "region. That is,\nthe latitude and longitude grids must have a "
       << "length >= 2.";
    throw runtime_error(os.str());
  }
  //
  if (nn != 2) {
    ostringstream os;
    os << "The data in *complex_refr_index_field* must have exactly two "
       << "pages. One page each\nfor the real and imaginary part of the "
       << "complex refractive index.";
    throw runtime_error(os.str());
  }

  const Vector& GFlat = complex_n_field.get_numeric_grid(gfield_latID);
  const Vector& GFlon = complex_n_field.get_numeric_grid(gfield_lonID);

  // Determine true geographical position
  Vector lat(1), lon(1);
  pos2true_latlon(
      lat[0], lon[0], atmosphere_dim, lat_grid, lat_true, lon_true, rtp_pos);

  // Ensure correct coverage of lon grid
  Vector lon_shifted;
  lon_shiftgrid(lon_shifted, GFlon, lon[0]);

  // Check if lat/lon we need are actually covered
  chk_if_in_range("rtp_pos.lat", lat[0], GFlat[0], GFlat[nlat - 1]);
  chk_if_in_range("rtp_pos.lon", lon[0], lon_shifted[0], lon_shifted[nlon - 1]);

  // Size and fills grids of *surface_complex_refr_index*
  surface_complex_refr_index.resize(nf, nt, 2);
  surface_complex_refr_index.set_grid_name(0, "Frequency");
  surface_complex_refr_index.set_grid(
      0, complex_n_field.get_numeric_grid(gfield_fID));
  surface_complex_refr_index.set_grid_name(1, "Temperature");
  surface_complex_refr_index.set_grid(
      1, complex_n_field.get_numeric_grid(gfield_tID));
  surface_complex_refr_index.set_grid_name(2, "Complex");
  surface_complex_refr_index.set_grid(2, {"real", "imaginary"});

  // Interpolate in lat and lon
  //
  GridPos gp_lat, gp_lon;
  gridpos(gp_lat, GFlat, lat[0]);
  gridpos(gp_lon, lon_shifted, lon[0]);
  Vector itw(4);
  interpweights(itw, gp_lat, gp_lon);
  //
  for (Index iv = 0; iv < nf; iv++) {
    for (Index it = 0; it < nt; it++) {
      surface_complex_refr_index.data(iv, it, 0) = interp(
          itw, complex_n_field.data(iv, it, 0, joker, joker), gp_lat, gp_lon);
      surface_complex_refr_index.data(iv, it, 1) = interp(
          itw, complex_n_field.data(iv, it, 1, joker, joker), gp_lat, gp_lon);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_reflectivityFromGriddedField6(Tensor3& surface_reflectivity,
                                           const Index& stokes_dim,
                                           const Vector& f_grid,
                                           const Index& atmosphere_dim,
                                           const Vector& lat_grid,
                                           const Vector& lat_true,
                                           const Vector& lon_true,
                                           const Vector& rtp_pos,
                                           const Vector& rtp_los,
                                           const GriddedField6& r_field,
                                           const Verbosity&) {
  // Basic checks and sizes
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  chk_latlon_true(atmosphere_dim, lat_grid, lat_true, lon_true);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);
  r_field.checksize_strict();
  chk_griddedfield_gridname(r_field, 0, "Frequency");
  chk_griddedfield_gridname(r_field, 1, "Stokes element");
  chk_griddedfield_gridname(r_field, 2, "Stokes element");
  chk_griddedfield_gridname(r_field, 3, "Incidence angle");
  chk_griddedfield_gridname(r_field, 4, "Latitude");
  chk_griddedfield_gridname(r_field, 5, "Longitude");
  //
  const Index nf_in = r_field.data.nvitrines();
  const Index ns2 = r_field.data.nshelves();
  const Index ns1 = r_field.data.nbooks();
  const Index nza = r_field.data.npages();
  const Index nlat = r_field.data.nrows();
  const Index nlon = r_field.data.ncols();
  //
  if (nlat < 2 || nlon < 2) {
    ostringstream os;
    os << "The data in *r_field* must span a geographical region. That is,\n"
       << "the latitude and longitude grids must have a length >= 2.";
    throw runtime_error(os.str());
  }
  //
  if (nza < 2) {
    ostringstream os;
    os << "The data in *r_field* must span a range of zenith angles. That\n"
       << "is the zenith angle grid must have a length >= 2.";
    throw runtime_error(os.str());
  }
  if (ns1 < stokes_dim || ns2 < stokes_dim || ns1 > 4 || ns2 > 4) {
    ostringstream os;
    os << "The \"Stokes dimensions\" must have a size that is >= "
       << "*stokes_dim* (but not exceeding 4).";
    throw runtime_error(os.str());
  }

  // Determine true geographical position
  Vector lat(1), lon(1);
  pos2true_latlon(
      lat[0], lon[0], atmosphere_dim, lat_grid, lat_true, lon_true, rtp_pos);

  // Ensure correct coverage of lon grid
  Vector lon_shifted;
  lon_shiftgrid(lon_shifted, r_field.get_numeric_grid(5), lon[0]);

  // Interpolate in lat and lon
  Tensor4 r_f_za(nf_in, stokes_dim, stokes_dim, nza);
  {
    chk_interpolation_grids(
        "Latitude interpolation", r_field.get_numeric_grid(4), lat[0]);
    chk_interpolation_grids("Longitude interpolation", lon_shifted, lon[0]);
    GridPos gp_lat, gp_lon;
    gridpos(gp_lat, r_field.get_numeric_grid(4), lat[0]);
    gridpos(gp_lon, lon_shifted, lon[0]);
    Vector itw(4);
    interpweights(itw, gp_lat, gp_lon);
    for (Index iv = 0; iv < nf_in; iv++) {
      for (Index iz = 0; iz < nza; iz++) {
        for (Index is1 = 0; is1 < stokes_dim; is1++) {
          for (Index is2 = 0; is2 < stokes_dim; is2++) {
            r_f_za(iv, is1, is2, iz) =
                interp(itw,
                       r_field.data(iv, is1, is2, iz, joker, joker),
                       gp_lat,
                       gp_lon);
          }
        }
      }
    }
  }

  // Interpolate in incidence angle, cubic if possible
  Tensor3 r_f(nf_in, stokes_dim, stokes_dim);
  Index order = 3;
  if (nza < 4) {
    order = 1;
  }
  {
    chk_interpolation_grids(
      "Incidence angle interpolation", r_field.get_numeric_grid(3), 180 - rtp_los[0]);
    const LagrangeInterpolation lag(0, 180 - rtp_los[0], r_field.get_numeric_grid(3), order);
    const auto itw = interpweights(lag);
    //
    for (Index i = 0; i < nf_in; i++) {
      for (Index is1 = 0; is1 < stokes_dim; is1++) {
        for (Index is2 = 0; is2 < stokes_dim; is2++) {
          r_f(i, is1, is2) = interp(r_f_za(i, is1, is2, joker), itw, lag);
        }
      }
    }
  }

  // Extract or interpolate in frequency
  //
  if (nf_in == 1) {
    surface_reflectivity = r_f;
  } else {
    chk_interpolation_grids(
        "Frequency interpolation", r_field.get_numeric_grid(0), f_grid);
    const Index nf_out = f_grid.nelem();
    surface_reflectivity.resize(nf_out, stokes_dim, stokes_dim);
    //
    ArrayOfGridPos gp(nf_out);
    Matrix itw(nf_out, 2);
    gridpos(gp, r_field.get_numeric_grid(0), f_grid);
    interpweights(itw, gp);
    for (Index is1 = 0; is1 < stokes_dim; is1++) {
      for (Index is2 = 0; is2 < stokes_dim; is2++) {
        interp(surface_reflectivity(joker, is1, is2),
               itw,
               r_f(joker, is1, is2),
               gp);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_scalar_reflectivityFromGriddedField4(
    Vector& surface_scalar_reflectivity,
    const Index& stokes_dim,
    const Vector& f_grid,
    const Index& atmosphere_dim,
    const Vector& lat_grid,
    const Vector& lat_true,
    const Vector& lon_true,
    const Vector& rtp_pos,
    const Vector& rtp_los,
    const GriddedField4& r_field,
    const Verbosity&) {
  // Basic checks and sizes
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_if_in_range("stokes_dim", stokes_dim, 1, 1);
  chk_latlon_true(atmosphere_dim, lat_grid, lat_true, lon_true);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  chk_rte_los(atmosphere_dim, rtp_los);
  r_field.checksize_strict();
  chk_griddedfield_gridname(r_field, 0, "Frequency");
  chk_griddedfield_gridname(r_field, 1, "Incidence angle");
  chk_griddedfield_gridname(r_field, 2, "Latitude");
  chk_griddedfield_gridname(r_field, 3, "Longitude");
  //
  const Index nf_in = r_field.data.nbooks();
  const Index nza = r_field.data.npages();
  const Index nlat = r_field.data.nrows();
  const Index nlon = r_field.data.ncols();
  //
  if (nlat < 2 || nlon < 2) {
    ostringstream os;
    os << "The data in *r_field* must span a geographical region. That is,\n"
       << "the latitude and longitude grids must have a length >= 2.";
    throw runtime_error(os.str());
  }
  //
  if (nza < 2) {
    ostringstream os;
    os << "The data in *r_field* must span a range of zenith angles. That\n"
       << "is the zenith angle grid must have a length >= 2.";
    throw runtime_error(os.str());
  }

  // Determine true geographical position
  Vector lat(1), lon(1);
  pos2true_latlon(
      lat[0], lon[0], atmosphere_dim, lat_grid, lat_true, lon_true, rtp_pos);

  // Ensure correct coverage of lon grid
  Vector lon_shifted;
  lon_shiftgrid(lon_shifted, r_field.get_numeric_grid(3), lon[0]);

  // Interpolate in lat and lon
  Matrix r_f_za(nf_in, nza);
  {
    chk_interpolation_grids(
        "Latitude interpolation", r_field.get_numeric_grid(2), lat[0]);
    chk_interpolation_grids("Longitude interpolation", lon_shifted, lon[0]);
    GridPos gp_lat, gp_lon;
    gridpos(gp_lat, r_field.get_numeric_grid(2), lat[0]);
    gridpos(gp_lon, lon_shifted, lon[0]);
    Vector itw(4);
    interpweights(itw, gp_lat, gp_lon);
    for (Index iv = 0; iv < nf_in; iv++) {
      for (Index iz = 0; iz < nza; iz++) {
        r_f_za(iv, iz) =
            interp(itw, r_field.data(iv, iz, joker, joker), gp_lat, gp_lon);
      }
    }
  }

  // Interpolate in incidence angle, cubic if possible
  Vector r_f(nf_in);
  Index order = 3;
  if (nza < 4) {
    order = 1;
  }
  {
    chk_interpolation_grids(
      "Incidence angle interpolation", r_field.get_numeric_grid(1), 180 - rtp_los[0]);
    const LagrangeInterpolation lag(0, 180 - rtp_los[0], r_field.get_numeric_grid(1), order);
    const auto itw=interpweights(lag);
    for (Index i = 0; i < nf_in; i++) {
      r_f[i] = interp(r_f_za(i, joker), itw, lag);
    }
  }

  // Extract or interpolate in frequency
  //
  if (nf_in == 1) {
    surface_scalar_reflectivity.resize(1);
    surface_scalar_reflectivity[0] = r_f[0];
  } else {
    chk_interpolation_grids(
        "Frequency interpolation", r_field.get_numeric_grid(0), f_grid);
    const Index nf_out = f_grid.nelem();
    surface_scalar_reflectivity.resize(nf_out);
    //
    ArrayOfGridPos gp(nf_out);
    Matrix itw(nf_out, 2);
    gridpos(gp, r_field.get_numeric_grid(0), f_grid);
    interpweights(itw, gp);
    interp(surface_scalar_reflectivity, itw, r_f, gp);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_scalar_reflectivityFromSurface_rmatrix(
    Vector& surface_scalar_reflectivity,
    const Tensor4& surface_rmatrix,
    const Verbosity&) {
  const Index nf = surface_rmatrix.npages();
  const Index nlos = surface_rmatrix.nbooks();

  surface_scalar_reflectivity.resize(nf);
  surface_scalar_reflectivity = 0;

  for (Index i = 0; i < nf; i++) {
    for (Index l = 0; l < nlos; l++) {
      surface_scalar_reflectivity[i] += surface_rmatrix(l, i, 0, 0);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_typeInterpTypeMask(ArrayOfIndex& surface_types,
                                Vector& surface_types_aux,
                                Vector& surface_types_weights,
                                const Index& atmosphere_dim,
                                const Vector& lat_grid,
                                const Vector& lat_true,
                                const Vector& lon_true,
                                const Vector& rtp_pos,
                                const GriddedField2& surface_type_mask,
                                const String& method,
                                const Verbosity&) {
  // Set expected order of grids
  Index gfield_latID = 0;
  Index gfield_lonID = 1;

  // Basic checks and sizes
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_latlon_true(atmosphere_dim, lat_grid, lat_true, lon_true);
  chk_rte_pos(atmosphere_dim, rtp_pos);
  surface_type_mask.checksize_strict();
  //
  chk_griddedfield_gridname(surface_type_mask, gfield_latID, "Latitude");
  chk_griddedfield_gridname(surface_type_mask, gfield_lonID, "Longitude");
  //
  const Index nlat = surface_type_mask.data.nrows();
  const Index nlon = surface_type_mask.data.ncols();
  //
  if (nlat < 2 || nlon < 2) {
    ostringstream os;
    os << "The data in *surface_type_mask* must span a geographical "
       << "region. Accordingly,\nthe latitude and longitude grids must "
       << "both have a length >= 2.";
    throw runtime_error(os.str());
  }

  const Vector& GFlat = surface_type_mask.get_numeric_grid(gfield_latID);
  const Vector& GFlon = surface_type_mask.get_numeric_grid(gfield_lonID);

  // Determine true geographical position
  Vector lat(1), lon(1);
  pos2true_latlon(
      lat[0], lon[0], atmosphere_dim, lat_grid, lat_true, lon_true, rtp_pos);

  // Ensure correct coverage of lon grid
  Vector lon_shifted;
  lon_shiftgrid(lon_shifted, GFlon, lon[0]);

  // Check if lat/lon we need are actually covered
  chk_if_in_range("rtp_pos.lat", lat[0], GFlat[0], GFlat[nlat - 1]);
  chk_if_in_range("rtp_pos.lon", lon[0], lon_shifted[0], lon_shifted[nlon - 1]);

  // Grid positions
  GridPos gp_lat, gp_lon;
  gridpos(gp_lat, GFlat, lat[0]);
  gridpos(gp_lon, lon_shifted, lon[0]);
  
  if (method == "nearest" ) {
    // Extract closest point
    Index ilat, ilon;
    if (gp_lat.fd[0] < 0.5) {
      ilat = gp_lat.idx;
    } else {
      ilat = gp_lat.idx + 1;
    }
    if (gp_lon.fd[0] < 0.5) {
      ilon = gp_lon.idx;
    } else {
      ilon = gp_lon.idx + 1;
    }
    //
    surface_types.resize(1);
    surface_types_aux.resize(1);
    surface_types_weights.resize(1);
    surface_types[0] = (Index)floor(surface_type_mask.data(ilat, ilon));
    surface_types_aux[0] = surface_type_mask.data(ilat, ilon) -
                                              Numeric(surface_types[0]);
    surface_types_weights[0] = 1.0;
  }
  
  else if (method == "linear" ) {
    // Determine types involved
    ArrayOfIndex types0(4), types(4);
    Index ntypes = 1;
    //
    types0[0] = (Index)floor(surface_type_mask.data(gp_lat.idx,gp_lon.idx));
    types0[1] = (Index)floor(surface_type_mask.data(gp_lat.idx,gp_lon.idx+1));
    types0[2] = (Index)floor(surface_type_mask.data(gp_lat.idx+1,gp_lon.idx));
    types0[3] = (Index)floor(surface_type_mask.data(gp_lat.idx+1,gp_lon.idx+1));
    //
    types[0] = types0[0];
    for (Index t=1; t<4; t++) {
      bool unique = true;
      for (Index n=0; n<ntypes && unique; n++) 
        if (types0[t] == types[n] ) { unique = false; }
      if (unique) {
        types[ntypes] = types0[t];
        ntypes += 1;
      }
    }
    //
    surface_types.resize(ntypes);
    surface_types_aux.resize(ntypes);
    surface_types_weights.resize(ntypes);
    
    // Interpolation weights
    Vector itw(4);
    interpweights(itw, gp_lat, gp_lon);

    // Determine weight for each type, and make an interpolation of aux value
    for (Index t=0; t<ntypes; t++ ) {
      Numeric wtot = 0, auxtot = 0;
      Numeric ntype = (Numeric)types[t];
      if (types[t] == types0[0]) { wtot += itw[0];
        auxtot += itw[0]*(surface_type_mask.data(gp_lat.idx,gp_lon.idx)-ntype);
      };
      if (types[t] == types0[1]) { wtot += itw[1];
        auxtot += itw[1]*(surface_type_mask.data(gp_lat.idx,gp_lon.idx+1)-ntype);
      };
      if (types[t] == types0[2]) { wtot += itw[2];
        auxtot += itw[2]*(surface_type_mask.data(gp_lat.idx+1,gp_lon.idx)-ntype);
      };
      if (types[t] == types0[3]) { wtot += itw[3];
        auxtot += itw[3]*(surface_type_mask.data(gp_lat.idx+1,gp_lon.idx+1)-ntype);
      };
      //
      surface_types[t] = types[t];
      surface_types_aux[t] = auxtot / wtot;
      surface_types_weights[t] = wtot;
    }
  }

  else {
    ostringstream os;
    os << "The allowed options for *method are: \"nearest\" and \"linear\",\n"
       << "but you have selected: " << method;
    throw runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void surface_rtpropCallAgendaX(Workspace& ws,
                               Numeric& surface_skin_t,
                               Matrix& surface_los,
                               Tensor4& surface_rmatrix,
                               Matrix& surface_emission,
                               const Vector& f_grid,
                               const Vector& rtp_pos,
                               const Vector& rtp_los,
                               const ArrayOfAgenda& surface_rtprop_agenda_array,
                               const ArrayOfIndex& surface_types,
                               const Vector& surface_types_aux,
                               const Vector& surface_types_weights,
                               const Verbosity&) {
  if (surface_types.nelem() != 1) {
    ostringstream os;
    os << "This method requires that *surface_types* have length 1.\n"
       << "If you are trying to apply a mixture model, you have to use\n"
       << "a setup based on *iySurfaceCallAgendaX* instead.";
    throw runtime_error(os.str());
  }
  if (surface_types[0] < 0)
    throw runtime_error("*surface_type* is not allowed to be negative.");
  if (surface_types[0] >= surface_rtprop_agenda_array.nelem()) {
    ostringstream os;
    os << "*surface_rtprop_agenda_array* has only "
       << surface_rtprop_agenda_array.nelem()
       << " elements,\n while you have selected element " << surface_types[0];
    throw runtime_error(os.str());
  }
  if (abs(surface_types_weights[0]-1)>1e-4)
    throw runtime_error("Sum of *surface_types_weights* deviates from 1.");  

  surface_rtprop_agenda_arrayExecute(ws,
                                     surface_skin_t,
                                     surface_emission,
                                     surface_los,
                                     surface_rmatrix,
                                     surface_types[0],
                                     f_grid,
                                     rtp_pos,
                                     rtp_los,
                                     surface_types_aux[0],
                                     surface_rtprop_agenda_array);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void SurfaceBlackbody(Matrix& surface_los,
                      Tensor4& surface_rmatrix,
                      ArrayOfTensor4& dsurface_rmatrix_dx,
                      Matrix& surface_emission,
                      ArrayOfMatrix& dsurface_emission_dx,
                      const Index& stokes_dim,
                      const Index& atmosphere_dim,
                      const Vector& lat_grid,
                      const Vector& lon_grid,
                      const Vector& f_grid,
                      const Vector& rtp_pos,
                      const Vector& rtp_los,
                      const Tensor3& surface_props_data,
                      const ArrayOfString& surface_props_names,
                      const ArrayOfString& dsurface_names,
                      const Index& jacobian_do,
                      const Verbosity& verbosity) {
  // Check surface_data
  surface_props_check(atmosphere_dim,
                      lat_grid,
                      lon_grid,
                      surface_props_data,
                      surface_props_names);

  // Interplation grid positions and weights
  ArrayOfGridPos gp_lat(1), gp_lon(1);
  Matrix itw;
  rte_pos2gridpos(
      gp_lat[0], gp_lon[0], atmosphere_dim, lat_grid, lon_grid, rtp_pos);
  interp_atmsurface_gp2itw(itw, atmosphere_dim, gp_lat, gp_lon);

  // Skin temperature
  Vector skin_t(1);
  surface_props_interp(skin_t,
                       "Skin temperature",
                       atmosphere_dim,
                       gp_lat,
                       gp_lon,
                       itw,
                       surface_props_data,
                       surface_props_names);

  surfaceBlackbody(surface_los,
                   surface_rmatrix,
                   surface_emission,
                   atmosphere_dim,
                   f_grid,
                   stokes_dim,
                   rtp_pos,
                   rtp_los,
                   skin_t[0],
                   verbosity);

  surface_rmatrix.resize(1, f_grid.nelem(), stokes_dim, stokes_dim);
  surface_rmatrix = 0.0;

  // Jacobian part
  if (jacobian_do) {
    dsurface_check(surface_props_names,
                   dsurface_names,
                   dsurface_rmatrix_dx,
                   dsurface_emission_dx);

    Index irq;

    // Skin temperature
    irq = find_first(dsurface_names, String("Skin temperature"));
    if (irq >= 0) {
      Matrix dbdt(f_grid.nelem(), 1);
      dplanck_dt(dbdt(joker, 0), f_grid, skin_t[0]);
      dsurface_emission_dx[irq] = dbdt;
      dsurface_rmatrix_dx[irq].resize(surface_rmatrix.nbooks(),
                                      surface_rmatrix.npages(),
                                      surface_rmatrix.nrows(),
                                      surface_rmatrix.ncols());
      dsurface_rmatrix_dx[irq] = 0;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void SurfaceDummy(ArrayOfTensor4& dsurface_rmatrix_dx,
                  ArrayOfMatrix& dsurface_emission_dx,
                  const Index& atmosphere_dim,
                  const Vector& lat_grid,
                  const Vector& lon_grid,
                  const Tensor3& surface_props_data,
                  const ArrayOfString& surface_props_names,
                  const ArrayOfString& dsurface_names,
                  const Index& jacobian_do,
                  const Verbosity&) {
  if (surface_props_names.nelem())
    throw runtime_error(
        "When calling this method, *surface_props_names* should be empty.");

  surface_props_check(atmosphere_dim,
                      lat_grid,
                      lon_grid,
                      surface_props_data,
                      surface_props_names);

  if (jacobian_do) {
    dsurface_check(surface_props_names,
                   dsurface_names,
                   dsurface_rmatrix_dx,
                   dsurface_emission_dx);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void SurfaceFastem(Matrix& surface_los,
                   Tensor4& surface_rmatrix,
                   ArrayOfTensor4& dsurface_rmatrix_dx,
                   Matrix& surface_emission,
                   ArrayOfMatrix& dsurface_emission_dx,
                   const Index& stokes_dim,
                   const Index& atmosphere_dim,
                   const Vector& lat_grid,
                   const Vector& lon_grid,
                   const Vector& f_grid,
                   const Vector& rtp_pos,
                   const Vector& rtp_los,
                   const Tensor3& surface_props_data,
                   const ArrayOfString& surface_props_names,
                   const ArrayOfString& dsurface_names,
                   const Index& jacobian_do,
                   const Vector& transmittance,
                   const Index& fastem_version,
                   const Verbosity& verbosity) {
  // Check surface_data
  surface_props_check(atmosphere_dim,
                      lat_grid,
                      lon_grid,
                      surface_props_data,
                      surface_props_names);

  // Interplation grid positions and weights
  ArrayOfGridPos gp_lat(1), gp_lon(1);
  Matrix itw;
  rte_pos2gridpos(
      gp_lat[0], gp_lon[0], atmosphere_dim, lat_grid, lon_grid, rtp_pos);
  interp_atmsurface_gp2itw(itw, atmosphere_dim, gp_lat, gp_lon);

  // Skin temperature
  Vector skin_t(1);
  surface_props_interp(skin_t,
                       "Water skin temperature",
                       atmosphere_dim,
                       gp_lat,
                       gp_lon,
                       itw,
                       surface_props_data,
                       surface_props_names);

  // Wind speed
  Vector wind_speed(1);
  surface_props_interp(wind_speed,
                       "Wind speed",
                       atmosphere_dim,
                       gp_lat,
                       gp_lon,
                       itw,
                       surface_props_data,
                       surface_props_names);

  // Wind direction
  Vector wind_direction(1);
  surface_props_interp(wind_direction,
                       "Wind direction",
                       atmosphere_dim,
                       gp_lat,
                       gp_lon,
                       itw,
                       surface_props_data,
                       surface_props_names);

  // Salinity
  Vector salinity(1);
  surface_props_interp(salinity,
                       "Salinity",
                       atmosphere_dim,
                       gp_lat,
                       gp_lon,
                       itw,
                       surface_props_data,
                       surface_props_names);

  // Call FASTEM
  surfaceFastem(surface_los,
                surface_rmatrix,
                surface_emission,
                atmosphere_dim,
                stokes_dim,
                f_grid,
                rtp_pos,
                rtp_los,
                skin_t[0],
                salinity[0],
                wind_speed[0],
                wind_direction[0],
                transmittance,
                fastem_version,
                verbosity);

  // Jacobian part
  if (jacobian_do) {
    dsurface_check(surface_props_names,
                   dsurface_names,
                   dsurface_rmatrix_dx,
                   dsurface_emission_dx);

    Index irq;

    // Skin temperature
    irq = find_first(dsurface_names, String("Water skin temperature"));
    if (irq >= 0) {
      const Numeric dd = 0.1;
      Matrix surface_los2;
      surfaceFastem(surface_los2,
                    dsurface_rmatrix_dx[irq],
                    dsurface_emission_dx[irq],
                    atmosphere_dim,
                    stokes_dim,
                    f_grid,
                    rtp_pos,
                    rtp_los,
                    skin_t[0] + dd,
                    salinity[0],
                    wind_speed[0],
                    wind_direction[0],
                    transmittance,
                    fastem_version,
                    verbosity);
      //
      dsurface_rmatrix_dx[irq] -= surface_rmatrix;
      dsurface_rmatrix_dx[irq] /= dd;
      dsurface_emission_dx[irq] -= surface_emission;
      dsurface_emission_dx[irq] /= dd;
    }

    // Wind speed
    irq = find_first(dsurface_names, String("Wind speed"));
    if (irq >= 0) {
      const Numeric dd = 0.1;
      Matrix surface_los2;
      surfaceFastem(surface_los2,
                    dsurface_rmatrix_dx[irq],
                    dsurface_emission_dx[irq],
                    atmosphere_dim,
                    stokes_dim,
                    f_grid,
                    rtp_pos,
                    rtp_los,
                    skin_t[0],
                    salinity[0],
                    wind_speed[0] + dd,
                    wind_direction[0],
                    transmittance,
                    fastem_version,
                    verbosity);
      //
      dsurface_rmatrix_dx[irq] -= surface_rmatrix;
      dsurface_rmatrix_dx[irq] /= dd;
      dsurface_emission_dx[irq] -= surface_emission;
      dsurface_emission_dx[irq] /= dd;
    }

    // Wind direction
    irq = find_first(dsurface_names, String("Wind direction"));
    if (irq >= 0) {
      const Numeric dd = 1;
      Matrix surface_los2;
      surfaceFastem(surface_los2,
                    dsurface_rmatrix_dx[irq],
                    dsurface_emission_dx[irq],
                    atmosphere_dim,
                    stokes_dim,
                    f_grid,
                    rtp_pos,
                    rtp_los,
                    skin_t[0],
                    salinity[0],
                    wind_speed[0],
                    wind_direction[0] + dd,
                    transmittance,
                    fastem_version,
                    verbosity);
      //
      dsurface_rmatrix_dx[irq] -= surface_rmatrix;
      dsurface_rmatrix_dx[irq] /= dd;
      dsurface_emission_dx[irq] -= surface_emission;
      dsurface_emission_dx[irq] /= dd;
    }

    // Salinity
    irq = find_first(dsurface_names, String("Salinity"));
    if (irq >= 0) {
      const Numeric dd = 0.0005;
      Matrix surface_los2;
      surfaceFastem(surface_los2,
                    dsurface_rmatrix_dx[irq],
                    dsurface_emission_dx[irq],
                    atmosphere_dim,
                    stokes_dim,
                    f_grid,
                    rtp_pos,
                    rtp_los,
                    skin_t[0],
                    salinity[0] + dd,
                    wind_speed[0],
                    wind_direction[0],
                    transmittance,
                    fastem_version,
                    verbosity);
      //
      dsurface_rmatrix_dx[irq] -= surface_rmatrix;
      dsurface_rmatrix_dx[irq] /= dd;
      dsurface_emission_dx[irq] -= surface_emission;
      dsurface_emission_dx[irq] /= dd;
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void SurfaceFlatScalarReflectivity(Matrix& surface_los,
                                   Tensor4& surface_rmatrix,
                                   ArrayOfTensor4& dsurface_rmatrix_dx,
                                   Matrix& surface_emission,
                                   ArrayOfMatrix& dsurface_emission_dx,
                                   const Index& stokes_dim,
                                   const Index& atmosphere_dim,
                                   const Vector& lat_grid,
                                   const Vector& lon_grid,
                                   const Vector& f_grid,
                                   const Vector& rtp_pos,
                                   const Vector& rtp_los,
                                   const Vector& specular_los,
                                   const Tensor3& surface_props_data,
                                   const ArrayOfString& surface_props_names,
                                   const ArrayOfString& dsurface_names,
                                   const Index& jacobian_do,
                                   const Vector& f_reflectivities,
                                   const Verbosity& verbosity) {
  // Check surface_data
  surface_props_check(atmosphere_dim,
                      lat_grid,
                      lon_grid,
                      surface_props_data,
                      surface_props_names);

  // Check of GINs
  chk_if_increasing("GIN *f_reflectivities*", f_reflectivities);

  // Sizes
  const Index nr = f_reflectivities.nelem();
  const Index nf = f_grid.nelem();

  // Interplation grid positions and weights
  ArrayOfGridPos gp_lat(1), gp_lon(1);
  Matrix itw;
  rte_pos2gridpos(
      gp_lat[0], gp_lon[0], atmosphere_dim, lat_grid, lon_grid, rtp_pos);
  interp_atmsurface_gp2itw(itw, atmosphere_dim, gp_lat, gp_lon);

  // Skin temperature
  Vector surface_skin_t(1);
  surface_props_interp(surface_skin_t,
                       "Skin temperature",
                       atmosphere_dim,
                       gp_lat,
                       gp_lon,
                       itw,
                       surface_props_data,
                       surface_props_names);

  // Reflectivities
  Vector reflectivities(nr);
  //
  for (Index i=0; i<nr; i++) {
    Vector rv(1);
    ostringstream sstr;
    sstr << "Scalar reflectivity " << i;
    surface_props_interp(rv,
                         sstr.str(),
                         atmosphere_dim,
                         gp_lat,
                         gp_lon,
                         itw,
                         surface_props_data,
                         surface_props_names);
    reflectivities[i] = rv[0];
  }

  // Create surface_scalar_reflectivity
  Vector surface_scalar_reflectivity(nf);
  ArrayOfGridPos gp(nf);
  if (nr==1) {
    gp4length1grid(gp);
  } else {
    gridpos(gp, f_reflectivities, f_grid, 1e99);
    jacobian_type_extrapol(gp);
  }
  Matrix itw2(nf, 2);
  interpweights(itw2, gp);
  interp( surface_scalar_reflectivity, itw2, reflectivities, gp);

  // Call non-Jacobian method
  surfaceFlatScalarReflectivity(surface_los,
                                surface_rmatrix,
                                surface_emission,
                                f_grid,
                                stokes_dim,
                                atmosphere_dim,
                                rtp_pos,
                                rtp_los,
                                specular_los,
                                surface_skin_t[0],
                                surface_scalar_reflectivity,
                                verbosity);
  
  // Jacobian part
  if (jacobian_do) {
    dsurface_check(surface_props_names,
                   dsurface_names,
                   dsurface_rmatrix_dx,
                   dsurface_emission_dx);

    Index irq;

    // Skin temperature
    irq = find_first(dsurface_names, String("Skin temperature"));
    if (irq >= 0) {
      const Numeric dd = 0.1;
      Matrix surface_los2;
      surfaceFlatScalarReflectivity(surface_los2,
                                    dsurface_rmatrix_dx[irq],
                                    dsurface_emission_dx[irq],
                                    f_grid,
                                    stokes_dim,
                                    atmosphere_dim,
                                    rtp_pos,
                                    rtp_los,
                                    specular_los,
                                    surface_skin_t[0]+dd,
                                    surface_scalar_reflectivity,
                                    verbosity);
      dsurface_rmatrix_dx[irq] -= surface_rmatrix;
      dsurface_rmatrix_dx[irq] /= dd;
      dsurface_emission_dx[irq] -= surface_emission;
      dsurface_emission_dx[irq] /= dd;
    }

    // Reflectivities
    for (Index i=0; i<nr; i++) {
      ostringstream sstr;
      sstr << "Scalar reflectivity " << i;
      irq = find_first(dsurface_names, String(sstr.str()));
      if (irq >= 0) {
        const Numeric dd = 1e-4;
        Vector rpert = reflectivities;
        rpert[i] += dd;
        interp( surface_scalar_reflectivity, itw2, rpert, gp);
        //
        Matrix surface_los2;
        surfaceFlatScalarReflectivity(surface_los2,
                                      dsurface_rmatrix_dx[irq],
                                      dsurface_emission_dx[irq],
                                      f_grid,
                                      stokes_dim,
                                      atmosphere_dim,
                                      rtp_pos,
                                      rtp_los,
                                      specular_los,
                                      surface_skin_t[0],
                                      surface_scalar_reflectivity,
                                      verbosity);
        dsurface_rmatrix_dx[irq] -= surface_rmatrix;
        dsurface_rmatrix_dx[irq] /= dd;
        dsurface_emission_dx[irq] -= surface_emission;
        dsurface_emission_dx[irq] /= dd;
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void SurfaceTessem(Matrix& surface_los,
                   Tensor4& surface_rmatrix,
                   ArrayOfTensor4& dsurface_rmatrix_dx,
                   Matrix& surface_emission,
                   ArrayOfMatrix& dsurface_emission_dx,
                   const Index& stokes_dim,
                   const Index& atmosphere_dim,
                   const Vector& lat_grid,
                   const Vector& lon_grid,
                   const Vector& f_grid,
                   const Vector& rtp_pos,
                   const Vector& rtp_los,
                   const TessemNN& net_h,
                   const TessemNN& net_v,
                   const Tensor3& surface_props_data,
                   const ArrayOfString& surface_props_names,
                   const ArrayOfString& dsurface_names,
                   const Index& jacobian_do,
                   const Verbosity& verbosity) {
  // Check surface_data
  surface_props_check(atmosphere_dim,
                      lat_grid,
                      lon_grid,
                      surface_props_data,
                      surface_props_names);

  // Interplation grid positions and weights
  ArrayOfGridPos gp_lat(1), gp_lon(1);
  Matrix itw;
  rte_pos2gridpos(
      gp_lat[0], gp_lon[0], atmosphere_dim, lat_grid, lon_grid, rtp_pos);
  interp_atmsurface_gp2itw(itw, atmosphere_dim, gp_lat, gp_lon);

  // Skin temperature
  Vector skin_t(1);
  surface_props_interp(skin_t,
                       "Water skin temperature",
                       atmosphere_dim,
                       gp_lat,
                       gp_lon,
                       itw,
                       surface_props_data,
                       surface_props_names);

  // Wind speed
  Vector wind_speed(1);
  surface_props_interp(wind_speed,
                       "Wind speed",
                       atmosphere_dim,
                       gp_lat,
                       gp_lon,
                       itw,
                       surface_props_data,
                       surface_props_names);

  // Salinity
  Vector salinity(1);
  surface_props_interp(salinity,
                       "Salinity",
                       atmosphere_dim,
                       gp_lat,
                       gp_lon,
                       itw,
                       surface_props_data,
                       surface_props_names);

  // Call TESSEM
  surfaceTessem(surface_los,
                surface_rmatrix,
                surface_emission,
                atmosphere_dim,
                stokes_dim,
                f_grid,
                rtp_pos,
                rtp_los,
                skin_t[0],
                net_h,
                net_v,
                salinity[0],
                wind_speed[0],
                verbosity);

  // Jacobian part
  if (jacobian_do) {
    dsurface_check(surface_props_names,
                   dsurface_names,
                   dsurface_rmatrix_dx,
                   dsurface_emission_dx);

    Index irq;

    // Skin temperature
    irq = find_first(dsurface_names, String("Water skin temperature"));
    if (irq >= 0) {
      const Numeric dd = 0.1;
      Matrix surface_los2;
      surfaceTessem(surface_los2,
                    dsurface_rmatrix_dx[irq],
                    dsurface_emission_dx[irq],
                    atmosphere_dim,
                    stokes_dim,
                    f_grid,
                    rtp_pos,
                    rtp_los,
                    skin_t[0] + dd,
                    net_h,
                    net_v,
                    salinity[0],
                    wind_speed[0],
                    verbosity);
      //
      dsurface_rmatrix_dx[irq] -= surface_rmatrix;
      dsurface_rmatrix_dx[irq] /= dd;
      dsurface_emission_dx[irq] -= surface_emission;
      dsurface_emission_dx[irq] /= dd;
    }

    // Wind speed
    irq = find_first(dsurface_names, String("Wind speed"));
    if (irq >= 0) {
      const Numeric dd = 0.1;
      Matrix surface_los2;
      surfaceTessem(surface_los2,
                    dsurface_rmatrix_dx[irq],
                    dsurface_emission_dx[irq],
                    atmosphere_dim,
                    stokes_dim,
                    f_grid,
                    rtp_pos,
                    rtp_los,
                    skin_t[0],
                    net_h,
                    net_v,
                    salinity[0],
                    wind_speed[0] + dd,
                    verbosity);
      //
      dsurface_rmatrix_dx[irq] -= surface_rmatrix;
      dsurface_rmatrix_dx[irq] /= dd;
      dsurface_emission_dx[irq] -= surface_emission;
      dsurface_emission_dx[irq] /= dd;
    }

    // Salinity
    irq = find_first(dsurface_names, String("Salinity"));
    if (irq >= 0) {
      const Numeric dd = 0.0005;
      Matrix surface_los2;
      surfaceTessem(surface_los2,
                    dsurface_rmatrix_dx[irq],
                    dsurface_emission_dx[irq],
                    atmosphere_dim,
                    stokes_dim,
                    f_grid,
                    rtp_pos,
                    rtp_los,
                    skin_t[0],
                    net_h,
                    net_v,
                    salinity[0] + dd,
                    wind_speed[0],
                    verbosity);
      //
      dsurface_rmatrix_dx[irq] -= surface_rmatrix;
      dsurface_rmatrix_dx[irq] /= dd;
      dsurface_emission_dx[irq] -= surface_emission;
      dsurface_emission_dx[irq] /= dd;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void transmittanceFromIy_aux(Vector& transmittance,
                             const ArrayOfString& iy_aux_vars,
                             const ArrayOfMatrix& iy_aux,
                             const Verbosity&) {
  Index ihit = -1;

  for (Index i = 0; i < iy_aux_vars.nelem(); i++) {
    if (iy_aux_vars[i] == "Optical depth") {
      ihit = i;
      break;
    }
  }

  if (ihit < 0)
    throw runtime_error("No element in *iy_aux* holds optical depths.");

  const Index n = iy_aux[ihit].nrows();

  transmittance.resize(n);

  for (Index i = 0; i < n; i++) {
    transmittance[i] = exp(-iy_aux[ihit](i, 0));
  }
}
