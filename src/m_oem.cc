/* Copyright (C) 2015
   Patrick Eriksson <patrick.eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the oem::GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   oem::GNU General Public License for more details.

   You should have received a copy of the oem::GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*===========================================================================
  === File description
  ===========================================================================*/

/*!
  \file   m_oem.cc
  \author Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2015-09-08

  \brief  Workspace functions related to making OEM inversions.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include "array.h"
#include "arts.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "jacobian.h"
#include "math_funcs.h"
#include "physics_funcs.h"
#include "rte.h"
#include "special_interp.h"
#include "surface.h"

#ifdef OEM_SUPPORT
#include "oem.h"
#endif

/* Workspace method: Doxygen documentation will be auto-generated */
void particle_bulkprop_fieldClip(Tensor4& particle_bulkprop_field,
                                 const ArrayOfString& particle_bulkprop_names,
                                 const String& bulkprop_name,
                                 const Numeric& limit_low,
                                 const Numeric& limit_high,
                                 const Verbosity&) {
  Index iq = -1;
  if (bulkprop_name == "ALL") {
  }

  else {
    for (Index i = 0; i < particle_bulkprop_names.nelem(); i++) {
      if (particle_bulkprop_names[i] == bulkprop_name) {
        iq = i;
        break;
      }
    }
    ARTS_USER_ERROR_IF (iq < 0,
      "Could not find ", bulkprop_name,
      " in particle_bulkprop_names.\n")
  }

  Tensor4Clip(particle_bulkprop_field, iq, limit_low, limit_high);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void vmr_fieldClip(Tensor4& vmr_field,
                   const ArrayOfArrayOfSpeciesTag& abs_species,
                   const String& species,
                   const Numeric& limit_low,
                   const Numeric& limit_high,
                   const Verbosity&) {
  Index iq = -1;
  if (species == "ALL") {
  }

  else {
    for (Index i = 0; i < abs_species.nelem(); i++) {
      if (abs_species[i].Species() == SpeciesTag(species).Spec()) {
        iq = i;
        break;
      }
    }
    ARTS_USER_ERROR_IF (iq < 0,
      "Could not find ", species, " in abs_species.\n")
  }

  Tensor4Clip(vmr_field, iq, limit_low, limit_high);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void xClip(Vector& x,
           const ArrayOfRetrievalQuantity& jacobian_quantities,
           const Index& ijq,
           const Numeric& limit_low,
           const Numeric& limit_high,
           const Verbosity&) {
  // Sizes
  const Index nq = jacobian_quantities.nelem();

  ARTS_USER_ERROR_IF (ijq < -1, "Argument *ijq* must be >= -1.");
  ARTS_USER_ERROR_IF (ijq >= nq,
      "Argument *ijq* is too high.\n"
      "You have selected index: ", ijq, "\n"
      "but the number of quantities is only: ", nq, "\n"
      "(Note that zero-based indexing is used)\n")

  // Jacobian indices
  ArrayOfArrayOfIndex ji;
  {
    bool any_affine;
    jac_ranges_indices(ji, any_affine, jacobian_quantities);
  }

  Index ifirst = 0, ilast = x.nelem() - 1;
  if (ijq > -1) {
    ifirst = ji[ijq][0];
    ilast = ji[ijq][1];
  }

  if (!std::isinf(limit_low)) {
    for (Index i = ifirst; i <= ilast; i++) {
      if (x[i] < limit_low) x[i] = limit_low;
    }
  }
  if (!std::isinf(limit_high)) {
    for (Index i = ifirst; i <= ilast; i++) {
      if (x[i] > limit_high) x[i] = limit_high;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void xaStandard(Workspace& ws,
                Vector& xa,
                const ArrayOfRetrievalQuantity& jacobian_quantities,
                const Index& atmfields_checked,
                const Index& atmgeom_checked,
                const Index& atmosphere_dim,
                const Vector& p_grid,
                const Vector& lat_grid,
                const Vector& lon_grid,
                const Tensor3& t_field,
                const Tensor4& vmr_field,
                const ArrayOfArrayOfSpeciesTag& abs_species,
                const Index& cloudbox_on,
                const Index& cloudbox_checked,
                const Tensor4& particle_bulkprop_field,
                const ArrayOfString& particle_bulkprop_names,
                const Tensor3& wind_u_field,
                const Tensor3& wind_v_field,
                const Tensor3& wind_w_field,
                const Tensor3& mag_u_field,
                const Tensor3& mag_v_field,
                const Tensor3& mag_w_field,
                const Tensor3& surface_props_data,
                const ArrayOfString& surface_props_names,
                const Agenda& water_p_eq_agenda,
                const Verbosity&) {
  // Basics
  //
  ARTS_USER_ERROR_IF (atmfields_checked != 1,
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");
  ARTS_USER_ERROR_IF (atmgeom_checked != 1,
        "The atmospheric geometry must be flagged to have "
        "passed a consistency check (atmgeom_checked=1).");
  ARTS_USER_ERROR_IF (cloudbox_checked != 1,
        "The cloudbox must be flagged to have "
        "passed a consistency check (cloudbox_checked=1).");

  // Jacobian indices
  ArrayOfArrayOfIndex ji;
  {
    bool any_affine;
    jac_ranges_indices(ji, any_affine, jacobian_quantities, true);
  }

  // Sizes
  const Index nq = jacobian_quantities.nelem();
  //
  xa.resize(ji[nq - 1][1] + 1);

  // Loop retrieval quantities and fill *xa*
  for (Index q = 0; q < jacobian_quantities.nelem(); q++) {
    // Index range of this retrieval quantity
    const Index np = ji[q][1] - ji[q][0] + 1;
    Range ind(ji[q][0], np);

    // Atmospheric temperatures
    if (jacobian_quantities[q] == Jacobian::Atm::Temperature) {
      // Here we need to interpolate *t_field*
      ArrayOfGridPos gp_p, gp_lat, gp_lon;
      get_gp_atmgrids_to_rq(gp_p,
                            gp_lat,
                            gp_lon,
                            jacobian_quantities[q],
                            atmosphere_dim,
                            p_grid,
                            lat_grid,
                            lon_grid);
      Tensor3 t_x;
      regrid_atmfield_by_gp(t_x, atmosphere_dim, t_field, gp_p, gp_lat, gp_lon);
      flat(xa[ind], t_x);
    }

    // Abs species
    else if (jacobian_quantities[q].Target().isSpeciesVMR()) {
      // Index position of species
      ArrayOfSpeciesTag atag(jacobian_quantities[q].Subtag());
      const Index isp = chk_contains("abs_species", abs_species, atag);

      if (jacobian_quantities[q].Mode() == "rel") {
        // This one is simple, just a vector of ones
        xa[ind] = 1;
      } else {
        // For all remaining options we need to interpolate *vmr_field*
        ArrayOfGridPos gp_p, gp_lat, gp_lon;
        get_gp_atmgrids_to_rq(gp_p,
                              gp_lat,
                              gp_lon,
                              jacobian_quantities[q],
                              atmosphere_dim,
                              p_grid,
                              lat_grid,
                              lon_grid);
        Tensor3 vmr_x;
        regrid_atmfield_by_gp(vmr_x,
                              atmosphere_dim,
                              vmr_field(isp, joker, joker, joker),
                              gp_p,
                              gp_lat,
                              gp_lon);

        if (jacobian_quantities[q].Mode() == "vmr") {
          flat(xa[ind], vmr_x);
        } else if (jacobian_quantities[q].Mode() == "nd") {
          // Here we need to also interpolate *t_field*
          Tensor3 t_x;
          regrid_atmfield_by_gp(
              t_x, atmosphere_dim, t_field, gp_p, gp_lat, gp_lon);
          // Calculate number density for species (vmr*nd_tot)
          Index i = 0;
          for (Index i3 = 0; i3 < vmr_x.ncols(); i3++) {
            for (Index i2 = 0; i2 < vmr_x.nrows(); i2++) {
              for (Index i1 = 0; i1 < vmr_x.npages(); i1++) {
                xa[ji[q][0] + i] =
                    vmr_x(i1, i2, i3) *
                    number_density(jacobian_quantities[q].Grids()[0][i1],
                                   t_x(i1, i2, i3));
                i += 1;
              }
            }
          }
        } else if (jacobian_quantities[q].Mode() == "rh") {
          // Here we need to also interpolate *t_field*
          Tensor3 t_x;
          regrid_atmfield_by_gp(
              t_x, atmosphere_dim, t_field, gp_p, gp_lat, gp_lon);
          Tensor3 water_p_eq;
          water_p_eq_agendaExecute(ws, water_p_eq, t_x, water_p_eq_agenda);
          // Calculate relative humidity (vmr*p/p_sat)
          Index i = 0;
          for (Index i3 = 0; i3 < vmr_x.ncols(); i3++) {
            for (Index i2 = 0; i2 < vmr_x.nrows(); i2++) {
              for (Index i1 = 0; i1 < vmr_x.npages(); i1++) {
                xa[ji[q][0] + i] = vmr_x(i1, i2, i3) *
                                   jacobian_quantities[q].Grids()[0][i1] /
                                   water_p_eq(i1, i2, i3);
                i += 1;
              }
            }
          }
        } else if (jacobian_quantities[q].Mode() == "q") {
          // Calculate specific humidity q, from mixing ratio r and
          // vapour pressure e, as
          // q = r(1+r); r = 0.622e/(p-e); e = vmr*p;
          Index i = 0;
          for (Index i3 = 0; i3 < vmr_x.ncols(); i3++) {
            for (Index i2 = 0; i2 < vmr_x.nrows(); i2++) {
              for (Index i1 = 0; i1 < vmr_x.npages(); i1++) {
                const Numeric e =
                    vmr_x(i1, i2, i3) * jacobian_quantities[q].Grids()[0][i1];
                const Numeric r =
                    0.622 * e / (jacobian_quantities[q].Grids()[0][i1] - e);
                xa[ji[q][0] + i] = r / (1 + r);
                i += 1;
              }
            }
          }
        } else {
          ARTS_ASSERT(0);
        }
      }
    }

    // Scattering species
    else if (jacobian_quantities[q] == Jacobian::Special::ScatteringString) {
      if (cloudbox_on) {
        ARTS_USER_ERROR_IF (particle_bulkprop_field.empty(),
              "One jacobian quantity belongs to the "
              "scattering species category, but *particle_bulkprop_field* "
              "is empty.");
        ARTS_USER_ERROR_IF (particle_bulkprop_field.nbooks() !=
                            particle_bulkprop_names.nelem(),
              "Mismatch in size between "
              "*particle_bulkprop_field* and *particle_bulkprop_names*.");

        const Index isp = find_first(particle_bulkprop_names,
                                     jacobian_quantities[q].SubSubtag());
        ARTS_USER_ERROR_IF (isp < 0,
            "Jacobian quantity with index ", q, " covers a "
            "scattering species, and the field quantity is set to \"",
            jacobian_quantities[q].SubSubtag(), "\", but this quantity "
            "could not found in *particle_bulkprop_names*.")

        ArrayOfGridPos gp_p, gp_lat, gp_lon;
        get_gp_atmgrids_to_rq(gp_p,
                              gp_lat,
                              gp_lon,
                              jacobian_quantities[q],
                              atmosphere_dim,
                              p_grid,
                              lat_grid,
                              lon_grid);
        Tensor3 pbp_x;
        regrid_atmfield_by_gp(pbp_x,
                              atmosphere_dim,
                              particle_bulkprop_field(isp, joker, joker, joker),
                              gp_p,
                              gp_lat,
                              gp_lon);
        flat(xa[ind], pbp_x);
      } else {
        xa[ind] = 0;
      }
    }

    // Wind
    else if (jacobian_quantities[q].Target().isWind()) {
      ConstTensor3View source_field(wind_u_field);
      if (jacobian_quantities[q] == Jacobian::Atm::WindV) {
        source_field = wind_v_field;
      } else if (jacobian_quantities[q] == Jacobian::Atm::WindW) {
        source_field = wind_w_field;
      }

      // Determine grid positions for interpolation from retrieval grids back
      // to atmospheric grids
      ArrayOfGridPos gp_p, gp_lat, gp_lon;
      get_gp_atmgrids_to_rq(gp_p,
                            gp_lat,
                            gp_lon,
                            jacobian_quantities[q],
                            atmosphere_dim,
                            p_grid,
                            lat_grid,
                            lon_grid);

      Tensor3 wind_x;
      regrid_atmfield_by_gp(
          wind_x, atmosphere_dim, source_field, gp_p, gp_lat, gp_lon);
      flat(xa[ind], wind_x);
    }

    // Magnetism
    else if (jacobian_quantities[q].Target().isMagnetic()) {
      if (jacobian_quantities[q] == Jacobian::Atm::MagneticMagnitude) {
        // Determine grid positions for interpolation from retrieval grids back
        // to atmospheric grids
        ArrayOfGridPos gp_p, gp_lat, gp_lon;
        get_gp_atmgrids_to_rq(gp_p,
                              gp_lat,
                              gp_lon,
                              jacobian_quantities[q],
                              atmosphere_dim,
                              p_grid,
                              lat_grid,
                              lon_grid);

        //all three component's hyoptenuse is the strength
        Tensor3 mag_u, mag_v, mag_w;
        regrid_atmfield_by_gp(
            mag_u, atmosphere_dim, mag_u_field, gp_p, gp_lat, gp_lon);
        regrid_atmfield_by_gp(
            mag_v, atmosphere_dim, mag_v_field, gp_p, gp_lat, gp_lon);
        regrid_atmfield_by_gp(
            mag_w, atmosphere_dim, mag_w_field, gp_p, gp_lat, gp_lon);

        Tensor3 mag_x(gp_p.nelem(), gp_lat.nelem(), gp_lon.nelem());
        for (Index i = 0; i < gp_p.nelem(); i++)
          for (Index j = 0; j < gp_lat.nelem(); j++)
            for (Index k = 0; k < gp_lon.nelem(); k++)
              mag_x(i, j, k) = std::hypot(mag_u(i, j, k), mag_v(i, j, k), mag_w(i, j, k));
        flat(xa[ind], mag_x);
      } else {
        ConstTensor3View source_field(mag_u_field);
        if (jacobian_quantities[q] == Jacobian::Atm::MagneticV) {
          source_field = mag_v_field;
        } else if (jacobian_quantities[q] == Jacobian::Atm::MagneticW) {
          source_field = mag_w_field;
        }

        // Determine grid positions for interpolation from retrieval grids back
        // to atmospheric grids
        ArrayOfGridPos gp_p, gp_lat, gp_lon;
        get_gp_atmgrids_to_rq(gp_p,
                              gp_lat,
                              gp_lon,
                              jacobian_quantities[q],
                              atmosphere_dim,
                              p_grid,
                              lat_grid,
                              lon_grid);

        Tensor3 mag_x;
        regrid_atmfield_by_gp(
            mag_x, atmosphere_dim, source_field, gp_p, gp_lat, gp_lon);
        flat(xa[ind], mag_x);
      }
    }

    // Surface
    else if (jacobian_quantities[q] == Jacobian::Special::SurfaceString) {
      surface_props_check(atmosphere_dim,
                          lat_grid,
                          lon_grid,
                          surface_props_data,
                          surface_props_names);
      ARTS_USER_ERROR_IF (surface_props_data.empty(),
            "One jacobian quantity belongs to the "
            "surface category, but *surface_props_data* is empty.");

      const Index isu =
          find_first(surface_props_names, jacobian_quantities[q].Subtag());
      ARTS_USER_ERROR_IF (isu < 0,
          "Jacobian quantity with index ", q, " covers a "
          "surface property, and the field Subtag is set to \"",
          jacobian_quantities[q].Subtag(), "\", but this quantity "
          "could not found in *surface_props_names*.")

      ArrayOfGridPos gp_lat, gp_lon;
      get_gp_atmsurf_to_rq(gp_lat,
                           gp_lon,
                           jacobian_quantities[q],
                           atmosphere_dim,
                           lat_grid,
                           lon_grid);
      Matrix surf_x;
      regrid_atmsurf_by_gp_oem(surf_x,
                               atmosphere_dim,
                               surface_props_data(isu, joker, joker),
                               gp_lat,
                               gp_lon);
      flat(xa[ind], surf_x);
    }

    // All variables having zero as a priori
    // ----------------------------------------------------------------------------
    else if (jacobian_quantities[q].Target().isPointing() ||
             jacobian_quantities[q].Target().isFrequency() ||
             jacobian_quantities[q] == Jacobian::Sensor::Polyfit ||
             jacobian_quantities[q] == Jacobian::Sensor::Sinefit) {
      xa[ind] = 0;
    }

    else {
      ARTS_USER_ERROR (
        "Found a retrieval quantity that is not yet handled by\n"
        "internal retrievals: ", jacobian_quantities[q], '\n')
    }
  }

  // Apply transformations
  transform_x(xa, jacobian_quantities);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void x2artsAtmAndSurf(Workspace& ws,
                      Tensor4& vmr_field,
                      Tensor3& t_field,
                      Tensor4& particle_bulkprop_field,
                      Tensor3& wind_u_field,
                      Tensor3& wind_v_field,
                      Tensor3& wind_w_field,
                      Tensor3& mag_u_field,
                      Tensor3& mag_v_field,
                      Tensor3& mag_w_field,
                      Tensor3& surface_props_data,
                      const ArrayOfRetrievalQuantity& jacobian_quantities,
                      const Vector& x,
                      const Index& atmfields_checked,
                      const Index& atmgeom_checked,
                      const Index& atmosphere_dim,
                      const Vector& p_grid,
                      const Vector& lat_grid,
                      const Vector& lon_grid,
                      const ArrayOfArrayOfSpeciesTag& abs_species,
                      const Index& cloudbox_on,
                      const Index& cloudbox_checked,
                      const ArrayOfString& particle_bulkprop_names,
                      const ArrayOfString& surface_props_names,
                      const Agenda& water_p_eq_agenda,
                      const Verbosity&) {
  // Basics
  //
  ARTS_USER_ERROR_IF (atmfields_checked != 1,
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");
  ARTS_USER_ERROR_IF (atmgeom_checked != 1,
        "The atmospheric geometry must be flagged to have "
        "passed a consistency check (atmgeom_checked=1).");
  ARTS_USER_ERROR_IF (cloudbox_checked != 1,
        "The cloudbox must be flagged to have "
        "passed a consistency check (cloudbox_checked=1).");

  // Revert transformation
  Vector x_t(x);
  transform_x_back(x_t, jacobian_quantities);

  // Main sizes
  const Index nq = jacobian_quantities.nelem();

  // Jacobian indices
  ArrayOfArrayOfIndex ji;
  {
    bool any_affine;
    jac_ranges_indices(ji, any_affine, jacobian_quantities, true);
  }

  // Check input
  ARTS_USER_ERROR_IF (x_t.nelem() != ji[nq - 1][1] + 1,
        "Length of *x* does not match length implied by "
        "*jacobian_quantities*.");

  // Note that when this method is called, vmr_field and other output variables
  // have original values, i.e. matching the a priori state.

  // Loop retrieval quantities
  for (Index q = 0; q < nq; q++) {
    // Index range of this retrieval quantity
    const Index np = ji[q][1] - ji[q][0] + 1;
    Range ind(ji[q][0], np);

    // Atmospheric temperatures
    // ----------------------------------------------------------------------------
    if (jacobian_quantities[q] == Jacobian::Atm::Temperature) {
      // Determine grid positions for interpolation from retrieval grids back
      // to atmospheric grids
      ArrayOfGridPos gp_p, gp_lat, gp_lon;
      Index n_p, n_lat, n_lon;
      get_gp_rq_to_atmgrids(gp_p,
                            gp_lat,
                            gp_lon,
                            n_p,
                            n_lat,
                            n_lon,
                            jacobian_quantities[q].Grids(),
                            atmosphere_dim,
                            p_grid,
                            lat_grid,
                            lon_grid);

      // Map values in x back to t_field
      Tensor3 t_x(n_p, n_lat, n_lon);
      reshape(t_x, x_t[ind]);
      regrid_atmfield_by_gp_oem(
          t_field, atmosphere_dim, t_x, gp_p, gp_lat, gp_lon);
    }

    // Abs species
    // ----------------------------------------------------------------------------
    else if (jacobian_quantities[q].Target().isSpeciesVMR()) {
      // Index position of species
      ArrayOfSpeciesTag atag(jacobian_quantities[q].Subtag());
      const Index isp = chk_contains("abs_species", abs_species, atag);

      // Map part of x to a full atmospheric field
      Tensor3 x_field(vmr_field.npages(), vmr_field.nrows(), vmr_field.ncols());
      {
        ArrayOfGridPos gp_p, gp_lat, gp_lon;
        Index n_p, n_lat, n_lon;
        get_gp_rq_to_atmgrids(gp_p,
                              gp_lat,
                              gp_lon,
                              n_p,
                              n_lat,
                              n_lon,
                              jacobian_quantities[q].Grids(),
                              atmosphere_dim,
                              p_grid,
                              lat_grid,
                              lon_grid);
        //
        Tensor3 t3_x(n_p, n_lat, n_lon);
        reshape(t3_x, x_t[ind]);
        regrid_atmfield_by_gp_oem(
            x_field, atmosphere_dim, t3_x, gp_p, gp_lat, gp_lon);
      }
      //
      if (jacobian_quantities[q].Mode() == "rel") {
        // vmr = vmr0 * x
        vmr_field(isp, joker, joker, joker) *= x_field;
      } else if (jacobian_quantities[q].Mode() == "vmr") {
        // vmr = x
        vmr_field(isp, joker, joker, joker) = x_field;
      } else if (jacobian_quantities[q].Mode() == "nd") {
        // vmr = nd / nd_tot
        for (Index i3 = 0; i3 < vmr_field.ncols(); i3++) {
          for (Index i2 = 0; i2 < vmr_field.nrows(); i2++) {
            for (Index i1 = 0; i1 < vmr_field.npages(); i1++) {
              vmr_field(isp, i1, i2, i3) =
                  x_field(i1, i2, i3) /
                  number_density(p_grid[i1], t_field(i1, i2, i3));
            }
          }
        }
      } else if (jacobian_quantities[q].Mode() == "rh") {
        // vmr = x * p_sat / p
        Tensor3 water_p_eq;
        water_p_eq_agendaExecute(ws, water_p_eq, t_field, water_p_eq_agenda);
        for (Index i3 = 0; i3 < vmr_field.ncols(); i3++) {
          for (Index i2 = 0; i2 < vmr_field.nrows(); i2++) {
            for (Index i1 = 0; i1 < vmr_field.npages(); i1++) {
              vmr_field(isp, i1, i2, i3) =
                  x_field(i1, i2, i3) * water_p_eq(i1, i2, i3) / p_grid[i1];
            }
          }
        }
      } else if (jacobian_quantities[q].Mode() == "q") {
        // We have that specific humidity q, mixing ratio r and
        // vapour pressure e, are related as
        // q = r(1+r); r = 0.622e/(p-e); e = vmr*p;
        // That is: vmr=e/p; e = rp/(0.622+r); r = q/(1-q)
        for (Index i3 = 0; i3 < vmr_field.ncols(); i3++) {
          for (Index i2 = 0; i2 < vmr_field.nrows(); i2++) {
            for (Index i1 = 0; i1 < vmr_field.npages(); i1++) {
              const Numeric r = x_field(i1, i2, i3) / (1 - x_field(i1, i2, i3));
              const Numeric e = r * p_grid[i1] / (0.622 + r);
              vmr_field(isp, i1, i2, i3) = e / p_grid[i1];
            }
          }
        }
      } else {
        ARTS_ASSERT(0);
      }
    }

    // Scattering species
    // ----------------------------------------------------------------------------
    else if (jacobian_quantities[q] == Jacobian::Special::ScatteringString) {
      // If no cloudbox, we assume that there is nothing to do
      if (cloudbox_on) {
        ARTS_USER_ERROR_IF (particle_bulkprop_field.empty(),
              "One jacobian quantity belongs to the "
              "scattering species category, but *particle_bulkprop_field* "
              "is empty.");
        ARTS_USER_ERROR_IF (particle_bulkprop_field.nbooks() !=
                            particle_bulkprop_names.nelem(),
              "Mismatch in size between "
              "*particle_bulkprop_field* and *particle_bulkprop_field*.");

        const Index isp = find_first(particle_bulkprop_names,
                                     jacobian_quantities[q].SubSubtag());
        ARTS_USER_ERROR_IF (isp < 0,
            "Jacobian quantity with index ", q, " covers a "
            "scattering species, and the field quantity is set to \"",
            jacobian_quantities[q].SubSubtag(), "\", but this quantity "
            "could not found in *particle_bulkprop_names*.")

        // Determine grid positions for interpolation from retrieval grids back
        // to atmospheric grids
        ArrayOfGridPos gp_p, gp_lat, gp_lon;
        Index n_p, n_lat, n_lon;
        get_gp_rq_to_atmgrids(gp_p,
                              gp_lat,
                              gp_lon,
                              n_p,
                              n_lat,
                              n_lon,
                              jacobian_quantities[q].Grids(),
                              atmosphere_dim,
                              p_grid,
                              lat_grid,
                              lon_grid);
        // Map x to particle_bulkprop_field
        Tensor3 pbfield_x(n_p, n_lat, n_lon);
        reshape(pbfield_x, x_t[ind]);
        Tensor3 pbfield;
        regrid_atmfield_by_gp_oem(
            pbfield, atmosphere_dim, pbfield_x, gp_p, gp_lat, gp_lon);
        particle_bulkprop_field(isp, joker, joker, joker) = pbfield;
      }
    }

    // Wind
    // ----------------------------------------------------------------------------
    else if (jacobian_quantities[q].Target().isWind()) {
      // Determine grid positions for interpolation from retrieval grids back
      // to atmospheric grids
      ArrayOfGridPos gp_p, gp_lat, gp_lon;
      Index n_p, n_lat, n_lon;
      get_gp_rq_to_atmgrids(gp_p,
                            gp_lat,
                            gp_lon,
                            n_p,
                            n_lat,
                            n_lon,
                            jacobian_quantities[q].Grids(),
                            atmosphere_dim,
                            p_grid,
                            lat_grid,
                            lon_grid);

      // TODO Could be done without copying.
      Tensor3 wind_x(n_p, n_lat, n_lon);
      reshape(wind_x, x_t[ind]);

      Tensor3View target_field(wind_u_field);

      Tensor3 wind_field(
          target_field.npages(), target_field.nrows(), target_field.ncols());
      regrid_atmfield_by_gp_oem(
          wind_field, atmosphere_dim, wind_x, gp_p, gp_lat, gp_lon);

      if (jacobian_quantities[q] == Jacobian::Atm::WindU) {
        wind_u_field = wind_field;
      } else if (jacobian_quantities[q] == Jacobian::Atm::WindV) {
        wind_v_field = wind_field;
      } else if (jacobian_quantities[q] == Jacobian::Atm::WindW) {
        wind_w_field = wind_field;
      }
    }

    // Magnetism
    // ----------------------------------------------------------------------------
    else if (jacobian_quantities[q].Target().isMagnetic()) {
      // Determine grid positions for interpolation from retrieval grids back
      // to atmospheric grids
      ArrayOfGridPos gp_p, gp_lat, gp_lon;
      Index n_p, n_lat, n_lon;
      get_gp_rq_to_atmgrids(gp_p,
                            gp_lat,
                            gp_lon,
                            n_p,
                            n_lat,
                            n_lon,
                            jacobian_quantities[q].Grids(),
                            atmosphere_dim,
                            p_grid,
                            lat_grid,
                            lon_grid);

      // TODO Could be done without copying.
      Tensor3 mag_x(n_p, n_lat, n_lon);
      reshape(mag_x, x_t[ind]);

      Tensor3View target_field(mag_u_field);

      Tensor3 mag_field(
          target_field.npages(), target_field.nrows(), target_field.ncols());
      regrid_atmfield_by_gp_oem(
          mag_field, atmosphere_dim, mag_x, gp_p, gp_lat, gp_lon);
      if (jacobian_quantities[q] == Jacobian::Atm::MagneticU) {
        mag_u_field = mag_field;
      } else if (jacobian_quantities[q] == Jacobian::Atm::MagneticV) {
        mag_v_field = mag_field;
      } else if (jacobian_quantities[q] == Jacobian::Atm::MagneticW) {
        mag_w_field = mag_field;
      } else if (jacobian_quantities[q] == Jacobian::Atm::MagneticMagnitude) {
        for (Index i = 0; i < n_p; i++) {
          for (Index j = 0; j < n_lat; j++) {
            for (Index k = 0; k < n_lon; k++) {
              const Numeric scale = mag_x(i, j, k) /
                std::hypot(mag_u_field(i, j, k), mag_v_field(i, j, k), mag_w_field(i, j, k));
              mag_u_field(i, j, k) *= scale;
              mag_v_field(i, j, k) *= scale;
              mag_w_field(i, j, k) *= scale;
            }
          }
        }
      } else {
        ARTS_USER_ERROR ( "Unsupported magnetism type");
      }
    }

    // Surface
    // ----------------------------------------------------------------------------
    else if (jacobian_quantities[q] == Jacobian::Special::SurfaceString) {
      surface_props_check(atmosphere_dim,
                          lat_grid,
                          lon_grid,
                          surface_props_data,
                          surface_props_names);
      ARTS_USER_ERROR_IF (surface_props_data.empty(),
            "One jacobian quantity belongs to the "
            "surface category, but *surface_props_data* is empty.");

      const Index isu =
          find_first(surface_props_names, jacobian_quantities[q].Subtag());
      ARTS_USER_ERROR_IF (isu < 0,
          "Jacobian quantity with index ", q, " covers a "
          "surface property, and the field Subtag is set to \"",
          jacobian_quantities[q].Subtag(), "\", but this quantity "
          "could not found in *surface_props_names*.")

      // Determine grid positions for interpolation from retrieval grids back
      // to atmospheric grids
      ArrayOfGridPos gp_lat, gp_lon;
      Index n_lat, n_lon;
      get_gp_rq_to_atmgrids(gp_lat,
                            gp_lon,
                            n_lat,
                            n_lon,
                            jacobian_quantities[q].Grids(),
                            atmosphere_dim,
                            lat_grid,
                            lon_grid);
      // Map values in x back to surface_props_data
      Matrix surf_x(n_lat, n_lon);
      reshape(surf_x, x_t[ind]);
      Matrix surf;
      regrid_atmsurf_by_gp_oem(surf, atmosphere_dim, surf_x, gp_lat, gp_lon);
      surface_props_data(isu, joker, joker) = surf;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void x2artsSensor(Workspace& ws,
                  Matrix& sensor_los,
                  Vector& f_backend,
                  Vector& y_baseline,
                  Sparse& sensor_response,
                  Vector& sensor_response_f,
                  ArrayOfIndex& sensor_response_pol,
                  Matrix& sensor_response_dlos,
                  Vector& sensor_response_f_grid,
                  ArrayOfIndex& sensor_response_pol_grid,
                  Matrix& sensor_response_dlos_grid,
                  Matrix& mblock_dlos_grid,
                  const ArrayOfRetrievalQuantity& jacobian_quantities,
                  const Vector& x,
                  const Agenda& sensor_response_agenda,
                  const Index& sensor_checked,
                  const ArrayOfTime& sensor_time,
                  const Verbosity&) {
  // Basics
  //
  ARTS_USER_ERROR_IF (sensor_checked != 1,
        "The sensor response must be flagged to have "
        "passed a consistency check (sensor_checked=1).");

  // Revert transformation
  Vector x_t(x);
  transform_x_back(x_t, jacobian_quantities);

  // Main sizes
  const Index nq = jacobian_quantities.nelem();

  // Jacobian indices
  ArrayOfArrayOfIndex ji;
  {
    bool any_affine;
    jac_ranges_indices(ji, any_affine, jacobian_quantities, true);
  }

  // Check input
  ARTS_USER_ERROR_IF (x_t.nelem() != ji[nq - 1][1] + 1,
        "Length of *x* does not match length implied by "
        "*jacobian_quantities*.");

  // Flag indicating that y_baseline is not set
  bool yb_set = false;

  // Shall sensor responses be calculed?
  bool do_sensor = false;

  // Loop retrieval quantities
  for (Index q = 0; q < nq; q++) {
    // Index range of this retrieval quantity
    const Index np = ji[q][1] - ji[q][0] + 1;

    // Pointing off-set
    // ----------------------------------------------------------------------------
    if (jacobian_quantities[q].Target().isPointing()) {
      // Handle pointing "jitter" seperately
      if (jacobian_quantities[q].Grids()[0][0] == -1) {
        ARTS_USER_ERROR_IF (sensor_los.nrows() != np,
              "Mismatch between pointing jacobian and *sensor_los* found.");
        // Simply add retrieved off-set(s) to za column of *sensor_los*
        for (Index i = 0; i < np; i++) {
          sensor_los(i, 0) += x_t[ji[q][0] + i];
        }
      }
      // Polynomial representation
      else {
        ARTS_USER_ERROR_IF (sensor_los.nrows() != sensor_time.nelem(),
              "Sizes of *sensor_los* and *sensor_time* do not match.");
        Vector w;
        for (Index c = 0; c < np; c++) {
          polynomial_basis_func(w, time_vector(sensor_time), c);
          for (Index i = 0; i < w.nelem(); i++) {
            sensor_los(i, 0) += w[i] * x_t[ji[q][0] + c];
          }
        }
      }
    }

    // Frequency shift or stretch
    // ----------------------------------------------------------------------------
    else if (jacobian_quantities[q].Target().isFrequency()) {
      if (jacobian_quantities[q] == Jacobian::Sensor::FrequencyShift) {
        ARTS_ASSERT(np == 1);
        if (x_t[ji[q][0]] != 0) {
          do_sensor = true;
          f_backend += x_t[ji[q][0]];
        }
      } else if (jacobian_quantities[q] == Jacobian::Sensor::FrequencyStretch) {
        ARTS_ASSERT(np == 1);
        if (x_t[ji[q][0]] != 0) {
          do_sensor = true;
          Vector w;
          polynomial_basis_func(w, f_backend, 1);
          for (Index i = 0; i < w.nelem(); i++) {
            f_backend[i] += w[i] * x_t[ji[q][0]];
          }
        }
      } else {
        ARTS_ASSERT(0);
      }
    }

    // Baseline fit: polynomial or sinusoidal
    // ----------------------------------------------------------------------------
    else if (jacobian_quantities[q] == Jacobian::Sensor::Polyfit ||
             jacobian_quantities[q] == Jacobian::Sensor::Sinefit) {
      if (!yb_set) {
        yb_set = true;
        Index y_size = sensor_los.nrows() * sensor_response_f_grid.nelem() *
                       sensor_response_pol_grid.nelem() *
                       sensor_response_dlos_grid.nrows();
        y_baseline.resize(y_size);
        y_baseline = 0;
      }

      for (Index mb = 0; mb < sensor_los.nrows(); ++mb) {
        calcBaselineFit(y_baseline,
                        x_t,
                        mb,
                        sensor_response,
                        sensor_response_pol_grid,
                        sensor_response_f_grid,
                        sensor_response_dlos_grid,
                        jacobian_quantities[q],
                        q,
                        ji);
      }
    }
  }

  // *y_baseline* not yet set?
  if (!yb_set) {
    y_baseline.resize(1);
    y_baseline[0] = 0;
  }

  // Recalculate sensor_response?
  if (do_sensor) {
    sensor_response_agendaExecute(ws,
                                  sensor_response,
                                  sensor_response_f,
                                  sensor_response_f_grid,
                                  sensor_response_pol,
                                  sensor_response_pol_grid,
                                  sensor_response_dlos,
                                  sensor_response_dlos_grid,
                                  mblock_dlos_grid,
                                  f_backend,
                                  sensor_response_agenda);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void x2artsSpectroscopy(const Verbosity&) {
  ARTS_USER_ERROR ( "Retrievals of spectroscopic variables not yet handled.");
}


#ifdef OEM_SUPPORT
/* Workspace method: Doxygen documentation will be auto-generated */
void OEM(Workspace& ws,
         Vector& x,
         Vector& yf,
         Matrix& jacobian,
         Matrix& dxdy,
         Vector& oem_diagnostics,
         Vector& lm_ga_history,
         ArrayOfString& errors,
         const Vector& xa,
         const CovarianceMatrix& covmat_sx,
         const Vector& y,
         const CovarianceMatrix& covmat_se,
         const ArrayOfRetrievalQuantity& jacobian_quantities,
         const Agenda& inversion_iterate_agenda,
         const String& method,
         const Numeric& max_start_cost,
         const Vector& x_norm,
         const Index& max_iter,
         const Numeric& stop_dx,
         const Vector& lm_ga_settings,
         const Index& clear_matrices,
         const Index& display_progress,
         const Verbosity&) {
  // Main sizes
  const Index n = covmat_sx.nrows();
  const Index m = y.nelem();

  // Checks
  covmat_sx.compute_inverse();
  covmat_se.compute_inverse();

  OEM_checks(ws,
             x,
             yf,
             jacobian,
             inversion_iterate_agenda,
             xa,
             covmat_sx,
             y,
             covmat_se,
             jacobian_quantities,
             method,
             x_norm,
             max_iter,
             stop_dx,
             lm_ga_settings,
             clear_matrices,
             display_progress);

  // Size diagnostic output and init with NaNs
  oem_diagnostics.resize(5);
  oem_diagnostics = NAN;
  //
  if (method == "ml" || method == "lm" || method == "ml_cg" ||
      method == "lm_cg") {
    lm_ga_history.resize(max_iter + 1);
    lm_ga_history = NAN;
  } else {
    lm_ga_history.resize(0);
  }

  // Check for start vector and precomputed yf, jacobian
  if (x.nelem() != n) {
    x = xa;
    yf.resize(0);
    jacobian.resize(0, 0);
  }

  // If no precomputed value given, we compute yf and jacobian to
  // compute initial cost (and use in the first OEM iteration).
  if (yf.nelem() == 0) {
    inversion_iterate_agendaExecute(
        ws, yf, jacobian, xa, 1, 0, inversion_iterate_agenda);
  }

  ARTS_USER_ERROR_IF (yf.nelem() not_eq y.nelem(),
      "Mismatch between simulated y and input y.\n"
      "Input y is size ", y.nelem(), " but simulated y is ",
      yf.nelem(), "\n"
      "Use your frequency grid vector and your sensor response matrix to match simulations with measurements.\n")

  // TODO: Get this from invlib log.
  // Start value of cost function
  Numeric cost_start = NAN;
  if (method == "ml" || method == "lm" || display_progress ||
      max_start_cost > 0) {
    Vector dy = y;
    dy -= yf;
    Vector sdy = y;
    mult_inv(sdy, covmat_se, dy);
    Vector dx = x;
    dx -= xa;
    Vector sdx = x;
    mult_inv(sdx, covmat_sx, dx);
    cost_start = dx * sdx + dy * sdy;
    cost_start /= static_cast<Numeric>(m);
  }
  oem_diagnostics[1] = cost_start;

  // Handle cases with too large start cost
  if (max_start_cost > 0 && cost_start > max_start_cost) {
    // Flag no inversion in oem_diagnostics, and let x to be undefined
    oem_diagnostics[0] = 99;
    //
    if (display_progress) {
      cout << "\n   No OEM inversion, too high start cost:\n"
           << "        Set limit : " << max_start_cost << endl
           << "      Found value : " << cost_start << endl
           << endl;
    }
  }
  // Otherwise do inversion
  else {
    bool apply_norm = false;
    oem::Matrix T{};
    if (x_norm.nelem() == n) {
      T.resize(n, n);
      T *= 0.0;
      T.diagonal() = x_norm;
      for (Index i = 0; i < n; i++) {
        T(i, i) = x_norm[i];
      }
      apply_norm = true;
    }

    oem::CovarianceMatrix Se(covmat_se), Sa(covmat_sx);
    oem::Vector xa_oem(xa), y_oem(y), x_oem(x);
    oem::AgendaWrapper aw(&ws,
                          (unsigned int)m,
                          (unsigned int)n,
                          jacobian,
                          yf,
                          &inversion_iterate_agenda);
    oem::OEM_STANDARD<oem::AgendaWrapper> oem(aw, xa_oem, Sa, Se);
    oem::OEM_MFORM<oem::AgendaWrapper> oem_m(aw, xa_oem, Sa, Se);
    int oem_verbosity = static_cast<int>(display_progress);

    int return_code = 0;

    try {
      if (method == "li") {
        oem::Std s(T, apply_norm);
        oem::GN gn(stop_dx, 1, s);  // Linear case, only one step.
        return_code = oem.compute<oem::GN, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history, true);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "li_m") {
        oem::Std s(T, apply_norm);
        oem::GN gn(stop_dx, 1, s);  // Linear case, only one step.
        return_code = oem_m.compute<oem::GN, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history, true);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "li_cg") {
        oem::CG cg(T, apply_norm, 1e-10, 0);
        oem::GN_CG gn(stop_dx, 1, cg);  // Linear case, only one step.
        return_code = oem.compute<oem::GN_CG, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history, true);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "li_cg_m") {
        oem::CG cg(T, apply_norm, 1e-10, 0);
        oem::GN_CG gn(stop_dx, 1, cg);  // Linear case, only one step.
        return_code = oem_m.compute<oem::GN_CG, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history, true);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "gn") {
        oem::Std s(T, apply_norm);
        oem::GN gn(stop_dx, (unsigned int)max_iter, s);
        return_code = oem.compute<oem::GN, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "gn_m") {
        oem::Std s(T, apply_norm);
        oem::GN gn(stop_dx, (unsigned int)max_iter, s);
        return_code = oem_m.compute<oem::GN, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "gn_cg") {
        oem::CG cg(T, apply_norm, 1e-10, 0);
        oem::GN_CG gn(stop_dx, (unsigned int)max_iter, cg);
        return_code = oem.compute<oem::GN_CG, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "gn_cg_m") {
        oem::CG cg(T, apply_norm, 1e-10, 0);
        oem::GN_CG gn(stop_dx, (unsigned int)max_iter, cg);
        return_code = oem_m.compute<oem::GN_CG, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if ((method == "lm") || (method == "ml")) {
        oem::Std s(T, apply_norm);
        Sparse diagonal = Sparse::diagonal(covmat_sx.inverse_diagonal());
        CovarianceMatrix SaDiag{};
        SaDiag.add_correlation_inverse(Block(Range(0, n),
                                             Range(0, n),
                                             std::make_pair(0, 0),
                                             make_shared<Sparse>(diagonal)));
        oem::CovarianceMatrix SaInvLM = inv(oem::CovarianceMatrix(SaDiag));
	oem::LM lm(SaInvLM, s);

        lm.set_tolerance(stop_dx);
        lm.set_maximum_iterations((unsigned int)max_iter);
        lm.set_lambda(lm_ga_settings[0]);
        lm.set_lambda_decrease(lm_ga_settings[1]);
        lm.set_lambda_increase(lm_ga_settings[2]);
        lm.set_lambda_maximum(lm_ga_settings[3]);
        lm.set_lambda_threshold(lm_ga_settings[4]);
        lm.set_lambda_constraint(lm_ga_settings[5]);

        return_code = oem.compute<oem::LM&, oem::ArtsLog>(
            x_oem, y_oem, lm, oem_verbosity, lm_ga_history);
        oem_diagnostics[0] = static_cast<Index>(return_code);
        if (lm.get_lambda() > lm.get_lambda_maximum()) {
          oem_diagnostics[0] = 2;
        }
      } else if ((method == "lm_cg") || (method == "ml_cg")) {
        oem::CG cg(T, apply_norm, 1e-10, 0);

        Sparse diagonal = Sparse::diagonal(covmat_sx.inverse_diagonal());
        CovarianceMatrix SaDiag{};
        SaDiag.add_correlation_inverse(Block(Range(0, n),
                                             Range(0, n),
                                             std::make_pair(0, 0),
                                             make_shared<Sparse>(diagonal)));
	oem::LM_CG lm(SaDiag, cg);

        lm.set_maximum_iterations((unsigned int)max_iter);
        lm.set_lambda(lm_ga_settings[0]);
        lm.set_lambda_decrease(lm_ga_settings[1]);
        lm.set_lambda_increase(lm_ga_settings[2]);
        lm.set_lambda_threshold(lm_ga_settings[3]);
        lm.set_lambda_maximum(lm_ga_settings[4]);

        return_code = oem.compute<oem::LM_CG&, oem::ArtsLog>(
            x_oem, y_oem, lm, oem_verbosity, lm_ga_history);
        oem_diagnostics[0] = static_cast<Index>(return_code);
        if (lm.get_lambda() > lm.get_lambda_maximum()) {
          oem_diagnostics[0] = 2;
        }
      }

      oem_diagnostics[2] = oem.cost / static_cast<Numeric>(m);
      oem_diagnostics[3] = oem.cost_y / static_cast<Numeric>(m);
      oem_diagnostics[4] = static_cast<Numeric>(oem.iterations);
    } catch (const std::exception& e) {
      oem_diagnostics[0] = 9;
      oem_diagnostics[2] = oem.cost;
      oem_diagnostics[3] = oem.cost_y;
      oem_diagnostics[4] = static_cast<Numeric>(oem.iterations);
      x_oem *= NAN;
      std::vector<std::string> sv = oem::handle_nested_exception(e);
      for (auto& s : sv) {
        std::stringstream ss{s};
        std::string t{};
        while (std::getline(ss, t)) {
          errors.push_back(t.c_str());
        }
      }
    } catch (...) {
      throw;
    }

    x = x_oem;
    yf = aw.get_measurement_vector();

    // Shall empty jacobian and dxdy be returned?
    if (clear_matrices) {
      jacobian.resize(0, 0);
      dxdy.resize(0, 0);
    } else if (oem_diagnostics[0] <= 2) {
      dxdy.resize(n, m);
      Matrix tmp1(n, m), tmp2(n, n), tmp3(n, n);
      mult_inv(tmp1, transpose(jacobian), covmat_se);
      mult(tmp2, tmp1, jacobian);
      add_inv(tmp2, covmat_sx);
      inv(tmp3, tmp2);
      mult(dxdy, tmp3, tmp1);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void covmat_soCalc(Matrix& covmat_so,
                   const Matrix& dxdy,
                   const CovarianceMatrix& covmat_se,
                   const Verbosity& /*v*/) {
  Index n(dxdy.nrows()), m(dxdy.ncols());
  Matrix tmp1(m, n);

  ARTS_USER_ERROR_IF ((m == 0) || (n == 0),
        "The gain matrix *dxdy* is required to compute the observation error covariance matrix.");
  ARTS_USER_ERROR_IF ((covmat_se.nrows() != m) || (covmat_se.ncols() != m),
        "The covariance matrix covmat_se has invalid dimensions.");

  covmat_so.resize(n, n);
  mult(tmp1, covmat_se, transpose(dxdy));
  mult(covmat_so, dxdy, tmp1);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void covmat_ssCalc(Matrix& covmat_ss,
                   const Matrix& avk,
                   const CovarianceMatrix& covmat_sx,
                   const Verbosity& /*v*/) {
  Index n(avk.ncols());
  Matrix tmp1(n, n), tmp2(n, n);

  ARTS_USER_ERROR_IF (n == 0,
        "The averaging kernel matrix *dxdy* is required to compute the smoothing error covariance matrix.");
  ARTS_USER_ERROR_IF ((covmat_sx.nrows() != n) || (covmat_sx.ncols() != n),
        "The covariance matrix *covmat_sx* invalid dimensions.");

  covmat_ss.resize(n, n);

  // Sign doesn't matter since we're dealing with a quadratic form.
  id_mat(tmp1);
  tmp1 -= avk;

  mult(tmp2, covmat_sx, transpose(tmp1));
  mult(covmat_ss, tmp1, tmp2);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixFromCovarianceMatrix(Matrix& S,
                                const CovarianceMatrix& Sc,
                                const Verbosity& /*v*/) {
  S = Matrix(Sc);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void avkCalc(Matrix& avk,
             const Matrix& dxdy,
             const Matrix& jacobian,
             const Verbosity& /*v*/) {
  Index m(jacobian.nrows()), n(jacobian.ncols());
  ARTS_USER_ERROR_IF ((m == 0) || (n == 0),
                      "The Jacobian matrix is empty.");
  ARTS_USER_ERROR_IF ((dxdy.nrows() != n) || (dxdy.ncols() != m),
      "Matrices have inconsistent sizes.\n"
      "  Size of gain matrix: ", dxdy.nrows(), " x ", dxdy.ncols(),
      "\n"
      "     Size of Jacobian: ", jacobian.nrows(), " x ",
      jacobian.ncols(), "\n")

  avk.resize(n, n);
  mult(avk, dxdy, jacobian);
}

#else

void covmat_soCalc(Matrix& /* covmat_so */,
                   const Matrix& /* dxdy */,
                   const CovarianceMatrix& /* covmat_se*/,
                   const Verbosity& /*v*/) {
  ARTS_USER_ERROR (
      "WSM is not available because ARTS was compiled without "
      "OEM support.");
}

void covmat_ssCalc(Matrix& /*covmat_ss*/,
                   const Matrix& /*avk*/,
                   const CovarianceMatrix& /*covmat_sx*/,
                   const Verbosity& /*v*/) {
  ARTS_USER_ERROR (
      "WSM is not available because ARTS was compiled without "
      "OEM support.");
}

void avkCalc(Matrix& /* avk */,
             const Matrix& /* dxdy */,
             const Matrix& /* jacobian */,
             const Verbosity& /*v*/) {
  ARTS_USER_ERROR (
      "WSM is not available because ARTS was compiled without "
      "OEM support.");
}

void OEM(Workspace&,
         Vector&,
         Vector&,
         Matrix&,
         Matrix&,
         Vector&,
         Vector&,
         ArrayOfString&,
         const Vector&,
         const CovarianceMatrix&,
         const Vector&,
         const CovarianceMatrix&,
         const Index&,
         const ArrayOfRetrievalQuantity&,
         const ArrayOfArrayOfIndex&,
         const Agenda&,
         const String&,
         const Numeric&,
         const Vector&,
         const Index&,
         const Numeric&,
         const Vector&,
         const Index&,
         const Index&,
         const Verbosity&) {
  ARTS_USER_ERROR (
      "WSM is not available because ARTS was compiled without "
      "OEM support.");
}
#endif

#if defined(OEM_SUPPORT) && 0

#include "agenda_wrapper_mpi.h"
#include "oem_mpi.h"

//
// Performs manipulations of workspace variables necessary for distributed
// retrievals with MPI:
//
//   - Splits up sensor positions evenly over processes
//   - Splits up inverse covariance matrices.
//
void MPI_Initialize(Matrix& sensor_los,
                    Matrix& sensor_pos,
                    Vector& sensor_time) {
  int initialized;

  MPI_Initialized(&initialized);
  if (!initialized) {
    MPI_Init(nullptr, nullptr);
  }

  int rank, nprocs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  int nmblock = (int)sensor_pos.nrows();
  int mblock_range = nmblock / nprocs;
  int mblock_start = mblock_range * rank;
  int remainder = nmblock % std::max(mblock_range, nprocs);

  //
  // Split up sensor positions.
  //

  if (rank < remainder) {
    mblock_range += 1;
    mblock_start += rank;
  } else {
    mblock_start += remainder;
  }

  if (nmblock > 0) {
    Range range = Range(mblock_start, mblock_range);

    Matrix tmp_m = sensor_los(range, joker);
    sensor_los = tmp_m;

    tmp_m = sensor_pos(range, joker);
    sensor_pos = tmp_m;

    Vector tmp_v = sensor_time[range];
    sensor_time = tmp_v;
  } else {
    sensor_los.resize(0, 0);
    sensor_pos.resize(0, 0);
    sensor_time.resize(0);
  }
}

void OEM_MPI(Workspace& ws,
             Vector& x,
             Vector& yf,
             Matrix& jacobian,
             Matrix& dxdy,
             Vector& oem_diagnostics,
             Vector& lm_ga_history,
             Matrix& sensor_los,
             Matrix& sensor_pos,
             Vector& sensor_time,
             CovarianceMatrix& covmat_sx,
             CovarianceMatrix& covmat_se,
             const Vector& xa,
             const Vector& y,
             const ArrayOfRetrievalQuantity& jacobian_quantities,
             const Agenda& inversion_iterate_agenda,
             const String& method,
             const Numeric& max_start_cost,
             const Vector& x_norm,
             const Index& max_iter,
             const Numeric& stop_dx,
             const Vector& lm_ga_settings,
             const Index& clear_matrices,
             const Index& display_progress,
             const Verbosity& /*v*/) {
  // Main sizes
  const Index n = covmat_sx.nrows();
  const Index m = y.nelem();

  // Check WSVs
  OEM_checks(ws,
             x,
             yf,
             jacobian,
             inversion_iterate_agenda,
             xa,
             covmat_sx,
             covmat_se,
             jacobian_quantities,
             method,
             x_norm,
             max_iter,
             stop_dx,
             lm_ga_settings,
             clear_matrices,
             display_progress);

  // Calculate spectrum and Jacobian for a priori state
  // Jacobian is also input to the agenda, and to flag this is this first
  // call, this WSV must be set to be empty.
  jacobian.resize(0, 0);

  // Initialize MPI environment.
  MPI_Initialize(sensor_los, sensor_pos, sensor_time);

  // Setup distributed matrices.
  MPICovarianceMatrix SeInvMPI(covmat_se);
  MPICovarianceMatrix SaInvMPI(covmat_sx);

  // Create temporary MPI vector from local results and use conversion to
  // standard vector to broadcast results to all processes.
  oem::Vector tmp;
  inversion_iterate_agendaExecute(
      ws, tmp, jacobian, xa, 1, inversion_iterate_agenda);
  yf = MPIVector(tmp);

  // Size diagnostic output and init with NaNs
  oem_diagnostics.resize(5);
  oem_diagnostics = NAN;
  //
  if (method == "ml" || method == "lm") {
    lm_ga_history.resize(max_iter);
    lm_ga_history = NAN;
  } else {
    lm_ga_history.resize(0);
  }

  // Start value of cost function. Covariance matrices are already distributed
  // over processes, so we need to use invlib matrix algebra.
  Numeric cost_start = NAN;
  if (method == "ml" || method == "lm" || display_progress ||
      max_start_cost > 0) {
    oem::Vector dy = y;
    dy -= yf;
    cost_start = dot(dy, SeInvMPI * dy);
  }
  oem_diagnostics[1] = cost_start;

  // Handle cases with too large start cost
  if (max_start_cost > 0 && cost_start > max_start_cost) {
    // Flag no inversion in oem_diagnostics, and let x to be undefined
    oem_diagnostics[0] = 99;
    //
    if (display_progress) {
      cout << "\n   No OEM inversion, too high start cost:\n"
           << "        Set limit : " << max_start_cost << endl
           << "      Found value : " << cost_start << endl
           << endl;
    }
  }

  // Otherwise do inversion
  else {
    // Size remaining output arguments
    x.resize(n);
    dxdy.resize(n, m);

    OEMVector xa_oem(xa), y_oem(y), x_oem;
    oem::AgendaWrapperMPI aw(&ws, &inversion_iterate_agenda, m, n);

    OEM_PS_PS_MPI<AgendaWrapperMPI> oem(aw, xa_oem, SaInvMPI, SeInvMPI);

    // Call selected method
    int return_value = 99;

    if (method == "li") {
      oem::CG cg(1e-12, 0);
      oem::GN_CG gn(stop_dx, (unsigned int)max_iter, cg);
      return_value = oem.compute<oem::GN_CG, invlib::MPILog>(
          x_oem, y_oem, gn, 2 * (int)display_progress);
    } else if (method == "gn") {
	    oem::CG cg(1e-12, 0);
      oem::GN_CG gn(stop_dx, (unsigned int)max_iter, cg);
      return_value = oem.compute<oem::GN_CG, invlib::MPILog>(
          x_oem, y_oem, gn, 2 * (int)display_progress);
    } else if ((method == "lm") || (method == "ml")) {
	    oem::CG cg(1e-12, 0);
      LM_CG_MPI lm(SaInvMPI, cg);

      lm.set_tolerance(stop_dx);
      lm.set_maximum_iterations((unsigned int)max_iter);
      lm.set_lambda(lm_ga_settings[0]);
      lm.set_lambda_decrease(lm_ga_settings[1]);
      lm.set_lambda_increase(lm_ga_settings[2]);
      lm.set_lambda_threshold(lm_ga_settings[3]);
      lm.set_lambda_maximum(lm_ga_settings[4]);

      return_value = oem.compute<oem::LM_CG_MPI, invlib::MPILog>(
          x_oem, y_oem, lm, 2 * (int)display_progress);
    }

    oem_diagnostics[0] = return_value;
    oem_diagnostics[2] = oem.cost;
    oem_diagnostics[3] = oem.cost_y;
    oem_diagnostics[4] = static_cast<Numeric>(oem.iterations);

    x = x_oem;
    // Shall empty jacobian and dxdy be returned?
    if (clear_matrices && (oem_diagnostics[0])) {
      jacobian.resize(0, 0);
      dxdy.resize(0, 0);
    }
  }
  MPI_Finalize();
}

#else

void OEM_MPI(Workspace&,
             Vector&,
             Vector&,
             Matrix&,
             Matrix&,
             Vector&,
             Vector&,
             Matrix&,
             Matrix&,
             Vector&,
             CovarianceMatrix&,
             CovarianceMatrix&,
             const Vector&,
             const Vector&,
             const Index&,
             const ArrayOfRetrievalQuantity&,
             const Agenda&,
             const String&,
             const Numeric&,
             const Vector&,
             const Index&,
             const Numeric&,
             const Vector&,
             const Index&,
             const Index&,
             const Verbosity&) {
  ARTS_USER_ERROR (
      "You have to compile ARTS with OEM support "
      " and enable MPI to use OEM_MPI.");
}

#endif  // OEM_SUPPORT && ENABLE_MPI
