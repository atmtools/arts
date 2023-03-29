/* Copyright (C) 2011-2017
   Jana Mendrok     <jana.mendrok@gmail.com>
   Daniel Kreyling  <daniel.kreyling@nict.go.jp>
   Manfred Brath    <manfred.brath@uni-hamburg.de>
   Patrick Eriksson <patrick.eriksson@chalmers.se>

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
  \file   m_microphysics.cc
  \author Jana Mendrok, Daniel Kreyling, Manfred Brath, Patrick Eriksson
  \date   2017-07-10

  \brief  Workspace functions related to particle micophysics (e.g. size
          distributions).

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h. */

/*===========================================================================
  === External declarations
  ===========================================================================*/
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "array.h"
#include "arts.h"
#include "arts_constants.h"
#include "auto_md.h"
#include "check_input.h"
#include "cloudbox.h"
#include "disort.h"
#include "file.h"
#include "interpolation.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "microphysics.h"
#include "optproperties.h"
#include "parameters.h"
#include "physics_funcs.h"
#include "psd.h"
#include "rte.h"
#include "sorting.h"
#include "special_interp.h"
#include "xml_io.h"

inline constexpr Numeric PI=Constant::pi;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void HydrotableCalc(Workspace& ws,
                    GriddedField4& hydrotable,
                    const ArrayOfAgenda& pnd_agenda_array,
                    const ArrayOfArrayOfString& pnd_agenda_array_input_names,
                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                    const Index& scat_data_checked,
                    const Vector& f_grid,
                    const Index& iss,
                    const Vector& T_grid,
                    const Vector& wc_grid,                    
                    const Verbosity&)
{
  // Sizes
  const Index nss = scat_data.nelem(); 
  const Index nf = f_grid.nelem();
  const Index nt = T_grid.nelem();
  const Index nw = wc_grid.nelem();

  ARTS_USER_ERROR_IF (pnd_agenda_array.nelem() != nss,
        "*scat_data* and *pnd_agenda_array* are inconsistent "
        "in size.");
  ARTS_USER_ERROR_IF (pnd_agenda_array_input_names.nelem() != nss,
        "*scat_data* and *pnd_agenda_array_input_names* are "
        "inconsistent in size.");
  ARTS_USER_ERROR_IF (pnd_agenda_array_input_names[iss].nelem() != 1,
        "This method requires one-moment PSDs, but *pnd_agenda_array_input_names* "
        "for the selected scattering species does not have length one.");
  ARTS_USER_ERROR_IF (!scat_data_checked,
                      "The scat_data must be flagged to have passed a "
                      "consistency check (scat_data_checked=1).");

  // Allocate *hydrotable*
  hydrotable.set_name("Table of particle optical properties");
  hydrotable.data.resize(4, nf, nt, nw);
  //
  hydrotable.set_grid_name(0, "Quantity");
  hydrotable.set_grid(0, {"Extinction [m-1]",
                          "Single scattering albedo [-]",
                          "Asymmetry parameter [-]",
                          "Radar reflectivity [m2]"});
  hydrotable.set_grid_name(1, "Frequency [Hz]");
  hydrotable.set_grid(1, f_grid);
  hydrotable.set_grid_name(2, "Temperature [K]");
  hydrotable.set_grid(2, T_grid);
  hydrotable.set_grid_name(3, "Particle content [kg/m3]");
  hydrotable.set_grid(3, wc_grid);

  // Scattering angle grid
  const Index nsa = 361;
  Vector sa_grid;
  nlinspace(sa_grid, 0, 180, nsa);

  // Local variables
  Matrix pnd_data;
  Tensor3 dpnd_data_dx;
  Matrix pnd_agenda_input(nt, 1);
  ArrayOfString dpnd_data_dx_names(0);
  Matrix ext(nf, nt);
  Matrix abs(nf, nt);
  Tensor3 pfun(nf, nt, nsa);
  const Numeric fourpi = 4.0 * PI;
  ArrayOfIndex cbox_limits = {0, nt-1};
  
  // Loop and fill table
  for (Index iw = 0; iw < nw; iw++) {
    // Call pnd_agenda
    pnd_agenda_input = wc_grid[iw];
    pnd_agenda_arrayExecute(ws,
                            pnd_data,
                            dpnd_data_dx,
                            iss,
                            T_grid,
                            pnd_agenda_input,
                            pnd_agenda_array_input_names[iss],
                            dpnd_data_dx_names,
                            pnd_agenda_array);

    // Calculate extinsion, absorbtion and phase function
    ext = 0.0;
    abs = 0.0;
    pfun = 0.0;
    ext_abs_pfun_from_tro(ext,
                          abs,
                          pfun,
                          scat_data[iss],
                          iss,
                          transpose(pnd_data),
                          cbox_limits,
                          T_grid,
                          sa_grid);
    
    // Fill the hydrotable for present particle content
    for (Index iv = 0; iv < nf; iv++) {
      for (Index it = 0; it < nt; it++) {
        hydrotable.data(0,iv,it,iw) = ext(iv,it);
        hydrotable.data(1,iv,it,iw) = 1.0 - (abs(iv,it) / ext(iv,it));
        hydrotable.data(2,iv,it,iw) = asymmetry_parameter(sa_grid,
                                                          pfun(iv,it,joker));
        hydrotable.data(3,iv,it,iw) = fourpi * pfun(iv,it,nsa-1);
      }
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void particle_massesFromMetaDataSingleCategory(
    Matrix& particle_masses,
    const ArrayOfArrayOfScatteringMetaData& scat_meta,
    const Verbosity&) {
  const Index np_total = TotalNumberOfElements(scat_meta);

  particle_masses.resize(np_total, 1);

  Index i_se_flat = 0;
  for (Index i_ss = 0; i_ss < scat_meta.nelem(); i_ss++) {
    for (Index i_se = 0; i_se < scat_meta[i_ss].nelem(); i_se++) {
      ARTS_USER_ERROR_IF (std::isnan(scat_meta[i_ss][i_se].mass) ||
          scat_meta[i_ss][i_se].mass <= 0 || scat_meta[i_ss][i_se].mass > 1.,
          "A presumably incorrect value found for "
          "scat_meta[", i_ss, "][", i_se, "].mass.\n"
          "The value is ", scat_meta[i_ss][i_se].mass)

      particle_masses(i_se_flat, 0) = scat_meta[i_ss][i_se].mass;

      i_se_flat++;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void particle_massesFromMetaData(  //WS Output:
    Matrix& particle_masses,
    // WS Input:
    const ArrayOfArrayOfScatteringMetaData& scat_meta,
    const Verbosity&) {
  // resize particle_masses to required diemsions and properly initialize values
  particle_masses.resize(TotalNumberOfElements(scat_meta), scat_meta.nelem());
  particle_masses = 0.;

  // calculate and set particle_masses
  Index i_se_flat = 0;
  for (Index i_ss = 0; i_ss < scat_meta.nelem(); i_ss++) {
    for (Index i_se = 0; i_se < scat_meta[i_ss].nelem(); i_se++) {
      ARTS_USER_ERROR_IF (std::isnan(scat_meta[i_ss][i_se].mass) ||
          scat_meta[i_ss][i_se].mass <= 0 || scat_meta[i_ss][i_se].mass > 1.,
          "A presumably incorrect value found for "
          "scat_meta[", i_ss, "][", i_se, "].mass.\n"
          "The value is ", scat_meta[i_ss][i_se].mass)

      particle_masses(i_se_flat, i_ss) = scat_meta[i_ss][i_se].mass;
      i_se_flat++;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void pndFromPsdBasic(Matrix& pnd_data,
                     Tensor3& dpnd_data_dx,
                     const Vector& pnd_size_grid,
                     const Matrix& psd_data,
                     const Vector& psd_size_grid,
                     const Tensor3& dpsd_data_dx,
                     const Index& quad_order,
                     const Verbosity&) {
  // Some sizes
  const Index np = psd_data.nrows();
  const Index ng = psd_size_grid.nelem();
  Index ndx = 0;
  const bool do_dx = !dpsd_data_dx.empty();

  // Checks
  ARTS_USER_ERROR_IF (ng < 2,
        "The method requires that length of *psd_size_grid* is >= 2.");
  ARTS_USER_ERROR_IF (ng != pnd_size_grid.nelem(),
        "So far, the method requires that *psd_size_grid* and "
        "*pnd_size_grid* have same length.");
  for (Index i = 0; i < ng; i++) {
    ARTS_USER_ERROR_IF (psd_size_grid[i] != pnd_size_grid[i],
          "So far, the method requires that *psd_size_grid* and "
          "*pnd_size_grid* are identical.");
  }
  ARTS_USER_ERROR_IF (psd_data.ncols() != ng,
        "Number of columns in *psd_data* and length of "
        "*psd_size_grid* must match.");

  pnd_data.resize(np, ng);
  if (do_dx) {
    ARTS_USER_ERROR_IF (dpsd_data_dx.ncols() != ng,
          "Number of columns in *dpsd_data_dx* and length of "
          "*psd_size_grid* must match.");
    ndx = dpsd_data_dx.npages();
    dpnd_data_dx.resize(ndx, np, ng);
  } else {
    dpnd_data_dx.resize(0, 0, 0);
  }

  // Get sorted version of psd_size_grid (and, since pnd_size_grid so far is
  // identical, of this as well implicitly)
  ArrayOfIndex intarr;
  Vector psd_size_grid_sorted(ng);
  get_sorted_indexes(intarr, psd_size_grid);
  for (Index i = 0; i < ng; i++)
    psd_size_grid_sorted[i] = psd_size_grid[intarr[i]];

  ARTS_USER_ERROR_IF (!is_increasing(psd_size_grid_sorted),
        "*psd_size_grid* is not allowed to contain "
        "duplicate values.");

  // Calculate pnd by intrgation of psd for given nodes/bins
  Vector quadweights(ng);
  bin_quadweights(quadweights, psd_size_grid_sorted, quad_order);

  for (Index i = 0; i < ng; i++) {
    for (Index ip = 0; ip < np; ip++) {
      pnd_data(ip, intarr[i]) = quadweights[i] * psd_data(ip, intarr[i]);
    }

    if (do_dx) {
      for (Index ip = 0; ip < np; ip++) {
        for (Index ix = 0; ix < ndx; ix++) {
          dpnd_data_dx(ix, ip, intarr[i]) =
              quadweights[i] * dpsd_data_dx(ix, ip, intarr[i]);
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void pndFromPsd(Matrix& pnd_data,
                Tensor3& dpnd_data_dx,
                const Vector& pnd_size_grid,
                const Matrix& psd_data,
                const Vector& psd_size_grid,
                const Tensor3& dpsd_data_dx,
                const ArrayOfArrayOfSingleScatteringData& scat_data,
                const Vector& f_grid,
                const Index& scat_data_checked,
                const Index& quad_order,
                const Index& scat_index,
                const Numeric& threshold_rsec,
                const Numeric& threshold_bext,
                const Numeric& threshold_rpnd,
                const Verbosity&) {
  // Some sizes
  const Index np = psd_data.nrows();
  const Index ng = psd_size_grid.nelem();
  Index ndx = 0;
  const bool do_dx = !dpsd_data_dx.empty();

  // Checks
  ARTS_USER_ERROR_IF (ng < 3,
        "The method requires that length of *psd_size_grid* is >= 3.");
  ARTS_USER_ERROR_IF (ng != pnd_size_grid.nelem(),
        "So far, the method requires that *psd_size_grid* and"
        " *pnd_size_grid* have the same length.");
  for (Index i = 0; i < ng; i++) {
    ARTS_USER_ERROR_IF (psd_size_grid[i] != pnd_size_grid[i],
          "So far, the method requires that *psd_size_grid* and"
          " *pnd_size_grid* are identical.");
  }
  ARTS_USER_ERROR_IF (psd_data.ncols() != ng,
        "Number of columns in *psd_data* and length of"
        " *psd_size_grid* must match.");

  pnd_data.resize(np, ng);
  if (do_dx) {
    ARTS_USER_ERROR_IF (dpsd_data_dx.ncols() != ng,
          "Number of columns in *dpsd_data_dx* and length of"
          " *psd_size_grid* must match.");
    ndx = dpsd_data_dx.npages();
    dpnd_data_dx.resize(ndx, np, ng);
  } else {
    dpnd_data_dx.resize(0, 0, 0);
  }

  ARTS_USER_ERROR_IF (!scat_data_checked,
        "*scat_data* must have passed a consistency check"
        " (scat_data_checked=1).\n"
        "Alternatively, use *pndFromPsdBasic*.");
  ARTS_USER_ERROR_IF (scat_index >= scat_data.nelem(),
        "*scat_index* exceeds the number of available"
        " scattering species.");
  ARTS_USER_ERROR_IF (scat_data[scat_index].nelem() != ng,
        "Number of scattering elements in this scattering"
        " species (*scat_index*) inconsistent with length of"
        " *pnd_size_grid*.");

  // Get sorted version of psd_size_grid (and, since pnd_size_grid so far is
  // identical, of this as well implicitly)
  ArrayOfIndex intarr;
  Vector psd_size_grid_sorted(ng);
  get_sorted_indexes(intarr, psd_size_grid);
  for (Index i = 0; i < ng; i++)
    psd_size_grid_sorted[i] = psd_size_grid[intarr[i]];

  ARTS_USER_ERROR_IF (!is_increasing(psd_size_grid_sorted),
        "*psd_size_grid* is not allowed to contain "
        "duplicate values.");

  // Calculate pnd by integration of psd for given nodes/bins
  Vector quadweights(ng);
  bin_quadweights(quadweights, psd_size_grid_sorted, quad_order);

  for (Index i = 0; i < ng;
       i++)  //loop over pnd_size_grid aka scattering elements
  {
    for (Index ip = 0; ip < np; ip++)  //loop over pressure levels
    {
      pnd_data(ip, intarr[i]) = quadweights[i] * psd_data(ip, intarr[i]);
    }

    if (do_dx) {
      for (Index ip = 0; ip < np; ip++) {
        for (Index ix = 0; ix < ndx; ix++) {
          dpnd_data_dx(ix, ip, intarr[i]) =
              quadweights[i] * dpsd_data_dx(ix, ip, intarr[i]);
        }
      }
    }
  }

  ArrayOfSingleScatteringData sds = scat_data[scat_index];
  Index fstart = 0;
  Index nf = f_grid.nelem();
  Matrix bulkext(np, nf, 0.);
  Vector ext(nf), ext_s0(nf), ext_s1(nf), ext_l0(nf), ext_l1(nf);

  // check that pnd_size_grid and scat_data cover (scalar) bulk extinction properly.
  // (formerly we used mass. but heavy particles might not necessarily
  // contribute much to optprops. same is valid for total particle number.)
  //
  // We do that by
  // (a) checking that psd*ext decreases at the edges, i.e. to increase the
  //     probability that the main extinction peak is covered and in order to
  //     avoid passing check (b) through making the edge bins tiny
  // (b) checking that the extinction contributions of the edge bins is small
  //     compared to the total bulk bin.
  // The tests are somewhat over-sensitive on the small particle side in the
  // sense that the checks are triggered easily by low mass density cases - which
  // have their main extinction contributed by small particles, but for which
  // the total extinction is quite low. Hence, there is also a total extinction
  // threshold, below which the checks are not applied. Furthermore, such low
  // densities do typically not occur throughout the whole atmospheric profile,
  // but higher density layers might exists, meaning the low density layers
  // contribute very little to the total optical properties, hence should not
  // dominate the check criteria (no need to be very strict with them if they
  // hardly contribute anyways). Therefore, checks are also not applied if the
  // edge bin number density at a given atmospheric level is low compared to the
  // maximum number of this scattering element in the profile.
  //
  // Technically, we need to cover all freqs separately (as optprops vary
  // strongly with freq). We indeed do that here.
  // FIXME: Shall a single freq or reduced freq-number option be implemented?
  // If so, the don't-apply-test criteria need to be checked again though (so it
  // does not happen that effectively all of the reduced-freq testing is skipped
  // due to them).
  //
  // Technically, we also need to use identical temperatures (and preferably the
  // ones used for deriving psd). But all scat elements can be on different
  // T-grids, ie we would need to interpolate. Instead we just use the first in
  // line.
  // FIXME: Maybe refine in future allowing to select between low, high, median
  // grid point or even specify one temp to use.
  //
  // And, technically ext might be directional dependent (for oriented
  // particles). Also this might be different for each scat element. However, at
  // least for now, we only check the 0-direction.
  // FIXME: Check further directions?

  for (Index ise = 0; ise < ng;
       ise++)  //loop over pnd_size_grid aka scattering elements
  {
    if (sds[intarr[ise]].ext_mat_data.nshelves() > 1)
      ext = sds[intarr[ise]].ext_mat_data(joker, 0, 0, 0, 0);
    else
      ext = sds[intarr[ise]].ext_mat_data(0, 0, 0, 0, 0);

    // keep the ext data of the edge bins
    if (ise == 0)
      ext_s0 = ext;
    else if (ise == 1)
      ext_s1 = ext;
    else if (ise == ng - 2)
      ext_l1 = ext;
    else if (ise == ng - 1)
      ext_l0 = ext;

    for (Index ip = 0; ip < np; ip++)  //loop over pressure levels
      if (abs(sum(pnd_data(ip, joker))) > 0.)
        for (Index f = fstart; f < (fstart + nf); f++)
          bulkext(ip, f) += pnd_data(ip, intarr[ise]) * ext[f];
  }

  Numeric max0 = 0, max1 = 0;
  for (Index ip = 0; ip < np; ip++)  //loop over pressure levels
  {
    max0 = max(abs(pnd_data(ip, intarr[0])), max0);
    max1 = max(abs(pnd_data(ip, intarr[ng - 1])), max1);
  }

  Numeric contrib;
  for (Index ip = 0; ip < np; ip++)  //loop over pressure levels
  {
    if (abs(sum(pnd_data(ip, joker))) > 0.) {
      for (Index f = fstart; f < (fstart + nf); f++) {
        /*        for( Index ise=0; ise<ng; ise++ )
        {
          if( sds[ise].ext_mat_data.nshelves()>1 )
            ext = sds[ise].ext_mat_data(joker,0,0,0,0);
          else
            ext = sds[ise].ext_mat_data(0,0,0,0,0);
          cout << "    se #" << ise << ": contrib = pnd*ext/bext = "
               << abs(pnd_data(ip,ise)) << "*" << ext[f] << "/"
               << abs(bulkext(ip,f)) << " = "
               << 1e2*abs(pnd_data(ip,ise))*ext[f]/abs(bulkext(ip,f))
               << "%\n";
        }*/

        // check that bin-width normalized extinction (or ext*psd) is
        // decreasing.
        if (abs(bulkext(ip, f)) > 1e-2 * threshold_bext) {
          ARTS_USER_ERROR_IF (abs(psd_data(ip, intarr[0])) > 0. and
              ext_s0[f] * abs(psd_data(ip, intarr[0])) >=
                  ext_s1[f] * abs(psd_data(ip, intarr[1])),
              "  Bin-width normalized extinction (ext*psd) not decreasing"
              " at small size edge\n"
              "  at atm level #", ip, " and freq point #", f, ".\n"
              "  ext_s0=", ext_s0[f],
              ", psd_s0=", abs(psd_data(ip, intarr[0])),
              ", ext_s0*psd_s0=", ext_s0[f] * abs(psd_data(ip, intarr[0])),
              "\n    LARGER EQUAL\n"
              "  ext_s1=", ext_s1[f],
              ", psd_s1=", abs(psd_data(ip, intarr[1])),
              ", ext_s1*psd_s1=", ext_s1[f] * abs(psd_data(ip, intarr[1])),
              "\n"
              "    Total bulk ext = ", abs(bulkext(ip, f)), "\n"
              "  Need to add smaller sized particles!\n")

          ARTS_USER_ERROR_IF (abs(psd_data(ip, intarr[ng - 1])) > 0. and
              ext_l0[f] * abs(psd_data(ip, intarr[ng - 1])) >=
                  ext_l1[f] * abs(psd_data(ip, intarr[ng - 2])),
              "Bin-width normalized extinction (ext*psd) not decreasing"
              " at large size edge\n"
              "at atm level #", ip, " and freq point #", f, ".\n"
              "  ext_l0=", ext_l0[f],
              ", psd_l0=", abs(psd_data(ip, intarr[ng - 1])),
              ", ext_l0*psd_l0=",
              ext_l0[f] * abs(psd_data(ip, intarr[ng - 1])),
              "\n    LARGER EQUAL\n"
              "  ext_l1=", ext_l1[f],
              ", psd_l1=", abs(psd_data(ip, intarr[ng - 2])),
              ", ext_l1*psd_l1=",
              ext_l1[f] * abs(psd_data(ip, intarr[ng - 2])), "\n"
              "    Total bulk ext = ", abs(bulkext(ip, f)), "\n"
              "  Need to add larger sized particles!\n")
        }

        // check that contribution of edge bins to total extinction is
        // sufficiently small
        if (abs(bulkext(ip, f)) > threshold_bext) {
          if (abs(pnd_data(ip, intarr[0])) > threshold_rpnd * max0) {
            contrib =
                abs(pnd_data(ip, intarr[0])) * ext_s0[f] / abs(bulkext(ip, f));
            //cout << "    small edge contrib = pnd*ext/bext = "
            //     << abs(pnd_data(ip,intarr[0])) << "*" << ext_s0[f] << "/"
            //     << abs(bulkext(ip,f)) << " = " << contrib << "\n";
            ARTS_USER_ERROR_IF (abs(pnd_data(ip, intarr[0])) * ext_s0[f] >
                threshold_rsec * abs(bulkext(ip, f)),
                "Contribution of edge bin to total extinction too high"
                " (", contrib * 1e2, "% of ", abs(bulkext(ip, f)),
                ") at small size edge\n"
                "at atm level #", ip, " and freq point #", f, ".\n"
                "  Need to add smaller sized particles or refine the size"
                " grid on the small size edge!\n")
          }
          if (abs(pnd_data(ip, intarr[ng - 1])) > threshold_rpnd * max1) {
            contrib = abs(pnd_data(ip, intarr[ng - 1])) * ext_l0[f] /
                      abs(bulkext(ip, f));
            //cout << "    large edge contrib = pnd*ext/bext = "
            //     << abs(pnd_data(ip,ng-1)) << "*" << ext_l0[f] << "/"
            //     << abs(bulkext(ip,f)) << " = " << contrib << "\n";
            ARTS_USER_ERROR_IF (abs(pnd_data(ip, intarr[ng - 1])) * ext_l0[f] >
                threshold_rsec * abs(bulkext(ip, f)),
                "Contribution of edge bin to total extinction too high"
                " (", contrib * 1e2, "% of ", abs(bulkext(ip, f)),
                ") at large size edge\n"
                "at atm level #", ip, " and freq point #", f, ".\n"
                "  Need to add larger sized particles or refine the size"
                " grid on the large size edge!\n")
          }
        }
      }
    }
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void pnd_fieldCalcFromParticleBulkProps(
    Workspace& ws,
    Tensor4& pnd_field,
    ArrayOfTensor4& dpnd_field_dx,
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const Tensor3& t_field,
    const Index& cloudbox_on,
    const ArrayOfIndex& cloudbox_limits,
    const ArrayOfString& scat_species,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const ArrayOfArrayOfScatteringMetaData& scat_meta,
    const Tensor4& particle_bulkprop_field,
    const ArrayOfString& particle_bulkprop_names,
    const ArrayOfAgenda& pnd_agenda_array,
    const ArrayOfArrayOfString& pnd_agenda_array_input_names,
    const Index& jacobian_do,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Verbosity&) {

  // Do nothing if cloudbox is inactive
  if (!cloudbox_on) {
    return;
  }

  ARTS_USER_ERROR_IF (particle_bulkprop_field.empty(),
                      "*particle_bulkprop_field* is empty.");

  // Checks (not totally complete, but should cover most mistakes)
  chk_if_in_range("atmosphere_dim", atmosphere_dim, 1, 3);
  chk_atm_grids(atmosphere_dim, p_grid, lat_grid, lon_grid);
  chk_atm_field("t_field", t_field, atmosphere_dim, p_grid, lat_grid, lon_grid);
  chk_atm_field("particle_bulkprop_field",
                particle_bulkprop_field,
                atmosphere_dim,
                particle_bulkprop_names.nelem(),
                p_grid,
                lat_grid,
                lon_grid);

  // Number of scattering species
  const Index nss = scat_data.nelem();

  ARTS_USER_ERROR_IF (particle_bulkprop_names.nelem() == 0 ||
      particle_bulkprop_names.nelem() != particle_bulkprop_field.nbooks(),
        "Number of fields in *particle_bulkprop_field*"
        " inconsistent with number of names in"
        " *particle_bulkprop_names*.");
  ARTS_USER_ERROR_IF (particle_bulkprop_field.nbooks() < nss,
        "At least one field per scattering species required"
        " in *particle_bulkprop_field*.");

  ARTS_USER_ERROR_IF (cloudbox_limits.nelem() != 2 * atmosphere_dim,
        "Length of *cloudbox_limits* incorrect with respect "
        "to *atmosphere_dim*.");
  ARTS_USER_ERROR_IF (cloudbox_limits[1] <= cloudbox_limits[0] || cloudbox_limits[0] < 0 ||
      cloudbox_limits[1] >= p_grid.nelem(),
                      "Invalid data in pressure part of *cloudbox_limits*.");
  if (atmosphere_dim > 1) {
    ARTS_USER_ERROR_IF (cloudbox_limits[3] <= cloudbox_limits[2] || cloudbox_limits[2] < 0 ||
        cloudbox_limits[3] >= lat_grid.nelem(),
          "Invalid data in latitude part of *cloudbox_limits*.");
    if (atmosphere_dim > 2) {
      ARTS_USER_ERROR_IF (cloudbox_limits[5] <= cloudbox_limits[4] || cloudbox_limits[4] < 0 ||
          cloudbox_limits[5] >= lon_grid.nelem(),
            "Invalid data in longitude part of *cloudbox_limits*.");
    }
  }

  ARTS_USER_ERROR_IF (nss < 1, "*scat_data* is empty!.");
  ARTS_USER_ERROR_IF (scat_species.nelem() != nss,
        "*scat_data* and *scat_species* are inconsistent in size.");
  ARTS_USER_ERROR_IF (scat_meta.nelem() != nss,
        "*scat_data* and *scat_meta* are inconsistent in size.");
  ARTS_USER_ERROR_IF (pnd_agenda_array.nelem() != nss,
        "*scat_data* and *pnd_agenda_array* are inconsistent "
        "in size.");
  ARTS_USER_ERROR_IF (pnd_agenda_array_input_names.nelem() != nss,
        "*scat_data* and *pnd_agenda_array_input_names* are "
        "inconsistent in size.");
  // Further checks of scat_data vs. scat_meta below

  // Effective lengths of cloudbox
  const Index np = cloudbox_limits[1] - cloudbox_limits[0] + 1;
  Index np_nonzero = np;                 // Can be changed below
  Index ip_offset = cloudbox_limits[0];  // Can be changed below
  Index pf_offset = 0;                   // Can be changed below
  Index nlat = 1;
  Index ilat_offset = 0;
  if (atmosphere_dim > 1) {
    nlat = cloudbox_limits[3] - cloudbox_limits[2] + 1;
    ilat_offset = cloudbox_limits[2];
  }
  Index nlon = 1;
  Index ilon_offset = 0;
  if (atmosphere_dim > 2) {
    nlon = cloudbox_limits[5] - cloudbox_limits[4] + 1;
    ilon_offset = cloudbox_limits[4];
  }

  // Check that *particle_bulkprop_field* contains zeros outside and at
  // cloudbox boundaries
  constexpr std::string_view estring =
      "*particle_bulkprop_field* allowed to contain"
      " non-zero values only inside the cloudbox.";
  // Pressure end ranges
  if (cloudbox_limits[0] != 0) {
    ARTS_USER_ERROR_IF (max(particle_bulkprop_field(
            joker, Range(0, ip_offset + 1), joker, joker)) > 0 ||
        min(particle_bulkprop_field(
            joker, Range(0, ip_offset + 1), joker, joker)) < 0,
                        estring);
    np_nonzero--;
    ip_offset++;
    pf_offset = 1;
  }
  if (cloudbox_limits[1] != p_grid.nelem() - 1) {
    const Index np_above = p_grid.nelem() + 1 - (np + ip_offset);
    ARTS_USER_ERROR_IF (max(particle_bulkprop_field(
            joker, Range(cloudbox_limits[1], np_above), joker, joker)) > 0 ||
        min(particle_bulkprop_field(
            joker, Range(cloudbox_limits[1], np_above), joker, joker)) < 0,
                        estring);
    np_nonzero--;
  }

  // Cumulative number of scattering elements
  ArrayOfIndex ncumse(nss + 1);
  ncumse[0] = 0;
  for (Index i = 0; i < nss; i++) {
    ARTS_USER_ERROR_IF (scat_data[i].nelem() != scat_meta[i].nelem(),
          "*scat_data* and *scat_meta* have inconsistent sizes.");
    ncumse[i + 1] = ncumse[i] + scat_data[i].nelem();
  }

  // Allocate output variables
  //
  pnd_field.resize(ncumse[nss], np, nlat, nlon);
  pnd_field = 0.0;
  //
  // Help variables for partial derivatives
  Index nq = 0;
  ArrayOfArrayOfIndex scatspecies_to_jq;
  //
  if (!jacobian_do) {
    dpnd_field_dx.resize(0);
  } else {
    nq = jacobian_quantities.nelem();
    dpnd_field_dx.resize(nq);
    scatspecies_to_jq.resize(nss);
    //
    for (Index iq = 0; iq < nq; iq++) {
      if (jacobian_quantities[iq] == Jacobian::Special::ScatteringString) {
        const Index ihit =
            find_first(scat_species, jacobian_quantities[iq].Subtag());
        ARTS_USER_ERROR_IF (ihit < 0,
            "Jacobian quantity with index ", iq, " refers to\n"
            "  ", jacobian_quantities[iq].Subtag(),
            "\nbut this species could not be found in *scat_species*.")
        scatspecies_to_jq[ihit].push_back(iq);
        dpnd_field_dx[iq].resize(ncumse[nss], np, nlat, nlon);
        dpnd_field_dx[iq] = 0.0;
      }
    }
  }

  // Extract data from pnd-agenda array
  for (Index is = 0; is < nss; is++) {
    // Index range with respect to pnd_field
    Range se_range(ncumse[is], ncumse[is + 1] - ncumse[is]);

    // Determine how pnd_agenda_array_input_names are related to input fields
    //
    const Index nin = pnd_agenda_array_input_names[is].nelem();
    ArrayOfIndex i_pbulkprop(nin);
    //
    for (Index i = 0; i < nin; i++) {
      i_pbulkprop[i] = find_first(particle_bulkprop_names,
                                  pnd_agenda_array_input_names[is][i]);
      ARTS_USER_ERROR_IF (i_pbulkprop[i] < 0,
          "Pnd-agenda with index ", is, " is set to require \"",
           pnd_agenda_array_input_names[is][i], "\",\nbut this quantity "
           "could not found in *particle_bulkprop_names*.\n"
           "(Note that temperature must be written as \"Temperature\")")
    }

    // Set *dpnd_data_dx_names*
    //
    Index ndx = 0;
    ArrayOfString dpnd_data_dx_names(0);
    //
    if (jacobian_do) {
      ndx = scatspecies_to_jq[is].nelem();
      dpnd_data_dx_names.resize(ndx);
      for (Index ix = 0; ix < ndx; ix++) {
        dpnd_data_dx_names[ix] =
            jacobian_quantities[scatspecies_to_jq[is][ix]].SubSubtag();
      }
    }

    // Loop lat/lon positions and call *pnd_agenda*
    for (Index ilon = 0; ilon < nlon; ilon++) {
      for (Index ilat = 0; ilat < nlat; ilat++) {

        // Nothing to do for lat/lon end points, if not 1D
        if ((nlat > 1 && (ilat == 0 || ilat == nlat - 1)) ||
            (nlon > 1 && (ilon == 0 || ilon == nlon - 1))) {
          continue;
        }

        Matrix pnd_agenda_input(np_nonzero, nin);
        //
        for (Index i = 0; i < nin; i++) {
          for (Index ip = 0; ip < np_nonzero; ip++) {
            pnd_agenda_input(ip, i) =
                particle_bulkprop_field(i_pbulkprop[i],
                                        ip_offset + ip,
                                        ilat_offset + ilat,
                                        ilon_offset + ilon);
          }
        }

        Vector pnd_agenda_input_t(np);
        //
        for (Index ip = 0; ip < np_nonzero; ip++) {
          pnd_agenda_input_t[ip] =
              t_field(ip_offset + ip, ilat_offset + ilat, ilon_offset + ilon);
        }

        // Call pnd-agenda array
        Matrix pnd_data;
        Tensor3 dpnd_data_dx;
        //
        pnd_agenda_arrayExecute(ws,
                                pnd_data,
                                dpnd_data_dx,
                                is,
                                pnd_agenda_input_t,
                                pnd_agenda_input,
                                pnd_agenda_array_input_names[is],
                                dpnd_data_dx_names,
                                pnd_agenda_array);

        // Copy to output variables
        for (Index ip = 0; ip < np_nonzero; ip++) {
          pnd_field(se_range, pf_offset+ip, ilat, ilon) = pnd_data(ip, joker);
        }
        for (Index ix = 0; ix < ndx; ix++) {
          for (Index ip = 0; ip < np_nonzero; ip++) {
            dpnd_field_dx[scatspecies_to_jq[is][ix]](se_range, pf_offset+ip, ilat, ilon) =
                dpnd_data_dx(ix, ip, joker);
          }
        }
      }
    }
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void ScatSpeciesSizeMassInfo(Vector& scat_species_x,
                             Numeric& scat_species_a,
                             Numeric& scat_species_b,
                             const ArrayOfArrayOfScatteringMetaData& scat_meta,
                             const Index& species_index,
                             const String& x_unit,
                             const Numeric& x_fit_start,
                             const Numeric& x_fit_end,
                             const Index& do_only_x,
                             const Verbosity&) {
  // Checks
  const Index nss = scat_meta.nelem();
  ARTS_USER_ERROR_IF (nss == 0, "*scat_meta* is empty!");
  ARTS_USER_ERROR_IF (nss < species_index + 1,
      "Selected scattering species index is ", species_index,
      " but this "
      "is not allowed since *scat_meta* has only ", nss, " elements.")
  //
  const Index nse = scat_meta[species_index].nelem();
  ARTS_USER_ERROR_IF (nse < 2,
        "The scattering species must have at least two "
        "elements to use this method.");

  // Extract particle masses
  //
  Vector mass(nse);
  //
  for (Index i = 0; i < nse; i++) {
    mass[i] = scat_meta[species_index][i].mass;
  }

  // Create size grid
  //
  scat_species_x.resize(nse);
  //
  if (x_unit == "dveq") {
    for (Index i = 0; i < nse; i++) {
      scat_species_x[i] = scat_meta[species_index][i].diameter_volume_equ;
    }
    //
    if (do_only_x) {
      scat_species_a = -1;
      scat_species_b = -1;
    } else
      derive_scat_species_a_and_b(scat_species_a,
                                  scat_species_b,
                                  scat_species_x,
                                  mass,
                                  x_fit_start,
                                  x_fit_end);
  }

  else if (x_unit == "dmax") {
    for (Index i = 0; i < nse; i++) {
      scat_species_x[i] = scat_meta[species_index][i].diameter_max;
    }
    //
    if (do_only_x) {
      scat_species_a = -1;
      scat_species_b = -1;
    } else
      derive_scat_species_a_and_b(scat_species_a,
                                  scat_species_b,
                                  scat_species_x,
                                  mass,
                                  x_fit_start,
                                  x_fit_end);
  }

  else if (x_unit == "area") {
    for (Index i = 0; i < nse; i++) {
      scat_species_x[i] =
          scat_meta[species_index][i].diameter_area_equ_aerodynamical;
    }
    //
    if (do_only_x) {
      scat_species_a = -1;
      scat_species_b = -1;
    } else
      derive_scat_species_a_and_b(scat_species_a,
                                  scat_species_b,
                                  scat_species_x,
                                  mass,
                                  x_fit_start,
                                  x_fit_end);
  }

  else if (x_unit == "mass") {
    scat_species_x = mass;
    //
    scat_species_a = 1;
    scat_species_b = 1;
  }

  else {
    ARTS_USER_ERROR ("You have selected the x_unit: ", x_unit, "\nwhile accepted "
                     "choices are: \"dveq\", \"dmax\", \"mass\" and \"area\"")
  }
}
