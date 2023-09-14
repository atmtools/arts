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
#include "arts_constants.h"
#include "atm.h"
#include "check_input.h"
#include "cloudbox.h"
#include "disort.h"
#include "file.h"
#include "interpolation.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "microphysics.h"
#include "optproperties.h"
#include "parameters.h"
#include "physics_funcs.h"
#include "psd.h"
#include "rte.h"
#include "sorting.h"
#include "special_interp.h"
#include "xml_io.h"

#include <workspace.h>

inline constexpr Numeric PI=Constant::pi;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void HydrotableCalc(const Workspace& ws,
                    GriddedField4& hydrotable,
                    const ArrayOfAgenda& pnd_agenda_array,
                    const ArrayOfArrayOfString& pnd_agenda_array_input_names,
                    const ArrayOfArrayOfSingleScatteringData& scat_data,
                    const Index& scat_data_checked,
                    const Vector& f_grid,
                    const Index& iss,
                    const Vector& T_grid,
                    const Vector& wc_grid)
{
  // Sizes
  const Index nss = scat_data.size(); 
  const Index nf = f_grid.nelem();
  const Index nt = T_grid.nelem();
  const Index nw = wc_grid.nelem();

  ARTS_USER_ERROR_IF (pnd_agenda_array.size() != static_cast<Size>(nss),
        "*scat_data* and *pnd_agenda_array* are inconsistent "
        "in size.");
  ARTS_USER_ERROR_IF (pnd_agenda_array_input_names.size() != static_cast<Size>(nss),
        "*scat_data* and *pnd_agenda_array_input_names* are "
        "inconsistent in size.");
  ARTS_USER_ERROR_IF (pnd_agenda_array_input_names[iss].size() != 1,
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
    const ArrayOfArrayOfScatteringMetaData& scat_meta) {
  const Index np_total = TotalNumberOfElements(scat_meta);

  particle_masses.resize(np_total, 1);

  Index i_se_flat = 0;
  for (Size i_ss = 0; i_ss < scat_meta.size(); i_ss++) {
    for (Size i_se = 0; i_se < scat_meta[i_ss].size(); i_se++) {
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
    const ArrayOfArrayOfScatteringMetaData& scat_meta) {
  // resize particle_masses to required diemsions and properly initialize values
  particle_masses.resize(TotalNumberOfElements(scat_meta), scat_meta.size());
  particle_masses = 0.;

  // calculate and set particle_masses
  Index i_se_flat = 0;
  for (Size i_ss = 0; i_ss < scat_meta.size(); i_ss++) {
    for (Size i_se = 0; i_se < scat_meta[i_ss].size(); i_se++) {
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
                     const Index& quad_order) {
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
                const Numeric& threshold_rpnd) {
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
  ARTS_USER_ERROR_IF (static_cast<Size>(scat_index) >= scat_data.size(),
        "*scat_index* exceeds the number of available"
        " scattering species.");
  ARTS_USER_ERROR_IF (scat_data[scat_index].size() != static_cast<Size>(ng),
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
    max0 = std::max(std::abs(pnd_data(ip, intarr[0])), max0);
    max1 = std::max(std::abs(pnd_data(ip, intarr[ng - 1])), max1);
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
void ScatSpeciesSizeMassInfo(Vector& scat_species_x,
                             Numeric& scat_species_a,
                             Numeric& scat_species_b,
                             const ArrayOfArrayOfScatteringMetaData& scat_meta,
                             const Index& species_index,
                             const String& x_unit,
                             const Numeric& x_fit_start,
                             const Numeric& x_fit_end,
                             const Index& do_only_x) {
  // Checks
  const Index nss = scat_meta.size();
  ARTS_USER_ERROR_IF (nss == 0, "*scat_meta* is empty!");
  ARTS_USER_ERROR_IF (nss < species_index + 1,
      "Selected scattering species index is ", species_index,
      " but this "
      "is not allowed since *scat_meta* has only ", nss, " elements.")
  //
  const Index nse = scat_meta[species_index].size();
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
