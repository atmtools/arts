
/* Copyright (C) 2000-2012
   Stefan Buehler   <sbuehler@ltu.se>
   Patrick Eriksson <patrick.eriksson@chalmers.se>
   Axel von Engeln  <engeln@uni-bremen.de>
   Thomas Kuhn      <tkuhn@uni-bremen.de>

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

//

/**
   \file   m_abs.cc

   Stuff related to the calculation of absorption coefficients.

   \author Stefan Buehler
   \date   2001-03-12
*/
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>

#include "absorption.h"
#include "absorptionlines.h"
#include "agenda_class.h"
#include "agenda_set.h"
#include "array.h"
#include "arts.h"
#include "arts_constants.h"
#include "arts_omp.h"
#include "artstime.h"
#include "auto_md.h"
#include "check_input.h"
#include "debug.h"
#include "depr.h"
#include "file.h"
#include "global_data.h"
#include "hitran_species.h"
#include "jacobian.h"
#include "lineshape.h"
#include "m_xml.h"
#include "math_funcs.h"
#include "matpack_data.h"
#include "messages.h"
#include "methods.h"
#include "montecarlo.h"
#include "optproperties.h"
#include "parameters.h"
#include "physics_funcs.h"
#include "rte.h"
#include "species_tags.h"
#include "xml_io.h"

#ifdef ENABLE_NETCDF
#include <netcdf.h>

#include "nc_io.h"
#endif

inline constexpr Numeric ELECTRON_CHARGE=-Constant::elementary_charge;
inline constexpr Numeric ELECTRON_MASS=Constant::electron_mass;
inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;
inline constexpr Numeric VACUUM_PERMITTIVITY=Constant::vacuum_permittivity;

/* Workspace method: Doxygen documentation will be auto-generated */
void AbsInputFromRteScalars(  // WS Output:
    Vector& abs_p,
    Vector& abs_t,
    Matrix& abs_vmrs,
    // WS Input:
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const Vector& rtp_vmr,
    const Verbosity&) {
  // Prepare abs_p:
  abs_p.resize(1);
  abs_p = rtp_pressure;

  // Prepare abs_t:
  abs_t.resize(1);
  abs_t = rtp_temperature;

  // Prepare abs_vmrs:
  abs_vmrs.resize(rtp_vmr.nelem(), 1);
  abs_vmrs = ExhaustiveMatrixView{rtp_vmr};
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesCreateFromLines(  // WS Output:
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    // WS Input:
    const ArrayOfAbsorptionLines& abs_lines,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Verbosity&) {
  // Size is set but inner size will now change from the original definition of species tags...
  abs_lines_per_species.resize(abs_species.nelem());

  // The inner arrays need to be emptied, because they may contain lines
  // from a previous calculation
  for (auto& lines : abs_lines_per_species) lines.resize(0);

#pragma omp parallel for schedule(dynamic) if (!arts_omp_in_parallel())
  for (Index ilines = 0; ilines < abs_lines.nelem(); ilines++) {
    AbsorptionLines lines = abs_lines[ilines];
    
    // Skip empty lines
    if (lines.NumLines() == 0) continue;

    // Loop all the tags
    for (Index i = 0; i < abs_species.nelem() and lines.NumLines(); i++) {
      for (auto& this_tag : abs_species[i]) {
        // Test isotopologue, we have to hit the end of the list for no isotopologue or the exact value
        if (not same_or_joker(this_tag.Isotopologue(), lines.Isotopologue()))
          continue;

        // If there is a frequency range, we have to check so that only selected lines are included
        if (this_tag.lower_freq >= 0 or this_tag.upper_freq >= 0) {
          const Numeric low = (this_tag.lower_freq >= 0)
                                  ? this_tag.lower_freq
                                  : std::numeric_limits<Numeric>::lowest();
          const Numeric upp = (this_tag.upper_freq >= 0)
                                  ? this_tag.upper_freq
                                  : std::numeric_limits<Numeric>::max();

          // Fill up a copy of the line record to match with the wished frequency criteria
          AbsorptionLines these_lines = lines;
          these_lines.lines.resize(0);
          for (Index k = lines.NumLines() - 1; k >= 0; k--)
            if (low <= lines.lines[k].F0 and upp >= lines.lines[k].F0)
              these_lines.AppendSingleLine(lines.PopLine(k));

          // Append these lines after sorting them if there are any of them
          if (these_lines.NumLines()) {
            these_lines.ReverseLines();
#pragma omp critical
            abs_lines_per_species[i].push_back(these_lines);
          }

          // If this means we have deleted all lines, then we leave
          if (lines.NumLines() == 0) goto leave_inner_loop;
        } else {
#pragma omp critical
          abs_lines_per_species[i].push_back(lines);
          goto leave_inner_loop;
        }
      }
    }
  leave_inner_loop : {}
  }

  abs_lines_per_species.shrink_to_fit();
  for (auto& spec_band : abs_lines_per_species)
    std::sort(spec_band.begin(), spec_band.end(), [](auto& a, auto& b) {
      return a.lines.size() and b.lines.size() and
             a.lines.front().F0 < b.lines.front().F0;
    });
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesDefineAllInScenario(  // WS Output:
    ArrayOfArrayOfSpeciesTag& tgs,
    Index& propmat_clearsky_agenda_checked,
    // Control Parameters:
    const String& basename,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  // Invalidate agenda check flags
  propmat_clearsky_agenda_checked = false;

  // We want to make lists of included and excluded species:
  ArrayOfString included(0), excluded(0);

  tgs.resize(0);

  for (Index i = 0; i < Index(Species::Species::FINAL); ++i) {
    const String specname = Species::toShortName(Species::Species(i));

    String filename = basename;
    if (basename.length() && basename[basename.length() - 1] != '/')
      filename += ".";
    filename += specname;

    try {
      find_xml_file(filename, verbosity);
      // Add to included list:
      included.push_back(specname);

      // Add this tag group to tgs:
      tgs.emplace_back(ArrayOfSpeciesTag(specname));
    } catch (const std::runtime_error& e) {
      // The file for the species could not be found.
      excluded.push_back(specname);
    }
  }

  // Some nice output:
  out2 << "  Included Species (" << included.nelem() << "):\n";
  for (Index i = 0; i < included.nelem(); ++i)
    out2 << "     " << included[i] << "\n";

  out2 << "  Excluded Species (" << excluded.nelem() << "):\n";
  for (Index i = 0; i < excluded.nelem(); ++i)
    out2 << "     " << excluded[i] << "\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesDefineAll(  // WS Output:
    ArrayOfArrayOfSpeciesTag& abs_species,
    Index& propmat_clearsky_agenda_checked,
    // Control Parameters:
    const Verbosity& verbosity) {
  // Species lookup data:

  // We want to make lists of all species
  ArrayOfString specs(0);
  for (Index i = 0; i < Index(Species::Species::FINAL); ++i) {
    if (Species::Species(i) not_eq Species::Species::Bath) {
      specs.emplace_back(Species::toShortName(Species::Species(i)));
    }
  }

  // Set the values
  abs_speciesSet(abs_species,
                 propmat_clearsky_agenda_checked,
                 specs,
                 verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void AbsInputFromAtmFields(  // WS Output:
    Vector& abs_p,
    Vector& abs_t,
    Matrix& abs_vmrs,
    // WS Input:
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Tensor3& t_field,
    const Tensor4& vmr_field,
    const Verbosity&) {
  // First, make sure that we really have a 1D atmosphere:
  ARTS_USER_ERROR_IF(1 != atmosphere_dim,
                     "Atmospheric dimension must be 1D, but atmosphere_dim is ",
                     atmosphere_dim,
                     ".")

  abs_p = p_grid;
  abs_t = t_field(joker, 0, 0);
  abs_vmrs = vmr_field(joker, joker, 0, 0);
}

//======================================================================
//             Methods related to continua
//======================================================================

/* Workspace method: Doxygen documentation will be auto-generated */
void nlte_sourceFromTemperatureAndSrcCoefPerSpecies(  // WS Output:
    StokesVector& nlte_source,
    ArrayOfStokesVector& dnlte_source_dx,
    // WS Input:
    const ArrayOfMatrix& src_coef_per_species,
    const ArrayOfMatrix& dsrc_coef_dx,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& f_grid,
    const Numeric& rtp_temperature,
    const Verbosity&) {
  // nlte_source has format
  // [ abs_species, f_grid, stokes_dim ].
  // src_coef_per_species has format ArrayOfMatrix (over species),
  // where for each species the matrix has format [f_grid, abs_p].

  Index n_species = src_coef_per_species.nelem();  // # species

  ARTS_USER_ERROR_IF(not n_species, "Must have at least one species.")

  Index n_f = src_coef_per_species[0].nrows();  // # frequencies

  // # pressures must be 1:
  ARTS_USER_ERROR_IF(1 not_eq src_coef_per_species[0].ncols(),
                     "Must have exactly one pressure.")

  // Check frequency dimension of propmat_clearsky
  ARTS_USER_ERROR_IF(nlte_source.NumberOfFrequencies() not_eq n_f,
                     "Frequency dimension of nlte_source does not\n"
                     "match abs_coef_per_species.")

  const Vector B = planck(f_grid, rtp_temperature);

  StokesVector sv(n_f, nlte_source.StokesDimensions());
  for (Index si = 0; si < n_species; ++si) {
    sv.Kjj() = src_coef_per_species[si](joker, 0);
    sv *= B;
    nlte_source.Kjj() += sv.Kjj();
  }

  // Jacobian
  for (Index ii = 0; ii < jacobian_quantities.nelem(); ii++) {
    const auto& deriv = jacobian_quantities[ii];

    if (not deriv.propmattype()) continue;

    if (deriv == Jacobian::Atm::Temperature) {
      const Vector dB = dplanck_dt(f_grid, rtp_temperature);

      for (Index si = 0; si < n_species; ++si) {
        sv.Kjj() = src_coef_per_species[si](joker, 0);
        sv *= dB;
        dnlte_source_dx[ii].Kjj() += sv.Kjj();
      }

      sv.Kjj() = dsrc_coef_dx[ii](joker, 0);
      sv *= B;
      dnlte_source_dx[ii].Kjj() += sv.Kjj();
    } else if (is_frequency_parameter(deriv)) {
      const Vector dB = dplanck_df(f_grid, rtp_temperature);

      for (Index si = 0; si < n_species; ++si) {
        sv.Kjj() = src_coef_per_species[si](joker, 0);
        sv *= dB;
        dnlte_source_dx[ii].Kjj() += sv.Kjj();
      }

      sv.Kjj() = dsrc_coef_dx[ii](joker, 0);
      sv *= B;
      dnlte_source_dx[ii].Kjj() += sv.Kjj();
    } else {
      sv.Kjj() = dsrc_coef_dx[ii](joker, 0);
      sv *= B;
      dnlte_source_dx[ii].Kjj() += sv.Kjj();
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyInit(  //WS Output
    PropagationMatrix& propmat_clearsky,
    StokesVector& nlte_source,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    ArrayOfStokesVector& dnlte_source_dx,
    //WS Input
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& f_grid,
    const Index& stokes_dim,
    const Index& propmat_clearsky_agenda_checked,
    const Verbosity&) {
  const Index nf = f_grid.nelem();
  const Index nq = jacobian_quantities.nelem();

  ARTS_USER_ERROR_IF(
      !propmat_clearsky_agenda_checked,
      "You must call *propmat_clearsky_agenda_checkedCalc* before calling this method.")

  ARTS_USER_ERROR_IF(not nf, "No frequencies");

  ARTS_USER_ERROR_IF(stokes_dim < 1 or stokes_dim > 4,
                     "stokes_dim not in [1, 2, 3, 4]");

  // Set size of propmat_clearsky or reset it's values
  if (propmat_clearsky.StokesDimensions() == stokes_dim and
      propmat_clearsky.NumberOfFrequencies() == nf) {
    propmat_clearsky.SetZero();
  } else {
    propmat_clearsky = PropagationMatrix(nf, stokes_dim);
  }

  // Set size of dpropmat_clearsky_dx or reset it's values
  if (dpropmat_clearsky_dx.nelem() not_eq nq) {
    dpropmat_clearsky_dx =
        ArrayOfPropagationMatrix(nq, PropagationMatrix(nf, stokes_dim));
  } else {
    for (auto& pm : dpropmat_clearsky_dx) {
      if (pm.StokesDimensions() == stokes_dim and
          pm.NumberOfFrequencies() == nf) {
        pm.SetZero();
      } else {
        pm = PropagationMatrix(nf, stokes_dim);
      }
    }
  }

  // Set size of nlte_source or reset it's values
  if (nlte_source.StokesDimensions() == stokes_dim and
      nlte_source.NumberOfFrequencies() == nf) {
    nlte_source.SetZero();
  } else {
    nlte_source = StokesVector(nf, stokes_dim);
  }

  // Set size of dnlte_source_dx or reset it's values
  if (dnlte_source_dx.nelem() not_eq nq) {
    dnlte_source_dx = ArrayOfStokesVector(nq, StokesVector(nf, stokes_dim));
  } else {
    for (auto& pm : dnlte_source_dx) {
      if (pm.StokesDimensions() == stokes_dim and
          pm.NumberOfFrequencies() == nf) {
        pm.SetZero();
      } else {
        pm = StokesVector(nf, stokes_dim);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddFaraday(
    PropagationMatrix& propmat_clearsky,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    const Index& stokes_dim,
    const Index& atmosphere_dim,
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& rtp_vmr,
    const Vector& rtp_los,
    const Vector& rtp_mag,
    const Verbosity&) {
  Index ife = -1;
  for (Index sp = 0; sp < abs_species.nelem() && ife < 0; sp++) {
    if (abs_species[sp].FreeElectrons()) {
      ife = sp;
    }
  }

  ARTS_USER_ERROR_IF(ife < 0,
                     "Free electrons not found in *abs_species* and "
                     "Faraday rotation can not be calculated.");

  // Allow early exit for lookup table calculations
  if (select_abs_species.nelem() and select_abs_species not_eq abs_species[ife]) return;

  // All the physical constants joined into one static constant:
  // (abs as e defined as negative)
  static const Numeric FRconst =
      abs(ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE /
          (8 * PI * PI * SPEED_OF_LIGHT * VACUUM_PERMITTIVITY * ELECTRON_MASS *
           ELECTRON_MASS));

  const bool do_magn_jac = do_magnetic_jacobian(jacobian_quantities);
  const Numeric dmag = magnetic_field_perturbation(jacobian_quantities);

  ARTS_USER_ERROR_IF(
      stokes_dim < 3,
      "To include Faraday rotation, stokes_dim >= 3 is required.")
  ARTS_USER_ERROR_IF(
      atmosphere_dim == 1 && rtp_los.nelem() < 1,
      "For applying propmat_clearskyAddFaraday, los needs to be specified\n"
      "(at least zenith angle component for atmosphere_dim==1),\n"
      "but it is not.\n")
  ARTS_USER_ERROR_IF(
      atmosphere_dim > 1 && rtp_los.nelem() < 2,
      "For applying propmat_clearskyAddFaraday, los needs to be specified\n"
      "(both zenith and azimuth angle components for atmosphere_dim>1),\n"
      "but it is not.\n")

  const Numeric ne = rtp_vmr[ife];

  if (ne != 0 && (rtp_mag[0] != 0 || rtp_mag[1] != 0 || rtp_mag[2] != 0)) {
    // Include remaining terms, beside /f^2
    const Numeric c1 =
        2 * FRconst *
        dotprod_with_los(
            rtp_los, rtp_mag[0], rtp_mag[1], rtp_mag[2], atmosphere_dim);

    Numeric dc1_u = 0.0, dc1_v = 0.0, dc1_w = 0.0;
    if (do_magn_jac) {
      dc1_u = (2 * FRconst *
                   dotprod_with_los(rtp_los,
                                    rtp_mag[0] + dmag,
                                    rtp_mag[1],
                                    rtp_mag[2],
                                    atmosphere_dim) -
               c1) /
              dmag;
      dc1_v = (2 * FRconst *
                   dotprod_with_los(rtp_los,
                                    rtp_mag[0],
                                    rtp_mag[1] + dmag,
                                    rtp_mag[2],
                                    atmosphere_dim) -
               c1) /
              dmag;
      dc1_w = (2 * FRconst *
                   dotprod_with_los(rtp_los,
                                    rtp_mag[0],
                                    rtp_mag[1],
                                    rtp_mag[2] + dmag,
                                    atmosphere_dim) -
               c1) /
              dmag;
    }

    for (Index iv = 0; iv < f_grid.nelem(); iv++) {
      const Numeric f2 = f_grid[iv] * f_grid[iv];
      const Numeric r = ne * c1 / f2;
      propmat_clearsky.AddFaraday(r, iv);

      // The Jacobian loop
      for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
        if (is_frequency_parameter(jacobian_quantities[iq]))
          dpropmat_clearsky_dx[iq].AddFaraday(-2.0 * ne * r / f_grid[iv], iv);
        else if (jacobian_quantities[iq] == Jacobian::Atm::MagneticU)
          dpropmat_clearsky_dx[iq].AddFaraday(ne * dc1_u / f2, iv);
        else if (jacobian_quantities[iq] == Jacobian::Atm::MagneticV)
          dpropmat_clearsky_dx[iq].AddFaraday(ne * dc1_v / f2, iv);
        else if (jacobian_quantities[iq] == Jacobian::Atm::MagneticW)
          dpropmat_clearsky_dx[iq].AddFaraday(ne * dc1_w / f2, iv);
        else if (jacobian_quantities[iq] == Jacobian::Atm::Electrons)
          dpropmat_clearsky_dx[iq].AddFaraday(r, iv);
        else if (jacobian_quantities[iq] == abs_species[ife])
          dpropmat_clearsky_dx[iq].AddFaraday(r, iv);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddParticles(
    // WS Output:
    PropagationMatrix& propmat_clearsky,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    // WS Input:
    const Index& stokes_dim,
    const Index& atmosphere_dim,
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& rtp_vmr,
    const Vector& rtp_los,
    const Numeric& rtp_temperature,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Index& scat_data_checked,
    const Index& use_abs_as_ext,
    // Verbosity object:
    const Verbosity& verbosity) {
  CREATE_OUT1;

  ARTS_USER_ERROR_IF(select_abs_species.nelem(), R"--(
  We do not yet support select_abs_species for lookup table calculations
  )--")

  // (i)yCalc only checks scat_data_checked if cloudbox is on. It is off here,
  // though, i.e. we need to check it here explicitly. (Also, cloudboxOff sets
  // scat_data_checked=0 as it does not check it and as we ususally don't need
  // scat_data for clearsky cases, hence don't want to check them by
  // scat_data_checkedCalc in that case. This approach seems to be the more
  // handy compared to cloudboxOff setting scat_data_checked=1 without checking
  // it assuming we won't use it anyways.)
  ARTS_USER_ERROR_IF(scat_data_checked != 1,
                     "The scat_data must be flagged to have "
                     "passed a consistency check (scat_data_checked=1).")

  const Index ns = TotalNumberOfElements(scat_data);
  Index np = 0;
  for (Index sp = 0; sp < abs_species.nelem(); sp++) {
    if (abs_species[sp].Particles()) {
      np++;
    }
  }

  ARTS_USER_ERROR_IF(
      np == 0,
      "For applying propmat_clearskyAddParticles, *abs_species* needs to"
      "contain species 'particles', but it does not.\n")

  ARTS_USER_ERROR_IF(
      ns != np,
      "Number of 'particles' entries in abs_species and of elements in\n"
      "*scat_data* needs to be identical. But you have ",
      np,
      " 'particles' entries\n"
      "and ",
      ns,
      " *scat_data* elements.\n")

  ARTS_USER_ERROR_IF(
      atmosphere_dim == 1 && rtp_los.nelem() < 1,
      "For applying *propmat_clearskyAddParticles*, *rtp_los* needs to be specified\n"
      "(at least zenith angle component for atmosphere_dim==1),\n"
      "but it is not.\n")
  ARTS_USER_ERROR_IF(
      atmosphere_dim > 1 && rtp_los.nelem() < 2,
      "For applying *propmat_clearskyAddParticles*, *rtp_los* needs to be specified\n"
      "(both zenith and azimuth angle components for atmosphere_dim>1),\n"
      "but it is not.\n")

  // Use for rescaling vmr of particulates
  Numeric rtp_vmr_sum = 0.0;

  // Tests and setup partial derivatives
  const bool do_jac_temperature = do_temperature_jacobian(jacobian_quantities);
  const bool do_jac_frequencies = do_frequency_jacobian(jacobian_quantities);
  const Numeric dT = temperature_perturbation(jacobian_quantities);

  const Index na = abs_species.nelem();
  Vector rtp_los_back;
  mirror_los(rtp_los_back, rtp_los, atmosphere_dim);

  // 170918 JM: along with transition to use of new-type (aka
  // pre-f_grid-interpolated) scat_data, freq perturbation switched off. Typical
  // clear-sky freq perturbations yield insignificant effects in particle
  // properties. Hence, this feature is neglected here.
  if (do_jac_frequencies) {
    out1 << "WARNING:\n"
         << "Frequency perturbation not available for absorbing particles.\n";
  }

  // creating temporary output containers
  ArrayOfArrayOfTensor5 ext_mat_Nse;
  ArrayOfArrayOfTensor4 abs_vec_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;

  // preparing input in format needed
  Vector T_array;
  if (do_jac_temperature) {
    T_array.resize(2);
    T_array = rtp_temperature;
    T_array[1] += dT;
  } else {
    T_array.resize(1);
    T_array = rtp_temperature;
  }
  Matrix dir_array(1, 2);
  dir_array(0, joker) = rtp_los_back;

  // ext/abs per scat element for all freqs at once
  opt_prop_NScatElems(ext_mat_Nse,
                      abs_vec_Nse,
                      ptypes_Nse,
                      t_ok,
                      scat_data,
                      stokes_dim,
                      T_array,
                      dir_array,
                      -1);

  const Index nf = abs_vec_Nse[0][0].nbooks();
  Tensor3 tmp(nf, stokes_dim, stokes_dim);

  // Internal computations necessary since it relies on zero start
  PropagationMatrix internal_propmat(propmat_clearsky.NumberOfFrequencies(),
                                     propmat_clearsky.StokesDimensions());

  // loop over the scat_data and link them with correct vmr_field entry according
  // to the position of the particle type entries in abs_species.
  Index sp = 0;
  Index i_se_flat = 0;
  for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++) {
    for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++) {
      // forward to next particle entry in abs_species
      while (sp < na && not abs_species[sp].Particles()) sp++;
      internal_propmat.SetZero();

      // running beyond number of abs_species entries when looking for next
      // particle entry. shouldn't happen, though.
      ARTS_ASSERT(sp < na);
      ARTS_USER_ERROR_IF(
          rtp_vmr[sp] < 0.,
          "Negative absorbing particle 'vmr' (aka number density)"
          " encountered:\n"
          "scat species #",
          i_ss,
          ", scat elem #",
          i_se,
          " (vmr_field entry #",
          sp,
          ")\n")

      if (rtp_vmr[sp] > 0.) {
        ARTS_USER_ERROR_IF(t_ok(i_se_flat, 0) < 0.,
                           "Temperature interpolation error:\n"
                           "scat species #",
                           i_ss,
                           ", scat elem #",
                           i_se,
                           "\n")
        if (use_abs_as_ext) {
          if (nf > 1)
            for (Index iv = 0; iv < f_grid.nelem(); iv++)
              internal_propmat.AddAbsorptionVectorAtPosition(
                  abs_vec_Nse[i_ss][i_se](iv, 0, 0, joker), iv);
          else
            for (Index iv = 0; iv < f_grid.nelem(); iv++)
              internal_propmat.AddAbsorptionVectorAtPosition(
                  abs_vec_Nse[i_ss][i_se](0, 0, 0, joker), iv);
        } else {
          if (nf > 1)
            for (Index iv = 0; iv < f_grid.nelem(); iv++)
              internal_propmat.SetAtPosition(
                  ext_mat_Nse[i_ss][i_se](iv, 0, 0, joker, joker), iv);
          else
            for (Index iv = 0; iv < f_grid.nelem(); iv++)
              internal_propmat.SetAtPosition(
                  ext_mat_Nse[i_ss][i_se](0, 0, 0, joker, joker), iv);
        }
        propmat_clearsky += rtp_vmr[sp] * internal_propmat;
      }

      // For temperature derivatives (so we don't need to check it in jac loop)
      if (do_jac_temperature) {
        ARTS_USER_ERROR_IF(
            t_ok(i_se_flat, 1) < 0.,
            "Temperature interpolation error (in perturbation):\n"
            "scat species #",
            i_ss,
            ", scat elem #",
            i_se,
            "\n")
      }

      // For number density derivatives
      if (jacobian_quantities.nelem()) rtp_vmr_sum += rtp_vmr[sp];

      for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
        const auto& deriv = jacobian_quantities[iq];

        if (not deriv.propmattype()) continue;

        if (deriv == Jacobian::Atm::Temperature) {
          if (use_abs_as_ext) {
            tmp(joker, joker, 0) = abs_vec_Nse[i_ss][i_se](joker, 1, 0, joker);
            tmp(joker, joker, 0) -= abs_vec_Nse[i_ss][i_se](joker, 0, 0, joker);
          } else {
            tmp = ext_mat_Nse[i_ss][i_se](joker, 1, 0, joker, joker);
            tmp -= ext_mat_Nse[i_ss][i_se](joker, 0, 0, joker, joker);
          }

          tmp *= rtp_vmr[sp];
          tmp /= dT;

          if (nf > 1)
            for (Index iv = 0; iv < f_grid.nelem(); iv++)
              if (use_abs_as_ext)
                dpropmat_clearsky_dx[iq].AddAbsorptionVectorAtPosition(
                    tmp(iv, joker, 0), iv);
              else
                dpropmat_clearsky_dx[iq].AddAtPosition(tmp(iv, joker, joker),
                                                       iv);
          else
            for (Index iv = 0; iv < f_grid.nelem(); iv++)
              if (use_abs_as_ext)
                dpropmat_clearsky_dx[iq].AddAbsorptionVectorAtPosition(
                    tmp(0, joker, 0), iv);
              else
                dpropmat_clearsky_dx[iq].AddAtPosition(tmp(0, joker, joker),
                                                       iv);
        }

        else if (deriv == Jacobian::Atm::Particulates) {
          for (Index iv = 0; iv < f_grid.nelem(); iv++)
            dpropmat_clearsky_dx[iq].AddAtPosition(internal_propmat, iv);
        }

        else if (deriv == abs_species[sp]) {
          dpropmat_clearsky_dx[iq] += internal_propmat;
        }
      }
      sp++;
      i_se_flat++;
    }
  }

  //checking that no further 'particle' entry left after all scat_data entries
  //are processes. this is basically not necessary. but checking it anyway to
  //really be safe. remove later, when more extensively tested.
  while (sp < na) {
    ARTS_ASSERT(abs_species[sp][0].Type() != Species::TagType::Particles);
    sp++;
  }

  if (rtp_vmr_sum != 0.0) {
    for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
      const auto& deriv = jacobian_quantities[iq];

      if (not deriv.propmattype()) continue;

      if (deriv == Jacobian::Atm::Particulates) {
        dpropmat_clearsky_dx[iq] /= rtp_vmr_sum;
      }
    }
  }
}

void sparse_f_gridFromFrequencyGrid(Vector& sparse_f_grid,
                                    const Vector& f_grid,
                                    const Numeric& sparse_df,
                                    const String& speedup_option,
                                    // Verbosity object:
                                    const Verbosity&) {
  // Return empty for nothing
  if (not f_grid.nelem()) {
    sparse_f_grid.resize(0);
    return;
  };

  switch (Options::toLblSpeedupOrThrow(speedup_option)) {
    case Options::LblSpeedup::LinearIndependent:
      sparse_f_grid = LineShape::linear_sparse_f_grid(f_grid, sparse_df);
      ARTS_ASSERT(LineShape::good_linear_sparse_f_grid(f_grid, sparse_f_grid))
      break;
    case Options::LblSpeedup::QuadraticIndependent:
      sparse_f_grid = LineShape::triple_sparse_f_grid(f_grid, sparse_df);
      break;
    case Options::LblSpeedup::None:
      sparse_f_grid.resize(0);
      break;
    case Options::LblSpeedup::FINAL: { /* Leave last */
    }
  }
}

Vector create_sparse_f_grid_internal(const Vector& f_grid,
                                     const Numeric& sparse_df,
                                     const String& speedup_option,
                                     // Verbosity object:
                                     const Verbosity& verbosity) {
  Vector sparse_f_grid;
  sparse_f_gridFromFrequencyGrid(
      sparse_f_grid, f_grid, sparse_df, speedup_option, verbosity);
  return sparse_f_grid;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddLines(  // Workspace reference:
    // WS Output:
    PropagationMatrix& propmat_clearsky,
    StokesVector& nlte_source,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    ArrayOfStokesVector& dnlte_source_dx,
    // WS Input:
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const EnergyLevelMap& rtp_nlte,
    const Vector& rtp_vmr,
    const Index& nlte_do,
    const Index& lbl_checked,
    // WS User Generic inputs
    const Numeric& sparse_df,
    const Numeric& sparse_lim,
    const String& speedup_option,
    const Index& robust,
    // Verbosity object:
    const Verbosity& verbosity) {
  // Size of problem
  const Index nf = f_grid.nelem();
  const Index nq = jacobian_quantities.nelem();
  const Index ns = abs_species.nelem();

  // Possible things that can go wrong in this code (excluding line parameters)
  ARTS_USER_ERROR_IF(not lbl_checked, "Must check LBL calculations")
  check_abs_species(abs_species);
  ARTS_USER_ERROR_IF(rtp_vmr.nelem() not_eq abs_species.nelem(),
                     "*rtp_vmr* must match *abs_species*")
  ARTS_USER_ERROR_IF(propmat_clearsky.NumberOfFrequencies() not_eq nf,
                     "*f_grid* must match *propmat_clearsky*")
  ARTS_USER_ERROR_IF(nlte_source.NumberOfFrequencies() not_eq nf,
                     "*f_grid* must match *nlte_source*")
  ARTS_USER_ERROR_IF(
      not nq and (nq not_eq dpropmat_clearsky_dx.nelem()),
      "*dpropmat_clearsky_dx* must match derived form of *jacobian_quantities*")
  ARTS_USER_ERROR_IF(
      not nq and bad_propmat(dpropmat_clearsky_dx, f_grid),
      "*dpropmat_clearsky_dx* must have frequency dim same as *f_grid*")
  ARTS_USER_ERROR_IF(
      nlte_do and (nq not_eq dnlte_source_dx.nelem()),
      "*dnlte_source_dx* must match derived form of *jacobian_quantities* when non-LTE is on")
  ARTS_USER_ERROR_IF(
      nlte_do and bad_propmat(dnlte_source_dx, f_grid),
      "*dnlte_source_dx* must have frequency dim same as *f_grid* when non-LTE is on")
  ARTS_USER_ERROR_IF(any_negative(f_grid),
                     "Negative frequency (at least one value).")
  ARTS_USER_ERROR_IF((any_cutoff(abs_lines_per_species) or speedup_option not_eq "None") and not is_increasing(f_grid),
                     "Must be sorted and increasing if any cutoff or speedup is used.")
  ARTS_USER_ERROR_IF(any_negative(rtp_vmr),
                     "Negative VMR (at least one value).")
  ARTS_USER_ERROR_IF(any_negative(rtp_nlte.value),
                     "Negative NLTE (at least one value).")
  ARTS_USER_ERROR_IF(rtp_temperature <= 0, "Non-positive temperature")
  ARTS_USER_ERROR_IF(rtp_pressure <= 0, "Non-positive pressure")
  ARTS_USER_ERROR_IF(
      sparse_lim > 0 and sparse_df > sparse_lim,
      "If sparse grids are to be used, the limit must be larger than the grid-spacing.\n"
      "The limit is ",
      sparse_lim,
      " Hz and the grid_spacing is ",
      sparse_df,
      " Hz")

  if (not nf) return;

  // Deal with sparse computational grid
  const Vector f_grid_sparse = create_sparse_f_grid_internal(
      f_grid, sparse_df, speedup_option, verbosity);
  const Options::LblSpeedup speedup_type =
      f_grid_sparse.nelem() ? Options::toLblSpeedupOrThrow(speedup_option)
                            : Options::LblSpeedup::None;
  ARTS_USER_ERROR_IF(
      sparse_lim <= 0 and speedup_type not_eq Options::LblSpeedup::None,
      "Must have a sparse limit if you set speedup_option")

  // Calculations data
  LineShape::ComputeData com(f_grid, jacobian_quantities, nlte_do);
  LineShape::ComputeData sparse_com(
      f_grid_sparse, jacobian_quantities, nlte_do);

  if (arts_omp_in_parallel()) {
    for (Index ispecies = 0; ispecies < ns; ispecies++) {
      if (select_abs_species.nelem() and
          select_abs_species not_eq abs_species[ispecies])
        continue;

      // Skip it if there are no species or there is Zeeman requested
      if (not abs_species[ispecies].nelem() or abs_species[ispecies].Zeeman() or
          not abs_lines_per_species[ispecies].nelem())
        continue;

      for (auto& band : abs_lines_per_species[ispecies]) {
        LineShape::compute(com,
                          sparse_com,
                          band,
                          jacobian_quantities,
                          rtp_nlte,
                          band.BroadeningSpeciesVMR(rtp_vmr, abs_species),
                          abs_species[ispecies],
                          rtp_vmr[ispecies],
                          isotopologue_ratios[band.Isotopologue()],
                          rtp_pressure,
                          rtp_temperature,
                          0,
                          sparse_lim,
                          Zeeman::Polarization::None,
                          speedup_type,
                          robust not_eq 0);
      }
    }
  } else {  // In parallel
    const Index nbands = [](auto& lines) {
      Index n = 0;
      for (auto& abs_lines : lines) n += abs_lines.nelem();
      return n;
    }(abs_lines_per_species);

    std::vector<LineShape::ComputeData> vcom(
        arts_omp_get_max_threads(),
        LineShape::ComputeData{
            f_grid, jacobian_quantities, static_cast<bool>(nlte_do)});
    std::vector<LineShape::ComputeData> vsparse_com(
        arts_omp_get_max_threads(),
        LineShape::ComputeData{
            f_grid_sparse, jacobian_quantities, static_cast<bool>(nlte_do)});

#pragma omp parallel for schedule(dynamic)
    for (Index i = 0; i < nbands; i++) {
      const auto [ispecies, iband] =
          flat_index(i, abs_species, abs_lines_per_species);
          
      if (select_abs_species.nelem() and
          select_abs_species not_eq abs_species[ispecies])
        continue;

      // Skip it if there are no species or there is Zeeman requested
      if (not abs_species[ispecies].nelem() or abs_species[ispecies].Zeeman() or
          not abs_lines_per_species[ispecies].nelem())
        continue;

      auto& band = abs_lines_per_species[ispecies][iband];
      LineShape::compute(vcom[arts_omp_get_thread_num()],
                         vsparse_com[arts_omp_get_thread_num()],
                         band,
                         jacobian_quantities,
                         rtp_nlte,
                         band.BroadeningSpeciesVMR(rtp_vmr, abs_species),
                         abs_species[ispecies],
                         rtp_vmr[ispecies],
                         isotopologue_ratios[band.Isotopologue()],
                         rtp_pressure,
                         rtp_temperature,
                         0,
                         sparse_lim,
                         Zeeman::Polarization::None,
                         speedup_type,
                         robust not_eq 0);
    }

    for (auto& pcom: vcom) com += pcom;
    for (auto& pcom: vsparse_com) sparse_com += pcom;
  }

  switch (speedup_type) {
    case Options::LblSpeedup::LinearIndependent:
      com.interp_add_even(sparse_com);
      break;
    case Options::LblSpeedup::QuadraticIndependent:
      com.interp_add_triplequad(sparse_com);
      break;
    case Options::LblSpeedup::None: /* Do nothing */
      break;
    case Options::LblSpeedup::FINAL: { /* Leave last */
    }
  }

  // Sum up the propagation matrix
  propmat_clearsky.Kjj() += com.F.real();

  // Sum up the Jacobian
  for (Index j = 0; j < nq; j++) {
    if (not jacobian_quantities[j].propmattype()) continue;
    dpropmat_clearsky_dx[j].Kjj() += com.dF.real()(joker, j);
  }

  if (nlte_do) {
    // Sum up the source vector
    nlte_source.Kjj() += com.N.real();

    // Sum up the Jacobian
    for (Index j = 0; j < nq; j++) {
      if (not jacobian_quantities[j].propmattype()) continue;
      dnlte_source_dx[j].Kjj() += com.dN.real()(joker, j);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyZero(PropagationMatrix& propmat_clearsky,
                          const Vector& f_grid,
                          const Index& stokes_dim,
                          const Verbosity&) {
  propmat_clearsky = PropagationMatrix(f_grid.nelem(), stokes_dim);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyForceNegativeToZero(PropagationMatrix& propmat_clearsky,
                                         const Verbosity&) {
  for (Index i = 0; i < propmat_clearsky.NumberOfFrequencies(); i++)
    if (propmat_clearsky.Kjj()[i] < 0.0) propmat_clearsky.SetAtPosition(0.0, i);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void isotopologue_ratiosInitFromBuiltin(
    SpeciesIsotopologueRatios& isotopologue_ratios, const Verbosity&) {
  isotopologue_ratios = Species::isotopologue_ratiosInitFromBuiltin();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void isotopologue_ratiosInitFromHitran(
    SpeciesIsotopologueRatios& isotopologue_ratios, const Verbosity&) {
  isotopologue_ratios = Hitran::isotopologue_ratios();
}

#ifdef ENABLE_NETCDF
/* Workspace method: Doxygen documentation will be auto-generated */
/* Included by Claudia Emde, 20100707 */
void WriteMolTau(  //WS Input
    const Vector& f_grid,
    const Tensor3& z_field,
    const Tensor7& propmat_clearsky_field,
    const Index& atmosphere_dim,
    //Keyword
    const String& filename,
    const Verbosity&) {
  int retval, ncid;
  int nlev_dimid, nlyr_dimid, nwvl_dimid, stokes_dimid, none_dimid;
  int dimids[4];
  int wvlmin_varid, wvlmax_varid, z_varid, wvl_varid, tau_varid;

  ARTS_USER_ERROR_IF(atmosphere_dim != 1,
                     "WriteMolTau can only be used for atmosphere_dim=1")
#pragma omp critical(netcdf__critical_region)
  {
    // Open file
    if ((retval = nc_create(filename.c_str(), NC_CLOBBER, &ncid)))
      nca_error(retval, "nc_create");

    // Define dimensions
    if ((retval = nc_def_dim(ncid, "nlev", (int)z_field.npages(), &nlev_dimid)))
      nca_error(retval, "nc_def_dim");

    if ((retval =
             nc_def_dim(ncid, "nlyr", (int)z_field.npages() - 1, &nlyr_dimid)))
      nca_error(retval, "nc_def_dim");

    if ((retval = nc_def_dim(ncid, "nwvl", (int)f_grid.nelem(), &nwvl_dimid)))
      nca_error(retval, "nc_def_dim");

    if ((retval = nc_def_dim(ncid, "none", 1, &none_dimid)))
      nca_error(retval, "nc_def_dim");

    if ((retval = nc_def_dim(ncid,
                             "nstk",
                             (int)propmat_clearsky_field.nbooks(),
                             &stokes_dimid)))
      nca_error(retval, "nc_def_dim");

    // Define variables
    if ((retval = nc_def_var(
             ncid, "wvlmin", NC_DOUBLE, 1, &none_dimid, &wvlmin_varid)))
      nca_error(retval, "nc_def_var wvlmin");

    if ((retval = nc_def_var(
             ncid, "wvlmax", NC_DOUBLE, 1, &none_dimid, &wvlmax_varid)))
      nca_error(retval, "nc_def_var wvlmax");

    if ((retval = nc_def_var(ncid, "z", NC_DOUBLE, 1, &nlev_dimid, &z_varid)))
      nca_error(retval, "nc_def_var z");

    if ((retval =
             nc_def_var(ncid, "wvl", NC_DOUBLE, 1, &nwvl_dimid, &wvl_varid)))
      nca_error(retval, "nc_def_var wvl");

    dimids[0] = nlyr_dimid;
    dimids[1] = nwvl_dimid;
    dimids[2] = stokes_dimid;
    dimids[3] = stokes_dimid;

    if ((retval =
             nc_def_var(ncid, "tau", NC_DOUBLE, 4, &dimids[0], &tau_varid)))
      nca_error(retval, "nc_def_var tau");

    // Units
    if ((retval = nc_put_att_text(ncid, wvlmin_varid, "units", 2, "nm")))
      nca_error(retval, "nc_put_att_text");

    if ((retval = nc_put_att_text(ncid, wvlmax_varid, "units", 2, "nm")))
      nca_error(retval, "nc_put_att_text");

    if ((retval = nc_put_att_text(ncid, z_varid, "units", 2, "km")))
      nca_error(retval, "nc_put_att_text");

    if ((retval = nc_put_att_text(ncid, wvl_varid, "units", 2, "nm")))
      nca_error(retval, "nc_put_att_text");

    if ((retval = nc_put_att_text(ncid, tau_varid, "units", 1, "-")))
      nca_error(retval, "nc_put_att_text");

    // End define mode. This tells netCDF we are done defining
    // metadata.
    if ((retval = nc_enddef(ncid))) nca_error(retval, "nc_enddef");

    // Assign data
    double wvlmin[1];
    wvlmin[0] = SPEED_OF_LIGHT / f_grid[f_grid.nelem() - 1] * 1e9;
    if ((retval = nc_put_var_double(ncid, wvlmin_varid, &wvlmin[0])))
      nca_error(retval, "nc_put_var");

    double wvlmax[1];
    wvlmax[0] = SPEED_OF_LIGHT / f_grid[0] * 1e9;
    if ((retval = nc_put_var_double(ncid, wvlmax_varid, &wvlmax[0])))
      nca_error(retval, "nc_put_var");

    double z[z_field.npages()];
    for (int iz = 0; iz < z_field.npages(); iz++)
      z[iz] = z_field(z_field.npages() - 1 - iz, 0, 0) * 1e-3;

    if ((retval = nc_put_var_double(ncid, z_varid, &z[0])))
      nca_error(retval, "nc_put_var");

    double wvl[f_grid.nelem()];
    for (int iv = 0; iv < f_grid.nelem(); iv++)
      wvl[iv] = SPEED_OF_LIGHT / f_grid[f_grid.nelem() - 1 - iv] * 1e9;

    if ((retval = nc_put_var_double(ncid, wvl_varid, &wvl[0])))
      nca_error(retval, "nc_put_var");

    const Index zfnp = z_field.npages() - 1;
    const Index fgne = f_grid.nelem();
    const Index amfnb = propmat_clearsky_field.nbooks();

    Tensor4 tau(zfnp, fgne, amfnb, amfnb, 0.);

    // Calculate average tau for layers
    for (int is = 0; is < propmat_clearsky_field.nlibraries(); is++)
      for (int iz = 0; iz < zfnp; iz++)
        for (int iv = 0; iv < fgne; iv++)
          for (int is1 = 0; is1 < amfnb; is1++)
            for (int is2 = 0; is2 < amfnb; is2++)
              // sum up all species
              tau(iz, iv, is1, is2) +=
                  0.5 *
                  (propmat_clearsky_field(is,
                                          f_grid.nelem() - 1 - iv,
                                          is1,
                                          is2,
                                          z_field.npages() - 1 - iz,
                                          0,
                                          0) +
                   propmat_clearsky_field(is,
                                          f_grid.nelem() - 1 - iv,
                                          is1,
                                          is2,
                                          z_field.npages() - 2 - iz,
                                          0,
                                          0)) *
                  (z_field(z_field.npages() - 1 - iz, 0, 0) -
                   z_field(z_field.npages() - 2 - iz, 0, 0));

    if ((retval = nc_put_var_double(ncid, tau_varid, tau.get_c_array())))
      nca_error(retval, "nc_put_var");

    // Close the file
    if ((retval = nc_close(ncid))) nca_error(retval, "nc_close");
  }
}

#else

void WriteMolTau(  //WS Input
    const Vector& f_grid _U_,
    const Tensor3& z_field _U_,
    const Tensor7& propmat_clearsky_field _U_,
    const Index& atmosphere_dim _U_,
    //Keyword
    const String& filename _U_,
    const Verbosity&) {
  ARTS_USER_ERROR_IF(true,
                     "The workspace method WriteMolTau is not available"
                     "because ARTS was compiled without NetCDF support.");
}

#endif /* ENABLE_NETCDF */

void propmat_clearsky_agendaAuto(// Workspace reference:
    Workspace& ws,
    // WS Output:
    Agenda& propmat_clearsky_agenda,
    Index& propmat_clearsky_agenda_checked,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    // WS Generic Input:
    const Numeric& H,
    const Numeric& T_extrapolfac,
    const Numeric& eta,
    const Numeric& extpolfac,
    const Numeric& force_p,
    const Numeric& force_t,
    const Index& ignore_errors,
    const Numeric& lines_sparse_df,
    const Numeric& lines_sparse_lim,
    const String& lines_speedup_option,
    const Index& manual_mag_field,
    const Index& no_negatives,
    const Numeric& theta,
    const Index& use_abs_as_ext,
    const Index& use_abs_lookup_ind,
    // Verbosity object:
    const Verbosity& verbosity) {
  using namespace AgendaManip;

  propmat_clearsky_agenda_checked = 0;  // In case of crash

  AgendaCreator agenda(ws, "propmat_clearsky_agenda");

  // Use bool because logic is easier
  const bool use_abs_lookup = static_cast<bool>(use_abs_lookup_ind);

  const SpeciesTagTypeStatus any_species(abs_species);
  const AbsorptionTagTypesStatus any_lines(abs_lines_per_species);

  // propmat_clearskyInit
  agenda.add("propmat_clearskyInit");

  // propmat_clearskyAddFromLookup
  if (use_abs_lookup) {
    agenda.add("propmat_clearskyAddFromLookup",
               SetWsv{"extpolfac", extpolfac},
               SetWsv{"no_negatives", no_negatives});
  }

  // propmat_clearskyAddLines
  if (not use_abs_lookup and any_species.Plain and
      (any_lines.population.LTE or any_lines.population.NLTE or
       any_lines.population.VibTemps)) {
    agenda.add("propmat_clearskyAddLines",
               SetWsv{"lines_sparse_df", lines_sparse_df},
               SetWsv{"lines_sparse_lim", lines_sparse_lim},
               SetWsv{"lines_speedup_option", lines_speedup_option},
               SetWsv{"no_negatives", no_negatives});
  }

  // propmat_clearskyAddZeeman
  if (any_species.Zeeman and
      (any_lines.population.LTE or any_lines.population.NLTE or
       any_lines.population.VibTemps)) {
    agenda.add("propmat_clearskyAddZeeman",
               SetWsv{"manual_mag_field", manual_mag_field},
               SetWsv{"H", H},
               SetWsv{"theta", theta},
               SetWsv{"eta", eta});
  }

  //propmat_clearskyAddHitranXsec
  if (not use_abs_lookup and any_species.XsecFit) {
    agenda.add("propmat_clearskyAddXsecFit",
               SetWsv{"force_p", force_p},
               SetWsv{"force_t", force_t});
  }

  //propmat_clearskyAddOnTheFlyLineMixing
  if (not use_abs_lookup and any_species.Plain and
      (any_lines.population.ByMakarovFullRelmat or
       any_lines.population.ByRovibLinearDipoleLineMixing)) {
    agenda.add("propmat_clearskyAddOnTheFlyLineMixing");
  }

  //propmat_clearskyAddOnTheFlyLineMixingWithZeeman
  if (any_species.Zeeman and
      (any_lines.population.ByMakarovFullRelmat or
       any_lines.population.ByRovibLinearDipoleLineMixing)) {
    agenda.add("propmat_clearskyAddOnTheFlyLineMixingWithZeeman");
  }

  //propmat_clearskyAddCIA
  if (not use_abs_lookup and any_species.Cia) {
    agenda.add("propmat_clearskyAddCIA",
               SetWsv{"T_extrapolfac", T_extrapolfac},
               SetWsv{"ignore_errors", ignore_errors});
  }

  //propmat_clearskyAddPredefined
  if (not use_abs_lookup and any_species.Predefined) {
    agenda.add("propmat_clearskyAddPredefined");
  }

  //propmat_clearskyAddParticles
  if (any_species.Particles) {
    agenda.add("propmat_clearskyAddParticles",
               SetWsv{"use_abs_as_ext", use_abs_as_ext});
  }

  //propmat_clearskyAddFaraday
  if (any_species.FreeElectrons) {
    agenda.add("propmat_clearskyAddFaraday");
  }

  // propmat_clearskyAddHitranLineMixingLines
  if (not use_abs_lookup and any_species.Plain and
      (any_lines.population.ByHITRANFullRelmat or
       any_lines.population.ByHITRANRosenkranzRelmat)) {
    agenda.add("propmat_clearskyAddHitranLineMixingLines");
  }

  // Extra check (should really never ever fail when species exist)
  propmat_clearsky_agenda = agenda.finalize();
  propmat_clearsky_agenda_checked = 1;

  CREATE_OUT3;
  if (out3.sufficient_priority()) {
    out3 << "propmat_clearsky_agendaAuto sets propmat_clearsky_agenda to:\n\n"
         << propmat_clearsky_agenda << '\n';
  }
}
