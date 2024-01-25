//

/**
   \file   m_abs.cc

   Stuff related to the calculation of absorption coefficients.

   \author Stefan Buehler
   \date   2001-03-12
*/
#include <workspace.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <utility>

#include "absorptionlines.h"
#include "array.h"
#include "arts_constants.h"
#include "arts_omp.h"
#include "atm.h"
#include "check_input.h"
#include "debug.h"
#include "file.h"
#include "hitran_species.h"
#include "jacobian.h"
#include "lbl_data.h"
#include "lineshape.h"
#include "math_funcs.h"
#include "matpack_concepts.h"
#include "matpack_data.h"
#include "nlte.h"
#include "optproperties.h"
#include "path_point.h"
#include "rte.h"
#include "species.h"
#include "species_tags.h"

#ifdef ENABLE_NETCDF
#include <netcdf.h>

#include "nc_io.h"
#endif

inline constexpr Numeric ELECTRON_CHARGE = -Constant::elementary_charge;
inline constexpr Numeric ELECTRON_MASS = Constant::electron_mass;
inline constexpr Numeric PI = Constant::pi;
inline constexpr Numeric SPEED_OF_LIGHT = Constant::speed_of_light;
inline constexpr Numeric VACUUM_PERMITTIVITY = Constant::vacuum_permittivity;

/* Workspace method: Doxygen documentation will be auto-generated */
void absorption_speciesDefineAllInScenario(  // WS Output:
    ArrayOfArrayOfSpeciesTag& tgs,
    // Control Parameters:
    const String& basename) {
  // We want to make lists of included and excluded species:
  ArrayOfString included(0), excluded(0);

  tgs.resize(0);

  for (Index i = 0; i < Index(Species::Species::FINAL); ++i) {
    const String specname{Species::toShortName(Species::Species(i))};

    String filename = basename;
    if (basename.length() && basename[basename.length() - 1] != '/')
      filename += ".";
    filename += specname;

    try {
      find_xml_file(filename);
      // Add to included list:
      included.push_back(specname);

      // Add this tag group to tgs:
      tgs.emplace_back(ArrayOfSpeciesTag(specname));
    } catch (const std::runtime_error& e) {
      // The file for the species could not be found.
      excluded.push_back(specname);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void absorption_speciesDefineAll(  // WS Output:
    ArrayOfArrayOfSpeciesTag& absorption_species) {
  // Species lookup data:

  // We want to make lists of all species
  ArrayOfString specs(0);
  for (Index i = 0; i < Index(Species::Species::FINAL); ++i) {
    if (Species::Species(i) not_eq Species::Species::Bath) {
      specs.emplace_back(Species::toShortName(Species::Species(i)));
    }
  }

  // Set the values
  absorption_speciesSet(absorption_species, specs);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void AbsInputFromAtmFields(  // WS Output:
    Vector& abs_p,
    Vector& abs_t,
    Matrix& abs_vmrs,
    // WS Input:
    const Vector& p_grid,
    const Tensor3& t_field,
    const Tensor4& vmr_field) {
  // First, make sure that we really have a 1D atmosphere:
  ARTS_USER_ERROR_IF(
      1 != 3, "Atmospheric dimension must be 1D, but 3 is ", 3, ".")

  abs_p = p_grid;
  abs_t = t_field(joker, 0, 0);
  abs_vmrs = vmr_field(joker, joker, 0, 0);
}

//======================================================================
//             Methods related to continua
//======================================================================

/* Workspace method: Doxygen documentation will be auto-generated */
void propagation_matrixInit(  //WS Output
    PropmatVector& propagation_matrix,
    StokvecVector& source_vector_nonlte,
    PropmatMatrix& propagation_matrix_jacobian,
    StokvecMatrix& source_vector_nonlte_jacobian,
    //WS Input
    const JacobianTargets& jacobian_targets,
    const AscendingGrid& frequency_grid) {
  const Index nf = frequency_grid.nelem();
  const Index nq = jacobian_targets.target_count();

  ARTS_USER_ERROR_IF(not nf, "No frequencies");

  // Set size of propagation_matrix and reset it's values
  propagation_matrix.resize(nf);
  propagation_matrix = 0.0;

  // Set size of source_vector_nonlte and reset it's values
  source_vector_nonlte.resize(nf);
  source_vector_nonlte = 0.0;

  // Set size of propagation_matrix_jacobian and reset it's values
  propagation_matrix_jacobian.resize(nq, nf);
  propagation_matrix_jacobian = 0.0;

  // Set size of source_vector_nonlte_jacobian and reset it's values
  source_vector_nonlte_jacobian.resize(nq, nf);
  source_vector_nonlte_jacobian = 0.0;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propagation_matrixAddFaraday(
    PropmatVector& propagation_matrix,
    PropmatMatrix& propagation_matrix_jacobian,
    const AscendingGrid& frequency_grid,
    const ArrayOfArrayOfSpeciesTag& absorption_species,
    const SpeciesEnum& select_absorption_species,
    const JacobianTargets& jacobian_targets,
    const AtmPoint& atm_point,
    const PropagationPathPoint& path_point) {
  Index ife = -1;
  for (Size sp = 0; sp < absorption_species.size() && ife < 0; sp++) {
    if (absorption_species[sp].FreeElectrons()) {
      ife = sp;
    }
  }

  const Vector rtp_los{path::mirror(path_point.los)};

  ARTS_USER_ERROR_IF(ife < 0,
                     "Free electrons not found in *absorption_species* and "
                     "Faraday rotation can not be calculated.");

  // Allow early exit for lookup table calculations
  if (select_absorption_species != SpeciesEnum::Bath and
      select_absorption_species != absorption_species[ife].Species())
    return;

  // All the physical constants joined into one static constant:
  // (abs as e defined as negative)
  static const Numeric FRconst =
      abs(ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE /
          (8 * PI * PI * SPEED_OF_LIGHT * VACUUM_PERMITTIVITY * ELECTRON_MASS *
           ELECTRON_MASS));

  const auto jacs = jacobian_targets.find_all<Jacobian::AtmTarget>(
      Atm::Key::mag_u,
      Atm::Key::mag_v,
      Atm::Key::mag_w,
      Atm::Key::wind_u,
      Atm::Key::wind_v,
      Atm::Key::wind_w,
      absorption_species[ife].Species());
  const Numeric dmag = field_perturbation(std::span{jacs.data(), 3});

  const Numeric ne = atm_point[absorption_species[ife].Species()];

  if (ne != 0 && not atm_point.zero_mag()) {
    // Include remaining terms, beside /f^2
    const Numeric c1 =
        2 * FRconst *
        dotprod_with_los(
            rtp_los, atm_point.mag[0], atm_point.mag[1], atm_point.mag[2]);

    std::array<Numeric, 3> dc1{0., 0., 0.};
    if (dmag != 0.0) {
      dc1[0] = (2 * FRconst *
                    dotprod_with_los(rtp_los,
                                     atm_point.mag[0] + dmag,
                                     atm_point.mag[1],
                                     atm_point.mag[2]) -
                c1) /
               dmag;
      dc1[1] = (2 * FRconst *
                    dotprod_with_los(rtp_los,
                                     atm_point.mag[0],
                                     atm_point.mag[1] + dmag,
                                     atm_point.mag[2]) -
                c1) /
               dmag;
      dc1[2] = (2 * FRconst *
                    dotprod_with_los(rtp_los,
                                     atm_point.mag[0],
                                     atm_point.mag[1],
                                     atm_point.mag[2] + dmag) -
                c1) /
               dmag;
    }

    for (Index iv = 0; iv < frequency_grid.nelem(); iv++) {
      const Numeric f2 = frequency_grid[iv] * frequency_grid[iv];
      const Numeric r = ne * c1 / f2;
      propagation_matrix[iv].U() += r;

      for (Size i = 0; i < 3; i++) {
        if (jacs[i].first) {
          propagation_matrix_jacobian(jacs[i].second->target_pos, iv).U() +=
              ne * dc1[i] / f2;
        }
      }

      for (Size i = 3; i < 6; i++) {
        if (jacs[i].first) {
          propagation_matrix_jacobian(jacs[i].second->target_pos, iv).U() +=
              -2.0 * ne * r / frequency_grid[iv];
        }
      }

      if (jacs[6].first) {
        propagation_matrix_jacobian(jacs[6].second->target_pos, iv).U() += r;
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propagation_matrixAddParticles(
    // WS Output:
    PropmatVector& propagation_matrix,
    PropmatMatrix& propagation_matrix_jacobian,
    // WS Input:
    const Vector& frequency_grid,
    const ArrayOfArrayOfSpeciesTag& absorption_species,
    const ArrayOfSpeciesTag& select_absorption_species,
    const JacobianTargets& jacobian_targets,
    const PropagationPathPoint& path_point,
    const AtmPoint& atm_point,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Index& scat_data_checked,
    const Index& use_abs_as_ext) {
  ARTS_USER_ERROR_IF(select_absorption_species.size(), R"--(
  We do not yet support select_absorption_species for lookup table calculations
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
  for (Size sp = 0; sp < absorption_species.size(); sp++) {
    if (absorption_species[sp].Particles()) {
      np++;
    }
  }

  ARTS_USER_ERROR_IF(
      np == 0,
      "For applying propagation_matrixAddParticles, *absorption_species* needs to"
      "contain species 'particles', but it does not.\n")

  ARTS_USER_ERROR_IF(
      ns != np,
      "Number of 'particles' entries in absorption_species and of elements in\n"
      "*scat_data* needs to be identical. But you have ",
      np,
      " 'particles' entries\n"
      "and ",
      ns,
      " *scat_data* elements.\n")

  const auto jac_temperature =
      jacobian_targets.find<Jacobian::AtmTarget>(Atm::Key::t);
  const Numeric dT = jac_temperature.first ? jac_temperature.second->d : 0.0;

  const Index na = absorption_species.size();
  const Vector rtp_los_back{path_point.los};

  // creating temporary output containers
  ArrayOfArrayOfTensor5 ext_mat_Nse;
  ArrayOfArrayOfTensor4 abs_vec_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;

  // preparing input in format needed
  Vector T_array;
  if (jac_temperature.first) {
    T_array.resize(2);
    T_array = atm_point.temperature;
    T_array[1] += dT;
  } else {
    T_array.resize(1);
    T_array = atm_point.temperature;
  }
  Matrix dir_array(1, 2);
  dir_array(0, joker) = rtp_los_back;

  // ext/abs per scat element for all freqs at once
  opt_prop_NScatElems(ext_mat_Nse,
                      abs_vec_Nse,
                      ptypes_Nse,
                      t_ok,
                      scat_data,
                      T_array,
                      dir_array,
                      -1);

  const Index nf = abs_vec_Nse[0][0].nbooks();
  Tensor3 tmp(nf, 4, 4);

  // Internal computations necessary since it relies on zero start
  PropmatVector internal_propmat(propagation_matrix.nelem());

  // loop over the scat_data and link them with correct vmr_field entry according
  // to the position of the particle type entries in absorption_species.
  Index sp = 0;
  Index i_se_flat = 0;
  for (Size i_ss = 0; i_ss < scat_data.size(); i_ss++) {
    for (Size i_se = 0; i_se < scat_data[i_ss].size(); i_se++) {
      // forward to next particle entry in absorption_species
      while (sp < na && not absorption_species[sp].Particles()) sp++;
      internal_propmat = 0.0;

      // running beyond number of absorption_species entries when looking for next
      // particle entry. shouldn't happen, though.
      ARTS_ASSERT(sp < na);
      ARTS_USER_ERROR_IF(
          atm_point[absorption_species[sp].Species()] < 0.,
          "Negative absorbing particle 'vmr' (aka number density)"
          " encountered:\n"
          "scat species #",
          i_ss,
          ", scat elem #",
          i_se,
          " (vmr_field entry #",
          sp,
          ")\n")

      if (atm_point[absorption_species[sp].Species()] > 0.) {
        ARTS_USER_ERROR_IF(t_ok(i_se_flat, 0) < 0.,
                           "Temperature interpolation error:\n"
                           "scat species #",
                           i_ss,
                           ", scat elem #",
                           i_se,
                           "\n")
        if (use_abs_as_ext) {
          for (Index iv = 0; iv < frequency_grid.nelem(); iv++) {
            internal_propmat[iv].A() += abs_vec_Nse[i_ss][i_se](iv, 0, 0, 0);
            internal_propmat[iv].B() += abs_vec_Nse[i_ss][i_se](iv, 0, 0, 1);
            internal_propmat[iv].C() += abs_vec_Nse[i_ss][i_se](iv, 0, 0, 2);
            internal_propmat[iv].D() += abs_vec_Nse[i_ss][i_se](iv, 0, 0, 3);
          }
        } else {
          for (Index iv = 0; iv < frequency_grid.nelem(); iv++) {
            internal_propmat[iv] = rtepack::to_propmat(
                ext_mat_Nse[i_ss][i_se](iv, 0, 0, joker, joker));
          }
        }

        const Numeric vmr = atm_point[absorption_species[sp].Species()];
        for (Index iv = 0; iv < frequency_grid.nelem(); iv++) {
          propagation_matrix[iv] += vmr * internal_propmat[iv];
        }
      }

      // For temperature derivatives (so we don't need to check it in jac loop)
      ARTS_USER_ERROR_IF(jac_temperature.first and t_ok(i_se_flat, 1) < 0.,
                         "Temperature interpolation error (in perturbation):\n"
                         "scat species #",
                         i_ss,
                         ", scat elem #",
                         i_se,
                         "\n")

      if (jac_temperature.first) {
        const auto iq = jac_temperature.second->target_pos;

        if (use_abs_as_ext) {
          tmp(joker, joker, 0) = abs_vec_Nse[i_ss][i_se](joker, 1, 0, joker);
          tmp(joker, joker, 0) -= abs_vec_Nse[i_ss][i_se](joker, 0, 0, joker);
        } else {
          tmp = ext_mat_Nse[i_ss][i_se](joker, 1, 0, joker, joker);
          tmp -= ext_mat_Nse[i_ss][i_se](joker, 0, 0, joker, joker);
        }

        tmp *= atm_point[absorption_species[sp].Species()];
        tmp /= dT;

        for (Index iv = 0; iv < frequency_grid.nelem(); iv++) {
          if (use_abs_as_ext) {
            propagation_matrix_jacobian(iq, iv).A() += tmp(iv, 0, 0);
            propagation_matrix_jacobian(iq, iv).B() += tmp(iv, 1, 0);
            propagation_matrix_jacobian(iq, iv).C() += tmp(iv, 2, 0);
            propagation_matrix_jacobian(iq, iv).D() += tmp(iv, 3, 0);
          } else {
            propagation_matrix_jacobian(iq, iv) +=
                rtepack::to_propmat(tmp(iv, joker, joker));
          }
        }
      }

      if (const auto jac_species = jacobian_targets.find<Jacobian::AtmTarget>(
              absorption_species[sp].Species());
          jac_species.first) {
        const auto iq = jac_species.second->target_pos;

        for (Index iv = 0; iv < frequency_grid.nelem(); iv++)
          propagation_matrix_jacobian(iq, iv) += internal_propmat[iv];
      }

      sp++;
      i_se_flat++;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propagation_matrixZero(PropmatVector& propagation_matrix,
                            const AscendingGrid& frequency_grid) {
  propagation_matrix = PropmatVector(frequency_grid.nelem());
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propagation_matrixForceNegativeToZero(PropmatVector& propagation_matrix) {
  for (Index i = 0; i < propagation_matrix.nelem(); i++)
    if (propagation_matrix[i].A() < 0.0) propagation_matrix[i] = 0.0;
  ;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void isotopologue_ratiosInitFromBuiltin(
    SpeciesIsotopologueRatios& isotopologue_ratios) {
  isotopologue_ratios = Species::isotopologue_ratiosInitFromBuiltin();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void isotopologue_ratiosInitFromHitran(
    SpeciesIsotopologueRatios& isotopologue_ratios) {
  isotopologue_ratios = Hitran::isotopologue_ratios();
}

#ifdef ENABLE_NETCDF
/* Workspace method: Doxygen documentation will be auto-generated */
/* Included by Claudia Emde, 20100707 */
void WriteMolTau(  //WS Input
    const Vector& frequency_grid,
    const Tensor3& z_field,
    const Tensor7& propagation_matrix_field,
    //Keyword
    const String& filename) {
  int retval, ncid;
  int nlev_dimid, nlyr_dimid, nwvl_dimid, stokes_dimid, none_dimid;
  int dimids[4];
  int wvlmin_varid, wvlmax_varid, z_varid, wvl_varid, tau_varid;

  ARTS_USER_ERROR_IF(3 != 1, "WriteMolTau can only be used for 3=1")
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

    if ((retval = nc_def_dim(
             ncid, "nwvl", (int)frequency_grid.nelem(), &nwvl_dimid)))
      nca_error(retval, "nc_def_dim");

    if ((retval = nc_def_dim(ncid, "none", 1, &none_dimid)))
      nca_error(retval, "nc_def_dim");

    if ((retval = nc_def_dim(ncid,
                             "nstk",
                             (int)propagation_matrix_field.nbooks(),
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
    wvlmin[0] =
        SPEED_OF_LIGHT / frequency_grid[frequency_grid.nelem() - 1] * 1e9;
    if ((retval = nc_put_var_double(ncid, wvlmin_varid, &wvlmin[0])))
      nca_error(retval, "nc_put_var");

    double wvlmax[1];
    wvlmax[0] = SPEED_OF_LIGHT / frequency_grid[0] * 1e9;
    if ((retval = nc_put_var_double(ncid, wvlmax_varid, &wvlmax[0])))
      nca_error(retval, "nc_put_var");

    double z[z_field.npages()];
    for (int iz = 0; iz < z_field.npages(); iz++)
      z[iz] = z_field(z_field.npages() - 1 - iz, 0, 0) * 1e-3;

    if ((retval = nc_put_var_double(ncid, z_varid, &z[0])))
      nca_error(retval, "nc_put_var");

    double wvl[frequency_grid.nelem()];
    for (int iv = 0; iv < frequency_grid.nelem(); iv++)
      wvl[iv] = SPEED_OF_LIGHT /
                frequency_grid[frequency_grid.nelem() - 1 - iv] * 1e9;

    if ((retval = nc_put_var_double(ncid, wvl_varid, &wvl[0])))
      nca_error(retval, "nc_put_var");

    const Index zfnp = z_field.npages() - 1;
    const Index fgne = frequency_grid.nelem();
    const Index amfnb = propagation_matrix_field.nbooks();

    Tensor4 tau(zfnp, fgne, amfnb, amfnb, 0.);

    // Calculate average tau for layers
    for (int is = 0; is < propagation_matrix_field.nlibraries(); is++)
      for (int iz = 0; iz < zfnp; iz++)
        for (int iv = 0; iv < fgne; iv++)
          for (int is1 = 0; is1 < amfnb; is1++)
            for (int is2 = 0; is2 < amfnb; is2++)
              // sum up all species
              tau(iz, iv, is1, is2) +=
                  0.5 *
                  (propagation_matrix_field(is,
                                            frequency_grid.nelem() - 1 - iv,
                                            is1,
                                            is2,
                                            z_field.npages() - 1 - iz,
                                            0,
                                            0) +
                   propagation_matrix_field(is,
                                            frequency_grid.nelem() - 1 - iv,
                                            is1,
                                            is2,
                                            z_field.npages() - 2 - iz,
                                            0,
                                            0)) *
                  (z_field(z_field.npages() - 1 - iz, 0, 0) -
                   z_field(z_field.npages() - 2 - iz, 0, 0));

    if ((retval = nc_put_var_double(ncid, tau_varid, tau.unsafe_data_handle())))
      nca_error(retval, "nc_put_var");

    // Close the file
    if ((retval = nc_close(ncid))) nca_error(retval, "nc_close");
  }
}

#else

void WriteMolTau(  //WS Input
    const Vector& frequency_grid [[maybe_unused]],
    const Tensor3& z_field [[maybe_unused]],
    const Tensor7& propagation_matrix_field [[maybe_unused]],
    //Keyword
    const String& filename [[maybe_unused]]) {
  ARTS_USER_ERROR_IF(true,
                     "The workspace method WriteMolTau is not available"
                     "because ARTS was compiled without NetCDF support.");
}

#endif /* ENABLE_NETCDF */

void propagation_matrix_agendaAuto(  // Workspace reference:
    Agenda& propagation_matrix_agenda,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& absorption_species,
    const AbsorptionBands& absorption_bands,
    // WS Generic Input:
    const Numeric& T_extrapolfac,
    const Numeric& extpolfac,
    const Numeric& force_p,
    const Numeric& force_t,
    const Index& ignore_errors,
    const Index& no_negatives,
    const Index& use_abs_lookup_ind) {
  AgendaCreator agenda("propagation_matrix_agenda");

  // Use bool because logic is easier
  const bool use_abs_lookup = static_cast<bool>(use_abs_lookup_ind);

  const SpeciesTagTypeStatus any_species(absorption_species);

  // propagation_matrixInit
  agenda.add("propagation_matrixInit");

  // propagation_matrixAddFromLookup
  if (use_abs_lookup) {
    agenda.add("propagation_matrixAddFromLookup",
               SetWsv{"extpolfac", extpolfac},
               SetWsv{"no_negatives", no_negatives});
  }

  // propagation_matrixAddLines
  if (not use_abs_lookup and absorption_bands.size()) {
    agenda.add("propagation_matrixAddLines");
  }

  //propagation_matrixAddHitranXsec
  if (not use_abs_lookup and any_species.XsecFit) {
    agenda.add("propagation_matrixAddXsecFit",
               SetWsv{"force_p", force_p},
               SetWsv{"force_t", force_t});
  }

  //propagation_matrixAddCIA
  if (not use_abs_lookup and any_species.Cia) {
    agenda.add("propagation_matrixAddCIA",
               SetWsv{"T_extrapolfac", T_extrapolfac},
               SetWsv{"ignore_errors", ignore_errors});
  }

  //propagation_matrixAddPredefined
  if (not use_abs_lookup and any_species.Predefined) {
    agenda.add("propagation_matrixAddPredefined");
  }

  //propagation_matrixAddFaraday
  if (any_species.FreeElectrons) {
    agenda.add("propagation_matrixAddFaraday");
  }

  // Extra check (should really never ever fail when species exist)
  propagation_matrix_agenda = std::move(agenda).finalize();
}
