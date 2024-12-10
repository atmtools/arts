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
#include "enumsAtmKey.h"
#include "file.h"
#include "hitran_species.h"
#include "jacobian.h"
#include "lbl_data.h"
#include "math_funcs.h"
#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_view.h"
#include "nlte.h"
#include "optproperties.h"
#include "path_point.h"
#include "sorted_grid.h"
#include "species.h"
#include "species_tags.h"

#ifdef ENABLE_NETCDF
#include <netcdf.h>

#include "nc_io.h"
#endif

void mirror_los(Vector& los_mirrored, const ConstVectorView& los) {
  los_mirrored.resize(2);
  los_mirrored[0] = 180 - los[0];
  los_mirrored[1] = los[1] + 180;
  if (los_mirrored[1] > 180) {
    los_mirrored[1] -= 360;
  }
}

Numeric dotprod_with_los(const ConstVectorView& los,
                         const Numeric& u,
                         const Numeric& v,
                         const Numeric& w) {
  // Strength of field
  const Numeric f = sqrt(u * u + v * v + w * w);

  // Zenith and azimuth angle for field (in radians)
  const Numeric za_f = acos(w / f);
  const Numeric aa_f = atan2(u, v);

  // Zenith and azimuth angle for photon direction (in radians)
  Vector los_p;
  mirror_los(los_p, los);
  const Numeric za_p = Conversion::deg2rad(1) * los_p[0];
  const Numeric aa_p = Conversion::deg2rad(1) * los_p[1];

  return f * (cos(za_f) * cos(za_p) + sin(za_f) * sin(za_p) * cos(aa_f - aa_p));
}

inline constexpr Numeric ELECTRON_CHARGE     = -Constant::elementary_charge;
inline constexpr Numeric ELECTRON_MASS       = Constant::electron_mass;
inline constexpr Numeric PI                  = Constant::pi;
inline constexpr Numeric SPEED_OF_LIGHT      = Constant::speed_of_light;
inline constexpr Numeric VACUUM_PERMITTIVITY = Constant::vacuum_permittivity;

/* Workspace method: Doxygen documentation will be auto-generated */
void absorption_speciesSet(  // WS Output:
    ArrayOfArrayOfSpeciesTag& absorption_species,
    // Control Parameters:
    const ArrayOfString& names) try {
  absorption_species.resize(names.size());

  //cout << "Names: " << names << "\n";

  // Each element of the array of Strings names defines one tag
  // group. Let's work through them one by one.
  for (Size i = 0; i < names.size(); ++i) {
    // This part has now been moved to array_species_tag_from_string.
    // Call this function.
    absorption_species[i] = ArrayOfSpeciesTag(names[i]);
  }
}
ARTS_METHOD_ERROR_CATCH

/* Workspace method: Doxygen documentation will be auto-generated */
void absorption_speciesDefineAllInScenario(  // WS Output:
    ArrayOfArrayOfSpeciesTag& tgs,
    // Control Parameters:
    const String& basename) {
  // We want to make lists of included and excluded species:
  ArrayOfString included(0), excluded(0);

  tgs.resize(0);

  for (Size i = 0; i < enumsize::SpeciesEnumSize; ++i) {
    const String specname{toString<1>(SpeciesEnum(i))};

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
  for (Size i = 0; i < enumsize::SpeciesEnumSize; ++i) {
    if (SpeciesEnum(i) not_eq SpeciesEnum::Bath) {
      specs.emplace_back(toString<1>(SpeciesEnum(i)));
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
  ARTS_USER_ERROR_IF(1 != 3, "Atmospheric dimension must be 1D, but 3 is 3")

  abs_p    = p_grid;
  abs_t    = t_field(joker, 0, 0);
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

  const Vector rtp_los{ExhaustiveConstVectorView{path::mirror(path_point.los)}};

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
      AtmKey::mag_u,
      AtmKey::mag_v,
      AtmKey::mag_w,
      AtmKey::wind_u,
      AtmKey::wind_v,
      AtmKey::wind_w,
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
      const Numeric f2            = frequency_grid[iv] * frequency_grid[iv];
      const Numeric r             = ne * c1 / f2;
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
  Index np       = 0;
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
      "*scat_data* needs to be identical. But you have {} 'particles' entries\nand {} *scat_data* elements.\n",
      np,
      ns)

  const auto jac_temperature =
      jacobian_targets.find<Jacobian::AtmTarget>(AtmKey::t);
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
    T_array     = atm_point.temperature;
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
  Index sp        = 0;
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
          "scat species #{}, scat elem #{} (vmr_field entry #{})\n",
          i_ss,
          i_se,
          sp)

      if (atm_point[absorption_species[sp].Species()] > 0.) {
        ARTS_USER_ERROR_IF(t_ok(i_se_flat, 0) < 0.,
                           "Temperature interpolation error:\n"
                           "scat species #{}, scat elem #{}\n",
                           i_ss,
                           i_se)
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
                         "scat species #{}, scat elem #{}\n",
                         i_ss,
                         i_se)

      if (jac_temperature.first) {
        const auto iq = jac_temperature.second->target_pos;

        if (use_abs_as_ext) {
          tmp(joker, joker, 0)  = abs_vec_Nse[i_ss][i_se](joker, 1, 0, joker);
          tmp(joker, joker, 0) -= abs_vec_Nse[i_ss][i_se](joker, 0, 0, joker);
        } else {
          tmp  = ext_mat_Nse[i_ss][i_se](joker, 1, 0, joker, joker);
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

void propagation_matrix_agendaAuto(
    Agenda& propagation_matrix_agenda,
    const ArrayOfArrayOfSpeciesTag& absorption_species,
    const AbsorptionBands& absorption_bands,
    const Index& use_absorption_lookup_table,
    const Numeric& T_extrapolfac,
    const Index& ignore_errors,
    const Index& no_negative_absorption,
    const Numeric& force_p,
    const Numeric& force_t,
    const Index& p_interp_order,
    const Index& t_interp_order,
    const Index& water_interp_order,
    const Index& f_interp_order,
    const Numeric& extpolfac) {
  AgendaCreator agenda("propagation_matrix_agenda");

  const SpeciesTagTypeStatus any_species(absorption_species);

  // propagation_matrixInit
  agenda.add("propagation_matrixInit");

  // propagation_matrixAddLines or propagation_matrixAddLookup
  if (use_absorption_lookup_table) {
    agenda.add("propagation_matrixAddLookup",
               SetWsv{"no_negative_absorption", no_negative_absorption},
               SetWsv{"p_interp_order", p_interp_order},
               SetWsv{"t_interp_order", t_interp_order},
               SetWsv{"water_interp_order", water_interp_order},
               SetWsv{"f_interp_order", f_interp_order},
               SetWsv{"extpolfac", extpolfac});
  } else if (absorption_bands.size()) {
    agenda.add("propagation_matrixAddLines",
               SetWsv{"no_negative_absorption", no_negative_absorption});
  }

  //propagation_matrixAddHitranXsec
  if (any_species.XsecFit) {
    agenda.add("propagation_matrixAddXsecFit",
               SetWsv{"force_p", force_p},
               SetWsv{"force_t", force_t});
  }

  //propagation_matrixAddCIA
  if (any_species.Cia) {
    agenda.add("propagation_matrixAddCIA",
               SetWsv{"T_extrapolfac", T_extrapolfac},
               SetWsv{"ignore_errors", ignore_errors});
  }

  //propagation_matrixAddPredefined
  if (any_species.Predefined) {
    agenda.add("propagation_matrixAddPredefined");
  }

  //propagation_matrixAddFaraday
  if (std::ranges::any_of(absorption_species,
                          [](auto& spec) { return spec.FreeElectrons(); })) {
    agenda.add("propagation_matrixAddFaraday");
  }

  agenda.add("propagation_matrix_jacobianWindFix");

  // Extra check (should really never ever fail when species exist)
  propagation_matrix_agenda = std::move(agenda).finalize(true);
}
