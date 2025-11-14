//

/**
   \file   m_abs.cc

   Stuff related to the calculation of absorption coefficients.

   \author Stefan Buehler
   \date   2001-03-12
*/
#include <arts_constants.h>
#include <debug.h>
#include <enumsAtmKey.h>
#include <file.h>
#include <hitran_species.h>
#include <jacobian.h>
#include <lbl_data.h>
#include <path_point.h>
#include <species_tags.h>
#include <workspace.h>

#include <algorithm>
#include <cmath>
#include <ranges>
#include <utility>

namespace {
void mirror_los(Vector& los_mirrored, const ConstVectorView& los) {
  ARTS_TIME_REPORT
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
  ARTS_TIME_REPORT
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
}  // namespace

inline constexpr Numeric ELECTRON_CHARGE     = -Constant::elementary_charge;
inline constexpr Numeric ELECTRON_MASS       = Constant::electron_mass;
inline constexpr Numeric PI                  = Constant::pi;
inline constexpr Numeric SPEED_OF_LIGHT      = Constant::speed_of_light;
inline constexpr Numeric VACUUM_PERMITTIVITY = Constant::vacuum_permittivity;

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesSet(ArrayOfSpeciesTag& abs_species,
                    const ArrayOfString& names) try {
  ARTS_TIME_REPORT

  abs_species.resize(names.size());

  const auto op = [](const String& name) { return SpeciesTag{name}; };

  stdr::transform(names, abs_species.begin(), op);
}
ARTS_METHOD_ERROR_CATCH

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesDefineAll(  // WS Output:
    ArrayOfSpeciesTag& abs_species) {
  ARTS_TIME_REPORT

  // Species lookup data:

  // We want to make lists of all species
  ArrayOfString specs(0);
  for (Size i = 0; i < enumsize::SpeciesEnumSize; ++i) {
    if (SpeciesEnum(i) not_eq SpeciesEnum::Bath) {
      specs.emplace_back(toString<1>(SpeciesEnum(i)));
    }
  }

  // Set the values
  abs_speciesSet(abs_species, specs);
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
    const JacobianTargets& jac_targets,
    const AscendingGrid& freq_grid) {
  ARTS_TIME_REPORT

  const Index nf = freq_grid.size();
  const Index nq = jac_targets.target_count();

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
void propagation_matrixAddFaraday(PropmatVector& propagation_matrix,
                                  PropmatMatrix& propagation_matrix_jacobian,
                                  const AscendingGrid& freq_grid,
                                  const SpeciesEnum& select_abs_species,
                                  const JacobianTargets& jac_targets,
                                  const AtmPoint& atm_point,
                                  const PropagationPathPoint& path_point) {
  ARTS_TIME_REPORT

  constexpr SpeciesEnum electrons_key = "free_electrons"_spec;

  if (select_abs_species != electrons_key or
      select_abs_species == SpeciesEnum::Bath) {
    // If the selected species is not free electrons, we do not add Faraday rotation.
    return;
  }

  const Vector rtp_los{path::mirror(path_point.los)};

  // All the physical constants joined into one static constant:
  // (abs as e defined as negative)
  constexpr Numeric FRconst =
      nonstd::abs(ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE /
                  (8 * PI * PI * SPEED_OF_LIGHT * VACUUM_PERMITTIVITY *
                   ELECTRON_MASS * ELECTRON_MASS));

  const auto end  = jac_targets.atm.end();
  const auto jacs = std::array{jac_targets.find(AtmKey::mag_u),
                               jac_targets.find(AtmKey::mag_v),
                               jac_targets.find(AtmKey::mag_w),
                               jac_targets.find(AtmKey::wind_u),
                               jac_targets.find(AtmKey::wind_v),
                               jac_targets.find(AtmKey::wind_w),
                               jac_targets.find(electrons_key)};

  const Numeric dmag = jacs[0] != end   ? jacs[0]->d
                       : jacs[1] != end ? jacs[1]->d
                       : jacs[2] != end ? jacs[2]->d
                                        : 0.0;

  const Numeric ne = atm_point[electrons_key];

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

    for (Size iv = 0; iv < freq_grid.size(); iv++) {
      const Numeric f2            = freq_grid[iv] * freq_grid[iv];
      const Numeric r             = ne * c1 / f2;
      propagation_matrix[iv].U() += r;

      for (Size i = 0; i < 3; i++) {
        if (jacs[i] != end) {
          propagation_matrix_jacobian[jacs[i]->target_pos, iv].U() +=
              ne * dc1[i] / f2;
        }
      }

      for (Size i = 3; i < 6; i++) {
        if (jacs[i] != end) {
          propagation_matrix_jacobian[jacs[i]->target_pos, iv].U() +=
              -2.0 * ne * r / freq_grid[iv];
        }
      }

      if (jacs[6] != end) {
        propagation_matrix_jacobian[jacs[6]->target_pos, iv].U() += r;
      }
    }
  }
}

void propagation_matrix_agendaAuto(Agenda& propagation_matrix_agenda,
                                   const ArrayOfSpeciesTag& abs_species,
                                   const AbsorptionBands& abs_bands,
                                   const Index& use_abs_lookup_data,
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
  ARTS_TIME_REPORT

  AgendaCreator agenda("propagation_matrix_agenda");

  const SpeciesTagTypeStatus any_species(abs_species);

  // propagation_matrixInit
  agenda.add("propagation_matrixInit");

  // propagation_matrixAddLines or propagation_matrixAddLookup
  if (use_abs_lookup_data) {
    agenda.add("propagation_matrixAddLookup",
               SetWsv{"no_negative_absorption", no_negative_absorption},
               SetWsv{"p_interp_order", p_interp_order},
               SetWsv{"t_interp_order", t_interp_order},
               SetWsv{"water_interp_order", water_interp_order},
               SetWsv{"f_interp_order", f_interp_order},
               SetWsv{"extpolfac", extpolfac});
  } else if (abs_bands.size()) {
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
  if (std::ranges::any_of(abs_species, [](auto& spec) {
        return spec.Spec() == "free_electrons"_spec;
      })) {
    agenda.add("propagation_matrixAddFaraday");
  }

  agenda.add("propagation_matrix_jacobianWindFix");

  // Extra check (should really never ever fail when species exist)
  propagation_matrix_agenda = std::move(agenda).finalize(true);
}
