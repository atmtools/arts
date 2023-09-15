/**
  @file   m_rte.cc
  @author Patrick Eriksson <patrick.eriksson@chalmers.se>
  @date   2002-05-11

  @brief  Workspace methods for solving clear sky radiative transfer.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts_constants.h"
#include "arts_omp.h"
#include "atm.h"
#include "atm_path.h"
#include <workspace.h>
#include "check_input.h"
#include "debug.h"
#include "geodetic.h"
#include "gridded_fields.h"
#include "jacobian.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpack_arrays.h"
#include "matpack_data.h"
#include "montecarlo.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "rtepack.h"
#include "special_interp.h"
#include "species_tags.h"
#include "sun.h"
#include "surf.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
#include <stdexcept>

inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;

/*===========================================================================
  === Workspace methods
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void iyApplyUnit(Matrix& iy,
                 ArrayOfMatrix& iy_aux,
                 const Vector& f_grid,
                 const ArrayOfString& iy_aux_vars,
                 const String& iy_unit) {
  ARTS_USER_ERROR_IF (iy_unit == "1",
    "No need to use this method with *iy_unit* = \"1\".");

  ARTS_USER_ERROR_IF (max(iy(joker, 0)) > 1e-3,
      "The spectrum matrix *iy* is required to have original radiance\n"
      "unit, but this seems not to be the case. This as a value above\n"
      "1e-3 is found in *iy*.")

  // Polarisation index variable
  ArrayOfIndex i_pol(4);
  for (Index is = 0; is < 4; is++) {
    i_pol[is] = is + 1;
  }

  apply_iy_unit(iy, iy_unit, f_grid, 1, i_pol);

  for (Size i = 0; i < iy_aux_vars.size(); i++) {
    if (iy_aux_vars[i] == "iy" || iy_aux_vars[i] == "Error" ||
        iy_aux_vars[i] == "Error (uncorrelated)") {
      apply_iy_unit(iy_aux[i], iy_unit, f_grid, 1, i_pol);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iyCalc(const Workspace& ws,
            Matrix& iy,
            ArrayOfMatrix& iy_aux,
            Ppath& ppath,
            Vector& geo_pos,
            const Index& atmfields_checked,
            const Index& atmgeom_checked,
            const ArrayOfString& iy_aux_vars,
            const Index& iy_id,
            const Index& cloudbox_on,
            const Index& cloudbox_checked,
            const Index& scat_data_checked,
            const Vector& f_grid,
            const AtmField& atm_field,
            const Vector& rte_pos,
            const Vector& rte_los,
            const Vector& rte_pos2,
            const String& iy_unit,
            const Agenda& iy_main_agenda) {
  // Basics
  //
  ARTS_USER_ERROR_IF (atmfields_checked != 1,
        "The atmospheric fields must be flagged to have\n"
        "passed a consistency check (atmfields_checked=1).");
  ARTS_USER_ERROR_IF (atmgeom_checked != 1,
        "The atmospheric geometry must be flagged to have\n"
        "passed a consistency check (atmgeom_checked=1).");
  ARTS_USER_ERROR_IF (cloudbox_checked != 1,
        "The cloudbox must be flagged to have\n"
        "passed a consistency check (cloudbox_checked=1).");
  if (cloudbox_on)
    ARTS_USER_ERROR_IF (scat_data_checked != 1,
          "The scattering data must be flagged to have\n"
          "passed a consistency check (scat_data_checked=1).");

  // iy_transmittance is just input and can be left empty for first call
  Tensor3 iy_transmittance(0, 0, 0);
  ArrayOfTensor3 diy_dx;

  iy_main_agendaExecute(ws,
                        iy,
                        iy_aux,
                        ppath,
                        diy_dx,
                        geo_pos,
                        1,
                        iy_transmittance,
                        iy_aux_vars,
                        iy_id,
                        iy_unit,
                        cloudbox_on,
                        0,
                        f_grid,
                        atm_field,
                        rte_pos,
                        rte_los,
                        rte_pos2,
                        iy_main_agenda);

  // Don't allow NaNs (should suffice to check first stokes element)
  for (Index i = 0; i < iy.nrows(); i++) {
    ARTS_USER_ERROR_IF (std::isnan(iy(i, 0)),
                        "One or several NaNs found in *iy*.");
  }
}

void ppvar_atmFromPath(ArrayOfAtmPoint &ppvar_atm, const Ppath &ppath,
                       const AtmField &atm_field) try {
  forward_atm_path(atm_path_resize(ppvar_atm, ppath), ppath, atm_field);
} ARTS_METHOD_ERROR_CATCH

void ppvar_fFromPath(ArrayOfVector &ppvar_f, const Vector &f_grid,
                     const Ppath &ppath, const ArrayOfAtmPoint &ppvar_atm,
                     const Numeric &rte_alonglos_v) try {
  forward_path_freq(path_freq_resize(ppvar_f, f_grid, ppvar_atm), f_grid, ppath,
                    ppvar_atm, rte_alonglos_v);
} ARTS_METHOD_ERROR_CATCH

void ppvar_radCalcEmission(
    ArrayOfStokvecVector &ppvar_rad, ArrayOfStokvecMatrix &ppvar_drad,
    const StokvecVector &background_rad, const ArrayOfStokvecVector &ppvar_src,
    const ArrayOfStokvecMatrix &ppvar_dsrc,
    const ArrayOfMuelmatVector &ppvar_tramat,
    const ArrayOfMuelmatVector &ppvar_cumtramat,
    const ArrayOfArrayOfMuelmatMatrix &ppvar_dtramat) try {
  const Size np = ppvar_src.size();

  ARTS_USER_ERROR_IF(np not_eq ppvar_dsrc.size(),
                     "ppvar_dsrc must have (np) elements")
  ARTS_USER_ERROR_IF(np not_eq ppvar_tramat.size(),
                     "ppvar_tramat must have (np) elements")
  ARTS_USER_ERROR_IF(np not_eq ppvar_cumtramat.size(),
                     "ppvar_cumtramat must have (np) elements")
  ARTS_USER_ERROR_IF(2 not_eq ppvar_dtramat.size() or
                         ppvar_dtramat.front().size() not_eq
                             ppvar_dtramat.back().size() or
                         ppvar_dtramat.front().size() not_eq np,
                     "ppvar_dtramat must (2 x np) elements")

  if (np == 0) {
    ppvar_rad.resize(0);
    ppvar_drad.resize(0);
    return;
  }

  const Index nq = ppvar_dsrc.front().nrows();
  const Index nf = ppvar_dsrc.front().ncols();

  const auto test_nf = [nf](auto &v) { return v.size() not_eq nf; };
  ARTS_USER_ERROR_IF(nf not_eq background_rad.size(),
                     "background_rad must have nf elements")
  ARTS_USER_ERROR_IF(std::any_of(ppvar_src.begin(), ppvar_src.end(), test_nf),
                     "ppvar_src must have (nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ppvar_tramat.begin(), ppvar_tramat.end(), test_nf),
      "ppvar_tramat must have (nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ppvar_cumtramat.begin(), ppvar_cumtramat.end(), test_nf),
      "ppvar_src must have (nf) inner elements")

  const auto test_nfnq = [nf, nq](auto &v) {
    return v.ncols() not_eq nf or v.nrows() not_eq nq;
  };
  ARTS_USER_ERROR_IF(
      std::any_of(ppvar_dsrc.begin(), ppvar_dsrc.end(), test_nfnq),
      "ppvar_dsrc must have (nq x nf) inner elements")
  ARTS_USER_ERROR_IF(std::any_of(ppvar_dtramat.front().begin(),
                                 ppvar_dtramat.front().end(), test_nfnq),
                     "ppvar_dtramat must have (nq x nf) inner elements")
  ARTS_USER_ERROR_IF(std::any_of(ppvar_dtramat.back().begin(),
                                 ppvar_dtramat.back().end(), test_nfnq),
                     "ppvar_dtramat must have (nq x nf) inner elements")

  ppvar_rad.resize(np, background_rad);
  ppvar_drad.resize(np, StokvecMatrix(nq, nf));
  for (Index ip = np - 2; ip >= 0; ip--) {
    ppvar_rad[ip] = ppvar_rad[ip + 1];
    two_level_linear_emission_step(
        ppvar_rad[ip], ppvar_drad[ip], ppvar_drad[ip + 1], ppvar_src[ip],
        ppvar_src[ip + 1], ppvar_dsrc[ip], ppvar_dsrc[ip + 1],
        ppvar_tramat[ip + 1], ppvar_cumtramat[ip], ppvar_dtramat[0][ip + 1],
        ppvar_dtramat[1][ip + 1]);
  }
}
ARTS_METHOD_ERROR_CATCH

void ppvar_radCalcTransmission(
    ArrayOfStokvecVector &ppvar_rad,
    ArrayOfStokvecMatrix &ppvar_drad,
    const ArrayOfMuelmatVector &ppvar_tramat,
    const ArrayOfMuelmatVector &ppvar_cumtramat,
    const ArrayOfArrayOfMuelmatMatrix &ppvar_dtramat) try {
  const Size np = ppvar_tramat.size();

  ARTS_USER_ERROR_IF(np not_eq ppvar_tramat.size(),
                     "ppvar_tramat must have (np) elements")
  ARTS_USER_ERROR_IF(np not_eq ppvar_cumtramat.size(),
                     "ppvar_cumtramat must have (np) elements")
  ARTS_USER_ERROR_IF(2 not_eq ppvar_dtramat.size() or
                         ppvar_dtramat.front().size() not_eq
                             ppvar_dtramat.back().size() or
                         ppvar_dtramat.front().size() not_eq np,
                     "ppvar_dtramat must (2 x np) elements")

  if (np == 0) {
    ppvar_rad.resize(0);
    ppvar_drad.resize(0);
    return;
  }

  const Index nq = ppvar_dtramat.front().front().nrows();
  const Index nf = ppvar_dtramat.front().front().ncols();

  const auto test_nf = [nf](auto &v) { return v.size() not_eq nf; };
  ARTS_USER_ERROR_IF(
      std::any_of(ppvar_tramat.begin(), ppvar_tramat.end(), test_nf),
      "ppvar_tramat must have (nf) inner elements")
  ARTS_USER_ERROR_IF(
      std::any_of(ppvar_cumtramat.begin(), ppvar_cumtramat.end(), test_nf),
      "ppvar_src must have (nf) inner elements")

  const auto test_nfnq = [nf, nq](auto &v) {
    return v.ncols() not_eq nf or v.nrows() not_eq nq;
  };
  ARTS_USER_ERROR_IF(std::any_of(ppvar_dtramat.front().begin(),
                                 ppvar_dtramat.front().end(), test_nfnq),
                     "ppvar_dtramat must have (nq x nf) inner elements")
  ARTS_USER_ERROR_IF(std::any_of(ppvar_dtramat.back().begin(),
                                 ppvar_dtramat.back().end(), test_nfnq),
                     "ppvar_dtramat must have (nq x nf) inner elements")

  ppvar_rad.resize(np, StokvecVector(nf, Stokvec{1, 0, 0, 0}));
  ppvar_drad.resize(np, StokvecMatrix(nq, nf));
  for (Index ip = np - 2; ip >= 0; ip--) {
    ppvar_rad[ip] = ppvar_rad[ip + 1];
    two_level_linear_transmission_step(
        ppvar_rad[ip], ppvar_drad[ip], ppvar_drad[ip + 1], ppvar_tramat[ip + 1],
        ppvar_cumtramat[ip], ppvar_dtramat[0][ip + 1],
        ppvar_dtramat[1][ip + 1]);
  }
}
ARTS_METHOD_ERROR_CATCH

void iyUnitConversion(Matrix &iy, ArrayOfTensor3 &diy_dx,
                      Tensor3 &ppvar_iy, const Vector &f_grid,
                      const Ppath &ppath,
                      const ArrayOfRetrievalQuantity &jacobian_quantities,
                      const String &iy_unit, const Index &jacobian_do,
                      const Index &iy_agenda_call1) {
  // Radiance unit conversions
  if (iy_agenda_call1) {
    rtmethods_unit_conversion(
        iy, diy_dx, ppvar_iy, f_grid, ppath, jacobian_quantities,
        jacobian_do ? do_analytical_jacobian<2>(jacobian_quantities) : 0,
        iy_unit);
  }
}

void iy_transmittance_backgroundFromRte(Tensor3 &iy_transmittance_background,
                                        const Tensor3 &iy_transmittance,
                                        const MuelmatVector &cumtramat,
                                        const Index &iy_agenda_call1) {
  if (iy_agenda_call1) {
    iy_transmittance_background = to_tensor3(cumtramat);
  } else {
    if (iy_transmittance_background.data_handle() ==
        iy_transmittance.data_handle()) {

    } else {
      iy_transmittance_mult(iy_transmittance_background, iy_transmittance,
                            to_tensor3(cumtramat));
    }
  }
}

void ppvar_propmatCalc(const Workspace &ws,
                       ArrayOfPropmatVector &ppvar_propmat,
                       ArrayOfStokvecVector &ppvar_nlte,
                       ArrayOfPropmatMatrix &ppvar_dpropmat,
                       ArrayOfStokvecMatrix &ppvar_dnlte,
                       const Agenda &propmat_clearsky_agenda,
                       const ArrayOfRetrievalQuantity &jacobian_quantities,
                       const ArrayOfVector &ppvar_f, const Ppath &ppath,
                       const ArrayOfAtmPoint &ppvar_atm,
                       const Index &jacobian_do) try {
  const Index j_analytical_do =
      jacobian_do ? do_analytical_jacobian<2>(jacobian_quantities) : 0;

  ArrayOfString fail_msg;
  bool do_abort = false;

  const Index np = ppath.np;
  if (np == 0) {
    ppvar_propmat.resize(0);
    ppvar_nlte.resize(0);
    ppvar_dpropmat.resize(0);
    ppvar_dnlte.resize(0);
    return;
  }

  ppvar_propmat.resize(np);
  ppvar_nlte.resize(np);
  ppvar_dpropmat.resize(np);
  ppvar_dnlte.resize(np);

  // Loop ppath points and determine radiative properties
#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Index ip = 0; ip < np; ip++) {
    if (do_abort)
      continue;
    try {
      get_stepwise_clearsky_propmat(
          ws, ppvar_propmat[ip], ppvar_nlte[ip], ppvar_dpropmat[ip],
          ppvar_dnlte[ip], propmat_clearsky_agenda, jacobian_quantities,
          ppvar_f[ip], Vector{ppath.los(ip, joker)}, ppvar_atm[ip],
          j_analytical_do);
    } catch (const std::runtime_error &e) {
#pragma omp critical(iyEmissionStandard_source)
      {
        do_abort = true;
        fail_msg.push_back(var_string("Runtime-error in propagation radiative "
                                      "properties calculation at index ",
                                      ip, ": \n", e.what()));
      }
    }
  }

  ARTS_USER_ERROR_IF(do_abort, "Error messages from failed cases:\n", fail_msg)
} ARTS_METHOD_ERROR_CATCH

void ppvar_srcFromPropmat(ArrayOfStokvecVector &ppvar_src,
                          ArrayOfStokvecMatrix &ppvar_dsrc,
                          const ArrayOfPropmatVector &ppvar_propmat,
                          const ArrayOfStokvecVector &ppvar_nlte,
                          const ArrayOfPropmatMatrix &ppvar_dpropmat,
                          const ArrayOfStokvecMatrix &ppvar_dnlte,
                          const ArrayOfVector &ppvar_f,
                          const ArrayOfAtmPoint &ppvar_atm,
                          const ArrayOfRetrievalQuantity &jacobian_quantities,
                          const Index &jacobian_do) try {
  const Index j_analytical_do =
      jacobian_do ? do_analytical_jacobian<2>(jacobian_quantities) : 0;

  ArrayOfString fail_msg;
  bool do_abort = false;

  const Index np = ppvar_atm.size();
  if (np == 0) {
    ppvar_src.resize(0);
    ppvar_dsrc.resize(0);
    return;
  }

  const Index nf = ppvar_propmat.front().size();
  const Index nq = j_analytical_do ? jacobian_quantities.size() : 0;

  ppvar_src.resize(np, StokvecVector(nf));
  ppvar_dsrc.resize(np, StokvecMatrix(nq, nf));

  Vector B(nf);
  Matrix dB(nq, nf);

  // Loop ppath points and determine radiative properties
#pragma omp parallel for if (!arts_omp_in_parallel())                          \
    firstprivate(B, dB)
  for (Index ip = 0; ip < np; ip++) {
    if (do_abort)
      continue;
    try {
      get_stepwise_blackbody_radiation(B, dB, ppvar_f[ip],
                                       ppvar_atm[ip].temperature,
                                       jacobian_quantities, j_analytical_do);

      rtepack::source::level_nlte(ppvar_src[ip], ppvar_dsrc[ip],
                                  ppvar_propmat[ip], ppvar_nlte[ip],
                                  ppvar_dpropmat[ip], ppvar_dnlte[ip], B, dB);
    } catch (const std::runtime_error &e) {
#pragma omp critical(iyEmissionStandard_source)
      {
        do_abort = true;
        fail_msg.push_back(
            var_string("Runtime-error in source calculation at index ", ip,
                       ": \n", e.what()));
      }
    }
  }
} ARTS_METHOD_ERROR_CATCH

void ppvar_tramatCalc(ArrayOfMuelmatVector &ppvar_tramat,
                      ArrayOfArrayOfMuelmatMatrix &ppvar_dtramat,
                      Vector &ppvar_distance,
                      ArrayOfArrayOfVector &ppvar_ddistance,
                      const ArrayOfPropmatVector &ppvar_propmat,
                      const ArrayOfPropmatMatrix &ppvar_dpropmat,
                      const Ppath &ppath, const ArrayOfAtmPoint &ppvar_atm,
                      const ArrayOfRetrievalQuantity &jacobian_quantities,
                      const Index &jacobian_do) try {
  ArrayOfString fail_msg;
  bool do_abort = false;

  const Index j_analytical_do =
      jacobian_do ? do_analytical_jacobian<2>(jacobian_quantities) : 0;

  // HSE variables
  const Index temperature_derivative_position =
      j_analytical_do ? pos<Jacobian::Atm::Temperature>(jacobian_quantities)
                      : -1;

  const Index np = ppath.np;

  if (np == 0) {
    ppvar_tramat.resize(0);
    ppvar_dtramat.resize(2, ArrayOfMuelmatMatrix{});
    ppvar_distance.resize(0);
    ppvar_ddistance.resize(2, ArrayOfVector{});
    return;
  }

  const Index nf = ppvar_propmat.front().size();
  const Index nq = j_analytical_do ? jacobian_quantities.size() : 0;

  ppvar_tramat.resize(np, MuelmatVector(nf));
  ppvar_dtramat.resize(2, ArrayOfMuelmatMatrix(np, MuelmatMatrix(nq, nf)));
  ppvar_distance.resize(np);
  ppvar_ddistance.resize(2, ArrayOfVector(np, Vector(nq, 0)));

#pragma omp parallel for if (!arts_omp_in_parallel())
  for (Index ip = 1; ip < np; ip++) {
    if (do_abort)
      continue;
    try {

      ppvar_distance[ip - 1] = ppath.lstep[ip - 1];
      if (temperature_derivative_position >= 0 and
          jacobian_quantities[temperature_derivative_position].Subtag() ==
              "HSE on") {
        ppvar_ddistance[0][ip][temperature_derivative_position] =
            ppath.lstep[ip - 1] / (2.0 * ppvar_atm[ip - 1].temperature);
        ppvar_ddistance[1][ip][temperature_derivative_position] =
            ppath.lstep[ip - 1] / (2.0 * ppvar_atm[ip].temperature);
      }

      two_level_exp(ppvar_tramat[ip], ppvar_dtramat[0][ip],
                    ppvar_dtramat[1][ip], ppvar_propmat[ip - 1],
                    ppvar_propmat[ip], ppvar_dpropmat[ip - 1],
                    ppvar_dpropmat[ip], ppvar_distance[ip - 1],
                    ppvar_ddistance[0][ip], ppvar_ddistance[1][ip]);
    } catch (const std::runtime_error &e) {
#pragma omp critical(iyEmissionStandard_transmission)
      {
        do_abort = true;
        fail_msg.push_back(
            var_string("Runtime-error in transmission calculation at index ",
                       ip, ": \n", e.what()));
      }
    }
  }

  ARTS_USER_ERROR_IF(do_abort, "Error messages from failed cases:\n", fail_msg)
} ARTS_METHOD_ERROR_CATCH

void iy_auxFromVars(ArrayOfMatrix &iy_aux, const ArrayOfString &iy_aux_vars,
                    const MuelmatVector &background_transmittance,
                    const Ppath &ppath) {
  const Index np = ppath.np;

  // Init iy_aux and fill where possible
  const Index naux = iy_aux_vars.size();
  iy_aux.resize(naux);

  if (np == 0) {
    return;
  }

  const Index nf = background_transmittance.size();

  //
  Index auxOptDepth = -1;
  //
  for (Index i = 0; i < naux; i++) {
    iy_aux[i].resize(nf, 4);
    iy_aux[i] = 0;

    if (iy_aux_vars[i] == "Optical depth")
      auxOptDepth = i;
    else {
      ARTS_USER_ERROR("The only allowed strings in *iy_aux_vars* are:\n"
                      "  \"Radiative background\"\n"
                      "  \"Optical depth\"\n"
                      "but you have selected: \"",
                      iy_aux_vars[i], "\"")
    }
  }

  // iy_aux: Optical depth
  if (auxOptDepth >= 0)
    for (Index iv = 0; iv < nf; iv++)
      iy_aux[auxOptDepth](iv, 0) = -std::log(background_transmittance[iv](0, 0));
}

void iyCopyPath(Matrix &iy, Tensor3 &ppvar_iy, Tensor4 &ppvar_trans_cumulat, Tensor4 &ppvar_trans_partial,
             ArrayOfTensor3 &diy_dpath, 
             const ArrayOfStokvecVector &ppvar_rad,
             const ArrayOfStokvecMatrix &ppvar_drad,
             const ArrayOfMuelmatVector &ppvar_cumtramat,
             const ArrayOfMuelmatVector &ppvar_tramat,
             const ArrayOfRetrievalQuantity &jacobian_quantities,
             const Index &jacobian_do) {
  const Index j_analytical_do =
      jacobian_do ? do_analytical_jacobian<2>(jacobian_quantities) : 0;

  const Index np = ppvar_rad.size();
  if (np == 0) {
    ppvar_trans_cumulat.resize(np, 0, 0, 0);
    ppvar_trans_partial.resize(np, 0, 0, 0);
    ppvar_iy.resize(0, 0, np);
    iy.resize(0, 0);
    return;
  }

  const Index nf = ppvar_rad.front().size();

  // Copy back to ARTS external style
  diy_dpath = j_analytical_do ? get_standard_diy_dpath(jacobian_quantities, np,
                                                       nf, false)
                              : ArrayOfTensor3(0);
  ppvar_trans_cumulat.resize(np, nf, 4, 4);
  ppvar_trans_partial.resize(np, nf, 4, 4);
  ppvar_iy.resize(nf, 4, np);
  iy = to_matrix(ppvar_rad.front());
  for (Size ip = 0; ip < ppvar_rad.size(); ip++) {
    ppvar_trans_cumulat(ip, joker, joker, joker) = to_tensor3(ppvar_cumtramat[ip]);
    ppvar_trans_partial(ip, joker, joker, joker) = to_tensor3(ppvar_tramat[ip]);
    ppvar_iy(joker, joker, ip) = to_matrix(ppvar_rad[ip]);
    if (j_analytical_do)
      FOR_ANALYTICAL_JACOBIANS_DO(diy_dpath[iq](ip, joker, joker) =
                                      to_matrix(ppvar_drad[ip][iq]););
  }
}

void diy_dxTransform(const Workspace &ws, ArrayOfTensor3 &diy_dx,
                     ArrayOfTensor3 &diy_dpath, const Ppath &ppath,
                     const ArrayOfAtmPoint &ppvar_atm,
                     const ArrayOfArrayOfSpeciesTag &abs_species,
                     const Tensor3 &iy_transmittance,
                     const Agenda &water_p_eq_agenda,
                     const ArrayOfRetrievalQuantity &jacobian_quantities,
                     const Index &jacobian_do, const Index &iy_agenda_call1) {
  const Index j_analytical_do =
      jacobian_do ? do_analytical_jacobian<2>(jacobian_quantities) : 0;

  const ArrayOfIndex jac_species_i =
      j_analytical_do ? get_pointers_for_analytical_species(jacobian_quantities,
                                                            abs_species)
                      : ArrayOfIndex(0);

  if (const Index np = ppath.np; np not_eq 0 and j_analytical_do) {
    const Index nf = diy_dpath.front().nrows();

    rtmethods_jacobian_finalisation(
        ws, diy_dx, diy_dpath, nf, np, ppath, ppvar_atm, iy_agenda_call1,
        iy_transmittance, water_p_eq_agenda, jacobian_quantities, abs_species,
        jac_species_i);
  }
}

void iyBackground(const Workspace &ws, Matrix &iy, ArrayOfTensor3 &diy_dx,
                   const Tensor3 &iy_transmittance,
                   const MuelmatVector &total_transmittance,
                   const SurfaceField &surface_field, const Vector &f_grid,
                   const Vector &rte_pos2, const Ppath &ppath,
                   const AtmField &atm_field,
                   const ArrayOfRetrievalQuantity &jacobian_quantities,
                   const Index &jacobian_do, const Index &cloudbox_on,
                   const Index &iy_id, const Index &iy_agenda_call1,
                   const Agenda &iy_main_agenda, const Agenda &iy_space_agenda,
                   const Agenda &iy_surface_agenda,
                   const Agenda &iy_cloudbox_agenda, const String &iy_unit) try {
  // iy_transmittance
  Tensor3 iy_transmittance_background;
  iy_transmittance_backgroundFromRte(iy_transmittance_background,
                                     iy_transmittance, total_transmittance,
                                     iy_agenda_call1);

  const Index nf = f_grid.size();
  const Index np = ppath.np;

  const Index j_analytical_do =
      jacobian_do ? do_analytical_jacobian<2>(jacobian_quantities) : 0;

  if (j_analytical_do and iy_agenda_call1)
    diy_dx =
        get_standard_starting_diy_dx(jacobian_quantities, np, nf, false);

  get_iy_of_background(
      ws, iy, diy_dx, iy_transmittance_background, iy_id, jacobian_do,
      jacobian_quantities, ppath, rte_pos2, atm_field, cloudbox_on, f_grid,
      iy_unit, surface_field, iy_main_agenda, iy_space_agenda,
      iy_surface_agenda, iy_cloudbox_agenda, iy_agenda_call1);
} ARTS_METHOD_ERROR_CATCH

void background_radFromMatrix(StokvecVector &background_rad, const Matrix &iy) {
  ARTS_USER_ERROR_IF(iy.ncols() > 5 or iy.ncols() < 1,
                     "Only for stokes dimensions [1, 4].")
  background_rad = rtepack::to_stokvec_vector(iy);
}

void background_transmittanceFromBack(
    MuelmatVector &background_transmittance,
    const ArrayOfMuelmatVector &ppvar_cumtramat) {
  ARTS_USER_ERROR_IF(ppvar_cumtramat.size() == 0, "Cannot extract from empty list.")
  background_transmittance = ppvar_cumtramat.back();
}

void background_transmittanceFromFront(
    MuelmatVector &background_transmittance,
    const ArrayOfMuelmatVector &ppvar_cumtramat) {
  ARTS_USER_ERROR_IF(ppvar_cumtramat.size() == 0, "Cannot extract from empty list.")
  background_transmittance = ppvar_cumtramat.front();
}

void ppvar_cumtramatForward(ArrayOfMuelmatVector &ppvar_cumtramat,
                            const ArrayOfMuelmatVector &ppvar_tramat) {
  ppvar_cumtramat = forward_cumulative_transmission(ppvar_tramat);
}

void ppvar_cumtramatReverse(ArrayOfMuelmatVector &ppvar_cumtramat,
                            const ArrayOfMuelmatVector &ppvar_tramat) {
  ppvar_cumtramat = reverse_cumulative_transmission(ppvar_tramat);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iyLoopFrequencies(const Workspace& ws,
                       Matrix& iy,
                       ArrayOfMatrix& iy_aux,
                       Ppath& ppath,
                       ArrayOfTensor3& diy_dx,
                       const ArrayOfString& iy_aux_vars,
                       const Index& iy_agenda_call1,
                       const Tensor3& iy_transmittance,
                       const Vector& rte_pos,
                       const Vector& rte_los,
                       const Vector& rte_pos2,
                       const Vector& f_grid,
                       const Agenda& iy_loop_freqs_agenda) {
  // Throw error if unsupported features are requested
  ARTS_USER_ERROR_IF (!iy_agenda_call1,
        "Recursive usage not possible (iy_agenda_call1 must be 1).");
  ARTS_USER_ERROR_IF (iy_transmittance.ncols(),
        "*iy_transmittance* must be empty.");

  const Index nf = f_grid.size();

  for (Index i = 0; i < nf; i++) {
    // Variables for 1 frequency
    Matrix iy1;
    ArrayOfMatrix iy_aux1;
    ArrayOfTensor3 diy_dx1;

    iy_loop_freqs_agendaExecute(ws,
                                iy1,
                                iy_aux1,
                                ppath,
                                diy_dx1,
                                iy_agenda_call1,
                                iy_transmittance,
                                iy_aux_vars,
                                0,
                                Vector(1, f_grid[i]),
                                rte_pos,
                                rte_los,
                                rte_pos2,
                                iy_loop_freqs_agenda);

    // After first frequency, give output its size
    if (i == 0) {
      iy.resize(nf, 4);
      //
      iy_aux.resize(iy_aux1.size());
      for (Size q = 0; q < iy_aux1.size(); q++) {
        iy_aux[q].resize(nf, 4);
      }
      //
      diy_dx.resize(diy_dx1.size());
      for (Size q = 0; q < diy_dx1.size(); q++) {
        diy_dx[q].resize(diy_dx1[q].npages(), nf, 4);
      }
    }

    // Copy to output variables
    iy(i, joker) = iy1(0, joker);
    for (Size q = 0; q < iy_aux1.size(); q++) {
      iy_aux[q](i, joker) = iy_aux1[q](0, joker);
    }
    for (Size q = 0; q < diy_dx1.size(); q++) {
      diy_dx[q](joker, i, joker) = diy_dx1[q](joker, 0, joker);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iyReplaceFromAux(Matrix& iy,
                      const ArrayOfMatrix& iy_aux,
                      const ArrayOfString& iy_aux_vars,
                      const Index& jacobian_do,
                      const String& aux_var) {
  ARTS_USER_ERROR_IF (iy_aux.size() != iy_aux_vars.size(),
        "*iy_aux* and *iy_aux_vars* must have the same "
        "number of elements.");

  ARTS_USER_ERROR_IF (jacobian_do,
        "This method can not provide any jacobians and "
        "*jacobian_do* must be 0.");

  bool ready = false;

  for (Size i = 0; i < iy_aux.size() && !ready; i++) {
    if (iy_aux_vars[i] == aux_var) {
      iy = iy_aux[i];
      ready = true;
    }
  }

  ARTS_USER_ERROR_IF (!ready,
        "The selected auxiliary variable to insert in *iy* "
        "is either not defined at all or is not set.");
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ppvar_optical_depthFromPpvar_trans_cumulat(
    Matrix& ppvar_optical_depth,
    const Tensor4& ppvar_trans_cumulat) {
  ppvar_optical_depth = ppvar_trans_cumulat(joker, joker, 0, 0);
  transform(ppvar_optical_depth, log, ppvar_optical_depth);
  ppvar_optical_depth *= -1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void yApplyUnit(Vector& y,
                Matrix& jacobian,
                const Vector& y_f,
                const ArrayOfIndex& y_pol,
                const String& iy_unit) {
  ARTS_USER_ERROR_IF (iy_unit == "1",
                      "No need to use this method with *iy_unit* = \"1\".");

  ARTS_USER_ERROR_IF (max(y) > 1e-3,
      "The spectrum vector *y* is required to have original radiance\n"
      "unit, but this seems not to be the case. This as a value above\n"
      "1e-3 is found in *y*.")

  // Is jacobian set?
  //
  const Index ny = y.size();
  //
  const bool do_j = jacobian.nrows() == ny;

  // Some jacobian quantities can not be handled
  ARTS_USER_ERROR_IF (do_j && max(jacobian) > 1e-3,
      "The method can not be used with jacobian quantities that are not\n"
      "obtained through radiative transfer calculations. One example on\n"
      "quantity that can not be handled is *jacobianAddPolyfit*.\n"
      "The maximum value of *jacobian* indicates that one or several\n"
      "such jacobian quantities are included.")

  // Planck-Tb
  //--------------------------------------------------------------------------
  if (iy_unit == "PlanckBT") {
    // Hard to use telescoping here as the data are sorted differently in y
    // and jacobian, than what is expected apply_iy_unit. Copy to temporary
    // variables instead.

    // Handle the elements in "frequency chunks"

    Index i0 = 0;  // Index of first element for present chunk
    //
    while (i0 < ny) {
      // Find number of values for this chunk
      Index n = 1;
      //
      while (i0 + n < ny && y_f[i0] == y_f[i0 + n]) {
        n++;
      }

      Matrix yv(1, n);
      ArrayOfIndex i_pol(n);
      bool any_quv = false;
      //
      for (Index i = 0; i < n; i++) {
        const Index ix = i0 + i;
        yv(0, i) = y[ix];
        i_pol[i] = y_pol[ix];
        if (i_pol[i] > 1 && i_pol[i] < 5) {
          any_quv = true;
        }
      }

      // Index of elements to convert
      Range ii(i0, n);

      if (do_j) {
        ARTS_USER_ERROR_IF (any_quv && i_pol[0] != 1,
            "The conversion to PlanckBT, of the Jacobian and "
            "errors for Q, U and V, requires that I (first Stokes "
            "element) is at hand and that the data are sorted in "
            "such way that I comes first for each frequency.")

        // Jacobian
        Tensor3 J(jacobian.ncols(), 1, n);
        J(joker, 0, joker) = transpose(jacobian(ii, joker));
        apply_iy_unit2(J, yv, iy_unit, ExhaustiveConstVectorView{y_f[i0]}, 1, i_pol);
        jacobian(ii, joker) = transpose(J(joker, 0, joker));
      }

      // y (must be done last)
      apply_iy_unit(yv, iy_unit, ExhaustiveConstVectorView{y_f[i0]}, 1, i_pol);
      y[ii] = yv(0, joker);

      i0 += n;
    }
  }

  // Other conversions
  //--------------------------------------------------------------------------
  else {
    // Here we take each element of y separately.

    Matrix yv(1, 1);
    ArrayOfIndex i_pol(1);

    for (Index i = 0; i < ny; i++) {
      yv(0, 0) = y[i];
      i_pol[0] = y_pol[i];

      // Jacobian
      if (do_j) {
        apply_iy_unit2(
            ExhaustiveTensor3View{ExhaustiveMatrixView{jacobian(i, joker)}}, yv, iy_unit, ExhaustiveConstVectorView{y_f[i]}, 1, i_pol);
      }

      // y (must be done last)
      apply_iy_unit(yv, iy_unit, ExhaustiveConstVectorView{y_f[i]}, 1, i_pol);
      y[i] = yv(0, 0);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void y_geo_seriesFromY_geo(Matrix& y_geo_series,
                           const Matrix& y_geo,
                           const Vector& sensor_response_f_grid)
{
  // Sizes
  const Index ly = y_geo.nrows();
  const Index nchannel = sensor_response_f_grid.size();
  const Index lseries = ly / nchannel;
  
  ARTS_USER_ERROR_IF (nchannel * lseries != ly,
    "Row size of *y_geo* not an even multiple of length of *sensor_response_f_grid*.")

  y_geo_series.resize(lseries, y_geo.ncols());

  Index i = 0;
  for (Index s=0; s<lseries; ++s) {
    y_geo_series(s, joker) = y_geo(i, joker);
    i += nchannel;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void y_geo_swathFromY_geo(Tensor3& y_geo_swath,
                          const Matrix& y_geo,
                          const Vector& sensor_response_f_grid,
                          const Index& npixel)
{
  // Sizes
  const Index ly = y_geo.nrows();
  const Index nchannel = sensor_response_f_grid.size();
  const Index nswath = ly / (nchannel * npixel);
  
  ARTS_USER_ERROR_IF (nchannel * npixel * nswath != ly,
    "Row size of *y_geo* does not match given *npixel* and *sensor_response_f_grid*.")

  y_geo_swath.resize(nswath, npixel, y_geo.ncols());

  Index i = 0;
  for (Index s=0; s<nswath; ++s) {
    for (Index p=0; p<npixel; ++p) {
      y_geo_swath(s, p, joker) = y_geo(i, joker);
      i += nchannel;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void y_seriesFromY(Matrix& y_series,
                   const Vector& y,
                   const Vector& y_f,
                   const Vector& sensor_response_f_grid,
                   const Index& safe)
{
  // Sizes
  const Index ly = y.size();
  const Index nchannel = sensor_response_f_grid.size();
  const Index lseries = ly / nchannel;
  
  ARTS_USER_ERROR_IF (nchannel * lseries != ly,
    "Length of *y* not an even multiple of length of *sensor_response_f_grid*.")

  y_series.resize(lseries, nchannel);

  Index i = 0;
  for (Index s=0; s<lseries; ++s) {
    for (Index c=0; c<nchannel; ++c) {
      if (safe && s > 0) {
        ARTS_USER_ERROR_IF (fabs(y_f[i] - y_f[i-nchannel]) > 1,
                          "At least one channel varies in frequency.")
      }
      y_series(s, c) = y[i++];
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void y_swathFromY(Tensor3& y_swath,
                  const Vector& y,
                  const Vector& y_f,
                  const Vector& sensor_response_f_grid,
                  const Index& npixel,
                  const Index& safe)
{
  // Sizes
  const Index ly = y.size();
  const Index nchannel = sensor_response_f_grid.size();
  const Index nswath = ly / (nchannel * npixel);
  
  ARTS_USER_ERROR_IF (nchannel * npixel * nswath != ly,
    "Length of *y* does not match given *npixel* and *sensor_response_f_grid*.")

  y_swath.resize(nswath, npixel, nchannel);

  Index i = 0;
  for (Index s=0; s<nswath; ++s) {
    for (Index p=0; p<npixel; ++p) {
      for (Index c=0; c<nchannel; ++c) {
        if (safe && (p > 0 || s > 0)) {
          ARTS_USER_ERROR_IF (fabs(y_f[i] - y_f[i-nchannel]) > 1,
                              "At least one channel varies in frequency.")
        }
        y_swath(s, p, c) = y[i++];
      }
    }
  }
}
