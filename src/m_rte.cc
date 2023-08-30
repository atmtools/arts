/**
  @file   m_rte.cc
  @author Patrick Eriksson <patrick.eriksson@chalmers.se>
  @date   2002-05-11

  @brief  Workspace methods for solving clear sky radiative transfer.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts.h"
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

  for (Index i = 0; i < iy_aux_vars.nelem(); i++) {
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
  const Index np = ppvar_src.nelem();

  ARTS_USER_ERROR_IF(np not_eq ppvar_dsrc.nelem(),
                     "ppvar_dsrc must have (np) elements")
  ARTS_USER_ERROR_IF(np not_eq ppvar_tramat.nelem(),
                     "ppvar_tramat must have (np) elements")
  ARTS_USER_ERROR_IF(np not_eq ppvar_cumtramat.nelem(),
                     "ppvar_cumtramat must have (np) elements")
  ARTS_USER_ERROR_IF(2 not_eq ppvar_dtramat.nelem() or
                         ppvar_dtramat.front().nelem() not_eq
                             ppvar_dtramat.back().nelem() or
                         ppvar_dtramat.front().nelem() not_eq np,
                     "ppvar_dtramat must (2 x np) elements")

  if (np == 0) {
    ppvar_rad.resize(0);
    ppvar_drad.resize(0);
    return;
  }

  const Index nq = ppvar_dsrc.front().nrows();
  const Index nf = ppvar_dsrc.front().ncols();

  const auto test_nf = [nf](auto &v) { return v.nelem() not_eq nf; };
  ARTS_USER_ERROR_IF(nf not_eq background_rad.nelem(),
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
  const Index np = ppvar_tramat.nelem();

  ARTS_USER_ERROR_IF(np not_eq ppvar_tramat.nelem(),
                     "ppvar_tramat must have (np) elements")
  ARTS_USER_ERROR_IF(np not_eq ppvar_cumtramat.nelem(),
                     "ppvar_cumtramat must have (np) elements")
  ARTS_USER_ERROR_IF(2 not_eq ppvar_dtramat.nelem() or
                         ppvar_dtramat.front().nelem() not_eq
                             ppvar_dtramat.back().nelem() or
                         ppvar_dtramat.front().nelem() not_eq np,
                     "ppvar_dtramat must (2 x np) elements")

  if (np == 0) {
    ppvar_rad.resize(0);
    ppvar_drad.resize(0);
    return;
  }

  const Index nq = ppvar_dtramat.front().front().nrows();
  const Index nf = ppvar_dtramat.front().front().ncols();

  const auto test_nf = [nf](auto &v) { return v.nelem() not_eq nf; };
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

  const Index np = ppvar_atm.nelem();
  if (np == 0) {
    ppvar_src.resize(0);
    ppvar_dsrc.resize(0);
    return;
  }

  const Index nf = ppvar_propmat.front().nelem();
  const Index nq = j_analytical_do ? jacobian_quantities.nelem() : 0;

  ppvar_src.resize(np, StokvecVector(nf));
  ppvar_dsrc.resize(np, StokvecMatrix(nq, nf));

  Vector B(nf);
  Matrix dB(nq, nf);

  const bool temperature_jacobian =
      j_analytical_do and do_temperature_jacobian(jacobian_quantities);

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

  const Index nf = ppvar_propmat.front().nelem();
  const Index nq = j_analytical_do ? jacobian_quantities.nelem() : 0;

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
                    const Ppath &ppath, const Index &iy_agenda_call1) {
  const Index np = ppath.np;

  // Init iy_aux and fill where possible
  const Index naux = iy_aux_vars.nelem();
  iy_aux.resize(naux);

  if (np == 0) {
    return;
  }

  const Index nf = background_transmittance.nelem();

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

  const Index np = ppvar_rad.nelem();
  if (np == 0) {
    ppvar_trans_cumulat.resize(np, 0, 0, 0);
    ppvar_trans_partial.resize(np, 0, 0, 0);
    ppvar_iy.resize(0, 0, np);
    iy.resize(0, 0);
    return;
  }

  const Index nf = ppvar_rad.front().nelem();

  // Copy back to ARTS external style
  diy_dpath = j_analytical_do ? get_standard_diy_dpath(jacobian_quantities, np,
                                                       nf, false)
                              : ArrayOfTensor3(0);
  ppvar_trans_cumulat.resize(np, nf, 4, 4);
  ppvar_trans_partial.resize(np, nf, 4, 4);
  ppvar_iy.resize(nf, 4, np);
  iy = to_matrix(ppvar_rad.front());
  for (Index ip = 0; ip < ppvar_rad.nelem(); ip++) {
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

  const Index nf = f_grid.nelem();
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

void RadiativePropertiesCalc(
    const Workspace &ws, ArrayOfPropmatVector &ppvar_propmat,
    ArrayOfPropmatMatrix &ppvar_dpropmat,
    ArrayOfStokvecVector &ppvar_src,
    ArrayOfStokvecMatrix &ppvar_dsrc,
    ArrayOfMuelmatVector &ppvar_tramat,
    ArrayOfArrayOfMuelmatMatrix &ppvar_dtramat,
    Vector &ppvar_distance, ArrayOfArrayOfVector &ppvar_ddistance,
    ArrayOfMuelmatVector &ppvar_cumtramat, const Ppath &ppath,
    const ArrayOfAtmPoint &ppvar_atm, const ArrayOfVector &ppvar_f,
    const Index &jacobian_do, const Agenda &ppvar_rtprop_agenda) {
  ppvar_rtprop_agendaExecute(
      ws, ppvar_propmat, ppvar_dpropmat, ppvar_src, ppvar_dsrc, ppvar_tramat,
      ppvar_dtramat, ppvar_distance, ppvar_ddistance, ppvar_cumtramat, ppath,
      ppvar_atm, ppvar_f, jacobian_do, ppvar_rtprop_agenda);
}

void RadiationBackgroundCalc(const Workspace &ws, StokvecVector &background_rad,
                             ArrayOfTensor3 &diy_dx, const Ppath &ppath,
                             const AtmField &atm_field, const Vector &f_grid,
                             const Tensor3 &iy_transmittance,
                             const MuelmatVector &background_transmittance,
                             const Index &jacobian_do,
                             const Agenda &rte_background_agenda) {
  rte_background_agendaExecute(
      ws, background_rad, diy_dx, ppath, atm_field, f_grid, iy_transmittance,
      background_transmittance, jacobian_do, rte_background_agenda);
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

  const Index nf = f_grid.nelem();

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
      iy_aux.resize(iy_aux1.nelem());
      for (Index q = 0; q < iy_aux1.nelem(); q++) {
        iy_aux[q].resize(nf, 4);
      }
      //
      diy_dx.resize(diy_dx1.nelem());
      for (Index q = 0; q < diy_dx1.nelem(); q++) {
        diy_dx[q].resize(diy_dx1[q].npages(), nf, 4);
      }
    }

    // Copy to output variables
    iy(i, joker) = iy1(0, joker);
    for (Index q = 0; q < iy_aux1.nelem(); q++) {
      iy_aux[q](i, joker) = iy_aux1[q](0, joker);
    }
    for (Index q = 0; q < diy_dx1.nelem(); q++) {
      diy_dx[q](joker, i, joker) = diy_dx1[q](joker, 0, joker);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iyMC(const Workspace& ws,
          Matrix& iy,
          ArrayOfMatrix& iy_aux,
          ArrayOfTensor3& diy_dx,
          const Index& iy_agenda_call1,
          const Tensor3& iy_transmittance,
          const Vector& rte_pos,
          const Vector& rte_los,
          const ArrayOfString& iy_aux_vars,
          const Index& jacobian_do,
          const AtmField& atm_field,
          const SurfaceField& surface_field,
          const Index& cloudbox_on,
          const ArrayOfIndex& cloudbox_limits,
          const Vector& f_grid,
          const ArrayOfArrayOfSingleScatteringData& scat_data,
          const Agenda& iy_space_agenda,
          const Agenda& surface_rtprop_agenda,
          const Agenda& propmat_clearsky_agenda,
          const Agenda& ppath_step_agenda,
          const Numeric& ppath_lmax,
          const Numeric& ppath_lraytrace,
          const Tensor4& pnd_field,
          const String& iy_unit,
          const Numeric& mc_std_err,
          const Index& mc_max_time,
          const Index& mc_max_iter,
          const Index& mc_min_iter,
          const Numeric& mc_taustep_limit,
          const Index& t_interp_order) {
  // Throw error if unsupported features are requested
  ARTS_USER_ERROR_IF (!cloudbox_on,
        "The cloudbox must be activated (cloudbox_on must be 1)");
  ARTS_USER_ERROR_IF (jacobian_do,
        "This method does not provide any jacobians (jacobian_do must be 0)");
  ARTS_USER_ERROR_IF (!iy_agenda_call1,
        "Recursive usage not possible (iy_agenda_call1 must be 1)");
  ARTS_USER_ERROR_IF (iy_transmittance.ncols(),
        "*iy_transmittance* must be empty");

  // Size output variables
  //
  const Index nf = f_grid.nelem();
  //
  iy.resize(nf, 4);
  diy_dx.resize(0);
  //
  //=== iy_aux part ===========================================================
  Index auxError = -1;
  {
    const Index naux = iy_aux_vars.nelem();
    iy_aux.resize(naux);
    //
    for (Index i = 0; i < naux; i++) {
      if (iy_aux_vars[i] == "Error (uncorrelated)") {
        auxError = i;
        iy_aux[i].resize(nf, 4);
      } else {
        ARTS_USER_ERROR (
          "In *iy_aux_vars* you have included: \"", iy_aux_vars[i],
          "\"\nThis choice is not recognised.")
      }
    }
  }
  //===========================================================================

  // Some MC variables are only local here
  //
  MCAntenna mc_antenna;
  mc_antenna.set_pencil_beam();

  // Pos and los must be matrices
  Matrix pos(1, 3), los(1, 2);
  //
  pos(0, joker) = rte_pos;
  los(0, joker) = rte_los;

  String fail_msg;
  bool failed = false;

  if (nf) {
#pragma omp parallel for if (!arts_omp_in_parallel() && nf > 1)
    for (Index f_index = 0; f_index < nf; f_index++) {
      if (failed) continue;

      try {
        // Seed reset for each loop. If not done, the errors
        // appear to be highly correlated.
        Index mc_seed;
        MCSetSeedFromTime(mc_seed);

        Vector y, mc_error;
        Index mc_iteration_count;
        Tensor3 mc_points;
        ArrayOfIndex mc_scat_order, mc_source_domain;

        MCGeneral(ws,
                  y,
                  mc_iteration_count,
                  mc_error,
                  mc_points,
                  mc_scat_order,
                  mc_source_domain,
                  mc_antenna,
                  f_grid,
                  f_index,
                  pos,
                  los,
                  ppath_step_agenda,
                  ppath_lmax,
                  ppath_lraytrace,
                  iy_space_agenda,
                  surface_rtprop_agenda,
                  propmat_clearsky_agenda,
                  surface_field,
                  atm_field,
                  cloudbox_on,
                  cloudbox_limits,
                  pnd_field,
                  scat_data,
                  1,
                  1,
                  1,
                  1,
                  iy_unit,
                  mc_seed,
                  mc_std_err,
                  mc_max_time,
                  mc_max_iter,
                  mc_min_iter,
                  mc_taustep_limit,
                  1,
                  t_interp_order);

        ARTS_ASSERT(y.nelem() == 4);

        iy(f_index, joker) = y;

        if (auxError >= 0) {
          iy_aux[auxError](f_index, joker) = mc_error;
        }
      } catch (const std::exception& e) {
        ostringstream os;
        os << "Error for f_index = " << f_index << " (" << f_grid[f_index]
           << ")" << endl
           << e.what();
#pragma omp critical(iyMC_fail)
        {
          failed = true;
          fail_msg = os.str();
        }
        continue;
      }
    }
  }

  ARTS_USER_ERROR_IF (failed, fail_msg);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iyReplaceFromAux(Matrix& iy,
                      const ArrayOfMatrix& iy_aux,
                      const ArrayOfString& iy_aux_vars,
                      const Index& jacobian_do,
                      const String& aux_var) {
  ARTS_USER_ERROR_IF (iy_aux.nelem() != iy_aux_vars.nelem(),
        "*iy_aux* and *iy_aux_vars* must have the same "
        "number of elements.");

  ARTS_USER_ERROR_IF (jacobian_do,
        "This method can not provide any jacobians and "
        "*jacobian_do* must be 0.");

  bool ready = false;

  for (Index i = 0; i < iy_aux.nelem() && !ready; i++) {
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
void yCalc(const Workspace& ws,
           Vector& y,
           Vector& y_f,
           ArrayOfIndex& y_pol,
           Matrix& y_pos,
           Matrix& y_los,
           ArrayOfVector& y_aux,
           Matrix& y_geo,
           Matrix& jacobian,
           const Index& atmgeom_checked,
           const Index& atmfields_checked,
           const AtmField& atm_field,
           const Index& cloudbox_on,
           const Index& cloudbox_checked,
           const Index& scat_data_checked,
           const Index& sensor_checked,
           const Vector& f_grid,
           const Matrix& sensor_pos,
           const Matrix& sensor_los,
           const Matrix& transmitter_pos,
           const Matrix& mblock_dlos,
           const Sparse& sensor_response,
           const Vector& sensor_response_f,
           const ArrayOfIndex& sensor_response_pol,
           const Matrix& sensor_response_dlos,
           const String& iy_unit,
           const Agenda& iy_main_agenda,
           const Agenda& jacobian_agenda,
           const Index& jacobian_do,
           const ArrayOfRetrievalQuantity& jacobian_quantities,
           const ArrayOfString& iy_aux_vars) {
  // Basics
  //
  //
  ARTS_USER_ERROR_IF (f_grid.empty(),
                      "The frequency grid is empty.");
  chk_if_increasing("f_grid", f_grid);
  ARTS_USER_ERROR_IF (f_grid[0] <= 0,
                      "All frequencies in *f_grid* must be > 0.");
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
    ARTS_USER_ERROR_IF (sensor_checked != 1,
        "The sensor variables must be flagged to have\n"
        "passed a consistency check (sensor_checked=1).");

  // Some sizes
  const Index nf = f_grid.nelem();
  const Index nlos = mblock_dlos.nrows();
  const Index n1y = sensor_response.nrows();
  const Index nmblock = sensor_pos.nrows();
  const Index niyb = nf * nlos * 4;

  //---------------------------------------------------------------------------
  // Allocations and resizing
  //---------------------------------------------------------------------------

  // Resize *y* and *y_XXX*
  //
  y.resize(nmblock * n1y);
  y_f.resize(nmblock * n1y);
  y_pol.resize(nmblock * n1y);
  y_pos.resize(nmblock * n1y, sensor_pos.ncols());
  y_los.resize(nmblock * n1y, sensor_los.ncols());
  y_geo.resize(nmblock * n1y, 5);
  y_geo = NAN;  // Will be replaced if relavant data are provided (*geo_pos*)

  // For y_aux we don't know the number of quantities, and we need to
  // store all output
  ArrayOfArrayOfVector iyb_aux_array(nmblock);

  // Jacobian variables
  //
  Index j_analytical_do = 0;
  ArrayOfArrayOfIndex jacobian_indices;
  //
  if (jacobian_do) {
    bool any_affine;
    jac_ranges_indices(jacobian_indices, any_affine, jacobian_quantities, true);
    jacobian.resize(nmblock * n1y,
                    jacobian_indices[jacobian_indices.nelem() - 1][1] + 1);
    jacobian = 0;
    //
    FOR_ANALYTICAL_JACOBIANS_DO2(j_analytical_do = 1;)
  } else {
    jacobian.resize(0, 0);
  }

  //---------------------------------------------------------------------------
  // The calculations
  //---------------------------------------------------------------------------

  String fail_msg;
  bool failed = false;

  if (nmblock >= arts_omp_get_max_threads() ||
      (nf <= nmblock && nmblock >= nlos)) {
#pragma omp parallel for
    for (Index mblock_index = 0; mblock_index < nmblock; mblock_index++) {
      // Skip remaining iterations if an error occurred
      if (failed) continue;

      yCalc_mblock_loop_body(failed,
                             fail_msg,
                             iyb_aux_array,
                             ws,
                             y,
                             y_f,
                             y_pol,
                             y_pos,
                             y_los,
                             y_geo,
                             jacobian,
                             atm_field,
                             cloudbox_on,
                             f_grid,
                             sensor_pos,
                             sensor_los,
                             transmitter_pos,
                             mblock_dlos,
                             sensor_response,
                             sensor_response_f,
                             sensor_response_pol,
                             sensor_response_dlos,
                             iy_unit,
                             iy_main_agenda,
                             jacobian_agenda,
                             jacobian_do,
                             jacobian_quantities,
                             jacobian_indices,
                             iy_aux_vars,
                             mblock_index,
                             n1y,
                             j_analytical_do);
    }  // End mblock loop
  } else {
    for (Index mblock_index = 0; mblock_index < nmblock; mblock_index++) {
      // Skip remaining iterations if an error occurred
      if (failed) continue;

      yCalc_mblock_loop_body(failed,
                             fail_msg,
                             iyb_aux_array,
                             ws,
                             y,
                             y_f,
                             y_pol,
                             y_pos,
                             y_los,
                             y_geo,
                             jacobian,
                             atm_field,
                             cloudbox_on,
                             f_grid,
                             sensor_pos,
                             sensor_los,
                             transmitter_pos,
                             mblock_dlos,
                             sensor_response,
                             sensor_response_f,
                             sensor_response_pol,
                             sensor_response_dlos,
                             iy_unit,
                             iy_main_agenda,
                             jacobian_agenda,
                             jacobian_do,
                             jacobian_quantities,
                             jacobian_indices,
                             iy_aux_vars,
                             mblock_index,
                             n1y,
                             j_analytical_do);
    }  // End mblock loop
  }

  // Rethrow exception if a runtime error occurred in the mblock loop
  ARTS_USER_ERROR_IF (failed, fail_msg);

  // Compile y_aux
  //
  const Index nq = iyb_aux_array[0].nelem();
  y_aux.resize(nq);
  //
  for (Index q = 0; q < nq; q++) {
    y_aux[q].resize(nmblock * n1y);
    //
    for (Index mblock_index = 0; mblock_index < nmblock; mblock_index++) {
      const Range rowind =
          get_rowindex_for_mblock(sensor_response, mblock_index);
      const Index row0 = rowind.offset;

      // The sensor response must be applied in a special way for
      // uncorrelated errors. Schematically: sqrt( H.^2 * y.^2 )
      if (iy_aux_vars[q] == "Error (uncorrelated)") {
        for (Index i = 0; i < n1y; i++) {
          const Index row = row0 + i;
          y_aux[q][row] = 0;
          for (Index j = 0; j < niyb; j++) {
            y_aux[q][row] +=
                pow(sensor_response(i, j) * iyb_aux_array[mblock_index][q][j],
                    (Numeric)2.0);
          }
          y_aux[q][row] = sqrt(y_aux[q][row]);
        }
      } else {
        mult(y_aux[q][rowind], sensor_response, iyb_aux_array[mblock_index][q]);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void yCalcAppend(const Workspace& ws,
                 Vector& y,
                 Vector& y_f,
                 ArrayOfIndex& y_pol,
                 Matrix& y_pos,
                 Matrix& y_los,
                 ArrayOfVector& y_aux,
                 Matrix& y_geo,
                 Matrix& jacobian,
                 ArrayOfRetrievalQuantity& jacobian_quantities,
                 const Index& atmfields_checked,
                 const Index& atmgeom_checked,
                 const AtmField& atm_field,
                 const Index& cloudbox_on,
                 const Index& cloudbox_checked,
                 const Index& scat_data_checked,
                 const Index& sensor_checked,
                 const Vector& f_grid,
                 const Matrix& sensor_pos,
                 const Matrix& sensor_los,
                 const Matrix& transmitter_pos,
                 const Matrix& mblock_dlos,
                 const Sparse& sensor_response,
                 const Vector& sensor_response_f,
                 const ArrayOfIndex& sensor_response_pol,
                 const Matrix& sensor_response_dlos,
                 const String& iy_unit,
                 const Agenda& iy_main_agenda,
                 const Agenda& jacobian_agenda,
                 const Index& jacobian_do,
                 const ArrayOfString& iy_aux_vars,
                 const ArrayOfRetrievalQuantity& jacobian_quantities_copy,
                 const Index& append_instrument_wfs) {
  // The jacobian indices of old and new part (without transformations)
  ArrayOfArrayOfIndex jacobian_indices, jacobian_indices_copy;
  {
    bool any_affine;
    jac_ranges_indices(
        jacobian_indices_copy, any_affine, jacobian_quantities_copy, true);
    jac_ranges_indices(jacobian_indices, any_affine, jacobian_quantities, true);
  }

  // Check consistency of data representing first measurement
  const Index n1 = y.nelem();
  Index nrq1 = 0;
  ARTS_USER_ERROR_IF (y.empty(),
                      "Input *y* is empty. Use *yCalc*");
  ARTS_USER_ERROR_IF (y_f.nelem() != n1,
                      "Lengths of input *y* and *y_f* are inconsistent.");
  ARTS_USER_ERROR_IF (y_pol.nelem() != n1,
                      "Lengths of input *y* and *y_pol* are inconsistent.");
  ARTS_USER_ERROR_IF (y_pos.nrows() != n1,
                      "Sizes of input *y* and *y_pos* are inconsistent.");
  ARTS_USER_ERROR_IF (y_los.nrows() != n1,
                      "Sizes of input *y* and *y_los* are inconsistent.");
  ARTS_USER_ERROR_IF (y_geo.nrows() != n1,
                      "Sizes of input *y* and *y_geo* are inconsistent.");
  if (jacobian_do) {
    nrq1 = jacobian_quantities_copy.nelem();
    ARTS_USER_ERROR_IF (jacobian.nrows() != n1,
                        "Sizes of *y* and *jacobian* are inconsistent.");
    ARTS_USER_ERROR_IF (jacobian.ncols() != jacobian_indices_copy[nrq1 - 1][1] + 1,
          "Size of input *jacobian* and size implied "
          "*jacobian_quantities_copy* are inconsistent.");
  }

  // Calculate new measurement
  //
  Vector y2, y_f2;
  Matrix y_pos2, y_los2, y_geo2, jacobian2;
  ArrayOfIndex y_pol2;
  ArrayOfVector y_aux2;
  //
  yCalc(ws,
        y2,
        y_f2,
        y_pol2,
        y_pos2,
        y_los2,
        y_aux2,
        y_geo2,
        jacobian2,
        atmfields_checked,
        atmgeom_checked,
        atm_field,
        cloudbox_on,
        cloudbox_checked,
        scat_data_checked,
        sensor_checked,
        f_grid,
        sensor_pos,
        sensor_los,
        transmitter_pos,
        mblock_dlos,
        sensor_response,
        sensor_response_f,
        sensor_response_pol,
        sensor_response_dlos,
        iy_unit,
        iy_main_agenda,
        jacobian_agenda,
        jacobian_do,
        jacobian_quantities,
        iy_aux_vars);

  // Consistency checks
  ARTS_USER_ERROR_IF (y_pos.ncols() != y_pos2.ncols(),
        "Different number of columns in *y_pos* between the measurements.");
  ARTS_USER_ERROR_IF (y_los.ncols() != y_los2.ncols(),
        "Different number of columns in *y_los* between the measurements.");

  // y and y_XXX
  //
  const Index n2 = y2.nelem();
  //
  {
    // Make copy of old measurement
    const Vector y1 = y, y_f1 = y_f;
    const Matrix y_pos1 = y_pos, y_los1 = y_los, y_geo1 = y_geo;
    const ArrayOfIndex y_pol1 = y_pol;
    const ArrayOfVector y_aux1 = y_aux;
    //
    y.resize(n1 + n2);
    y[Range(0, n1)] = y1;
    y[Range(n1, n2)] = y2;
    //
    y_f.resize(n1 + n2);
    y_f[Range(0, n1)] = y_f1;
    y_f[Range(n1, n2)] = y_f2;
    //
    y_pos.resize(n1 + n2, y_pos1.ncols());
    y_pos(Range(0, n1), joker) = y_pos1;
    y_pos(Range(n1, n2), joker) = y_pos2;
    //
    y_los.resize(n1 + n2, y_los1.ncols());
    y_los(Range(0, n1), joker) = y_los1;
    y_los(Range(n1, n2), joker) = y_los2;
    //
    y_geo.resize(n1 + n2, y_geo1.ncols());
    y_geo(Range(0, n1), joker) = y_geo1;
    y_geo(Range(n1, n2), joker) = y_geo2;
    //
    y_pol.resize(n1 + n2);
    for (Index i = 0; i < n1; i++) {
      y_pol[i] = y_pol1[i];
    }
    for (Index i = 0; i < n2; i++) {
      y_pol[n1 + i] = y_pol2[i];
    }

    // y_aux
    const Index na1 = y_aux1.nelem();
    const Index na2 = y_aux2.nelem();
    const Index na = max(na1, na2);
    //
    y_aux.resize(na);
    //
    for (Index a = 0; a < na; a++) {
      y_aux[a].resize(n1 + n2);
      if (a < na1) {
        y_aux[a][Range(0, n1)] = y_aux1[a];
      } else {
        y_aux[a][Range(0, n1)] = 0;
      }
      if (a < na2) {
        y_aux[a][Range(n1, n2)] = y_aux2[a];
      } else {
        y_aux[a][Range(n1, n2)] = 0;
      }
    }
  }

  // Jacobian and friends
  if (jacobian_do) {
    // Put in old jacobian_quantities and jacobian_indices as first part in
    // new version of these variables
    ArrayOfRetrievalQuantity jacobian_quantities2 = jacobian_quantities;
    ArrayOfArrayOfIndex jacobian_indices2 = jacobian_indices;
    //
    jacobian_quantities = jacobian_quantities_copy;
    jacobian_indices = jacobian_indices_copy;

    // Loop new jacobian_quantities to determine how new jacobian data shall
    // be inserted
    //
    const Index nrq2 = jacobian_quantities2.nelem();
    ArrayOfIndex map_table(nrq2);
    //
    for (Index q2 = 0; q2 < nrq2; q2++) {
      Index pos = -1;

      // Compare to old quantities, to determine if append shall be
      // considered. Some special checks performed here, grids checked later
      if (jacobian_quantities2[q2].Target().isSpeciesVMR() ||
          jacobian_quantities2[q2] == Jacobian::Atm::Temperature ||
          jacobian_quantities2[q2] == Jacobian::Special::ScatteringString ||
          jacobian_quantities2[q2].Target().isWind() ||
          jacobian_quantities2[q2] == Jacobian::Special::SurfaceString ||
          append_instrument_wfs) {
        for (Index q1 = 0; q1 < nrq1; q1++ && pos < 0) {  // FIXME: What is with this "&& pos < 0"
          if (jacobian_quantities2[q2].Target().sameTargetType(jacobian_quantities_copy[q1].Target())) {
            if (jacobian_quantities2[q2].Target().isSpeciesVMR()) {
              if (jacobian_quantities2[q2].Subtag() ==
                  jacobian_quantities_copy[q1].Subtag()) {
                if (jacobian_quantities2[q2].Mode() ==
                    jacobian_quantities_copy[q1].Mode()) {
                  pos = q1;
                } else {
                  ARTS_USER_ERROR (
                    "Jacobians for ", jacobian_quantities2[q2],
                    " shall be appended.\nThis requires "
                    "that the same retrieval unit is used "
                    "but it seems that this requirement is "
                    "not met.")
                }
              }
            } else if (jacobian_quantities2[q2] == Jacobian::Atm::Temperature) {
              if (jacobian_quantities2[q2].Subtag() ==
                  jacobian_quantities_copy[q1].Subtag()) {
                pos = q1;
              } else {
                ARTS_USER_ERROR (
                  "Jacobians for ", jacobian_quantities2[q2],
                  " shall be appended.\nThis requires "
                  "that HSE is either ON or OFF for both "
                  "parts but it seems that this requirement "
                  "is not met.")
              }
            } else if (jacobian_quantities[q2] == Jacobian::Special::ScatteringString) {
              if ((jacobian_quantities2[q2].Subtag() ==
                    jacobian_quantities_copy[q1].Subtag()) &&
                  (jacobian_quantities2[q2].SubSubtag() ==
                    jacobian_quantities_copy[q1].SubSubtag())) {
                pos = q1;
              }
            } else if (jacobian_quantities2[q2].Subtag() == jacobian_quantities_copy[q1].Subtag()) {
              pos = q1;
            }
          }
        }
      }

      // New quantity
      if (pos < 0) {
        map_table[q2] = jacobian_quantities.nelem();
        jacobian_quantities.push_back(jacobian_quantities2[q2]);
        ArrayOfIndex indices(2);
        indices[0] = jacobian_indices[jacobian_indices.nelem() - 1][1] + 1;
        indices[1] =
            indices[0] + jacobian_indices2[q2][1] - jacobian_indices2[q2][0];
        jacobian_indices.push_back(indices);
      }
      // Existing quantity
      else {
        map_table[q2] = pos;
        // Check if grids are equal
        ArrayOfVector grids1 = jacobian_quantities_copy[pos].Grids();
        ArrayOfVector grids2 = jacobian_quantities2[q2].Grids();
        bool any_wrong = false;
        if (grids1.nelem() != grids2.nelem()) {
          any_wrong = true;
        } else {
          for (Index g = 0; g < grids1.nelem(); g++) {
            if (grids1[g].nelem() != grids2[g].nelem()) {
              any_wrong = true;
            } else {
              for (Index e = 0; e < grids1[g].nelem(); e++) {
                const Numeric v1 = grids1[g][e];
                const Numeric v2 = grids2[g][e];
                if ((v1 == 0 && abs(v2) > 1e-9) || abs(v1 - v2) / v1 > 1e-6) {
                  any_wrong = true;
                }
              }
            }
          }
        }
        if (any_wrong) {
          ARTS_USER_ERROR (
            "Jacobians for ", jacobian_quantities2[q2],
            " shall be appended.\nThis requires that the "
            "same grids are used for both measurements,\nbut "
            "it seems that this requirement is not met.")
        }
      }
    }

    // Create and fill *jacobian*
    //
    const Index nrq = jacobian_quantities.nelem();
    const Matrix jacobian1 = jacobian;
    //
    jacobian.resize(n1 + n2, jacobian_indices[nrq - 1][1] + 1);
    jacobian = 0;
    //
    // Put in old part in top-left corner
    jacobian(Range(0, n1), Range(0, jacobian_indices_copy[nrq1 - 1][1] + 1)) =
        jacobian1;
    // New parts
    for (Index q2 = 0; q2 < nrq2; q2++) {
      jacobian(Range(n1, n2),
               Range(jacobian_indices[map_table[q2]][0],
                     jacobian_indices[map_table[q2]][1] -
                         jacobian_indices[map_table[q2]][0] + 1)) =
          jacobian2(
              joker,
              Range(jacobian_indices2[q2][0],
                    jacobian_indices2[q2][1] - jacobian_indices2[q2][0] + 1));
    }
  }
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
  const Index ny = y.nelem();
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
  const Index nchannel = sensor_response_f_grid.nelem();
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
  const Index nchannel = sensor_response_f_grid.nelem();
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
  const Index ly = y.nelem();
  const Index nchannel = sensor_response_f_grid.nelem();
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
  const Index ly = y.nelem();
  const Index nchannel = sensor_response_f_grid.nelem();
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
