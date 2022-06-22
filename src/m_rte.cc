/* Copyright (C) 2002-2012
   Patrick Eriksson <patrick.eriksson@chalmers.se>
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
#include "arts_omp.h"
#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "jacobian.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "montecarlo.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"
#include "transmissionmatrix.h"
#include <cmath>
#include <stdexcept>

extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT;
extern const Index GFIELD4_FIELD_NAMES;
extern const Index GFIELD4_P_GRID;
extern const Index GFIELD4_LAT_GRID;
extern const Index GFIELD4_LON_GRID;

/*===========================================================================
  === Workspace methods 
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void iyApplyUnit(Matrix& iy,
                 ArrayOfMatrix& iy_aux,
                 const Index& stokes_dim,
                 const Vector& f_grid,
                 const ArrayOfString& iy_aux_vars,
                 const String& iy_unit,
                 const Verbosity&) {
  ARTS_USER_ERROR_IF (iy_unit == "1",
    "No need to use this method with *iy_unit* = \"1\".");

  ARTS_USER_ERROR_IF (max(iy(joker, 0)) > 1e-3,
      "The spectrum matrix *iy* is required to have original radiance\n"
      "unit, but this seems not to be the case. This as a value above\n"
      "1e-3 is found in *iy*.")

  // Polarisation index variable
  ArrayOfIndex i_pol(stokes_dim);
  for (Index is = 0; is < stokes_dim; is++) {
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
void iyCalc(Workspace& ws,
            Matrix& iy,
            ArrayOfMatrix& iy_aux,
            Ppath& ppath,
            const Index& atmfields_checked,
            const Index& atmgeom_checked,
            const ArrayOfString& iy_aux_vars,
            const Index& iy_id,
            const Index& cloudbox_on,
            const Index& cloudbox_checked,
            const Index& scat_data_checked,
            const Vector& f_grid,            
            const EnergyLevelMap& nlte_field,
            const Vector& rte_pos,
            const Vector& rte_los,
            const Vector& rte_pos2,
            const String& iy_unit,
            const Agenda& iy_main_agenda,
            const Verbosity&) {
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
                        1,
                        iy_transmittance,
                        iy_aux_vars,
                        iy_id,
                        iy_unit,
                        cloudbox_on,
                        0,
                        f_grid,
                        nlte_field,
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

/* Workspace method: Doxygen documentation will be auto-generated */
void iyEmissionHybrid(Workspace& ws,
              Matrix& iy,
              ArrayOfMatrix& iy_aux,
              ArrayOfTensor3& diy_dx,
              Vector& ppvar_p,
              Vector& ppvar_t,
              EnergyLevelMap& ppvar_nlte,
              Matrix& ppvar_vmr,
              Matrix& ppvar_wind,
              Matrix& ppvar_mag,
              Matrix& ppvar_pnd,
              Matrix& ppvar_f,
              Tensor3& ppvar_iy,
              Tensor4& ppvar_trans_cumulat,
              Tensor4& ppvar_trans_partial,
              const Index& iy_id,
              const Index& stokes_dim,
              const Vector& f_grid,
              const Index& atmosphere_dim,
              const Vector& p_grid,
              const Tensor3& t_field,
              const EnergyLevelMap& nlte_field,
              const Tensor4& vmr_field,
              const ArrayOfArrayOfSpeciesTag& abs_species,
              const Tensor3& wind_u_field,
              const Tensor3& wind_v_field,
              const Tensor3& wind_w_field,
              const Tensor3& mag_u_field,
              const Tensor3& mag_v_field,
              const Tensor3& mag_w_field,
              const Index& cloudbox_on,
              const ArrayOfIndex& cloudbox_limits,
              const Tensor4& pnd_field,
              const ArrayOfTensor4& dpnd_field_dx,
              const ArrayOfString& scat_species,
              const ArrayOfArrayOfSingleScatteringData& scat_data,
              const String& iy_unit,
              const ArrayOfString& iy_aux_vars,
              const Index& jacobian_do,
              const ArrayOfRetrievalQuantity& jacobian_quantities,
              const Agenda& propmat_clearsky_agenda,
              const Agenda& water_p_eq_agenda,
              const String& rt_integration_option,               
              const Agenda& iy_main_agenda,
              const Agenda& iy_space_agenda,
              const Agenda& iy_surface_agenda,
              const Agenda& iy_cloudbox_agenda,
              const Index& iy_agenda_call1,
              const Tensor3& iy_transmittance,
              const Ppath& ppath,
              const Vector& rte_pos2,
              const Numeric& rte_alonglos_v,
              const Tensor3& surface_props_data,
              const Tensor7& cloudbox_field,
              const Vector& za_grid,
              const Index& Naa,
              const Index& t_interp_order,
              const Verbosity& verbosity) {
  // If cloudbox off, switch to use clearsky method
  if (!cloudbox_on) {
    Tensor4 dummy;
    iyEmissionStandard(ws,
                       iy,
                       iy_aux,
                       diy_dx,
                       ppvar_p,
                       ppvar_t,
                       ppvar_nlte,
                       ppvar_vmr,
                       ppvar_wind,
                       ppvar_mag,
                       ppvar_f,
                       ppvar_iy,
                       ppvar_trans_cumulat,
                       ppvar_trans_partial,
                       iy_id,
                       stokes_dim,
                       f_grid,
                       atmosphere_dim,
                       p_grid,
                       t_field,
                       nlte_field,
                       vmr_field,
                       abs_species,
                       wind_u_field,
                       wind_v_field,
                       wind_w_field,
                       mag_u_field,
                       mag_v_field,
                       mag_w_field,
                       cloudbox_on,
                       iy_unit,
                       iy_aux_vars,
                       jacobian_do,
                       jacobian_quantities,
                       ppath,
                       rte_pos2,
                       propmat_clearsky_agenda,
                       water_p_eq_agenda,
                       rt_integration_option,
                       iy_main_agenda,
                       iy_space_agenda,
                       iy_surface_agenda,
                       iy_cloudbox_agenda,
                       iy_agenda_call1,
                       iy_transmittance,
                       rte_alonglos_v,
                       surface_props_data,
                       verbosity);
    return;
  }
  //  Init Jacobian quantities?
  const Index j_analytical_do = jacobian_do ? do_analytical_jacobian<1>(jacobian_quantities) : 0;

  // Some basic sizes
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index nq = j_analytical_do ? jacobian_quantities.nelem() : 0;

  // Radiative background index
  const Index rbi = ppath_what_background(ppath);

  // Throw error if unsupported features are requested
  if (atmosphere_dim != 1)
    throw runtime_error(
        "With cloudbox on, this method handles only 1D calculations.");
  if (Naa < 3) throw runtime_error("Naa must be > 2.");
  if (jacobian_do)
    if (dpnd_field_dx.nelem() != jacobian_quantities.nelem())
      throw runtime_error(
          "*dpnd_field_dx* not properly initialized:\n"
          "Number of elements in dpnd_field_dx must be equal number of jacobian"
          " quantities.\n(Note: jacobians have to be defined BEFORE *pnd_field*"
          " is calculated/set.");
  if (rbi < 1 || rbi > 9)
    throw runtime_error(
        "ppath.background is invalid. Check your "
        "calculation of *ppath*?");
  if (rbi == 3 || rbi == 4)
    throw runtime_error(
        "The propagation path ends inside or at boundary of "
        "the cloudbox.\nFor this method, *ppath* must be "
        "calculated in this way:\n   ppathCalc( cloudbox_on = 0 ).");
  // iy_aux_vars checked below
  // Checks of i_field
  if (cloudbox_field.ncols() != stokes_dim)
    throw runtime_error(
        "Obtained *cloudbox_field* number of Stokes elements inconsistent with "
        "*stokes_dim*.");
  if (cloudbox_field.nrows() != 1)
    throw runtime_error(
        "Obtained *cloudbox_field* has wrong number of azimuth angles.");
  if (cloudbox_field.npages() != za_grid.nelem())
    throw runtime_error(
        "Obtained *cloudbox_field* number of zenith angles inconsistent with "
        "*za_grid*.");
  if (cloudbox_field.nbooks() != 1)
    throw runtime_error(
        "Obtained *cloudbox_field* has wrong number of longitude points.");
  if (cloudbox_field.nshelves() != 1)
    throw runtime_error(
        "Obtained *cloudbox_field* has wrong number of latitude points.");
  if (cloudbox_field.nvitrines() != cloudbox_limits[1] - cloudbox_limits[0] + 1)
    throw runtime_error(
        "Obtained *cloudbox_field* number of pressure points inconsistent with "
        "*cloudbox_limits*.");
  if (cloudbox_field.nlibraries() != nf)
    throw runtime_error(
        "Obtained *cloudbox_field* number of frequency points inconsistent with "
        "*f_grid*.");
  
  // Set diy_dpath if we are doing are doing jacobian calculations
  ArrayOfTensor3 diy_dpath = j_analytical_do ? get_standard_diy_dpath(jacobian_quantities, np, nf, ns, false) : ArrayOfTensor3(0);
  
  // Set the species pointers if we are doing jacobian
  const ArrayOfIndex jac_species_i = j_analytical_do ? get_pointers_for_analytical_species(jacobian_quantities, abs_species) : ArrayOfIndex(0);
  
  // Start diy_dx out if we are doing the first run and are doing jacobian calculations
  if (j_analytical_do and iy_agenda_call1) diy_dx = get_standard_starting_diy_dx(jacobian_quantities, np, nf, ns, false);
  
  // Checks that the scattering species are treated correctly if their derivatives are needed (we can here discard the Array)
  if (j_analytical_do and iy_agenda_call1) get_pointers_for_scat_species(jacobian_quantities, scat_species, cloudbox_on);

  // Init iy_aux and fill where possible
  const Index naux = iy_aux_vars.nelem();
  iy_aux.resize(naux);
  //
  for (Index i = 0; i < naux; i++) {
    iy_aux[i].resize(nf, ns);

    if (iy_aux_vars[i] == "Optical depth") { /*pass*/
    }                                        // Filled below
    else if (iy_aux_vars[i] == "Radiative background")
      iy_aux[i] = (Numeric)min((Index)2, rbi - 1);
    else {
      ostringstream os;
      os << "The only allowed strings in *iy_aux_vars* are:\n"
         << "  \"Radiative background\"\n"
         << "  \"Optical depth\"\n"
         << "but you have selected: \"" << iy_aux_vars[i] << "\"";
      throw runtime_error(os.str());
    }
  }

  // Get atmospheric and radiative variables along the propagation path
  ppvar_trans_cumulat.resize(np, nf, ns, ns);
  ppvar_trans_partial.resize(np, nf, ns, ns);
  ppvar_iy.resize(nf, ns, np);

  ArrayOfTransmissionMatrix lyr_tra(np, TransmissionMatrix(nf, ns));
  ArrayOfRadiationVector lvl_rad(np, RadiationVector(nf, ns));
  ArrayOfArrayOfRadiationVector dlvl_rad(
      np, ArrayOfRadiationVector(nq, RadiationVector(nf, ns)));
  ArrayOfRadiationVector src_rad(np, RadiationVector(nf, ns));
  ArrayOfArrayOfRadiationVector dsrc_rad(
      np, ArrayOfRadiationVector(nq, RadiationVector(nf, ns)));

  ArrayOfArrayOfTransmissionMatrix dlyr_tra_above(
      np, ArrayOfTransmissionMatrix(nq, TransmissionMatrix(nf, ns)));
  ArrayOfArrayOfTransmissionMatrix dlyr_tra_below(
      np, ArrayOfTransmissionMatrix(nq, TransmissionMatrix(nf, ns)));

  ArrayOfIndex clear2cloudy;
  //
  if (np == 1 && rbi == 1) {  // i.e. ppath is totally outside the atmosphere:
    ppvar_p.resize(0);
    ppvar_t.resize(0);
    ppvar_vmr.resize(0, 0);
    ppvar_wind.resize(0, 0);
    ppvar_mag.resize(0, 0);
    ppvar_f.resize(0, 0);
    ppvar_trans_cumulat = 0;
    ppvar_trans_partial = 0;
    for (Index iv = 0; iv < nf; iv++) {
      for (Index is = 0; is < ns; is++) {
        ppvar_trans_cumulat(0,iv,is,is) = 1;
        ppvar_trans_partial(0,iv,is,is) = 1;
      }
    }
  } else {
    // Basic atmospheric variables
    get_ppath_atmvars(ppvar_p,
                      ppvar_t,
                      ppvar_nlte,
                      ppvar_vmr,
                      ppvar_wind,
                      ppvar_mag,
                      ppath,
                      atmosphere_dim,
                      p_grid,
                      t_field,
                      nlte_field,
                      vmr_field,
                      wind_u_field,
                      wind_v_field,
                      wind_w_field,
                      mag_u_field,
                      mag_v_field,
                      mag_w_field);

    get_ppath_f(
        ppvar_f, ppath, f_grid, atmosphere_dim, rte_alonglos_v, ppvar_wind);

    // here, the cloudbox is on, ie we don't need to check and branch this here
    // anymore.
    ArrayOfMatrix ppvar_dpnd_dx;
    //
    get_ppath_cloudvars(clear2cloudy,
                        ppvar_pnd,
                        ppvar_dpnd_dx,
                        ppath,
                        atmosphere_dim,
                        cloudbox_limits,
                        pnd_field,
                        dpnd_field_dx);

    // Size radiative variables always used
    Vector B(nf);
    PropagationMatrix K_this(nf, ns), K_past(nf, ns), Kp(nf, ns);
    StokesVector a(nf, ns), S(nf, ns), Sp(nf, ns);
    ArrayOfIndex lte(np);

    // Init variables only used if analytical jacobians done
    Vector dB_dT(0);
    ArrayOfPropagationMatrix dK_this_dx(nq), dK_past_dx(nq), dKp_dx(nq);
    ArrayOfStokesVector da_dx(nq), dS_dx(nq), dSp_dx(nq);

    // HSE variables
    Index temperature_derivative_position = -1;
    bool do_hse = false;

    if (j_analytical_do) {
      dB_dT.resize(nf);
      FOR_ANALYTICAL_JACOBIANS_DO(
          dK_this_dx[iq] = PropagationMatrix(nf, ns);
          dK_past_dx[iq] = PropagationMatrix(nf, ns);
          dKp_dx[iq] = PropagationMatrix(nf, ns);
          da_dx[iq] = StokesVector(nf, ns);
          dS_dx[iq] = StokesVector(nf, ns);
          dSp_dx[iq] = StokesVector(nf, ns);
          if (jacobian_quantities[iq] == Jacobian::Atm::Temperature) {
            temperature_derivative_position = iq;
            do_hse = jacobian_quantities[iq].Subtag() == "HSE on";
          })
    }
    const bool temperature_jacobian =
        j_analytical_do and do_temperature_jacobian(jacobian_quantities);

    // Loop ppath points and determine radiative properties
    for (Index ip = 0; ip < np; ip++) {
      get_stepwise_blackbody_radiation(
          B, dB_dT, ppvar_f(joker, ip), ppvar_t[ip], temperature_jacobian);

      get_stepwise_clearsky_propmat(ws,
                                    K_this,
                                    S,
                                    lte[ip],
                                    dK_this_dx,
                                    dS_dx,
                                    propmat_clearsky_agenda,
                                    jacobian_quantities,
                                    ppvar_f(joker, ip),
                                    ppvar_mag(joker, ip),
                                    ppath.los(ip, joker),
                                    ppvar_nlte[ip],
                                    ppvar_vmr(joker, ip),
                                    ppvar_t[ip],
                                    ppvar_p[ip],
                                    j_analytical_do);

      if (j_analytical_do)
        adapt_stepwise_partial_derivatives(dK_this_dx,
                                           dS_dx,
                                           jacobian_quantities,
                                           ppvar_f(joker, ip),
                                           ppath.los(ip, joker),
                                           lte[ip],
                                           atmosphere_dim,
                                           j_analytical_do);

      if (clear2cloudy[ip] + 1) {
        get_stepwise_scattersky_propmat(a,
                                        Kp,
                                        da_dx,
                                        dKp_dx,
                                        jacobian_quantities,
                                        ppvar_pnd(joker, Range(ip, 1)),
                                        ppvar_dpnd_dx,
                                        ip,
                                        scat_data,
                                        ppath.los(ip, joker),
                                        ppvar_t[Range(ip, 1)],
                                        atmosphere_dim,
                                        jacobian_do);
        a += K_this;
        K_this += Kp;

        if (j_analytical_do)
          FOR_ANALYTICAL_JACOBIANS_DO(da_dx[iq] += dK_this_dx[iq];
                                      dK_this_dx[iq] += dKp_dx[iq];)

        Vector aa_grid;
        nlinspace(aa_grid, 0, 360, Naa);
        //
        get_stepwise_scattersky_source(Sp,
                                       dSp_dx,
                                       jacobian_quantities,
                                       ppvar_pnd(joker, ip),
                                       ppvar_dpnd_dx,
                                       ip,
                                       scat_data,
                                       cloudbox_field,
                                       za_grid,
                                       aa_grid,
                                       ppath.los(Range(ip, 1), joker),
                                       ppath.gp_p[ip],
                                       ppvar_t[Range(ip, 1)],
                                       atmosphere_dim,
                                       jacobian_do,
                                       t_interp_order);
        S += Sp;

        if (j_analytical_do)
          FOR_ANALYTICAL_JACOBIANS_DO(dS_dx[iq] += dSp_dx[iq];)
      } else {  // no particles present at this level
        a = K_this;
        if (j_analytical_do)
          FOR_ANALYTICAL_JACOBIANS_DO(da_dx[iq] = dK_this_dx[iq];)
      }

      if (ip not_eq 0) {
        const Numeric dr_dT_past =
            do_hse ? ppath.lstep[ip - 1] / (2.0 * ppvar_t[ip - 1]) : 0;
        const Numeric dr_dT_this =
            do_hse ? ppath.lstep[ip - 1] / (2.0 * ppvar_t[ip]) : 0;
        stepwise_transmission(lyr_tra[ip],
                              dlyr_tra_above[ip],
                              dlyr_tra_below[ip],
                              K_past,
                              K_this,
                              dK_past_dx,
                              dK_this_dx,
                              ppath.lstep[ip - 1],
                              dr_dT_past,
                              dr_dT_this,
                              temperature_derivative_position);
      }

      stepwise_source(src_rad[ip],
                      dsrc_rad[ip],
                      K_this,
                      a,
                      S,
                      dK_this_dx,
                      da_dx,
                      dS_dx,
                      B,
                      dB_dT,
                      jacobian_quantities,
                      jacobian_do);

      swap(K_past, K_this);
      swap(dK_past_dx, dK_this_dx);
    }
  }

  const ArrayOfTransmissionMatrix tot_tra =
      cumulative_transmission(lyr_tra, CumulativeTransmission::Forward);

  // iy_transmittance
  Tensor3 iy_trans_new;
  if (iy_agenda_call1)
    iy_trans_new = tot_tra[np - 1];
  else
    iy_transmittance_mult(iy_trans_new, iy_transmittance, tot_tra[np - 1]);

  // Copy transmission to iy_aux
  for (Index i = 0; i < naux; i++)
    if (iy_aux_vars[i] == "Optical depth")
      for (Index iv = 0; iv < nf; iv++)
        iy_aux[i](iv, joker) = -log(ppvar_trans_cumulat(np - 1, iv, 0, 0));

  // Radiative background
  get_iy_of_background(ws,
                       iy,
                       diy_dx,
                       iy_trans_new,
                       iy_id,
                       jacobian_do,
                       jacobian_quantities,
                       ppath,
                       rte_pos2,
                       atmosphere_dim,
                       nlte_field,
                       cloudbox_on,
                       stokes_dim,
                       f_grid,
                       iy_unit,
                       surface_props_data,
                       iy_main_agenda,
                       iy_space_agenda,
                       iy_surface_agenda,
                       iy_cloudbox_agenda,
                       iy_agenda_call1,
                       verbosity);

  lvl_rad[np - 1] = iy;

  // Radiative transfer calculations
  for (Index ip = np - 2; ip >= 0; ip--) {
    lvl_rad[ip] = lvl_rad[ip + 1];
    update_radiation_vector(lvl_rad[ip],
                            dlvl_rad[ip],
                            dlvl_rad[ip + 1],
                            src_rad[ip],
                            src_rad[ip + 1],
                            dsrc_rad[ip],
                            dsrc_rad[ip + 1],
                            lyr_tra[ip + 1],
                            tot_tra[ip],
                            dlyr_tra_above[ip + 1],
                            dlyr_tra_below[ip + 1],
                            PropagationMatrix(),
                            PropagationMatrix(),
                            ArrayOfPropagationMatrix(),
                            ArrayOfPropagationMatrix(),
                            Numeric(),
                            Vector(),
                            Vector(),
                            0,
                            0,
                            RadiativeTransferSolver::Emission);
  }

  // Copy back to ARTS external style
  iy = lvl_rad[0];
  for (Index ip = 0; ip < lvl_rad.nelem(); ip++) {
    ppvar_trans_cumulat(ip, joker, joker, joker) = tot_tra[ip];
    ppvar_trans_partial(ip, joker, joker, joker) = lyr_tra[ip];
    ppvar_iy(joker, joker, ip) = lvl_rad[ip];
    if (j_analytical_do)
      FOR_ANALYTICAL_JACOBIANS_DO(diy_dpath[iq](ip, joker, joker) =
                                      dlvl_rad[ip][iq];);
  }

  // Finalize analytical Jacobians
  if (j_analytical_do)
    rtmethods_jacobian_finalisation(ws,
                                    diy_dx,
                                    diy_dpath,
                                    ns,
                                    nf,
                                    np,
                                    atmosphere_dim,
                                    ppath,
                                    ppvar_p,
                                    ppvar_t,
                                    ppvar_vmr,
                                    iy_agenda_call1,
                                    iy_transmittance,
                                    water_p_eq_agenda,
                                    jacobian_quantities,
                                    jac_species_i);

  // Unit conversions
  if (iy_agenda_call1)
    rtmethods_unit_conversion(iy,
                              diy_dx,
                              ppvar_iy,
                              ns,
                              np,
                              f_grid,
                              ppath,
                              jacobian_quantities,
                              j_analytical_do,
                              iy_unit);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iyEmissionStandard(
    Workspace& ws,
    Matrix& iy,
    ArrayOfMatrix& iy_aux,
    ArrayOfTensor3& diy_dx,
    Vector& ppvar_p,
    Vector& ppvar_t,
    EnergyLevelMap& ppvar_nlte,
    Matrix& ppvar_vmr,
    Matrix& ppvar_wind,
    Matrix& ppvar_mag,
    Matrix& ppvar_f,
    Tensor3& ppvar_iy,
    Tensor4& ppvar_trans_cumulat,
    Tensor4& ppvar_trans_partial,
    const Index& iy_id,
    const Index& stokes_dim,
    const Vector& f_grid,
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Tensor3& t_field,
    const EnergyLevelMap& nlte_field,
    const Tensor4& vmr_field,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Tensor3& wind_u_field,
    const Tensor3& wind_v_field,
    const Tensor3& wind_w_field,
    const Tensor3& mag_u_field,
    const Tensor3& mag_v_field,
    const Tensor3& mag_w_field,
    const Index& cloudbox_on,
    const String& iy_unit,
    const ArrayOfString& iy_aux_vars,
    const Index& jacobian_do,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Ppath& ppath,
    const Vector& rte_pos2,
    const Agenda& propmat_clearsky_agenda,
    const Agenda& water_p_eq_agenda,
    const String& rt_integration_option,
    const Agenda& iy_main_agenda,
    const Agenda& iy_space_agenda,
    const Agenda& iy_surface_agenda,
    const Agenda& iy_cloudbox_agenda,
    const Index& iy_agenda_call1,
    const Tensor3& iy_transmittance,
    const Numeric& rte_alonglos_v,
    const Tensor3& surface_props_data,
    const Verbosity& verbosity) {
  //  Init Jacobian quantities?
  const Index j_analytical_do = jacobian_do ? do_analytical_jacobian<2>(jacobian_quantities) : 0;
  
  // Some basic sizes
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index nq = j_analytical_do ? jacobian_quantities.nelem() : 0;

  // Radiative background index
  const Index rbi = ppath_what_background(ppath);

  // Checks of input
  ARTS_USER_ERROR_IF (rbi < 1 || rbi > 9,
        "ppath.background is invalid. Check your "
        "calculation of *ppath*?");
  ARTS_USER_ERROR_IF (!iy_agenda_call1 && np == 1 && rbi == 2,
        "A secondary propagation path starting at the "
        "surface and is going directly into the surface "
        "is found. This is not allowed.");

  // Set diy_dpath if we are doing are doing jacobian calculations
  ArrayOfTensor3 diy_dpath = j_analytical_do ? get_standard_diy_dpath(jacobian_quantities, np, nf, ns, false) : ArrayOfTensor3(0);
  
  // Set the species pointers if we are doing jacobian
  const ArrayOfIndex jac_species_i = j_analytical_do ? get_pointers_for_analytical_species(jacobian_quantities, abs_species) : ArrayOfIndex(0);
  
  // Start diy_dx out if we are doing the first run and are doing jacobian calculations
  if (j_analytical_do and iy_agenda_call1) diy_dx = get_standard_starting_diy_dx(jacobian_quantities, np, nf, ns, false);

  // Init iy_aux and fill where possible
  const Index naux = iy_aux_vars.nelem();
  iy_aux.resize(naux);
  //
  Index auxOptDepth = -1;
  //
  for (Index i = 0; i < naux; i++) {
    iy_aux[i].resize(nf, ns);
    iy_aux[i] = 0;

    if (iy_aux_vars[i] == "Radiative background")
      iy_aux[i](joker, 0) = (Numeric)min((Index)2, rbi - 1);
    else if (iy_aux_vars[i] == "Optical depth")
      auxOptDepth = i;
    else {
      ARTS_USER_ERROR (
        "The only allowed strings in *iy_aux_vars* are:\n"
        "  \"Radiative background\"\n"
        "  \"Optical depth\"\n"
        "but you have selected: \"", iy_aux_vars[i], "\"")
    }
  }

  // Get atmospheric and radiative variables along the propagation path
  ppvar_trans_cumulat.resize(np, nf, ns, ns);
  ppvar_trans_partial.resize(np, nf, ns, ns);
  ppvar_iy.resize(nf, ns, np);

  ArrayOfRadiationVector lvl_rad(np, RadiationVector(nf, ns));
  ArrayOfArrayOfRadiationVector dlvl_rad(
      np, ArrayOfRadiationVector(nq, RadiationVector(nf, ns)));

  ArrayOfRadiationVector src_rad(np, RadiationVector(nf, ns));
  ArrayOfArrayOfRadiationVector dsrc_rad(
      np, ArrayOfRadiationVector(nq, RadiationVector(nf, ns)));

  ArrayOfTransmissionMatrix lyr_tra(np, TransmissionMatrix(nf, ns));
  ArrayOfArrayOfTransmissionMatrix dlyr_tra_above(
      np, ArrayOfTransmissionMatrix(nq, TransmissionMatrix(nf, ns)));
  ArrayOfArrayOfTransmissionMatrix dlyr_tra_below(
      np, ArrayOfTransmissionMatrix(nq, TransmissionMatrix(nf, ns)));
  
  ArrayOfPropagationMatrix K(np, PropagationMatrix(nf, ns));
  ArrayOfArrayOfPropagationMatrix dK_dx(np);
  Vector r(np);
  ArrayOfVector dr_below(np, Vector(nq, 0));
  ArrayOfVector dr_above(np, Vector(nq, 0));

  if (np == 1 && rbi == 1) {  // i.e. ppath is totally outside the atmosphere:
    ppvar_p.resize(0);
    ppvar_t.resize(0);
    ppvar_vmr.resize(0, 0);
    ppvar_wind.resize(0, 0);
    ppvar_mag.resize(0, 0);
    ppvar_f.resize(0, 0);
    ppvar_trans_cumulat = 0;
    ppvar_trans_partial = 0;
    for (Index iv = 0; iv < nf; iv++) {
      for (Index is = 0; is < ns; is++) {
        ppvar_trans_cumulat(0,iv,is,is) = 1;
        ppvar_trans_partial(0,iv,is,is) = 1;
      }
    }
    
  } else {
    // Basic atmospheric variables
    get_ppath_atmvars(ppvar_p,
                      ppvar_t,
                      ppvar_nlte,
                      ppvar_vmr,
                      ppvar_wind,
                      ppvar_mag,
                      ppath,
                      atmosphere_dim,
                      p_grid,
                      t_field,
                      nlte_field,
                      vmr_field,
                      wind_u_field,
                      wind_v_field,
                      wind_w_field,
                      mag_u_field,
                      mag_v_field,
                      mag_w_field);

    get_ppath_f(
        ppvar_f, ppath, f_grid, atmosphere_dim, rte_alonglos_v, ppvar_wind);
    
    const bool temperature_jacobian =
        j_analytical_do and do_temperature_jacobian(jacobian_quantities);

    // Size radiative variables always used
    Vector B(nf);
    StokesVector a(nf, ns), S(nf, ns);

    // Init variables only used if analytical jacobians done
    Vector dB_dT(temperature_jacobian ? nf : 0);
    ArrayOfStokesVector da_dx(nq), dS_dx(nq);

    // HSE variables
    Index temperature_derivative_position = -1;
    bool do_hse = false;

    if (j_analytical_do) {
      for (Index ip = 0; ip < np; ip++) {
        dK_dx[ip].resize(nq);
        FOR_ANALYTICAL_JACOBIANS_DO(dK_dx[ip][iq] = PropagationMatrix(nf, ns);)
      }
      FOR_ANALYTICAL_JACOBIANS_DO(
          da_dx[iq] = StokesVector(nf, ns); dS_dx[iq] = StokesVector(nf, ns);
          if (jacobian_quantities[iq] == Jacobian::Atm::Temperature) {
            temperature_derivative_position = iq;
            do_hse = jacobian_quantities[iq].Subtag() == "HSE on";
          })
    }

    ArrayOfString fail_msg;
    bool do_abort = false;

    // Loop ppath points and determine radiative properties
#pragma omp parallel for if (!arts_omp_in_parallel()) \
    firstprivate(ws, a, B, dB_dT, S, da_dx, dS_dx)
    for (Index ip = 0; ip < np; ip++) {
      if (do_abort) continue;
      try {
        get_stepwise_blackbody_radiation(
            B, dB_dT, ppvar_f(joker, ip), ppvar_t[ip], temperature_jacobian);

        Index lte;
        get_stepwise_clearsky_propmat(ws,
                                      K[ip],
                                      S,
                                      lte,
                                      dK_dx[ip],
                                      dS_dx,
                                      propmat_clearsky_agenda,
                                      jacobian_quantities,
                                      ppvar_f(joker, ip),
                                      ppvar_mag(joker, ip),
                                      ppath.los(ip, joker),
                                      ppvar_nlte[ip],
                                      ppvar_vmr(joker, ip),
                                      ppvar_t[ip],
                                      ppvar_p[ip],
                                      j_analytical_do);

        if (j_analytical_do)
          adapt_stepwise_partial_derivatives(dK_dx[ip],
                                             dS_dx,
                                             jacobian_quantities,
                                             ppvar_f(joker, ip),
                                             ppath.los(ip, joker),
                                             lte,
                                             atmosphere_dim,
                                             j_analytical_do);

        // Here absorption equals extinction
        a = K[ip];
        if (j_analytical_do)
          FOR_ANALYTICAL_JACOBIANS_DO(da_dx[iq] = dK_dx[ip][iq];);

        stepwise_source(src_rad[ip],
                        dsrc_rad[ip],
                        K[ip],
                        a,
                        S,
                        dK_dx[ip],
                        da_dx,
                        dS_dx,
                        B,
                        dB_dT,
                        jacobian_quantities,
                        jacobian_do);
      } catch (const std::runtime_error& e) {
        ostringstream os;
        os << "Runtime-error in source calculation at index " << ip
           << ": \n";
        os << e.what();
#pragma omp critical(iyEmissionStandard_source)
        {
          do_abort = true;
          fail_msg.push_back(os.str());
        }
      }
    }

#pragma omp parallel for if (!arts_omp_in_parallel())
    for (Index ip = 1; ip < np; ip++) {
      if (do_abort) continue;
      try {
        const Numeric dr_dT_past =
            do_hse ? ppath.lstep[ip - 1] / (2.0 * ppvar_t[ip - 1]) : 0;
        const Numeric dr_dT_this =
            do_hse ? ppath.lstep[ip - 1] / (2.0 * ppvar_t[ip]) : 0;
        stepwise_transmission(lyr_tra[ip],
                              dlyr_tra_above[ip],
                              dlyr_tra_below[ip],
                              K[ip - 1],
                              K[ip],
                              dK_dx[ip - 1],
                              dK_dx[ip],
                              ppath.lstep[ip - 1],
                              dr_dT_past,
                              dr_dT_this,
                              temperature_derivative_position);
        
        r[ip - 1] = ppath.lstep[ip - 1];
        if (temperature_derivative_position >= 0){
          dr_below[ip][temperature_derivative_position] = dr_dT_past;
          dr_above[ip][temperature_derivative_position] = dr_dT_this;
        }
      } catch (const std::runtime_error& e) {
        ostringstream os;
        os << "Runtime-error in transmission calculation at index " << ip
           << ": \n";
        os << e.what();
#pragma omp critical(iyEmissionStandard_transmission)
        {
          do_abort = true;
          fail_msg.push_back(os.str());
        }
      }
    }

    ARTS_USER_ERROR_IF (do_abort,
      "Error messages from failed cases:\n", fail_msg)
  }

  const ArrayOfTransmissionMatrix tot_tra =
      cumulative_transmission(lyr_tra, CumulativeTransmission::Forward);

  // iy_transmittance
  Tensor3 iy_trans_new;
  if (iy_agenda_call1)
    iy_trans_new = tot_tra[np - 1];
  else
    iy_transmittance_mult(iy_trans_new, iy_transmittance, tot_tra[np - 1]);

  // iy_aux: Optical depth
  if (auxOptDepth >= 0)
    for (Index iv = 0; iv < nf; iv++)
      iy_aux[auxOptDepth](iv, 0) = -std::log(tot_tra[np - 1](iv, 0, 0));

  // Radiative background
  get_iy_of_background(ws,
                       iy,
                       diy_dx,
                       iy_trans_new,
                       iy_id,
                       jacobian_do,
                       jacobian_quantities,
                       ppath,
                       rte_pos2,
                       atmosphere_dim,
                       nlte_field,
                       cloudbox_on,
                       stokes_dim,
                       f_grid,
                       iy_unit,
                       surface_props_data,
                       iy_main_agenda,
                       iy_space_agenda,
                       iy_surface_agenda,
                       iy_cloudbox_agenda,
                       iy_agenda_call1,
                       verbosity);

  lvl_rad[np - 1] = iy;

  // Radiative transfer calculations
  if (rt_integration_option == "first order" || rt_integration_option == "default") {
    for (Index ip = np - 2; ip >= 0; ip--) {
      lvl_rad[ip] = lvl_rad[ip + 1];
      update_radiation_vector(lvl_rad[ip],
                              dlvl_rad[ip],
                              dlvl_rad[ip + 1],
                              src_rad[ip],
                              src_rad[ip + 1],
                              dsrc_rad[ip],
                              dsrc_rad[ip + 1],
                              lyr_tra[ip + 1],
                              tot_tra[ip],
                              dlyr_tra_above[ip + 1],
                              dlyr_tra_below[ip + 1],
                              PropagationMatrix(),
                              PropagationMatrix(),
                              ArrayOfPropagationMatrix(),
                              ArrayOfPropagationMatrix(),
                              Numeric(),
                              Vector(),
                              Vector(),
                              0,
                              0,
                              RadiativeTransferSolver::Emission);
    }
  } else if (rt_integration_option == "second order") {
    for (Index ip = np - 2; ip >= 0; ip--) {
      lvl_rad[ip] = lvl_rad[ip + 1];
      update_radiation_vector(lvl_rad[ip],
                              dlvl_rad[ip],
                              dlvl_rad[ip + 1],
                              src_rad[ip],
                              src_rad[ip + 1],
                              dsrc_rad[ip],
                              dsrc_rad[ip + 1],
                              lyr_tra[ip + 1],
                              tot_tra[ip],
                              dlyr_tra_above[ip + 1],
                              dlyr_tra_below[ip + 1],
                              K[ip],
                              K[ip + 1],
                              dK_dx[ip + 1],
                              dK_dx[ip + 1],
                              r[ip],
                              dr_above[ip + 1],
                              dr_below[ip + 1],
                              0,
                              0,
                              RadiativeTransferSolver::LinearWeightedEmission);
    }
  } else {
    ARTS_USER_ERROR ( "Only allowed choices for *integration order* are "
                        "1 and 2.");    
  }
  
  // Copy back to ARTS external style
  iy = lvl_rad[0];
  for (Index ip = 0; ip < lvl_rad.nelem(); ip++) {
    ppvar_trans_cumulat(ip, joker, joker, joker) = tot_tra[ip];
    ppvar_trans_partial(ip, joker, joker, joker) = lyr_tra[ip];
    ppvar_iy(joker, joker, ip) = lvl_rad[ip];
    if (j_analytical_do)
      FOR_ANALYTICAL_JACOBIANS_DO(diy_dpath[iq](ip, joker, joker) =
                                      dlvl_rad[ip][iq];);
  }

  // Finalize analytical Jacobians
  if (j_analytical_do) {
    rtmethods_jacobian_finalisation(ws,
                                    diy_dx,
                                    diy_dpath,
                                    ns,
                                    nf,
                                    np,
                                    atmosphere_dim,
                                    ppath,
                                    ppvar_p,
                                    ppvar_t,
                                    ppvar_vmr,
                                    iy_agenda_call1,
                                    iy_transmittance,
                                    water_p_eq_agenda,
                                    jacobian_quantities,
                                    jac_species_i);
  }

  // Radiance unit conversions
  if (iy_agenda_call1) {
    rtmethods_unit_conversion(iy,
                              diy_dx,
                              ppvar_iy,
                              ns,
                              np,
                              f_grid,
                              ppath,
                              jacobian_quantities,
                              j_analytical_do,
                              iy_unit);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iyIndependentBeamApproximation(Workspace& ws,
                                    Matrix& iy,
                                    ArrayOfMatrix& iy_aux,
                                    Ppath& ppath,
                                    ArrayOfTensor3& diy_dx,
                                    GriddedField4& atm_fields_compact,
                                    const Index& iy_id,
                                    const Vector& f_grid,
                                    const Index& atmosphere_dim,
                                    const Vector& p_grid,
                                    const Vector& lat_grid,
                                    const Vector& lon_grid,
                                    const Vector& lat_true,
                                    const Vector& lon_true,
                                    const Tensor3& t_field,
                                    const Tensor3& z_field,
                                    const Tensor4& vmr_field,
                                    const EnergyLevelMap& nlte_field,
                                    const Tensor3& wind_u_field,
                                    const Tensor3& wind_v_field,
                                    const Tensor3& wind_w_field,
                                    const Tensor3& mag_u_field,
                                    const Tensor3& mag_v_field,
                                    const Tensor3& mag_w_field,
                                    const Index& cloudbox_on,
                                    const ArrayOfIndex& cloudbox_limits,
                                    const Tensor4& pnd_field,
                                    const Matrix& particle_masses,
                                    const Agenda& ppath_agenda,
                                    const Numeric& ppath_lmax,
                                    const Numeric& ppath_lraytrace,
                                    const Index& iy_agenda_call1,
                                    const String& iy_unit,
                                    const Tensor3& iy_transmittance,
                                    const Vector& rte_pos,
                                    const Vector& rte_los,
                                    const Vector& rte_pos2,
                                    const Index& jacobian_do,
                                    const ArrayOfString& iy_aux_vars,
                                    const Agenda& iy_independent_beam_approx_agenda,
                                    const Index& return_atm1d,
                                    const Index& skip_vmr,
                                    const Index& skip_pnd,
                                    const Index& return_masses,
                                    const Verbosity&) {
  // Throw error if unsupported features are requested
  ARTS_USER_ERROR_IF (jacobian_do,
        "Jacobians not provided by the method, *jacobian_do* "
        "must be 0.");
  ARTS_USER_ERROR_IF (!nlte_field.value.empty(),
        "This method does not yet support non-empty *nlte_field*.");
  ARTS_USER_ERROR_IF (!wind_u_field.empty(),
        "This method does not yet support non-empty *wind_u_field*.");
  ARTS_USER_ERROR_IF (!wind_v_field.empty(),
        "This method does not yet support non-empty *wind_v_field*.");
  ARTS_USER_ERROR_IF (!wind_w_field.empty(),
        "This method does not yet support non-empty *wind_w_field*.");
  ARTS_USER_ERROR_IF (!mag_u_field.empty(),
        "This method does not yet support non-empty *mag_u_field*.");
  ARTS_USER_ERROR_IF (!mag_v_field.empty(),
        "This method does not yet support non-empty *mag_v_field*.");
  ARTS_USER_ERROR_IF (!mag_w_field.empty(),
        "This method does not yet support non-empty *mag_w_field*.");
  //
  if (return_masses) {
    ARTS_USER_ERROR_IF (pnd_field.nbooks() != particle_masses.nrows(),
          "Sizes of *pnd_field* and *particle_masses* "
          "are inconsistent.");
  }
  chk_latlon_true(atmosphere_dim,
                  lat_grid,
                  lat_true,
                  lon_true);
  
  // Note that input 1D atmospheres are handled exactly as 2D and 3D, to make
  // the function totally general. And 1D must be handled for iterative calls.

  // Determine propagation path (with cloudbox deactivated) and check
  // that is OK for ICA
  //
  ppath_agendaExecute(ws,
                      ppath,
                      ppath_lmax,
                      ppath_lraytrace,
                      rte_pos,
                      rte_los,
                      rte_pos2,
                      0,
                      0,
                      f_grid,
                      ppath_agenda);
  //
  error_if_limb_ppath(ppath);

  // If scattering and sensor inside atmosphere, we need a pseudo-ppath that
  // samples altitudes not covered by main ppath. We make this second path
  // strictly vertical.
  //
  Ppath ppath2;
  //
  if (cloudbox_on && ppath.end_lstep == 0) {
    Vector los_tmp = rte_los;
    if (abs(rte_los[0]) < 90) {
      los_tmp[0] = 180;
    } else {
      los_tmp[0] = 0;
    }
    //
    ppath_agendaExecute(ws,
                        ppath2,
                        ppath_lmax,
                        ppath_lraytrace,
                        rte_pos,
                        los_tmp,
                        rte_pos2,
                        0,
                        0,
                        f_grid,
                        ppath_agenda);
  } else {
    ppath2.np = 1;
  }

  // Merge grid positions, and sort from bottom to top of atmosphere
  const Index np = ppath.np + ppath2.np - 1;
  ArrayOfGridPos gp_p(np), gp_lat(np), gp_lon(np);
  if (ppath.np > 1 &&
      ppath.pos(0, 0) >
          ppath.pos(1, 0)) {  // Ppath is sorted in downward direction
    // Copy ppath in reversed order
    for (Index i = 0; i < ppath.np; i++) {
      const Index ip = ppath.np - i - 1;
      gp_p[i] = ppath.gp_p[ip];
      if (atmosphere_dim > 1) {
        gp_lat[i] = ppath.gp_lat[ip];
        if (atmosphere_dim == 3) {
          gp_lon[i] = ppath.gp_lon[ip];
        }
      }
    }
    // Append ppath2, but skipping element [0]
    for (Index i = ppath.np; i < np; i++) {
      const Index ip = i - ppath.np + 1;
      gp_p[i] = ppath2.gp_p[ip];
      if (atmosphere_dim > 1) {
        gp_lat[i] = ppath2.gp_lat[ip];
        if (atmosphere_dim == 3) {
          gp_lon[i] = ppath2.gp_lon[ip];
        }
      }
    }
  } else {
    // Copy ppath2 in reversed order, but skipping element [0]
    for (Index i = 0; i < ppath2.np - 1; i++) {
      const Index ip = ppath2.np - i - 1;
      gp_p[i] = ppath2.gp_p[ip];
      if (atmosphere_dim > 1) {
        gp_lat[i] = ppath2.gp_lat[ip];
        if (atmosphere_dim == 3) {
          gp_lon[i] = ppath2.gp_lon[ip];
        }
      }
    }
    // Append ppath
    for (Index i = ppath2.np - 1; i < np; i++) {
      const Index ip = i - ppath2.np + 1;
      gp_p[i] = ppath.gp_p[ip];
      if (atmosphere_dim > 1) {
        gp_lat[i] = ppath.gp_lat[ip];
        if (atmosphere_dim == 3) {
          gp_lon[i] = ppath.gp_lon[ip];
        }
      }
    }
  }

  // 1D version of p_grid
  Matrix itw;
  Vector p1(np);
  ArrayOfGridPos gp0(0), gp1(1);
  interp_atmfield_gp2itw(itw, 1, gp_p, gp0, gp0);
  itw2p(p1, p_grid, gp_p, itw);

  // 1D version of lat and lon variables
  Vector lat1(0), lon1(0);
  Vector lat_true1(1), lon_true1(1);
  if (atmosphere_dim == 3) {
    gp1[0] = gp_lat[0];
    interp_atmfield_gp2itw(itw, 1, gp1, gp0, gp0);
    interp(lat_true1, itw, lat_grid, gp1);
    gp1[0] = gp_lon[0];
    interp_atmfield_gp2itw(itw, 1, gp1, gp0, gp0);
    interp(lon_true1, itw, lon_grid, gp1);
  } else if (atmosphere_dim == 2) {
    gp1[0] = gp_lat[0];
    interp_atmfield_gp2itw(itw, 1, gp1, gp0, gp0);
    interp(lat_true1, itw, lat_true, gp1);
    interp(lon_true1, itw, lon_true, gp1);
  } else {
    lat_true1[0] = lat_true[0];
    lon_true1[0] = lon_true[0];
  }

  // 2D/3D interpolation weights
  interp_atmfield_gp2itw(itw, atmosphere_dim, gp_p, gp_lat, gp_lon);

  // 1D temperature field
  Tensor3 t1(np, 1, 1);
  interp_atmfield_by_itw(
      t1(joker, 0, 0), atmosphere_dim, t_field, gp_p, gp_lat, gp_lon, itw);

  // 1D altitude field
  Tensor3 z1(np, 1, 1);
  interp_atmfield_by_itw(
      z1(joker, 0, 0), atmosphere_dim, z_field, gp_p, gp_lat, gp_lon, itw);

  // 1D VMR field
  Tensor4 vmr1(vmr_field.nbooks(), np, 1, 1);
  for (Index is = 0; is < vmr_field.nbooks(); is++) {
    interp_atmfield_by_itw(vmr1(is, joker, 0, 0),
                           atmosphere_dim,
                           vmr_field(is, joker, joker, joker),
                           gp_p,
                           gp_lat,
                           gp_lon,
                           itw);
  }

  // 1D surface altitude
  Matrix zsurf1(1, 1);
  zsurf1(0, 0) = z1(0, 0, 0);

  // 1D version of rte_pos/los
  Vector pos1(1);
  pos1[0] = rte_pos[0];
  Vector los1(1);
  los1[0] = abs(rte_los[0]);
  Vector pos2(0);
  if (rte_pos2.nelem()) {
    pos2 = rte_pos2[Range(0, rte_pos2.nelem())];
  }

  // Cloudbox variables
  //
  Index cbox_on1 = cloudbox_on;
  ArrayOfIndex cbox_lims1(0);
  Tensor4 pnd1(0, 0, 0, 0);
  //
  if (cloudbox_on) {
    // Determine what p1-levels that are inside cloudbox
    Index ifirst = np;
    Index ilast = -1;
    for (Index i = 0; i < np; i++) {
      if (is_gp_inside_cloudbox(gp_p[i],
                                gp_lat[i],
                                gp_lon[i],
                                cloudbox_limits,
                                true,
                                atmosphere_dim)) {
        if (i < ifirst) {
          ifirst = i;
        }
        ilast = i;
      }
    }

    // If no hit, deactive cloudbox
    if (ifirst == np) {
      cbox_on1 = 0;
    }

    // Otherwise set 1D cloud variables
    else {
      // We can enter the cloudbox from the side, and we need to add 1
      // level on each side to be safe (if not limits already at edges)
      //
      const Index extra_bot = ifirst == 0 ? 0 : 1;
      const Index extra_top = ilast == np - 1 ? 0 : 1;
      //
      cbox_lims1.resize(2);
      cbox_lims1[0] = ifirst - extra_bot;
      cbox_lims1[1] = ilast + extra_top;

      // pnd_field
      //
      pnd1.resize(pnd_field.nbooks(), cbox_lims1[1] - cbox_lims1[0] + 1, 1, 1);
      pnd1 = 0;   // As lowermost and uppermost level not always filled
      //
      itw.resize(1, Index(pow(2.0, Numeric(atmosphere_dim))));
      //
      for (Index i = extra_bot; i < pnd1.npages() - extra_top; i++) {
        const Index i0 = cbox_lims1[0] + i;
        ArrayOfGridPos gpc_p(1), gpc_lat(1), gpc_lon(1);
        interp_cloudfield_gp2itw(itw(0, joker),
                                 gpc_p[0],
                                 gpc_lat[0],
                                 gpc_lon[0],
                                 gp_p[i0],
                                 gp_lat[i0],
                                 gp_lon[i0],
                                 atmosphere_dim,
                                 cloudbox_limits);
        for (Index p = 0; p < pnd_field.nbooks(); p++) {
          interp_atmfield_by_itw(pnd1(p, i, 0, 0),
                                 atmosphere_dim,
                                 pnd_field(p, joker, joker, joker),
                                 gpc_p,
                                 gpc_lat,
                                 gpc_lon,
                                 itw);
        }
      }
    }
  }

  // Call sub agenda
  //
  {
    const Index adim1 = 1;
    const Numeric lmax1 = -1;
    Ppath ppath1d;
    //
    iy_independent_beam_approx_agendaExecute(ws,
                                             iy,
                                             iy_aux,
                                             ppath1d,
                                             diy_dx,
                                             iy_agenda_call1,
                                             iy_unit,
                                             iy_transmittance,
                                             iy_aux_vars,
                                             iy_id,
                                             adim1,
                                             p1,
                                             lat1,
                                             lon1,
                                             lat_true1,
                                             lon_true1,
                                             t1,
                                             z1,
                                             vmr1,
                                             zsurf1,
                                             lmax1,
                                             ppath_lraytrace,
                                             cbox_on1,
                                             cbox_lims1,
                                             pnd1,
                                             jacobian_do,
                                             f_grid,
                                             pos1,
                                             los1,
                                             pos2,
                                             iy_independent_beam_approx_agenda);
  }

  // Fill *atm_fields_compact*?
  if (return_atm1d) {
    // Sizes and allocate memory
    const Index nvmr = skip_vmr ? 0 : vmr1.nbooks();
    const Index npnd = skip_pnd ? 0 : pnd1.nbooks();
    const Index nmass = return_masses ? particle_masses.ncols() : 0;
    const Index ntot = 2 + nvmr + npnd + nmass;
    ArrayOfString field_names(ntot);
    atm_fields_compact.resize(ntot, np, 1, 1);

    // Altitudes
    field_names[0] = "Geometric altitudes";
    atm_fields_compact.data(0, joker, 0, 0) = z1(joker, 0, 0);

    // Temperature
    field_names[1] = "Temperature";
    atm_fields_compact.data(1, joker, 0, 0) = t1(joker, 0, 0);

    // VMRs
    if (nvmr) {
      for (Index i = 0; i < nvmr; i++) {
        const Index iout = 2 + i;
        ostringstream sstr;
        sstr << "VMR species " << i;
        field_names[iout] = sstr.str();
        atm_fields_compact.data(iout, joker, 0, 0) = vmr1(i, joker, 0, 0);
      }
    }

    // PNDs
    if (npnd) {
      for (Index i = 0; i < npnd; i++) {
        const Index iout = 2 + nvmr + i;
        ostringstream sstr;
        sstr << "Scattering element " << i;
        field_names[iout] = sstr.str();
        atm_fields_compact.data(iout, joker, 0, 0) = 0;
        atm_fields_compact.data(
            iout, Range(cbox_lims1[0], pnd1.npages()), 0, 0) =
            pnd1(i, joker, 0, 0);
      }
    }

    // Masses
    if (nmass) {
      for (Index i = 0; i < nmass; i++) {
        const Index iout = 2 + nvmr + npnd + i;
        ostringstream sstr;
        sstr << "Mass category " << i;
        field_names[iout] = sstr.str();
        atm_fields_compact.data(iout, joker, 0, 0) = 0;
        for (Index ip = cbox_lims1[0]; ip < pnd1.npages(); ip++) {
          for (Index is = 0; is < pnd1.nbooks(); is++) {
            atm_fields_compact.data(iout, ip, 0, 0) +=
                particle_masses(is, i) * pnd1(is, ip, 0, 0);
          }
        }
      }
    }

    // Finally, set grids and names
    //
    atm_fields_compact.set_name(
        "Data created by *iyIndependentBeamApproximation*");
    //
    atm_fields_compact.set_grid_name(GFIELD4_FIELD_NAMES,
                                     "Atmospheric quantity");
    atm_fields_compact.set_grid(GFIELD4_FIELD_NAMES, field_names);
    atm_fields_compact.set_grid_name(GFIELD4_P_GRID, "Pressure");
    atm_fields_compact.set_grid(GFIELD4_P_GRID, p1);
    atm_fields_compact.set_grid_name(GFIELD4_LAT_GRID, "Latitude");
    atm_fields_compact.set_grid(GFIELD4_LAT_GRID, lat_true1);
    atm_fields_compact.set_grid_name(GFIELD4_LON_GRID, "Longitude");
    atm_fields_compact.set_grid(GFIELD4_LON_GRID, lon_true1);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iyLoopFrequencies(Workspace& ws,
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
                       const Index& stokes_dim,
                       const Vector& f_grid,
                       const Agenda& iy_loop_freqs_agenda,
                       const Verbosity&) {
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
      iy.resize(nf, stokes_dim);
      //
      iy_aux.resize(iy_aux1.nelem());
      for (Index q = 0; q < iy_aux1.nelem(); q++) {
        iy_aux[q].resize(nf, stokes_dim);
      }
      //
      diy_dx.resize(diy_dx1.nelem());
      for (Index q = 0; q < diy_dx1.nelem(); q++) {
        diy_dx[q].resize(diy_dx1[q].npages(), nf, stokes_dim);
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
void iyMC(Workspace& ws,
          Matrix& iy,
          ArrayOfMatrix& iy_aux,
          ArrayOfTensor3& diy_dx,
          const Index& iy_agenda_call1,
          const Tensor3& iy_transmittance,
          const Vector& rte_pos,
          const Vector& rte_los,
          const ArrayOfString& iy_aux_vars,
          const Index& jacobian_do,
          const Index& atmosphere_dim,
          const Vector& p_grid,
          const Vector& lat_grid,
          const Vector& lon_grid,
          const Tensor3& z_field,
          const Tensor3& t_field,
          const Tensor4& vmr_field,
          const Vector& refellipsoid,
          const Matrix& z_surface,
          const Index& cloudbox_on,
          const ArrayOfIndex& cloudbox_limits,
          const Index& stokes_dim,
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
          const Index& t_interp_order,
          const Verbosity& verbosity) {
  // Throw error if unsupported features are requested
  ARTS_USER_ERROR_IF (atmosphere_dim != 3,
        "Only 3D atmospheres are allowed (atmosphere_dim must be 3)");
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
  iy.resize(nf, stokes_dim);
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
        iy_aux[i].resize(nf, stokes_dim);
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

  if (nf)
#pragma omp parallel for if (!arts_omp_in_parallel() && nf > 1) firstprivate(ws)
    for (Index f_index = 0; f_index < nf; f_index++) {
      if (failed) continue;

      try {
        // Seed reset for each loop. If not done, the errors
        // appear to be highly correlated.
        Index mc_seed;
        MCSetSeedFromTime(mc_seed, verbosity);

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
                  stokes_dim,
                  atmosphere_dim,
                  ppath_step_agenda,
                  ppath_lmax,
                  ppath_lraytrace,
                  iy_space_agenda,
                  surface_rtprop_agenda,
                  propmat_clearsky_agenda,
                  p_grid,
                  lat_grid,
                  lon_grid,
                  z_field,
                  refellipsoid,
                  z_surface,
                  t_field,
                  vmr_field,
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
                  t_interp_order,
                  verbosity);

        ARTS_ASSERT(y.nelem() == stokes_dim);

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

    ARTS_USER_ERROR_IF (failed, fail_msg);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iyReplaceFromAux(Matrix& iy,
                      const ArrayOfMatrix& iy_aux,
                      const ArrayOfString& iy_aux_vars,
                      const Index& jacobian_do,
                      const String& aux_var,
                      const Verbosity&) {
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
    const Tensor4& ppvar_trans_cumulat,
    const Verbosity&) {
  ppvar_optical_depth = ppvar_trans_cumulat(joker, joker, 0, 0);
  transform(ppvar_optical_depth, log, ppvar_optical_depth);
  ppvar_optical_depth *= -1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void yCalc(Workspace& ws,
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
           const Index& atmosphere_dim,
           const EnergyLevelMap& nlte_field,
           const Index& cloudbox_on,
           const Index& cloudbox_checked,
           const Index& scat_data_checked,
           const Index& sensor_checked,
           const Index& stokes_dim,
           const Vector& f_grid,
           const Matrix& sensor_pos,
           const Matrix& sensor_los,
           const Matrix& transmitter_pos,
           const Matrix& mblock_dlos_grid,
           const Sparse& sensor_response,
           const Vector& sensor_response_f,
           const ArrayOfIndex& sensor_response_pol,
           const Matrix& sensor_response_dlos,
           const String& iy_unit,
           const Agenda& iy_main_agenda,
           const Agenda& geo_pos_agenda,
           const Agenda& jacobian_agenda,
           const Index& jacobian_do,
           const ArrayOfRetrievalQuantity& jacobian_quantities,
           const ArrayOfString& iy_aux_vars,
           const Verbosity& verbosity) {
  CREATE_OUT3;

  // Basics
  //
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
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
  const Index nlos = mblock_dlos_grid.nrows();
  const Index n1y = sensor_response.nrows();
  const Index nmblock = sensor_pos.nrows();
  const Index niyb = nf * nlos * stokes_dim;

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
    out3 << "  Parallelizing mblock loop (" << nmblock << " iterations)\n";

#pragma omp parallel for firstprivate(ws)
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
                             atmosphere_dim,
                             nlte_field,
                             cloudbox_on,
                             stokes_dim,
                             f_grid,
                             sensor_pos,
                             sensor_los,
                             transmitter_pos,
                             mblock_dlos_grid,
                             sensor_response,
                             sensor_response_f,
                             sensor_response_pol,
                             sensor_response_dlos,
                             iy_unit,
                             iy_main_agenda,
                             geo_pos_agenda,
                             jacobian_agenda,
                             jacobian_do,
                             jacobian_quantities,
                             jacobian_indices,
                             iy_aux_vars,
                             verbosity,
                             mblock_index,
                             n1y,
                             j_analytical_do);
    }  // End mblock loop
  } else {
    out3 << "  Not parallelizing mblock loop (" << nmblock << " iterations)\n";

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
                             atmosphere_dim,
                             nlte_field,
                             cloudbox_on,
                             stokes_dim,
                             f_grid,
                             sensor_pos,
                             sensor_los,
                             transmitter_pos,
                             mblock_dlos_grid,
                             sensor_response,
                             sensor_response_f,
                             sensor_response_pol,
                             sensor_response_dlos,
                             iy_unit,
                             iy_main_agenda,
                             geo_pos_agenda,
                             jacobian_agenda,
                             jacobian_do,
                             jacobian_quantities,
                             jacobian_indices,
                             iy_aux_vars,
                             verbosity,
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
      const Index row0 = rowind.get_start();

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
void yCalcAppend(Workspace& ws,
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
                 const Index& atmosphere_dim,
                 const EnergyLevelMap& nlte_field,
                 const Index& cloudbox_on,
                 const Index& cloudbox_checked,
                 const Index& scat_data_checked,
                 const Index& sensor_checked,
                 const Index& stokes_dim,
                 const Vector& f_grid,
                 const Matrix& sensor_pos,
                 const Matrix& sensor_los,
                 const Matrix& transmitter_pos,
                 const Matrix& mblock_dlos_grid,
                 const Sparse& sensor_response,
                 const Vector& sensor_response_f,
                 const ArrayOfIndex& sensor_response_pol,
                 const Matrix& sensor_response_dlos,
                 const String& iy_unit,
                 const Agenda& iy_main_agenda,
                 const Agenda& geo_pos_agenda,
                 const Agenda& jacobian_agenda,
                 const Index& jacobian_do,
                 const ArrayOfString& iy_aux_vars,
                 const ArrayOfRetrievalQuantity& jacobian_quantities_copy,
                 const Index& append_instrument_wfs,
                 const Verbosity& verbosity) {
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
        atmosphere_dim,
        nlte_field,
        cloudbox_on,
        cloudbox_checked,
        scat_data_checked,
        sensor_checked,
        stokes_dim,
        f_grid,
        sensor_pos,
        sensor_los,
        transmitter_pos,
        mblock_dlos_grid,
        sensor_response,
        sensor_response_f,
        sensor_response_pol,
        sensor_response_dlos,
        iy_unit,
        iy_main_agenda,
        geo_pos_agenda,
        jacobian_agenda,
        jacobian_do,
        jacobian_quantities,
        iy_aux_vars,
        verbosity);

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
                const String& iy_unit,
                const Verbosity&) {
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
        apply_iy_unit2(J, yv, iy_unit, y_f[i0], 1, i_pol);
        jacobian(ii, joker) = transpose(J(joker, 0, joker));
      }

      // y (must be done last)
      apply_iy_unit(yv, iy_unit, y_f[i0], 1, i_pol);
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
            MatrixView(jacobian(i, joker)), yv, iy_unit, y_f[i], 1, i_pol);
      }

      // y (must be done last)
      apply_iy_unit(yv, iy_unit, y_f[i], 1, i_pol);
      y[i] = yv(0, 0);
    }
  }
}
