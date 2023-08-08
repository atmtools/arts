/*===========================================================================
  ===  File description
  ===========================================================================*/

/**
  @file   m_transmitter.cc
  @author Patrick Eriksson <patrick.eriksson@chalmers.se>
  @date   2012-10-31

  @brief  Workspace functions related to transmitters and radiative transfer
  for transmitted signals.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include <workspace.h>
#include "geodetic.h"
#include "jacobian.h"
#include "lin_alg.h"
#include "logic.h"
#include "math_funcs.h"
#include "matpack_complex.h"
#include "rte.h"
#include "rtepack.h"
#include "rtepack_multitype.h"
#include "sensor.h"
#include <cmath>
#include <stdexcept>

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric PI=Constant::pi;
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);
inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;

/* Workspace method: Doxygen documentation will be auto-generated */
void iyTransmissionStandard(const Workspace& ws,
                            Matrix& iy,
                            ArrayOfMatrix& iy_aux,
                            ArrayOfTensor3& diy_dx,
                            ArrayOfAtmPoint& ppvar_atm,
                            Matrix& ppvar_pnd,
                            ArrayOfVector& ppvar_f,
                            Tensor3& ppvar_iy,
                            Tensor4& ppvar_trans_cumulat,
                            Tensor4& ppvar_trans_partial,
                            const Vector& f_grid,
                            const ArrayOfArrayOfSpeciesTag& abs_species,
                            const AtmField& atm_field,
                            const Index& cloudbox_on,
                            const ArrayOfIndex& cloudbox_limits,
                            const Index& gas_scattering_do,
                            const Tensor4& pnd_field,
                            const ArrayOfTensor4& dpnd_field_dx,
                            const ArrayOfString& scat_species,
                            const ArrayOfArrayOfSingleScatteringData& scat_data,
                            const ArrayOfString& iy_aux_vars,
                            const Index& jacobian_do,
                            const ArrayOfRetrievalQuantity& jacobian_quantities,
                            const Ppath& ppath,
                            const Matrix& iy_transmitter,
                            const Agenda& propmat_clearsky_agenda,
                            const Agenda& water_p_eq_agenda,
                            const Agenda& gas_scattering_agenda,
                            const Index& iy_agenda_call1,
                            const Tensor3& iy_transmittance,
                            const Numeric& rte_alonglos_v) {
  //  Init Jacobian quantities?
  const Index j_analytical_do = jacobian_do ? do_analytical_jacobian<1>(jacobian_quantities) : 0;
    
  // Some basic sizes
  const Index nf = f_grid.nelem();
  const Index np = ppath.np;
  const Index nq = j_analytical_do ? jacobian_quantities.nelem() : 0;

  // Radiative background index
  const Index rbi = 0;;//ppath_what_background(ppath);
ARTS_USER_ERROR("ERROR")

  // Checks of input
  // Throw error if unsupported features are requested
  if (!iy_agenda_call1)
    throw runtime_error(
        "Recursive usage not possible (iy_agenda_call1 must be 1)");
  if (!iy_transmittance.empty())
    throw runtime_error("*iy_transmittance* must be empty");
  if (rbi < 1 || rbi > 9)
    throw runtime_error(
        "ppath.background is invalid. Check your "
        "calculation of *ppath*?");
  if (jacobian_do) {
    if (dpnd_field_dx.nelem() != jacobian_quantities.nelem())
      throw runtime_error(
          "*dpnd_field_dx* not properly initialized:\n"
          "Number of elements in dpnd_field_dx must be equal number of jacobian"
          " quantities.\n(Note: jacobians have to be defined BEFORE *pnd_field*"
          " is calculated/set.");
  }

  ARTS_USER_ERROR_IF(jacobian_quantities.nelem() && gas_scattering_do, R"--(
Jacobian calculation are not supported when gas scattering or suns are included.
This feature will be added in a future version.
)--");

  // iy_aux_vars checked below

  // Transmitted signal
//  iy_transmitter=iy_transmitter;
  if (iy_transmitter.ncols() != 4 || iy_transmitter.nrows() != nf) {
    ostringstream os;
    os << "The size of *iy_transmitter* returned from *iy_transmitter_agenda* is\n"
       << "not correct:\n"
       << "  expected size = [" << nf << "," << 4 << "]\n"
       << "  size of iy_transmitter    = [" << iy_transmitter.nrows() << "," << iy_transmitter.ncols() << "]\n";
    throw runtime_error(os.str());
  }
  
  // Set diy_dpath if we are doing are doing jacobian calculations
  ArrayOfTensor3 diy_dpath = j_analytical_do ? get_standard_diy_dpath(jacobian_quantities, np, nf, false) : ArrayOfTensor3(0);
  
  // Set the species pointers if we are doing jacobian
  const ArrayOfIndex jac_species_i = j_analytical_do ? get_pointers_for_analytical_species(jacobian_quantities, abs_species) : ArrayOfIndex(0);
  
  // Start diy_dx out if we are doing the first run and are doing jacobian calculations
  if (j_analytical_do and iy_agenda_call1) diy_dx = get_standard_starting_diy_dx(jacobian_quantities, np, nf, false);
  
  // Checks that the scattering species are treated correctly if their derivatives are needed (we can here discard the Array)
  if (j_analytical_do and iy_agenda_call1) get_pointers_for_scat_species(jacobian_quantities, scat_species, cloudbox_on);

  // Init iy_aux and fill where possible
  const Index naux = iy_aux_vars.nelem();
  iy_aux.resize(naux);
  //
  Index auxOptDepth = -1;
  //
  for (Index i = 0; i < naux; i++) {
    iy_aux[i].resize(nf, 4);
    iy_aux[i] = 0;

    if (iy_aux_vars[i] == "Radiative background")
      iy_aux[i](joker, 0) = (Numeric)min((Index)2, rbi - 1);
    else if (iy_aux_vars[i] == "Optical depth")
      auxOptDepth = i;
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
  //
  ppvar_trans_cumulat.resize(np, nf, 4, 4);
  ppvar_trans_partial.resize(np, nf, 4, 4);

  ArrayOfMuelmatVector lyr_tra(np, MuelmatVector(nf));
  ArrayOfStokvecVector lvl_rad(np, StokvecVector(nf));
  ArrayOfStokvecMatrix dlvl_rad(np, StokvecMatrix(nq, nf));

  ArrayOfMuelmatMatrix dlyr_tra_above(np, MuelmatMatrix(nq, nf));
  ArrayOfMuelmatMatrix dlyr_tra_below(dlyr_tra_above);

  ArrayOfIndex clear2cloudy;
  //
  if (np == 1 && rbi == 1) {  // i.e. ppath is totally outside the atmosphere:
    ppvar_pnd.resize(0, 0);
    ppvar_f.resize(0, 0);
    ppvar_iy.resize(0, 0, 0);
    ppvar_trans_cumulat = 0;
    ppvar_trans_partial = 0;
    for (Index iv = 0; iv < nf; iv++) {
      for (Index is = 0; is < 4; is++) {
        ppvar_trans_cumulat(0,iv,is,is) = 1;
        ppvar_trans_partial(0,iv,is,is) = 1;        
      }
    }
    
  } else {
    // ppvar_iy
    ppvar_iy.resize(nf, 4, np);
    ppvar_iy(joker, joker, np - 1) = iy_transmitter;

    // pnd_field
    ArrayOfMatrix ppvar_dpnd_dx(0);
    //
    if (cloudbox_on)
      get_ppath_cloudvars(clear2cloudy,
                          ppvar_pnd,
                          ppvar_dpnd_dx,
                          ppath,
                          cloudbox_limits,
                          pnd_field,
                          dpnd_field_dx);
    else {
      clear2cloudy.resize(np);
      for (Index ip = 0; ip < np; ip++) clear2cloudy[ip] = -1;
    }




    // Size radiative variables always used
    PropmatVector K_this(nf), K_past(nf), Kp(nf);
    StokvecVector a(nf), S(nf), Sp(nf);

    // size gas scattering variables
    MuelmatVector gas_scattering_mat;
    Vector sca_fct_dummy;
    PropmatVector K_sca(gas_scattering_do * nf);

    Vector gas_scattering_los_in, gas_scattering_los_out;

    // Init variables only used if analytical jacobians done
    Vector dr1(nq, 0.0), dr2(nq, 0.0);
    Matrix dB_dT(nq, nf);
    PropmatMatrix dK_this_dx(nq, nf), dK_past_dx(nq, nf), dKp_dx(nq, nf);
    StokvecMatrix da_dx(nq, nf), dS_dx(nq, nf), dSp_dx(nq, nf);
    //
    Index temperature_derivative_position = -1;
    bool do_hse = false;

    if (j_analytical_do) {
      FOR_ANALYTICAL_JACOBIANS_DO(if (jacobian_quantities[iq] == Jacobian::Atm::Temperature) {
                                    temperature_derivative_position = iq;
                                    do_hse = jacobian_quantities[iq].Subtag() ==
                                             "HSE on";
                                  })
    }

    // Loop ppath points and determine radiative properties
    for (Index ip = 0; ip < np; ip++) {
      get_stepwise_clearsky_propmat(ws,
                                    K_this,
                                    S,
                                    dK_this_dx,
                                    dS_dx,
                                    propmat_clearsky_agenda,
                                    jacobian_quantities,
                                    ppvar_f[ip],
                                    Vector{ppath.los(ip, joker)},
                                    ppvar_atm[ip],
                                    j_analytical_do);

      // get Extinction from gas scattering
      if (gas_scattering_do) {
        gas_scattering_agendaExecute(ws,
                                     K_sca,
                                     gas_scattering_mat,
                                     sca_fct_dummy,
                                     ppvar_f[ip],  // FIXME: WAS THIS A BUG?
                                     ppvar_atm[ip],
                                     gas_scattering_los_in,
                                     gas_scattering_los_out,
                                     0,
                                     gas_scattering_agenda);

        K_this += K_sca;
      }

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
                                        ExhaustiveVectorView{ppvar_atm[ip].temperature},
                                        jacobian_do);
        K_this += Kp;

        if (j_analytical_do) {
          FOR_ANALYTICAL_JACOBIANS_DO(dK_this_dx[iq] += dKp_dx[iq];)
        }
      }

      // Transmission
      if (ip not_eq 0) {
        const Numeric dr_dT_past =
            do_hse ? ppath.lstep[ip - 1] / (2.0 * ppvar_atm[ip - 1].temperature) : 0;
        const Numeric dr_dT_this =
            do_hse ? ppath.lstep[ip - 1] / (2.0 * ppvar_atm[ip].temperature) : 0;
        dr1[temperature_derivative_position] = dr_dT_past;
        dr2[temperature_derivative_position] = dr_dT_this;
        two_level_exp(lyr_tra[ip], dlyr_tra_above[ip], dlyr_tra_below[ip],
                      K_past, K_this, dK_past_dx, dK_this_dx,
                      ppath.lstep[ip - 1], dr1, dr2);
      }

      swap(K_past, K_this);
      swap(dK_past_dx, dK_this_dx);
    }
  }

  const auto tot_tra = forward_cumulative_transmission(lyr_tra);

  // iy_aux: Optical depth
  if (auxOptDepth >= 0) {
    for (Index iv = 0; iv < nf; iv++)
      iy_aux[auxOptDepth](iv, 0) = -std::log(tot_tra[np - 1][iv](0, 0));
  }

  lvl_rad[np - 1] = rtepack::to_stokvec_vector(iy_transmitter);

  // Radiative transfer calculations
  for (Index ip = np - 2; ip >= 0; ip--) {
    lvl_rad[ip] = lvl_rad[ip + 1];
    two_level_linear_transmission_step(
        lvl_rad[ip], dlvl_rad[ip], dlvl_rad[ip + 1], lyr_tra[ip + 1],
        tot_tra[ip], dlyr_tra_above[ip + 1], dlyr_tra_below[ip + 1]);
  }

  // Copy back to ARTS external style
  iy = to_matrix(lvl_rad[0]);
  //
  if (np > 1) {
    for (Index ip = 0; ip < lvl_rad.nelem(); ip++) {
      ppvar_trans_cumulat(ip, joker, joker, joker) = to_tensor3(tot_tra[ip]);
      ppvar_trans_partial(ip, joker, joker, joker) = to_tensor3(lyr_tra[ip]);
      ppvar_iy(joker, joker, ip) = to_matrix(lvl_rad[ip]);
      if (j_analytical_do)
        FOR_ANALYTICAL_JACOBIANS_DO(diy_dpath[iq](ip, joker, joker) =
                                        to_matrix(dlvl_rad[ip][iq]););
    }
  }

  // Finalize analytical Jacobians
  if (j_analytical_do) {
    rtmethods_jacobian_finalisation(ws,
                                    diy_dx,
                                    diy_dpath,
                                    nf,
                                    np,
                                    ppath,
                                    ppvar_atm,
                                    iy_agenda_call1,
                                    iy_transmittance,
                                    water_p_eq_agenda,
                                    jacobian_quantities,
                                    abs_species,
                                    jac_species_i);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iy_transmitterMultiplePol(Matrix& iy_transmitter,
                               const Vector& f_grid,
                               const ArrayOfIndex& instrument_pol) {
  const Index nf = f_grid.nelem();

  if (instrument_pol.nelem() != nf)
    throw runtime_error(
        "The length of *f_grid* and the number of elements "
        "in *instrument_pol* must be equal.");

  iy_transmitter.resize(nf, 4);

  for (Index i = 0; i < nf; i++) {
    stokes2pol(iy_transmitter(i, joker), instrument_pol[i], 1);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iy_transmitterSinglePol(Matrix& iy_transmitter,
                             const Vector& f_grid,
                             const ArrayOfIndex& instrument_pol) {
  const Index nf = f_grid.nelem();

  if (instrument_pol.nelem() != 1)
    throw runtime_error(
        "The number of elements in *instrument_pol* must be 1.");

  iy_transmitter.resize(nf, 4);

  stokes2pol(iy_transmitter(0, joker), instrument_pol[0], 1);

  for (Index i = 1; i < nf; i++) {
    iy_transmitter(i, joker) = iy_transmitter(0, joker);
  }
}
