/* Copyright (C) 2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler(at)ltu.se>

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
  @file   m_cloudradar.cc
  @author Patrick Eriksson <patrick.eriksson@chalmers.se>
  @date   2010-10-31

  @brief  Workspace functions related to simulation of cloud radars.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>

#include "arts.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "logic.h"
#include "messages.h"
#include "montecarlo.h"
#include "propagationmatrix.h"
#include "rte.h"
#include "sensor.h"

extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric LOG10_EULER_NUMBER;

// Index of grids inside radar inversion tables
const Index GFIELD3_DB_GRID = 1;
const Index GFIELD3_T_GRID = 2;


/* Workspace method: Doxygen documentation will be auto-generated */
void iyRadarSingleScat(Workspace& ws,
                       Matrix& iy,
                       ArrayOfMatrix& iy_aux,
                       ArrayOfTensor3& diy_dx,
                       Vector& ppvar_p,
                       Vector& ppvar_t,
                       Matrix& ppvar_vmr,
                       Matrix& ppvar_wind,
                       Matrix& ppvar_mag,
                       Matrix& ppvar_pnd,
                       Matrix& ppvar_f,
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
                       const Index& scat_data_checked,
                       const ArrayOfString& iy_aux_vars,
                       const Index& jacobian_do,
                       const ArrayOfRetrievalQuantity& jacobian_quantities,
                       const Ppath& ppath,
                       const Agenda& propmat_clearsky_agenda,
                       const Agenda& water_p_eq_agenda,
                       const Agenda& iy_transmitter_agenda,
                       const Numeric& rte_alonglos_v,
                       const Index& trans_in_jacobian,
                       const Numeric& pext_scaling,
                       const Index& t_interp_order,
                       const Verbosity& verbosity [[maybe_unused]]) {
  //  Init Jacobian quantities?
  const Index j_analytical_do = jacobian_do ?
    do_analytical_jacobian<1>(jacobian_quantities) : 0;
  
  // Some basic sizes
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index ne = pnd_field.nbooks();
  const Index nq = j_analytical_do ? jacobian_quantities.nelem() : 0;
  const Index naux = iy_aux_vars.nelem();

  // Radiative background index
  const Index rbi = ppath_what_background(ppath);

  // Checks of input
  // Throw error if unsupported features are requested
  ARTS_USER_ERROR_IF (rbi < 1 || rbi > 9,
        "ppath.background is invalid. Check your "
        "calculation of *ppath*?");
  ARTS_USER_ERROR_IF (rbi == 3 || rbi == 4,
        "The propagation path ends inside or at boundary of "
        "the cloudbox.\nFor this method, *ppath* must be "
        "calculated in this way:\n   ppathCalc( cloudbox_on = 0 ).");
  if (cloudbox_on) {
    ARTS_USER_ERROR_IF (scat_data_checked != 1,
          "The scattering data must be flagged to have\n"
          "passed a consistency check (scat_data_checked=1).");
    ARTS_USER_ERROR_IF (ne != TotalNumberOfElements(scat_data),
          "*pnd_field* and *scat_data* inconsistent regarding total number of"
          " scattering elements.");
  }
  if (jacobian_do) {
    // FIXME: These needs to be tested properly
    ARTS_USER_ERROR ( "Jacobian calculations *iyActiveSingleScat* need "
                              "revision before safe to use.");
    
    ARTS_USER_ERROR_IF (dpnd_field_dx.nelem() != jacobian_quantities.nelem(),
          "*dpnd_field_dx* not properly initialized:\n"
          "Number of elements in dpnd_field_dx must be equal number of jacobian"
          " quantities.\n(Note: jacobians have to be defined BEFORE *pnd_field*"
          " is calculated/set.");
  }
  // iy_aux_vars checked below
  chk_if_in_range("pext_scaling", pext_scaling, 0, 2);

  // Transmitted signal
  //
  Matrix iy0;
  //
  iy_transmitter_agendaExecute(ws,
                               iy0,
                               f_grid,
                               ppath.pos(np - 1, Range(0, atmosphere_dim)),
                               ppath.los(np - 1, joker),
                               iy_transmitter_agenda);
  //
  ARTS_USER_ERROR_IF (iy0.ncols() != ns || iy0.nrows() != nf,
    "The size of *iy* returned from *iy_transmitter_agenda* is\n"
    "not correct:\n"
    "  expected size = [", nf, ",", stokes_dim, "]\n"
    "  size of iy    = [", iy0.nrows(), ",", iy0.ncols(), "]\n")
  for (Index iv = 0; iv < nf; iv++)
    ARTS_USER_ERROR_IF (iy0(iv, 0) != 1,
          "The *iy* returned from *iy_transmitter_agenda* "
          "must have the value 1 in the first column.");

  // Set diy_dpath if we are doing are doing jacobian calculations
  ArrayOfTensor3 diy_dpath = j_analytical_do ?
    get_standard_diy_dpath(jacobian_quantities, np, nf, ns, true) :
    ArrayOfTensor3(0);
  
  // Set the species pointers if we are doing jacobian
  const ArrayOfIndex jac_species_i = j_analytical_do ?
    get_pointers_for_analytical_species(jacobian_quantities, abs_species) :
    ArrayOfIndex(0);
  
  // Start diy_dx out if we are doing the first run and are doing jacobian
  // calculations
  diy_dx = j_analytical_do ?
    get_standard_starting_diy_dx(jacobian_quantities, np, nf, ns, true) :
    ArrayOfTensor3(0);
  
  // Checks that the scattering species are treated correctly if their
  // derivatives are needed
  const ArrayOfIndex jac_scat_i = j_analytical_do ?
    get_pointers_for_scat_species(jacobian_quantities, scat_species, cloudbox_on) :
    ArrayOfIndex(0);

  // Init iy_aux
  iy_aux.resize(naux);
  //
  Index auxBackScat = -1;
  Index auxAbSpAtte = -1;
  Index auxPartAtte = -1;
  //
  for (Index i = 0; i < naux; i++) {
    iy_aux[i].resize(nf * np, ns);

    if (iy_aux_vars[i] == "Radiative background") {
      iy_aux[i] = 0;
      iy_aux[i](joker, 0) = (Numeric)min((Index)2, rbi - 1);
    } else if( iy_aux_vars[i] == "Backscattering" ) {
      iy_aux[i]   = 0;
      auxBackScat = i;
    } else if( iy_aux_vars[i] == "Abs species extinction" ) {
      iy_aux[i]   = 0;
      auxAbSpAtte = i; 
    } else if( iy_aux_vars[i] == "Particle extinction" ) {
      iy_aux[i]   = 0;
      auxPartAtte = i;
    } else {
      ARTS_USER_ERROR (
        "The only allowed strings in *iy_aux_vars* are:\n"
        "  \"Radiative background\"\n"
        "  \"Backscattering\"\n"
        //         "  \"Optical depth\"\n"
        //         "  \"Particle extinction\"\n"
        "but you have selected: \"", iy_aux_vars[i], "\"")
    }
  }

  // Get atmospheric and radiative variables along the propagation path
  //
  ArrayOfRadiationVector lvl_rad(np, RadiationVector(nf, ns));
  ArrayOfArrayOfArrayOfRadiationVector dlvl_rad(np,
      ArrayOfArrayOfRadiationVector(np,
        ArrayOfRadiationVector(nq, RadiationVector(nf, ns))));
  //
  ArrayOfTransmissionMatrix lyr_tra(np, TransmissionMatrix(nf, ns));
  ArrayOfArrayOfTransmissionMatrix dlyr_tra_above(np,
      ArrayOfTransmissionMatrix(nq, TransmissionMatrix(nf, ns)));
  ArrayOfArrayOfTransmissionMatrix dlyr_tra_below(np,
      ArrayOfTransmissionMatrix(nq, TransmissionMatrix(nf, ns)));

  // Consider to make Pe an ArrayOfArrayOfPropagationMatrix!
  Tensor5 Pe(ne, np, nf, ns, ns, 0);

  // Various 
  ArrayOfMatrix ppvar_dpnd_dx(0);
  ArrayOfIndex clear2cloudy;
  //  Matrix scalar_ext(np,nf,0);  // Only used for iy_aux
  const Index nf_ssd = scat_data[0][0].pha_mat_data.nlibraries();
  const Index duplicate_freqs = ((nf == nf_ssd) ? 0 : 1);
  Tensor6 pha_mat_1se(nf_ssd, 1, 1, 1, ns, ns);
  Vector t_ok(1), t_array(1);
  Matrix mlos_sca(1, 2), mlos_inc(1, 2);
  Index ptype;
  EnergyLevelMap ppvar_nlte;

  if (np == 1 && rbi == 1) {  // i.e. ppath is totally outside the atmosphere:
    ppvar_p.resize(0);
    ppvar_t.resize(0);
    ppvar_vmr.resize(0, 0);
    ppvar_wind.resize(0, 0);
    ppvar_mag.resize(0, 0);
    ppvar_pnd.resize(0, 0);
    ppvar_f.resize(0, 0);
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

    // pnd_field
    if (cloudbox_on)
      get_ppath_cloudvars(clear2cloudy,
                          ppvar_pnd,
                          ppvar_dpnd_dx,
                          ppath,
                          atmosphere_dim,
                          cloudbox_limits,
                          pnd_field,
                          dpnd_field_dx);
    else {
      clear2cloudy.resize(np);
      for (Index ip = 0; ip < np; ip++)
        clear2cloudy[ip] = -1;
    }

    // Size radiative variables always used
    PropagationMatrix K_this(nf, ns), K_past(nf, ns), Kp(nf, ns);
    StokesVector a(nf, ns), S(nf, ns);
    ArrayOfIndex lte(np);

    // Init variables only used if transmission part of jacobian
    Vector dB_dT(0);
    ArrayOfPropagationMatrix dK_this_dx(0), dK_past_dx(0), dKp_dx(0);
    ArrayOfStokesVector da_dx(0), dS_dx(0);

    // HSE variables
    Index temperature_derivative_position = -1;
    bool do_hse = false;

    if (trans_in_jacobian && j_analytical_do) {
      dK_this_dx.resize(nq);
      dK_past_dx.resize(nq);
      dKp_dx.resize(nq);
      da_dx.resize(nq);
      dS_dx.resize(nq);
      dB_dT.resize(nf);
      FOR_ANALYTICAL_JACOBIANS_DO(
          dK_this_dx[iq] = PropagationMatrix(nf, ns);
          dK_past_dx[iq] = PropagationMatrix(nf, ns);
          dKp_dx[iq] = PropagationMatrix(nf, ns);
          da_dx[iq] = StokesVector(nf, ns);
          dS_dx[iq] = StokesVector(nf, ns);
          if (jacobian_quantities[iq] == Jacobian::Atm::Temperature) {
            temperature_derivative_position = iq;
            do_hse = jacobian_quantities[iq].Subtag() == "HSE on";
          })
    }

    // Loop ppath points and determine radiative properties
    for (Index ip = 0; ip < np; ip++) {
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
                                    trans_in_jacobian && j_analytical_do);

      if (trans_in_jacobian && j_analytical_do)
        adapt_stepwise_partial_derivatives(
            dK_this_dx,
            dS_dx,
            jacobian_quantities,
            ppvar_f(joker, ip),
            ppath.los(ip, joker),
            lte[ip],
            atmosphere_dim,
            trans_in_jacobian && j_analytical_do);

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
                                        trans_in_jacobian && jacobian_do);

        if (abs(pext_scaling - 1) > 1e-6) {
          Kp *= pext_scaling;
          if (trans_in_jacobian && j_analytical_do) {
            FOR_ANALYTICAL_JACOBIANS_DO(dKp_dx[iq] *= pext_scaling;)
          }
        }

        // Some iy_aux quantities
        if (auxAbSpAtte >= 0) {
          for (Index iv = 0; iv < nf; iv++) {
            iy_aux[auxAbSpAtte](iv*np+ip, joker) = K_this.Kjj()[iv];
          }
        }
        if (auxPartAtte >= 0) {
          for (Index iv = 0; iv < nf; iv++) {
            iy_aux[auxPartAtte](iv*np+ip, joker) = Kp.Kjj()[iv];
          }
        }
        
        K_this += Kp;

        if (trans_in_jacobian && j_analytical_do)
          FOR_ANALYTICAL_JACOBIANS_DO(dK_this_dx[iq] += dKp_dx[iq];);

        // Get back-scattering per particle, where relevant
        {
          // Direction of outgoing scattered radiation (which is reverse to LOS).
          Vector los_sca;
          mirror_los(los_sca, ppath.los(ip, joker), atmosphere_dim);
          mlos_sca(0, joker) = los_sca;

          // Obtain a length-2 vector for incoming direction
          Vector los_inc;
          if (atmosphere_dim == 3) {
            los_inc = ppath.los(ip, joker);
          } else { // Mirror back to get a correct 3D LOS
            mirror_los(los_inc, los_sca, 3);
          }
          mlos_inc(0, joker) = los_inc;

          Index i_se_flat = 0;
          for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++)
            for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++) {
              // determine whether we have some valid pnd for this
              // scatelem (in pnd or dpnd)
              Index val_pnd = 0;
              if (ppvar_pnd(i_se_flat, ip) != 0)
                val_pnd = 1;
              else if (j_analytical_do)
                for (Index iq = 0; iq < nq && !val_pnd; iq++)
                  if (jac_scat_i[iq] >= 0)
                    if (ppvar_dpnd_dx[iq](i_se_flat, ip) != 0) {
                      val_pnd = 1;
                      break;
                    }
              if (val_pnd) {
                pha_mat_1ScatElem(pha_mat_1se,
                                  ptype,
                                  t_ok,
                                  scat_data[i_ss][i_se],
                                  ppvar_t[Range(ip, 1)],
                                  mlos_sca,
                                  mlos_inc,
                                  0,
                                  t_interp_order);
                if (t_ok[0] not_eq 0)
                  if (duplicate_freqs)
                    for (Index iv = 0; iv < nf; iv++)
                      Pe(i_se_flat, ip, iv, joker, joker) =
                          pha_mat_1se(0, 0, 0, 0, joker, joker);
                  else
                    Pe(i_se_flat, ip, joker, joker, joker) =
                        pha_mat_1se(joker, 0, 0, 0, joker, joker);
                else {
                  ARTS_USER_ERROR (
                    "Interpolation error for (flat-array) scattering"
                    " element #", i_se_flat, "\n"
                    "at location/temperature point #", ip, "\n")
                }
              }
              i_se_flat++;
            } // for i_ss

        }  // local scope
      }    // clear2cloudy

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

      swap(K_past, K_this);
      swap(dK_past_dx, dK_this_dx);
    }
  }

  const ArrayOfTransmissionMatrix tot_tra_forward =
      cumulative_transmission(lyr_tra, CumulativeTransmission::Reverse);
  const ArrayOfTransmissionMatrix tot_tra_reflect =
      cumulative_transmission(lyr_tra, CumulativeTransmission::Forward);
  const ArrayOfTransmissionMatrix reflect_matrix =
      bulk_backscatter(Pe, ppvar_pnd);
  const ArrayOfArrayOfTransmissionMatrix dreflect_matrix =
      bulk_backscatter_derivative(Pe, ppvar_dpnd_dx);

  lvl_rad[0] = iy0;
  RadiationVector rad_inc = RadiationVector(nf, ns);
  rad_inc = iy0;
  set_backscatter_radiation_vector(lvl_rad,
                                   dlvl_rad,
                                   rad_inc,
                                   lyr_tra,
                                   tot_tra_forward,
                                   tot_tra_reflect,
                                   reflect_matrix,
                                   dlyr_tra_above,
                                   dlyr_tra_below,
                                   dreflect_matrix,
                                   BackscatterSolver::CommutativeTransmission);


  // Fill iy and diy_dpath
  //
  iy.resize(nf * np, ns);  // iv*np + ip is the desired output order...
  //
  for (Index ip = 0; ip < np; ip++) {
    for (Index iv = 0; iv < nf; iv++) {
      for (Index is = 0; is < stokes_dim; is++) {
        iy(iv*np+ip, is) = lvl_rad[ip](iv, is);
        if (j_analytical_do) {
          FOR_ANALYTICAL_JACOBIANS_DO(for (Index ip2 = 0; ip2 < np; ip2++)
                                          diy_dpath[iq](ip, iv*np+ip2, is) =
                                              dlvl_rad[ip][ip2][iq](iv, is););
        }
      }
    }
  }

  // Finalize analytical Jacobian
  if (j_analytical_do) {
    const Index iy_agenda_call1 = 1;
    const Tensor3 iy_transmittance(0, 0, 0);

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

  // Remaining part of iy_aux
  if (auxBackScat >= 0) {
    for (Index ip = 0; ip < np; ip++) {
      for (Index iv = 0; iv < nf; iv++) {
        VectorView stokesvec = VectorView(iy_aux[auxBackScat](iv*np+ip, joker));
        MapToEigen(stokesvec).noalias() = reflect_matrix[ip].Mat(iv) * rad_inc.Vec(iv);
      }
    }
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void yRadar(Workspace& ws,
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
            const String& iy_unit_radar,
            const ArrayOfString& iy_aux_vars,
            const Index& stokes_dim,
            const Vector& f_grid,
            const Index& atmosphere_dim,
            const Index& cloudbox_on,
            const Index& cloudbox_checked,
            const Matrix& sensor_pos,
            const Matrix& sensor_los,
            const Index& sensor_checked,
            const Index& jacobian_do,
            const ArrayOfRetrievalQuantity& jacobian_quantities,
            const Agenda& iy_radar_agenda,
            const Agenda& geo_pos_agenda,
            const ArrayOfArrayOfIndex& instrument_pol_array,
            const Vector& range_bins,
            const Numeric& ze_tref,
            const Numeric& k2,
            const Numeric& dbze_min,
            const Verbosity&) {
  // Important sizes
  const Index npos = sensor_pos.nrows();
  const Index nbins = range_bins.nelem() - 1;
  const Index nf = f_grid.nelem();
  const Index naux = iy_aux_vars.nelem();

  //---------------------------------------------------------------------------
  // Input checks
  //---------------------------------------------------------------------------

  // Basics
  //
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  //
  ARTS_USER_ERROR_IF (f_grid.empty(), "The frequency grid is empty.");
  chk_if_increasing("f_grid", f_grid);
  ARTS_USER_ERROR_IF (f_grid[0] <= 0, "All frequencies in *f_grid* must be > 0.");
  //
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  ARTS_USER_ERROR_IF (atmfields_checked != 1,
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");
  ARTS_USER_ERROR_IF (atmgeom_checked != 1,
        "The atmospheric geometry must be flagged to have "
        "passed a consistency check (atmgeom_checked=1).");
  ARTS_USER_ERROR_IF (cloudbox_checked != 1,
        "The cloudbox must be flagged to have "
        "passed a consistency check (cloudbox_checked=1).");
  ARTS_USER_ERROR_IF (sensor_checked != 1,
        "The sensor variables must be flagged to have "
        "passed a consistency check (sensor_checked=1).");

  // Method specific variables
  bool is_z = max(range_bins) > 1;
  ARTS_USER_ERROR_IF (!is_increasing(range_bins),
        "The vector *range_bins* must contain strictly "
        "increasing values.");
  ARTS_USER_ERROR_IF (!is_z && min(range_bins) < 0,
        "The vector *range_bins* is not allowed to contain "
        "negative times.");
  ARTS_USER_ERROR_IF (instrument_pol_array.nelem() != nf,
        "The main length of *instrument_pol_array* must match "
        "the number of frequencies.");

  // iy_unit_radar and variables to handle conversion to Ze and dBZe
  Vector cfac(nf, 1.0);
  Numeric ze_min = 0;
  const Numeric jfac = 10 * log10(exp(1.0));
  if (iy_unit_radar == "1") {
  } else if (iy_unit_radar == "Ze") {
    ze_cfac(cfac, f_grid, ze_tref, k2);
  } else if (iy_unit_radar == "dBZe") {
    ze_cfac(cfac, f_grid, ze_tref, k2);
    ze_min = pow(10.0, dbze_min / 10);
  } else {
    ARTS_USER_ERROR (
        "For this method, *iy_unit_radar* must be set to \"1\", "
        "\"Ze\" or \"dBZe\".");
  }

  //---------------------------------------------------------------------------
  // Various  initializations
  //---------------------------------------------------------------------------

  // Variables to handle conversion from Stokes to instrument_pol
  ArrayOfIndex npolcum(nf + 1);
  npolcum[0] = 0;
  ArrayOfArrayOfVector W(nf);
  for (Index i = 0; i < nf; i++) {
    const Index ni = instrument_pol_array[i].nelem();
    npolcum[i + 1] = npolcum[i] + ni;
    W[i].resize(ni);
    for (Index j = 0; j < ni; j++) {
      W[i][j].resize(stokes_dim);
      stokes2pol(W[i][j], stokes_dim, instrument_pol_array[i][j], 0.5);
    }
  }

  // Resize and init *y* and *y_XXX*
  //
  const Index ntot = npos * npolcum[nf] * nbins;
  y.resize(ntot);
  y = NAN;
  y_f.resize(ntot);
  y_pol.resize(ntot);
  y_pos.resize(ntot, sensor_pos.ncols());
  y_los.resize(ntot, sensor_los.ncols());
  y_geo.resize(ntot, 5);
  y_geo = NAN;  // Will be replaced if relavant data are provided

  // y_aux
  y_aux.resize(naux);
  for (Index i = 0; i < naux; i++) {
    y_aux[i].resize(ntot);
    y_aux[i] = NAN;
  }

  // Jacobian variables
  //
  Index j_analytical_do = 0;
  Index njq = 0;
  ArrayOfArrayOfIndex jacobian_indices;
  //
  if (jacobian_do) {
    bool any_affine;
    jac_ranges_indices(jacobian_indices, any_affine, jacobian_quantities, true);

    jacobian.resize(ntot,
                    jacobian_indices[jacobian_indices.nelem() - 1][1] + 1);
    jacobian = 0;
    njq = jacobian_quantities.nelem();
    //
    FOR_ANALYTICAL_JACOBIANS_DO(j_analytical_do = 1;)
  } else {
    jacobian.resize(0, 0);
  }

  //---------------------------------------------------------------------------
  // The calculations
  //---------------------------------------------------------------------------

  // Loop positions
  for (Index p = 0; p < npos; p++) {
    // RT part
    ArrayOfTensor3 diy_dx;
    Matrix iy;
    Ppath ppath;
    ArrayOfMatrix iy_aux;
    const Index iy_id = (Index)1e6 * p;
    //
    iy_radar_agendaExecute(ws,
                          iy,
                          iy_aux,
                          ppath,
                          diy_dx,
                          iy_aux_vars,
                          iy_id,
                          cloudbox_on,
                          jacobian_do,
                          sensor_pos(p, joker),
                          sensor_los(p, joker),
                          iy_radar_agenda);

    // Check if path and size OK
    const Index np = ppath.np;
    ARTS_USER_ERROR_IF (np == 1,
          "A path consisting of a single point found. "
          "This is not allowed.");
    error_if_limb_ppath(ppath);
    ARTS_USER_ERROR_IF (iy.nrows() != nf * np,
          "The size of *iy* returned from *iy_radar_agenda* "
          "is not correct (for this method).");

    // Range of ppath, in altitude or time
    Vector range(np);
    if (is_z) {
      range = ppath.pos(joker, 0);
    } else {  // Calculate round-trip time
      range[0] = 2 * ppath.end_lstep / SPEED_OF_LIGHT;
      for (Index i = 1; i < np; i++) {
        range[i] = range[i - 1] + ppath.lstep[i - 1] *
                                      (ppath.ngroup[i - 1] + ppath.ngroup[i]) /
                                      SPEED_OF_LIGHT;
      }
    }
    const Numeric range_end1 = min(range[0], range[np - 1]);
    const Numeric range_end2 = max(range[0], range[np - 1]);

    // Loop radar bins
    for (Index b = 0; b < nbins; b++) {
      if (!(range_bins[b] >= range_end2 ||     // Otherwise bin totally outside
            range_bins[b + 1] <= range_end1))  // range of ppath
      {
        // Bin limits
        Numeric blim1 = max(range_bins[b], range_end1);
        Numeric blim2 = min(range_bins[b + 1], range_end2);

        // Determine weight vector to obtain mean inside bin
        Vector hbin(np);
        integration_bin_by_vecmult(hbin, range, blim1, blim2);
        // The function above handles integration over the bin, while we
        // want the average, so divide weights with bin width
        hbin /= (blim2 - blim1);

        for (Index iv = 0; iv < nf; iv++) {
          // Pick out part of iy for frequency
          Matrix I = iy(Range(iv * np, np), joker);
          ArrayOfTensor3 dI(njq);
          if (j_analytical_do) {
            FOR_ANALYTICAL_JACOBIANS_DO(
                dI[iq] = diy_dx[iq](joker, Range(iv * np, np), joker);)
          }
          ArrayOfMatrix A(naux);
          for (Index a = 0; a < naux; a++) {
            A[a] = iy_aux[a](Range(iv * np, np), joker);
          }

          // Variables to hold data for one freq and one pol
          Vector refl(np);
          ArrayOfMatrix drefl(njq);
          if (j_analytical_do) {
            FOR_ANALYTICAL_JACOBIANS_DO(drefl[iq].resize(dI[iq].npages(), np);)
          }
          ArrayOfVector auxvar(naux);
          for (Index a = 0; a < naux; a++) {
            auxvar[a].resize(np);
          }

          for (Index ip = 0; ip < instrument_pol_array[iv].nelem(); ip++) {
            // Apply weights on each Stokes element
            mult(refl, I, W[iv][ip]);
            if (j_analytical_do) {
              FOR_ANALYTICAL_JACOBIANS_DO(for (Index k = 0;
                                               k < drefl[iq].nrows();
                                               k++) {
                mult(drefl[iq](k, joker), dI[iq](k, joker, joker), W[iv][ip]);
              })
            }
            for (Index a = 0; a < naux; a++) {
              if (iy_aux_vars[a] == "Backscattering") {
                mult(auxvar[a], A[a], W[iv][ip]);
              } else {
                for (Index j = 0; j < np; j++) {
                  auxvar[a][j] = A[a](j, 0);
                }
              }
            }

            // Apply bin weight vector to get final values.
            Index iout = nbins * (p * npolcum[nf] + npolcum[iv] + ip) + b;
            y[iout] = cfac[iv] * (hbin * refl);
            //
            if (j_analytical_do) {
              FOR_ANALYTICAL_JACOBIANS_DO(
                  for (Index k = 0; k < drefl[iq].nrows(); k++) {
                    jacobian(iout, jacobian_indices[iq][0] + k) =
                        cfac[iv] * (hbin * drefl[iq](k, joker));
                    if (iy_unit_radar == "dBZe") {
                      jacobian(iout, jacobian_indices[iq][0] + k) *=
                          jfac / max(y[iout], ze_min);
                    }
                  })
            }

            if (iy_unit_radar == "dBZe") {
              y[iout] = y[iout] <= ze_min ? dbze_min : 10 * log10(y[iout]);
            }

            // Same for aux variables
            for (Index a = 0; a < naux; a++) {
              if (iy_aux_vars[a] == "Backscattering") {
                y_aux[a][iout] = cfac[iv] * (hbin * auxvar[a]);
                if (iy_unit_radar == "dBZe") {
                  y_aux[a][iout] = y_aux[a][iout] <= ze_min
                                       ? dbze_min
                                       : 10 * log10(y_aux[a][iout]);
                }
              } else {
                y_aux[a][iout] = hbin * auxvar[a];
              }
            }
          }
        }  // Frequency
      }
    }

    // Other aux variables
    //
    Vector geo_pos;
    geo_pos_agendaExecute(ws, geo_pos, ppath, geo_pos_agenda);
    ARTS_USER_ERROR_IF (geo_pos.nelem() && geo_pos.nelem() != atmosphere_dim,
          "Wrong size of *geo_pos* obtained from "
          "*geo_pos_agenda*.\nThe length of *geo_pos* must "
          "be zero or equal to *atmosphere_dim*.");
    //
    for (Index b = 0; b < nbins; b++) {
      for (Index iv = 0; iv < nf; iv++) {
        for (Index ip = 0; ip < instrument_pol_array[iv].nelem(); ip++) {
          const Index iout = nbins * (p * npolcum[nf] + npolcum[iv] + ip) + b;
          y_f[iout] = f_grid[iv];
          y_pol[iout] = instrument_pol_array[iv][ip];
          y_pos(iout, joker) = sensor_pos(p, joker);
          y_los(iout, joker) = sensor_los(p, joker);
          if (geo_pos.nelem()) {
            y_geo(iout, joker) = geo_pos;
          }
        }
      }
    }
  }
}



/* Workspace method: Doxygen documentation will be auto-generated */
void particle_bulkpropRadarOnionPeeling(
    Workspace& ws,
    Tensor4& particle_bulkprop_field,
    ArrayOfString& particle_bulkprop_names,
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Vector& lat_grid,
    const Vector& lon_grid,
    const Tensor3& t_field,
    const Tensor3& z_field,
    const Tensor4& vmr_field,
    const Matrix& z_surface,
    const Index& atmfields_checked,
    const Index& atmgeom_checked,
    const Vector& f_grid,
    const Agenda& propmat_clearsky_agenda,
    const ArrayOfString& scat_species,
    const ArrayOfGriddedField3& invtable,
    const Matrix& incangles,
    const Tensor3& dBZe,
    const Numeric& dbze_noise,
    const Matrix& h_clutter,
    const Index& fill_clutter,
    const Numeric& t_phase,
    const Numeric& wc_max,
    const Numeric& wc_clip,
    const Index& do_atten_abs,
    const Index& do_atten_hyd,
    const Numeric& atten_hyd_scaling,
    const Numeric& atten_hyd_max,
    const Verbosity&)
{
  const Index np = t_field.npages();
  const Index nlat = t_field.nrows();
  const Index nlon = t_field.ncols();

  // Checks of input
  ARTS_USER_ERROR_IF (atmfields_checked != 1,
        "The atmospheric fields must be flagged to have\n"
        "passed a consistency check (atmfields_checked=1).");
  ARTS_USER_ERROR_IF (atmgeom_checked != 1,
        "The atmospheric geometry must be flagged to have\n"
        "passed a consistency check (atmgeom_checked=1).");
  ARTS_USER_ERROR_IF (scat_species.nelem() != 2, "Length of *scat_data* must be two.");
  ARTS_USER_ERROR_IF (invtable.nelem() != 2, "Length of *invtable* must be two.");
  ARTS_USER_ERROR_IF (dbze_noise < invtable[0].get_numeric_grid(GFIELD3_DB_GRID)[0],
                      "*dbze_noise* not covered by invtable[0]." );
  ARTS_USER_ERROR_IF (dbze_noise < invtable[1].get_numeric_grid(GFIELD3_DB_GRID)[0],
                      "*dbze_noise* not covered by invtable[1]." );    
  chk_atm_field("GIN reflectivities", dBZe, atmosphere_dim, p_grid,
                lat_grid, lon_grid);
  chk_atm_surface("GIN incangles", incangles, atmosphere_dim, lat_grid,
                  lon_grid);
  chk_if_in_range("atten_hyd_scaling", atten_hyd_scaling, 0, 2);

  // h_clutter needs special attention as it is allowed to have size 1x1.
  // And to make the code simpler below, we always switch to matrix
  Matrix hclutterm;
  if (h_clutter.nrows() > 1  || h_clutter.ncols() > 1) {
    chk_atm_surface("GIN h_clutter", h_clutter, atmosphere_dim, lat_grid, lon_grid);
    hclutterm = h_clutter;
  } else {
    hclutterm.resize(nlat, nlon);
    hclutterm = h_clutter(0,0);
  }
  
  // Init output
  particle_bulkprop_field.resize(2, np, nlat, nlon);
  particle_bulkprop_field = 0;
  particle_bulkprop_names = scat_species;

  // We apply huge extrapolation in dbze, to handle values outside the table
  // range. Extrapolation should be OK dBZe and log(IWC) have close to linear
  // relationship. 
  const Numeric extrap_fac = 100;
  
  Workspace l_ws(ws);
  Agenda l_propmat_clearsky_agenda(propmat_clearsky_agenda);
  ArrayOfString fail_msg;

  // Loop all profiles
#pragma omp parallel for if (!arts_omp_in_parallel() && nlat + nlon > 2) \
    firstprivate(l_ws, l_propmat_clearsky_agenda) collapse(2)
  for (Index ilat = 0; ilat < nlat; ilat++) {
    for (Index ilon = 0; ilon < nlon; ilon++) {
      if (fail_msg.nelem() != 0) continue;
      try {
        // Check if any significant reflectivity at all in profile
        if (max(dBZe(joker, ilat, ilon)) > dbze_noise) {
          Numeric dbze_corr_abs = 0;  // Corrections for 2-way attenuation
          Numeric dbze_corr_hyd = 0;
          Numeric k_part_above = 0, k_abs_above = 0;
          // Take abs below if incangle wrongly given as viewing angle
          const Numeric hfac = abs(1 / cos(DEG2RAD * incangles(ilat, ilon)));

          for (Index ip = np - 1; ip >= 0; ip--) {
            // Above clutter zone
            if (z_field(ip, ilat, ilon) >= z_surface(ilat, ilon) +
                                           hclutterm(ilat, ilon)) {
              // Local dBZe, roughly corrected with attenuation for previous point
              Numeric dbze =
                  dBZe(ip, ilat, ilon) + dbze_corr_abs + dbze_corr_hyd;
                                
              // Local temperature
              const Numeric t = t_field(ip, ilat, ilon);

              // Prepare for interpolation of invtable, when needed
              Index phase = -1, it = -1;  // Dummy values
              if (dbze > dbze_noise) {
                // Liquid or ice?
                phase = t >= t_phase ? 0 : 1;
                
                // Find closest temperature
                const Index tlast =
                    invtable[phase].get_numeric_grid(GFIELD3_T_GRID).nelem() -
                    1;
                if (t <= invtable[phase].get_numeric_grid(GFIELD3_T_GRID)[0])
                  it = 0;
                else if (t >= invtable[phase].get_numeric_grid(
                                  GFIELD3_T_GRID)[tlast])
                  it = tlast;
                else {
                  GridPos gp;
                  gridpos(
                      gp, invtable[phase].get_numeric_grid(GFIELD3_T_GRID), t);
                  it = gp.fd[0] < 0.5 ? gp.idx : gp.idx + 1;
                }
              }

              // Calculate attenuation from previous point
              if (ip < np - 1) {
                // Attenuation due particles
                if (dBZe(ip, ilat, ilon) > dbze_noise && do_atten_hyd &&
                    dbze_corr_hyd < atten_hyd_max) {
                  // Get extinction
                  ARTS_USER_ERROR_IF ( dbze - 20 >
                    last(invtable[phase].get_numeric_grid(GFIELD3_DB_GRID)),
                    "Max 20 dB extrapolation allowed. Found a value\n"
                    "exceeding this limit.");                  
                  GridPos gp;
                  gridpos(gp,
                          invtable[phase].get_numeric_grid(GFIELD3_DB_GRID),
                          dbze,
                          extrap_fac);
                  Vector itw(2);
                  interpweights(itw, gp);
                  Numeric k_this =
                      interp(itw, invtable[phase].data(1, joker, it), gp);
                  // Optical thickness
                  Numeric tau =
                      0.5 * hfac * (k_part_above + k_this) *
                      (z_field(ip + 1, ilat, ilon) - z_field(ip, ilat, ilon));
                  // This equals -10*log10(exp(-tau)^2)
                  dbze_corr_hyd +=
                      20 * LOG10_EULER_NUMBER * (atten_hyd_scaling * tau);
                  // Make sure we don't pass max value
                  dbze_corr_hyd = min(dbze_corr_hyd, atten_hyd_max);
                  // k_this can now be shifted to be "above value"
                  k_part_above = k_this;
                } else {  // To handle case: dBZe(ip,ilat,ilon) <= dbze_noise
                  k_part_above = 0;
                }

                // Attenuation due to abs_species
                if (do_atten_abs) {
                  // Calculate local attenuation
                  Numeric k_this;
                  PropagationMatrix propmat;
                  ArrayOfPropagationMatrix partial_dummy;
                  StokesVector nlte_dummy;
                  ArrayOfStokesVector partial_nlte_dummy;
                  EnergyLevelMap rtp_nlte_local_dummy;
                  propmat_clearsky_agendaExecute(
                      l_ws,
                      propmat,
                      nlte_dummy,
                      partial_dummy,
                      partial_nlte_dummy,
                      ArrayOfRetrievalQuantity(0),
                      {},
                      f_grid,
                      Vector(3, 0),
                      Vector(0),
                      p_grid[ip],
                      t_field(ip, ilat, ilon),
                      rtp_nlte_local_dummy,
                      vmr_field(joker, ip, ilat, ilon),
                      l_propmat_clearsky_agenda);
                  k_this = propmat.Kjj()[0];
                  // Optical thickness
                  Numeric tau =
                      0.5 * hfac * (k_abs_above + k_this) *
                      (z_field(ip + 1, ilat, ilon) - z_field(ip, ilat, ilon));
                  // This equals -10*log10(exp(-tau)^2)
                  dbze_corr_abs += 20 * LOG10_EULER_NUMBER * tau;
                  // k_this can now be shifted to be "above value"
                  k_abs_above = k_this;
                }
              }

              // Invert
              if (dBZe(ip, ilat, ilon) > dbze_noise) {
                // Correct reflectivity with (updated) attenuation
                dbze = dBZe(ip, ilat, ilon) + dbze_corr_abs + dbze_corr_hyd;

                // Interpolate inversion table (table holds log10 of water content)
                ARTS_USER_ERROR_IF ( dbze - 20 >
                    last(invtable[phase].get_numeric_grid(GFIELD3_DB_GRID)),
                    "Max 20 dB extrapolation allowed. Found a value\n"
                    "exceeding this limit.");                  
                GridPos gp;
                gridpos(gp,
                        invtable[phase].get_numeric_grid(GFIELD3_DB_GRID),
                        dbze,
                        extrap_fac);
                Vector itw(2);
                interpweights(itw, gp);
                Numeric wc =
                    pow(10.0,
                        (interp(itw, invtable[phase].data(0, joker, it), gp)));
                // Apply max and clip values
                if (wc > wc_max)
                  wc = 0;
                else if (wc > wc_clip)
                  wc = wc_clip;
                particle_bulkprop_field(phase, ip, ilat, ilon) = wc;
              }

              // In clutter zone or below surface
            } else {
              if (fill_clutter) {
                particle_bulkprop_field(0, ip, ilat, ilon) =
                    particle_bulkprop_field(0, ip + 1, ilat, ilon);
                particle_bulkprop_field(1, ip, ilat, ilon) =
                    particle_bulkprop_field(1, ip + 1, ilat, ilon);
              }
            }

          }  // presuure
        }    // R > 0
      } catch (const std::exception& e) {
        ostringstream os;
        os << "Run-time error at ilat: " << ilat << " ilon: " << ilon << ":\n"
           << e.what();
#pragma omp critical(particle_bulkpropRadarOnionPeeling_latlon)
        fail_msg.push_back(os.str());
      }
    }  // lon
  }    // lat

  if (fail_msg.nelem()) {
    ostringstream os;
    for (const auto& msg : fail_msg) os << msg << '\n';
    ARTS_USER_ERROR(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void RadarOnionPeelingTableCalc(
    Workspace& ws,
    ArrayOfGriddedField3& invtable,
    const Vector& f_grid,
    const ArrayOfString& scat_species,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const ArrayOfArrayOfScatteringMetaData& scat_meta,
    const ArrayOfAgenda& pnd_agenda_array,
    const ArrayOfArrayOfString& pnd_agenda_array_input_names,
    const Index& i_species,
    const Vector& dbze_grid,
    const Vector& t_grid,
    const Numeric& wc_min,
    const Numeric& wc_max,
    const Numeric& ze_tref,
    const Numeric& k2,
    const Verbosity& verbosity)
{
  // Some index and sizes
  const Index nss = scat_data.nelem();
  const Index ndb = dbze_grid.nelem();
  const Index nt = t_grid.nelem();
  const Index iss = i_species;

  // Check input
  ARTS_USER_ERROR_IF (scat_data[0][0].T_grid[0] < 220,
                      "First element in T_grid of scattering species is "
                      "below 220 K.\nThat seems to be too low for liquid.\n"
                      "Please note that the onion peeling function expects "
                      "rain\nas the first scattering element.");
  ARTS_USER_ERROR_IF (f_grid.nelem() != 1,
                      "This method requires that *f_grid* has length 1.");
  ARTS_USER_ERROR_IF (i_species < 0 || i_species > 1,
                      "*i_species* must either be 0 or 1.");
  ARTS_USER_ERROR_IF (nss != 2,
                      "*scat_data* must contain data for exactly two "
                      "scattering species.");
  ARTS_USER_ERROR_IF (scat_species.nelem() != nss,
        "*scat_data* and *scat_species* are inconsistent in size.");
  ARTS_USER_ERROR_IF (scat_meta.nelem() != nss,
        "*scat_data* and *scat_meta* are inconsistent in size.");
  ARTS_USER_ERROR_IF (scat_data[iss].nelem() != scat_meta[iss].nelem(),
                      "*scat_data* and *scat_meta* have inconsistent sizes.");
  ARTS_USER_ERROR_IF (scat_data[iss][0].f_grid.nelem() != 1,
                      "*scat_data* should just contain one frequency.");
  ARTS_USER_ERROR_IF (pnd_agenda_array_input_names[iss].nelem() != 1,
                      "The PSD applied must be of 1-moment type.");
  
  // Allocate
  if (invtable.empty())
    invtable.resize(2);
  invtable[iss].set_name("Radar inversion table created by *RadarOnionPeelingTableCalc*");
  invtable[iss].resize(2, ndb, nt);
  invtable[iss].data = 0;
  invtable[iss].set_grid_name(0, "Radiative properties");
  invtable[iss].set_grid(0, {"Log of water content","Extinction"});
  invtable[iss].set_grid_name(1, "Radar reflectivity");
  invtable[iss].set_grid(1, dbze_grid);
  invtable[iss].set_grid_name(2, "Temperature");
  invtable[iss].set_grid(2, t_grid);  

  // Determine back-scattering and extinction on t_grid, at each size
  const Index nse = scat_data[iss].nelem();
  Matrix b(nt,nse), e(nt,nse);
  {
    Matrix itw(nt,2);
    ArrayOfGridPos gp(nt);
    for (Index i=0; i<nse; i++) {
      ARTS_USER_ERROR_IF (scat_data[iss][i].ptype != PTYPE_TOTAL_RND,
                          "So far only TRO is handled by this method.");
      const Index ib = scat_data[iss][i].za_grid.nelem() - 1;
      ARTS_USER_ERROR_IF (scat_data[iss][i].za_grid[ib] < 179.999,
          "All za_grid in scat_data must end with 180.\n"
          "This is not the case for scat_data[", iss, "][",
          i, "] (0-based)")
      gridpos(gp, scat_data[iss][i].T_grid, t_grid, 1);
      interpweights(itw, gp);
      interp(e(joker,i), itw, scat_data[iss][i].ext_mat_data(0,joker,0,0,0), gp);
      interp(b(joker,i), itw, scat_data[iss][i].pha_mat_data(0,joker,ib,0,0,0,0), gp);
    }
  }

  // Create test grid for water content
  Vector wc_grid;
  VectorLogSpace(wc_grid, wc_min, wc_max, 0.04, verbosity);
  const Index nwc = wc_grid.nelem();
  
  // Calculate dBZe and extinction for wc_grid
  //
  Tensor3 D(2, nwc, nt, 0);
  //
  const Vector& pnd_agenda_input_t = t_grid;
  Matrix pnd_agenda_input(nt, 1);
  ArrayOfString dpnd_data_dx_names(0);
  //
  Vector cfac(1);
  ze_cfac(cfac, f_grid, ze_tref, k2);
  //
  for (Index w=0; w<nwc; w++) {
    // Get pnd
    pnd_agenda_input = wc_grid[w];
    Matrix pnd_data;
    Tensor3 dpnd_data_dx;    
    pnd_agenda_arrayExecute(ws,
                            pnd_data,
                            dpnd_data_dx,
                            iss,
                            pnd_agenda_input_t,
                            pnd_agenda_input,
                            pnd_agenda_array_input_names[iss],
                            dpnd_data_dx_names,
                            pnd_agenda_array);
    
    // Sum up to get bulk back-scattering and extinction
    for (Index t=0; t<nt; t++) {
      for (Index i=0; i<nse; i++) {
        ARTS_USER_ERROR_IF (b(t,i) < 0,
          "A negative back-scattering found for scat_species ", iss,
          ",\ntemperature ", t_grid[t], "K and element ", i)
        ARTS_USER_ERROR_IF (e(t,i) < 0,
          "A negative extinction found for scat_species ", iss, 
          ",\ntemperature ", t_grid[t], "K and element ", i)
        ARTS_USER_ERROR_IF (pnd_data(t,i) < 0,
          "A negative PSD value found for scat_species ", iss, 
          ",\ntemperature ", t_grid[t], "K and ", wc_grid[w],
          " kg/m3")
        D(0,w,t) += pnd_data(t,i) * b(t,i);
        D(1,w,t) += pnd_data(t,i) * e(t,i);
      }
      // Convert to dBZe
      D(0,w,t) = 10 * log10(cfac[0] * D(0,w,t));      
    }
  }

  // Get water content and extinction as a function of dBZe by interpolation
  Matrix itw(ndb,2);
  ArrayOfGridPos gp(ndb);
  // Water content interpolated in log
  Vector wc_log(nwc);
  transform(wc_log, log10, wc_grid);
  for (Index t=0; t<nt; t++) {
    if (!is_increasing(D(0,joker,t))) {
      for (Index w=0; w<nwc; w++) {
        cout << wc_grid[w] << " " << D(0,w,t) << endl;
      }
      ARTS_USER_ERROR (
        "A case found of non-increasing dBZe.\n"
        "Found for scat_species ", iss, " and ", t_grid[t], "K.")
    }
    if (D(0,0,t) > dbze_grid[0]) {
      for (Index w=0; w<nwc; w++) {
        cout << wc_grid[w] << " " << D(0,w,t) << endl;
      }
      ARTS_USER_ERROR (
        "A case found where start of dbze_grid not covered.\n"
        "Found for scat_species ", iss, " and ", t_grid[t], "K.")      
    }
    if (D(0,nwc-1,t) < dbze_grid[ndb-1]) {
      for (Index w=0; w<nwc; w++) {
        cout << wc_grid[w] << " " << D(0,w,t) << endl;
      }
      ARTS_USER_ERROR (
        "A case found where end of dbze_grid not covered.\n"
        "Found for scat_species ", iss, " and ", t_grid[t], "K.")
    }
    //
    gridpos(gp, D(0,joker,t), dbze_grid);
    interpweights(itw, gp);
    interp(invtable[iss].data(0,joker,t), itw, wc_log, gp);
    interp(invtable[iss].data(1,joker,t), itw, D(1,joker,t), gp);
  }
}
