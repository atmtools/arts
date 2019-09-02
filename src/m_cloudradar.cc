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

/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   m_cloudradar.cc
  \author Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2010-10-31

  \brief  Workspace functions related to simulation of cloud radars.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "auto_md.h"
#include "logic.h"
#include "messages.h"
#include "montecarlo.h"
#include "propagationmatrix.h"
#include "rte.h"
#include "sensor.h"

extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT;
extern const String SCATSPECIES_MAINTAG;

/* Workspace method: Doxygen documentation will be auto-generated */
void iyActiveSingleScat(Workspace& ws,
                        Matrix& iy,
                        ArrayOfMatrix& iy_aux,
                        ArrayOfTensor3& diy_dx,
                        Vector& ppvar_p,
                        Vector& ppvar_t,
                        Matrix& ppvar_nlte,
                        Matrix& ppvar_vmr,
                        Matrix& ppvar_wind,
                        Matrix& ppvar_mag,
                        Matrix& ppvar_pnd,
                        Matrix& ppvar_f,
                        Tensor4& ppvar_trans_cumulat,
                        const Index& stokes_dim,
                        const Vector& f_grid,
                        const Index& atmosphere_dim,
                        const Vector& p_grid,
                        const Tensor3& t_field,
                        const Tensor4& nlte_field,
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
                        const Index& iy_agenda_call1,
                        const Tensor3& iy_transmission,
                        const Numeric& rte_alonglos_v,
                        const Index& trans_in_jacobian,
                        const Numeric& pext_scaling,
                        const Index& t_interp_order,
                        const Verbosity& verbosity _U_) {
  // Some basic sizes
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index ne = pnd_field.nbooks();

  // Radiative background index
  const Index rbi = ppath_what_background(ppath);

  // Checks of input
  // Throw error if unsupported features are requested
  if (!iy_agenda_call1)
    throw runtime_error(
        "Recursive usage not possible (iy_agenda_call1 must be 1)");
  if (!iy_transmission.empty())
    throw runtime_error("*iy_transmission* must be empty");
  if (rbi < 1 || rbi > 9)
    throw runtime_error(
        "ppath.background is invalid. Check your "
        "calculation of *ppath*?");
  if (rbi == 3 || rbi == 4)
    throw runtime_error(
        "The propagation path ends inside or at boundary of "
        "the cloudbox.\nFor this method, *ppath* must be "
        "calculated in this way:\n   ppathCalc( cloudbox_on = 0 ).");
  if (cloudbox_on) {
    if (scat_data_checked != 1)
      throw runtime_error(
          "The scattering data must be flagged to have\n"
          "passed a consistency check (scat_data_checked=1).");
    if (ne != TotalNumberOfElements(scat_data))
      throw runtime_error(
          "*pnd_field* and *scat_data* inconsistent regarding total number of"
          " scattering elements.");
  }
  if (jacobian_do) {
    if (dpnd_field_dx.nelem() != jacobian_quantities.nelem())
      throw runtime_error(
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
  if (iy0.ncols() != ns || iy0.nrows() != nf) {
    ostringstream os;
    os << "The size of *iy* returned from *iy_transmitter_agenda* is\n"
       << "not correct:\n"
       << "  expected size = [" << nf << "," << stokes_dim << "]\n"
       << "  size of iy    = [" << iy0.nrows() << "," << iy0.ncols() << "]\n";
    throw runtime_error(os.str());
  }
  for (Index iv = 0; iv < nf; iv++) {
    if (iy0(iv, 0) != 1)
      throw runtime_error(
          "The *iy* returned from *iy_transmitter_agenda* "
          "must have the value 1 in the first column.");
  }

  //  Init Jacobian quantities
  Index j_analytical_do = 0;
  if (jacobian_do) {
    FOR_ANALYTICAL_JACOBIANS_DO(j_analytical_do = 1;)
  }
  //
  const Index nq = j_analytical_do ? jacobian_quantities.nelem() : 0;
  ArrayOfTensor3 diy_dpath(nq);
  ArrayOfIndex jac_species_i(nq), jac_scat_i(nq), jac_is_t(nq), jac_wind_i(nq);
  ArrayOfIndex jac_mag_i(nq), jac_other(nq);

  if (j_analytical_do) {
    rtmethods_jacobian_init(jac_species_i,
                            jac_scat_i,
                            jac_is_t,
                            jac_wind_i,
                            jac_mag_i,
                            jac_other,
                            diy_dx,
                            diy_dpath,
                            ns,
                            nf,
                            np,
                            nq,
                            abs_species,
                            cloudbox_on,
                            scat_species,
                            dpnd_field_dx,
                            jacobian_quantities,
                            iy_agenda_call1,
                            true);
  }

  // Init iy_aux
  const Index naux = iy_aux_vars.nelem();
  iy_aux.resize(naux);
  Index auxBackScat = -1;
  Index auxOptDepth = -1;
  Index auxPartAtte = -1;
  //
  for (Index i = 0; i < naux; i++) {
    iy_aux[i].resize(nf * np, ns);
    iy_aux[i] = 0;

    if (iy_aux_vars[i] == "Radiative background") {
      iy_aux[i](joker, 0) = (Numeric)min((Index)2, rbi - 1);
    } else if (iy_aux_vars[i] == "Backscattering") {
      iy_aux[i] = 0;
      auxBackScat = i;
    } else if (iy_aux_vars[i] == "Optical depth") {
      auxOptDepth = i;
    } else if (iy_aux_vars[i] == "Particle extinction") {
      iy_aux[i](Range(0, nf, np), joker) = 0;
      auxPartAtte = i;
    } else {
      ostringstream os;
      os << "The only allowed strings in *iy_aux_vars* are:\n"
         << "  \"Radiative background\"\n"
         << "  \"Backscattering\"\n"
         << "  \"Optical depth\"\n"
         << "  \"Particle extinction\"\n"
         << "but you have selected: \"" << iy_aux_vars[i] << "\"";
      throw runtime_error(os.str());
    }
  }

  // Get atmospheric and radiative variables along the propagation path
  //
  ppvar_trans_cumulat.resize(np, nf, ns, ns);
  Tensor3 J(np, nf, ns);
  Tensor4 trans_partial(np, nf, ns, ns);
  Tensor5 dtrans_partial_dx_above(0, 0, 0, 0, 0);
  Tensor5 dtrans_partial_dx_below(0, 0, 0, 0, 0);
  Tensor5 Pe(ne, np, nf, ns, ns, 0);
  ArrayOfMatrix ppvar_dpnd_dx(0);
  ArrayOfIndex clear2cloudy;
  Matrix scalar_ext(np, nf, 0);  // Only used for iy_aux
  Index nf_ssd = scat_data[0][0].pha_mat_data.nlibraries();
  Index duplicate_freqs = ((nf == nf_ssd) ? 0 : 1);
  Tensor6 pha_mat_1se(nf_ssd, 1, 1, 1, ns, ns);
  Vector t_ok(1), t_array(1);
  Matrix pdir(1, 2), idir(1, 2);
  Index ptype;

  //
  if (np == 1 && rbi == 1)  // i.e. ppath is totally outside the atmosphere:
  {
    ppvar_p.resize(0);
    ppvar_t.resize(0);
    ppvar_vmr.resize(0, 0);
    ppvar_nlte.resize(0, 0);
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
    if (cloudbox_on) {
      get_ppath_cloudvars(clear2cloudy,
                          ppvar_pnd,
                          ppvar_dpnd_dx,
                          ppath,
                          atmosphere_dim,
                          cloudbox_limits,
                          pnd_field,
                          dpnd_field_dx);
    } else {
      clear2cloudy.resize(np);
      for (Index ip = 0; ip < np; ip++) {
        clear2cloudy[ip] = -1;
      }
    }

    // Size radiative variables always used
    PropagationMatrix K_this(nf, ns), K_past(nf, ns), Kp(nf, ns);
    StokesVector a(nf, ns), S(nf, ns), Sp(nf, ns);
    ArrayOfIndex lte(np);

    // Init variables only used if transmission part of jacobian
    Vector dB_dT(0);
    ArrayOfPropagationMatrix dK_this_dx(0), dK_past_dx(0), dKp_dx(0);
    ArrayOfStokesVector da_dx(0), dS_dx(0), dSp_dx(0);
    //
    if (trans_in_jacobian && j_analytical_do) {
      dtrans_partial_dx_above.resize(np, nq, nf, ns, ns);
      dtrans_partial_dx_below.resize(np, nq, nf, ns, ns);
      dtrans_partial_dx_above = 0;  // Here all values are set to zero, to
      dtrans_partial_dx_below = 0;  // allow looping over non-defined values
      dK_this_dx.resize(nq);
      dK_past_dx.resize(nq);
      dKp_dx.resize(nq);
      da_dx.resize(nq);
      dS_dx.resize(nq);
      dSp_dx.resize(nq);
      dB_dT.resize(nf);
      FOR_ANALYTICAL_JACOBIANS_DO(dK_this_dx[iq] = PropagationMatrix(nf, ns);
                                  dK_past_dx[iq] = PropagationMatrix(nf, ns);
                                  dKp_dx[iq] = PropagationMatrix(nf, ns);
                                  da_dx[iq] = StokesVector(nf, ns);
                                  dS_dx[iq] = StokesVector(nf, ns);
                                  dSp_dx[iq] = StokesVector(nf, ns);)
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
                                    ppvar_nlte(joker, ip),
                                    ppvar_vmr(joker, ip),
                                    ppvar_t[ip],
                                    ppvar_p[ip],
                                    jac_species_i,
                                    trans_in_jacobian && j_analytical_do);

      if (trans_in_jacobian && j_analytical_do) {
        adapt_stepwise_partial_derivatives(
            dK_this_dx,
            dS_dx,
            jacobian_quantities,
            ppvar_f(joker, ip),
            ppath.los(ip, joker),
            ppvar_vmr(joker, ip),
            ppvar_t[ip],
            ppvar_p[ip],
            jac_species_i,
            jac_wind_i,
            lte[ip],
            atmosphere_dim,
            trans_in_jacobian && j_analytical_do);
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
                                        ppvar_t[Range(ip, 1)],
                                        atmosphere_dim,
                                        trans_in_jacobian && jacobian_do);

        if (abs(pext_scaling - 1) > 1e-6) {
          Kp *= pext_scaling;
          if (trans_in_jacobian && j_analytical_do) {
            FOR_ANALYTICAL_JACOBIANS_DO(dKp_dx[iq] *= pext_scaling;)
          }
        }

        K_this += Kp;

        if (auxPartAtte >= 0) {
          scalar_ext(ip, joker) = Kp.Kjj();
        }

        if (trans_in_jacobian && j_analytical_do) {
          FOR_ANALYTICAL_JACOBIANS_DO(dK_this_dx[iq] += dKp_dx[iq];)
        }

        // Get back-scattering per particle, where relevant
        {
          // Direction of outgoing scattered radiation (which is reverse to LOS).
          Vector los_sca;
          mirror_los(los_sca, ppath.los(ip, joker), atmosphere_dim);
          pdir(0, joker) = los_sca;

          // Obtain a length-2 vector for incoming direction
          Vector los_inc;
          if (atmosphere_dim == 3) {
            los_inc = ppath.los(ip, joker);
          } else  // Mirror back to get a correct 3D LOS
          {
            mirror_los(los_inc, los_sca, 3);
          }
          idir(0, joker) = los_inc;

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
                                  pdir,
                                  idir,
                                  0,
                                  t_interp_order);
                if (t_ok[0] != 0)
                  if (duplicate_freqs)
                    for (Index iv = 0; iv < nf; iv++)
                      Pe(i_se_flat, ip, iv, joker, joker) =
                          pha_mat_1se(0, 0, 0, 0, joker, joker);
                  else
                    Pe(i_se_flat, ip, joker, joker, joker) =
                        pha_mat_1se(joker, 0, 0, 0, joker, joker);
                else {
                  ostringstream os;
                  os << "Interpolation error for (flat-array) scattering"
                     << " element #" << i_se_flat << "\n"
                     << "at location/temperature point #" << ip << "\n";
                  throw runtime_error(os.str());
                }
              }
              i_se_flat++;
            }

        }  // local scope
      }    // clear2cloudy

      if (trans_in_jacobian && j_analytical_do) {
        get_stepwise_transmission_matrix(
            ppvar_trans_cumulat(ip, joker, joker, joker),
            trans_partial(ip, joker, joker, joker),
            dtrans_partial_dx_above(ip, joker, joker, joker, joker),
            dtrans_partial_dx_below(ip, joker, joker, joker, joker),
            (ip > 0) ? ppvar_trans_cumulat(ip - 1, joker, joker, joker)
                     : Tensor3(0, 0, 0),
            K_past,
            K_this,
            dK_past_dx,
            dK_this_dx,
            (ip > 0) ? ppath.lstep[ip - 1] : Numeric(1.0),
            ip == 0);
      } else {
        get_stepwise_transmission_matrix(
            ppvar_trans_cumulat(ip, joker, joker, joker),
            trans_partial(ip, joker, joker, joker),
            Tensor4(0, 0, 0, 0),
            Tensor4(0, 0, 0, 0),
            (ip > 0) ? ppvar_trans_cumulat(ip - 1, joker, joker, joker)
                     : Tensor3(0, 0, 0),
            K_past,
            K_this,
            dK_past_dx,
            dK_this_dx,
            (ip > 0) ? ppath.lstep[ip - 1] : Numeric(1.0),
            ip == 0);
      }

      swap(K_past, K_this);
      swap(dK_past_dx, dK_this_dx);
    }
  }

  // Transmission for reversed direction
  //Tensor3 tr_rev( nf, ns, ns );
  //for( Index iv=0; iv<nf; iv++ )
  //  { id_mat( tr_rev(iv,joker,joker) ); }

  // Size iy and set to zero
  iy.resize(nf * np, ns);
  iy = 0;

  // If neither cloudbox or aux variables, nothing more to do
  if (!cloudbox_on && !naux) {
    return;
  }

  // Radiative transfer calculations
  for (Index ip = 0; ip < np; ip++) {
    // Radar return only possible if inside cloudbox
    if (clear2cloudy[ip] >= 0) {
      for (Index iv = 0; iv < nf; iv++) {
        // Calculate bulk back-scattering
        Matrix P(ns, ns, 0);
        for (Index i = 0; i < ne; i++) {
          if (ppvar_pnd(i, ip) != 0) {
            for (Index is1 = 0; is1 < ns; is1++) {
              for (Index is2 = 0; is2 < ns; is2++) {
                P(is1, is2) += ppvar_pnd(i, ip) * Pe(i, ip, iv, is1, is2);
              }
            }
          }
        }

        // Combine iy0, double transmission and scattering matrix
        Vector iy1(ns), iy2(ns);
        const Index iout = iv * np + ip;
        //mult( iy1, tr_rev(iv,joker,joker), iy0(iv,joker) );
        mult(iy1, ppvar_trans_cumulat(ip, iv, joker, joker), iy0(iv, joker));
        mult(iy2, P, iy1);
        mult(iy(iout, joker), ppvar_trans_cumulat(ip, iv, joker, joker), iy2);

        if (auxBackScat >= 0) {
          mult(iy_aux[auxBackScat](iout, joker), P, iy0(iv, joker));
        }
        if (auxPartAtte >= 0 && ip > 0) {
          iy_aux[auxPartAtte](iout, 0) =
              iy_aux[auxPartAtte](iout - 1, 0) +
              ppath.lstep[ip - 1] *
                  (scalar_ext(ip - 1, iv) + scalar_ext(ip - 1, iv));
        }

        // Jacobians
        if (j_analytical_do) {
          // All this is to cover impact on back-scattering of
          // scattering species
          for (Index iq = 0; iq < nq; iq++) {
            if (jacobian_quantities[iq].Analytical()) {
              if (jac_scat_i[iq] >= 0) {
                // Change of scattering scattering matrix
                P = 0.0;
                for (Index i = 0; i < ne; i++) {
                  if (ppvar_dpnd_dx[iq](i, ip) != 0) {
                    for (Index is1 = 0; is1 < ns; is1++) {
                      for (Index is2 = 0; is2 < ns; is2++) {
                        P(is1, is2) +=
                            ppvar_dpnd_dx[iq](i, ip) * Pe(i, ip, iv, is1, is2);
                      }
                    }
                  }
                }

                // Apply transmissions as above
                Vector iy_tmp(ns);
                mult(iy_tmp, P, iy1);
                mult(diy_dpath[iq](ip, iout, joker),
                     ppvar_trans_cumulat(ip, iv, joker, joker),
                     iy_tmp);
              }
            }
          }
        }  // j_analytical_do
      }    // for iv
    }      // if cloudy
    else {
      if (auxPartAtte >= 0 && ip > 0) {
        for (Index iv = 0; iv < nf; iv++) {
          const Index iout = iv * np + ip;
          iy_aux[auxPartAtte](iout, 0) = iy_aux[auxPartAtte](iout - nf, 0);
        }
      }
    }

    if (auxOptDepth >= 0) {
      for (Index iv = 0; iv < nf; iv++) {
        const Index iout = iv * np + ip;
        iy_aux[auxOptDepth](iout, 0) =
            -2 * log(ppvar_trans_cumulat(ip, iv, 0, 0));
      }
    }

    // Update tr_rev
    /*
      if( ip<np-1 )
        {
          for( Index iv=0; iv<nf; iv++ )
            {          
              Matrix tmp = tr_rev(iv,joker,joker);
              mult( tr_rev(iv,joker,joker), trans_partial(ip+1,iv,joker,joker),
                    tmp );
            }
        }
      */
  }  // for ip

  // Finalize analytical Jacobian
  if (j_analytical_do) {
    // Add impact on transmission, if selected to be included
    if (trans_in_jacobian) {
      Vector jterm(stokes_dim);
      FOR_ANALYTICAL_JACOBIANS_DO(for (Index iv = 0; iv < nf; iv++) {
        for (Index ip = 0; ip < np; ip++)  // Index of x-point
        {
          for (Index j = ip; j < np; j++)  // Index of back-scattering
          {
            // Nothing to do for j==0
            if (j == 0) {
              continue;
            }

            const Index iout = iv * np + j;
            // Impact on transmission towards the sensor
            if (ip > 0) {
              mult(jterm,
                   dtrans_partial_dx_above(ip, iq, iv, joker, joker),
                   iy(iout, joker));
              Matrix tr_inv(ns, ns);
              inv(tr_inv, trans_partial(ip, iv, joker, joker));
              Vector jnew(ns);
              mult(jnew, tr_inv, jterm);
              jnew *= 2;
              diy_dpath[iq](ip, iout, joker) += jnew;
            }
            // Impact on transmission away from the sensor
            if (j > ip && j < np - 1) {
              mult(jterm,
                   dtrans_partial_dx_below(ip, iq, iv, joker, joker),
                   iy(iout, joker));
              Matrix tr_inv(ns, ns);
              inv(tr_inv, trans_partial(ip + 1, iv, joker, joker));
              Vector jnew(ns);
              mult(jnew, tr_inv, jterm);
              jnew *= 2;
              diy_dpath[iq](ip, iout, joker) += jnew;
            }
          }
        }
      })
    }

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
                                    iy_transmission,
                                    water_p_eq_agenda,
                                    jacobian_quantities,
                                    jac_species_i,
                                    jac_is_t);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void iyActiveSingleScat2(Workspace& ws,
                         Matrix& iy,
                         ArrayOfMatrix& iy_aux,
                         ArrayOfTensor3& diy_dx,
                         Vector& ppvar_p,
                         Vector& ppvar_t,
                         Matrix& ppvar_nlte,
                         Matrix& ppvar_vmr,
                         Matrix& ppvar_wind,
                         Matrix& ppvar_mag,
                         Matrix& ppvar_pnd,
                         Matrix& ppvar_f,
                         Tensor4& ppvar_trans_cumulat,
                         const Index& stokes_dim,
                         const Vector& f_grid,
                         const Index& atmosphere_dim,
                         const Vector& p_grid,
                         const Tensor3& t_field,
                         const Tensor4& nlte_field,
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
                         const Index& iy_agenda_call1,
                         const Tensor3& iy_transmission,
                         const Numeric& rte_alonglos_v,
                         const Index& trans_in_jacobian,
                         const Numeric& pext_scaling,
                         const Index& t_interp_order,
                         const Verbosity& verbosity _U_) {
  // Some basic sizes
  const Index nf = f_grid.nelem();
  const Index ns = stokes_dim;
  const Index np = ppath.np;
  const Index ne = pnd_field.nbooks();

  // Radiative background index
  const Index rbi = ppath_what_background(ppath);

  // Checks of input
  // Throw error if unsupported features are requested
  if (!iy_agenda_call1)
    throw runtime_error(
        "Recursive usage not possible (iy_agenda_call1 must be 1)");
  if (!iy_transmission.empty())
    throw runtime_error("*iy_transmission* must be empty");
  if (rbi < 1 || rbi > 9)
    throw runtime_error(
        "ppath.background is invalid. Check your "
        "calculation of *ppath*?");
  if (rbi == 3 || rbi == 4)
    throw runtime_error(
        "The propagation path ends inside or at boundary of "
        "the cloudbox.\nFor this method, *ppath* must be "
        "calculated in this way:\n   ppathCalc( cloudbox_on = 0 ).");
  if (cloudbox_on) {
    if (scat_data_checked != 1)
      throw runtime_error(
          "The scattering data must be flagged to have\n"
          "passed a consistency check (scat_data_checked=1).");
    if (ne != TotalNumberOfElements(scat_data))
      throw runtime_error(
          "*pnd_field* and *scat_data* inconsistent regarding total number of"
          " scattering elements.");
  }
  if (jacobian_do)
    if (dpnd_field_dx.nelem() != jacobian_quantities.nelem())
      throw runtime_error(
          "*dpnd_field_dx* not properly initialized:\n"
          "Number of elements in dpnd_field_dx must be equal number of jacobian"
          " quantities.\n(Note: jacobians have to be defined BEFORE *pnd_field*"
          " is calculated/set.");
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

  if (iy0.ncols() != ns || iy0.nrows() != nf) {
    ostringstream os;
    os << "The size of *iy* returned from *iy_transmitter_agenda* is\n"
       << "not correct:\n"
       << "  expected size = [" << nf << "," << stokes_dim << "]\n"
       << "  size of iy    = [" << iy0.nrows() << "," << iy0.ncols() << "]\n";
    throw runtime_error(os.str());
  }
  for (Index iv = 0; iv < nf; iv++)
    if (iy0(iv, 0) != 1)
      throw runtime_error(
          "The *iy* returned from *iy_transmitter_agenda* "
          "must have the value 1 in the first column.");

  //  Init Jacobian quantities
  Index j_analytical_do = 0;
  if (jacobian_do) FOR_ANALYTICAL_JACOBIANS_DO(j_analytical_do = 1;);
  //
  const Index nq = j_analytical_do ? jacobian_quantities.nelem() : 0;
  ArrayOfTensor3 diy_dpath(nq);
  ArrayOfIndex jac_species_i(nq), jac_scat_i(nq), jac_is_t(nq), jac_wind_i(nq);
  ArrayOfIndex jac_mag_i(nq), jac_other(nq);

  if (j_analytical_do)
    rtmethods_jacobian_init(jac_species_i,
                            jac_scat_i,
                            jac_is_t,
                            jac_wind_i,
                            jac_mag_i,
                            jac_other,
                            diy_dx,
                            diy_dpath,
                            ns,
                            nf,
                            np,
                            nq,
                            abs_species,
                            cloudbox_on,
                            scat_species,
                            dpnd_field_dx,
                            jacobian_quantities,
                            iy_agenda_call1,
                            true);

  // Init iy_aux
  const Index naux = iy_aux_vars.nelem();
  iy_aux.resize(naux);

  // Not implemented in this function yet
  //  Index auxBackScat = -1;
  //  Index auxOptDepth = -1;
  //  Index auxPartAtte = -1;
  //
  for (Index i = 0; i < naux; i++) {
    iy_aux[i].resize(nf * np, ns);
    iy_aux[i] = 0;

    if (iy_aux_vars[i] == "Radiative background")
      iy_aux[i](joker, 0) = (Numeric)min((Index)2, rbi - 1);
    //    else if( iy_aux_vars[i] == "Backscattering" ) {
    //      iy_aux[i]   = 0;
    //      auxBackScat = i;
    //    }
    //    else if( iy_aux_vars[i] == "Optical depth" )
    //      auxOptDepth = i;
    //    else if( iy_aux_vars[i] == "Particle extinction" ) {
    //      iy_aux[i](Range(0,nf,np),joker) = 0;
    //      auxPartAtte = i;
    //    }
    else {
      ostringstream os;
      os << "The only allowed strings in *iy_aux_vars* are:\n"
         << "  \"Radiative background\"\n"
         //      << "  \"Backscattering\"\n"
         //      << "  \"Optical depth\"\n"
         //      << "  \"Particle extinction\"\n"
         << "but you have selected: \"" << iy_aux_vars[i] << "\"";
      throw runtime_error(os.str());
    }
  }

  // Get atmospheric and radiative variables along the propagation path
  ppvar_trans_cumulat.resize(np, nf, ns, ns);

  ArrayOfRadiationVector lvl_rad(np, RadiationVector(nf, ns));
  ArrayOfArrayOfRadiationVector dlvl_rad(
      np, ArrayOfRadiationVector(nq, RadiationVector(nf, ns)));

  ArrayOfTransmissionMatrix lyr_tra(np, TransmissionMatrix(nf, ns));
  ArrayOfArrayOfTransmissionMatrix dlyr_tra_above(
      np, ArrayOfTransmissionMatrix(nq, TransmissionMatrix(nf, ns)));
  ArrayOfArrayOfTransmissionMatrix dlyr_tra_below(
      np, ArrayOfTransmissionMatrix(nq, TransmissionMatrix(nf, ns)));

  // TEMPORARY VARIABLE, THIS COULD BE AN ArrayOfArrayOfPropagationMatrix most likely, to speedup calculations
  Tensor5 Pe(ne, np, nf, ns, ns, 0);

  ArrayOfMatrix ppvar_dpnd_dx(0);
  ArrayOfIndex clear2cloudy;
  //  Matrix scalar_ext(np,nf,0);  // Only used for iy_aux
  Index nf_ssd = scat_data[0][0].pha_mat_data.nlibraries();
  Index duplicate_freqs = ((nf == nf_ssd) ? 0 : 1);
  Tensor6 pha_mat_1se(nf_ssd, 1, 1, 1, ns, ns);
  Vector t_ok(1), t_array(1);
  Matrix pdir(1, 2), idir(1, 2);
  Index ptype;

  if (np == 1 && rbi == 1) {  // i.e. ppath is totally outside the atmosphere:
    ppvar_p.resize(0);
    ppvar_t.resize(0);
    ppvar_vmr.resize(0, 0);
    ppvar_nlte.resize(0, 0);
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
      for (Index ip = 0; ip < np; ip++) clear2cloudy[ip] = -1;
    }

    // Size radiative variables always used
    PropagationMatrix K_this(nf, ns), K_past(nf, ns), Kp(nf, ns);
    StokesVector a(nf, ns), S(nf, ns), Sp(nf, ns);
    ArrayOfIndex lte(np);

    // Init variables only used if transmission part of jacobian
    Vector dB_dT(0);
    ArrayOfPropagationMatrix dK_this_dx(0), dK_past_dx(0), dKp_dx(0);
    ArrayOfStokesVector da_dx(0), dS_dx(0), dSp_dx(0);

    // HSE variables
    Index temperature_derivative_position = -1;
    bool do_hse = false;

    if (trans_in_jacobian && j_analytical_do) {
      dK_this_dx.resize(nq);
      dK_past_dx.resize(nq);
      dKp_dx.resize(nq);
      da_dx.resize(nq);
      dS_dx.resize(nq);
      dSp_dx.resize(nq);
      dB_dT.resize(nf);
      FOR_ANALYTICAL_JACOBIANS_DO(dK_this_dx[iq] = PropagationMatrix(nf, ns);
                                  dK_past_dx[iq] = PropagationMatrix(nf, ns);
                                  dKp_dx[iq] = PropagationMatrix(nf, ns);
                                  da_dx[iq] = StokesVector(nf, ns);
                                  dS_dx[iq] = StokesVector(nf, ns);
                                  dSp_dx[iq] = StokesVector(nf, ns);
                                  if (jacobian_quantities[iq].IsTemperature()) {
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
                                    lte[ip],
                                    dK_this_dx,
                                    dS_dx,
                                    propmat_clearsky_agenda,
                                    jacobian_quantities,
                                    ppvar_f(joker, ip),
                                    ppvar_mag(joker, ip),
                                    ppath.los(ip, joker),
                                    ppvar_nlte(joker, ip),
                                    ppvar_vmr(joker, ip),
                                    ppvar_t[ip],
                                    ppvar_p[ip],
                                    jac_species_i,
                                    trans_in_jacobian && j_analytical_do);

      if (trans_in_jacobian && j_analytical_do)
        adapt_stepwise_partial_derivatives(
            dK_this_dx,
            dS_dx,
            jacobian_quantities,
            ppvar_f(joker, ip),
            ppath.los(ip, joker),
            ppvar_vmr(joker, ip),
            ppvar_t[ip],
            ppvar_p[ip],
            jac_species_i,
            jac_wind_i,
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

        K_this += Kp;

        //        if( auxPartAtte >= 0 )
        //          scalar_ext(ip,joker) = Kp.Kjj();

        if (trans_in_jacobian && j_analytical_do)
          FOR_ANALYTICAL_JACOBIANS_DO(dK_this_dx[iq] += dKp_dx[iq];);

        // Get back-scattering per particle, where relevant
        {
          // Direction of outgoing scattered radiation (which is reverse to LOS).
          Vector los_sca;
          mirror_los(los_sca, ppath.los(ip, joker), atmosphere_dim);
          pdir(0, joker) = los_sca;

          // Obtain a length-2 vector for incoming direction
          Vector los_inc;
          if (atmosphere_dim == 3) {
            los_inc = ppath.los(ip, joker);
          } else  // Mirror back to get a correct 3D LOS
          {
            mirror_los(los_inc, los_sca, 3);
          }
          idir(0, joker) = los_inc;

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
                                  pdir,
                                  idir,
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
                  ostringstream os;
                  os << "Interpolation error for (flat-array) scattering"
                     << " element #" << i_se_flat << "\n"
                     << "at location/temperature point #" << ip << "\n";
                  throw runtime_error(os.str());
                }
              }
              i_se_flat++;
            }

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
      cumulative_transmission(lyr_tra, CumulativeTransmission::ForwardReverse);
  const ArrayOfTransmissionMatrix tot_tra_reflect =
      cumulative_transmission(lyr_tra, CumulativeTransmission::Reflect);
  const ArrayOfTransmissionMatrix reflect_matrix =
      cumulative_backscatter(Pe, ppvar_pnd);
  const ArrayOfArrayOfTransmissionMatrix dreflect_matrix =
      cumulative_backscatter_derivative(Pe, ppvar_dpnd_dx);

  lvl_rad[0] = iy0;
  set_backscatter_radiation_vector(lvl_rad,
                                   dlvl_rad,
                                   lyr_tra,
                                   tot_tra_forward,
                                   tot_tra_reflect,
                                   reflect_matrix,
                                   dlyr_tra_above,
                                   dlyr_tra_below,
                                   dreflect_matrix,
                                   BackscatterSolver::Full);

  // Size iy and set to zero
  iy.resize(nf * np, ns);  // iv*np + ip is the desired output order...
  for (Index ip = 0; ip < np; ip++) {
    for (Index iv = 0; iv < nf; iv++) {
      for (Index is = 0; is < stokes_dim; is++) {
        iy(iv * np + ip, is) = lvl_rad[ip](iv, is);
        FOR_ANALYTICAL_JACOBIANS_DO(diy_dpath[iq](ip, iv * np + ip, is) =
                                        dlvl_rad[ip][iq](iv, is););
      }
    }
  }

  // FIXME: Add the aux-variables back
  // FIXME: Add the diy_dpath

  // Finalize analytical Jacobian
  if (j_analytical_do and false)  // FIXME: remove false
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
                                    iy_transmission,
                                    water_p_eq_agenda,
                                    jacobian_quantities,
                                    jac_species_i,
                                    jac_is_t);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void yActive(Workspace& ws,
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
             const String& iy_unit,
             const ArrayOfString& iy_aux_vars,
             const Index& stokes_dim,
             const Vector& f_grid,
             const Index& atmosphere_dim,
             const Tensor3& t_field,
             const Tensor3& z_field,
             const Tensor4& vmr_field,
             const Tensor4& nlte_field,
             const Index& cloudbox_on,
             const Index& cloudbox_checked,
             const Matrix& sensor_pos,
             const Matrix& sensor_los,
             const Index& sensor_checked,
             const Index& jacobian_do,
             const ArrayOfRetrievalQuantity& jacobian_quantities,
             const Agenda& iy_main_agenda,
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
  if (f_grid.empty()) {
    throw runtime_error("The frequency grid is empty.");
  }
  chk_if_increasing("f_grid", f_grid);
  if (f_grid[0] <= 0) {
    throw runtime_error("All frequencies in *f_grid* must be > 0.");
  }
  //
  chk_if_in_range("stokes_dim", stokes_dim, 1, 4);
  if (atmfields_checked != 1)
    throw runtime_error(
        "The atmospheric fields must be flagged to have "
        "passed a consistency check (atmfields_checked=1).");
  if (atmgeom_checked != 1)
    throw runtime_error(
        "The atmospheric geometry must be flagged to have "
        "passed a consistency check (atmgeom_checked=1).");
  if (cloudbox_checked != 1)
    throw runtime_error(
        "The cloudbox must be flagged to have "
        "passed a consistency check (cloudbox_checked=1).");
  if (sensor_checked != 1)
    throw runtime_error(
        "The sensor variables must be flagged to have "
        "passed a consistency check (sensor_checked=1).");

  // Method specific variables
  bool is_z = max(range_bins) > 1;
  if (!is_increasing(range_bins))
    throw runtime_error(
        "The vector *range_bins* must contain strictly "
        "increasing values.");
  if (!is_z && min(range_bins) < 0)
    throw runtime_error(
        "The vector *range_bins* is not allowed to contain "
        "negative times.");
  if (instrument_pol_array.nelem() != nf)
    throw runtime_error(
        "The main length of *instrument_pol_array* must match "
        "the number of frequencies.");

  // iy_unit and variables to handle conversion to Ze and dBZe
  Vector cfac(nf, 1.0);
  Numeric ze_min = 0;
  const Numeric jfac = 10 * log10(exp(1.0));
  if (iy_unit == "1") {
  } else if (iy_unit == "Ze") {
    ze_cfac(cfac, f_grid, ze_tref, k2);
  } else if (iy_unit == "dBZe") {
    ze_cfac(cfac, f_grid, ze_tref, k2);
    ze_min = pow(10.0, dbze_min / 10);
  } else {
    throw runtime_error(
        "For this method, *iy_unit* must be set to \"1\", "
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
    Tensor3 iy_transmission(0, 0, 0);
    ArrayOfTensor3 diy_dx;
    Vector rte_pos2(0);
    Matrix iy;
    Ppath ppath;
    ArrayOfMatrix iy_aux;
    const Index iy_id = (Index)1e6 * p;
    //
    iy_main_agendaExecute(ws,
                          iy,
                          iy_aux,
                          ppath,
                          diy_dx,
                          1,
                          iy_unit,
                          iy_transmission,
                          iy_aux_vars,
                          iy_id,
                          cloudbox_on,
                          jacobian_do,
                          t_field,
                          z_field,
                          vmr_field,
                          nlte_field,
                          f_grid,
                          sensor_pos(p, joker),
                          sensor_los(p, joker),
                          rte_pos2,
                          iy_main_agenda);

    // Check if path and size OK
    const Index np = ppath.np;
    if (np == 1)
      throw runtime_error(
          "A path consisting of a single point found. "
          "This is not allowed.");
    error_if_limb_ppath(ppath);
    if (iy.nrows() != nf * np)
      throw runtime_error(
          "The size of *iy* returned from *iy_main_agenda* "
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
                    if (iy_unit == "dBZe") {
                      jacobian(iout, jacobian_indices[iq][0] + k) *=
                          jfac / max(y[iout], ze_min);
                    }
                  })
            }

            if (iy_unit == "dBZe") {
              y[iout] = y[iout] <= ze_min ? dbze_min : 10 * log10(y[iout]);
            }

            // Same for aux variables
            for (Index a = 0; a < naux; a++) {
              if (iy_aux_vars[a] == "Backscattering") {
                y_aux[a][iout] = cfac[iv] * (hbin * auxvar[a]);
                if (iy_unit == "dBZe") {
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
    if (geo_pos.nelem() && geo_pos.nelem() != atmosphere_dim) {
      throw runtime_error(
          "Wrong size of *geo_pos* obtained from "
          "*geo_pos_agenda*.\nThe length of *geo_pos* must "
          "be zero or equal to *atmosphere_dim*.");
    }
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
