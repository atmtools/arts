/* Copyright (C) 2015
   Richard Larsson
                            
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
 * @file m_radiation_field.cc
 * @author Richard Larsson
 * @date 2015-09-13
 * 
 * @brief Radiation field calculations for the user
 */

#include "absorption.h"
#include "arts.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "linefunctions.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "propmat_field.h"
#include "radiation_field.h"

void line_irradianceCalcForSingleSpeciesNonOverlappingLinesPseudo2D(
    Workspace& ws,
    Matrix& line_irradiance,
    Tensor3& line_transmission,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfArrayOfLineRecord& abs_lines_per_species,
    const Tensor4& nlte_field,
    const Tensor4& vmr_field,
    const Tensor3& t_field,
    const Tensor3& z_field,
    const Vector& p_grid,
    const Vector& refellipsoid,
    const Tensor3& surface_props_data,
    const Agenda& ppath_agenda,
    const Agenda& iy_main_agenda,
    const Agenda& iy_space_agenda,
    const Agenda& iy_surface_agenda,
    const Agenda& iy_cloudbox_agenda,
    const Agenda& propmat_clearsky_agenda,
    const Numeric& df,
    const Index& nz,
    const Index& nf,
    const Numeric& r,
    const Verbosity& verbosity) {
  if (abs_lines_per_species.nelem() not_eq 1)
    throw std::runtime_error("Only for one species...");
  if (nf % 2 not_eq 1)
    throw std::runtime_error("Must hit line center, nf % 2 must be 1.");
  const Index nl = abs_lines_per_species[0].nelem();
  const Index np = p_grid.nelem();

  // Compute variables
  ArrayOfTensor3 diy_dx;
  FieldOfPropagationMatrix propmat_field;
  FieldOfStokesVector absorption_field, additional_source_field;
  ArrayOfPpath ppath_field;

  // Check that the lines and nf is correct
  Vector f_grid(nf * nl);
  for (Index il = 0; il < nl; il++) {
    const Numeric d = abs_lines_per_species[0][il].F() * df;
    nlinspace(f_grid[Range(il * nf, nf)],
              abs_lines_per_species[0][il].F() - d,
              abs_lines_per_species[0][il].F() + d,
              nf);
  }

  ppath_fieldFromDownUpLimbGeoms(ws,
                                 ppath_field,
                                 ppath_agenda,
                                 -1,
                                 1e99,
                                 1,
                                 z_field,
                                 f_grid,
                                 0,
                                 1,
                                 0,
                                 Vector(1, 0),
                                 Vector(1, 0),
                                 Vector(0),
                                 refellipsoid,
                                 1,
                                 nz,
                                 verbosity);

  ArrayOfArrayOfIndex sorted_index;
  ArrayOfVector cos_zenith_angles;
  sorted_index_of_ppath_field(sorted_index, cos_zenith_angles, ppath_field);

  for (Index ip = 0; ip < np; ip++)
    error_in_integrate(
        "Your lineshape integration does normalize.  Increase nf and decrease df until it does.",
        test_integrate_zenith(cos_zenith_angles[ip], sorted_index[ip]));

  field_of_propagation(ws,
                       propmat_field,
                       absorption_field,
                       additional_source_field,
                       1,
                       f_grid,
                       p_grid,
                       z_field,
                       t_field,
                       nlte_field,
                       vmr_field,
                       ArrayOfRetrievalQuantity(0),
                       propmat_clearsky_agenda);

  Array<Array<Eigen::VectorXcd>> lineshapes(
      nl, Array<Eigen::VectorXcd>(np, Eigen::VectorXcd(nf * nl)));
  for (Index il = 0; il < nl; il++)
    for (Index ip = 0; ip < np; ip++)
      Linefunctions::set_lineshape(lineshapes[il][ip],
                                   MapToEigen(f_grid),
                                   abs_lines_per_species[0][il],
                                   Vector(1, vmr_field(0, ip, 0, 0)),
                                   t_field(ip, 0, 0),
                                   p_grid[ip],
                                   0.0,
                                   0.0,
                                   abs_species);
  for (auto& aols : lineshapes)
    for (auto& ls : aols)
      error_in_integrate(
          "Your lineshape integration does normalize.  Increase nf and decrease df until it does.",
          test_integrate_convolved(ls, f_grid));

  // Counting the path index so we can make the big loop parallel
  ArrayOfArrayOfIndex counted_path_index(ppath_field.nelem());
  ArrayOfIndex counter(np, 0);
  for (Index i = 0; i < ppath_field.nelem(); i++) {
    const Ppath& path = ppath_field[i];
    counted_path_index[i].resize(path.np);
    for (Index ip_path = 0; ip_path < path.np; ip_path++) {
      const Index ip_grid = grid_index_from_gp(path.gp_p[ip_path]);
      counted_path_index[i][ip_path] = counter[ip_grid];
      counter[ip_grid]++;
    }
  }

  line_irradiance = Matrix(nl, np, 0.0);
  line_transmission = Tensor3(1, nl, np, 0.0);

  ArrayOfMatrix line_radiance(np);
  for (Index i = 0; i < np; i++)
    line_radiance[i].resize(sorted_index[i].nelem(), nl);

  Workspace l_ws(ws);
  Agenda l_iy_main_agenda(iy_main_agenda);
  Agenda l_iy_space_agenda(iy_space_agenda);
  Agenda l_iy_surface_agenda(iy_surface_agenda);
  Agenda l_iy_cloudbox_agenda(iy_cloudbox_agenda);

#pragma omp parallel for if (not arts_omp_in_parallel())               \
    schedule(guided) default(shared) firstprivate(l_ws,                \
                                                  l_iy_main_agenda,    \
                                                  l_iy_space_agenda,   \
                                                  l_iy_surface_agenda, \
                                                  l_iy_cloudbox_agenda)
  for (Index i = 0; i < ppath_field.nelem(); i++) {
    const Ppath& path = ppath_field[i];

    thread_local ArrayOfRadiationVector lvl_rad;
    thread_local ArrayOfRadiationVector src_rad;
    thread_local ArrayOfTransmissionMatrix lyr_tra;
    thread_local ArrayOfTransmissionMatrix tot_tra;

    emission_from_propmat_field(l_ws,
                                lvl_rad,
                                src_rad,
                                lyr_tra,
                                tot_tra,
                                propmat_field,
                                absorption_field,
                                additional_source_field,
                                f_grid,
                                t_field,
                                path,
                                l_iy_main_agenda,
                                l_iy_space_agenda,
                                l_iy_surface_agenda,
                                l_iy_cloudbox_agenda,
                                surface_props_data,
                                verbosity);

    for (Index ip_path = 0; ip_path < path.np; ip_path++) {
      const Index ip_grid = grid_index_from_gp(path.gp_p[ip_path]);
      for (Index il = 0; il < nl; il++)
        line_radiance[ip_grid](counted_path_index[i][ip_path], il) =
            integrate_convolved(
                lvl_rad[ip_path], lineshapes[il][ip_grid], f_grid);
    }
  }

  for (Index ip = 0; ip < np; ip++)
    for (Index il = 0; il < nl; il++)
      line_irradiance(il, ip) = integrate_zenith(line_radiance[ip](joker, il),
                                                 cos_zenith_angles[ip],
                                                 sorted_index[ip]);

  if (r > 0) {
    const FieldOfTransmissionMatrix transmat_field =
        transmat_field_calc_from_propmat_field(propmat_field, r);
    for (Index ip = 0; ip < np; ip++)
      for (Index il = 0; il < nl; il++)
        line_transmission(0, il, ip) = integrate_convolved(
            transmat_field(ip, 0, 0), lineshapes[il][ip], f_grid);
  }
}

void line_irradianceCalcForSingleSpeciesNonOverlappingLines(
    Workspace& ws,
    Matrix& line_irradiance,
    Tensor3& line_transmission,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfArrayOfLineRecord& abs_lines_per_species,
    const Tensor4& nlte_field,
    const Tensor4& vmr_field,
    const Tensor3& t_field,
    const Tensor3& z_field,
    const Vector& p_grid,
    const Index& atmosphere_dim,
    const Tensor3& surface_props_data,
    const Agenda& iy_space_agenda,
    const Agenda& iy_surface_agenda,
    const Agenda& iy_cloudbox_agenda,
    const Agenda& propmat_clearsky_agenda,
    const Agenda& water_p_eq_agenda,
    const Numeric& df,
    const Index& nz,
    const Index& nf,
    const Verbosity& verbosity) {
  extern const Numeric DEG2RAD;
  const Index ns = abs_species.nelem();
  if (ns not_eq 1)
    throw std::runtime_error("Only for a single species in this test version");
  const Index nl = abs_lines_per_species[0].nelem();
  const Index np = p_grid.nelem();
  if (nz % 2)
    throw std::runtime_error(
        "Cannot hit the tangent point in this test version; nz must be even");
  if (nf % 2 not_eq 1)
    throw std::runtime_error(
        "Must include central frequency to test validity; nf must be uneven.");

  // Find Zenith angles and weigh them by their area
  Vector za_grid;
  nlinspace(za_grid, 0.0, 180.0, nz);
  Vector wzad(nz - 1);
  for (Index iz = 0; iz < nz - 1; iz++)
    wzad[iz] = cos(DEG2RAD * za_grid[iz]) - cos(DEG2RAD * za_grid[iz + 1]);

  line_irradiance = Matrix(nl, np, 0.0);
  line_transmission = Tensor3(1, nl, np, 0.0);

  Tensor7 doit_i_field;
  Tensor3 trans_field;
  for (Index il = 0; il < nl; il++) {
    const Numeric d = abs_lines_per_species[0][il].F() * df;
    Vector f_grid;
    nlinspace(f_grid,
              abs_lines_per_species[0][il].F() - d,
              abs_lines_per_species[0][il].F() + d,
              nf);
    doit_i_fieldClearskyPlaneParallel(ws,
                                      doit_i_field,
                                      trans_field,
                                      propmat_clearsky_agenda,
                                      water_p_eq_agenda,
                                      iy_space_agenda,
                                      iy_surface_agenda,
                                      iy_cloudbox_agenda,
                                      1,
                                      f_grid,
                                      atmosphere_dim,
                                      p_grid,
                                      z_field,
                                      t_field,
                                      nlte_field,
                                      vmr_field,
                                      abs_species,
                                      Tensor3(0, 0, 0),
                                      Tensor3(0, 0, 0),
                                      Tensor3(0, 0, 0),
                                      Tensor3(0, 0, 0),
                                      Tensor3(0, 0, 0),
                                      Tensor3(0, 0, 0),
                                      z_field(0, joker, joker),
                                      1e99,
                                      0.0,
                                      surface_props_data,
                                      za_grid,
                                      0,
                                      verbosity);

    // Integrate over the sphere
    for (Index ip = 0; ip < np; ip++) {
      Eigen::VectorXcd F(nf);
      Linefunctions::set_lineshape(F,
                                   MapToEigen(f_grid),
                                   abs_lines_per_species[0][il],
                                   Vector(1, vmr_field(0, ip, 0, 0)),
                                   t_field(ip, 0, 0),
                                   p_grid[ip],
                                   0.0,
                                   0.0,
                                   abs_species);
      Numeric sx = 0;
      for (Index iv = 0; iv < nf - 1; iv++) {
        const Numeric intF =
            (F[iv].real() + F[iv + 1].real()) * (f_grid[iv + 1] - f_grid[iv]);
        for (Index iz = 0; iz < nz - 1; iz++) {
          const Numeric x = wzad[iz] * 0.25 * intF;
          sx += x;
          line_irradiance(il, ip) +=
              (doit_i_field(iv, ip, 0, 0, iz, 0, 0) +
               doit_i_field(iv, ip, 0, 0, iz + 1, 0, 0) +
               doit_i_field(iv + 1, ip, 0, 0, iz, 0, 0) +
               doit_i_field(iv + 1, ip, 0, 0, iz + 1, 0, 0)) *
              0.25 * x;

          const Numeric exp_mtau =
              0.25 *
              (trans_field(iv, ip, iz) + trans_field(iv, ip, iz + 1) +
               trans_field(iv + 1, ip, iz) + trans_field(iv + 1, ip, iz + 1));
          if (abs(1 - exp_mtau) < 1e-6) { /* do nothing */
          } else
            line_transmission(0, il, ip) +=
                (1 + (1 - exp_mtau) / log(exp_mtau)) * x;
        }
      }

      if (abs(sx - 1) > 1e-3) {
        ostringstream os;
        os << "Integrated iy normalizes to " << sx << " instead of 1\n";
        os << "This means your frequency grid spanning " << f_grid[0] << " to "
           << f_grid[nf - 1] << " is not good enough\n";
        if (sx < 1)
          os << "Please consider increasing df by about the inverse of the missing normalizations (a factor about "
             << 1 / sx << " is approximately enough)\n";
        else
          os << "Please consider increasing nf so that the frequency grid better represents the lineshape\n";
        throw std::runtime_error(os.str());
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void doit_i_fieldClearskyPlaneParallel(
    Workspace& ws,
    Tensor7& doit_i_field,
    Tensor3& trans_field,
    const Agenda& propmat_clearsky_agenda,
    const Agenda& water_p_eq_agenda,
    const Agenda& iy_space_agenda,
    const Agenda& iy_surface_agenda,
    const Agenda& iy_cloudbox_agenda,
    const Index& stokes_dim,
    const Vector& f_grid,
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Tensor3& z_field,
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
    const Matrix& z_surface,
    const Numeric& ppath_lmax,
    const Numeric& rte_alonglos_v,
    const Tensor3& surface_props_data,
    const Vector& scat_za_grid,
    const Index& use_parallel_iy,
    const Verbosity& verbosity) {
  // Check input
  if (atmosphere_dim != 1)
    throw runtime_error("This method only works for atmosphere_dim = 1.");
  chk_if_increasing("scat_za_grid", scat_za_grid);
  if (abs(z_surface(0, 0) - z_field(0, 0, 0)) > 1e-3)
    throw runtime_error(
        "The surface must be placed exactly at the first pressure level.");

  // Sizes
  const Index nl = p_grid.nelem();
  const Index nf = f_grid.nelem();
  const Index nza = scat_za_grid.nelem();

  // Init doit_i_field and trans_field
  doit_i_field.resize(nf, nl, 1, 1, nza, 1, stokes_dim);
  trans_field.resize(nf, nl, nza);

  // De-activate cloudbox
  const Index cloudbox_on = 0, ppath_inside_cloudbox_do = 0;
  const ArrayOfIndex cloudbox_limits(0);

  // Various input variables
  const String iy_unit = "1";
  const ArrayOfString iy_aux_vars(0);
  const Vector rte_pos2(0);
  const Index iy_agenda_call1 = 1;
  const Tensor3 iy_transmission(0, 0, 0);
  const Index jacobian_do = 0;
  const ArrayOfRetrievalQuantity jacobian_quantities(0);
  // Create one altitde just above TOA
  const Numeric z_space = z_field(nl - 1, 0, 0) + 10;

  Workspace l_ws(ws);
  ArrayOfString fail_msg;
  bool failed = false;

  // Define iy_main_agenda to be consistent with the methods applied inside
  // this method. This definition of iy_main_agenda will be used to when
  // calculating the the radiation reflected by the surface
  Agenda iy_main_agenda;
  iy_main_agenda.append("ppathPlaneParallel", TokVal());
  iy_main_agenda.append("iyEmissionStandard", TokVal());
  iy_main_agenda.set_name("iy_main_agenda");
  iy_main_agenda.check(ws, verbosity);

  // Loop zenith angles
  //
  if (nza)
#pragma omp parallel for if (!arts_omp_in_parallel() && nza > 1 && \
                             !use_parallel_iy) firstprivate(l_ws)
    for (Index i = 0; i < nza; i++) {
      if (failed) continue;
      try {
        // Define output variables
        Ppath ppath;
        Vector ppvar_p, ppvar_t;
        Matrix iy, ppvar_nlte, ppvar_vmr, ppvar_wind, ppvar_mag, ppvar_f;
        Tensor3 ppvar_iy;
        Tensor4 ppvar_trans_cumulat, ppvar_trans_partial;
        ArrayOfMatrix iy_aux;
        ArrayOfTensor3 diy_dx;

        Index iy_id = i;
        Vector rte_los(1, scat_za_grid[i]);
        Vector rte_pos(1, scat_za_grid[i] < 90 ? z_surface(0, 0) : z_space);

        ppathPlaneParallel(ppath,
                           atmosphere_dim,
                           z_field,
                           z_surface,
                           cloudbox_on,
                           cloudbox_limits,
                           ppath_inside_cloudbox_do,
                           rte_pos,
                           rte_los,
                           ppath_lmax,
                           verbosity);

        if (use_parallel_iy) {
          iyEmissionStandardParallel(l_ws,
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
                                     iy_main_agenda,
                                     iy_space_agenda,
                                     iy_surface_agenda,
                                     iy_cloudbox_agenda,
                                     iy_agenda_call1,
                                     iy_transmission,
                                     rte_alonglos_v,
                                     surface_props_data,
                                     verbosity);
        } else {
          iyEmissionStandard(l_ws,
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
                             iy_main_agenda,
                             iy_space_agenda,
                             iy_surface_agenda,
                             iy_cloudbox_agenda,
                             iy_agenda_call1,
                             iy_transmission,
                             rte_alonglos_v,
                             surface_props_data,
                             verbosity);
        }

        assert(iy.ncols() == stokes_dim);

        // First and last points are most easily handled separately
        if (scat_za_grid[i] < 90) {
          doit_i_field(joker, 0, 0, 0, i, 0, joker) = ppvar_iy(joker, joker, 0);
          doit_i_field(joker, nl - 1, 0, 0, i, 0, joker) =
              ppvar_iy(joker, joker, ppath.np - 1);
          trans_field(joker, 0, i) = ppvar_trans_partial(0, joker, 0, 0);
          trans_field(joker, nl - 1, i) =
              ppvar_trans_partial(ppath.np - 1, joker, 0, 0);
        } else {
          doit_i_field(joker, nl - 1, 0, 0, i, 0, joker) =
              ppvar_iy(joker, joker, 0);
          doit_i_field(joker, 0, 0, 0, i, 0, joker) =
              ppvar_iy(joker, joker, ppath.np - 1);
          trans_field(joker, nl - 1, i) = ppvar_trans_partial(0, joker, 0, 0);
          trans_field(joker, 0, i) =
              ppvar_trans_partial(ppath.np - 1, joker, 0, 0);
        }

        // Remaining points
        for (Index p = 1; p < ppath.np - 1; p++) {
          // We just store values at pressure levels
          if (ppath.gp_p[p].fd[0] < 1e-2) {
            doit_i_field(joker, ppath.gp_p[p].idx, 0, 0, i, 0, joker) =
                ppvar_iy(joker, joker, p);
            trans_field(joker, ppath.gp_p[p].idx, i) =
                ppvar_trans_partial(p, joker, 0, 0);
          }
        }
      } catch (const std::exception& e) {
#pragma omp critical(planep_setabort)
        failed = true;

        ostringstream os;
        os << "Run-time error at nza #" << i << ": \n" << e.what();
#pragma omp critical(planep_push_fail_msg)
        fail_msg.push_back(os.str());
      }
    }

  if (fail_msg.nelem()) {
    ostringstream os;
    for (auto& msg : fail_msg) os << msg << '\n';
    throw runtime_error(os.str());
  }
}
