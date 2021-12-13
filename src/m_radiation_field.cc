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
#include "lineshape.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "propmat_field.h"
#include "radiation_field.h"

void line_irradianceCalcForSingleSpeciesNonOverlappingLinesPseudo2D(
    Workspace& ws,
    Matrix& line_irradiance,
    Tensor3& line_transmission,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const EnergyLevelMap& nlte_field,
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
    const Verbosity& verbosity)
{
  ARTS_USER_ERROR_IF (abs_lines_per_species.nelem() not_eq 1,
                      "Only for one species...");
  ARTS_USER_ERROR_IF (nf % 2 not_eq 1,
                      "Must hit line center, nf % 2 must be 1.");
  const Index nl = nelem(abs_lines_per_species);
  const Index np = p_grid.nelem();

  // Compute variables
  ArrayOfTensor3 diy_dx;
  FieldOfPropagationMatrix propmat_field;
  FieldOfStokesVector absorption_field, additional_source_field;
  ArrayOfPpath ppath_field;

  // Check that the lines and nf is correct
  Vector f_grid(nf * nl);
  Index il=0;
  for (auto& lines: abs_lines_per_species) {
    for (auto& band: lines) {
      for (Index k=0; k<band.NumLines(); k++) {
        nlinspace(f_grid[Range(il * nf, nf)], band.lines[k].F0 * (1 - df), band.lines[k].F0 * (1 + df), nf);
        il++;
      }
    }
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
  for (Index ip=0; ip<np; ip++) {
    il=0;
    for (auto& lines: abs_lines_per_species) {
      for (auto& band: lines) {
        const Numeric DC = band.DopplerConstant(t_field(ip, 0, 0));
        const Vector vmrs = band.BroadeningSpeciesVMR(vmr_field(joker, ip, 0, 0), abs_species);
        for (Index k=0; k<band.NumLines(); k++) {
          const auto X = band.ShapeParameters(k, t_field(ip, 0, 0), p_grid[ip], vmrs);
          LineShape::Calculator ls(band.lineshapetype, band.lines[k].F0, X, DC, 0, false);
          for (Index iv=0; iv<f_grid.nelem(); iv++) {
            lineshapes[il][ip][iv] = ls(f_grid[iv]);
          }
          il++;
        }
      }
    }
  }
  for (auto& aols : lineshapes)
    for (auto& ls : aols)
      error_in_integrate(
          "Your lineshape integration does not normalize.  Increase nf and decrease df until it does.",
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
                                                  l_iy_cloudbox_agenda,\
                                                  il)
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
                                nlte_field,
                                path,
                                l_iy_main_agenda,
                                l_iy_space_agenda,
                                l_iy_surface_agenda,
                                l_iy_cloudbox_agenda,
                                surface_props_data,
                                verbosity);
    
    for (Index ip_path = 0; ip_path < path.np; ip_path++) {
      const Index ip_grid = grid_index_from_gp(path.gp_p[ip_path]);
      for (il = 0; il < nl; il++)
        line_radiance[ip_grid](counted_path_index[i][ip_path], il) =
            integrate_convolved(
                lvl_rad[ip_path], lineshapes[il][ip_grid], f_grid);
    }
  }

  for (Index ip = 0; ip < np; ip++) {
    for (il = 0; il < nl; il++) {
      line_irradiance(il, ip) = integrate_zenith(line_radiance[ip](joker, il),
                                                 cos_zenith_angles[ip],
                                                 sorted_index[ip]);
    }
  }

  if (r > 0) {
    const FieldOfTransmissionMatrix transmat_field =
        transmat_field_calc_from_propmat_field(propmat_field, r);
    for (Index ip = 0; ip < np; ip++)
      for (il = 0; il < nl; il++)
        line_transmission(0, il, ip) = integrate_convolved(
            transmat_field(ip, 0, 0), lineshapes[il][ip], f_grid);
  }
}

