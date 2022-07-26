/* Copyright (C) 2018
   Manfred Brath  <manfred.brath@uni-hamburg.de>

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
#include <iostream>
#include <stdexcept>
#include "absorption.h"
#include "agenda_class.h"
#include "auto_md.h"
#include "check_input.h"
#include "legendre.h"
#include "math_funcs.h"
#include "matpackVII.h"
#include "messages.h"
#include "sorting.h"
#include "surface.h"
#include "workspace_ng.h"
#include "check_input.h"
#include "global_data.h"

/*!
  \file   m_fluxes.cc
  \author Manfred Brath  <manfred.brath@uni-hamburg.de>
  \date   2018-06-22

  \brief  Workspace functions related to simulation of radiation fluxes.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

extern const Numeric PI;
extern const Numeric DEG2RAD;

/*===========================================================================
  === The functions
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void AngularGridsSetFluxCalc(Vector& za_grid,
                             Vector& aa_grid,
                             Vector& za_grid_weights,
                             // Keywords:
                             const Index& N_za_grid,
                             const Index& N_aa_grid,
                             const String& za_grid_type,
                             const Verbosity&) {
  // Azimuth angle grid
  if (N_aa_grid > 1)
    nlinspace(aa_grid, 0, 360, N_aa_grid);
  else if (N_aa_grid < 1) {
    ostringstream os;
    os << "N_aa_grid must be > 0 (even for 1D).";
    throw std::runtime_error(os.str());
  } else {
    aa_grid.resize(1);
    aa_grid[0] = 0.;
  }

  if (N_za_grid % 2 == 1) {
    ostringstream os;
    os << "N_za_grid must be even.";
    throw runtime_error(os.str());
  }

  Index nph = N_za_grid / 2;

  //calculate zenith angle grid
  za_grid.resize(N_za_grid);
  za_grid = 0.;
  za_grid_weights.resize(N_za_grid);
  za_grid_weights = 0;

  if (za_grid_type == "double_gauss") {
    Vector x;
    Vector w;
    Vector xtemp;
    Vector wtemp;
    //Numeric theta;

    //calculate legendre weights and evaluation points
    gsl_integration_glfixed_table_alloc(xtemp, wtemp, nph);

    x.resize(nph);
    w.resize(nph);

    // reorder and expand weights and abscissa vectors
    // transform abscissa vector from cos(theta)-space to theta-space
    // and adjust the domain and weights
    if (nph % 2 == 1) {
      x[xtemp.nelem() - 1] = acos((xtemp[0] + 1) / 2) / DEG2RAD;
      w[wtemp.nelem() - 1] = wtemp[0] / 2;

      for (Index i = 0; i < xtemp.nelem() - 1; i++) {
        x[i] = acos((xtemp[xtemp.nelem() - 1 - i] + 1) / 2.) / DEG2RAD;
        x[xtemp.nelem() + i] = acos(1 - (xtemp[i + 1] + 1) / 2.) / DEG2RAD;

        w[i] = wtemp[wtemp.nelem() - 1 - i] / 2;
        w[wtemp.nelem() + i] = wtemp[i + 1] / 2;
      }
    } else {
      for (Index i = 0; i < xtemp.nelem(); i++) {
        x[i] = acos((xtemp[xtemp.nelem() - 1 - i] + 1) / 2.) / DEG2RAD;
        x[xtemp.nelem() + i] = acos(1 - (xtemp[i] + 1) / 2.) / DEG2RAD;

        w[i] = wtemp[wtemp.nelem() - 1 - i] / 2;
        w[wtemp.nelem() + i] = wtemp[i] / 2;
      }
    }

    for (Index i = 0; i < nph; i++) {
      //set the angles
      //theta=x[i];//acos((x[i]+1)/2)/DEG2RAD;
      za_grid[i] = x[i];
      za_grid[za_grid.nelem() - 1 - i] = 180 - x[i];

      // set the weights to the right component
      za_grid_weights[i] = w[i];
      za_grid_weights[za_grid_weights.nelem() - 1 - i] = w[i];
    }

  } else if (za_grid_type == "linear") {
    Vector x;
    Vector w;
    calculate_weights_linear(x, w, nph);

    for (Index i = 0; i < N_za_grid; i++) {
      za_grid[i] = (x[i] + 1) * 90.;

      // set the weights to the right component
      // by adjusting the domain, we also have to adjust the weights
      za_grid_weights[i] = w[i] * sin(za_grid[i] * DEG2RAD);
    }
  } else if (za_grid_type == "linear_mu") {
    Vector x;
    Vector w;

    //calculate weights in cos(theta) space
    calculate_weights_linear(x, w, nph);

    //allocate
    Vector za_grid_temp;
    za_grid_temp.resize(x.nelem());

    for (Index i = 0; i < N_za_grid; i++) {
      za_grid_temp[i] = acos(x[i]) / DEG2RAD;
    }

    //#sort weights and theta in increasing direction of za_grid
    za_grid = za_grid_temp[Range(x.nelem() - 1, x.nelem(), -1)];
    za_grid_weights = w[Range(x.nelem() - 1, x.nelem(), -1)];

  } else {
    ostringstream os;
    os << "The selected grid type is not implemented";
    throw std::runtime_error(os.str());
  }

  //be sure that the first and the last angle are within the closed interval
  //between 0 and 180 deg, because ARTS is picky if the angles are due to numerics
  // slightly below and above,respectively.
  if (za_grid[0] < 0) {
    za_grid[0] = 0.;
  }

  if (za_grid[za_grid.nelem() - 1] > 180) {
    za_grid[za_grid.nelem() - 1] = 180.;
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void heating_ratesFromIrradiance(Tensor3& heating_rates,
                                 const Vector& p_grid,
                                 const Tensor4& irradiance_field,
                                 const Tensor3& specific_heat_capacity,
                                 const Numeric& g0,
                                 const Verbosity&) {
  //allocate
  heating_rates.resize(irradiance_field.nbooks(),
                       irradiance_field.npages(),
                       irradiance_field.nrows());
  heating_rates = 0;

  // allocate some auxiliary variables
  Numeric net_flux_b;  //net_flux bottom
  Numeric net_flux_c;  //net_flux center
  Numeric net_flux_t;  //net_flux top
  Index idx;

  // calculate heating rates, we skip the upper and lower boundary here, because to achieve the same
  // second order accuracy for the derivation of the net flux at the edged, we use
  // a differentiation based on polynomial interpolation
  for (Index b = 1; b < irradiance_field.nbooks() - 1; b++) {
    for (Index p = 0; p < irradiance_field.npages(); p++) {
      for (Index r = 0; r < irradiance_field.nrows(); r++) {
        net_flux_b = (irradiance_field(b - 1, p, r, 0) +
                      irradiance_field(b - 1, p, r, 1));
        net_flux_t = (irradiance_field(b + 1, p, r, 0) +
                      irradiance_field(b + 1, p, r, 1));

        heating_rates(b, p, r) = (net_flux_t - net_flux_b) /
                                 (p_grid[b + 1] - p_grid[b - 1]) * g0 /
                                 specific_heat_capacity(b, p, r);
      }
    }
  }

  idx = irradiance_field.nbooks();

  // now calculate the heating rates for the upper and lower boundary
  for (Index p = 0; p < irradiance_field.npages(); p++) {
    for (Index r = 0; r < irradiance_field.nrows(); r++) {
      // lower boundary
      net_flux_b =
          (irradiance_field(0, p, r, 0) + irradiance_field(0, p, r, 1));
      net_flux_c =
          (irradiance_field(1, p, r, 0) + irradiance_field(1, p, r, 1));
      net_flux_t =
          (irradiance_field(2, p, r, 0) + irradiance_field(0, p, r, 1));

      heating_rates(0, p, r) = (-3 * net_flux_b + 4 * net_flux_c - net_flux_t) /
                               (p_grid[2] - p_grid[0]) * g0 /
                               specific_heat_capacity(0, p, r);

      // lower boundary
      net_flux_t = (irradiance_field(idx - 1, p, r, 0) +
                    irradiance_field(idx - 1, p, r, 1));
      net_flux_c = (irradiance_field(idx - 2, p, r, 0) +
                    irradiance_field(idx - 2, p, r, 1));
      net_flux_b = (irradiance_field(idx - 3, p, r, 0) +
                    irradiance_field(idx - 3, p, r, 1));

      heating_rates(idx - 1, p, r) =
          -(-3 * net_flux_t + 4 * net_flux_c - net_flux_b) /
          (p_grid[2] - p_grid[0]) * g0 / specific_heat_capacity(0, p, r);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void irradiance_fieldFromRadiance(Tensor4& irradiance_field,
                                  const Tensor5& radiance_field,
                                  const Vector& za_grid,
                                  const Vector& aa_grid,
                                  const Vector& za_grid_weights,
                                  const Verbosity&) {
  // Number of zenith angles.
  const Index N_scat_za = za_grid.nelem();
  const Index N_scat_aa = aa_grid.nelem();

  Tensor4 radiance_field_aa_integrated;

  //azimuth integration
  if (N_scat_aa == 1)  //1D case no azimuth dependency
  {
    radiance_field_aa_integrated =
        radiance_field(joker, joker, joker, joker, 0);
    radiance_field_aa_integrated *= 2 * PI;

  } else  //general case with azimuth dependency
  {
    radiance_field_aa_integrated.resize(radiance_field.nshelves(),
                                        radiance_field.nbooks(),
                                        radiance_field.npages(),
                                        radiance_field.nrows());
    radiance_field_aa_integrated = 0.;

    for (Index b = 0; b < radiance_field_aa_integrated.nbooks(); b++) {
      for (Index p = 0; p < radiance_field_aa_integrated.npages(); p++) {
        for (Index r = 0; r < radiance_field_aa_integrated.nrows(); r++) {
          for (Index c = 0; c < radiance_field_aa_integrated.ncols(); c++) {
            for (Index i = 0; i < N_scat_aa - 1; i++) {
              radiance_field_aa_integrated(b, p, r, c) +=
                  (radiance_field(b, p, r, c, i) +
                   radiance_field(b, p, r, c, i + 1)) /
                  2. * abs(aa_grid[i + 1] - aa_grid[i]) * DEG2RAD;
            }
          }
        }
      }
    }
  }

  //allocate
  irradiance_field.resize(radiance_field.nshelves(),
                          radiance_field.nbooks(),
                          radiance_field.npages(),
                          2);
  irradiance_field = 0;

  // zenith angle integration

  for (Index b = 0; b < irradiance_field.nbooks(); b++) {
    for (Index p = 0; p < irradiance_field.npages(); p++) {
      for (Index r = 0; r < irradiance_field.nrows(); r++) {
        for (Index i = 0; i < N_scat_za; i++) {
          if (za_grid[i] <= 90.) {
            irradiance_field(b, p, r, 0) +=
                radiance_field_aa_integrated(b, p, r, i) *
                cos(za_grid[i] * DEG2RAD) * (-1.) * za_grid_weights[i];
          } else {
            irradiance_field(b, p, r, 1) +=
                radiance_field_aa_integrated(b, p, r, i) *
                cos(za_grid[i] * DEG2RAD) * (-1.) * za_grid_weights[i];
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void RadiationFieldSpectralIntegrate(Tensor4& radiation_field,
                                     const Vector& f_grid,
                                     const Tensor5& spectral_radiation_field,
                                     const Verbosity&) {
  if (f_grid.nelem() != spectral_radiation_field.nshelves()) {
    throw runtime_error(
        "The length of f_grid does not match with\n"
        " the first dimension of the spectral_radiation_field");
  }

  //allocate
  radiation_field.resize(spectral_radiation_field.nbooks(),
                         spectral_radiation_field.npages(),
                         spectral_radiation_field.nrows(),
                         spectral_radiation_field.ncols());
  radiation_field = 0;

  // frequency integration
  for (Index i = 0; i < spectral_radiation_field.nshelves() - 1; i++) {
    const Numeric df = f_grid[i + 1] - f_grid[i];

    for (Index b = 0; b < radiation_field.nbooks(); b++) {
      for (Index p = 0; p < radiation_field.npages(); p++) {
        for (Index r = 0; r < radiation_field.nrows(); r++) {
          for (Index c = 0; c < radiation_field.ncols(); c++) {
            radiation_field(b, p, r, c) +=
                (spectral_radiation_field(i + 1, b, p, r, c) +
                 spectral_radiation_field(i, b, p, r, c)) /
                2 * df;
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void RadiationFieldSpectralIntegrate(Tensor5& radiation_field,
                                     const Vector& f_grid,
                                     const Tensor7& spectral_radiation_field,
                                     const Verbosity&) {
  if (f_grid.nelem() != spectral_radiation_field.nlibraries()) {
    throw runtime_error(
        "The length of f_grid does not match with\n"
        " the first dimension of the spectral_radiation_field");
  }

  //allocate
  radiation_field.resize(spectral_radiation_field.nvitrines(),
                         spectral_radiation_field.nshelves(),
                         spectral_radiation_field.nbooks(),
                         spectral_radiation_field.npages(),
                         spectral_radiation_field.nrows());
  radiation_field = 0;

  // frequency integration
  for (Index i = 0; i < spectral_radiation_field.nlibraries() - 1; i++) {
    const Numeric df = f_grid[i + 1] - f_grid[i];

    for (Index s = 0; s < radiation_field.nshelves(); s++) {
      for (Index b = 0; b < radiation_field.nbooks(); b++) {
        for (Index p = 0; p < radiation_field.npages(); p++) {
          for (Index r = 0; r < radiation_field.nrows(); r++) {
            for (Index c = 0; c < radiation_field.ncols(); c++) {
              radiation_field(s, b, p, r, c) +=
                  (spectral_radiation_field(i + 1, s, b, p, r, c, 0) +
                   spectral_radiation_field(i, s, b, p, r, c, 0)) /
                  2 * df;
            }
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void spectral_irradiance_fieldFromSpectralRadianceField(
    Tensor5& spectral_irradiance_field,
    const Tensor7& spectral_radiance_field,
    const Vector& za_grid,
    const Vector& aa_grid,
    const Vector& za_grid_weights,
    const Verbosity&) {
  // Number of zenith angles.
  const Index N_scat_za = za_grid.nelem();
  const Index N_scat_aa = aa_grid.nelem();

  Tensor5 iy_field_aa_integrated;

  //azimuth integration
  if (N_scat_aa == 1)  //1D case no azimuth dependency
  {
    iy_field_aa_integrated =
        spectral_radiance_field(joker, joker, joker, joker, joker, 0, 0);
    iy_field_aa_integrated *= 2 * PI;

  } else  //general case with azimuth dependency
  {
    iy_field_aa_integrated.resize(spectral_radiance_field.nlibraries(),
                                  spectral_radiance_field.nvitrines(),
                                  spectral_radiance_field.nshelves(),
                                  spectral_radiance_field.nbooks(),
                                  spectral_radiance_field.npages());
    iy_field_aa_integrated = 0.;

    for (Index s = 0; s < iy_field_aa_integrated.nshelves(); s++) {
      for (Index b = 0; b < iy_field_aa_integrated.nbooks(); b++) {
        for (Index p = 0; p < iy_field_aa_integrated.npages(); p++) {
          for (Index r = 0; r < iy_field_aa_integrated.nrows(); r++) {
            for (Index c = 0; c < iy_field_aa_integrated.ncols(); c++) {
              for (Index i = 0; i < N_scat_aa - 1; i++) {
                iy_field_aa_integrated(s, b, p, r, c) +=
                    (spectral_radiance_field(s, b, p, r, c, i, 0) +
                     spectral_radiance_field(s, b, p, r, c, i + 1, 0)) /
                    2. * abs(aa_grid[i + 1] - aa_grid[i]) * DEG2RAD;
              }
            }
          }
        }
      }
    }
  }

  //allocate
  spectral_irradiance_field.resize(spectral_radiance_field.nlibraries(),
                                   spectral_radiance_field.nvitrines(),
                                   spectral_radiance_field.nshelves(),
                                   spectral_radiance_field.nbooks(),
                                   2);
  spectral_irradiance_field = 0;

  // zenith angle integration
  for (Index s = 0; s < spectral_irradiance_field.nshelves(); s++) {
    for (Index b = 0; b < spectral_irradiance_field.nbooks(); b++) {
      for (Index p = 0; p < spectral_irradiance_field.npages(); p++) {
        for (Index r = 0; r < spectral_irradiance_field.nrows(); r++) {
          for (Index i = 0; i < N_scat_za; i++) {
            if (za_grid[i] <= 90.) {
              spectral_irradiance_field(s, b, p, r, 0) +=
                  iy_field_aa_integrated(s, b, p, r, i) *
                  cos(za_grid[i] * DEG2RAD) * (-1.) * za_grid_weights[i];
            } else {
              spectral_irradiance_field(s, b, p, r, 1) +=
                  iy_field_aa_integrated(s, b, p, r, i) *
                  cos(za_grid[i] * DEG2RAD) * (-1.) * za_grid_weights[i];
            }
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void spectral_radiance_fieldClearskyPlaneParallel(
    Workspace& ws,
    Tensor7& spectral_radiance_field,
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
    const EnergyLevelMap& nlte_field,
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
    const String& rt_integration_option,
    const Tensor3& surface_props_data,
    const Vector& za_grid,
    const Index& use_parallel_za  [[maybe_unused]],
    const Verbosity& verbosity) {
  // Check input
  if (atmosphere_dim != 1)
    throw runtime_error("This method only works for atmosphere_dim = 1.");

  // Sizes
  const Index nl = p_grid.nelem();
  const Index nf = f_grid.nelem();
  const Index nza = za_grid.nelem();

  // Init spectral_radiance_field and trans_field
  spectral_radiance_field.resize(nf, nl, 1, 1, nza, 1, stokes_dim);
  trans_field.resize(nf, nl, nza);

  // De-activate cloudbox
  const Index cloudbox_on = 0, ppath_inside_cloudbox_do = 0;
  const ArrayOfIndex cloudbox_limits(0);

  // Various variables
  const String iy_unit = "1";
  const ArrayOfString iy_aux_vars(0);
  const Vector rte_pos2(0);
  const Index iy_agenda_call1 = 1;
  const Tensor3 iy_transmittance(0, 0, 0);
  const Index jacobian_do = 0;
  const ArrayOfRetrievalQuantity jacobian_quantities(0);
  // Create one altitude just above TOA
  const Numeric z_space = z_field(nl - 1, 0, 0) + 10;

  ArrayOfString fail_msg;
  bool failed = false;

  // Define iy_main_agenda to be consistent with the assumptions of
  // this method. This definition of iy_main_agenda will be used to when
  // calculating the the radiation reflected by the surface
  Agenda iy_main_agenda(ws.original_workspace);
  iy_main_agenda.append("ppathPlaneParallel", TokVal());
  iy_main_agenda.append("iyEmissionStandard", TokVal());
  iy_main_agenda.push_back(MRecord(global_data::MdMap.at("VectorSet"),
                                   {ws.WsvMap_ptr->at("geo_pos")},
                                   {},
                                   Vector{},
                                   Agenda{ws.original_workspace}));
  iy_main_agenda.set_name("iy_main_agenda");
  iy_main_agenda.check(ws, verbosity);

  // Index in p_grid where field at surface shall be placed
  const Index i0 = index_of_zsurface(z_surface(0, 0), z_field(joker, 0, 0));

  // Loop zenith angles
  //
  if (nza)
#pragma omp parallel for if (!arts_omp_in_parallel() && nza > 1 && \
                             use_parallel_za) firstprivate(ws)
    for (Index i = 0; i < nza; i++) {
      if (failed) continue;
      try {
        // Define output variables
        Ppath ppath;
        Vector ppvar_p, ppvar_t;
        Matrix iy, ppvar_vmr, ppvar_wind, ppvar_mag, ppvar_f;
        EnergyLevelMap ppvar_nlte;
        Tensor3 ppvar_iy;
        Tensor4 ppvar_trans_cumulat, ppvar_trans_partial;
        ArrayOfMatrix iy_aux;
        ArrayOfTensor3 diy_dx;

        Index iy_id = i;
        Vector rte_los(1, za_grid[i]);
        Vector rte_pos(1, za_grid[i] < 90 ? z_surface(0, 0) : z_space);

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
        ARTS_ASSERT(ppath.gp_p[ppath.np - 1].idx == i0 ||
               ppath.gp_p[ppath.np - 1].idx == nl - 2);

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
        ARTS_ASSERT(iy.ncols() == stokes_dim);

        // First and last points are most easily handled separately
        if (za_grid[i] < 90) {
          spectral_radiance_field(joker, i0, 0, 0, i, 0, joker) =
              ppvar_iy(joker, joker, 0);
          spectral_radiance_field(joker, nl - 1, 0, 0, i, 0, joker) =
              ppvar_iy(joker, joker, ppath.np - 1);
          trans_field(joker, 0, i) = ppvar_trans_partial(0, joker, 0, 0);
          trans_field(joker, nl - 1, i) =
              ppvar_trans_partial(ppath.np - 1, joker, 0, 0);
        } else {
          spectral_radiance_field(joker, nl - 1, 0, 0, i, 0, joker) =
              ppvar_iy(joker, joker, 0);
          spectral_radiance_field(joker, i0, 0, 0, i, 0, joker) =
              ppvar_iy(joker, joker, ppath.np - 1);
          trans_field(joker, nl - 1, i) = ppvar_trans_partial(0, joker, 0, 0);
          trans_field(joker, 0, i) =
              ppvar_trans_partial(ppath.np - 1, joker, 0, 0);
        }

        // Remaining points
        for (Index p = 1; p < ppath.np - 1; p++) {
          // We just store values at pressure levels
          if (ppath.gp_p[p].fd[0] < 1e-2) {
            spectral_radiance_field(
                joker, ppath.gp_p[p].idx, 0, 0, i, 0, joker) =
                ppvar_iy(joker, joker, p);
            trans_field(joker, ppath.gp_p[p].idx, i) =
                ppvar_trans_partial(p, joker, 0, 0);
          }
        }

        // We don't want undefined values to possibly affect an interpolation,
        // and to be safe we set the fields for underground levels to equal the
        // ones at the surface
        for (Index p = 0; p < i0; p++) {
          spectral_radiance_field(joker, p, 0, 0, i, 0, joker) =
              spectral_radiance_field(joker, i0, 0, 0, i, 0, joker);
          trans_field(joker, p, i) = trans_field(joker, i0, i);
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

/* Workspace method: Doxygen documentation will be auto-generated */
void spectral_radiance_fieldCopyCloudboxField(
    Tensor7& spectral_radiance_field,
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Index& cloudbox_on,
    const ArrayOfIndex& cloudbox_limits,
    const Tensor7& cloudbox_field,
    const Verbosity&) {
  if (atmosphere_dim > 1)
    throw runtime_error("This method can only be used for 1D calculations.\n");
  if (!cloudbox_on)
    throw runtime_error("Cloudbox is off. This is not handled by this method.");
  if (cloudbox_limits[0] || cloudbox_limits[1] != p_grid.nelem() - 1)
    throw runtime_error(
        "The cloudbox must cover all pressure levels "
        "to use this method.");

  // If all OK, it is just to copy
  spectral_radiance_field = cloudbox_field;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void spectral_radiance_fieldExpandCloudboxField(
    Workspace& ws,
    Tensor7& spectral_radiance_field,
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
    const EnergyLevelMap& nlte_field,
    const Tensor4& vmr_field,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const Tensor3& wind_u_field,
    const Tensor3& wind_v_field,
    const Tensor3& wind_w_field,
    const Tensor3& mag_u_field,
    const Tensor3& mag_v_field,
    const Tensor3& mag_w_field,
    const Matrix& z_surface,
    const Index& cloudbox_on,
    const ArrayOfIndex& cloudbox_limits,
    const Tensor7& cloudbox_field,
    const Numeric& ppath_lmax,
    const Numeric& rte_alonglos_v,
    const String& rt_integration_option,
    const Tensor3& surface_props_data,
    const Vector& za_grid,
    const Index& use_parallel_za  [[maybe_unused]],
    const Verbosity& verbosity) {
  // Check input
  if (atmosphere_dim != 1)
    throw runtime_error("This method only works for atmosphere_dim = 1.");
  if (!cloudbox_on)
    throw runtime_error("No ned to use this method with cloudbox=0.");
  if (cloudbox_limits[0])
    throw runtime_error(
        "The first element of *cloudbox_limits* must be zero "
        "to use this method.");

  // Sizes
  const Index nl = p_grid.nelem();
  const Index nf = f_grid.nelem();
  const Index nza = za_grid.nelem();

  // Init spectral_radiance_field
  spectral_radiance_field.resize(nf, nl, 1, 1, nza, 1, stokes_dim);

  // and copy the part taken from cloudbox_field
  spectral_radiance_field(joker,
                          Range(0, cloudbox_limits[1] + 1),
                          joker,
                          joker,
                          joker,
                          joker,
                          joker) = cloudbox_field;

  // Various variables
  const Index ppath_inside_cloudbox_do = 0;
  const String iy_unit = "1";
  const ArrayOfString iy_aux_vars(0);
  const Vector rte_pos2(0);
  const Index iy_agenda_call1 = 1;
  const Tensor3 iy_transmittance(0, 0, 0);
  const Index jacobian_do = 0;
  const ArrayOfRetrievalQuantity jacobian_quantities(0);
  // Create one altitude just above TOA
  const Numeric z_space = z_field(nl - 1, 0, 0) + 10;

  ArrayOfString fail_msg;
  bool failed = false;

  // Define iy_main_agenda to be consistent with the assumptions of
  // this method (but the agenda will not be used).
  Agenda iy_main_agenda{ws.original_workspace};
  iy_main_agenda.append("ppathPlaneParallel", TokVal());
  iy_main_agenda.append("iyEmissionStandard", TokVal());
  iy_main_agenda.push_back(MRecord(global_data::MdMap.at("VectorSet"),
                                   {ws.WsvMap_ptr->at("geo_pos")},
                                   {},
                                   Vector{},
                                   Agenda{ws.original_workspace}));
  iy_main_agenda.set_name("iy_main_agenda");
  iy_main_agenda.check(ws, verbosity);

  // Variables related to the top of the cloudbox
  const Index i0 = cloudbox_limits[1];
  const Numeric z_top = z_field(i0 + 1, 0, 0);  // Note i0+1

  // Loop zenith angles
  //
  if (nza)
#pragma omp parallel for if (!arts_omp_in_parallel() && nza > 1 && \
                             use_parallel_za) firstprivate(ws)
    for (Index i = 0; i < nza; i++) {
      if (failed) continue;
      try {
        // Define output variables
        Ppath ppath;
        Vector ppvar_p, ppvar_t;
        Matrix iy, ppvar_vmr, ppvar_wind, ppvar_mag, ppvar_f;
        EnergyLevelMap ppvar_nlte;
        Tensor3 ppvar_iy;
        Tensor4 ppvar_trans_cumulat, ppvar_trans_partial;
        ArrayOfMatrix iy_aux;
        ArrayOfTensor3 diy_dx;

        Index iy_id = i;
        Vector rte_los(1, za_grid[i]);
        Vector rte_pos(1, za_grid[i] < 90 ? z_top : z_space);

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
        ARTS_ASSERT(ppath.gp_p[ppath.np - 1].idx == i0 ||
               ppath.gp_p[ppath.np - 1].idx == nl - 2);

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
        ARTS_ASSERT(iy.ncols() == stokes_dim);

        // First and last points are most easily handled separately
        // But field at top cloudbox already known from copying above
        if (za_grid[i] < 90) {
          spectral_radiance_field(joker, i0 + 1, 0, 0, i, 0, joker) =
              ppvar_iy(joker, joker, 0);
          spectral_radiance_field(joker, nl - 1, 0, 0, i, 0, joker) =
              ppvar_iy(joker, joker, ppath.np - 1);
        } else {
          spectral_radiance_field(joker, nl - 1, 0, 0, i, 0, joker) =
              ppvar_iy(joker, joker, 0);
        }

        // Remaining points
        for (Index p = 1; p < ppath.np - 1; p++) {
          // We just store values at pressure levels
          if (ppath.gp_p[p].fd[0] < 1e-2) {
            spectral_radiance_field(
                joker, ppath.gp_p[p].idx, 0, 0, i, 0, joker) =
                ppvar_iy(joker, joker, p);
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
