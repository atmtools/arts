/*===========================================================================
  ===  File description
  ===========================================================================*/

#include <workspace.h>

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "arts_constants.h"
#include "arts_conversions.h"
#include "atm.h"
#include "gsl_gauss_legendre.h"
#include "math_funcs.h"
#include "matpack_data.h"

/*!
  \file   m_fluxes.cc
  \author Manfred Brath  <manfred.brath@uni-hamburg.de>
  \date   2018-06-22

  \brief  Workspace functions related to simulation of radiation fluxes.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

inline constexpr Numeric PI = Constant::pi;
inline constexpr Numeric DEG2RAD = Conversion::deg2rad(1);

/*===========================================================================
  === The functions
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void AngularGridsSetFluxCalc(AscendingGrid& za_grid,
                             AscendingGrid& aa_grid,
                             Vector& za_grid_weights,
                             // Keywords:
                             const Index& N_za_grid,
                             const Index& N_aa_grid,
                             const String& za_grid_type) {
  // Azimuth angle grid
  if (N_aa_grid > 1){
    Vector aa_grid_temp;
    nlinspace(aa_grid_temp, -180, 180, N_aa_grid);
    aa_grid = std::move(aa_grid_temp);
  }else if (N_aa_grid < 1) {
    std::ostringstream os;
    os << "N_aa_grid must be > 0 (even for 1D).";
    throw std::runtime_error(os.str());
  } else {
    aa_grid.unsafe_resize(1);
    *aa_grid.unsafe_begin() = 0.;
  }

  if (N_za_grid % 2 == 1) {
    std::ostringstream os;
    os << "N_za_grid must be even.";
    throw std::runtime_error(os.str());
  }

  Index nph = N_za_grid / 2;

  //calculate zenith angle grid
  Vector za_grid_temp2(N_za_grid, 0.);
  za_grid_weights.resize(N_za_grid);
  za_grid_weights = 0;

  if (za_grid_type == "double_gauss") {
    Vector x;
    Vector w;
    Vector xtemp;
    Vector wtemp;
    //Numeric theta;

    //calculate legendre weights and evaluation points
    GSL::Integration::GaussLegendre(xtemp, wtemp, nph);

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
      za_grid_temp2[i] = x[i];
      za_grid_temp2[za_grid_temp2.nelem() - 1 - i] = 180 - x[i];

      // set the weights to the right component
      za_grid_weights[i] = w[i];
      za_grid_weights[za_grid_weights.nelem() - 1 - i] = w[i];
    }

  } else if (za_grid_type == "linear") {
    Vector x;
    Vector w;
    calculate_weights_linear(x, w, nph);

    for (Index i = 0; i < N_za_grid; i++) {
      za_grid_temp2[i] = (x[i] + 1) * 90.;

      // set the weights to the right component
      // by adjusting the domain, we also have to adjust the weights
      za_grid_weights[i] = w[i] * sin(za_grid_temp2[i] * DEG2RAD);
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
    za_grid_temp2 = reverse(za_grid_temp);
    za_grid_weights = reverse(w);

  } else {
    std::ostringstream os;
    os << "The selected grid type is not implemented";
    throw std::runtime_error(os.str());
  }

  //be sure that the first and the last angle are within the closed interval
  //between 0 and 180 deg, because ARTS is picky if the angles are due to numerics
  // slightly below and above,respectively.
  if (za_grid_temp2[0] < 0) {
    za_grid_temp2[0] = 0.;
  }

  if (za_grid_temp2[za_grid_temp2.nelem() - 1] > 180) {
    za_grid_temp2[za_grid_temp2.nelem() - 1] = 180.;
  }

  za_grid = std::move(za_grid_temp2);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void heating_ratesFromIrradiance(Tensor3& heating_rates,
                                 const ArrayOfAtmPoint& ppvar_atm,
                                 const Tensor4& irradiance_field,
                                 const Numeric& g0,
                                 const Tensor3& specific_heat_capacity) {
  const Vector p_grid = [&] {
    Vector out(ppvar_atm.size());
    std::transform(
        ppvar_atm.begin(), ppvar_atm.end(), out.begin(), [](auto& point) {
          return point.pressure;
        });
    return out;
  }();

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

      // upper boundary
      net_flux_t = (irradiance_field(idx - 1, p, r, 0) +
                    irradiance_field(idx - 1, p, r, 1));
      net_flux_c = (irradiance_field(idx - 2, p, r, 0) +
                    irradiance_field(idx - 2, p, r, 1));
      net_flux_b = (irradiance_field(idx - 3, p, r, 0) +
                    irradiance_field(idx - 3, p, r, 1));

      heating_rates(idx - 1, p, r) =
          -(-3 * net_flux_t + 4 * net_flux_c - net_flux_b) /
          (p_grid[idx - 1] - p_grid[idx - 3]) * g0 /
          specific_heat_capacity(0, p, r);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void irradiance_fieldFromRadiance(Tensor4& irradiance_field,
                                  const AscendingGrid& za_grid,
                                  const AscendingGrid& aa_grid,
                                  const Vector& za_grid_weights,
                                  const Tensor5& radiance_field) {
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
                                     const AscendingGrid& f_grid,
                                     const Tensor5& spectral_radiation_field) {
  if (f_grid.nelem() != spectral_radiation_field.nshelves()) {
    throw std::runtime_error(
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
                                     const AscendingGrid& f_grid,
                                     const Tensor7& spectral_radiation_field) {
  if (f_grid.nelem() != spectral_radiation_field.nlibraries()) {
    throw std::runtime_error(
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
    const AscendingGrid& za_grid,
    const AscendingGrid& aa_grid,
    const Vector& za_grid_weights) {
  // Number of zenith angles.
  const Index N_scat_za = spectral_radiance_field.npages();
  const Index N_scat_aa = spectral_radiance_field.nrows();

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
void spectral_radiance_fieldCopyCloudboxField(
    Tensor7& spectral_radiance_field,
    const Vector& p_grid,
    const Index& cloudbox_on,
    const ArrayOfIndex& cloudbox_limits,
    const Tensor7& cloudbox_field) {
  throw std::runtime_error(
      "This method can only be used for 1D calculations.\n");
  if (!cloudbox_on)
    throw std::runtime_error(
        "Cloudbox is off. This is not handled by this method.");
  if (cloudbox_limits[0] || cloudbox_limits[1] != p_grid.nelem() - 1)
    throw std::runtime_error(
        "The cloudbox must cover all pressure levels "
        "to use this method.");

  // If all OK, it is just to copy
  spectral_radiance_field = cloudbox_field;
}
