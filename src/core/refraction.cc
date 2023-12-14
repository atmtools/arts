/*===========================================================================
  === File description
  ===========================================================================*/

/*!
  \file   refraction.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2003-01-17
  
  \brief  Functions releated to calculation of refractive index.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "refraction.h"

#include <cmath>

#include "arts_constants.h"
#include "arts_conversions.h"
#include "check_input.h"
#include "matpack_complex.h"

inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric RAD2DEG=Conversion::rad2deg(1);
inline constexpr Numeric TEMP_0_C=Constant::temperature_at_0c;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! complex_n_water_liebe93
/*! 
  Complex refractive index of liquid water according to Liebe 1993.

  The method treats liquid water without salt. Thus, not valid below 10 GHz.
  Upper frequency limit not known, here set to 1000 GHz. Model parameters taken
  from Atmlab function epswater93 (by C. Maetzler), which refer to Liebe 1993
  without closer specifications.
 
  Temperature must be between 0 and 100 degrees Celsius.

  The output matrix has two columns, where column 0 is real part and column 1
  is imaginary part. And rows matches f_grid.
   
   \param   complex_n   Out: Complex refractive index.        
   \param   f_grid      As the WSV with the same name.
   \param   t           Temperature

   \author Patrick Eriksson
   \date   2003-08-15
*/
void complex_n_water_liebe93(Matrix& complex_n,
                             const Vector& f_grid,
                             const Numeric& t) {
  chk_if_in_range("t", t, TEMP_0_C - 40, TEMP_0_C + 100);
  chk_if_in_range("min of f_grid", min(f_grid), 10e9, 1000e9);
  chk_if_in_range("max of f_grid", max(f_grid), 10e9, 1000e9);

  const Index nf = f_grid.nelem();

  complex_n.resize(nf, 2);

  // Implementation following epswater93.m (by C. Mätzler), part of Atmlab,
  // but numeric values strictly following the paper version (146, not 146.4)
  const Numeric theta = 1 - 300 / t;
  const Numeric e0 = 77.66 - 103.3 * theta;
  const Numeric e1 = 0.0671 * e0;
  const Numeric f1 = 20.2 + 146 * theta + 316 * theta * theta;
  const Numeric e2 = 3.52;
  const Numeric f2 = 39.8 * f1;

  for (Index iv = 0; iv < nf; iv++) {
    const Complex ifGHz(0.0, f_grid[iv] / 1e9);

    Complex n = sqrt(e2 + (e1 - e2) / (Numeric(1.0) - ifGHz / f2) +
                     (e0 - e1) / (Numeric(1.0) - ifGHz / f1));

    complex_n(iv, 0) = n.real();
    complex_n(iv, 1) = n.imag();
  }
}

//! complex_n_ice_matzler06
/*! 
  Complex refractive index of water ice according to Matzler 2006 (equivalent to
  Warren 2008).

  The method treats pure water ice (no impurities like salt). Valid from 10 MHz
  up to 3 THz. Thus, not valid below 10 GHz. Follows the atmlab implementation,
  including some relaxation of upper temperature limit to 280K.

  The output matrix has two columns, where column 0 is real part and column 1
  is imaginary part; rows match f_grid.
   
   \param   complex_n   Out: Complex refractive index.        
   \param   f_grid      As the WSV with the same name.
   \param   t           Temperature

   \author Jana Mendrok
   \date   2016-03-21
*/
void complex_n_ice_matzler06(Matrix& complex_n,
                             const Vector& f_grid,
                             const Numeric& t) {
  chk_if_in_range("t", t, 20., 280.);
  chk_if_in_range("min of f_grid", min(f_grid), 10e6, 3000e9);
  chk_if_in_range("max of f_grid", max(f_grid), 10e6, 3000e9);

  const Index nf = f_grid.nelem();

  complex_n.resize(nf, 2);

  // some parametrization constants
  const Numeric B1 = 0.0207;
  const Numeric B2 = 1.16e-11;
  const Numeric b = 335.;

  const Numeric deltabeta = exp(-9.963 + 0.0372 * (t - 273));
  const Numeric ebdt = exp(b / t);
  const Numeric betam = (B1 / t) * ebdt / ((ebdt - 1.) * (ebdt - 1.));

  const Numeric theta = 300. / t - 1;
  const Numeric alfa = (0.00504 + 0.0062 * theta) * exp(-22.1 * theta);
  const Numeric reps = 3.1884 + 9.1e-4 * (t - 273);

  for (Index iv = 0; iv < nf; iv++) {
    Numeric f = f_grid[iv] / 1e9;
    Numeric beta = betam + B2 * f * f + deltabeta;
    Numeric ieps = alfa / f + beta * f;

    Complex eps(reps, ieps);
    Complex n = sqrt(eps);
    complex_n(iv, 0) = n.real();
    complex_n(iv, 1) = n.imag();
  }
}

void refractive_index_water_and_steam_VisNIR(Numeric& n,
                                             const Index& only_valid_range,
                                             const Numeric& frequency,
                                             const Numeric& temperature,
                                             const Numeric& density) {
  //convert frequency to wavelength
  const Numeric wavelength = Conversion::freq2wavelen(frequency) * 1e6;  // [µm]

  //Reference values
  const Numeric T_star = 273.15;      // [K]
  const Numeric rho_star = 1000;      // [kg/m^3]
  const Numeric lambda_star = 0.589;  // [µm]

  //check input
  bool T_ok = (temperature > T_star - 12) && (temperature < T_star + 500);
  bool rho_ok = (density > 0) && (density < 1060);
  bool wvl_ok = (wavelength > 0.2) && (wavelength < 1.9);

  if (only_valid_range) {
    ARTS_USER_ERROR_IF(!T_ok,
                       "Refractive index is calculated outside range of "
                       "validity \n",
                       "Temperature must be between ",
                       T_star - 12,
                       "K and ",
                       T_star + 500,
                       "K\n"
                       "Desired temperature: ",
                       temperature,
                       "K \n")

    ARTS_USER_ERROR_IF(!rho_ok,
                       "Refractive index is calculated outside range of "
                       "validity \n",
                       "Density must be between ",
                       0,
                       "kg m^-3 and ",
                       1060,
                       "kg m^-3\n"
                       "Desired density: ",
                       density,
                       "kg m^-3 \n")

    Numeric frq_upper_limit = Conversion::wavelen2freq(0.2 * 1e-6);
    Numeric frq_lower_limit = Conversion::wavelen2freq(1.9 * 1e-6);

    ARTS_USER_ERROR_IF(!wvl_ok,
                       "Refractive index is calculated outside range of "
                       "validity \n",
                       "Frequency must be between ",
                       frq_lower_limit,
                       "Hz and ",
                       frq_upper_limit,
                       "Hz\n"
                       "Desired density: ",
                       frequency,
                       "Hz \n")
  }

  //coefficients
  const Numeric a0 = 0.244257733;
  const Numeric a1 = 9.74634476e-3;
  const Numeric a2 = -3.73234996e-3;
  const Numeric a3 = 2.68678472e-4;
  const Numeric a4 = 1.58920570e-3;
  const Numeric a5 = 2.45934259e-3;
  const Numeric a6 = 0.900704920;
  const Numeric a7 = -1.66626219e-2;

  const Numeric lambda_uv = 0.2292020;
  const Numeric lambda_ir = 5.432937;

  //normalize input variables
  Numeric T_bar = temperature / T_star;
  Numeric rho_bar = density / rho_star;
  Numeric lambda_bar = wavelength / lambda_star;

  //set up right hands side of eq A1
  Numeric rhs = a0;
  rhs += a1 * rho_bar;
  rhs += a2 * T_bar;
  rhs += a3 * lambda_bar * lambda_bar * T_bar;
  rhs += a4 / (lambda_bar * lambda_bar);
  rhs += a5 / (lambda_bar * lambda_bar - lambda_uv * lambda_uv);
  rhs += a6 / (lambda_bar * lambda_bar - lambda_ir * lambda_ir);
  rhs += a7 * rho_bar * rho_bar;

  Complex a;
  a = rho_bar * rhs;

  //solve Eq A1 (see paper in documentation after n
  Complex n_cmplx = sqrt(-2 * a - 1);
  n_cmplx /= sqrt(a - 1);

  //refractive index
  n = - real(n_cmplx);
}
