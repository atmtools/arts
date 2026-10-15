/**
 * @file   physics_funcs.cc
 * @author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
 * @date   2002-05-08
 *
 * @brief  This file contains the code of functions of physical character.
 *
 *  Modified by Claudia Emde (2002-05-28).
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "physics_funcs.h"

#include <debug.h>

#include <cmath>
#include <tuple>

inline constexpr Numeric BOLTZMAN_CONST = Constant::boltzmann_constant;
inline constexpr Numeric DEG2RAD        = Conversion::deg2rad(1);
inline constexpr Numeric PLANCK_CONST   = Constant::planck_constant;
inline constexpr Numeric SPEED_OF_LIGHT = Constant::speed_of_light;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/** barometric_heightformula
 *
 *  Barometric heightformula for isothermal earth atmosphere.
 *
 * @param[in] p  Atmospheric pressure at starting level [Pa].
 * @param[in] dh Vertical displacement to starting pressure level [m].
 *
 * @return p1 Pressure in displacement level [Pa].
 *
 * @author Daniel Kreyling
 * @date 2011-01-20
 */
Numeric barometric_heightformula(  //output is p1
    //input
    const Numeric& p,
    const Numeric& dh)

{
  /* taken from: Seite „Barometrische Höhenformel“. In: Wikipedia,
 * Die freie Enzyklopädie. Bearbeitungsstand: 3. April 2011, 20:28 UTC.
 * URL: http://de.wikipedia.org/w/index.php?title=Barometrische_H%C3%B6henformel&oldid=87257486
 * (Abgerufen: 15. April 2011, 15:41 UTC)
 */

  //barometric height formula
  Numeric M = 0.02896;  //mean molar mass of air [kg mol^-1]
  Numeric g = 9.807;    //earth acceleration [kg m s^-1]
  Numeric R = 8.314;    //universal gas constant [J K^−1 mol^−1]
  Numeric T = 253;      //median tropospheric reference temperature [K]

  // calculation
  Numeric p1 = p * exp(-(-dh) / (R * T / (M * g)));

  return p1;
}

/** dinvplanckdI
 *
 * Calculates the derivative of inverse-Planck with respect to intensity.
 *
 * @param[in]  i  Radiance.
 * @param[in]  f  Frequency.
 *
 * @return     The derivative.
 *
 * @author Patrick Eriksson
 * @date   2010-10-26
 */
Numeric dinvplanckdI(const Numeric& i, const Numeric& f) {
  constexpr Numeric a    = PLANCK_CONST / BOLTZMAN_CONST;
  constexpr Numeric b    = 2 * PLANCK_CONST / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
  const Numeric     d    = b * f * f * f / i;
  const Numeric     binv = a * f / log1p(d);

  return binv * binv / (a * f * i * (1 / d + 1));
}

Numeric refractive_index_water_visible_nir_harvey98(const Numeric frequency,
                                                    const Numeric temperature,
                                                    const Numeric density,
                                                    const bool    check_validity) {
  constexpr Numeric reference_temperature = 273.15;  // K
  constexpr Numeric reference_density     = 1000.0;  // kg/m3
  constexpr Numeric reference_wavelength  = 0.589;   // micrometres

  const Numeric wavelength = Conversion::freq2wavelen(frequency) * 1e6;

  if (check_validity) {
    ARTS_USER_ERROR_IF(not(temperature > 261.15 and temperature < 773.15),
                       "Harvey98 water refractive index is valid for 261.15 K < temperature < 773.15 K. "
                       "Got {} K.",
                       temperature)
    ARTS_USER_ERROR_IF(not(density > 0.0 and density < 1060.0),
                       "Harvey98 water refractive index is valid for 0 kg/m3 < density < 1060 kg/m3. "
                       "Got {} kg/m3.",
                       density)
    ARTS_USER_ERROR_IF(not(wavelength > 0.2 and wavelength < 1.9),
                       "Harvey98 water refractive index is valid for 0.2 micrometres < wavelength < "
                       "1.9 micrometres. Got frequency {} Hz (wavelength {} micrometres).",
                       frequency,
                       wavelength)
  }

  constexpr Numeric a0 = 0.244257733;
  constexpr Numeric a1 = 9.74634476e-3;
  constexpr Numeric a2 = -3.73234996e-3;
  constexpr Numeric a3 = 2.68678472e-4;
  constexpr Numeric a4 = 1.58920570e-3;
  constexpr Numeric a5 = 2.45934259e-3;
  constexpr Numeric a6 = 0.900704920;
  constexpr Numeric a7 = -1.66626219e-2;

  constexpr Numeric lambda_uv = 0.2292020;
  constexpr Numeric lambda_ir = 5.432937;

  const Numeric reduced_temperature = temperature / reference_temperature;
  const Numeric reduced_density     = density / reference_density;
  const Numeric reduced_wavelength  = wavelength / reference_wavelength;
  const Numeric wavelength_squared  = Math::pow2(reduced_wavelength);

  Numeric rhs  = a0;
  rhs         += a1 * reduced_density;
  rhs         += a2 * reduced_temperature;
  rhs         += a3 * wavelength_squared * reduced_temperature;
  rhs         += a4 / wavelength_squared;
  rhs         += a5 / (wavelength_squared - Math::pow2(lambda_uv));
  rhs         += a6 / (wavelength_squared - Math::pow2(lambda_ir));
  rhs         += a7 * Math::pow2(reduced_density);

  const Numeric lorentz_lorenz = reduced_density * rhs;
  return std::sqrt((1.0 + 2.0 * lorentz_lorenz) / (1.0 - lorentz_lorenz));
}

/** fresnel
 *
 * Calculates complex AMPLITUDE reflection coeffcients for a specular
 *   reflection.
 *
 *  The properties of the two involved media are given as the complex
 *  refractive index, n. A dielectric constant, eps, is converted as
 *  n = sqrt( eps ). The power reflection coefficient, r, for one
 *  polarisation is r = abs(R)^2.
 *
 *  Snell's law is applied with the complex indices, n1 sin(theta1) =
 *  n2 sin(theta2), so the transmitted cosine is complex for an absorbing
 *  medium and beyond total reflection.  Its root is that of a wave
 *  decaying into the reflecting medium (non-negative real part, and
 *  non-negative imaginary part for a zero real part).
 *
 *  @param[out]  Rv    Reflection coefficient for vertical polarisation.
 *  @param[out]  Rh    Reflection coefficient for vertical polarisation.
 *  @param[in]   n1    Refractive index of medium where radiation propagates.
 *  @param[in]   n2    Refractive index of reflecting medium.
 *  @param[in]   theta Propagation angle from normal of radiation to be.
 *                     reflected [deg]
 *
 *  @author Patrick Eriksson
 *  @date   2004-09-21
 */
void fresnel(Complex& Rv, Complex& Rh, const Complex& n1, const Complex& n2, const Numeric& theta) {
  std::tie(Rv, Rh) = fresnel(n1, n2, theta);
}

std::pair<Complex, Complex> fresnel(const Complex& n1, const Complex& n2, const Numeric& theta) {
  const Numeric theta1    = DEG2RAD * theta;
  const Numeric costheta1 = std::cos(theta1);
  const Complex sintheta2 = n1 * std::sin(theta1) / n2;
  Complex       costheta2 = std::sqrt(1.0 - sintheta2 * sintheta2);
  if (costheta2.real() < 0.0 or (costheta2.real() == 0.0 and costheta2.imag() < 0.0)) costheta2 = -costheta2;

  const Complex a = n2 * costheta1;
  const Complex b = n1 * costheta2;
  const Complex c = n1 * costheta1;
  const Complex d = n2 * costheta2;

  return {(a - b) / (a + b), (c - d) / (c + d)};
}

/** invplanck
 *
 * Converts a radiance to Planck brightness temperature.
 *
 * @param[in]  i   Radiance.
 * @param[in]  f  Frequency.
 *
 * @return     Planck brightness temperature.
 *
 * @author Patrick Eriksson
 * @date   2002-08-11
*/
Numeric invplanck(const Numeric& i, const Numeric& f) {
  constexpr Numeric a = PLANCK_CONST / BOLTZMAN_CONST;
  constexpr Numeric b = 2 * PLANCK_CONST / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);

  return (a * f) / log1p((b * f * f * f) / i);
}

/** invrayjean
 *
 * Converts a radiance to Rayleigh-Jean brightness temperature.
 *
 * @param[in]  i  Radiance.
 * @param[in]  f  Frequency.
 *
 * @return     RJ brightness temperature.
 *
 * @author Patrick Eriksson
 * @date   2000-09-28
 */
Numeric invrayjean(const Numeric& i, const Numeric& f) {
  constexpr Numeric a = SPEED_OF_LIGHT * SPEED_OF_LIGHT / (2 * BOLTZMAN_CONST);

  return (a * i) / (f * f);
}

/** planck
 *
 * Calculates the Planck function for a single temperature.

 * Note that this expression gives the intensity for both polarisations.
 *
 *  @param[in]  f  Frequency.
 *  @param[in]  t  Temperature.
 *
 *  @return     Blackbody radiation.
 *
 *  @author Patrick Eriksson
 *  @date   2000-04-08
 */
Numeric planck(const Numeric& f, const Numeric& t) {
  constexpr Numeric a = 2 * Constant::h / Math::pow2(Constant::c);
  constexpr Numeric b = Constant::h / Constant::k;

  return a * Math::pow3(f) / std::expm1((b * f) / t);
}

/** planck
 *
 * Calculates the Planck function for a single temperature and a vector of
 * frequencies.
 *
 * Note that this expression gives the intensity for both polarisations.
 *
 * @param[in]  f  Frequency.
 * @param[in]  t  Temperature.
 *
 * @return     Blackbody radiation.
 *
 * @author Patrick Eriksson
 * @date   2015-12-15
 */
void planck(StridedVectorView b, const ConstVectorView& f, const Numeric& t) {
  ARTS_USER_ERROR_IF(b.size() not_eq f.size(), "Vector size mismatch: frequency dim is bad")
  for (Size i = 0; i < f.size(); i++) b[i] = planck(f[i], t);
}

/** planck
 *
 * Calculates the Planck function for a single temperature and a vector of
 * frequencies.
 *
 * Note that this expression gives the intensity for both polarisations.
 *
 * @param[in]  f  Frequency.
 * @param[in]  t  Temperature.
 *
 * @return     Blackbody radiation.
 *
 * @author Patrick Eriksson
 * @date   2015-12-15
 */
Vector planck(const StridedConstVectorView& f, const Numeric& t) {
  Vector b(f.size());
  for (Size i = 0; i < f.size(); i++) b[i] = planck(f[i], t);
  return b;
}

/** dplanck_dt
 *
 * Calculates the temperature derivative of the Planck function
 * for a single temperature and frequency.
 *
 * @param[in]  f  Frequency.
 * @param[in]  t  Temperature.
 *
 * @return     Blackbody radiation temperature derivative.
 *
 * @author Richard Larsson
 * @date   2015-09-15
 */
Numeric dplanck_dt(const Numeric& f, const Numeric& t) {
  constexpr Numeric a = 2 * Constant::h / Math::pow2(Constant::c);
  constexpr Numeric b = Constant::h / Constant::k;

  // nb. expm1(x) should be more accurate than exp(x) - 1, so use it
  const Numeric inv_exp_t_m1 = 1.0 / std::expm1(b * f / t);

  return a * b * Math::pow4(f) * inv_exp_t_m1 * (1 + inv_exp_t_m1) / Math::pow2(t);
}

/** dplanck_dt
 * 
 * Calculates the Planck function temperature derivative for a single
 * temperature and a vector of frequencies.
 *
 * @param[in]  f  Frequency.
 * @param[in]  t  Temperature.
 *
 * @return     Blackbody radiation temperature derivative.
 *
 * @author Richard Larsson
 * @date   2019-10-11
 */
void dplanck_dt(VectorView dbdt, const ConstVectorView& f, const Numeric& t) {
  ARTS_USER_ERROR_IF(dbdt.size() not_eq f.size(), "Vector size mismatch: frequency dim is bad")
  for (Size i = 0; i < f.size(); i++) dbdt[i] = dplanck_dt(f[i], t);
}

/** dplanck_dt
 * 
 * Calculates the Planck function temperature derivative for a single
 * temperature and a vector of frequencies.
 *
 * @param[in]  f  Frequency.
 * @param[in]  t  Temperature.
 *
 * @return     Blackbody radiation temperature derivative.
 *
 * @author Richard Larsson
 * @date   2019-10-11
 */
Vector dplanck_dt(const ConstVectorView& f, const Numeric& t) {
  Vector dbdt(f.size());
  for (Size i = 0; i < f.size(); i++) dbdt[i] = dplanck_dt(f[i], t);
  return dbdt;
}

/** dplanck_df
 *
 * Calculates the frequency derivative of the Planck function
 * for a single temperature and frequency.
 *
 * @param[in]  f  Frequency.
 * @param[in]  t  Temperature.
 *
 * @return     Blackbody radiation frequency derivative.
 *
 * @author Richard Larsson
 * @date   2015-09-15
 */
Numeric dplanck_df(const Numeric& f, const Numeric& t) {
  constexpr Numeric a = 2 * Constant::h / Math::pow2(Constant::c);
  constexpr Numeric b = Constant::h / Constant::k;

  const Numeric inv_exp_t_m1 = 1.0 / std::expm1(b * f / t);

  return a * Math::pow2(f) * (3.0 - (b * f / t) * (1 + inv_exp_t_m1)) * inv_exp_t_m1;
}

/** dplanck_df
 * 
 * Calculates the frequency derivative of the Planck function
 * for a single temperature and frequency.
 *
 * @param[in]  f  Frequency.
 * @param[in]  t  Temperature.
 *
 * @return     Blackbody radiation frequency derivative.
 *
 * @author Richard Larsson
 * @date   2015-09-15
 */
Vector dplanck_df(const ConstVectorView& f, const Numeric& t) {
  Vector dbdf(f.size());
  for (Size i = 0; i < f.size(); i++) dbdf[i] = dplanck_df(f[i], t);
  return dbdf;
}

/** rayjean
 *
 * Converts a Rayleigh-Jean brightness temperature to radiance
 *
 * @param[in]  tb  RJ brightness temperature.
 * @param[in]  f   Frequency.
 *
 * @return     Radiance.
 *
 * @author Patrick Eriksson
 * @date   2011-07-13
 */
Numeric rayjean(const Numeric& f, const Numeric& tb) {
  constexpr Numeric a = SPEED_OF_LIGHT * SPEED_OF_LIGHT / (2 * BOLTZMAN_CONST);

  return (f * f) / (a * tb);
}
