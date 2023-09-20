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
#include "arts_constants.h"
#include "arts_conversions.h"
#include "mystring.h"
#include "physics_funcs.h"
#include <cmath>
#include <stdexcept>

inline constexpr Numeric BOLTZMAN_CONST=Constant::boltzmann_constant;
inline constexpr Numeric DEG2RAD=Conversion::deg2rad(1);
inline constexpr Numeric PLANCK_CONST=Constant::planck_constant;
inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;

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
  ARTS_USER_ERROR_IF (i <= 0, "Non-positive radiance")
  ARTS_USER_ERROR_IF (f <= 0, "Non-positive frequency")

  static const Numeric a = PLANCK_CONST / BOLTZMAN_CONST;
  static const Numeric b = 2 * PLANCK_CONST / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
  const Numeric d = b * f * f * f / i;
  const Numeric binv = a * f / log1p(d);

  return binv * binv / (a * f * i * (1 / d + 1));
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
 *  @param[out]  Rv    Reflection coefficient for vertical polarisation.
 *  @param[out]  Rh    Reflection coefficient for vertical polarisation.
 *  @param[in]   n1    Refractive index of medium where radiation propagates.
 *  @param[in]   n2    Refractive index of reflecting medium.
 *  @param[in]   theta Propagation angle from normal of radiation to be.
 *                     reflected
 *
 *  @author Patrick Eriksson
 *  @date   2004-09-21
 */
void fresnel(Complex& Rv,
             Complex& Rh,
             const Complex& n1,
             const Complex& n2,
             const Numeric& theta) {
  const Numeric theta1 = DEG2RAD * theta;
  const Numeric costheta1 = cos(theta1);
  const Numeric costheta2 = cos(asin(n1.real() * sin(theta1) / n2.real()));

  Complex a, b;
  a = n2 * costheta1;
  b = n1 * costheta2;
  Rv = (a - b) / (a + b);
  a = n1 * costheta1;
  b = n2 * costheta2;
  Rh = (a - b) / (a + b);
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
  ARTS_USER_ERROR_IF (i <= 0, "Non-positive radiance")
  ARTS_USER_ERROR_IF (f < 0, "Non-positive frequency")

  static const Numeric a = PLANCK_CONST / BOLTZMAN_CONST;
  static const Numeric b = 2 * PLANCK_CONST / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);

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
  //   ARTS_USER_ERROR_IF (i <  0, "Negative radiance")
  ARTS_USER_ERROR_IF (f <= 0, "Non-positive frequency")

  static const Numeric a =
      SPEED_OF_LIGHT * SPEED_OF_LIGHT / (2 * BOLTZMAN_CONST);

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
  ARTS_USER_ERROR_IF (t <= 0, "Non-positive temperature")
  ARTS_USER_ERROR_IF (f <= 0, "Non-positive frequency")
  
  constexpr Numeric a = 2 * Constant::h / Math::pow2(Constant::c);;
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
void planck(VectorView b, const ConstVectorView& f, const Numeric& t) {
  ARTS_USER_ERROR_IF (b.nelem() not_eq f.nelem(),
                      "Vector size mismatch: frequency dim is bad")
  for (Index i = 0; i < f.nelem(); i++) b[i] = planck(f[i], t);
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
Vector planck(const ConstVectorView& f, const Numeric& t) {
  Vector b(f.nelem());
  for (Index i = 0; i < f.nelem(); i++) b[i] = planck(f[i], t);
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
  ARTS_USER_ERROR_IF (t <= 0, "Non-positive temperature")
  ARTS_USER_ERROR_IF (f <= 0, "Non-positive frequency")

  constexpr Numeric a = 2 * Constant::h / Math::pow2(Constant::c);;
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
  ARTS_USER_ERROR_IF (dbdt.nelem() not_eq f.nelem(),
                      "Vector size mismatch: frequency dim is bad")
  for (Index i = 0; i < f.nelem(); i++) dbdt[i] = dplanck_dt(f[i], t);
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
  Vector dbdt(f.nelem());
  for (Index i = 0; i < f.nelem(); i++) dbdt[i] = dplanck_dt(f[i], t);
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
Numeric dplanck_df(const Numeric& f, const Numeric& t)  {
  ARTS_USER_ERROR_IF (t <= 0, "Non-positive temperature")
  ARTS_USER_ERROR_IF (f <= 0, "Non-positive frequency")
  
  constexpr Numeric a = 2 * Constant::h / Math::pow2(Constant::c);;
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
Vector dplanck_df(const ConstVectorView& f, const Numeric& t)  {
  Vector dbdf(f.nelem());
  for (Index i = 0; i < f.nelem(); i++) dbdf[i] = dplanck_df(f[i], t);
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
  ARTS_USER_ERROR_IF (tb <= 0, "Non-positive temperature")
  ARTS_USER_ERROR_IF (f < 0, "Negative frequency")

  static const Numeric a =
      SPEED_OF_LIGHT * SPEED_OF_LIGHT / (2 * BOLTZMAN_CONST);

  return (f * f) / (a * tb);
}
