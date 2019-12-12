/* Copyright (C) 2002-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>

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
#include <cmath>
#include <stdexcept>
#include "messages.h"
#include "mystring.h"
#include "physics_funcs.h"

extern const Numeric BOLTZMAN_CONST;
extern const Numeric DEG2RAD;
extern const Numeric PLANCK_CONST;
extern const Numeric SPEED_OF_LIGHT;

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
Numeric dinvplanckdI(const Numeric& i, const Numeric& f) try {
  if (i <= 0) throw "Non-positive radiance";
  if (f <= 0) throw "Non-positive frequency";

  static const Numeric a = PLANCK_CONST / BOLTZMAN_CONST;
  static const Numeric b = 2 * PLANCK_CONST / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
  const Numeric d = b * f * f * f / i;
  const Numeric binv = a * f / log(d + 1);

  return binv * binv / (a * f * i * (1 / d + 1));
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *dinvplanckdI* internal function:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
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
Numeric invplanck(const Numeric& i, const Numeric& f) try {
  if (i <= 0) throw "Non-positive radiance";
  if (f < 0) throw "Non-positive frequency";

  static const Numeric a = PLANCK_CONST / BOLTZMAN_CONST;
  static const Numeric b = 2 * PLANCK_CONST / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);

  return (a * f) / log((b * f * f * f) / i + 1.0);
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *invplanck* internal function:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
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
Numeric invrayjean(const Numeric& i, const Numeric& f) try {
  //   if(i <  0) throw "Negative radiance";
  if (f <= 0) throw "Non-positive frequency";

  static const Numeric a =
      SPEED_OF_LIGHT * SPEED_OF_LIGHT / (2 * BOLTZMAN_CONST);

  return (a * i) / (f * f);
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *invrayjean* internal function:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
}

/** number_density
 *
 * Calculates the atmospheric number density.
 *
 * @param[in]  p  Pressure.
 * @param[in]  t  Temperature.
 *
 * @return     Number density.
 *
 * @author Patrick Eriksson
 * @date   2000-04-08
*/
Numeric number_density(const Numeric& p, const Numeric& t) try {
  if (p < 0) throw "Negative pressure";
  if (t <= 0) throw "Non-positive temperature";

  return p / (t * BOLTZMAN_CONST);
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *number_density* internal function:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
}

/** dnumber_density_dT
 *
 * Calculates the atmospheric number density derivative with temperature.
 *
 * @param[in]  p  Pressure.
 * @param[in]  t  Temperature.
 *
 * @return     Number density.
 *
 * @author Richard Larsson
 * @date   2015-09-22
 */
Numeric dnumber_density_dt(const Numeric& p, const Numeric& t) try {
  if (p < 0) throw "Negative pressure";
  if (t <= 0) throw "Non-positive temperature";

  return -p / (t * BOLTZMAN_CONST * t);
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *dnumber_density_dt* internal function:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
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
Numeric planck(const Numeric& f, const Numeric& t) try {
  if (t <= 0) throw "Non-positive temperature";
  if (f <= 0) throw "Non-positive frequency";

  static const Numeric a = 2 * PLANCK_CONST / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
  static const Numeric b = PLANCK_CONST / BOLTZMAN_CONST;

  return (a * f * f * f) / (exp((b * f) / t) - 1.0);
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *planck* internal function:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
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
void planck(VectorView b, ConstVectorView f, const Numeric& t) try {
  if (b.nelem() not_eq f.nelem())
    throw "Vector size mismatch: frequency dim is bad";

  for (Index i = 0; i < f.nelem(); i++) b[i] = planck(f[i], t);
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *planck* internal function:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *planck* internal function:\n";
  os << e.what();

  Index n = 0;
  for (auto& F : f)
    if (F <= 0) n++;
  if (n)
    os << '\t' << "You have " << n
       << " frequency grid points that reports a negative frequency!\n";

  throw std::runtime_error(os.str());
}

/** dplanck_dt
 *
 * Calculates the temperature derivative of the Planck function
 * for a single temperature and frequency.
 *
 * @param[in]  f  Frequency.
 * @param[in]  t  Temperature.
 *
 * @return     Blackbody radiation.
 *
 * @author Richard Larsson
 * @date   2015-09-15
 */
Numeric dplanck_dt(const Numeric& f, const Numeric& t) try {
  if (t <= 0) throw "Non-positive temperature";
  if (f <= 0) throw "Non-positive frequency";

  static const Numeric a = 2 * PLANCK_CONST / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
  static const Numeric b = PLANCK_CONST / BOLTZMAN_CONST;

  const Numeric exp_t = exp(b * f / t);
  const Numeric exp_t_m1 = exp_t - 1.0;
  const Numeric f2 = f * f;

  return a * b * f2 * f2 * exp_t / (t * t * exp_t_m1 * exp_t_m1);
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *dplanck_dt* internal function:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
}

/** dplanck_dt
 * 
 * Calculates the Planck function temperature derivative for a single
 * temperature and a vector of frequencies.
 *
 * @param[in]  f  Frequency.
 * @param[in]  t  Temperature.
 *
 * @return     Blackbody radiation.
 *
 * @author Richard Larsson
 * @date   2019-10-11
 */
void dplanck_dt(VectorView dbdt, ConstVectorView f, const Numeric& t) try {
  if (dbdt.nelem() not_eq f.nelem())
    throw "Vector size mismatch: frequency dim is bad";
  
  for (Index i = 0; i < f.nelem(); i++) dbdt[i] = dplanck_dt(f[i], t);
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *planck* internal function:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *planck* internal function:\n";
  os << e.what();
  
  Index n = 0;
  for (auto& F : f)
    if (F <= 0) n++;
  if (n)
    os << '\t' << "You have " << n
    << " frequency grid points that reports a negative frequency!\n";
    
  throw std::runtime_error(os.str());
}

/** dplanck_df
 *
 * Calculates the frequency derivative of the Planck function
 * for a single temperature and frequency.
 *
 * @param[in]  f  Frequency.
 * @param[in]  t  Temperature.
 *
 * @return     Blackbody radiation.
 *
 * @author Richard Larsson
 * @date   2015-09-15
 */
Numeric dplanck_df(const Numeric& f, const Numeric& t) try {
  if (t <= 0) throw "Non-positive temperature";
  if (f <= 0) throw "Non-positive frequency";

  static const Numeric a = 2 * PLANCK_CONST / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
  static const Numeric b = PLANCK_CONST / BOLTZMAN_CONST;

  const Numeric exp_t = exp(b * f / t);
  const Numeric exp_t_m1 = exp_t - 1.0;

  return -(a * f * f * (3.0 * t - 3.0 * t * exp_t + b * f * exp_t)) /
         (t * exp_t_m1 * exp_t_m1);
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *dplanck_df* internal function:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
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
Numeric rayjean(const Numeric& f, const Numeric& tb) try {
  if (tb <= 0) throw "Non-positive temperature";
  if (f < 0) throw "Negative frequency";

  static const Numeric a =
      SPEED_OF_LIGHT * SPEED_OF_LIGHT / (2 * BOLTZMAN_CONST);

  return (f * f) / (a * tb);
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *rayjean* internal function:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
}
