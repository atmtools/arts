/* Copyright (C) 2019
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/*!
 * \file   constants.h
 * \brief  Constants of physical expressions as constexpr
 * 
 * Following May 2019 reported SI unit updates, several previously
 * important constant are proper constants and not derived constants.
 * 
 * This file contains the intended May 2019 constants update meaning
 * that the Planck constant, Boltzmann constant, elementary charge,
 * Avogadro number, the definition of a luminosity, the frequency of
 * Cesium, and the speed of light have been given constant values.
 * 
 * All other values are derived.  To be internally consistent in ARTS,
 * the adaptation of these rules below take the approach that only the
 * fine structure constant and the Rydberg constant are 'measurable'
 * quantities.  All other quantities are derived by their formal
 * relationships.  As an example, the resting mass of an electron is
 * defined as 2 h R_inf / c alpha^2.  This means that they should be
 * constant in their internal relationships to the level permitted by
 * the select floating point accuracy.
 * 
 * Note 1: The constants below contains both named and a few convenience
 * variables
 * 
 * Note 2: Derived constants of convenience, such as the Doppler broadening
 * constant, should also go into here if they have a clear name.
 * 
 * Note 3: The PrintPhysicalConstants convenience function is part of the 
 * global ARTS namespace.  Please update this if and when new constants
 * are added to this namespace.
 * 
 * \author Richard Larsson
 * \date   2019-04-01
 */

#ifndef CONSTANTS_IN_ARTS_H
#define CONSTANTS_IN_ARTS_H

#include <cmath>
#include "matpack.h"

/** Namespace containing several constants, physical and mathematical **/
namespace Constant {
/** power of two */
template <class T>
constexpr T pow2(T x) {
  return x * x;
}

/** power of three */
template <class T>
constexpr T pow3(T x) {
  return pow2(x) * x;
}

/** power of four */
template <class T>
constexpr T pow4(T x) {
  return pow2(pow2(x));
}

/** The following mathematical constants are generated in python Decimal package by the code:
    * 
    import decimal as d
    pi_str = '3.141592653589793238462643383279502884197169399375105820974944592307816406286'
    d.getcontext().prec = len(pi_str) - 1
    pi = d.Decimal(pi_str)
    one = d.Decimal('1')
    two = d.Decimal('2')
    ten = d.Decimal('10')
    print('pi =', pi)
    print('inv_pi =', 1/pi)
    print('two_pi =', two*pi)
    print('inv_two_pi =', 1/(two*pi))
    print('sqrt_pi =', pi.sqrt())
    print('inv_sqrt_pi =', 1/pi.sqrt())
    print('euler =', one.exp())
    print('inv_euler =', 1/one.exp())
    print('log10_euler =', one.exp().log10())
    print('ln_10 =', ten.log())
    print('sqrt_2 =', two.sqrt())
    print('inv_sqrt_2 =', 1/two.sqrt())
    print('ln_2 = ', two.ln())
    print('inv_ln_2 =', 1/two.ln())
    print('sqrt_ln_2 =', two.ln().sqrt())
    print('inv_sqrt_ln_2 =', 1/two.ln().sqrt())
    *
  To improve the numerical accuracy further, insert larger pi string */

/** Pi, related to circles
   From: //www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html 2019-04-01
   **/
static constexpr Numeric pi =
    3.141592653589793238462643383279502884197169399375105820974944592307816406286;

/** Inverse of pi */
static constexpr Numeric inv_pi =
    0.3183098861837906715377675267450287240689192914809128974953346881177935952685;

/** Two times pi **/
static constexpr Numeric two_pi =
    6.283185307179586476925286766559005768394338798750211641949889184615632812572;

/** Inverse of two pi **/
static constexpr Numeric inv_two_pi =
    0.1591549430918953357688837633725143620344596457404564487476673440588967976342;

/** Square root of pi */
static constexpr Numeric sqrt_pi =
    1.772453850905516027298167483341145182797549456122387128213807789852911284591;

/** Inverse of the square root of pi */
static constexpr Numeric inv_sqrt_pi =
    0.5641895835477562869480794515607725858440506293289988568440857217106424684415;

/** Euler's number */
static constexpr Numeric euler =
    2.718281828459045235360287471352662497757247093699959574966967627724076630354;

/** Inverse of Euler's number */
static constexpr Numeric inv_euler =
    0.3678794411714423215955237701614608674458111310317678345078368016974614957448;

/** Ten's logarithm of Euler's number */
static constexpr Numeric log10_euler =
    0.4342944819032518276511289189166050822943970058036665661144537831658646492089;

/** Natural logarithm of 10 */
static constexpr Numeric ln_10 =
    2.302585092994045684017991454684364207601101488628772976033327900967572609677;

/** Square root of 2 */
static constexpr Numeric sqrt_2 =
    1.414213562373095048801688724209698078569671875376948073176679737990732478462;

/** Inverse of the square root of 2 */
static constexpr Numeric inv_sqrt_2 =
    0.7071067811865475244008443621048490392848359376884740365883398689953662392311;

/** Natural logarithm of 2 */
static constexpr Numeric ln_2 =
    0.6931471805599453094172321214581765680755001343602552541206800094933936219697;

/** Inverse of the natural logarithm of 2 */
static constexpr Numeric inv_ln_2 =
    1.442695040888963407359924681001892137426645954152985934135449406931109219181;

/** Square root of natural logarithm of 2 */
static constexpr Numeric sqrt_ln_2 =
    0.8325546111576977563531646448952010476305888522644407291668291172340794351973;

/** Inverse of the square root of the natural logarithm of 2 */
static constexpr Numeric inv_sqrt_ln_2 =
    1.201122408786449794857803286095221722566764028068699423868879896733837175546;

/** Cesium-133 Unperturbed ground-state hyperfine transition frequency [Hz]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
   **/
static constexpr Numeric Delta_nu_Cs = 9192631770;

/** Speed of light [m/s]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
    **/
static constexpr Numeric speed_of_light = 299792458;

/** Speed of light convenience name [m/s] **/
static constexpr Numeric c = speed_of_light;

/** Planck constant [J s]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
    **/
static constexpr Numeric planck_constant = 6.62607015e-34;

/** Planck constant convenience name [J s] **/
static constexpr Numeric h = planck_constant;

/** Reduced planck constant [J s] **/
static constexpr Numeric reduced_planck_constant = h * inv_two_pi;

/** Reduced planck constant convenience name [J s] **/
static constexpr Numeric h_bar = reduced_planck_constant;

/** Elementary charge [C]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
    **/
static constexpr Numeric elementary_charge = 1.602176634e-19;

/** Elementary charge convenience name [C] **/
static constexpr Numeric e = elementary_charge;

/** Boltzmann constant [J/K]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
    **/
static constexpr Numeric boltzmann_constant = 1.380649e-23;

/** Boltzmann constant convenience name [J/K] **/
static constexpr Numeric k = boltzmann_constant;

/** Avogadro constant [1/mol]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
    **/
static constexpr Numeric avogadro_constant = 6.02214076e23;

/** Avogadro constant convenience name [1/mol] **/
static constexpr Numeric NA = avogadro_constant;

/** Luminous efficacy of monochromatic 540 THz radiation [lm / W]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
    **/
static constexpr Numeric K_cd = 683;

/** Fine structure constant [-]
   From: https://physics.nist.gov/cgi-bin/cuu/Value?alph
   Date: 2019-06-18
   Reported error: (11)
    **/
static constexpr Numeric fine_structure_constant = 7.2973525693e-3;

/** Fine structure constant convenience name [-] **/
static constexpr Numeric alpha = fine_structure_constant;

/** Rydberg constant [1/m]
    From: https://physics.nist.gov/cgi-bin/cuu/Value?ryd
    Date: 2016-06-18
    Reported error: (21)
    **/
static constexpr Numeric rydberg_constant = 10973731.568160;

/** Rydberg constant convenience name [1/m] **/
static constexpr Numeric R_inf = rydberg_constant;

/** Magnetic constant [H/m] **/
static constexpr Numeric magnetic_constant = 2 * h * alpha / (c * pow2(e));

/** Magnetic constant convenience name [H/m] **/
static constexpr Numeric mu_0 = magnetic_constant;

/** Vacuum permittivity [F/m] **/
static constexpr Numeric vacuum_permittivity = pow2(e) / (2 * h * c * alpha);

/** Vacuum permittivity convenience name [F/m] **/
static constexpr Numeric epsilon_0 = vacuum_permittivity;

/** Mass of resting electron [kg] **/
static constexpr Numeric electron_mass = 2 * h * R_inf / (c * pow2(alpha));

/** Mass of resting electron convenience name [kg] **/
static constexpr Numeric m_e = electron_mass;

/** Mass ratio of electrons to protons [-]
    From: https://physics.nist.gov/cgi-bin/cuu/Value?mpsme
    Date: 2020-01-08
    Reported error: (11)
    **/
static constexpr Numeric mass_ratio_electrons_per_proton = 1'836.152'673'43;

/** Mass of a proton [kg] */
static constexpr Numeric proton_mass = electron_mass * mass_ratio_electrons_per_proton;

/** Mass ratio of electrons to protons [-]
    From: https://physics.nist.gov/cgi-bin/cuu/Value?mnsme
    Date: 2020-01-08
    Reported error: (89)
    **/
static constexpr Numeric mass_ratio_electrons_per_neutron = 1'838.683'661'73;

/** Mass of a neutron [kg] */
static constexpr Numeric neutron_mass = electron_mass * mass_ratio_electrons_per_neutron;

/** Bohr magneton [J/T] **/
static constexpr Numeric bohr_magneton = e * h_bar / (2 * m_e);

/** Ideal gas constant [J/mol K] **/
static constexpr Numeric ideal_gas_constant = k * NA;

/** Ideal gas constant convenience name [J/mol K] **/
static constexpr Numeric R = ideal_gas_constant;

/** Doppler broadening constant squared [kg/T]^2 **/
static constexpr Numeric doppler_broadening_const_squared = 2000 * R / pow2(c);
};  // namespace Constant

/** Namespace containing several practical unit conversions, physical and mathematical **/
namespace Conversion {
using namespace Constant;

/** Conversion constant degrees to radians and back.  Use conversion formulae instead of pure constant if possible. NOTE:  No constexpr cos etal in ARTS yet. **/
static constexpr Numeric DEG2RAD = pi / 180;
static constexpr Numeric RAD2DEG = 1 / DEG2RAD;
/** Converts degrees to radians  */
template <class T>
constexpr Numeric deg2rad(T x) {
  return x * DEG2RAD;
}

/** Converts radians to degrees  */
template <class T>
constexpr Numeric rad2deg(T x) {
  return x * RAD2DEG;
}

/** Returns the cosine of the deg2rad of the input  */
template <class T>
Numeric cosd(T x) {
  return std::cos(deg2rad(x));
}

/** Returns the sine of the deg2rad of the input  */
template <class T>
Numeric sind(T x) {
  return std::sin(deg2rad(x));
}

/** Returns the tangent of the deg2rad of the input  */
template <class T>
Numeric tand(T x) {
  return std::tan(deg2rad(x));
}

/** Returns rad2deg of the arc-cosine of the input  */
template <class T>
Numeric acosd(T x) {
  return rad2deg(std::acos(x));
}

/** Returns rad2deg of the arc-sine of the input  */
template <class T>
Numeric asind(T x) {
  return rad2deg(std::asin(x));
}

/** Returns rad2deg of the arc-tangent of the input  */
template <class T>
Numeric atand(T x) {
  return rad2deg(std::atan(x));
}

/** Returns rad2deg of the arc-tangent of inputs #T1/#T2  */
template <class T1, class T2>
Numeric atan2d(T1 y, T2 x) {
  return rad2deg(std::atan2(y, x));
}

/** Conversion constant Kayser wavenumber to frequency and back.  Use conversion formulae instead of pure constant if possible. **/
static constexpr Numeric KAYCM2FREQ = 100 * c;
static constexpr Numeric FREQ2KAYCM = 1 / KAYCM2FREQ;
template <class T>
constexpr Numeric kaycm2freq(T x) {
  return x * KAYCM2FREQ;
}
template <class T>
constexpr Numeric freq2kaycm(T x) {
  return x * FREQ2KAYCM;
}

/** Conversion constant Angular wavenumber to frequency and back.  Use conversion formulae instead of pure constant if possible. **/
static constexpr Numeric ANGCM2FREQ = KAYCM2FREQ * inv_two_pi;
static constexpr Numeric FREQ2ANGCM = 1 / ANGCM2FREQ;
template <class T>
constexpr Numeric angcm2freq(T x) {
  return x * ANGCM2FREQ;
}
template <class T>
constexpr Numeric freq2angcm(T x) {
  return x * FREQ2ANGCM;
}

/** Conversion constant Angular frequency to frequency and back **/
template <class T>
constexpr Numeric angfreq2freq(T x) {
  return x * inv_two_pi;
}
template <class T>
constexpr Numeric freq2angfreq(T x) {
  return x * two_pi;
}

/** Conversion wavelength to frequency and back. **/
template <class T>
constexpr Numeric wavelen2freq(T x) {
  return c / x;
}
template <class T>
constexpr Numeric freq2wavelen(T x) {
  return c / x;
}

/** Conversion constant 1 atmosphere to 1 Pascal and back.  Use conversion formulae instead of pure constant if possible. **/
static constexpr Numeric ATM2PA = 101325;
static constexpr Numeric PA2ATM = 1 / ATM2PA;
template <class T>
constexpr Numeric atm2pa(T x) {
  return x * ATM2PA;
}
template <class T>
constexpr Numeric pa2atm(T x) {
  return x * PA2ATM;
}

/** Conversion constant 1 torr to 1 Pascal and back.  Use conversion formulae instead of pure constant if possible. **/
static constexpr Numeric TORR2PA = ATM2PA / 760;
static constexpr Numeric PA2TORR = 1 / TORR2PA;
template <class T>
constexpr Numeric torr2pa(T x) {
  return x * TORR2PA;
}
template <class T>
constexpr Numeric pa2torr(T x) {
  return x * PA2TORR;
}

/** Conversion constant Celsius to Kelvin and back.  Use conversion formulae instead of pure constant if possible. **/
static constexpr Numeric CEL2KEL = 273.15;
template <class T>
constexpr Numeric celsius2kelvin(T x) {
  return x + CEL2KEL;
}
template <class T>
constexpr Numeric kelvin2celsius(T x) {
  return x - CEL2KEL;
}

/** Conversion from cm-1 per molecule per cm^2 to Hz per molecule per m^2 **/
static constexpr Numeric HITRAN2ARTS_LS = KAYCM2FREQ * 1e-4;
static constexpr Numeric ARTS2HITRAN_LS = 1 / HITRAN2ARTS_LS;
template <class T>
constexpr T hitran2arts_linestrength(T x) {
  return x * HITRAN2ARTS_LS;
}
template <class T>
constexpr T arts2hitran_linestrength(T x) {
  return x * ARTS2HITRAN_LS;
}

/** Conversion from cm-1 per atmosphere to Hz per Pascal **/
static constexpr Numeric HITRAN2ARTS_GAMMA = KAYCM2FREQ / ATM2PA;
static constexpr Numeric ARTS2HITRAN_GAMMA = 1 / HITRAN2ARTS_GAMMA;
template <class T>
constexpr T hitran2arts_broadening(T x) {
  return x * HITRAN2ARTS_GAMMA;
}
template <class T>
constexpr T arts2hitran_broadening(T x) {
  return x * ARTS2HITRAN_GAMMA;
}

/** Conversion from cm-1 to Joule **/
static constexpr Numeric HITRAN2ARTS_ENERGY = h * KAYCM2FREQ;
static constexpr Numeric ARTS2HITRAN_ENERGY = 1 / HITRAN2ARTS_ENERGY;
template <class T>
constexpr T hitran2arts_energy(T x) {
  return x * HITRAN2ARTS_ENERGY;
}
template <class T>
constexpr T arts2hitran_energy(T x) {
  return x * ARTS2HITRAN_ENERGY;
}
};  // namespace Conversion

#endif
