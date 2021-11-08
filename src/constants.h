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

#include "enums.h"
#include "matpack.h"

/** Namespace containing several constants, physical and mathematical **/
namespace Constant {
/** power of two */
template <class T>
constexpr auto pow2(T x) noexcept -> decltype(x * x) {
  return x * x;
}

/** power of three */
template <class T>
constexpr auto pow3(T x) noexcept -> decltype(pow2(x) * x) {
  return pow2(x) * x;
}

/** power of four */
template <class T>
constexpr auto pow4(T x) noexcept -> decltype(pow2(pow2(x))) {
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
    print('ln_10 =', ten.ln())
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

/** Unified atomic mass unit [kg] 
    From: https://physics.nist.gov/cgi-bin/cuu/Value?ukg
    Date: 2020-02-18
    Reported error: (50)
    **/
static constexpr Numeric unified_atomic_mass_unit = 1.66053906660e-27;

/** Unified atomic mass unit convenience name [kg] **/
static constexpr Numeric m_u = unified_atomic_mass_unit;

/** Mass ratio of electrons to protons [-]
    From: https://physics.nist.gov/cgi-bin/cuu/Value?mpsme
    Date: 2020-01-08
    Reported error: (11)
    **/
static constexpr Numeric mass_ratio_electrons_per_proton = 1'836.152'673'43;

/** Mass of a proton [kg] */
static constexpr Numeric proton_mass =
    electron_mass * mass_ratio_electrons_per_proton;

/** Mass ratio of electrons to protons [-]
    From: https://physics.nist.gov/cgi-bin/cuu/Value?mnsme
    Date: 2020-01-08
    Reported error: (89)
    **/
static constexpr Numeric mass_ratio_electrons_per_neutron = 1'838.683'661'73;

/** Mass of a neutron [kg] */
static constexpr Numeric neutron_mass =
    electron_mass * mass_ratio_electrons_per_neutron;

/** Bohr magneton [J/T] **/
static constexpr Numeric bohr_magneton = e * h_bar / (2 * m_e);

/** Ideal gas constant [J/mol K] **/
static constexpr Numeric ideal_gas_constant = k * NA;

/** Ideal gas constant convenience name [J/mol K] **/
static constexpr Numeric R = ideal_gas_constant;

/** Doppler broadening constant squared [kg/T]^2 **/
static constexpr Numeric doppler_broadening_const_squared = 2'000 * R / pow2(c);

/** One degree in radians */
static constexpr Numeric one_degree_in_radians = pi / 180;
};  // namespace Constant

/** Namespace containing several practical unit conversions, physical and mathematical **/
namespace Conversion {
using namespace Constant;

/** Converts degrees to radians */
template <class T>
constexpr auto deg2rad(T x) noexcept -> decltype(x * one_degree_in_radians) {
  return x * one_degree_in_radians;
}

/** Converts radians to degrees */
template <class T>
constexpr auto rad2deg(T x) noexcept -> decltype(x / one_degree_in_radians) {
  return x / one_degree_in_radians;
}

/** Returns the cosine of the deg2rad of the input */
template <class T>
auto cosd(T x) noexcept -> decltype(std::cos(deg2rad(x))) {
  return std::cos(deg2rad(x));
}

/** Returns the sine of the deg2rad of the input */
template <class T>
auto sind(T x) noexcept -> decltype(std::sin(deg2rad(x))) {
  return std::sin(deg2rad(x));
}

/** Returns the tangent of the deg2rad of the input */
template <class T>
auto tand(T x) noexcept -> decltype(std::tan(deg2rad(x))) {
  return std::tan(deg2rad(x));
}

/** Returns rad2deg of the arc-cosine of the input */
template <class T>
auto acosd(T x) noexcept -> decltype(rad2deg(std::acos(x))) {
  return rad2deg(std::acos(x));
}

/** Returns rad2deg of the arc-sine of the input */
template <class T>
auto asind(T x) noexcept -> decltype(rad2deg(std::asin(x))) {
  return rad2deg(std::asin(x));
}

/** Returns rad2deg of the arc-tangent of the input */
template <class T>
auto atand(T x) noexcept -> decltype(rad2deg(std::atan(x))) {
  return rad2deg(std::atan(x));
}

/** Returns rad2deg of the arc-tangent of inputs #T1/#T2  */
template <class T1, class T2>
auto atan2d(T1 y, T2 x) noexcept -> decltype(rad2deg(std::atan2(y, x))) {
  return rad2deg(std::atan2(y, x));
}

/** Conversion from Kayser wavenumber to Hz */
template <class T>
constexpr auto kaycm2freq(T x) noexcept -> decltype(x * (100 * c)) {
  return x * (100 * c);
}

/** Conversion from Hz to Kayser wavenumber */
template <class T>
constexpr auto freq2kaycm(T x) noexcept -> decltype(x / (100 * c)) {
  return x / (100 * c);
}

/** Conversion from Angular wavenumber to Hz */
template <class T>
constexpr auto angcm2freq(T x) noexcept -> decltype(kaycm2freq(inv_two_pi)) {
  return x * kaycm2freq(inv_two_pi);
}

/** Conversion from Hz to Angular wavenumber */
template <class T>
constexpr auto freq2angcm(T x) noexcept
    -> decltype(x / kaycm2freq(inv_two_pi)) {
  return x / kaycm2freq(inv_two_pi);
}

/** Conversion from Angular Hz to Hz */
template <class T>
constexpr auto angfreq2freq(T x) noexcept -> decltype(x * inv_two_pi) {
  return x * inv_two_pi;
}

/** Conversion from Hz to Angular Hz */
template <class T>
constexpr auto freq2angfreq(T x) noexcept -> decltype(x * two_pi) {
  return x * two_pi;
}

/** Conversion from wavelength to Hz */
template <class T>
constexpr auto wavelen2freq(T x) noexcept -> decltype(c / x) {
  return c / x;
}

/** Conversion from Hz to wavelength */
template <class T>
constexpr auto freq2wavelen(T x) noexcept -> decltype(c / x) {
  return c / x;
}

/** Conversion from wavelength to Hz */
template <class T>
constexpr auto hz2ghz(T x) noexcept -> decltype(x * 1e-9) {
  return x * 1e-9;
}

/** Conversion from Hz to wavelength */
template <class T>
constexpr auto ghz2hz(T x) noexcept -> decltype(x * 1e9) {
  return x * 1e9;
}

/** Conversion from Atm to Pa */
template <class T>
constexpr auto atm2pa(T x) noexcept -> decltype(x * 101'325.0) {
  return x * 101'325.0;
}

/** Conversion from Pa to Atm */
template <class T>
constexpr auto pa2atm(T x) noexcept -> decltype(x / 101'325.0) {
  return x / 101'325.0;
}

/** Conversion from bar to Pa */
template <class T>
constexpr auto bar2pa(T x) noexcept -> decltype(x * 1e5) {
  return x * 1e5;
}

/** Conversion from Pa to bar */
template <class T>
constexpr auto pa2bar(T x) noexcept -> decltype(x * 1e-5) {
  return x * 1e-5;
}

/** Conversion from Torr to Pa */
template <class T>
constexpr auto torr2pa(T x) noexcept -> decltype(x * atm2pa(1.0 / 760.0)) {
  return x * atm2pa(1.0 / 760.0);
}

/** Conversion from Pa to Torr */
template <class T>
constexpr auto pa2torr(T x) noexcept -> decltype(x / atm2pa(1.0 / 760.0)) {
  return x / atm2pa(1.0 / 760.0);
}

/** Conversion from MHz/Torr to Hz/Pa */
template <class T>
constexpr auto mhz_per_torr2hz_per_pa(T x) noexcept
    -> decltype(x * pa2torr(1e6)) {
  return x * pa2torr(1e6);
}

/** Conversion from C to K */
template <class T>
constexpr auto celsius2kelvin(T x) noexcept -> decltype(x + 273.15) {
  return x + 273.15;
}

/** Conversion from K to C */
template <class T>
constexpr auto kelvin2celsius(T x) noexcept -> decltype(x - 273.15) {
  return x - 273.15;
}

/** Conversion from cm-1 per molecule per cm^2 to Hz per molecule per m^2 **/
template <class T>
constexpr auto kaycm_per_cmsquared2hz_per_msquared(T x) noexcept
    -> decltype(x * kaycm2freq(1e-4)) {
  return x * kaycm2freq(1e-4);
}

/** Conversion from Hz per molecule per m^2 to cm-1 per molecule per cm^2 **/
template <class T>
constexpr auto hz_per_msquared2kaycm_per_cmsquared(T x) noexcept
    -> decltype(x * freq2kaycm(1e4)) {
  return x * freq2kaycm(1e4);
}

/** Conversion from cm-1 per atmosphere to Hz per Pascal **/
template <class T>
constexpr auto kaycm_per_atm2hz_per_pa(T x) noexcept
    -> decltype(x * kaycm2freq(pa2atm(1))) {
  return x * kaycm2freq(pa2atm(1));
}

/** Conversion from Hz per Pascal to cm-1 per atmosphere **/
template <class T>
constexpr auto hz_per_pa2kaycm_per_atm(T x) noexcept
    -> decltype(x * freq2kaycm(atm2pa(1))) {
  return x * freq2kaycm(atm2pa(1));
}

/** Conversion from cm-1 to Joule **/
template <class T>
constexpr auto kaycm2joule(T x) noexcept -> decltype(x * kaycm2freq(h)) {
  return x * kaycm2freq(h);
}

/** Conversion from MHz to Joule **/
template <class T>
constexpr auto hz2joule(T x) noexcept -> decltype(x * h) {
  return x * h;
}

/** Conversion from MHz to Joule **/
template <class T>
constexpr auto mhz2joule(T x) noexcept -> decltype(hz2joule(x) * 1e6) {
  return hz2joule(x) * 1e6;
}

/** Conversion from Kelvin to Joule **/
template <class T>
constexpr auto kelvin2joule(T x) noexcept -> decltype(x * k) {
  return x * k;
}

/** Conversion from Hz to Joule **/
template <class T>
constexpr auto joule2hz(T x) noexcept -> decltype(x / h) {
  return x / h;
}

/** Conversion from Joule to cm-1 **/
template <class T>
constexpr auto joule2kaycm(T x) noexcept -> decltype(x / kaycm2freq(h)) {
  return x / kaycm2freq(h);
}

/** Conversion from Å to meter **/
template <class T>
constexpr auto angstrom2meter(T x) noexcept -> decltype(x * 1e-10) {
  return x * 1e-10;
}

/** Conversion from meter to Å **/
template <class T>
constexpr auto meter2angstrom(T x) noexcept -> decltype(x * 1e10) {
  return x * 1e10;
}
};  // namespace Conversion

namespace Options {
/** Keep time options available to switch over them */
ENUMCLASS(
    TimeStep, char, hour, hours, h, minute, minutes, min, second, seconds, s)

/** Possible line shape coefficients */
ENUMCLASS(LineShapeCoeff, char, X0, X1, X2, X3)

/** Possible Jacobian for Wind and Magnetic field */
ENUMCLASS(WindMagJacobian, char, u, v, w, strength)

/** Possible Jacobian for basic line parameters */
ENUMCLASS(BasicCatParamJacobian, char, LineStrength, LineCenter)

/** Possible HITRAN types of catalogs */
ENUMCLASS(
    HitranType,
    char,
    Pre2004,         // 2004 version changed the .par-length
    From2004To2012,  // New par length but old isotopologues order
    Post2012,        // 2012 version changed the order of isotopologues
    Online  // Onine expects a modern .par line followed by Upper then Lower quantum numbers
)

/** Possible AddLines Speedups */
ENUMCLASS(LblSpeedup, char, None, QuadraticIndependent, LinearIndependent)

ENUMCLASS(SortingOption, char, ByFrequency, ByEinstein)
}  // namespace Options

#endif
