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
 * \file   arts_constants.h
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

#include <numbers>

#include "arts_constexpr_math.h"
#include "matpack.h"

/** Namespace containing several constants, physical and mathematical **/
namespace Constant {
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
inline constexpr Numeric pi = std::numbers::pi;

/** Inverse of pi */
inline constexpr Numeric inv_pi = std::numbers::inv_pi;

/** Two times pi **/
inline constexpr Numeric two_pi = 2 * pi;

/** Inverse of two pi **/
inline constexpr Numeric inv_two_pi = 0.5 * inv_pi;

/** Square root of pi */
inline constexpr Numeric sqrt_pi = 1.0/std::numbers::inv_sqrtpi;

/** Inverse of the square root of pi */
inline constexpr Numeric inv_sqrt_pi = std::numbers::inv_sqrtpi;

/** Euler's number */
inline constexpr Numeric euler = std::numbers::e;

/** Inverse of Euler's number */
inline constexpr Numeric inv_euler = 1.0 / euler;

/** Ten's logarithm of Euler's number */
inline constexpr Numeric log10_euler = std::numbers::log10e;

/** Natural logarithm of 10 */
inline constexpr Numeric ln_10 = std::numbers::ln10;

/** Square root of 2 */
inline constexpr Numeric sqrt_2 = std::numbers::sqrt2;

/** Inverse of the square root of 2 */
inline constexpr Numeric inv_sqrt_2 = 1.0 / sqrt_2;

/** Natural logarithm of 2 */
inline constexpr Numeric ln_2 = std::numbers::ln2;

/** Inverse of the natural logarithm of 2 */
inline constexpr Numeric inv_ln_2 = 1.0 / ln_2;

/** Square root of natural logarithm of 2 */
inline constexpr Numeric sqrt_ln_2 =
    0.8325546111576977563531646448952010476305888522644407291668291172340794351973;

/** Inverse of the square root of the natural logarithm of 2 */
inline constexpr Numeric inv_sqrt_ln_2 =
    1.201122408786449794857803286095221722566764028068699423868879896733837175546;

/** Cesium-133 Unperturbed ground-state hyperfine transition frequency [Hz]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
   **/
inline constexpr Numeric Delta_nu_Cs = 9192631770;

/** Speed of light [m/s]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
    **/
inline constexpr Numeric speed_of_light = 299792458;

/** Speed of light convenience name [m/s] **/
inline constexpr Numeric c = speed_of_light;

/** Planck constant [J s]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
    **/
inline constexpr Numeric planck_constant = 6.62607015e-34;

/** Planck constant convenience name [J s] **/
inline constexpr Numeric h = planck_constant;

/** Reduced planck constant [J s] **/
inline constexpr Numeric reduced_planck_constant = h * inv_two_pi;

/** Reduced planck constant convenience name [J s] **/
inline constexpr Numeric h_bar = reduced_planck_constant;

/** Elementary charge [C]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
    **/
inline constexpr Numeric elementary_charge = 1.602176634e-19;

/** Elementary charge convenience name [C] **/
inline constexpr Numeric e = elementary_charge;

/** Boltzmann constant [J/K]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
    **/
inline constexpr Numeric boltzmann_constant = 1.380649e-23;

/** Boltzmann constant convenience name [J/K] **/
inline constexpr Numeric k = boltzmann_constant;

/** Avogadro constant [1/mol]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
    **/
inline constexpr Numeric avogadro_constant = 6.02214076e23;

/** Avogadro constant convenience name [1/mol] **/
inline constexpr Numeric NA = avogadro_constant;

/** Luminous efficacy of monochromatic 540 THz radiation [lm / W]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units
   Date: 2019-04-01
    **/
inline constexpr Numeric K_cd = 683;

/** Fine structure constant [-]
   From: https://physics.nist.gov/cgi-bin/cuu/Value?alph
   Date: 2019-06-18
   Reported error: (11)
    **/
inline constexpr Numeric fine_structure_constant = 7.2973525693e-3;

/** Fine structure constant convenience name [-] **/
inline constexpr Numeric alpha = fine_structure_constant;

/** Rydberg constant [1/m]
    From: https://physics.nist.gov/cgi-bin/cuu/Value?ryd
    Date: 2016-06-18
    Reported error: (21)
    **/
inline constexpr Numeric rydberg_constant = 10973731.568160;

/** Rydberg constant convenience name [1/m] **/
inline constexpr Numeric R_inf = rydberg_constant;

/** Magnetic constant [H/m] **/
inline constexpr Numeric magnetic_constant = 2 * h * alpha / (c * Math::pow2(e));

/** Magnetic constant convenience name [H/m] **/
inline constexpr Numeric mu_0 = magnetic_constant;

/** Vacuum permittivity [F/m] **/
inline constexpr Numeric vacuum_permittivity = Math::pow2(e) / (2 * h * c * alpha);

/** Vacuum permittivity convenience name [F/m] **/
inline constexpr Numeric epsilon_0 = vacuum_permittivity;

/** Mass of resting electron [kg] **/
inline constexpr Numeric electron_mass = 2 * h * R_inf / (c * Math::pow2(alpha));

/** Mass of resting electron convenience name [kg] **/
inline constexpr Numeric m_e = electron_mass;

/** Unified atomic mass unit [kg] 
    From: https://physics.nist.gov/cgi-bin/cuu/Value?ukg
    Date: 2020-02-18
    Reported error: (50)
    **/
inline constexpr Numeric unified_atomic_mass_unit = 1.66053906660e-27;

/** Unified atomic mass unit convenience name [kg] **/
inline constexpr Numeric m_u = unified_atomic_mass_unit;

/** Mass ratio of electrons to protons [-]
    From: https://physics.nist.gov/cgi-bin/cuu/Value?mpsme
    Date: 2020-01-08
    Reported error: (11)
    **/
inline constexpr Numeric mass_ratio_electrons_per_proton = 1'836.152'673'43;

/** Mass of a proton [kg] */
inline constexpr Numeric proton_mass =
    electron_mass * mass_ratio_electrons_per_proton;

/** Mass ratio of electrons to protons [-]
    From: https://physics.nist.gov/cgi-bin/cuu/Value?mnsme
    Date: 2020-01-08
    Reported error: (89)
    **/
inline constexpr Numeric mass_ratio_electrons_per_neutron = 1'838.683'661'73;

/** Mass of a neutron [kg] */
inline constexpr Numeric neutron_mass =
    electron_mass * mass_ratio_electrons_per_neutron;

/** Bohr magneton [J/T] **/
inline constexpr Numeric bohr_magneton = e * h_bar / (2 * m_e);

/** Ideal gas constant [J/mol K] **/
inline constexpr Numeric ideal_gas_constant = k * NA;

/** Ideal gas constant convenience name [J/mol K] **/
inline constexpr Numeric R = ideal_gas_constant;

/** Doppler broadening constant squared [kg/T]^2 **/
inline constexpr Numeric doppler_broadening_const_squared = 2'000 * R / Math::pow2(c);

/** One degree in radians */
inline constexpr Numeric one_degree_in_radians = pi / 180;

/** Stefan-Boltzmann constant [W/(K^4*m^2)] */
inline constexpr Numeric stefan_boltzmann_constant =
    Math::pow2(pi) * Math::pow4(k) / (Math::pow3(h_bar) * Math::pow2(c)) / 60.;

/** Stefan-Boltzmann constant convenience name [W/(K^4*m^2)] */
inline constexpr Numeric sigma = stefan_boltzmann_constant;

/** Global constant, Density of water ice at 0C [kg/m3]
    source: http://en.wikipedia.org/wiki/Ice
    \author Jana Mendrok 
    \date   2014-10-15
*/
inline constexpr Numeric density_of_ice_at_0c = 0.9167e3;

/** Global constant, Density of liquid water +4C [kg/m3]
    source: http://en.wikipedia.org/wiki/Water
    \author Jana Mendrok 
    \date   2014-10-15
*/
inline constexpr Numeric denity_of_water_at_4c = 1e3;

/** Global constant, Planck temperature for cosmic background radiation [K]
    \author Patrick Eriksson 
    \date   08.04.2000
*/
inline constexpr Numeric cosmic_microwave_background_temperature = 2.735;

/** Global constant, the radius of the Earth [m]
    \author Patrick Eriksson 
    \date   08.04.2000
*/
inline constexpr Numeric earth_radius = 6.3781e6;

/** Global constant, Temperature in Celsius of 0 Kelvin.
    \author Axel von Engeln
    \date   2000-12-19
*/
inline constexpr Numeric temperature_at_0c = 273.15;
};  // namespace Constant

#endif
