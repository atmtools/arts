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


/** Namespace containing several constants, physical and mathematical **/
namespace Constant {
  template <class T> constexpr T pow2(T x) {return x*x;}
  template <class T> constexpr T pow3(T x) {return x*x*x;}
  template <class T> constexpr Numeric inv(T x) {return 1e0/x;}
  template <class T> constexpr Numeric inv_pow2(T x) {return inv(pow2(x));}
  template <class T> constexpr Numeric inv_pow3(T x) {return inv(pow3(x));}
  
  /** The following mathematical constants are generated in python Decimal package by the code:
    * 
    import decimal as d
    pi_str = '3.141592653589793238462643383279502884197169399375105820974944592307816406286'
    d.getcontext().prec = len(pi_str) - 1
    pi = d.Decimal(pi_str)
    two = d.Decimal('2')
    print('pi =', pi)
    print('inv_pi =', 1/pi)
    print('two_pi =', two*pi)
    print('inv_two_pi =', 1/two_pi)
    print('sqrt_pi =', pi.sqrt())
    print('inv_sqrt_pi =', 1/pi.sqrt())
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
  constexpr Numeric pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
  
  /** Inverse of pi */
  constexpr Numeric inv_pi = 0.3183098861837906715377675267450287240689192914809128974953346881177935952685;
  
  /** Two times pi **/
  constexpr Numeric two_pi = 6.283185307179586476925286766559005768394338798750211641949889184615632812572;
  
  /** Inverse of two pi **/
  constexpr Numeric inv_two_pi = 0.1591549430918953357688837633725143620344596457404564487476673440588967976342;
  
  /** Square root of pi */
  constexpr Numeric sqrt_pi = 1.772453850905516027298167483341145182797549456122387128213807789852911284591;
  
  /** Inverse of the square root of pi */
  constexpr Numeric inv_sqrt_pi = 0.5641895835477562869480794515607725858440506293289988568440857217106424684415;
  
  /** Square root of 2 */
  constexpr Numeric sqrt_2 = 1.414213562373095048801688724209698078569671875376948073176679737990732478462;
  
  /** Inverse of the square root of 2 */
  constexpr Numeric inv_sqrt_2 = 0.7071067811865475244008443621048490392848359376884740365883398689953662392311;
  
  /** Natural logarithm of 2 */
  constexpr Numeric ln_2 = 0.6931471805599453094172321214581765680755001343602552541206800094933936219697;
  
  /** Inverse of the natural logarithm of 2 */
  constexpr Numeric inv_ln_2 = 1.442695040888963407359924681001892137426645954152985934135449406931109219181;
  
  /** Square root of natural logarithm of 2 */
  constexpr Numeric sqrt_ln_2 = 0.8325546111576977563531646448952010476305888522644407291668291172340794351973;
  
  /** Inverse of the square root of the natural logarithm of 2 */
  constexpr Numeric inv_sqrt_ln_2 = 1.201122408786449794857803286095221722566764028068699423868879896733837175546;
  
  /** Cesium-133 Unperturbed ground-state hyperfine transition frequency [Hz]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units 2019-04-01
   **/
  constexpr Numeric Delta_nu_Cs = 9192631770;
  
  /** Speed of light [m/s]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units 2019-04-01
    **/
  constexpr Numeric speed_of_light = 299792458;
  
  /** Speed of light convenience name [m/s] **/
  constexpr Numeric c = speed_of_light;
  
  /** Planck constant [J s]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units 2019-04-01
    **/
  constexpr Numeric planck_constant = 6.62607015e-34;
  
  /** Planck constant convenience name [J s] **/
  constexpr Numeric h = planck_constant;
  
  /** Reduced planck constant [J s] **/
  constexpr Numeric reduced_planck_constant = h * inv_two_pi;
  
  /** Reduced planck constant convenience name [J s] **/
  constexpr Numeric h_bar = reduced_planck_constant;
  
  /** Elementary charge [C]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units 2019-04-01
    **/
  constexpr Numeric elementary_charge = 1.602176634e-19;
  
  /** Elementary charge convenience name [C] **/
  constexpr Numeric e = elementary_charge;
  
  /** Boltzmann constant [J/K]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units 2019-04-01
    **/
  constexpr Numeric boltzmann_constant = 1.380649e-23;
  
  /** Boltzmann constant convenience name [J/K] **/
  constexpr Numeric k = boltzmann_constant;
  
  /** Avogadro constant [1/mol]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units 2019-04-01
    **/
  constexpr Numeric avogadro_constant = 6.02214076e23;
  
  /** Avogadro constant convenience name [1/mol] **/
  constexpr Numeric NA = avogadro_constant;
  
  /** Luminous efficacy of monochromatic 540 THz radiation [lm / W]
   From: https://en.wikipedia.org/wiki/2019_redefinition_of_SI_base_units 2019-04-01
    **/
  constexpr Numeric K_cd = 683;
  
  /** Fine structure constant [-]
   From: https://physics.nist.gov/cgi-bin/cuu/Value?alph 2019-04-01
   Reported error: (17)
   
   NOTE: Value from before update to new SI.  Double check in May.
    **/
  constexpr Numeric fine_structure_constant = 7.2973525664e-3;
  
  /** Fine structure constant convenience name [-] **/
  constexpr Numeric alpha = fine_structure_constant;
  
  /** Rydberg constant [1/m]
    From: https://physics.nist.gov/cgi-bin/cuu/Value?ryd 2019-04-01
    Reported error: (65)
    
    NOTE: Value from before update to new SI.  Double check in May.
    **/
  constexpr Numeric rydberg_constant = 10973731.568508;
  
  /** Rydberg constant convenience name [1/m] **/
  constexpr Numeric R_inf = rydberg_constant;
  
  /** Magnetic constant [H/m] **/
  constexpr Numeric magnetic_constant = 2 * h * alpha / (c * pow2(e));
  
  /** Magnetic constant convenience name [H/m] **/
  constexpr Numeric mu_0 = magnetic_constant;
  
  /** Vacuum permittivity [F/m] **/
  constexpr Numeric vacuum_permittivity = pow2(e) / (2 * h * c * alpha);
  
  /** Vacuum permittivity convenience name [F/m] **/
  constexpr Numeric epsilon_0 = vacuum_permittivity;
  
  /** Mass of resting electron [kg] **/
  constexpr Numeric electron_mass = 2 * h * R_inf / (c * pow2(alpha));
  
  /** Mass of resting electron convenience name [kg] **/
  constexpr Numeric m_e = electron_mass;
  
  /** Bohr magneton [J/T] **/
  constexpr Numeric bohr_magneton = e * h_bar / (2 * m_e);
  
  /** Ideal gas constant [J/mol K] **/
  constexpr Numeric ideal_gas_constant = k * NA;
  
  /** Ideal gas constant convenience name [J/mol K] **/
  constexpr Numeric R = ideal_gas_constant;
  
  /** Doppler broadening constant squared [kg/T]^2 **/
  constexpr Numeric doppler_broadening_const_squared = 2000 * R / pow2(c);
};


/** Namespace containing several practical unit conversions, physical and mathematical **/
namespace Conversion {
  using namespace Constant;
  
  /** Conversion constant degrees to radians and back.  Use conversion formulae instead of pure constant if possible. NOTE:  No constexpr cos etal in ARTS yet. **/
  constexpr Numeric DEG2RAD = pi/180;
  constexpr Numeric RAD2DEG = 1/DEG2RAD;
  template <class T> constexpr Numeric deg2rad(T x) {return x*DEG2RAD;}
  template <class T> constexpr Numeric rad2deg(T x) {return x*RAD2DEG;}
  template <class T> Numeric cosd(T x) {return std::cos(deg2rad(x));}
  template <class T> Numeric sind(T x) {return std::sin(deg2rad(x));}
  template <class T> Numeric tand(T x) {return std::tan(deg2rad(x));}
  template <class T> Numeric acosd(T x) {return rad2deg(std::acos(x));}
  template <class T> Numeric asind(T x) {return rad2deg(std::asin(x));}
  template <class T> Numeric atand(T x) {return rad2deg(std::atan(x));}
  template <class T1, class T2> Numeric atand2(T1 x, T2 y) {return rad2deg(std::atan2(x, y));}
  
  /** Conversion constant Kayser wavenumber to frequency and back.  Use conversion formulae instead of pure constant if possible. **/
  constexpr Numeric KAYCM2FREQ = 100*c;
  constexpr Numeric FREQ2KAYCM = 1/KAYCM2FREQ;
  template <class T> constexpr Numeric kaycm2freq(T x) {return x*KAYCM2FREQ;}
  template <class T> constexpr Numeric freq2kaycm(T x) {return x*FREQ2KAYCM;}
  
  /** Conversion constant Angular wavenumber to frequency and back.  Use conversion formulae instead of pure constant if possible. **/
  constexpr Numeric ANGCM2FREQ = KAYCM2FREQ * inv_two_pi;
  constexpr Numeric FREQ2ANGCM = 1/ANGCM2FREQ;
  template <class T> constexpr Numeric angcm2freq(T x) {return x*ANGCM2FREQ;}
  template <class T> constexpr Numeric freq2angcm(T x) {return x*FREQ2ANGCM;}
  
  /** Conversion constant Angular frequency to frequency and back **/
  template <class T> constexpr Numeric angfreq2freq(T x) {return x*inv_two_pi;}
  template <class T> constexpr Numeric freq2angfreq(T x) {return x*two_pi;}
  
  /** Conversion wavelength to frequency and back. **/
  template <class T> constexpr Numeric wavelen2freq(T x) {return c/x;}
  template <class T> constexpr Numeric freq2wavelen(T x) {return c/x;}
  
  /** Conversion constant 1 atmosphere to 1 Pascal and back.  Use conversion formulae instead of pure constant if possible. **/
  constexpr Numeric ATM2PA = 101325;
  constexpr Numeric PA2ATM = 1/ATM2PA;
  template <class T> constexpr Numeric atm2pa(T x) {return x*ATM2PA;}
  template <class T> constexpr Numeric pa2atm(T x) {return x*PA2ATM;}
  
  /** Conversion constant 1 torr to 1 Pascal and back.  Use conversion formulae instead of pure constant if possible. **/
  constexpr Numeric TORR2PA = ATM2PA/760;
  constexpr Numeric PA2TORR = 1/TORR2PA;
  template <class T> constexpr Numeric torr2pa(T x) {return x*TORR2PA;}
  template <class T> constexpr Numeric pa2torr(T x) {return x*PA2TORR;}
  
  /** Conversion constant Celsius to Kelvin and back.  Use conversion formulae instead of pure constant if possible. **/
  constexpr Numeric CEL2KEL = 273.15;
  template <class T> constexpr Numeric celsius2kelvin(T x) {return x + CEL2KEL;}
  template <class T> constexpr Numeric kelvin2celsius(T x) {return x - CEL2KEL;}
  
  /** Conversion from cm-1 per molecule per cm^2 to Hz per molecule per m^2 **/
  constexpr Numeric HITRAN2ARTS_LS = KAYCM2FREQ * 1e-4;
  constexpr Numeric ARTS2HITRAN_LS = 1 / HITRAN2ARTS_LS;
  template <class T> constexpr Numeric hitran2arts_linestrength(T x) {return x*HITRAN2ARTS_LS;}
  template <class T> constexpr Numeric arts2hitran_linestrength(T x) {return x*ARTS2HITRAN_LS;}
  
  /** Conversion from cm-1 per atmosphere to Hz per Pascal **/
  constexpr Numeric HITRAN2ARTS_GAMMA = KAYCM2FREQ / ATM2PA;
  constexpr Numeric ARTS2HITRAN_GAMMA = 1 / HITRAN2ARTS_GAMMA;
  template <class T> constexpr Numeric hitran2arts_broadening(T x) {return x*HITRAN2ARTS_GAMMA;}
  template <class T> constexpr Numeric arts2hitran_broadening(T x) {return x*ARTS2HITRAN_GAMMA;}
  
  /** Conversion from cm-1 to Joule **/
  constexpr Numeric HITRAN2ARTS_ENERGY = h * KAYCM2FREQ;
  constexpr Numeric ARTS2HITRAN_ENERGY = 1 / HITRAN2ARTS_ENERGY;
  template <class T> constexpr Numeric hitran2arts_energy(T x) {return x*HITRAN2ARTS_ENERGY;}
  template <class T> constexpr Numeric arts2hitran_energy(T x) {return x*ARTS2HITRAN_ENERGY;}
};

#endif
