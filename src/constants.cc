/* Copyright (C) 2000-2008
   Stefan Buehler <sbuehler@ltu.se>
   Patrick Eriksson <patrick@rss.chalmers.se>

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



////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
   \file   constants.cc

   This file contains global constants. You can use them anywhere 
   by declaring them as in the following example:

   extern const Numeric PI;

   See for example the
   <a href="http://physics.nist.gov/cuu/index.html">National Institute of Standards and Technology (NIST) </a>  
   home page for the values of specific constants.<br>

   \author Patrick Eriksson
   \date 2000-09-14 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "arts.h"
#include "matpackI.h"


////////////////////////////////////////////////////////////////////////////
// The constants
////////////////////////////////////////////////////////////////////////////

/** Global constant, the radius of the Earth [m]
    \author Patrick Eriksson 
    \date   08.04.2000
*/
extern const Numeric EARTH_RADIUS   = 6.378e6;

/** Global constant, conversion from radians to degrees

    Multiply your value in radians by this constant to get 
    the value in degrees.

    \author Patrick Eriksson 
    \date   08.04.2000
*/
extern const Numeric RAD2DEG        = 57.29577951308232;

/** Global constant, conversion from degrees to radians

    Multiply your value in degrees by this constant to get 
    the value in radians.

    \author Patrick Eriksson 
    \date   08.04.2000
*/
extern const Numeric DEG2RAD        = 0.01745329251994;

/** Global constant, the Planck constant [Js]
    \author Patrick Eriksson 
    \date   08.04.2000
*/
extern const Numeric PLANCK_CONST   = 6.626180e-34;

/** Global constant, spped of light in vaccum [m/s]
    \author Patrick Eriksson 
    \date   08.04.2000
*/
extern const Numeric SPEED_OF_LIGHT = 2.99792458e8;

/** Global constant, the Boltzmann constant [J/K]
    \author Patrick Eriksson 
    \date   08.04.2000
*/
extern const Numeric BOLTZMAN_CONST = 1.380662e-23;

/** Global constant, the Avogadro's number [molec/kg]
    \author Patrick Eriksson 
    \date   08.04.2000
*/
extern const Numeric AVOGADROS_NUMB = 6.0220450e26;

/** Atomic mass unit, 12th of a C^12_6 atom.
    \author Oliver Lemke
    \date   27.09.2000
*/
extern const Numeric ATOMIC_MASS_UNIT = 1.6606E-27;

/** Global constant, Planck temperature for cosmic background radiation [K]
    \author Patrick Eriksson 
    \date   08.04.2000
*/
extern const Numeric COSMIC_BG_TEMP = 2.735;

/** Global constant, Planck temperature for solar radiation [K]
    \author Patrick Eriksson 
    \date   08.04.2000
*/
extern const Numeric SUN_TEMP = 6000.0;

/** Global constant, e (Euler's number)
    \author Thomas Kuhn
    \date   08.11.2001
*/
// reference: http://thesapps.com/Doug/exp/
extern const Numeric EULER_NUMBER   = 2.7182818284590452;

/** Global constant, log10(Euler's number)
    \author Thomas Kuhn
    \date   08.11.2001
*/
extern const Numeric LOG10_EULER_NUMBER = 0.43429448190325176;

/** Global constant, ln(10)
    \author Thomas Kuhn
    \date   08.11.2001
*/
extern const Numeric NAT_LOG_TEN = 2.3025850929940459;

/** Global constant, ln(2)
    \author Patrick Eriksson 
    \date   08.04.2000
*/
extern const Numeric NAT_LOG_2      = 0.69314718055994;

/** Global constant, sqrt(ln(2))
    \author Axel von Engeln
    \date   25.09.2000
*/
extern const Numeric SQRT_NAT_LOG_2      = 0.832554611;

/** Global constant, pi
    \author Patrick Eriksson 
    \date   08.04.2000
*/
extern const Numeric PI             = 3.14159265358979;

/** Global constant, converts atm to Pa.

    Multiply your value in atm by this constant to get the value in Pa.

    \author Stefan Buehler 
    \date   2000-06-17
*/
extern const Numeric ATM2PA         = 1.01325e5;

/** Global constant, converts torr to Pa.

    Multiply your value in torr by this constant to get the value in Pa.

    \author Axel von Engeln
    \date   2000-10-31
*/
extern const Numeric TORR2PA        = 133.3227;

/** Global constant, Temperature in Celsius of 0 Kelvin.

    \author Axel von Engeln
    \date   2000-12-19
*/
extern const Numeric TEMP_0_C =  273.15;

/** Global constant, Standard pressure in Pa.

    \author Axel von Engeln
    \date   2000-12-19
*/
extern const Numeric PRES_STAND = 101300.25;

/** Global constant,  Loschmidt constant [m^-3].

    \author Axel von Engeln
    \date   2000-12-19
*/
extern const Numeric LOSCHMIDT_CONST = 2.686763E25;

/** Global constant,   Earth gravitational constant [m^3/s^2].

    \author Carlos Jimenez
    \date   2001-04-20
*/
extern const Numeric EARTH_GRAV_CONST = 3.98601E14;

/** Global constant, converts Hz to cm-1.

    Multiply your value in Hz by this constant to get the value in cm-1.

    \author Patrick Eriksson
    \date   2003-09-07
*/
extern const Numeric HZ2CM =  0.01 / SPEED_OF_LIGHT;

/** Global constant, Index of the frequency grid in GField1.

    \author Patrick Eriksson
    \date   2008-07-02
*/
extern const Index GFIELD1_F_GRID = 0;

/** Global constant, Index of the pressure grid in GField3.

    \author Oliver Lemke
    \date   2008-06-24
*/
extern const Index GFIELD3_P_GRID = 0;

/** Global constant, Index of the latitude grid in GField3.

    \author Oliver Lemke
    \date   2008-06-24
*/
extern const Index GFIELD3_LAT_GRID = 1;

/** Global constant, Index of the longitude grid in GField3.

    \author Oliver Lemke
    \date   2008-06-24
*/
extern const Index GFIELD3_LON_GRID = 2;

/** Global constant, Index of the field names in GField4.

    \author Oliver Lemke
    \date   2008-06-25
*/
extern const Index GFIELD4_FIELD_NAMES = 0;

/** Global constant, Index of incidence angles in GField4.

    \author Patrick Eriksson
    \date   2008-09-20
*/
extern const Index GFIELD4_IA_GRID = 0;

/** Global constant, Index of the pressure grid in GField4.

    \author Oliver Lemke
    \date   2008-06-25
*/
extern const Index GFIELD4_P_GRID = 1;

/** Global constant, Index of the frequency grid in GField4.

    \author Patrick Eriksson
    \date   2008-07-01
*/
extern const Index GFIELD4_F_GRID = 1;

/** Global constant, Index of the latitude grid in GField4.

    \author Oliver Lemke
    \date   2008-06-25
*/
extern const Index GFIELD4_LAT_GRID = 2;

/** Global constant, Index of the zenith angle grid in GField4.

    \author Patrick Eriksson
    \date   2008-07-01
*/
extern const Index GFIELD4_ZA_GRID = 2;

/** Global constant, Index of the longitude grid in GField4.

    \author Oliver Lemke
    \date   2008-06-25
*/
extern const Index GFIELD4_LON_GRID = 3;

/** Global constant, Index of the azimuth angle grid in GField4.

    \author Patrick Eriksson
    \date   2008-07-01
*/
extern const Index GFIELD4_AA_GRID = 3;

/** Define the global joker objekt. This is used by Matpack to specify joker ranges.  
    
    \author Stefan Buehler
    \date   2001-09-10
*/
Joker joker;

