/*-----------------------------------------------------------------------
FILE:      constants.cc

INCLUDES:  This file contains global constants. You can use them anywhere 
           by declaring them as in the following example:
             extern const Numeric PI;

CONSTANTS: EARTH_RADIUS
           RAD2DEG
           DEG2RAD
           PLANCK_CONST
           SPEED_OF_LIGHT
           BOLTZMAN_CONST
           AVOGADROS_NUMB
           COSMIC_BG_TEMP
           SUN_TEMP
           NAT_LOG_2
           PI
	   WAVENUM2Hz
	   ATM2HPA

HISTORY:   08.04.2000 Created by Patrick Eriksson.
           22.05.2000 Stefan Buehler: 
	   	      - Put # around units.
		      - Added ATM2HPA
-----------------------------------------------------------------------*/

#include "arts.h"


/** Global constant, the radius of the Earth [#m#]
    @author Patrick Eriksson 08.04.2000 */
extern const Numeric EARTH_RADIUS   = 6.378e6;

/** Global constant, conversion from radians to degrees
    @author Patrick Eriksson 08.04.2000 */
extern const Numeric RAD2DEG        = 57.29577951308232;

/** Global constant, conversion from degrees to radians
    @author Patrick Eriksson 08.04.2000 */
extern const Numeric DEG2RAD        = 0.01745329251994;

/** Global constant, the Planck constant [#Js#]
    @author Patrick Eriksson 08.04.2000 */
extern const Numeric PLANCK_CONST   = 6.626180e-34;

/** Global constant, spped of light in vaccum [#m/s#]
    @author Patrick Eriksson 08.04.2000 */
extern const Numeric SPEED_OF_LIGHT = 2.99792458e8;

/** Global constant, the Boltzmann constant [#J/K#]
    @author Patrick Eriksson 08.04.2000 */
extern const Numeric BOLTZMAN_CONST = 1.380662e-23;

/** Global constant, the Avogadro's number [#molec/kg#]
    @author Patrick Eriksson 08.04.2000 */
extern const Numeric AVOGADROS_NUMB = 6.0220450e26;

/** Global constant, Planck temperature for cosmic background radiation [#K#]
    @author Patrick Eriksson 08.04.2000 */
extern const Numeric COSMIC_BG_TEMP = 2.735;

/** Global constant, Planck temperature for solar radiation [#K#]
    @author Patrick Eriksson 12.04.2000 */
extern const Numeric SUN_TEMP = 6000.0;

/** Global constant, ln(2)
    @author Patrick Eriksson 08.04.2000 */
extern const Numeric NAT_LOG_2      = 0.69314718055994;

/** Global constant, pi
    @author Patrick Eriksson 08.04.2000 */
extern const Numeric PI             = 3.14159265358979;

/** Global constant, converts atm to hPa. Multiply the your value in
    atm by this constant to get the value in hPa.

    @author Stefan Buehler 22.05.2000  */
extern const Numeric ATM2HPA        = 1.01325e3;
