/* This file contains global constants. You can use them anywhere by
   declaring them as in the following example:
   extern Numeric PI; */

#include "arts.h"

// FIXME: Somebody should document the units here...

/// The radius of Earth in Meters.
Numeric EARTH_RADIUS   = 6.378e6;
Numeric RAD2DEG        = 57.29577951308232;
Numeric DEG2RAD        = 0.01745329251994;
Numeric PLANCK_CONST   = 6.626180e-34;
Numeric SPEED_OF_LIGHT = 2.99792458e8;
Numeric BOLTZMAN_CONST = 1.380662e-23;
Numeric AVOGADROS_NUMB = 6.0220450e26;
Numeric COSMIC_BG_TEMP = 2.735;
Numeric NAT_LOG_2      = 0.69314718055994;
Numeric PI             = 3.14159265358979 ;

// This does not work:
//Numeric INF             9.99e99                // [**better solution ??]

