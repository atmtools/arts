/**
   The ARTS running version number.
   Please always increase this before you do a `cvs
   commit', however minor your change may be.
 */

#include <string>
#include "config.h"

#define SUBVERSION "5"

string subversion = SUBVERSION;
string full_name  = static_cast<string>(PACKAGE) + "-" + VERSION + "." + SUBVERSION;

