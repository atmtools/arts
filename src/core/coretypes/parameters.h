/**
   \file   parameters.h

   This file contains header information for the dealing with command
   line parameters.

   \author Stefan Buehler
   \date   2001-07-24
*/

#ifndef parameters_h
#define parameters_h

#include <mystring.h>

struct Parameters {
  /** List of paths to search for include files. */
  ArrayOfString includepath{};
  /** List of paths to search for data files. */
  ArrayOfString datapath{};
};

void parse_path_from_environment(const String& envvar, ArrayOfString& paths);

#endif  // parameters_h
